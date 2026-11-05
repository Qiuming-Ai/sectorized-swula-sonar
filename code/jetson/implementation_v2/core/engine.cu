#include "engine.h"

#include <cuda_runtime.h>
#include <cufft.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace sonar_v2 {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr int kBeamStreams = 8;
constexpr int kThreads = 256;

void checkCuda(cudaError_t code, const char* operation) {
    if (code != cudaSuccess) {
        throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(code));
    }
}

void checkCufft(cufftResult code, const char* operation) {
    if (code != CUFFT_SUCCESS) {
        throw std::runtime_error(std::string(operation) + ": cuFFT code " +
                                 std::to_string(static_cast<int>(code)));
    }
}

__device__ __forceinline__ cufftComplex cmul(cufftComplex a, cufftComplex b) {
    return make_cuFloatComplex(a.x * b.x - a.y * b.y,
                              a.x * b.y + a.y * b.x);
}

__global__ void copyAndPadRecord(const cufftComplex* raw, cufftComplex* spectrum,
                                int channels, int rawLength, int fftLength) {
    const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    const std::size_t count = static_cast<std::size_t>(channels) * fftLength;
    if (index >= count) return;
    const int channel = static_cast<int>(index / fftLength);
    const int sample = static_cast<int>(index % fftLength);
    spectrum[index] = sample < rawLength ? raw[static_cast<std::size_t>(channel) * rawLength + sample]
                                        : make_cuFloatComplex(0.0f, 0.0f);
}

__global__ void applyFrontend(cufftComplex* spectrum, const cufftComplex* frontend,
                              std::size_t count, int fftLength) {
    const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index < count) spectrum[index] = cmul(spectrum[index], frontend[index % fftLength]);
}

// Complex even-length scipy.signal.resample convention: the aliased native
// Nyquist bin is divided equally between the positive and negative fine bins.
__global__ void constructFineMatchedSpectrum(const cufftComplex* frontendSpectrum,
                                             const cufftComplex* templateSpectrum,
                                             cufftComplex* fineSpectrum,
                                             int channels, int fftLength, int fineLength) {
    const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    const std::size_t count = static_cast<std::size_t>(channels) * fftLength;
    if (index >= count) return;
    const int channel = static_cast<int>(index / fftLength);
    const int bin = static_cast<int>(index % fftLength);
    cufftComplex q = templateSpectrum[bin];
    q.y = -q.y;
    cufftComplex value = cmul(frontendSpectrum[index], q);
    const std::size_t offset = static_cast<std::size_t>(channel) * fineLength;
    const int half = fftLength / 2;
    if (fineLength == fftLength) {
        fineSpectrum[offset + bin] = value;
    } else if (bin < half) {
        fineSpectrum[offset + bin] = value;
    } else if (bin == half) {
        value.x *= 0.5f;
        value.y *= 0.5f;
        fineSpectrum[offset + half] = value;
        fineSpectrum[offset + fineLength - half] = value;
    } else {
        fineSpectrum[offset + fineLength - fftLength + bin] = value;
    }
}

__global__ void gatherNativeMatchedFilter(const cufftComplex* fineTime,
                                          cufftComplex* nativeTime,
                                          std::size_t count, int fftLength,
                                          int fineLength, int oversampling) {
    const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= count) return;
    const std::size_t channel = index / fftLength;
    const std::size_t sample = index % fftLength;
    cufftComplex value = fineTime[channel * fineLength + sample * oversampling];
    const float normalization = 1.0f / static_cast<float>(fftLength);
    value.x *= normalization;
    value.y *= normalization;
    nativeTime[index] = value;
}

__global__ void sectorDas(const cufftComplex* fineTime, cufftComplex* output,
                          const double* rangeBaseIndex, const double* delayIndex,
                          const cufftComplex* weightCarrier,
                          const cufftComplex* rangeCarrier,
                          const int* angleIndices, int angleCount,
                          int channels, int rangeCount, int angleColumns,
                          int fineLength, int fftLength) {
    const int index = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (index >= angleCount * rangeCount) return;
    const int localAngle = index / rangeCount;
    const int rangeRow = index % rangeCount;
    const int angle = angleIndices[localAngle];
    float sumReal = 0.0f;
    float sumImag = 0.0f;
    for (int channel = 0; channel < channels; ++channel) {
        const std::size_t geometryIndex = static_cast<std::size_t>(angle) * channels + channel;
        const double sample = rangeBaseIndex[rangeRow] + delayIndex[geometryIndex];
        // Constructor rejects grids outside the guarded interpolation domain.
        // This safety branch is not a silent substitute for a valid fixture.
        if (sample < 0.0 || sample > static_cast<double>(fineLength - 1)) continue;
        const int lower = static_cast<int>(floor(sample));
        const int upper = lower + 1 < fineLength ? lower + 1 : lower;
        const float fraction = static_cast<float>(sample - static_cast<double>(lower));
        const std::size_t offset = static_cast<std::size_t>(channel) * fineLength;
        const cufftComplex left = fineTime[offset + lower];
        const cufftComplex right = fineTime[offset + upper];
        const cufftComplex interpolated = make_cuFloatComplex(
            left.x + fraction * (right.x - left.x),
            left.y + fraction * (right.y - left.y));
        const cufftComplex value = cmul(interpolated, weightCarrier[geometryIndex]);
        sumReal += value.x;
        sumImag += value.y;
    }
    // cuFFT inverse is unnormalized. Dividing by Nfft (NOT Nfine) matches
    // native IFFT followed by scipy.signal.resample's Nfine/Nfft scaling.
    const float normalization = 1.0f / static_cast<float>(fftLength);
    cufftComplex sum = make_cuFloatComplex(sumReal * normalization, sumImag * normalization);
    output[static_cast<std::size_t>(rangeRow) * angleColumns + angle] = cmul(sum, rangeCarrier[rangeRow]);
}

unsigned int blocks(std::size_t count) {
    return static_cast<unsigned int>((count + kThreads - 1) / kThreads);
}

void validateConfig(const Config& cfg) {
    if (cfg.M < 2 || cfg.Nraw < 1 || cfg.Nfft < cfg.Nraw || cfg.Nfft % 2 != 0 ||
        cfg.OS < 1 || cfg.Nr < 1 || cfg.Na < 1 || cfg.K < 1 ||
        !std::isfinite(cfg.fs) || !std::isfinite(cfg.fc) ||
        !std::isfinite(cfg.c) || !std::isfinite(cfg.gate) ||
        cfg.fs <= 0.0 || cfg.fc <= 0.0 || cfg.c <= 0.0) {
        throw std::invalid_argument("Invalid reconstruction dimensions or physical parameters");
    }
    if (cfg.Nfft > std::numeric_limits<int>::max() / cfg.OS) {
        throw std::invalid_argument("Oversampled FFT length exceeds cuFFT integer dimensions");
    }
    if (cfg.Na > std::numeric_limits<int>::max() / cfg.Nr) {
        throw std::invalid_argument("Image element count exceeds kernel integer dimensions");
    }
    if (cfg.positions.size() != static_cast<std::size_t>(cfg.M) ||
        cfg.angles.size() != static_cast<std::size_t>(cfg.Na) ||
        cfg.ranges.size() != static_cast<std::size_t>(cfg.Nr) ||
        cfg.sector_map.size() != static_cast<std::size_t>(cfg.Na) ||
        cfg.G.size() != static_cast<std::size_t>(cfg.Nfft) ||
        cfg.Q.size() != static_cast<std::size_t>(cfg.K) * cfg.Nfft) {
        throw std::invalid_argument("Config array sizes do not agree with dimensions");
    }
    for (const int sector : cfg.sector_map) {
        if (sector < 0 || sector >= cfg.K) throw std::invalid_argument("Out-of-range sector index");
    }
    for (const double angle : cfg.angles) {
        if (!std::isfinite(angle) || angle < -50.0 * kPi / 180.0 - 1e-12 ||
            angle > 50.0 * kPi / 180.0 + 1e-12) {
            throw std::invalid_argument("Angle is nonfinite or outside the declared +-50 degree sector");
        }
    }
    for (const double position : cfg.positions) {
        if (!std::isfinite(position)) throw std::invalid_argument("Nonfinite channel position");
    }
    for (const double range : cfg.ranges) {
        if (!std::isfinite(range)) throw std::invalid_argument("Nonfinite reconstruction range");
    }
    for (const auto value : cfg.G) {
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
            throw std::invalid_argument("Nonfinite receive frontend spectrum");
        }
    }
    for (const auto value : cfg.Q) {
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
            throw std::invalid_argument("Nonfinite effective matched-filter template spectrum");
        }
    }
}

} // namespace

struct Engine::Impl {
    Config config;
    bool cached;
    bool diagnostic;
    int fineLength;
    MemoryReport memory;
    Diagnostics stage;
    cudaStream_t computeStream = nullptr;
    std::array<cudaStream_t, kBeamStreams> beamStreams{};
    cudaEvent_t inverseReady = nullptr;
    std::array<cudaEvent_t, kBeamStreams> beamComplete{};
    cufftHandle forwardPlan = 0;
    cufftHandle inversePlan = 0;
    cufftComplex* rawDevice = nullptr;
    cufftComplex* frontendDevice = nullptr;
    cufftComplex* fineDevice = nullptr;
    cufftComplex* outputDevice = nullptr;
    cufftComplex* frontendTransfer = nullptr;
    cufftComplex* templateTransfer = nullptr;
    cufftComplex* weightCarrierDevice = nullptr;
    cufftComplex* rangeCarrierDevice = nullptr;
    double* rangeIndexDevice = nullptr;
    double* delayIndexDevice = nullptr;
    int* angleIndicesDevice = nullptr;
    cufftComplex* nativeDiagnosticDevice = nullptr;
    void* workspaceDevice = nullptr;
    std::vector<std::pair<void*, std::size_t>> allocations;
    std::vector<int> sectorOffsets;

    template<class T>
    void allocate(T*& pointer, std::size_t count) {
        const std::size_t bytes = count * sizeof(T);
        if (bytes == 0) return;
        checkCuda(cudaMalloc(reinterpret_cast<void**>(&pointer), bytes), "cudaMalloc");
        allocations.emplace_back(pointer, bytes);
        memory.explicit_live_bytes += bytes;
        memory.explicit_peak_bytes = std::max(memory.explicit_peak_bytes, memory.explicit_live_bytes);
        ++memory.allocation_count;
    }

    void allocateWorkspace(std::size_t bytes) {
        if (bytes == 0) return;
        checkCuda(cudaMalloc(&workspaceDevice, bytes), "cudaMalloc shared FFT workspace");
        allocations.emplace_back(workspaceDevice, bytes);
        memory.explicit_live_bytes += bytes;
        memory.explicit_peak_bytes = std::max(memory.explicit_peak_bytes, memory.explicit_live_bytes);
        ++memory.allocation_count;
    }

    template<class T>
    void upload(T* destination, const std::vector<T>& source) {
        checkCuda(cudaMemcpyAsync(destination, source.data(), source.size() * sizeof(T),
                                  cudaMemcpyHostToDevice, computeStream), "static H2D upload");
    }

    Impl(const Config& cfg, bool useCache, bool collectDiagnostics)
        : config(cfg), cached(useCache), diagnostic(collectDiagnostics), fineLength(0) {
        validateConfig(config);
        fineLength = config.Nfft * config.OS;
        static_assert(sizeof(std::complex<float>) == sizeof(cufftComplex),
                      "This platform must use two adjacent floats for complex<float>");
        try {
            checkCuda(cudaMemGetInfo(&memory.cuda_free_bytes_before, &memory.cuda_total_bytes),
                      "cudaMemGetInfo before setup");
            checkCuda(cudaStreamCreateWithFlags(&computeStream, cudaStreamNonBlocking),
                      "create compute stream");
            for (int stream = 0; stream < kBeamStreams; ++stream) {
                checkCuda(cudaStreamCreateWithFlags(&beamStreams[stream], cudaStreamNonBlocking),
                          "create DAS stream");
                checkCuda(cudaEventCreateWithFlags(&beamComplete[stream], cudaEventDisableTiming),
                          "create DAS completion event");
            }
            checkCuda(cudaEventCreateWithFlags(&inverseReady, cudaEventDisableTiming),
                      "create inverse-ready event");

            checkCufft(cufftCreate(&forwardPlan), "create forward plan");
            checkCufft(cufftSetAutoAllocation(forwardPlan, 0), "disable forward auto-allocation");
            int forwardLength[1] = {config.Nfft};
            checkCufft(cufftMakePlanMany(forwardPlan, 1, forwardLength, nullptr, 1, config.Nfft,
                                       nullptr, 1, config.Nfft, CUFFT_C2C, config.M,
                                       &memory.forward_plan_workspace_bytes), "plan batched forward FFT");
            checkCufft(cufftGetSize(forwardPlan, &memory.forward_plan_workspace_bytes),
                       "query forward FFT workspace");
            checkCufft(cufftSetStream(forwardPlan, computeStream), "set forward stream");
            checkCufft(cufftCreate(&inversePlan), "create inverse plan");
            checkCufft(cufftSetAutoAllocation(inversePlan, 0), "disable inverse auto-allocation");
            int inverseLength[1] = {fineLength};
            checkCufft(cufftMakePlanMany(inversePlan, 1, inverseLength, nullptr, 1, fineLength,
                                       nullptr, 1, fineLength, CUFFT_C2C, config.M,
                                       &memory.inverse_plan_workspace_bytes), "plan batched inverse FFT");
            checkCufft(cufftGetSize(inversePlan, &memory.inverse_plan_workspace_bytes),
                       "query inverse FFT workspace");
            checkCufft(cufftSetStream(inversePlan, computeStream), "set inverse stream");
            memory.shared_plan_workspace_bytes = std::max(memory.forward_plan_workspace_bytes,
                                                          memory.inverse_plan_workspace_bytes);

            allocate(rawDevice, static_cast<std::size_t>(config.M) * config.Nraw);
            allocate(frontendDevice, static_cast<std::size_t>(config.M) * config.Nfft);
            allocate(fineDevice, static_cast<std::size_t>(config.M) * fineLength);
            allocate(outputDevice, static_cast<std::size_t>(config.Nr) * config.Na);
            allocate(frontendTransfer, static_cast<std::size_t>(config.Nfft));
            allocate(templateTransfer, static_cast<std::size_t>(config.K) * config.Nfft);
            allocate(weightCarrierDevice, static_cast<std::size_t>(config.Na) * config.M);
            allocate(rangeCarrierDevice, static_cast<std::size_t>(config.Nr));
            allocate(rangeIndexDevice, static_cast<std::size_t>(config.Nr));
            allocate(delayIndexDevice, static_cast<std::size_t>(config.Na) * config.M);
            allocate(angleIndicesDevice, static_cast<std::size_t>(config.Na));
            if (diagnostic) allocate(nativeDiagnosticDevice, static_cast<std::size_t>(config.M) * config.Nfft);
            memory.non_workspace_bytes = memory.explicit_live_bytes;
            allocateWorkspace(memory.shared_plan_workspace_bytes);
            checkCufft(cufftSetWorkArea(forwardPlan, workspaceDevice), "set forward workspace");
            checkCufft(cufftSetWorkArea(inversePlan, workspaceDevice), "set inverse workspace");

            std::vector<cufftComplex> frontend(config.G.size());
            std::vector<cufftComplex> templates(config.Q.size());
            for (std::size_t i = 0; i < frontend.size(); ++i) frontend[i] = make_cuFloatComplex(config.G[i].real(), config.G[i].imag());
            for (std::size_t i = 0; i < templates.size(); ++i) templates[i] = make_cuFloatComplex(config.Q[i].real(), config.Q[i].imag());
            std::vector<double> rangeIndices(config.Nr);
            std::vector<double> delayIndices(static_cast<std::size_t>(config.Na) * config.M);
            std::vector<cufftComplex> rangeCarriers(config.Nr);
            std::vector<cufftComplex> weightCarriers(delayIndices.size());
            std::vector<double> weights(config.M);
            double weightSum = 0.0;
            for (int m = 0; m < config.M; ++m) {
                weights[m] = 0.54 - 0.46 * std::cos(2.0 * kPi * m / (config.M - 1));
                weightSum += weights[m];
            }
            for (int r = 0; r < config.Nr; ++r) {
                const double delay = 2.0 * (config.ranges[r] - config.gate) / config.c;
                rangeIndices[r] = delay * config.fs * config.OS;
                const double phase = 2.0 * kPi * config.fc * delay;
                rangeCarriers[r] = make_cuFloatComplex(static_cast<float>(std::cos(phase)),
                                                      static_cast<float>(std::sin(phase)));
            }
            for (int a = 0; a < config.Na; ++a) {
                const double sine = std::sin(config.angles[a]);
                for (int m = 0; m < config.M; ++m) {
                    const std::size_t index = static_cast<std::size_t>(a) * config.M + m;
                    const double delay = -config.positions[m] * sine / config.c;
                    delayIndices[index] = delay * config.fs * config.OS;
                    const double phase = 2.0 * kPi * config.fc * delay;
                    const double weight = weights[m] / weightSum;
                    weightCarriers[index] = make_cuFloatComplex(static_cast<float>(weight * std::cos(phase)),
                                                               static_cast<float>(weight * std::sin(phase)));
                }
            }
            const auto rangeBounds = std::minmax_element(rangeIndices.begin(), rangeIndices.end());
            for (const double delayIndex : delayIndices) {
                if (*rangeBounds.first + delayIndex < 0.0 ||
                    *rangeBounds.second + delayIndex >= static_cast<double>(fineLength - 1)) {
                    throw std::invalid_argument("Image grid lies outside guarded linear interpolation record");
                }
            }
            std::vector<int> angleIndices;
            angleIndices.reserve(config.Na);
            sectorOffsets.push_back(0);
            for (int k = 0; k < config.K; ++k) {
                for (int a = 0; a < config.Na; ++a) {
                    if (config.sector_map[a] == k) angleIndices.push_back(a);
                }
                sectorOffsets.push_back(static_cast<int>(angleIndices.size()));
            }
            upload(frontendTransfer, frontend);
            upload(templateTransfer, templates);
            upload(rangeIndexDevice, rangeIndices);
            upload(delayIndexDevice, delayIndices);
            upload(rangeCarrierDevice, rangeCarriers);
            upload(weightCarrierDevice, weightCarriers);
            upload(angleIndicesDevice, angleIndices);
            checkCuda(cudaStreamSynchronize(computeStream), "finish static setup uploads");
            checkCuda(cudaMemGetInfo(&memory.cuda_free_bytes_after_setup, &memory.cuda_total_bytes),
                      "cudaMemGetInfo after setup");
            if (diagnostic) {
                stage.frontend_spectrum.resize(static_cast<std::size_t>(config.M) * config.Nfft);
                stage.matched_filter_native.resize(static_cast<std::size_t>(config.K) * config.M * config.Nfft);
            }
        } catch (...) {
            cleanup();
            throw;
        }
    }

    ~Impl() { cleanup(); }

    void cleanup() noexcept {
        if (computeStream) cudaStreamSynchronize(computeStream);
        for (const auto stream : beamStreams) if (stream) cudaStreamSynchronize(stream);
        if (forwardPlan) cufftDestroy(forwardPlan);
        if (inversePlan) cufftDestroy(inversePlan);
        forwardPlan = inversePlan = 0;
        for (auto iterator = allocations.rbegin(); iterator != allocations.rend(); ++iterator) {
            cudaFree(iterator->first);
            memory.explicit_live_bytes -= iterator->second;
        }
        allocations.clear();
        if (inverseReady) cudaEventDestroy(inverseReady);
        inverseReady = nullptr;
        for (int i = 0; i < kBeamStreams; ++i) {
            if (beamComplete[i]) cudaEventDestroy(beamComplete[i]);
            if (beamStreams[i]) cudaStreamDestroy(beamStreams[i]);
            beamComplete[i] = nullptr;
            beamStreams[i] = nullptr;
        }
        if (computeStream) cudaStreamDestroy(computeStream);
        computeStream = nullptr;
    }

    void run(const std::vector<std::complex<float>>& raw,
             std::vector<std::complex<float>>& output) {
        const std::size_t rawCount = static_cast<std::size_t>(config.M) * config.Nraw;
        const std::size_t nativeCount = static_cast<std::size_t>(config.M) * config.Nfft;
        const std::size_t fineCount = static_cast<std::size_t>(config.M) * fineLength;
        const std::size_t outputCount = static_cast<std::size_t>(config.Nr) * config.Na;
        if (raw.size() != rawCount) throw std::invalid_argument("Raw frame size mismatch");
        output.resize(outputCount);
        checkCuda(cudaMemcpyAsync(rawDevice, raw.data(), rawCount * sizeof(cufftComplex),
                                  cudaMemcpyHostToDevice, computeStream), "fresh raw-frame H2D");
        copyAndPadRecord<<<blocks(nativeCount), kThreads, 0, computeStream>>>(rawDevice, frontendDevice,
                                                                            config.M, config.Nraw, config.Nfft);
        checkCuda(cudaGetLastError(), "launch record padding");
        checkCufft(cufftExecC2C(forwardPlan, frontendDevice, frontendDevice, CUFFT_FORWARD), "execute forward FFT");
        applyFrontend<<<blocks(nativeCount), kThreads, 0, computeStream>>>(frontendDevice, frontendTransfer,
                                                                         nativeCount, config.Nfft);
        checkCuda(cudaGetLastError(), "launch one common receive frontend");
        if (diagnostic) {
            checkCuda(cudaMemcpyAsync(stage.frontend_spectrum.data(), frontendDevice,
                                      nativeCount * sizeof(cufftComplex), cudaMemcpyDeviceToHost,
                                      computeStream), "diagnostic frontend-spectrum D2H");
        }

        for (int k = 0; k < config.K; ++k) {
            const int begin = sectorOffsets[k];
            const int end = sectorOffsets[k + 1];
            const int angleCount = end - begin;
            if (angleCount == 0 && !diagnostic) continue;
            checkCuda(cudaMemsetAsync(fineDevice, 0, fineCount * sizeof(cufftComplex), computeStream),
                      "clear oversampled MF spectrum");
            constructFineMatchedSpectrum<<<blocks(nativeCount), kThreads, 0, computeStream>>>(
                frontendDevice, templateTransfer + static_cast<std::size_t>(k) * config.Nfft,
                fineDevice, config.M, config.Nfft, fineLength);
            checkCuda(cudaGetLastError(), "launch MF spectrum and complex Nyquist split");
            checkCufft(cufftExecC2C(inversePlan, fineDevice, fineDevice, CUFFT_INVERSE), "execute oversampled inverse FFT");
            if (diagnostic) {
                gatherNativeMatchedFilter<<<blocks(nativeCount), kThreads, 0, computeStream>>>(
                    fineDevice, nativeDiagnosticDevice, nativeCount, config.Nfft, fineLength, config.OS);
                checkCuda(cudaGetLastError(), "launch native MF diagnostic decimation");
                checkCuda(cudaMemcpyAsync(stage.matched_filter_native.data() + static_cast<std::size_t>(k) * nativeCount,
                                          nativeDiagnosticDevice, nativeCount * sizeof(cufftComplex),
                                          cudaMemcpyDeviceToHost, computeStream), "diagnostic MF-native D2H");
            }
            checkCuda(cudaEventRecord(inverseReady, computeStream), "record inverse-ready event");
            std::array<bool, kBeamStreams> used{};
            if (cached) {
                const int chunk = (angleCount + kBeamStreams - 1) / kBeamStreams;
                if (chunk > 0) {
                    for (int stream = 0; stream < kBeamStreams; ++stream) {
                        const int localBegin = stream * chunk;
                        if (localBegin >= angleCount) break;
                        const int localCount = std::min(chunk, angleCount - localBegin);
                        checkCuda(cudaStreamWaitEvent(beamStreams[stream], inverseReady, 0), "DAS wait for inverse FFT");
                        sectorDas<<<blocks(static_cast<std::size_t>(localCount) * config.Nr), kThreads, 0, beamStreams[stream]>>>(
                            fineDevice, outputDevice, rangeIndexDevice, delayIndexDevice,
                            weightCarrierDevice, rangeCarrierDevice, angleIndicesDevice + begin + localBegin,
                            localCount, config.M, config.Nr, config.Na, fineLength, config.Nfft);
                        checkCuda(cudaGetLastError(), "launch cached sector-angle block DAS");
                        used[stream] = true;
                    }
                }
            } else {
                for (int localAngle = 0; localAngle < angleCount; ++localAngle) {
                    const int stream = localAngle % kBeamStreams;
                    if (!used[stream]) {
                        checkCuda(cudaStreamWaitEvent(beamStreams[stream], inverseReady, 0), "DAS wait for inverse FFT");
                        used[stream] = true;
                    }
                    sectorDas<<<blocks(static_cast<std::size_t>(config.Nr)), kThreads, 0, beamStreams[stream]>>>(
                        fineDevice, outputDevice, rangeIndexDevice, delayIndexDevice,
                        weightCarrierDevice, rangeCarrierDevice, angleIndicesDevice + begin + localAngle,
                        1, config.M, config.Nr, config.Na, fineLength, config.Nfft);
                    checkCuda(cudaGetLastError(), "launch uncached single-angle DAS");
                }
            }
            // The fine buffer is shared between waveforms. Each next sector's
            // memset/IFFT and final output D2H must wait for every DAS consumer.
            for (int stream = 0; stream < kBeamStreams; ++stream) {
                if (!used[stream]) continue;
                checkCuda(cudaEventRecord(beamComplete[stream], beamStreams[stream]), "record DAS completion");
                checkCuda(cudaStreamWaitEvent(computeStream, beamComplete[stream], 0), "protect shared fine buffer");
            }
        }
        checkCuda(cudaMemcpyAsync(output.data(), outputDevice, outputCount * sizeof(cufftComplex),
                                  cudaMemcpyDeviceToHost, computeStream), "complex reconstruction D2H");
        checkCuda(cudaStreamSynchronize(computeStream), "complete reconstruction and output transfer");
    }
};

Engine::Engine(const Config& config, bool cached, bool diagnostic)
    : impl_(std::make_unique<Impl>(config, cached, diagnostic)) {}
Engine::~Engine() = default;
Engine::Engine(Engine&&) noexcept = default;
Engine& Engine::operator=(Engine&&) noexcept = default;

void Engine::process(const std::vector<std::complex<float>>& raw,
                     std::vector<std::complex<float>>& output) {
    if (!impl_) throw std::logic_error("Cannot process with a moved-from Engine");
    impl_->run(raw, output);
}

const Diagnostics& Engine::diagnostics() const {
    if (!impl_ || !impl_->diagnostic) throw std::logic_error("Diagnostics were not enabled for this Engine");
    return impl_->stage;
}

MemoryReport Engine::memoryReport() const {
    if (!impl_) throw std::logic_error("Cannot query a moved-from Engine");
    return impl_->memory;
}

std::string Engine::implementationName() const {
    if (!impl_) throw std::logic_error("Cannot query a moved-from Engine");
    return impl_->cached ? "persistent_plans_buffers_batched_sector_DAS" : "per_frame_plans_buffers_single_angle_DAS";
}

} // namespace sonar_v2
