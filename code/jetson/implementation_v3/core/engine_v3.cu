#include "engine_v3.h"

#include <cuda_runtime.h>
#include <cufft.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace sonar_v3 {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr int kBeamStreams = 8;
constexpr int kThreads = 256;

void cudaCheck(cudaError_t code, const char* operation) {
    if (code != cudaSuccess) {
        throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(code));
    }
}

void fftCheck(cufftResult code, const char* operation) {
    if (code != CUFFT_SUCCESS) {
        throw std::runtime_error(std::string(operation) + ": cuFFT code " +
                                 std::to_string(static_cast<int>(code)));
    }
}

__device__ __forceinline__ cufftComplex multiply(cufftComplex a, cufftComplex b) {
    return make_cuFloatComplex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

__global__ void padRaw(const cufftComplex* raw, cufftComplex* spectrum,
                       int channels, int rawLength, int fftLength) {
    const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= static_cast<std::size_t>(channels) * fftLength) return;
    const int channel = static_cast<int>(index / fftLength);
    const int sample = static_cast<int>(index % fftLength);
    spectrum[index] = sample < rawLength ? raw[static_cast<std::size_t>(channel) * rawLength + sample]
                                        : make_cuFloatComplex(0.0f, 0.0f);
}

__global__ void frontend(cufftComplex* spectrum, const cufftComplex* response,
                         std::size_t count, int fftLength) {
    const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index < count) spectrum[index] = multiply(spectrum[index], response[index % fftLength]);
}

__global__ void constructMf(const cufftComplex* spectrum, const cufftComplex* q,
                            cufftComplex* matched, std::size_t count, int fftLength) {
    const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= count) return;
    cufftComplex value = q[index % fftLength];
    value.y = -value.y;
    matched[index] = multiply(spectrum[index], value);
}

__global__ void normalize(cufftComplex* data, std::size_t count, float scale) {
    const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index < count) {
        data[index].x *= scale;
        data[index].y *= scale;
    }
}

__global__ void mixDown(cufftComplex* data, const cufftComplex* mixing,
                        std::size_t count, int fftLength) {
    const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index < count) data[index] = multiply(data[index], mixing[index % fftLength]);
}

// For native IFFT-normalized time x[n], decimated x[D*n] has low-rate DFT
// (1/D) sum_{d=0}^{D-1} X[k+d*Nlow]. The AA response is multiplied before
// alias folding. The ensuing low-rate inverse must be normalized by Nlow.
__global__ void filterFold(const cufftComplex* spectrum, const cufftComplex* response,
                           cufftComplex* low, int channels, int nativeLength,
                           int lowLength, int decimation) {
    const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= static_cast<std::size_t>(channels) * lowLength) return;
    const int channel = static_cast<int>(index / lowLength);
    const int bin = static_cast<int>(index % lowLength);
    cufftComplex sum = make_cuFloatComplex(0.0f, 0.0f);
    for (int alias = 0; alias < decimation; ++alias) {
        const int nativeBin = bin + alias * lowLength;
        const cufftComplex value = multiply(
            spectrum[static_cast<std::size_t>(channel) * nativeLength + nativeBin], response[nativeBin]);
        sum.x += value.x;
        sum.y += value.y;
    }
    const float scale = 1.0f / static_cast<float>(decimation);
    low[index] = make_cuFloatComplex(sum.x * scale, sum.y * scale);
}

__device__ __forceinline__ cufftComplex interpolateLinear(const cufftComplex* time,
                                                          std::size_t offset,
                                                          int lower, float fraction) {
    const cufftComplex left = time[offset + lower];
    const cufftComplex right = time[offset + lower + 1];
    return make_cuFloatComplex(left.x + fraction * (right.x - left.x),
                              left.y + fraction * (right.y - left.y));
}

__device__ __forceinline__ cufftComplex interpolateSinc(const cufftComplex* time,
                                                        std::size_t offset,
                                                        int lower, float fraction,
                                                        const float* lut, int taps,
                                                        int phases) {
    const float phase = fraction * phases;
    int phaseIndex = static_cast<int>(floorf(phase));
    phaseIndex = phaseIndex < phases ? phaseIndex : phases - 1;
    const float phaseFraction = phase - phaseIndex;
    const int first = lower - taps / 2 + 1;
    const std::size_t row = static_cast<std::size_t>(phaseIndex) * taps;
    cufftComplex sum = make_cuFloatComplex(0.0f, 0.0f);
    for (int tap = 0; tap < taps; ++tap) {
        const float left = lut[row + tap];
        const float right = lut[row + taps + tap];
        const float weight = left + phaseFraction * (right - left);
        const cufftComplex sample = time[offset + first + tap];
        sum.x += weight * sample.x;
        sum.y += weight * sample.y;
    }
    return sum;
}

__global__ void lowRateDas(const cufftComplex* low, cufftComplex* output,
                           const double* rangeIndex, const double* delayIndex,
                           const cufftComplex* weightCarrier,
                           const cufftComplex* rangeCarrier,
                           const int* angleIndices, int angleCount,
                           int channels, int rangeCount, int angleColumns,
                           int lowLength, const float* lut, int taps, int phases) {
    const int index = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (index >= angleCount * rangeCount) return;
    const int localAngle = index / rangeCount;
    const int rangeRow = index % rangeCount;
    const int angle = angleIndices[localAngle];
    cufftComplex sum = make_cuFloatComplex(0.0f, 0.0f);
    for (int channel = 0; channel < channels; ++channel) {
        const std::size_t geometry = static_cast<std::size_t>(angle) * channels + channel;
        const double sample = rangeIndex[rangeRow] + delayIndex[geometry];
        const int lower = static_cast<int>(floor(sample));
        const float fraction = static_cast<float>(sample - static_cast<double>(lower));
        const std::size_t offset = static_cast<std::size_t>(channel) * lowLength;
        // Configuration is validated against both the linear and sinc record
        // guard; no extrapolation or wraparound is silently used in DAS.
        const cufftComplex value = taps == 0 ?
            interpolateLinear(low, offset, lower, fraction) :
            interpolateSinc(low, offset, lower, fraction, lut, taps, phases);
        const cufftComplex weighted = multiply(value, weightCarrier[geometry]);
        sum.x += weighted.x;
        sum.y += weighted.y;
    }
    output[static_cast<std::size_t>(rangeRow) * angleColumns + angle] =
        multiply(sum, rangeCarrier[rangeRow]);
}

unsigned int blocks(std::size_t count) {
    return static_cast<unsigned int>((count + kThreads - 1) / kThreads);
}

void validate(const Config& cfg) {
    if (cfg.M < 2 || cfg.Nraw < 1 || cfg.Nfft < cfg.Nraw || cfg.Nfft % 2 != 0 ||
        cfg.Nr < 1 || cfg.Na < 1 || cfg.K < 1 || cfg.D < 1 || cfg.Nfft % cfg.D != 0 ||
        cfg.Nfft / cfg.D < 2 || cfg.interpolation_taps < 0 || cfg.interpolation_taps % 2 != 0 ||
        cfg.lut_phases < 1 || !std::isfinite(cfg.fs) || !std::isfinite(cfg.fc) ||
        !std::isfinite(cfg.c) || !std::isfinite(cfg.gate) ||
        cfg.fs <= 0.0 || cfg.fc <= 0.0 || cfg.c <= 0.0 ||
        cfg.Na > std::numeric_limits<int>::max() / cfg.Nr) {
        throw std::invalid_argument("Invalid v3 dimensions or physical parameters");
    }
    if (cfg.positions.size() != static_cast<std::size_t>(cfg.M) ||
        cfg.angles.size() != static_cast<std::size_t>(cfg.Na) ||
        cfg.ranges.size() != static_cast<std::size_t>(cfg.Nr) ||
        cfg.sector_map.size() != static_cast<std::size_t>(cfg.Na) ||
        cfg.G.size() != static_cast<std::size_t>(cfg.Nfft) ||
        cfg.Q.size() != static_cast<std::size_t>(cfg.K) * cfg.Nfft ||
        cfg.H.size() != static_cast<std::size_t>(cfg.Nfft) ||
        cfg.mix.size() != static_cast<std::size_t>(cfg.K) * cfg.Nfft ||
        cfg.offsets.size() != static_cast<std::size_t>(cfg.K) ||
        cfg.interpolation_lut.size() != static_cast<std::size_t>(cfg.lut_phases + 1) * cfg.interpolation_taps) {
        throw std::invalid_argument("v3 config array lengths do not match dimensions");
    }
    for (int sector : cfg.sector_map) {
        if (sector < 0 || sector >= cfg.K) throw std::invalid_argument("Invalid sector assignment");
    }
    for (double angle : cfg.angles) {
        if (!std::isfinite(angle) || std::abs(angle) > 50.0 * kPi / 180.0 + 1e-12)
            throw std::invalid_argument("Nonfinite or outside-sector reconstruction angle");
    }
    for (double position : cfg.positions) if (!std::isfinite(position)) throw std::invalid_argument("Nonfinite position");
    for (double range : cfg.ranges) if (!std::isfinite(range)) throw std::invalid_argument("Nonfinite range");
    for (double offset : cfg.offsets) if (!std::isfinite(offset)) throw std::invalid_argument("Nonfinite offset");
    for (const auto* values : {&cfg.G, &cfg.Q, &cfg.H, &cfg.mix}) {
        for (const auto value : *values) {
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
                throw std::invalid_argument("Nonfinite static response");
        }
    }
    for (float weight : cfg.interpolation_lut) if (!std::isfinite(weight)) throw std::invalid_argument("Nonfinite sinc LUT");
}

} // namespace

struct Engine::Impl {
    Config config;
    bool diagnostic;
    int lowLength;
    MemoryReport memory;
    Diagnostics stage;
    cudaStream_t computeStream = nullptr;
    std::array<cudaStream_t, kBeamStreams> beamStreams{};
    cudaEvent_t lowReady = nullptr;
    std::array<cudaEvent_t, kBeamStreams> beamComplete{};
    cufftHandle forwardPlan = 0;
    cufftHandle nativeInversePlan = 0;
    cufftHandle lowInversePlan = 0;
    cufftComplex* rawDevice = nullptr;
    cufftComplex* frontendDevice = nullptr;
    cufftComplex* matchedDevice = nullptr;
    cufftComplex* lowDevice = nullptr;
    cufftComplex* outputDevice = nullptr;
    cufftComplex* frontendTransfer = nullptr;
    cufftComplex* templateTransfer = nullptr;
    cufftComplex* filterTransfer = nullptr;
    cufftComplex* mixingTransfer = nullptr;
    cufftComplex* weightCarrierDevice = nullptr;
    cufftComplex* rangeCarrierDevice = nullptr;
    double* rangeIndexDevice = nullptr;
    double* delayIndexDevice = nullptr;
    int* angleIndicesDevice = nullptr;
    float* interpolationDevice = nullptr;
    void* workspaceDevice = nullptr;
    std::vector<std::pair<void*, std::size_t>> allocations;
    std::vector<int> sectorOffsets;

    template<class T>
    void allocate(T*& pointer, std::size_t count) {
        const std::size_t bytes = count * sizeof(T);
        if (bytes == 0) return;
        cudaCheck(cudaMalloc(reinterpret_cast<void**>(&pointer), bytes), "cudaMalloc");
        allocations.emplace_back(pointer, bytes);
        memory.explicit_live_bytes += bytes;
        memory.explicit_peak_bytes = std::max(memory.explicit_peak_bytes, memory.explicit_live_bytes);
        ++memory.allocation_count;
    }

    template<class T>
    void upload(T* destination, const std::vector<T>& values) {
        if (values.empty()) return;
        cudaCheck(cudaMemcpyAsync(destination, values.data(), values.size() * sizeof(T),
                                  cudaMemcpyHostToDevice, computeStream), "static configuration upload");
    }

    void createPlan(cufftHandle& handle, int length, std::size_t& workspace, const char* operation) {
        fftCheck(cufftCreate(&handle), "create FFT plan");
        fftCheck(cufftSetAutoAllocation(handle, 0), "disable implicit FFT allocation");
        int dimensions[1] = {length};
        fftCheck(cufftMakePlanMany(handle, 1, dimensions, nullptr, 1, length,
                                  nullptr, 1, length, CUFFT_C2C, config.M, &workspace), operation);
        fftCheck(cufftGetSize(handle, &workspace), "query FFT workspace");
        fftCheck(cufftSetStream(handle, computeStream), "set FFT compute stream");
    }

    Impl(const Config& cfg, bool collectDiagnostics)
        : config(cfg), diagnostic(collectDiagnostics), lowLength(0) {
        validate(config);
        lowLength = config.Nfft / config.D;
        static_assert(sizeof(std::complex<float>) == sizeof(cufftComplex), "Complex layout mismatch");
        try {
            cudaCheck(cudaMemGetInfo(&memory.cuda_free_bytes_before, &memory.cuda_total_bytes), "memory before setup");
            cudaCheck(cudaStreamCreateWithFlags(&computeStream, cudaStreamNonBlocking), "create compute stream");
            for (int stream = 0; stream < kBeamStreams; ++stream) {
                cudaCheck(cudaStreamCreateWithFlags(&beamStreams[stream], cudaStreamNonBlocking), "create DAS stream");
                cudaCheck(cudaEventCreateWithFlags(&beamComplete[stream], cudaEventDisableTiming), "create DAS event");
            }
            cudaCheck(cudaEventCreateWithFlags(&lowReady, cudaEventDisableTiming), "create decimation event");
            createPlan(forwardPlan, config.Nfft, memory.forward_plan_workspace_bytes, "plan native forward FFT");
            createPlan(nativeInversePlan, config.Nfft, memory.native_inverse_plan_workspace_bytes, "plan native MF inverse FFT");
            createPlan(lowInversePlan, lowLength, memory.low_inverse_plan_workspace_bytes, "plan low-rate inverse FFT");
            memory.inverse_plan_workspace_bytes = std::max(memory.native_inverse_plan_workspace_bytes,
                                                           memory.low_inverse_plan_workspace_bytes);
            memory.shared_plan_workspace_bytes = std::max(memory.forward_plan_workspace_bytes,
                                                          memory.inverse_plan_workspace_bytes);

            allocate(rawDevice, static_cast<std::size_t>(config.M) * config.Nraw);
            allocate(frontendDevice, static_cast<std::size_t>(config.M) * config.Nfft);
            allocate(matchedDevice, static_cast<std::size_t>(config.M) * config.Nfft);
            allocate(lowDevice, static_cast<std::size_t>(config.M) * lowLength);
            allocate(outputDevice, static_cast<std::size_t>(config.Nr) * config.Na);
            allocate(frontendTransfer, static_cast<std::size_t>(config.Nfft));
            allocate(templateTransfer, static_cast<std::size_t>(config.K) * config.Nfft);
            allocate(filterTransfer, static_cast<std::size_t>(config.Nfft));
            allocate(mixingTransfer, static_cast<std::size_t>(config.K) * config.Nfft);
            allocate(weightCarrierDevice, static_cast<std::size_t>(config.Na) * config.M);
            allocate(rangeCarrierDevice, static_cast<std::size_t>(config.K) * config.Nr);
            allocate(rangeIndexDevice, static_cast<std::size_t>(config.Nr));
            allocate(delayIndexDevice, static_cast<std::size_t>(config.Na) * config.M);
            allocate(angleIndicesDevice, static_cast<std::size_t>(config.Na));
            allocate(interpolationDevice, config.interpolation_lut.size());
            memory.non_workspace_bytes = memory.explicit_live_bytes;
            if (memory.shared_plan_workspace_bytes > 0) {
                cudaCheck(cudaMalloc(&workspaceDevice, memory.shared_plan_workspace_bytes), "allocate shared FFT workspace");
                allocations.emplace_back(workspaceDevice, memory.shared_plan_workspace_bytes);
                memory.explicit_live_bytes += memory.shared_plan_workspace_bytes;
                memory.explicit_peak_bytes = memory.explicit_live_bytes;
                ++memory.allocation_count;
            }
            fftCheck(cufftSetWorkArea(forwardPlan, workspaceDevice), "set forward FFT workspace");
            fftCheck(cufftSetWorkArea(nativeInversePlan, workspaceDevice), "set native inverse FFT workspace");
            fftCheck(cufftSetWorkArea(lowInversePlan, workspaceDevice), "set low inverse FFT workspace");

            const auto toCuda = [](const std::vector<std::complex<float>>& values) {
                std::vector<cufftComplex> converted(values.size());
                for (std::size_t i = 0; i < values.size(); ++i)
                    converted[i] = make_cuFloatComplex(values[i].real(), values[i].imag());
                return converted;
            };
            const auto g = toCuda(config.G);
            const auto q = toCuda(config.Q);
            const auto h = toCuda(config.H);
            const auto mix = toCuda(config.mix);
            std::vector<double> rangeIndex(config.Nr);
            std::vector<double> delayIndex(static_cast<std::size_t>(config.Na) * config.M);
            std::vector<cufftComplex> rangeCarrier(static_cast<std::size_t>(config.K) * config.Nr);
            std::vector<cufftComplex> weightCarrier(delayIndex.size());
            std::vector<double> weights(config.M);
            double weightSum = 0.0;
            for (int m = 0; m < config.M; ++m) {
                weights[m] = 0.54 - 0.46 * std::cos(2.0 * kPi * m / (config.M - 1));
                weightSum += weights[m];
            }
            const double lowFs = config.fs / config.D;
            for (int r = 0; r < config.Nr; ++r) {
                const double delay = 2.0 * (config.ranges[r] - config.gate) / config.c;
                rangeIndex[r] = delay * lowFs;
                for (int k = 0; k < config.K; ++k) {
                    const double phase = 2.0 * kPi * (config.fc + config.offsets[k]) * delay;
                    rangeCarrier[static_cast<std::size_t>(k) * config.Nr + r] =
                        make_cuFloatComplex(static_cast<float>(std::cos(phase)), static_cast<float>(std::sin(phase)));
                }
            }
            for (int a = 0; a < config.Na; ++a) {
                const double sine = std::sin(config.angles[a]);
                const int k = config.sector_map[a];
                for (int m = 0; m < config.M; ++m) {
                    const std::size_t index = static_cast<std::size_t>(a) * config.M + m;
                    const double delay = -config.positions[m] * sine / config.c;
                    delayIndex[index] = delay * lowFs;
                    const double phase = 2.0 * kPi * (config.fc + config.offsets[k]) * delay;
                    const double weight = weights[m] / weightSum;
                    weightCarrier[index] = make_cuFloatComplex(static_cast<float>(weight * std::cos(phase)),
                                                              static_cast<float>(weight * std::sin(phase)));
                }
            }
            const auto bounds = std::minmax_element(rangeIndex.begin(), rangeIndex.end());
            const int guard = config.interpolation_taps / 2;
            for (double delay : delayIndex) {
                const double minimum = *bounds.first + delay;
                const double maximum = *bounds.second + delay;
                if ((config.interpolation_taps == 0 && (minimum < 0.0 || maximum >= lowLength - 1.0)) ||
                    (config.interpolation_taps > 0 &&
                     (std::floor(minimum) - guard + 1 < 0.0 || std::floor(maximum) + guard >= lowLength))) {
                    throw std::invalid_argument("Image grid outside guarded low-rate interpolation record");
                }
            }
            std::vector<int> angleIndices;
            angleIndices.reserve(config.Na);
            sectorOffsets.push_back(0);
            for (int k = 0; k < config.K; ++k) {
                for (int a = 0; a < config.Na; ++a)
                    if (config.sector_map[a] == k) angleIndices.push_back(a);
                sectorOffsets.push_back(static_cast<int>(angleIndices.size()));
            }
            upload(frontendTransfer, g);
            upload(templateTransfer, q);
            upload(filterTransfer, h);
            upload(mixingTransfer, mix);
            upload(rangeIndexDevice, rangeIndex);
            upload(delayIndexDevice, delayIndex);
            upload(rangeCarrierDevice, rangeCarrier);
            upload(weightCarrierDevice, weightCarrier);
            upload(angleIndicesDevice, angleIndices);
            upload(interpolationDevice, config.interpolation_lut);
            cudaCheck(cudaStreamSynchronize(computeStream), "complete static setup");
            cudaCheck(cudaMemGetInfo(&memory.cuda_free_bytes_after_setup, &memory.cuda_total_bytes), "memory after setup");
            if (diagnostic) {
                stage.frontend_spectrum.resize(static_cast<std::size_t>(config.M) * config.Nfft);
                stage.matched_filter_native.resize(static_cast<std::size_t>(config.K) * config.M * config.Nfft);
                stage.downsampled_matched_filter.resize(static_cast<std::size_t>(config.K) * config.M * lowLength);
            }
        } catch (...) {
            cleanup();
            throw;
        }
    }

    ~Impl() { cleanup(); }

    void cleanup() noexcept {
        if (computeStream) cudaStreamSynchronize(computeStream);
        for (auto stream : beamStreams) if (stream) cudaStreamSynchronize(stream);
        if (forwardPlan) cufftDestroy(forwardPlan);
        if (nativeInversePlan) cufftDestroy(nativeInversePlan);
        if (lowInversePlan) cufftDestroy(lowInversePlan);
        forwardPlan = nativeInversePlan = lowInversePlan = 0;
        for (auto iterator = allocations.rbegin(); iterator != allocations.rend(); ++iterator) {
            cudaFree(iterator->first);
            memory.explicit_live_bytes -= iterator->second;
        }
        allocations.clear();
        if (lowReady) cudaEventDestroy(lowReady);
        lowReady = nullptr;
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
        const std::size_t lowCount = static_cast<std::size_t>(config.M) * lowLength;
        const std::size_t outputCount = static_cast<std::size_t>(config.Nr) * config.Na;
        if (raw.size() != rawCount) throw std::invalid_argument("Raw frame size mismatch");
        output.resize(outputCount);
        cudaCheck(cudaMemcpyAsync(rawDevice, raw.data(), rawCount * sizeof(cufftComplex),
                                  cudaMemcpyHostToDevice, computeStream), "fresh raw frame H2D");
        padRaw<<<blocks(nativeCount), kThreads, 0, computeStream>>>(rawDevice, frontendDevice,
                                                                  config.M, config.Nraw, config.Nfft);
        cudaCheck(cudaGetLastError(), "pad raw frame");
        fftCheck(cufftExecC2C(forwardPlan, frontendDevice, frontendDevice, CUFFT_FORWARD), "raw forward FFT");
        frontend<<<blocks(nativeCount), kThreads, 0, computeStream>>>(frontendDevice, frontendTransfer,
                                                                    nativeCount, config.Nfft);
        cudaCheck(cudaGetLastError(), "single common receive frontend");
        if (diagnostic)
            cudaCheck(cudaMemcpyAsync(stage.frontend_spectrum.data(), frontendDevice,
                                      nativeCount * sizeof(cufftComplex), cudaMemcpyDeviceToHost,
                                      computeStream), "diagnostic frontend D2H");

        for (int k = 0; k < config.K; ++k) {
            const int begin = sectorOffsets[k];
            const int angleCount = sectorOffsets[k + 1] - begin;
            if (angleCount == 0 && !diagnostic) continue;
            constructMf<<<blocks(nativeCount), kThreads, 0, computeStream>>>(
                frontendDevice, templateTransfer + static_cast<std::size_t>(k) * config.Nfft,
                matchedDevice, nativeCount, config.Nfft);
            cudaCheck(cudaGetLastError(), "construct native matched spectrum");
            fftCheck(cufftExecC2C(nativeInversePlan, matchedDevice, matchedDevice, CUFFT_INVERSE), "native matched inverse FFT");
            normalize<<<blocks(nativeCount), kThreads, 0, computeStream>>>(matchedDevice, nativeCount,
                                                                          1.0f / config.Nfft);
            cudaCheck(cudaGetLastError(), "normalize native matched time");
            if (diagnostic)
                cudaCheck(cudaMemcpyAsync(stage.matched_filter_native.data() + static_cast<std::size_t>(k) * nativeCount,
                                          matchedDevice, nativeCount * sizeof(cufftComplex), cudaMemcpyDeviceToHost,
                                          computeStream), "diagnostic native MF D2H");
            mixDown<<<blocks(nativeCount), kThreads, 0, computeStream>>>(
                matchedDevice, mixingTransfer + static_cast<std::size_t>(k) * config.Nfft,
                nativeCount, config.Nfft);
            cudaCheck(cudaGetLastError(), "subband offset demodulation");
            fftCheck(cufftExecC2C(forwardPlan, matchedDevice, matchedDevice, CUFFT_FORWARD), "demodulated forward FFT");
            filterFold<<<blocks(lowCount), kThreads, 0, computeStream>>>(
                matchedDevice, filterTransfer, lowDevice, config.M, config.Nfft, lowLength, config.D);
            cudaCheck(cudaGetLastError(), "anti-alias filtering and alias-fold decimation");
            fftCheck(cufftExecC2C(lowInversePlan, lowDevice, lowDevice, CUFFT_INVERSE), "low-rate inverse FFT");
            normalize<<<blocks(lowCount), kThreads, 0, computeStream>>>(lowDevice, lowCount,
                                                                       1.0f / lowLength);
            cudaCheck(cudaGetLastError(), "normalize low-rate time");
            if (diagnostic)
                cudaCheck(cudaMemcpyAsync(stage.downsampled_matched_filter.data() + static_cast<std::size_t>(k) * lowCount,
                                          lowDevice, lowCount * sizeof(cufftComplex), cudaMemcpyDeviceToHost,
                                          computeStream), "diagnostic low-rate MF D2H");
            cudaCheck(cudaEventRecord(lowReady, computeStream), "record low-rate ready event");
            std::array<bool, kBeamStreams> used{};
            const int chunk = (angleCount + kBeamStreams - 1) / kBeamStreams;
            for (int stream = 0; stream < kBeamStreams && chunk > 0; ++stream) {
                const int localBegin = stream * chunk;
                if (localBegin >= angleCount) break;
                const int localCount = std::min(chunk, angleCount - localBegin);
                cudaCheck(cudaStreamWaitEvent(beamStreams[stream], lowReady, 0), "DAS waits for low-rate record");
                lowRateDas<<<blocks(static_cast<std::size_t>(localCount) * config.Nr), kThreads, 0, beamStreams[stream]>>>(
                    lowDevice, outputDevice, rangeIndexDevice, delayIndexDevice,
                    weightCarrierDevice, rangeCarrierDevice + static_cast<std::size_t>(k) * config.Nr,
                    angleIndicesDevice + begin + localBegin, localCount, config.M, config.Nr,
                    config.Na, lowLength, interpolationDevice, config.interpolation_taps, config.lut_phases);
                cudaCheck(cudaGetLastError(), "low-rate fractional-delay DAS");
                used[stream] = true;
            }
            // A single low buffer is reused between sectors; its producer and
            // final D2H explicitly wait for every consumer on all DAS streams.
            for (int stream = 0; stream < kBeamStreams; ++stream) {
                if (!used[stream]) continue;
                cudaCheck(cudaEventRecord(beamComplete[stream], beamStreams[stream]), "record DAS completion");
                cudaCheck(cudaStreamWaitEvent(computeStream, beamComplete[stream], 0), "protect shared low-rate buffer");
            }
        }
        cudaCheck(cudaMemcpyAsync(output.data(), outputDevice, outputCount * sizeof(cufftComplex),
                                  cudaMemcpyDeviceToHost, computeStream), "complex output D2H");
        cudaCheck(cudaStreamSynchronize(computeStream), "complete frame reconstruction");
    }
};

Engine::Engine(const Config& config, bool diagnostic)
    : impl_(std::make_unique<Impl>(config, diagnostic)) {}
Engine::~Engine() = default;
Engine::Engine(Engine&&) noexcept = default;
Engine& Engine::operator=(Engine&&) noexcept = default;

void Engine::process(const std::vector<std::complex<float>>& raw,
                     std::vector<std::complex<float>>& output) {
    if (!impl_) throw std::logic_error("Moved-from v3 engine");
    impl_->run(raw, output);
}

const Diagnostics& Engine::diagnostics() const {
    if (!impl_ || !impl_->diagnostic) throw std::logic_error("v3 diagnostics disabled");
    return impl_->stage;
}

MemoryReport Engine::memoryReport() const {
    if (!impl_) throw std::logic_error("Moved-from v3 engine");
    return impl_->memory;
}

std::string Engine::implementationName() const {
    if (!impl_) throw std::logic_error("Moved-from v3 engine");
    return "cached_lowrate_D" + std::to_string(impl_->config.D) +
        (impl_->config.interpolation_taps == 0 ? "_linear" :
         "_sinc" + std::to_string(impl_->config.interpolation_taps));
}

} // namespace sonar_v3
