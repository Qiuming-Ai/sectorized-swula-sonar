#pragma once

#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace sonar_v2 {

// All spectra and records use channel-major, complex-baseband storage.
// Q is the non-conjugated FFT of the unit-energy effective Tx/Rx template.
struct Config {
    int M = 0;
    int Nraw = 8192;
    int Nfft = 16384;
    int OS = 8;
    int Nr = 321;
    int Na = 401;
    int K = 1;
    double fs = 400000.0;
    double fc = 450000.0;
    double c = 1500.0;
    double gate = 28.5;
    std::vector<double> positions;
    std::vector<double> angles;  // radians, inside the evaluation sector
    std::vector<double> ranges;
    std::vector<int> sector_map; // Na entries in [0, K), supplied by the fixture
    std::vector<std::complex<float>> G; // Nfft, single common receive frontend
    std::vector<std::complex<float>> Q; // K * Nfft
};

struct MemoryReport {
    std::size_t explicit_live_bytes = 0;
    std::size_t explicit_peak_bytes = 0;
    std::size_t non_workspace_bytes = 0;
    std::size_t forward_plan_workspace_bytes = 0;
    std::size_t inverse_plan_workspace_bytes = 0;
    std::size_t shared_plan_workspace_bytes = 0;
    std::size_t allocation_count = 0;
    std::size_t cuda_free_bytes_before = 0;
    std::size_t cuda_free_bytes_after_setup = 0;
    std::size_t cuda_total_bytes = 0;
    // This report is explicit device-allocation accounting, NOT module/process
    // peak memory. CUDA context, cuFFT internals, driver and host memory are not
    // implied to be included by explicit_peak_bytes.
};

struct Diagnostics {
    // Only produced when diagnostic=true, not in the performance path.
    std::vector<std::complex<float>> frontend_spectrum; // M * Nfft
    std::vector<std::complex<float>> matched_filter_native; // K * M * Nfft
};

// cached=true batches each assigned-sector angle block in CUDA DAS kernels.
// cached=false launches one angle per DAS kernel. Neither caches frame results.
// For the allocation/plan uncached comparator, construct and destroy one Engine
// for each frame in the outer harness; include that lifecycle in wall timing.
class Engine {
public:
    Engine(const Config& config, bool cached, bool diagnostic = false);
    ~Engine();
    Engine(Engine&&) noexcept;
    Engine& operator=(Engine&&) noexcept;
    Engine(const Engine&) = delete;
    Engine& operator=(const Engine&) = delete;

    // Synchronously completes H2D -> reconstruction -> D2H for a fresh raw frame.
    // Output is complex, range-major [Nr, Na], without normalization or dB mapping.
    void process(const std::vector<std::complex<float>>& raw,
                 std::vector<std::complex<float>>& output);
    const Diagnostics& diagnostics() const;
    MemoryReport memoryReport() const;
    std::string implementationName() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace sonar_v2
