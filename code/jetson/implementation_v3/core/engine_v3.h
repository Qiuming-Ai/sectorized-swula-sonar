#pragma once

#include "engine.h"

#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace sonar_v3 {

// Static metadata is inherited from the frozen v2 full-rate FDM fixtures.
// No frame values, matched-filter outputs or final images are cached.
struct Config : sonar_v2::Config {
    int D = 4;
    int interpolation_taps = 0; // 0: linear; otherwise even, finite sinc taps.
    int lut_phases = 1024;
    std::vector<double> offsets; // K baseband offsets, in Hz.
    std::vector<std::complex<float>> H; // Nfft circular zero-phase AA response.
    std::vector<std::complex<float>> mix; // K*Nfft exp(-j*2*pi*offset*n/fs).
    // Phase-major [(lut_phases+1), interpolation_taps], normalized each row.
    // j = -taps/2+1,...,taps/2 relative to floor(sample). Adjacent LUT rows
    // are linearly blended. The independent fixture generator defines rows.
    std::vector<float> interpolation_lut;
};

struct MemoryReport : sonar_v2::MemoryReport {
    std::size_t native_inverse_plan_workspace_bytes = 0;
    std::size_t low_inverse_plan_workspace_bytes = 0;
};

struct Diagnostics {
    std::vector<std::complex<float>> frontend_spectrum; // M*Nfft
    std::vector<std::complex<float>> matched_filter_native; // K*M*Nfft
    std::vector<std::complex<float>> downsampled_matched_filter; // K*M*(Nfft/D)
};

class Engine {
public:
    explicit Engine(const Config& config, bool diagnostic = false);
    ~Engine();
    Engine(Engine&&) noexcept;
    Engine& operator=(Engine&&) noexcept;
    Engine(const Engine&) = delete;
    Engine& operator=(const Engine&) = delete;

    // H2D -> common frontend -> native MF -> offset demodulation -> AA /
    // decimation -> low-rate fractional-delay DAS -> complex D2H. Synchronous.
    void process(const std::vector<std::complex<float>>& raw,
                 std::vector<std::complex<float>>& output);
    const Diagnostics& diagnostics() const;
    MemoryReport memoryReport() const;
    std::string implementationName() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace sonar_v3
