#include "engine.h"
#include "engine_v3.h"
#include <cuda_runtime.h>
#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <time.h>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;
using C32 = std::complex<float>;
using C64 = std::complex<double>;
constexpr int kWarmups = 50, kMeasurements = 200, kFixtureCount = 20;
constexpr std::int64_t kSecond = 1000000000LL, kFrameTimeout = 30 * kSecond;
constexpr double kPi = 3.141592653589793238462643383279502884;

std::int64_t nowNs() {
    timespec t{};
    if (clock_gettime(CLOCK_MONOTONIC, &t) != 0)
        throw std::runtime_error("CLOCK_MONOTONIC failure: " + std::string(std::strerror(errno)));
    return std::int64_t(t.tv_sec) * kSecond + t.tv_nsec;
}
double milliseconds(std::int64_t end, std::int64_t start) { return double(end - start) / 1e6; }
void cudaCheck(cudaError_t e, const char* operation) {
    if (e != cudaSuccess) throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(e));
}
std::string jsonString(const std::string& s) {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : s) {
        if (c == '"' || c == '\\') out << '\\' << char(c);
        else if (c == '\n') out << "\\n";
        else if (c == '\r') out << "\\r";
        else if (c == '\t') out << "\\t";
        else if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec << std::setfill(' ');
        else out << char(c);
    }
    out << '"'; return out.str();
}
std::string jsonNumber(double value) {
    if (!std::isfinite(value)) return "null";
    std::ostringstream out; out << std::setprecision(17) << value; return out.str();
}
std::string csvString(const std::string& s) {
    std::string result = "\"";
    for (char c : s) { if (c == '"') result += '"'; result += c; }
    return result + '"';
}
void event(const std::string& name, std::int64_t ns, const std::string& additional = "") {
    std::cout << "{\"event\":" << jsonString(name) << ",\"monotonic_ns\":" << ns
              << ",\"pid\":" << getpid() << additional << "}\n" << std::flush;
}
std::ofstream newOutput(const fs::path& path, bool binary = false) {
    if (fs::exists(path)) throw std::runtime_error("Refusing output overwrite: " + path.string());
    std::ofstream f(path, std::ios::out | (binary ? std::ios::binary : std::ios::openmode(0)));
    if (!f) throw std::runtime_error("Cannot create output: " + path.string());
    f << std::setprecision(17); return f;
}
template<class T> std::vector<T> readBinary(const fs::path& path, std::size_t count) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f || f.tellg() != std::streamoff(count * sizeof(T)))
        throw std::runtime_error("Missing/wrong-sized input: " + path.string());
    std::vector<T> values(count); f.seekg(0);
    f.read(reinterpret_cast<char*>(values.data()), std::streamsize(count * sizeof(T)));
    if (!f) throw std::runtime_error("Input read failure: " + path.string());
    return values;
}
template<class T> void writeBinary(const fs::path& path, const std::vector<T>& values) {
    auto f = newOutput(path, true);
    f.write(reinterpret_cast<const char*>(values.data()), std::streamsize(values.size() * sizeof(T)));
    if (!f) throw std::runtime_error("Output write failure: " + path.string());
}
std::map<std::string, std::string> readMeta(const fs::path& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Missing metadata: " + path.string());
    std::map<std::string, std::string> values; std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        auto equals = line.find('=');
        if (equals == std::string::npos || equals == 0) throw std::runtime_error("Invalid metadata line");
        if (!values.emplace(line.substr(0, equals), line.substr(equals + 1)).second)
            throw std::runtime_error("Duplicate metadata key");
    }
    return values;
}
void approximately(double actual, double expected, const char* label, double tolerance = 1e-9) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance)
        throw std::runtime_error(std::string("Frozen parameter mismatch: ") + label);
}
sonar_v2::Config loadBase(const fs::path& bank, const std::string& method) {
    auto m = readMeta(bank / "meta.txt");
    auto ni = [&](const char* k) { return std::stoi(m.at(k)); };
    auto nd = [&](const char* k) { return std::stod(m.at(k)); };
    if (ni("schema_version") != 2) throw std::runtime_error("Expected v2 input schema");
    sonar_v2::Config c;
    c.M = ni("channels"); c.Nraw = ni("record_samples"); c.Nfft = ni("fft_samples");
    c.OS = ni("oversampling"); c.Nr = ni("range_rows"); c.Na = ni("angle_columns");
    c.K = ni("waveforms"); c.fs = nd("fs"); c.fc = nd("fc"); c.c = nd("c"); c.gate = nd("gate_range");
    if (c.M != (method == "dense" ? 128 : 32) || c.K != (method == "dense" ? 1 : 6) ||
        c.Nraw != 8192 || c.Nfft != 16384 || c.OS != 8 || c.Nr != 321 || c.Na != 401)
        throw std::runtime_error("Frozen base dimensions do not match selected method");
    approximately(c.fs, 400000., "fs"); approximately(c.fc, 450000., "fc");
    approximately(c.c, 1500., "sound_speed"); approximately(c.gate, 28.5, "gate");
    c.positions = readBinary<double>(bank / "positions.bin", c.M);
    c.angles = readBinary<double>(bank / "angles.bin", c.Na);
    c.ranges = readBinary<double>(bank / "ranges.bin", c.Nr);
    c.sector_map = readBinary<int>(bank / "sector_map.bin", c.Na);
    c.G = readBinary<C32>(bank / "G.bin", c.Nfft);
    c.Q = readBinary<C32>(bank / "Q.bin", std::size_t(c.K) * c.Nfft);
    const double pitch = (method == "dense" ? 0.5 : 2.0) * c.c / c.fc;
    for (int channel = 0; channel < c.M; ++channel)
        approximately(c.positions[channel], (channel - (c.M - 1) / 2.0) * pitch, "Rx coordinate");
    for (int angle = 0; angle < c.Na; ++angle) {
        approximately(c.angles[angle], (-50.0 + angle * .25) * kPi / 180., "angle grid");
        if (c.sector_map[angle] < 0 || c.sector_map[angle] >= c.K)
            throw std::runtime_error("Invalid sector assignment");
    }
    for (int range = 0; range < c.Nr; ++range)
        approximately(c.ranges[range], 30.0 + range * .0025, "range grid");
    return c;
}
sonar_v3::Config loadFdm(const sonar_v2::Config& b, const fs::path& bank) {
    auto m = readMeta(bank / "meta.txt");
    if (std::stoi(m.at("schema_version")) != 3 || std::stoi(m.at("decimation")) != 4 ||
        std::stoi(m.at("interpolation_taps")) != 16 || std::stoi(m.at("lut_phases")) != 1024 ||
        std::stoi(m.at("fir_taps")) != 129 || std::stoi(m.at("low_samples")) != 4096 ||
        std::stoi(m.at("native_fft_samples")) != 16384)
        throw std::runtime_error("Only frozen D4 sinc16 FDM bank is allowed");
    approximately(std::stod(m.at("fs_low")), 100000., "FDM fs_low");
    approximately(std::stod(m.at("cutoff_hz")), 30000., "FDM AA cutoff");
    approximately(std::stod(m.at("kaiser_beta")), 8.6, "FDM Kaiser beta");
    approximately(std::stod(m.at("fir_delay_compensated_samples")), 64., "FDM AA delay compensation");
    sonar_v3::Config c; static_cast<sonar_v2::Config&>(c) = b;
    c.D = 4; c.interpolation_taps = 16; c.lut_phases = 1024;
    c.offsets = readBinary<double>(bank / "offsets.bin", c.K);
    constexpr std::array<double, 6> offsets{-24000., -8000., 40000., 24000., 8000., -38000.};
    for (int k = 0; k < c.K; ++k) approximately(c.offsets[k], offsets[k], "FDM offset");
    c.H = readBinary<C32>(bank / "H.bin", c.Nfft);
    c.mix = readBinary<C32>(bank / "mix.bin", std::size_t(c.K) * c.Nfft);
    c.interpolation_lut = readBinary<float>(bank / "lut.bin", std::size_t(c.lut_phases + 1) * c.interpolation_taps);
    return c;
}
std::string noiseId(int i) {
    std::ostringstream out; out << "noise" << std::setw(2) << std::setfill('0') << i; return out.str();
}
std::vector<std::string> verificationCases() {
    std::vector<std::string> result{"paper_clean", "zero", "boundary_points", "offgrid", "complex_impulse"};
    for (int i = 0; i < kFixtureCount; ++i) result.push_back(noiseId(i));
    std::sort(result.begin(), result.end()); return result;
}
long statusKiB(const std::string& key) {
    std::ifstream f("/proc/self/status"); std::string line;
    while (std::getline(f, line)) if (line.rfind(key + ":", 0) == 0) {
        std::istringstream in(line.substr(key.size() + 1)); long value = -1; in >> value; return value;
    }
    return -1;
}
struct HostSnapshot { long rss = -1, highWater = -1; std::int64_t ns = 0; };
HostSnapshot hostSnapshot() { return {statusKiB("VmRSS"), statusKiB("VmHWM"), nowNs()}; }

class Processor {
public:
    Processor(const sonar_v2::Config& base, const sonar_v3::Config* fdm) {
        if (fdm) v3 = std::make_unique<sonar_v3::Engine>(*fdm, false);
        else v2 = std::make_unique<sonar_v2::Engine>(base, true, false);
    }
    void process(const std::vector<C32>& raw, std::vector<C32>& output) {
        if (v3) v3->process(raw, output); else v2->process(raw, output);
    }
    sonar_v3::MemoryReport memory() const {
        if (v3) return v3->memoryReport();
        sonar_v3::MemoryReport result;
        static_cast<sonar_v2::MemoryReport&>(result) = v2->memoryReport(); return result;
    }
    std::string name() const { return v3 ? v3->implementationName() : v2->implementationName(); }
private:
    std::unique_ptr<sonar_v2::Engine> v2;
    std::unique_ptr<sonar_v3::Engine> v3;
};
struct Comparison {
    double relative = 0., maximum = 0., correlation = 0.;
    bool zeroReference = false, correlationNA = true, finite = true, shape = true, passed = false;
};
Comparison compare(const std::vector<C32>& gpu, const std::vector<C64>& cpu) {
    Comparison result; result.shape = gpu.size() == cpu.size();
    if (!result.shape || cpu.empty()) return result;
    long double numerator = 0., denominator = 0., sumCpu = 0., sumGpu = 0.;
    std::vector<double> magnitudeCpu(cpu.size()), magnitudeGpu(cpu.size());
    for (std::size_t i = 0; i < cpu.size(); ++i) {
        const C64 actual(gpu[i]);
        if (!std::isfinite(actual.real()) || !std::isfinite(actual.imag()) ||
            !std::isfinite(cpu[i].real()) || !std::isfinite(cpu[i].imag())) {
            result.finite = false; continue;
        }
        magnitudeCpu[i] = std::abs(cpu[i]); magnitudeGpu[i] = std::abs(actual);
        if (!std::isfinite(magnitudeCpu[i]) || !std::isfinite(magnitudeGpu[i])) result.finite = false;
        const double error = std::abs(actual - cpu[i]); result.maximum = std::max(result.maximum, error);
        numerator += std::norm(actual - cpu[i]); denominator += std::norm(cpu[i]);
        sumCpu += magnitudeCpu[i]; sumGpu += magnitudeGpu[i];
    }
    result.zeroReference = denominator == 0.;
    result.relative = denominator > 0. ? std::sqrt(double(numerator / denominator)) : result.maximum;
    if (result.finite) {
        long double covariance = 0., varianceCpu = 0., varianceGpu = 0.;
        const long double meanCpu = sumCpu / cpu.size(), meanGpu = sumGpu / cpu.size();
        for (std::size_t i = 0; i < cpu.size(); ++i) {
            const long double x = magnitudeCpu[i] - meanCpu, y = magnitudeGpu[i] - meanGpu;
            covariance += x * y; varianceCpu += x * x; varianceGpu += y * y;
        }
        const auto cpuBounds = std::minmax_element(magnitudeCpu.begin(), magnitudeCpu.end());
        const auto gpuBounds = std::minmax_element(magnitudeGpu.begin(), magnitudeGpu.end());
        result.correlationNA = *cpuBounds.first == *cpuBounds.second;
        if (!result.correlationNA) {
            // Constant candidate with nonconstant reference is an undefined
            // correlation and a FAILED gate, not a reference-constant N/A.
            result.correlation = std::numeric_limits<double>::quiet_NaN();
            if (*gpuBounds.first != *gpuBounds.second && varianceCpu > 0. && varianceGpu > 0.) {
                result.correlation = double(covariance / std::sqrt(varianceCpu * varianceGpu));
                if (result.correlation > 1. && result.correlation <= 1. + 1e-12) result.correlation = 1.;
                if (result.correlation < -1. && result.correlation >= -1. - 1e-12) result.correlation = -1.;
            }
        }
    }
    result.passed = result.finite && result.shape && std::isfinite(result.relative) &&
        result.relative <= (result.zeroReference ? 1e-7 : 1e-4) &&
        (result.correlationNA || (std::isfinite(result.correlation) && result.correlation >= .999 && result.correlation <= 1.));
    return result;
}
constexpr std::array<std::size_t, 14> checkpointIndices{
    80 * 401 + 40, 80 * 401 + 168, 80 * 401 + 296,
    160 * 401 + 104, 160 * 401 + 232, 160 * 401 + 360,
    240 * 401 + 40, 240 * 401 + 168, 240 * 401 + 296,
    0, 400, 320 * 401, 321 * 401 - 1, 160 * 401 + 200
};
struct Frame {
    std::string id; std::vector<C32> raw;
    std::array<C64, checkpointIndices.size()> checkpoints{};
};
struct LightCheck { bool passed = false, finite = true, shape = true; double relative = 0., maximum = 0.; };
LightCheck lightCheck(const std::vector<C32>& output, const Frame& frame) {
    LightCheck check; check.shape = output.size() == 321 * 401;
    if (!check.shape) return check;
    long double numerator = 0., denominator = 0.;
    for (std::size_t i = 0; i < checkpointIndices.size(); ++i) {
        const C64 value(output[checkpointIndices[i]]), reference = frame.checkpoints[i];
        check.finite &= std::isfinite(value.real()) && std::isfinite(value.imag());
        const C64 error = value - reference;
        numerator += std::norm(error); denominator += std::norm(reference);
        check.maximum = std::max(check.maximum, std::abs(error));
    }
    check.relative = denominator > 0. ? std::sqrt(double(numerator / denominator)) : check.maximum;
    check.passed = check.shape && check.finite && std::isfinite(check.relative) &&
        check.relative <= (denominator > 0. ? 1e-4 : 1e-7);
    return check;
}
struct Row {
    std::string phase, fixture, status, message; int index = 0;
    std::int64_t begin = 0, ready = 0, checked = 0; LightCheck check;
};
void writeRows(const fs::path& path, const std::vector<Row>& rows, const std::string& method, int session) {
    auto f = newOutput(path);
    f << "session,method,phase,block_id,fixture_id,t_start_ns,t_result_ready_ns,t_light_check_end_ns,wall_ms,light_check_ms,shape_pass,checkpoint_finite,checkpoint_relative_l2,checkpoint_max_abs,status,message\n";
    for (const auto& r : rows)
        f << session << ',' << method << ',' << r.phase << ',' << r.index << ',' << r.fixture << ','
          << r.begin << ',' << r.ready << ',' << r.checked << ',' << milliseconds(r.ready, r.begin) << ','
          << milliseconds(r.checked, r.ready) << ',' << r.check.shape << ',' << r.check.finite << ','
          << r.check.relative << ',' << r.check.maximum << ',' << r.status << ',' << csvString(r.message) << '\n';
    if (!f) throw std::runtime_error("Timing output failure");
}
bool stopRequested(const fs::path& stopFile) { return !stopFile.empty() && fs::exists(stopFile); }
bool runFrame(Processor& processor, const Frame& frame, std::vector<C32>& image,
              std::vector<Row>& rows, const std::string& phase, int index, const fs::path& stopFile) {
    if (stopRequested(stopFile)) return false;
    Row row; row.phase = phase; row.index = index; row.fixture = frame.id; row.begin = nowNs();
    try {
        processor.process(frame.raw, image); row.ready = nowNs();
        row.check = lightCheck(image, frame); row.checked = nowNs();
        row.status = row.ready - row.begin > kFrameTimeout ? "timeout" :
            (!row.check.passed ? "light_check_failed" : (stopRequested(stopFile) ? "external_stop" : "ok"));
    } catch (const std::exception& e) {
        row.ready = nowNs(); row.checked = row.ready; row.status = "compute_error"; row.message = e.what();
    }
    const bool success = row.status == "ok"; rows.push_back(std::move(row)); return success;
}
void writeMemory(const fs::path& path, const sonar_v3::MemoryReport& m, const HostSnapshot& context,
                 const HostSnapshot& live, const std::string& implementation, const std::string& phase) {
    auto f = newOutput(path);
    f << "{\n\"implementation\":" << jsonString(implementation)
      << ",\n\"snapshot_phase\":" << jsonString(phase)
      << ",\n\"explicit_live_bytes\":" << m.explicit_live_bytes << ",\n\"explicit_peak_bytes\":" << m.explicit_peak_bytes
      << ",\n\"non_workspace_bytes\":" << m.non_workspace_bytes
      << ",\n\"forward_plan_workspace_bytes\":" << m.forward_plan_workspace_bytes
      << ",\n\"inverse_plan_workspace_bytes\":" << m.inverse_plan_workspace_bytes
      << ",\n\"native_inverse_plan_workspace_bytes\":" << m.native_inverse_plan_workspace_bytes
      << ",\n\"low_inverse_plan_workspace_bytes\":" << m.low_inverse_plan_workspace_bytes
      << ",\n\"shared_plan_workspace_bytes\":" << m.shared_plan_workspace_bytes
      << ",\n\"allocation_count\":" << m.allocation_count
      << ",\n\"cuda_free_bytes_before_setup\":" << m.cuda_free_bytes_before
      << ",\n\"cuda_free_bytes_after_setup\":" << m.cuda_free_bytes_after_setup
      << ",\n\"cuda_total_bytes\":" << m.cuda_total_bytes
      << ",\n\"host_context_baseline\":{\"monotonic_ns\":" << context.ns << ",\"RSS_KiB\":" << context.rss << ",\"VmHWM_KiB\":" << context.highWater << "}"
      << ",\n\"host_live_phase\":{\"monotonic_ns\":" << live.ns << ",\"RSS_KiB\":" << live.rss << ",\"VmHWM_KiB\":" << live.highWater << "}"
      << ",\n\"context_baseline_includes_preloaded_host_inputs_and_config\":true"
      << ",\n\"not_total_module_peak\":true,\n\"unified_memory_quantities_must_not_be_summed\":true\n}\n";
}
void bufferTable(const fs::path& path, const sonar_v2::Config& c, const sonar_v3::Config* d,
                 const sonar_v3::MemoryReport& memory) {
    auto f = newOutput(path); f << "buffer,shape,dtype,bytes,scope\n"; std::size_t total = 0;
    auto add = [&](const char* name, const std::string& shape, const char* dtype, std::size_t count, std::size_t elementBytes) {
        const std::size_t bytes = count * elementBytes; total += bytes;
        f << name << ',' << csvString(shape) << ',' << dtype << ',' << bytes << ",source_derived_simultaneously_live_device_allocation\n";
    };
    auto shape2 = [](int a, int b) { return std::to_string(a) + "x" + std::to_string(b); };
    const auto mn = std::size_t(c.M) * c.Nfft, geometry = std::size_t(c.Na) * c.M;
    add("raw", shape2(c.M, c.Nraw), "complex_fp32", std::size_t(c.M) * c.Nraw, 8);
    add("frontend_spectrum", shape2(c.M, c.Nfft), "complex_fp32", mn, 8);
    if (d) {
        add("matched_native", shape2(c.M, c.Nfft), "complex_fp32", mn, 8);
        add("low_record", shape2(c.M, c.Nfft / d->D), "complex_fp32", mn / d->D, 8);
    } else add("fine_record", shape2(c.M, c.Nfft * c.OS), "complex_fp32", mn * c.OS, 8);
    add("output", shape2(c.Nr, c.Na), "complex_fp32", std::size_t(c.Nr) * c.Na, 8);
    add("frontend_transfer", std::to_string(c.Nfft), "complex_fp32", c.Nfft, 8);
    add("template_transfer", shape2(c.K, c.Nfft), "complex_fp32", std::size_t(c.K) * c.Nfft, 8);
    if (d) {
        add("AA_response", std::to_string(c.Nfft), "complex_fp32", c.Nfft, 8);
        add("mixing", shape2(c.K, c.Nfft), "complex_fp32", std::size_t(c.K) * c.Nfft, 8);
    }
    add("weight_carrier", shape2(c.Na, c.M), "complex_fp32", geometry, 8);
    add("range_carrier", d ? shape2(c.K, c.Nr) : std::to_string(c.Nr), "complex_fp32", std::size_t(c.Nr) * (d ? c.K : 1), 8);
    add("range_index", std::to_string(c.Nr), "fp64", c.Nr, 8);
    add("delay_index", shape2(c.Na, c.M), "fp64", geometry, 8);
    add("angle_indices", std::to_string(c.Na), "int32", c.Na, 4);
    if (d) add("fractional_delay_LUT", shape2(d->lut_phases + 1, d->interpolation_taps), "fp32", d->interpolation_lut.size(), 4);
    if (total != memory.non_workspace_bytes) throw std::runtime_error("Device buffer ledger does not match Engine report");
    f << "shared_fft_workspace," << csvString("max(plan workspace requirements)") << ",uint8," << memory.shared_plan_workspace_bytes << ",single_explicit_shared_workspace_allocation\n";
    if (total + memory.shared_plan_workspace_bytes != memory.explicit_live_bytes ||
        memory.explicit_live_bytes != memory.explicit_peak_bytes)
        throw std::runtime_error("Explicit device high-water ledger mismatch");
}

int main(int argc, char** argv) {
    fs::path activeOut; bool ownsOut = false;
    try {
        const std::array<std::string, 8> keys{"mode", "method", "data", "fdm-bank", "permit", "session", "out", "stop-file"};
        std::map<std::string, std::string> args;
        for (int i = 1; i < argc; ++i) {
            const std::string option(argv[i]);
            if (option.rfind("--", 0) != 0 || i + 1 >= argc) throw std::runtime_error("Arguments require --key value pairs");
            const auto key = option.substr(2);
            if (std::find(keys.begin(), keys.end(), key) == keys.end() || !args.emplace(key, argv[++i]).second)
                throw std::runtime_error("Unknown or duplicate option: " + option);
        }
        const std::string mode = args.at("mode"), method = args.at("method");
        if (mode != "verify" && mode != "latency" && mode != "memory" && mode != "power") throw std::runtime_error("Unknown mode");
        if (method != "dense" && method != "cdm" && method != "fdm") throw std::runtime_error("Unknown method");
        const int session = std::stoi(args.at("session"));
        if (session < 0 || session > 2) throw std::runtime_error("Session must be 0, 1 or 2");
        auto permit = readMeta(args.at("permit"));
        if (permit.at("selected_engineering") != "1" || permit.at("equivalence_claimed") != "0" ||
            (mode != "verify" && permit.at("implementation_verified") != "1"))
            throw std::runtime_error("Formal engineering/implementation permit required");
        if (permit.at("implementation_verified") != "0" && permit.at("implementation_verified") != "1")
            throw std::runtime_error("Invalid implementation_verified field");
        const fs::path data = args.at("data"), fdmBank = args.at("fdm-bank");
        const fs::path stopFile = args.count("stop-file") ? fs::path(args.at("stop-file")) : fs::path();
        activeOut = args.at("out");
        if (activeOut.empty() || (fs::exists(activeOut) && (!fs::is_directory(activeOut) || !fs::is_empty(activeOut))))
            throw std::runtime_error("Output must be new or empty; refusing overwrite");
        fs::create_directories(activeOut);
        if (!fs::create_directory(activeOut / ".runner_lock")) throw std::runtime_error("Output already claimed");
        ownsOut = true;
        const fs::path methodRoot = data / method;
        const auto base = loadBase(methodRoot / "static", method);
        std::unique_ptr<sonar_v3::Config> fdm;
        if (method == "fdm") fdm = std::make_unique<sonar_v3::Config>(loadFdm(base, fdmBank));
        const std::size_t rawCount = std::size_t(base.M) * base.Nraw, imageCount = std::size_t(base.Nr) * base.Na;
        auto referencePath = [&](const std::string& id) {
            return method == "fdm" ? fdmBank / "cases" / id / "reference.bin" : methodRoot / id / "reference.bin";
        };
        std::vector<Frame> frames;
        if (mode != "verify") {
            frames.reserve(kFixtureCount);
            for (int i = 0; i < kFixtureCount; ++i) {
                Frame frame; frame.id = noiseId(i); frame.raw = readBinary<C32>(methodRoot / frame.id / "raw.bin", rawCount);
                auto reference = readBinary<C64>(referencePath(frame.id), imageCount);
                for (std::size_t j = 0; j < checkpointIndices.size(); ++j) frame.checkpoints[j] = reference[checkpointIndices[j]];
                frames.push_back(std::move(frame));
            }
        }
        std::vector<C32> image(imageCount);
        const auto contextBegin = nowNs();
        cudaCheck(cudaFree(nullptr), "initialize CUDA context");
        cudaCheck(cudaDeviceSynchronize(), "quiescent CUDA context");
        const auto contextReady = nowNs();
        const auto contextHost = hostSnapshot();
        const auto initializeBegin = nowNs();
        auto processor = std::make_unique<Processor>(base, fdm.get());
        const auto initializeReady = nowNs();
        const std::string implementation = processor->name();
        event("initialized", initializeReady, ",\"method\":" + jsonString(method) + ",\"mode\":" + jsonString(mode) + ",\"session\":" + std::to_string(session));
        bool success = true; std::vector<Row> rows; std::string lastId; bool haveImage = false;
        sonar_v3::MemoryReport memoryReport; HostSnapshot liveHost;
        if (mode == "verify") {
            auto csv = newOutput(activeOut / "verification.csv");
            csv << "method,fixture,range_rows,angle_columns,complex_relative_l2,complex_max_abs,magnitude_correlation,correlation_na,reference_is_zero,finite,shape_pass,overall_pass\n";
            int count = 0;
            for (const auto& id : verificationCases()) {
                if (stopRequested(stopFile)) { success = false; break; }
                auto raw = readBinary<C32>(methodRoot / id / "raw.bin", rawCount);
                auto cpu = readBinary<C64>(referencePath(id), imageCount);
                processor->process(raw, image); const auto result = compare(image, cpu);
                writeBinary(activeOut / (id + ".bin"), image);
                csv << method << ',' << id << ',' << base.Nr << ',' << base.Na << ',';
                if (!result.zeroReference) csv << result.relative;
                csv << ',' << result.maximum << ',';
                if (!result.correlationNA) csv << result.correlation;
                csv << ',' << result.correlationNA << ',' << result.zeroReference << ',' << result.finite << ',' << result.shape << ',' << result.passed << '\n';
                auto f = newOutput(activeOut / (id + ".json"));
                f << "{\"method\":" << jsonString(method) << ",\"fixture\":" << jsonString(id)
                  << ",\"shape\":[" << base.Nr << ',' << base.Na << "],\"finite\":" << (result.finite ? "true" : "false")
                  << ",\"shape_pass\":" << (result.shape ? "true" : "false")
                  << ",\"complex_relative_l2\":" << (result.finite && !result.zeroReference ? jsonNumber(result.relative) : "null")
                  << ",\"complex_max_abs\":" << (result.finite ? jsonNumber(result.maximum) : "null")
                  << ",\"magnitude_correlation\":" << (!result.correlationNA && result.finite ? jsonNumber(result.correlation) : "null")
                  << ",\"correlation_na\":" << (result.correlationNA ? "true" : "false")
                  << ",\"reference_is_zero\":" << (result.zeroReference ? "true" : "false")
                  << ",\"overall_pass\":" << (result.passed ? "true" : "false") << "}\n";
                success &= result.passed; ++count;
            }
            auto f = newOutput(activeOut / "verification_summary.json"); success &= count == 25;
            f << "{\"method\":" << jsonString(method) << ",\"fixture_count\":" << count
              << ",\"expected_fixture_count\":25,\"all_passed\":" << (success ? "true" : "false")
              << ",\"equivalence_to_full_rate_FDM_claimed\":false,\"reference_complex_l2_limit\":0.0001,\"zero_max_abs_limit\":1e-7}\n";
        } else if (mode == "latency" || mode == "memory") {
            rows.reserve(kWarmups + kMeasurements + 1);
            for (int i = 0; i < kWarmups && success; ++i)
                success = runFrame(*processor, frames[i % kFixtureCount], image, rows, "warmup", i, stopFile);
            for (int i = 0; i < kMeasurements && success; ++i)
                success = runFrame(*processor, frames[i % kFixtureCount], image, rows,
                                   mode == "latency" ? "measure" : "memory_probe", i, stopFile);
            if (!rows.empty()) { lastId = rows.back().fixture; haveImage = rows.back().status != "compute_error"; }
            if (mode == "memory" && success) {
                success = runFrame(*processor, frames.back(), image, rows, "untimed_memory_snapshot", 0, stopFile);
                lastId = frames.back().id; haveImage = !rows.empty() && rows.back().status != "compute_error";
                memoryReport = processor->memory(); liveHost = hostSnapshot();
                event("memory_snapshot", liveHost.ns);
            }
            writeRows(activeOut / "latency_raw.csv", rows, method, session);
        } else {
            rows.reserve(60000);
            const auto idleBegin = nowNs(); event("idle_begin", idleBegin);
            while (nowNs() - idleBegin < 30 * kSecond && !stopRequested(stopFile))
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            if (stopRequested(stopFile)) success = false;
            const auto warmupBegin = nowNs(); event("warmup_begin", warmupBegin); int warmupCount = 0;
            while (success && nowNs() - warmupBegin < 60 * kSecond) {
                success = runFrame(*processor, frames[warmupCount % kFixtureCount], image, rows, "power_warmup", warmupCount, stopFile);
                if (success) ++warmupCount;
            }
            const auto t0 = nowNs();
            event("window_start", t0, ",\"t0_ns\":" + std::to_string(t0) + ",\"target_submission_ns\":" + std::to_string(120 * kSecond));
            int attempted = 0, completed = 0;
            while (success && nowNs() - t0 < 120 * kSecond) {
                const auto previousRows = rows.size();
                success = runFrame(*processor, frames[attempted % kFixtureCount], image, rows, "power_measure", attempted, stopFile);
                if (rows.size() > previousRows) ++attempted;
                if (success) ++completed;
            }
            const auto t1 = nowNs();
            event("window_end", t1, ",\"t0_ns\":" + std::to_string(t0) + ",\"t1_ns\":" + std::to_string(t1)
                  + ",\"completed_blocks\":" + std::to_string(completed) + ",\"workload_valid\":" + (success && completed > 0 ? "true" : "false"));
            if (!rows.empty()) { lastId = rows.back().fixture; haveImage = rows.back().status != "compute_error"; }
            writeRows(activeOut / "power_blocks.csv", rows, method, session);
            auto f = newOutput(activeOut / "window.json");
            f << "{\"method\":" << jsonString(method) << ",\"session\":" << session << ",\"t0_ns\":" << t0 << ",\"t1_ns\":" << t1
              << ",\"target_submission_ns\":" << 120 * kSecond << ",\"actual_window_ns\":" << t1 - t0
              << ",\"warmup_completed_blocks\":" << warmupCount << ",\"attempted_blocks\":" << attempted << ",\"completed_blocks\":" << completed
              << ",\"workload_valid\":" << (success && completed > 0 ? "true" : "false")
              << ",\"includes_light_checks_and_stop_file_polling\":true,\"includes_final_engine_destruction\":false,\"power_integration_performed_by_parent_collector\":true}\n";
        }
        if (mode != "verify") {
            if (mode != "memory" || !success) { memoryReport = processor->memory(); liveHost = hostSnapshot(); }
            writeMemory(activeOut / "memory.json", memoryReport, contextHost, liveHost, implementation,
                        success && haveImage ? "untimed_live_engine_after_complete_frame" : "failed_or_incomplete_live_engine_snapshot");
            bufferTable(activeOut / "buffers.csv", base, fdm.get(), memoryReport);
            if (haveImage) writeBinary(activeOut / ("last_" + lastId + ".bin"), image);
        }
        const auto finalizeBegin = nowNs(); processor.reset(); const auto finalizeReady = nowNs();
        auto lifecycle = newOutput(activeOut / "lifecycle.json");
        lifecycle << "{\"method\":" << jsonString(method) << ",\"mode\":" << jsonString(mode) << ",\"session\":" << session
                  << ",\"implementation\":" << jsonString(implementation)
                  << ",\"context_begin_ns\":" << contextBegin << ",\"context_ready_ns\":" << contextReady
                  << ",\"initialize_begin_ns\":" << initializeBegin << ",\"initialize_ready_ns\":" << initializeReady
                  << ",\"finalize_begin_ns\":" << finalizeBegin << ",\"finalize_ready_ns\":" << finalizeReady
                  << ",\"preloaded_noise_fixtures\":" << frames.size() << ",\"warmups_latency_memory\":50,\"measurements_latency_memory\":200"
                  << ",\"main_frame_wall_endpoint\":\"process_return_after_complex_D2H_and_stream_sync\""
                  << ",\"light_check_outside_frame_wall\":true,\"success\":" << (success ? "true" : "false") << "}\n";
        event("completed", nowNs(), ",\"success\":" + std::string(success ? "true" : "false"));
        return success ? 0 : 2;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        if (ownsOut && !fs::exists(activeOut / "failure.json")) {
            try { auto f = newOutput(activeOut / "failure.json"); f << "{\"failed\":true,\"message\":" << jsonString(e.what()) << ",\"monotonic_ns\":" << nowNs() << "}\n"; }
            catch (...) {}
        }
        return 1;
    }
}
