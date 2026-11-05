// stl_raytrace_arrivals_cuda.cu (amplitude semantics, glint-suppressed, energy-conserving BRDF)
// -----------------------------------------------------------------------------
// MEX-CUDA: STL 场景递归光线追踪（PRT + WRT，至多 3 次弹跳）→ “等效到达散射点”
//
// 传播（幅度语义):
//   A ∝ [1/(4π r_tx)] * [1/(4π r_rx)] * exp(-alpha_amp * L_total) * Gain_mat * BRDF_amp * Gtx * Grx
// 若 amp_type=1 输出为强度：I = A^2。
// BRDF 可选：
//   brdf_model=0 : 旧版（镜面 cos^m + Kirchhoff ρ + 漫反）幅度混合（非能量归一）
//   brdf_model=1 : 能量守恒（功率归一的 Phong + Lambert，经 ρ 分配），幅度=√(功率BRDF)
// 末端散射面可启用法向“坡度抖动”；可限制 BRDF 幅度上限；可给 Tx/Rx 小圆锥抖动。
// pos 为真实命中点；多跳“等效两段差”只在 extra_delay 给出；不做坐标外推。
// 方向性：p<0 表示各向同性（G=1）；p=0 半空间门限；p>0 余弦^p。
// -----------------------------------------------------------------------------

#include "mex.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdint.h>
#include <string>
#include <cuda_runtime.h>

#define CHECK_CUDA(stmt) do{ cudaError_t e=(stmt); if(e!=cudaSuccess){ mexErrMsgIdAndTxt("stl:cuda", cudaGetErrorString(e)); }}while(0)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ---------------------------- 基础向量 ----------------------------
struct float3x { float x,y,z; __host__ __device__ float3x(){} __host__ __device__ float3x(float a,float b,float c):x(a),y(b),z(c){} };
__host__ __device__ inline float3x make_f3(float a,float b,float c){ return float3x(a,b,c); }
__host__ __device__ inline float  dot3(const float3x& a,const float3x& b){ return a.x*b.x + a.y*b.y + a.z*b.z; }
__host__ __device__ inline float3x cross3(const float3x& a,const float3x& b){ return make_f3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); }
__host__ __device__ inline float  norm3(const float3x& a){ return sqrtf(dot3(a,a)); }
__host__ __device__ inline float3x add3(const float3x& a,const float3x& b){ return make_f3(a.x+b.x,a.y+b.y,a.z+b.z); }
__host__ __device__ inline float3x sub3(const float3x& a,const float3x& b){ return make_f3(a.x-b.x,a.y-b.y,a.z-b.z); }
__host__ __device__ inline float3x mul3(const float3x& a,float s){ return make_f3(a.x*s,a.y*s,a.z*s); }
__host__ __device__ inline float3x normalize3(const float3x& a){ float n = norm3(a); return (n>0.f)? make_f3(a.x/n,a.y/n,a.z/n) : make_f3(0,0,0); }

// RNG: LCG + uniform + normal（Box-Muller）
__host__ __device__ inline uint32_t lcg(uint64_t &state){ state = state*6364136223846793005ULL + 1ULL; return (uint32_t)(state >> 32); }
__host__ __device__ inline float randu(uint64_t &state){ return ( (lcg(state) & 0x00FFFFFF) / float(0x01000000) ); }
__host__ __device__ inline float randn(uint64_t& state){
    float u1 = fmaxf(randu(state), 1e-7f);
    float u2 = fmaxf(randu(state), 1e-7f);
    float r = sqrtf(-2.0f * logf(u1));
    float th= 2.0f*(float)M_PI * u2;
    return r * cosf(th);
}

// 射线-三角交 Möller–Trumbore
__host__ __device__ inline bool rayTriMT(const float3x& orig, const float3x& dir,
                                         const float3x& v0, const float3x& v1, const float3x& v2,
                                         float& t, float& u, float& v, float eps=1e-6f){
    float3x e1 = sub3(v1,v0), e2 = sub3(v2,v0);
    float3x p  = cross3(dir, e2);
    float det  = dot3(e1, p);
    if (fabsf(det) < eps) return false;
    float inv = 1.0f / det;
    float3x s = sub3(orig, v0);
    u = dot3(s, p) * inv; if (u < -eps || u > 1.f+eps) return false;
    float3x q = cross3(s, e1);
    v = dot3(dir, q) * inv; if (v < -eps || u+v > 1.f+eps) return false;
    t = dot3(e2, q) * inv;
    return t > eps;
}

// 反射向量
__host__ __device__ inline float3x reflect(const float3x& v_in, const float3x& n_raw){
    float3x n = n_raw; if (dot3(v_in, n) > 0.0f) n = mul3(n, -1.0f);
    return sub3(v_in, mul3(n, 2.0f*dot3(v_in,n)));
}

// Orthonormal basis
__host__ __device__ inline void orthonormal_basis(const float3x& n, float3x& t1, float3x& t2){
    if (fabsf(n.z) < 0.999f) t1 = normalize3( cross3(make_f3(0,0,1), n) );
    else                     t1 = normalize3( cross3(make_f3(0,1,0), n) );
    t2 = cross3(n, t1);
}

// 法向坡度抖动（仅用于 BRDF 评估，不改几何）
__host__ __device__ inline float3x jitter_normal(const float3x& n, float slope_std, uint64_t& st){
    if (slope_std <= 1e-9f) return n;
    float3x t1,t2; orthonormal_basis(n, t1, t2);
    float s1 = slope_std * randn(st);
    float s2 = slope_std * randn(st);
    float3x nj = normalize3( add3( n, add3(mul3(t1,s1), mul3(t2,s2)) ) );
    return nj;
}

// 圆锥方向抖动（半角 angle_rad）
__host__ __device__ inline float3x cone_dither(const float3x& dir, float angle_rad, uint64_t& st){
    if (angle_rad <= 1e-7f) return dir;
    float3x z = normalize3(dir), x, y; orthonormal_basis(z, x, y);
    float u = randu(st), v = randu(st);
    float cos_t = 1.0f - u * (1.0f - cosf(angle_rad));
    float sin_t = sqrtf(fmaxf(0.f, 1.0f - cos_t*cos_t));
    float phi   = 2.0f*(float)M_PI * v;
    float3x local = add3( add3( mul3(x, cosf(phi)*sin_t ), mul3(y, sinf(phi)*sin_t) ), mul3(z, cos_t) );
    return normalize3(local);
}

// Thorp：f[kHz] → alpha_amp [Np/m]（幅度域）
__host__ __device__ inline float thorp_alpha_np_per_m(float f_khz){
    float f2=f_khz*f_khz;
    float a_db_km = 0.11f*(f2/(1.0f+f2)) + 44.0f*(f2/(4100.0f+f2)) + 2.75e-4f*f2 + 0.003f;
    float a_db_m  = a_db_km * 1e-3f;
    return a_db_m / 8.6858896f; // dB/m → Np/m（幅度）
}

// ---------------------------- 场景结构 ----------------------------
struct Tri {
    float3x v0,v1,v2;
    float3x n;
    float   refl;   // reflectivity（语义由 reflect_is_power 决定）
    float   scat;   // scattering_coeff (0~1)
    float   f0_hz;  // center_frequency (Hz)
    float   trans;  // 预留
};

struct ParamsGPU {
    // 声学与几何
    float3x S;  float c;  float eps;

    // 弹跳/采样
    int   max_bounces;
    int   samples_per_tri;
    unsigned long long seed;

    // Rx 指向 + 抖动
    int   rx_has_per_channel;
    float3x rx_dir_common;
    int   rx_cos_power;   // <0: 各向同性；=0: 半空间门限；>0: cos^p
    float rx_dither_rad;

    // Tx 指向 + 抖动
    float3x tx_dir;
    int   tx_cos_power;   // <0: 各向同性；=0: 半空间门限；>0: cos^p
    float tx_dither_rad;

    // BRDF / 语义
    int   brdf_model;       // 0=旧版幅度混合, 1=能量守恒(推荐)
    float brdf_cos_power;   // 镜面指数 m
    float rough_sigma;      // 粗糙 rms 高度 σ[m]
    float brdf_amp_cap;     // >0 则对 BRDF 幅度做上限
    int   amp_type;         // 0=幅度, 1=强度
    int   reflect_is_power; // 1=功率系数, 0=幅度系数

    // 微表面法向抖动
    int   normal_jitter_on;
    float normal_jitter_slope; // rms 坡度

    // 多跳附加衰减（幅度域）
    float order2_extra;
    float order3_extra;
};

// 输出结构
struct Arrival { float3x pos; float amp; float extra_dt; float path_len; int order; int rx_idx; };

// 线段遮挡
__device__ bool isSegmentOccluded(const float3x& orig, const float3x& dest,
                                  const Tri* tris, int nTri,
                                  int tri_ignore_a, int tri_ignore_b, float eps){
    float3x dir = sub3(dest, orig);
    float   L   = norm3(dir);
    if (L <= eps) return false;
    dir = mul3(dir, 1.0f/L);
    for(int i=0;i<nTri;i++){
        if (i==tri_ignore_a || i==tri_ignore_b) continue;
        float t,u,v;
        if (rayTriMT(orig, dir, tris[i].v0, tris[i].v1, tris[i].v2, t,u,v, eps)){
            if (t > eps && t < (L - eps)) return true;
        }
    }
    return false;
}

// 指向性（p<0 → 各向同性：恒为 1）
__host__ __device__ inline float gain_cos_power(float cosang, int p){
    if (p < 0) return 1.f;                     // 各向同性
    if (p == 0) return (cosang>0.f)? 1.f:0.f;  // 半空间门限
    float ct = fmaxf(0.0f, cosang);
    return powf(ct, (float)p);                 // 余弦^p
}

// BRDF: 旧版（幅度混合）
__host__ __device__ inline float brdf_old_amp(
    const float3x& n_eff, const float3x& S, const float3x& P,
    float S_coeff, float cos_power_m, float rough_sigma, float lambda)
{
    float3x toL = normalize3(sub3(S, P));
    float cosT  = fmaxf(0.0f, dot3(n_eff, toL));
    float spec  = (cos_power_m>0.f) ? powf(cosT, cos_power_m) : cosT;
    float a     = 4.0f * (float)M_PI * rough_sigma * cosT / fmaxf(lambda, 1e-9f);
    float rho   = expf(-a*a);
    float brdf  = sqrtf(fmaxf(1.0f - S_coeff, 0.0f)) * (rho * spec)
                + sqrtf(fmaxf(S_coeff, 0.0f))       * ((1.0f - rho) * cosT);
    return brdf;
}

// BRDF: 能量守恒（功率归一的 Phong + Lambert，经 ρ 分配），幅度=√(功率BRDF)
__host__ __device__ inline float brdf_energy_conserving_amp(
    const float3x& n_eff, const float3x& S, const float3x& P,
    float S_coeff, float cos_power_m, float rough_sigma, float lambda)
{
    float3x toL   = normalize3(sub3(S, P));
    float cosT    = fmaxf(0.0f, dot3(n_eff, toL));
    float a       = 4.0f * (float)M_PI * rough_sigma * cosT / fmaxf(lambda, 1e-9f);
    float rho     = expf(-a*a);                       // 相干保留
    float specPow = ((cos_power_m + 2.0f) / (2.0f * (float)M_PI)) * powf(cosT, cos_power_m); // 归一化 Phong 功率
    float diffPow = (1.0f/(float)M_PI) * cosT;        // Lambert 功率
    float Pbrdf   = fmaxf(0.f,
                     (1.0f - S_coeff) * rho * specPow
                   +    S_coeff       * (1.0f - rho) * diffPow );
    return sqrtf(Pbrdf); // 幅度
}

__host__ __device__ inline float eval_brdf_amp(
    int brdf_model, const float3x& n_eff, const float3x& S, const float3x& P,
    float S_coeff, float cos_power_m, float rough_sigma, float lambda)
{
    if (brdf_model==1) return brdf_energy_conserving_amp(n_eff,S,P,S_coeff,cos_power_m,rough_sigma,lambda);
    return brdf_old_amp(n_eff,S,P,S_coeff,cos_power_m,rough_sigma,lambda);
}

// -------------------- 关键修改：传播“幅度”语义 --------------------
__host__ __device__ inline float two_way_amp(float r_tx, float r_rx, float L_total, float f0_hz){
    // const float inv4pi = 1.0f / (4.0f * (float)M_PI);
    const float inv4pi = 1.0f;
    // 几何（幅度）：两程相乘 → 1/(4π r_tx) * 1/(4π r_rx)
    float geom_amp = powf((inv4pi / fmaxf(r_tx, 1e-6f)) * (inv4pi / fmaxf(r_rx, 1e-6f)),0.65f);
    // 介质（幅度）：alpha_amp[Np/m]，幅度单程 e^{-alpha L}（这里 L_total 为两程总长）
    float alpha_amp = thorp_alpha_np_per_m(f0_hz * 1e-3f);
    float med_amp   = expf(-alpha_amp * fmaxf(L_total, 0.0f));
    return geom_amp * med_amp;
}

__host__ __device__ inline float reflect_to_amp(float refl, int reflect_is_power){
    if (reflect_is_power) return sqrtf(fmaxf(refl, 0.0f));
    return refl;
}

// ---------------------------- Kernel ----------------------------
__global__ void kernelTraceArrivals(
    const Tri* __restrict__ tris, int nTri,
    const float3x* __restrict__ R, const float3x* __restrict__ rx_dirs, int nR,
    const ParamsGPU prm,
    Arrival* __restrict__ out, int max_out, int* __restrict__ counter)
{
    long long total = (long long)nTri * (long long)prm.samples_per_tri * (long long)nR;
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    if ((long long)tid >= total) return;

    int triIdx   = tid / (prm.samples_per_tri * nR);
    int rem      = tid % (prm.samples_per_tri * nR);
    int sampId   = rem / nR;
    int rIdx     = rem % nR;

    const Tri &T = tris[triIdx];

    // RNG 状态（样本稳定随机）
    uint64_t st = prm.seed
        ^ (0x9e3779b97f4a7c15ULL * ((unsigned long long)triIdx*1315423911ULL
                                  + (unsigned long long)sampId*2654435761ULL
                                  + (unsigned long long)rIdx*8895235929ULL));

    // Rx 方向（可抖动）
    float3x rxDir = prm.rx_has_per_channel ? rx_dirs[rIdx] : prm.rx_dir_common;
    rxDir = normalize3(rxDir);
    if (prm.rx_dither_rad > 1e-7f) rxDir = cone_dither(rxDir, prm.rx_dither_rad, st);

    // Tx 方向（可抖动）
    float3x txDir = normalize3(prm.tx_dir);
    if (prm.tx_dither_rad > 1e-7f) txDir = cone_dither(txDir, prm.tx_dither_rad, st);

    // 样本点：质心 or 重心随机
    float3x P;
    if (prm.samples_per_tri <= 1 || sampId==0){
        P = mul3(add3(add3(T.v0, T.v1), T.v2), 1.0f/3.0f);
    } else {
        float u = randu(st), v = randu(st); if (u+v>1.f){ u=1.f-u; v=1.f-v; }
        float w = 1.f - u - v;
        P = add3( add3( mul3(T.v0, w), mul3(T.v1, u) ), mul3(T.v2, v) );
    }

    // 末端法向（可做微表面坡度抖动）
    float R_amp_P = reflect_to_amp(T.refl, prm.reflect_is_power);
    float lambdaP = prm.c / fmaxf(T.f0_hz, 1e-6f);

    // ---------- 一次弹跳（PRT 背散） ----------
    bool occSP = isSegmentOccluded(prm.S, P, tris, nTri, triIdx, -1, prm.eps);
    if (!occSP){
        bool occPR = isSegmentOccluded(P, R[rIdx], tris, nTri, triIdx, -1, prm.eps);
        if (!occPR){
            float3x dirSP = normalize3(sub3(P,prm.S));    // S->P
            float3x dirRP = normalize3(sub3(P, R[rIdx])); // R->P

            float Grx = gain_cos_power( dot3(rxDir, dirRP), prm.rx_cos_power );
            if (Grx > 0.f){
                float Gtx = gain_cos_power( dot3(txDir, dirSP), prm.tx_cos_power );

                float rSP = norm3(sub3(P,prm.S));
                float rPR = norm3(sub3(R[rIdx],P));
                float L1  = rSP + rPR;

                float3x n_eff = T.n;
                if (prm.normal_jitter_on) n_eff = jitter_normal(n_eff, prm.normal_jitter_slope, st);

                float brdf = eval_brdf_amp(prm.brdf_model, n_eff, prm.S, P,
                                           T.scat, prm.brdf_cos_power, prm.rough_sigma, lambdaP);
                if (prm.brdf_amp_cap > 1e-9f) brdf = fminf(brdf, prm.brdf_amp_cap);

                float base = two_way_amp(rSP, rPR, L1, T.f0_hz);   // 幅度
                float amp  = R_amp_P * brdf * base * Gtx * Grx;    // 幅度

                if (prm.amp_type==1) amp = amp*amp; // 强度（I=A^2）

                int idx = atomicAdd(counter, 1);
                if (idx < max_out){
                    out[idx].pos      = P;
                    out[idx].amp      = amp;
                    out[idx].extra_dt = 0.0f;
                    out[idx].path_len = L1;
                    out[idx].order    = 1;
                    out[idx].rx_idx   = rIdx + 1;
                }
            }
        }
    }

    if (prm.max_bounces < 2) return;

    // ---------- 二次弹跳（P 镜面 → Q 背散） ----------
    float3x dir_in = normalize3(sub3(P, prm.S));   // S->P
    float3x vref   = reflect(dir_in, T.n);         // 在 P 处镜面
    float tMin=1e30f; int tri2=-1; float tt=0,uu=0,vv=0;
    for(int i=0;i<nTri;i++){
        if (i==triIdx) continue;
        float t,u,v;
        if(rayTriMT(P, vref, tris[i].v0,tris[i].v1,tris[i].v2, t,u,v, prm.eps)){
            if (t< tMin){ tMin=t; tri2=i; tt=t; uu=u; vv=v; }
        }
    }
    if (tri2>=0){
        const Tri &T2 = tris[tri2];
        float3x Q = add3(P, mul3(vref, tt));
        bool occPQ = isSegmentOccluded(P, Q, tris, nTri, triIdx, tri2, prm.eps);
        if (!occPQ){
            bool occQR = isSegmentOccluded(Q, R[rIdx], tris, nTri, tri2, -1, prm.eps);
            if (!occQR){
                float3x dirRQ = normalize3(sub3(Q, R[rIdx]));   // R->Q
                float3x dirSQ = normalize3(sub3(Q, prm.S));     // S->Q
                float Grx = gain_cos_power( dot3(rxDir, dirRQ), prm.rx_cos_power );
                if (Grx > 0.f){
                    float Gtx = gain_cos_power( dot3(txDir, dirSQ), prm.tx_cos_power );

                    float rSP = norm3(sub3(P,prm.S));
                    float rPQ = norm3(sub3(Q,P));
                    float rQR = norm3(sub3(R[rIdx],Q));
                    float L2  = rSP + rPQ + rQR;

                    float3x n2 = T2.n;
                    if (prm.normal_jitter_on) n2 = jitter_normal(n2, prm.normal_jitter_slope, st);

                    float lambdaQ = prm.c / fmaxf(T2.f0_hz, 1e-6f);
                    float brdfQ   = eval_brdf_amp(prm.brdf_model, n2, prm.S, Q,
                                                  T2.scat, prm.brdf_cos_power, prm.rough_sigma, lambdaQ);
                    if (prm.brdf_amp_cap > 1e-9f) brdfQ = fminf(brdfQ, prm.brdf_amp_cap);

                    // 注意 two_way_amp 仍按“直连”SQ/RQ 的两段距离给几何，但吸收用全 L2
                    float base    = two_way_amp(norm3(sub3(Q,prm.S)), norm3(sub3(R[rIdx],Q)), L2, T2.f0_hz);
                    float amp2    = R_amp_P * brdfQ * base * Gtx * Grx * prm.order2_extra; // 幅度
                    if (prm.amp_type==1) amp2 = amp2*amp2;

                    float extra = ( L2 - ( norm3(sub3(Q,prm.S)) + norm3(sub3(R[rIdx],Q)) ) ) / prm.c;

                    int idx2 = atomicAdd(counter, 1);
                    if (idx2 < max_out){
                        out[idx2].pos      = Q;
                        out[idx2].amp      = amp2;
                        out[idx2].extra_dt = extra;
                        out[idx2].path_len = L2;
                        out[idx2].order    = 2;
                        out[idx2].rx_idx   = rIdx + 1;
                    }

                    if (prm.max_bounces >= 3){
                        // ---------- 三次弹跳：P 镜面 → Q 镜面 → W 背散 ----------
                        float3x dir_in2 = normalize3(sub3(Q, P));
                        float3x vref2   = reflect(dir_in2, T2.n);
                        float tMin2=1e30f; int tri3=-1; float t2=0,u2=0,v2=0;
                        for(int j=0;j<nTri;j++){
                            if (j==tri2) continue;
                            float t,u,v;
                            if(rayTriMT(Q, vref2, tris[j].v0,tris[j].v1,tris[j].v2, t,u,v, prm.eps)){
                                if (t< tMin2){ tMin2=t; tri3=j; t2=t; u2=u; v2=v; }
                            }
                        }
                        if (tri3>=0){
                            const Tri &T3 = tris[tri3];
                            float3x W = add3(Q, mul3(vref2, t2));
                            bool occQW = isSegmentOccluded(Q, W, tris, nTri, tri2, tri3, prm.eps);
                            if (!occQW){
                                bool occWR = isSegmentOccluded(W, R[rIdx], tris, nTri, tri3, -1, prm.eps);
                                if (!occWR){
                                    float3x dirRW = normalize3(sub3(W, R[rIdx])); // R->W
                                    float3x dirSW = normalize3(sub3(W, prm.S));   // S->W
                                    float Grx3 = gain_cos_power( dot3(rxDir, dirRW), prm.rx_cos_power );
                                    if (Grx3 > 0.f){
                                        float Gtx3 = gain_cos_power( dot3(txDir, dirSW), prm.tx_cos_power );

                                        float rSP = norm3(sub3(P,prm.S));
                                        float rPQ = norm3(sub3(Q,P));
                                        float rQW = norm3(sub3(W,Q));
                                        float rWR = norm3(sub3(R[rIdx],W));
                                        float L3  = rSP + rPQ + rQW + rWR;

                                        float3x n3 = T3.n;
                                        if (prm.normal_jitter_on) n3 = jitter_normal(n3, prm.normal_jitter_slope, st);

                                        float lambdaW = prm.c / fmaxf(T3.f0_hz, 1e-6f);
                                        float brdfW   = eval_brdf_amp(prm.brdf_model, n3, prm.S, W,
                                                                      T3.scat, prm.brdf_cos_power, prm.rough_sigma, lambdaW);
                                        if (prm.brdf_amp_cap > 1e-9f) brdfW = fminf(brdfW, prm.brdf_amp_cap);

                                        float base3   = two_way_amp(norm3(sub3(W,prm.S)), norm3(sub3(R[rIdx],W)), L3, T3.f0_hz);
                                        float R_amp_2 = reflect_to_amp(T2.refl, prm.reflect_is_power);
                                        float amp3    = R_amp_P * R_amp_2 * brdfW * base3 * Gtx3 * Grx3 * prm.order3_extra; // 幅度
                                        if (prm.amp_type==1) amp3 = amp3*amp3;

                                        float extra3 = ( L3 - ( norm3(sub3(W,prm.S)) + norm3(sub3(R[rIdx],W)) ) ) / prm.c;

                                        int idx3 = atomicAdd(counter, 1);
                                        if (idx3 < max_out){
                                            out[idx3].pos      = W;
                                            out[idx3].amp      = amp3;
                                            out[idx3].extra_dt = extra3;
                                            out[idx3].path_len = L3;
                                            out[idx3].order    = 3;
                                            out[idx3].rx_idx   = rIdx + 1;
                                        }
                                    }
                                }
                            }
                        }
                    } // max_bounces>=3
                }
            }
        }
    }
}

// ---------------------------- Host 侧工具 ----------------------------
static double getFieldOrDefaultDouble(const mxArray* S, const char* name, double defv){
    if (!S) return defv;
    mxArray* f = mxGetField(S, 0, name);
    if(!f){ std::string low(name); std::transform(low.begin(),low.end(),low.begin(),::tolower); f=mxGetField(S,0,low.c_str()); }
    return f? mxGetScalar(f):defv;
}
static int getFieldOrDefaultInt(const mxArray* S, const char* name, int defv){
    if (!S) return defv;
    mxArray* f = mxGetField(S, 0, name);
    if(!f){ std::string low(name); std::transform(low.begin(),low.end(),low.begin(),::tolower); f=mxGetField(S,0,low.c_str()); }
    return f? (int)mxGetScalar(f):defv;
}
static unsigned long long getFieldOrDefaultULL(const mxArray* S, const char* name, unsigned long long defv){
    if (!S) return defv;
    mxArray* f = mxGetField(S, 0, name);
    if(!f){ std::string low(name); std::transform(low.begin(),low.end(),low.begin(),::tolower); f=mxGetField(S,0,low.c_str()); }
    return f? (unsigned long long)mxGetScalar(f):defv;
}
static void getVec3OrDefault(const mxArray* S, const char* name, double defx, double defy, double defz, double* out3){
    if (!S){ out3[0]=defx; out3[1]=defy; out3[2]=defz; return; }
    mxArray* f = mxGetField(S,0,name);
    if(!f){ std::string low(name); std::transform(low.begin(),low.end(),low.begin(),::tolower); f=mxGetField(S,0,low.c_str()); }
    if (!f || mxIsEmpty(f)){ out3[0]=defx; out3[1]=defy; out3[2]=defz; return; }
    if (mxGetNumberOfElements(f)==3){ double* p = mxGetPr(f); out3[0]=p[0]; out3[1]=p[1]; out3[2]=p[2]; }
    else { out3[0]=defx; out3[1]=defy; out3[2]=defz; }
}

// ---------------------------- MEX 入口 ----------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if (nrhs < 4) mexErrMsgIdAndTxt("stl:args", "Usage: out = stl_raytrace_arrivals_cuda(models_cell, S, R, params)");

    const mxArray* models_cell = prhs[0];
    const mxArray* S_in        = prhs[1];
    const mxArray* R_in        = prhs[2];
    const mxArray* params_in   = prhs[3];

    if (!mxIsCell(models_cell)) mexErrMsgIdAndTxt("stl:args", "models_cell must be a cell array");
    if (mxGetM(S_in)*mxGetN(S_in)!=3) mexErrMsgIdAndTxt("stl:args", "S must be [1x3]");
    if (mxGetN(R_in)!=3) mexErrMsgIdAndTxt("stl:args", "R must be [Nr x 3]");

    // S, R
    double* Sd = mxGetPr(S_in);
    float3x S = make_f3((float)Sd[0], (float)Sd[1], (float)Sd[2]);

    size_t nR = mxGetM(R_in); double* Rd = mxGetPr(R_in);
    std::vector<float3x> h_R(nR);
    for(size_t i=0;i<nR;i++) h_R[i] = make_f3((float)Rd[i], (float)Rd[i+nR], (float)Rd[i+2*nR]);

    // params
    ParamsGPU prm;
    prm.S   = S;
    prm.c   = (float)getFieldOrDefaultDouble(params_in, "c", 1500.0);
    prm.eps = (float)getFieldOrDefaultDouble(params_in, "eps", 1e-5);
    prm.seed= getFieldOrDefaultULL(params_in, "seed", 1234567ULL);

    prm.max_bounces     = std::min(3, std::max(1, getFieldOrDefaultInt(params_in, "max_bounces", 3)));
    prm.samples_per_tri = std::max(1, getFieldOrDefaultInt(params_in, "samples_per_tri", 1));

    // Rx 指向 + 抖动
    double rxdir0[3]; getVec3OrDefault(params_in, "rx_dir", 0,0,1, rxdir0);
    prm.rx_dir_common = normalize3( make_f3((float)rxdir0[0],(float)rxdir0[1],(float)rxdir0[2]) );
    prm.rx_cos_power  = getFieldOrDefaultInt(params_in, "rx_cos_power", 0); // -1 → 各向性
    double rx_dith_deg= getFieldOrDefaultDouble(params_in, "rx_dither_deg", 0.0);
    prm.rx_dither_rad = (float)(rx_dith_deg * M_PI / 180.0);

    // Tx 指向 + 抖动
    double txdir0[3]; getVec3OrDefault(params_in, "tx_dir", 0,0,1, txdir0);
    prm.tx_dir       = normalize3( make_f3((float)txdir0[0],(float)txdir0[1],(float)txdir0[2]) );
    prm.tx_cos_power = getFieldOrDefaultInt(params_in, "tx_cos_power", 0); // -1 → 各向性
    double tx_dith_deg= getFieldOrDefaultDouble(params_in, "tx_dither_deg", 0.0);
    prm.tx_dither_rad = (float)(tx_dith_deg * M_PI / 180.0);

    // BRDF / 语义
    prm.brdf_model      = getFieldOrDefaultInt(params_in, "brdf_model", 1); // 默认能量守恒
    prm.brdf_cos_power  = (float)getFieldOrDefaultDouble(params_in, "brdf_cos_power", 8.0);
    prm.rough_sigma     = (float)getFieldOrDefaultDouble(params_in, "rough_sigma", 0.0);
    prm.brdf_amp_cap    = (float)getFieldOrDefaultDouble(params_in, "brdf_amp_cap", 0.0); // 0=禁用
    prm.amp_type        = getFieldOrDefaultInt(params_in, "amp_type", 0);
    prm.reflect_is_power= getFieldOrDefaultInt(params_in, "reflect_is_power", 1);

    // 法向抖动
    prm.normal_jitter_on    = getFieldOrDefaultInt(params_in, "normal_jitter_on", 0);
    prm.normal_jitter_slope = (float)getFieldOrDefaultDouble(params_in, "normal_jitter_slope", 0.0);

    // 多跳附加衰减（幅度）
    prm.order2_extra     = (float)getFieldOrDefaultDouble(params_in, "order2_extra", 1.0);
    prm.order3_extra     = (float)getFieldOrDefaultDouble(params_in, "order3_extra", 1.0);

    // rx_dir 每通道
    std::vector<float3x> h_rxdir(nR);
    prm.rx_has_per_channel = 0;
    if (params_in){
        mxArray* f = mxGetField(params_in, 0, "rx_dir"); if(!f){ f=mxGetField(params_in,0,"RX_DIR"); }
        if (f && !mxIsEmpty(f)){
            if (mxGetNumberOfElements(f)==3){
                double* p = mxGetPr(f);
                prm.rx_dir_common = normalize3( make_f3((float)p[0],(float)p[1],(float)p[2]) );
                prm.rx_has_per_channel = 0;
            } else if (mxGetN(f)==3 && mxGetM(f)==nR){
                double* p = mxGetPr(f);
                for(size_t i=0;i<nR;i++) h_rxdir[i] = normalize3( make_f3((float)p[i], (float)p[i+nR], (float)p[i+2*nR]) );
                prm.rx_has_per_channel = 1;
            }
        }
    }

    // 展平三角网格（读取每模型参数）
    size_t nModels = mxGetNumberOfElements(models_cell);
    std::vector<Tri> h_tris; h_tris.reserve(1024);
    for(size_t m=0;m<nModels;m++){
        const mxArray* mdl = mxGetCell(models_cell, (mwIndex)m); if (!mdl) continue;
        auto field=[&](const char* name)->const mxArray*{
            const mxArray* f=mxGetField(mdl,0,name);
            if(!f){ std::string low(name); std::transform(low.begin(),low.end(),low.begin(),::tolower); f=mxGetField(mdl,0,low.c_str()); }
            return f;
        };
        mxArray* P = (mxArray*)field("Points");
        mxArray* C = (mxArray*)field("ConnectivityList");
        if (!P || !C) mexErrMsgIdAndTxt("stl:args", "Each model must have Points and ConnectivityList");

        double refl = field("reflectivity")? mxGetScalar(field("reflectivity")) : 0.7;
        double scat = field("scattering_coeff")? mxGetScalar(field("scattering_coeff")) : 0.5;
        double f0hz = field("center_frequency")? mxGetScalar(field("center_frequency")) : 450e3;
        double trans= field("transmittance")? mxGetScalar(field("transmittance")) : 1.0;

        size_t nP = mxGetM(P); if (mxGetN(P)!=3) mexErrMsgIdAndTxt("stl:args", "Points must be [N x 3]");
        double* Pd = mxGetPr(P);
        size_t nT = mxGetM(C); double* Cd = mxGetPr(C);

        for(size_t t=0;t<nT;t++){
            int i0=(int)Cd[t]-1, i1=(int)Cd[t+nT]-1, i2=(int)Cd[t+2*nT]-1;
            if (i0<0||i0>=(int)nP||i1<0||i1>=(int)nP||i2<0||i2>=(int)nP) mexErrMsgIdAndTxt("stl:args", "Connectivity index OOR");
            float3x v0 = make_f3((float)Pd[i0], (float)Pd[i0+nP], (float)Pd[i0+2*nP]);
            float3x v1 = make_f3((float)Pd[i1], (float)Pd[i1+nP], (float)Pd[i1+2*nP]);
            float3x v2 = make_f3((float)Pd[i2], (float)Pd[i2+nP], (float)Pd[i2+2*nP]);
            float3x n  = normalize3( cross3( sub3(v1,v0), sub3(v2,v0) ) );
            Tri T; T.v0=v0; T.v1=v1; T.v2=v2; T.n=n; T.refl=(float)refl; T.scat=(float)scat; T.f0_hz=(float)f0hz; T.trans=(float)trans;
            h_tris.push_back(T);
        }
    }
    int nTri = (int)h_tris.size(); if (nTri==0) mexErrMsgIdAndTxt("stl:args", "No triangles");

    int max_arrivals = getFieldOrDefaultInt(params_in, "max_arrivals",
        (int)(nR*(size_t)nTri*std::max(1,prm.max_bounces)/8 + 4096));
    max_arrivals = std::max(max_arrivals, 4096);

    // GPU 内存
    Tri* d_tris=nullptr; CHECK_CUDA(cudaMalloc(&d_tris, nTri*sizeof(Tri)));
    CHECK_CUDA(cudaMemcpy(d_tris, h_tris.data(), nTri*sizeof(Tri), cudaMemcpyHostToDevice));
    float3x* d_R=nullptr; CHECK_CUDA(cudaMalloc(&d_R, nR*sizeof(float3x)));
    CHECK_CUDA(cudaMemcpy(d_R, h_R.data(), nR*sizeof(float3x), cudaMemcpyHostToDevice));
    float3x* d_rxdir=nullptr;
    if (prm.rx_has_per_channel){
        CHECK_CUDA(cudaMalloc(&d_rxdir, nR*sizeof(float3x)));
        CHECK_CUDA(cudaMemcpy(d_rxdir, h_rxdir.data(), nR*sizeof(float3x), cudaMemcpyHostToDevice));
    }

    Arrival* d_out=nullptr; CHECK_CUDA(cudaMalloc(&d_out, max_arrivals*sizeof(Arrival)));
    int* d_count=nullptr;   CHECK_CUDA(cudaMalloc(&d_count, sizeof(int)));
    CHECK_CUDA(cudaMemset(d_count, 0, sizeof(int)));

    // Kernel
    long long total = (long long)nTri * (long long)prm.samples_per_tri * (long long)nR;
    int block=256, grid=(int)((total + block - 1)/block);
    kernelTraceArrivals<<<grid, block>>>(d_tris, nTri, d_R, (prm.rx_has_per_channel? d_rxdir:nullptr), (int)nR, prm, d_out, max_arrivals, d_count);
    CHECK_CUDA(cudaDeviceSynchronize());
    CHECK_CUDA(cudaGetLastError());

    int h_count=0; CHECK_CUDA(cudaMemcpy(&h_count, d_count, sizeof(int), cudaMemcpyDeviceToHost));
    h_count = std::min(h_count, max_arrivals);
    std::vector<Arrival> h_out(h_count);
    if (h_count>0) CHECK_CUDA(cudaMemcpy(h_out.data(), d_out, h_count*sizeof(Arrival), cudaMemcpyDeviceToHost));

    // 释放 GPU
    cudaFree(d_tris); cudaFree(d_R); if (d_rxdir) cudaFree(d_rxdir); cudaFree(d_out); cudaFree(d_count);

    // ------------------ 组装单结构体输出 ------------------
    const char* fnames[] = {"pos","amp","extra_delay","path_len","order","rx_idx"};
    plhs[0] = mxCreateStructMatrix(1,1, 6, fnames);

    mxArray* Apos = mxCreateDoubleMatrix(h_count, 3, mxREAL);
    mxArray* Aamp = mxCreateDoubleMatrix(h_count, 1, mxREAL);
    mxArray* Adt  = mxCreateDoubleMatrix(h_count, 1, mxREAL);
    mxArray* Alen = mxCreateDoubleMatrix(h_count, 1, mxREAL);
    mxArray* Aord = mxCreateNumericMatrix(h_count, 1, mxINT32_CLASS, mxREAL);
    mxArray* Arxi = mxCreateNumericMatrix(h_count, 1, mxINT32_CLASS, mxREAL);

    double* Ppos = mxGetPr(Apos);
    double* Pamp = mxGetPr(Aamp);
    double* Pdt  = mxGetPr(Adt);
    double* Plen = mxGetPr(Alen);
    int32_t* Pord= (int32_t*)mxGetData(Aord);
    int32_t* Prxi= (int32_t*)mxGetData(Arxi);

    for(int i=0;i<h_count;i++){
        Ppos[i]            = (double)h_out[i].pos.x;
        Ppos[i + h_count]  = (double)h_out[i].pos.y;
        Ppos[i + 2*h_count]= (double)h_out[i].pos.z;
        Pamp[i] = (double)h_out[i].amp;        // amp: 幅度 或 强度（取决于 amp_type）
        Pdt[i]  = (double)h_out[i].extra_dt;
        Plen[i] = (double)h_out[i].path_len;
        Pord[i] = (int32_t)h_out[i].order;
        Prxi[i] = (int32_t)h_out[i].rx_idx;
    }

    mxSetField(plhs[0], 0, "pos",         Apos);
    mxSetField(plhs[0], 0, "amp",         Aamp);
    mxSetField(plhs[0], 0, "extra_delay", Adt);
    mxSetField(plhs[0], 0, "path_len",    Alen);
    mxSetField(plhs[0], 0, "order",       Aord);
    mxSetField(plhs[0], 0, "rx_idx",      Arxi);
}
