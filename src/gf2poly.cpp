#include "qkdpp/gf2poly.hpp"

#include <algorithm>
#include <cstring>
#include <stdexcept>

#if defined(_MSC_VER)
#include <intrin.h>          // _BitScanForward64 / __cpuid（MSVC 内建，必须显式包含）
#endif
#if defined(__PCLMUL__)
#include <immintrin.h>
#endif

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#define QKDPP_X86 1
#else
#define QKDPP_X86 0
#endif

namespace qkdpp {

namespace {

// ---------------------------------------------------------------------------
// 编译器可移植的内建封装
//
// 本项目以 GCC/Clang 为主，但保留 MSVC 目标（CPU 探测分支即为此而加），
// 因此不能无条件使用 GCC 专有内建——否则 MSVC 下整文件编译失败，
// 连软件回退路径都无法实例化。
// ---------------------------------------------------------------------------

/// 64 位尾零计数（前置条件：x != 0）。
inline int ctz64(u64 x) noexcept {
#if defined(_MSC_VER)
    unsigned long idx = 0;
    _BitScanForward64(&idx, x);      // x != 0 时必然返回非零
    return static_cast<int>(idx);
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(x);
#else
    int n = 0;
    while ((x & 1ull) == 0) { x >>= 1; ++n; }
    return n;
#endif
}

// ---------------------------------------------------------------------------
// CLMUL 能力探测
//
// 两个不同的量必须分开：
//   * 编译期 `__PCLMUL__`  —— 本翻译单元是否被允许发射 PCLMULQDQ 指令。
//     为 0 时执行该指令会直接编译失败，硬件路径根本不存在。
//   * 运行时 CPUID       —— 当前 CPU 是否真的支持 PCLMULQDQ。
//     `-mpclmul` 编译的二进制跑在不支持的 CPU 上会 SIGILL。
//
// 只有两者同时成立才能走硬件路径。仅返回编译期常量的 `gf2_has_clmul()`
// 会造成「报 yes 却崩」或「报 no 却本可加速」两种误判。
// ---------------------------------------------------------------------------
constexpr bool kCompiledWithClmul =
#if defined(__PCLMUL__)
    true;
#else
    false;
#endif

/// 运行时 CPUID 探测 PCLMULQDQ（CPUID leaf 1，ECX bit 1）。
/// 非 x86 平台恒为 false。
bool cpu_supports_clmul() noexcept {
#if QKDPP_X86 && defined(_MSC_VER)
    int regs[4] = {0, 0, 0, 0};
    __cpuid(regs, 1);
    return (regs[2] & (1 << 1)) != 0;   // ECX bit 1 = PCLMULQDQ
#elif QKDPP_X86 && (defined(__GNUC__) || defined(__clang__))
    // __builtin_cpu_supports 由编译器在 init 阶段经 CPUID 填充，可安全调用。
    return __builtin_cpu_supports("pclmul") != 0;
#else
    return false;
#endif
}

/// 仅当「编译期允许」且「运行时支持」时启用硬件路径。
/// 用函数内 static 保证 CPUID 只查一次。
bool clmul_enabled() noexcept {
    static const bool enabled = kCompiledWithClmul && cpu_supports_clmul();
    return enabled;
}

// ---------------------------------------------------------------------------
// 64×64 → 128 位 GF(2) 乘法（无进位）
// ---------------------------------------------------------------------------

inline void clmul64_sw(u64 a, u64 b, u64& lo, u64& hi) noexcept {
    lo = 0;
    hi = 0;
    while (b) {
        const int i = ctz64(b);
        if (i == 0) {
            lo ^= a;
        } else {
            lo ^= (a << i);
            hi ^= (a >> (64 - i));
        }
        b &= (b - 1);
    }
}

inline void clmul64_hw(u64 a, u64 b, u64& lo, u64& hi) noexcept {
#if defined(__PCLMUL__)
    const __m128i va = _mm_set_epi64x(0, static_cast<long long>(a));
    const __m128i vb = _mm_set_epi64x(0, static_cast<long long>(b));
    const __m128i r  = _mm_clmulepi64_si128(va, vb, 0x00);
    // 仅用 SSE2 指令，避免对 SSE4.1 的依赖
    lo = static_cast<u64>(_mm_cvtsi128_si64(r));
    hi = static_cast<u64>(_mm_cvtsi128_si64(_mm_srli_si128(r, 8)));
#else
    (void)a; (void)b; (void)lo; (void)hi;
#endif
}

inline void clmul64(u64 a, u64 b, u64& lo, u64& hi) noexcept {
    if (clmul_enabled())
        clmul64_hw(a, b, lo, hi);
    else
        clmul64_sw(a, b, lo, hi);
}

// ---------------------------------------------------------------------------
// 位向量辅助运算（就地，避免频繁分配）
// ---------------------------------------------------------------------------

inline void trim(BitVec& v) { v.normalize_tail(); }

void xor_into(BitVec& dst, const BitVec& src, std::size_t at_bit) {
    if (src.empty()) return;
    const std::size_t wo = at_bit >> 6;
    const unsigned    bo = static_cast<unsigned>(at_bit & 63u);
    std::vector<u64>& D  = detail::BitVecRawAccess::words(dst);
    const std::vector<u64>& S = src.words();
    const std::size_t n = D.size();

    // 越界丢位本身不一定是 bug：Karatsuba 会把 a0/a1、b0/b1 零填充到公共
    // 宽度 q 再相乘，z1 的高位天然为 0，超出目标长度的部分就该丢掉。真正的
    // bug 是「丢掉的是非零位」——那说明调用方的尺寸/偏移算错了。
    //
    // 这里刻意**不做** assert-only 的检查：assert 在 -DNDEBUG 的发布构建里
    // 被完全移除，静默截断会在正式构建中发生且毫无反馈。改为无条件运行期
    // 判断，代价只是越界分支上的几次比较（正常路径的 `wo+k < n` 本就为真）。
    if (bo == 0) {
        for (std::size_t k = 0; k < S.size(); ++k) {
            if (wo + k < n) {
                D[wo + k] ^= S[k];
            } else if (S[k] != 0) {
                throw std::logic_error(
                    "xor_into: non-zero source word truncated (destination too small)");
            }
        }
    } else {
        const unsigned lsh = 64u - bo;
        for (std::size_t k = 0; k < S.size(); ++k) {
            const u64 low_part  = S[k] << bo;
            const u64 high_part = S[k] >> lsh;
            if (wo + k < n) {
                D[wo + k] ^= low_part;
            } else if (low_part != 0) {
                throw std::logic_error(
                    "xor_into: non-zero source bits truncated (destination too small)");
            }
            if (wo + k + 1 < n) {
                D[wo + k + 1] ^= high_part;
            } else if (high_part != 0) {
                throw std::logic_error(
                    "xor_into: non-zero source bits truncated (destination too small)");
            }
        }
    }
    trim(dst);
}

// 抽取 [start, start+len) 位（按位，非按字）
BitVec extract_bits(const BitVec& v, std::size_t start, std::size_t len) {
    if (len == 0) return BitVec(0, 0);
    BitVec out(len, 0);
    const std::size_t w0 = start >> 6;
    const unsigned    b0 = static_cast<unsigned>(start & 63u);
    std::vector<u64>& O  = detail::BitVecRawAccess::words(out);
    const std::vector<u64>& S = v.words();

    if (b0 == 0) {
        for (std::size_t k = 0; k < O.size(); ++k)
            if (w0 + k < S.size()) O[k] = S[w0 + k];
    } else {
        const unsigned rsh = b0;
        const unsigned lsh = 64u - b0;
        for (std::size_t k = 0; k < O.size(); ++k) {
            const u64 lo = (w0 + k     < S.size()) ? S[w0 + k]     : 0ull;
            const u64 hi = (w0 + k + 1 < S.size()) ? S[w0 + k + 1] : 0ull;
            O[k] = (lo >> rsh) | (hi << lsh);
        }
    }
    trim(out);
    return out;
}

// 64 位字级的 schoolbook 乘法（基元用 CLMUL）
BitVec base_mul(const BitVec& a, const BitVec& b) {
    const std::size_t na = a.size(), nb = b.size();
    if (na == 0 || nb == 0) return BitVec(0, 0);

    const std::size_t wa = a.word_count(), wb = b.word_count();
    BitVec            out(na + nb - 1, 0);
    std::vector<u64>& O = detail::BitVecRawAccess::words(out);

    const std::vector<u64>& A = a.words();
    const std::vector<u64>& B = b.words();

    for (std::size_t i = 0; i < wa; ++i) {
        if (A[i] == 0) continue;
        for (std::size_t j = 0; j < wb; ++j) {
            if (B[j] == 0) continue;
            u64 lo = 0, hi = 0;
            clmul64(A[i], B[j], lo, hi);
            if (i + j     < O.size()) O[i + j]     ^= lo;
            if (i + j + 1 < O.size()) O[i + j + 1] ^= hi;
        }
    }
    trim(out);
    return out;
}

// Karatsuba 阈值（按位计）
constexpr std::size_t kKaratsubaThreshold = 4096;

BitVec poly_mul_rec(const BitVec& a, const BitVec& b) {
    const std::size_t na = a.size(), nb = b.size();
    if (na == 0 || nb == 0) return BitVec(0, 0);
    if (na < kKaratsubaThreshold || nb < kKaratsubaThreshold) return base_mul(a, b);

    const std::size_t h = (std::max(na, nb) + 1) / 2;

    BitVec a0 = extract_bits(a, 0, std::min(h, na));
    BitVec a1 = (na > h) ? extract_bits(a, h, na - h) : BitVec(0, 0);
    BitVec b0 = extract_bits(b, 0, std::min(h, nb));
    BitVec b1 = (nb > h) ? extract_bits(b, h, nb - h) : BitVec(0, 0);

    const std::size_t q = std::max(std::max(a0.size(), a1.size()),
                                   std::max(b0.size(), b1.size()));

    BitVec as(q, 0), bs(q, 0);
    xor_into(as, a0, 0);
    xor_into(as, a1, 0);
    xor_into(bs, b0, 0);
    xor_into(bs, b1, 0);

    BitVec z0 = poly_mul_rec(a0, b0);
    BitVec z2 = poly_mul_rec(a1, b1);
    BitVec z1 = poly_mul_rec(as, bs);

    xor_into(z1, z0, 0);
    xor_into(z1, z2, 0);

    const std::size_t total = na + nb - 1;
    BitVec            out(total, 0);
    xor_into(out, z0, 0);
    if (h < total) xor_into(out, z1, h);
    if (2 * h < total) xor_into(out, z2, 2 * h);
    trim(out);
    return out;
}

}  // namespace

bool gf2_has_clmul() noexcept { return clmul_enabled(); }

BitVec gf2_polymul(const BitVec& a, const BitVec& b) {
    if (a.empty() || b.empty()) return BitVec(0, 0);
    return poly_mul_rec(a, b);
}

BitVec toeplitz_multiply_fast(const BitVec& x, const BitVec& seed, std::size_t m) {
    const std::size_t n = x.size();
    if (m == 0 || n == 0) return BitVec(m, 0);
    // `n + m - 1` wraps when m is near SIZE_MAX, which would bypass this check
    // and let the convolution read past the end of the seed.  Compare in the
    // subtraction form instead.
    const std::size_t needed = n - 1u;   // seed must hold (n-1) + m bits
    if (seed.size() < needed || m > seed.size() - needed)
        throw std::invalid_argument("toeplitz_multiply_fast: seed shorter than n+m-1");

    // y[i] = (seed * x)[i + n - 1]
    const BitVec conv = gf2_polymul(seed, x);
    return extract_bits(conv, n - 1, m);
}

}  // namespace qkdpp