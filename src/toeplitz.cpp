#include "qkdpp/toeplitz.hpp"

#include "qkdpp/entropy.hpp"

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <vector>

// ---- 平台级 CSPRNG 头文件 ------------------------------------------------
#if defined(_WIN32)
#include <windows.h>
#include <bcrypt.h>
#if defined(_MSC_VER)
#pragma comment(lib, "bcrypt.lib")   // MSVC 自动链接；MinGW 用 -lbcrypt
#endif
#elif defined(__linux__)
#include <sys/random.h>
#include <cerrno>
#else
#include <cerrno>
#include <cstdio>
#endif

namespace qkdpp {

namespace {

/// §IV-A 的输出长度：m = n·(1 − h2(q)) − n·h2(q)·f_ec，钳制到 [0, n]。
///
/// 这是全局**唯一**的 m 计算实现。所有重载与内部路径都必须经由此函数，
/// 否则「种子按 m1 生成、实际按 m2 使用」的错配会导致抛异常或长度不符。
///
/// 前置条件：qber 与 ec_efficiency 必须已由调用方确认非 NaN（本函数是
/// 内部工具，不做参数校验，只保证自己不产生 UB）。
inline std::size_t output_length_for(std::size_t n, double qber,
                                     double ec_efficiency) noexcept {
    if (n == 0) return 0;
    const double h2  = binary_entropy(clamp(qber, 0.0, 1.0));
    if (std::isnan(h2)) return 0;
    const double leak = static_cast<double>(n) * h2 * ec_efficiency;
    double m_d = static_cast<double>(n) * (1.0 - h2) - leak;
    if (!(m_d > 0.0)) return 0;         // 同时吞掉 NaN 与负值
    std::size_t m = static_cast<std::size_t>(m_d);
    if (m > n) m = n;
    return m;
}

/// 统一的入参校验：qber 与 ec_efficiency 必须为有限值且在有效区间内。
/// 两个重载共用，保证相同非法输入在两条路径上行为一致（都抛异常，
/// 且都在任何内存分配之前抛出）。
inline void validate_params(double qber, double ec_efficiency) {
    if (std::isnan(qber) || std::isnan(ec_efficiency))
        throw std::invalid_argument(
            "privacy_amplify: qber and ec_efficiency must not be NaN");
    if (!std::isfinite(ec_efficiency) || ec_efficiency <= 0.0)
        throw std::invalid_argument(
            "privacy_amplify: ec_efficiency must be a finite positive value "
            "(<= 0 would silently skip privacy amplification)");
    if (!std::isfinite(qber) || qber < 0.0 || qber > 1.0)
        throw std::invalid_argument(
            "privacy_amplify: qber must be a finite value in [0, 1]");
}

/// 从字节源填充位向量。要求源在整个过程中可靠（不可用即抛异常）。
/// \param fill 把 [0, nbytes) 个随机字节写入 buf，失败抛 std::runtime_error。
template <class Filler>
BitVec seed_from_filler(std::size_t seed_bits, Filler&& fill) {
    BitVec out(seed_bits, 0);
    if (seed_bits == 0) return out;

    const std::size_t nbytes = (seed_bits + 7u) / 8u;
    std::vector<unsigned char> buf(nbytes);
    fill(buf.data(), nbytes);

    for (std::size_t i = 0; i < seed_bits; ++i)
        out.set(i, (buf[i >> 3] >> (i & 7u)) & 1u);
    return out;
}

}  // namespace

// 公开的算长辅助：调用方据此为 CSPRNG 种子定长，避免自行复制公式导致
// 「种子过短」异常或长度错配。
std::size_t output_bits(std::size_t n, double qber, double ec_efficiency) {
    validate_params(qber, ec_efficiency);
    return output_length_for(n, qber, ec_efficiency);
}

BitVec toeplitz_reference(const BitVec& x, const BitVec& seed, std::size_t m) {
    const std::size_t n = x.size();
    // 减法形式比较：`n + m - 1` 在 m 接近 SIZE_MAX 时回绕，会绕过校验。
    if (m != 0 && n != 0) {
        const std::size_t needed = n - 1u;
        if (seed.size() < needed || m > seed.size() - needed)
            throw std::invalid_argument("toeplitz_reference: seed shorter than n+m-1");
    }

    BitVec y(m, 0);
    for (std::size_t i = 0; i < m; ++i) {
        int acc = 0;
        for (std::size_t j = 0; j < n; ++j) {
            if (x.at(j)) acc ^= seed.at(i + n - 1 - j);
        }
        y.set(i, acc);
    }
    return y;
}

// ---------------------------------------------------------------------------
// 数学推导（用于核对实现）
//
//   y[i] = XOR_{j=0}^{n-1} seed[i - j + n - 1] · x[j]
//   令 t = i + n - 1 - j  ⇒  j = i + n - 1 - t，t 从 i 到 i+n-1
//   y[i] = XOR_{t=0}^{n-1} seed[i + t] · x[n - 1 - t]
//   定义 xr[t] = x[n - 1 - t]（x 的位反转）
//   ⇒ y[i] = parity( seed[i .. i+n-1]  AND  xr[0 .. n-1] )
//
//   即：seed 上长度为 n 的滑窗 与 反转后的输入 做 GF(2) 内积。
//   位打包后每次处理 64 位，滑窗跨字用移位拼接。
// ---------------------------------------------------------------------------
BitVec toeplitz_multiply(const BitVec& x, const BitVec& seed, std::size_t m) {
    const std::size_t n = x.size();
    if (m == 0 || n == 0) return BitVec(m, 0);
    {
        const std::size_t needed = n - 1u;
        if (seed.size() < needed || m > seed.size() - needed)
            throw std::invalid_argument("toeplitz_multiply: seed shorter than n+m-1");
    }

    BitVec xr(n, 0);
    for (std::size_t j = 0; j < n; ++j)
        xr.set(j, x.at(n - 1 - j));

    const std::size_t xw       = xr.word_count();
    const std::size_t sw_total = seed.word_count();
    const u64* const  S        = seed.words().data();
    const u64* const  X        = xr.words().data();

    BitVec y(m, 0);
    for (std::size_t i = 0; i < m; ++i) {
        const std::size_t wi  = i >> 6;
        const unsigned    off = static_cast<unsigned>(i & 63u);
        int               acc = 0;

        if (off == 0) {
            for (std::size_t k = 0; k < xw; ++k) {
                const std::size_t a  = wi + k;
                const u64         sw = (a < sw_total) ? S[a] : 0ull;
                acc ^= BitVec::popcount(sw & X[k]) & 1;
            }
        } else {
            const unsigned lsh = 64u - off;   // off ∈ [1,63] ⇒ lsh ∈ [1,63]，无 UB
            for (std::size_t k = 0; k < xw; ++k) {
                const std::size_t a  = wi + k;
                const u64         lo = (a     < sw_total) ? S[a]     : 0ull;
                const u64         hi = (a + 1 < sw_total) ? S[a + 1] : 0ull;
                const u64         sw = (lo >> off) | (hi << lsh);
                acc ^= BitVec::popcount(sw & X[k]) & 1;
            }
        }
        y.set(i, acc);
    }
    return y;
}

// ---------------------------------------------------------------------------
// 通用内核：给定已校验的参数与 m，用调用方提供的完整种子做压缩。
// 两个重载都汇聚到这里，杜绝 m 的重复计算。
//
// 关于 leak_ec：它在**所有**返回路径上都被赋值。此前 m==0 的提前返回会
// 让 leak_ec 停留在默认值 0，与确定性重载显式设置 leak_ec 的行为不一致。
// 现在统一在函数开头计算，两条重载对同一输入给出完全相同的 PAResult。
// ---------------------------------------------------------------------------
namespace {

PAResult make_result(const BitVec& sifted, double qber, double ec_efficiency,
                     std::size_t m) {
    PAResult r;
    r.input_bits  = sifted.size();
    r.h2          = binary_entropy(clamp(qber, 0.0, 1.0));
    r.output_bits = m;
    r.leak_ec     = static_cast<double>(sifted.size()) * r.h2 * ec_efficiency;
    return r;
}

PAResult amplify_with_seed(const BitVec& sifted, double qber, double ec_efficiency,
                           std::size_t m, const BitVec& seed) {
    PAResult r = make_result(sifted, qber, ec_efficiency, m);

    const std::size_t n = sifted.size();
    if (n == 0 || m == 0) return r;

    // 种子必须覆盖 n+m-1 位；比较用减法形式避免 n+m-1 回绕。
    const std::size_t needed = n - 1u;
    if (seed.size() < needed || m > seed.size() - needed)
        throw std::invalid_argument(
            "privacy_amplify: Toeplitz seed must be at least n+m-1 bits "
            "(generate it from a CSPRNG, see random_seed_from_os)");

    r.key = toeplitz_multiply(sifted, seed, m);
    return r;
}

}  // namespace

PAResult privacy_amplify_deterministic(const BitVec& sifted, double qber, u64 seed,
                                       double ec_efficiency) {
    // 与 BitVec 重载完全一致的入口校验：相同非法输入在两条路径上行为相同，
    // 且都在任何分配之前抛出。
    validate_params(qber, ec_efficiency);

    const std::size_t n = sifted.size();
    if (n == 0) return make_result(sifted, qber, ec_efficiency, 0);

    const std::size_t m = output_length_for(n, qber, ec_efficiency);
    if (m == 0) return make_result(sifted, qber, ec_efficiency, 0);

    // ⚠️ 单 u64 种子仅可用于可复现实验（密钥空间 ≤ 2^64）。
    // m 与下游 amplify_with_seed 用的是**同一个** output_length_for 结果，
    // 不存在「按 m1 生成种子、按 m2 使用」的错配。
    const BitVec seedvec = BitVec::random(n + m - 1, seed);
    return amplify_with_seed(sifted, qber, ec_efficiency, m, seedvec);
}

PAResult privacy_amplify(const BitVec& sifted, double qber, const BitVec& seed,
                         double ec_efficiency) {
    validate_params(qber, ec_efficiency);

    const std::size_t n = sifted.size();
    if (n == 0) return make_result(sifted, qber, ec_efficiency, 0);

    const std::size_t m = output_length_for(n, qber, ec_efficiency);
    return amplify_with_seed(sifted, qber, ec_efficiency, m, seed);
}

// ---------------------------------------------------------------------------
// 平台级 CSPRNG
// ---------------------------------------------------------------------------
namespace {

#if defined(_WIN32)

BitVec seed_from_os_impl(std::size_t seed_bits) {
    return seed_from_filler(seed_bits, [](unsigned char* buf, std::size_t nbytes) {
        // BCryptGenRandom 是 Windows 上受支持的进程级 CSPRNG。
        if (BCryptGenRandom(nullptr, buf, static_cast<ULONG>(nbytes),
                            BCRYPT_USE_SYSTEM_PREFERRED_RNG) != 0)
            throw std::runtime_error(
                "random_seed_from_os: BCryptGenRandom failed (no OS entropy available)");
    });
}

#elif defined(__linux__)

BitVec seed_from_os_impl(std::size_t seed_bits) {
    return seed_from_filler(seed_bits, [](unsigned char* buf, std::size_t nbytes) {
        std::size_t got = 0;
        while (got < nbytes) {
            const ssize_t r = ::getrandom(buf + got, nbytes - got, 0);
            if (r < 0) {
                if (errno == EINTR) continue;
                throw std::runtime_error(
                    "random_seed_from_os: getrandom() failed (no OS entropy available)");
            }
            // getrandom(2) 在请求长度非 0 时理论上不返回 0；但若返回 0，
            // 既非错误也不前进，while 条件恒真会陷入无限自旋。显式处理。
            if (r == 0) {
                throw std::runtime_error(
                    "random_seed_from_os: getrandom() returned 0 (would spin forever)");
            }
            got += static_cast<std::size_t>(r);
        }
    });
}

#else  // POSIX 兜底：/dev/urandom

BitVec seed_from_os_impl(std::size_t seed_bits) {
    return seed_from_filler(seed_bits, [](unsigned char* buf, std::size_t nbytes) {
        // RAII：任何后续插桩若抛异常也不会泄漏 FILE*。
        using FilePtr = std::unique_ptr<std::FILE, int (*)(std::FILE*)>;
        std::FILE* raw = std::fopen("/dev/urandom", "rb");
        if (raw == nullptr)
            throw std::runtime_error(
                "random_seed_from_os: cannot open /dev/urandom (no OS entropy available)");
        FilePtr fp(raw, &std::fclose);

        // fread 允许短读（尤其对字符设备在信号中断或大缓冲时），必须循环
        // 累加直到读满，否则正常情况会被误判为致命错误。
        std::size_t got = 0;
        while (got < nbytes) {
            const std::size_t r = std::fread(buf + got, 1, nbytes - got, fp.get());
            if (r == 0) {
                if (std::ferror(fp.get()))
                    throw std::runtime_error(
                        "random_seed_from_os: read error on /dev/urandom");
                throw std::runtime_error(
                    "random_seed_from_os: short read from /dev/urandom "
                    "(EOF before buffer filled)");
            }
            got += r;
        }
    });
}

#endif

}  // namespace

BitVec random_seed_from_os(std::size_t seed_bits) {
    return seed_from_os_impl(seed_bits);
}

BitVec random_seed_csprng(std::size_t seed_bits) {
    return random_seed_from_os(seed_bits);
}

}  // namespace qkdpp