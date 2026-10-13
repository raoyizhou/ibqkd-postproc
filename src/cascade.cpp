#include "qkdpp/cascade.hpp"

#include "qkdpp/entropy.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace qkdpp {

namespace {

inline u64 splitmix64(u64& state) noexcept {
    state += 0x9E3779B97F4A7C15ull;
    u64 z = state;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

/// 计算位向量在 [lo, hi) 区间内的奇偶（popcount mod 2），按 64 位字并行。
inline int parity_of_range(const BitVec& v, std::size_t lo, std::size_t hi) {
    if (hi <= lo) return 0;
    const std::vector<u64>& W    = v.words();
    const std::size_t       wlo  = lo >> 6;
    const std::size_t       whi  = (hi - 1) >> 6;
    int                     acc  = 0;

    if (wlo == whi) {
        u64 mask = (~u64{0}) << (lo & 63u);
        const unsigned top = static_cast<unsigned>((hi - 1) & 63u);
        if (top < 63u) mask &= ((u64{1} << (top + 1)) - 1u);
        return BitVec::popcount(W[wlo] & mask) & 1;
    }

    {
        u64 mask = (~u64{0}) << (lo & 63u);
        acc ^= BitVec::popcount(W[wlo] & mask) & 1;
    }
    for (std::size_t w = wlo + 1; w < whi; ++w)
        acc ^= BitVec::popcount(W[w]) & 1;
    {
        const unsigned top  = static_cast<unsigned>((hi - 1) & 63u);
        const u64      mask = (top >= 63u) ? ~u64{0} : ((u64{1} << (top + 1)) - 1u);
        acc ^= BitVec::popcount(W[whi] & mask) & 1;
    }
    return acc;
}

}  // namespace

Cascade::Cascade(double q_est, std::size_t passes, u64 seed, std::size_t block_size)
    : q_est_(q_est), passes_(passes == 0 ? 1 : passes), seed_(seed), block_size_(block_size) {
    // q_est 是 k1 ≈ 0.73/q 的分母：0 会得到 +inf，NaN 会得到 NaN，两者经
    // `static_cast<std::size_t>` 都是未定义行为，且会在很晚才表现为静默错误。
    // > 1 的 QBER 无物理意义，同样拒绝。与 entropy.hpp / privacy_amplify 的
    // 「NaN 不静默降级」策略保持一致。
    if (std::isnan(q_est) || !std::isfinite(q_est) || !(q_est > 0.0) || q_est > 1.0)
        throw std::invalid_argument(
            "Cascade: q_est must be a finite value in (0, 1]");
}

double Cascade::theoretical_leak(std::size_t n, double q, double ec_efficiency) {
    return static_cast<double>(n) * binary_entropy(clamp(q, 0.0, 1.0)) * ec_efficiency;
}

CascadeResult Cascade::correct(BitVec& bob_key, const BitVec& alice_key) const {
    const std::size_t n = bob_key.size();
    if (alice_key.size() != n)
        throw std::invalid_argument("Cascade::correct key size mismatch");

    CascadeResult res;
    if (n == 0) return res;

    // 首轮块大小：Brassard–Salvail 建议 k1 ≈ 0.73/q
    std::size_t block = block_size_;
    if (block == 0) {
        // 构造函数已保证 q_est_ ∈ (0,1]，此处再兜一层：`std::max(NaN, x)`
        // 会返回 NaN，而 `std::max(x, NaN)` 返回 x —— 顺序依赖且不可靠。
        // 用显式比较保证得到的是有限正数。
        const double q = (std::isfinite(q_est_) && q_est_ > 1e-4) ? q_est_ : 1e-4;
        block          = static_cast<std::size_t>(0.73 / q);
    }
    // 钳制顺序必须是 min(max(4, block), n)：先取下界再取上界。
    // 反之在 n < 4 时会得到 block = 4 > n，循环条件 `block <= n` 一轮都不执行，
    // 残余错误被静默带入隐私放大阶段。
    block = std::min<std::size_t>(std::max<std::size_t>(4u, block), n);
    if (block < 2) {
        // 单个比特无法构成可校验的块，必须显式报告而非静默返回
        throw std::invalid_argument(
            "Cascade::correct: key too short for block-level parity checks (need n >= 2)");
    }

    u64                      st = seed_;
    std::vector<std::size_t> perm(n);
    BitVec                   ap(n, 0), bp(n, 0);

    for (std::size_t pass = 0; pass < passes_ && block <= n; ++pass) {
        ++res.rounds;

        // Fisher–Yates 确定性置换
        for (std::size_t i = 0; i < n; ++i) perm[i] = i;
        for (std::size_t i = n; i > 1; --i) {
            const std::size_t j = static_cast<std::size_t>(splitmix64(st) % i);
            std::swap(perm[i - 1], perm[j]);
        }

        // 物化置换后的视图，使分块奇偶校验可按字并行
        for (std::size_t i = 0; i < n; ++i) {
            ap.set(i, alice_key.at(perm[i]));
            bp.set(i, bob_key.at(perm[i]));
        }

        for (std::size_t start = 0; start < n; start += block) {
            const std::size_t end = std::min(start + block, n);
            if (end - start < 2) continue;

            ++res.parity_bits_sent;
            if (parity_of_range(ap, start, end) == parity_of_range(bp, start, end))
                continue;

            // 二分定位块内错误位
            std::size_t lo = start, hi = end;
            while (hi - lo > 1) {
                const std::size_t mid = lo + (hi - lo) / 2;
                ++res.parity_bits_sent;
                if (parity_of_range(ap, lo, mid) != parity_of_range(bp, lo, mid))
                    hi = mid;
                else
                    lo = mid;
            }
            bob_key.flip(perm[lo]);
            // 同步更新 bp 视图：bp 是 bob_key 的置换镜像，翻转后必须立即
            // 跟上，否则同一轮内 bp 相对 bob_key 变为过期镜像。当前各块区间
            // 互不重叠、也未对已纠正块做复核，所以尚未出错；但 Cascade 的
            // 关键步骤正是「纠正后对该块再分割/复核」（用于捕获块内多重错误），
            // 一旦补上就会读到过期奇偶并把正确比特翻回错误。此处提前对齐。
            bp.flip(lo);
            ++res.corrections;
        }

        block *= 2;
    }

    res.residual_errors = bob_key.hamming_distance(alice_key);
    return res;
}

}  // namespace qkdpp
