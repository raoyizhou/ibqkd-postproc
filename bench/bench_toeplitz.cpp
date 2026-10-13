// bench/bench_toeplitz.cpp — Toeplitz 隐私放大吞吐基准
//
// 对照基准：Python 参考实现（MP-QKD-TP sim/netsquid/qkd_protocol.py
// 的 QKDSession._privacy_amplification）对每个输出位执行一次 numpy
// 切片归约，n ≈ 7.6e6 时需约 3×10^6 次 numpy 调用。
//
#include "qkdpp/bitvec.hpp"
#include "qkdpp/entropy.hpp"
#include "qkdpp/gf2poly.hpp"
#include "qkdpp/toeplitz.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>

using namespace qkdpp;

/// 消费整个结果，制造对**所有**输出字的真实数据依赖。
/// 若只读单个比特（如 `y.at(r % m)`），CMake 开启 LTO/内联后编译器可能把
/// 其余位乃至大量乘加当死代码删除，得到虚高吞吐；而 sec==0 的守卫只能识别
/// 「整体被删」，识别不了这种「部分被删」。
static inline u64 sink_of(const BitVec& v) noexcept {
    u64 h = 0;
    for (u64 w : v.words()) h ^= w;
    return h;
}

/// 交叉校验：两条路径必须逐位一致。
///
/// `toeplitz_multiply`（滑动窗）与 `toeplitz_multiply_fast`（Karatsuba+CLMUL）
/// 在 gf2poly.hpp / README 中都承诺「逐位一致」，而本文件是唯一同时持有两者
/// 输出的场所。若快速路径存在正确性缺陷，漂亮的 Mbit/s 数据会完全将其掩盖，
/// 因此必须在计时**之外**先做一次一致性验证。
/// @return 一致返回 true；不一致则打印首个差异并返回 false
static bool cross_check(std::size_t n, std::size_t m) {
    BitVec x    = BitVec::random(n, 20261011u + n);
    BitVec seed = BitVec::random(n + m - 1, 424242u + m);

    const BitVec slow = toeplitz_multiply(x, seed, m);
    const BitVec fast = toeplitz_multiply_fast(x, seed, m);

    if (slow.size() != fast.size()) {
        std::printf("  [MISMATCH] size: slow=%zu fast=%zu (n=%zu m=%zu)\n",
                    slow.size(), fast.size(), n, m);
        return false;
    }
    for (std::size_t i = 0; i < slow.size(); ++i) {
        if (slow.at(i) != fast.at(i)) {
            std::printf("  [MISMATCH] bit %zu: slow=%d fast=%d (n=%zu m=%zu)\n",
                        i, slow.at(i), fast.at(i), n, m);
            return false;
        }
    }
    return true;
}

/// use_fast=false → 滑动窗 O(m·n/64)；true → Karatsuba+CLMUL O(N^1.585)
static void bench_one(std::size_t n, std::size_t m, int reps, bool use_fast) {
    if (reps <= 0) {
        std::printf("  n=%9zu  m=%9zu   -- skipped (reps<=0)\n", n, m);
        return;
    }
    if (m == 0) {
        std::printf("  n=%9zu  m=%9zu   -- skipped (m==0)\n", n, m);
        return;
    }
    BitVec x    = BitVec::random(n, 20261011u + n);
    BitVec seed = BitVec::random(n + m - 1, 424242u + m);

    BitVec warm = use_fast ? toeplitz_multiply_fast(x, seed, m)
                           : toeplitz_multiply(x, seed, m);
    volatile u64 sink0 = sink_of(warm);
    (void)sink0;

    volatile u64 sink = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (int r = 0; r < reps; ++r) {
        BitVec y = use_fast ? toeplitz_multiply_fast(x, seed, m)
                            : toeplitz_multiply(x, seed, m);
        sink ^= sink_of(y);
    }
    const auto t1 = std::chrono::steady_clock::now();
    (void)sink;

    const double total = std::chrono::duration<double>(t1 - t0).count();
    const double sec   = total / static_cast<double>(reps);

    // 计时器分辨率不足时 elapsed 可能为 0（或被优化掉），此时吞吐为 inf；
    // 直接打印 inf 会污染结果表而不报错，因此显式标注不可测。
    if (!(sec > 0.0) || !std::isfinite(sec)) {
        std::printf("  n=%9zu  m=%9zu   %10.6f s/run   throughput: n/a (timer resolution)\n",
                    n, m, sec);
        return;
    }
    const double mbit_per_s = (static_cast<double>(m) / 1e6) / sec;

    std::printf("  n=%9zu  m=%9zu   %10.6f s/run   %9.3f Mbit/s\n",
                n, m, sec, mbit_per_s);
}

/// 每个用例独立捕获异常，使坏用例不影响其余结果表。
/// 否则某条用例抛出（种子过短 invalid_argument、大 n 分配失败 bad_alloc…）
/// 会越过 main 直接 std::terminate，把已算出的行连同诊断一起丢掉 —— 这与
/// main 里 `std::setvbuf(..., _IONBF, 0)` 表达的「崩溃时保留进度」意图矛盾。
static void safe_bench(std::size_t n, std::size_t m, int reps, bool use_fast) {
    try {
        bench_one(n, m, reps, use_fast);
    } catch (const std::exception& e) {
        std::printf("  n=%9zu  m=%9zu   FAILED: %s\n", n, m, e.what());
    } catch (...) {
        std::printf("  n=%9zu  m=%9zu   FAILED: unknown non-std:: exception\n", n, m);
    }
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);   // 无缓冲，避免崩溃时丢失进度
    std::printf("=== libqkd-postproc Toeplitz benchmark ===\n");
    std::printf("hw CLMUL: %s | QBER=3%% f_ec=1.1 -> m/n = 0.5917\n\n",
                gf2_has_clmul() ? "yes" : "no");

    const double ratio = 1.0 - binary_entropy(0.03) * 2.1;   // ≈ 0.5917

    std::printf("[A] sliding-window  O(m*n/64)  -- v1 baseline\n");
    safe_bench(50000, static_cast<std::size_t>(50000 * ratio), 20, false);

    std::printf("\n[B] Karatsuba + CLMUL  O(N^1.585)  -- v2\n");
    // 单个结构体数组：避免 ns/rp 两个隐式等长并行数组不同步导致 rp[k] 越界。
    struct Case { std::size_t n; int reps; };
    const Case cases[] = {
        {50000, 200}, {200000, 50}, {1000000, 20}, {4000000, 5}, {8000000, 2},
    };
    for (const Case& c : cases) {
        safe_bench(c.n, static_cast<std::size_t>(static_cast<double>(c.n) * ratio),
                   c.reps, true);
    }

    // ---------------------------------------------------------------------
    // 正确性交叉校验：必须在吞吐测量之后做（避免污染预热与计时），但必须
    // 影响退出码 —— 否则「快速路径有缺陷」会被漂亮的 Mbit/s 完全掩盖。
    // ---------------------------------------------------------------------
    std::printf("\n[C] cross-check fast vs sliding-window (bit-exact)\n");
    bool all_ok = true;
    const struct { std::size_t n, m; } xcases[] = {
        {64, 37}, {1000, 592}, {50000, 29585}, {200000, 118340},
    };
    for (const auto& c : xcases) {
        const bool ok = cross_check(c.n, c.m);
        std::printf("  n=%9zu  m=%9zu   %s\n", c.n, c.m,
                    ok ? "consistent" : "INCONSISTENT");
        all_ok = all_ok && ok;
    }

    std::printf("\n%s\n", all_ok ? "done." : "done (WITH MISMATCHES — see above).");
    return all_ok ? 0 : 1;
}
