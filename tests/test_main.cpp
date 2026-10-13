// tests/test_main.cpp — 无外部依赖的自测驱动
#include "qkdpp/bitvec.hpp"
#include "qkdpp/cascade.hpp"
#include "qkdpp/entropy.hpp"
#include "qkdpp/gf2poly.hpp"
#include "qkdpp/toeplitz.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace qkdpp;

static int g_checks = 0;
static int g_fails  = 0;

#define CHECK(cond)                                                              \
    do {                                                                         \
        ++g_checks;                                                              \
        if (!(cond)) {                                                           \
            ++g_fails;                                                           \
            std::printf("  [FAIL] line %d: %s\n", __LINE__, #cond);              \
        }                                                                        \
    } while (0)

/// 抑制 `random_seed_csprng` 的 deprecation 警告（本文件刻意测试旧名兼容性）
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

static inline u64 lcg(u64& s) noexcept {
    s = s * 6364136223846793005ull + 1442695040888963407ull;
    return s >> 33;
}

// ---------------------------------------------------------------------------
static void test_bitvec() {
    std::printf("test_bitvec ...\n");
    BitVec v(130, 0);
    CHECK(v.size() == 130);
    CHECK(v.word_count() == 3);

    v.set(0, 1); v.set(63, 1); v.set(64, 1); v.set(129, 1);
    CHECK(v.at(0) == 1); CHECK(v.at(63) == 1);
    CHECK(v.at(64) == 1); CHECK(v.at(129) == 1);
    CHECK(v.at(1) == 0); CHECK(v.at(128) == 0);

    v.flip(1);
    CHECK(v.at(1) == 1);
    v.flip(1);
    CHECK(v.at(1) == 0);

    // 全 1 构造时尾部无效位必须为 0
    BitVec ones(70, 1);
    CHECK(ones.at(69) == 1);
    CHECK(ones.words().back() >> 6 == 0);   // 位 70..127 必须为 0

    // set_word_masked / xor_word_masked 带边界检查，且写入末字时自动掩码
    BitVec mw(64, 0);
    mw.set_word_masked(0, ~u64{0});
    CHECK(mw.at(63) == 1);
    bool ref_threw = false;
    try { mw.set_word_masked(1, 0); }
    catch (const std::out_of_range&) { ref_threw = true; }
    CHECK(ref_threw);
    ref_threw = false;
    try { mw.xor_word_masked(1, 0); }
    catch (const std::out_of_range&) { ref_threw = true; }
    CHECK(ref_threw);

    // 写入末字时自动掩码：70 位向量的末字写全 1 后，无效高位必须仍为 0
    BitVec masked(70, 0);
    masked.set_word_masked(1, ~u64{0});
    CHECK(masked.at(69) == 1);
    CHECK((masked.words().back() >> 6) == 0);
    CHECK(masked.check_invariants());

    // to_hex 按 nbits 裁剪：长度 = ceil(nbits/4)，且能区分不同比特长度
    {
        BitVec h4(4, 0);   h4.set(1, 1); h4.set(3, 1);   // 0b1010 = 0xA
        BitVec h64(64, 0); h64.set(1, 1); h64.set(3, 1);
        CHECK(h4.to_hex() == "a");
        CHECK(h64.to_hex() == "000000000000000a");
        CHECK(h4.to_hex() != h64.to_hex());
        BitVec h70(70, 0); h70.set(69, 1);
        CHECK(h70.to_hex().size() == 18);                // ceil(70/4)
        CHECK(h70.to_hex() == "200000000000000000");     // 位 69 = 2^69
        CHECK(BitVec(0, 0).to_hex().empty());
    }

    // Hamming 距离
    BitVec a(64, 0), b(64, 0);
    a.set(5, 1); a.set(40, 1);
    CHECK(a.hamming_distance(b) == 2);
    CHECK(b.hamming_distance(a) == 2);
    CHECK(a.hamming_distance(a) == 0);

    // slice：含逐位对照（含非 8 字节对齐的起始偏移）
    BitVec src(200, 0);
    for (std::size_t i = 0; i < 200; ++i) src.set(i, static_cast<int>(i & 1u));
    {
        BitVec s = BitVec::slice(src, 10, 50);
        CHECK(s.size() == 50);
        CHECK(s.at(0) == src.at(10));
        CHECK(s.at(49) == src.at(59));
        for (std::size_t i = 0; i < 50; ++i) CHECK(s.at(i) == src.at(10 + i));
    }
    // 跨字、字对齐、以及各种非对齐起点的按字实现在此接受逐位验证
    for (std::size_t start : {std::size_t(0), std::size_t(1), std::size_t(63),
                              std::size_t(64), std::size_t(65), std::size_t(127)}) {
        for (std::size_t len : {std::size_t(1), std::size_t(63), std::size_t(64),
                                std::size_t(65), std::size_t(100)}) {
            if (start + len > 200) continue;
            BitVec s = BitVec::slice(src, start, len);
            CHECK(s.size() == len);
            std::size_t diff = 0;
            for (std::size_t i = 0; i < len; ++i)
                if (s.at(i) != src.at(start + i)) ++diff;
            CHECK(diff == 0);
            // 尾部无效位必须为 0（按字实现容易漏掉掩码）
            const std::size_t rem = len % 64u;
            if (rem) CHECK((s.words().back() >> rem) == 0);
        }
    }

    // random 的确定性与尾部清理
    BitVec r1 = BitVec::random(100, 1234);
    BitVec r2 = BitVec::random(100, 1234);
    CHECK(r1.hamming_distance(r2) == 0);
    CHECK(r1.words().back() >> 36 == 0);     // 100 位：末字仅低 36 位有效
}

// ---------------------------------------------------------------------------
static void test_toeplitz_equivalence() {
    std::printf("test_toeplitz_equivalence ...\n");
    const std::size_t ns[] = {1, 2, 7, 63, 64, 65, 127, 200, 513};
    const std::size_t ms[] = {1, 5, 63, 64, 100, 256};

    for (std::size_t n : ns) {
        for (std::size_t m : ms) {
            BitVec x    = BitVec::random(n, 1234567u + n * 31 + m);
            BitVec seed = BitVec::random(n + m - 1, 7654321u + n * 17 + m * 5);

            BitVec y_ref = toeplitz_reference(x, seed, m);
            BitVec y_opt = toeplitz_multiply(x, seed, m);

            CHECK(y_ref.size() == m);
            CHECK(y_opt.size() == m);
            CHECK(y_ref.hamming_distance(y_opt) == 0);
        }
    }
}

// ---------------------------------------------------------------------------
static void test_toeplitz_linearity() {
    std::printf("test_toeplitz_linearity ...\n");
    const std::size_t n = 300, m = 128;
    BitVec x1 = BitVec::random(n, 111);
    BitVec x2 = BitVec::random(n, 222);
    BitVec seed = BitVec::random(n + m - 1, 333);

    BitVec x12(n, 0);
    for (std::size_t i = 0; i < n; ++i) x12.set(i, x1.at(i) ^ x2.at(i));

    BitVec y1  = toeplitz_multiply(x1, seed, m);
    BitVec y2  = toeplitz_multiply(x2, seed, m);
    BitVec y12 = toeplitz_multiply(x12, seed, m);

    std::size_t diff = 0;
    for (std::size_t i = 0; i < m; ++i)
        if (y12.at(i) != (y1.at(i) ^ y2.at(i))) ++diff;
    CHECK(diff == 0);   // T(x1 ⊕ x2) == T x1 ⊕ T x2
}

// ---------------------------------------------------------------------------
static void test_privacy_amplify() {
    std::printf("test_privacy_amplify ...\n");
    const std::size_t n = 10000;
    const double      q = 0.03;

    BitVec   sifted = BitVec::random(n, 4242);
    PAResult r      = privacy_amplify_deterministic(sifted, q, 99, 1.1);

    // 不逐字复刻实现公式（那样是同义反复：公式错也跟着错），改为断言
    // 「隐私放大必须满足的信息论必要条件」——它们独立于实现的具体舍入方式。
    CHECK(r.input_bits == n);
    CHECK(r.output_bits > 0);
    CHECK(r.output_bits <= n);                 // m 钳制在 [0, n]
    CHECK(r.key.size() == r.output_bits);

    // 必要条件 1：m 不得超过 Shannon 极限 n·(1 − h2(q))（无限 f_ec 时的上界）
    const double h2 = binary_entropy(q);
    const double shannon_cap = static_cast<double>(n) * (1.0 - h2);
    CHECK(static_cast<double>(r.output_bits) <= shannon_cap + 1.0);

    // 必要条件 2：f_ec = 1.1 > 1，故 m 必须严格小于 Shannon 上限
    CHECK(static_cast<double>(r.output_bits) < shannon_cap);

    // 必要条件 3：h2 与 leak_ec 必须自洽（同一 h2 定义）
    CHECK(std::abs(r.h2 - h2) < 1e-12);
    CHECK(std::abs(r.leak_ec - static_cast<double>(n) * h2 * 1.1) < 1e-6);

    // 必要条件 4：公开的 output_bits() 辅助与内部计算必须一致
    CHECK(output_bits(n, q, 1.1) == r.output_bits);

    // 相同输入与种子必须可复现
    PAResult r2 = privacy_amplify_deterministic(sifted, q, 99, 1.1);
    CHECK(r.key.hamming_distance(r2.key) == 0);

    // 不同种子必须给出不同密钥（雪崩）
    PAResult r3 = privacy_amplify_deterministic(sifted, q, 100, 1.1);
    CHECK(r.key.hamming_distance(r3.key) > r.key.size() / 4);

    // 随 q 单调：QBER 越大，可提取的密钥越短
    std::size_t m_low  = privacy_amplify_deterministic(sifted, 0.01, 99, 1.1).output_bits;
    std::size_t m_high = privacy_amplify_deterministic(sifted, 0.08, 99, 1.1).output_bits;
    CHECK(m_low >= m_high);

    std::printf("    n=%zu q=%.2f -> m=%zu (Shannon cap %.1f)\n",
                n, q, r.output_bits, shannon_cap);
}

// ---------------------------------------------------------------------------
// GF(2) 多项式乘法：小尺寸逐位对照 + 大尺寸不变量校验
static BitVec naive_conv(const BitVec& a, const BitVec& b) {
    if (a.empty() || b.empty()) return BitVec(0, 0);
    BitVec c(a.size() + b.size() - 1, 0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (!a.at(i)) continue;
        for (std::size_t j = 0; j < b.size(); ++j)
            if (b.at(j)) c.flip(i + j);
    }
    return c;
}

static void test_gf2_polymul() {
    std::printf("test_gf2_polymul ...\n");
    std::printf("    hardware CLMUL: %s\n", gf2_has_clmul() ? "yes" : "no (software fallback)");

    // 小尺寸：与朴素卷积逐位对照（覆盖 base_mul 与 Karatsuba 两侧阈值）
    const std::size_t small[] = {1, 63, 64, 65, 127, 1023, 4095, 4096, 4097};
    for (std::size_t na : small) {
        for (std::size_t nb : {std::size_t(1), std::size_t(64),
                               std::size_t(1023), std::size_t(4096)}) {
            BitVec a = BitVec::random(na, 0xA1u + na * 7 + nb);
            BitVec b = BitVec::random(nb, 0xB2u + na * 13 + nb * 3);

            BitVec c_naive = naive_conv(a, b);
            BitVec c_fast  = gf2_polymul(a, b);

            CHECK(c_naive.size() == c_fast.size());
            if (c_naive.size() == c_fast.size())
                CHECK(c_naive.hamming_distance(c_fast) == 0);
        }
    }

    // 大尺寸：朴素卷积 O(na·nb) 在此已不可行（20000×12345 ≈ 2.5e8 内层循环），
    // 改用与参考实现无关的不变量 —— 与 seed 无关的代数恒等式与长度契约。
    for (std::size_t na : {std::size_t(9000), std::size_t(20000)}) {
        BitVec a = BitVec::random(na, 0xC3u + na);
        BitVec b = BitVec::random(12345, 0xD4u + na);

        BitVec c = gf2_polymul(a, b);
        CHECK(c.size() == na + b.size() - 1);

        // 乘 1（多项式 x^0）为恒等：a * (1) == a 补齐到长度 na
        BitVec one(1, 0); one.set(0, 1);
        BitVec a_one = gf2_polymul(a, one);
        CHECK(a_one.size() == na);
        if (a_one.size() == na) CHECK(a_one.hamming_distance(a) == 0);

        // 线性性：(a⊕a') * b == a*b ⊕ a'*b
        BitVec a2 = BitVec::random(na, 0xE5u + na);
        BitVec ax(na, 0);
        for (std::size_t i = 0; i < na; ++i) ax.set(i, a.at(i) ^ a2.at(i));
        BitVec c_a2 = gf2_polymul(a2, b);
        BitVec c_ax = gf2_polymul(ax, b);
        BitVec c_lin = gf2_polymul(a, b);
        std::size_t diff = 0;
        for (std::size_t i = 0; i < c_lin.size(); ++i)
            if (c_ax.at(i) != (c_lin.at(i) ^ c_a2.at(i))) ++diff;
        CHECK(diff == 0);
    }
}

static void test_toeplitz_fast() {
    std::printf("test_toeplitz_fast ...\n");
    // 小/中尺寸与参考实现逐位对照；大尺寸样本量收敛以避免 O(n·m) 超时
    const std::size_t ns[] = {1, 63, 64, 65, 200, 4096, 8192};
    const std::size_t ms[] = {1, 5, 64, 1000, 4096};

    for (std::size_t n : ns) {
        for (std::size_t m : ms) {
            BitVec x    = BitVec::random(n, 555u + n * 11 + m);
            BitVec seed = BitVec::random(n + m - 1, 777u + n * 19 + m * 7);

            BitVec y_ref  = toeplitz_reference(x, seed, m);
            BitVec y_fast = toeplitz_multiply_fast(x, seed, m);

            CHECK(y_ref.size() == y_fast.size());
            if (y_ref.size() == y_fast.size())
                CHECK(y_ref.hamming_distance(y_fast) == 0);
        }
    }

    // 大尺寸：只校验长度契约 + 与滑窗实现的抽样一致性（逐位 O(n·m) 太慢）
    {
        const std::size_t n = 20000, m = 12000;
        BitVec x    = BitVec::random(n, 555u + n * 11 + m);
        BitVec seed = BitVec::random(n + m - 1, 777u + n * 19 + m * 7);

        BitVec y_fast = toeplitz_multiply_fast(x, seed, m);
        CHECK(y_fast.size() == m);

        // 抽样 64 个输出位，用参考实现逐点复核
        std::size_t mismatch = 0;
        for (std::size_t s = 0; s < 64; ++s) {
            const std::size_t i = (s * 131u) % m;
            int acc = 0;
            for (std::size_t j = 0; j < n; ++j)
                if (x.at(j)) acc ^= seed.at(i + n - 1 - j);
            if (y_fast.at(i) != acc) ++mismatch;
        }
        CHECK(mismatch == 0);
    }
}

// ---------------------------------------------------------------------------
static void test_cascade() {
    std::printf("test_cascade ...\n");
    const std::size_t n = 8192;
    const double      q = 0.03;

    BitVec alice = BitVec::random(n, 2024);
    BitVec bob   = alice;

    // 注入约 q·n 个不重复错误位
    const std::size_t nerr = static_cast<std::size_t>(n * q);
    u64 st = 13579;
    std::vector<char> used(n, 0);
    std::size_t placed = 0;
    while (placed < nerr) {
        const std::size_t pos = static_cast<std::size_t>(lcg(st) % n);
        if (used[pos]) continue;
        used[pos] = 1;
        bob.flip(pos);
        ++placed;
    }

    CHECK(bob.hamming_distance(alice) == nerr);

    Cascade       cas(/*q_est=*/q, /*passes=*/6, /*seed=*/777);
    CascadeResult res = cas.correct(bob, alice);

    const double theo = Cascade::theoretical_leak(n, q, 1.1);
    std::printf("    injected=%zu  corrected=%zu  residual=%zu\n",
                nerr, res.corrections, res.residual_errors);
    std::printf("    parity_sent=%zu bits (%.4f n)   theoretical_leak=%.1f bits (%.4f n)   ratio=%.2f\n",
                res.parity_bits_sent,
                static_cast<double>(res.parity_bits_sent) / static_cast<double>(n),
                theo, theo / static_cast<double>(n),
                static_cast<double>(res.parity_bits_sent) / theo);

    CHECK(res.corrections >= 1);
    CHECK(res.residual_errors == 0);          // 多轮后应完全收敛
    CHECK(res.parity_bits_sent > 0);
    CHECK(res.parity_bits_sent < n);          // 泄露必须显著小于原始比特数
    CHECK(static_cast<double>(res.parity_bits_sent) < 2.0 * theo);  // 且贴近理论下界
}

// ---------------------------------------------------------------------------
// 边界与溢出回归（OCR 修复项）
// ---------------------------------------------------------------------------
static void test_size_overflow_guards() {
    std::printf("test_size_overflow_guards ...\n");

    // slice 的 `start + len` 回绕不得绕过越界检查
    BitVec v(128, 0);
    bool threw = false;
    try { (void)BitVec::slice(v, 64, static_cast<std::size_t>(-1)); }
    catch (const std::out_of_range&) { threw = true; }
    CHECK(threw);

    threw = false;
    try { (void)BitVec::slice(v, 129, 1); }
    catch (const std::out_of_range&) { threw = true; }
    CHECK(threw);

    // 合法切片仍应工作
    BitVec ok = BitVec::slice(v, 64, 64);
    CHECK(ok.size() == 64);

    // 极端位数不得让 word_count 回绕（否则 at/set 越界）
    BitVec empty(0, 0);
    CHECK(empty.word_count() == 0);
    BitVec one(1, 0);
    CHECK(one.word_count() == 1);
    BitVec exact(64, 0);
    CHECK(exact.word_count() == 1);
    BitVec over(65, 0);
    CHECK(over.word_count() == 2);

    // seed 长度校验在 m 极大时不得回绕
    BitVec x(64, 0), seed(64, 0);
    threw = false;
    try { (void)toeplitz_multiply(x, seed, static_cast<std::size_t>(-1)); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);

    threw = false;
    try { (void)toeplitz_multiply_fast(x, seed, static_cast<std::size_t>(-1)); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
}

static void test_toeplitz_explicit_seed() {
    std::printf("test_toeplitz_explicit_seed ...\n");
    const std::size_t n = 2000;
    const double      q = 0.03;

    BitVec sifted = BitVec::random(n, 777);
    // 生产路径：调用方提供完整种子向量
    BitVec full_seed = random_seed_from_os(n + 2000);
    CHECK(full_seed.size() == n + 2000);

    PAResult r = privacy_amplify(sifted, q, full_seed, 1.1);
    CHECK(r.output_bits > 0);
    CHECK(r.key.size() == r.output_bits);

    // 同一份种子 -> 同一结果（确定性）
    PAResult r2 = privacy_amplify(sifted, q, full_seed, 1.1);
    CHECK(r.key.hamming_distance(r2.key) == 0);

    // OS 熵源接口行为：长度契约 + 位未被全部置 1（粗略非退化检查）。
    // 不做「两次调用必然不同」的断言 —— 那依赖平台的熵保证，可能偶发失败。
    BitVec s2 = random_seed_from_os(256);
    CHECK(s2.size() == 256);
    BitVec zero256(256, 0);
    BitVec one256(256, 1);
    CHECK(s2.hamming_distance(zero256) < 256);   // 不能全 0
    CHECK(s2.hamming_distance(one256) < 256);    // 不能全 1

    // 旧名仍可用（兼容），并给出编译期警告
    BitVec s3 = random_seed_csprng(64);
    CHECK(s3.size() == 64);

    // 过短的种子必须被拒绝，而不是静默使用
    BitVec too_short(16, 0);
    bool threw = false;
    try { (void)privacy_amplify(sifted, q, too_short, 1.1); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);

    // 0 位种子
    BitVec zero_seed(0, 0);
    threw = false;
    try { (void)privacy_amplify(sifted, q, zero_seed, 1.1); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
}

static void test_entropy_nan_propagation() {
    std::printf("test_entropy_nan_propagation ...\n");
    const double nan = std::nan("");

    // NaN 必须被传播，绝不能被静默当成 h2=0（那等于不扣减任何泄露）
    CHECK(std::isnan(binary_entropy(nan)));
    CHECK(std::isnan(clamp(nan, 0.0, 1.0)));

    // 端点连续延拓为 0（这是唯一允许返回 0 的两点）
    CHECK(binary_entropy(0.0) == 0.0);
    CHECK(binary_entropy(1.0) == 0.0);

    // ⚠️ 越界 QBER 必须 NaN，绝不能静默映射为 h2=0 —— 后者意味着
    //    隐私放大量 m = n(1-h2) - n·h2·f_ec 不扣减任何泄露，密钥率被高估。
    //    QBER 估计出错（负数 / >1 / ±inf）恰好落在这条路径上。
    CHECK(std::isnan(binary_entropy(-1.0)));
    CHECK(std::isnan(binary_entropy(2.0)));
    CHECK(std::isnan(binary_entropy(-std::numeric_limits<double>::infinity())));
    CHECK(std::isnan(binary_entropy(std::numeric_limits<double>::infinity())));

    CHECK(std::abs(binary_entropy(0.5) - 1.0) < 1e-12);
    CHECK(clamp(0.5, 0.0, 1.0) == 0.5);
    CHECK(clamp(-1.0, 0.0, 1.0) == 0.0);
    CHECK(clamp(2.0, 0.0, 1.0) == 1.0);

    // lo > hi 属编程错误，以 NaN 暴露而非静默纠偏
    CHECK(std::isnan(clamp(0.5, 1.0, 0.0)));

    // privacy_amplify 必须拒绝 NaN qber，而不是产出错误密钥
    BitVec sifted = BitVec::random(256, 1);
    bool threw = false;
    try { (void)privacy_amplify_deterministic(sifted, nan, 7, 1.1); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);

    // NaN ec_efficiency 同样必须被拒绝（否则 m_d 变 NaN → cast UB）
    threw = false;
    try { (void)privacy_amplify_deterministic(sifted, 0.03, 7, nan); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);

    // ec_efficiency <= 0 会让 m 不降反升、静默跳过隐私放大 → 必须拒绝
    for (double bad : {0.0, -1.0}) {
        threw = false;
        try { (void)privacy_amplify_deterministic(sifted, 0.03, 7, bad); }
        catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
    }

    // qber 越界同样拒绝
    for (double bad_q : {-0.1, 1.5}) {
        threw = false;
        try { (void)privacy_amplify_deterministic(sifted, bad_q, 7, 1.1); }
        catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
    }

    // 空输入不得触发 n+m-1 下溢导致的巨额分配
    BitVec empty(0, 0);
    PAResult re = privacy_amplify_deterministic(empty, 0.03, 7, 1.1);
    CHECK(re.input_bits == 0);
    CHECK(re.output_bits == 0);
    CHECK(re.key.size() == 0);

    // 两条重载对同一输入的 PAResult 必须逐字段一致（含 m==0 时的 leak_ec）。
    // 此前 m==0 的提前返回会让确定性重载的 leak_ec 停在默认 0，
    // 而 BitVec 重载显式赋值，两者行为不一致。
    {
        const double q_hi = 0.5;   // h2(0.5)=1 ⇒ m=0，走提前返回路径
        PAResult da = privacy_amplify_deterministic(sifted, q_hi, 7, 1.1);
        CHECK(da.output_bits == 0);
        // leak_ec 必须仍按 §IV-A 赋值，而不是停在默认 0
        CHECK(da.leak_ec > 0.0);
        CHECK(std::abs(da.leak_ec -
                       static_cast<double>(sifted.size()) * binary_entropy(q_hi) * 1.1)
              < 1e-6);

        BitVec seed_for_full(sifted.size() + 64, 0);
        PAResult db = privacy_amplify(sifted, q_hi, seed_for_full, 1.1);
        CHECK(db.output_bits == da.output_bits);
        CHECK(db.input_bits == da.input_bits);
        CHECK(std::abs(db.leak_ec - da.leak_ec) < 1e-12);
        CHECK(std::abs(db.h2 - da.h2) < 1e-12);
    }
}

static void test_cascade_short_key() {
    std::printf("test_cascade_short_key ...\n");

    // n < 2 无法构成可校验的块：必须抛异常而不是静默返回残余错误
    Cascade c(0.03, 4, 1);
    BitVec a1(1, 0), b1(1, 0);
    bool threw = false;
    try { (void)c.correct(b1, a1); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);

    // n = 0 是合法的空输入，返回零结果
    Cascade c0(0.03, 4, 1);
    BitVec a0(0, 0), b0(0, 0);
    CascadeResult r0 = c0.correct(b0, a0);
    CHECK(r0.rounds == 0);
    CHECK(r0.residual_errors == 0);

    // 小 n（>=2）必须真正执行至少一轮并纠正错误
    // 块大小钳制顺序修正后，n=8 也应至少跑一轮
    Cascade c8(0.25, 4, 5, /*block_size=*/4);
    BitVec a8(8, 0), b8(8, 0);
    a8.set(0, 1); a8.set(3, 1);
    b8.set(3, 1);                       // bob 的 bit0 错了
    CascadeResult r8 = c8.correct(b8, a8);
    CHECK(r8.rounds >= 1);
    CHECK(r8.residual_errors == 0);

    // q_est 契约：0 / NaN / inf / > 1 必须在构造期拒绝
    for (double bad : {0.0, -0.1, std::nan(""), std::numeric_limits<double>::infinity(), 1.5}) {
        bool q_threw = false;
        try { Cascade cbad(bad, 4, 1); (void)cbad; }
        catch (const std::invalid_argument&) { q_threw = true; }
        CHECK(q_threw);
    }
    // 边界值 1.0 合法：精确捕获 std::exception，既能表达「不该抛」的意图，
    // 又不会像 catch (...) 那样把非 std:: 异常（如 SEH、bad_alloc 之外的
    // 实现定义异常）一并吞掉而掩盖真实问题。
    {
        bool ok = true;
        try { Cascade c1(1.0, 4, 1); (void)c1; }
        catch (const std::exception& e) {
            ok = false;
            std::printf("    unexpected exception for q_est=1.0: %s\n", e.what());
        }
        CHECK(ok);
    }
}

// ---------------------------------------------------------------------------
int main() {
    std::printf("=== libqkd-postproc unit tests ===\n\n");

    test_bitvec();
    test_toeplitz_equivalence();
    test_toeplitz_linearity();
    test_gf2_polymul();
    test_toeplitz_fast();
    test_privacy_amplify();
    test_size_overflow_guards();
    test_toeplitz_explicit_seed();
    test_entropy_nan_propagation();
    test_cascade_short_key();
    test_cascade();

    std::printf("\n--- %d checks, %d failures ---\n", g_checks, g_fails);
    if (g_fails == 0) {
        std::printf("ALL TESTS PASSED\n");
        return 0;
    }
    std::printf("TESTS FAILED\n");
    return 1;
}

#pragma GCC diagnostic pop