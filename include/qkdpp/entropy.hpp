// include/qkdpp/entropy.hpp
//
// 信息论工具：与 MP-QKD-TP 论文 §IV-A / §V-B 的隐私放大量一致。
//
#pragma once

#include <cmath>
#include <limits>

namespace qkdpp {

/// 二元熵 h2(q) = -q·log2(q) - (1-q)·log2(1-q)，单位 bit/bit。
///
/// 契约（错误一律以 NaN 表达，恰与本文件 `clamp` 一致）：
///   - `q` 为 NaN          → 返回 NaN（原样传播）。
///   - `q == 0` 或 `q == 1` → 返回 0（h2 在端点连续延拓为 0，**仅端点**）。
///   - `q < 0`、`q > 1`、`±inf` → 返回 NaN。这些是非法 QBER，绝不能
///     静默映射为 h2 = 0。
///   - 其余（0 < q < 1）   → 返回正常的二元熵。
///
/// ⚠️ 为什么越界必须是 NaN 而不是 0：h2 出现在隐私放大量
/// `m = n·(1 − h2) − n·h2·f_ec` 中，h2 = 0 意味着**不扣减任何泄露**，
/// 于是密钥率被高估。QBER 估计出错（负数或 > 1）恰恰会落到越界分支，
/// 若返回 0 就构成安全路径上的静默降级。早期实现正是如此，已修正。
///
/// `[[nodiscard]]` 是刻意的：本函数以返回值（而非异常）传递错误，若调用方
/// 丢弃返回值则错误不可观测。同时也不应把结果与 NaN 一并静默吞掉 —— 上层
/// （如 privacy_amplify）必须显式 `std::isnan()` 检查后再继续计算。
[[nodiscard]] inline double binary_entropy(double q) noexcept {
    if (std::isnan(q)) return q;                 // 传播，绝不静默降级为 0
    if (q == 0.0 || q == 1.0) return 0.0;        // 端点连续延拓
    // 越界（q<0 / q>1）与 ±inf 一律 NaN：见上方「为什么越界必须是 NaN」。
    if (!(q > 0.0) || !(q < 1.0))
        return std::numeric_limits<double>::quiet_NaN();
    return -q * std::log2(q) - (1.0 - q) * std::log2(1.0 - q);
}

/// 将 x 限制到 [lo, hi]。
///
/// 契约（与实现严格一致）：
///   - `std::isnan(x)`  → 返回 NaN（原样传播，不做纠偏）。
///   - `lo`/`hi` 任一为 NaN，或 `lo > hi` → 返回 quiet_NaN。这属调用方
///     编程错误，但本函数 noexcept 且不抛异常，故以 NaN 作为**唯一的**
///     错误信号。
///   - 其余情况返回钳制后的有限值。
///
/// ⚠️ 调用方**必须**检查返回值：NaN 不会中断执行，只会在隐私放大量计算里
/// 持续传播。`[[nodiscard]]` 是编译期的第一道防线。
///
/// 这里刻意**不**用 `assert` 拦截 `lo > hi`：那会让「非法输入 → NaN」这条
/// 对外契约在 debug 构建下变成 abort，两种行为无法共存，也让契约不可测。
/// 统一为返回值语义后，debug 与 release 行为一致。
[[nodiscard]] inline double clamp(double x, double lo, double hi) noexcept {
    if (std::isnan(x)) return x;                 // 传播 NaN
    // lo > hi 属于调用方编程错误：不做静默纠偏，返回 NaN 以便尽早暴露。
    if (std::isnan(lo) || std::isnan(hi) || lo > hi)
        return std::numeric_limits<double>::quiet_NaN();
    return x < lo ? lo : (x > hi ? hi : x);
}

}  // namespace qkdpp