// include/qkdpp/toeplitz.hpp
//
// Toeplitz 矩阵隐私放大（privacy amplification）。
// 对应论文 §IV-A "privacy-amplified via Toeplitz-matrix hashing"
// 与 §V-B "privacy amplification"。
//
#pragma once

#include "qkdpp/bitvec.hpp"

#include <cstddef>

namespace qkdpp {

/// 参考实现：逐位三重循环 y[i] = XOR_j T[i][j]·x[j]，T[i][j] = seed[i-j+n-1]。
/// 仅用于正确性对照，复杂度 O(m·n)——**不要**在真实密钥规模（n ≥ 10^6）上调用。
///
/// ⚠️ 安全前提：Toeplitz 矩阵的种子向量必须是**密码学安全的均匀随机比特**
/// 且对窃听者保密。本模块只做位运算，不负责生成种子——调用方必须从
/// CSPRNG（系统熵源 / 经认证的 DRBG，可用下方 `random_seed_from_os()`）
/// 生成完整的 n+m-1 位种子。
///
/// @pre  当 `m > 0` 且 `x.size() > 0` 时 `seed.size() >= x.size() + m - 1`。
/// @throws std::invalid_argument 种子过短。
/// @return 长度 m 的向量；`m == 0` 或 `x.empty()` 时返回 m 位全 0。
BitVec toeplitz_reference(const BitVec& x, const BitVec& seed, std::size_t m);

/// 位打包实现：先反转 x，再对每个输出位取 seed 的 n 位滑窗做 AND + popcount。
/// 每步处理 64 位，复杂度 O(m·n/64)。
///
/// @pre / @throws / @return 同 `toeplitz_reference`。
BitVec toeplitz_multiply(const BitVec& x, const BitVec& seed, std::size_t m);

/// 从操作系统熵源生成 `seed_bits` 位的 Toeplitz 种子（不可复现）。
///
/// 优先使用平台级 CSPRNG：Windows `BCryptGenRandom`、Linux `getrandom(2)`、
/// 其余 POSIX 回落到 `/dev/urandom`。三者都不可用时**抛异常**，而不是静默
/// 退回 `std::random_device`——后者的密码学强度由实现定义，在部分工具链上
/// 会退化为确定性伪随机，而 Toeplitz 的安全性完全依赖种子不可预测。
///
/// @throws std::runtime_error 所有系统熵源均不可用。
BitVec random_seed_from_os(std::size_t seed_bits);

/// @deprecated 语义易误导的旧名（暗示已保证密码学安全）。
/// 现等价于 `random_seed_from_os()`，保留仅为兼容既有调用方。
[[deprecated("use random_seed_from_os(); this name overstates the guarantee")]]
BitVec random_seed_csprng(std::size_t seed_bits);

/// 按 §IV-A 计算输出长度：m = n·(1 − h2(q)) − n·h2(q)·f_ec，钳制到 [0, n]。
///
/// 生产路径需要按 `n + m - 1` 为 CSPRNG 种子定长，若调用方自行复制该公式，
/// 一旦公式/舍入/钳制细节与内部实现不一致就会触发「种子过短」异常或长度
/// 错配。公开此函数即可让调用方在调用前可靠地确定所需种子位数。
///
/// @pre  `0 <= qber <= 1`、`ec_efficiency > 0`，且两者均为有限值。
/// @throws std::invalid_argument 参数非有限或超出区间。
[[nodiscard]] std::size_t output_bits(std::size_t n, double qber,
                                      double ec_efficiency);

/// 隐私放大一次调用的完整结果。
///
/// 各重载对同一组输入返回的 PAResult 完全一致——包括 `n == 0` 或 `m == 0`
/// 时 `leak_ec` 仍按 §IV-A 赋值，而非停留在默认值 0。
struct PAResult {
    BitVec      key;          ///< 最终密钥
    std::size_t input_bits  = 0;  ///< 输入 sifted 比特数 n
    std::size_t output_bits = 0;  ///< 输出密钥比特数 m
    double      h2          = 0.0;///< h2(QBER)
    double      leak_ec     = 0.0;///< 纠错泄露估计（比特）
};

/// 可复现实验路径：用单个 u64 `seed` 派生整段 Toeplitz 种子（splitmix64）。
///
/// ⚠️ 密钥空间上限为 2^64，且种子泄漏即矩阵泄漏——**仅限实验**。函数名带
/// `_deterministic` 后缀即为此意；生产路径请用下面的 BitVec 种子重载。
/// 没有默认参数是刻意的：任何 3 实参调用都会因缺少第 4 参而编译失败，
/// 从而无法把生产路径静默降级到 2^64 密钥空间。
///
/// @pre  `0 <= qber <= 1`、`ec_efficiency > 0`，且两者均为有限值。
/// @throws std::invalid_argument qber 或 ec_efficiency 非有限值 / 超出区间。
[[nodiscard]] PAResult privacy_amplify_deterministic(const BitVec& sifted, double qber,
                                                     u64 seed, double ec_efficiency);

/// 生产路径重载：接受调用方提供的完整 Toeplitz 种子向量。
///
/// @pre  `seed.size() >= n + m - 1`，由 CSPRNG 生成并保持保密；
///       `0 <= qber <= 1`、`ec_efficiency > 0`，且两者均为有限值。
/// @throws std::invalid_argument 种子过短，或 qber/ec_efficiency 非有限值 / 越界。
[[nodiscard]] PAResult privacy_amplify(const BitVec& sifted, double qber,
                                       const BitVec& seed, double ec_efficiency);

}  // namespace qkdpp