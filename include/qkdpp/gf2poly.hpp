// include/qkdpp/gf2poly.hpp
//
// GF(2) 多项式乘法 —— Toeplitz 隐私放大的算法内核。
//
// 关键恒等式：Toeplitz 矩阵乘向量等价于一次多项式卷积
//   y[i] = (seed * x)[i + n - 1],     i = 0 .. m-1
// 因此把 O(m·n/64) 的滑动窗算法替换为 O(N^1.585) 的 Karatsuba，
// 基元用 x86 PCLMULQDQ（carry-less multiply）加速。
//
#pragma once

#include "qkdpp/bitvec.hpp"

#include <cstddef>

namespace qkdpp {

/// GF(2) 多项式乘法（无进位）。
/// @return 长度 `a.size() + b.size() - 1` 的结果；任一输入为空时返回空向量。
BitVec gf2_polymul(const BitVec& a, const BitVec& b);

/// 当前进程**实际**使用的乘法基元是否为硬件 CLMUL。
///
/// 这是运行时判定，等价于「编译期启用了 `__PCLMUL__`」且「当前 CPU 经 CPUID
/// 确认支持 PCLMULQDQ」两个条件同时成立。因此：
///   - 以 `-mpclmul` 编译但运行在无 PCLMULQDQ 的老 CPU 上 → 返回 false，
///     并自动回落到软件实现，**不会**触发 SIGILL；
///   - 未开启 `-mpclmul` 的构建即使 CPU 支持也返回 false（硬件路径不存在）。
bool gf2_has_clmul() noexcept;

/// Toeplitz 矩阵乘向量（多项式卷积实现）。
///
/// @pre 当 `m > 0` 且 `x.size() > 0` 时，须满足
///      `seed.size() >= x.size() + m - 1`。该式用减法形式比较，
///      避免 `n + m - 1` 在 m 接近 SIZE_MAX 时回绕绕过校验。
/// @throws std::invalid_argument 前置条件不满足（种子过短）。
/// @return 长度 m 的向量；`m == 0` 或 `x.empty()` 时返回 m 位全 0。
BitVec toeplitz_multiply_fast(const BitVec& x, const BitVec& seed, std::size_t m);

}  // namespace qkdpp