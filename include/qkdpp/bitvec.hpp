// include/qkdpp/bitvec.hpp
//
// 紧凑位向量：以 64 位字打包存储 {0,1} 序列。
// 用于 QKD 后处理链的密钥比特流（sifted / corrected / final key）。
//
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace qkdpp {

using u64 = std::uint64_t;

namespace detail {
/// 前置声明：BitVec 的 private 缓冲区仅对本结构开放（见类内 friend 声明）。
/// 完整定义在 BitVec 之后，以便改名/收敛 mutable 访问时只需改一处。
struct BitVecRawAccess;
}  // namespace detail

class BitVec {
public:
    /// 类不变量（所有成员函数均须维持）：
    ///   I1. `w_.size() == nwords_of(nbits_)`，即字数恰为上取整 (nbits_/64)。
    ///   I2. 尾部无效位（下标 >= nbits_ 的高位）恒为 0。
    /// I2 不是装饰性的：parity_and / hamming_distance / to_hex 都直接遍历
    /// 整个字，若尾部混入 1 会静默给出错误结果。
    BitVec() = default;

    /// 构造 nbits 位，fill=1 时全置位（尾部无效位强制清零）
    explicit BitVec(std::size_t nbits, int fill = 0);

    std::size_t size() const noexcept { return nbits_; }
    std::size_t word_count() const noexcept { return w_.size(); }
    bool        empty() const noexcept { return nbits_ == 0; }

    /// @throws std::out_of_range i >= size()
    int at(std::size_t i) const;
    /// @throws std::out_of_range i >= size()
    void set(std::size_t i, int v);
    /// @throws std::out_of_range i >= size()
    void flip(std::size_t i);

    const std::vector<u64>& words() const noexcept { return w_; }

    /// GF(2) 内积：popcount(A & B) mod 2
    /// @pre   `other.size() == size()`
    /// @throws std::invalid_argument 长度不等
    int parity_and(const BitVec& other) const;

    /// Hamming 距离
    /// @pre   `other.size() == size()`
    /// @throws std::invalid_argument 长度不等
    std::size_t hamming_distance(const BitVec& other) const;

    /// 取窗口 [start, start+len) 生成新位向量。
    /// @throws std::out_of_range 当 start > size() 或 len > size() - start。
    static BitVec slice(const BitVec& src, std::size_t start, std::size_t len);

    /// 确定性伪随机位流（splitmix64）。
    ///
    /// ⚠️ **非密码学安全，禁止用于任何安全材料**。splitmix64 是纯确定性
    /// 函数：给定 seed 即可完全重建输出，且输出位可由 seed 经简单代数推出。
    /// 因此**绝不可**用它生成 Toeplitz 隐私放大种子或密钥——
    /// 那会直接击穿隐私放大的安全性。专用接口见：
    ///   - `random_seed_from_os()`（toeplitz.hpp）用于生产种子
    ///   - `privacy_amplify_deterministic()` 内部才用它，且函数名已声明是实验路径
    /// 本函数仅用于可复现实验（测试向量、基准输入）。
    static BitVec random(std::size_t nbits, u64 seed);

    /// 十六进制摘要，高位对齐、**按 nbits 裁剪**到 ceil(nbits/4) 个字符。
    /// 与 Python 侧 golden 值比对时不会带上尾部填充产生的多余前导零，
    /// 也能区分不同比特长度（如 4 位 0xA 与 64 位 0xA 结果不同）。
    std::string to_hex() const;

    static int popcount(u64 x) noexcept;

    // ---- 就地运算的低层接口 ------------------------------------------------
    //
    // 下面两个接口服务于本库内部的高性能例程（gf2poly / toeplitz 的批量
    // xor 与掩码）。它们刻意**没有**把 vector 的可变性整体外泄：既不返回
    // 容器引用，也不允许改变长度，因此 I1 不可能被破坏。

    /// 受控的单字写入：越界抛 std::out_of_range；写入后自动把末字的无效高位
    /// 清零（写入的若是末字），因此 **I2 不会被破坏**。
    /// @throws std::out_of_range idx >= word_count()
    void set_word_masked(std::size_t idx, u64 value);

    /// 就地 XOR 一个字（同上，末字自动掩码）。
    /// @throws std::out_of_range idx >= word_count()
    void xor_word_masked(std::size_t idx, u64 value);

    /// 把末字中下标 >= nbits_ 的无效位清零，恢复不变量 I2。
    /// 供本库内部例程在批量就地写之后调用（O(1)）。
    void normalize_tail() noexcept;

    /// 校验类不变量 I1 / I2（O(字数) 只读遍历）。
    /// **在 release 构建下依然可用** —— 它不做断言，只返回判定结果，
    /// 便于调用方（含 test / bench / 外部绑定）随时兜底检查。
    [[nodiscard]] bool check_invariants() const noexcept;

private:
    /// 暴露给 gf2poly 等例程的连续字缓冲。命名即告诫：只允许读写字内容，
    /// 不得 resize/clear/增删，也不得写入下表 >= nbits_ 的高位（那会破坏
    /// I2，需要时用 normalize_tail() 收尾）。需要改变长度时请新建 BitVec。
    /// 只有 qkdpp::detail::BitVecRawAccess 可访问（见下）。
    friend struct detail::BitVecRawAccess;

    std::size_t      nbits_ = 0;
    std::vector<u64> w_;
};

/// 内部访问器：把「连续字缓冲的可变引用」收敛到单一白名单入口。
///
/// ⚠️ **仅供库内部例程使用**（`bitvec.cpp` / `gf2poly.cpp` 的批量位运算内核）。
/// 之所以放在 `qkdpp::detail` 下而非公共命名空间，是因为它返回的是
/// `std::vector<u64>&` —— 一旦调用方持有它，就可以 `resize/clear` 破坏不变量
/// I1，或写入末字的高位破坏 I2，从而静默污染 parity_and / hamming_distance /
/// to_hex 的结果。放进 `detail` 命名空间是「请勿使用」的显式信号。
///
/// 使用约束（调用方自行保证）：
///   - 只读写 `[0, words().size())` 区间内的字，**不得**改长度；
///   - 不得写入下标 >= `size()` 的尾部高位；若确实写了，须在例程末尾调用
///     `normalize_tail()` 收尾，或依赖 `check_invariants()` 做调试期校验。
namespace detail {
struct BitVecRawAccess {
    static std::vector<u64>& words(BitVec& v) noexcept { return v.w_; }
};
}  // namespace detail

}  // namespace qkdpp