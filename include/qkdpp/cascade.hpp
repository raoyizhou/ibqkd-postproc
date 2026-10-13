// include/qkdpp/cascade.hpp
//
// Cascade 纠错协议（Brassard–Salvail 1994）。
// 对应论文 §IV-A / §V-B 中"error-corrected via the Cascade protocol"
// —— 仿真代码 qkd_protocol.py 中仅有空配置类，此处为参考实现。
//
#pragma once

#include "qkdpp/bitvec.hpp"

#include <cstddef>

namespace qkdpp {

struct CascadeResult {
    std::size_t corrections      = 0;  ///< 翻转的比特总数
    std::size_t parity_bits_sent = 0;  ///< 传输的奇偶校验位总数（泄露量）
    std::size_t rounds           = 0;  ///< 实际执行的轮数
    std::size_t residual_errors  = 0;  ///< 结束后与 Alice 的残余 Hamming 距离
};

class Cascade {
public:
    /// @param q_est       预期 QBER，用于选取首轮块大小 k1 ≈ 0.73/q。
    ///                    必须为有限值且 ∈ (0, 1]；否则抛 std::invalid_argument。
    ///                    传 0 或 NaN 会让 k1 变成 inf/NaN，再经
    ///                    `static_cast<std::size_t>` 落入未定义行为，故在此拒绝。
    /// @param passes      轮数（经典建议 4）；传 0 视为 1
    /// @param seed        置换的确定性随机种子
    /// @param block_size  首轮块大小。**非原样生效**：实际值会被钳制为
    ///                    `min(max(4, block_size), n)`，即传 1–3 会被提升到 4，
    ///                    传大于 n 会被截断为 n；只有钳制后仍 < 2（即 n < 2）
    ///                    才抛 std::invalid_argument。传 0 表示由 q_est 自动推导
    ///                    （k1 ≈ 0.73/q_est，同样经过上述钳制）。
    /// @throws std::invalid_argument 当 q_est 非有限、<= 0 或 > 1
    explicit Cascade(double q_est = 0.03, std::size_t passes = 4, u64 seed = 0,
                     std::size_t block_size = 0);

    /// 就地纠正 bob_key 使其向 alice_key 收敛。
    ///
    /// 前置条件与契约：
    ///   - `bob_key.size() == alice_key.size()`，否则抛 std::invalid_argument。
    ///   - `n == 0` 合法：直接返回全零 CascadeResult。
    ///   - `n == 1` 抛 std::invalid_argument —— 单个比特构不成可校验的块，
    ///     静默返回会让残余错误流入隐私放大阶段。
    ///   - **未收敛不抛异常**：本函数只要执行完就正常返回，残余错误以
    ///     `CascadeResult::residual_errors` 暴露。调用方**必须**检查该字段
    ///     （期望为 0）后再把 bob_key 交给隐私放大，否则可能用错误密钥
    ///     生成「看似成功」的最终密钥。`[[nodiscard]]` 即为强制这一点：
    ///     `cas.correct(bob, alice);` 这种丢弃返回值的写法会触发编译告警。
    ///   - 副作用：`bob_key` 被就地修改；`alice_key` 只读。
    /// @throws std::invalid_argument 长度不匹配，或长度 < 2 且非 0。
    [[nodiscard]] CascadeResult correct(BitVec& bob_key, const BitVec& alice_key) const;

    /// 理论纠错泄露：n·h2(q)·f_ec（与实测 parity_bits_sent 对照）
    static double theoretical_leak(std::size_t n, double q, double ec_efficiency = 1.1);

private:
    double      q_est_;
    std::size_t passes_;
    u64         seed_;
    std::size_t block_size_;
};

}  // namespace qkdpp