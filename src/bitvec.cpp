#include "qkdpp/bitvec.hpp"

#include <algorithm>
#include <stdexcept>

namespace qkdpp {

namespace {
/// 上取整除以 64。用减法形式避免 `nbits + 63` 在 nbits 接近 SIZE_MAX 时回绕，
/// 否则会算出错误的字数，令 at/set/flip 的 `w_[i >> 6]` 越界（UB）。
inline std::size_t nwords_of(std::size_t nbits) noexcept {
    if (nbits == 0) return 0;
    return 1u + (nbits - 1u) / 64u;
}

/// splitmix64：确定性、可复现，仅用于实验随机位流
inline u64 splitmix64(u64& state) noexcept {
    state += 0x9E3779B97F4A7C15ull;
    u64 z = state;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}
}  // namespace

int BitVec::popcount(u64 x) noexcept {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);
#else
    int c = 0;
    while (x) { x &= (x - 1); ++c; }
    return c;
#endif
}

BitVec::BitVec(std::size_t nbits, int fill) : nbits_(nbits), w_(nwords_of(nbits), 0) {
    if (fill) {
        std::fill(w_.begin(), w_.end(), ~u64{0});
        normalize_tail();
    }
}

void BitVec::normalize_tail() noexcept {
    const std::size_t rem = nbits_ % 64u;
    if (rem && !w_.empty()) w_.back() &= ((u64{1} << rem) - 1u);
}

bool BitVec::check_invariants() const noexcept {
    if (w_.size() != nwords_of(nbits_)) return false;   // I1
    const std::size_t rem = nbits_ % 64u;
    if (rem && !w_.empty() && (w_.back() >> rem) != 0) return false;  // I2
    return true;
}

void BitVec::set_word_masked(std::size_t idx, u64 value) {
    if (idx >= w_.size()) throw std::out_of_range("BitVec::set_word_masked out of range");
    w_[idx] = value;
    // 写入的若是末字，立刻把无效高位清掉，维持 I2。
    if (idx + 1 == w_.size()) normalize_tail();
}

void BitVec::xor_word_masked(std::size_t idx, u64 value) {
    if (idx >= w_.size()) throw std::out_of_range("BitVec::xor_word_masked out of range");
    w_[idx] ^= value;
    if (idx + 1 == w_.size()) normalize_tail();
}

int BitVec::at(std::size_t i) const {
    if (i >= nbits_) throw std::out_of_range("BitVec::at out of range");
    return static_cast<int>((w_[i >> 6] >> (i & 63u)) & 1u);
}

void BitVec::set(std::size_t i, int v) {
    if (i >= nbits_) throw std::out_of_range("BitVec::set out of range");
    u64&      word = w_[i >> 6];
    const u64 mask = u64{1} << (i & 63u);
    if (v) word |= mask;
    else   word &= ~mask;
}

void BitVec::flip(std::size_t i) {
    if (i >= nbits_) throw std::out_of_range("BitVec::flip out of range");
    w_[i >> 6] ^= (u64{1} << (i & 63u));
}

int BitVec::parity_and(const BitVec& other) const {
    if (nbits_ != other.nbits_)
        throw std::invalid_argument("BitVec::parity_and size mismatch");
    int acc = 0;
    for (std::size_t k = 0; k < w_.size(); ++k)
        acc ^= (popcount(w_[k] & other.w_[k]) & 1);
    return acc;
}

std::size_t BitVec::hamming_distance(const BitVec& other) const {
    if (nbits_ != other.nbits_)
        throw std::invalid_argument("BitVec::hamming_distance size mismatch");
    std::size_t d = 0;
    for (std::size_t k = 0; k < w_.size(); ++k)
        d += static_cast<std::size_t>(popcount(w_[k] ^ other.w_[k]));
    return d;
}

BitVec BitVec::slice(const BitVec& src, std::size_t start, std::size_t len) {
    // Subtraction form: `start + len` can wrap to a small value when len is
    // huge, which would defeat the bounds check and cause an oversized
    // allocation followed by out-of-range reads.
    if (start > src.size() || len > src.size() - start)
        throw std::out_of_range("BitVec::slice out of range");

    BitVec out(len, 0);
    if (len == 0) return out;

    // 按 64 位字批量搬运：逐位 set/at 对 n=10^6 的密钥约有 64 倍常数开销。
    const std::size_t        s0 = start >> 6;
    const unsigned           sh = static_cast<unsigned>(start & 63u);
    const std::vector<u64>&  S  = src.w_;
    std::vector<u64>&        O  = out.w_;

    if (sh == 0) {
        for (std::size_t k = 0; k < O.size(); ++k)
            if (s0 + k < S.size()) O[k] = S[s0 + k];
    } else {
        // sh ∈ [1,63] ⇒ lsh ∈ [1,63]，不存在 `>> 64` / `<< 64` 这类 UB；
        // 越界读取以 0 兜底（src 尾部无效位本就为 0）。
        const unsigned lsh = 64u - sh;
        for (std::size_t k = 0; k < O.size(); ++k) {
            const u64 lo = (s0 + k     < S.size()) ? S[s0 + k]     : 0ull;
            const u64 hi = (s0 + k + 1 < S.size()) ? S[s0 + k + 1] : 0ull;
            O[k] = (lo >> sh) | (hi << lsh);
        }
    }
    out.normalize_tail();
    return out;
}

BitVec BitVec::random(std::size_t nbits, u64 seed) {
    BitVec out(nbits, 0);
    std::vector<u64>& O = detail::BitVecRawAccess::words(out);
    u64 state = seed;
    for (std::size_t k = 0; k < O.size(); ++k) O[k] = splitmix64(state);
    out.normalize_tail();
    return out;
}

std::string BitVec::to_hex() const {
    static const char* kHex = "0123456789abcdef";
    if (nbits_ == 0) return std::string();

    // 先按整字数拼出完整十六进制串，再**从头部裁掉填充零**。
    // 按 nbits 裁剪的意义：高位对齐、长度恰为 ceil(nbits/4)。若按整字数
    // 输出，nbits 非 64 倍数时会多出前导填充零，且无法区分不同比特长度
    // （4 位 0xA 与 64 位 0xA 会得到同一字符串）。
    std::string full;
    full.reserve(w_.size() * 16);
    for (auto it = w_.rbegin(); it != w_.rend(); ++it) {
        const u64 v = *it;
        for (int nib = 15; nib >= 0; --nib)
            full.push_back(kHex[(v >> (nib * 4)) & 0xFu]);
    }

    const std::size_t nhex    = (nbits_ + 3u) / 4u;
    const std::size_t skipped = full.size() - nhex;   // 头部填充零的个数
    return full.substr(skipped);
}

}  // namespace qkdpp