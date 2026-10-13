# libqkd-postproc

QKD 后处理链的 C++ 参考实现 —— MP-QKD-TP 协议栈（旗舰②）。

对应论文 *MP-QKD-TP* 的：

| 论文位置 | 原文陈述 | 本仓库对应 |
|---|---|---|
| §IV-A | "error-corrected via the **Cascade protocol**, and privacy-amplified via **Toeplitz-matrix hashing**" | `cascade.*`, `toeplitz.*` |
| §V-B | "After screening, trusted paths undergo **Cascade error correction and privacy amplification**" | 同上 |
| §IV-A 式 | `m = n·(1 − h₂(q)) − n·h₂(q)·f_ec` | `privacy_amplify()` |

> **背景**：论文将 Cascade 与 Toeplitz 隐私放大作为已实现的管线环节陈述，但仿真代码
> `sim/netsquid/qkd_protocol.py` 中 Cascade 仅有空的 `CascadeConfig` 配置类（用随机比特顶替），
> Toeplitz 则是逐位 Python 循环。本仓库为二者的可运行参考实现。

---

## 构建

要求：C++17 编译器、CMake ≥ 3.16。

```bash
cmake -S . -B build -G "MinGW Makefiles" \
      -DCMAKE_CXX_COMPILER="<path>/g++.exe" \
      -DCMAKE_MAKE_PROGRAM="<path>/mingw32-make.exe"
cmake --build build
```

无 CMake 时可直接编译（Windows 需链接 `bcrypt` 以使用系统熵源）：

```bash
g++ -std=c++17 -O3 -mpopcnt -mpclmul -Iinclude \
    src/bitvec.cpp src/toeplitz.cpp src/cascade.cpp src/gf2poly.cpp \
    tests/test_main.cpp -lbcrypt -o build/bin/qkdpp_tests.exe
```

`-mpclmul` 只决定「可否发射 PCLMULQDQ 指令」；是否**真正启用**由 `src/gf2poly.cpp`
的运行时 CPUID 探测决定（两者同时成立才走硬件路径）。因此在以 `-mpclmul` 编译、
却运行于不支持 PCLMULQDQ 的 CPU 上时，会自动回落到软件实现而**不会** SIGILL。
`gf2_has_clmul()` 返回的是这个「实际在用」的运行时判定结果。

---

## 模块

| 文件 | 职责 |
|---|---|
| `bitvec.*` | 紧凑位向量（uint64 打包），支持 GF(2) 内积、Hamming 距离、批量窗口切片 |
| `entropy.hpp` | 二元熵 h₂(q) 与 `clamp`，均以 `[[nodiscard]]` + NaN 传播表达错误 |
| `toeplitz.*` | Toeplitz 隐私放大：`toeplitz_reference`（朴素教学用）、`toeplitz_multiply`（滑动窗 O(m·n/64)）、`privacy_amplify`、`random_seed_from_os` |
| `gf2poly.*` | GF(2) 多项式乘法内核：Karatsuba + PCLMULQDQ（含运行时 CPUID 探测），`toeplitz_multiply_fast` O(N^1.585) |
| `cascade.*` | Cascade 纠错（Fisher–Yates 置换 + 分块奇偶校验 + 二分定位） |

### API 安全性设计

隐私放大路径按「实验」与「生产」明确分离，避免误用：

| 接口 | 用途 | 密钥空间 |
|---|---|---|
| `privacy_amplify_deterministic(sifted, qber, u64_seed, f_ec)` | **仅实验**，可复现 | ≤ 2^64（`splitmix64` 派生） |
| `privacy_amplify(sifted, qber, BitVec_seed, f_ec)` | **生产**，种子由调用方从 CSPRNG 提供 | 完整 n+m−1 位 |
| `random_seed_from_os(bits)` | 平台级 CSPRNG：Windows `BCryptGenRandom` / Linux `getrandom` / POSIX `/dev/urandom` | — |

三条硬约束（编译器强制，非注释约定）：

1. **无默认参数**：两个 `privacy_amplify*` 的 `ec_efficiency` 都必须显式给出。
   旧版 `privacy_amplify(sifted, qber, u64)` 这样的 3 实参调用会静默落到确定性路径，
   把密钥空间从完整种子向量降级到 2^64。去掉默认参数后这类调用直接编译失败。
2. **实验路径改名**：`privacy_amplify` → `privacy_amplify_deterministic`，函数名即告知语义。
3. **入参校验前置**：`qber` / `ec_efficiency` 的有限性与取值范围在**任何内存分配之前**校验，
   两个重载行为一致。`ec_efficiency <= 0` 会让 m 不降反升、静默跳过隐私放大，故直接拒绝。

### 关键恒等式

Toeplitz 矩阵乘向量等价于一次多项式卷积：

```
y[i] = (seed * x)[i + n − 1],   i = 0 … m−1
```

因此 `toeplitz_multiply_fast` 用 Karatsuba 多项式乘法把复杂度从 O(m·n/64) 降到 O(N^1.585)。

---

## 验证结果

```
=== libqkd-postproc unit tests ===
test_bitvec ...
test_toeplitz_equivalence ...
test_toeplitz_linearity ...
test_gf2_polymul ...
    hardware CLMUL: yes
test_toeplitz_fast ...
test_privacy_amplify ...
    n=10000 q=0.03 -> m=5917 (Shannon cap 8056.1)
test_size_overflow_guards ...
test_toeplitz_explicit_seed ...
test_entropy_nan_propagation ...
test_cascade_short_key ...
test_cascade ...
    injected=245  corrected=245  residual=0
    parity_sent=2035 bits (0.2484 n)   theoretical_leak=1751.7 bits (0.2138 n)   ratio=1.16

--- 555 checks, 0 failures ---
ALL TESTS PASSED
```

覆盖点：

- `slice` 的按字批量实现与逐位参考在**所有**字对齐/非对齐起点下逐位一致，且尾部无效位为 0
- `toeplitz_multiply` / `toeplitz_multiply_fast` 与 `toeplitz_reference` **逐位一致**
  （小尺寸全量；n=20000, m=12000 抽样 64 点复核）
- `gf2_polymul` 小尺寸与朴素卷积逐位一致（跨越 4096 Karatsuba 阈值）；
  大尺寸改用代数不变量（乘 1 恒等、线性性）——避免 O(na·nb) 的 2.5×10⁸ 次内层循环
- Toeplitz 线性性 `T(x₁ ⊕ x₂) = Tx₁ ⊕ Tx₂`
- 相同输入+种子可复现；不同种子雪崩；m 随 q 单调递减
- **二元熵的错误语义**：`h2` 只在 `q == 0` / `q == 1` 两端点返回 0；
  `q < 0`、`q > 1`、`±inf` 一律返回 NaN。这是安全关键 —— `h2` 出现在
  `m = n(1−h2) − n·h2·f_ec` 中，`h2 = 0` 意味着**不扣减任何泄露**、
  密钥率被高估，而 QBER 估计出错恰好会落到越界分支。
- `privacy_amplify*` 拒绝 NaN `qber`、NaN `ec_efficiency`、
  `ec_efficiency <= 0`、越界 `qber`
- `Cascade` 构造函数拒绝 `q_est ∈ {0, NaN, ±inf, >1}`；`n=1` 抛异常；`n=0` 返回零结果
- Cascade 注入 245 个错误**全部纠正、残余为 0**，泄露 0.2484n 贴近理论下界 0.2138n（比值 1.16）

## 性能

QBER = 3%、f_ec = 1.1 ⇒ m/n = 0.5917。

| n (输入比特) | m (输出比特) | [A] 滑动窗 O(m·n/64) | [B] Karatsuba+CLMUL | 加速比 |
|---|---|---|---|---|
| 50 000 | 29 588 | 0.0287 s | 0.000595 s | **48×** |
| 200 000 | 118 355 | — | 0.00545 s | — |
| 1 000 000 | 591 777 | — | 0.0755 s | — |
| 4 000 000 | 2 367 108 | — | 0.688 s | — |
| 8 000 000 | 4 734 216 | — | 2.077 s | — |

> n = 8×10⁶ 对应 SKR ≈ 7.6 Mbps 累积一秒的 sifted key 量，用本实现处理耗时约 2.1 s。
> 作为对比，Python 参考实现（`qkd_protocol.py::_privacy_amplification`）对 n = 10⁶
> 需执行约 5.9×10⁵ 次 numpy 切片归约，量级为数百秒。

`bench/bench_toeplitz.cpp` 在吞吐测量**之外**还做一次**逐位交叉校验**
（`toeplitz_multiply` vs `toeplitz_multiply_fast`），不一致则以退出码 1 报告。
这是必要的：若快速路径存在正确性缺陷，漂亮的 Mbit/s 数据会完全将其掩盖，
而 bench 是唯一同时持有两条路径输出的场所。实测 4 组尺寸（n = 64 … 200000）
全部 `consistent`。

---

## 已知限制

1. **Cascade 未实现"级联"（cascading）**：当前为多轮独立分块 + 二分定位
   （即 "Cascade without cascading"）。在测试工况下已能完全收敛，但对高 QBER 场景
   的泄露效率不如完整版。完整版需维护跨轮块重叠关系图。
2. **`bitvec.hpp` 无可变长 resize**：Karatsuba 通过 `extract_bits` 生成子多项式，
   带来额外分配。后续可用 slab 分配器优化。
3. **未做 FFT/NTT**：n > 10⁷ 时可考虑把 Karatsuba 换成 GF(2) FFT（additive FFT），
   复杂度可进一步降到 O(N log N)。
4. **非恒定时间**：当前实现面向吞吐与正确性，未做侧信道防护。
   生产环境若需抵抗时序攻击，应对分支与表访问做遮蔽。

## 后续路线

- [ ] Cascade 完整级联（提升高 QBER 下的泄露效率）
- [ ] 与 Python 侧的 golden 值交叉验证（`to_hex()` 已为此预留）
- [ ] pybind11 绑定，接回 MP-QKD-TP 仿真管线
- [ ] GF(2) additive FFT（面向 n > 10⁷）
- [ ] 恒定时间实现（面向侧信道对抗场景）
