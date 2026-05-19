# TBTD Bug 定位分析

## 现状

- OBTD 全对，BIG 单体公式正确
- TBTD 失败模式：
  - XX 块 same-state Jab=Jcd K-奇：错 (PN 块此情形对)
  - XX 块 Jab≠Jcd K≥2：错
  - PN 块 Jab≠Jcd：错
  - cross-state 全部混合错（含本征矢相位与公式差异）

## XX 块 HC 分支机制（src/bdenslib5.f90 line 388-477）

BIGSTICK `couple_2bdensXX` 内层 loop 限制 `cpair <= dpair`（上三角 m-pair）：

- **Direct branch (cpair ≤ dpair)**: 用 `v(itbme) = dmatpp(itbme)`，填 `x2bden(indx_direct)%v(Jab, Jcd, Jtot)`. 第三 CG = `cleb(2*Jab, 2*m, 2*Jcd, -2*m, 2*Jtot, 0)`，时反相位 `(-1)^{Jcd-m}`。
- **HC branch (cpair < dpair)**: 用 `vhc(itbme) = dmatpphc(itbme)`，填 **transposed** 入口 `x2bden(indx_HC)%v(Jcd, Jab, Jtot)`. 第三 CG = `cleb(2*Jcd, -2*m, 2*Jab, 2*m, 2*Jtot, 0)`，时反相位 `(-1)^{Jab-m}`。

由 CG m-反演 $\langle J_0,-m,J_0,m|K,0\rangle = (-1)^{2J_0-K}\langle J_0,m,J_0,-m|K,0\rangle$，对 J0 整数：

$$\text{HC's 3rd CG} = (-1)^K \cdot \text{Direct's 3rd CG}$$

## 相同-orbital Jab=Jcd 入口分析

对 same-orbital (a=b=c=d 同 orbit)，`indx_direct = indx_HC`，且 Jab=Jcd 时 `%v(Jab,Jcd)=%v(Jcd,Jab)`。**两分支填同一入口**。

对实波函数 same-state diag: $dmatpp = dmatpphc$（Hermitian 对称）。

对每个 (cpair=p1 < dpair=p2):
- Direct + HC = $v \cdot C_1 + vhc \cdot C_1 \cdot (-1)^K = C_1 \cdot v \cdot [1 + (-1)^K]$
- K-偶: $2 v C_1$ ✓ (匹配 KSHELL 该 m-pair 对的贡献 = $v(p_2,p_1) C_1 + v(p_1,p_2) C_1 = 2 v C_1$)
- K-奇: $0$ ✗ **BIG 丢失 cpair<dpair 部分的贡献**，只剩 cpair=dpair 对角 m-pair 部分

KSHELL 通过 `pair_sum` 中 `do ij` `do kl` 双独立求和，自然遍历所有 (ij, kl) 组合（包括 ij<kl 和 ij>kl），不需要 HC 分支补全。

## KSHELL pair 存储约定 (operator_mscheme.f90:789)

```fortran
if (k1==k2 .and. m1inv(m1)>=m2inv(m2)) cycle
```

对 same orbital (k1=k2) 仅存 m1<m2 的 ordered pair。对 different orbitals 不限制。**KSHELL 双独立 sum 把每个 (ij, kl) 组合算一次，自动覆盖 ij<kl, ij=kl, ij>kl 全部。**

## 实验性数值证据

state 2 same-state, (a=b=c=d=orbit 2 (d5/2), Jab=Jcd=2):
- K=0: ratio = 1.000 ✓
- K=1: BIG=0.0951, KSH=0.000063 (≈0 by Hermiticity) ✗
- K=2: ratio = 1.000 ✓
- K=3: BIG=0.0508, KSH=0.0319, ratio=1.59 ✗
- K=4: ratio = 1.000 ✓

**模式**：偶 K 完美对；奇 K BIG 与 KSH 都不一致。具体来说：
- K=1 BIG 非零、KSH ≈ 0：BIG 有了不该有的贡献？或 KSH 物理对而 BIG 错？
- K=3 BIG 与 KSH 都非零但比例 ≈ 1.59：两边都不为零，可能都对（差归一化）或都错（不同方式）。

## 待回答的关键问题

1. **物理上 K-奇 same-orbital Jab=Jcd same-state diag 是否必为零？** 我的 Hermitian 论证得"是"，但 KSHELL K=3 输出非零（0.0319）说明论证可能漏一个相位。需要严格检查 $(O^K)^\dagger$ 的形式。
2. **BIG HC 分支的 $(-1)^K$ 抵消是 bug 还是 feature？** 它在 K-偶时给出正确双倍贡献（与 KSHELL 一致），但 K-奇时却消去——这暗示 BIG HC 公式是按"K-偶应双倍 / K-奇应消去"约定设计的，可能基于上文 Hermiticity 论证。如果论证对，BIG 是对的，KSHELL 多了不该有的贡献。
3. **谁对**：需要构造解析可验证情形（如最简模型 18F 或 18O sd-shell 双体配对）手算 TBTD，对比 BIG / KSHELL。

## 修复方向（待 #2 / #3 决定后）

**如果 BIG 错**（K-奇应有贡献而 BIG 因 HC 抵消把它消掉）：
- HC 第三 CG 改为 `cleb(2*Jab, 2*m, 2*Jcd, -2*m, 2*Jtot, 0)`（与 direct 同），phase 改 `(-1)^{Jcd-m}`。
- 同时需检查 different-orbital 情况（HC 填 TRANSPOSED 入口）是否仍正确。

**如果 KSHELL 错**（同 orbital K-奇 不该非零）：
- BIG 当前公式对，KSHELL `get_cpld_tbtd:pair_sum` 在某些情形下应对 (ij, kl) 求和加额外限制（如 ij≤kl + Hermitian 反对称化），但目前是全 (ij, kl) 双独立 sum。

## 下一步

不预设谁对：构造一个 18F 或 4-particle 简单 sd-shell 模型，手算 TBTD 验证。或者用 BIG/KSHELL 同时跑 N=2 single-orbit 双体配对模型对比。

## 18O 解析判决（已完成）

取 18O (Z=0, N=2)，BIG/KSH 能谱完美对齐（-11.93/-9.93/-8.40/-7.57/-7.34）。看 state 2+ 同态对角 (a=b=c=d=nd5/2 即 orbit 5) Jab=Jcd=2 入口对 K=0..4：

| K | BIG M=1 | KSHELL | 解析（纯 (d5/2)² J=2 假设）|
|---|---------|--------|---|
| 0 | 0.6056080 | 0.6056080 | 1.000 |
| 1 | 0.3584212 | 0.6056080 | 1.000 |
| 2 | 0.6056080 | 0.6056080 | 1.000 |
| 3 | 0.3584212 | 0.6056080 | 1.000 |
| 4 | 0.6056081 | 0.6056080 | 1.000 |

解析推导：单 J=2 pair 态 $|i\rangle = -\hat{A}^\dagger(aa,2,0)|\text{core}\rangle/\sqrt{2}$。算符 $T^K_q = [\hat{A}^\dagger \otimes \tilde{A}]^K_q$ 作用：$T^K_q|2 M_i\rangle = \pm c_d c_c \cdot CG(2,M_i+q,2,-M_i|K,q)|2,M_i+q\rangle$，由 Wigner-Eckart 得 reduced ME，再按 BIGSTICK manual Eq. 5.23 取 TBTD：

$$\text{TBTD} = \frac{\langle 2\|T^K\|2\rangle}{[K]\sqrt{(1+\delta_{aa})^2}} = 1 \text{ for all } K$$

（恒等 1 因为纯单 pair；c_d c_c = 2 与 [K]·2 在每个 K 抵消）。

KSHELL 数值 0.6056（全 K 同）= 解析×缩放因子 0.6056（18O 2+ 并非完全纯单 pair，是 (d5/2)²+(d3/2)(s1/2)+s1/2² 等混合，配对幅度平方≈0.6056）。

**结论：KSHELL 物理对。BIG K-奇错误**。

之前的 Hermiticity 论证错——$\langle 2 0|T^K_0|2 0\rangle = CG(2,0,2,0|K,0)\cdot D=0$ 对 K 奇仅说明 M=0,q=0 矩阵元零，但 reduced ME 非零（其他 M-choice 可见之）。

## Bug 确认与修复方向

`src/bdenslib5.f90` `couple_2bdensXX` HC 分支（line 459-477）：当前用 `cleb(2*Jcd, -2*m, 2*Jab, 2*m, 2*Jtot, 0)` 与 phase `(-1)^{Jab-m}`。对 same-orbital Jab=Jcd 入口，相对 direct 多一个 $(-1)^{J_{tot}}$ 因子，K-奇时让 off-diagonal m-pair 贡献抵消。

**修复策略 A**：HC 分支用与 direct 一致的 CG `cleb(2*Jab, 2*m, 2*Jcd, -2*m, 2*Jtot, 0)` 与 phase `(-1)^{Jcd-m}`，但只对 vhc 使用。需小心 Jab≠Jcd 与 different-orbital 情况下 HC 填 transposed 入口的正确性。

**修复策略 B**：保留 HC 当前 CG，但在最终累加时给 HC 项额外乘 $(-1)^{Jab+Jcd-Jtot}$ 相位以抵消相位反演。等价于策略 A 但只改一个相位常数。

需要修改后重跑 Mg24 + 18O 全部对比 BIG vs KSHELL 验证：
- 18O 同态 K-奇 应回到 0.6056080
- Mg24 各失败类别（XX K-奇、XX Jab≠Jcd、PN Jab≠Jcd、cross）的影响。
- PN 块的 bug 可能是不同来源（PN 无 HC 分支），需另行排查。

## 修复实施记录（XX HC 分支，已完成）

`src/bdenslib5.f90` `couple_2bdensXX` HC 分支（line ~469-470 修改后）：
```fortran
* (-1)**(Jab -m)   & ! PHASE FROM time-reversal of second pair
* (-1)**(Jab + Jcd - Jtot)  & ! BUGFIX: compensate (-1)^{Jab+Jcd-Jtot} from CG m-reversal in HC's 3rd CG vs direct
```

## 修复后验证结果

### 18O (纯 nn 块)
- state 2+ self-diag (5,5,5,5) Jab=Jcd=2 K=0..4 全部 = 0.6056080，匹配 KSHELL ✓
- state 2+ → 2+ 全 125 个 common 条目 max diff = 3e-7 ✓

### Mg24 (含 pp/nn 和 pn 块)
完美对上 (max < 1e-5):
- SAME pp/nn Jab=Jcd 全 K：max 5.8e-6 ✓ (修复前 K-奇 max 0.25)
- SAME pp/nn Jab≠Jcd 全 K：max 8.1e-6 ✓ (修复前 K≥2 max 0.32)
- SAME pn Jab=Jcd 全 K：max 5.6e-6 ✓ (修复前已对)
- (1,1)=0+→0+ 整组：max 2e-7 ✓

仍未对上 (修复未触及):
- SAME pn Jab≠Jcd：max 0.253 ✗ (PN 块 bug，与 HC 无关)
- CROSS pp/nn Jab=Jcd：max 0.154 ✗
- CROSS pp/nn Jab≠Jcd：max 0.203 ✗
- CROSS pn Jab=Jcd：max 0.112 ✗
- CROSS pn Jab≠Jcd：max 0.193 ✗

### 总体改善
Mg24 max diff: 0.61 → 0.25, mean: 0.026 → 0.015.
12720 条目中 5616 个 (44.2%) 现在 |diff| < 1e-4。
3768 条目 (29.6%) 是修复直接解决的 SAME-state XX 类别。

## PN 块 Jab≠Jcd Bug — 已定位与修复

**Bug 在 `print_out_2bdens` 中 PN 块 indx 计算与 `couple_2bdensPN` 不一致**：

- `couple_2bdensPN` line 770（compute）: `indx = pairblocksize*(pair1-pcref-1) + pair2-pcref + pcstart`，其中 `pair1 = destruction orbital pair, pair2 = creation orbital pair` — **(destruction, creation) layout**
- `print_out_2bdens` line 1170（print, 修复前）: `indx = (abcouple - ...)*parblocksize + cdcouple - ...`，其中 `abcouple = creation, cdcouple = destruction` — **(creation, destruction) layout**

二者互为转置！Print 读到的是 compute 写到"转置位置"的数据。对 Jab=Jcd 入口不显现（同一 slot），对 Jab≠Jcd 入口表现为 BIG[orig] = KSH[transp]，BIG[transp] = KSH[orig]。

对比 XX：`couple_2bdensXX` line 303 用 `indx = (dcouple-pcref-1)*coupleblocksize + ccouple-...` (destruction, creation)；XX print 也用 `(cdcouple, abcouple)` (destruction, creation) — 一致 ✓。所以 XX 没此问题。

**修复**（已实施）：把 PN print 的 indx 公式中 `abcouple ↔ cdcouple` 交换，对齐 compute 的 (destruction, creation) 约定。

修复后效果（Mg24）：
- SAME pn Jab≠Jcd: max 0.253 → **5.7e-6** ✓
- CROSS pn Jab=Jcd: max 0.112 → **5.8e-6** ✓
- CROSS pn Jab≠Jcd: max 0.193 → **4.6e-6** ✓

## 当前状态总结

完成两个修复后，Mg24 状态（含 PN print 修复）：

| 类别 | 状态 | max diff | n |
|---|---|---|---|
| SAME pp/nn Jab=Jcd | ✓ | 5.8e-6 | 836 |
| SAME pp/nn Jab≠Jcd | ✓ | 8.1e-6 | 1352 |
| SAME pn Jab=Jcd | ✓ | 5.6e-6 | 1580 |
| SAME pn Jab≠Jcd | ✓ | 5.7e-6 | 2904 |
| CROSS pn Jab=Jcd | ✓ | 5.8e-6 | 2048 |
| CROSS pn Jab≠Jcd | ✓ | 4.6e-6 | 4568 |
| **CROSS pp/nn Jab=Jcd** | ✗ | 0.154 | 1040 |
| **CROSS pp/nn Jab≠Jcd** | ✗ | 0.203 | 2200 |

总体: max 0.61 → **0.20**, mean 0.026 → **0.0045**, **83.1% 条目** |diff| < 1e-4.

## 剩余 CROSS pp/nn Bug（深入分析）

PN 块 print indx bug 修复后（详见上节），PN cross 已对（max 5.8e-6 ✓）。**剩余仅 XX cross-state pp/nn**，max 0.20。

### 详细推导（HC 修复对 cross-state 也应正确）

对 orig 入口 (a,b,c,d,Jab,Jcd) 的贡献分解：
- 每对 (P_ab m-pair index, P_cd m-pair index) 组合唯一访问一次 iter
- If P_ab < P_cd → iter (cpair=P_ab, dpair=P_cd)：direct 贡献 orig 入口
- If P_cd < P_ab → iter (cpair=P_cd, dpair=P_ab)：HC 贡献 orig 入口

每组合 m-scheme 算符是 $O_1 = a^+_{m_a}a^+_{m_b}a_{m_d}a_{m_c}$：
- direct 用 $v(\text{iter1}, f→i) = ⟨f|O_1|i⟩$
- HC 用 $v_{hc}(\text{iter2}, f→i) = ⟨f|\text{backward of iter2}|i⟩$。Forward op of iter2 = $a^+_{m_c}a^+_{m_d}a_{m_b}a_{m_a} = O_1^\dagger$。Backward = $O_1$。所以 $v_{hc}(\text{iter2}) = ⟨f|O_1|i⟩$。

数值上 $v(\text{iter1}) = v_{hc}(\text{iter2})$。CG 因子在 my-fix 后等于 direct CG 因子（公式推导：HC 第三 CG 经 m-反演 = $(-1)^{Jab+Jcd-K}$ × direct 第三 CG；my-fix 乘 $(-1)^{Jab+Jcd-K}$ 抵消）。

所以两分支对每个 (P_ab, P_cd) 组合贡献相同的 $⟨f|O_1|i⟩ × \text{CG}$，cover 全部 m-channels 各一次 — **理论上应得到正确 reduced ME**。

### 但实证 0.20 残差 — 矛盾来自何处？

Top 条目呈 Hermitian 反对称模式：`BIG_orig - BIG_transp ≈ KSH_orig - KSH_transp`，但 `BIG_orig + BIG_transp ≠ KSH_orig + KSH_transp`。即 BIG 和 KSH 仅在 **对称组合**（orig + transp）上有差异。

可能性：
1. **Lanczos phase**: state 2 vs state 3 在 M=1 file 中相位歧义。已用 OBTD 比例确认 state 1: +1, state 2: -1。Lanczos 相位是 per-eigenvector global，不该 per-matrix-element 不同。已排除。
2. **KSHELL 卷积公式**（用 OBTD 乘积合成 TBTD）在 cross-state 下可能与 BIG 直接公式有 ±1 不同的对称化处理 — KSHELL 可能错。
3. **BIG 内部其他相位**（phasekl, phaseij, hfactor）在 HC 分支中与 direct 不对称（已 verify 相同）— 排除。
4. **算符约定差异**：BIG manual ρ = -⟨f‖[A^+ × Ã]^K‖i⟩/norm；KSHELL/Brown ρ = +⟨…⟩/norm。在 cross-state 下，∓ 符号与 orbital-transpose 对称性互动 — 可能源。

### 实测数据（Mg24, 后所有已知修复）

```
=== Per (ini, fin) state pair ===
(1->1): max 7.8e-7    (2->2): max 5.5e-6    (3->3): max 5.4e-6  [diag]
(1->2): max 4.0e-6    (1->3): max 4.6e-6    (2->1): max 6.5e-6  [0+↔2+ XX]
(2->3): max 0.20254   (3->2): max 0.20254                       [2+↔2+' XX]

=== By block ===
所有 PN 同/跨态:        max ~5e-6 ✓
XX SAME state pp/nn:    max ~8e-6 ✓
XX CROSS pp/nn K=0..4:  max 0.155-0.203 ✗
```

约 20% 条目（XX cross-state pp/nn）残留 0.20 max，所有其他类别 ≤ 5e-6 ✓。

### 下一步候选

- 用 sympy 在 sd-shell 内做 Mg24 一个具体 cross 入口的 KSHELL `pair_sum_decompose` 公式 vs BIG 直接公式数值对算，硬判一谁对
- 或写小型独立 m-scheme 直接积分代码，对 single cross-state 入口求 reduced ME，与 BIG/KSH 比较
- 考查 BIG manual ρ 定义里的负号与 KSHELL/Brown 公式的符号一致性，看 cross-state 在两边 ∓ 因子处理是否对称破缺
