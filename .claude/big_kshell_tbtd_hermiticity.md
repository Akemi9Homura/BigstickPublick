# BIG vs KSH TBTD: Hermiticity 与 per-entry 全核验证

**重要更正（2026-05-20）**：本文件之前版本记录"BIGSTICK TBTD per-entry 在 cross-state Hermitian-反对称分量上还有 0.15 bug"是**错的**。那个结论来自一个错误的比较方法（没按 (ini,fin) route M-projection）。现在正确做法 + commit 5670a1c 的 "full non-Hermitian same-species density path" 修复后，**所有测过的核 (Mg24/USDB, Mg24/IMSRG, Na21/USDB, O18/USDB) BIG 与 KSH TBTD 都逐条对得上 ~1e-5**。

## Hermiticity 关系（已验证）

对实波函数 + 实算符：
$$\rho^{fi}_K(ab,J_{ab};cd,J_{cd}) = (-1)^{J_f - J_i + J_{cd} - J_{ab}} \, \rho^{if}_K(cd,J_{cd};ab,J_{ab})$$

整数 J 下相对的 $(-1)^{2J_{ab}}=1$，公式等价 $(-1)^{J_f-J_i+J_{ab}+J_{cd}}$；半整数 J 下两形式仍等价。**相位不含 K**。

KSH 与 BIG 都严格满足。Mg24/IMSRG 5态包含 3+（Jf-Ji 半整数倍）情况下也都验过。

## per-entry 比较：正确做法

每个 (ini_KSH, fin_KSH) 对必须 route 到 reduce CG 非零的 M 选择，否则 BIGSTICK 输出 -999 哨兵或 -706 scaled 哨兵。KSHELL 内部 m-up 机制（`*** called mup ***`）会自动升 ket 解决，但 BIGSTICK 不会。

| 案例 | 全 0+, 2+ M=0 OK? | 解决 |
|---|---|---|
| 0+ ↔ 任意 | 是（M=0 reduce CG 全非零）| 用 M=0 run |
| 2+ ↔ 2+ 偶 K | M=0 OK | 用 M=0 run |
| 2+ ↔ 2+ 奇 K | M=0 出 −999（reduce CG=0）| **必须用 M=1 (2Jz=2) run** |

奇 A（如 Na21 4 个低态 3/2+, 5/2+, 7/2+, 9/2+）单一 M=1/2 (2Jz=1) run 没有 reduce CG 零，直接全部对得上。

## 当前验证状态（5670a1c 修复后）

| 核 + H | states | M routing | n common | max\|Δ\| | mean\|Δ\| |
|---|---|---|---|---|---|
| Mg24/USDB | 0+, 2+_1, 2+_2 | M=0 含 0+；M=1 (2+↔2+) | 16528 | 8.1e-6 | 3.8e-7 |
| Mg24/IMSRG | 0+, 2+, 2+, 4+, 3+ | M=0 含 0+；M=1 其他 | 63042 | 1.5e-5 | 2.8e-7 |
| Na21/USDB | 3/2+, 5/2+, 7/2+, 9/2+ | 单 M=1/2 | 53189 | 2.6e-5 | 3.4e-7 |

精度 1e-5 ~ 几 e-5 = KSHELL transit log 7 位打印精度。

## commit 5670a1c 关键修复

`src/bdenslib4.f90` `applyhPPbundled_den` / `applyhNNbundled_den` 的 'b' (backward) 分支：当 `dens2bflag = .true.` 时整段 `cycle` 跳过，forward 分支同时把 `dmatpp` 与 `dmatpphc` 都纳入 OMP `reduction`，并把 backward 的 ME 直接累加到 `dmatpp`（非 HC 矩阵）。

`src/bdenslib5.f90` `couple_2bdensXX`：dens2bflag 模式下 `pairblocksize = XX2(it)%block(m,par)`（全 block 而非 triangular），循环 `cpair = 1, npairXX(it)` 全 m-pair 遍历，**不限 cpair ≤ dpair**；HC 分支整段 `cycle` 跳过。

物理含义：dens2bflag 模式把 same-species 密度按完全非 Hermitian 的 full block 算，每个 (cpair, dpair) m-pair 组合恰好访问一次，避免了 triangular indexing + HC 分支带来的 cross-state 相位错。

## 比较脚本

- `runs/mg24_usdb_density_compare/my_independent_compare.py` — 3 state, route M=0/M=1
- `runs/mg24_imsrg/compare_tbtd_imsrg_5_omp4.py` — 5 state IMSRG, route + sign search
- `runs/na21_usdb/compare_tbtd_na21.py` — odd-A, 单 M=1/2

## 启示

1. **不要相信旧的"per-entry 反对称分量 bug"结论** — 那是错的比较方法 + 没用 commit 5670a1c 的二进制
2. KSHELL m-up 让 KSHELL 一份 wfn 文件能输出所有 (Ji, Jf) 跃迁；BIGSTICK 没此机制，所以必须按 (ini,fin) route 选合适的 M
3. 同 M 选择正确后，BIGSTICK TBTD per-entry 与 KSHELL 一致到 transit log 打印精度，三核（偶 A 0+/2+/4+/3+ 含 USDB & IMSRG；奇 A Na21 USDB）共测得 max ≤ 3e-5
