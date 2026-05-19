# Mg24/USDB BIGSTICK ↔ KSHELL 密度比较

## 设置

- BIGSTICK: `runs/mg24_usdb_density_compare/`，编译用 `make/gfortran-openmp`
- KSHELL: `/home/mengziyan/kshell/run_mg24_usdb_density_compare/`
- 单粒子空间 sd-shell（`docs` 与 KSHELL `usdb.snt`），相同的轨道映射:
  - 1 = p 0d3/2，2 = p 0d5/2，3 = p 1s1/2
  - 4 = n 0d3/2，5 = n 0d5/2，6 = n 1s1/2

## 重跑命令链（独立重跑前清掉旧文件）

```
cd runs/mg24_usdb_density_compare
OMP_NUM_THREADS=1 bigstick-openmp.x < input.bigstick_mg24_usdb_m0_3       # lanczos + 简单密度，M=0
OMP_NUM_THREADS=1 bigstick-openmp.x < input.bigstick_mg24_usdb_m1_2       # lanczos + 简单密度，M=1
OMP_NUM_THREADS=1 bigstick-openmp.x < input.bigstick_mg24_usdb_m0_3_2full # 完整 TBTD M=0
OMP_NUM_THREADS=1 bigstick-openmp.x < input.bigstick_mg24_usdb_m1_2_2full # 完整 TBTD M=1
```

KSHELL：

```
cd /home/mengziyan/kshell/run_mg24_usdb_density_compare
./kshell.exe kshell_mg24_usdb.input > log_Mg24_usdb_m0p.txt
./transit.exe kshell_mg24_usdb_all_density.input > log_Mg24_usdb_all_density.txt
```

KSHELL transit log 含 `*** called mup ***`，表示触发自动 m-up（升 ket，详见 [[kshell_mup]]）。

## 能谱（重跑后完美对齐）

| 态 | J | BIG M=0 | BIG M=1 | KSHELL |
|---|---|---|---|---|
| 1 | 0+ | -87.10445 | — | -87.104449 |
| 2 | 2+ first | -85.60215 | -85.60215 | -85.602148 |
| 3 | 2+ second | -82.98830 | -82.98830 | -82.988299 |

## M 选取规则（不可跳 -999；详见 [[bigstick_m_selection]]）

对每个 (initial, final) 跃迁对挑一个让 BIG 整段无 -999 的 M：

| KSH (ini, fin) | 用 M | BIG 态映射 |
|---|---|---|
| (1,1) | 0 | (1,1) |
| (1,2) (1,3) | 0 | 同 KSH |
| (2,1) (3,1) | 0 | 同 KSH |
| (2,2) | 1 | (1,1) |
| (2,3) | 1 | (1,2) |
| (3,2) | 1 | (2,1) |
| (3,3) | 1 | (2,2) |

验证：使用 M=0 的段和 M=1 的段在用到的 (ini, fin) 块中各自 sentinel = 0。

## 比较脚本（独立写、不依赖旧的 compare_tbtd_*.py）

- `runs/mg24_usdb_density_compare/my_independent_compare.py` — TBTD 比较
- `runs/mg24_usdb_density_compare/my_obtd_compare.py` — OBTD 比较（BIG isospin → KSHELL pn 转换：ρp = (T0+T1)/√2，ρn = (T0−T1)/√2）

## OBTD 结论 — 完全对上

去掉本征矢相位歧义后，286 个常见条目全部 max < 1e-5 量级。最大差 < 5e-6。
**BIG 的 OBTD reduce、CG 约定、isospin → pn 转换均正确。**

## 本征矢相位歧义（Lanczos 自由度）

OBTD ratio 直接揭示：

- BIG M=0 文件：所有 3 个态相位与 KSHELL 同号（ratio = +1）。
- BIG M=1 文件：state 1 与 KSHELL state 2 同号；state 2 与 KSHELL state 3 反号。

→ 比较 TBTD 前对 BIG M=1 中 state-2 的项整体乘 −1（仅修正 Lanczos 相位歧义，非物理）。

## TBTD 结论（共 12720 条 common；max diff = 0.32；mean = 0.019）

去除相位歧义后：

| 类别 | n | max diff | 状态 |
|---|---|---|---|
| (1,1) 0+→0+ | 260 | 2e-7 | ✓ |
| pp/nn same-state Jab=Jcd K∈{0,2,4} | 556 | 6e-6 | ✓ |
| pp/nn same-state Jab=Jcd K∈{1,3} | 280 | 0.25 | ✗ |
| pp/nn same-state Jab≠Jcd K=1 | 240 | 5e-6 | ✓ |
| pp/nn same-state Jab≠Jcd K∈{2,3,4} | 1112 | 0.32 | ✗ |
| pn same-state Jab=Jcd 全 K | 1580 | 6e-6 | ✓ |
| pn same-state Jab≠Jcd 全 K | 1424 | 0.25 | ✗ |
| 所有 cross-state (2,3)/(3,2) 等 | 7528 | 0.31 | ✗ |

正常对上 ~2636 条，失败 ~10084 条。

## Bug 定位（src/bdenslib5.f90）

1. **`couple_2bdensXX`**: K-奇 (Jab=Jcd 时) 或 Jab≠Jcd 任何 K≥2。涉及 direct/HC 双分支、(-1)^{Jab+Jcd-Jtot} 相对相位、pair-swap phase。
2. **`couple_2bdensPN`**: Jab≠Jcd 全 K。该子程序仅有 direct 分支，需查 proton-neutron 对在 Jab≠Jcd 块的对子排列与第三 CG（`cleb(2*Jab, 2*m, 2*Jcd, -2*m, 2*Jtot, 0)`）相关相位。
3. Cross-state 残差需进一步分析：相位修正后 (2,3)/(3,2) 仍 max 0.31，说明耦合公式 bug 也影响 cross-state（不仅相位差）。

## 历史失败实验（已回退、勿重试）

- 将 XX 块从三角 Hermitian 索引切到全非 Hermitian 块索引 — 破坏已对上的 rank-0
- 关闭 XX 块的 Hermitian-conjugate 耦合 — 破坏已对上的 rank-0
- 把 `dens2bflag` 分支重新打开 — 同上
- 合并 M=0 + M=1 后跳 -999 — 概念上违反"−999 表示 M 不可用"原则；正确做法是按 (ini,fin) 选 M

## 已确认正确的 BIGSTICK 部分

- 能谱（与 KSHELL 完全一致）
- OBTD 全部条目（reduce、isospin→pn 转换）
- TBTD 中 PN 块 Jab=Jcd 全 K
- TBTD 中 XX 块 Jab=Jcd 偶 K
- TBTD 中 XX 块 Jab≠Jcd K=1
- `src/bdenslib4.f90` 中的 PN OpenMP race fix（用 `reduction(+:dmatpn)` 而非 `dmatpnhc`）
