# BIGSTICK (Public branch)

大规模壳模型对角化代码（Fortran），Calvin W. Johnson 等。本仓库当前重点工作是排查 BIGSTICK 的两体跃迁密度 (TBTD) 与 KSHELL 输出在某些耦合分支上的不一致。

## 目录结构

- `src/` — Fortran 主代码
- `docs/` — 手册与教程 (`BIGSTICK_Manual.pdf` 中 Applications 节有 OBTD/TBTD 公式)
- `bin/` — 编译后的可执行文件 (`bigstick-openmp.x`)
- `make/` — Makefile
- `examples/`、`test/`、`util/` — 例子、单元测试与工具
- `runs/mg24_usdb_density_compare/` — 当前 Mg24/USDB 的 BIG↔KSHELL 比较工作目录（输入、密度输出、独立比较脚本）
- 平级仓库 `../kshell/` — KSHELL，参考实现 (`src/transit.F90`, `src/operator_mscheme.f90`, `Brown..._density.pdf`)

## 主要模块（密度部分）

- 单体密度耦合 `coupled_densities` in `src/bdensities.f90`
- 两体密度 m-scheme 累加 `src/bdenslib4.f90`（pp/nn 在 `applyhPPbundled_den`、`applyhNNbundled_den`；pn 在 `applyhPNbundled_den`；含 OpenMP 归约）
- 两体密度 J-coupling `src/bdenslib5.f90`（XX 在 `couple_2bdensXX`；PN 在 `couple_2bdensPN`；输出 `print_out_2bdens`）
- 其他: `bdenslib1..3,6.f90`, `bdensmathlib.f90`, `bdensities.f90`

## 当前工作背景与发现

Mg24/USDB 三个最低态 (0+, 2+, 2+) 与 KSHELL 完整比对（OBTD + TBTD）。详细发现与定位见：

- 比较方法、相位约定、Bug 模式 → @.claude/mg24_usdb_density_compare.md
- 公式与代码对照 → @.claude/tbtd_formulas.md
- TBTD bug 机制分析（XX HC 分支 K-奇问题）→ @.claude/bug_localization.md
- KSHELL 自动 m-up 机制 → @.claude/kshell_mup.md
- BIGSTICK -999 哨兵与 M 选择规则 → @.claude/bigstick_m_selection.md

进度详见 `TODO.md`。
