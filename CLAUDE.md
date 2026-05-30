# BIGSTICK (Public branch)

大规模壳模型对角化代码（Fortran），Calvin W. Johnson 等。本仓库当前主要用于 Li8 等核的 NO2B no-core 壳模型（NCSM）能谱计算；历史上还做过 BIG↔KSHELL 两体跃迁密度（TBTD）验证（见末尾）。

## 目录结构

- `src/` — Fortran 主代码
- `docs/` — 手册与教程（`BIGSTICK_Manual.pdf`）
- `bin/` — 编译后的可执行文件（`bigstick-openmp.x`）
- `make/` — Makefile（用 `make -C make gfortran-openmp` 编译；本环境无 `ifort`）
- `examples/`、`test/`、`util/` — 例子、单元测试与工具
- `runs/` — 运行工作目录（见下方 Li8 约定）
- 平级仓库 `../kshell/` — KSHELL 参考实现

## NO2B 二进制相互作用输入（format code `no2`）

- NO2B 二进制 reader 从平级 `BigstickPublick-3N` 移植而来，实现在 `src/binput_tbme.f90`。
- 输入：平级 `normal-order/src/no2b` 程序写出的 `no2b_*.bin` minipack 二进制。
- reader 以 unformatted stream 读：header（hw, Emax, orbit 数, 单体数, 两体数）、orbit 表（n,l,2j,2tz）、零体项、单体 ME（折进有效两体）、归一化两体 ME（填进 `ppme`/`pnme`/`nnme`）。
- BIGSTICK 交互输入用法：① 输入 `no2` ② 输入 `.bin` 文件名 ③ 按提示给振子频率与质心强度。
- **解析坑**：`tbme_menu` 只看文件名前三字符，`no2b_*.bin` 开头的 `no2` 会被当成另一个格式码。每个运行目录里把源二进制软链成 `li8_hw<hw>_no2b.bin` 这类名字，再把这个名字喂给 BIGSTICK。

## NO2B NCSM 运行约定（通用）

适用于任意核的 NO2B no-core 壳模型能谱扫描。下面以 Li8（Z=3, N=5）+ `EM1.8_2.0`/`DNNLO_go394` 核力为例，换核只改 `<nucleus>`、`<Z> <N>` 与 NO2B 源文件。

### 目录与命名

- 目录布局：`runs/<nucleus>/<interaction>_<emax>_<e3max>/<hw_Nmax_beta>/{input,output}`
  - 例：`runs/Li8/EM1.8_2.0_emax12_e3max16/hw16_Nmax2_beta0`
- 路径/输入名/输出名/Slurm 日志一律用**完整** interaction 标签（如 `EM1.8_2.0_emax12_e3max16`），不要简写成 `EM18` 或省略 `emax12_e3max16`。
- 核力标签照 NO2B 数据目录名取（如 `EM1.8_2.0`、`DNNLO_go394`），路径/文件名保持一致即可。
- 每格 output basename 建议带全参数，如 `<Nuc>_<标签>_hw<hw>_Nmax<N>_beta<b>_auto12_no2b_<nkeep>states`。

### 核力怎么读（NO2B `no2`）

- NO2B 源：`/tns/mengziyan/tools/no2b/normal-order/data/<INT>/<Nuc>/no2b_<Nuc>_<标签>_hw<hw>_emax<emax>_e3max<e3max>.bin`（emax12 文件每个 ~1.9G）。
- BIGSTICK 用 `no2` 格式码读它（reader 见上一节）。**每格运行目录里把源二进制软链成 `<nuc>_hw<hw>_no2b.bin` 这类名字再喂给 BIGSTICK**，避免 `no2b_*.bin` 文件名被当成格式码（前三字符 `no2`）。
- `auto` 始终设 `12`（匹配 emax12 文件）；用多体 `Nmax` 做实际截断，**不要**拿更小的 `auto` 当截断手段。

### 输入 deck 模板（菜单 `n` = 普通 Lanczos）

```
n                       # 普通 Lanczos
<output basename>       # 输出文件名前缀
auto                    # 自动单粒子空间
12                      # max Nprincipal（emax12 文件固定 12）
<Z> <N>                 # valence 质子数 中子数（Li8: 3 5）
0                       # 2Jz（偶 A 取 0；要奇 K 跃迁/奇 A 见下）
+                       # parity
y                       # 是否做 W 截断：是
<Nmax>                  # W 截断 = Nmax
no2                     # 相互作用格式码
<linked .bin name>      # 软链后的 NO2B 文件名
<hw> <beta_cm>          # 振子频率 与 Lawson 质心项强度
end
ld                      # Lanczos（大 Nmax 换 td，见内存节）
<nkeep> <maxiter>       # 保留态数 与 最大迭代步数
```

- 振子频率经验值 `hbar_omega ≈ 45 A^(-1/3) − 25 A^(-2/3)`（A=8 → ~16.25 MeV），实际取最接近的现成 NO2B 文件 `hw`。
- 每格配一个 wrapper（`run_*.slurm`）：里面做 `ln -sf` 软链 + `cd output` + `"$EXE" < input.in`，把 `cd` 和重定向放进脚本，避免在命令行 "cd && 重定向" 触发权限询问。

### 参数扫描（hw × Nmax × beta_cm）

- **hw**（振子频率，MeV）：典型 `8, 12, 16, 20, 24, 28`。物理最优在 ~16 附近；高频（hw≥24）Lanczos 收敛更慢。
- **Nmax**（W 截断，多体激发上限）：`2, 4, 6, 8, 10(, 12)`。Nmax 越大越接近无截断极限，能量单调下降；计算量/内存随 Nmax 急增（见内存节）。
- **beta_cm**（Lawson 质心项强度）：典型 `0, 5, 10`。把质心激发（spurious CM）态推到高能、净化谱。
  - **beta=0**：低激发谱里混有 CM 污染态，不可信。
  - **beta=5 与 beta=10**：谱几乎一致（CM 态已推走），取它们看物理谱。基态能量对 beta 几乎不敏感（~0.1 MeV）。
- 逐格运行/收敛状态记进各核自己的 `<Nuc>_RUN_STATUS.md`（含 Current Slurm jobs、Aggregate、Convergence margins）。

### 收敛与 maxiter

- **收敛判据**：BIGSTICK `.out` 末尾 `(energy convergence X > criterion 0.00100)`，`X < 0.00100` 即收敛（X = 保留态中最差那个的能量变化）。撞到 maxiter 上限而 X 仍 ≥ 0.001 = **没收敛**，需加大 maxiter 重跑。
- **判 DONE 看非空**：`.res`/`.wfn` 在作业开始时是 0 字节占位，完成才写入；扫状态要 `[ -s file ]` 才算完成，否则是还在跑。
- **maxiter 经验**（下列数据 = Li8/emax12, EM1.8_2.0, hw16, beta5, 4 态实测；换核/emax/hw 会变）：Nmax 2/4/6/8/10 ≈ 110/150/190/220/250 步，约 +30/档且放缓；高频（hw28）需要 400~500 步。

### 内存标度与 `ld` vs `td`

下列绝对值是 Li8/emax12 实测，仅作量级参考；标度规律（单矢量 ∝ 基矢维数、随 Nmax 每档约 ×6–8）通用。

- 单 Lanczos 矢量 ∝ 基矢维数（real4 → 4 字节/分量），随 Nmax 每档约 ×6–8。Li8/emax12 hw16 beta5 实测：基矢 Nmax6/8/10/12 = 1.6e6/1.8e7/1.4e8/9.46e8；单矢量 ~6.6M/71M/575M/3.79G；jumps ~1G/11G/96G/**623G**。
- **`ld`（标准全重正交，菜单 `ld`）**：所有 Lanczos 矢量都留内存，峰值 ≈ 步数 × 单矢量 + jumps。Nmax10 峰值 ~330G（矢量 ~231G@400 步 + jumps ~96G），`c128m1024`（~1TB）够用。Nmax≤10 都用 `ld`。输入两行：`ld` / `<nkeep> <maxiter>`。
- **`td`（thick-restart，菜单 `td`，手册 Lanczos 节）**：维度关系 `Nkeep < Nthick < Niter`——每轮建 `Niter+1` 维 Krylov 子空间、截到 `Nthick` 维、再加回，重复到收敛或用满最大重启数。峰值矢量内存只跟 **Niter** 有关（保留 ~Niter+1 个矢量），与总迭代步数**解耦**。
  - 输入**三行**（漏第三行会读到 EOF 崩，已踩坑）：
    ```
    td
    <nkeep> <niter>      # 保留态数; 重启前迭代数（= 内存保留矢量数，手册建议尽量取大）
    <max_restarts>       # 最大重启次数（td 默认收敛会自动停，这是上限）
    ```
  - `td` 自动取 `Nthick = max(3·nkeep, nkeep+5)`；手册经验 `Nthick > ~3·nkeep`、`Niter` 尽量大。
  - 选 `Niter` 的内存账：`(Niter+1) × 单矢量 + jumps < 节点 RAM`。
- **Nmax12 实测（Li8/emax12, c128m1024）**：基矢 9.46e8、单矢量 3.79G、**jumps 623G**（`RAM for jumps in storage 637633 Mb`，`jumps built` 成功）。`ld` 跑 ~300 步光矢量 ~1.1TB → OOM，**必须 `td`**。`c128m512`（~503G）连 jumps（623G）都装不下 → 实测 OOM-killed；Nmax12 只能上 `c128m1024`。一个能跑的配置：`td` / `3 60` / `30`（窗口 ~61 矢量 ~231G + jumps 623G ≈ 854G，留 ~150G 余量）。
- 建 jumps 前打印的 `Estimated time per mat-vec multiply`（Nmax12 ~2.8e13 ns ≈ 7.7h）是**极度悲观的预估、与真实性能无关**：Nmax10 也印同款（~45 min/步）但实际只 ~36 秒/步。判可行性看 `RAM for jumps in storage` 真实数字 + 是否 `jumps built`，别看这个估计。

### Slurm

- 单节点、单 MPI rank；`OMP_NUM_THREADS` 32（小 Nmax）/ 64（大 Nmax）。
- 分区按空闲选 `c128m512` / `c128m1024`；Nmax10 用 `c128m1024` + ~512G。提交前 `sinfo` / `squeue` 看资源。
- **`c128m512` 内存上限 ~503G**（515000 MB）：请 `--mem=512G` 会一直 PENDING；Nmax≤8 实测只要 ~40–60G，256G 足够。
- wrapper 脚本里做 `ln -sf` 软链 + `cd output` + `bigstick < input.in`，避免 "cd && 重定向" 触发权限询问。

## 主要模块（密度部分，历史 TBTD 验证）

- 单体密度耦合 `coupled_densities` in `src/bdensities.f90`
- 两体密度 m-scheme 累加 `src/bdenslib4.f90`；J-coupling `src/bdenslib5.f90`（XX `couple_2bdensXX`、PN `couple_2bdensPN`、输出 `print_out_2bdens`）
- BIG↔KSHELL TBTD per-entry 已逐条验证（Mg24/USDB、Mg24/IMSRG、Na21/USDB，max ~1e-5）；关键修复 commit `5670a1c`（dens2bflag 下 full non-Hermitian same-species density path，配套 `bdenslib4/5`、`bjumplib_master/weld`、`bparallel_opbundles`）。判正确性只看与 KSHELL/独立参考的直接吻合，不看差异是否变小。

## TODO — Nmax12 待在外部平台跑（本机 c128m1024 已排满）

目标：把 Li8 `EM1.8_2.0` 最优参数（beta5）推到 **Nmax12**、**最低 3 态**、**thick-restart Lanczos**，hw = **12, 16, 20** 三个频率。本机 1 TB 节点排满，需换到另一台有 ~1 TB 内存节点的平台跑。下面是另一平台直接照搬即可的全部信息。

### 要带过去的文件

- BIGSTICK OpenMP 可执行：本仓库 `make -C make gfortran-openmp` → `bin/bigstick-openmp.x`。
- 三个 NO2B 二进制（~1.9G/个），从本机拷：
  - `no2b_Li8_EM1.8_2.0_hw12_emax12_e3max16.bin`
  - `no2b_Li8_EM1.8_2.0_hw16_emax12_e3max16.bin`
  - `no2b_Li8_EM1.8_2.0_hw20_emax12_e3max16.bin`
  - 本机源目录：`/tns/mengziyan/tools/no2b/normal-order/data/EM1.8_2.0/Li8/`
- 每个运行目录里把对应 bin 软链成 `li8_hw<hw>_no2b.bin` 再喂给 BIGSTICK（`no2b_*` 原名会被当格式码）。

### 输入 deck（hw16 例；换 hw 只改 ① 频率行 `<hw> 5` ② 软链文件名 ③ output basename 里的 hw）

```
n
Li8_EM1.8_2.0_emax12_e3max16_hw16_Nmax12_beta5_auto12_no2b_3states
auto
12
3 5
0
+
y
12
no2
li8_hw16_no2b.bin
16 5
end
td
3 60
30
```

`td 3 60 30` = thick-restart：保留 3 态、重启前 60 步（= 内存里保留的矢量窗口）、最多 30 次重启（`td` 默认收敛会自动停）。**三行缺一不可**，漏第三行会读 EOF 崩。细节见上面「内存标度与 `ld` vs `td`」。

### 资源（关键，决定能不能跑）

- 实测 Nmax12 内存（Li8/emax12）：基矢 9.46e8、单矢量 3.79G、**jumps 623G**。峰值 ≈ 窗口 ~61 矢量(~231G) + jumps(623G) ≈ **854G**。
- **必须 ~1 TB 内存节点**（≥ ~900G 可用）；~512G 节点连 jumps（623G）都装不下，会 OOM。
- 单节点、单 MPI rank、OpenMP；`OMP_NUM_THREADS=128`（占满核）。
- 建 jumps 前打印的 `Estimated time per mat-vec ~7.7h` 是假警报，忽略（见上节）；jumps 建好后每步 matvec 实际几十秒量级。

### 跑完验证

- 收敛：`.out` 末尾 `(energy convergence X > criterion 0.00100)`，`X < 0.00100` 即收敛。
- `.res` 非空才算跑完（开跑时是 0 字节占位）。基态 = state 1（Li8 基态 2+）。
- 三个 hw 的基态能量应比 Nmax10 更低（能量随 Nmax 单调下降）。

## 验证（编译/格式）

- `git diff --check`
- `make -C make gfortran-openmp`（本环境无 `ifort`，不用默认 `make openmp`）

进度详见 `TODO.md`。
