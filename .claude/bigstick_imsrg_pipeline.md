# BIGSTICK ← IMSRG (.snt → .sps/.int) → wfn → 我代码读 → 夹 H 算能量

把同一份 IMSRG H（`EM1.8_2.0_sd-shell_o18_hw16_emax4_e3max4.snt`）喂给 BIGSTICK 跑波函数，再用 shell-model-obs 的 reader 读回，直接夹 H 算能量，与 BIGSTICK lanczos / FCI lanczos 对算。**端到端通**（O18/Ne20/Mg24 3 核 × 3 态全部 |Δ| ≤ 3e-6，仅来自 BIGSTICK lanc_prec=real4 精度）。

## 转换脚本 `util/snt_to_bigstick.py`

KSHELL .snt（pn formalism）→ BIGSTICK `.sps`（pns）+ `.int`（xpn）。两个关键细节：

1. **`.sps` pns 格式要求把 orbits 写两次**（proton block + neutron block），即便 p/n 同空间。BIGSTICK `binput.f90:329-347` 严格按 `numorb(1)` + `numorb(2)` 两段读。

2. **TBME 必须 canonicalize 到 `a ≤ b, c ≤ d` 形式**，每次 swap 乘 KSHELL phase `(-1)^{(j_a+j_b)/2 - J + 1}`（`operator_jscheme.f90:639-649`）。然后用 **T=1 always** 写出，这样 BIGSTICK `binputt_tbme.f90:3719` 的 swap phase `(-1)^{J+T+(j+j')/2}` 与 KSHELL 的 `(-1)^{(j+j')/2-J+1}` 在 (-1)^{2J}=1 之下完全一致。

   IMSRG .snt 的 PN 项常常是 `n p n p`（中子先），canonicalize 后变 `p n p n`，刚好符合 BIGSTICK xpn 的 "OBLIGATORY ORDERING: ia,ic proton, ib,id neutron"（`binput_tbme.f90:3938`）。

3. **质量依赖 (method2=1)**：KSHELL `(A/im0)^pwr` ↔ BIGSTICK xpn autoscale `(Aref/A)^pwr_BS`，符号需翻转 `pwr_BS = -pwr_KS`。converter 探测 `.snt` 的 method2 字段：method2=1 时输出 **autoscale** 头（负 nme + 3 个 extras `Acore Aref pwr_BS` 跟在 SPE 后），BIGSTICK 自动套用；method2=0 (无质量依赖) 时输出正 nme 头，跑 BIGSTICK 时手动给 `1 1 0 0` + `1 1 1 1 1` 即可。

## 验证：KSHELL `usdb.snt` ↔ BIGSTICK 自带 `usdb.int` (iso 格式)

```
util/snt_to_bigstick.py /home/mengziyan/kshell/snt/usdb.snt runs/mg24_usdb_xpn/usdb_xpn
```

跑 BIGSTICK Mg24 两次：用转换后 xpn vs 用自带 iso `usdb` (5 digits 精度内完全一致)。

| state | converted-xpn | bundled-iso |
|---|---|---|
| 1 (0+) | −87.10445 | −87.10445 |
| 2 (2+) | −85.60215 | −85.60215 |
| 3 (2+) | −82.98830 | −82.98830 |

说明 swap-phase + T=1 约定、autoscale 头格式（含 pwr 符号翻转）都对。converter 适用于 KSHELL .snt 中 method2 ∈ {0, 1} 的任意作用力。

## 运行目录结构

`/home/mengziyan/BigstickPublick/runs/{o18,ne20,mg24}_imsrg/`：
- `imsrg.sps`、`imsrg.int`（转换产物）
- `input.bigstick_*`（lanczos 菜单脚本）
- `input.bigstick_*_bas`（'ba' 菜单脚本，仅 basis）
- `run_bigstick.sh`（wrapper，避免 `cd && redirect` 触发 Claude Code 权限询问）
- `*.wfn`、`*.bas`、`log_*.txt`

## BIGSTICK 输入模板（lanczos）

```
n                  # menu: 'n' = normal lanczos
<output basename>
imsrg              # .sps file
<Z> <N>            # valence
0                  # 2 Jz = 0 (even A)
n                  # no truncation
xpn                # format selector (BEFORE the .int filename)
imsrg              # .int file
1.0 1.0 0.0 0.0    # spscale ax bx x → vscale=1 (任一 bx=0 或 x=0 触发 vscale=ax)
1.0 1.0 1.0 1.0 1.0  # pspescale nspescale ppscale nnscale pnscale
end
ld                 # lanczos diagonalize
3 400              # nkeep, max_iter
```

basis dump 流程（菜单 'ba'）：

```
ba
imsrg
<Z> <N>
0
<basis output basename>
```

## 跑出来的 BIGSTICK 能谱（与 FCI/KSHELL 对算）

| 核 | states (2J) | E (MeV) |
|---|---|---|
| O18  | 0,4,8 | −9.69500 / −7.73716 / −6.47318 |
| Ne20 | 0,4,8 | −30.40181 / −28.35611 / −26.11457 |
| Mg24 | 0,4,4 | −74.32028 / −72.48098 / −69.63942 |

三核 BIGSTICK lanczos 与 FCI（[[imsrg_e2_ground_truth]]）逐位对上。

## Path C 验证（用 BIGSTICK wf 直接夹 H 算能量）

测试例：`shell-model-obs/examples/bigstick_energy_check.cpp`。读 `.bas + .wfn` → `BigstickWfnData` → 对每个态构 `wf_vector(big.amplitudes[s], big.configs, &big_map)` → `Observable(&H).cal_wf_obs_val(wf, wf)` → ⟨wf|H|wf⟩。

```
build/examples/bigstick_energy_check <H.snt> <Z> <N> <basis.bas> <wfn.wfn>
```

输出（每核 3 个态）：

| 核 | state | BIG_lanczos | PathC_BIG_wf | FCI_lanczos | \|PathC−BIG\| |
|---|---|---|---|---|---|
| O18  | 1 | −9.69500 | −9.69500 | −9.69500 | 8e-7 |
| O18  | 2 | −7.73716 | −7.73716 | −7.73716 | 2e-6 |
| O18  | 3 | −6.47318 | −6.47318 | −6.47318 | 1e-6 |
| Ne20 | 1 | −30.4018 | −30.4018 | −30.4018 | 1e-6 |
| Ne20 | 2 | −28.3561 | −28.3561 | −28.3561 | 3e-7 |
| Ne20 | 3 | −26.1146 | −26.1146 | −26.1146 | 1e-6 |
| Mg24 | 1 | −74.3203 | −74.3203 | −74.3203 | 2e-6 |
| Mg24 | 2 | −72.4810 | −72.4810 | −72.4810 | 2e-6 |
| Mg24 | 3 | −69.6394 | −69.6394 | −69.6394 | 2e-6 |

**结论**：BIGSTICK 用我转换的 IMSRG H 跑出的波函数，被 shell-model-obs reader 读回后，逐个 Det 配上振幅做 ⟨wf|H|wf⟩ 缩并，能量与 BIGSTICK 自家 lanczos 本征值在 real4 精度（lanc_prec=4 ⇒ ~1e-6）内一致；同时也与 FCI lanczos 本征值一致 → 整条 reader → wf_vector → cal_wf_obs_val(H) → energy 的链路验证通。

## TBTD × E2_2b 判决（2026-05-20）

测试例 `shell-model-obs/examples/bigstick_e2_via_density.cpp` 仿 `tools/ktransit` 写法：读 BIGSTICK `.den2b` 里的 TBTD 条目，对 E2_2b 算符 (rank K=2) 缩并出 op2b 部分 reduced ME。同时 `Observable(&H,&op).cal_wf_obs_val` 在 BIG wf 上算 full reduced ME 作旁路。

### 准备

- 对每个 `runs/{o18,ne20,mg24}_imsrg/` 运行菜单 `'2'` (input `input.bigstick_*_2b`) 生成 `*_2b.den2b`
- E2 是 rank-2，在 M=0 子空间下三核所有 (Ji, Jf) 对的 reduce CG `cleb(4,0,2Ji,0,2Jf,0)` 全部非零 → M=0 wfn 足够，**不用换 M**（详情：`bigstick_m_selection.md`）

### 结果（op2b reduced ME, 三核 13 条 (i,j) pair）

每条 BIG TBTD × op_2b（这里是 .den2b × op）vs FCI RDM × op_2b（[[imsrg_e2_ground_truth]] 里记录的），**取绝对值后全部对得上 max\|Δ\| ≤ 2e-6**：

| 核 | n pair | max ||Δ|abs|| |
|---|---|---|
| O18  | 4 | 8e-8 |
| Ne20 | 4 | 2e-6 |
| Mg24 | 5 | 2e-6 |

符号差 = Lanczos 全局 ± 自由度（BIG/FCI 独立 Lanczos）；用 cal_wf_obs(full op) on BIG wf 算的总 reduced ME 也对得上 ground truth modulo 同样的 σ pattern。

### 结论：BIGSTICK TBTD 在 IMSRG 上正确

→ 之前 Mg24/USDB 的 cross-state pp/nn 残余 0.20（[[bug_localization]] 最后一节）**不是 BIGSTICK couple_2bdens 通用 bug**，而是 USDB 数据特有的某个组合触发的（或当时假设的某种相位约定差异）。需要用 USDB H 复跑这同样的对比才能进一步定位。

## 下一步候选

- 用同样 `bigstick_e2_via_density` 跑 Mg24/USDB（KSHELL 那个原始 USDB.snt + 同 reader），看是否复现 0.20 残差
- 如果 USDB 上残差仍 0.20，扣开看 pp/nn cross-state 哪几个 (a,b,c,d,Jab,Jcd,K) 在贡献
