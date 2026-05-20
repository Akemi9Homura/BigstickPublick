# TODO

## Phase 1 — TBTD vs KSHELL 一致性排查（✅ 完成）

- [x] 在 Mg24/USDB 上对齐能谱 (BIG M=0, BIG M=1, KSHELL)
- [x] 独立写 OBTD 比较脚本（含 BIG isospin → KSHELL pn 转换）；确认 OBTD 全部对上
- [x] 独立写 TBTD 比较脚本，严格按 (ini,fin) 选 M、不跳 -999；定位 bug 类别
- [x] 用 OBTD ratio 确认 BIG M=1 第二个 2+ 与 KSHELL 反相（Lanczos 相位歧义），并在 TBTD 比较中扣除
- [x] 抽 BIGSTICK 手册和 KSHELL Brown PDF 的 OBTD/TBTD 公式（manual Eq.5.23、Brown Eq.22.70），证明定义等价
- [x] 修复 `src/bdenslib4.f90` PN backward 分支 OpenMP race（`reduction(+:dmatpn)`）
- [x] 修复 `src/bdenslib5.f90` `couple_2bdensXX` HC 分支（commit f3c063c：加 `(-1)^{Jab+Jcd-Jtot}`）
- [x] 修复 `print_out_2bdens` PN print indx 公式（abcouple↔cdcouple）
- [x] **关键修复 commit 5670a1c**：`dens2bflag` 模式下走 full non-Hermitian same-species density path（XX 分支不再用 triangular + HC 分支，改为 full block 全 m-pair 遍历）
- [x] Mg24/USDB 全核 BIG↔KSHELL per-entry 验证：16528 common，max |Δ|=8.1e-6
- [x] Mg24/IMSRG 5 态（含 0+,2+,2+,4+,3+）BIG↔KSHELL per-entry 验证：63042 common，max |Δ|=1.5e-5
- [x] 奇 A Na21/USDB BIG↔KSHELL per-entry 验证：53189 common，max |Δ|=2.6e-5

## Phase 2 — Ground truth 交叉验证（✅ 完成）

- [x] FCI 端逐条对算 KSHELL transit log 中 Mg24/USDB OBTD/TBTD（max|Δ|=1.7e-5）
- [x] 验证 ktransit(KSHELL 密度 × E2 op) ≡ FCI(RDM × E2 op)（4 核 O18/Ne20/Mg24/Si28 + E2_2b）
- [x] Path C 验证：`Observable::cal_wf_obs(fci, fci, bra, ket)` 直接夹 m-scheme 算符（O18/Ne20/Mg24 三核，三方一致）
- [x] KSHELL .snt → BIGSTICK .sps/.int 转换器 `util/snt_to_bigstick.py`（USDB 与自带文件对算 5 位精度内一致）
- [x] BIGSTICK wfn reader（shell-model-obs/libfci）+ ⟨wf|H|wf⟩ 验证：O18/Ne20/Mg24 三核 lanc_prec=real4 精度
- [x] BIGSTICK `.den2b` × E2_2b 算符 vs FCI 同算（4 核三方全对）
- [x] **反 Hermitian V ground truth 测试**：用 V_E2_antiHerm（rank-2 反 Hermitian rank-2 op，pp d5/2² J=0↔J=2 单 entry）撬开 Hermitian 算符看不见的 ρ_A 子空间。Pre-fix（f3c063c）反 Hermitian diff 量级 0.01~0.3；post-fix（5670a1c）diff ~1e-7 ≡ wf 直接夹。在 shell-model-obs 加 `antiHermitian` load flag + `bigstick_antiherm_test.cpp`

## Phase 3 — 后续

- [ ] 在更重核（如 Si28、S32）做 BIG↔KSHELL per-entry 验证作为更广回归
- [ ] 把目前的 5670a1c 修复正式审定，commit message 由 "Tentative" 改为 "Confirmed"，并加 regression test 到 `test/`
- [ ] 同步 AGENTS.md（5670a1c 修复 + 三核 per-entry 验证 + 反 Hermitian V 判决）
- [ ] 决定是否给 BIGSTICK 也加自动 m-up 功能（KSHELL 那种），避免用户手动 route 多个 M run
- [ ] 清理 /tmp/bigstick_prefix worktree（验证完成后可移除）
