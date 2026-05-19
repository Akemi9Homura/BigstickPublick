# TODO

## Phase 1 — TBTD vs KSHELL 一致性排查

- [x] 在 Mg24/USDB 上对齐能谱 (BIG M=0, BIG M=1, KSHELL)
- [x] 独立写 OBTD 比较脚本（含 BIG isospin → KSHELL pn 转换）；确认 OBTD 全部对上
- [x] 独立写 TBTD 比较脚本，严格按 (ini,fin) 选 M、不跳 -999；定位 bug 类别
- [x] 用 OBTD ratio 确认 BIG M=1 第二个 2+ 与 KSHELL 反相（Lanczos 相位歧义），并在 TBTD 比较中扣除
- [x] 抽 BIGSTICK 手册和 KSHELL Brown PDF 的 OBTD/TBTD 公式（manual Eq.5.23、Brown Eq.22.70），证明定义等价
- [x] 对比 KSHELL `get_cpld_tbtd:pair_sum` 与 BIGSTICK `couple_2bdensXX` 的累加因子，发现 HC 分支第三 CG 多一 `(-1)^{Jab+Jcd-Jtot}` 相位
- [x] 用 18O 解析（纯 (d5/2)² J=2 配对模型）判决 KSHELL 对 BIGSTICK 错
- [x] 修复 `src/bdenslib5.f90` `couple_2bdensXX` HC 分支：加 `(-1)^{Jab+Jcd-Jtot}` 因子
- [x] 验证修复：18O state 2+ self-diag K=0..4 全 0.6056080 ✓；Mg24 SAME pp/nn 全 K 全 Jab=Jcd 与 Jab≠Jcd 都达到 max < 1e-5 ✓
- [x] 追 `print_out_2bdens` PN print indx 公式：发现 (creation, destruction) 与 compute 的 (destruction, creation) 反，已交换 abcouple↔cdcouple 修复
- [x] 验证 PN 修复：Mg24 PN 块所有类别 max < 6e-6 ✓
- [ ] 追 XX cross-state pp/nn 残差：当前 max 0.20，呈 Hermitian-反对称模式（BIG 与 KSH 对差一致但对和不一致），需分析 cross-state 下 HC 公式是否需调整
  - [ ] **用自写的直接夹算符代码做 ground truth**：直接读波函数文件，按 ⟨f|O^K|i⟩ 直接夹两体算符算 reduced ME；判决某 Mg24 cross 入口 BIG/KSHELL 谁对（密度矩阵 × 算符 缩并 必须等于波函数直接夹算符）
  - [ ] 比较 BIG `direct branch + HC branch` 求和 与 KSHELL `pair_sum (a^+_p p_p product over m's)` 的 m-scheme 累加结构
- [ ] 最终验证：Mg24 max diff 目标 < 1e-5 全类别

## Phase 2 — 修复后回归测试

- [ ] 在另一对核（如 Ne20, O18 或 Si28）上重做一遍 BIG↔KSHELL 比较
- [ ] 加单元/回归测试到 `test/`，覆盖 XX 和 PN 在 Jab=Jcd, Jab≠Jcd, K 偶/奇 的各分支

## Phase 3 — 文档与 Codex 同步

- [ ] AGENTS.md 同步当前结论（见 `/同步codex` 技能）
- [ ] BIGSTICK 手册若有相关章节需要更新，记下要改的部分
