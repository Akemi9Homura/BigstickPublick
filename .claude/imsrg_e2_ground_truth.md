# IMSRG-H + E2_2b 在 O18/Ne20/Mg24 上的 ground truth

将来要把 BIGSTICK 用同一组 IMSRG H 跑出的波函数喂回 shell-model-obs，做"BIGSTICK wf 直接夹 op" vs 下表算 ground truth 比对。

## 来源

- IMSRG H：`shell-model-obs/temp/EM1.8_2.0_sd-shell_o18_hw16_emax4_e3max4.snt`
- E2_2b 算符：`shell-model-obs/temp/EM1.8_2.0_sd-shell_o18_hw16_emax4_e3max4E2_2b.snt`
- 验证由三条独立路径完成（A=ktransit, B=FCI RDM×op, C=FCI wf-sandwich），结果一致。
- 见 [[ktransit_fci_verification]]。

数值取 Path B (RDM × op) 与 Path C (wf-sandwich) 的共同值，符号约定固定为 FCI Lanczos 给的符号（state 1 = 0+_1，state 2 = 2+_1，state 3 = 2+_2 或 4+_1，按能量顺序）。

## 能谱 (MeV)

| state | 2J | O18 | Ne20 | Mg24 |
|---|---|---|---|---|
| 1 | 0 | −9.69500 | −30.40181 | −74.32028 |
| 2 | 4 | −7.73716 | −28.35611 | −72.48098 |
| 3 | 8 (O18/Ne20) / 4 (Mg24) | −6.47318 | −26.11457 | −69.63942 |

(KSHELL Lanczos 收敛值；FCI 各自打印 6 位精度内一致)

## E2 reduced ME ⟨bra‖E2_2b‖ket⟩（取 Path B/C 共同值；单位由 op 文件决定）

### O18 (states: 0+, 2+, 4+)

| (bra, ket) | sum (= op1b + op2b) | 说明 |
|---|---|---|
| (1, 2) | −1.77212 | 0+ → 2+ |
| (1, 3) | — | 不耦合 (rank 2 不连 0↔4) |
| (2, 2) | −0.125296 | 2+ diag |
| (2, 3) | −2.15739 | 2+ → 4+ |
| (3, 3) | −2.19679 | 4+ diag |

### Ne20 (states: 0+, 2+, 4+)

| (bra, ket) | sum |
|---|---|
| (1, 2) | +9.88618 |
| (2, 2) | −11.38080 |
| (2, 3) | +14.33480 |
| (3, 3) | −14.92340 |

### Mg24 (states: 0+, 2+_1, 2+_2)

| (bra, ket) | sum |
|---|---|
| (1, 2) | +13.10590 |
| (1, 3) | −3.25477 |
| (2, 2) | −14.77890 |
| (2, 3) | −6.01743 |
| (3, 3) | +14.98590 |

## 对称性

- `(i, j) sum == (j, i) sum` (在 Path B/C 输出中)
- 沿对角的诊断值 (i, i) 与 BIG 自家 OBTD/TBTD 缩并的对算才是 BIGSTICK bug 的最终判决

## 使用方法

将来若用 BIGSTICK wf reader + `cal_wf_obs_val(wf_bra, wf_ket)`：
- H 标量直接 sandwich → 与上表能量对（无 Wigner-Eckart）
- E2 非标量 → sandwich 后按 `(-1)^{(J2bra-Mbra)/2} × wigner_3j(...)` 反推 reduced ME，与上表 sum 对。注意 BIGSTICK 的 Lanczos 全局符号可能与 FCI 不同，cross-state 项允许 ±1 全局相位偏差。
