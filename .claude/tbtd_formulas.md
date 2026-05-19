# TBTD 公式对照（BIGSTICK manual vs KSHELL Brown PDF）

## KSHELL（Brown PDF 第 22 章）

- 对子算符（含归一化因子 $N_{12} = 1/\sqrt{1+\delta_{ab}}$）：
  $$A^+(k_\alpha k_\beta J_o M_o) = -N_{12}[a^+_{k_\alpha} \otimes a^+_{k_\beta}]^{J_o}_{M_o}$$
  $$\tilde{A}(k_\gamma k_\delta J_o' M_o') = N_{12}[\tilde{a}_{k_\gamma} \otimes \tilde{a}_{k_\delta}]^{J_o'}_{M_o'}$$
- 时间反演 $\tilde{a}_{k,m} = (-1)^{j+m}a_{k,-m}$。
- TBTD：
  $$\mathrm{TBTD}(abkJ_oJ_o'\lambda) = \frac{\langle n\omega_a J_a \| [A^+(k_\alpha k_\beta J_o) \otimes \tilde{A}(k_\gamma k_\delta J_o')]^\lambda \| n\omega_b J_b\rangle}{\sqrt{2\lambda+1}}$$
- 求和的 pair 限制 $k_\alpha \le k_\beta$，$k_\gamma \le k_\delta$（避免重复计数）。

## BIGSTICK（manual A.4-A.6, 5.23）

- 对子算符**不含** $N_{12}$：
  $$\hat{A}^\dagger_{JM}(ab) = (\hat{a}^\dagger_a \otimes \hat{b}^\dagger_b)_{JM}$$
  $$\tilde{A}_{JM}(cd) = -(\tilde{c}_c \otimes \tilde{d}_d)_{JM}$$
- 时间反演同上 $\tilde{c}_{j,m} = (-1)^{j+m}c_{j,-m}$。
- TBTD（含外因子 $1/\sqrt{(1+\delta_{ab})(1+\delta_{cd})}$）：
  $$\rho^{fi}_J(ab,J_{ab}; cd,J_{cd}) = \frac{1}{[J]}\frac{\langle J_f \| [\hat{A}^\dagger_{J_{ab}}(ab) \otimes \tilde{A}_{J_{cd}}(cd)]_J \| J_i\rangle}{\sqrt{(1+\delta_{ab})(1+\delta_{cd})}}$$

## 关键观察：两个定义在公式层面**等价**

BIGSTICK $\hat{A}^\dagger$ = $-\sqrt{1+\delta_{ab}}\cdot A^+_{\rm KS}$；BIGSTICK $\tilde{A}$ = $\sqrt{1+\delta_{cd}}\cdot \tilde{A}_{\rm KS}$（负号、相位都一致）。外面除的 $\sqrt{(1+\delta_{ab})(1+\delta_{cd})}$ 抵消两边 $\sqrt{1+\delta}$ 后剩下：
$$\rho^{fi}_J(BS) = \frac{1}{[J]}\langle J_f \|[A^+_{\rm KS} \otimes \tilde{A}_{\rm KS}]_J\|J_i\rangle = \mathrm{TBTD}(\rm KS)$$

所以**两者公式相同**。差异必然来自代码实现。

## Wigner-Eckart 约定（两边都是 Edmonds）

Manual Eq. A.1：
$$\langle J_f M_f|\hat{O}_{KM}|J_i M_i\rangle = [J_f]^{-1}(J_iM_i,KM|J_fM_f)\,\langle J_f\|\hat{O}_K\|J_i\rangle$$

反过来 $\langle J_f\|\hat{O}_K\|J_i\rangle = [J_f]\,\frac{\langle J_fM_f|\hat{O}_{KM}|J_iM_i\rangle}{(J_iM_i,KM|J_fM_f)}$。

KSHELL `transit.F90` 中应用 reduce：
```fortran
c = dcg(jr, mr, jl, -ml, jj*2, mr-ml)
c = (-1d0)**((jr-mr)/2) / c
TBTD = c * x   ! x = get_cpld_tbtd
```
即除以 $\langle J_r M_r, J_l\,-M_l | jj, M_r-M_l\rangle$ 再乘 $(-1)^{J_r-M_r}$。对 $M_r=M_l=M$：reduce 因子 = $(-1)^{J_r-M}/\langle J_rM, J_l\,-M|K,0\rangle$。

BIGSTICK `src/bdenslib5.f90` (couple_2bdensXX 直接分支)：
```fortran
cleb4reduce = cleb(jtot*2, 0, JJi, Jz, JJf, Jz)  ! = <K 0, J_i M | J_f M>
... * sqrt(JJf+1.0) * (-1)**(Jtot + (JJf-JJi)/2) / cleb4reduce / sqrt(2*Jtot+1)
```
即除以 $\langle K\,0, J_i M | J_f M\rangle$，乘 $\sqrt{2J_f+1}\,(-1)^{K+J_f-J_i}/\sqrt{2K+1}$。

**两个 reduce 因子在数学上等价**（之前手工推导验证过，差异是 CG 行交换、3j 列对易、$(-1)^{2K}=1$ 等恒等式）。所以 reduce 不是 bug 源。

## 当前定位

- 公式两边一致；OBTD 完全对得上 → BIGSTICK 单体公式实现正确。
- TBTD bug **只能在**：
  1. m-scheme 累加阶段 (`src/bdenslib4.f90` 中 `dmatpp/dmatnn/dmatpn` 的填充)，或
  2. m→J 耦合阶段 (`src/bdenslib5.f90` `couple_2bdensXX/PN` 的 hfactor、CG、相位、HC 分支)。

具体怀疑点（写到 [[mg24_usdb_density_compare]] 中的 bug 类别）：
- XX direct vs HC 分支的第三 CG `cleb(2*Jab, 2*m, 2*Jcd, -2*m, ...)` vs `cleb(2*Jcd, -2*m, 2*Jab, 2*m, ...)` 差 $(-1)^{J_{ab}+J_{cd}-J_{\rm tot}}$。
- hfactor = $(1+\delta_{ab})(1+\delta_{cd})$ 是否与 BIGSTICK 自身的 TBTD 定义中 $1/\sqrt{(1+\delta_{ab})(1+\delta_{cd})}$ 一致——hfactor 应该是 $\sqrt{(1+\delta_{ab})(1+\delta_{cd})}$ 还是 $(1+\delta_{ab})(1+\delta_{cd})$？这取决于 `v(itbme)` 存的是什么归一化的 m-scheme ME。
- PN 块单 direct 分支：`Jab` 与 `Jcd` 在 `indx = pairblocksize*(pair1-pcref-1)+pair2-pcref+pcstart` 中的角色（`pair1` 是 destruction = cd, `pair2` 是 creation = ab）。对 Jab≠Jcd 的非对角对子块需要核对。
