# KSHELL 自动 m-up 机制

KSHELL transit 在算 OBTD/TBTD 前，会检查每个 (left=final, right=initial) 配对的 reduce CG 是否为零：

```fortran
c = dcg(jr, mr, jl, -ml, jj*2, mr-ml)
if (abs(c) < thd_mup) then
   if (jr >= mtotr+2) is_mup(i,j) = .true.
   cycle
end if
```

被标记的 (i,j) 走 `calc_mup` 分支，**升的是 ket (right)**：

```fortran
! src/transit.F90:1042
! |ini J, M+1> = 1/N * J+ |ini J, M>
mupr = mtotr + 2   ! ket 的 doubled M 升 2，即 Jz +1
```

bra 维持原 `mtotl`，重新算 reduce 用 `dcg(jr, mupr, op%irank*2, mtotl - mupr, jl, mtotl)`。所以 KSHELL 输出的"真值"对应了 **不同 M-projection** 但**相同 reduced ME**。

这导致 BIGSTICK 必须用合适的 M（让 reduce CG ≠ 0）才能跟 KSHELL 比，否则 BIGSTICK 那个 (i,f) 的整段不可用（输出 -999 哨兵）。详见 [[bigstick_m_selection]]。

提示信息：`*** called mup ***` 出现在 transit log 中一次（对每次 transit 调用最多一次，因为只用了一个 `mupr = mtotr + 2` 的 ket）。
