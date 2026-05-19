# BIGSTICK -999 哨兵与 M 选择规则

BIGSTICK 把 m-scheme 累加完的密度乘以 reduce 因子 `sqrt(JJf+1.0) * (-1)^{Jtot+(JJf-JJi)/2} / cleb(2*Jtot, 0, JJi, Jz, JJf, Jz)`（在 `src/bdenslib5.f90` 三个 `couple_2bdens*` 子程序中）。当 `cleb4reduce = cleb(2*Jtot, 0, JJi, Jz, JJf, Jz)` 接近 0 时，密度矩阵元被设为 **-999.** 哨兵，并跳过本条：

```fortran
cleb4reduce = cleb(jtot*2, 0, JJi, Jz, JJf, Jz)
if (abs(cleb4reduce) < 0.00001) then
   x2bden(indx)%v(Jab, Jcd, Jtot) = -999.
   cycle
end if
```

邮件（`打印邮件.pdf`，Calvin W. Johnson 2024-06-30）：**遇到这种情况要换 M (2M=2) 重跑该跃迁，不要把 -999 项跳过当成可用结果。**

## Mg24 / USDB 的 M 选取实例

- (0+,0+)、(0+,2+)、(2+,0+) 在 M=0 下 reduce CG 全部非零 → 用 M=0
- (2+,2+)、(2+,2+'）等 2+↔2+ 跃迁：M=0 下 K=1, 3 (奇 K) 出 −999（因 CG(1,0,2,0|2,0)=0、CG(3,0,2,0|2,0)=0）→ 整段不可用 → 换 M=1 (2Jz=2)
- 在 M=1 下，对每个 K 检查 reduce CG `cleb(2*K, 0, 4, 2, 4, 2)` 均非零，整段可用

## 验证段是否可用

```bash
# 列出每个 (ini, fin) 段的 -999 个数
awk '/^ !# Ini state/{wi=1;next}
     /^ !# Fin state/{wf=1;next}
     wi && NF>=3 {ini=$1; wi=0; next}
     wf && NF>=3 {fin=$1; wf=0; next}
     NF>=8 { for(i=9;i<=NF;i++) if($i+0 < -900) s[ini"->"fin]++ }
     END { for(k in s) print k, s[k] }' mg24_usdb_m0_3_2full.den2b
```

正确选 M 后，目标段计数应为 0。

## 切勿

- 合并 M=0 和 M=1 的输出"补全"−999；这违反原则（同一 ME 在不同 M 下经过不同的 reduce 因子，约定一致才能合并，且实测合并后 K 偶/奇分布跟单独跑的差不一致）。
- 把 K 偶部分留 M=0、K 奇部分换 M=1 用于同一跃迁；同一 (ini,fin) 的所有 K 必须来自同一个 M，否则隐含约定混乱。
