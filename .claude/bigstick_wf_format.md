# BIGSTICK 波函数解码规格（`.wfn` + `.bas`）

不用改 BIGSTICK。生成 `.bas`（basis listing）一次，加上原本就有的 `.wfn`（amplitude），就足够在 shell-model-obs 里端到端解出 (Slater determinant bitset, amplitudes)。

## 生成 `.bas` 文件（一次）

BIGSTICK 自带菜单 `'b'`（binary）/ `'ba'`（ASCII）：在 sps/Z/N/Jz 设置好后写出 sector + SD 列表。Mg24/USDB M=0 的输入：

```
ba
sd
4 4
0
mg24_basis_m0
```

M=1 同样：

```
ba
sd
4 4
2
mg24_basis_m1
```

输出 `mg24_basis_m0.bas` ≈ 97 KB（Mg24/USDB），均为 ASCII，1186 行左右。**写 `.bas` 不需要 interaction，不跑 Lanczos**，只构造 basis。OpenMP/串行都可以，但 `'ba'` 在 MPI 下会拒绝。

## `.bas` (ASCII) 格式

源码：`src/boutputlib.f90:basis_out4postprocessing` 在 `baseASCII=.true.` 分支（line 1043-1247 附近）。

### Header（line 1–39 in Mg24 example）

```
wfn_magic   = 31415926        int
wfn_version = 1                int
np(1) np(2)                    int int        ! valence Z, N
isoflag                        T/F
numorb(1) numorb(2)            int int
for it=1,2; for n=1..numorb(it):
    nr  j  l  par  w           5 × int        ! orbqn (j is doubled)
nsps(1) nsps(2)                int int
for it=1,2; for n=1..nsps(it):
    nr  j  m  l  w  par  orb  group   8 × int  ! spsqn (j,m doubled)
jz  cparity  maxWtot           int char int   ! e.g. "0 + 0"
allsameparity allsameW spinless   3 × logical
dimbasis                       int            ! e.g. 28503
```

### Haiku 半空间单粒子表（line 40–64）

```
nhsps array                    5 × int  ! 实为 nhsps(-2:2)；e.g. "6 6 0 6 6"
                               ! nhsps(0) 没用 =0；其余 = ith=-2,-1,1,2 的 half-space sps 数
for each ith ∈ {-1, +1, -2, +2}:  ! 顺序：proton-left, proton-right, neutron-left, neutron-right
    for i = 1 .. nhsps(ith):
        nstate  nr  l  j  m  w  par     7 × int   ! hspsqn (j,m doubled)
```

`hspsqn` 是把 spsqn 按 m 符号拆成左半 (m<0) / 右半 (m>0) 两块，proton 部分占 nstate=1..(nhsps(-1)+nhsps(1))，neutron 接在后面。

### Proton sectors（line 65 开始）

```
nsectors(1) nxSD(1)            int int   ! e.g. "13  495" → 13 proton sectors, 495 proton SDs

for ps = 1 .. nsectors(1):
    ps  jzX  parX  Wx          4 × int            ! sector quantum numbers (jzX doubled)
    xsdstart  xsdend  nxsd     3 × int            ! proton SD index range in this sector
    ncsectors                  int                ! 此 ps 对应的 neutron 共轭 sector 数
    csector(1..ncsectors)      ncsectors × int    ! 共轭 neutron sector 列表
    basisstart basisend        2 × int            ! 当前 'b'/'ba' 调用点尚未填入，常见 0 0
    for SD in [xsdstart, xsdend]:
        ip  pstart(ip)  occupied_sps(1..np(1))    1+1+np(1) × int
                                                  ! occupied_sps 是 hspsqn 表里的 nstate
```

### Neutron sectors（紧跟 proton sectors 之后）

```
nsectors(2) nxSD(2)            int int   ! e.g. "13  495"

for ns = 1 .. nsectors(2):
    ns  jzX  parX  Wx
    xsdstart  xsdend  nxsd
    ncsectors
    csector(1..ncsectors)               ! 共轭 proton sector 列表
    basisstart basisend
    for SD in [xsdstart, xsdend]:
        in  nstart(in)  occupied_sps(1..np(2))    ! occupied_sps 是 hspsqn 表里的 nstate
                                                   ! 但注意 nstate ∈ [nhsps(-1)+nhsps(1)+1, totalsps]
```

## Basis index 公式（最关键）

源码：`boutputlib.f90:1038` 等价。

```
ibasis = nstart(in) + pstart(ip)
```

其中：

- 对每个 proton sector ps，它的 ncsectors 列表给出"哪些 neutron sectors 与之共轭"
- 对每个共轭对 (ps, ns)，遍历 ip ∈ proton-SD-list(ps)、in ∈ neutron-SD-list(ns)
- `pstart(ip)` 与 `nstart(in)` 在 `.bas` 文件里直接给出（不需要计算）
- 这样 `ibasis` 跑遍 [1, dimbasis]

Mg24 sector 1 验证：

```
ps=1, jzX=-12, parX=1; 3 proton SDs ip=1,2,3 with pstart = 0, 3, 6
   csector = {1}; neutron sector 1 has 3 SDs in=1,2,3 with nstart = 1, 2, 3
ibasis = nstart(in) + pstart(ip):
   (ip=1, in=1) → 0+1 = 1
   (ip=1, in=2) → 0+2 = 2
   (ip=1, in=3) → 0+3 = 3
   (ip=2, in=1) → 3+1 = 4
   (ip=2, in=2) → 3+2 = 5
   ...
   (ip=3, in=3) → 6+3 = 9
```

ibasis 顺序：**外层 ip，内层 in**（与 `boutputlib.f90:1037` 的双 do 循环 nladd/nradd × pladd/pradd 一致）。

## `.wfn` (binary stream)

源码：`src/bwfnlib.f90:write_wfn_header` + `wfn_writeeigenvec`。`open(... access='stream', form='unformatted')`，纯字节流。

### Header

```
wfn_magic = 31415926                   int4
wfn_version = 1                         int4
0                                       int4    ! placeholder
offset_to_vec_data                      int4    ! 第二轮 inquire 后回填
0  0  0  0                              4 × int4
np(1) np(2)                             2 × int4
isoflag                                 int4 logical (0/1)
numorb(1) numorb(2)                     2 × int4
for it=1,2; for n=1..numorb(it):
    nr  j  l  par  w                    5 × int4
nsps(1) nsps(2)                         2 × int4
for it=1,2; for n=1..nsps(it):
    nr  j  m  l  w  par  orb  group     8 × int4
jz  ICHAR(cparity)  maxWtot             3 × int4
allsameparity allsameW spinless         3 × int4 logical
dimbasis                                basis_prec   ! int4 或 int8 by compile
```

### Vector data

```
nkeep                                   int4
for i = 1 .. nkeep:
    wfn_vecmagic = 27182818              int4
    i                                    int4
    e  xj  xt2                           3 × real4
    v(1 .. dimbasis)                     dimbasis × lanc_prec
```

`lanc_prec`、`basis_prec` 由编译期决定。Mg24/USDB 经验：`lanc_prec = 4`、`basis_prec = 4`（从 `m0_3.wfn` 文件大小 343080 B、dimbasis=28503、3 个 nkeep 反推）。

## 在 shell-model-obs 里落地

`src/libbigstick/`（与 `libkshell/` 并列）：

| 文件 | 内容 |
|---|---|
| `BigstickBas.{hpp,cpp}` | 解析 `.bas` ASCII：sps 表、hsps 表、sector 表、SD 表 + `pstart/nstart` |
| `BigstickWfn.{hpp,cpp}` | 解析 `.wfn` header + vector blocks |
| `BigstickBasis.{hpp,cpp}` | 把 `.bas` + `.wfn` 组合：对每个 ibasis 给 (proton occupation set, neutron occupation set) → bitset → amplitude |
| `tests/test_bigstick_mg24.cpp` | 验证 ‖ψ‖²=1、维度=28503、能量与 KSHELL 对应态吻合 |

## sps 编号映射注意

`.bas` 给的 occupied_sps 是 `hspsqn` 表的 nstate（1..24，含 proton-left/right + neutron-left/right）。要映射到 shell-model-obs 的 `bitset<NMO>`，需查 hspsqn 行的 `(nr, l, 2j, 2m)` 量子数，再对应到 shell-model-obs 自己的 (orbit, m) bit 编号。Mg24 sps 表（`.bas` 行 41-64）：

```
proton-left (m<0):   1..6   = (p 0d3/2 m=-1/2,1/2), (p 0d5/2 m=-1/2,1/2), (p 1s1/2 m=-1/2,1/2)  按 m 升序
                     ↑ 注意：hspsqn 的 left 半也含 m=-1/2 与 m=+1/2 (m<0 strict 还是 ≤0 看 BIGSTICK 约定，需要复核)
```

精确映射需要 verify 一次，建议读 reader 时 print 出来对照。

## 后续

完成 reader 后，先在 shell-model-obs 里复现 ‖ψ‖²=1、⟨ψ_i|ψ_j⟩=δ_ij。通过后再加 ⟨f|O|i⟩ 直接夹算符，与 BIGSTICK/KSHELL 输出的 TBTD 缩并对算符的结果对算（ground truth）。
