# KSHELL 波函数格式与解码规格（给 shell-model-obs 的 reader）

要从 KSHELL 取出 (determinant, amplitude) 列表做 ground truth，必须**同时**读 `.snt + .ptn + .wav` 三个文件，复刻 KSHELL 的 partition 枚举顺序。不能只读 `.wav`：它只是一段扁平的 amplitude 字节数组，没有任何 determinant 标签。

## 文件 1: `.snt` (single-particle space)

源码：`kshell/src/model_space.f90:read_sps`。

```
# 跳过开头所有以 ! 或 # 起头的注释
n_jorb_p  n_jorb_n  n_core_p  n_core_n          ! 第一行数据
{ index, n, l, j, tz }                          ! proton orbits (tz=-1) ×  n_jorb_p
{ index, n, l, j, tz }                          ! neutron orbits (tz=+1) × n_jorb_n
```

约定：j 是 doubled (即 j=3 表示 3/2)，tz=-1 是 proton、tz=+1 是 neutron。proton orbits 必须在前。

## 文件 2: `.ptn` (partition)

源码：`kshell/src/partition.F90:init_partition` line 90-128。

```
n_proton  n_neutron  parity
n_id_p   n_id_n
{ id, occ_orb1, occ_orb2, ..., occ_orb_n_jorb_p }   ×  n_id_p
{ id, occ_orb1, ..., occ_orb_n_jorb_n }             ×  n_id_n
n_pn_combination
{ idp, idn }                                         ×  n_pn_combination
```

每个 partition id = 一种轨道占据分配（不区分 Jz）。combined list = 物理空间允许的 (proton-partition, neutron-partition) 组合。

## 文件 3: `.wav` (波函数，二进制 stream)

源码：`kshell/src/wavefunction.F90:load_wf` (non-MPI 分支 line 771-803)、`bp_io.F90:bp_save_wf`。

Fortran `open(..., form='unformatted', access='stream')` — 无 record marker，纯字节流。

```
header:
   integer(4)  n_eigen                    ! 本征态数
   integer(4)  mtotal                     ! = 2 Jz
   real(8)     eval(n_eigen)              ! 能量
   integer(4)  jj(n_eigen)                ! 2 J
body:
   real(kwf)   psi_1[ndim]                ! ndim = ptn%ndim
   real(kwf)   psi_2[ndim]
   ...
   real(kwf)   psi_n_eigen[ndim]
```

`kwf` 编译期常数（4 或 8）。从 `(filesize - header) / n_eigen / ndim` 反推；Mg24/USDB Mg24_usdb_m0p.wav (684116 B, n_eigen=3, ndim=28503) → kwf=4 (real4)。

## KSHELL m-scheme bit 编号（必须 1:1 复刻）

源码：`partition.F90:init_mbit_orb` line 547-613。

- **proton 与 neutron 各用一个独立整数 mask `mbp`, `mbn`**，不共享位
- 每个 ipn 内 **bit 0 跳过不用**（line 589 `i = 1`）
- proton orbit k 占 bit `[1 + sum_{k'<k}(j_{k'}+1)] .. [+j_k+1]`
- 每个 orbit 内部，bit 由低到高对应 m = -j, -j+2, ..., +j

orbit-level 子 mask（`mbit_orb(k, n, ipn)`）= 所有 (j+1)-bit 的 popcount=n 整数，**按数值 0..2^(j+1)-1 升序遍历筛 popcount=n 后入数组**。

## Partition × Mz 内 mb 列表的枚举顺序

源码：`partition.F90:set_ptn_mbit_arr` line 400-416。

对 (proton 或 neutron 一边的) partition id 给定 nocc 数组：

```
for i = 0 .. (∏_k C(j_k+1, n_k)) - 1:
   j = i;  mb = 0;  mm = 0
   for k = 1 .. n_jorb(ipn):
      n_choices = mbit_orb(k, n_k, ipn)%n   = C(j_k+1, n_k)
      pick      = j mod n_choices
      j         = j / n_choices
      mb       |= mbit_orb(k, n_k, ipn)%mbit[pick]   ! 已经 ishft 过到 orbit's bit start
      mm       += m_of(pick)
   bucket  mb  →  mbs[mm]
```

mixed-radix，**外层 = orbit n_jorb，内层 = orbit 1**。每个 (id, mm) bucket 得 `mbit[1..n]` 数组。

## Sorted partition (pidpnM_pid_srt) 顺序

源码：`partition.F90` line 208-226。

```
for k = 1 .. n_pn_combination:                    ! combined list 的输入顺序
   i = pidpn_pid(1, k);  j = pidpn_pid(2, k)
   for mp = id_p(i)%min_m  step 2  to  id_p(i)%max_m:
       mn = mtotal - mp
       if mn ∉ [id_n(j)%min_m, id_n(j)%max_m]: skip
       emit (i, j, mp)
```

`min_m / max_m` 初值 = `±max_m_nocc(nocc, ipn)` (line 143)，然后被同 combined-list 的 partner partition 钳到 `mtotal − partner_max_m` (line 161-173)。这步要按精确顺序复刻。

## `.wav` body 内 amplitude 顺序（最终关键）

```
idx = 0
for n = 1 .. n_pidpnM (sorted):
    (i, j, mp) = pidpnM_pid_srt(:, n)
    mn = mtotal - mp
    for in = 1 .. id_n(j)%mz(mn)%n:          ! 外层 = neutron
       mbn = id_n(j)%mz(mn)%mbit[in]
       for ip = 1 .. id_p(i)%mz(mp)%n:       ! 内层 = proton
          mbp = id_p(i)%mz(mp)%mbit[ip]
          idx += 1
          → psi_state[idx] 对应 determinant (mbp, mbn)
```

`ndim = Σ id_n(j)%mz(mn)%n × id_p(i)%mz(mp)%n` over sorted partitions = `.wav` 文件 body 中每个本征矢的元素数。

## 在 shell-model-obs 里落地建议

子模块：`src/libkshell/`

| 文件 | 内容 |
|---|---|
| `KshellSnt.{hpp,cpp}` | 解析 `.snt` → `OrbitTable` + m-scheme bit layout |
| `KshellPtn.{hpp,cpp}` | 解析 `.ptn` → proton / neutron partition lists + combined list |
| `KshellMbitEnumerator.{hpp,cpp}` | `init_mbit_orb` + `set_ptn_mbit_arr` + min/max_m 钳制 + sorted partition |
| `KshellWav.{hpp,cpp}` | 读 `.wav` header + body，按 sorted 顺序绑定 amplitude ↔ (mbp, mbn) |
| `tests/test_kshell_mg24.cpp` | Mg24/USDB 验证 ‖ψ‖²=1、⟨ψ_i|ψ_j⟩=δ_ij、ndim 一致 |

预估代码量 ~600-800 行 C++，partition 枚举核心 ~200 行。

## 主要坑

- **bit 0 跳过**：KSHELL bit 位是 1-based offset；shell-model-obs 的 `bitset<NMO>` 是 0-based。映射 `your_bit_index = (kshell_bit - 1) + orbit_offset_in_your_layout`
- **proton/neutron 分两个 mask**：KSHELL 是 `(mbp, mbn)` 双整数；shell-model-obs 可能是合并的单 bitset
- **`kwf` 4 vs 8**：从文件大小反推或用户指定
- **partition `min_m/max_m` 钳制**：必须严格按 partition.F90:161-173 的逻辑做，否则枚举顺序错位
