#!/usr/bin/env python3
"""Convert KSHELL/IMSRG .snt (pn formalism) -> BIGSTICK .sps (pns) + .int (xpn).

Usage:
    snt_to_bigstick.py input.snt out_basename

Writes:
    out_basename.sps         BIGSTICK pns-formalism single-particle space
    out_basename.int         BIGSTICK xpn-format interaction (auto-scaling if input
                             .snt has TBME mass dependence method2=1)

Notes:
  - Requires numorb_p == numorb_n with identical orbit sets in identical order
    (BIGSTICK pns format constraint; nearly always true for sd/pf/etc.).
  - TBMEs are canonicalised to (a<=b, c<=d) using the KSHELL swap phase
    (-1)^{(j_a+j_b)/2 - J + 1}.  T field is then always set to 1, so that
    BIGSTICK's internal swap phase (-1)^{J+T+(j+j')/2} matches in (-1)^{2J}=1.
  - If the .snt TBME header specifies method2=1 (mass dependence (mass/im0)^pwr),
    the converter emits BIGSTICK xpn autoscale header (negative nme + Acore, Aref,
    pwr_BS = -pwr_KS).  BIGSTICK then re-applies the same scaling at runtime.
  - method2=0 (no mass dependence) emits a positive-nme header with TBME values as-is;
    when running BIGSTICK enter `1.0 1.0 0.0 0.0` then `1 1 1 1 1` for scaling.
  - method2=10 (no-core kinetic-energy term) is not supported.
"""
import sys

if len(sys.argv) < 3:
    print(__doc__)
    sys.exit(1)
snt, base = sys.argv[1], sys.argv[2]

# ---------- strip comments, keep token-rich lines (preserve all fields incl method2) ----------
raw_lines = []
for raw in open(snt):
    s = raw.split('!')[0].split('#')[0].strip()
    if s:
        raw_lines.append(s)
it = iter(raw_lines)

# ---------- header: numorb_p, numorb_n, ncore_p, ncore_n ----------
hdr = next(it).split()
nP, nN, ncP, ncN = (int(x) for x in hdr[:4])
assert nP == nN, f"BIGSTICK pns requires numorb_p == numorb_n; got {nP} vs {nN}"

p_orbs, n_orbs = [], []
for _ in range(nP):
    parts = next(it).split()
    _, n, l, j2, _ = (int(x) for x in parts[:5])
    p_orbs.append((n, l, j2))
for _ in range(nN):
    parts = next(it).split()
    _, n, l, j2, _ = (int(x) for x in parts[:5])
    n_orbs.append((n, l, j2))
assert p_orbs == n_orbs, "proton & neutron orbit sets must match in same order"

# ---------- SPE section: n_spe  method1  [hw]  then  a b val ----------
spe_hdr = next(it).split()
nSPE = int(spe_hdr[0])
spe = [0.0] * (nP + nN)
spe_off = []
for _ in range(nSPE):
    parts = next(it).split()
    a, b, val = int(parts[0]), int(parts[1]), float(parts[2])
    if a == b:
        spe[a - 1] = val
    else:
        spe_off.append((a, b, val))
if spe_off:
    print(f"WARNING: {len(spe_off)} off-diagonal SPE entries ignored "
          f"(BIGSTICK xpn header only takes diagonal SPE).", file=sys.stderr)

# ---------- TBME header: n_tbme  method2  [im0 pwr]  or  [hw] ----------
tbme_hdr = next(it).split()
nTBME = int(tbme_hdr[0])
method2 = int(tbme_hdr[1]) if len(tbme_hdr) >= 2 else 0
im0 = pwr = None
if method2 == 1:
    im0 = int(tbme_hdr[2])
    pwr = float(tbme_hdr[3])
elif method2 == 10:
    raise NotImplementedError("method2=10 (no-core hw/A kinetic term) not supported")
elif method2 != 0:
    raise NotImplementedError(f"method2={method2} not supported")

all_orbs = p_orbs + n_orbs
j2_of = [0] + [orb[2] for orb in all_orbs]  # 1-based: j2_of[idx] = 2*j

def swap_with_phase(a, b, J):
    """KSHELL canonical swap: isign = (-1)^{(j_a+j_b)/2 - J + 1}.  Return (a, b, phase)."""
    if a > b:
        a, b = b, a
        return a, b, ((-1) ** ((j2_of[a] + j2_of[b]) // 2 - J + 1))
    return a, b, 1

tbmes = []
for _ in range(nTBME):
    parts = next(it).split()
    a, b, c, d, J, val = (int(parts[0]), int(parts[1]), int(parts[2]),
                          int(parts[3]), int(parts[4]), float(parts[5]))
    a, b, ph1 = swap_with_phase(a, b, J)
    c, d, ph2 = swap_with_phase(c, d, J)
    tbmes.append((a, b, c, d, J, val * ph1 * ph2))

# ---------- write .sps : pns formalism, proton block then neutron block ----------
with open(base + ".sps", "w") as f:
    f.write("pns\n")
    f.write(f"{nP}\n")
    for (n, l, j2) in p_orbs:
        f.write(f"  {n}.0  {l}.0  {j2/2:.1f}  2\n")
    for (n, l, j2) in n_orbs:
        f.write(f"  {n}.0  {l}.0  {j2/2:.1f}  2\n")

# ---------- write .int : BIGSTICK xpn format ----------
#  - autoscale header (-nme + 3 extras: Acore, Aref, pwr_BS):
#      KSHELL    : v *= (A/im0)**pwr_KS
#      BIGSTICK  : v *= (Aref/A)**pwr_BS = (Aref/A)**(-pwr_KS)
#    so pwr_BS = -pwr_KS, Aref = im0, Acore = ncP + ncN.
def emit_int_autoscale(f):
    A_core = ncP + ncN
    A_ref = im0
    pwr_BS = -pwr
    nme_signed = -len(tbmes)
    f.write(f"{nme_signed}  "
            + "  ".join(f"{v:.8f}" for v in spe)
            + f"  {A_core}  {A_ref}  {pwr_BS:.6f}\n")

def emit_int_noscale(f):
    f.write(f"{len(tbmes)}  " + "  ".join(f"{v:.8f}" for v in spe) + "\n")

with open(base + ".int", "w") as f:
    f.write("!XPN converted from KSHELL .snt by util/snt_to_bigstick.py\n")
    f.write(f"!     proton orbits 1..{nP}, neutron orbits {nP + 1}..{nP + nN}\n")
    if method2 == 1:
        f.write(f"!     autoscale enabled: TBME *= ({im0}/A)^{-pwr:.4f}  (Acore={ncP + ncN})\n")
        emit_int_autoscale(f)
    else:
        f.write("!     no mass dependence; BIGSTICK will prompt for scaling (use 1 1 0 0)\n")
        emit_int_noscale(f)
    for (a, b, c, d, J, val) in tbmes:
        f.write(f"  {a:3d} {b:3d} {c:3d} {d:3d}  {J:3d}    1  {val:.8f}\n")

print(f"wrote {base}.sps and {base}.int")
print(f"  nP={nP} nN={nN} ncoreP={ncP} ncoreN={ncN} nSPE={nSPE} nTBME={nTBME}",
      f"method2={method2}",
      (f"im0={im0} pwr={pwr}" if method2 == 1 else ""))
