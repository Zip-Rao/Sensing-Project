"""Root-cause check: the code's G3 (diagonal single-time kernel integral)
vs the genuine cubic Taylor coefficient of p_diff(Delta).

Both are pure-numpy and exact for a 2-level system. We:
  (A) fit p_diff(Delta) over a wide range  -> TRUE cubic coeff (defines
      G3 via p_diff = G1*D + (1/6)*G3_true*D^3)
  (B) replicate frequency.py's path: build the *diagonal* omega-kernel
      k3(t) by perturbing the detuning at ONE time point t_j with a narrow
      Gaussian and taking the 5-point FD third derivative w.r.t. amplitude;
      integrate -> G3_diag = int k3^diag dt   (this is what the code plugs
      into the cubic-Newton inversion).
Then compare.
"""
import numpy as np

TWO_PI = 2 * np.pi
T_PI2 = 10.0
OMEGA = (np.pi / 2.0) / T_PI2
DT = 0.5
SIGMA_T = 2.0 * DT

a = np.array([[0.0, 1.0], [0.0, 0.0]], dtype=complex)
adag = a.conj().T
N_op = adag @ a                       # = diag(0,1); sigma_z/2 + const
sz = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
psi0 = np.array([1.0, 0.0], dtype=complex)
P_e = np.array([[0.0, 0.0], [0.0, 1.0]], dtype=complex)


def H_ctrl(phase):
    return OMEGA / 2.0 * (a * np.exp(1j * phase) + adag * np.exp(-1j * phase))


def expm2(M):
    w, V = np.linalg.eig(M)
    return V @ np.diag(np.exp(w)) @ np.linalg.inv(V)


# Time grid spanning the two pi/2 pulses (Y then X), tau=0.
tgrid = np.arange(0.0, 2 * T_PI2, DT)


def pulse_phase(t):
    """phase1=pi/2 on pulse 1 (Y), phase2=0 on pulse 2 (X)."""
    return np.pi / 2.0 if t < T_PI2 else 0.0


def propagate(delta, vz_amp=0.0, vz_center=None):
    """Step-evolve |0> through both pulses with constant detuning `delta`
    plus an optional narrow Gaussian detuning stimulus (the VZ omega probe)
    centred at vz_center.  Return p_e at the end."""
    psi = psi0.copy()
    for t in tgrid:
        d_eff = delta
        if vz_amp != 0.0 and vz_center is not None:
            g = np.exp(-0.5 * ((t - vz_center) / SIGMA_T) ** 2)
            norm = 1.0 / (SIGMA_T * np.sqrt(2 * np.pi))
            d_eff = delta + vz_amp * norm * g
        H = d_eff * N_op + H_ctrl(pulse_phase(t))
        U = expm2(-1j * H * DT)
        psi = U @ psi
    return float(np.real(np.conj(psi) @ (P_e @ psi)))


def p_diff(delta, **kw):
    pX = propagate(delta, **kw)
    # second pulse phase pi  -> -X. emulate by flipping pulse-2 phase.
    global pulse_phase
    return None  # placeholder, replaced below


# --- p_diff with selectable second-pulse phase --------------------------
def propagate_seq(delta, phase2, vz_amp=0.0, vz_center=None):
    psi = psi0.copy()
    for t in tgrid:
        ph = np.pi / 2.0 if t < T_PI2 else phase2
        d_eff = delta
        if vz_amp != 0.0 and vz_center is not None:
            g = np.exp(-0.5 * ((t - vz_center) / SIGMA_T) ** 2)
            norm = 1.0 / (SIGMA_T * np.sqrt(2 * np.pi))
            d_eff = delta + vz_amp * norm * g
        H = d_eff * N_op + H_ctrl(ph)
        U = expm2(-1j * H * DT)
        psi = U @ psi
    return float(np.real(np.conj(psi) @ (P_e @ psi)))


def pdiff(delta, vz_amp=0.0, vz_center=None):
    pX = propagate_seq(delta, 0.0, vz_amp, vz_center)
    pmX = propagate_seq(delta, np.pi, vz_amp, vz_center)
    return (pX - pmX) / 2.0


# ===== (A) TRUE cubic Taylor coefficient ================================
d = np.linspace(-0.3, 0.3, 41)
d = d[d != 0]
pd = np.array([pdiff(x) for x in d])
A = np.column_stack([d, d**3, d**5, d**7])
coef, *_ = np.linalg.lstsq(A, pd, rcond=None)
G1_taylor = coef[0]
G3_taylor = 6.0 * coef[1]          # since p_diff = G1*D + (1/6)G3*D^3
print("=== (A) genuine Taylor coefficients of p_diff(Delta) ===")
print(f"  G1_taylor = {G1_taylor:.4f}")
print(f"  G3_taylor = {G3_taylor:.4f}   (=6*cubic-fit-coeff)")
print(f"  G3/G1     = {G3_taylor/G1_taylor:.3f} ns^2")

# ===== (B) diagonal single-time omega kernel (the code's G3) ============
# Replicate _extract_kn_omega: at each t_j, 5 VZ amplitudes {-2h..2h},
# 5-point stencils for k1 (1st deriv) and k3 (3rd deriv) of p_diff wrt amp.
h = 0.005
k1_diag = np.zeros(len(tgrid))
k3_diag = np.zeros(len(tgrid))
for i, tj in enumerate(tgrid):
    pv = np.array([pdiff(0.0, vz_amp=m * h, vz_center=tj)
                   for m in (-2, -1, 0, 1, 2)])
    fm2, fm1, f0, fp1, fp2 = pv
    k1_diag[i] = (fm2 - 8*fm1 + 8*fp1 - fp2) / (12 * h)
    k3_diag[i] = (-fm2 + 2*fm1 - 2*fp1 + fp2) / (2 * h**3)

G1_diag = np.trapezoid(k1_diag, tgrid)
G3_diag = np.trapezoid(k3_diag, tgrid)
print("\n=== (B) diagonal single-time kernel integrals (code's path) ===")
print(f"  G1_diag = {G1_diag:.4f}   (should match G1_taylor)")
print(f"  G3_diag = {G3_diag:.4f}   (code feeds THIS into cubic Newton)")
print(f"  G3_diag/G1_diag = {G3_diag/G1_diag:.4f}   (the doc's '|G3/G1|=1')")

print("\n=== ratio of the two 'G3' objects ===")
print(f"  G3_taylor / G3_diag = {G3_taylor/G3_diag:.2f} ns^2")
print(f"  sqrt(|that|)        = {np.sqrt(abs(G3_taylor/G3_diag)):.2f} ns "
      f"(~ effective pulse time; 2*T_pi2={2*T_PI2} ns)")
