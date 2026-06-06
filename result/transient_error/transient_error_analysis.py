"""Transient-method error analysis (faithful 2-level reproduction)."""
import numpy as np

TWO_PI = 2 * np.pi
EC = 0.2 * TWO_PI
EJ0 = 10.0 * TWO_PI
T_PI2 = 10.0
OMEGA = (np.pi / 2.0) / T_PI2

a = np.array([[0.0, 1.0], [0.0, 0.0]], dtype=complex)
adag = a.conj().T
N_op = adag @ a
psi0 = np.array([1.0, 0.0], dtype=complex)


def omega_q(flux):
    EJ = EJ0 * abs(np.cos(np.pi * flux))
    return np.sqrt(8 * EJ * EC) - EC


def H_ctrl(phase):
    return OMEGA / 2.0 * (a * np.exp(1j * phase) + adag * np.exp(-1j * phase))


def expm_2x2(M):
    w, V = np.linalg.eig(M)
    return V @ np.diag(np.exp(w)) @ np.linalg.inv(V)


def p_e_sequence(delta, phase2):
    H1 = delta * N_op + H_ctrl(np.pi / 2.0)
    H2 = delta * N_op + H_ctrl(phase2)
    U1 = expm_2x2(-1j * H1 * T_PI2)
    U2 = expm_2x2(-1j * H2 * T_PI2)
    psi = U2 @ (U1 @ psi0)
    return abs(psi[1]) ** 2


def p_diff(delta):
    return (p_e_sequence(delta, 0.0) - p_e_sequence(delta, np.pi)) / 2.0


def extract_G1_G3(h=1e-3):
    d = np.array([-3, -2, -1, 1, 2, 3]) * h
    pd = np.array([p_diff(x) for x in d])
    A = np.column_stack([d, d**3])
    coef, *_ = np.linalg.lstsq(A, pd, rcond=None)
    return coef[0], 6.0 * coef[1]


def transient_measure(delta_true, G1, G3, order=1):
    pd = p_diff(delta_true)
    dw = pd / G1
    if order >= 3 and abs(G3) > 1e-10:
        for _ in range(50):
            f = G1 * dw + (G3 / 6.0) * dw**3 - pd
            fp = G1 + (G3 / 2.0) * dw**2
            if abs(fp) < 1e-15:
                break
            dw_new = dw - f / fp
            if abs(dw_new - dw) < 1e-14 * max(abs(dw), 1e-14):
                dw = dw_new
                break
            dw = dw_new
    return dw


def run_flux_point(flux_bias):
    omega_d = omega_q(flux_bias)
    G1, G3 = extract_G1_G3()
    flux_scan = np.linspace(0, 0.01, 15)
    rows = []
    for fa in flux_scan:
        f_true = omega_q(flux_bias + fa)
        delta_true = f_true - omega_d
        dw1 = transient_measure(delta_true, G1, G3, order=1)
        dw3 = transient_measure(delta_true, G1, G3, order=3)
        err_pred = (G3 / (6.0 * G1)) * delta_true**3
        rows.append((fa, flux_bias + fa, delta_true, f_true,
                     omega_d + dw1, omega_d + dw3, err_pred))
    return G1, G3, omega_d, np.array(rows)


if __name__ == "__main__":
    for bias in (0.0, 0.9):
        G1, G3, wd, R = run_flux_point(bias)
        print(f"\n===== flux bias = {bias} =====")
        print(f"omega_d={wd/TWO_PI:.6f} GHz; G1={G1:.4f} G3={G3:.4f} "
              f"G3/G1={G3/G1:.4f}")
        delta = R[:, 2]; f_true = R[:, 3]
        err_o1 = (R[:, 4] - f_true) / TWO_PI * 1e3
        err_o3 = (R[:, 5] - f_true) / TWO_PI * 1e3
        err_pred = R[:, 6] / TWO_PI * 1e3
        print(f"{'fa':>7}{'Delta(MHz)':>12}{'err_o1':>11}{'err_pred':>11}"
              f"{'err_o3':>12}")
        for i in range(len(R)):
            print(f"{R[i,0]:7.4f}{delta[i]/TWO_PI*1e3:12.3f}"
                  f"{err_o1[i]:11.4f}{err_pred[i]:11.4f}{err_o3[i]:12.6f}")
        mask = np.abs(err_pred) > 1e-6
        if mask.any():
            ratio = err_o1[mask] / err_pred[mask]
            print(f"  mean(err_o1/err_pred)={ratio.mean():.4f} "
                  f"std={ratio.std():.4f}")
        rms1 = np.sqrt(np.mean(err_o1**2)); rms3 = np.sqrt(np.mean(err_o3**2))
        print(f"  RMS: o1={rms1:.4f} MHz o3={rms3:.6f} MHz "
              f"reduction x{rms1/max(rms3,1e-12):.1f}")
