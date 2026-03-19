"""
Propagación Acústica en el Océano
==================================
Implementación Python del código Fortran acustica.f90.

Métodos implementados:
  - Trazado de rayos : RK45 adaptativo  (scipy.integrate.solve_ivp)
  - Modos normales   : Método de disparo con algoritmo de Numerov
                       + brentq (método de Brent, scipy.optimize)

Perfil de velocidad: Munk  c0=1500 m/s, eps=0.00737, zM=1300 m
Condiciones de frontera: ψ(0)=0, ψ(D)=0
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")   # sin GUI; cambiar a "TkAgg" si se quiere ventana
import matplotlib.pyplot as plt

# ===========================================================
#  PERFIL DE VELOCIDAD DE MUNK  (ec. 33 del Fortran)
# ===========================================================
C0  = 1500.0   # m/s
EPS = 0.00737
ZM  = 1300.0   # m  (eje del canal SOFAR)

def c_munk(z):
    """Velocidad del sonido [m/s] — perfil de Munk (vectorizable)."""
    g = 2.0 * (np.asarray(z) - ZM) / ZM
    return C0 * (1.0 + EPS * (g - 1.0 + np.exp(-g)))

def dc_munk(z):
    """Derivada dc/dz del perfil de Munk (vectorizable)."""
    g = 2.0 * (np.asarray(z) - ZM) / ZM
    return (2.0 * C0 / ZM) * EPS * (1.0 - np.exp(-g))


# ===========================================================
#  TRAZADO DE RAYOS — RK45 ADAPTATIVO
# ===========================================================
# Estado Y = [x, z, px, pz]
#   px = cos θ / c(z),   pz = sin θ / c(z)
# Ecuaciones (ec. 7-8):
#   dx/ds = c·px,  dz/ds = c·pz
#   dpx/ds = 0,    dpz/ds = −dc/dz / c²

def _derivs_ray(s, Y):
    x, z, px, pz = Y
    cs = float(c_munk(z))
    return [cs * px, cs * pz, 0.0, -float(dc_munk(z)) / cs**2]


def _reflect(Y, D):
    """Reflexión en superficie (z≤0) o fondo (z≥D): invierte pz."""
    x, z, px, pz = Y
    cs = float(c_munk(max(z, 0.0)))
    # Conserva magnitud de los slowness, invierte componente vertical
    return [x, np.clip(z, 0.0, D), px, -pz]


def trazar_rayo(theta_deg, z0, D, r_max, rtol=1e-8, atol=1e-10):
    """
    Integra una trayectoria de rayo con RK45 adaptativo.

    Parámetros
    ----------
    theta_deg : ángulo de emisión [°] (+ hacia abajo)
    z0        : profundidad de la fuente [m]
    D         : profundidad del fondo [m]
    r_max     : rango máximo [m]

    Retorna
    -------
    xs, zs : arrays (rango, profundidad) a lo largo del rayo
    """
    theta = np.radians(theta_deg)
    cs0   = float(c_munk(z0))
    Y     = np.array([0.0, z0,
                      np.cos(-theta) / cs0,
                      np.sin(-theta) / cs0])

    xs, zs = [Y[0]], [Y[1]]
    s      = 0.0
    s_max  = 5.0 * r_max

    def hit_surface(s, Y): return Y[1]
    def hit_bottom(s, Y):  return Y[1] - D
    hit_surface.terminal  = True;  hit_surface.direction = -1
    hit_bottom.terminal   = True;  hit_bottom.direction  =  1

    while s < s_max and Y[0] < r_max:
        sol = solve_ivp(
            _derivs_ray, [s, s_max], Y,
            method="RK45",
            events=[hit_surface, hit_bottom],
            rtol=rtol, atol=atol,
            max_step=50.0,
        )
        xs.extend(sol.y[0, 1:].tolist())
        zs.extend(sol.y[1, 1:].tolist())

        if sol.status == 1:          # rebote
            Y = np.array(_reflect(sol.y[:, -1], D))
            s = sol.t[-1]
        else:
            break

    return np.array(xs), np.array(zs)


# ===========================================================
#  MODOS NORMALES — MÉTODO DE DISPARO CON NUMEROV
# ===========================================================
# Ecuación de Helmholtz 1-D (ec. 28-29):
#   ψ'' + Q(z) ψ = 0,   Q(z) = (ω/c(z))² − k²
#
# Método de Numerov (O(h⁴)) para y'' = f(x)·y:
#   ψ_{n+1} = [2ψ_n(1 + 5h²Qn/12) − ψ_{n-1}(1 − h²Q_{n-1}/12)]
#              / (1 − h²Q_{n+1}/12)
#
# Condiciones iniciales: ψ(0) = 0,  ψ(h) = h  (disparo desde z=0)
# La función de disparo F(k) = ψ(D; k) se anula en los eigenvalores.

def _numerov_psi_D(k, omega, z_arr, Q_cache=None):
    """
    Integra ψ desde z=0 a z=D con Numerov; retorna ψ(D).
    Si Q_cache es None, calcula Q(z) = (ω/c(z))² − k²  internamente.
    """
    N  = len(z_arr)
    h  = z_arr[1] - z_arr[0]
    h2 = h * h

    # Q(z; k) = (ω/c)² − k²
    Q = (omega / c_munk(z_arr))**2 - k**2   # array vectorizado

    # Arranque
    psi_prev = 0.0        # ψ(z=0) = 0
    psi_curr = h          # ψ(z=h) = h  (derivada inicial = 1, escala arbitraria)

    f_prev = Q[0]
    f_curr = Q[1]

    for i in range(2, N):
        f_next = Q[i]
        denom  = 1.0 - h2 * f_next / 12.0
        psi_next = (2.0 * psi_curr * (1.0 + 5.0 * h2 * f_curr / 12.0)
                    - psi_prev * (1.0 - h2 * f_prev / 12.0)) / denom
        psi_prev = psi_curr
        psi_curr = psi_next
        f_prev   = f_curr
        f_curr   = f_next

    return psi_curr   # ψ(D)


def buscar_eigenvalores(omega, D, N_pts=1000, N_scan=3500, tol=1e-12):
    """
    Encuentra los eigenvalores kn en el rango físico usando:
      1. Barrido de F(k) = ψ(D; k)  para detectar cambios de signo  (ZBRAK)
      2. brentq de scipy para refinar cada raíz               (ZBRENT)

    Parámetros
    ----------
    N_pts  : puntos de z para Numerov  (equiv. NPT del Fortran)
    N_scan : subintervalos del barrido  (equiv. N1=3500 del Fortran)
    """
    # Rango de búsqueda: igual que Fortran (x1 = w/vson(1300), x2 = w/vson(5000))
    k1 = omega / float(c_munk(ZM))       # k mínimo (velocidad máxima en SOFAR)
    k2 = omega / float(c_munk(5000.0))   # k máximo (velocidad en el fondo)

    z_arr = np.linspace(0.0, D, N_pts + 1)

    # Barrido grueso
    ks   = np.linspace(k1, k2, N_scan + 1)
    vals = np.array([_numerov_psi_D(k, omega, z_arr) for k in ks])

    eigenvalores = []
    for i in range(N_scan):
        if vals[i] * vals[i + 1] < 0.0:
            kn = brentq(
                _numerov_psi_D, ks[i], ks[i + 1],
                args=(omega, z_arr),
                xtol=tol,
            )
            eigenvalores.append(kn)

    return np.array(eigenvalores), z_arr


def calcular_modo(kn, omega, z_arr):
    """
    Integra ψ(z) para el eigenvalor kn con Numerov y normaliza
    ∫|ψ|² dz = 1.  Retorna (z_arr, psi_normalizada).
    """
    N  = len(z_arr)
    h  = z_arr[1] - z_arr[0]
    h2 = h * h
    Q  = (omega / c_munk(z_arr))**2 - kn**2

    psi = np.zeros(N)
    psi[0] = 0.0
    psi[1] = h

    for i in range(2, N):
        denom   = 1.0 - h2 * Q[i] / 12.0
        psi[i]  = (2.0 * psi[i-1] * (1.0 + 5.0 * h2 * Q[i-1] / 12.0)
                   - psi[i-2] * (1.0 - h2 * Q[i-2] / 12.0)) / denom

    norm = np.sqrt(np.trapezoid(psi**2, z_arr))
    if norm > 0:
        psi /= norm
    return psi


# ===========================================================
#  PÉRDIDA POR TRANSMISIÓN  (ec. 20 del Fortran)
# ===========================================================

def perdida_transmision(r_arr, psi_src, psi_rec, kn_vals, omega):
    """
    TL [dB] para un array de rangos usando la suma modal vectorizada.

    p(r) = sqrt(2π/r) · Σ_n ψ_src,n · ψ_rec,n · exp(i·kn·r) / sqrt(kn)
    TL   = −20 log10|p(r)|
    """
    r_arr = np.asarray(r_arr)
    # Matriz (n_modos × n_rangos) de contribuciones
    A     = (psi_src * psi_rec / np.sqrt(kn_vals))[:, None]  # (M,1)
    phase = np.exp(1j * kn_vals[:, None] * r_arr[None, :])   # (M,R)
    p     = np.sqrt(2.0 * np.pi / r_arr) * np.sum(A * phase, axis=0)
    return -20.0 * np.log10(np.abs(p) + 1e-30)


# ===========================================================
#  PROGRAMA PRINCIPAL
# ===========================================================

def main():
    # ---- Parámetros ----
    D             = 5000.0   # profundidad del fondo [m]
    z0            = 200.0    # profundidad de la fuente [m]
    freq          = 50.0     # frecuencia [Hz]
    omega         = 2.0 * np.pi * freq
    r_max         = 100e3    # rango máximo [m]
    theta_min     = -10      # ángulo mínimo de emisión [°]
    theta_max     =  10      # ángulo máximo de emisión [°]
    N_pts         = 1000     # puntos de z (equiv. NPT)
    N_scan        = 3500     # subintervalos del barrido (equiv. N1)

    print("=" * 60)
    print("  PROPAGACIÓN ACÚSTICA — Python / RK45 + Shooting (Numerov)")
    print("=" * 60)
    print(f"  D={D} m | z_src={z0} m | f={freq} Hz | r_max={r_max/1e3} km")
    print()

    # =========================================================
    # 1. TRAZADO DE RAYOS — RK45 adaptativo
    # =========================================================
    print("[1] Trazado de rayos (RK45 adaptativo)...")
    fig_r, ax_r = plt.subplots(figsize=(13, 5))

    for theta in range(theta_min, theta_max + 1):
        if theta == 0:
            continue
        xs, zs = trazar_rayo(theta, z0, D, r_max)
        ax_r.plot(xs / 1e3, zs, lw=0.7, alpha=0.85)

    ax_r.set_xlim(0, r_max / 1e3)
    ax_r.set_ylim(D, 0)
    ax_r.set_xlabel("Rango [km]")
    ax_r.set_ylabel("Profundidad [m]")
    ax_r.set_title(f"Trazado de Rayos — RK45 adaptativo  (f = {freq} Hz)")
    ax_r.axhline(0, color="steelblue", lw=1.5, label="Superficie")
    ax_r.axhline(D, color="saddlebrown", lw=1.5, label="Fondo")
    ax_r.legend(fontsize=8)
    fig_r.tight_layout()
    fig_r.savefig("rayos_rk45.png", dpi=150)
    print("   → rayos_rk45.png")

    # =========================================================
    # 2. MODOS NORMALES — Método de disparo (Numerov + Brent)
    # =========================================================
    print(f"\n[2] Buscando eigenvalores kn (Numerov + Brent, N_scan={N_scan})...")
    kn_vals, z_arr = buscar_eigenvalores(omega, D, N_pts=N_pts,
                                         N_scan=N_scan, tol=1e-12)
    M = len(kn_vals)
    print(f"   → {M} modos encontrados")
    print(f"\n   {'Modo':>5}  {'kn [1/m]':>14}  {'cp [m/s]':>12}")
    print("   " + "-" * 36)
    for j, kn in enumerate(kn_vals):
        print(f"   {j+1:5d}  {kn:14.8f}  {omega/kn:12.4f}")

    with open("modos_shooting.txt", "w") as fout:
        fout.write("# Modo   kn[1/m]           cp[m/s]\n")
        for j, kn in enumerate(kn_vals):
            fout.write(f"{j+1:4d}  {kn:.10f}  {omega/kn:.6f}\n")
    print("   → modos_shooting.txt")

    # =========================================================
    # 3. FUNCIONES MODALES — Numerov para cada modo
    # =========================================================
    print("\n[3] Calculando funciones modales...")
    psi_modes = np.zeros((M, len(z_arr)))   # psi_modes[j, i] = ψ_j(z_i)
    for j, kn in enumerate(kn_vals):
        psi_modes[j] = calcular_modo(kn, omega, z_arr)

    # Interpolación en z_src y en profundidades de receptor (cada 100 m)
    psi_src    = np.array([np.interp(z0, z_arr, psi_modes[j]) for j in range(M)])
    z_rec_arr  = np.arange(100.0, D, 100.0)
    psi_rec_mat = np.array([
        [np.interp(zr, z_arr, psi_modes[j]) for j in range(M)]
        for zr in z_rec_arr
    ])   # shape: (n_rec, M)

    # Figura de los primeros 6 modos
    fig_m, ax_m = plt.subplots(figsize=(7, 8))
    for j in range(min(6, M)):
        ax_m.plot(psi_modes[j], z_arr, label=f"Modo {j+1}  k={kn_vals[j]:.5f}")
    ax_m.axhline(z0, color="r", ls="--", lw=1, label=f"z_src={z0} m")
    ax_m.set_ylim(D, 0)
    ax_m.set_xlabel("ψ(z) [normalizado]")
    ax_m.set_ylabel("Profundidad [m]")
    ax_m.set_title(f"Funciones modales — Shooting/Numerov  (f = {freq} Hz)")
    ax_m.legend(fontsize=8)
    fig_m.tight_layout()
    fig_m.savefig("modos_shooting.png", dpi=150)
    print("   → modos_shooting.png")

    # =========================================================
    # 4. PÉRDIDA POR TRANSMISIÓN — suma modal vectorizada
    # =========================================================
    print("\n[4] Calculando pérdida por transmisión...")
    r_arr = np.arange(1.0, r_max + 1.0, 1.0)

    fig_tl, ax_tl = plt.subplots(figsize=(13, 5))
    with open("tl_modos.txt", "w") as fout:
        fout.write("# Rango[km]  Profundidad[m]  TL[dB]\n")
        for ir, zr in enumerate(z_rec_arr):
            TL = perdida_transmision(r_arr, psi_src, psi_rec_mat[ir],
                                     kn_vals, omega)
            # Solo cada 200 m para la figura (menos curvas)
            if ir % 2 == 0:
                ax_tl.plot(r_arr / 1e3, TL, lw=0.6, alpha=0.75,
                           label=f"z={zr:.0f} m")
            # Guarda igual que el Fortran: cada 100 m en rango
            for i in range(99, len(r_arr), 100):
                fout.write(f"{r_arr[i]/1e3:.3f}  {-zr:.1f}  {TL[i]:.4f}\n")
            fout.write("\n")

    ax_tl.set_xlim(0, r_max / 1e3)
    ax_tl.set_xlabel("Rango [km]")
    ax_tl.set_ylabel("TL [dB]")
    ax_tl.set_title(f"Pérdida por Transmisión — Modos Normales  (f = {freq} Hz)")
    ax_tl.legend(fontsize=6, ncol=3)
    fig_tl.tight_layout()
    fig_tl.savefig("tl_modos.png", dpi=150)
    print("   → tl_modos.png")
    print("   → tl_modos.txt")

    print("\nListo.")
    plt.close("all")


if __name__ == "__main__":
    main()
