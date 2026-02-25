"""
generar_figuras.py
Genera las figuras de los tres pozos potenciales con sus niveles de energía.
Uso: python3 generar_figuras.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Estilo global ────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'axes.linewidth': 1.2,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'legend.framealpha': 0.9,
})

# ── Parámetros de la malla ────────────────────────────────────────────────────
n = 2000
x = np.linspace(-4.0, 4.0, n + 1)
h = x[1] - x[0]

# ── Potenciales ───────────────────────────────────────────────────────────────
def V_caja(x):
    return np.where(np.abs(x) < 1.5, 0.0, 25.0)

def V_pozo_V(x):
    return 4.0 * np.abs(x)

def V_doble(x):
    return 0.5 * x**4 - 5.0 * x**2

# ── Numerov shooting ──────────────────────────────────────────────────────────
def solve_numerov(E, V_arr):
    k2 = 2.0 * (E - V_arr)
    p = np.zeros(n + 1)
    p[0] = 0.0
    p[1] = 1e-5
    h2_12 = h**2 / 12.0
    for j in range(1, n):
        p[j + 1] = (
            2.0 * (1.0 - 5.0 * h2_12 * k2[j]) * p[j]
            - (1.0 + h2_12 * k2[j - 1]) * p[j - 1]
        ) / (1.0 + h2_12 * k2[j + 1])
        if abs(p[j + 1]) > 1e10:
            p[j + 1] = 1e10   # tal como está en el código Fortran
    return p

# ── Bisección para encontrar eigenvalor ──────────────────────────────────────
def biseccion(E1, E2, V_arr):
    for _ in range(60):
        Em = (E1 + E2) / 2.0
        pf1 = solve_numerov(E1, V_arr)[n]
        pfm = solve_numerov(Em, V_arr)[n]
        if pf1 * pfm < 0:
            E2 = Em
        else:
            E1 = Em
    return (E1 + E2) / 2.0

# ── Buscar niveles ────────────────────────────────────────────────────────────
def encontrar_niveles(V_arr, E_start=-15.0, E_end=20.0, dE=0.2, max_niveles=8):
    niveles = []
    E1 = E_start
    while E1 < E_end and len(niveles) < max_niveles:
        E2 = E1 + dE
        pf1 = solve_numerov(E1, V_arr)[n]
        pf2 = solve_numerov(E2, V_arr)[n]
        if pf1 * pf2 < 0.0:
            em = biseccion(E1, E2, V_arr)
            niveles.append(em)
        E1 = E2
    return np.array(niveles)

# ── Funciones de onda normalizadas para graficar ──────────────────────────────
def get_psi(E, V_arr, escala=1.5):
    psi = solve_numerov(E, V_arr)
    # normalizar al rango [-1, 1]
    pmax = np.max(np.abs(psi))
    if pmax > 0:
        psi = psi / pmax
    return psi * escala   # escala para visualización

# ═════════════════════════════════════════════════════════════════════════════
# FIGURA 1: Pozo Cuadrado (Caja)
# ═════════════════════════════════════════════════════════════════════════════
print("Generando figura 1: Pozo cuadrado...")
V1 = V_caja(x)
niveles_caja = encontrar_niveles(V1)
print(f"  Niveles caja: {niveles_caja}")

fig, ax = plt.subplots(figsize=(7, 5))

# Sombreado interior del pozo
ax.fill_between(x, -2, 0, where=np.abs(x) < 1.5, color='lightblue', alpha=0.35, label='Interior del pozo')

# Potencial
ax.plot(x, V1, color='#1a237e', lw=2.5, label=r'$V(x)$')
ax.axhline(0, color='gray', lw=0.6, ls='--', alpha=0.5)

# Niveles de energía con funciones de onda
colores = plt.cm.Reds(np.linspace(0.45, 0.9, len(niveles_caja)))
for idx, E in enumerate(niveles_caja):
    psi = get_psi(E, V1, escala=min(1.2, (niveles_caja[-1] - niveles_caja[0]) / 8))
    ax.axhline(E, color=colores[idx], lw=1.3, ls='--', alpha=0.7)
    ax.plot(x, psi + E, color=colores[idx], lw=1.5)
    ax.text(3.25, E + 0.25, rf'$E_{{{idx+1}}}={E:.3f}$', color=colores[idx],
            fontsize=8.5, va='bottom')

ax.set_xlim(-4.2, 4.2)
ax.set_ylim(-2, 27)
ax.set_xlabel(r'Posición $x$ (u.a.)', fontsize=12)
ax.set_ylabel(r'Energía (u.a.)', fontsize=12)
ax.set_title('Pozo Cuadrado: $V(x) = 0$ para $|x| < 1.5$,\n$V(x) = V_0 = 25$ en otro caso', fontsize=11)
ax.legend(loc='upper center', fontsize=9)
ax.grid(True, which='both', ls=':', alpha=0.4)
plt.tight_layout()
plt.savefig('fig_caja.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Guardado: fig_caja.png")

# ═════════════════════════════════════════════════════════════════════════════
# FIGURA 2: Pozo en V (Lineal)
# ═════════════════════════════════════════════════════════════════════════════
print("Generando figura 2: Pozo en V...")
V2 = V_pozo_V(x)
niveles_V = encontrar_niveles(V2)
print(f"  Niveles V: {niveles_V}")

fig, ax = plt.subplots(figsize=(7, 5))

# Sombreado interior del pozo (debajo del nivel más alto encontrado)
ax.fill_between(x, V2, niveles_V[-1] + 2, color='lightyellow', alpha=0.4)

# Potencial
ax.plot(x, V2, color='#1b5e20', lw=2.5, label=r'$V(x) = 4|x|$')
ax.axhline(0, color='gray', lw=0.6, ls='--', alpha=0.5)

# Valores analíticos (funciones Airy): E_n = 2 * z_n
# zeros de Ai'(-z): 1.0188, 3.2482, 4.8201, 6.1633 (estados pares)
# zeros de Ai(-z):  2.3381, 4.0879, 5.5206, 6.7867  (estados impares)
zeros_par  = [1.01879, 3.24820, 4.82010, 6.16331]
zeros_impar= [2.33811, 4.08795, 5.52056, 6.78671]
todos_z = sorted(zeros_par[:4] + zeros_impar[:4])
E_analitico = [2 * z for z in todos_z]

colores = plt.cm.Greens(np.linspace(0.45, 0.9, len(niveles_V)))
for idx, E in enumerate(niveles_V):
    psi = get_psi(E, V2, escala=0.8)
    ax.axhline(E, color=colores[idx], lw=1.3, ls='--', alpha=0.7)
    ax.plot(x, psi + E, color=colores[idx], lw=1.5)
    ax.text(3.25, E + 0.25, rf'$E_{{{idx+1}}}={E:.3f}$', color=colores[idx],
            fontsize=8.5, va='bottom')

ax.set_xlim(-4.2, 4.2)
ax.set_ylim(-1, 17)
ax.set_xlabel(r'Posición $x$ (u.a.)', fontsize=12)
ax.set_ylabel(r'Energía (u.a.)', fontsize=12)
ax.set_title(r'Pozo en V: $V(x) = 4|x|$', fontsize=12)
ax.legend(loc='upper center', fontsize=9)
ax.grid(True, which='both', ls=':', alpha=0.4)
plt.tight_layout()
plt.savefig('fig_pozo_V.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Guardado: fig_pozo_V.png")

# ═════════════════════════════════════════════════════════════════════════════
# FIGURA 3: Doble Pozo Cuártico
# ═════════════════════════════════════════════════════════════════════════════
print("Generando figura 3: Doble pozo...")
V3 = V_doble(x)
niveles_dp = encontrar_niveles(V3)
print(f"  Niveles doble pozo: {niveles_dp}")

fig, ax = plt.subplots(figsize=(7, 5))

# Sombreado del interior de los dos pozos (debajo de E=0)
ax.fill_between(x, V3, 0, where=V3 < 0, color='lavender', alpha=0.5)
# Sombreado de la región barrera
ax.fill_between(x, 0, color='mistyrose', alpha=0.3, where=(V3 >= 0) & (np.abs(x) < 2.5))

# Potencial (recortar para visualización)
V3_plot = np.clip(V3, -15, 20)
ax.plot(x, V3_plot, color='#4a148c', lw=2.5, label=r'$V(x) = \frac{1}{2}x^4 - 5x^2$')

# Líneas de referencia
ax.axhline(0, color='firebrick', lw=1.0, ls=':', alpha=0.8, label='Tope de barrera ($E=0$)')
ax.axhline(-12.5, color='steelblue', lw=1.0, ls=':', alpha=0.8, label=r'Mínimo del pozo ($E=-12.5$)')

colores = plt.cm.Purples(np.linspace(0.45, 0.9, len(niveles_dp)))
for idx, E in enumerate(niveles_dp):
    psi = get_psi(E, V3, escala=0.8)
    ax.axhline(E, color=colores[idx], lw=1.3, ls='--', alpha=0.7)
    ax.plot(x, psi + E, color=colores[idx], lw=1.5)
    ax.text(3.25, E + 0.25, rf'$E_{{{idx+1}}}={E:.3f}$', color=colores[idx],
            fontsize=8.5, va='bottom')

ax.set_xlim(-4.2, 4.2)
ax.set_ylim(-14, 17)
ax.set_xlabel(r'Posición $x$ (u.a.)', fontsize=12)
ax.set_ylabel(r'Energía (u.a.)', fontsize=12)
ax.set_title(r'Doble Pozo: $V(x) = \frac{1}{2}x^4 - 5x^2$', fontsize=12)
ax.legend(loc='upper center', fontsize=9)
ax.grid(True, which='both', ls=':', alpha=0.4)
plt.tight_layout()
plt.savefig('fig_doble.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Guardado: fig_doble.png")

# ═════════════════════════════════════════════════════════════════════════════
# FIGURA 4: Comparación de los 3 pozos (panel 1×3)
# ═════════════════════════════════════════════════════════════════════════════
print("Generando figura 4: Panel comparativo...")
fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=False)

datos = [
    (V_caja(x),   niveles_caja, r'$V(x)=V_0\,\Theta(|x|-a)$', 'Pozo Cuadrado',
     '#1a237e', 'lightblue', (-2, 27)),
    (V_pozo_V(x), niveles_V,    r'$V(x) = 4|x|$',             'Pozo en V (Lineal)',
     '#1b5e20', 'lightyellow', (-1, 17)),
    (np.clip(V_doble(x), -15, 20), niveles_dp, r'$V(x) = \frac{1}{2}x^4 - 5x^2$', 'Doble Pozo Cuártico',
     '#4a148c', 'lavender', (-14, 17)),
]
cmaps = [plt.cm.Reds, plt.cm.Greens, plt.cm.Purples]

for ax, (V_plot, niv, formula, titulo, color, fill_c, ylim), cm in zip(axes, datos, cmaps):
    ax.plot(x, V_plot, color=color, lw=2.5)
    ax.axhline(0, color='gray', lw=0.6, ls='--', alpha=0.5)
    cols = cm(np.linspace(0.45, 0.9, len(niv)))
    for idx, E in enumerate(niv):
        ax.axhline(E, color=cols[idx], lw=1.5, ls='--', alpha=0.8)
        ax.text(3.5, E + 0.4, f'$E_{{{idx+1}}}$', color=cols[idx], fontsize=8)
    ax.set_title(f'{titulo}\n{formula}', fontsize=10, pad=4)
    ax.set_xlabel(r'$x$ (u.a.)', fontsize=11)
    ax.set_ylabel(r'$E$ (u.a.)', fontsize=11)
    ax.set_xlim(-4.2, 4.2)
    ax.set_ylim(*ylim)
    ax.grid(True, ls=':', alpha=0.4)

plt.suptitle('Comparación de Pozos Potenciales — Método de Numerov', fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig('fig_comparacion.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Guardado: fig_comparacion.png")
print("\nTodas las figuras generadas correctamente.")
