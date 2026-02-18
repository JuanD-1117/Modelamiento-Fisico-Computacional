"""
Taller 1 - Modelamiento Físico Computacional 2026-1
Modelo de Rashevsky: Dinámica de Inconformistas
Implementación Python: Euler Explícito, Taylor Orden 2, Trapecio Implícito
"""

import numpy as np
import matplotlib.pyplot as plt
import time

# ============================================================
# PARÁMETROS DEL MODELO
# ============================================================
p0    = 0.01
b     = 0.02
d     = 0.015
r     = 0.1
k     = r * b       # k = rb = 0.002  (d se cancela en la EDO)
h     = 1.0
t_max = 50
steps = int(t_max / h)
t_vals = np.linspace(0, t_max, steps + 1)

print("=" * 70)
print("TALLER 1 - Modelo de Rashevsky (Inconformistas)")
print(f"Parámetros: b={b}, d={d}, r={r}, k=rb={k}, h={h}, t_max={t_max}")
print("=" * 70)

# ============================================================
# SOLUCIÓN EXACTA:  p(t) = 1 - (1-p0)*exp(-k*t)
# ============================================================
p_exact = 1 - (1 - p0) * np.exp(-k * t_vals)

# ============================================================
# MÉTODO 1: Euler Explícito
#   p[n+1] = p[n] + h*k*(1 - p[n])
# ============================================================
p_euler = np.zeros(steps + 1)
p_euler[0] = p0
for i in range(steps):
    p_euler[i+1] = p_euler[i] + h * k * (1.0 - p_euler[i])

# ============================================================
# MÉTODO 2: Taylor Orden 2
#   p'  = k(1-p)
#   p'' = -k^2(1-p)
#   factor = h*k - (h^2*k^2)/2
#   p[n+1] = p[n] + (1-p[n]) * factor
# ============================================================
p_taylor = np.zeros(steps + 1)
p_taylor[0] = p0
factor_t2 = h * k - 0.5 * (h**2) * (k**2)
for i in range(steps):
    p_taylor[i+1] = p_taylor[i] + (1.0 - p_taylor[i]) * factor_t2

# ============================================================
# MÉTODO 3: Trapecio Implícito (Crank-Nicolson)
#   p[n+1] = p[n] + (h/2)*[k(1-p[n]) + k(1-p[n+1])]
#   Despejando: p[n+1] = (p[n]*(1 - hk/2) + hk) / (1 + hk/2)
# ============================================================
p_trap = np.zeros(steps + 1)
p_trap[0] = p0
hk2    = h * k / 2.0
a_trap = 1.0 - hk2
b_trap = h * k
c_trap = 1.0 + hk2
for i in range(steps):
    p_trap[i+1] = (p_trap[i] * a_trap + b_trap) / c_trap

# ============================================================
# ERRORES ABSOLUTOS
# ============================================================
err_euler  = np.abs(p_exact - p_euler)
err_taylor = np.abs(p_exact - p_taylor)
err_trap   = np.abs(p_exact - p_trap)

# ============================================================
# TABLA DE VALORES (cada 5 años)
# ============================================================
print("\nTabla de valores comparativa (cada 5 años):")
header = f"{'t':>4} | {'Exacta':>12} | {'Euler':>12} | {'Taylor 2':>12} | {'Trapecio':>12} | {'Err Euler':>10} | {'Err T2':>10} | {'Err Trap':>10}"
print(header)
print("-" * len(header))
for i in range(0, steps + 1, 5):
    t = t_vals[i]
    print(f"{t:>4.0f} | {p_exact[i]:>12.8f} | {p_euler[i]:>12.8f} | "
          f"{p_taylor[i]:>12.8f} | {p_trap[i]:>12.8f} | "
          f"{err_euler[i]:>10.3e} | {err_taylor[i]:>10.3e} | {err_trap[i]:>10.3e}")

print("\nResumen p(50):")
print(f"  Exacta   : {p_exact[-1]:.10f}")
print(f"  Euler    : {p_euler[-1]:.10f}  Error: {err_euler[-1]:.4e}")
print(f"  Taylor 2 : {p_taylor[-1]:.10f}  Error: {err_taylor[-1]:.4e}")
print(f"  Trapecio : {p_trap[-1]:.10f}  Error: {err_trap[-1]:.4e}")

# ============================================================
# GRÁFICA COMPARATIVA
# ============================================================
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 9), sharex=True)

ax1.plot(t_vals, p_exact,  'k-',   lw=2.5, label='Solución exacta', zorder=5)
ax1.plot(t_vals, p_euler,  'b--o', lw=1.5, ms=4, alpha=0.85, label='Euler explícito')
ax1.plot(t_vals, p_taylor, 'r--s', lw=1.5, ms=4, alpha=0.85, label='Taylor Orden 2')
ax1.plot(t_vals, p_trap,   'g--^', lw=1.5, ms=4, alpha=0.85, label='Trapecio implícito')
ax1.set_ylabel("$p(t)$ — Proporción de inconformistas", fontsize=12)
ax1.set_title("Modelo de Rashevsky: Comparación de métodos numéricos\n"
              r"$\frac{dp}{dt}=rb(1-p),\quad p(0)=0.01,\quad k=rb=0.002,\quad h=1$",
              fontsize=12)
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3)

ax2.semilogy(t_vals, err_euler  + 1e-16, 'b--o', lw=1.5, ms=4, label='Error Euler')
ax2.semilogy(t_vals, err_taylor + 1e-16, 'r--s', lw=1.5, ms=4, label='Error Taylor 2')
ax2.semilogy(t_vals, err_trap   + 1e-16, 'g--^', lw=1.5, ms=4, label='Error Trapecio')
ax2.set_xlabel("Tiempo $t$ (años)", fontsize=12)
ax2.set_ylabel("Error absoluto $|p_n - p(t_n)|$ (log)", fontsize=12)
ax2.set_title("Evolución del error absoluto por método", fontsize=12)
ax2.legend(fontsize=11)
ax2.grid(True, which='both', alpha=0.3)

plt.tight_layout()
plt.savefig("taller1_comparacion.png", dpi=150, bbox_inches='tight')
print("\nGráfica guardada: taller1_comparacion.png")
plt.show()

# ============================================================
# BENCHMARK: 10 millones de iteraciones por método
# ============================================================
ITERS = 10_000_000
print(f"\n{'=' * 60}")
print(f"BENCHMARK Python ({ITERS:,} iter × {steps} pasos)")
print(f"{'=' * 60}")

# --- Euler Explícito ---
start = time.perf_counter()
for _ in range(ITERS):
    p = p0
    for i in range(steps):
        p = p + h * k * (1.0 - p)
t1 = time.perf_counter() - start
print(f"Euler Explícito  : {t1:.4f} s  | p(50) = {p:.8f}")

# --- Taylor Orden 2 ---
ft2 = factor_t2
start = time.perf_counter()
for _ in range(ITERS):
    p = p0
    for i in range(steps):
        p = p + (1.0 - p) * ft2
t2 = time.perf_counter() - start
print(f"Taylor Orden 2   : {t2:.4f} s  | p(50) = {p:.8f}")

# --- Trapecio Implícito ---
start = time.perf_counter()
for _ in range(ITERS):
    p = p0
    for i in range(steps):
        p = (p * a_trap + b_trap) / c_trap
t3 = time.perf_counter() - start
print(f"Trapecio Impl.   : {t3:.4f} s  | p(50) = {p:.8f}")

print(f"{'=' * 60}")
print(f"TIEMPOS PYTHON [s]: Euler={t1:.4f}  Taylor2={t2:.4f}  Trapecio={t3:.4f}")
