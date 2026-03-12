"""
Solución numérica: Ecuaciones de Poisson y Laplace
Método de Diferencias Finitas — Iteración de Jacobi
Python serial (numpy)

Ejecutar: python poisson_serial.py [M] [N]   (default M=N=50)
"""

import sys
import time
import numpy as np

TOL      = 1e-6
MAX_ITER = 400000


def solve_poisson(case, M=50, N=50):
    """
    Resuelve la ecuación de Poisson/Laplace para el caso dado
    usando diferencias finitas con iteración de Jacobi.
    Retorna: (T_numerica, X, Y, iteraciones)
    """
    domains = {
        1: (0.0, 2.0, 0.0, 1.0),
        2: (1.0, 2.0, 0.0, 1.0),
        3: (1.0, 2.0, 0.0, 2.0),
        4: (1.0, 2.0, 1.0, 2.0),
    }
    x0, xf, y0, yf = domains[case]

    h  = (xf - x0) / M
    k  = (yf - y0) / N
    h2 = h * h
    k2 = k * k
    denom = 2.0 * (h2 + k2)

    # Grilla: X[i,j] = x[j], Y[i,j] = y[i]  (igual que Fortran)
    x = np.linspace(x0, xf, M + 1)
    y = np.linspace(y0, yf, N + 1)
    X, Y = np.meshgrid(x, y)   # shape (N+1, M+1)

    # Término fuente
    if case == 1:
        source = (X**2 + Y**2) * np.exp(X * Y)
    elif case == 2:
        source = np.zeros_like(X)
    elif case == 3:
        source = np.full_like(X, 4.0)
    else:  # case 4
        source = X / Y + Y / X

    # Inicializar T
    T = np.zeros((N + 1, M + 1))

    # Condiciones de contorno
    if case == 1:
        T[:, 0] = 1.0                       # V(x0,y) = 1
        T[:, M] = np.exp(2.0 * y)           # V(xf,y) = e^(2y)
        T[0, :] = 1.0                       # V(x,y0) = 1
        T[N, :] = np.exp(x)                 # V(x,yf) = e^x
    elif case == 2:
        T[:, 0] = np.log(y**2 + 1.0)       # V(1,y) = ln(y²+1)
        T[:, M] = np.log(y**2 + 4.0)       # V(2,y) = ln(y²+4)
        T[0, :] = 2.0 * np.log(x)          # V(x,0) = 2ln(x)
        T[N, :] = np.log(x**2 + 1.0)       # V(x,1) = ln(x²+1)
    elif case == 3:
        T[:, 0] = (x0 - y)**2              # V(1,y) = (1-y)²
        T[:, M] = (xf - y)**2              # V(2,y) = (2-y)²
        T[0, :] = (x - y0)**2              # V(x,0) = x²
        T[N, :] = (x - yf)**2              # V(x,2) = (x-2)²
    else:  # case 4
        T[:, 0] = y * np.log(y)            # V(1,y) = y·ln(y)
        T[:, M] = 2.0 * y * np.log(2.0 * y)  # V(2,y) = 2y·ln(2y)
        T[0, :] = x * np.log(x)            # V(x,1) = x·ln(x)
        T[N, :] = 2.0 * x * np.log(2.0 * x)  # V(x,2) = 2x·ln(2x)

    # Iteración de Jacobi
    iters = MAX_ITER
    for it in range(1, MAX_ITER + 1):
        T_new = T.copy()
        T_new[1:-1, 1:-1] = (
            (T[2:,  1:-1] + T[:-2, 1:-1]) * k2
          + (T[1:-1, 2:]  + T[1:-1, :-2]) * h2
          - source[1:-1, 1:-1] * h2 * k2
        ) / denom

        delta = np.max(np.abs(T_new[1:-1, 1:-1] - T[1:-1, 1:-1]))
        T = T_new
        if delta < TOL:
            iters = it
            break

    return T, X, Y, iters


def analytical(case, X, Y):
    if case == 1:
        return np.exp(X * Y)
    elif case == 2:
        return np.log(X**2 + Y**2)
    elif case == 3:
        return (X - Y)**2
    else:  # case 4
        return X * Y * np.log(X * Y)


def main():
    M = int(sys.argv[1]) if len(sys.argv) >= 2 else 50
    N = int(sys.argv[2]) if len(sys.argv) >= 3 else 50

    case_names = {
        1: "Poisson  f=(x²+y²)·e^(xy)  [0,2]x[0,1]",
        2: "Laplace  f=0                [1,2]x[0,1]",
        3: "Poisson  f=4                [1,2]x[0,2]",
        4: "Poisson  f=x/y+y/x         [1,2]x[1,2]",
    }

    print("=" * 60)
    print(f"  Jacobi FDM — Poisson/Laplace  (malla {M}x{N}, tol=1e-6)")
    print("  Python serial (numpy)")
    print("=" * 60)

    times = []
    for case in range(1, 5):
        t_start = time.perf_counter()
        T, X, Y, iters = solve_poisson(case, M, N)
        t_end = time.perf_counter()
        elapsed = t_end - t_start
        times.append(elapsed)

        T_ana    = analytical(case, X, Y)
        err_max  = np.max(np.abs(T - T_ana))
        err_mean = np.mean(np.abs(T - T_ana))

        print(f"\n----- Ejercicio {case}: {case_names[case]} -----")
        print(f"  Convergió en {iters} iteraciones")
        print(f"  Error máximo  : {err_max:.4e}")
        print(f"  Error promedio: {err_mean:.4e}")
        print(f"  >> Tiempo: {elapsed:.6f} s")

    print("\n" + "=" * 60)
    print("  RESUMEN DE TIEMPOS (perf_counter)")
    print("=" * 60)
    for i, t in enumerate(times, 1):
        print(f"  Ejercicio {i}: {t:.6f} s")
    print(f"  Total      : {sum(times):.6f} s")
    print("=" * 60)


if __name__ == "__main__":
    main()
