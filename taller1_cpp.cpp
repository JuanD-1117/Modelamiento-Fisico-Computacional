/*
 * Taller 1 - Modelamiento Físico Computacional 2026-1
 * Modelo de Rashevsky: Dinámica de Inconformistas
 * Implementación C++: Euler Explícito, Taylor Orden 2, Trapecio Implícito
 */

#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <string>

using namespace std;
using namespace chrono;

int main() {
    // --------------------------------------------------------
    // Parámetros del modelo
    // --------------------------------------------------------
    const double p0 = 0.01;
    const double b  = 0.02;
    const double d  = 0.015;   // solo informativo, se cancela en dp/dt
    const double r  = 0.1;
    const double k  = r * b;   // k = 0.002
    const double h  = 1.0;
    const int    steps = 50;

    // Constantes precomputadas (invariantes del bucle)
    const double hk            = h * k;
    const double factor_t2     = hk - 0.5 * h * h * k * k;
    const double hk2           = hk / 2.0;
    const double a_trap        = 1.0 - hk2;
    const double b_trap        = hk;
    const double c_trap        = 1.0 + hk2;

    const long long ITERS = 10000000LL;

    cout << string(60, '=') << "\n";
    cout << "Taller 1 - Benchmark C++\n";
    cout << "Parametros: b=" << b << " d=" << d << " r=" << r
         << " k=" << k << " h=" << h << "\n";
    cout << string(60, '=') << "\n";
    cout << fixed << setprecision(8);

    // --------------------------------------------------------
    // Solución de una trayectoria (para verificar valores)
    // --------------------------------------------------------
    cout << "\n--- Verificacion: solucion en una trayectoria ---\n";
    cout << setw(5) << "t"
         << setw(14) << "Exacta"
         << setw(14) << "Euler"
         << setw(14) << "Taylor2"
         << setw(14) << "Trapecio" << "\n";
    cout << string(61, '-') << "\n";

    double pe = p0, pt = p0, pz = p0;
    for (int i = 0; i <= steps; i++) {
        if (i % 5 == 0) {
            double t_exact = 1.0 - (1.0 - p0) * exp(-k * i * h);
            cout << setw(5) << i
                 << setw(14) << t_exact
                 << setw(14) << pe
                 << setw(14) << pt
                 << setw(14) << pz << "\n";
        }
        if (i < steps) {
            pe = pe + hk * (1.0 - pe);
            pt = pt + (1.0 - pt) * factor_t2;
            pz = (pz * a_trap + b_trap) / c_trap;
        }
    }

    // --------------------------------------------------------
    // BENCHMARK: 10^7 iteraciones
    // --------------------------------------------------------
    double p_final;
    cout << "\n" << string(60, '=') << "\n";
    cout << "BENCHMARK (" << ITERS << " iteraciones x " << steps << " pasos)\n";
    cout << string(60, '=') << "\n";
    cout << fixed << setprecision(6);

    // --- Euler Explícito ---
    auto t_start = high_resolution_clock::now();
    for (long long j = 0; j < ITERS; j++) {
        double p = p0;
        for (int i = 0; i < steps; i++)
            p = p + hk * (1.0 - p);
        if (j == ITERS - 1) p_final = p;
    }
    auto t_end = high_resolution_clock::now();
    duration<double> dur1 = t_end - t_start;
    cout << "Euler Explicito  : " << setw(10) << dur1.count()
         << " s  | p(50) = " << setprecision(8) << p_final << "\n";

    // --- Taylor Orden 2 ---
    cout << setprecision(6);
    t_start = high_resolution_clock::now();
    for (long long j = 0; j < ITERS; j++) {
        double p = p0;
        for (int i = 0; i < steps; i++)
            p = p + (1.0 - p) * factor_t2;
        if (j == ITERS - 1) p_final = p;
    }
    t_end = high_resolution_clock::now();
    duration<double> dur2 = t_end - t_start;
    cout << "Taylor Orden 2   : " << setw(10) << dur2.count()
         << " s  | p(50) = " << setprecision(8) << p_final << "\n";

    // --- Trapecio Implícito ---
    cout << setprecision(6);
    t_start = high_resolution_clock::now();
    for (long long j = 0; j < ITERS; j++) {
        double p = p0;
        for (int i = 0; i < steps; i++)
            p = (p * a_trap + b_trap) / c_trap;
        if (j == ITERS - 1) p_final = p;
    }
    t_end = high_resolution_clock::now();
    duration<double> dur3 = t_end - t_start;
    cout << "Trapecio Impl.   : " << setw(10) << dur3.count()
         << " s  | p(50) = " << setprecision(8) << p_final << "\n";

    cout << string(60, '=') << "\n";
    cout << setprecision(4);
    cout << "TIEMPOS C++ [s]: Euler=" << dur1.count()
         << "  Taylor2=" << dur2.count()
         << "  Trapecio=" << dur3.count() << "\n";

    return 0;
}
