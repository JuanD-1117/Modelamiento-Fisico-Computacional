// --- ARCHIVO: poisson_serial.cpp ---
// Solución numérica: Ecuaciones de Poisson y Laplace
// Método de Diferencias Finitas — Iteración de Jacobi
// Compilar: g++ -O2 -o poisson_serial poisson_serial.cpp
// Ejecutar: ./poisson_serial [M] [N]   (default M=N=50)

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <chrono>

// --- VARIABLES GLOBALES DE CONFIGURACIÓN ---
int case_selector = 1;
double x_ini, x_fin, y_ini, y_fin;
static constexpr double EPS = 1e-12;
static constexpr int MAX_ITERATIONS = 400000;
static constexpr double TOLERANCE = 1e-6;

// --- FUNCIONES DEL PROBLEMA FÍSICO ---
double source_term(double x, double y) {
    switch (case_selector) {
        case 1: return (x*x + y*y) * std::exp(x * y);
        case 2: return 0.0;
        case 3: return 4.0;
        case 4: return x / y + y / x;
        default: return 0.0;
    }
}

double boundary_condition(double x, double y) {
    switch (case_selector) {
        case 1:
            if (std::abs(x - x_ini) < EPS) return 1.0;
            if (std::abs(x - x_fin) < EPS) return std::exp(2.0 * y);
            if (std::abs(y - y_ini) < EPS) return 1.0;
            if (std::abs(y - y_fin) < EPS) return std::exp(x);
            break;
        case 2:
            if (std::abs(x - x_ini) < EPS) return std::log(y*y + 1.0);
            if (std::abs(x - x_fin) < EPS) return std::log(y*y + 4.0);
            if (std::abs(y - y_ini) < EPS) return 2.0 * std::log(x);
            if (std::abs(y - y_fin) < EPS) return std::log(x*x + 1.0);
            break;
        case 3:
            if (std::abs(x - x_ini) < EPS) return (1.0 - y) * (1.0 - y);
            if (std::abs(x - x_fin) < EPS) return (2.0 - y) * (2.0 - y);
            if (std::abs(y - y_ini) < EPS) return x * x;
            if (std::abs(y - y_fin) < EPS) return (x - 2.0) * (x - 2.0);
            break;
        case 4:
            if (std::abs(x - x_ini) < EPS) return y * std::log(y);
            if (std::abs(x - x_fin) < EPS) return 2.0 * y * std::log(2.0 * y);
            if (std::abs(y - y_ini) < EPS) return x * std::log(x);
            if (std::abs(y - y_fin) < EPS) return 2.0 * x * std::log(2.0 * x);
            break;
    }
    return 0.0;
}

// --- FUNCIONES DE CÁLCULO ---
void initialize_grid(int M, int N, std::vector<std::vector<double>> &V, double h, double k) {
    V.assign(M + 1, std::vector<double>(N + 1, 0.0));
    for (int i = 0; i <= M; ++i) {
        for (int j = 0; j <= N; ++j) {
            if (i == 0 || i == M || j == 0 || j == N) {
                double x = x_ini + i * h;
                double y = y_ini + j * k;
                V[i][j] = boundary_condition(x, y);
            }
        }
    }
}

void solve_poisson_fdm(int M, int N, std::vector<std::vector<double>> &V, double h, double k) {
    std::vector<std::vector<double>> V_old = V;
    double delta = TOLERANCE + 1.0;
    int iterations = 0;
    const double h2k2 = h*h*k*k;
    const double denom = 2.0 * (h*h + k*k);

    while (delta > TOLERANCE && iterations < MAX_ITERATIONS) {
        delta = 0.0;
        V_old = V;

        for (int i = 1; i < M; ++i) {
            for (int j = 1; j < N; ++j) {
                double x = x_ini + i * h;
                double y = y_ini + j * k;
                double f = source_term(x, y);

                double numer = (V_old[i+1][j] + V_old[i-1][j]) * (k*k)
                             + (V_old[i][j+1] + V_old[i][j-1]) * (h*h) - f * h2k2;
                
                V[i][j] = numer / denom;

                double diff = std::abs(V[i][j] - V_old[i][j]);
                if (diff > delta) {
                    delta = diff;
                }
            }
        }
        iterations++;
    }
}

// --- FUNCIÓN MAIN ---
int main(int argc, char *argv[]) {
    // Tamaño de malla opcional por argumento (default 50x50)
    const int M = (argc >= 2) ? std::atoi(argv[1]) : 50;
    const int N = (argc >= 3) ? std::atoi(argv[2]) : 50;

    struct CaseConfig {
        int id;
        double x0, xf, y0, yf;
        const char* name;
    };
    CaseConfig cases[] = {
        {1, 0.0, 2.0, 0.0, 1.0, "Poisson f=(x^2+y^2)e^(xy)  [0,2]x[0,1]"},
        {2, 1.0, 2.0, 0.0, 1.0, "Laplace f=0                 [1,2]x[0,1]"},
        {3, 1.0, 2.0, 0.0, 2.0, "Poisson f=4                 [1,2]x[0,2]"},
        {4, 1.0, 2.0, 1.0, 2.0, "Poisson f=x/y+y/x           [1,2]x[1,2]"},
    };

    std::cout << "============================================================\n";
    std::cout << "  Jacobi FDM — Poisson/Laplace  (malla " << M << "x" << N
              << ", tol=1e-6)\n";
    std::cout << "  C++ serial\n";
    std::cout << "============================================================\n";

    double total_time = 0.0;
    double times[4];

    for (int c = 0; c < 4; ++c) {
        case_selector = cases[c].id;
        x_ini = cases[c].x0;  x_fin = cases[c].xf;
        y_ini = cases[c].y0;  y_fin = cases[c].yf;

        double h = (x_fin - x_ini) / M;
        double k = (y_fin - y_ini) / N;
        std::vector<std::vector<double>> V;

        auto t_start = std::chrono::high_resolution_clock::now();
        initialize_grid(M, N, V, h, k);
        solve_poisson_fdm(M, N, V, h, k);
        auto t_end = std::chrono::high_resolution_clock::now();

        double elapsed = std::chrono::duration<double>(t_end - t_start).count();
        times[c] = elapsed;
        total_time += elapsed;

        std::cout << "\n----- Ejercicio " << cases[c].id << ": "
                  << cases[c].name << " -----\n";
        std::cout << "  >> Tiempo: " << elapsed << " s\n";
    }

    std::cout << "\n============================================================\n";
    std::cout << "  RESUMEN DE TIEMPOS (wall time)\n";
    std::cout << "============================================================\n";
    for (int c = 0; c < 4; ++c) {
        std::cout << "  Ejercicio " << (c+1) << ": " << times[c] << " s\n";
    }
    std::cout << "  Total      : " << total_time << " s\n";
    std::cout << "============================================================\n";

    return 0;
}