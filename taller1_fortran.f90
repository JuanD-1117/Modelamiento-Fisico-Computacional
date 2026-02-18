! =============================================================
! Taller 1 - Modelamiento Físico Computacional 2026-1
! Modelo de Rashevsky: Dinámica de Inconformistas
! Implementación Fortran 90: Euler Explícito, Taylor Orden 2,
!                            Trapecio Implícito
! =============================================================
program taller1
    implicit none

    ! ---- Parámetros del modelo ----
    real(8), parameter :: p0   = 0.01d0
    real(8), parameter :: b    = 0.02d0
    real(8), parameter :: d    = 0.015d0   ! informativo, se cancela
    real(8), parameter :: r    = 0.1d0
    real(8), parameter :: k    = r * b     ! k = 0.002
    real(8), parameter :: h    = 1.0d0
    integer, parameter :: steps = 50

    ! ---- Constantes precomputadas ----
    real(8), parameter :: hk         = h * k
    real(8), parameter :: factor_t2  = hk - 0.5d0 * h * h * k * k
    real(8), parameter :: hk2        = hk / 2.0d0
    real(8), parameter :: a_trap     = 1.0d0 - hk2
    real(8), parameter :: b_trap     = hk
    real(8), parameter :: c_trap     = 1.0d0 + hk2

    integer(8), parameter :: ITERS = 10000000_8

    integer          :: i
    integer(8)       :: j
    real(8)          :: pe, pt, pz, t_exact
    real(8)          :: start_time, end_time

    ! ===========================================================
    write(*,'(A)') repeat("=", 60)
    write(*,'(A)') "Taller 1 - Benchmark Fortran"
    write(*,'(A,F6.4,A,F6.4,A,F4.2,A,F6.4)') &
        "Parametros: b=",b,"  d=",d,"  r=",r,"  k=",k
    write(*,'(A)') repeat("=", 60)

    ! ===========================================================
    ! Verificación: solución en una trayectoria
    ! ===========================================================
    write(*,'(/,A)') "--- Verificacion: solucion en una trayectoria ---"
    write(*,'(A5,4A14)') "t", "Exacta", "Euler", "Taylor2", "Trapecio"
    write(*,'(A)') repeat("-", 61)

    pe = p0;  pt = p0;  pz = p0
    do i = 0, steps
        if (mod(i,5) == 0) then
            t_exact = 1.0d0 - (1.0d0 - p0) * exp(-k * dble(i) * h)
            write(*,'(I5, 4F14.8)') i, t_exact, pe, pt, pz
        end if
        if (i < steps) then
            pe = pe + hk * (1.0d0 - pe)
            pt = pt + (1.0d0 - pt) * factor_t2
            pz = (pz * a_trap + b_trap) / c_trap
        end if
    end do

    ! ===========================================================
    ! BENCHMARK: 10^7 iteraciones
    ! ===========================================================
    write(*,'(/,A)') repeat("=", 60)
    write(*,'(A,I0,A,I0,A)') "BENCHMARK (", ITERS, " iteraciones x ", steps, " pasos)"
    write(*,'(A)') repeat("=", 60)

    ! --- Euler Explícito ---
    call cpu_time(start_time)
    do j = 1_8, ITERS
        pe = p0
        do i = 1, steps
            pe = pe + hk * (1.0d0 - pe)
        end do
    end do
    call cpu_time(end_time)
    write(*,'(A,F10.4,A,F12.8)') &
        "Euler Explicito  : ", end_time-start_time, " s  | p(50) = ", pe

    ! --- Taylor Orden 2 ---
    call cpu_time(start_time)
    do j = 1_8, ITERS
        pt = p0
        do i = 1, steps
            pt = pt + (1.0d0 - pt) * factor_t2
        end do
    end do
    call cpu_time(end_time)
    write(*,'(A,F10.4,A,F12.8)') &
        "Taylor Orden 2   : ", end_time-start_time, " s  | p(50) = ", pt

    ! --- Trapecio Implícito ---
    call cpu_time(start_time)
    do j = 1_8, ITERS
        pz = p0
        do i = 1, steps
            pz = (pz * a_trap + b_trap) / c_trap
        end do
    end do
    call cpu_time(end_time)
    write(*,'(A,F10.4,A,F12.8)') &
        "Trapecio Impl.   : ", end_time-start_time, " s  | p(50) = ", pz

    write(*,'(A)') repeat("=", 60)

end program taller1
