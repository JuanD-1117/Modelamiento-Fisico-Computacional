program pozo_multipotencial
    implicit none
    
    integer, parameter :: n = 2000
    integer :: tipo_pozo ! 1: Caja, 2: V, 3: Onda (Doble pozo)
    double precision :: x(0:n), E1, E2, dE, h, x_min, x_max
    double precision :: energies(10)
    integer :: i, count, j
    
    ! --- CONFIGURACIÓN ---
    tipo_pozo = 3  ! <--- 1= Caja cuadrada 2=Pozo en V 3=Pozo doble
    x_min = -4.0d0
    x_max = 4.0d0
    h = (x_max - x_min) / dble(n)
    ! ---------------------

    ! 1. Generar y guardar el potencial
    open(10, file='potencial.csv', status='replace')
    write(10, '(A)') "x,V"
    do i = 0, n
        x(i) = x_min + i * h
        write(10, '(F12.6, ",", F12.6)') x(i), V_func(x(i), tipo_pozo)
    end do
    close(10)

    ! 2. Buscar niveles de energía
    count = 0
    dE = 0.2d0
    E1 = -15.0d0 ! Rango amplio para captar diferentes fondos
    
    print *, "Analizando pozo tipo ", tipo_pozo, "..."
    do while (E1 < 20.0d0 .and. count < 8)
        E2 = E1 + dE
        if (solve_numerov(E1, n, h, x, tipo_pozo) * solve_numerov(E2, n, h, x, tipo_pozo) < 0.d0) then
            count = count + 1
            energies(count) = biseccion(E1, E2, n, h, x, tipo_pozo)
            print *, "Nivel ", count, " encontrado en E =", energies(count)
        end if
        E1 = E2
    end do

    ! 3. Guardar niveles
    open(20, file='niveles.csv', status='replace')
    write(20, '(A)') "E"
    do j = 1, count
        write(20, '(F12.6)') energies(j)
    end do
    close(20)

contains

    function V_func(pos, selector)
        double precision :: pos, V_func
        integer :: selector
        select case (selector)
        case(1) ! Caja Cuadrada
            if (abs(pos) < 1.5d0) then; V_func = 0.d0; else; V_func = 25.d0; end if
        case(2) ! Pozo en V (Lineal)
            V_func = 4.0d0 * abs(pos)
        case(3) ! Pozo tipo "Onda" (Doble pozo)
            V_func = 0.5d0*pos**4 - 5.0d0*pos**2
        end select
    end function

    function solve_numerov(E, n, h, x, selector) result(psif)
        integer :: n, j, selector
        double precision :: E, h, x(0:n), psif, p(0:n), k2(0:n), h2_12
        h2_12 = (h**2)/12.d0
        do j = 0, n
            k2(j) = 2.0d0 * (E - V_func(x(j), selector))
        end do
        p(0) = 0.d0; p(1) = 1e-5
        do j = 1, n-1
            p(j+1) = (2.d0*(1.d0-5.d0*h2_12*k2(j))*p(j) - (1.d0+h2_12*k2(j-1))*p(j-1))/(1.d0+h2_12*k2(j+1))
            if (abs(p(j+1)) > 1e10) p(j+1) = 1e10 ! Evitar overflow
        end do
        psif = p(n)
    end function

    function biseccion(E1, E2, n, h, x, selector) result(Em)
        double precision :: E1, E2, Em
        integer :: n, k, selector
        double precision :: h, x(0:n)
        do k = 1, 40
            Em = (E1 + E2) / 2.d0
            if (solve_numerov(E1, n, h, x, selector) * solve_numerov(Em, n, h, x, selector) < 0.d0) then
                E2 = Em
            else
                E1 = Em
            end if
        end do
    end function
end program