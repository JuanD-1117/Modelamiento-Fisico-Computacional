!=============================================================================!
!  SOLUCIÓN NUMÉRICA: ECUACIONES DE POISSON Y LAPLACE                        !
!  Método de Diferencias Finitas — Iteración de Jacobi                        !
!                                                                             !
!  Compilar: gfortran -O2 -o poisson_laplace poisson_laplace.f90             !
!  Ejecutar: ./poisson_laplace [M] [N]   (default M=N=50)                    !
!=============================================================================!

program poisson_laplace
  implicit none
  integer      :: M, N
  real(8)      :: t1, t2, times(4)
  character(20):: arg

  ! Leer tamaño de malla desde argumentos (default 50)
  M = 50;  N = 50
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg);  read(arg,*) M
  end if
  if (command_argument_count() >= 2) then
    call get_command_argument(2, arg);  read(arg,*) N
  end if

  write(*,*) "============================================================"
  write(*,'(A,I0,A,I0,A)') &
    "  Jacobi FDM — Poisson/Laplace  (malla ", M, "x", N, ", tol=1e-6)"
  write(*,*) "  Fortran serial"
  write(*,*) "============================================================"

  write(*,*)
  write(*,*) "----- Ejercicio 1: Poisson  f=(x²+y²)e^(xy)  [0,2]x[0,1] -----"
  call cpu_time(t1);  call ejercicio1(M, N);  call cpu_time(t2)
  times(1) = t2 - t1
  write(*,'(A,F12.6,A)') "  >> Tiempo: ", times(1), " s"

  write(*,*)
  write(*,*) "----- Ejercicio 2: Laplace  f=0               [1,2]x[0,1] -----"
  call cpu_time(t1);  call ejercicio2(M, N);  call cpu_time(t2)
  times(2) = t2 - t1
  write(*,'(A,F12.6,A)') "  >> Tiempo: ", times(2), " s"

  write(*,*)
  write(*,*) "----- Ejercicio 3: Poisson  f=4               [1,2]x[0,2] -----"
  call cpu_time(t1);  call ejercicio3(M, N);  call cpu_time(t2)
  times(3) = t2 - t1
  write(*,'(A,F12.6,A)') "  >> Tiempo: ", times(3), " s"

  write(*,*)
  write(*,*) "----- Ejercicio 4: Poisson  f=x/y+y/x        [1,2]x[1,2] -----"
  call cpu_time(t1);  call ejercicio4(M, N);  call cpu_time(t2)
  times(4) = t2 - t1
  write(*,'(A,F12.6,A)') "  >> Tiempo: ", times(4), " s"

  write(*,*)
  write(*,*) "============================================================"
  write(*,*) "  RESUMEN DE TIEMPOS (CPU time)"
  write(*,*) "============================================================"
  write(*,'(A,F12.6,A)') "  Ejercicio 1: ", times(1), " s"
  write(*,'(A,F12.6,A)') "  Ejercicio 2: ", times(2), " s"
  write(*,'(A,F12.6,A)') "  Ejercicio 3: ", times(3), " s"
  write(*,'(A,F12.6,A)') "  Ejercicio 4: ", times(4), " s"
  write(*,'(A,F12.6,A)') "  Total      : ", sum(times), " s"
  write(*,*) "============================================================"

end program poisson_laplace


!=============================================================================!
! SUBRUTINA AUXILIAR: Imprime estadísticas de error                           !
!=============================================================================!
subroutine print_stats(T_num, T_ana, M, N, label)
  implicit none
  integer,  intent(in) :: M, N
  real(8),  intent(in) :: T_num(0:N,0:M), T_ana(0:N,0:M)
  character(len=*), intent(in) :: label
  real(8) :: err_max, err_mean, err
  integer :: i, j

  err_max  = 0.0d0
  err_mean = 0.0d0
  do i = 0, N
    do j = 0, M
      err = abs(T_num(i,j) - T_ana(i,j))
      if (err > err_max) err_max = err
      err_mean = err_mean + err
    end do
  end do
  err_mean = err_mean / real((M+1)*(N+1), 8)

  write(*,'(A,A)')       "  Resultado: ", label
  write(*,'(A,ES12.4)')  "  Error máximo    : ", err_max
  write(*,'(A,ES12.4)')  "  Error promedio  : ", err_mean
end subroutine print_stats


!=============================================================================!
! EJERCICIO 1                                                                 !
! ∇²V = (x²+y²)·e^(xy)   en [0,2]x[0,1]                                    !
! Analítica: V(x,y) = e^(xy)                                                 !
!=============================================================================!
subroutine ejercicio1(M, N)
  implicit none
  integer, intent(in) :: M, N
  real(8), parameter :: x0=0.0d0, xf=2.0d0, y0=0.0d0, yf=1.0d0
  real(8), parameter :: tol = 1.0d-6
  real(8) :: h, k, h2, k2, denom
  real(8), allocatable :: X(:,:), Y(:,:), T(:,:), T_new(:,:), T_ana(:,:), source(:,:)
  integer :: i, j, it
  real(8) :: delta

  allocate(X(0:N,0:M), Y(0:N,0:M), T(0:N,0:M), T_new(0:N,0:M), &
           T_ana(0:N,0:M), source(0:N,0:M))

  h = (xf - x0) / M;  k = (yf - y0) / N
  h2 = h*h;  k2 = k*k;  denom = 2.0d0*(h2 + k2)

  do i = 0, N
    do j = 0, M
      X(i,j) = x0 + j*h;  Y(i,j) = y0 + i*k
    end do
  end do

  T = 0.0d0
  do i = 0, N
    T(i, 0) = 1.0d0;  T(i, M) = exp(2.0d0 * Y(i,M))
  end do
  do j = 0, M
    T(0, j) = 1.0d0;  T(N, j) = exp(X(N,j))
  end do

  do i = 0, N
    do j = 0, M
      source(i,j) = (X(i,j)**2 + Y(i,j)**2) * exp(X(i,j)*Y(i,j))
    end do
  end do

  T_new = T
  do it = 1, 400000
    delta = 0.0d0
    do i = 1, N-1
      do j = 1, M-1
        T_new(i,j) = ( (T(i+1,j) + T(i-1,j))*k2 &
                     + (T(i,j+1) + T(i,j-1))*h2  &
                     - source(i,j)*h2*k2          ) / denom
        delta = max(delta, abs(T_new(i,j) - T(i,j)))
      end do
    end do
    T = T_new
    if (delta < tol) then
      write(*,'(A,I6,A)') "  Convergió en ", it, " iteraciones."
      exit
    end if
  end do

  do i = 0, N
    do j = 0, M
      T_ana(i,j) = exp(X(i,j)*Y(i,j))
    end do
  end do
  call print_stats(T, T_ana, M, N, "Ejercicio 1 — Poisson")
  deallocate(X, Y, T, T_new, T_ana, source)
end subroutine ejercicio1


!=============================================================================!
! EJERCICIO 2                                                                 !
! ∇²V = 0   (Laplace)  en [1,2]x[0,1]                                       !
! Analítica: V(x,y) = ln(x²+y²)                                              !
!=============================================================================!
subroutine ejercicio2(M, N)
  implicit none
  integer, intent(in) :: M, N
  real(8), parameter :: x0=1.0d0, xf=2.0d0, y0=0.0d0, yf=1.0d0
  real(8), parameter :: tol = 1.0d-6
  real(8) :: h, k, h2, k2, denom
  real(8), allocatable :: X(:,:), Y(:,:), T(:,:), T_new(:,:), T_ana(:,:)
  integer :: i, j, it
  real(8) :: delta

  allocate(X(0:N,0:M), Y(0:N,0:M), T(0:N,0:M), T_new(0:N,0:M), T_ana(0:N,0:M))

  h = (xf - x0) / M;  k = (yf - y0) / N
  h2 = h*h;  k2 = k*k;  denom = 2.0d0*(h2 + k2)

  do i = 0, N
    do j = 0, M
      X(i,j) = x0 + j*h;  Y(i,j) = y0 + i*k
    end do
  end do

  T = 0.0d0
  do i = 0, N
    T(i, 0) = log(Y(i,0)**2 + 1.0d0)
    T(i, M) = log(Y(i,M)**2 + 4.0d0)
  end do
  do j = 0, M
    T(0, j) = 2.0d0*log(X(0,j))
    T(N, j) = log(X(N,j)**2 + 1.0d0)
  end do

  T_new = T
  do it = 1, 400000
    delta = 0.0d0
    do i = 1, N-1
      do j = 1, M-1
        T_new(i,j) = ( (T(i+1,j) + T(i-1,j))*k2 &
                     + (T(i,j+1) + T(i,j-1))*h2  ) / denom
        delta = max(delta, abs(T_new(i,j) - T(i,j)))
      end do
    end do
    T = T_new
    if (delta < tol) then
      write(*,'(A,I6,A)') "  Convergió en ", it, " iteraciones."
      exit
    end if
  end do

  do i = 0, N
    do j = 0, M
      T_ana(i,j) = log(X(i,j)**2 + Y(i,j)**2)
    end do
  end do
  call print_stats(T, T_ana, M, N, "Ejercicio 2 — Laplace")
  deallocate(X, Y, T, T_new, T_ana)
end subroutine ejercicio2


!=============================================================================!
! EJERCICIO 3                                                                 !
! ∇²V = 4   en [1,2]x[0,2]                                                  !
! Analítica: V(x,y) = (x-y)²                                                 !
!=============================================================================!
subroutine ejercicio3(M, N)
  implicit none
  integer, intent(in) :: M, N
  real(8), parameter :: x0=1.0d0, xf=2.0d0, y0=0.0d0, yf=2.0d0
  real(8), parameter :: tol = 1.0d-6
  real(8) :: h, k, h2, k2, denom
  real(8), allocatable :: X(:,:), Y(:,:), T(:,:), T_new(:,:), T_ana(:,:), source(:,:)
  integer :: i, j, it
  real(8) :: delta

  allocate(X(0:N,0:M), Y(0:N,0:M), T(0:N,0:M), T_new(0:N,0:M), &
           T_ana(0:N,0:M), source(0:N,0:M))

  h = (xf - x0) / M;  k = (yf - y0) / N
  h2 = h*h;  k2 = k*k;  denom = 2.0d0*(h2 + k2)

  do i = 0, N
    do j = 0, M
      X(i,j) = x0 + j*h;  Y(i,j) = y0 + i*k
    end do
  end do

  T = 0.0d0
  do i = 0, N
    T(i, 0) = (x0 - Y(i,0))**2
    T(i, M) = (xf - Y(i,M))**2
  end do
  do j = 0, M
    T(0, j) = (X(0,j) - y0)**2
    T(N, j) = (X(N,j) - yf)**2
  end do

  source = 4.0d0

  T_new = T
  do it = 1, 400000
    delta = 0.0d0
    do i = 1, N-1
      do j = 1, M-1
        T_new(i,j) = ( (T(i+1,j) + T(i-1,j))*k2 &
                     + (T(i,j+1) + T(i,j-1))*h2  &
                     - source(i,j)*h2*k2          ) / denom
        delta = max(delta, abs(T_new(i,j) - T(i,j)))
      end do
    end do
    T = T_new
    if (delta < tol) then
      write(*,'(A,I6,A)') "  Convergió en ", it, " iteraciones."
      exit
    end if
  end do

  do i = 0, N
    do j = 0, M
      T_ana(i,j) = (X(i,j) - Y(i,j))**2
    end do
  end do
  call print_stats(T, T_ana, M, N, "Ejercicio 3 — Poisson f=4")
  deallocate(X, Y, T, T_new, T_ana, source)
end subroutine ejercicio3


!=============================================================================!
! EJERCICIO 4                                                                 !
! ∇²V = x/y + y/x   en [1,2]x[1,2]                                          !
! Analítica: V(x,y) = xy·ln(xy)                                              !
!=============================================================================!
subroutine ejercicio4(M, N)
  implicit none
  integer, intent(in) :: M, N
  real(8), parameter :: x0=1.0d0, xf=2.0d0, y0=1.0d0, yf=2.0d0
  real(8), parameter :: tol = 1.0d-6
  real(8) :: h, k, h2, k2, denom
  real(8), allocatable :: X(:,:), Y(:,:), T(:,:), T_new(:,:), T_ana(:,:), source(:,:)
  integer :: i, j, it
  real(8) :: delta

  allocate(X(0:N,0:M), Y(0:N,0:M), T(0:N,0:M), T_new(0:N,0:M), &
           T_ana(0:N,0:M), source(0:N,0:M))

  h = (xf - x0) / M;  k = (yf - y0) / N
  h2 = h*h;  k2 = k*k;  denom = 2.0d0*(h2 + k2)

  do i = 0, N
    do j = 0, M
      X(i,j) = x0 + j*h;  Y(i,j) = y0 + i*k
    end do
  end do

  T = 0.0d0
  do i = 0, N
    T(i, 0) = Y(i,0) * log(Y(i,0))
    T(i, M) = 2.0d0*Y(i,M)*log(2.0d0*Y(i,M))
  end do
  do j = 0, M
    T(0, j) = X(0,j) * log(X(0,j))
    T(N, j) = X(N,j) * log(4.0d0*X(N,j)**2)
  end do

  do i = 0, N
    do j = 0, M
      source(i,j) = X(i,j)/Y(i,j) + Y(i,j)/X(i,j)
    end do
  end do

  T_new = T
  do it = 1, 400000
    delta = 0.0d0
    do i = 1, N-1
      do j = 1, M-1
        T_new(i,j) = ( (T(i+1,j) + T(i-1,j))*k2 &
                     + (T(i,j+1) + T(i,j-1))*h2  &
                     - source(i,j)*h2*k2          ) / denom
        delta = max(delta, abs(T_new(i,j) - T(i,j)))
      end do
    end do
    T = T_new
    if (delta < tol) then
      write(*,'(A,I6,A)') "  Convergió en ", it, " iteraciones."
      exit
    end if
  end do

  do i = 0, N
    do j = 0, M
      T_ana(i,j) = X(i,j)*Y(i,j)*log(X(i,j)*Y(i,j))
    end do
  end do
  call print_stats(T, T_ana, M, N, "Ejercicio 4 — Poisson f=x/y+y/x")
  deallocate(X, Y, T, T_new, T_ana, source)
end subroutine ejercicio4
