PROGRAM ACUSTICA
C===========================================================
C===========================================================
C         PROPAGACION ACUSTICA EN EL OCEANO
C===========================================================
C===========================================================
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER MODOS*10, TL*10, RAYOS*10
      PARAMETER (N=4, NPT=1000, NBMAX=1000, N1=3500)
C
      COMPLEX*16 SUM,CI
C
      COMMON /DET/W,HZ,NN
C
      DIMENSION Y(N), DYDX(N), YOUT(N),
     &          XB1(NBMAX), XB2(NBMAX), PSIF(NPT),
     &          PSI(100,NPT), C(NPT),
     &          DIAG(NPT), SUPERD(NPT), SUBD(NPT),
     &          VK(NPT), VKM1(NPT)
      REAL*8 MODO(NBMAX)
      EXTERNAL DETERM
      PI=4.*DATAN(1.D0)
      CI=(0.D0,1.D0)
C===========================================================
      WRITE (6,*) 'PROFUNDIDAD DEL FONDO (MTS.) ? '
      READ  (5,*) D
      WRITE (6,*) 'PROFUNDIDAD DE LA FUENTE (MTS.)? '
      READ  (5,*) Z0
      WRITE (6,*) 'CUAL ES LA FRECUENCIA DE LA SENAL? '
      READ  (5,*) W
      WRITE (6,*) 'CUAL ES EL CONO DE EMISION DE LOS RAYOS '
      WRITE (6,*) 'TETA MINIMO EN GRADOS '
      READ  (5,*) NTMIN
      WRITE (6,*) 'TETA MAXIMO EN GRADOS '
      READ  (5,*) NTMAX
      WRITE (6,*) 'HASTA QUE RANGO DESEA DIBUJAR (KM.)? '
      READ  (5,*) PRR
      PRR=1000.*PRR
      WRITE (6,*) 'NOMBRE DEL ARCHIVO PARA RAYOS '
      READ  (5,*) RAYOS
      OPEN  (97,FILE=RAYOS)
      WRITE (6,*) 'ARCHIVO PARA LOS DATOS DE LOS MODOS '
      READ  (5,*) MODOS
      OPEN  (98,FILE=MODOS)
      WRITE (6,*)'NOMBRE DEL ARCHIVO PARA TL '
      READ  (5,*) TL
      OPEN  (99,FILE=TL)
C===========================================================
C     TRAZADO DE RAYOS
      HS=.5                   ! PASO PARA S
                              ! EN LA INTEGRACION
      DO NAN=NTMIN,NTMAX
        PHI=FLOAT(NAN)
        S=0.D0                ! ARGUMENTO INICIAL
C                             ! CONDICIONES INICIALES
        CS=VSON(Z0)
        Y(1)=0.D0
        Y(2)=Z0
        Y(3)=DCOS(-PHI*PI/180.D0)/CS
        Y(4)=DSIN(-PHI*PI/180.D0)/CS
        CALL DERIVS(Y,DYDX)
C                             ! CICLO SOBRE EL INTERVALO S
        DO L=1,5*INT(PRR)
          S=S+HS
          CALL RK4(Y,DYDX,N,HS,YOUT)
          DO I=1,N
C                             ! CONDICION REBOTE EN LA
C                             ! SUPERFICIE O FONDO
            IF (YOUT(2).LE.1D-003.OR.YOUT(2).GE.D) THEN
              CS=VSON(Y(2))
              Y(3)=DCOS(-DATAN(YOUT(4)/YOUT(3)))/CS
              Y(4)=DSIN(-DATAN(YOUT(4)/YOUT(3)))/CS
            ELSE
              Y(I)=YOUT(I)
            END IF
          END DO
          IF(MOD(L,100).EQ.0)
     &      WRITE (97,'(I4,2F16.6)') NAN, YOUT(1),YOUT(2)
          IF(YOUT(1).GE.PRR) GOTO 100
          CALL DERIVS(Y,DYDX)
        END DO
  100   WRITE (97,'()')
      END DO
C===========================================================
C     SOLUCION DE LA ECUACION DE HELMHOLTZ MEDIANTE
C     EL METODO DE MODOS NORMALES
C     PASO
      HZ=D/FLOAT(NPT)
C                             ! RANGOS DE BUSQUEDA
C                             ! DE LOS MODOS
      X1=W/VSON(1300.D0)
      X2=W/VSON(5000.D0)
      NB=NBMAX
      NN=NPT
      CALL ZBRAK(DETERM,X1,X2,N1,XB1,XB2,NB)
      JJ=0
      DO I=1,NB
        TOL=(1.0D-12)*(XB1(I)+XB2(I))/2.D0
        RAIZ=ZBRENT(DETERM,XB1(I),XB2(I),TOL)
        JJ=JJ+1
        MODO(JJ)=RAIZ
        WRITE (98,*) JJ,MODO(JJ), W/MODO(JJ)
      END DO
C===========================================================
C     CALCULO DE LAS FUNCIONES NORMALIZADAS
      DO I=1,NPT
        C(I)=VSON(FLOAT(I)*HZ)
      END DO
      DO J=1,JJ
        VP=MODO(J)
        DO K=1,NPT
          DIAG(K)=-2.D0+HZ*HZ*(W/C(K))*(W/C(K))-HZ*HZ*VP*VP
        END DO
        DO K=1,NPT-1
          SUPERD(K)=1.D0
        END DO
        DO K=2,NPT
          SUBD(K)=1.D0
        END DO
        SUBD(NPT)=2.D0
C                             ! VECTOR DE ARRANQUE
        DO K=1,NPT
          VK(K)=0.D0
        END DO
        VK(1)=1.
C                             ! BUSQUEDA DE LA FUNCION
C                             ! NORMALIZADA
        DO KK=1, 10
          CALL TRIDAG(SUBD,DIAG,SUPERD,VK,VKM1,NPT)
          CALL NORMA(NPT,VKM1,VK)
        END DO
        PSIF(J)=VK(INT(Z0/HZ)+1)
        DO IP=1, INT(D/100.)
          PROF_DET=(IP-1)*100
          PSI(IP,J)=VK(INT(PROF_DET/HZ)+1)
        END DO
      END DO
C===========================================================
C     DETERMINACION DE LA PERDIDA POR TRANSMISION (EC 20)
      DO IP=1,INT(D/100)
        DO I=1,INT(PRR)
          R=FLOAT(I)
          SUM=0.D0
          DO J=1,JJ
            SUM=SUM+PSIF(J)*PSI(IP,J)*(CDEXP(CI*MODO(J)*R))/
     &                                  DSQRT(MODO(J))
          END DO
          PERD=-20.D0*DLOG10(CDABS(SUM*(DSQRT(2.D0*PI/R))))
          IF(MOD(I,100).EQ.0) WRITE (99,*) R/1000.,-100.*IP,PERD
        END DO
        WRITE (99, '()')
      END DO
C===========================================================
      CLOSE (97)
      CLOSE (98)
      CLOSE (99)
      END
C===========================================================
      SUBROUTINE DERIVS(Y,DYDX)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION Y(*),DYDX(*)
C     ECUACIONES DIFERENCIALES PARA EL TRAZADO
C     DE RAYOS (EC 7,8)
      CS=VSON(Y(2))
      DYDX(1)=CS*Y(3)
      DYDX(2)=CS*Y(4)
      DYDX(3)=0.D0
      DYDX(4)=-VDER(Y(2))/(CS*CS)
      RETURN
      END
C===========================================================
      FUNCTION VSON(Z)
      IMPLICIT REAL*8(A-H,O-Z)
C     PERFIL DE VELOCIDAD DE MUNK (EC 33)
      PARAMETER (V0=1500.D0, EPS=0.00737D0, ZM=1300.D0)
      G=2.D0*(Z-ZM)/ZM
      VSON=V0*(1.D0+EPS*(G-1.D0+DEXP(-G)))
      RETURN
      END
C===========================================================
      FUNCTION VDER(Z)
      IMPLICIT REAL*8(A-H,O-Z)
C     DERIVADA DEL PERFIL DE VELOCIDAD DE MUNK
      PARAMETER (V0=1500.D0, EPS=0.00737D0, ZM=1300.D0)
C     DERIVADA DE LA VELOCIDAD DEL SONIDO
      G=2.D0*(Z-ZM)/ZM
      VDER=(2.D0*V0/ZM)*EPS*(1.D0-DEXP(-G))
      RETURN
      END
C===========================================================
      FUNCTION DETERM(X)
      IMPLICIT REAL*8(A-H,O-Z)
C     CALCULO DEL DETERMINANTE (EC 30)
      DIMENSION P(0:2000), C(2000)
      COMMON /DET/ W, H, N
C     CONSTRUCCION DEL VECTOR VELOCIDAD
      DO I=1,N
        C(I)=VSON(H*FLOAT(I))
      END DO
      P(0)=0.D0
      P(1)=1.D0
      DO I=2,N
        P(I)=(-2.D0+H*H*(W/C(I-1))*(W/C(I-1))-H*H*X*X)*
     &       P(I-1)-P(I-2)
      END DO
      DETERM=((-2.D0+H*H*(W/C(N))*(W/C(N))-H*H*X*X)*
     &       P(N)-2.D0*P(N-1))
      RETURN
      END
C===========================================================
      SUBROUTINE NORMA(N,VN,V)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VN(N),V(N)
C     NORMALIZACION DEL VECTOR
      SUM=0.
      DO I=1,N
        SUM=SUM+VN(I)*VN(I)
      END DO
      DO I=1,N
        V(I)=VN(I)/DSQRT(SUM)
      END DO
      RETURN
      END
C===========================================================
C     SUBRUTINA RK4 - INTEGRADOR RUNGE-KUTTA ORDEN 4
C     Numerical Recipes in Fortran, Cap. 16
C     Integra un sistema de N ODEs un paso H
C     Y(N)    : valores actuales
C     DYDX(N) : derivadas en el punto actual
C     N       : numero de ecuaciones
C     H       : paso de integracion
C     YOUT(N) : solucion al paso siguiente
C===========================================================
      SUBROUTINE RK4(Y,DYDX,N,H,YOUT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(N),DYDX(N),YOUT(N)
      DIMENSION YT(4),DYT(4),DYM(4)
      HH=H*0.5D0
      H6=H/6.D0
C     PRIMER PASO
      DO I=1,N
        YT(I)=Y(I)+HH*DYDX(I)
      END DO
      CALL DERIVS(YT,DYT)
C     SEGUNDO PASO
      DO I=1,N
        YT(I)=Y(I)+HH*DYT(I)
      END DO
      CALL DERIVS(YT,DYM)
C     TERCER PASO
      DO I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
      END DO
      CALL DERIVS(YT,DYT)
C     COMBINACION DE LOS CUATRO PASOS
      DO I=1,N
        YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.D0*DYM(I))
      END DO
      RETURN
      END
C===========================================================
C     SUBRUTINA ZBRAK - BUSQUEDA DE BRACKETS DE RAICES
C     Numerical Recipes in Fortran, Cap. 9
C     Divide el intervalo [X1,X2] en N subintervalos y
C     detecta cambios de signo de FUNC (posibles raices)
C     XB1, XB2 : extremos de los brackets encontrados
C     NB (entrada): maximo de brackets buscados
C     NB (salida) : numero de brackets encontrados
C===========================================================
      SUBROUTINE ZBRAK(FUNC,X1,X2,N,XB1,XB2,NB)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XB1(NB),XB2(NB)
      NBB=NB
      NB=0
      DX=(X2-X1)/FLOAT(N)
      X=X1
      FP=FUNC(X)
      DO I=1,N
        X=X+DX
        FC=FUNC(X)
        IF(FC*FP.LT.0.D0) THEN
          NB=NB+1
          XB1(NB)=X-DX
          XB2(NB)=X
          IF(NB.EQ.NBB) RETURN
        END IF
        FP=FC
      END DO
      RETURN
      END
C===========================================================
C     FUNCION ZBRENT - REFINAMIENTO DE RAICES (METODO DE BRENT)
C     Numerical Recipes in Fortran, Cap. 9
C     Encuentra la raiz de FUNC en [X1,X2] con tolerancia TOL
C     Combina biseccion, secante e interpolacion inversa
C     Garantiza convergencia si hay cambio de signo en [X1,X2]
C===========================================================
      FUNCTION ZBRENT(FUNC,X1,X2,TOL)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (ITMAX=100, EPS=3.D-15)
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF(FB*FA.GT.0.D0) THEN
        WRITE(6,*) 'ZBRENT: la raiz no esta entre X1 y X2'
        ZBRENT=0.D0
        RETURN
      END IF
      FC=FB
      DO ITER=1,ITMAX
        IF(FB*FC.GT.0.D0) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        END IF
        IF(DABS(FC).LT.DABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        END IF
        TOL1=2.D0*EPS*DABS(B)+0.5D0*TOL
        XM=0.5D0*(C-B)
        IF(DABS(XM).LE.TOL1 .OR. FB.EQ.0.D0) THEN
          ZBRENT=B
          RETURN
        END IF
        IF(DABS(E).GE.TOL1 .AND. DABS(FA).GT.DABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
C           INTERPOLACION LINEAL (SECANTE)
            P=2.D0*XM*S
            Q=1.D0-S
          ELSE
C           INTERPOLACION INVERSA CUADRATICA
            Q=FA/FC
            R=FB/FC
            P=S*(2.D0*XM*Q*(Q-R)-(B-A)*(R-1.D0))
            Q=(Q-1.D0)*(R-1.D0)*(S-1.D0)
          END IF
          IF(P.GT.0.D0) THEN
            Q=-Q
          ELSE
            P=-P
          END IF
          IF(2.D0*P .LT. MIN(3.D0*XM*Q-DABS(TOL1*Q),
     &                        DABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          END IF
        ELSE
C         BISECCION
          D=XM
          E=D
        END IF
        A=B
        FA=FB
        IF(DABS(D).GT.TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        END IF
        FB=FUNC(B)
      END DO
      WRITE(6,*) 'ZBRENT: supero el maximo de iteraciones'
      ZBRENT=B
      RETURN
      END
C===========================================================
C     SUBRUTINA TRIDAG - SOLVER TRIDIAGONAL (ALGORITMO DE THOMAS)
C     Numerical Recipes in Fortran, Cap. 2
C     Resuelve el sistema tridiagonal: A*u = r
C     A(J) : subdiagonal  (j=2..N)
C     D(J) : diagonal principal (j=1..N)
C     C(J) : superdiagonal (j=1..N-1)
C     R(J) : vector del lado derecho
C     U(J) : solucion
C===========================================================
      SUBROUTINE TRIDAG(A,D,C,R,U,N)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=2000)
      DIMENSION A(N),D(N),C(N),R(N),U(N)
      DIMENSION GAM(NMAX)
      IF(D(1).EQ.0.D0) THEN
        WRITE(6,*) 'TRIDAG: error, D(1)=0'
        RETURN
      END IF
      BET=D(1)
      U(1)=R(1)/BET
C     DESCOMPOSICION Y SUSTITUCION HACIA ADELANTE
      DO J=2,N
        GAM(J)=C(J-1)/BET
        BET=D(J)-A(J)*GAM(J)
        IF(BET.EQ.0.D0) THEN
          WRITE(6,*) 'TRIDAG: fallo en el paso ',J
          RETURN
        END IF
        U(J)=(R(J)-A(J)*U(J-1))/BET
      END DO
C     SUSTITUCION HACIA ATRAS
      DO J=N-1,1,-1
        U(J)=U(J)-GAM(J+1)*U(J+1)
      END DO
      RETURN
      END
C===========================================================
