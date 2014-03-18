SUBROUTINE f_Delta(Delta,m,Ef,g)
USE sys_par
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::Delta,m,Ef
REAL(KIND=double),INTENT(OUT)::g
REAL(KIND=double),DIMENSION(2)::k,wk
REAL(KIND=double)::Mk,e2,tmp
COMPLEX(KIND=double)::tk
INTEGER::h1,h2,i
g=0
DO h1=1,N
  DO h2=1,N
     k=h1/REAL(N)*b1+h2/REAL(N)*b2
     CALL writePar(m,k,Mk,tk)
	 Mk=ABS(Mk)
	 e2=Mk*Mk+ABS(tk)*ABS(tk)
	 tmp=SQRT(e2*Ef*Ef+Mk*Mk*Delta*Delta)
     wk(1)=SQRT(e2+Delta*Delta+Ef*Ef-2*tmp)
     wk(2)=SQRT(e2+Delta*Delta+Ef*Ef+2*tmp)
	 
	 DO i=1,2
	 g=g+(Delta+(2*i-3)*Mk*Mk/tmp)/wk(i)
     END DO
  END DO
END DO
g=g/REAL(4*Ns)
END SUBROUTINE f_Delta
