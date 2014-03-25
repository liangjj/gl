SUBROUTINE findEf(m,Delta,Ne,Ef)
USE sys_par
USE mkl95_lapack
USE mkl95_precision
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::m,Delta,Ne
REAL(KIND=double),INTENT(OUT)::Ef
REAL(KIND=double)::Nab,lower,upper
REAL(KIND=double),DIMENSION(2)::k
REAL(KIND=double),DIMENSION(4)::Ek,nrk
INTEGER::h1,h2,i,j
COMPLEX(KIND=double),DIMENSION(4,4)::Hk
EF=0.1
upper=1.
lower=0
outmost: DO
Nab=0
DO h1=1,N                              !calculate particle number
  DO h2=1,N
     k=h1/REAL(N)*b1+h2/REAL(N)*b2
     CALL writehk(Ef,m,Delta,k,Hk)
     CALL heev(Hk,Ek,'V')
     DO i=1,4
        IF(Ek(i)>0) THEN !Bogoliubov quasi-particles fermi energy is 0
           nrk(i)=0.
        ELSE
           nrk(i)=1.
        END IF
        Nab=Nab+ABS(hk(1,i))*ABS(hk(1,i))*nrk(i)+ABS(hk(2,i))*ABS(hk(2,i))*nrk(i)
		Nab=Nab-ABS(hk(3,i))*ABS(hk(3,i))*nrk(i)-ABS(hk(4,i))*ABS(hk(4,i))*nrk(i)
     END DO
  END DO
END DO
Nab=Nab/REAL(Ns)
write(*,*) "Nab,Ef",Nab,Ef
IF(ABS(Nab-Ne)<eps) THEN
  EXIT outmost
ELSE IF(Nab<Ne) THEN
  lower=Ef
  Ef=(EF+upper)/2
ELSE
  upper=Ef
  Ef=(Ef+lower)/2
END IF

END DO outmost
END SUBROUTINE findEf



