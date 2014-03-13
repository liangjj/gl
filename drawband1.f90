SUBROUTINE drawbandfull(Ef,m,Delta)
USE sys_par
USE mkl95_lapack
USE mkl95_precision
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::Ef,m
COMPLEX(KIND=double),INTENT(IN)::Delta
REAL(KIND=double),DIMENSION(2)::k
COMPLEX(KIND=double),DIMENSION(4,4)::Hk
INTEGER::h1,h2
REAL(double),DIMENSION(4)::Ek
OPEN(11,FILE="fullband.txt")
DO h1=1,N                              !calculate particle number
  DO h2=1,N
     k=h1/REAL(N)*b1+h2/REAL(N)*b2
     CALL writehk(Ef,m,Delta,k,Hk)
     CALL heev(Hk,Ek,'V')
	 WRITE(11,"(5F10.5)") k,Ek(2:3),ef
  END DO
END DO
CLOSE(11)
END SUBROUTINE drawbandfull
