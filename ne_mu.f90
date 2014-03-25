SUBROUTINE ne_mu(m,Ef,Delta)
USE sys_par
USE mkl95_lapack
USE mkl95_precision
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::m,Ef,Delta
REAL(KIND=double)::Ne,mu
REAL(KIND=double),DIMENSION(2)::k
REAL(KIND=double),DIMENSION(4)::Ek,nrk
INTEGER::h1,h2,i,j
COMPLEX(KIND=double),DIMENSION(4,4)::Hk
!OPEN(19,FILE="mu_ne.txt")
!DO mu=0,0.3,0.01
mu=0.3
Ne=0
DO h1=1,N                              !calculate particle number
  DO h2=1,N
     k=h1/REAL(N)*b1+h2/REAL(N)*b2
     CALL writehk(Ef,m,Delta,k,Hk)
     CALL heev(Hk,Ek,'V')
     DO i=1,4
        IF(Ek(i)>mu) THEN !Bogoliubov quasi-particles fermi energy is 0
           nrk(i)=0.
        ELSE
           nrk(i)=1.
        END IF
        Ne=Ne+nrk(i)
     END DO
  END DO
END DO
Ne=Ne/REAL(Ns)
write(*,*) Ne
END SUBROUTINE ne_mu


