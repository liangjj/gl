SUBROUTINE drawbandfull(ef,m_lit,delta_a,delta_b)
USE sys_par
USE mkl95_lapack
USE mkl95_precision
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::ef,m_lit
COMPLEX(KIND=double),INTENT(IN)::delta_a,delta_b
REAL(KIND=double),DIMENSION(2)::k_vec
COMPLEX(KIND=double),DIMENSION(4,4)::hk
INTEGER::h1,h2
REAL(double),DIMENSION(4)::ek
OPEN(11,FILE="fullband.txt")
DO h1=1,n_lattice                              !calculate particle number
  DO h2=1,n_lattice
     k_vec=h1/REAL(n_lattice)*a1star+h2/REAL(n_lattice)*a2star
     CALL writehk(ef,m_lit,delta_a,delta_b,k_vec,hk)
     CALL heev(hk,ek,'V')
	 WRITE(11,"(5F10.5)") k_vec,ek(2:3),ef
  END DO
END DO
CLOSE(11)
END SUBROUTINE drawbandfull
