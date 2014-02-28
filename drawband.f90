SUBROUTINE drawband(ef,m_lit,delta_a,delta_b)
USE sys_par
USE mkl95_lapack
USE mkl95_precision
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::ef,m_lit
COMPLEX(KIND=double),INTENT(IN)::delta_a,delta_b
REAL(KIND=double),DIMENSION(2)::k_vec
COMPLEX(KIND=double),DIMENSION(4,4)::hk
REAL(double)::cyclelimit,k,kx,ky
REAL(double),DIMENSION(4)::ek
cyclelimit=4*PI/3/SQRT(3.)
OPEN(20,FILE="energyband.txt")
DO kx=0,cyclelimit,0.01
	ky=0
	k_vec(1)=kx
	k_vec(2)=ky
		CALL writehk(ef,m_lit,delta_a,delta_b,k_vec,hk)
		CALL heev(hk,ek,'V')
	WRITE(20,"(5F10.5)")kx,ek
END DO

DO kx=cyclelimit,cyclelimit/2,-0.005
	ky=-SQRT(3.)*kx+4*Pi/3
	k_vec(1)=kx
	k_vec(2)=ky
	k=SQRT(ky*ky+(kx-cyclelimit)*(kx-cyclelimit))+cyclelimit
		CALL writehk(ef,m_lit,delta_a,delta_b,k_vec,hk)
		CALL heev(hk,ek,'V')
	WRITE(20,"(5F10.5)")k,ek
END DO

DO kx=cyclelimit/2,0,-0.01
	ky=2*PI/3
	k_vec(1)=kx
	k_vec(2)=ky
	k=5./2.*cyclelimit-kx
	CALL writehk(ef,m_lit,delta_a,delta_b,k_vec,hk)
	CALL heev(hk,ek,'V')
	WRITE(20,"(5F10.5)")k,ek
END DO

DO ky=2*PI/3,0,-0.01
	kx=0
	k_vec(1)=kx
	k_vec(2)=ky
	k=5./2.*cyclelimit+2*PI/3-ky
	CALL writehk(ef,m_lit,delta_a,delta_b,k_vec,hk)
	CALL heev(hk,ek,'V')
	WRITE(20,"(5F10.5)")k,ek
END DO
CLOSE(20)
END SUBROUTINE drawband
