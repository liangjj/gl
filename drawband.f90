SUBROUTINE drawband(Ef,m,Delta)
USE sys_par
USE mkl95_lapack
USE mkl95_precision
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::Ef,m
COMPLEX(KIND=double),INTENT(IN)::Delta
REAL(KIND=double),DIMENSION(2)::k
COMPLEX(KIND=double),DIMENSION(4,4)::Hk
REAL(double)::cyclelimit,x,kx,ky
REAL(double),DIMENSION(4)::Ek
cyclelimit=4*PI/3/SQRT(3.)
OPEN(20,FILE="energyband.txt")
DO kx=0,cyclelimit,0.01
	ky=0
	k(1)=kx
	k(2)=ky
        CALL writehk(Ef,m,Delta,k,Hk)
		CALL heev(Hk,Ek,'V')
	WRITE(20,"(6F10.5)")kx,Ek,Ef
END DO

DO kx=cyclelimit,cyclelimit/2,-0.005
	ky=-SQRT(3.)*kx+4*Pi/3
	k(1)=kx
	k(2)=ky
	x=SQRT(ky*ky+(kx-cyclelimit)*(kx-cyclelimit))+cyclelimit
		CALL writehk(Ef,m,Delta,k,Hk)
		CALL heev(Hk,Ek,'V')
	WRITE(20,"(6F10.5)")x,Ek,Ef
END DO

DO kx=cyclelimit/2,0,-0.01
	ky=2*PI/3
	k(1)=kx
	k(2)=ky
	x=5./2.*cyclelimit-kx
	CALL writehk(Ef,m,Delta,k,Hk)
	CALL heev(Hk,Ek,'V')
	WRITE(20,"(6F10.5)")x,Ek,Ef
END DO

DO ky=2*PI/3,0,-0.01
	kx=0
	k(1)=kx
	k(2)=ky
	x=5./2.*cyclelimit+2*PI/3-ky
	CALL writehk(Ef,m,Delta,k,Hk)
	CALL heev(Hk,Ek,'V')
	WRITE(20,"(6F10.5)")x,Ek,Ef
END DO
CLOSE(20)
END SUBROUTINE drawband
