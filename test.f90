!MODULE sys_par
!IMPLICIT NONE
!INTEGER,PARAMETER:: double=8
!INTEGER,DIMENSION(2,3),PARAMETER::mat=(/1,2,3,4,5,6/)
!REAL(KIND=double),PARAMETER::pi=3.1415
!REAL(KIND=double),DIMENSION(2),PARAMETER::a2star=(/0.*5,pi/)
!END MODULE sys_par

PROGRAM MAIN
!USE sys_par
IMPLICIT NONE
integer,parameter::double=8
REAL(KIND=8),parameter::pi=3.141592653
!REAL(KIND=8),DIMENSION(2)::a,b
COMPLEX(KIND=8)::x
real(kind=8) a,b
a=1
b=1
x=cmplx(a,b,double)
 write(*,*)x
END PROGRAM MAIN


