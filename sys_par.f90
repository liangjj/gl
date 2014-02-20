MODULE sys_par
IMPLICIT NONE
SAVE
INTEGER,PARAMETER::n_lattice=100
INTEGER,PARAMETER::n_site=n_lattice*n_lattice
INTEGER,PARAMETER::double=8
REAL(KIND=double),PARAMETER::eps=1.e-6
REAL(KIND=double),PARAMETER::pi=3.141592653
REAL(KIND=double),DIMENSION(2,3),PARAMETER::delta=(/ 0.5*SQRT(3.),0.5,  -0.5*SQRT(3.),0.5,  0,-1. /)
REAL(KIND=double),DIMENSION(2),PARAMETER::a1star=(/ 2*pi/SQRT(3.),-2*pi/3. /) !reciprotical vector
REAL(KIND=double),DIMENSION(2),PARAMETER::a2star=(/ 0.*pi,4*pi/3. /)           
REAL(KIND=double),DIMENSION(2,3),PARAMETER::d_vec=(/ SQRT(3.),0,-SQRT(3.)/2.,-1.5,-SQRT(3.)/2.,1.5 /)
END MODULE sys_par

PROGRAM MAIN
USE sys_par
IMPLICIT NONE
write(*,*) delta
END PROGRAM
