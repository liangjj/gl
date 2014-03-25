MODULE sys_par
IMPLICIT NONE
SAVE
INTEGER,PARAMETER::double=8
INTEGER(KIND=DOUBLE),PARAMETER::N=1000
INTEGER(KIND=DOUBLE),PARAMETER::NS=N*N
REAL(KIND=double),PARAMETER::eps=1.e-5
REAL(KIND=double),PARAMETER::eps1=5.e-5
REAL(KIND=double),PARAMETER::Pi=3.141592653589793238462643383279502884197169399
REAL(KIND=double),DIMENSION(2,3),PARAMETER::a=(/ 0.5*SQRT(3.),0.5,  -0.5*SQRT(3.),0.5,  0,-1. /)
REAL(KIND=double),DIMENSION(2),PARAMETER::b1=(/ 2*pi/SQRT(3.),-2*pi/3. /) !reciprotical vector
REAL(KIND=double),DIMENSION(2),PARAMETER::b2=(/ 0.*pi,4*pi/3. /)           
REAL(KIND=double),DIMENSION(2,3),PARAMETER::d=(/ SQRT(3.),0,-SQRT(3.)/2.,-1.5,-SQRT(3.)/2.,1.5 /)
END MODULE sys_par

PROGRAM graphene_light
USE sys_par
IMPLICIT NONE
REAL(KIND=double)::m,Ef
REAL(KIND=double)::Delta,g,Ne
Ne=0.02
m=0.04
!CALL ne_mu(m,Ef,Delta)
!CALL drawband(Ef,m,Delta)
OPEN(15,FILE="m_0.txt")
DO Delta=0.00001,0.001,0.0001
  CALL findEf_Newton(m,Delta,Ne,Ef)
  CALL f_Delta(Delta,m,Ef,g)
  WRITE(15,"(2F10.5)") 1./g,Delta
  WRITE(*,*)"U,Delta,Ne=0.2,m=0.04",1./g,Delta
END DO

DO Delta=0.001,0.1,0.001
  CALL findEf_Newton(m,Delta,Ne,Ef)
  CALL f_Delta(Delta,m,Ef,g)
  WRITE(15,"(2F10.5)") 1./g,Delta
  WRITE(*,*)"U,Delta,Ne=0.2,m=0.04",1./g,Delta
END DO

CLOSE(15)
END PROGRAM graphene_light


SUBROUTINE findNe(U,m,Ef,Delta,Ne)
USE sys_par
USE mkl95_lapack
USE mkl95_precision
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::U,Ef,m,Delta
REAL(KIND=double),INTENT(OUT)::Ne
REAL(KIND=double),DIMENSION(2)::k
REAL(KIND=double),DIMENSION(4)::Ek,nrk
INTEGER::h1,h2,i,j
COMPLEX(KIND=double),DIMENSION(4,4)::Hk
Ne=0
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
        Ne=Ne+ABS(hk(1,i))*ABS(hk(1,i))*nrk(i)+ABS(hk(2,i))*ABS(hk(2,i))*nrk(i)
		Ne=Ne-ABS(hk(3,i))*ABS(hk(3,i))*nrk(i)-ABS(hk(4,i))*ABS(hk(4,i))*nrk(i)
     END DO
  END DO
END DO
Ne=Ne/REAL(Ns)
END SUBROUTINE findNe


SUBROUTINE writePar(m,k,Mk,tk)
USE sys_par
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::m
REAL(KIND=double),DIMENSION(2),INTENT(IN)::k
REAL(KIND=double),INTENT(OUT)::Mk
COMPLEX(KIND=double),INTENT(OUT)::tk
REAL(KIND=double)::tk_re,tk_im,tmp
INTEGER::i
tk_re=0
tk_im=0
tk=0
mk=0
DO i=1,3
  mk=mk+2*m/3./SQRT(3.)*SIN(DOT_PRODUCT(k,d(:,i)))
  tmp=DOT_PRODUCT(k,a(:,i))
  tk_re=tk_re+COS(tmp)
  tk_im=tk_im+SIN(tmp)
END DO
tk=CMPLX(tk_re,tk_im,double)
END SUBROUTINE writePar


SUBROUTINE writehk(Ef,m,Delta,k,Hk)
USE sys_par
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::Ef,m
REAL(KIND=double),INTENT(IN)::Delta
REAL(KIND=double),DIMENSION(2),INTENT(IN)::k
COMPLEX(KIND=double),DIMENSION(4,4),INTENT(OUT)::Hk
COMPLEX(KIND=double)::tk
REAL(KIND=double)::Mk
CALL writePar(m,k,Mk,tk)
 hk=0
 hk(1,1)=mk-ef
 hk(2,2)=-mk-ef
 hk(3,3)=mk+ef
 hk(4,4)=-mk+ef
 hk(1,2)=-tk
 hk(1,3)=-Delta
 hk(2,4)=-Delta
 hk(3,4)=tk
 hk(2,1)=-CONJG(tk)
 hk(3,1)=-Delta
 hk(4,2)=-Delta
 hk(4,3)=CONJG(tk)
END SUBROUTINE writehk


INCLUDE "drawband.f90"   !draw band around BZ boundary
INCLUDE "drawband1.f90"  !draw the full band, 3D
INCLUDE "f_Delta.f90"    !calculate mu!=0 delta vs U,with and without m
INCLUDE "g_Delta.f90"    !calculate mu=0 delta vs U, with and without m
INCLUDE "ne_mu.f90"      !find Ne vs mu with no interaction
INCLUDE "findEf.f90"     !find Ef self consistantly with interaction 
INCLUDE "findEf_Newton.f90"     !find Ef self consistantly with interaction using Newton method

SUBROUTINE findDelta(U,m,Ef,Delta)
USE sys_par
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::U,Ef,m
REAL(KIND=double),INTENT(OUT)::Delta
REAL(KIND=double),DIMENSION(2)::k,wk
REAL(KIND=double)::Mk,tmp,tmp1,upper,lower,e2
COMPLEX(KIND=double)::tk
INTEGER::h1,h2,i
upper=3
lower=0
outmost: DO
tmp1=0
DO h1=1,N
  DO h2=1,N
     k=h1/REAL(N)*b1+h2/REAL(N)*b2
     CALL writePar(m,k,Mk,tk)
	 e2=Mk*Mk+ABS(tk)*ABS(tk)
	 tmp=SQRT(e2*Ef*Ef+Mk*Mk*Delta*Delta)
     wk(1)=SQRT(e2+Ef*Ef+Delta*Delta-2*tmp)
     wk(2)=SQRT(e2+Ef*Ef+Delta*Delta+2*tmp)
	 IF(m<eps) THEN
	    tmp=0
	 ELSE 
	    tmp=Mk*Mk/SQRT(e2*Ef*Ef+Mk*Mk*Delta*Delta)
		write(*,*) tmp
	 END IF
     DO i=1,2
	 tmp1=tmp1+(1+(2*i-3)*tmp)/wk(i)
     END DO
  END DO
END DO
tmp1=tmp1/REAL(4*Ns)
!write(*,*)"tmp1,Delta",tmp1,Delta
IF(ABS(Delta/U-Delta*tmp1)<eps) THEN 
	EXIT outmost
ELSE IF(tmp1>1/U) THEN           !the right handside is decress function of Delta 
	lower=Delta
    Delta=(Delta+upper)/2
ELSE 
	upper=Delta
	Delta=(lower+Delta)/2
END IF
END DO outmost
END SUBROUTINE findDelta

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! DESCRIPTION:
! This program deal with "negative U Hubbard" model, with circular light shed on it. Originally, the  
! negative U hubbard model will lead to the S-wave pairing superconducting states at some particular  
! U value. After shed light on graphene, this U value might change. Our goal is to find out how the   
! light will affect the phase diagram of this model. We will first use an effective model to represent
! the circular light. 
!-----------------------------------------------------------------------------------------------------
! AUTHOR:
! Luming   |   hiluming@gmail.com
!-----------------------------------------------------------------------------------------------------
! VERSION 1.0
! graphene negative hubbard U model, without sheding light on it.
! DATE:
! Thu Feb 20 14:15:37 CST 2014
!-----------------------------------------------------------------------------------------------------
! VERSION 1.1
! add subroutine to self-consistent dertermine ef
! DATE:
! Wed Feb 26 11:07:41 CST 2014
!-----------------------------------------------------------------------------------------------------
! VERSION 1.2
! add subroutine drawband. The above versions don't work fine
! DATE:
! Wed Feb 26 22:32:33 CST 2014
!-----------------------------------------------------------------------------------------------------
! VERSION 1.3
! Using Nambu representation, to make the matrix 4*4, also figure out the main bug---that is ---
! forget initialize the matrix. The bug is fixed. So this version works fine
! DATE:
! Wed Feb 26 11:07:41 CST 2014
!-----------------------------------------------------------------------------------------------------
! VERSION 1.4
! add subroutine drawbandfull(drawband1.f90). modify all the delta_norm to delta_a and delta_b,which 
! is complex number. Also modify the method to compute particle number, using qusi-particle number 
! for counting. 
! DATE:
! Fri Feb 28 23:19:26 CST 2014
!-----------------------------------------------------------------------------------------------------
! VERSION 1.5
! Change the variable name to the simple ones.
! Solve the gap consistently using gap equation
! Change the qusiparticle "Fermi Energy" to zero
! DATE:
! Thu Mar 13 21:48:08 CST 2014
!-----------------------------------------------------------------------------------------------------
! VERSION 1.6
! Add g_Delta.f90 and f_Delta.f90. Which are the function to obtain U using delta, half-filled and doped
! respectively. Another trick is to use the Delta limits to 0 instead of setting Delta to 0 directly. 
! This avoid the divergent problem. 
! DATE:
! Tue Mar 18 19:24:13 CST 2014
!-----------------------------------------------------------------------------------------------------
! VERSION 1.7
! Add findEf.f90 and findEf_Newton.f90 to find Ef with a given Ne selfconsistently. Also add a simple 
! program ne_mu to find the relation bewteen Ne and Mu without interaction. Change two places that 
! wrongly define Delta as a complex number.
! DATE:
!Tue Mar 25 22:16:09 CST 2014
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
