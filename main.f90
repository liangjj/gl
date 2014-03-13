MODULE sys_par
IMPLICIT NONE
SAVE
INTEGER,PARAMETER::N=100
INTEGER,PARAMETER::NS=N*N
INTEGER,PARAMETER::double=8
REAL(KIND=double),PARAMETER::eps=1.e-5
REAL(KIND=double),PARAMETER::eps1=5.e-5
REAL(KIND=double),PARAMETER::Pi=3.141592653
REAL(KIND=double),DIMENSION(2,3),PARAMETER::a=(/ 0.5*SQRT(3.),0.5,  -0.5*SQRT(3.),0.5,  0,-1. /)
REAL(KIND=double),DIMENSION(2),PARAMETER::b1=(/ 2*pi/SQRT(3.),-2*pi/3. /) !reciprotical vector
REAL(KIND=double),DIMENSION(2),PARAMETER::b2=(/ 0.*pi,4*pi/3. /)           
REAL(KIND=double),DIMENSION(2,3),PARAMETER::d=(/ SQRT(3.),0,-SQRT(3.)/2.,-1.5,-SQRT(3.)/2.,1.5 /)
END MODULE sys_par

PROGRAM graphene_light
USE sys_par
IMPLICIT NONE
REAL(KIND=double)::U,m,Ne,Ef
REAL(KIND=double)::Delta
Ne=0.1
U=2.8
m=0
!OPEN(11,FILE="pd_0.5p")
!DO m=0,1.01,0.05
!DO U=2.,3.01,0.1
  CALL meanfield(U,Ne,m,Delta,Ef)
!  WRITE(10,"(3F10.5)")m,U,Delta
  CALL drawband(Ef,m,Delta)
  WRITE(*,"(4F10.5)")m,U,Delta,Ef
!END DO
!END DO
!CLOSE(11)
END PROGRAM graphene_light


SUBROUTINE meanfield(U,Ne,m,Delta,Ef)
USE sys_par
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::U,Ne,m
REAL(KIND=double),INTENT(OUT)::Ef,Delta
REAL(KIND=double)::Delta1
Delta1=0.5
Ef=0.3
DO
CALL findef(U,Ne,m,Delta,Ef)
CALL findDelta(U,Ef,m,Delta1)
write(*,*)Delta1,Delta,Ef
IF(ABS(Delta-Delta1)<eps) EXIT
Delta=Delta1                      !else, update delta
END DO
END SUBROUTINE meanfield


SUBROUTINE findEf(U,Ne,m,Delta,Ef)
USE sys_par
USE mkl95_lapack
USE mkl95_precision
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::U,Ne,m,Delta
REAL(KIND=double),INTENT(OUT)::Ef
REAL(KIND=double),DIMENSION(2)::k
REAL(KIND=double)::Nab,upper,lower
REAL(KIND=double),DIMENSION(4)::Ek,nrk
INTEGER::h1,h2,i,j
COMPLEX(KIND=double),DIMENSION(4,4)::Hk
upper=3
lower=-3
outmost:DO
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
write(*,*) "Nab,Ef,Delta",Nab,Ef,Delta
IF(ABS(Nab-Ne)<eps1) THEN 
	EXIT outmost
ELSE IF(Nab<Ne) THEN           !ef is too low 
	lower=ef
    Ef=(Ef+upper)/2
ELSE 
	upper=Ef
	Ef=(lower+Ef)/2
END IF
END DO outmost
END SUBROUTINE findEf


SUBROUTINE findDelta(U,Ef,m,Delta)
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
	 tmp=Mk*Mk/SQRT(e2*Ef*Ef+Mk*Mk*Delta*Delta)
     wk(1)=SQRT(e2+Ef*Ef+Delta*Delta-2*tmp)
     wk(2)=SQRT(e2+Ef*Ef+Delta*Delta+2*tmp)
     DO i=1,2
	 tmp1=tmp1+(1+(2*i-3)*tmp)/wk(i)
     END DO
  END DO
END DO
tmp1=tmp1/REAL(4*Ns)
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
COMPLEX(KIND=double),INTENT(IN)::Delta
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


INCLUDE "drawband.f90"
INCLUDE "drawband1.f90"


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
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
