MODULE sys_par
IMPLICIT NONE
SAVE
INTEGER,PARAMETER::n_lattice=100
INTEGER,PARAMETER::n_site=n_lattice*n_lattice
INTEGER,PARAMETER::double=8
REAL(KIND=double),PARAMETER::eps=1.e-5
REAL(KIND=double),PARAMETER::eps1=1.e-3
REAL(KIND=double),PARAMETER::pi=3.141592653
REAL(KIND=double),DIMENSION(2,3),PARAMETER::delta=(/ 0.5*SQRT(3.),0.5,  -0.5*SQRT(3.),0.5,  0,-1. /)
REAL(KIND=double),DIMENSION(2),PARAMETER::a1star=(/ 2*pi/SQRT(3.),-2*pi/3. /) !reciprotical vector
REAL(KIND=double),DIMENSION(2),PARAMETER::a2star=(/ 0.*pi,4*pi/3. /)           
REAL(KIND=double),DIMENSION(2,3),PARAMETER::d_vec=(/ SQRT(3.),0,-SQRT(3.)/2.,-1.5,-SQRT(3.)/2.,1.5 /)
END MODULE sys_par

PROGRAM graphene_light
USE sys_par
IMPLICIT NONE
REAL(KIND=double)::u_hub,m_lit,n_dope,ef
COMPLEX(KIND=double)::delta_a,delta_b
n_dope=-0.01
!u_hub=3
!m_lit=0
OPEN(11,FILE="pd_0.5p")
DO m_lit=0,1.01,0.05
DO u_hub=0.,4.01,0.05
  CALL meanfield(u_hub,n_dope,m_lit,delta_a,delta_b,ef)
  WRITE(10,"(3F10.5)")m_lit,u_hub,ABS(delta_a)
!  CALL drawbandfull(ef,m_lit,delta_a,delta_b)
  WRITE(*,"(4F10.5)")m_lit,u_hub,ABS(delta_a),ABS(delta_b),ef
END DO
END DO
CLOSE(11)
END PROGRAM graphene_light




SUBROUTINE meanfield(u_hub,n_dope,m_lit,delta_a,delta_b,ef)
USE sys_par
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::u_hub,n_dope,m_lit
REAL(KIND=double),INTENT(OUT)::ef
COMPLEX(KIND=double),INTENT(OUT)::delta_a,delta_b
COMPLEX(KIND=double)::delta_a1,delta_b1
delta_a=0.5
delta_b=0.4
DO
CALL findef(u_hub,n_dope,m_lit,delta_a,delta_b,ef)
CALL findnext(u_hub,ef,m_lit,delta_a,delta_b,delta_a1,delta_b1)
IF(ABS(ABS(delta_a)-ABS(delta_a1))<eps.AND.ABS(ABS(delta_b)-ABS(delta_b1))<eps) EXIT
delta_a=delta_a1                      !else, update delta_a
delta_b=delta_b1
END DO
END SUBROUTINE meanfield


SUBROUTINE findef(u_hub,n_dope,m_lit,delta_a,delta_b,ef)
USE sys_par
USE mkl95_lapack
USE mkl95_precision
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::u_hub,n_dope,m_lit
COMPLEX(KIND=double),INTENT(IN)::delta_a,delta_b
REAL(KIND=double),INTENT(OUT)::ef
REAL(KIND=double),DIMENSION(2)::k_vec
REAL(KIND=double)::nab,upper,lower
REAL(KIND=double),DIMENSION(4)::ek,nrk
INTEGER::h1,h2,i,j
COMPLEX(KIND=double),DIMENSION(4,4)::hk
upper=3
lower=-3
outmost:DO
nab=0
DO h1=1,n_lattice                              !calculate particle number
  DO h2=1,n_lattice
     k_vec=h1/REAL(n_lattice)*a1star+h2/REAL(n_lattice)*a2star
     CALL writehk(ef,m_lit,delta_a,delta_b,k_vec,hk)
     CALL heev(hk,ek,'V')
     DO i=1,4
        IF(ek(i)>ef) THEN
           nrk(i)=0.
        ELSE
           nrk(i)=1.
        END IF
	 END DO
	 DO i=1,4
	 nab=nab+nrk(i)
       ! nab=nab+ABS(hk(1,i))*ABS(hk(1,i))*nrk(i)+ABS(hk(2,i))*ABS(hk(2,i))*nrk(i)
		!nab=nab+ABS(hk(3,i))*ABS(hk(3,i))*(1-nrk(i))+ABS(hk(4,i))*ABS(hk(4,i))*(1-nrk(i))
     END DO
  END DO
END DO
nab=nab/REAL(n_site)
write(*,*) nab,ef,abs(delta_a),abs(delta_b)
IF(ABS(nab-2.-n_dope)<eps1) THEN 
	EXIT outmost
ELSE IF(nab<2+n_dope) THEN           !ef is too low 
	lower=ef
    ef=(ef+upper)/2
ELSE 
	upper=ef
	ef=(lower+ef)/2
END IF
END DO outmost
END SUBROUTINE findef


SUBROUTINE findnext(u_hub,ef,m_lit,delta_a,delta_b,delta_a1,delta_b1)
USE sys_par
USE mkl95_lapack
USE mkl95_precision
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::u_hub,ef,m_lit
COMPLEX(KIND=double),INTENT(IN)::delta_a,delta_b
COMPLEX(KIND=double),INTENT(OUT)::delta_a1,delta_b1
REAL(KIND=double),DIMENSION(2)::k_vec
REAL(KIND=double),DIMENSION(4)::ek,nrk
INTEGER::h1,h2,i
COMPLEX(KIND=double),DIMENSION(4,4)::hk
delta_a1=0
delta_b1=0
DO h1=1,n_lattice
  DO h2=1,n_lattice
     k_vec=h1/REAL(n_lattice)*a1star+h2/REAL(n_lattice)*a2star
     CALL writehk(ef,m_lit,delta_a,delta_b,k_vec,hk)
     CALL heev(hk,ek,'V')
     DO i=1,4
        IF(ek(i)>ef) THEN
           nrk(i)=0
        ELSE
           nrk(i)=1
        END IF
        delta_a1=delta_a1+u_hub*CONJG(hk(3,i))*hk(1,i)*nrk(i)
        delta_b1=delta_b1+u_hub*CONJG(hk(4,i))*hk(2,i)*nrk(i)
     END DO
  END DO
END DO
delta_a1=delta_a1/REAL(n_site)
delta_b1=delta_b1/REAL(n_site)
END SUBROUTINE findnext


SUBROUTINE writehk(ef,m_lit,delta_a,delta_b,k_vec,hk)
USE sys_par
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::ef,m_lit
COMPLEX(KIND=double),INTENT(IN)::delta_a,delta_b
REAL(KIND=double),DIMENSION(2),INTENT(IN)::k_vec
COMPLEX(KIND=double),DIMENSION(4,4),INTENT(OUT)::hk
COMPLEX(KIND=double)::tk
REAL(KIND=double)::tk_re,tk_im,mk,tmp
INTEGER::i,j
tk_re=0
tk_im=0
tk=0
mk=0
DO i=1,3
  mk=mk+2*m_lit/3./SQRT(3.)*SIN(DOT_PRODUCT(k_vec,d_vec(:,i)))
  tmp=DOT_PRODUCT(k_vec,delta(:,i))
  tk_re=tk_re+COS(tmp)
  tk_im=tk_im+SIN(tmp)
END DO
tk=CMPLX(tk_re,tk_im,double)
 hk=0
 hk(1,1)=mk-ef
 hk(2,2)=-mk-ef
 hk(3,3)=mk+ef
 hk(4,4)=-mk+ef
 hk(1,2)=-tk
 hk(1,3)=-delta_a
 hk(2,4)=-delta_b
 hk(3,4)=tk
 hk(2,1)=-CONJG(tk)
 hk(3,1)=-CONJG(delta_a)
 hk(4,2)=-CONJG(delta_b)
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
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
