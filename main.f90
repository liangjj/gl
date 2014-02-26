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
! add subroutine drawband. The above version doesn't work fine
! DATE:
! Wed Feb 26 22:32:33 CST 2014
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
MODULE sys_par
IMPLICIT NONE
SAVE
INTEGER,PARAMETER::n_lattice=100
INTEGER,PARAMETER::n_site=n_lattice*n_lattice
INTEGER,PARAMETER::double=8
REAL(KIND=double),PARAMETER::eps=1.e-5
REAL(KIND=double),PARAMETER::eps1=1.e-4
REAL(KIND=double),PARAMETER::pi=3.141592653
REAL(KIND=double),DIMENSION(2,3),PARAMETER::delta=(/ 0.5*SQRT(3.),0.5,  -0.5*SQRT(3.),0.5,  0,-1. /)
REAL(KIND=double),DIMENSION(2),PARAMETER::a1star=(/ 2*pi/SQRT(3.),-2*pi/3. /) !reciprotical vector
REAL(KIND=double),DIMENSION(2),PARAMETER::a2star=(/ 0.*pi,4*pi/3. /)           
REAL(KIND=double),DIMENSION(2,3),PARAMETER::d_vec=(/ SQRT(3.),0,-SQRT(3.)/2.,-1.5,-SQRT(3.)/2.,1.5 /)
END MODULE sys_par

PROGRAM graphene_light
USE sys_par
IMPLICIT NONE
REAL(KIND=double)::u_hub,m_lit,n_dope,delta_norm,ef
n_dope=0.
u_hub=0
m_lit=0
ef=0
OPEN(10,FILE="pd3")
!DO m_lit=0,1.01,0.1
!DO u_hub=3.,6.01,0.1
!  CALL meanfield(u_hub,n_dope,m_lit,delta_norm)
  delta_norm=0
  CALL drawband(u_hub,ef,m_lit,delta_norm)
!  WRITE(10,"(3F10.5)")m_lit,u_hub,delta_norm
!  WRITE(*,"(3F10.5)")m_lit,u_hub,delta_norm
!END DO
!END DO
CLOSE(10)
END PROGRAM graphene_light


INCLUDE "drawband.f90"


SUBROUTINE meanfield(u_hub,n_dope,m_lit,delta_norm)
USE sys_par
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::u_hub,n_dope,m_lit
REAL(KIND=double),INTENT(OUT)::delta_norm
REAL(KIND=double)::delta_norm1,ef
delta_norm=0
DO
CALL findef(u_hub,n_dope,m_lit,delta_norm,ef)
CALL findnext(u_hub,ef,m_lit,delta_norm,delta_norm1)
IF(ABS(delta_norm-delta_norm1)<eps) EXIT
delta_norm=delta_norm1                      !else, update delta_norm
!write(*,*) ef,delta_norm
END DO
END SUBROUTINE meanfield


SUBROUTINE findef(u_hub,n_dope,m_lit,delta_norm,ef)
USE sys_par
USE mkl95_lapack
USE mkl95_precision
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::u_hub,n_dope,m_lit,delta_norm
REAL(KIND=double),INTENT(OUT)::ef
REAL(KIND=double),DIMENSION(2)::k_vec
REAL(KIND=double)::nab,upper,lower
REAL(KIND=double),DIMENSION(8)::ek,nrk
INTEGER::h1,h2,i,j
COMPLEX(KIND=double),DIMENSION(8,8)::hk
nab=0
ef=0
upper=3
lower=-3
outmost:DO 
DO h1=1,n_lattice                              !calculate particle number
  DO h2=1,n_lattice
     k_vec=h1/REAL(n_lattice)*a1star+h2/REAL(n_lattice)*a2star
     CALL writehk(u_hub,ef,m_lit,delta_norm,k_vec,hk)
     CALL heev(hk,ek,'V')
     DO i=1,8
        IF(ek(i)>ef) THEN
           nrk(i)=0
        ELSE
           nrk(i)=1
        END IF
		DO j=1,4
        	nab=nab+CONJG(hk(j,i))*hk(j,i)*nrk(i)
		END DO
     END DO
  END DO
END DO
nab=nab/REAL(n_site)

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


SUBROUTINE findnext(u_hub,ef,m_lit,delta_norm,delta_norm1)
USE sys_par
USE mkl95_lapack
USE mkl95_precision
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::u_hub,ef,m_lit,delta_norm
REAL(KIND=double),INTENT(OUT)::delta_norm1
REAL(KIND=double),DIMENSION(2)::k_vec
REAL(KIND=double),DIMENSION(8)::ek,nrk
COMPLEX(KIND=double)::delta_a
INTEGER::h1,h2,i
COMPLEX(KIND=double),DIMENSION(8,8)::hk
delta_a=0
DO h1=1,n_lattice
  DO h2=1,n_lattice
     k_vec=h1/REAL(n_lattice)*a1star+h2/REAL(n_lattice)*a2star
     CALL writehk(u_hub,ef,m_lit,delta_norm,k_vec,hk)
     CALL heev(hk,ek,'V')
     DO i=1,8
        IF(ek(i)>ef) THEN
           nrk(i)=0
        ELSE
           nrk(i)=1
        END IF
        delta_a=delta_a+CONJG(hk(6,i))*hk(1,i)*nrk(i)
     END DO
  END DO
END DO
delta_norm1=ABS(delta_a)/REAL(n_site)
END SUBROUTINE findnext


SUBROUTINE writehk(u_hub,ef,m_lit,delta_norm,k_vec,hk)
USE sys_par
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::u_hub,ef,m_lit,delta_norm
REAL(KIND=double),DIMENSION(2),INTENT(IN)::k_vec
COMPLEX(KIND=double),DIMENSION(8,8),INTENT(OUT)::hk
COMPLEX(KIND=double)::tk
REAL(KIND=double)::tk_re,tk_im,mk
INTEGER::i,j
tk_re=0
tk_im=0
mk=0
DO i=1,3
  mk=mk+2*m_lit/3./SQRT(3.)*SIN(DOT_PRODUCT(k_vec,d_vec(:,i)))
  tk_re=tk_re+COS(DOT_PRODUCT(k_vec,delta(:,i)))
  tk_im=tk_im+SIN(DOT_PRODUCT(k_vec,delta(:,i)))
END DO
tk=CMPLX(tk_re,tk_im,double)
DO i=1,2
 hk(i,i)=mk-ef
 hk(i+2,i+2)=-mk-ef
 hk(i+4,i+4)=mk+ef
 hk(i+6,i+6)=-mk+ef
 hk(i,i+2)=-tk
 hk(i+4,i+6)=tk
END DO
hk(2,5)=u_hub*delta_norm
hk(1,6)=-u_hub*delta_norm
hk(3:4,7:8)=hk(1:2,5:6)
DO i=1,8
 DO j=1,8
  IF(i<j) hk(j,i)=CONJG(hk(i,j))
 END DO
END DO
hk=0.5*hk
END SUBROUTINE writehk


