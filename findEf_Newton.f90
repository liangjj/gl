SUBROUTINE findEf_Newton(m,Delta,Ne,Ef)
USE sys_par
USE mkl95_lapack
USE mkl95_precision
IMPLICIT NONE
REAL(KIND=double),INTENT(IN)::m,Delta,Ne
REAL(KIND=double),INTENT(OUT)::Ef
REAL(KIND=double)::Nab,fk_1,fk,uk_1,uk,uk1
REAL(KIND=double),DIMENSION(2)::k
REAL(KIND=double),DIMENSION(4)::Ek,nrk
INTEGER::h1,h2,i,j
COMPLEX(KIND=double),DIMENSION(4,4)::Hk
uk=0.15
uk_1=0.1
fk_1=-0.002
outmost: DO
Nab=0
DO h1=1,N                              !calculate particle number
  DO h2=1,N
     k=h1/REAL(N)*b1+h2/REAL(N)*b2
     CALL writehk(uk,m,Delta,k,Hk)     !using uk to obtain fk
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
fk=Nab-Ne
uk1=uk-fk/((fk-fk_1)/(uk-uk_1))         !formula to obtain  uk1
!write(*,*) "Ef,fk",uk1,fk
IF(ABS(fk)<eps) THEN
  EXIT outmost
ELSE 
  fk_1=fk
  uk_1=uk
  uk=uk1
END IF
END DO outmost
Ef=uk1
END SUBROUTINE findEf_Newton

