!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE getEl_DIC_M_S(IELEM,ELCOR,NODEL,NDIME,&
		   NGAUSS,NORDER,K_W,Omega,&
		   ELM,ELK,PTLD_Phi,PTLD_Vlcty,&
		   NANGL,NODEF,BETA,ENR_Type)
           !Description: This subroutine builds elementary 
           		!matrices for Mass, Stiffness. 
           		!For other purposes, this could 
           		!be updated to compute decompositions
           		!for example, DIC for u0,g0 etc.
           !Arguments out
           		!ELM: Elementary Mass 
           		!ELK: Elementary stiffness 

IMPLICIT NONE
double precision ELCOR
integer NODEL,NDIME,NORDER
DIMENSION ELCOR(NODEL,NDIME)

INTEGER NGAUSS, IGAUSS, JGAUSS,I,J, IELEM
DOUBLE PRECISION WG(NGAUSS), XG(NGAUSS),&
                 XI, ETA, WTX, WTY, WTXY,&
                 X,Y
                 
double precision SF(NODEL), SFDL(NDIME,NODEL)
double precision JCBT(NDIME,NDIME),&
                 JCBTI(NDIME,NDIME), DETJ


double precision n(NDIME)
integer EDG_TYPE

!definitions for DIC,M,S
double precision SFDG(NDIME,NODEL),&
                 r_mod,K_W,Omega
double complex	 ELM(NODEF,NODEF),&
		 ELK(NODEF,NODEF),&
		 FZ_Phi,FZ_Vlcty,&
		 PTLD_Phi(NODEF),&
		 PTLD_Vlcty(NODEF)
		 
!definitions for PUFEM
integer NANGL,NODEF
double complex SFN(NODEF),SFDGN(NDIME,NODEF)
    !definitions for ENR_Type
    integer ENR_type		! 1=ENR2, 0=ENR1
    !changes for ENR_Type
    double precision BETA(NANGL-ENR_type)

  
  !Initialize elementary matrices and vectors
  ELM = dcmplx(0.0D0,0.0D0)
  ELK = dcmplx(0.0D0,0.0D0)
  PTLD_Phi = dcmplx(0.0D0,0.0D0)
  PTLD_Vlcty = dcmplx(0.0D0,0.0D0)
  
  !Get integration points.
  call GAULEG(NGAUSS, XG, NGAUSS, WG, NGAUSS)
  
  !Integrate inside domain.
    DO IGAUSS = 1,NGAUSS
      XI = XG(IGAUSS)
      WTX = WG(IGAUSS)
      DO JGAUSS = 1,NGAUSS
        ETA = XG(JGAUSS)
        WTY = WG(JGAUSS)
        WTXY = WTX*WTY
        !Get shape functions
        call getSF(XI, ETA, SF, NODEL, SFDL, NDIME,NORDER)
        !Get Jacobians.
          !Get JCBT = transpose(Jacobian)
          CALL MATMUL(SFDL, NDIME, NODEL, ELCOR, NODEL, NDIME, JCBT,&
          NDIME, NDIME, NDIME, NDIME, NODEL)
          !Invert JCBT
          CALL invJCBT(JCBT, NDIME, JCBTI, DETJ)


        !Update weights
        WTXY = WTXY*DETJ
        
        !Get global partial derivatives, i.e. (d/dx,d/dy) 
        call MATMUL(JCBTI,NDIME,NDIME,SFDL,NDIME,NODEL,SFDG,&
        	    NDIME,NODEL,NDIME,NODEL,NDIME)
        
        !Get (X,Y) from (XI,ETA)
        X = 0.0D0
        Y = 0.0D0
        DO I=1,NODEL
          X = X + SF(I)*ELCOR(I,1)
          Y = Y + SF(I)*ELCOR(I,2)
        ENDDO        
        
        !get Enriched shape functions and derivatives
          !get SFN
          call getSFN(X,Y,SF,NODEL,K_W,NANGL,BETA,SFN,ENR_type)
          !get SFDGN
          call getSFDGN(X, Y,SF,NODEL,SFDG,NDIME,K_W,& 
                    NANGL,BETA,SFDGN,ENR_type)
          !write(*,*)'X = ',x,'Y = ',y
          !Do I = 1,NODEF
          !   write(*,*)'SFDGN (1,',I,') =',SFDGN(1,I)
          !  write(*,*)'SFDGN (2,',I,') =',SFDGN(2,I)
          !endDO

        !Build Mass and Stiffness
        DO I = 1,NODEF
          DO J = 1,NODEF
          ELM(J,I) = ELM(J,I) + SFN(J)*SFN(I)*WTXY
          ELK(J,I) = ELK(J,I) + (SFDGN(1,J)*SFDGN(1,I) + &
                                 SFDGN(2,J)*SFDGN(2,I))*WTXY
          endDO
        endDO

        !Build RHS for DIC
          !Calculate load          
      !     r_mod = sqrt(X*X + Y*Y)
      !     FZ_Phi = CDEXP(DCMPLX(0.0D0,1.0D0)*(K_W*r_mod&
			! 		 - Omega*(0.0D0)))
      !     FZ_Vlcty = CDEXP(DCMPLX(0.0D0,1.0D0)*(K_W*r_mod&
			! 		 - Omega*(0.0D0)))&
      ! *(DCMPLX(0.0D0,-1.0D0)*Omega)
      FZ_Phi = dcmplx(0.0D0,0.0D0)
      FZ_Vlcty = dcmplx(0.0D0,0.0D0)
					 
          !Calculate inner product of load with weight functions
          DO J = 1,NODEF
            PTLD_Phi(J) = PTLD_Phi(J) + FZ_Phi*SFN(J)*WTXY
            PTLD_Vlcty(J) = PTLD_Vlcty(J) + FZ_Vlcty*SFN(J)*WTXY
          endDO
        
      endDO
    endDO
 
  
end SUBROUTINE getEl_DIC_M_S
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
