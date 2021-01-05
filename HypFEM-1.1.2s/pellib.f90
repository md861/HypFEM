!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE getPTLD(IELEM,ELCOR,NODEL,NDIME,NGAUSS,EFACT,&
                   NMARK,ELTYPE,NORDER,K_W,Omega,T,Cp,PTLD,&
                   NANGL,NODEF,BETA,ENR_Type)

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


integer NMARK,ELTYPE,EFACT(NMARK,ELTYPE)
double precision n(NDIME)
integer EDG_TYPE

!definitions for getRHS
double precision r_mod,K_W,Omega,T,Cp
double complex FZ,GZ,PTLD(NODEF),&
               dUdx,dUdy,Uxy
!definitions for PUFEM
integer NANGL,NODEF
  !definitions for ENR_Type
  integer ENR_type		! 1=ENR2, 0=ENR1
  !changes for ENR_Type
  double precision BETA(NANGL-ENR_type)
double complex SFN(NODEF)

!difinitions for sources
double precision ea,eb,ox,oy,sprd,limx,limy

  !Initialize elementary vectors
  PTLD = dcmplx(0.0D0,0.0D0)
  
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

        
        !Get (X,Y) from (XI,ETA)
        X = 0.0D0
        Y = 0.0D0
        DO I=1,NODEL
          X = X + SF(I)*ELCOR(I,1)
          Y = Y + SF(I)*ELCOR(I,2)
        ENDDO        
        
        !get Enriched shape functions
          !get SFN
          call getSFN(X,Y,SF,NODEL,K_W,NANGL,BETA,SFN,ENR_type)
        
        !Compute load
        FZ = dcmplx(0.0D0,0.0D0)
        ! !parameters for bench_geo
        !   ea = 0.09
        !   eb = 0.14
        !   ox = 0.56
        !   oy = 0.67
        !parameters for hex_screw
          ea = 0.1
          eb = 0.1
          ox = 0.9
          oy = 0.9
          sprd = 10
          limx = ((y-oy)/(eb))*((y-oy)/(eb))
          limx = ea*ea*(1 - limx)
          limy = ((x-ox)/(ea))*((x-ox)/(ea))
          limy = eb*eb*(1-limy)
          if(( (x-ox)*(x-ox) .lt. limx).and.( (y-oy)*(y-oy) .lt. limy))then
            r_mod = sqrt((x-ox)*(x-ox) + (y-oy)*(y-oy))
            FZ = exp(-r_mod*r_mod*sprd*sprd)*dsin(Omega*T)
          endif
        
        !Compute inner product of load with weight functions
        DO J = 1,NODEF
            PTLD(J) = PTLD(J) + Cp*Cp*FZ*SFN(J)*WTXY            
        endDO
        
      endDO
    endDO
      
    
  
  ! !Integrate over edges.
  !   !edge 1
  !   if(EFACT(1,1).EQ. 1)THEN   
  !   EDG_TYPE = 1
  !   !write(*,*)'-----------------'
  !   !write(*,*) 'ELEMENT = ',IELEM,' EDGE = 1'
  !     Do IGAUSS = 1,NGAUSS
  !       XI = XG(IGAUSS)
  !       WTX = WG(IGAUSS)
  !       ETA = -1.0
  !       WTY = 1
  !       WTXY = WTX*WTY
  !       !Get shape functions
  !       call getSF(XI, ETA, SF, NODEL, SFDL, NDIME,NORDER)
  !       !Get Jacobians.
  !         !Get JCBT = transpose(Jacobian)
  !         CALL MATMUL(SFDL, NDIME, NODEL, ELCOR, NODEL, NDIME, JCBT,&
  !         NDIME, NDIME, NDIME, NDIME, NODEL)      
        
  !       !Update weights
  !         !get Determinant for line integral = dsqrt((dx/dXI)^2 + (dy/dXI)^2)
  !         DETJ = dsqrt((JCBT(1,1)**2) + (JCBT(1,2)**2))
  !         !Update weight
  !         WTXY = WTXY*DETJ 
  !         !write(*,*) 'DETJ = ',DETJ
        
  !       !Get (X,Y) from (XI,ETA)
  !       X = 0.0D0
  !       Y = 0.0D0
  !       DO I=1,NODEL
  !         X = X + SF(I)*ELCOR(I,1)
  !         Y = Y + SF(I)*ELCOR(I,2)
  !       ENDDO 
        
  !       !get Enriched shape functions
  !         !get SFN
  !         call getSFN(X,Y,SF,NODEL,K_W,NANGL,BETA,SFN,ENR_type)
        
  !       !Get normal to edge
  !       call getNormVec(EDG_TYPE,JCBT,NDIME,NODEL,X,Y,ELCOR,n,NORDER)
  !       !write(*,*)'n = ',n  
        
  !       !Compute load
  !       r_mod = sqrt(X*X + Y*Y)
  !       GZ = dcmplx(0.0D0,0.0D0)
  !       Uxy  = CDEXP(DCMPLX(0.0D0,1)*(K_W*r_mod&
	! 		-Omega*T))
  !       dUdx = DCMPLX(0.0D0, K_W)*(X/r_mod)*Uxy
  !       dUdy = DCMPLX(0.0D0, K_W)*(Y/r_mod)*Uxy
  !       GZ = dUdx*n(1) + dUdy*n(2)
        
  !       !Compute inner product of load with weight functions
  !       DO J = 1,NODEF
  !           PTLD(J) = PTLD(J) + Cp*Cp*GZ*SFN(J)*WTXY            
  !       endDO
  !     endDO
  !     !write(*,*)'-----------------'
  !   endif
        
        
  !   !edge 2
  !   if(EFACT(1,2).EQ. 1)THEN  
  !   EDG_TYPE = 0 
  !   !write(*,*)'-----------------'
  !   !write(*,*) 'ELEMENT = ',IELEM,' EDGE = 2'
  !     Do IGAUSS = 1,NGAUSS
  !       ETA = XG(IGAUSS)
  !       WTY = WG(IGAUSS)
  !       XI = 1.0
  !       WTX = 1
  !       WTXY = WTX*WTY
  !       !Get shape functions
  !       call getSF(XI, ETA, SF, NODEL, SFDL, NDIME,NORDER)
  !       !Get Jacobians.
  !         !Get JCBT = transpose(Jacobian)
  !         CALL MATMUL(SFDL, NDIME, NODEL, ELCOR, NODEL, NDIME, JCBT,&
  !         NDIME, NDIME, NDIME, NDIME, NODEL)          
        
        
  !       !Update weights
  !         !get Determinant for line integral = dsqrt((dx/dETA)^2 + (dy/dETA)^2)
  !         DETJ = dsqrt((JCBT(2,1)**2) + (JCBT(2,2)**2))
  !         !Update weight
  !         WTXY = WTXY*DETJ  
  !         !write(*,*) 'DETJ = ',DETJ

        
  !       !Get (X,Y) from (XI,ETA)
  !       X = 0.0D0
  !       Y = 0.0D0
  !       DO I=1,NODEL
  !         X = X + SF(I)*ELCOR(I,1)
  !         Y = Y + SF(I)*ELCOR(I,2)
  !       ENDDO 
        
  !       !get Enriched shape functions
  !         !get SFN
  !         call getSFN(X,Y,SF,NODEL,K_W,NANGL,BETA,SFN,ENR_type)
        
  !       !Get normal to edge
  !       call getNormVec(EDG_TYPE,JCBT,NDIME,NODEL,X,Y,ELCOR,n,NORDER)
  !       !write(*,*)'n = ',n  
        
  !       !Compute load
  !       r_mod = sqrt(X*X + Y*Y)
  !       GZ = dcmplx(0.0D0,0.0D0)
  !       Uxy  = CDEXP(DCMPLX(0.0D0,1)*(K_W*r_mod&
	! 		-Omega*T))
  !       dUdx = DCMPLX(0.0D0, K_W)*(X/r_mod)*Uxy
  !       dUdy = DCMPLX(0.0D0, K_W)*(Y/r_mod)*Uxy
  !       GZ = dUdx*n(1) + dUdy*n(2)
        
  !       !Compute inner product of load with weight functions
  !       DO J = 1,NODEF
  !           PTLD(J) = PTLD(J) + Cp*Cp*GZ*SFN(J)*WTXY            
  !       endDO
  !     endDO
  !   !write(*,*)'-----------------'
  !   endif
    
    
  !   !edge 3
  !   if(EFACT(1,3).EQ. 1)THEN   
  !   EDG_TYPE = 1
  !   !write(*,*)'-----------------'
  !   !write(*,*) 'ELEMENT = ',IELEM,' EDGE = 3'
  !     Do IGAUSS = 1,NGAUSS
  !       XI = XG(IGAUSS)
  !       WTX = WG(IGAUSS)
  !       ETA = 1.0
  !       WTY = 1
  !       WTXY = WTX*WTY
  !       !Get shape functions
  !       call getSF(XI, ETA, SF, NODEL, SFDL, NDIME,NORDER)
        
  !       !Get Jacobians.
  !         !Get JCBT = transpose(Jacobian)
  !         CALL MATMUL(SFDL, NDIME, NODEL, ELCOR, NODEL, NDIME, JCBT,&
  !         NDIME, NDIME, NDIME, NDIME, NODEL)
        
  !       !Update weights
  !         !get Determinant for line integral = dsqrt((dx/dXI)^2 + (dy/dXI)^2)
  !         DETJ = dsqrt((JCBT(1,1)**2) + (JCBT(1,2)**2))
  !         !Update weight
  !         WTXY = WTXY*DETJ        
  !         !write(*,*) 'DETJ = ',DETJ
        
  !       !Get (X,Y) from (XI,ETA)
  !       X = 0.0D0
  !       Y = 0.0D0
  !       DO I=1,NODEL
  !         X = X + SF(I)*ELCOR(I,1)
  !         Y = Y + SF(I)*ELCOR(I,2)
  !       ENDDO 
        
  !       !get Enriched shape functions
  !         !get SFN
  !         call getSFN(X,Y,SF,NODEL,K_W,NANGL,BETA,SFN,ENR_type)
        
  !       !Get normal to edge  
  !       call getNormVec(EDG_TYPE,JCBT,NDIME,NODEL,X,Y,ELCOR,n,NORDER)
  !       !write(*,*)'n = ',n  
          
  !       !Compute load
  !       r_mod = sqrt(X*X + Y*Y)
  !       GZ = dcmplx(0.0D0,0.0D0)
  !       Uxy  = CDEXP(DCMPLX(0.0D0,1)*(K_W*r_mod&
	! 		-Omega*T))
  !       dUdx = DCMPLX(0.0D0, K_W)*(X/r_mod)*Uxy
  !       dUdy = DCMPLX(0.0D0, K_W)*(Y/r_mod)*Uxy
  !       GZ = dUdx*n(1) + dUdy*n(2)
        
  !       !Compute inner product of load with weight functions
  !       DO J = 1,NODEF
  !           PTLD(J) = PTLD(J) + Cp*Cp*GZ*SFN(J)*WTXY            
  !       endDO
  !     endDO
  !   !write(*,*)'-----------------'
  !   endif
    
  !   !edge 4
  !   if(EFACT(1,4).EQ. 1)THEN   
  !   EDG_TYPE = 0
  !   !write(*,*)'-----------------'
  !   !write(*,*) 'ELEMENT = ',IELEM,' EDGE = 4'
  !     Do IGAUSS = 1,NGAUSS
  !       ETA = XG(IGAUSS)
  !       WTY = WG(IGAUSS)
  !       XI = -1.0
  !       WTX = 1
  !       WTXY = WTX*WTY
  !       !Get shape functions
  !       call getSF(XI, ETA, SF, NODEL, SFDL, NDIME,NORDER)
  !       !Get Jacobians.
  !         !Get JCBT = transpose(Jacobian)
  !         CALL MATMUL(SFDL, NDIME, NODEL, ELCOR, NODEL, NDIME, JCBT,&
  !         NDIME, NDIME, NDIME, NDIME, NODEL)
          
        
  !       !Update weights
  !         !get Determinant for line integral = dsqrt((dx/dETA)^2 + (dy/dETA)^2)
  !         DETJ = dsqrt((JCBT(2,1)**2) + (JCBT(2,2)**2))
  !         !Update weight
  !         WTXY = WTXY*DETJ  
  !         !write(*,*) 'DETJ = ',DETJ
          
        
  !       !Get (X,Y) from (XI,ETA)
  !       X = 0.0D0
  !       Y = 0.0D0
  !       DO I=1,NODEL
  !         X = X + SF(I)*ELCOR(I,1)
  !         Y = Y + SF(I)*ELCOR(I,2)
  !       ENDDO 
        
  !       !get Enriched shape functions
  !         !get SFN
  !         call getSFN(X,Y,SF,NODEL,K_W,NANGL,BETA,SFN,ENR_type)

  !       !Get normal to edge
  !       call getNormVec(EDG_TYPE,JCBT,NDIME,NODEL,X,Y,ELCOR,n,NORDER)
  !       !write(*,*)'n = ',n  
        
  !       !Compute load
  !       r_mod = sqrt(X*X + Y*Y)
  !       GZ = dcmplx(0.0D0,0.0D0)
  !       Uxy  = CDEXP(DCMPLX(0.0D0,1)*(K_W*r_mod&
	! 		-Omega*T))
  !       dUdx = DCMPLX(0.0D0, K_W)*(X/r_mod)*Uxy
  !       dUdy = DCMPLX(0.0D0, K_W)*(Y/r_mod)*Uxy
  !       GZ = dUdx*n(1) + dUdy*n(2)
        
  !       !Compute inner product of load with weight functions
  !       DO J = 1,NODEF
  !           PTLD(J) = PTLD(J) + Cp*Cp*GZ*SFN(J)*WTXY            
  !       endDO
  !     endDO
  !   !write(*,*)'-----------------'
  !   endif
 
end SUBROUTINE getPTLD
