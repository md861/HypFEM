!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getRHS(TOTELS,NODEL,NORDER,NGAUSS,ofst_ELDAT,ELDAT,&
                 NDIME,COORD,TOTCOORD,NMARK,ELTYPE,FACT,GLDF,NANGL,&
                 IStep,dt,Cp,K_W,Omega,RHS,TOTDOF,&
                 NODEF,BETA,ENR_Type,&
                 c_RKM)

           !Description: This subroutine controls the RHS for  
           		!the time loop. This could be modified
           		!to be used as getSLOPE, wherein, 
           		!it could produce the slope for 
           		!ODE in time, useful for Higher Order Time (HOT) integration.
           !Arguments out
           		!RHS:	The global vector that stores the right hand side 
           			!for a given problem, integrating over domain 
           			!and the boundaries.
IMPLICIT NONE

integer TOTELS,NODEL,NORDER,NDIME,NGAUSS,&
        ofst_ELDAT,TOTCOORD,NMARK,ELTYPE,&
        NANGL,IStep,TOTDOF
integer ELNODDAT(NODEL),ELNDS(NODEL),&
        ELDAT(TOTELS, ofst_ELDAT + NODEL)
double precision COORD(TOTCOORD,NDIME),&
	ELCOR(NODEL,NDIME),dt,T,K_W,Omega,&
	Cp
integer EFACT(NMARK,ELTYPE),&
        FACT(TOTELS,NMARK,ELTYPE),&
        GLDF(TOTELS,NODEL),&
        ELDF(NODEL*NANGL)
integer I,J,IELEM

double complex PTLD(NODEF),RHS(TOTDOF)

!definitions for PUFEM
integer NODEF
  !definitions for ENR_Type
  integer ENR_type		! 1=ENR2, 0=ENR1
  !changes for ENR_Type
  double precision BETA(NANGL-ENR_Type)
  
!definitions for Explicit RKM
double precision c_RKM


T = IStep*dt + (c_RKM*dt)

!Initialize matrices and vectors
RHS = dcmplx(0.0D0,0.0D0)

  !Assemble
  Do IELEM = 1,TOTELS
    !get element co-ordinates
    call GetELNODE(ELNODDAT,NODEL,NORDER)
    ELNDS = ELDAT(IELEM, ofst_ELDAT + ELNODDAT)
    DO I = 1,NODEL
      DO J = 1,NDIME
        ELCOR(I,J) = COORD(ELNDS(I),J)
      endDO
    endDO
    
    !get element edge-factors
    DO I = 1,NMARK
      EFACT(I,1:ELTYPE)=FACT(IELEM,I,1:ELTYPE)
    endDO
        
    !Build local system    
    call getPTLD(IELEM,ELCOR,NODEL,NDIME,NGAUSS,EFACT,&
                   NMARK,ELTYPE,NORDER,K_W,Omega,T,Cp,PTLD,&
                   NANGL,NODEF,BETA,ENR_type)
    
    !Build global system
      !get indices that map from ELememnt Degree of Freedom {ELDF} to GLobal Degree of Freedom {GLDF table}
      call getELDF(IELEM,GLDF,TOTELS,NODEL,ELNODDAT,NORDER,ELDF,NANGL)
      !transfer ELementary matrices/vectors to Global matrices/vectors
      DO I = 1,NODEF
        RHS(ELDF(I)) = RHS(ELDF(I)) + PTLD(I)
      endDO
  endDO
end subroutine getRHS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
