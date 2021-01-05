use module_data_read

implicit none
!Definitions
  integer NDIME,TOTELS,NODEL,TOTCOORD,NODEDGE,NORDER,NMARK,NGAUSS
  integer, allocatable :: ELDAT(:,:)
  integer ofst_ELDAT,ELTYPE
  double precision, allocatable :: COORD(:,:),ELCOR(:,:)
  integer, allocatable :: ELEDGDAT(:,:), ELNODDAT(:), ELNDS(:)
  integer, allocatable :: FACT(:,:,:),EFACT(:,:)  
  integer IStep, NPOINTPLOT, IELEM, I, J,counter
  
  integer TOTNOD
  integer, allocatable ::  GLDF(:,:),MapCoord2Coef(:),ELDF(:)
  integer NANGL
  !definitions for book-keeping 
  character (len=52) char_name
  integer n_record_IStep
  !definitions for DIC,M,S
  integer TOTDOF
  double precision K_W,Omega,PI
  double complex, allocatable :: ELM(:,:),ELK(:,:),&!Elementary Mass & Stiffness
  				 GM(:,:),GS(:,:),&  !Global Mass & Stiffness
  				 InvGM(:,:),&
  				 PTLD_Phi(:),PTLD_Vlcty(:),&
  				 RHS_Phi(:),RHS_Vlcty(:),&
  				 Phi(:),Vlcty(:)
  !definitions for getRHS
  double precision Cp	!Wave speed
  double precision dt
  double complex, allocatable :: RHS(:)
  !definitions for Implicit formulation
    !none
    
  integer n_istep
  
  !definitions for other postprocessing
  double precision error1,error2
  
  !definitions for PUFEM
  integer NODEF 		!NODEF : NODEL*NANGL
  double precision, allocatable :: BETA(:)
    !definitions for ENR_Type
    integer ENR_type		! 1=ENR2, 0=ENR1

  !definitions for Explicit (RKM)
  integer  RKM_Order		!This variable holds the order of RKM
  double precision, allocatable :: a_RKM(:,:),&    !The coefficients for RKM
  				   c_RKM(:),&	   !The nodes for RKM
  				   b_RKM(:)	   !The weights for RKM
  double complex, allocatable :: PrevPhi(:),&                                 	   
  				 PrevVlcty(:),&
  				 temp_vec(:),&     !Holds RKM vector e.g. {vn + dt*sum(a,k)}
  				 Amat(:,:),&	   !A = -c*c*[M]^(-1)[S]
  				 AY(:),&	   !AY= [A]*[RKM vector]
  				 R_tilde(:),&  	   !=[M]^(-1)*(c*c*<BC,w> + <f,w>)|tn + c(i)*dt
  				 k_u(:,:),k_v(:,:) !vectors to store k values for RKM
  				 
  !definitions for Explicit (RKM2)
  double precision alpha
  
  !definitions for Implicit (RKM)
  double complex, allocatable :: InvMat(:,:,:),&
                                 temp_mat(:,:)
  
  
    
  !Get input data  
    call getData(NDIME,TOTELS,NODEL,TOTCOORD,NODEDGE,NORDER,NMARK,&
               ELDAT,ofst_ELDAT,COORD,ELEDGDAT,FACT,ELTYPE,GLDF,&
               TOTNOD,MapCoord2Coef)  
               
  !book-keeping
  char_name = "case_default"
  call system("mkdir "//char_name)	!Used for storing all results
  call system("mkdir PLOTS")		!Used for plots
  open(1002, file='logfile.txt', status="old", position="append", action="write")
  
               
  !Set-up problem related parameters---------------------------------Starts
  PI = 4.0D0*DATAN(1.0D0)
  K_W = 10.0D0*PI
    !Changes/Addition in parameters for PUFEM
    NGAUSS = 6!NGAUSS + CEILING(dsqrt(1.0D0/TOTELS)*K_W*10.0D0/(2.0D0*PI))
    NANGL = 1	!Must be replaced to 1 + n_Beta if ENR_type = 1 (ENR2)
    NODEF = NODEL*NANGL
    !Changes for ENR_Type
    ENR_type = 1! 1=ENR2, 0=ENR1
    DO I = 1,NANGL - ENR_type
      if(I.EQ.1)THEN
        allocate (BETA(NANGL-ENR_type))
      endif
      BETA(I) = 2.0d0*PI*dble(i)/(nangl-ENR_type)
      write(*,*)'BETA(',I,') = ',BETA(I)
    endDO
  TOTDOF = TOTNOD*NANGL 
  write(*,*)'NGAUSS = ',NGAUSS,' This must be (h/Lambda)*10'
  write(1002,*)'NGAUSS = ',NGAUSS  
  write(1002,*)'NANGL = ',NANGL,' is 1 for FEM'   
  write(1002,*)'TOTDOF = ',TOTDOF
  Omega = 0.5D0
  Cp = Omega/K_W
  dt = 0.01 
  n_istep = 6000
  !plot related parameters
  NPOINTPLOT = 30	!# of integration points for plots
  n_record_IStep = CEILING(dble(n_istep)/100.0D0)
    !Parameters for RKM
    RKM_Order = 2
    alpha = 2.0D0/3.0D0 
  !Set-up problem related parameters---------------------------------Stops
  
  !allocate variables--------------------------------------------Starts
    allocate(ELNODDAT(NODEL),ELNDS(NODEL),ELCOR(NODEL,NDIME))
    allocate(EFACT(NMARK,ELTYPE))
    allocate(ELDF(NODEF))
    !allocate for DIC,M,S
    allocate (ELM(NODEF,NODEF),ELK(NODEF,NODEF))
    allocate (PTLD_Phi(NODEF),PTLD_Vlcty(NODEF))
    allocate (GM(TOTDOF,TOTDOF),InvGM(TOTDOF,TOTDOF))
    allocate (GS(TOTDOF,TOTDOF))
    allocate (RHS_Phi(TOTDOF),RHS_Vlcty(TOTDOF))
    allocate (Phi(TOTDOF),Vlcty(TOTDOF))
    !allocate for getRHS
    allocate (RHS(TOTDOF))
    !allocate for Implicit formulation
      !none
    !allocate for Explicit (RKM)
    allocate (Amat(TOTDOF,TOTDOF))
    allocate (k_u(RKM_Order,TOTDOF),k_v(RKM_Order,TOTDOF))
    allocate (a_RKM(RKM_Order,RKM_Order))
    allocate (c_RKM(RKM_Order))
    allocate (b_RKM(RKM_Order))
    allocate (AY(TOTDOF),R_tilde(TOTDOF))
    allocate (temp_vec(TOTDOF))
    allocate (PrevVlcty(TOTDOF))
    allocate (PrevPhi(TOTDOF))
    !allocate for Implicit (RKM)
    allocate (InvMat(RKM_Order,TOTDOF,TOTDOF))
    allocate (temp_mat(TOTDOF,TOTDOF))
  !allocate variables--------------------------------------------Stops
  
  !~~~~~~~~~~~~~~~~~~~~DIC,M,S~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!Starts
  write(*,*)'~~~~~~~~~~~~~~~~~~~'
  write(*,*)'Processing DIC,M,S'
  write(*,*)'~~~~~~~~~~~~~~~~~~~'  
  !Initialize matrices and vectors
  GM = dcmplx(0.0D0,0.0D0)
  GS = dcmplx(0.0D0,0.0D0)
  RHS_Phi = dcmplx(0.0D0,0.0D0)
  RHS_Vlcty = dcmplx(0.0D0,0.0D0)
  Phi = dcmplx(0.0D0,0.0D0)
  Vlcty = dcmplx(0.0D0,0.0D0)
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
        
    !Build local system    
    call getEl_DIC_M_S(IELEM,ELCOR,NODEL,NDIME,&
		   NGAUSS,NORDER,K_W,Omega,&
		   ELM,ELK,PTLD_Phi,PTLD_Vlcty,&
		   NANGL,NODEF,BETA,ENR_Type)
		       
    !Build global system
      !get indices that map from ELememnt Degree of Freedom {ELDF} to GLobal Degree of Freedom {GLDF table}
      call getELDF(IELEM,GLDF,TOTELS,NODEL,ELNODDAT,NORDER,ELDF,NANGL)
      !transfer ELementary matrices/vectors to Global matrices/vectors
      DO I = 1,NODEF
        DO J = 1,NODEF
          GM(ELDF(I),ELDF(J)) = GM(ELDF(I),ELDF(J)) + ELM(I,J)
          GS(ELDF(I),ELDF(J)) = GS(ELDF(I),ELDF(J)) + ELK(I,J)
        endDO
        RHS_Phi(ELDF(I)) = RHS_Phi(ELDF(I)) + PTLD_Phi(I)
        RHS_Vlcty(ELDF(I)) = RHS_Vlcty(ELDF(I)) + PTLD_Vlcty(I)
      endDO    
  endDO  
  !Solve for Decomposition of Phi0,Vlcty0
    !Invert Mass
    call inv(GM,TOTDOF,InvGM,TOTDOF)
    !Multiply inverted Mass with RHS to obtain Phi0,Vlcty0
    CALL MATMULCPLX(InvGM,TOTDOF,TOTDOF,RHS_Phi,TOTDOF,1,Phi,TOTDOF,1,TOTDOF,1,TOTDOF)
    CALL MATMULCPLX(InvGM,TOTDOF,TOTDOF,RHS_Vlcty,TOTDOF,1,Vlcty,TOTDOF,1,TOTDOF,1,TOTDOF)
    !get Phi(-1) from Phi0,Vlcty0: Phi(-1) = Phi0 - dt*Vlcty0
    !PrevPhi = Phi - (dt*Vlcty)	!This is useful for implicit formulation

  !~~~~~~~~~~~~~~~~~~~~DIC,M,S~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!Stops
  write(*,*)'~~~~~~~~~~~~~~~~~~~'
  write(*,*)'Matrices for explicit-wave'
  write(*,*)'~~~~~~~~~~~~~~~~~~~'  
  
  !Update Stiffness for wave problem
  GS = (Cp*Cp)*GS
  
  
  !Compute A = -c*c*[M]^(-1)[S] == -[invGM]*[GS]
  call MATMULCPLX(invGM,TOTDOF,TOTDOF,&
	    (dcmplx(-1.0D0,0.0D0)*GS),TOTDOF,TOTDOF,&
            Amat,TOTDOF,TOTDOF,TOTDOF,TOTDOF,TOTDOF)

  !Initialize 
    !for Implicit RKM-1
    !a_RKM(1,1) = 1
    !c_RKM(1) = 1
    !b_RKM(1) = 1
    !for Implicit RKM-2 {Guass-Legendre}
    a_RKM(1,1) = 0.25
    a_RKM(1,2) = 0.25 - (dsqrt(3.0D0)/6.0D0)
    a_RKM(2,1) = 0.25 + (dsqrt(3.0D0)/6.0D0)
    a_RKM(2,2) = 0.25
    c_RKM(1) = 0.5 - (dsqrt(3.0D0)/6.0D0)
    c_RKM(2) = 0.5 + (dsqrt(3.0D0)/6.0D0)
    b_RKM(1) = 0.5
    b_RKM(2) = 0.5
    !compute inv(1-(alpha^2)*A)
    Do I = 1,RKM_Order
      temp_mat = (-1.0D0)*(dt*dt*a_RKM(I,I)*a_RKM(I,I))*Amat
      Do J = 1,TOTDOF
      temp_mat(J,J) = temp_mat(J,J) + 1.0D0
      endDO
      call inv(temp_mat,TOTDOF,InvMat(I,:,:),TOTDOF)
    endDO
  !~~~~~~~~~~~~~~~~~~~~TimeLoop~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!Starts
  Do Istep = 1,n_istep
    write(*,*)'~~~~~~~~~~~~~~~~~~~'
    write(*,*)'Time-step ', IStep
    write(*,*)'~~~~~~~~~~~~~~~~~~~'     
    !Save history of Phi and Vlcty
    PrevPhi = Phi
    PrevVlcty  = Vlcty
    !Compute k vectors for RKM
      Do I = 1,RKM_Order
        !get explicit k vectors
          !get k_u
            !update temp_vec   
            temp_vec = dcmplx(0.0D0,0.0D0)
            Do J = 1,I-1
              temp_vec = temp_vec + a_RKM(I,J)*k_v(J,:)
            endDO         
            k_u(I,:) = PrevVlcty + dt*(temp_vec)
        
          !get k_v
            !update temp_vec
            temp_vec = dcmplx(0.0D0,0.0D0)
            Do J = 1,I-1
              temp_vec = temp_vec + a_RKM(I,J)*k_u(J,:)
            endDO
            temp_vec = PrevPhi + dt*(temp_vec)
            !calculate (c*c*<BC,w> + c*c*<f,w>)  at t = tn + c(i)*dt  
            call getRHS(TOTELS,NODEL,NORDER,NGAUSS,ofst_ELDAT,ELDAT,&
                     NDIME,COORD,TOTCOORD,NMARK,ELTYPE,FACT,GLDF,NANGL,&
                     IStep,dt,Cp,K_W,Omega,RHS,TOTDOF,&
                     NODEF,BETA,ENR_Type,&
                     c_RKM(I))
            !calculate [M]^(-1)*(c*c*<BC,w> + c*c*<f,w>) == R_tilde
            call MATMULCPLX(invGM,TOTDOF,TOTDOF,&
                        RHS,TOTDOF,1,&
                        R_tilde,TOTDOF,1,TOTDOF,1,TOTDOF)
            !calculate [A]*{Un} == AY
            call MATMULCPLX(Amat,TOTDOF,TOTDOF,&
                        temp_vec,TOTDOF,1,&
                        AY,TOTDOF,1,TOTDOF,1,TOTDOF) 
          
            k_v(I,:) = AY + R_tilde  
        !get implicit k vectors from explicit k vectors
            !calculate [A]k_u == temp_vec
            call MATMULCPLX(Amat,TOTDOF,TOTDOF,&
                        k_u(I,:),TOTDOF,1,&
                        temp_vec,TOTDOF,1,TOTDOF,1,TOTDOF)         
            
            !get k_v = [inv(1-(alpha^2)*A)]*{k_v - (alpha)[A]k_u}
            temp_vec = k_v(I,:) + (dt*a_RKM(I,I))*temp_vec   
            call MATMULCPLX(InvMat(I,:,:),TOTDOF,TOTDOF,&
                        temp_vec,TOTDOF,1,&
                        k_v(I,:),TOTDOF,1,TOTDOF,1,TOTDOF) 
            !get k_u = k_u - (alpha)k_v
            k_u(I,:) = k_u(I,:) + (dt*a_RKM(I,I))*k_v(I,:)
      endDO     
       
        
    !Comnpute solution by adding up k vectors from RKM        
      !compute Phi
        temp_vec = dcmplx(0.0D0,0.0D0)
        Do I = 1,RKM_Order
          temp_vec = temp_vec + b_RKM(I)*k_u(I,:)
        endDO
        Phi = PrevPhi + dt*(temp_vec)
      !compute Vlcty
        temp_vec = dcmplx(0.0D0,0.0D0)
        Do I = 1,RKM_Order
          temp_vec = temp_vec + b_RKM(I)*k_v(I,:)
        endDO
        Vlcty = PrevVlcty + dt*(temp_vec)              


  
  
    !Post Processing      
      ! !Errors
      ! call getLn_norm(IStep,TOTELS,ELDAT,COORD,ofst_ELDAT,NODEL,&
      !              TOTCOORD,NDIME,NORDER,NANGL,TOTDOF,Phi,GLDF,&
      !              ELCOR,NGAUSS,&
      !              K_W,Omega,dt,Cp,error1,error2,&
      !              NODEF,BETA,ENR_type)
      ! write(*,*)'L2 error % =',error2*100
      ! write(*,*)'L1 error % =',error1*100
      !Plots            
      if(((IStep .EQ. 1).OR.(mod(IStep,n_record_IStep).EQ. 0)).OR.(IStep .EQ. n_istep))Then
       call plotPara(IStep,NPOINTPLOT,TOTELS,ELDAT,COORD,ofst_ELDAT,NODEL,&
                  TOTCOORD,NDIME,NORDER,NANGL,TOTDOF,Phi,Vlcty,GLDF,&
                  K_W,Omega,dt,&
                  NODEF,BETA,ENR_type)
      endif                   
      !Exports
      !open(1020, file = 'realMASS')
      !open(1021, file = 'imagMASS')
      !Do I = 1,TOTDOF
      !  DO J = 1,TOTDOF
      !    write(1020,1102,advance='no')dreal(GM(I,J))
      !	write(1021,1102,advance='no')dimag(GM(I,J))
      !	1102 format(f16.7)
      !  endDO
      !  write(1020,*)
      !  write(1021,*)
      !endDO
      !close(1020)
      !close(1021)  
  endDO !{Time}
  !~~~~~~~~~~~~~~~~~~~~TimeLoop~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!Stops
  
  !book-keeping	
  call system("cp dat "//char_name)
  call system("cp -r PLOTS "//char_name)
  call system("cp error_1_data "//char_name)
  call system("cp error_2_data "//char_name)
  call system("cp logfile.txt "//char_name)
  !call system("rm dat")
  call system("rm -r PLOTS")
  call system("rm error_1_data")
  call system("rm error_2_data")
  call system("rm logfile.txt")
                   
END
