!#################################################################################################
!*************************************************************************************************
!-------------------------------------------------------------------------------------------------
!                                                 Copyright (C) 2013 Lei HAN & Zhanpeng ZHUANG
!                                                 PROGRAM NAME : masnum_circ_model            
!                                                 PRESENT VERSION : 2013-09-05                
!                                                                                                 
! --- GENERAL DESCRIPTION : 
!     The MASNUM circulation model single-thread version              
! --- DEPEND ON: Netcdf library
! --- Model structure:
!     Part 1: module of array definitions
!             -- contains array-initializing subroutines
!     Part 2: main program
!             -- contains subroutines called by the main program
!
!-------------------------------------------------------------------------------------------------
!*************************************************************************************************
!-------------------------------------------------------------------------------------------------

!*************************************************************************************************
!-------------------------------------------------------------------------------------------------
!
! MASNUM model Part 1: Module of array definitions
! --- Description:
!
!-------------------------------------------------------------------------------------------------
!*************************************************************************************************

module def_common_mod
  
  use time_mod
  
  implicit none
  
  public init_common_array, read_parameters
  public get_2D_field, get_3D_field
  public bilinear_4points_HL,  trilinear_8points_HL
  public test_interp_field_2d, test_interp_field_3d
  
  public
  
  !---------------------------------------------------------------------------------------------
  ! Array size setup:
  !---------------------------------------------------------------------------------------------
  integer,parameter :: precision = 4          ! single=4, double=8
  integer :: i,j,k,m,n,im,jm,ks,kb,l
  integer :: imm1,jmm1,ksm1,kbm1
  
  integer :: n_month

  !---------------------------------------------------------------------------------------------
  ! Constants:
  !---------------------------------------------------------------------------------------------
  real(kind=precision) :: Ra, PI, GRAV, UMOL
  real(kind=precision) :: HORCON
  real(kind=precision) :: RHO_REF
  real(kind=precision) :: Boltzmann
  real(kind=precision) :: TPRNI
  real(kind=precision) :: SMALL
  real(kind=precision) :: ZERO, ONE
  
  real(kind=precision) :: Land_Depth, Min_Depth, Max_Depth
  integer Num_topo_smooth
  integer Use_restart_date
  integer sst_bound
  
  !---------------------------------------------------------------------------------------------
  ! Spherical coordinate 
  !---------------------------------------------------------------------------------------------
  real(kind=precision) :: lon0, lat0, inv_dlon, inv_dlat
  real(kind=precision), allocatable, target :: vlon(:), vlon_u(:)
  real(kind=precision), allocatable, target :: vlat(:), vlat_v(:)

  !---------------------------------------------------------------------------------------------
  ! Time-related variables:
  !---------------------------------------------------------------------------------------------
  real(kind=precision) :: dte, dti, isplit
  integer :: iend, iint, iext, iint0
  integer :: Nint, step_per_day
  real(kind=precision) :: period, period_obc, period_tide
  real(kind=precision) :: ramp, ramp_obc, ramp_tide
  
  integer :: ymdhms(6), ymdhms_prev(6)
  integer :: ymdhms_start(6), ymdhms_finish(6)
  character(len=8) date_char, date_char_prev
  real(kind=8) :: days_now, days_begin, model_days, yearly_day

  !---------------------------------------------------------------------------------------------
  ! Controlling parameters
  !---------------------------------------------------------------------------------------------
  integer :: Baropg_method
  integer :: Steps_upd_advdiff
  integer :: Mode
  integer :: Cold_start
  integer :: IF_save_restart
  integer :: Load_Bv
  integer :: Load_OB_tide
  integer :: Regional_run
  integer :: ChangJiang_diluting
  character(len=50) ext_data_dir

  integer :: wav_freq          ! frequency of getting the wave-induced mixing (hours)
  
  !---------------------------------------------------------------------------------------------
  ! Vertical Coordinate
  !---------------------------------------------------------------------------------------------
  real(kind=precision) :: kl1, kl2
  real(kind=precision), allocatable :: z(:), zz(:), zt(:), dz(:), dzz(:)

  !---------------------------------------------------------------------------------------------
  ! Coriolis parameter
  !---------------------------------------------------------------------------------------------
  real(kind=precision), allocatable :: Cor_t(:,:), Cor_u(:,:), Cor_v(:,:)

  !---------------------------------------------------------------------------------------------
  ! 2D Common variables
  !---------------------------------------------------------------------------------------------
  real(kind=precision), allocatable, target :: H_c(:,:), H_u(:,:), H_v(:,:)

  !---------------------------------------------------------------------------------------------
  ! External mode
  !---------------------------------------------------------------------------------------------
  real(kind=precision), allocatable :: etb_c(:,:)
  real(kind=precision), allocatable, target :: etf_c(:,:)
  real(kind=precision), allocatable :: et_mean(:,:)

  !---------------------------------------------------------------------------------------------
  ! Internal mode
  !---------------------------------------------------------------------------------------------
  real(kind=precision), allocatable,target :: zetf_c(:,:), zetb_c(:,:)             
  real(kind=precision), allocatable, target :: hb_inc(:,:), hb_inu(:,:), hb_inv(:,:)
  real(kind=precision), allocatable, target :: hm_inc(:,:), hm_inu(:,:), hm_inv(:,:)
  real(kind=precision), allocatable, target :: hf_inc(:,:), hf_inu(:,:), hf_inv(:,:)
  real(kind=precision), allocatable :: utf(:,:), vtf(:,:)                   
  
  real(kind=precision), allocatable, target :: fsm(:,:), dum(:,:), dvm(:,:), fsm3d(:,:,:)                   
  real(kind=precision), allocatable, target :: dum01(:,:), dvm01(:,:)                              
  real(kind=precision), allocatable :: dx_t(:,:), dy_t(:,:), dx_u(:,:), dy_u(:,:)          
  real(kind=precision), allocatable :: dx_v(:,:), dy_v(:,:), dx_u_down(:,:), dy_v_left(:,:)
  real(kind=precision), allocatable :: art(:,:), aru(:,:), arv(:,:)                        
  real(kind=precision), allocatable :: int_vm_u(:,:), int_vm_v(:,:)                        
  real(kind=precision), allocatable :: Am_vm(:,:)                               
  real(kind=precision), allocatable, target :: uab(:,:), ua(:,:)                         
  real(kind=precision), allocatable, target :: vab(:,:), va(:,:) 
  real(kind=precision), allocatable, target :: uaf(:,:), vaf(:,:)                         
  real(kind=precision), allocatable :: adv_ua(:,:), adv_va(:,:)                            
  real(kind=precision), allocatable :: diff_ua(:,:), diff_va(:,:)                          
  
  real(kind=precision) :: cbcmin, cbcmax, z0b, kappa
  real(kind=precision), allocatable :: cbc(:,:)                                          
  real(kind=precision), allocatable, target :: wusurf(:,:), wvsurf(:,:), wtsurf(:,:), wssurf(:,:)
  real(kind=precision), allocatable :: wubot(:,:) , wvbot(:,:)                           
  real(kind=precision), allocatable :: fric_ub(:,:), fric_vb(:,:)                        
  
  real(kind=precision), allocatable :: psi(:,:), ml_dep(:,:)

  !---------------------------------------------------------------------------------------------
  ! Spatial moothing factor
  !---------------------------------------------------------------------------------------------
  real(kind=precision) :: smooth_factor
    
  !=============================================================================================
  ! 3D Common variables
  !=============================================================================================
  !---------------------------------------------------------------------------------------------
  ! Prognostic variables
  !---------------------------------------------------------------------------------------------
  real(kind=precision), allocatable, target :: ub(:,:,:), uf(:,:,:)
  real(kind=precision), allocatable, target :: vb(:,:,:), vf(:,:,:) 
  real(kind=precision), allocatable :: wf(:,:,:)           
  real(kind=precision), allocatable, target :: tb(:,:,:), tf(:,:,:)
  real(kind=precision), allocatable, target :: sb(:,:,:), sf(:,:,:)

  real(kind=precision), allocatable :: um(:,:,:), vm(:,:,:)
  real(kind=precision), allocatable :: tm(:,:,:), sm(:,:,:) 
    
  !---------------------------------------------------------------------------------------------
  ! ADV terms
  !---------------------------------------------------------------------------------------------
  real(kind=precision), allocatable :: adv_diff_u(:,:,:), adv_diff_v(:,:,:), adv_t(:,:,:)
    
  !---------------------------------------------------------------------------------------------
  ! Diffusive terms
  !---------------------------------------------------------------------------------------------
  real(kind=precision), allocatable :: diff_t(:,:,:)
  real(kind=precision), allocatable :: Am(:,:,:)                       
  real(kind=precision), allocatable :: Km(:,:,:), Kh(:,:,:), Kq(:,:,:)
  real(kind=precision), allocatable, target :: bv_var(:,:,:), bti_var(:,:,:)
  real(kind=precision), allocatable, target :: bsm1(:,:,:), bsm2(:,:,:)
  real(kind=precision), allocatable, target :: bsmt(:,:,:), bsms(:,:,:)
  
  !---------------------------------------------------------------------------------------------
  ! Coriolis terms
  !---------------------------------------------------------------------------------------------
  real(kind=precision), allocatable :: fcor_u(:,:,:), fcor_v(:,:,:)
  
  !---------------------------------------------------------------------------------------------
  ! Pressure terms
  !---------------------------------------------------------------------------------------------
  real(kind=precision), allocatable :: epg_int_x(:,:), epg_int_y(:,:)
  real(kind=precision), allocatable :: ipg_x(:,:,:),  ipg_y(:,:,:)   

  !---------------------------------------------------------------------------------------------
  ! Thermodynamic variables
  !---------------------------------------------------------------------------------------------
  real(kind=precision), allocatable, target :: rho(:,:,:), rclim(:,:,:)
  real(kind=precision), allocatable, target :: tclim(:,:,:), sclim(:,:,:)
  real(kind=precision), allocatable :: N2(:,:,:)

  !---------------------------------------------------------------------------------------------
  ! POM Profq and Advq variabes
  !---------------------------------------------------------------------------------------------
  real(kind=precision), allocatable, target :: q2b(:,:,:), q2lb(:,:,:)
  real(kind=precision), allocatable, target :: q2f(:,:,:), q2lf(:,:,:)

  !---------------------------------------------------------------------------------------------
  ! Variables for debug use
  !---------------------------------------------------------------------------------------------
  integer, allocatable :: ix(:),jx(:),kx(:)
  
  !---------------------------------------------------------------------------------------------
  ! Surface flux variables 
  !---------------------------------------------------------------------------------------------
  real(kind=precision), allocatable, target :: netheat(:,:)
  real(kind=precision), allocatable, target :: dqdsst(:,:), sst_dq(:,:), sst_OI(:,:)

  !=============================================================================================
  ! Open boundary variables for regional case 
  !=============================================================================================
  !---------------------------------------------------------------------------------------------
  ! Open boundary variables
  !---------------------------------------------------------------------------------------------
  integer :: nmax, numboundary
  integer :: ib,jb
  !!!integer, parameter :: nmax = im*2+jm
  integer, allocatable :: ij_obc(:,:)
  integer, allocatable :: bctype(:)
  real(kind=precision) :: rfe, rfw, rfs, rfn   ! Special variables in cntide case (irregular OB)
  
  !---------------------------------------------------------------------------------------------
  ! Open boundary water level input                        
  !---------------------------------------------------------------------------------------------
  integer, parameter :: numtide = 1 
  real(kind=precision) ::  Omega_tide(4)
    
  !---------------------------------------------------------------------------------------------
  ! OBC variables
  !---------------------------------------------------------------------------------------------
  real(kind=precision), allocatable, target :: et_obc(:), ua_obc(:), va_obc(:)  
  real(kind=precision), allocatable, target :: ub_obc(:,:), vb_obc(:,:)       
  real(kind=precision), allocatable, target :: tb_obc(:,:), sb_obc(:,:)       
  
  !---------------------------------------------------------------------------------------------
  ! Open boundary tide harmonic constants
  !---------------------------------------------------------------------------------------------
  real(kind=precision), allocatable :: Amplitude(:,:), Phase(:,:)
  
  !---------------------------------------------------------------------------------------------
  ! Bv variables
  !---------------------------------------------------------------------------------------------
  integer, parameter :: i0_bv=109, i1_bv=199, j0_bv=139, j1_bv=247
  integer, parameter :: ks_bv=16
  integer, parameter :: im_bv = i1_bv-i0_bv+1
  integer, parameter :: jm_bv = j1_bv-j0_bv+1
  integer, allocatable :: Bv_all(:,:,:,:)
  real(kind=precision), allocatable :: Bv_CN(:,:,:)
  
  !---------------------------------------------------------------------------------------------
  ! Changjiang flux
  !---------------------------------------------------------------------------------------------
  real(kind=precision) :: chj_q(12), time_chj_q
  integer, parameter :: chj_i=136, chj_j=189
  
  !-----------------------------------------------------------------------------------------------
  ! Model execution time statistic
  !-----------------------------------------------------------------------------------------------
  real(kind=precision) :: cpu_start, cpu_stop
  real(kind=precision) :: cpu_realtime

  !-----------------------------------------------------------------------------------------------
  real(kind=precision) :: total_KE             ! Energy variables
  integer :: LL                             ! OBC
  real(kind=precision), allocatable, target :: courant(:,:)
                                            ! Maximum courant number 
  real(kind=precision), allocatable :: tempx(:,:,:), tempy(:,:,:)
                                            ! Temporary variables

  !---------------------------------------------------------------------------------------------
  ! Restart file
  !---------------------------------------------------------------------------------------------
  character(len=100) hot_file_load

  !---------------------------------------------------------------------------------------------
  ! Record time series files
  !---------------------------------------------------------------------------------------------
  integer :: pts,lys,cols
  parameter(lys=5)
  integer ly_id(lys)

  real(kind=precision) :: x1, x2, y1, y2

  !---------------------------------------------------------------------------------------------
  ! Import and export data
  !---------------------------------------------------------------------------------------------
  real(kind=precision) :: missing_value
  
  ! monthly mean data output
  integer IF_export_clim
  integer n_clim
  real(kind=precision), allocatable :: ssh_mon(:,:)
  real(kind=precision), allocatable :: tb_mon(:,:,:), sb_mon(:,:,:)
  real(kind=precision), allocatable :: ub_mon(:,:,:), vb_mon(:,:,:)
  
  ! instant file
  integer IF_export_instant
  integer nth_field, interval_export_instant
  character(len=200) ncfile_instant

  ! OBC data
  integer IF_import_obc
  integer nth_obc_in
  character(len=200) ncfile_obc_in
  
  ! OBC export parameters
  integer ipts
  real(kind=precision) :: inv_dlon_reg, inv_dlat_reg

  !---------------------------------------------------------------------------------------------
  ! Input file variable
  !---------------------------------------------------------------------------------------------
  ! wind speed and wind stress
  real(kind=precision), allocatable, target :: windx(:,:), windy(:,:)  
  real(kind=precision), allocatable, target :: uflx(:,:), vflx(:,:)
  ! heat flux at the surface
  real(kind=precision), allocatable, target :: uswrf(:,:), dswrf(:,:)  
  real(kind=precision), allocatable, target :: ulwrf(:,:), dlwrf(:,:)  
  real(kind=precision), allocatable, target :: shtfl(:,:), lhtfl(:,:)
  real(kind=precision), allocatable :: jloc_gauss_u(:), jloc_gauss_v(:)
  
  !---------------------------------------------------------------------------------------------
  ! Data assimilation variables
  !---------------------------------------------------------------------------------------------
  integer IF_EAKF_argo
  character(len=100) EAKF_adjusted_file
  character(len=100) EAKF_dir
  character(len=2) ens_char
  integer n_param
  
  integer IF_record_series
  
  !---------------------------------------------------------------------------------------------
  ! Open boundary points coordinates
  !---------------------------------------------------------------------------------------------
  real(kind=precision), allocatable :: reg_xy(:,:)
  
  !---------------------------------------------------------------------------------------------
  ! SST OI Assimilation
  !---------------------------------------------------------------------------------------------
  integer IF_sst_OI_assimilate
  real(kind=precision) :: sst_coeff   ! SST OI assimilation
  
  real(kind=precision),allocatable :: taux(:,:,:),tauy(:,:,:),dqdsst_pom(:,:,:)
  real(kind=precision),allocatable :: netheat_pom(:,:,:),sst_pom(:,:,:)
  
  real(kind=precision),allocatable :: annual_psi(:,:)
  

  character(len=200) glb_dir, init_dir, wavemix_dir, bti_dir

  !---------------------------------------------------------------------------------------------
  ! Namelist variables to be loaded from the parameters file
  !---------------------------------------------------------------------------------------------
  namelist /ctrlparams/                                     &
          ens_char,                                         &
          im, jm, ks,                                       &
          cold_start, IF_save_restart, hot_file_load,       &
          period, period_tide, period_obc,                  &
          lon0, lat0, inv_dlon, inv_dlat,                   &
          Load_OB_tide, Changjiang_diluting,                &
          dte, isplit,                                      &
          mode, Baropg_method, Steps_upd_advdiff,           &
          ymdhms_start, ymdhms_finish,                      &
          ext_data_dir, Use_restart_date, IF_export_clim,   &
          IF_export_instant, interval_export_instant,       &
          IF_import_obc, IF_EAKF_argo, IF_record_series,    &
          IF_sst_OI_assimilate, EAKF_dir, Nint, glb_dir,    &
          init_dir, wavemix_dir, bti_dir,                   &
          wav_freq
  
  contains
  
  !---------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------
  
  subroutine read_parameters                ! Load basic parameters
    implicit none
    
    character(3) windname
    character(2) sstmonth

    !---------------------------------------------------------------------------------------------
    ! Initial parameter values
    !---------------------------------------------------------------------------------------------
    Omega_tide = (/ 28.98410424, 30.0, 15.04106864, 13.94303559 /)
                                            ! In deg/hour
    chj_q = (/ 11008,11903,16825,25254,33345,40342, &
               52183,44065,39315,32952,21817,13413  /)
                                            ! Changjiang diluting monthly fluxes
    
    
    open(1,file='MASNUM_CtrlPARAMS_SETUP_cmp2.50.txt',status='old')
    read(1,nml=ctrlparams)
    close(1)
    
    print*, 'starting date = ', ymdhms_start(1:3)
    
    imm1=im-1
    jmm1=jm-1
    ksm1=ks-1
    kb=ks
    kbm1=kb-1

    nmax = im*2+jm*2                          ! Added by Zhuang Zhanpeng
        
    !---------------------------------------------------------------------------------------------
    ! Define integration time steps and durations
    !---------------------------------------------------------------------------------------------
    dti = dte * isplit                      ! Time step of the internal mode
    iint0 = 0
    
    days_begin = datenum( ymdhms_start )
    days_now   = days_begin
    
    ymdhms    = datevec(days_now)
    date_char = datestr(days_now)    
    
    iend = ( datenum( ymdhms_finish ) - days_begin ) * 86400.0 / dti
    
    !---------------------------------------------------------------------------------------------
    ! Define constants
    !---------------------------------------------------------------------------------------------
    Ra = 6370.e3                            ! Mean radius of the Earth (m)
    PI = 3.141592654
    GRAV = 9.806
    Boltzmann = 1.3806503e-23               ! m2*kg/s2/K  Boltsmann Constant
    Omega_tide = Omega_tide * 24.
    ZERO = 0.0    
    ONE = 1.0
    
    ramp = 1.0
    ramp_obc = 1.0
    ramp_tide = 1.0
    
    cbcmax = 1.
    cbcmin = .0025
    z0b = .01
    kappa = 0.4	
    
    umol = 2.0e-5
    rho_ref = 1000.
    horcon = 0.2
    small = 1.e-10
    TPRNI = 1.0
    
    sst_coeff = dti / 86400.0
    
    Land_Depth = 1.0
    Min_Depth = 10.0
    		
    missing_value = -999.
    
    n_clim = 0
    
    nth_field = 0
		
    sst_bound = 0
		
    !---------------------------------------------------------------------------------------------
    ! Screen display the parameters' value setups
    !---------------------------------------------------------------------------------------------
    print*, ''
    print*, '=========  Control Parameters Setup in this run =========='
    print*, 'cold_start = ',cold_start
    print*, 'MODE = ', mode
    print*, 'IF_save_restart = ',IF_save_restart
    print*, 'period = ',period
    print*, 'hot_file_load = ',hot_file_load
    print*, '[lon0 lat0] = ', lon0,lat0
    print*, 'inv_dlon(dlat) = ',inv_dlon,inv_dlat
    print*, 'dte = ',dte
    print*, 'isplit = ',isplit
    print*, 'Baropg method = ',Baropg_method  
    print*, '=========  Control Parameters Setup =========='
    print*, ''    

    !---------------------------------------------------------------------------------------------
    ! Define the grid interval : dx & dy 
    ! -- lon0 : leftmost line of the grids U(i=1)
    ! -- lat0 : lowest line of the grids   V(j=1)
    ! -- [lon,lat]: defined at the centers of the T-cells
    !---------------------------------------------------------------------------------------------
    allocate(vlon(im), vlon_u(im), vlat(jm), vlat_v(jm))
      vlon = ZERO
      vlat = ZERO
      vlon_u = ZERO
      vlat_v = ZERO
    
    do i=1,im
   	   vlon(i) = lon0 + 1./inv_dlon*(float(i)-0.5)
    enddo		  
    do j=1,jm
   	   vlat(j) = lat0 + (float(j)-0.5)/inv_dlat
    enddo   
    vlon_u = vlon - 1./inv_dlon/2.0           ! Longitudes at UV-point
    vlat_v = vlat - 1./inv_dlat/2.0           ! Latitudes  at UV-point

    !---------------------------------------------------------------------------------------------
    ! Auto-judge whether the case is global or regional
    !---------------------------------------------------------------------------------------------
    if (vlon(im)-360.0 > vlon(1)) then
      regional_run = 0
      print*, '------------------------'
      print*, '| This is a Global run |'
      print*, '------------------------'
    else
      regional_run = 1
      print*, '--------------------------'
      print*, '| This is a Regional run |'
      print*, '--------------------------'
    endif  
    
        
    !---------------------------------------------------------------------------------------------
    ! Load depth   
    !---------------------------------------------------------------------------------------------
    allocate(H_c(im,jm), H_u(im,jm), H_v(im,jm))
    
      H_c = Land_Depth
      H_u = Land_Depth
      H_v = Land_Depth

    allocate(fsm(im,jm))
      fsm = ONE


    if (regional_run) then

      H_c = Land_Depth
      H_u = Land_Depth
      H_v = Land_Depth

      H_c = 0.
      call get_2D_field('Gebco08','Global')
      	      
    else
    	      
      H_c = 0.
      call get_2D_field('Gebco08','Global') 
      
      
    endif
    
    !---------------------------------------------------------------------------------------------
    ! Special treatment of topography in different cases
    !---------------------------------------------------------------------------------------------
    if (regional_run==0) then
    	
      smooth_factor = 0.025
      Max_Depth = 5500.0
      Num_topo_smooth = 0.4/smooth_factor
      
      ! irregular sole water point hard to be eliminated by the current algorithm
      ! very bad treatment!
      if (lat0==-79) H_c(248,182) = Land_Depth
      
      do j=1,jm
    	  do i=1,im
      
          !----------------------------------------------------------------------------------------
          ! Mask the Black Sea
          !----------------------------------------------------------------------------------------
          if (vlon(i)>26.69 .and. vlon(i)<43.69 .and. vlat(j)>40.39 .and. vlat(j)<49.01) then
          	H_c(i,j) = Land_Depth
          endif
          
          !----------------------------------------------------------------------------------------
          ! Mask the Caspian Sea
          !----------------------------------------------------------------------------------------
          if (vlon(i)>44.01 .and. vlon(i)<64.01 .and. vlat(j)>34.49 .and. vlat(j)<50.01) then
          	H_c(i,j) = Land_Depth
          endif
          
          !----------------------------------------------------------------------------------------
          ! Mask the Red Sea
          !----------------------------------------------------------------------------------------
          if (vlon(i)>30.55 .and. vlon(i)<43.69 .and. vlat(j)>12.27 .and. vlat(j)<30.69) then
          	H_c(i,j) = Land_Depth
          endif
          
          !----------------------------------------------------------------------------------------
          ! Mask the Arctic Ocean
          !----------------------------------------------------------------------------------------
          if (vlat(j)>=65) then
          	H_c(i,j) = Land_Depth
          endif
          
          !----------------------------------------------------------------------------------------
          ! Turn on the Qiongzhou Strait of the Hainan Island
          !----------------------------------------------------------------------------------------
          if (vlon(i)==110.5 .and. vlat(j)==20) then
          	H_c(i,j) = 44.
          endif
    
        enddo  
      enddo
    
    elseif (regional_run==1) then
    	
      smooth_factor = 0.025
      Max_Depth = 4000.0
      Num_topo_smooth = 0.4/smooth_factor
    	
!      H_c(1,:) = 0.                         ! Close the western boundary
!      H_c(153,76)=0.0                       ! Reduece the ponit with too shallow water in the ...
!      H_c(173, 6)=0.0                       ! Phillipines
!      do i=128,153                          ! Remove the area in the south Indonesia
!        do j=1,4                            ! ...
!          H_c(i,j)=0.0                      ! ...
!        enddo                               ! ...
!      enddo                                 ! ...          
      !!!!!!!!
      ! set the points land on the northern and southern bound
      H_c(:,jm) = 0.0
      H_c(:, 1) = 0.0

      !!!!!!!!
      ! Mediterranean Sea is set to land
      do i=1,im
        do j=1,jm
          if(vlon(i)<80 .and. vlat(j)>28.5) then
            H_c(i,j) = 0.0
          endif
        enddo
      enddo

      !!!!!!!!
      ! Crimson Sea is set to land
      do i=1,im
        do j=1,jm
          if(vlon(i)<=56 .and. vlon(i)>48 .and. vlat(j)>23 .and. vlat(j)< 29) then
            H_c(i,j) = 0.0
          endif
        enddo
      enddo

      if (ChangJiang_diluting) then         ! Turn on Changjiang Diluting water
        H_c(chj_i,chj_j) = 10.              ! ...
        H_c(chj_i-1,chj_j) = 10.            ! ...
      endif                                 ! ...

      !  land or sea
      !do j=1,jm
      !	do i=1,im
      !		if (vlon(i)>113.39 .and. vlon(i)<115.21 .and. vlat(j)>15.23 .and. vlat(j)<16.48) then
      !			if (H_c(i,j)<=Land_Depth) then
      !  			H_c(i,j) = Min_Depth    ! set to be water
        			!H_c(i,j) = Land_Depth    ! set to be land
      !  		endif
      !		endif
      !	enddo
      !enddo
      
    endif
    
    
    !---------------------------------------------------------------------------------------------
    ! Use box to knock out areas of "dead water"
    !---------------------------------------------------------------------------------------------
    do j=1,jm
      do i=1,im
        if (H_c(i,j)<=Land_Depth) H_c(i,j) = zero    ! Set land depth to be zero for knock out
      enddo
    enddo

    print*,'knock out sole water point within the land',minval(abs(H_c)),maxval(abs(H_c))
    do n=0,10
    do m=0,10
      do j=2,jmm1-n
      do i=2,imm1-m
        x1 = sum( H_c( i-1  , j:j+n ) )
        x2 = sum( H_c( i+m+1, j:j+n ) )
        y1 = sum( H_c( i:i+m, j-1   ) )
        y2 = sum( H_c( i:i+m, j+n+1 ) )
        if ( x1+x2+y1+y2==0 ) then
          H_c( i:i+m, j:j+n ) = ZERO
        endif
      enddo
      enddo
    enddo
    enddo
    print*,'knock out finish'


    
    return
  end subroutine read_parameters
  
!=================================================================================================
  
  subroutine init_common_array
  
    implicit none
    
    !---------------------------------------------------------------------------------------------
    ! Load parameters from setup file
    !---------------------------------------------------------------------------------------------
    call read_parameters
    
    !---------------------------------------------------------------------------------------------
    ! Vertical Coordinate
    !---------------------------------------------------------------------------------------------
    allocate(z(ks), zz(ks), zt(ks), dz(ks), dzz(ks))
      z   = ZERO
      zz  = ZERO
      zt  = ZERO
      dz  = ZERO
      dzz = ZERO
    
    !---------------------------------------------------------------------------------------------
    ! Coriolis parameter
    !---------------------------------------------------------------------------------------------
    allocate(Cor_t(im,jm), Cor_u(im,jm), Cor_v(im,jm))
      Cor_t = ZERO
      Cor_u = ZERO
      Cor_v = ZERO
    
    !---------------------------------------------------------------------------------------------
    ! 2D Common variables
    !---------------------------------------------------------------------------------------------
    
    !---------------------------------------------------------------------------------------------
    ! External mode
    !---------------------------------------------------------------------------------------------
    allocate(etb_c(im,jm), etf_c(im,jm))
    allocate(et_mean(im,jm))
      etb_c = ZERO
      etf_c = ZERO
      et_mean = ZERO
    
    !---------------------------------------------------------------------------------------------
    ! Internal mode
    !---------------------------------------------------------------------------------------------
    allocate(zetf_c(im,jm), zetb_c(im,jm))
    allocate(hb_inc(im,jm), hb_inu(im,jm), hb_inv(im,jm))
    allocate(hm_inc(im,jm), hm_inu(im,jm), hm_inv(im,jm))
    allocate(hf_inc(im,jm), hf_inu(im,jm), hf_inv(im,jm))
    allocate(utf(im,jm), vtf(im,jm))
      zetb_c = ZERO
      zetf_c = ZERO
      hb_inc = ZERO
      hb_inu = ZERO
      hb_inv = ZERO
      hm_inc = ZERO
      hm_inu = ZERO
      hm_inv = ZERO
      hf_inc = ZERO
      hf_inu = ZERO
      hf_inv = ZERO
      utf = ZERO
      vtf = ZERO
   
    allocate(dum(im,jm), dvm(im,jm))
    allocate(dum01(im,jm), dvm01(im,jm))
    allocate(dx_t(im,jm), dy_t(im,jm), dx_u(im,jm), dy_u(im,jm))
    allocate(dx_v(im,jm), dy_v(im,jm), dx_u_down(im,jm), dy_v_left(im,jm))
    allocate(art(im,jm), aru(im,jm), arv(im,jm))
    allocate(int_vm_u(im,jm), int_vm_v(im,jm))
    allocate(Am_vm(im,jm))
    allocate(uab(im,jm), ua(im,jm), uaf(im,jm))
    allocate(vab(im,jm), va(im,jm), vaf(im,jm))
    allocate(adv_ua(im,jm), adv_va(im,jm))
    allocate(diff_ua(im,jm), diff_va(im,jm))
      fsm = ZERO
      dum = ZERO
      dvm = ZERO
      dum01 = ZERO
      dvm01 = ZERO
      dx_t = ZERO
      dy_t = ZERO
      dx_u = ZERO
      dy_u = ZERO
      dx_v = ZERO
      dy_v = ZERO
      dx_u_down = ZERO 
      dy_v_left = ZERO
      art = ZERO
      aru = ZERO
      arv = ZERO
      int_vm_u = ZERO 
      int_vm_v = ZERO
      Am_vm = 500.
      uab = ZERO
      vab = ZERO
      uaf = ZERO
      vaf = ZERO
      ua  = ZERO
      va  = ZERO
      adv_ua = ZERO
      adv_va = ZERO
      diff_ua = ZERO
      diff_va = ZERO
      
    allocate(cbc(im,jm))
    allocate(wusurf(im,jm), wvsurf(im,jm), wtsurf(im,jm), wssurf(im,jm))
    allocate(wubot(im,jm) , wvbot(im,jm))
    allocate(fric_ub(im,jm), fric_vb(im,jm))
      cbc = ZERO
      wusurf = ZERO
      wvsurf = ZERO
      wtsurf = ZERO
      wssurf = ZERO
      wubot = ZERO
      wvbot = ZERO
      fric_ub = ZERO
      fric_vb = ZERO
      
    allocate(psi(im,jm), ml_dep(im,jm))
      psi = ZERO
      ml_dep = ZERO
        
    !=============================================================================================
    ! 3D Common variables
    !=============================================================================================
    !---------------------------------------------------------------------------------------------
    ! Prognostic variables
    !---------------------------------------------------------------------------------------------
    allocate(ub(ks,im,jm), vb(ks,im,jm))
    allocate(wf(ks,im,jm))  
    allocate(tb(ks,im,jm), sb(ks,im,jm))
      ub = ZERO
      vb = ZERO
      wf = ZERO
      tb = ZERO
      sb = ZERO
        
    !---------------------------------------------------------------------------------------------
    ! Mixing coefficients
    !---------------------------------------------------------------------------------------------
    allocate(Km(ks,im,jm), Kh(ks,im,jm), Kq(ks,im,jm))
      Km = 2.0e-4
      Kh = 2.0e-4
      Kq = 2.0e-4
    
    !---------------------------------------------------------------------------------------------
    ! Pressure terms
    !---------------------------------------------------------------------------------------------
    allocate(epg_int_x(im,jm), epg_int_y(im,jm))
      epg_int_x = ZERO
      epg_int_y = ZERO
    
    !---------------------------------------------------------------------------------------------
    ! Thermodynamic variables
    !---------------------------------------------------------------------------------------------
    allocate(rho(ks,im,jm), rclim(ks,im,jm))
    allocate(tclim(ks,im,jm), sclim(ks,im,jm))
      rho   = ZERO
      rclim = ZERO
      tclim = ZERO
      sclim = ZERO
    
    !---------------------------------------------------------------------------------------------
    ! POM Profq and Advq variabes
    !---------------------------------------------------------------------------------------------
    allocate(q2b(kb,im,jm), q2lb(kb,im,jm))
      q2b  = 1.0e-8
      q2lb = 1.0e-8
    
    !---------------------------------------------------------------------------------------------
    ! Variables for debug use
    !---------------------------------------------------------------------------------------------
    allocate(ix(30),jx(30),kx(30))
      ix = ZERO
      jx = ZERO
      kx = ZERO
    
    !---------------------------------------------------------------------------------------------
    ! Surface flux variables 
    !---------------------------------------------------------------------------------------------
    allocate(netheat(im,jm), dqdsst(im,jm), sst_dq(im,jm), sst_OI(im,jm))
      netheat = ZERO
      dqdsst  = ZERO
      sst_dq  = ZERO
      sst_OI  = ZERO
    
    !=============================================================================================
    ! Open boundary variables for regional case 
    !=============================================================================================
    !---------------------------------------------------------------------------------------------
    ! Open boundary variables
    !---------------------------------------------------------------------------------------------
    allocate(ij_obc(2,nmax))
    allocate(bctype(nmax))
      ij_obc = ZERO
      bctype = ZERO
        
    !---------------------------------------------------------------------------------------------
    ! OBC variables
    !---------------------------------------------------------------------------------------------
    allocate(et_obc(nmax), ua_obc(nmax), va_obc(nmax))
    allocate(ub_obc(ks,nmax), vb_obc(ks,nmax))
    allocate(tb_obc(ks,nmax), sb_obc(ks,nmax))  
      et_obc = ZERO
      ua_obc = ZERO
      va_obc = ZERO
      ub_obc = ZERO
      vb_obc = ZERO
      tb_obc = ZERO
      sb_obc = ZERO
    
    !---------------------------------------------------------------------------------------------
    ! Open boundary tide harmonic constants
    !---------------------------------------------------------------------------------------------
    allocate(Amplitude(nmax,numtide), Phase(nmax,numtide))
      Amplitude = ZERO
      Phase = ZERO
    
    !---------------------------------------------------------------------------------------------
    ! Bv variables
    !---------------------------------------------------------------------------------------------
    allocate(Bv_all(ks_bv,im_bv,jm_bv,36))
    allocate(Bv_CN(ks_bv,im_bv,jm_bv))
      Bv_all = ZERO
      Bv_CN = ZERO
      
    !-----------------------------------------------------------------------------------------------
    ! allocate the momery spaces
    !-----------------------------------------------------------------------------------------------
    allocate(courant(im,jm))
      courant = ZERO

    !-----------------------------------------------------------------------------------------------
    ! allocate monthly mean output variables
    !-----------------------------------------------------------------------------------------------
    if (IF_export_clim) then
      allocate(ssh_mon(im,jm))
      allocate(tb_mon(ks,im,jm), sb_mon(ks,im,jm))
      allocate(ub_mon(ks,im,jm), vb_mon(ks,im,jm))
      ssh_mon = ZERO
      tb_mon  = ZERO
      sb_mon  = ZERO
      ub_mon  = ZERO
      vb_mon  = ZERO
    endif
    
    allocate(annual_psi(im,jm))
    annual_psi = ZERO
    
    allocate(bv_var(ks,im,jm), bti_var(ks,im,jm))
      bv_var = ZERO
      bti_var = ZERO

    allocate(bsm1(ks,im,jm), bsm2(ks,im,jm))
      bsm1 = ZERO
      bsm2 = ZERO

    allocate(bsmt(ks,im,jm), bsms(ks,im,jm))
      bsmt = ZERO
      bsms = ZERO
    
    allocate(um(ks,im,jm), vm(ks,im,jm), tm(ks,im,jm), sm(ks,im,jm))
      um = ZERO
      vm = ZERO
      tm = ZERO
      sm = ZERO

    return
  end subroutine init_common_array


!=================================================================================================
  
  subroutine get_2D_field(dataset, var_name, time_res)
  
    use netcdf_mod
    use time_mod
    
    implicit none
    
    integer ncid
    character(len=*) :: dataset
    character(len=*), optional :: var_name, time_res
    
    integer i1,i3
    integer j1,j3
    real(kind=precision) :: x1,x3,x
    real(kind=precision) :: y1,y3,y
    real(kind=precision) :: fin(4),fx
    
    real(kind=8) :: dtime_hrs, dtime_ref
    real(kind=8) :: rec_num
    
    real(kind=8) :: temp_i, temp_j
    integer :: time_len, lon_len, lat_len, irec
    integer, allocatable :: iloc(:), jloc(:)

    character(len=200) file_2D_nm, dir_var
    
    real(kind=precision), allocatable :: time_2d(:), lon_2d(:), lat_2d(:)
    real(kind=precision), allocatable :: var_2d(:,:), var_3d(:,:,:), lsmask(:,:)
    real(kind=precision), allocatable :: time_dt(:), lon_dt(:), lat_dt(:)
    real(kind=precision), allocatable, target :: ones_2d(:,:)
    real(kind=precision), pointer :: var_2d_interp(:,:)
    real(kind=precision), pointer :: vlon_interp(:), vlat_interp(:), mask(:,:) 
    
    real(kind=precision) :: null_value, scale_factor, add_offset
    real(kind=precision) :: x1_offset, x3_offset
    real(kind=precision) :: maxlat, minlat
    real(kind=precision) :: Land_Default
        
    logical existed
    
    !---------------------------------------------------------------------------------------------
    ! End of definition part
    !---------------------------------------------------------------------------------------------

    allocate( vlon_interp(im), vlat_interp(jm) )
    
    dtime_hrs = days_now*24.0
    allocate( iloc(im), jloc(jm) )
     
    add_offset = ZERO
    scale_factor = ONE
    Land_Default = ZERO

    !---------------------------------------------------------------------------------------------
    ! 1) North Pacific topography extracted from Gebco-08 (0.5'x0.5') by Zhao Y.X.
    !---------------------------------------------------------------------------------------------
    if (dataset == 'Gebco08') then
      
      var_2d_interp => H_c
      vlon_interp   => vlon
      vlat_interp   => vlat
      mask => fsm
      
      if (var_name=='North_Pacific') then
        file_2D_nm = trim(adjustl(ext_data_dir))//'/Gebco08/gridone_np.nc'
      elseif (var_name=='Global') then
        file_2D_nm = trim(adjustl(ext_data_dir))//'/Gebco08/gridone_global_int.nc'
      endif
      print*,'loading terrain file: ', file_2D_nm
      
      call open_nc(ncid,file_2D_nm,'r')
      lon_len  = get_dimension_len(ncid,'lon')
      lat_len  = get_dimension_len(ncid,'lat')
      
      allocate( lon_2d(lon_len), lat_2d(lat_len) )
      allocate( var_2d(lon_len,lat_len) )
      
      call readnc(ncid, 'lon', lon_2d)
      call readnc(ncid, 'lat', lat_2d)
      call readnc(ncid, 'H_c', var_2d)
      
      call get_attribute(ncid, 'missing_value',  null_value, 'H_c')
            
      print*,'max val', maxval(var_2d), minval(var_2d)

    !---------------------------------------------------------------------------------------------
    ! 2) NCEP Reanalysis 1 wind-10m : 2.5x2.5 degrees
    !---------------------------------------------------------------------------------------------
    elseif (dataset == 'NCEP1-wind') then
      
      if (var_name=='uwnd') then      
        allocate( windx(im,jm) )
        var_2d_interp => windx
        vlon_interp   => vlon_u
        vlat_interp   => vlat
        mask          => dum
      elseif (var_name=='vwnd') then      
        allocate( windy(im,jm) )
        var_2d_interp => windy
        vlon_interp   => vlon
        vlat_interp   => vlat_v
        mask          => dvm
      endif      
      
      dtime_ref = datenum(1,1,1,0,0,0) * 24.0
      dtime_hrs = dtime_hrs - dtime_ref
      
      file_2D_nm = trim(adjustl(ext_data_dir))//'/NCEP Wind/'//var_name &
                                              //'.sig995.'//date_char(1:4)//'.nc'
      print*,'loading 2D ext. data file: ', file_2D_nm, 'at ',date_char(1:8)
      
      call open_nc(ncid,file_2D_nm,'r')
      lon_len  = get_dimension_len(ncid,'lon')
      lat_len  = get_dimension_len(ncid,'lat')
      time_len = get_dimension_len(ncid,'time')
      
      allocate( lon_2d(lon_len), lat_2d(lat_len), time_2d(time_len) )
      allocate( var_2d(lon_len,lat_len) )
            
      call readnc(ncid, 'time',time_2d)

      rec_num = (  dtime_hrs - time_2d(1) )   &
              / ( time_2d(2) - time_2d(1) )
      rec_num = max(rec_num, 0.0)
      irec = int( rec_num + 1 )

      call readnc(ncid, 'lon' , lon_2d)
      call readnc(ncid, 'lat' , lat_2d)
      call readnc(ncid, var_name, var_2d, irec)

      call get_attribute(ncid, 'missing_value',  null_value, var_name)
      call get_attribute(ncid, 'scale_factor', scale_factor, var_name)
      call get_attribute(ncid, 'add_offset',     add_offset, var_name)

      print*,'add_offset',add_offset,null_value, scale_factor

      call close_nc(ncid)

    !---------------------------------------------------------------------------------------------
    ! 3) QuickScat wind-10m : wind stress (positive downward) or wind speed
    !---------------------------------------------------------------------------------------------
    elseif (dataset == 'QS-wind') then
    
      print*, 'var-name ',var_name
      if (var_name=='surface_downward_eastward_stress') then
        var_2d_interp => wusurf
        vlon_interp   => vlon_u
        vlat_interp   => vlat
        mask          => dum
      elseif (var_name=='surface_downward_northward_stress') then
        var_2d_interp => wvsurf
        vlon_interp   => vlon
        vlat_interp   => vlat_v
        mask          => dvm
      elseif (var_name=='eastward_wind') then
        allocate(windx(im,jm))
        var_2d_interp => windx
        vlon_interp   => vlon_u
        vlat_interp   => vlat
        mask          => dum
      elseif (var_name=='northward_wind') then
        allocate(windy(im,jm))
        var_2d_interp => windy
        vlon_interp   => vlon
        vlat_interp   => vlat_v
        mask          => dvm
      endif      

      dtime_ref = datenum(1900,1,1,0,0,0) * 24.0
      dtime_hrs = dtime_hrs - dtime_ref
      
      file_2D_nm = trim(adjustl(ext_data_dir))//'/QuickScat-MWF/'//date_char(1:8)   &
                                              //'00_daily-ifremer-L3-MWF.nc'
      print*,'loading 2D ext. data file: ', trim(adjustl(file_2D_nm)),          &
             ' at ',date_char(1:8)
      
      call open_nc(ncid,file_2D_nm,'r')
      lon_len  = get_dimension_len(ncid,'longitude')
      lat_len  = get_dimension_len(ncid,'latitude')
      
      allocate( lon_2d(lon_len), lat_2d(lat_len) )
      allocate( var_2d(lon_len,lat_len) )
            
      call readnc(ncid, 'longitude' , lon_2d)
      call readnc(ncid, 'latitude'  , lat_2d)
      call readnc(ncid, var_name    , var_2d)

      call get_attribute(ncid, '_FillValue'  ,   null_value, var_name)
      call get_attribute(ncid, 'scale_factor', scale_factor, var_name)
      call get_attribute(ncid, 'add_offset'  ,   add_offset, var_name)
            
      print*,'add_offset', add_offset, null_value, scale_factor

      call close_nc(ncid)

    !---------------------------------------------------------------------------------------------
    ! 4) NCEP Reanalysis 2 : Heat and momentum fluxes at Gaussian grid
    !---------------------------------------------------------------------------------------------
    elseif (dataset == 'NCEP2:surface_fluxes') then
      
      if (var_name=='uflx') then
      	
        if(.not.allocated(uflx)) allocate( uflx(im,jm) )
        var_2d_interp => uflx
        vlon_interp   => vlon_u
        vlat_interp   => vlat
        mask          => dum
        
      elseif (var_name=='vflx') then
        
        if(.not.allocated(vflx)) allocate( vflx(im,jm) )
        var_2d_interp => vflx
        vlon_interp   => vlon
        vlat_interp   => vlat_v
        mask          => dvm
     
      else
      	
        !------------------------------------------------------------------------------------------
        ! six Heat flux variables
        !------------------------------------------------------------------------------------------
      	if (var_name=='uswrf') then      
          allocate( uswrf(im,jm) )
          var_2d_interp => uswrf
        elseif (var_name=='ulwrf') then
          allocate( ulwrf(im,jm) )
          var_2d_interp => ulwrf
        elseif (var_name=='dswrf') then
          allocate( dswrf(im,jm) )
          var_2d_interp => dswrf
        elseif (var_name=='dlwrf') then
          allocate( dlwrf(im,jm) )
          var_2d_interp => dlwrf
        elseif (var_name=='shtfl') then
          allocate( shtfl(im,jm) )
          var_2d_interp => shtfl
        elseif (var_name=='lhtfl') then
          allocate( lhtfl(im,jm) )
          var_2d_interp => lhtfl
        endif

        vlon_interp   => vlon
        vlat_interp   => vlat
        mask          => fsm
        
      endif
      
      !---------------------------------------------------------------------------------------------
      ! Define time(in hours) relative to the reference date of NCEP-2 
      !---------------------------------------------------------------------------------------------
      dtime_ref = datenum(1800,1,1,0,0,0) * 24.0
      dtime_hrs = dtime_hrs - dtime_ref
      
      !---------------------------------------------------------------------------------------------
      ! Load basic data information 
      !---------------------------------------------------------------------------------------------
      if (time_res=='4xdaily') then
        dir_var = '/NCEP-2/'//var_name//'.sfc.gauss.'
        file_2D_nm = trim(adjustl(ext_data_dir))//trim(adjustl(dir_var))//date_char(1:4)//'n.nc'
      elseif (time_res=='monthly') then
        file_2D_nm = trim(adjustl(ext_data_dir))//'/NCEP-2 monthly/'//var_name//'.sfc.mon.mean.nc'
      elseif (time_res=='clim') then
        file_2D_nm = trim(adjustl(ext_data_dir))//'/NCEP-2clim/'//var_name//'.sfc.mon.ltm.nc'
      endif
      
      print*,'loading 2D ext. data file: ', trim(file_2D_nm), 'at ',date_char(1:8)
      
      call open_nc(ncid,file_2D_nm,'r')

      lon_len  = get_dimension_len(ncid, 'lon' )
      lat_len  = get_dimension_len(ncid, 'lat' )
      time_len = get_dimension_len(ncid, 'time')

      allocate( lon_2d(lon_len), lat_2d(lat_len), time_2d(time_len) )
      allocate( var_2d(lon_len,lat_len) )
            
      call readnc(ncid, 'time',time_2d)

      !---------------------------------------------------------------------------------------------
      ! Determine which time slide is gonna be extracted
      !---------------------------------------------------------------------------------------------
      allocate( time_dt(time_len) )
      time_dt = time_2d - dtime_hrs
      irec = maxloc(time_dt, dim=1, mask=time_dt<=0)
      
      call readnc(ncid, 'lon' , lon_2d)
      call readnc(ncid, 'lat' , lat_2d)

      !---------------------------------------------------------------------------------------------
      ! Clim data is a normal 3D NC variable, while other data is a time-unlimited NC variable
      !---------------------------------------------------------------------------------------------
      if (time_res=='clim') then
      	if(.not.allocated(var_3d)) allocate( var_3d(lon_len,lat_len,time_len) )
        call readnc(ncid, var_name, var_3d)
        do j=1,lat_len
        	do i=1,lon_len
        		var_2d(i,j) = var_3d(i,j,ymdhms(2))
        	enddo
        enddo
        print*,'loading Climtological NCEP-2 data:',var_name
      else    	
        call readnc(ncid, var_name, var_2d, irec)
        print*,'irec',irec,lat_len,lon_len
      endif

      call get_attribute(ncid, 'missing_value',  null_value, var_name)
      call get_attribute(ncid, 'scale_factor', scale_factor, var_name)
      call get_attribute(ncid, 'add_offset',     add_offset, var_name)

      call close_nc(ncid)

      print*,'add_offset',add_offset,null_value, scale_factor
      
      !---------------------------------------------------------------------------------------------
      ! Load data's land-sea mask
      !---------------------------------------------------------------------------------------------
      if(.not.allocated(lsmask)) allocate(lsmask(lon_len,lat_len))
      file_2D_nm = trim(adjustl(ext_data_dir))//'/NCEP-2/lsmask.ncep2.gauss.nc'
      call open_nc(ncid,file_2D_nm,'r')
      call readnc(ncid, 'lsmask',lsmask)
      call close_nc(ncid)

      !---------------------------------------------------------------------------------------------
      ! Set Land value to be null
      !---------------------------------------------------------------------------------------------
      do j=1,lat_len
      	do i=1,lon_len
      	  if (lsmask(i,j)/=0) var_2d(i,j)=null_value
        enddo
      enddo

    !---------------------------------------------------------------------------------------------
    ! 5) monthly climtological dq/dsst
    !---------------------------------------------------------------------------------------------
    elseif (dataset == 'dqdsst') then
      
      if (allocated(ones_2d)) deallocate(ones_2d)
      allocate(ones_2d(im,jm))
      ones_2d = ONE
      
      var_2d_interp => dqdsst
      vlon_interp   => vlon
      vlat_interp   => vlat
      mask          => ones_2d
      
      file_2D_nm = trim(adjustl(ext_data_dir))//'/NCEP-2clim/dqdsst_mon_fill_land.nc'
      print*,'loading 2D ext. data file: ', trim(adjustl(file_2D_nm))
      
      call open_nc(ncid,file_2D_nm,'r')
      lon_len  = get_dimension_len(ncid,'lon')
      lat_len  = get_dimension_len(ncid,'lat')
      
      if(allocated(var_3d)) deallocate(var_3d)
     	allocate( var_3d(lon_len,lat_len,12) )
      allocate( var_2d(lon_len,lat_len) )
      allocate( lon_2d(lon_len), lat_2d(lat_len) )
     	      
      call readnc(ncid, 'lon'    , lon_2d)
      call readnc(ncid, 'lat'    , lat_2d)
      call readnc(ncid, 'dqdsst' , var_3d)
      
      do j=1,lat_len
      	do i=1,lon_len
      		var_2d(i,j) = var_3d(i,j,ymdhms(2))
       	enddo
      enddo
      
      print*,'dqdsst factor',add_offset,null_value, scale_factor

      call close_nc(ncid)


    !---------------------------------------------------------------------------------------------
    ! 5) monthly climtological dq/dsst
    !---------------------------------------------------------------------------------------------
    elseif (dataset == 'tmi_amsre_SST') then
      
      var_2d_interp => sst_OI
      vlon_interp   => vlon
      vlat_interp   => vlat
      mask          => fsm
      
      
      file_2D_nm = trim(adjustl(ext_data_dir))//'/SST_satellite/tmi_amsre/' &
      //date_char(1:8)//'-REMSS-L4LRfnd-GLOB-v01-fv03-tmi_amsre'
      
      ! Two kinds of suffixes of TMI-AMSRE files are present.
      ! Try both, and turn back if neither found.
      inquire( file=trim(adjustl(file_2D_nm))//'.nc', exist=existed )
      if (.not.existed) file_2D_nm = trim(adjustl(file_2D_nm))//'_rt'
      file_2D_nm = trim(adjustl(file_2D_nm))//'.nc'
      inquire( file=trim(adjustl(file_2D_nm)), exist=existed )
      if (.not.existed) then
      	print*,'tmi_amsre_SST file:',trim(adjustl(file_2D_nm)),' not found, sst_OI NOT updated'
      	return
      endif
      
      print*,'loading TMI-AMSRE L4 SST file: ', trim(adjustl(file_2D_nm))
            
      call open_nc(ncid,file_2D_nm,'r')
      lon_len  = get_dimension_len(ncid,'lon')
      lat_len  = get_dimension_len(ncid,'lat')
      
      if(allocated(lsmask)) deallocate(lsmask)
      allocate( lsmask(lon_len,lat_len) )
      allocate( var_2d(lon_len,lat_len) )
      allocate( lon_2d(lon_len), lat_2d(lat_len) )
     	      
      call readnc(ncid, 'lon'    , lon_2d)
      call readnc(ncid, 'lat'    , lat_2d)
      call readnc(ncid, var_name , var_2d)
      call readnc(ncid, 'mask'   , lsmask)
      
      call get_attribute(ncid, '_FillValue'  ,   null_value, var_name)
      call get_attribute(ncid, 'scale_factor', scale_factor, var_name)
      call get_attribute(ncid, 'add_offset',     add_offset, var_name)
      
      add_offset = add_offset - 273.16   ! transfer Kelvin to Centigrade
      
      print*,'add_offset',add_offset, null_value, scale_factor
      
      Land_Default = null_value
      
      !---------------------------------------------------------------------------------------------
      ! Set Land value to be null
      !---------------------------------------------------------------------------------------------
      do j=1,lat_len
      	do i=1,lon_len
      	  if (lsmask(i,j)/=33) var_2d(i,j)=null_value   ! 33:Ocean 65:Coast
        enddo
      enddo

      !---------------------------------------------------------------------------------------------
      ! rearrange the longitude and SST from format [-180 180] to format [0 360]
      !---------------------------------------------------------------------------------------------
      m = lon_len/2

      do i=1,lon_len
      	if (lon_2d(i)<0) lon_2d(i)=lon_2d(i)+360.0
      enddo
      
      lon_2d = (/lon_2d(m+1:lon_len), lon_2d(1:m)/)
      do j=1,lat_len
        var_2d(:,j) = (/var_2d(m+1:lon_len,j), var_2d(1:m,j)/)
      enddo
            
      call close_nc(ncid)

    !---------------------------------------------------------------------------------------------
    ! 5) monthly climtological surface flux from POM
    !---------------------------------------------------------------------------------------------
    elseif (dataset == 'pomsurface') then
      
      if (var_name=='wndx') then
      	
        if(.not.allocated(uflx)) allocate( uflx(im,jm) )
        var_2d_interp => uflx
        vlon_interp   => vlon_u
        vlat_interp   => vlat
        mask          => dum
        
      elseif (var_name=='wndy') then
        
        if(.not.allocated(vflx)) allocate( vflx(im,jm) )
        var_2d_interp => vflx
        vlon_interp   => vlon
        vlat_interp   => vlat_v
        mask          => dvm
     
      else
      	
        !------------------------------------------------------------------------------------------
        ! six Heat flux variables
        !------------------------------------------------------------------------------------------
      	if (var_name=='dqdsst') then      
!          allocate( dqdsst(im,jm) )
          dqdsst = ZERO
          var_2d_interp => dqdsst
        elseif (var_name=='sst') then
!          allocate( sst_dq(im,jm) )
          sst_dq = ZERO
          var_2d_interp => sst_dq
        elseif (var_name=='netheat') then
!          allocate( netheat(im,jm) )
          netheat = ZERO
          var_2d_interp => netheat
        endif

        vlon_interp   => vlon
        vlat_interp   => vlat
        mask          => fsm
        
      endif

      
      !---------------------------------------------------------------------------------------------
      ! Define time(in hours) relative to the reference date of NCEP-2 
      !---------------------------------------------------------------------------------------------
      dtime_ref = datenum(1800,1,1,0,0,0) * 24.0
      dtime_hrs = dtime_hrs - dtime_ref
      
      !---------------------------------------------------------------------------------------------
      ! Load basic data information 
      !---------------------------------------------------------------------------------------------
      if (time_res=='4xdaily') then
        dir_var = '/NCEP-2/'//var_name//'.sfc.gauss.'
        file_2D_nm = trim(adjustl(ext_data_dir))//trim(adjustl(dir_var))//date_char(1:4)//'n.nc'
      elseif (time_res=='monthly') then
        file_2D_nm = trim(adjustl(ext_data_dir))//'/NCEP-2 monthly/'//var_name//'.sfc.mon.mean.nc'
      elseif (time_res=='clim') then
        file_2D_nm = trim(adjustl(ext_data_dir))//'/surface.nc'
      endif
      
      print*,'loading 2D ext. data file: ', trim(file_2D_nm), 'at ',date_char(1:8)
      
      call open_nc(ncid,file_2D_nm,'r')

      lon_len  = get_dimension_len(ncid, 'lon' )
      lat_len  = get_dimension_len(ncid, 'lat' )
      time_len = get_dimension_len(ncid, 'time')

      allocate( lon_2d(lon_len), lat_2d(lat_len), time_2d(time_len) )
      allocate( var_2d(lon_len,lat_len) )
            
      call readnc(ncid, 'time',time_2d)

      !---------------------------------------------------------------------------------------------
      ! Determine which time slide is gonna be extracted
      !---------------------------------------------------------------------------------------------
      allocate( time_dt(time_len) )
      time_dt = time_2d - dtime_hrs
      irec = maxloc(time_dt, dim=1, mask=time_dt<=0)
      
      call readnc(ncid, 'lon' , lon_2d)
      call readnc(ncid, 'lat' , lat_2d)

      !---------------------------------------------------------------------------------------------
      ! Clim data is a normal 3D NC variable, while other data is a time-unlimited NC variable
      !---------------------------------------------------------------------------------------------
      if (time_res=='clim') then
      	if(.not.allocated(var_3d)) allocate( var_3d(lon_len,lat_len,time_len) )
        call readnc(ncid, var_name, var_3d)
        do j=1,lat_len
        	do i=1,lon_len
        		var_2d(i,j) = var_3d(i,j,ymdhms(2))
        	enddo
        enddo
        print*,'loading Climtological NCEP-2 data:',var_name
      else    	
        call readnc(ncid, var_name, var_2d, irec)
        print*,'irec',irec,lat_len,lon_len
      endif

      call get_attribute(ncid, 'null_value',  null_value, var_name)
!      call get_attribute(ncid, 'scale_factor', scale_factor, var_name)
!      call get_attribute(ncid, 'add_offset',     add_offset, var_name)

      call close_nc(ncid)

      print*,'add_offset',add_offset,null_value, scale_factor

    !---------------------------------------------------------------------------------------------
    ! End of the dataset selections
    !---------------------------------------------------------------------------------------------
    endif
    
    !---------------------------------------------------------------------------------------------
    ! Interpolation part
    !---------------------------------------------------------------------------------------------  
    
      !---------------------------------------------------------------------------------------------
      ! 1 : Find the nearest lower-left point
      !---------------------------------------------------------------------------------------------
      allocate( lon_dt(lon_len), lat_dt(lat_len) )
     	
      do i=1,im                           ! only Ascending longitude is considered
        x1_offset = zero
        if (vlon_interp(i) > 360) x1_offset = 360
        if (vlon_interp(i) <=  0) x1_offset = -360
        lon_dt = lon_2d - vlon_interp(i) + x1_offset
        iloc(i) = maxloc(lon_dt, dim=1, mask=lon_dt<=0)
      enddo

     	if (lat_2d(2)>lat_2d(1)) then       ! Ascending latitude
       	do j=1,jm
          lat_dt = lat_2d - vlat_interp(j)
          jloc(j) = maxloc(lat_dt, dim=1, mask=lat_dt<=0)
          if (jloc(j)>lat_len .or. jloc(j)<1)	jloc(j)=1   ! abnormal value
          jloc(j) = min( jloc(j), lat_len-1 )
        enddo
      else                                ! Descending latitude (NCEP custom)
      	do j=1,jm
          lat_dt = lat_2d - vlat_interp(j)
          jloc(j) = minloc(lat_dt, dim=1, mask=lat_dt>=0)
          if (jloc(j)>lat_len .or. jloc(j)<1)	jloc(j)=1   ! abnormal value
          jloc(j) = min( jloc(j), lat_len-1 )
        enddo
      endif
      
      maxlat = maxval(lat_2d)
      minlat = minval(lat_2d)
      
      !-------------------------------------------------------------------------------------------
      ! 2 : Interpolate to the model grids:
      ! Comment: null_value--"external field"   Land_default--"model variable"
      !------------------------------------------------------------------------------------------            
      var_2d_interp = Land_Default
      
      do i=1,im
         
        x1_offset = zero
        x3_offset = zero
         
        !-------------------------------------------------------------------------------------------
        ! Define the "squre box" which the output point (x,y) locates in.
        ! 1) the longitude direction
        !-------------------------------------------------------------------------------------------
        i1 = iloc(i)
        i3 = i1 + 1

        if (i1==lon_len) then
          i3 = 1
          x3_offset = 360.  
        elseif(i1>lon_len .or. i1<1) then 
        !-------------------------------------------------------------------------------------------
        ! abnormal value means lon_dt<=0 is false, i.e., point exceeds western bound.
        !-------------------------------------------------------------------------------------------
          i1 = lon_len
          i3 = 1
          x1_offset = -360.
        endif

        x  = vlon_interp(i)
        if(x > 360) x = x-360
        if(x <=  0) x = x+360
        
        x1 = lon_2d(i1) + x1_offset
        x3 = lon_2d(i3) + x3_offset
        
        do j=1,jm
        	
          if (mask(i,j)==zero) cycle
        	
          !-------------------------------------------------------------------------------------------
          ! Define the "squre box" which the output point (x,y) locates in.
          ! 2) the latitude direction
          !-------------------------------------------------------------------------------------------
          j1 = jloc(j)
          j3 = j1 + 1

          y1 = lat_2d(j1)
          y3 = lat_2d(j3)
          y  = vlat_interp(j)
          
          y  = max(y,minlat)
          y  = min(y,maxlat)
                                        
          fin(1) = var_2d(i1,j1)
          fin(2) = var_2d(i1,j3)
          fin(3) = var_2d(i3,j3)
          fin(4) = var_2d(i3,j1)
           
          !-------------------------------------------------------------------------------------------
          ! Bilinear interpolation
          !-------------------------------------------------------------------------------------------
          call bilinear_4points_HL(x,y,x1,x3,y1,y3,fin,fx,null_value)
          
          if ( fx/=null_value ) then
            var_2d_interp(i,j) = add_offset + scale_factor * fx
          endif
          
        enddo
      enddo
      
     deallocate( var_2d )
     if ( allocated(var_3d) ) deallocate( var_3d )
     
     nullify( var_2d_interp, vlon_interp, vlat_interp )
    
!    call test_interp_field_2d('H_c',zero)
    
    return
  end subroutine get_2D_field

!=================================================================================================

  
  subroutine get_3D_field(dataset, var_name, var_name2)
  
    use netcdf_mod
    use time_mod
    
    implicit none
    
    integer ncid, ncid2
    character(len=*) dataset
    character(len=*), optional :: var_name, var_name2
    
    integer i1,i3
    integer j1,j3
    integer k1,k2
    real(kind=precision) :: x1,x3,x
    real(kind=precision) :: y1,y3,y
    real(kind=precision) :: z1,z2,z_mets
    real(kind=precision) :: fin(8),fx, null_value
    
    real(kind=8) :: dtime_hrs, dtime_ref
    real(kind=8) :: rec_num
    
    real(kind=8) :: temp_i, temp_j
    integer :: time_len, lon_len, lat_len, dep_len, irec
    integer, allocatable :: iloc(:), jloc(:)

    character(len=200) dir, ncfile, ncfile2
    character(len=2) tt
    
    real(kind=precision), allocatable :: time_3d(:), lon_3d(:), lat_3d(:), dep_3d(:), dep_3dt(:)
    real(kind=precision), allocatable :: layer(:), depth(:,:)
    real(kind=precision), allocatable :: time_dt(:), lon_dt(:), lat_dt(:)    
    real(kind=precision), allocatable :: var_4d(:,:,:,:), var_3d(:,:,:)
    real(kind=precision), allocatable, target :: var_sst(:,:,:)
    
    real(kind=precision), pointer :: var_3d_interp(:,:,:)
    real(kind=precision), pointer :: vlon_interp(:), vlat_interp(:)
    real(kind=precision), pointer :: mask(:,:)
    
    real(kind=precision), allocatable :: x1_offset(:)
    
    real(kind=precision), allocatable :: rmean_1d(:)
    real(kind=precision) :: arsea, rtot, atot, yi, zi 
            
    allocate( vlon_interp(im), vlat_interp(jm) )
    allocate( iloc(im), jloc(jm) )
    

    !----------------------------------------------------------------------------------------------
    ! PHC climatological annual 1x1 degrees
    !---------------------------------------------------------------------------------------------- 
    if (dataset == 'PHC-annual' ) then
      
      vlon_interp   => vlon
      vlat_interp   => vlat
      mask => fsm
      
      if ( var_name == 't_an' ) then
        var_3d_interp => tclim
        ncfile = trim(adjustl(ext_data_dir))//'/Levitus/temperature_annual_phcfilled.nc'
        print*,'loading initial temperature file: ', ncfile
      elseif ( var_name == 's_an' ) then
        var_3d_interp => sclim
        ncfile = trim(adjustl(ext_data_dir))//'/Levitus/salinity_annual_phcfilled.nc'
        print*,'loading initial salinity file: ', ncfile
      elseif ( var_name == 'rhoz' ) then
        ncfile = trim(adjustl(ext_data_dir))//'/Levitus/rhoz_annual_phcfilled.nc'
        print*,'loading Climatological rhoz file: ', ncfile
      endif
      
      call open_nc(ncid, trim(adjustl(ncfile)), 'r')
      lon_len  = get_dimension_len(ncid, 'lon')
      lat_len  = get_dimension_len(ncid, 'lat')
      dep_len  = get_dimension_len(ncid, 'depth')
            
      allocate( lon_3d(lon_len), lat_3d(lat_len), dep_3d(dep_len) )
      allocate( var_3d(lon_len, lat_len, dep_len) )
      
      call readnc(ncid, 'lon'    , lon_3d)
      call readnc(ncid, 'lat'    , lat_3d)
      call readnc(ncid, 'depth'  , dep_3d)
      call readnc(ncid, var_name , var_3d)
                               
      print*,'max val', maxval(var_3d), minval(var_3d)
      call close_nc(ncid)
      
    !---------------------------------------------------------------------------------------------
    ! Levitus climatological monthly 1x1 degrees
    !---------------------------------------------------------------------------------------------  
    elseif (dataset == 'WOA09-monthly' ) then
      
      vlon_interp   => vlon
      vlat_interp   => vlat
      mask => fsm
      
      if ( var_name == 't_an' ) then
        if (var_name2 /= 'sst' ) then
        	var_3d_interp => tclim
        else
        	print*,'SST going to load...'
        	allocate(var_sst(ks,im,jm))
        	var_3d_interp => var_sst
        endif
        	
        ncfile = trim(adjustl(ext_data_dir))//'/Levitus/temperature_monthly_levitusfilled.nc'
        print*,'loading initial temperature file: ', trim(ncfile)
      elseif ( var_name == 's_an' ) then
        var_3d_interp => sclim
        ncfile = trim(adjustl(ext_data_dir))//'/Levitus/salinity_monthly_levitusfilled.nc'
        print*,'loading initial salinity file: ',  trim(ncfile)
      endif

      call open_nc(ncid, trim(adjustl(ncfile)), 'r')
      lon_len  = get_dimension_len(ncid, 'lon')
      lat_len  = get_dimension_len(ncid, 'lat')
      dep_len  = get_dimension_len(ncid, 'depth')
      time_len = get_dimension_len(ncid, 'time')
            
      allocate( lon_3d(lon_len), lat_3d(lat_len), dep_3d(dep_len) )      
      allocate( var_3d(lon_len, lat_len, dep_len) )
      allocate( var_4d(lon_len, lat_len, dep_len, time_len) )
      
      call readnc(ncid, 'lon'    , lon_3d)
      call readnc(ncid, 'lat'    , lat_3d)
      call readnc(ncid, 'depth'  , dep_3d)
      call readnc(ncid, var_name , var_4d)

      !--------------------------------------------------------------------------------------------
      ! Extract the climatological data of the month
      !--------------------------------------------------------------------------------------------
      do k=1,dep_len
      	do j=1,lat_len
      		do i=1,lon_len
      			var_3d(i,j,k) = var_4d(i,j,k,ymdhms(2))
      		enddo
      	enddo
      enddo
      
      deallocate(var_4d)
                   
      print*,'max val', maxval(var_3d), minval(var_3d)
      call close_nc(ncid)      

    elseif (dataset == 'global_res' ) then
      
      vlon_interp   => vlon
      vlat_interp   => vlat
      
      if(var_name == 'tb') then
        var_3d_interp => tclim
        print*,'loading initial temperature file: '
      elseif(var_name == 'sb') then
        var_3d_interp => sclim
        print*,'loading initial sallinity file: '
      elseif(var_name == 'rmean') then
        var_3d_interp => rclim
        print*,'loading initial density file: '
      endif
      
      ncfile = trim(adjustl(ext_data_dir))//'/initial.nc'
      print*, trim(adjustl(ncfile))

      call open_nc(ncid, trim(adjustl(ncfile)), 'r')
      lon_len  = get_dimension_len(ncid, 'lon')
      lat_len  = get_dimension_len(ncid, 'lat')
      dep_len  = get_dimension_len(ncid, 'layer')

      allocate( var_3d(lon_len, lat_len, dep_len) )
      
      call readnc(ncid, var_name , var_3d)

      call close_nc(ncid)

      do i=1,im
        do j=1,jm
          do k=1,ks
            if(var_3d(i,j+23,ks-k+1) < -90) then
              var_3d_interp(k,i,j) = 0.0
            else
              var_3d_interp(k,i,j) = 0.5 *     &
                (var_3d(i,j+23,ks-k+1) + var_3d(i,j+23,ks-k)) * fsm(i,j)
            endif
          enddo
        enddo
      enddo

      print*,'max val', maxval(var_3d_interp), minval(var_3d_interp)

      return

    !---------------------------------------------------------------------------------------------
    ! wave mixing data loading
    !---------------------------------------------------------------------------------------------  
    elseif (dataset == 'bv' ) then
    	
      vlon_interp   => vlon
      vlat_interp   => vlat
      mask => fsm
      var_3d_interp => bv_var
      
      ncfile=trim(adjustl(ext_data_dir))//'/Bv/monthly/2002'//date_char(5:6)//'_bvv.nc'

      print*,'loading wave mixing file (Bv): ', trim(ncfile)
      
      call open_nc(ncid, trim(adjustl(ncfile)), 'r')
      
      lon_len  = get_dimension_len(ncid, 'lon')
      lat_len  = get_dimension_len(ncid, 'lat')
      dep_len  = get_dimension_len(ncid, 'dep')
            
      allocate( lon_3d(lon_len), lat_3d(lat_len), dep_3d(dep_len) )      
      allocate( var_3d(lon_len, lat_len, dep_len) )
      
      call readnc(ncid, 'lon'  , lon_3d)
      call readnc(ncid, 'lat'  , lat_3d)
      call readnc(ncid, 'dep'  , dep_3d)
      call readnc(ncid, 'bt2'  , var_3d)
      
      !call get_attribute(ncid, 'missing_value',  null_value, 'bv')
      null_value = 999
      
      dep_3d = -dep_3d

      call close_nc(ncid)
      
      do i=1,lon_len
        do j=1,lat_len
          do k=1,dep_len
            if(abs(var_3d(i,j,k)) > 1e5) then
              var_3d(i,j,k) = null_value
            endif
          enddo
        enddo
      enddo
      
    !---------------------------------------------------------------------------------------------
    ! Interpolation part
    !---------------------------------------------------------------------------------------------  
      !---------------------------------------------------------------------------------------------
      ! Find the nearest lower-left point
      ! only ascending longitude and latitude are considered
      !---------------------------------------------------------------------------------------------
      allocate( lon_dt(lon_len), lat_dt(lat_len) )
      allocate( x1_offset(im) )
   	  x1_offset = zero
   	  
      do i=1,im-3
        !------------------------------------------------------------------------------------------
        ! Move the western out-of-range point to the eastern bound
        !------------------------------------------------------------------------------------------
        if (vlon_interp(i) < lon_3d(1)) x1_offset(i) = -360
        if (vlon_interp(i) > 360) x1_offset(i) = 360
        if (vlon_interp(i) <=  0) x1_offset(i) = -360

        lon_dt = lon_3d - vlon_interp(i) + x1_offset(i)
        iloc(i) = maxloc(lon_dt, dim=1, mask=lon_dt<=0)        	
        if (iloc(i)>=lon_len) then
          print*,'Warning : Get 3D field exceed eastern bound. (bv) i=',i,vlon_interp(i),iloc(i),lon_len
          stop
        endif
      enddo
      !print*,'iloc get3d',iloc

     	do j=2,jmm1
        lat_dt = lat_3d - vlat_interp(j)
        jloc(j) = maxloc(lat_dt, dim=1, mask=lat_dt<=0)
        if (jloc(j)>lat_len .or. jloc(j)<1)	jloc(j)=1   ! abnormal value
        jloc(j) = min( jloc(j), lat_len-1 )
      enddo
      !print*,'jloc get3d',jloc
                
      !---------------------------------------------------------------------------------------------
      ! 3D Interpolation Part
      !---------------------------------------------------------------------------------------------
      do i=1,im-3
      	
        i1 = iloc(i)
        i3 = i1 + 1 

        x  = vlon_interp(i) - x1_offset(i)
        x1 = lon_3d(i1)
        x3 = lon_3d(i3)
      	
        do j=2,jmm1
        	if (mask(i,j)==0) cycle
          
          !-------------------------------------------------------------------------------------------
          ! Define the "squre box" which the output point (x,y) locates in
          !-------------------------------------------------------------------------------------------
          j1 = jloc(j)
          j3 = j1 + 1
      
          y  = vlat_interp(j)
          y1 = lat_3d(j1)
          y3 = lat_3d(j3)

          do k=ks,1,-1
          	
            z_mets = - H_c(i,j) * zz(k)
            if (z_mets > 500) cycle
      
            !-------------------------------------------------------------------------------------
            ! Rules of vertical interpolation :
            ! For  annual climatological (down to 5500m), if z>5500, interpolate to the bottom;
            ! For monthly climatological (down to 1500m), if z>1500, do not interpolate, use
            ! the already loaded annual climatological data instead.
            !-------------------------------------------------------------------------------------
            do m=1,dep_len-1
              if ( z_mets>=dep_3d(m) .and. z_mets<dep_3d(m+1) ) then
                k1 = m
                k2 = m+1
                exit
              elseif ( z_mets<dep_3d(1) ) then
                k1 = 1
                k2 = 1
                exit
              elseif ( z_mets==dep_3d(dep_len) ) then
                k1 = dep_len
                k2 = dep_len
                exit
              elseif ( z_mets>dep_3d(dep_len)              &
                      .and. dataset == 'PHC-annual' ) then
                k1 = dep_len
                k2 = dep_len
                exit
              elseif ( z_mets>dep_3d(dep_len)              &
                      .and. dataset == 'WOA09-monthly' ) then
                z_mets = -1
                exit
              endif
            enddo

            if (z_mets==-1) cycle    ! Skip value settings for depths deeper than 1500m in loading
                                     ! monthly data 
            
            z1 = dep_3d(k1)
            z2 = dep_3d(k2)

            fin(1) = var_3d(i1,j1,k1)
            fin(2) = var_3d(i1,j3,k1)
            fin(3) = var_3d(i3,j3,k1)
            fin(4) = var_3d(i3,j1,k1)
      
            fin(5) = var_3d(i1,j1,k2)
            fin(6) = var_3d(i1,j3,k2)
            fin(7) = var_3d(i3,j3,k2)
            fin(8) = var_3d(i3,j1,k2)

            call trilinear_8points_HL(x,y,z_mets,x1,x3,y1,y3,z1,z2,fin,fx,null_value)
            if (fx == null_value) then
              fx = ZERO
            endif
      
            var_3d_interp(k,i,j) = fx

          enddo
        enddo
      enddo   	

      print*,'var_3d_interp:',maxval(var_3d_interp),minval(var_3d_interp)
      
      deallocate( var_3d )
      nullify( var_3d_interp, vlon_interp, vlat_interp )
    
      return

    !---------------------------------------------------------------------------------------------
    ! wave mixing data loading
    !---------------------------------------------------------------------------------------------  
    elseif (dataset == 'bti' ) then
    	
      vlon_interp   => vlon
      vlat_interp   => vlat
      mask => fsm
      var_3d_interp => bti_var
      
      write(tt,'(i2.2)') ymdhms(2)

      ncfile=trim(adjustl(bti_dir))//'/01_bti_'//tt//'.nc'

      print*,'loading internal wave mixing file (Bti): ', trim(ncfile)
      
      call open_nc(ncid, trim(adjustl(ncfile)), 'r')
      
      lon_len  = get_dimension_len(ncid, 'lon')
      lat_len  = get_dimension_len(ncid, 'lat')
      dep_len  = get_dimension_len(ncid, 'dep')
            
      allocate( lon_3d(lon_len), lat_3d(lat_len), dep_3d(dep_len) )      
      allocate( var_3d(lon_len, lat_len, dep_len) )
      
      call readnc(ncid, 'lon'  , lon_3d)
      call readnc(ncid, 'lat'  , lat_3d)
      call readnc(ncid, 'dep'  , dep_3d)

      call readnc(ncid, 'bti'   , var_3d)

      !call get_attribute(ncid, 'missing_value',  null_value, 'bv')
      null_value = 999.0
      
      !dep_3d = -dep_3d

      call close_nc(ncid)
      
      do i=1,lon_len
        do j=1,lat_len
          do k=1,dep_len
            if(abs(var_3d(i,j,k)) > 1e5) then
              var_3d(i,j,k) = null_value
            endif
          enddo
        enddo
      enddo

    !---------------------------------------------------------------------------------------------
    ! Interpolation part
    !---------------------------------------------------------------------------------------------  
      !---------------------------------------------------------------------------------------------
      ! Find the nearest lower-left point
      ! only ascending longitude and latitude are considered
      !---------------------------------------------------------------------------------------------
      allocate( lon_dt(lon_len), lat_dt(lat_len) )
      allocate( x1_offset(im) )
        x1_offset = zero
   	  
      do i=2,imm1
        !------------------------------------------------------------------------------------------
        ! Move the western out-of-range point to the eastern bound
        !------------------------------------------------------------------------------------------
        if (vlon_interp(i) < lon_3d(1)) x1_offset(i) = -360
        if (vlon_interp(i) > 360) x1_offset(i) = 360
        if (vlon_interp(i) <=  0) x1_offset(i) = -360

        lon_dt = lon_3d - vlon_interp(i) + x1_offset(i)
        iloc(i) = maxloc(lon_dt, dim=1, mask=lon_dt<=0)

        if (iloc(i)>=lon_len) then

          print*,'Warning : Get 3D field exceed eastern bound. (bti) i=',i,vlon_interp(i),iloc(i),lon_len
!          stop
          !!!!!!!!
          if (abs(vlon_interp(i)) < 5) then
            iloc(i) = 1
          elseif (abs(vlon_interp(i)) > 355) then
            iloc(i) = lon_len - 1
          else
            print*,'can not find the neighbour!!!'
            stop
          endif
          x1_offset(i) = 0

          print*,'get the value from neighbour point'
!          write(*,'(15G15.5)') vlon_interp(i),iloc(i),lon_3d(iloc(i)),lon_3d(iloc(i)+1)

        endif
      enddo
      !print*,'iloc get3d',iloc

      do j=2,jmm1
        lat_dt = lat_3d - vlat_interp(j)
        jloc(j) = maxloc(lat_dt, dim=1, mask=lat_dt<=0)
        if (jloc(j)>lat_len .or. jloc(j)<1)	jloc(j)=1   ! abnormal value
        jloc(j) = min( jloc(j), lat_len-1 )
      enddo
      !print*,'jloc get3d',jloc
                
      !---------------------------------------------------------------------------------------------
      ! 3D Interpolation Part
      !---------------------------------------------------------------------------------------------
      do i=2,imm1
      	
        i1 = iloc(i)
        i3 = i1 + 1 

        x  = vlon_interp(i) - x1_offset(i)
        x1 = lon_3d(i1)
        x3 = lon_3d(i3)
      	
        do j=2,jmm1
          if (mask(i,j)==0) cycle

          
          !-------------------------------------------------------------------------------------------
          ! Define the "squre box" which the output point (x,y) locates in
          !-------------------------------------------------------------------------------------------
          j1 = jloc(j)
          j3 = j1 + 1
      
          y  = vlat_interp(j)
          y1 = lat_3d(j1)
          y3 = lat_3d(j3)

          do k=ks,1,-1
          	
            z_mets = - H_c(i,j) * zz(k)
            if (z_mets > 1500) cycle

            !-------------------------------------------------------------------------------------
            ! Rules of vertical interpolation :
            ! For  annual climatological (down to 5500m), if z>5500, interpolate to the bottom;
            ! For monthly climatological (down to 1500m), if z>1500, do not interpolate, use
            ! the already loaded annual climatological data instead.
            !-------------------------------------------------------------------------------------
            do m=1,dep_len-1
              if ( z_mets>=dep_3d(m) .and. z_mets<dep_3d(m+1) ) then
                k1 = m
                k2 = m+1
                exit
              elseif ( z_mets<dep_3d(1) ) then
                k1 = 1
                k2 = 1
                exit
              elseif ( z_mets==dep_3d(dep_len) ) then
                k1 = dep_len
                k2 = dep_len
                exit
              elseif ( z_mets>dep_3d(dep_len)              &
                      .and. dataset == 'PHC-annual' ) then
                k1 = dep_len
                k2 = dep_len
                exit
              elseif ( z_mets>dep_3d(dep_len)              &
                      .and. dataset == 'WOA09-monthly' ) then
                z_mets = -1
                exit
              endif
            enddo

            if (z_mets==-1) cycle    ! Skip value settings for depths deeper than 1500m in loading
                                     ! monthly data 
            
            z1 = dep_3d(k1)
            z2 = dep_3d(k2)

            fin(1) = var_3d(i1,j1,k1)
            fin(2) = var_3d(i1,j3,k1)
            fin(3) = var_3d(i3,j3,k1)
            fin(4) = var_3d(i3,j1,k1)
      
            fin(5) = var_3d(i1,j1,k2)
            fin(6) = var_3d(i1,j3,k2)
            fin(7) = var_3d(i3,j3,k2)
            fin(8) = var_3d(i3,j1,k2)

            call trilinear_8points_HL(x,y,z_mets,x1,x3,y1,y3,z1,z2,fin,fx,null_value)

            if (fx == null_value) then
              fx = ZERO
            endif

            var_3d_interp(k,i,j) = fx

          enddo
        enddo
      enddo

!      do k=1,ks
!        var_3d_interp(k,imm1,:) = var_3d_interp(k,im-3,:)
!        var_3d_interp(k,im-2,:) = var_3d_interp(k,im-3,:)
!      enddo

      print*,'var_3d_interp:',maxval(var_3d_interp),minval(var_3d_interp)
      
      deallocate( var_3d )
      nullify( var_3d_interp, vlon_interp, vlat_interp )

      return

    !---------------------------------------------------------------------------------------------
    ! wave mixing data loading
    !---------------------------------------------------------------------------------------------  
    elseif (dataset == 'wave_mix' ) then

      vlon_interp   => vlon
      vlat_interp   => vlat
      mask => fsm

      if(var_name=='bv_wtv2') then
        bv_var = ZERO
        var_3d_interp => bv_var
      elseif(var_name=='bsm1') then
        bsm1 = ZERO
        var_3d_interp => bsm1
      elseif(var_name=='bsm2') then
        bsm2 = ZERO
        var_3d_interp => bsm2
      endif

      ncfile=trim(adjustl(wavemix_dir))//'/MSR_allmix_'//'2018'//date_char(5:8)//'.nc'

      irec = int(ymdhms(4) / wav_freq) + 1

      print*,'loading wave mixing file: ', trim(ncfile)
      print*, var_name, irec
      
      call open_nc(ncid, trim(adjustl(ncfile)), 'r')
      
      lon_len  = get_dimension_len(ncid, 'lon')
      lat_len  = get_dimension_len(ncid, 'lat')
      dep_len  = get_dimension_len(ncid, 'dep')
            
      allocate( lon_3d(lon_len), lat_3d(lat_len), dep_3d(dep_len) )      
      allocate( var_3d(lon_len, lat_len, dep_len) )
      
      call readnc(ncid, 'lon'  , lon_3d)
      call readnc(ncid, 'lat'  , lat_3d)
      call readnc(ncid, 'dep'  , dep_3d)
      call readnc(ncid, var_name  , var_3d, irec)

      dep_3d = -dep_3d

      null_value = 999.

      do i=1,lon_len
        do j=1,lat_len
          if(var_3d(i,j,1)>1e10) then
            var_3d(i,j,:) = null_value
          endif
        enddo
      enddo

!      print*,'max val', maxval(var_3d), minval(var_3d)
      call close_nc(ncid)

    !----------------------------------------------------------------------------------------------
    ! global analysis results by MASNUM circulation model
    !---------------------------------------------------------------------------------------------- 
    elseif (dataset == 'glo_res_MAS' ) then
           
      ncfile = trim(adjustl(init_dir))

      if ( var_name == 't' ) then
        vlon_interp   => vlon
        vlat_interp   => vlat
        mask => fsm
        var_3d_interp => tb
        print*,'loading initial temperature file: ', ncfile
      elseif ( var_name == 's' ) then
        vlon_interp   => vlon
        vlat_interp   => vlat
        mask => fsm
        var_3d_interp => sb
        print*,'loading initial salinity file: ', ncfile
      elseif ( var_name == 'u' ) then
        vlon_interp   => vlon_u
        vlat_interp   => vlat
        mask => dum
        var_3d_interp => ub
        print*,'loading initial u-velocity file: ', ncfile
      elseif ( var_name == 'v' ) then
        vlon_interp   => vlon
        vlat_interp   => vlat_v
        mask => dvm
        var_3d_interp => vb
        print*,'loading initial v-velocity file: ', ncfile
      endif
      
      call open_nc(ncid, trim(adjustl(ncfile)), 'r')
      lon_len  = get_dimension_len(ncid, 'lon')
      lat_len  = get_dimension_len(ncid, 'lat')
      dep_len  = get_dimension_len(ncid, 'dep')
            
      allocate( lon_3d(lon_len), lat_3d(lat_len))
      allocate( var_3d(lon_len, lat_len, dep_len), dep_3d(dep_len) )
      
      call readnc(ncid, 'lon'  , lon_3d)
      call readnc(ncid, 'lat'  , lat_3d)
      call readnc(ncid, 'dep'  , dep_3d)
      call readnc(ncid, var_name , var_3d)

      null_value = 999.0

      print*,'max val', maxval(var_3d), minval(var_3d)
      call close_nc(ncid)

    !---------------------------------------------------------------------------------------------
    ! Interpolation part
    !---------------------------------------------------------------------------------------------  

      !---------------------------------------------------------------------------------------------
      ! Find the nearest lower-left point
      ! only ascending longitude and latitude are considered
      !---------------------------------------------------------------------------------------------
      allocate( lon_dt(lon_len), lat_dt(lat_len) )
        allocate( x1_offset(im) )
        x1_offset = zero
   	  
        do i=1,im
        !------------------------------------------------------------------------------------------
        ! Move the western out-of-range point to the eastern bound
        !------------------------------------------------------------------------------------------
        if (vlon_interp(i) < lon_3d(1)) x1_offset(i) = -360
        if (vlon_interp(i) > 360) x1_offset(i) = 360
        if (vlon_interp(i) <=  0) x1_offset(i) = -360

        lon_dt = lon_3d - vlon_interp(i) + x1_offset(i)
        iloc(i) = maxloc(lon_dt, dim=1, mask=lon_dt<=0)        	
        if (iloc(i)>=lon_len) then
          print*,'Warning : Get 3D field exceed eastern bound. (glo_res_MAS) i=',i
          print*,vlon_interp(i),iloc(i),lon_len
          stop
        endif
      enddo
      !print*,'iloc get3d',iloc

      do j=1,jm
        lat_dt = lat_3d - vlat_interp(j)
        jloc(j) = maxloc(lat_dt, dim=1, mask=lat_dt<=0)
        if (jloc(j)>lat_len .or. jloc(j)<1)	jloc(j)=1   ! abnormal value
        jloc(j) = min( jloc(j), lat_len-1 )
      enddo
      !print*,'jloc get3d',jloc
      
      !---------------------------------------------------------------------------------------------
      ! 3D Interpolation Part
      !---------------------------------------------------------------------------------------------
      do i=1,im
      	
        i1 = iloc(i)
        i3 = i1 + 1 

        x  = vlon_interp(i) - x1_offset(i)
        x1 = lon_3d(i1)
        x3 = lon_3d(i3)
      	
        do j=1,jm
          if (mask(i,j)==0) cycle

          !-------------------------------------------------------------------------------------------
          ! Define the "squre box" which the output point (x,y) locates in
          !-------------------------------------------------------------------------------------------
          j1 = jloc(j)
          j3 = j1 + 1
      
          y  = vlat_interp(j)
          y1 = lat_3d(j1)
          y3 = lat_3d(j3)

          do k=ksm1,1,-1
          	
            z_mets = - H_c(i,j) * zz(k)
      
            !-------------------------------------------------------------------------------------
            ! Rules of vertical interpolation :
            ! For  annual climatological (down to 5500m), if z>5500, interpolate to the bottom;
            ! For monthly climatological (down to 1500m), if z>1500, do not interpolate, use
            ! the already loaded annual climatological data instead.
            !-------------------------------------------------------------------------------------
            do m=1,dep_len-1
              if ( z_mets<=dep_3d(m+1) .and. z_mets>dep_3d(m) ) then
                k1 = m
                k2 = m+1
                exit
              elseif ( z_mets<dep_3d(1) ) then
                k1 = 1
                k2 = 1
                exit
              elseif ( z_mets>=dep_3d(dep_len) ) then
                k1 = dep_len
                k2 = dep_len
                exit
              endif
            enddo

            if (z_mets==-1) cycle    ! Skip value settings for depths deeper than 1500m in loading
                                     ! monthly data 

            z1 = dep_3d(k1)
            z2 = dep_3d(k2)

            fin(1) = var_3d(i1,j1,k1)
            fin(2) = var_3d(i1,j3,k1)
            fin(3) = var_3d(i3,j3,k1)
            fin(4) = var_3d(i3,j1,k1)
      
            fin(5) = var_3d(i1,j1,k2)
            fin(6) = var_3d(i1,j3,k2)
            fin(7) = var_3d(i3,j3,k2)
            fin(8) = var_3d(i3,j1,k2)

            call trilinear_8points_HL(x,y,z_mets,x1,x3,y1,y3,z1,z2,fin,fx,null_value)

            if (fx == null_value) then
              fx = ZERO
            endif

            var_3d_interp(k,i,j) = fx
          
          enddo

        enddo
      enddo   	

      deallocate( lon_3d, lat_3d )
      deallocate( var_3d, dep_3d )

      print*,'var_3d_interp:',maxval(var_3d_interp),minval(var_3d_interp)
    
    return

    endif        
    
    !---------------------------------------------------------------------------------------------
    ! Interpolation part
    !---------------------------------------------------------------------------------------------  
      !---------------------------------------------------------------------------------------------
      ! Find the nearest lower-left point
      ! only ascending longitude and latitude are considered
      !---------------------------------------------------------------------------------------------
      allocate( lon_dt(lon_len), lat_dt(lat_len) )
        allocate( x1_offset(im) )
        x1_offset = zero
   	  
      do i=1,im
        !------------------------------------------------------------------------------------------
        ! Move the western out-of-range point to the eastern bound
        !------------------------------------------------------------------------------------------
        if (vlon_interp(i)<lon_3d(1)) x1_offset(i) = -360
        if (vlon_interp(i) > 360) x1_offset(i) = 360
        if (vlon_interp(i) <=  0) x1_offset(i) = -360

        lon_dt = lon_3d - vlon_interp(i) + x1_offset(i)
        iloc(i) = maxloc(lon_dt, dim=1, mask=lon_dt<=0)        	
        if (iloc(i)>=lon_len) then
          print*,'Warning : Get 3D field exceed eastern bound. i=',i,vlon_interp(i),iloc(i),lon_len
          stop
        endif
      enddo
      !print*,'iloc get3d',iloc

     	do j=1,jm
        lat_dt = lat_3d - vlat_interp(j)
        jloc(j) = maxloc(lat_dt, dim=1, mask=lat_dt<=0)
        if (jloc(j)>lat_len .or. jloc(j)<1)	jloc(j)=1   ! abnormal value
        jloc(j) = min( jloc(j), lat_len-1 )
      enddo
      !print*,'jloc get3d',jloc
                
      !---------------------------------------------------------------------------------------------
      ! Define rmean : 1) compute the horizontally averaged density (rmean_1d)
      !---------------------------------------------------------------------------------------------
      if (var_name=='rhoz') then
      	
        !---------------------------------------------------------------------------------------------
        ! 1) compute the horizontally averaged density (rmean_1d)
        !---------------------------------------------------------------------------------------------
      	allocate(rmean_1d(dep_len))
    	  do k=1,dep_len
    	  	rtot = ZERO
    	  	atot = ZERO
    	  	do j=1,jm
    	    do i=1,im
    	  		arsea = mask(i,j)*art(i,j)
    	      rtot  = rtot + var_3d(iloc(i),jloc(j),k)*arsea
    	      atot  = atot + arsea
    	    enddo
    	    enddo
    	    rmean_1d(k) = rtot / atot
    	  enddo
    	  print*,'rmean-1d :'
    	  print*, rmean_1d
        
        !---------------------------------------------------------------------------------------------
        ! 2) Interpolate rmean_1d to rmean-3d array (rclim)
        !---------------------------------------------------------------------------------------------
        do j=1,jm
        do i=1,im
        	do k=1,ksm1
        	  zi = -H_c(i,j)*zz(k)
        	  call interp1d(dep_len,dep_3d,rmean_1d,zi,yi,RHO_REF-1000)
        	  rclim(k,i,j) = yi
          enddo
          rclim(ks,i,j) = rclim(ksm1,i,j)
        enddo
        enddo
                        
        deallocate(rmean_1d)
        return
        
      endif
      

      !---------------------------------------------------------------------------------------------
      ! 3D Interpolation Part
      !---------------------------------------------------------------------------------------------
      do i=1,im
      	
        i1 = iloc(i)
        i3 = i1 + 1 

        x  = vlon_interp(i) - x1_offset(i)
        x1 = lon_3d(i1)
        x3 = lon_3d(i3)
      	
        do j=1,jm
        	if (mask(i,j)==0) cycle
          
          !-------------------------------------------------------------------------------------------
          ! Define the "squre box" which the output point (x,y) locates in
          !-------------------------------------------------------------------------------------------
          j1 = jloc(j)
          j3 = j1 + 1
      
          y  = vlat_interp(j)
          y1 = lat_3d(j1)
          y3 = lat_3d(j3)
          
          do k=ks,1,-1
          	
          	z_mets = - H_c(i,j) * zz(k)
      
            !-------------------------------------------------------------------------------------
            ! Rules of vertical interpolation :
            ! For  annual climatological (down to 5500m), if z>5500, interpolate to the bottom;
            ! For monthly climatological (down to 1500m), if z>1500, do not interpolate, use
            ! the already loaded annual climatological data instead.
            !-------------------------------------------------------------------------------------
          	do m=1,dep_len-1
          	  if ( z_mets>=dep_3d(m) .and. z_mets<dep_3d(m+1) ) then
          	  	k1 = m
          	  	k2 = m+1
          	  	exit
          	  elseif ( z_mets<dep_3d(1) ) then
          	  	k1 = 1
          	  	k2 = 1
          	  	exit
          	  elseif ( z_mets==dep_3d(dep_len) ) then
          	  	k1 = dep_len
          	  	k2 = dep_len
          	  	exit
          	  elseif ( z_mets>dep_3d(dep_len)              &
          	          .and. dataset == 'PHC-annual' ) then
          	    k1 = dep_len
          	    k2 = dep_len
          	    exit
          	  elseif ( z_mets>dep_3d(dep_len)              &
          	          .and. dataset == 'WOA09-monthly' ) then
          	    z_mets = -1
          	    exit
          	  endif
            enddo
            
            if (z_mets==-1) cycle    ! Skip value settings for depths deeper than 1500m in loading
                                     ! monthly data 
            
            z1 = dep_3d(k1)
            z2 = dep_3d(k2)
      
            fin(1) = var_3d(i1,j1,k1)
            fin(2) = var_3d(i1,j3,k1)
            fin(3) = var_3d(i3,j3,k1)
            fin(4) = var_3d(i3,j1,k1)
      
            fin(5) = var_3d(i1,j1,k2)
            fin(6) = var_3d(i1,j3,k2)
            fin(7) = var_3d(i3,j3,k2)
            fin(8) = var_3d(i3,j1,k2)
            
                        
            call trilinear_8points_HL(x,y,z_mets,x1,x3,y1,y3,z1,z2,fin,fx,null_value)
      
            if (fx == null_value) then
              fx = ZERO
            endif

            var_3d_interp(k,i,j) = fx
          
          enddo
          
        enddo
      enddo   	
      
      print*,'var_3d_interp:',maxval(var_3d_interp),minval(var_3d_interp)

    if (var_name2 == 'sst' ) then
      do j=1,jm
      	do i=1,im
      		sst_dq(i,j) = var_3d_interp(ksm1,i,j)
      	enddo
      enddo
      deallocate(var_sst)
    endif
      
    deallocate( var_3d )
    nullify( var_3d_interp, vlon_interp, vlat_interp )
    
    return
  end subroutine get_3D_field


!=================================================================================================

  
  subroutine test_interp_field_2d(var_name, null_value)
    
    use netcdf_mod
    use time_mod
    implicit none
    
    integer ncid
    real(kind=precision), optional :: null_value
    character(len=*) var_name
    character(len=100) file_name_check
    real(kind=precision), pointer :: var_2d(:,:)
        
    if (var_name=='H_c') then
      var_2d => H_c
    elseif (var_name=='fsm') then
      var_2d => fsm
    elseif (var_name=='windx') then
      var_2d => windx
    elseif (var_name=='uswrf') then
      var_2d => uswrf
    elseif (var_name=='lhtfl') then
      var_2d => lhtfl
    elseif (var_name=='uflx') then
      var_2d => uflx
    elseif (var_name=='vflx') then
      var_2d => vflx
    elseif (var_name=='wtsurf') then
      var_2d => wtsurf
    elseif (var_name=='courant') then
      var_2d => courant
    elseif (var_name=='netheat') then
      var_2d => netheat
    elseif (var_name=='dqdsst') then
      var_2d => dqdsst
    elseif (var_name=='sst_OI') then
      var_2d => sst_OI
    elseif (var_name=='zetf_c') then
      var_2d => zetf_c
    elseif (var_name=='ua') then
      var_2d => ua
    elseif (var_name=='va') then
      var_2d => va
    endif

    file_name_check = 'interp_'//var_name//'_test.nc'
    call open_nc(ncid, file_name_check, 'c')
    print*,'output file for check: ',file_name_check,ncid
    
    call dimension_define(ncid, 'lon', im, 'lon', 5)
    call dimension_define(ncid, 'lat', jm, 'lat', 5)
  
    call variable_define(ncid, var_name, 5, ['lon','lat'])
    if (present(null_value)) then
      call set_attribute(ncid, 'missing_value', null_value, var_name)
    endif
    
    call end_define(ncid)

    call writenc(ncid,'lon',vlon)
    call writenc(ncid,'lat',vlat)
    call writenc(ncid, var_name, var_2d)

    call close_nc(ncid)
        
    nullify(var_2d)
  
    return
  end subroutine test_interp_field_2d
  
!=================================================================================================
  
  subroutine test_interp_field_3d(var_name, null_value)
    
    use netcdf_mod
    use time_mod
    implicit none
    
    integer ncid
    real(kind=precision), optional :: null_value
    character(len=*) var_name
    character(len=100) file_name_check
    real(kind=precision), pointer :: var_3d(:,:,:)
        
    if (var_name=='tb') then
      var_3d => tb
    elseif (var_name=='sb') then
      var_3d => sb
    elseif (var_name=='rho') then
      var_3d => rho
    elseif (var_name=='rclim') then
      var_3d => rclim
    endif

    file_name_check = 'interp3d_'//var_name//'_test.nc'
    call open_nc(ncid, file_name_check, 'c')
    print*,'output file for check: ',file_name_check,ncid
    
    call dimension_define(ncid, 'lon',   im, 'lon',   5)
    call dimension_define(ncid, 'lat',   jm, 'lat',   5)
    call dimension_define(ncid, 'layer', ks, 'layer', 5)
  
    call variable_define(ncid, var_name, 5, ['layer','lon','lat'])
    if (present(null_value)) then
      call set_attribute(ncid, 'missing_value', null_value, var_name)
    endif
    
    call end_define(ncid)

    call writenc(ncid,'lon',vlon)
    call writenc(ncid,'lat',vlat)
    call writenc(ncid,'layer',zz)
    call writenc(ncid, var_name, var_3d)

    call close_nc(ncid)
        
    nullify(var_3d)
  
    return
  end subroutine test_interp_field_3d
  
  
!=======================================================================

  subroutine interp1d(klen,xo,yo,xi,yi,null_value)
    
    implicit none
    
    integer, intent(in) :: klen
    real(kind=precision), intent(in)  :: xo(klen), yo(klen), xi
    real(kind=precision), intent(in)  :: null_value
    real(kind=precision), intent(out) :: yi

    real(kind=precision) :: dx(klen)
    integer :: kloc

    dx = xo - xi
    
    
    if (xo(2)>xo(1)) then                   ! Ascending 
      kloc = maxloc(dx, dim=1, mask=dx<=0)
      if (xi>=xo(klen)) kloc = -1
      if (dx(1)>0) kloc = 0
    elseif (xo(2)<xo(1)) then               ! Descending 
      kloc = minloc(dx, dim=1, mask=dx>=0)
      if (xi<=xo(klen)) kloc = -1
      if (dx(1)<0) kloc = 0
    endif
    !if (null_value==-999)    print*,'kloc',kloc
    
    if (kloc==klen+1) kloc=0  ! sometimes kloc=ks_len+1 instead of 0 if out of range  
    
    if (kloc==0) then           ! xi below the range of xo
    	yi = null_value
    elseif (kloc==-1) then      ! xi above the range of xo
    	yi = yo(klen)
    else
    	yi = ( yo(kloc+1) - yo(kloc) ) / ( xo(kloc+1) - xo(kloc) )   &
   	     * ( xi - xo(kloc) ) + yo(kloc)
   	endif

    !if (null_value==-999) then
    !print*,'xo',xo
    !print*,'yo',yo
    !print*,'xi yi',xi,yi
    !print*,'dx',dx
    !endif

    
    return
  end subroutine interp1d
  


!======================================================================  


  subroutine bilinear_4points_HL(x,y,x1,x3,y1,y3,                 &
                                 fin,fx,null_value)
    !**********************************************
    !----------------------------------------------
    !
    !     x1 < x < x2, x4 < x < x3
    !     y1 < y < y4, y2 < y < y3
    !
    !          (2)---(3)
    !           |     |
    !           |     |
    !          (1)---(4)
    !
    !----------------------------------------------
    !**********************************************
    
    implicit none
    
    real(kind=precision), intent(in)  :: x1,x3,x
    real(kind=precision), intent(in)  :: y1,y3,y
    real(kind=precision), intent(in)  :: fin(4)
    real(kind=precision), intent(out) :: fx
    real(kind=precision), intent(in), optional  :: null_value
    integer :: ni
    
    real(kind=precision) :: t0x,t0y
    real(kind=precision) :: wt_square(4), wa
    
    t0x = abs(x-x1)/abs(x3-x1)
    t0y = abs(y-y1)/abs(y3-y1)
    
    wt_square(1) = (1. - t0x) * (1. - t0y)
    wt_square(2) = (1. - t0x) *       t0y
    wt_square(3) =       t0x  *       t0y
    wt_square(4) =       t0x  * (1. - t0y)
    
    ! watch out the missing value data
    wa = 0.0
    fx = 0.0
    if (.not.present(null_value)) then
    	fx = fx + sum(fin*wt_square)
      return
    else
    	do ni=1,4
        if ( fin(ni)/=null_value ) then
      	  wa = wa + wt_square(ni)
      	  fx = fx + fin(ni)*wt_square(ni)
        endif
      enddo
    endif
    
    if (minval(wt_square)<0) then
    	print*,'bi:',fin,fx
    	print*,'bi2',wt_square,wa
    	print*,'bi3',t0x,t0y
    	print*,'x',x,x1,x3
    	print*,'y',y,y1,y3
    endif

    
    ! re-normalize the weight
    if (wa==0.) then
      fx = null_value
    else
      fx = fx / wa
    endif
    
  return

  end subroutine bilinear_4points_HL


!=================================================================================================
  
  subroutine trilinear_8points_HL(x,y,zm,x1,x3,y1,y3,z1,z2,fin,      &   
                                      fx,null_value)
                                               
    !**********************************************
    !----------------------------------------------
    !
    !          (2)---(3)        z1 square(1,2,3,4)
    !           |     |
    !           |     |
    !          (1)---(4)        z2 square(5,6,7,8)
    !
    ! z1=z2 is allowed
    !----------------------------------------------
    !**********************************************
    
    implicit none
    
    real(kind=precision), intent(in)  :: x1,x3,x
    real(kind=precision), intent(in)  :: y1,y3,y
    real(kind=precision), intent(in)  :: z1,z2,zm
    real(kind=precision), intent(in)  :: fin(8)
    real(kind=precision), intent(out) :: fx
    real(kind=precision), intent(in), optional  :: null_value
    
    real(kind=precision) :: t0x,t0y,t0z
    real(kind=precision) :: wt_cube(8), wa
    integer :: ni
    
    ! cube length normalization
    t0x = abs(x-x1)/abs(x3-x1)
    t0y = abs(y-y1)/abs(y3-y1)
    if (z2==z1) then
    	t0z = 0.0
    else
    	t0z = abs(zm-z1)/abs(z2-z1)
    endif
    
    wt_cube(1) = (1. - t0x) * (1. - t0y) * (1. - t0z)
    wt_cube(2) = (1. - t0x) *       t0y  * (1. - t0z)
    wt_cube(3) =       t0x  *       t0y  * (1. - t0z)
    wt_cube(4) =       t0x  * (1. - t0y) * (1. - t0z)
    
    wt_cube(5) = (1. - t0x) * (1. - t0y) * t0z
    wt_cube(6) = (1. - t0x) *       t0y  * t0z
    wt_cube(7) =       t0x  *       t0y  * t0z
    wt_cube(8) =       t0x  * (1. - t0y) * t0z
    
    ! watch out the missing value data
    wa = 0.0
    fx = 0.0
    
    if (.not.present(null_value)) then
      fx = fx + sum(fin*wt_cube)
      return
    else
      do ni=1,8
        if ( fin(ni)/=null_value ) then
        	wa = wa + wt_cube(ni)
        	fx = fx + fin(ni)*wt_cube(ni)
        endif
      enddo
    endif
      
    ! re-normalize the weight
    if (wa==0.) then
      fx = null_value
    else
      fx = fx / wa
    endif
    
    return

  end subroutine trilinear_8points_HL

  
!=================================================================================================



end module def_common_mod  

!*************************************************************************************************
!-------------------------------------------------------------------------------------------------
!
! MASNUM model Part 2: Main Program  
! -- Originally developed by Lei HAN (c) 2012.August
! -- Code standardized by Lei HAN and Zhanpeng ZHUANG in September, 2013
! -- Description:
!
!-------------------------------------------------------------------------------------------------
!*************************************************************************************************

Program masnum_circ_model_main

  use def_common_mod
  
  implicit none
  
  include 'netcdf.inc'
                
  integer ii,jj,kk,kij(3)
  
    ii=260
    jj=55
    kk=15
  
  !-----------------------------------------------------------------------------------------------
  ! Model initialization
  !-----------------------------------------------------------------------------------------------
  call init_model_sub
  print*,'init done'

  call cpu_time(cpu_start)                  ! Program execution timing
      
  open(741,file='psi_annual.txt')
  open(742,file='psi_instant.txt')

  !===============================================================================================
  ! Internal mode integration
  !===============================================================================================
  
  Internal_Loop: do iint = 1 , iend

    if(ymdhms(1)==2015 .and. iint>1) call save_monthly

    if (IF_record_series==1) call record_time_series_sub

    !---------------------------------------------------------------------------------------------
    ! Update the time-related variables
    !---------------------------------------------------------------------------------------------
    call upd_time_sub           
    
    if (mode/=2) then
    !---------------------------------------------------------------------------------------------

    if (Baropg_method==1) then
      call IPG_linear_sub              ! Update the internal pressure gradient force
    elseif (Baropg_method==0) then
      allocate(ipg_x(ks,im,jm), ipg_y(ks,im,jm))
      ipg_x = ZERO
      ipg_y = ZERO
      deallocate(rho)
    endif    
    
    
    !---------------------------------------------------------------------------------------------
    ! Update horizontal mixing coefficient
    !---------------------------------------------------------------------------------------------
    call Am_Smagorinsky_sub
   
    !---------------------------------------------------------------------------------------------
    ! Update advection and diffusion for the next iteration
    !---------------------------------------------------------------------------------------------
    call diff_uv_sub                        ! Diffusion terms
    call adv_uv_sub                         ! Horizontal advection term
    !---------------------------------------------------------------------------------------------
    
    call vertical_mean_sub                  ! Vertically integrate the 3-D arrays
    
    endif
        
    !---------------------------------------------------------------------------------------------
      
    call compute_psi_sub
    annual_psi = annual_psi + psi  
      
      
    call ext_mode_sub                       ! Execute the external mode
    if (mode==2) cycle

    call upd_uv_sub

    !---------------------------------------------------------------------------------------------
    ! Adjust uf vf to adapt to external-mode flux
    !---------------------------------------------------------------------------------------------
    call adjust_uv_sub
     
    !---------------------------------------------------------------------------------------------
    ! Update the vertical velocity wf
    !---------------------------------------------------------------------------------------------
    call upd_w_sub
    
    !---------------------------------------------------------------------------------------------
    ! Solve the tracers' prognostic equations : Temperature (T)
    !---------------------------------------------------------------------------------------------
    call bsm_mix
    call upd_ts_sub
    
    !---------------------------------------------------------------------------------------------
    ! Update the density
    !---------------------------------------------------------------------------------------------
    call compute_dens_POM98
       
    !---------------------------------------------------------------------------------------------
    ! POM turbulence scheme
    !---------------------------------------------------------------------------------------------
    call advq_POM98('q2')
    call advq_POM98('q2l')

    call bv_mix
    call bti_mix

    call profq_POM98

    if(iint == 1) call save_data
    if(ymdhms(1)==2015 .and. (ymdhms(3)==5 .and. ymdhms(4)==12) ) call save_data

    call finalize_INTloop_sub
    
  end do Internal_Loop
  
  ! output annual psi and instantaneous psi
 	annual_psi = annual_psi / iend * 1e-6
 	psi = psi * 1e-6
	write(741,'(<im>e15.5)') ((annual_psi(i,j),i=1,im),j=1,jm)
	close(741)
	write(742,'(<im>e15.5)') ((psi(i,j),i=1,im),j=1,jm)
	close(742)


  !-----------------------------------------------------------------------------------------------
  ! Program execution time estimation
  !-----------------------------------------------------------------------------------------------
  call cpu_time(cpu_stop)
  
  print*, 'Main program executed', (cpu_stop - cpu_start) / 3600, ' hours'

contains
  
  !*************************************************************************************************
  !-----------------------------------------------------------------------------------------------
  !
  ! MASNUM model Part 3: Subroutines to be called in the main program
  ! --- Components:
  !     1) init_model_sub
  !
  !-----------------------------------------------------------------------------------------------
  !*************************************************************************************************
  
  !-----------------------------------------------------------------------------------------------
  ! 1)
  ! This init_model_sub subroutine defines the computational meshes, set initial values and 
  ! prepare the external data to spin up the circulation model.
  !
  ! -- Originally developed by Lei HAN (c) August, 2012
  ! -- Code standardized by Lei HAN and Zhanpeng ZHUANG in September, 2013
  !-----------------------------------------------------------------------------------------------
  subroutine init_model_sub
    implicit none
    
    !---------------------------------------------------------------------------------------------
    ! computational grid
    !---------------------------------------------------------------------------------------------
    integer :: kdz(12) 
    real(kind=precision) :: delz, z1
    character(len=2) level, mon_char

    integer maxcr(2),ki
  
    data kdz/1,1,2,4,8,16,32,64,128,256,512,1024/


    !---------------------------------------------------------------------------------------------
    ! allocate the momery spaces of global variable arraies
    !---------------------------------------------------------------------------------------------
    call init_common_array                             
                                                        

    !---------------------------------------------------------------------------------------------
    ! Generate mask and define open boundary points
    !---------------------------------------------------------------------------------------------
    call mask_define_sub
  
    print*,'0.1',maxval(H_c),minval(H_c)
    
    
    

    !---------------------------------------------------------------------------------------------
    ! Set values of H_u and H_v with H_c
    !---------------------------------------------------------------------------------------------
    do j=1,jm
      do i=2,im
        H_u(i,j) = 0.5 * (H_c(i,j)+H_c(i-1,j))
                                            ! Calculate H_u
      enddo
    enddo
    do j=2,jm
      do i=1,im
        H_v(i,j) = 0.5 * (H_c(i,j)+H_c(i,j-1))
                                            ! Calculate H_v
      enddo
    enddo
    H_v(:,1) = H_v(:,2)                     ! Define marginal H_v
    if (regional_run) then                  ! Define marginal H_u
      H_u(1,:) = H_u(2,:)                   ! ...
    else                                    ! ...
      H_u(1,:) = 0.5 * (H_c(1,:)+H_c(im-2,:))
    endif                                   ! ...
    
    !---------------------------------------------------------------------------------------------
    ! Set initial values of the dynamic depths
    !---------------------------------------------------------------------------------------------
    hb_inc = H_c
    hb_inu = H_u
    hb_inv = H_v
    hf_inc = H_c
    hf_inu = H_u
    hf_inv = H_v
    hm_inc = H_c
    hm_inu = H_u
    hm_inv = H_v
    
    !---------------------------------------------------------------------------------------------
    ! Define the vertical coordinate
    !---------------------------------------------------------------------------------------------
    kl1=8
    kl2=kb-2  
    z(1)=0.0  
    do k=2,kl1
      z(k)=z(k-1)+float(kdz(k-1))           ! kl1 shold not exceed 13
                                            ! exponential increase within k=1:kl1
    enddo
    delz=z(kl1)-z(kl1-1)                    ! =kdz(kl1-1)
    do k=kl1+1,kl2
      z(k)=z(k-1)+delz                      ! linear increase within k=kl1+1:kl2
    enddo
    do k=kl2+1,ks
      dz(k)=float(kdz(ks-k+1))*delz/float(kdz(ks-kl2))
      z(k)=z(k-1)+dz(k)                     ! exponential decrease within k=kl2+1:ks
    end do
    do k=1,ks
      z(k)=-z(k)/z(ks)                      ! Normalize
    end do

    zt = z
    do k=1,ks
      z(k) = zt(ks+1-k)                     ! Reverse the vertical
    enddo
    do k=1,ks-1
      zz(k)=0.5*(z(k)+z(k+1))               ! Central levels of the vertical grid
    end do
    zz(ks) = 0.0                            

    dz = 0.0                                ! Set vertical grid interval
    do k=1,ks-1                             ! ...
      dz(k)=z(k+1)-z(k)                     ! ...
      dzz(k)=zz(k+1)-zz(k)                  ! ...
    end do
    dz(ks)=0.0
    dzz(ks)=0.0

    !---------------------------------------------------------------------------------------------
    ! Compute the dx & dy of T-cell  
    !---------------------------------------------------------------------------------------------
    do j=1,jm
      do i=1,im
        dy_t(i,j) = Ra / inv_dlat * pi / 180.
        dx_t(i,j) = Ra / inv_dlon * pi / 180. * cosd( vlat(j) ) 
      enddo
    enddo

    !---------------------------------------------------------------------------------------------
    ! Compute the dx & dy of U-cell & V-cell
    !---------------------------------------------------------------------------------------------
    do j=1,jm
      do i=2,im
        dx_u(i,j) = ( dx_t(i,j) + dx_t(i-1,j) ) / 2.0
        dy_u(i,j) = ( dy_t(i,j) + dy_t(i-1,j) ) / 2.0
      enddo
    enddo
    do j=2,jm
      do i=1,im
        dx_v(i,j) = ( dx_t(i,j) + dx_t(i,j-1) ) / 2.0
        dy_v(i,j) = ( dy_t(i,j) + dy_t(i,j-1) ) / 2.0
      enddo
    enddo
    
    dx_v(:,1) = dx_v(:,2)                   ! Define marginal dx_v
    dy_v(:,1) = dy_v(:,2)                   ! Define marginal dy_v
    if (regional_run) then                  ! Define marginal dx_u and dy_u
        dx_u(1,:) = dx_u(2,:)
        dy_u(1,:) = dy_u(2,:)
    else
        dx_u(1,:) = dx_u(imm1,:)            ! Cyclic dx_u in global run
        dy_u(1,:) = dy_u(imm1,:)            ! Cyclic dy_u in global run
    endif
    
    do j=2,jm
      do i=1,im
        dx_u_down(i,j) = ( dx_u(i,j) + dx_u(i,j-1) ) / 2.0
                                            ! dx of the lower side of U-cell 
      enddo
    enddo    
    do j=1,jm
      do i=2,im
        dy_v_left(i,j) = ( dy_v(i,j) + dy_v(i-1,j) ) / 2.0
                                            ! dy of the left side of V-cell 
      enddo
    enddo
    do i=1,im
      dx_u_down(i,1) = dx_u_down(i,2)         ! Define marginal dx_u_down
    enddo
    if (regional_run) then
      do j=1,jm
        dy_v_left(1,j) = dy_v_left(2,j)       ! Define marginal dy_v_down
      enddo
    else   
      do j=1,jm 
        dy_v_left(1,j) = dy_v_left(imm1,j)    ! Cyclic dy_v_left in global run
      enddo
    endif
       
    !---------------------------------------------------------------------------------------------
    ! Compute areas of T-cell, U-cell, and V-cell
    !---------------------------------------------------------------------------------------------
    art = dx_t*dy_t
    aru = dx_u*dy_u
    arv = dx_v*dy_v

    !-------------------------------------------------------------------------------------------
    ! Compute the local Courant number (Stability check)
    !-------------------------------------------------------------------------------------------
    do i=1,im
      do j=1,jm
        if (fsm(i,j)==1) then
          courant(i,j) = sqrt(GRAV*hb_inc(i,j)) * dte / dx_t(i,j)
        endif
      enddo
    enddo
    
!    call test_interp_field_2d('courant')
    
    maxcr = maxloc(courant)
    i = maxcr(1)
    j = maxcr(2)
    print*, '---------------------------------------------------------------------------'
    print*, 'Max-Courant Number = ', maxval(courant)
    print*, 'Max-Courant Number index :', i, j
    print*, 'Max-Courant Number [lon lat H]:', vlon(i),vlat(j), H_c(i,j)
    print*, 'dte, isplit = ',dte, isplit
    print*, '---------------------------------------------------------------------------'
    if (maxval(courant)>=1) then
    	print*,'Stability WARNING(*) : MAX Courant number exceeds unity. Please check'
      stop
    endif
        
    !---------------------------------------------------------------------------------------------
    ! Compute the Coriolis parameter centered at T-cell, U-cell and V-cell
    !---------------------------------------------------------------------------------------------
    do j=1,jm
      do i=1,im
        Cor_t(i,j) = 2.0*2.*PI/86400.*sind(vlat(j))
        Cor_u(i,j) = Cor_t(i,j)
        Cor_v(i,j) = 2.0*2.*PI/86400.*sind(vlat_v(j))
      enddo
    enddo
        
    !---------------------------------------------------------------------------------------------
    ! Set the bottom friction coefficient
    !---------------------------------------------------------------------------------------------
    do j=1,jmm1
      do i=1,imm1
        cbc(i,j)=( KAPPA/log((1.e0+zz(1))*H_c(i,j)/z0b) )**2  
        cbc(i,j)=max(cbcmin,cbc(i,j))    
        cbc(i,j)=min(cbcmax,cbc(i,j))    
      end do
    end do

    print*,'********month=',ymdhms(2)
!    call get_3d_field('global_res','tb')
!    call get_3d_field('global_res','sb')
    call get_3d_field('PHC-annual','t_an')
    call get_3d_field('PHC-annual','s_an')

    do k=1,ks
      do j=1,jm
        do i=1,im
          tb(k,i,j) = tclim(k,i,j)
          sb(k,i,j) = sclim(k,i,j)
        enddo
      enddo
    enddo

!=================================================    

    deallocate(rho)
    allocate(tf(ks,im,jm),sf(ks,im,jm))
    tf = tb
    sf = sb
    call compute_dens_POM98
    deallocate(tf,sf)
    
    call get_3d_field('glo_res_MAS','t')
    call get_3d_field('glo_res_MAS','s')
    call get_3d_field('glo_res_MAS','u')
    call get_3d_field('glo_res_MAS','v')

    do i=1,im
      do j=1,100
        if(vb(ksm1,i,j) > 0.5) then
          vb(1:ksm1,i,j) = 0.5
        endif
        if(vb(ksm1,i,j) < -0.5) then
          vb(1:ksm1,i,j) = -0.5
        endif
      enddo
    enddo

    print*,'read init T',maxval(tb),minval(tb)
    print*,'read init s',maxval(sb),minval(sb)
    print*,'read init u',maxval(ub),minval(ub)
    print*,'read init v',maxval(vb),minval(vb)

    !---------------------------------------------------------------------------------------------
    ! Load the initial static background fields (T & S) 
    !---------------------------------------------------------------------------------------------
!    call get_3d_field('global_res','rmean')   
    call get_3d_field('PHC-annual','rhoz')    
    
    
    print*,'5.5', maxval(rclim),minval(rclim)
    
    
!    call test_interp_field_3d('rclim')

!  	call test_interp_field_3d('tb',zero)
!   	call test_interp_field_3d('sb',zero)
!   	call test_interp_field_3d('rho')
   	!stop

    !---------------------------------------------------------------------------------------------
    ! Load Open Boundary Data for regional run
    !---------------------------------------------------------------------------------------------
    !if (regional_run) call load_OBC_data_sub
    if (regional_run) call upd_open_bound
    print*,'8'
    !---------------------------------------------------------------------------------------------
    ! Load wave mixing data (Bv)
    !---------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------
    ! Load restart data
    !---------------------------------------------------------------------------------------------
    if (cold_start==0) then
    	
      open(1,file=trim(adjustl(hot_file_load)),form='unformatted')
      read(1)                                        &
         iint0, days_now                             &
        ,ua,va,uab,vab,etb_c,zetb_c,wubot,wvbot      &
        ,tb,sb,ub,vb,wf,q2b,q2lb,Km,Kh,Kq,rho        
      close(1)
      
      print*, '------------------------------------------------------------------------'
      print*, 'This is a MASNUM restart run from YY/MM/DD/HH/MM/SS: ', datestr(days_now)
      print*, '------------------------------------------------------------------------'
      
      !-------------------------------------------------------------------------------------------
      ! Redefine the total integration steps
      !-------------------------------------------------------------------------------------------
      if (Use_restart_date==1) then
      	days_begin = days_now
        ymdhms = datevec(days_now)
        iend = ( datenum( ymdhms_finish ) - days_begin ) * 86400.0 / dti 
        print*,'restart file date used.'
      else
      	days_now = days_begin
      	print*,'restart file date neglected, starting date in the control parameter file is used'
      endif
            
      !-------------------------------------------------------------------------------------------
      ! Update the surface elevations and dynamic depths
      !-------------------------------------------------------------------------------------------
      etf_c  = etb_c
      hb_inc = H_c  + zetb_c
      do j=2,jm
        do i=2,im
          hb_inu(i,j) = ( hb_inc(i,j)+hb_inc(i-1,j) ) / 2.0
          hb_inv(i,j) = ( hb_inc(i,j)+hb_inc(i,j-1) ) / 2.0
        enddo
      enddo
      if (regional_run) then
        hb_inu(1,:) = hb_inu(2,:)           ! Define marginal hb
        hb_inu(:,1) = hb_inu(:,2)
        hb_inv(1,:) = hb_inv(2,:)
        hb_inv(:,1) = hb_inv(:,2)
      else
        hb_inu(1,:) = hb_inu(imm1,:)        ! Cyclic in global run
        hb_inv(1,:) = hb_inv(imm1,:)        
      endif
      hf_inc = hb_inc                       ! Set initial value of hf
      hf_inu = hb_inu
      hf_inv = hb_inv
      
    endif

    !---------------------------------------------------------------------------------------------
    ! Open file to record time series of the total kinetic energy
    !---------------------------------------------------------------------------------------------
    if (cold_start==1) then
    	open(2012,file='total_energy_masnumV1.8.txt',status='replace')
    elseif (cold_start==0) then
    	open(2012,file='total_energy_masnumV1.8.txt',access='append')
    endif

    print*,'9'

    print*,'get_2d_field call================'

!    call get_2d_field('dqdsst')
!  	call test_interp_field_2d('dqdsst')

!    call get_2D_field('tmi_amsre_SST','analysed_sst')
!  	call test_interp_field_2d('sst_OI',-32768.)
    
    !stop
    
    return
  end subroutine init_model_sub
!=================================================================================================

!-----------------------------------------------------------------------------------------------
! 2)
! This ext_mode_sub subroutine solve the external equations.
!
! -- Originally developed by Lei HAN (c) August, 2012
! -- Code standardized by Lei HAN and Zhanpeng ZHUANG in September, 2013
!-----------------------------------------------------------------------------------------------
  subroutine ext_mode_sub
    implicit none 
    
    real(kind=precision), allocatable :: epg_ext_x(:,:), epg_ext_y(:,:)
    real(kind=precision), allocatable :: fcor_ua(:,:), fcor_va(:,:)
    real(kind=precision), allocatable :: flux_ua(:,:), flux_va(:,:)
    real(kind=precision), allocatable :: split_discrep_x(:,:), split_discrep_y(:,:)
    
    !---------------------------------------------------------------------------------------------
    ! allocate the momery spaces
    !---------------------------------------------------------------------------------------------
    allocate(epg_ext_x(im,jm), epg_ext_y(im,jm))
    allocate(fcor_ua(im,jm), fcor_va(im,jm))
    allocate(flux_ua(im,jm), flux_va(im,jm))
    allocate(split_discrep_x(im,jm), split_discrep_y(im,jm))
      epg_ext_x = ZERO
      epg_ext_y = ZERO
      fcor_ua = ZERO
      fcor_va = ZERO
      flux_ua = ZERO
      flux_va = ZERO
      split_discrep_x = ZERO
      split_discrep_y = ZERO
    
    !---------------------------------------------------------------------------------------------
    ! Compute 2-D advection and diffusion terms
    !---------------------------------------------------------------------------------------------    
    call adv_diff_ext_sub

    !---------------------------------------------------------------------------------------------
    ! Compute the difference between the vertically averaged 3-D adv&diff terms and the 2-D
    ! adv&diff terms. This quantity is used in the external loop by assuming that the difference 
    ! between the 2D and vertical-mean-3D adv&diff terms is held invariant within the entire 
    ! internal time step
    !---------------------------------------------------------------------------------------------        
    do j=1,jm
      do i=1,im
        split_discrep_x(i,j) = int_vm_u(i,j) - ( diff_ua(i,j) - adv_ua(i,j) )
        split_discrep_y(i,j) = int_vm_v(i,j) - ( diff_va(i,j) - adv_va(i,j) )
      enddo
    enddo

    
    et_mean = 0.                            ! Clear the arrays
    utf     = 0.
    vtf     = 0.
    
    !---------------------------------------------------------------------------------------------
    ! Set tidal forcing open bounadry condition
    !---------------------------------------------------------------------------------------------    
    if (load_OB_tide) call obcond_sub('SET:et-tide')
    
    !=============================================================================================
    ! Start integration of external mode
    !=============================================================================================
    
    External_Loop: DO iext = 1, isplit      
      
      !-------------------------------------------------------------------------------------------
      ! Compute the volume flux around the T-cell
      !-------------------------------------------------------------------------------------------      
      do j=1,jm
        do i=1,im
          flux_ua(i,j) = uab(i,j) * dy_u(i,j)
          flux_va(i,j) = vab(i,j) * dx_v(i,j) 
        enddo
      enddo
      
      utf = utf + uab                       ! The time-integrated volume flux within the entire
      vtf = vtf + vab                       ! internal time step, for adjusting the 3-D velocity
      
      !-------------------------------------------------------------------------------------------
      ! Update the external surface elevation (eta)
      !-------------------------------------------------------------------------------------------      
      do j=2,jmm1
        do i=2,imm1
          if (fsm(i,j).ne.1) cycle
          etf_c(i,j) = flux_ua(i,j) - flux_ua(i+1,j)   &
                     + flux_va(i,j) - flux_va(i,j+1)
          etf_c(i,j) = etb_c(i,j)                      &
                     + etf_c(i,j) * dte / art(i,j)
        enddo
      enddo
      
      etf_c = etf_c * fsm
      
      if (regional_run==0) call obcond_sub('CYC:et')
      
      call smooth_spatial_sub('etf')

      if (regional_run==0) call obcond_sub('CYC:et')

                                            ! Cyclic eta OBC in the global run 
            
      !-------------------------------------------------------------------------------------------
      ! Compute the external velocity components (ua & va)
      !-------------------------------------------------------------------------------------------      
      do j=2,jm
        do i=2,im
          ua(i,j) = uab(i,j)                                             &
                  / (H_u(i,j)*4. + etf_c(i,j)+etf_c(i-1,j)               &
                                 + etb_c(i,j)+etb_c(i-1,j) ) * 4.
          va(i,j) = vab(i,j)                                             &
                  / (H_v(i,j)*4. + etf_c(i,j)+etf_c(i,j-1)               &
                                 + etb_c(i,j)+etb_c(i,j-1) ) * 4.
        enddo
      enddo
      
      !-------------------------------------------------------------------------------------------
      ! OBC on ua & va (both regional & global)
      !-------------------------------------------------------------------------------------------      
      if (regional_run) then
        do n=1,numboundary
          ib = ij_obc(1,n)
          jb = ij_obc(2,n)     
          ua(ib,jb) = uab(ib,jb) / H_u(ib,jb)
          va(ib,jb) = vab(ib,jb) / H_v(ib,jb)
        enddo
        ua = ua * dum01
        va = va * dvm01
        etf_c( 1,:) = etf_c(   2,:)
        etf_c(im,:) = etf_c(imm1,:)
        etf_c(:, 1) = etf_c(:,   2)
        etf_c(:,jm) = etf_c(:,jmm1)        
        
        etf_c = etf_c * fsm
  
        do i=1,im
          ua(i,jm) = uab(i,jm) / H_u(i,jm)
          va(i,jm) = vab(i,jm) / H_v(i,jm)
          ua(i,2 ) = uab(i,2 ) / H_u(i,2 )
          va(i,2 ) = vab(i,2 ) / H_v(i,2 )
        enddo
        do j=1,jm
          ua(im,j) = uab(im,j) / H_u(im,j)
          va(im,j) = vab(im,j) / H_v(im,j)
          ua(2 ,j) = uab(2 ,j) / H_u(2 ,j)
          va(2 ,j) = vab(2 ,j) / H_v(2 ,j)
        enddo
        ua = ua * dum01
        va = va * dvm01
      else
        call obcond_sub('CYC:ua')
      end if
      
      et_mean = et_mean + etf_c             ! The time-averaged surface elevation throughout the 
                                            ! entire internal time step
                                            
      !-------------------------------------------------------------------------------------------
      ! Update the 2D advection and diffusion terms (not on every external time step)
      !-------------------------------------------------------------------------------------------      
      if ( mod(iext,steps_upd_advdiff)==0 ) then
        call adv_diff_ext_sub
      endif
      
      !-------------------------------------------------------------------------------------------
      ! Compute uaf (ua*h_u)
      !-------------------------------------------------------------------------------------------            
      do j=2,jmm1
        do i=2,imm1
          if (dum(i,j)/=1) cycle
          epg_ext_x(i,j) = - GRAV * ( etf_c(i,j)-etf_c(i-1,j) )       &
              * ( H_c(i,j)+H_c(i-1,j)+etf_c(i,j)+etf_c(i-1,j) ) / 2.0 / dx_u(i,j)      
          fcor_ua(i,j) = 0.25 * &
             ( Cor_t(i  ,j) * (va(i  ,j)+va(i  ,j+1)) * (H_c(i  ,j)+etf_c(i  ,j))   &
             + Cor_t(i-1,j) * (va(i-1,j)+va(i-1,j+1)) * (H_c(i-1,j)+etf_c(i-1,j)) )
          uaf(i,j) = epg_ext_x(i,j) + fcor_ua(i,j)                           &
                   + diff_ua(i,j)   - adv_ua(i,j) + split_discrep_x(i,j)     &
                   + ( wubot(i,j)-wusurf(i,j) )
          !if (i==166.and.j==50) then
          !	print*,'epg_ext',epg_ext_x(i,j),fcor_ua(i,j)
          !	print*,'diff_ua',diff_ua(i,j),adv_ua(i,j)
          !	print*,'split_discrep',split_discrep_x(i,j)
          !	print*,'wusurf wubot',wusurf(i,j),wubot(i,j)
          !	print*,'total:',uaf(i,j)
          !  	print*,'--------------'
          !endif
          uaf(i,j) = uab(i,j) + dte * uaf(i,j)
        enddo
      enddo
      
      !-------------------------------------------------------------------------------------------
      ! Compute vaf (va*h_v)
      !-------------------------------------------------------------------------------------------            
      do j=2,jmm1
        do i=2,imm1
          if (dvm(i,j)/=1) cycle
          epg_ext_y(i,j) = - GRAV * ( etf_c(i,j)-etf_c(i,j-1) )       &
           * ( H_c(i,j)+H_c(i,j-1)+etf_c(i,j)+etf_c(i,j-1) ) / 2.0 / dy_v(i,j)
          fcor_va(i,j) = - 0.25 * &
             ( Cor_t(i,j  ) * (ua(i,j  )+ua(i+1,j  )) * (H_c(i,j  )+etf_c(i,j  )) &
             + Cor_t(i,j-1) * (ua(i,j-1)+ua(i+1,j-1)) * (H_c(i,j-1)+etf_c(i,j-1)) )
          vaf(i,j) = epg_ext_y(i,j) + fcor_va(i,j)                           &
                   + diff_va(i,j)   - adv_va(i,j) + split_discrep_y(i,j)     &
                   + ( wvbot(i,j)-wvsurf(i,j) )
          vaf(i,j) = vab(i,j) + dte * vaf(i,j)
        enddo
      enddo  
      
      !-------------------------------------------------------------------------------------------
      ! OBC on uaf & vaf
      !-------------------------------------------------------------------------------------------            
      if (regional_run) then
        if (load_OB_tide) then
          call obcond_sub('tide:uaH')
        else
          call obcond_sub('FLA:uaH')
        endif     
      else
        call obcond_sub('CYC:uaH')
      endif
      
      !-------------------------------------------------------------------------------------------
      ! Dual spatial smoothing : uaf and vaf, respectively
      !-------------------------------------------------------------------------------------------            
      !call smooth_spatial_sub('uaf')
      !call smooth_spatial_sub('vaf')
                  
      !-------------------------------------------------------------------------------------------
      ! Step forward for the next iteration
      !-------------------------------------------------------------------------------------------                        
      etb_c = etf_c
      uab   = uaf
      vab   = vaf
    
    end do External_Loop
    
    et_mean = et_mean / isplit              ! For computing the internal EPG

    !---------------------------------------------------------------------------------------------
    ! Update the dynamic depths of the internal mode 
    !---------------------------------------------------------------------------------------------                            
    zetf_c = etf_c
    hf_inc = H_c + zetf_c
    do j=2,jm
      do i=2,im
        hf_inu(i,j) = ( hf_inc(i,j)+hf_inc(i-1,j) ) / 2.0
        hf_inv(i,j) = ( hf_inc(i,j)+hf_inc(i,j-1) ) / 2.0
      enddo
    enddo
    
    !---------------------------------------------------------------------------------------------
    ! Set the marginal value of the dynamic depths 
    !---------------------------------------------------------------------------------------------                                
    if (regional_run) then
      hf_inu(1,:) = hf_inu(2,:)
      hf_inu(:,1) = hf_inu(:,2)
      hf_inv(1,:) = hf_inv(2,:)
      hf_inv(:,1) = hf_inv(:,2)
    else
      hf_inu(1,:) = hf_inu(imm1,:)
      hf_inv(1,:) = hf_inv(imm1,:)
    endif
    
    !---------------------------------------------------------------------------------------------
    ! Set the dynamic depths of the internal mode at the middle time step 
    !---------------------------------------------------------------------------------------------                                    
    hm_inc = ( hb_inc+hf_inc ) / 2.0
    hm_inu = ( hb_inu+hf_inu ) / 2.0
    hm_inv = ( hb_inv+hf_inv ) / 2.0
    
    !---------------------------------------------------------------------------------------------
    ! Deallocate the momery spaces
    !---------------------------------------------------------------------------------------------
    deallocate(epg_ext_x, epg_ext_y)
    deallocate(fcor_ua, fcor_va)
    deallocate(flux_ua, flux_va)
    deallocate(split_discrep_x, split_discrep_y)
    
    return
  end subroutine ext_mode_sub


!-------------------------------------------------------------------------------------------------
! 3)
! This mask_define_sub subroutine define the mask.
!
! -- Originally developed by Lei HAN (c) August, 2012
! -- Code standardized by Zhanpeng ZHUANG in September, 2013
!-------------------------------------------------------------------------------------------------
  subroutine mask_define_sub
    implicit none
  
    integer :: count, ite
    logical existed
    real(kind=precision) :: x1, x2, y1, y2
    
    !---------------------------------------------------------------------------------------------
    ! Limit the water depth range before smoothing
    !---------------------------------------------------------------------------------------------
    do j=1,jm
      do i=1,im
        
        !------------------------------------------------------------------------------------------
        ! Land Depth
        !------------------------------------------------------------------------------------------
        if (H_c(i,j)<Land_Depth) H_c(i,j) = Land_Depth
        
        !------------------------------------------------------------------------------------------
        ! Mininum water depth
        !------------------------------------------------------------------------------------------
        if (H_c(i,j)>Land_Depth .and. H_c(i,j)<Min_Depth) H_c(i,j) = Min_Depth
        
        !------------------------------------------------------------------------------------------
        ! Maximum water depth
        !------------------------------------------------------------------------------------------
        if (H_c(i,j)>Max_Depth) H_c(i,j) = Max_Depth
        
      enddo
    enddo    

    !---------------------------------------------------------------------------------------------
    ! Generate mask(fsm) with H_c
    !---------------------------------------------------------------------------------------------
    fsm = 1.0
    do j=1,jm
      do i=1,im
        if(H_c(i,j)<=Land_Depth) fsm(i,j)=0          ! Land point
      end do
    end do


    !---------------------------------------------------------------------------------------------
    ! Set the information of the open boundary points
    !---------------------------------------------------------------------------------------------
    if (regional_run) then
    	
      do ite=1,2
    	  
        !-------------------------------------------------------------------------------------------
        ! Defaulted square domain (i.e., West Pac run)
        !-------------------------------------------------------------------------------------------
        count = 1
        ij_obc = ZERO
        bctype = ZERO
        
        !-------------------------------------------------------------------------------------------
        ! North
        !-------------------------------------------------------------------------------------------
        do i=1,im
          j=jm
          if (fsm(i,j)==1) then
            ij_obc(1,count)=i
            ij_obc(2,count)=j
            bctype(count)=0
            count = count + 1
          endif
        enddo

        !-------------------------------------------------------------------------------------------
        ! East
        !-------------------------------------------------------------------------------------------
        do j=1,jm
          i=im
          if (fsm(i,j)==1) then
            ij_obc(1,count)=i
            ij_obc(2,count)=j
            bctype(count)=3
            count = count + 1
          endif
        enddo    
        
        !-------------------------------------------------------------------------------------------
        ! South
        !-------------------------------------------------------------------------------------------
        do i=1,im
          j=1
          if (fsm(i,j)==1) then
            ij_obc(1,count)=i
            ij_obc(2,count)=j
            bctype(count)=6
            count = count + 1
          endif
        enddo
        
        !-------------------------------------------------------------------------------------------
        ! West
        !-------------------------------------------------------------------------------------------
        do j=1,jm
          i=1
          if (fsm(i,j)==1) then
            ij_obc(1,count)=i
            ij_obc(2,count)=j
            bctype(count)=9
            count = count + 1
          endif
        enddo
        
        numboundary = count - 1
                
        !-------------------------------------------------------------------------------------------
        ! Define the accurate coordinates of open boundary points of U,V,T, respectively
        !-------------------------------------------------------------------------------------------
        if (allocated(reg_xy)) deallocate(reg_xy)
        allocate( reg_xy(6,numboundary) )

        do n=1,numboundary
          ib = ij_obc(1,n)
          jb = ij_obc(2,n)
         
          if ( bctype(n) <= 3 ) then         ! right boundary: north and east
          
            reg_xy(:,n) =               &
              (/ vlon(ib)  , vlat(jb),  &
                 vlon_u(ib), vlat(jb),  &
                 vlon(ib)  , vlat_v(jb) /)
         
          elseif ( bctype(n) == 6 ) then     ! left boundary: south
                   
            reg_xy(:,n) =                 &
              (/ vlon(ib)  , vlat(jb),    &
                 vlon_u(ib), vlat(jb),    &
                 vlon(ib)  , vlat_v(jb+1) /)
        
          elseif ( bctype(n) == 9 ) then     ! left boundary: west
                   
            reg_xy(:,n) =                 &
              (/ vlon(ib)    , vlat(jb),  &
                 vlon_u(ib+1), vlat(jb),  &
                 vlon(ib)    , vlat_v(jb) /)
         
          endif
         
        enddo
        
        !------------------------------------------------------------------------------------------
        ! mark the point unable to be interpolated from the more expanded coaser model land
        !------------------------------------------------------------------------------------------
        if (ite==1) call upd_open_bound(1)   
        
      enddo
        
    endif

    !---------------------------------------------------------------------------------------------
    ! Smooth the topography and limit the slope of the topo
    !---------------------------------------------------------------------------------------------
    do l=1,Num_topo_smooth
    	call smooth_spatial_sub('H_c')
    	print*,'TOPO smooth ',Num_topo_smooth,' iterations.'
    enddo

    call slpmax    

    do j=1,jm
      do i=1,im        
        H_c(i,j) = max( H_c(i,j), Land_Depth ) 
      enddo
    enddo

    !---------------------------------------------------------------------------------------------
    ! Cyclic condition for Global Run
    !---------------------------------------------------------------------------------------------
    if (.not.regional_run) then
      print*,'Topo cyclic'
      do j=1,jm
        H_c( 1,j) = H_c(imm1,j)
        H_c(im,j) = H_c(   2,j)
      enddo
    endif
    
    
!    call test_interp_field_2d('fsm',one)
!    call test_interp_field_2d('H_c',one)
        
    !---------------------------------------------------------------------------------------------
    ! Set velocity mask: dum&dvm according to fsm
    !---------------------------------------------------------------------------------------------
    do j=1,jm
      do i=2,im
        dum(i,j) = fsm(i,j) * fsm(i-1,j)
      enddo
    enddo

    do j=2,jm
      do i=1,im
        dvm(i,j) = fsm(i,j) * fsm(i,j-1)
      enddo
    enddo

    !---------------------------------------------------------------------------------------------
    ! Further define MASKs for the open boundary points
    !---------------------------------------------------------------------------------------------
    
    if (regional_run) then  
      do j=1,jm    
        dum(1,j) = dum(2,j)
      enddo
      do i=1,im
        dvm(i,1) = dvm(i,2)
      enddo
      
      do n=1,numboundary
        ib = ij_obc(1,n)
        jb = ij_obc(2,n)
        if (bctype(n)==6) then              ! South
          dvm(ib,jb+1) = - dvm(ib,jb+1)
        elseif (bctype(n)==9) then          ! West
          dum(ib+1,jb) = - dum(ib+1,jb)
        else                                ! North and East
          dum(ib,jb) = - dum(ib,jb)
          dvm(ib,jb) = - dvm(ib,jb)
        endif
      enddo
      
    else
      
      do j=1,jm
        dum(im,j) = dum(   2,j)               ! cyclic
        dum( 1,j) = dum(imm1,j)  
        dvm(im,j) = dvm(   2,j)               ! cyclic
        dvm( 1,j) = dvm(imm1,j)
      enddo
      
    endif
    
    dum01 = abs(dum)
    dvm01 = abs(dvm)

    return
  end subroutine mask_define_sub
!=================================================================================================
  
!-------------------------------------------------------------------------------------------------
! 4)
! This upd_date_sub subroutine update date/time : yyyy/mm/dd/hh and Sun angle.
!
! -- Originally developed by Lei HAN (c) August, 2012
! -- Code standardized by Zhanpeng ZHUANG in September, 2013
!-------------------------------------------------------------------------------------------------
  subroutine upd_time_sub
    implicit none

    real(kind=precision) :: hour_passed, hour_total, hour_left
    character(len=100) :: work_dir
    integer :: kij(3)
  
    !---------------------------------------------------------------------------------------------
    ! Check the elapsed time and estimate the remaining time
    !---------------------------------------------------------------------------------------------
    if (mod(iint,10)==1 .and. 0) then               
      call getcwd(work_dir)
      print*,'Current work dir = ',trim(work_dir)
      print*, 'Current YYYY-MM-DD-HH-MM-SS: ', datestr(days_now)
      call cpu_time(cpu_realtime)           
      hour_passed = (cpu_realtime - cpu_start) / 3600.
      hour_total  = hour_passed * iend / iint
      hour_left   = hour_total - hour_passed
      print*, '------ Time passed away : ', hour_passed, ' hours'
      print*, '------ Estimated time to finish : ', hour_left, ' hours'
    endif
    
    !-------------------------------------------------------------------------------------------
    ! Update the time
    !-------------------------------------------------------------------------------------------
    days_now = days_now + dti/86400.0
    
    ymdhms_prev = ymdhms
    date_char_prev = date_char
    
    ymdhms    = datevec(days_now)
    date_char = datestr(days_now)
    
    model_days = days_now - days_begin      ! model's running days since cold/hot start
    
    yearly_day = days_now - datenum(ymdhms(1),1,1,0,0,0)
    
    !---------------------------------------------------------------------------------------------
    ! Set ramp for cold start
    !---------------------------------------------------------------------------------------------
    if (cold_start==1) then
      ramp = min( model_days / period, 1.0 )
      if (regional_run) then
        ramp_obc = min( model_days / period_obc, 1.0)
        ramp_tide = ramp
      endif
    endif
    
!    print*, '-- iint = ', iint, '  Total:', iend, '  Complete = ', float(iint)/float(iend)*100,'%  --'
!    print*, '   model days= ', int(model_days*100.)/100.
!    print*, '   date:', ymdhms(1:5)
    if(mod(iint,10)==1) write(*,'(15G15.5)') iint,ymdhms(1:5),sum(ub)
    
    
    !---------------------------------------------------------------------------------------------
    ! Monthly updated quantities from external files
    !---------------------------------------------------------------------------------------------
    if (ymdhms(1)<2020) then

      if (iint==1 .or. ymdhms(2)/=ymdhms_prev(2)) then
 	  
        !---------------------------------------------------------------------------------------------
        ! Update monthly climatological T and S
        !---------------------------------------------------------------------------------------------
        call get_2d_field('pomsurface','sst','clim')
        print*,'sst_dq',maxval(sst_dq),minval(sst_dq)
!        call test_interp_field_2d('sst_dq')
    	    	
        !---------------------------------------------------------------------------------------------
        ! Load surface fluxes
        !---------------------------------------------------------------------------------------------
        call get_2d_field('pomsurface','dqdsst','clim')
        print*,'dqdsst',maxval(dqdsst),minval(dqdsst)
!        call test_interp_field_2d('dqdsst')
    	
        call get_2d_field('pomsurface','wndx','clim')
        call get_2d_field('pomsurface','wndy','clim')
        call get_2d_field('pomsurface','netheat','clim')

      
!        call test_interp_field_2d('uflx',zero)
!        call test_interp_field_2d('vflx',zero)
      
        netheat = - netheat
      
!        call test_interp_field_2d('netheat',zero)
          
!        deallocate(uswrf, ulwrf, dswrf, dlwrf, shtfl, lhtfl)
      
      endif

      !---------------------------------------------------------------------------------------------
      ! Update surface heat and momentum fluxes: wtsurf, wusurf, wvsurf
      !---------------------------------------------------------------------------------------------
      do i=1,im
        do j=1,jm
          wtsurf(i,j) = ( netheat(i,j) + dqdsst(i,j)*(sst_dq(i,j)-tb(ksm1,i,j)) ) / 4.1876e6 
        enddo
      enddo

    else

      if (iint==1 .or. ymdhms(3)/=ymdhms_prev(3)) then

        call get_2d_field('NCEP2:surface_fluxes','uflx','4xdaily')
        call get_2d_field('NCEP2:surface_fluxes','vflx','4xdaily')
        
        call get_2d_field('NCEP2:surface_fluxes','uswrf','4xdaily')
        call get_2d_field('NCEP2:surface_fluxes','ulwrf','4xdaily')
        call get_2d_field('NCEP2:surface_fluxes','dswrf','4xdaily')
        call get_2d_field('NCEP2:surface_fluxes','dlwrf','4xdaily')
        call get_2d_field('NCEP2:surface_fluxes','shtfl','4xdaily')
        call get_2d_field('NCEP2:surface_fluxes','lhtfl','4xdaily')

        netheat = uswrf + ulwrf - dswrf - dlwrf + shtfl + lhtfl
          
        deallocate(uswrf, ulwrf, dswrf, dlwrf, shtfl, lhtfl)

      endif

      !---------------------------------------------------------------------------------------------
      ! Update surface heat and momentum fluxes: wtsurf, wusurf, wvsurf
      !---------------------------------------------------------------------------------------------
      do j=1,jm
        do i=1,im
          wtsurf(i,j) = ( netheat(i,j) ) / 4.1876e6 
        enddo
      enddo

    endif

    do j=1,jm
      do i=1,im
        wusurf(i,j) = -uflx(i,j) * ramp/1025. * dum01(i,j)
        wvsurf(i,j) = -vflx(i,j) * ramp/1025. * dvm01(i,j)
      enddo
    enddo

   
    !---------------------------------------------------------------------------------------------
    ! Update Open boundary data
    !---------------------------------------------------------------------------------------------
    if ( regional_run== 1 ) then
      if ( iint==1 .or. ymdhms(3)/=ymdhms_prev(3) ) then
!      	if (mod(int(yearly_day),5)==0) then
      	  call upd_open_bound
!        endif
      endif
    endif    
   
    !---------------------------------------------------------------------------------------------
    ! If turn on Changjiang diluting water
    !---------------------------------------------------------------------------------------------
    if (ChangJiang_diluting) then
       time_chj_q = chj_q(ymdhms(2)) * ramp_tide
    endif
    
    !---------------------------------------------------------------------------------------------
    ! Output instantaneous data to files : daily field
    !---------------------------------------------------------------------------------------------
    if ( IF_export_instant .and. nth_field==0 ) then
      nth_field = 1
    endif

    if ( nth_field/=0 .and. ymdhms(4)/=ymdhms_prev(4)   &
        .and. mod(ymdhms(4),interval_export_instant)==0 ) then
      print*,'  **** outputing instant file No.', nth_field
!      call output_instant_data
    endif

    !---------------------------------------------------------------------------------------------
    ! Output data to files : monthly climatological data
    !---------------------------------------------------------------------------------------------
    if ( IF_export_clim .and. ymdhms(3)/=ymdhms_prev(3)   &
        .and. iint/=1 ) then
!      call output_clim_data
    endif
    
    if ( IF_EAKF_argo .and. ymdhms(3)/=ymdhms_prev(3) ) then
!    	call output_instant_ts
    endif
            
    return
  end subroutine upd_time_sub
!=================================================================================================

!-------------------------------------------------------------------------------------------------
! 5)
! This IPG_linear_sub subroutine calculate baroclinic pressure gradient force using POM.
!
! -- Originally developed by Lei HAN (c) August, 2012
! -- Code standardized by Zhanpeng ZHUANG in September, 2013
!-------------------------------------------------------------------------------------------------
  subroutine IPG_linear_sub
    implicit none

    real(kind=precision), allocatable :: drhox(:,:,:), drhoy(:,:,:)
    integer :: ki
    
        
    rho = rho - rclim

    !---------------------------------------------------------------------------------------------
    ! Calculate x-component of baroclinic pressure gradient:
    !---------------------------------------------------------------------------------------------
    allocate(ipg_x(ks,im,jm), drhox(ks,im,jm))
      ipg_x = ZERO
      drhox = ZERO

    do j=2,jmm1
      do i=2,imm1
        drhox(1,i,j) = GRAV * (-zz(ks-1)) * hb_inu(i,j)         &
                            * (rho(ks-1,i,j) - rho(ks-1,i-1,j))
      enddo
    enddo

    do j=2,jmm1
      do i=2,imm1
        if (dum(i,j)/=1) cycle
        do k=2,kbm1    
          drhox(k,i,j) = drhox(k-1,i,j)                                       &
            + GRAV * 0.50 * ( zz(ks-k+1)-zz(ks-k) )                           &
                   * hb_inu(i,j) * ( rho(ks-k  ,i,j) - rho(ks-k  ,i-1,j)      &
                                   + rho(ks-k+1,i,j) - rho(ks-k+1,i-1,j) )    &
            + GRAV * 0.25 * ( zz(ks-k+1)+zz(ks-k) )                           &
                   * ( hb_inc(i,j) - hb_inc(i-1,j) )                          &
                   * ( rho(ks-k  ,i,j) + rho(ks-k  ,i-1,j)                    &
                     - rho(ks-k+1,i,j) - rho(ks-k+1,i-1,j) )
        enddo
      enddo
    end do

    do j=2,jmm1
      do i=2,imm1
        if (dum(i,j)/=1) cycle
        do k=1,kbm1
          drhox(k,i,j) = hb_inu(i,j) * drhox(k,i,j) * dy_u(i,j)
        enddo
      end do
    end do

    do j=1,jm
      do i=1,im
        if (dum(i,j)/=1) cycle
        do k=1,ksm1
          ki=kb-k
          ipg_x(k,i,j) = -drhox(ki,i,j) / rho_ref / aru(i,j) * ramp
        enddo
      enddo
    enddo

    
    if (0) then
    print*,'ipg-x',ipg_x(1,i+1,j)/ramp,ramp,aru(i,j),hb_inu(i,j)
    print*,'t',tb(:,i,j)
    print*,'t',tb(:,i+1,j)
    print*,'s',sb(:,i,j)
    print*,'s',sb(:,i+1,j)
    print*,'rho(165,50)',rho(:,i,j)
    print*,'rho(166,50)',rho(:,i+1,j)
    print*,'rclim(165,50)',rclim(:,i,j)
    print*,'rclim(166,50)',rclim(:,i+1,j)
    print*,'zz',zz
    print*,'H(165,50)',H_c(165:166,50)
   ! stop
  endif
    

    
    deallocate(drhox)

    !---------------------------------------------------------------------------------------------
    ! Calculate y-component of baroclinic pressure gradient:
    !---------------------------------------------------------------------------------------------
    allocate(ipg_y(ks,im,jm), drhoy(ks,im,jm))
      ipg_y = ZERO
      drhoy = ZERO
    
    do j=2,jmm1
      do i=2,imm1
        drhoy(1,i,j) = GRAV * (-zz(ks-1)) * hb_inv(i,j)         &
                            * (rho(ks-1,i,j) - rho(ks-1,i,j-1))
      enddo
    enddo

    do j=2,jmm1
      do i=2,imm1  
        if (dvm(i,j)/=1) cycle
        do k=2,kbm1
          drhoy(k,i,j) = drhoy(k-1,i,j)                                       &
            + GRAV * 0.50 * ( zz(ks-k+1)-zz(ks-k) )                           &
                   * hb_inv(i,j) * ( rho(ks-k  ,i,j) - rho(ks-k  ,i,j-1)      &
                                   + rho(ks-k+1,i,j) - rho(ks-k+1,i,j-1) )    &
            + GRAV * 0.25 * ( zz(ks-k+1)+zz(ks-k) )                           &
                   * ( hb_inc(i,j) - hb_inc(i,j-1) )                          &
                   * ( rho(ks-k  ,i,j) + rho(ks-k  ,i,j-1)                    &
                     - rho(ks-k+1,i,j) - rho(ks-k+1,i,j-1) )
        end do
      end do
    end do

    do j=2,jmm1
      do i=2,imm1
        if (dvm(i,j)/=1) cycle
        do k=1,kbm1
          drhoy(k,i,j) = hb_inv(i,j) * drhoy(k,i,j) * dx_v(i,j)
        enddo
      end do
    end do

    do j=1,jm
      do i=1,im
        if (dvm(i,j)/=1) cycle
        do k=1,ksm1
          ki=kb-k
          ipg_y(k,i,j) = -drhoy(ki,i,j) / rho_ref / arv(i,j) * ramp
        enddo
      enddo
    enddo
    

    deallocate(drhoy)
    deallocate(rho)
    
    
    return
  end subroutine IPG_linear_sub
!=================================================================================================

!-------------------------------------------------------------------------------------------------
! 6)
! This subroutine update IPGs with cubic spline method.
!
! -- Originally developed by Lei HAN (c) August, 2012
! -- Code standardized by Zhanpeng ZHUANG in September, 2013
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
! 7)
! This Am_Smagorinsky_sub subroutine calculate horizontal mixing coefficient with Smagorinsky.
!
! -- Originally developed by Lei HAN (c) August, 2012
! -- Code standardized by Zhanpeng ZHUANG in September, 2013
!-------------------------------------------------------------------------------------------------
  subroutine Am_Smagorinsky_sub
    
    implicit none 

    real(kind=precision) :: dudx, dudy, dvdx, dvdy
    real(kind=precision) :: uv_grad
  
    allocate(Am(ks,im,jm))
      Am = ZERO

    do j=2,jmm1
      do i=2,imm1
        if (fsm(i,j)/=1) cycle
        do k=1,ksm1
          dudx = (ub(k,i+1,j) - ub(k,i,j)) / dx_t(i,j)
          dvdy = (vb(k,i,j+1) - vb(k,i,j)) / dy_t(i,j)
          dudy = (ub(k,i,j+1) + ub(k,i+1,j+1) - ub(k,i,j-1) - ub(k,i+1,j-1))            &
               / (dy_v(i,j) + dy_v(i,j+1))
          dvdx = (vb(k,i+1,j) + vb(k,i+1,j+1) - vb(k,i-1,j) - vb(k,i-1,j+1))            &
               / (dx_u(i,j) + dx_u(i+1,j))
     
          uv_grad = dvdx + dudy
          uv_grad = dudx*dudx + dvdy*dvdy + uv_grad*uv_grad*0.125
          uv_grad = sqrt(uv_grad) 
          
          Am(k,i,j) = horcon * art(i,j) * uv_grad
        enddo
      enddo
    enddo

    if (regional_run) then
      !-------------------------------------------------------------------------------------------
      ! Define marginal area
      !-------------------------------------------------------------------------------------------
      do j=1,jm
        do k=1,ks
          Am(k, 1, j) = Am(k,   2, j)
          Am(k,im, j) = Am(k,imm1, j)
        enddo
      enddo
                
      do i=1,im
        do k=1,ks
          Am(k, i,  1) = Am(k, i,   2)
          Am(k, i, jm) = Am(k, i,jmm1)
        enddo
      enddo
      
    else
      !-------------------------------------------------------------------------------------------
      ! Limit the minimum value to 3000 according to Xia's POM code
      !-------------------------------------------------------------------------------------------
      Am = max(Am,3000.)
      !-------------------------------------------------------------------------------------------
      ! Cyclic : Quasi Globe
      !-------------------------------------------------------------------------------------------
      do j=1,jm
        do k=1,ks
          Am(k, 1, j) = Am(k,imm1, j)
          Am(k,im, j) = Am(k,   2, j)
        enddo
      enddo
      
    endif
    
    return
  end subroutine Am_Smagorinsky_sub
!=================================================================================================

!-------------------------------------------------------------------------------------------------
! 8)
! This diff_uv_sub subroutine calculate momentum diffusion(x component).
!
! -- Originally developed by Lei HAN (c) August, 2012
! -- Code standardized by Zhanpeng ZHUANG in September, 2013
!-------------------------------------------------------------------------------------------------
  subroutine diff_uv_sub
    
    implicit none 

    real(kind=precision), allocatable :: diff_u(:,:,:), diff_v(:,:,:)
    real(kind=precision), allocatable :: diff_x(:,:,:), diff_y(:,:,:)
    real(kind=precision), allocatable :: px(:,:,:)
    real(kind=precision) :: p1, p2

    !---------------------------------------------------------------------------------------------
    ! allocate the momery spaces
    !---------------------------------------------------------------------------------------------
    allocate(adv_diff_u(ks,im,jm), adv_diff_v(ks,im,jm))
    allocate(diff_x(ks,im,jm), diff_y(ks,im,jm))
    allocate(px(ks,im,jm))
      adv_diff_u = ZERO
      adv_diff_v = ZERO
      diff_x = ZERO
      diff_y = ZERO
      px =ZERO

    !---------------------------------------------------------------------------------------------  
    ! Compute the diffusion term of U-component : diff_u in the unit of acceleratoin (m/s**2)
    !---------------------------------------------------------------------------------------------  

    do j=2,jm
      do i=2,im
        if (fsm(i,j) + fsm(i-1,j) + fsm(i,j-1) + fsm(i-1,j-1)==0) cycle
        do k=1,ks
          px(k,i,j) = ( Am(k,i  ,j) + Am(k,i-1,j-1)                               &
                      + Am(k,i-1,j) + Am(k,i  ,j-1) + 4 * umol) * 0.25            &
                    * ( (ub(k,i,j) - ub(k,i,j-1)) / dy_v_left(i,j)                &
                      + (vb(k,i,j) - vb(k,i-1,j)) / dx_u_down(i,j) )
          diff_y(k,i,j) = px(k,i,j) * dx_u_down(i,j)                              &
                        * (hf_inu(i,j) + hf_inu(i,j-1)) / 2.0
        enddo
      enddo
    enddo

    do j=2,jm
      do i=2,im
        if (fsm(i-1,j)==0) cycle
        do k=1,ks
          p1 = ( Am(k,i-1,j) + umol ) * 2.0                                &
             * ( ub(k,i,j) - ub(k,i-1,j) ) / dx_t(i-1,j)
          diff_x(k,i,j) = p1 * dy_t(i-1,j) * hf_inc(i-1,j)
        enddo
      enddo
    enddo

    allocate(diff_u(ks,im,jm))
      diff_u = ZERO

    do j=2,jmm1
      do i=2,imm1
        if (dum(i,j)/=1) cycle
        do k=1,ksm1
          diff_u(k,i,j) = - diff_x(k,i,j) + diff_x(k,i+1,j  )       &
                          - diff_y(k,i,j) + diff_y(k,i  ,j+1)
          diff_u(k,i,j) = diff_u(k,i,j) / aru(i,j)
        enddo
      enddo
    enddo

    adv_diff_u = diff_u
    deallocate(diff_u)
    
    !---------------------------------------------------------------------------------------------  
    ! Compute the diffusion term of V-component : diff_v in the unit of acceleratoin (m/s**2)
    !---------------------------------------------------------------------------------------------  
    diff_x = 0.
    diff_y = 0.

    do j=2,jm
      do i=2,im
        if (fsm(i,j) + fsm(i-1,j) + fsm(i,j-1) + fsm(i-1,j-1)==0) cycle     
        do k=1,ks
          diff_x(k,i,j) = px(k,i,j) * dy_v_left(i,j)                              &
                        * (hf_inu(i,j) + hf_inu(i,j-1)) / 2.0
        enddo
      enddo
    enddo

    do j=2,jm
      do i=2,im
        if (fsm(i,j-1)==0) cycle
        do k=1,ks
          p2 = ( Am(k,i,j-1) + umol ) * 2.0                                &
             * ( vb(k,i,j) - vb(k,i,j-1) ) / dy_t(i,j-1)
          diff_y(k,i,j) = p2 * dx_t(i,j-1) * hf_inc(i,j-1)
        enddo
      enddo
    enddo
    
    allocate(diff_v(ks,im,jm))
      diff_v = ZERO

    do j=2,jmm1
      do i=2,imm1
        if (dvm(i,j)==0) cycle
        do k=1,ksm1
          diff_v(k,i,j) = - diff_x(k,i,j) + diff_x(k,i+1,j  )                     &
                          - diff_y(k,i,j) + diff_y(k,i  ,j+1)  
          diff_v(k,i,j) = diff_v(k,i,j) / arv(i,j)
        enddo
      enddo
    enddo

    adv_diff_v = diff_v
    deallocate(diff_v)

    !---------------------------------------------------------------------------------------------
    ! Deallocate the momery spaces
    !---------------------------------------------------------------------------------------------
    deallocate(diff_x, diff_y)
    deallocate(px)
        
    return
  end subroutine diff_uv_sub
!=================================================================================================

!-------------------------------------------------------------------------------------------------
! 9)
! This adv_uv_sub subroutine calculate momentum upwind advection.
!
! -- Originally developed by Lei HAN (c) August, 2012
! -- Code standardized by Zhanpeng ZHUANG in September, 2013
!-------------------------------------------------------------------------------------------------
  subroutine adv_uv_sub
    
    implicit none 

    real(kind=precision), allocatable :: adv_u(:,:,:), adv_v(:,:,:)
    real(kind=precision), allocatable :: adv_x(:,:,:), adv_y(:,:,:)
    real(kind=precision), allocatable :: adv_z(:)
    real(kind=precision) :: uface, vface, wface
    integer :: iup, jup, kup

    !---------------------------------------------------------------------------------------------
    ! allocate the momery spaces
    !---------------------------------------------------------------------------------------------
    allocate(adv_u(ks,im,jm))
    allocate(adv_x(ks,im,jm), adv_y(ks,im,jm), adv_z(ks))
      adv_u = ZERO 
      adv_x = ZERO
      adv_y = ZERO
      adv_z = ZERO

    !---------------------------------------------------------------------------------------------
    ! Compute the advection term of U-component : adv_u in the unit of acceleratoin (m/s**2)
    !---------------------------------------------------------------------------------------------
    do j=2,jm
      do i=2,im
        if (fsm(i-1,j)==0) cycle
        do k=1,ksm1  
          uface = ( ub(k,i,j)*hf_inu(i,j) + ub(k,i-1,j)*hf_inu(i-1,j) ) * 0.5
          adv_x(k,i,j) = uface * ( ub(k,i,j)+ub(k,i-1,j) ) * 0.5 * dy_t(i-1,j)
        enddo
      enddo
    enddo
    
    do j=2,jm
      do i=2,im
        if (dvm01(i,j) + dvm01(i-1,j)==0) cycle
        if (dum01(i,j) + dum01(i,j-1)==0) cycle
        do k=1,ksm1  
          vface = ( vb(k,i,j)*hf_inv(i,j) + vb(k,i-1,j)*hf_inv(i-1,j) ) * 0.5 
          adv_y(k,i,j) = vface * (ub(k,i,j) + ub(k,i,j-1)) * 0.5 * dx_u_down(i,j)
        enddo
      enddo
    enddo

    do j=2,jmm1
      do i=2,imm1
        if (dum(i,j)==0) cycle
        do k=1,ksm1
          adv_u(k,i,j) = - adv_x(k,i,j) + adv_x(k,i+1,j  )    &
                         - adv_y(k,i,j) + adv_y(k,i  ,j+1)
          adv_u(k,i,j) = adv_u(k,i,j) / aru(i,j)
        enddo
      enddo
    enddo
    
    !---------------------------------------------------------------------------------------------
    ! Compute the vertical component of momentum convection flux 
    !---------------------------------------------------------------------------------------------
    do j=2,jmm1
      do i=2,imm1
        if (dum(i,j)/=1) cycle
        do k=2,ksm1
          adv_z(k) = (wf(k,i,j) + wf(k,i-1,j))              &
                   * (ub(k,i,j) + ub(k-1,i,j)) * 0.25
        enddo
        do k=1,ksm1
          adv_u(k,i,j) = adv_u(k,i,j)                       & 
                       + ( -adv_z(k) + adv_z(k+1) ) / dz(k)
        enddo        
      enddo
    enddo
        
    adv_diff_u = adv_diff_u - adv_u
    deallocate(adv_u)
    
    !---------------------------------------------------------------------------------------------
    ! Compute the advection term of V-component : adv_v in the unit of acceleratoin (m/s**2)  
    !---------------------------------------------------------------------------------------------
    adv_x = ZERO
    adv_y = ZERO
    adv_z = ZERO
    
    allocate(adv_v(ks,im,jm))
      adv_v = ZERO
    
    do j=2,jm
      do i=2,im
        if (dum01(i,j) + dum01(i,j-1)==0) cycle
        if (dvm01(i,j) + dvm01(i-1,j)==0) cycle
        do k=1,ksm1
          uface = ( ub(k,i,j)*hf_inu(i,j) + ub(k,i,j-1)*hf_inu(i,j-1) ) * 0.5
          adv_x(k,i,j) = uface * ( vb(k,i,j)+vb(k,i-1,j) ) * 0.5 * dy_v_left(i,j)
        enddo
      enddo
    enddo

    do j=2,jm
      do i=2,im
        if (fsm(i,j-1)==0) cycle
        do k=1,ksm1
          vface = ( vb(k,i,j)*hf_inv(i,j) + vb(k,i,j-1)*hf_inv(i,j-1) ) * 0.5 
          adv_y(k,i,j) = vface * ( vb(k,i,j)+vb(k,i,j-1) ) * 0.5 * dx_t(i,j-1)
        enddo
      enddo
    enddo

    do j=2,jmm1
      do i=2,imm1
        if (dvm(i,j)/=1) cycle
        do k=1,ksm1
          adv_v(k,i,j) = - adv_x(k,i,j) + adv_x(k,i+1,j  )    &
                         - adv_y(k,i,j) + adv_y(k,i  ,j+1)
          adv_v(k,i,j) = adv_v(k,i,j) / arv(i,j)  
        enddo
      enddo
    enddo  

    !---------------------------------------------------------------------------------------------
    ! Compute the vertical component of momentum convection flux 
    !---------------------------------------------------------------------------------------------
    do j=2,jmm1
      do i=2,imm1
        if (dvm(i,j)/=1) cycle
        do k=2,ksm1
          adv_z(k) = (wf(k,i,j) + wf(k,i,j-1))              &
                   * (vb(k,i,j) + vb(k-1,i,j)) * 0.25 
        enddo
        do k=1,ksm1
          adv_v(k,i,j) = adv_v(k,i,j)                       &
                       + ( -adv_z(k) + adv_z(k+1) ) / dz(k) 
        enddo
      enddo
    enddo
    
    adv_diff_v = adv_diff_v - adv_v
    deallocate(adv_v)
    
    !---------------------------------------------------------------------------------------------
    ! Deallocate the momery spaces
    !---------------------------------------------------------------------------------------------
    deallocate(adv_x, adv_y, adv_z)
    deallocate(wf)
    
    return
  end subroutine adv_uv_sub
!=================================================================================================

!-------------------------------------------------------------------------------------------------
! 10)
! This vertical_mean_sub subroutine calculate vertical integration of 3d variables.
!
! -- Originally developed by Lei HAN (c) August, 2012
! -- Code standardized by Zhanpeng ZHUANG in September, 2013
!-------------------------------------------------------------------------------------------------
  subroutine vertical_mean_sub
    
    implicit none 

    real(kind=precision), allocatable :: adv_diff_vm_u(:,:), adv_diff_vm_v(:,:)
    real(kind=precision), allocatable :: IPG_vm_x(:,:), IPG_vm_y(:,:)

    !---------------------------------------------------------------------------------------------
    ! allocate the momery spaces
    !---------------------------------------------------------------------------------------------
    allocate(adv_diff_vm_u(im,jm), adv_diff_vm_v(im,jm))
    allocate(IPG_vm_x(im,jm), IPG_vm_y(im,jm))
      adv_diff_vm_u = ZERO
      adv_diff_vm_v = ZERO
      IPG_vm_x = ZERO
      IPG_vm_y = ZERO

    !---------------------------------------------------------------------------------------------
    ! Vertically integrate Am
    !---------------------------------------------------------------------------------------------
    if (iint/=1) then
    Am_vm = 0.
    do j=1,jm
      do i=1,im
        do k=1,ksm1
          Am_vm(i,j) = Am_vm(i,j) + Am(k,i,j) * dz(k)
        enddo
      enddo
    enddo
    endif

    !---------------------------------------------------------------------------------------------
    ! Vertically integrate IPG
    !---------------------------------------------------------------------------------------------
    do j=1,jm
      do i=1,im
        do k=1,ksm1
          IPG_vm_x(i,j) = IPG_vm_x(i,j) + ipg_x(k,i,j)*dz(k)
          IPG_vm_y(i,j) = IPG_vm_y(i,j) + ipg_y(k,i,j)*dz(k)
        enddo
      enddo
    enddo

    !---------------------------------------------------------------------------------------------
    ! Vertical average of 3D advection and diffusion terms:
    !---------------------------------------------------------------------------------------------  
    do j=1,jm
      do i=1,im
        do k=1,ksm1
          adv_diff_vm_u(i,j) = adv_diff_vm_u(i,j) + adv_diff_u(k,i,j) * dz(k)
          adv_diff_vm_v(i,j) = adv_diff_vm_v(i,j) + adv_diff_v(k,i,j) * dz(k)
        enddo
      enddo
    enddo
  
    !---------------------------------------------------------------------------------------------
    ! The total vertically-integrated 3-D forces applied in the external momentum equations
    !---------------------------------------------------------------------------------------------
    int_vm_u =  IPG_vm_x + adv_diff_vm_u
    int_vm_v =  IPG_vm_y + adv_diff_vm_v
        
    !---------------------------------------------------------------------------------------------
    ! Deallocate the momery spaces
    !---------------------------------------------------------------------------------------------
    deallocate(adv_diff_vm_u, adv_diff_vm_v)
    deallocate(IPG_vm_x, IPG_vm_y)
    
    return
  end subroutine vertical_mean_sub
!=================================================================================================

!-------------------------------------------------------------------------------------------------
! 11)
! This adv_diff_ext_sub subroutine calculate 2D advection & diffusion terms.
!
! -- Originally developed by Lei HAN (c) August, 2012
! -- Code standardized by Zhanpeng ZHUANG in September, 2013
!-------------------------------------------------------------------------------------------------
  subroutine adv_diff_ext_sub
    
    implicit none
  
    integer :: iup, jup
    real(kind=precision) :: uface, vface
    real(kind=precision), allocatable :: adv_uax(:,:), adv_uay(:,:)
    real(kind=precision), allocatable :: adv_vax(:,:), adv_vay(:,:)
  
    real(kind=precision), allocatable :: p11_ua(:,:), p21_ua(:,:)
    real(kind=precision), allocatable :: p12_va(:,:), p22_va(:,:)
    real(kind=precision) :: p21_temp
    real(kind=precision), allocatable :: hf_e(:,:)

    !---------------------------------------------------------------------------------------------
    ! allocate the momery spaces
    !---------------------------------------------------------------------------------------------
    allocate(adv_uax(im,jm), adv_uay(im,jm))
    allocate(adv_vax(im,jm), adv_vay(im,jm))
    allocate(p11_ua(im,jm), p21_ua(im,jm))
    allocate(p12_va(im,jm), p22_va(im,jm))
    allocate(hf_e(im,jm))
      adv_uax = ZERO
      adv_uay = ZERO
      adv_vax = ZERO
      adv_vay = ZERO
      p11_ua  = ZERO
      p21_ua  = ZERO
      p12_va  = ZERO
      p22_va  = ZERO
      hf_e    = ZERO

    do j=2,jm
      do i=2,im
        hf_e(i,j) = 0.25 * ( H_c(i,j) + H_c(i-1,j) + H_c(i,j-1) + H_c(i-1,j-1)             &
                         + etf_c(i,j) + etf_c(i-1,j) + etf_c(i,j-1) + etf_c(i-1,j-1))
      enddo
    enddo
    
    do i=1,im
      hf_e(i,1) = hf_e(i,2)
    enddo
    if (regional_run) then
      do j=1,jm
        hf_e(1,j) = hf_e(2,j)
      enddo
    else
      do j=1,jm
        hf_e(1,j) = hf_e(imm1,j)  ! cyclic
      enddo
    endif

    !---------------------------------------------------------------------------------------------
    ! Compute 2-D advective flux
    !---------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------  
    ! Compute viscous stress tensor and momentum diffusion
    !---------------------------------------------------------------------------------------------
    do j=2,jmm1
      do i=2,im
        if (fsm(i-1,j)==0) cycle
        uface = (uab(i,j) + uab(i-1,j)) * 0.5
        adv_uax(i,j) = uface * (ua(i,j) + ua(i-1,j)) * 0.5 * dy_t(i-1,j)        
        p11_ua(i,j) = ( Am_vm(i-1,j) + umol ) * 2.0 * ( ua(i,j) - ua(i-1,j) )            &
                    * ( H_c(i-1,j) + etf_c(i-1,j) ) * ( dy_t(i-1,j)/dx_t(i-1,j) )
      enddo
    enddo

    do j=2,jm
      do i=2,imm1
        if (fsm(i,j-1)==0) cycle
        vface = (vab(i,j) + vab(i,j-1)) * 0.5
        adv_vay(i,j) = vface * (va(i,j) + va(i,j-1)) * 0.5 * dx_t(i,j-1)        
        p22_va(i,j) = ( Am_vm(i,j-1) + umol ) * 2.0 * ( va(i,j) - va(i,j-1) )               &
                    * ( H_c(i,j-1) + etf_c(i,j-1) ) * ( dx_t(i,j-1) / dy_t(i,j-1))
      enddo
    enddo

    do j=2,jm
      do i=2,im
        if (dvm01(i,j)+dvm01(i-1,j)==0) cycle
        if (dum01(i,j)+dum01(i,j-1)==0) cycle
        vface = (vab(i,j) + vab(i-1,j)) * 0.5
        uface = (uab(i,j) + uab(i,j-1)) * 0.5
        adv_uay(i,j) = vface * ( ua(i,j) + ua(i,j-1) ) * 0.5 * dx_u_down(i,j)
        adv_vax(i,j) = uface * ( va(i,j) + va(i-1,j) ) * 0.5 * dy_v_left(i,j)
      enddo
    enddo

    do j=2,jm
      do i=2,im
        if (fsm(i,j) + fsm(i-1,j) + fsm(i,j-1) + fsm(i-1,j-1)==0) cycle
        p21_temp = ( Am_vm(i  ,j) + Am_vm(i-1,j-1)                         &
                   + Am_vm(i-1,j) + Am_vm(i  ,j-1) + 4*umol ) * 0.25       &
                 * ( (ua(i,j) - ua(i,j-1)) / dy_v_left(i,j)              &
                   + (va(i,j) - va(i-1,j)) / dx_u_down(i,j) )            &
                 * hf_e(i,j)
        p21_ua(i,j) = p21_temp * dx_u_down(i,j)
        p12_va(i,j) = p21_temp * dy_v_left(i,j)
      enddo
    enddo

    do j=2,jmm1
      do i=2,imm1
        adv_ua(i,j)  = - adv_uax(i,j) + adv_uax(i+1,j)    &
                       - adv_uay(i,j) + adv_uay(i,j+1)
        adv_va(i,j)  = - adv_vax(i,j) + adv_vax(i+1,j)    &
                       - adv_vay(i,j) + adv_vay(i,j+1)    
        diff_ua(i,j) = - p11_ua(i,j) + p11_ua(i+1,j)      &
                       - p21_ua(i,j) + p21_ua(i,j+1)
        diff_va(i,j) = - p12_va(i,j) + p12_va(i+1,j)      &
                       - p22_va(i,j) + p22_va(i,j+1)
        adv_ua(i,j)  = adv_ua(i,j)  / aru(i,j)
        adv_va(i,j)  = adv_va(i,j)  / arv(i,j)
        diff_ua(i,j) = diff_ua(i,j) / aru(i,j)
        diff_va(i,j) = diff_va(i,j) / arv(i,j)
      enddo
    enddo
    
    !---------------------------------------------------------------------------------------------
    ! 2D mode Bottom momentum flux
    !---------------------------------------------------------------------------------------------
    if (mode==2) then
      do j=2,jmm1
        do i=2,imm1     
          wubot(i,j) = -0.5 * ua(i,j) * (cbc(i,j) + cbc(i-1,j))                   &  
                     * sqrt( ua(i,j)**2                                           &
                           + (0.25*( va(i  ,j) + va(i  ,j+1)                      &
                                   + va(i-1,j) + va(i-1,j+1)))**2 )
          
          wvbot(i,j) = -0.5 * va(i,j) * (cbc(i,j) + cbc(i,j-1))                   &
                     * sqrt( va(i,j)**2                                           &
                           + (0.25*( ua(i,j  ) + ua(i+1,j  )                      &
                                   + ua(i,j-1) + ua(i+1,j-1)))**2 ) 
        enddo
      enddo
    endif

    !---------------------------------------------------------------------------------------------
    ! Deallocate the momery spaces
    !---------------------------------------------------------------------------------------------
    deallocate(adv_uax, adv_uay)
    deallocate(adv_vax, adv_vay)
    deallocate(p11_ua, p21_ua)
    deallocate(p12_va, p22_va)
    deallocate(hf_e)
    
    return
  end subroutine adv_diff_ext_sub
!=================================================================================================

!-------------------------------------------------------------------------------------------------
! 12)
! This solve_tridiag_sub subroutine use Thomas algorithm solve tridiagonal system.
!
! -- Originally developed by Lei HAN (c) August, 2012
! -- Code standardized by Zhanpeng ZHUANG in September, 2013
!-------------------------------------------------------------------------------------------------
  subroutine solve_tridiag_sub(var_name)
    
    implicit none 

    character(len=1), intent(in) :: var_name
    real(kind=precision), pointer :: wxsurf(:,:)
    real(kind=precision), pointer :: xb(:,:,:)
    real(kind=precision), pointer :: xf(:,:,:)

    real(kind=precision), pointer :: mask(:,:)
    real(kind=precision), pointer :: hb_inx(:,:),hm_inx(:,:)

    real(kind=precision), allocatable :: Rh(:,:)

    !---------------------------------------------------------------------------------------------
    ! Thomas algorithm variables
    !---------------------------------------------------------------------------------------------
    real(kind=precision) :: a_temp
    real(kind=precision), allocatable :: a_u(:), b_u(:), c_u(:), Fu(:)
    real(kind=precision), allocatable :: aa0(:), cc0(:)
  
    real(kind=precision), allocatable :: rhs(:,:,:), Kx_temp(:,:,:)
    real(kind=precision), allocatable :: bot_temp(:,:)

    !---------------------------------------------------------------------------------------------
    ! allocate the momery spaces
    !---------------------------------------------------------------------------------------------
    allocate(Rh(im,jm))
    allocate(a_u(ks-2), b_u(ksm1), c_u(ks-2), Fu(ksm1))
    allocate(aa0(ksm1), cc0(ks-2))
    allocate(rhs(ks,im,jm), Kx_temp(ks,im,jm))
    allocate(bot_temp(im,jm))
      Rh  = ZERO
      a_u = ZERO
      b_u = ZERO
      c_u = ZERO
      Fu  = ZERO
      aa0 = ZERO
      cc0 = ZERO
      rhs = ZERO
      Kx_temp  = ZERO
      bot_temp = ZERO

    !---------------------------------------------------------------------------------------------
    ! Define respective relating variables
    !---------------------------------------------------------------------------------------------
    bot_temp = 0.

    if (var_name=='u') then
      wxsurf => wusurf                      ! Added by Zhuang Zhanpeng
      mask => dum
      hb_inx => hb_inu
      hm_inx => hm_inu
      xb => ub
      allocate(uf(ks,im,jm))
      xf => uf
      Rh = hf_inu / hb_inu
      do j=2,jmm1
        do i=2,imm1
          do k=1,ksm1
            rhs(k,i,j) = adv_diff_u(k,i,j) + fcor_u(k,i,j)         &
                       + epg_int_x(i,j)    + ipg_x(k,i,j)
            
            if(0.and.i==166.and.j==50) then
            	print*,'--------------'
            	print*,'u,k=',i,j,k,ub(k,i,j),ua(i,j)
            	print*,'adv_u  fcor_u',adv_diff_u(k,i,j),fcor_u(k,i,j)
            	print*,'epg_int ipg_x',epg_int_x(i,j),ipg_x(k,i,j)
            	print*,'rhs',rhs(k,i,j)
            endif
            
            Kx_temp(k,i,j) = Km(k+1,i,j) + Km(k+1,i-1,j) + 2.0 * umol
          enddo
        enddo
      enddo
      bot_temp = fric_ub * dti / hb_inu / dz(1)
      
    elseif (var_name=='v') then
      wxsurf => wvsurf                      ! Added by Zhuang Zhanpeng
      mask => dvm
      hb_inx => hb_inv
      hm_inx => hm_inv
      xb => vb
      allocate(vf(ks,im,jm))
      xf => vf
      Rh = hf_inv / hb_inv
      do j=2,jmm1
        do i=2,imm1
          do k=1,ksm1
            rhs(k,i,j) = adv_diff_v(k,i,j) + fcor_v(k,i,j)         &
                       + epg_int_y(i,j)    + ipg_y(k,i,j)            
            Kx_temp(k,i,j) = Km(k+1,i,j) + Km(k+1,i,j-1) + 2.0 * umol
          enddo
        enddo
      enddo
      bot_temp = fric_vb * dti / hb_inv / dz(1)

    else                                    ! For tracer variables such as T&S
      mask => fsm
      hb_inx => hb_inc
      hm_inx => hm_inc
      Rh = hf_inc / hb_inc
      if (var_name=='t') then
        wxsurf => wtsurf                    ! Added by Zhuang Zhanpeng
      	xb => tb
      	allocate(tf(ks,im,jm))
      	xf => tf
      elseif (var_name=='s') then
        wxsurf => wssurf                    ! Added by Zhuang Zhanpeng
      	xb => sb
      	allocate(sf(ks,im,jm))
      	xf => sf
      endif
      
      do j=1,jm
        do i=1,im
          do k=1,ksm1
            rhs(k,i,j) = - adv_t(k,i,j) + diff_t(k,i,j)
            Kx_temp(k,i,j) = (Kh(k+1,i,j) + umol) * 2.0
          enddo
        enddo
      enddo
      
    endif
    
    xf = ZERO
  
    !---------------------------------------------------------------------------------------------
    ! Solve the tridiagonal matrix with the Thomas algorithm
    !---------------------------------------------------------------------------------------------
    do j=2,jmm1
      do i=2,imm1
        if (mask(i,j)/=1) cycle
        
        !-----------------------------------------------------------------------------------------
        ! Compute the tridiagonal coefficients
        !-----------------------------------------------------------------------------------------
        do k=1,ks-2
          a_temp = dti / hb_inx(i,j) / hm_inx(i,j) * Kx_temp(k,i,j) / (dz(k) + dz(k+1))
          a_u(k) = -a_temp / dz(k+1)
          c_u(k) = -a_temp / dz(k)
        enddo
        do k=2,ks-2  
          b_u(k) = Rh(i,j) - a_u(k-1) - c_u(k)     
        enddo
          b_u(1) = Rh(i,j) - c_u(1) + bot_temp(i,j)
          b_u(ksm1) = Rh(i,j) - a_u(ks-2)
        do k=1,ksm1
          Fu(k) = rhs(k,i,j) / hb_inx(i,j) * dti + xb(k,i,j)
        enddo
        
          Fu(ksm1) = Fu(ksm1) - wxsurf(i,j) * dti / dz(ksm1) / hb_inx(i,j)

        !-----------------------------------------------------------------------------------------
        ! Thomas algorithm solving tridiagonal
        !-----------------------------------------------------------------------------------------
        aa0(1) = a_u(1)
        do m=ks-2,1,-1
          aa0(m+1) = a_u(m)
        enddo
        cc0(1) = c_u(1) / b_u(1)
        Fu(1) = Fu(1) / b_u(1)
        do m=2,ks-2
          b_u(m) = b_u(m) - aa0(m) * cc0(m-1)
          cc0(m) = c_u(m) / b_u(m)
          Fu(m) = (Fu(m) - aa0(m) * Fu(m-1)) / b_u(m)
        enddo
        Fu(ksm1)=(Fu(ksm1) - aa0(ksm1) * Fu(ks-2)) / (b_u(ksm1) - aa0(ksm1) * cc0(ks-2))
        do m=ks-2,1,-1
          Fu(m) = Fu(m) - cc0(m) * Fu(m+1)
        enddo
        do k=1,ksm1
          xf(k,i,j) = Fu(k)
        enddo  
        
      enddo
    enddo
  
    nullify(wxsurf, xb, xf, mask, hb_inx, hm_inx)

    !---------------------------------------------------------------------------------------------
    ! Deallocate the momery spaces
    !---------------------------------------------------------------------------------------------
    deallocate(Rh)
    deallocate(a_u, b_u, c_u, Fu)
    deallocate(aa0, cc0)
    deallocate(rhs, Kx_temp)
    deallocate(bot_temp)
    
    return
  end subroutine solve_tridiag_sub
!=================================================================================================

!-------------------------------------------------------------------------------------------------
! 13)
! This adjust_uv_sub subroutine adjust uf vf to adapt to external mode flux.
!
! -- Originally developed by Lei HAN (c) August, 2012
! -- Code standardized by Zhanpeng ZHUANG in September, 2013
!-------------------------------------------------------------------------------------------------
  subroutine adjust_uv_sub
    
    implicit none

    real(kind=precision), allocatable :: tps_u(:,:), tps_v(:,:)

    !---------------------------------------------------------------------------------------------
    ! allocate the momery spaces
    !---------------------------------------------------------------------------------------------
    allocate(tps_u(im,jm), tps_v(im,jm))
      tps_u = ZERO
      tps_v = ZERO
  
    utf = utf / isplit
    vtf = vtf / isplit

    tps_u = 0.0
    do j=1,jm
      do i=1,im
        do k=1,ksm1
          tps_u(i,j) = tps_u(i,j) + uf(k,i,j) * dz(k)
        enddo
      enddo
    enddo
    do j=1,jm
      do i=1,im
        do k=1,ksm1
          uf(k,i,j) = uf(k,i,j) - tps_u(i,j) + utf(i,j)/hf_inu(i,j)
        enddo
      enddo
    enddo
  
    tps_v = 0.0
    do j=1,jm
      do i=1,im
        do k=1,ksm1
          tps_v(i,j) = tps_v(i,j) + vf(k,i,j) * dz(k)
        enddo
      enddo
    enddo
    do j=1,jm
      do i=1,im
        do k=1,ksm1
          vf(k,i,j) = vf(k,i,j) - tps_v(i,j) + vtf(i,j)/hf_inv(i,j)
        enddo
      enddo
    enddo

    do j=1,jm  
      do i=1,im
        do k=1,ks
          uf(k,i,j) = uf(k,i,j) * dum01(i,j)
          vf(k,i,j) = vf(k,i,j) * dvm01(i,j)
        enddo
      enddo
    enddo

    !---------------------------------------------------------------------------------------------
    ! Deallocate the momery spaces
    !---------------------------------------------------------------------------------------------
    deallocate(tps_u, tps_v)
    
    return
  end subroutine adjust_uv_sub
!=================================================================================================

!-------------------------------------------------------------------------------------------------
! 14)
! This upd_w_sub subroutine update the vertical velocity.
!
! -- Originally developed by Lei HAN (c) August, 2012
! -- Code standardized by Zhanpeng ZHUANG in September, 2013
!-------------------------------------------------------------------------------------------------
  subroutine upd_w_sub
    
    implicit none 

    real(kind=precision) :: tmp, tmp_zet
    real(kind=precision), allocatable :: w_up(:), w_dn(:), div_uv(:)
    real(kind=precision), allocatable :: flowx(:,:,:), flowy(:,:,:)

    !---------------------------------------------------------------------------------------------
    ! allocate the momery spaces
    !---------------------------------------------------------------------------------------------
    allocate(w_up(ks), w_dn(ks), div_uv(ks))
    allocate(flowx(ks,im,jm), flowy(ks,im,jm))
    allocate(wf(ks,im,jm))
      w_up = ZERO
      w_dn = ZERO
      div_uv = ZERO
      flowx = ZERO
      flowy = ZERO
      wf = ZERO

    do j=1,jm
      do i=1,im
        do k=1,ksm1
          flowx(k,i,j) = uf(k,i,j) * dy_u(i,j) * hf_inu(i,j)
          flowy(k,i,j) = vf(k,i,j) * dx_v(i,j) * hf_inv(i,j)
        enddo
      enddo
    enddo
  
    do j=2,jmm1
      do i=2,imm1
        if (fsm(i,j)==0) cycle     
        tmp_zet = (zetf_c(i,j) - zetb_c(i,j)) / dti
        do k=1,ksm1
          div_uv(k) = flowx(k,i,j) - flowx(k,i+1,j)    & 
                    + flowy(k,i,j) - flowy(k,i,j+1)
        enddo
        
        div_uv = div_uv / art(i,j)
        do k=ksm1,1,-1                      ! downward
          tmp = tmp_zet - div_uv(k)
          w_dn(k) = w_dn(k+1) + tmp * dz(k)
        enddo
        
        do k=1,ks
          wf(k,i,j) = w_dn(k)
        enddo
      enddo
    enddo

    if (regional_run) then
      call obcond_sub('GRD:wf')
    else
      call obcond_sub('CYC:wf')
    endif
  
    !---------------------------------------------------------------------------------------------
    ! Deallocate the momery spaces
    !---------------------------------------------------------------------------------------------
    deallocate(w_up, w_dn, div_uv)
    deallocate(flowx, flowy)
    
    return
  end subroutine upd_w_sub
!=================================================================================================

!-------------------------------------------------------------------------------------------------
! 14.5)
! This upd_ts_sub subroutine update the vertical velocity.
!
! -- Originally developed by Lei HAN (c) August, 2012
!-------------------------------------------------------------------------------------------------
  subroutine upd_ts_sub
    
    implicit none 
    
    real(kind=precision) :: sst_ice, dsst_obs, dt_prof
        
    if (mode==4) then
      allocate(tf(ks,im,jm),sf(ks,im,jm))
      deallocate(tb,sb)
      tf = tclim
      sf = sclim
      return
    endif

    allocate(adv_t(ks,im,jm),diff_t(ks,im,jm))

    !---------------------------------------------------------------------------------------------
    ! Solve the tracers' prognostic equations : Temperature (T)
    !---------------------------------------------------------------------------------------------
    adv_t  = ZERO
    diff_t = ZERO
    call adv_ts_sub('t')              ! Diffusion terms
    call diff_ts_sub('t')             ! Advection terms
    !adv_t  = ZERO
    !diff_t = ZERO
    if(sst_bound) then
      tb(ksm1,:,:) = sst_dq
      wtsurf = ZERO
    endif

    !!!!!!!!!!!!!!!!!!!!!
    ! add bsmt from MASNUM wave model
    do i=2,imm1
      do j=2,jmm1
        if(fsm(i,j) /= 1) cycle
        do k=1,ks
          diff_t(k,i,j) = diff_t(k,i,j)                                            &
                        - (tb(k,i,j) - tb(k,i-1,j)) / dx_t(i,j) * bsm1(k,i,j)      &
                        - (tb(k,i,j) - tb(k,i,j-1)) / dy_t(i,j) * bsm2(k,i,j)
          bsmt(k,i,j) = - (tb(k,i,j) - tb(k,i-1,j)) / dx_t(i,j) * bsm1(k,i,j)      &
                        - (tb(k,i,j) - tb(k,i,j-1)) / dy_t(i,j) * bsm2(k,i,j)
        enddo
      enddo
    enddo
    !!!!!!!!!!!!!!!!!!!!!
    
    call solve_tridiag_sub('t')    

    !---------------------------------------------------------------------------------------------
    ! Solve the tracers' prognostic equations : Salinity (S)
    !---------------------------------------------------------------------------------------------
    if (mode/=5) then
      adv_t  = ZERO
      diff_t = ZERO
      call adv_ts_sub('s')                     ! Advection terms
      call diff_ts_sub('s')                    ! Diffusion terms
    endif         

    !!!!!!!!!!!!!!!!!!!!!
    ! add bsms from MASNUM wave model
    do i=2,imm1
      do j=2,jmm1
        if(fsm(i,j) /= 1) cycle
        do k=1,ks
          diff_t(k,i,j) = diff_t(k,i,j)                                            &
                        - (sb(k,i,j) - sb(k,i-1,j)) / dx_t(i,j) * bsm1(k,i,j)      &
                        - (sb(k,i,j) - sb(k,i,j-1)) / dy_t(i,j) * bsm2(k,i,j)
          bsms(k,i,j) = - (sb(k,i,j) - sb(k,i-1,j)) / dx_t(i,j) * bsm1(k,i,j)      &
                        - (sb(k,i,j) - sb(k,i,j-1)) / dy_t(i,j) * bsm2(k,i,j)
        enddo
      enddo
    enddo
    !!!!!!!!!!!!!!!!!!!!!
     
    call solve_tridiag_sub('s')


    deallocate(adv_t,diff_t)


    !---------------------------------------------------------------------------------------------
    ! Update TS's open boundary value
    !---------------------------------------------------------------------------------------------
    if (regional_run) then
      call obcond_sub('ADV:TS')
!      call obcond_sub('SET:TS')
    else
      call obcond_sub('CYC:TS')
    endif
    
    deallocate(tb,sb)

    if (mode==5) sf = 35.
    
    !---------------------------------------------------------------------------------------------
    ! Impose the freezing temperature limit on SST
    !---------------------------------------------------------------------------------------------
    
    do j=1,jm
    	do i=1,im
    		if (fsm(i,j)==1) then
    			sst_ice = -0.0137 - ( 0.05199 + 7.225e-5*sf(ksm1,i,j) )*sf(ksm1,i,j)
    			do k=1,ks
       			tf(k,i,j) = max(sst_ice,tf(k,i,j))
       		enddo
       	endif
      enddo
    enddo

    !---------------------------------------------------------------------------------------------
    ! SST's OI (Optimum Interpolation) data assimilation from TMI-AMSRE's daily satellite data
    !---------------------------------------------------------------------------------------------
    if ( IF_sst_OI_assimilate ) then
    	
    	if ( iint==1 .or. ymdhms(3)/=ymdhms_prev(3) ) then
    		call get_2D_field('tmi_amsre_SST','analysed_sst')
      endif
    	
    	do j=1,jm
    	  do i=1,im
    		  
    		  if (fsm(i,j)==1 .and. sst_OI(i,j)>-99) then
    		  	dsst_obs = sst_coeff * ( sst_OI(i,j)-tf(ksm1,i,j) )
    		  	do k=ks-2,1,-1
    		  		dt_prof = tf(ksm1,i,j) - tf(k,i,j)
    		  		if (dt_prof>=0.0 .and. dt_prof<1.0) then
    		  			tf(k,i,j) = tf(k,i,j) + dsst_obs * ( 1.0-dt_prof )
    		  		endif
    		  	enddo
    		  	tf(ksm1,i,j) = tf(ksm1,i,j) + dsst_obs
    		  endif
    		  
    		enddo
    	enddo
    	
    	print*,'SST OI process finished. SST_OBS range:',   &
    	minval(sst_OI,mask=sst_OI>-99), maxval(sst_OI)
    	
    endif
    		  		



    
    return
  end subroutine upd_ts_sub
!=================================================================================================

!-------------------------------------------------------------------------------------------------
! 15)
! This obcond_sub subroutine updates the open boundary.
!
! -- Originally developed by Lei HAN (c) August, 2012
! -- Code standardized by Zhanpeng ZHUANG in September, 2013
!-------------------------------------------------------------------------------------------------
  subroutine obcond_sub(xobc)
    
    implicit none

    character(len=*),intent(in) :: xobc
  
    real(kind=precision) :: Elevation
    real(kind=precision) :: u1, v1
    real(kind=precision) :: e1
    real(kind=precision) :: wm                           ! ADV:TS as in POM
  
    !---------------------------------------------------------------------------------------------
    ! POM-uf
    !---------------------------------------------------------------------------------------------
    integer :: ibm1, ibp1, jbm1, jbp1
    real(kind=precision) :: ga, Hmax
  
    select case (xobc) 
    
  !===============================================================================================
  ! Global open boundary conditions
  !===============================================================================================
        
      !-------------------------------------------------------------------------------------------
      ! Periodic condition of external mode surface elevation for global run
      !-------------------------------------------------------------------------------------------
      case ('CYC:et')
        etf_c( 1,:) = etf_c(imm1,:)
        etf_c(im,:) = etf_c(   2,:)
        etf_c = etf_c * fsm                 ! mask
      return
      
      !-------------------------------------------------------------------------------------------
      ! Periodic condition of external mode ua,va for global run
      !-------------------------------------------------------------------------------------------
      case ('CYC:ua')
        ua( 1,:) = ua(imm1,:)
        ua(im,:) = ua(   2,:)
        va( 1,:) = va(imm1,:)
        va(im,:) = va(   2,:)
        ua = ua * dum01                     ! mask
        va = va * dvm01
      return
      
      !-------------------------------------------------------------------------------------------
      ! Periodic condition of external mode uah,vah for global run
      !-------------------------------------------------------------------------------------------
      case ('CYC:uaH')
        uaf( 1,:) = uaf(imm1,:)
        uaf(im,:) = uaf(   2,:)
        vaf( 1,:) = vaf(imm1,:)
        vaf(im,:) = vaf(   2,:)
        uaf = uaf * dum01                   ! mask
        vaf = vaf * dvm01
      return
      
      !-------------------------------------------------------------------------------------------
      ! Periodic condition of internal mode uf for global run
      !-------------------------------------------------------------------------------------------
      case ('CYC:uf')
        uf(:,1 ,:) = uf(:,imm1,:)
        uf(:,im,:) = uf(:,   2,:)
        do j=1,jm
          do i=1,im
            uf(:,i,j) = uf(:,i,j) * dum01(i,j)
          enddo
        enddo
      return
      
      !-------------------------------------------------------------------------------------------
      ! Periodic condition of internal mode vf for global run
      !-------------------------------------------------------------------------------------------
      case ('CYC:vf')
        vf(:, 1,:) = vf(:,imm1,:)
        vf(:,im,:) = vf(:,   2,:)
        do j=1,jm
          do i=1,im
            vf(:,i,j) = vf(:,i,j) * dvm01(i,j)
          enddo
        enddo
      return
      
      !-------------------------------------------------------------------------------------------
      ! Periodic condition of internal mode T,S for global run
      !-------------------------------------------------------------------------------------------
      case ('CYC:TS')
        tf(:,1 ,:) = tf(:,imm1,:)
        tf(:,im,:) = tf(:,   2,:)
        sf(:,1 ,:) = sf(:,imm1,:)
        sf(:,im,:) = sf(:,   2,:)
        do j=1,jm
          do i=1,im
            tf(:,i,j) = tf(:,i,j) * fsm(i,j)
            sf(:,i,j) = sf(:,i,j) * fsm(i,j)
          enddo
        enddo
      return
      
      !-------------------------------------------------------------------------------------------
      ! Periodic condition of vertical velocity for global run
      !-------------------------------------------------------------------------------------------
      case ('CYC:wf')
        wf(:,1 ,:) = wf(:,imm1,:)
        wf(:,im,:) = wf(:,   2,:)
        do j=1,jm
          do i=1,im
            wf(:,i,j) = wf(:,i,j) * fsm(i,j)
          enddo
        enddo
      return

      !-------------------------------------------------------------------------------------------
      ! Periodic condition of turbulence variables q2,q2l for global run
      !-------------------------------------------------------------------------------------------
      case ('CYC:q2l')
        q2f(:, 1,:) = q2f(:,imm1,:)
        q2f(:,im,:) = q2f(:,   2,:)    
        q2lf(:, 1,:) = q2lf(:,imm1,:)
        q2lf(:,im,:) = q2lf(:,   2,:)
        do j=1,jm
          do i=1,im
            q2f( :,i,j) = q2f( :,i,j) * fsm(i,j) + SMALL
            q2lf(:,i,j) = q2lf(:,i,j) * fsm(i,j) + SMALL
          enddo
        enddo
      return

  !===============================================================================================
  ! Regional OB conditions
  !===============================================================================================
  
      !-------------------------------------------------------------------------------------------
      ! Prescribed elevation from larger-domain simulations
      !-------------------------------------------------------------------------------------------
      case ('SET:et-circ') 
        do n=1,numboundary
          ib = ij_obc(1,n)
          jb = ij_obc(2,n)
          etf_c(ib,jb) = et_obc(n)
        enddo
        do j=1,jm
          do i=1,im
            etf_c(i,j) = etf_c(i,j) * fsm(i,j)
          enddo
        enddo
      return
    
      !-------------------------------------------------------------------------------------------
      ! Prescribed elevation for tidal forcing boundary
      !-------------------------------------------------------------------------------------------
      case ('SET:et-tide')
        do n=1,numboundary
          Elevation = 0.
          do j=1,numtide
            Elevation = Elevation     & 
                      + Amplitude(j,n) * cosd( Omega_tide(j)*days_now-Phase(j,n) )
          enddo
          ib = ij_obc(1,n)
          jb = ij_obc(2,n)
          et_obc(n) = et_obc(n) + Elevation * ramp_tide  
                                            ! circulation + tide
        enddo
      return

      !-------------------------------------------------------------------------------------------
      ! zero-Gradient OBC for vertical velocity   
      !-------------------------------------------------------------------------------------------
      case ('GRD:wf') 
        do n=1,numboundary
          ib = ij_obc(1,n)
          jb = ij_obc(2,n)
          select case (bctype(n))
            case (3)                        ! Eastern
              wf(:,ib,jb) = wf(:,ib-1,jb)    
            case (0)                        ! Northern
              wf(:,ib,jb) = wf(:,ib,jb-1)    
            case (9)                        ! Western
              wf(:,ib,jb) = wf(:,ib+1,jb)    
            case (6)                        ! Southern
              wf(:,ib,jb) = wf(:,ib,jb+1)
          end select
        enddo  
    
        do j=1,jm
          do i=1,im
            wf(:,i,j) = wf(:,i,j) * fsm(i,j)
          enddo
        enddo
      return
  
      !-------------------------------------------------------------------------------------------
      ! Radiation OBC for external mode uah,vah (Flather, 1976)
      !-------------------------------------------------------------------------------------------
      case ('FLA:uaH')
        rfe= 1.
        rfw= 1.
        rfs= 1.
        rfn= 1.
      
        do n=1,numboundary
          ib = ij_obc(1,n)
          jb = ij_obc(2,n)
          select case (bctype(n))
            case (3)                        ! Eastern Open boundary
              uaf(ib,jb) = ua_obc(n)                                   &
                         + rfe * sqrt( GRAV / H_c(ib-1,jb) )           &
                               * ( etf_c(ib-1,jb) - et_obc(n) )
              uaf(ib,jb) = uaf(ib,jb) * ramp_obc * H_u(ib,jb)
              vaf(ib,jb) = 0.
            case (0)                        ! Northern Open boundary
              vaf(ib,jb) = va_obc(n)                                   &
                         + rfn * sqrt( GRAV / H_c(ib,jb-1) )           &
                               * ( etf_c(ib,jb-1) - et_obc(n) )
              vaf(ib,jb) = vaf(ib,jb) * ramp_obc * H_v(ib,jb)
              uaf(ib,jb) = 0.
            case (9)                        ! Western Open boundary
              uaf(ib+1,jb) = ua_obc(n)                                 &
                           - rfw * sqrt( GRAV / H_c(ib+1,jb) )         &
                                 * (etf_c(ib+1,jb) - et_obc(n))
              uaf(ib+1,jb) = uaf(ib+1,jb) * ramp_obc * H_u(ib+1,jb)
              uaf(ib,jb) = uaf(ib+1,jb)
              vaf(ib,jb) = 0.
            case (6)                        ! Southern Open boundary
              vaf(ib,jb+1) = va_obc(n)                                 &
                           - rfs * sqrt( GRAV / H_c(ib,jb+1) )         &
                                 * ( etf_c(ib,jb+1) - et_obc(n) )
              vaf(ib,jb+1) = vaf(ib,jb+1) * ramp_obc * H_v(ib,jb+1)
              vaf(ib,jb) = vaf(ib,jb+1)
              uaf(ib,jb) = 0.  
          end select
        enddo
        
        uaf = uaf * dum01
        vaf = vaf * dvm01
        
        if (ChangJiang_diluting) then
          uaf(chj_i,chj_j) = time_chj_q / dy_t(chj_i,chj_j)
        endif
      return

      !-------------------------------------------------------------------------------------------
      ! Prescribed T&S from larger-domain simulations
      !-------------------------------------------------------------------------------------------
      case ('SET:TS')
        do n=1,numboundary
          ib = ij_obc(1,n)
          jb = ij_obc(2,n)
          do k=1,ks
            tf(k,ib,jb) = tb_obc(k,n)
            sf(k,ib,jb) = sb_obc(k,n)
          enddo
        enddo
        do j=1,jm
          do i=1,im
            tf(:,i,j) = tf(:,i,j) * fsm(i,j)
            sf(:,i,j) = sf(:,i,j) * fsm(i,j)
          enddo
        enddo
      return
        
      !-------------------------------------------------------------------------------------------
      ! Upstream T&S OBC (Mellor, 1996)
      !-------------------------------------------------------------------------------------------
      case ('ADV:TS')
        do n=1,numboundary
          ib = ij_obc(1,n)
          jb = ij_obc(2,n)
          ibm1 = ib
          jbm1 = jb
          
          
          
          select case (bctype(n))
            case (3)                        ! Eastern Open boundary
              do k=1,ksm1
                u1 = uf(k,ib,jb)*dti/dx_u(ib,jb)
                if (u1>0) then
                  tf(k,ib,jb) = tb(k,ib,jb) - u1 * (tb(k,ib,jb) - tb(k,ib-1,jb))
                  sf(k,ib,jb) = sb(k,ib,jb) - u1 * (sb(k,ib,jb) - sb(k,ib-1,jb))
                else
                  tf(k,ib,jb) = tb(k,ib,jb)                              &
                              - u1 * (tb_obc(k,n) - tb(k,ib,jb))
                  sf(k,ib,jb) = sb(k,ib,jb)                              &
                              - u1 * (sb_obc(k,n) - sb(k,ib,jb))
                endif
              enddo
            case (0)                        ! Northern Open boundary
              do k=1,ksm1
                v1 = vf(k,ib,jb) * dti / dy_v(ib,jb)
                if (v1>0) then
                  tf(k,ib,jb) = tb(k,ib,jb)                              &
                              - v1 * (tb(k,ib,jb) - tb(k,ib,jb-1))
                  sf(k,ib,jb) = sb(k,ib,jb)                              &
                              - v1 * (sb(k,ib,jb) - sb(k,ib,jb-1))
                else
                  tf(k,ib,jb) = tb(k,ib,jb)                              &
                              - v1 * (tb_obc(k,n) - tb(k,ib,jb))
                  sf(k,ib,jb) = sb(k,ib,jb)                              &
                              - v1 * (sb_obc(k,n) - sb(k,ib,jb))
                endif
              enddo
            case (9)                        ! Western Open boundary
              do k=1,ksm1
                u1 = uf(k,ib+1,jb) * dti / dx_u(ib+1,jb)
                if (u1<0) then
                  tf(k,ib,jb) = tb(k,ib,jb)                              & 
                              - u1 * (tb(k,ib+1,jb) - tb(k,ib,jb))
                  sf(k,ib,jb) = sb(k,ib,jb)                              &
                              - u1 * (sb(k,ib+1,jb) - sb(k,ib,jb))
                else
                  tf(k,ib,jb) = tb(k,ib,jb)                              &
                              - u1 * (tb(k,ib,jb) - tb_obc(k,n))
                  sf(k,ib,jb) = sb(k,ib,jb)                              &
                              - u1 * (sb(k,ib,jb) - sb_obc(k,n))
                endif
              enddo
            case (6)                        ! Southern Open boundary
              do k=1,ksm1
                v1 = vf(k,ib,jb+1) * dti / dy_v(ib,jb+1)
                if (v1<0) then
                  tf(k,ib,jb) = tb(k,ib,jb)                              &
                              - v1 * (tb(k,ib,jb+1) - tb(k,ib,jb))
                  sf(k,ib,jb) = sb(k,ib,jb)                              &
                              - v1 * (sb(k,ib,jb+1) - sb(k,ib,jb))
                else
                  tf(k,ib,jb) = tb(k,ib,jb)                              &
                              - v1 * (tb(k,ib,jb) - tb_obc(k,n))
                  sf(k,ib,jb) = sb(k,ib,jb)                              &
                              - v1 * (sb(k,ib,jb) - sb_obc(k,n))
                endif      
                
                
              enddo
          end select

          if (bctype(n)==6) tf(:,ib,jb) = tb_obc(:,n) 
                                            ! SET South-T OBC
        enddo   

        do j=1,jm
          do i=1,im
            tf(:,i,j) = tf(:,i,j) * fsm(i,j)
            sf(:,i,j) = sf(:,i,j) * fsm(i,j)
          enddo
        enddo

        if (ChangJiang_diluting) sf(:,chj_i,chj_j) = 1.0
                                            ! Changjiang diluting water
      return
  
      !-------------------------------------------------------------------------------------------
      ! Radiation OBC of uf (POM)
      !-------------------------------------------------------------------------------------------
      case ('POM:uf')
        Hmax = maxval(H_c)
        do n=1,numboundary
          ib = ij_obc(1,n)
          jb = ij_obc(2,n)
          ga = sqrt(H_c(ib,jb) / Hmax)
          ibm1 = max(ib-1, 1)
          ibp1 = min(ib+1,im)
          jbm1 = max(jb-1, 1)
          jbp1 = min(jb+1,jm)

          select case (bctype(n))
            case (3)                        ! Eastern Open boundary
              do k=1,ksm1
                uf(k,ib,jb) =                                                                    &
                      ga  * ( 0.25*ub(k,ib-1,jbm1) + 0.5*ub(k,ib-1,jb) + 0.25*ub(k,ib-1,jbp1) )  &
                + (1.-ga) * ( 0.25*ub(k,ib  ,jbm1) + 0.5*ub(k,ib  ,jb) + 0.25*ub(k,ib  ,jbp1) )  
              enddo
              
              if (ib==307.and.jb==47) then
                write(16,'(<ksm1>e15.5)') (uf(k,ib,jb),k=1,ksm1) 
              endif


            case (0)                        ! Northern Open boundary
              do k=1,ksm1
                uf(k,ib,jb) = 0.
              enddo

            case (9)                        ! Western Open boundary
              do k=1,ksm1
                uf(k,ib+1,jb) =                                                                  &       
                      ga  * ( 0.25*ub(k,ib+2,jbm1) + 0.5*ub(k,ib+2,jb) + 0.25*ub(k,ib+2,jbp1) )  &
                + (1.-ga) * ( 0.25*ub(k,ib+1,jbm1) + 0.5*ub(k,ib+1,jb) + 0.25*ub(k,ib+1,jbp1) )
                uf(k,ib,jb) = uf(k,ib+1,jb)
              enddo

            case (6)                        ! Southern Open boundary
              do k=1,ksm1
                !---------------------------------------------------------------------------------
                ! Force south --- From Xia CX WestPac run
                !---------------------------------------------------------------------------------
                uf(k,ib,jb+1) = ub_obc(k,n)
                uf(k,ib,jb  ) = uf(k,ib,jb+1)
              enddo
          end select
        enddo
        
        do j=1,jm
          do i=1,im
            uf(:,i,j) = uf(:,i,j) * dum01(i,j)
          enddo
        enddo
        
      return
      
      case ('POM:vf')
        Hmax = maxval(H_c) 
        do n=1,numboundary
          ib = ij_obc(1,n)
          jb = ij_obc(2,n)
          ga = sqrt(H_c(ib,jb) / Hmax)
          ibm1 = max(ib-1, 1)
          ibp1 = min(ib+1,im)
          jbm1 = max(jb-1, 1)
          jbp1 = min(jb+1,jm)

          select case (bctype(n))
            case (3)                        ! Eastern Open boundary
              do k=1,ksm1
                vf(k,ib,jb) = 0.
              enddo

            case (0)                        ! Northern Open boundary
              do k=1,ksm1
                vf(k,ib,jb) =                                                                    &
                      ga  * ( 0.25*vb(k,ibm1,jb-1) + 0.5*vb(k,ib,jb-1) + 0.25*vb(k,ibp1,jb-1))   &
                + (1.-ga) * ( 0.25*vb(k,ibm1,  jb) + 0.5*vb(k,ib,  jb) + 0.25*vb(k,ibp1,  jb))  
              enddo

            case (9)                        ! Western Open boundary
              do k=1,ksm1
                vf(k,ib,jb) = 0.
              enddo

            case (6)                        ! Southern Open boundary
              do k=1,ksm1
                !---------------------------------------------------------------------------------
                ! Force south --- From Xia CX WestPac run
                !---------------------------------------------------------------------------------
                vf(k,ib,jb+1) = vb_obc(k,n)
                vf(k,ib,jb  ) = vf(k,ib,jb+1)
              enddo
          end select
        enddo
        
        do j=1,jm
          do i=1,im
            vf(:,i,j) = vf(:,i,j) * dvm01(i,j)
          enddo
        enddo
        
      return
      
    end select
  end subroutine obcond_sub
!=================================================================================================

  !-----------------------------------------------------------------------------------------------
  ! 16)
  ! This adv_ts_sub subroutine computes the advection term of the tracers' equation. 
  !
  ! -- Originally developed by Lei HAN (c) August, 2012
  ! -- Code standardized by Lei HAN in September, 2013
  !-----------------------------------------------------------------------------------------------  
  subroutine adv_ts_sub(var_name)
    
    implicit none 

    character(len=1), intent(in) :: var_name
    real(kind=precision), pointer :: rb(:,:,:)
    
    real(kind=precision), allocatable :: adv_x(:,:,:), adv_y(:,:,:), adv_zt(:,:,:)
    real(kind=precision) :: uface, vface, wface
    integer :: iup, jup, kup

    !---------------------------------------------------------------------------------------------
    ! allocate the momery spaces
    !---------------------------------------------------------------------------------------------
    allocate(adv_x(ks,im,jm), adv_y(ks,im,jm), adv_zt(ks,im,jm))
      adv_x  = ZERO
      adv_y  = ZERO
      adv_zt = ZERO

    if (var_name=='t') then
    	rb => tb
    elseif (var_name=='s') then
    	rb => sb
    endif
        
    do j=2,jmm1
      do i=2,imm1
        if (fsm(i,j)==0) cycle
        do k=2,ks
          wface = wf(k,i,j)
          adv_zt(k,i,j) = wface * (rb(k,i,j)+rb(k-1,i,j))/2.0
        enddo
        adv_zt(1,i,j) = wf(1,i,j) * rb(1,i,j)
      enddo
    enddo
    do j=2,jm
      do i=2,im
        if (dum(i,j)==0) cycle
        do k=1,ksm1
          uface = uf(k,i,j)
          adv_x(k,i,j) = uface * (rb(k,i,j)+rb(k,i-1,j))*0.5      &
                               * hf_inu(i,j) * dy_u(i,j)
        enddo
      enddo
    enddo
    do j=2,jm
      do i=2,im
        if (dvm(i,j)==0) cycle
        do k=1,ksm1
          vface = vf(k,i,j)
          adv_y(k,i,j) = vface * (rb(k,i,j)+rb(k,i,j-1))*0.5      &
                               * hf_inv(i,j) * dx_v(i,j)
        enddo
      enddo
    enddo

    do j=2,jmm1
      do i=2,imm1
        if (fsm(i,j)==0) cycle
        do k=1,ksm1
          adv_t(k,i,j) = -adv_x(k,i,j) + adv_x(k,i+1,j)     &
                         -adv_y(k,i,j) + adv_y(k,i,j+1)  
          adv_t(k,i,j) =  adv_t(k,i,j) / art(i,j)           &
                       + ( -adv_zt(k,i,j) + adv_zt(k+1,i,j) ) / dz(k)
        enddo
      enddo
    enddo
    
    nullify(rb)

    !---------------------------------------------------------------------------------------------
    ! Deallocate the momery spaces
    !---------------------------------------------------------------------------------------------
    deallocate(adv_x, adv_y, adv_zt)
    
    return
  end subroutine adv_ts_sub

!====================================================================

  !-----------------------------------------------------------------------------------------------
  ! 17)
  ! This diff_ts_sub subroutine computes the diffusion term of the tracers' equation. 
  !
  ! -- Originally developed by Lei HAN (c) August, 2012
  ! -- Code standardized by Lei HAN in September, 2013
  !-----------------------------------------------------------------------------------------------  
  subroutine diff_ts_sub(var_name)
    
    implicit none
    
    character(len=1), intent(in) :: var_name
    real(kind=precision), pointer :: rb(:,:,:), xclim(:,:,:)
    
    real(kind=precision), allocatable :: diff_x(:,:,:), diff_y(:,:,:)
    real(kind=precision), allocatable :: Ah(:,:,:)

    real(kind=precision) :: p1, p2
    
    !---------------------------------------------------------------------------------------------
    ! allocate the momery spaces
    !---------------------------------------------------------------------------------------------
    allocate(diff_x(ks,im,jm), diff_y(ks,im,jm))
    allocate(Ah(ks,im,jm))
      diff_x = ZERO
      diff_y = ZERO
      Ah = ZERO
    
    Ah = Am * TPRNI
    
    if (var_name=='t') then
    	rb => tb
    	xclim => tclim
    elseif (var_name=='s') then
    	rb => sb
    	xclim => sclim
    endif
    
    do i=1,im                               ! Added by Zhuang Zhanpeng
      do j=1,jm
        do k=1,ks
          rb(k,i,j) = rb(k,i,j) - xclim(k,i,j)
                                            ! Assuming no diffusion for the climatological 
                                            ! background field
        enddo
      enddo
    enddo

    do j=2,jm
      do i=2,im
        if (dum(i,j)==0) cycle
        do k=1,ksm1
          p1 = ( Ah(k,i,j)+Ah(k,i-1,j) + umol*2. ) * 0.5  &
             * ( rb(k,i,j)-rb(k,i-1,j) ) / dx_u(i,j)
          diff_x(k,i,j) = p1 * dy_u(i,j) * hf_inu(i,j) * dum01(i,j)
        enddo
      enddo
    enddo
    
    do j=2,jm
      do i=2,im
        if (dvm(i,j)==0) cycle
        do k=1,ksm1
          p2 = ( Ah(k,i,j)+Ah(k,i,j-1) + umol*2. ) * 0.5 &
             * ( rb(k,i,j)-rb(k,i,j-1) ) / dy_v(i,j)
          diff_y(k,i,j) = p2 * dx_v(i,j) * hf_inv(i,j) * dvm01(i,j)
        enddo
      enddo
    enddo
    
    do j=2,jmm1
      do i=2,imm1
        if (fsm(i,j)==0) cycle  
        do k=1,ksm1
          diff_t(k,i,j) = - diff_x(k,i,j) + diff_x(k,i+1,j)   &
                          - diff_y(k,i,j) + diff_y(k,i,j+1) 
          diff_t(k,i,j) = diff_t(k,i,j) / art(i,j)
        enddo
      enddo
    enddo

    do i=1,im                               ! Added by Zhuang Zhanpeng
      do j=1,jm
        do k=1,ks
          rb(k,i,j) = rb(k,i,j) + xclim(k,i,j)
        enddo
      enddo
    enddo
    
    nullify(rb,xclim)

    !---------------------------------------------------------------------------------------------
    ! Deallocate the momery spaces
    !---------------------------------------------------------------------------------------------
    deallocate(diff_x, diff_y)
    deallocate(Ah)
    
    return
  end subroutine diff_ts_sub

!=================================================================================================

  !-----------------------------------------------------------------------------------------------
  ! 18)
  ! This teos_sub subroutine computes the in-situ density and the static stability 
  ! (or Brunt-Vaisala frequency) with the uptodate Equation of State of the sea water 
  ! (see Jackett, McDougall, Feistel, and Wright, 2006). 
  !
  ! -- Originally developed by Lei HAN (c) August, 2012
  ! -- Code standardized by Lei HAN in September, 2013
  !-----------------------------------------------------------------------------------------------  
  subroutine teos_sub
    
    implicit none
  
    real(kind=precision) :: A0, B0
    real(kind=precision) :: A(11), B(12)
    real(kind=precision) :: pr(ks)
    real(kind=precision) :: t1, s1, p1, t2, t3
    real(kind=precision) :: ss, sx, p2, p3, pt2
    real(kind=precision) :: Pn, Pd
    real(kind=precision) :: dpdt, dpds, drds
    real(kind=precision) :: dtdz, dsdz
    real(kind=precision) :: dzm
    real(kind=precision), allocatable :: abet(:,:,:), beta(:,:,:)
  
    allocate(abet(ks,im,jm), beta(ks,im,jm))
  
    !-----------------------------------------------------------------------------------------------      
    ! Define Coefficients  
    !-----------------------------------------------------------------------------------------------  
    data A / 7.3471625860981584    , &
            -5.3211231792841769e-2 , &
             3.6492439109814549e-4 , &
             2.5880571023991390    , &
            -6.7168282786692355e-3 , &
             1.9203202055760151e-3 , &
             1.1798263740430364e-2 , &
             9.8920219266399117e-8 , &
             4.6996642771754730e-6 , &
            -2.5862187075154352e-8 , &
            -3.2921414007960662e-12  /
    data B / 7.2815210113327091e-3 , &
            -4.4787265461983921e-5 , &
             3.3851002965802430e-7 , &
             1.3651202389758572e-10, &
             1.7632126669040377e-3 , &
            -8.8066583251206474e-6 , &
            -1.8832689434804897e-10, &
             5.7463776745432097e-6 , &
             1.4716275472242334e-9 , &
             6.7103246285651894e-6 , &
            -2.4461698007024582e-17, &
            -9.1534417604289062e-18  /            
    A0 = 9.9984085444849347e2
    B0 = 1.0
    
    !-----------------------------------------------------------------------------------------------      
    ! Update in-situ density, thermal & salinity expansion coefficients
    !-----------------------------------------------------------------------------------------------  
    rclim = rclim + 1000.
    do j=1,jm
      do i=1,im
        if (fsm(i,j)==0) cycle
        do k=ksm1,1,-1
          if (k==ksm1) then
            pr(k) = rclim(k,i,j) * hf_inc(i,j) * dz(k) / 2.0
          else
            pr(k) = pr(k+1) + ( rclim(k+1,i,j) + rclim(k,i,j) )   &
                              * hf_inc(i,j) * dzz(k) / 2.0
          endif
          
          !-----------------------------------------------------------------------------------------------      
          ! Calculate TSP terms
          !-----------------------------------------------------------------------------------------------  
          t1   =  tf(k,i,j)
          s1   =  sf(k,i,j)
          p1   =  pr(k)  * GRAV * 1.0e-4
          t2   =  t1 * t1
          t3   =  t2 * t1
    
          ss  =  sqrt(s1)
          sx   =  ss * s1
          
          pt2  =  p1 * t2
          p2   =  p1 * p1
          p3   =  p1 * p2
       
          !-----------------------------------------------------------------------------------------------      
          ! Calculate the factors of EOS
          !-----------------------------------------------------------------------------------------------       
          Pn = A0                                &
             + t1 * ( A(1) + A(2)*t1 + A(3)*t2 ) &
             + s1 * ( A(4) + A(5)*t1 + A(6)*s1 ) &
             + p1 * ( A(7) + A(8)*t2 + A(9)*s1 ) &
             + p2 * ( A(10) + A(11)*t2 )
          Pd = B0                                            &
             + t1 * ( B(1) + B(2)*t1 + B(3)*t2 + B(4)*t3 )   &
             + s1 * ( B(5) + B(6)*t1 + B(7)*t3 )             &
             + sx * ( B(8) + B(9)*t2 )                       &
             + p1 * ( B(10)+ t1*( B(11)*pt2 + B(12)*p2 ) )
          rho(k,i,j) = Pn / Pd
    
          !-----------------------------------------------------------------------------------------------      
          ! Calculate the static stablity
          !-----------------------------------------------------------------------------------------------       
          dpdt = A(1) + ( A(5)-B(6) )*s1 - B(12)*p3                         &
               + 2. * t1 * ( A(2) + A(8)*p1 + A(11)*p2 - B(2) - B(9)*sx )   &
               + 3. * t2 * ( A(3) - B(3) - B(7)*s1 - B(11)*p2 )             &
               - 4. * B(4) * t3
          dpds = A(4) - B(5) + t1 * ( A(5) - B(6) - B(7)*t2 )               &
               + A(9)*p1 + 2.*A(6)*s1 - 1.5*ss*( B(8)+B(9)*t2 )
          
          drds = dpds / Pd / Pd
          abet(k,i,j) = - dpdt / dpds
          beta(k,i,j) =   drds / rho(k,i,j)
        enddo
      enddo
    enddo
    
    do j=1,jm
      do i=1,im
        rho(ks,i,j) = rho(ksm1,i,j)
      enddo
    enddo
    
    !-----------------------------------------------------------------------------------------------      
    ! Compute the static stability N2
    !-----------------------------------------------------------------------------------------------          
    do j=1,jm
      do i=1,im
        if (fsm(i,j)==0) cycle
        do k=1,ksm1
          dzm  = hf_inc(i,j) * dzz(k)
          if (k==ksm1) then
            dtdz = ( tf(k,i,j)-tf(k-1,i,j) ) / dzm
            dsdz = ( sf(k,i,j)-sf(k-1,i,j) ) / dzm
          else     
            dtdz = ( tf(k+1,i,j)-tf(k,i,j) ) / dzm
            dsdz = ( sf(k+1,i,j)-sf(k,i,j) ) / dzm
          endif
          N2(k,i,j) = abet(k,i,j) * dtdz - dsdz
          N2(k,i,j) = N2(k,i,j) * GRAV * beta(k,i,j)
        enddo
      enddo
    enddo
    
    rclim = rclim - 1000.                   ! Transform to be the density excess

    return 
  end subroutine teos_sub

!====================================================================

  !-----------------------------------------------------------------------------------------------
  ! 19)
  ! This advq_POM98 subroutine computes the adv&diff terms of the turbulence sub-model equations. 
  ! This subroutine is entirely transplanted from the POM98 only that the time stepping algorithm
  ! is reformed to be two-level scheme.
  !
  ! -- Originally developed by POM98
  ! -- Code reformed and standardized by Lei HAN in September, 2013
  !-----------------------------------------------------------------------------------------------  
  subroutine advq_POM98(var_name)
    
    implicit none
    
    character(len=*), intent(in) :: var_name
    real(kind=precision), pointer :: qb(:,:,:), qf(:,:,:)
    real(kind=precision), allocatable :: xflux(:,:,:),yflux(:,:,:)

    !---------------------------------------------------------------------------------------------
    ! allocate the momery spaces
    !---------------------------------------------------------------------------------------------
    allocate(xflux(kb,im,jm),yflux(kb,im,jm))
      xflux = ZERO
      yflux = ZERO
      
    if (var_name=='q2') then
    	qb => q2b
    	allocate(q2f(ks,im,jm))
    	qf => q2f
    elseif (var_name=='q2l') then
    	qb => q2lb
    	allocate(q2lf(ks,im,jm))
    	qf => q2lf
    endif    	
    
    qf = ZERO
        
    !-----------------------------------------------------------------------------------------------
    ! Do horizontal advection
    !-----------------------------------------------------------------------------------------------
    do j=2,jm
      do i=2,im
        if (dum(i,j)==0) cycle
        do k=2,kbm1
          xflux(k,i,j) = .25 * ( qb(k,i,j)+qb(k,i-1,j) )          &
                             * hm_inu(i,j)                        &
                             * ( uf(ks-k,i,j)+uf(ks-k+1,i,j) )
        end do
      end do
    end do
    
    do j=2,jm
      do i=2,im
        if (dvm(i,j)==0) cycle
        do k=2,kbm1
          yflux(k,i,j) = .25 * ( qb(k,i,j)+qb(k,i,j-1) )          &
                             * hm_inv(i,j)                        &
                             * ( vf(ks-k,i,j)+vf(ks-k+1,i,j) )
        end do
      end do
    end do

    !-----------------------------------------------------------------------------------------------
    ! Do horizontal diffusion
    !-----------------------------------------------------------------------------------------------
    do j=2,jm
      do i=2,im
        if (dum(i,j)==0) cycle
        do k=2,kbm1
          xflux(k,i,j) = xflux(k,i,j)                             &
                -.25 * ( Am(ks-k  ,i,j) + Am(ks-k  ,i-1,j)        &
                       + Am(ks-k+1,i,j) + Am(ks-k+1,i-1,j) )      &
                     * hm_inu(i,j)                                &
                     * ( qb(k,i,j)-qb(k,i-1,j))                   &
                     / dx_u(i,j)
          xflux(k,i,j) = dy_u(i,j)*xflux(k,i,j)
        end do
      end do
    end do

    do j=2,jm
      do i=2,im
        if (dvm(i,j)==0) cycle
        do k=2,kbm1
          yflux(k,i,j) = yflux(k,i,j)                             &
                 -.25 * ( Am(ks-k  ,i,j) + Am(ks-k  ,i,j-1)       &
                        + Am(ks-k+1,i,j) + Am(ks-k+1,i,j-1) )     &
                      * hm_inv(i,j)                               &
                      * ( qb(k,i,j)-qb(k,i,j-1) )                 &
                      / dy_v(i,j)
          yflux(k,i,j) = dx_v(i,j)*yflux(k,i,j)
        end do
      end do
    end do

    !-----------------------------------------------------------------------------------------------
    ! Do vertical advection, add flux terms, then step forward in time:
    !-----------------------------------------------------------------------------------------------
    do j=2,jmm1
      do i=2,imm1
        if (fsm(i,j)==0) cycle
        do k=2,kbm1
          qf(k,i,j) = ( wf(ks-k+2,i,j) * qb(k-1,i,j)         &
                       -wf(ks-k  ,i,j) * qb(k+1,i,j) )       &
                    * art(i,j) / (dz(ks-k)+dz(ks-k+1))       &
                    + xflux(k,i+1,j) - xflux(k,i,j)          &
                    + yflux(k,i,j+1) - yflux(k,i,j)

          qf(k,i,j) = ( hb_inc(i,j)*art(i,j)*qb(k,i,j) - dti*qf(k,i,j)) &
                    / ( hf_inc(i,j)*art(i,j) )
        end do
      end do
    end do
        
    nullify(qb,qf)

    !---------------------------------------------------------------------------------------------
    ! Deallocate the momery spaces
    !---------------------------------------------------------------------------------------------
    deallocate(xflux,yflux)
    
    if (var_name=='q2l') deallocate(Am)
    
    return
  end subroutine advq_POM98

!=================================================================================================

  !-----------------------------------------------------------------------------------------------
  ! 21)
  ! This compute_dens_POM98 subroutine computes the in-situ density (adopted from POM98 code).
  !
  ! -- Originally developed by POM98
  ! -- Code standardized by Lei HAN in September, 2013
  !-----------------------------------------------------------------------------------------------  
  subroutine compute_dens_POM98
    
    implicit none
    
    real(kind=precision), pointer :: dt(:,:)
    real(kind=precision) :: tr,tbias,sr,sbias
    real(kind=precision) :: tr2,tr3,tr4
    real(kind=precision) :: p,rhor,cr
    
    allocate(rho(ks,im,jm))
      rho = rho_ref
            
    !---------------------------------------------------------------------------------------------  
    !    If using 32 bit precision, it is recommended that
    !    TR, SR, P, RHOR , CR be made double precision.
    !    and the E's in the constants should be changed 
    !    to D's.
    !
    !    THIS SUBROUTINE COMPUTES DENSITY- 1.000
    !    T = POTENTIAL TEMPERATURE
    ! ( See: Mellor, 1991, J. Atmos. Oceanic Tech., 609-611)
    !---------------------------------------------------------------------------------------------   
    
    dt => hf_inc
    tbias=0.
    sbias=0.
        
    DO J=1,JM
      DO I=1,IM
        if (fsm(i,j)==0) cycle
        DO K=1,KBM1
          TR=tf(ks-k,I,J)+TBIAS
          SR=sf(ks-k,I,J)+SBIAS
          TR2=TR*TR
          TR3=TR2*TR
          TR4=TR3*TR
          ! Approximate pressure in units of bars
          P=-GRAV*1.025*ZZ(ks-K)*DT(I,J)*0.01
          RHOR = -0.157406  + 6.793952E-2*TR              &
                 - 9.095290E-3*TR2 + 1.001685E-4*TR3      &
                 - 1.120083E-6*TR4 + 6.536332E-9*TR4*TR
          RHOR = RHOR + (0.824493 - 4.0899E-3*TR          &
                + 7.6438E-5*TR2 - 8.2467E-7*TR3           & 
                + 5.3875E-9*TR4) * SR                     &
                + (-5.72466E-3 + 1.0227E-4*TR             &
                - 1.6546E-6*TR2) * ABS(SR)**1.5           &
                + 4.8314E-4 * SR*SR
          CR=1449.1+.0821*P+4.55*TR-.045*TR2              &
                   +1.34*(SR-35.)
          RHOR=RHOR + 1.E5*P/(CR*CR)*(1.-2.*P/(CR*CR))
          rho(Ks-k,I,J)=RHOR*1.E-3*FSM(I,J)
        end do
      enddo
    enddo
    
    do j=1,jm
      do i=1,im
        rho(ks,i,j) = rho(ksm1,i,j)
      enddo
    enddo
    
    !---------------------------------------------------------------------------------------------  
    ! Transform (rho-1000)/1000 to be the (rho-1000)(density excess) [Added by Lei Han] 
    !---------------------------------------------------------------------------------------------  
    do j=1,jm
      do i=1,im
        do k=1,ks
          rho(k,i,j) = rho(k,i,j) * 1000. * fsm(i,j)
        enddo
      enddo
    enddo
      
    nullify(dt)
        
    return
  END subroutine compute_dens_POM98

!=================================================================================================

  !-----------------------------------------------------------------------------------------------
  ! 22)
  ! This profq_POM98 subroutine solves the turbulence sub-model and yields the vertical mixing
  ! coefficients to the circulation model. This subroutine is entirely transplanted from the POM98
  ! only that the time stepping algorithm is reformed to be two-level scheme.
  !
  ! -- Originally developed by POM
  ! -- Code reformed and standardized by Lei HAN in September, 2013
  !-----------------------------------------------------------------------------------------------  

  subroutine profq_POM98
    
    implicit none

    real(kind=precision), allocatable :: GH(:,:,:),SM(:,:,:),SH(:,:,:)
    real(kind=precision), allocatable :: PROD(:,:,:),KN(:,:,:),BOYGR(:,:,:)
    real(kind=precision), allocatable :: DH(:,:),CC(:,:,:),STF(:,:,:)
    real(kind=precision), allocatable :: dtef(:,:,:)
    real(kind=precision), allocatable :: gg(:,:,:),ee(:,:,:)
    real(kind=precision), allocatable :: a(:,:,:),c(:,:,:)
    real(kind=precision), allocatable :: l(:,:,:)               

    real(kind=precision) :: tp,sp,p
    real(kind=precision) :: coef1,coef2,coef3,coef4,coef5
    real(kind=precision) :: dt2,sq,sef
    real(kind=precision) :: tbias,sbias,const1,ghc
    integer :: ki
    real(kind=precision) :: a1,b1,a2,b2,c1,e1,e2,e3

    DATA A1,B1,A2,B2,C1/0.92,16.6,0.74,10.1,0.08/
    DATA E1/1.8/,E2/1.33/,E3/1.0/
    DATA SQ/0.20/,SEF/1./

    !---------------------------------------------------------------------------------------------
    ! allocate the momery spaces
    !---------------------------------------------------------------------------------------------
    allocate(GH(kb,IM,JM),SM(kb,IM,JM),SH(kb,IM,JM))
    allocate(PROD(kb,IM,JM),KN(kb,IM,JM),BOYGR(kb,IM,JM))
    allocate(DH(IM,JM),CC(kb,IM,JM),STF(kb,IM,JM))
    allocate(dtef(kb,im,jm))
    allocate(gg(kb,im,jm),ee(kb,im,jm))
    allocate(a(kb,im,jm),c(kb,im,jm))
    allocate(l(kb,im,jm))
      GH = ZERO 
      SM = ZERO
      SH = ZERO
      PROD = ZERO
      KN = ZERO
      BOYGR = ZERO
      DH = ZERO
      CC = ZERO
      STF = ZERO
      dtef = ZERO
      gg = ZERO
      ee = ZERO
      a = ZERO
      c = ZERO
      l = ZERO
  
    small = 1.0e-8
    tbias = 0.
    sbias = 0.
    dt2 = dti
    rho = rho/1000.
    l = 1.0
    
    DO J=1,JM
      DO I=1,IM
        DH(I,J)=hf_inc(i,j)
      enddo
    enddo
    
    DO J=1,JM
      DO I=1,IM
        if (fsm(i,j)==0) cycle
        DO K=2,KBM1
          A(k,I,J) = -DT2*(KQ(ks-k+1,I,J)+KQ(ks-K  ,I,J)+2.*UMOL)*0.5   &
                   /(DZZ(ks-K)*DZ(ks-K)*DH(I,J)*DH(I,J))
          C(k,I,J) = -DT2*(KQ(ks-k+2,I,J)+KQ(ks-k+1,I,J)+2.*UMOL)*0.5   &
                   /(DZZ(ks-K)*DZ(ks-K+1)*DH(I,J)*DH(I,J))
        enddo
      enddo
    enddo
    
    !-----------------------------------------------------------------------------------------------                                                            
    !          THE FOLLOWING SECTION SOLVES THE EQUATION                    
    !          DT2*(KQ*Q2')' - Q2*(2.*DT2*DTEF+1.) = -Q2B                    
    !-----------------------------------------------------------------------------------------------
    !-----  SURFACE AND BOTTOM B.C.S ----------------------------------
    CONST1=16.6**.6666667*SEF
    DO J=1,JMM1
      DO I=1,IMM1
        EE(1,I,J)=0.
        GG(1,I,J)  =SQRT( (.5*(WUSURF(I,J)+WUSURF(I+1,J)))**2            &
                         +(.5*(WVSURF(I,J)+WVSURF(I,J+1)))**2 )*CONST1
        q2f(kb,I,J)=SQRT( (.5*(WUBOT(I,J)+WUBOT(I+1,J)))**2              &
                         +(.5*(WVBOT(I,J)+WVBOT(I,J+1)))**2 )*CONST1
      enddo
    enddo
    !----- Calculate speed of sound squared ----------------------------
    DO J=1,JM
      DO I=1,IM
        if (fsm(i,j)==0) cycle
        DO K=1,KBM1
          TP=Tf(ks-k,I,J)+TBIAS
          SP=Sf(ks-k,I,J)+SBIAS
          !----- Calculate pressure in units of decibars ---------------
          P=-GRAV*1.025*ZZ(ks-K)*DH(I,J)*.1 
          CC(k,I,J)=1449.1+.00821*P+4.55*TP -.045*TP**2         &
                   +1.34*(SP- 35.0)  
          CC(k,I,J)=CC(k,I,J)/SQRT((1.-.01642*P/CC(k,I,J))      &
                   *(1.-0.40*P/CC(k,I,J)**2))
        enddo
      enddo
    enddo
    !----- Calculate buoyancy gradient ---------------------------------
    DO J=1,JM
      DO I=1,IM
        if (fsm(i,j)==0) cycle
        DO K=2,KBM1
          Q2B(k,I,J) =ABS(Q2B(k,I,J))
          Q2LB(k,I,J)=ABS(Q2LB(k,I,J))
          BOYGR(k,I,J)=GRAV*((RHO(ks-K+1,I,J)-RHO(ks-K,I,J))    &
                      /(DZZ(ks-K)*DH(I,J)))                     &
                      +GRAV**2*2.*1.025/(CC(k-1,I,J)**2+CC(k,I,J)**2)
        enddo
      enddo
    enddo
    
    
    DO J=1,JM
      DO I=1,IM
        if (fsm(i,j)==0) cycle
        DO K=2,KBM1
          L(k,I,J)=Q2LB(k,I,J)/Q2B(k,I,J)
          GH(k,I,J)=L(k,I,J)**2/Q2B(k,I,J)*BOYGR(k,I,J)
          GH(k,I,J)=MIN(GH(k,I,J),.028)
        enddo
      enddo
    enddo
    DO J=1,JM
      DO I=1,IM
        L(1,I,J)  = 0.   
        L(kb,I,J) = 0.   
        GH(1,I,J) = 0.  
        GH(kb,I,J)= 0.  
      enddo
    enddo
    
    !------ CALC. T.K.E. PRODUCTION -----------------------------------
    prod = 0.
    DO J=2,JMM1
      DO I=1,IMM1
        if (fsm(i,j)==0) cycle
        DO K=2,KBM1
          PROD(k,I,J)=KM(ks-k+1,I,J)*.25*SEF                      &
              *( (Uf(ks-k,I  ,J)-Uf(ks-k+1,I  ,J)                 &
                 +Uf(ks-k,I+1,J)-Uf(ks-k+1,I+1,J))**2             &
                +(Vf(ks-k,I,J  )-Vf(ks-k+1,I,J  )                 &
                 +Vf(ks-k,I,J+1)-Vf(ks-k+1,I,J+1))**2 )           &
              /(DZZ(ks-K)*DH(I,J))**2
        enddo
      enddo
    enddo
    
    DO J=2,JMM1
      DO I=1,IMM1
        if (fsm(i,j)==0) cycle
        DO K=2,KBM1
          PROD(k,I,J)=PROD(k,I,J)+KH(ks-k+1,I,J)*BOYGR(k,I,J)
        enddo
      enddo
    enddo
    
    GHC=-6.0
    DO J=1,JM
      DO I=1,IM
        if (fsm(i,j)==0) cycle
        DO K=1,KB
          STF(k,I,J)=1.
          IF(GH(k,I,J)<0. ) STF(k,I,J)=1.0-0.9*(GH(k,I,J)/GHC)**1.5
          IF(GH(k,I,J)<GHC) STF(k,I,J)=0.1
          DTEF(k,I,J)=Q2B(k,I,J)*SQRT(Q2B(k,I,J))/(B1*Q2LB(k,I,J)+SMALL)   &
                     *STF(k,I,J)
        enddo
      enddo
    enddo
    
    DO J=1,JM
      DO I=1,IM
        if (fsm(i,j)==0) cycle
        DO K=2,KBM1
          GG(k,I,J)=1./(A(k,I,J)+C(k,I,J)*(1.-EE(k-1,I,J))         &
                   -(2.*DT2*DTEF(k,I,J)+1.) )
          EE(k,I,J)=A(k,I,J)*GG(k,I,J)
          GG(k,I,J)=(-2.*DT2*PROD(k,I,J)                           &
                   +C(k,I,J)*GG(k-1,I,J)-q2F(k,I,J))*GG(k,I,J)
        enddo
      enddo
    enddo
    
    DO J=1,JM
      DO I=1,IM
        if (fsm(i,j)==0) cycle
        DO K=1,KBM1
          KI=KB-K
          q2F(ki,I,J)=EE(ki,I,J)*q2F(ki+1,I,J)+GG(ki,I,J)
        enddo
      enddo
    enddo
    
    !---------------------------------------------------------------------------------------------
    !          THE FOLLOWING SECTION SOLVES THE EQUATION                     *
    !          DT2(KQ*Q2L')' - Q2L*(DT2*DTEF+1.) = -Q2LB                     *
    !---------------------------------------------------------------------------------------------
    DO J=1,JM
      DO I=1,IM
        IF(q2F(2,I,J)<SMALL) q2F(2,I,J)=SMALL
        EE(2,I,J)=0.
        GG(2,I,J)=-KAPPA*Z(ks-2+1)*DH(I,J)*q2F(2,I,J)
        q2lF(kb,I,J)=0.                                        
      enddo
    enddo
    
    DO J=1,JM
      DO I=1,IM
        if (fsm(i,j)==0) cycle
        DO K=3,KBM1 
          DTEF(k,I,J) = DTEF(k,I,J)*(1.+E2*((1./ABS(Z(ks-K+1)-Z(ks-1+1))+ &
                      + 1./ABS(Z(ks-K+1)-Z(ks-KB+1))) *L(k,I,J)/(DH(I,J)*KAPPA))**2)
          GG(k,I,J) = 1./(A(k,I,J)+C(k,I,J)*(1.-EE(k-1,I,J)) &
                    - (DT2*DTEF(k,I,J)+1.))
          EE(k,I,J) = A(k,I,J)*GG(k,I,J)
          GG(k,I,J) = (DT2*(-PROD(k,I,J)  &
                       *L(k,I,J)*E1)+C(k,I,J)*GG(k-1,I,J)-q2lF(k,I,J))*GG(k,I,J)
        enddo
      enddo
    enddo
    
    DO J=1,JM
      DO I=1,IM
        if (fsm(i,j)==0) cycle
        DO K=1,KB-2
          KI=KB-K
          q2lF(ki,I,J)=EE(ki,I,J)*q2lF(ki+1,I,J)+GG(ki,I,J)
        enddo
      enddo
    enddo    
    
    DO J=1,JM
      DO I=1,IM
        if (fsm(i,j)==0) cycle
        DO K=2,KBM1
          IF(q2F(k,I,J)>SMALL.AND.q2lF(k,I,J)>SMALL) cycle
          q2F(k,I,J)  = SMALL
          q2lF(k,I,J) = SMALL
        enddo
      enddo
    enddo
        
    !---------------------------------------------------------------------------------------------
    !      THE FOLLOWING SECTION SOLVES FOR KM AND KH                        
    !---------------------------------------------------------------------------------------------
    COEF4=18.*A1*A1+9.*A1*A2
    COEF5=9.*A1*A2
    ! NOTE THAT SM,SH LIMIT TO INFINITY WHEN GH APPROACHES 0.0288
    DO J=1,JM
      DO I=1,IM
        if (fsm(i,j)==0) cycle
        DO K=1,KB
          COEF1=A2*(1.-6.*A1/B1*STF(k,I,J))
          COEF2=3.*A2*B2/STF(k,I,J)+18.*A1*A2
          COEF3=A1*(1.-3.*C1-6.*A1/B1*STF(k,I,J))
          SH(k,I,J)=COEF1/(1.-COEF2*GH(k,I,J))
          SM(k,I,J)=COEF3+SH(k,I,J)*COEF4*GH(k,I,J)
          SM(k,I,J)=SM(k,I,J)/(1.-COEF5*GH(k,I,J))
        enddo
      enddo
    enddo

    DO J=1,JM
      DO I=1,IM
        if (fsm(i,j)==0) cycle
        DO K=1,KB
          KN(k,I,J)=L(k,I,J)*SQRT(ABS(Q2b(k,I,J)))
          KQ(Ks-k+1,I,J) = ( KN(k,I,J)*.41*SM(k,I,J)    &
                            +KQ(Ks-k+1,I,J) )*.5
          KM(ks-k+1,I,J)=(KN(k,I,J)*SM(k,I,J)+KM(Ks-k+1,I,J))*.5
          KH(Ks-k+1,I,J)=(KN(k,I,J)*SH(k,I,J)+KH(ks-k+1,I,J))*.5
        enddo
      enddo
    enddo
         
    small = 1.0e-10
    
    if (regional_run.eq.0) call obcond_sub('CYC:q2l')
    
    q2b  = q2f
    q2lb = q2lf
    
    rho = rho*1000.

    !---------------------------------------------------------------------------------------------
    ! Incorporate Bv into Km
    !---------------------------------------------------------------------------------------------
    if (load_Bv) then
      Km(ks-ks_bv+1:ks, i0_bv:i1_bv, j0_bv:j1_bv) =       &
      Km(ks-ks_bv+1:ks, i0_bv:i1_bv, j0_bv:j1_bv) + Bv_CN
      Kh(ks-ks_bv+1:ks, i0_bv:i1_bv, j0_bv:j1_bv) =       &
      Kh(ks-ks_bv+1:ks, i0_bv:i1_bv, j0_bv:j1_bv) + Bv_CN
    endif

    !---------------------------------------------------------------------------------------------
    ! Deallocate the momery spaces
    !---------------------------------------------------------------------------------------------
    deallocate(GH,SM,SH)
    deallocate(PROD,KN,BOYGR)
    deallocate(DH,CC,STF)
    deallocate(dtef)
    deallocate(gg,ee)
    deallocate(a,c)
    deallocate(l)
    deallocate(q2f, q2lf)
    
    RETURN   
  END subroutine profq_POM98
  
!=================================================================================================

  !-----------------------------------------------------------------------------------------------
  ! 23)
  ! This compute_psi_sub subroutine computes the barotropic stream function by assuming the zero
  ! stream function set at the NorthWest of the computational domain.
  !
  ! -- Originally developed by Lei HAN (c) August, 2012
  ! -- Code standardized by Lei HAN in September, 2013
  !-----------------------------------------------------------------------------------------------  
  subroutine compute_psi_sub
    
    implicit none
 
    psi = 0.
    !---------------------------------------------------------------------------------------------
    ! Sweep southward
    !---------------------------------------------------------------------------------------------
    do j=jmm1,2,-1
      do i=2,im
        psi(i,j) = psi(i,j+1) + uab(i,j)*dy_u(i,j)
      enddo
    enddo     
  
    !---------------------------------------------------------------------------------------------
    ! Sweep eastward
    !---------------------------------------------------------------------------------------------
    do j=2,jm
      do i=2,imm1
        psi(i+1,j) = psi(i,j) + vab(i,j)*dx_v(i,j)
      enddo
    enddo
  
    !---------------------------------------------------------------------------------------------
    ! Cyclic
    !---------------------------------------------------------------------------------------------
    if (regional_run==0) then
      psi( 1,:) = psi(imm1,:)  
      psi(im,:) = psi(   2,:)  
    endif

    return
  end subroutine compute_psi_sub
  
!=================================================================================================


  !-----------------------------------------------------------------------------------------------
  ! 27)
  ! 
  ! -- Originally developed by Lei HAN (c) August, 2012
  ! -- Code standardized by Lei HAN and Zhanpeng ZHUANG in September, 2013
  !-----------------------------------------------------------------------------------------------  
  subroutine upd_uv_sub
    implicit none
    
    !---------------------------------------------------------------------------------------------
    ! Compute the Internal-mode EPGs
    !---------------------------------------------------------------------------------------------
    do j=2,jmm1
      do i=2,imm1
        epg_int_x(i,j) = -GRAV * (et_mean(i,j) - et_mean(i-1,j))   &
                       * hm_inu(i,j) / dx_u(i,j)
        epg_int_y(i,j) = -GRAV * (et_mean(i,j) - et_mean(i,j-1))   &
                       * hm_inv(i,j) / dy_v(i,j)
      enddo
    enddo
   
    !---------------------------------------------------------------------------------------------
    ! Compute the bottom stress
    !---------------------------------------------------------------------------------------------
    do j=2,jmm1
      do i=2,imm1
        fric_ub(i,j) = (cbc(i,j) + cbc(i-1,j)) * 0.5                                    &
                     * sqrt( ub(1,i,j)**2 + ( 0.25*(vb(1,i  ,j) + vb(1,i  ,j+1)         &
                                                   +vb(1,i-1,j) + vb(1,i-1,j+1)) )**2 )
        fric_vb(i,j) = (cbc(i,j) + cbc(i,j-1)) * 0.5                                    &
                     * sqrt( vb(1,i,j)**2 + ( 0.25*(ub(1,i,j  ) + ub(1,i+1,j  )         &
                                                   +ub(1,i,j-1) + ub(1,i+1,j-1)) )**2 ) 
      enddo
    enddo
    
    if (regional_run==0) then
      !-------------------------------------------------------------------------------------------
      ! Cyclic fric_ub fric_vb for wubot wvbot in the global run
      !-------------------------------------------------------------------------------------------
      fric_ub( 1,:) = fric_ub(imm1,:)
      fric_ub(im,:) = fric_ub(   2,:)
      fric_vb( 1,:) = fric_vb(imm1,:)
      fric_vb(im,:) = fric_vb(   2,:)
    endif   
   
    !---------------------------------------------------------------------------------------------
    ! Solve U,V equations alternatively with Thomas algorithm to achieve an indirect central 
    ! difference for the Coriolis term
    !---------------------------------------------------------------------------------------------
    allocate(fcor_u(ks,im,jm), fcor_v(ks,im,jm))
      fcor_u = ZERO
      fcor_v = ZERO
          
    if (mod(iint+iint0,2)==1) then          ! Update uf first and then vf at odd step
    	
      do j=2,jmm1
        do i=2,imm1
          if (dum(i,j)/=1) cycle
          do k=1,ksm1
            fcor_u(k,i,j) = 0.25 *                                               &
            (Cor_t(i-1,j) * hm_inc(i-1,j) * (vb(k,i-1,j) + vb(k,i-1,j+1)) +      &
             Cor_t(i  ,j) * hm_inc(i  ,j) * (vb(k,i  ,j) + vb(k,i  ,j+1)))       
          enddo
        enddo
      enddo
      
      call solve_tridiag_sub('u')      
      
      if(regional_run==0) then
        call obcond_sub('CYC:uf')         ! Update boundary condition for the global run
      elseif(regional_run==1) then
        call obcond_sub('POM:uf') 
      endif
   
      do j=2,jmm1
        do i=2,imm1
          if (dvm(i,j)/=1) cycle
          do k=1,ksm1
            fcor_v(k,i,j) = -0.25 *                                              &
            (Cor_t(i,j  ) * hm_inc(i,j  ) * (uf(k,i,j  ) + uf(k,i+1,j  )) +      &
             Cor_t(i,j-1) * hm_inc(i,j-1) * (uf(k,i,j-1) + uf(k,i+1,j-1)))
          enddo
        enddo  
      enddo
            
      call solve_tridiag_sub('v')
      
      if(regional_run==0) then
        call obcond_sub('CYC:vf')         ! Update boundary condition for the global run
      elseif(regional_run==1) then
        call obcond_sub('POM:vf') 
      endif
    
    !---------------------------------------------------------------------------------------------
    else                                    ! Update vf first and then uf at even step
    	
      do j=2,jmm1
        do i=2,imm1
          if (dvm(i,j)/=1) cycle
          do k=1,ksm1
            fcor_v(k,i,j) = - 0.25 *                                             &
            (Cor_t(i,j  ) * hm_inc(i,j  ) * (ub(k,i,j  ) + ub(k,i+1,j  )) +      &
             Cor_t(i,j-1) * hm_inc(i,j-1) * (ub(k,i,j-1) + ub(k,i+1,j-1)))
          enddo
        enddo
      enddo
            
      call solve_tridiag_sub('v')
      
      if(regional_run==0) then
        call obcond_sub('CYC:vf')         ! Update boundary condition for the global run
      elseif(regional_run==1) then
        call obcond_sub('POM:vf') 
      endif
            
      do j=2,jmm1
        do i=2,imm1
          if (dum(i,j)/=1) cycle
          do k=1,ksm1
            fcor_u(k,i,j) = 0.25 *                                               &
            (Cor_t(i-1,j) * hm_inc(i-1,j) * (vf(k,i-1,j) + vf(k,i-1,j+1)) +      & 
             Cor_t(i  ,j) * hm_inc(i  ,j) * (vf(k,i  ,j) + vf(k,i  ,j+1)))
          enddo
        enddo
      enddo

      call solve_tridiag_sub('u')
      
      if(regional_run==0) then
        call obcond_sub('CYC:uf')         ! Update boundary condition for the global run
      elseif(regional_run==1) then
        call obcond_sub('POM:uf') 
      endif
      
    endif
    
    !---------------------------------------------------------------------------------------------
    ! Update the bottom friction
    !---------------------------------------------------------------------------------------------
    do j=1,jm
      do i=1,im
        wubot(i,j) = -fric_ub(i,j) * uf(1,i,j)
        wvbot(i,j) = -fric_vb(i,j) * vf(1,i,j)
      enddo
    enddo
   

    deallocate(fcor_u, fcor_v)              ! Release the memory
    deallocate(ub, vb)
    deallocate(ipg_x, ipg_y)
    deallocate(adv_diff_u, adv_diff_v)
    
    return
  end subroutine upd_uv_sub  
!=================================================================================================

  !-----------------------------------------------------------------------------------------------
  ! 28)
  ! 
  ! -- Originally developed by Lei HAN (c) August, 2012
  ! -- Code standardized by Lei HAN and Zhanpeng ZHUANG in September, 2013
  !-----------------------------------------------------------------------------------------------  
  subroutine finalize_INTloop_sub
    implicit none  
  
    allocate(ub(ks,im,jm),vb(ks,im,jm))
    allocate(tb(ks,im,jm),sb(ks,im,jm))
      ub = ZERO
      vb = ZERO
      tb = ZERO
      sb = ZERO
    
    !---------------------------------------------------------------------------------------------
    ! Step forward the model variables for the next iteration
    !---------------------------------------------------------------------------------------------
    zetb_c = zetf_c
    hb_inc = hf_inc
    hb_inu = hf_inu
    hb_inv = hf_inv


    tb = tf
    sb = sf
    ub = uf
    vb = vf

    
    !---------------------------------------------------------------------------------------------
    ! Output yearly restart file for future use
    !---------------------------------------------------------------------------------------------
    if ( IF_save_restart==1 .and. ymdhms(1)/=ymdhms_prev(1) ) then
      open(1,file='restart_'//date_char(1:8)//'.bin'            &
           ,form='unformatted', status='replace')
      write(1)                                       &
         iint, days_now                              &
        ,ua,va,uab,vab,etb_c,zetb_c,wubot,wvbot      &
        ,tb,sb,ub,vb,wf,q2b,q2lb,Km,Kh,Kq,rho        
      close(1)
      print*,'||||| Yearly restart file saved at iint ', iint  &
          ,' or days = ', days_now, ' ||||'
    endif
    
    !---------------------------------------------------------------------------------------------
    ! Output monthly restart file in case this run stops
    !---------------------------------------------------------------------------------------------
!    if ( IF_save_restart==1 .and. ymdhms(2)/=ymdhms_prev(2) ) then
!      open(1,file='restart_realtime_monthly.bin'            &
!           ,form='unformatted', status='replace')
!      write(1)                                       &
!         iint, days_now                              &
!        ,ua,va,uab,vab,etb_c,zetb_c,wubot,wvbot      &
!        ,tb,sb,ub,vb,wf,q2b,q2lb,Km,Kh,Kq,rho        
!      close(1)
!      print*,'||||| RT restart file saved at iint ', iint  &
!          ,' or days = ', days_now, ' ||||'
!    endif
      
   
    !---------------------------------------------------------------------------------------------
    ! Sum up the whole KE within the computational domain
    !---------------------------------------------------------------------------------------------
    if ( ymdhms(3) /= ymdhms_prev(3) ) then
      total_KE = 0.
      do j=2,jmm1
        do i=2,imm1
          if (fsm(i,j)==0) cycle
          do k=1,ksm1
            total_KE = total_KE                                        &
                     + 0.125 * (1000. + rho(k,i,j))                    &  
                             * art(i,j) * hf_inc(i,j) * dz(k)          &
                             * ( (ub(k,i,j) + ub(k,i+1,j))**2.0        &
                                +(vb(k,i,j) + vb(k,i,j+1))**2.0 )      
          enddo
        enddo
      enddo
      write(2012,*) iint,days_now,total_KE 
    endif         
   
       
    !---------------------------------------------------------------------------------------------
    ! Internal mode integration step completed
    !---------------------------------------------------------------------------------------------
    if (mod(iint,500).eq.0) then
      print*, 'maxval(uf)', maxloc(uf), maxval(uf)
      print*, 'minval(uf)', minloc(uf), minval(uf)
      print*, 'maxval(vf)', maxloc(vf), maxval(vf)
      print*, 'minval(vf)', minloc(vf), minval(vf)
      print*, 'maxval(tf)', maxloc(tb), maxval(tb)
      print*, 'minval(tf)', minloc(tb), minval(tb)
      print*, 'maxval(sst)', maxloc(sst_dq), maxval(sst_dq)
      print*, 'minval(sst)', minloc(sst_dq), minval(sst_dq)
      print*, 'maxval(sf)', maxloc(sb), maxval(sb)
      print*, 'minval(sf)', minloc(sb), minval(sb)
      print*, 'maxval(ua)', maxloc(ua), maxval(ua)
      print*, 'minval(ua)', minloc(ua), minval(ua)
      print*, 'maxval(va)', maxloc(va), maxval(va)
      print*, 'minval(va)', minloc(va), minval(va)
      print*, 'maxval(wf)', maxloc(wf), maxval(wf)
      print*, 'minval(wf)', minloc(wf), minval(wf)
      print*, 'maxval(rho)',maxloc(rho),maxval(rho)
      print*, 'minval(rho)',minloc(rho),minval(rho)
    
    endif
    
!      write(17,*) '----------------------------------'
!      write(17,*) 'iint=',iint,'day=',date_char(1:8)
!      write(17,*) 'maxval(uf)', maxloc(uf), maxval(uf)
!      write(17,*) 'minval(uf)', minloc(uf), minval(uf)
!      write(17,*) 'maxval(vf)', maxloc(vf), maxval(vf)
!      write(17,*) 'minval(vf)', minloc(vf), minval(vf)
!      write(17,*) 'maxval(tf)', maxloc(tb), maxval(tb)
!      write(17,*) 'minval(tf)', minloc(tb), minval(tb)
!      write(17,*) 'maxval(sst)', maxloc(sst_dq), maxval(sst_dq)
!      write(17,*) 'minval(sst)', minloc(sst_dq), minval(sst_dq)
!      write(17,*) 'maxval(sf)', maxloc(sb), maxval(sb)
!      write(17,*) 'minval(sf)', minloc(sb), minval(sb)
!      write(17,*) 'maxval(ua)', maxloc(ua), maxval(ua)
!      write(17,*) 'minval(ua)', minloc(ua), minval(ua)
!      write(17,*) 'maxval(va)', maxloc(va), maxval(va)
!      write(17,*) 'minval(va)', minloc(va), minval(va)
!      write(17,*) 'maxval(wf)', maxloc(wf), maxval(wf)
!      write(17,*) 'minval(wf)', minloc(wf), minval(wf)
!      write(17,*) 'maxval(rho)',maxloc(rho),maxval(rho)
!      write(17,*) 'minval(rho)',minloc(rho),minval(rho)
    
    
    
    deallocate(uf,vf)
    deallocate(tf,sf)
     
    return
  end subroutine
!=================================================================================================

  subroutine smooth_spatial_sub(var_name)
    implicit none
  
    character(len=*), intent(in) :: var_name
    real(kind=precision), pointer :: filtered_var(:,:)
    real(kind=precision), pointer :: mask(:,:)
    real(kind=precision), allocatable :: tmp_var(:,:)
    real(kind=precision) :: waters_around,tmp
    integer dual_smooth_index
    
    allocate( tmp_var(im,jm) )
    
    if (var_name=='etf') then
    	filtered_var => etf_c
    	mask => fsm
    elseif (var_name=='uaf') then
    	filtered_var => uaf
    	mask => dum01
    elseif (var_name=='vaf') then
    	filtered_var => vaf
    	mask => dvm01
    elseif (var_name=='H_c') then
    	filtered_var => H_c
    	mask => fsm
    endif    	
        
    do dual_smooth_index=1,-1,-2
    	if (var_name=='H_c' .and. dual_smooth_index==-1) exit
    	
      tmp_var = filtered_var * mask
      if (regional_run) then
        tmp_var( 1,:) = tmp_var(   2,:)
        tmp_var(im,:) = tmp_var(imm1,:)
        tmp_var(:, 1) = tmp_var(:,   2)
        tmp_var(:,jm) = tmp_var(:,jmm1)
        tmp_var = tmp_var * mask
      endif
      do j=2,jmm1
        do i=2,imm1
          waters_around = mask(i+1,j)+mask(i-1,j)+mask(i,j+1)+mask(i,j-1)
          if (waters_around <= 0) cycle
          tmp = tmp_var(i+1,j)+tmp_var(i,j+1)   &
              + tmp_var(i-1,j)+tmp_var(i,j-1)-waters_around*tmp_var(i,j)
          filtered_var(i,j) = ( tmp_var(i,j) + dual_smooth_index*smooth_factor*tmp ) &
                            * mask(i,j)
        enddo
      enddo
    enddo
    
    nullify(filtered_var, mask)
    
    return
  end subroutine
!=================================================================================================

  subroutine record_time_series_sub
    implicit none
    
    real(kind=precision) :: lon_sp(90), lat_sp(90)
    real(kind=precision) :: lon_dt(im), lat_dt(jm)
                	    
    if (iint==1) then
      if (ks==16) then
        ly_id = (/1,4,8,11,15/)
      elseif (ks==21) then
        ly_id = (/1,6,11,16,20/)
      endif
	    if (regional_run) then
	      pts = 30
	      open(1,file='MP30.txt',action='read')
	      print*,'Time series point information loaded from <MP30_v2.txt>'
	    else
	      pts = 20
	      open(1,file='MP30.txt',action='read')
	      print*,'Time series point information loaded from <MP20_v2.txt>'
	    endif
!	      read(1,'(2f10.4)') (lon_sp(k),lat_sp(k),k=1,pts)

      	read(1,'(2i5)') (ix(k),jx(k),k=1,30)
      	do n=1,30
      		print*,ix(n),jx(n)
      	enddo

	      close(1)
	    cols = pts * lys

    !---------------------------------------------------------------------------------------------
    ! Find corresponding index points
    !---------------------------------------------------------------------------------------------
    !do i=1,pts
    !  lon_dt = vlon - lon_sp(i)
    !  ix(i) = maxloc(lon_dt, dim=1, mask=lon_dt<=0)
    !  lat_dt = vlat - lat_sp(i)
    !  jx(i) = maxloc(lat_dt, dim=1, mask=lat_dt<=0)
    !enddo

	    open(22,file=    'res/sdx_series_V2.15.txt')
	    open(23,file=    'res/sdy_series_V2.15.txt')	
	    open(21,file=    'res/temp_series_V2.15.txt')
	    open(39,file=      'res/ua_series_V2.15.txt')
	    open(52,file=     'res/etf_series_V2.15.txt')
	    open(53,file=     'res/zet_series_V2.15.txt')
	    open(54,file=      'res/va_series_V2.15.txt')
	    open(79,file=      'res/vf_series_V2.15.txt')
	    open(88,file=      'res/uf_series_V2.15.txt')
	    open(90,file='res/salinity_series_V2.15.txt')
	    open(103,file=     'res/wf_series_V2.15.txt')
	    open(104,file=     'res/Km_series_V2.15.txt')
	    open(102,file=     'res/Kh_series_V2.15.txt')
	    open(105,file=   'res/advt_series_V2.15.txt')
	    open(106,file=  'res/difft_series_V2.15.txt')
	    open(107,file=     'res/Sm_series_V2.15.txt')
	    open(108,file=    'res/q2b_series_V2.15.txt')
	    open(109,file=   'res/q2lb_series_V2.15.txt')
	    open(110,file=     'res/Gh_series_V2.15.txt')
	    open(111,file=  'res/wusurf_series_V2.15.txt')
	    open(112,file=  'res/wvsurf_series_V2.15.txt')
	    open(113,file=  'res/advu2_series_V2.15.txt')
    endif
    	
    write(52,'(<pts>e15.5)') (etf_c(ix(m),jx(m)),m=1,pts)
    write(39,'(<pts>e15.5)') (ua(ix(m),jx(m)),m=1,pts)
    write(54,'(<pts>e15.5)') (va(ix(m),jx(m)),m=1,pts)
    write(53,'(<pts>e15.5)') (zetf_c(ix(m),jx(m)),m=1,pts)
	  write(21,'(<cols>e15.5)')                                &
      ((tb(ly_id(k),ix(m),jx(m)),k=1,lys),m=1,pts)
    write(90,'(<cols>e15.5)')                                &
      ((sb(ly_id(k),ix(m),jx(m)),k=1,lys),m=1,pts)
    write(79,'(<cols>e15.5)')                                &
      ((vb(ly_id(k),ix(m),jx(m)),k=1,lys),m=1,pts)
    write(88,'(<cols>e15.5)')                                &
      ((ub(ly_id(k),ix(m),jx(m)),k=1,lys),m=1,pts)
    write(103,'(<cols>e15.5)')                               &
      ((wf(ly_id(k),ix(m),jx(m)),k=1,lys),m=1,pts)
    write(104,'(<cols>e15.5)')                               &
      ((Km(ly_id(k),ix(m),jx(m)),k=1,lys),m=1,pts)
    write(102,'(<cols>e15.5)')                               &
      ((Kh(ly_id(k),ix(m),jx(m)),k=1,lys),m=1,pts)
    write(108,'(<cols>e15.5)')                               &
      ((q2b(ly_id(k),ix(m),jx(m)),k=1,lys),m=1,pts)
    write(109,'(<cols>e15.5)')                               &
      ((q2lb(ly_id(k),ix(m),jx(m)),k=1,lys),m=1,pts)
      
    write(111,'(<pts>e15.5)') (wusurf(ix(m),jx(m)),m=1,pts)
    write(112,'(<pts>e15.5)') (wvsurf(ix(m),jx(m)),m=1,pts)

    
    return
  end subroutine record_time_series_sub


!=================================================================================================

  
  subroutine output_clim_data
  
    use netcdf_mod
    implicit none
    
    integer ncid
    character(len=200) ncfile
          
    !---------------------------------------------------------------------------------------------
    ! Output the data to netcdf file
    !---------------------------------------------------------------------------------------------
    if ( ymdhms(2) /= ymdhms_prev(2) ) then
      
      !---------------------------------------------------------------------------------------------
      ! Average the accumulated sum
      !---------------------------------------------------------------------------------------------
      ssh_mon = ssh_mon / n_clim
      tb_mon  = tb_mon  / n_clim
      sb_mon  = sb_mon  / n_clim
      ub_mon  = ub_mon  / n_clim
      vb_mon  = vb_mon  / n_clim
      
      !---------------------------------------------------------------------------------------------
      ! Define missing value
      !---------------------------------------------------------------------------------------------
      do j=1,jm
        do i=1,im
          if (fsm(i,j)==0) then
            ssh_mon(i,j) = missing_value
            do k=1,ks
              tb_mon(k,i,j) = missing_value
              sb_mon(k,i,j) = missing_value
            enddo
          endif
          if (dum(i,j)==0) then
            do k=1,ks
              ub_mon(k,i,j) = missing_value
            enddo
          endif
          if (dvm(i,j)==0) then
            do k=1,ks
              vb_mon(k,i,j) = missing_value
            enddo
          endif
        enddo
      enddo  
                
      ncfile = 'monthly_clim_'//date_char_prev(1:6)//'_masnum.nc'
      print*,'Outputing NC monthly file: ',trim(adjustl(ncfile))

      call open_nc(ncid,trim(adjustl(ncfile)),'c')
      
      call dimension_define(ncid, 'sigma', ks, 'sigma', 5)
      call dimension_define(ncid, 'lon', im, 'lon', 5)
      call dimension_define(ncid, 'lat', jm, 'lat', 5)
      
      call set_attribute(ncid, 'units', 'sigma layers',  'sigma')
      call set_attribute(ncid, 'units', 'degrees_north', 'lat')
      call set_attribute(ncid, 'units', 'degrees_east',  'lon')

      call variable_define(ncid,'fsm',5, ['lon','lat'])
      call variable_define(ncid,'dum',5, ['lon','lat'])
      call variable_define(ncid,'dvm',5, ['lon','lat'])
    
      call variable_define(ncid,'depth',5, ['lon','lat'])
      call variable_define(ncid,'ssh'  ,5, ['lon','lat'])
      call variable_define(ncid,'temp' ,5, ['sigma','lon','lat'])
      call variable_define(ncid,'salt' ,5, ['sigma','lon','lat'])
      call variable_define(ncid,'u'    ,5, ['sigma','lon','lat'])
      call variable_define(ncid,'v'    ,5, ['sigma','lon','lat'])
      
      call set_attribute(ncid, 'missing_value', missing_value, 'ssh')
      call set_attribute(ncid, 'missing_value', missing_value, 'temp')
      call set_attribute(ncid, 'missing_value', missing_value, 'salt')
      call set_attribute(ncid, 'missing_value', missing_value, 'u')
      call set_attribute(ncid, 'missing_value', missing_value, 'v')
      
      call end_define(ncid)

      call writenc(ncid,'sigma',zz)
      call writenc(ncid,'lon',vlon)
      call writenc(ncid,'lat',vlat)

      call writenc(ncid,'fsm',fsm) 
      call writenc(ncid,'dum',dum01) 
      call writenc(ncid,'dvm',dvm01) 

      call writenc(ncid,'depth',H_c) 
      call writenc(ncid,'ssh',ssh_mon) 
      call writenc(ncid,'temp',tb_mon) 
      call writenc(ncid,'salt',sb_mon) 
      call writenc(ncid,'u',ub_mon) 
      call writenc(ncid,'v',vb_mon) 
          
      call close_nc(ncid)
      
      ssh_mon = ZERO
      tb_mon  = ZERO
      sb_mon  = ZERO
      ub_mon  = ZERO
      vb_mon  = ZERO
      n_clim  = ZERO
      
    endif
    
    !---------------------------------------------------------------------------------------------
    ! Accumulate the monthly data
    !---------------------------------------------------------------------------------------------
    n_clim = n_clim + 1
    do j=1,jm
      do i=1,im
        ssh_mon(i,j) = ssh_mon(i,j) + zetf_c(i,j)
        do k=1,ks
          tb_mon(k,i,j) = tb_mon(k,i,j) + tb(k,i,j)
          sb_mon(k,i,j) = sb_mon(k,i,j) + sb(k,i,j)
          ub_mon(k,i,j) = ub_mon(k,i,j) + ub(k,i,j)
          vb_mon(k,i,j) = vb_mon(k,i,j) + vb(k,i,j)
        enddo
      enddo
    enddo
    
  
  return
  end subroutine output_clim_data
  
!=================================================================================================

  
  subroutine output_instant_data
  
    use netcdf_mod
    use time_mod
    
    implicit none

    integer :: ncid
    integer :: null_value
    character(len=200) ncfile_instant
    real(kind=4) :: scale_factor, add_offset
    integer, allocatable :: var2d(:,:), var3d(:,:,:)
    
    scale_factor = 0.001
    add_offset = 20
    null_value = 32767
    
    allocate(var2d(im,jm),var3d(ks,im,jm))
    
    !---------------------------------------------------------------------------------------------
    ! Output latest 7days' data into netcdf file
    !---------------------------------------------------------------------------------------------
    if (nth_field==1) then
            
      ncfile_instant = 'instant_output/daily_'//date_char(1:8)//'_masnum.nc'
      print*,'outputing NC instant file: ',trim(adjustl(ncfile_instant))
      
      call open_nc(ncid,trim(adjustl(ncfile_instant)),'c')

      call dimension_define(ncid, 'sigma', ks, 'sigma', 5)
      call dimension_define(ncid, 'lon', im, 'lon', 5)
      call dimension_define(ncid, 'lat', jm, 'lat', 5)
      call dimension_define(ncid, 'time', 0, 'time', 6)
      
      call set_attribute(ncid, 'units', 'sigma layers',  'sigma')
      call set_attribute(ncid, 'units', 'degrees_north', 'lat')
      call set_attribute(ncid, 'units', 'degrees_east',  'lon')
      call set_attribute(ncid, 'units', 'days since 0-1-1 00:00', 'time')
      call set_attribute(ncid, 'starting_date', date_char(1:8))
        
      call variable_define(ncid,'fsm', 3, ['lon','lat'])
      call variable_define(ncid,'dum', 3, ['lon','lat'])
      call variable_define(ncid,'dvm', 3, ['lon','lat'])
      
      call variable_define(ncid, 'Hc', 5, ['lon','lat'])
      call variable_define(ncid, 'Hu', 5, ['lon','lat'])
      call variable_define(ncid, 'Hv', 5, ['lon','lat'])
      
      call variable_define(ncid,'ssh'  , 3, ['lon','lat','time'])
      call variable_define(ncid,'ua'   , 3, ['lon','lat','time'])
      call variable_define(ncid,'va'   , 3, ['lon','lat','time'])
      call variable_define(ncid,'temp' , 3, ['sigma','lon','lat','time'])
      call variable_define(ncid,'salt' , 3, ['sigma','lon','lat','time'])
      call variable_define(ncid,'u'    , 3, ['sigma','lon','lat','time'])
      call variable_define(ncid,'v'    , 3, ['sigma','lon','lat','time'])

      call set_attribute(ncid,'scale_factor',scale_factor,'ssh')
      call set_attribute(ncid,'scale_factor',scale_factor,'ua')
      call set_attribute(ncid,'scale_factor',scale_factor,'va')
      call set_attribute(ncid,'scale_factor',scale_factor,'u')
      call set_attribute(ncid,'scale_factor',scale_factor,'v')
      call set_attribute(ncid,'scale_factor',scale_factor,'temp')
      call set_attribute(ncid,'scale_factor',scale_factor,'salt')

      call set_attribute(ncid,'missing_value',null_value,'ssh')
      call set_attribute(ncid,'missing_value',null_value,'ua')
      call set_attribute(ncid,'missing_value',null_value,'va')
      call set_attribute(ncid,'missing_value',null_value,'u')
      call set_attribute(ncid,'missing_value',null_value,'v')
      call set_attribute(ncid,'missing_value',null_value,'temp')
      call set_attribute(ncid,'missing_value',null_value,'salt')

      call set_attribute(ncid,'add_offset',add_offset,'temp')
      call set_attribute(ncid,'add_offset',add_offset,'salt')

      call end_define(ncid)
    
      call writenc(ncid,'sigma' , zz   )
      call writenc(ncid,'lon'   , vlon )
      call writenc(ncid,'lat'   , vlat )
      call writenc(ncid,'fsm'   , fsm  )
      call writenc(ncid,'dum'   , dum01)
      call writenc(ncid,'dvm'   , dvm01)
      call writenc(ncid,'Hc'    , H_c  )
      call writenc(ncid,'Hu'    , H_u  )
      call writenc(ncid,'Hv'    , H_v  )
      
    else
    
      call open_nc(ncid,trim(adjustl(ncfile_instant)),'w')
    
    endif
      
    call writenc(ncid, 'time', days_now, nth_field) 

    !---------------------------------------------------------------------------------------------
    ! Export as short INT (2D variables)
    !---------------------------------------------------------------------------------------------    
    do j=1,jm
    	do i=1,im
    		var2d(i,j) = zetf_c(i,j)/scale_factor
    		if (fsm(i,j)==0) var2d(i,j) = null_value
      enddo
    enddo
    call writenc(ncid, 'ssh', var2d, nth_field) 

    do j=1,jm
    	do i=1,im
    		var2d(i,j) = ua(i,j)/scale_factor
    		if (dum01(i,j)==0) var2d(i,j) = null_value
      enddo
    enddo
    call writenc(ncid, 'ua', var2d, nth_field) 

    do j=1,jm
    	do i=1,im
    		var2d(i,j) = va(i,j)/scale_factor
    		if (dvm01(i,j)==0) var2d(i,j) = null_value
      enddo
    enddo
    call writenc(ncid, 'va', var2d, nth_field) 

    !---------------------------------------------------------------------------------------------
    ! Export as short INT (3D variables)
    !---------------------------------------------------------------------------------------------    
    do j=1,jm
    	do i=1,im
    		do k=1,ks
    		  var3d(k,i,j) = ub(k,i,j)/scale_factor
        enddo
      		if (dum01(i,j)==0) var3d(:,i,j) = null_value
      enddo
    enddo
    call writenc(ncid, 'u', var3d, nth_field)     

    do j=1,jm
    	do i=1,im
    		do k=1,ks
    		  var3d(k,i,j) = vb(k,i,j)/scale_factor
        enddo
          if (dvm01(i,j)==0) var3d(:,i,j) = null_value
      enddo
    enddo
    call writenc(ncid, 'v', var3d, nth_field) 

    do j=1,jm
    	do i=1,im
    		do k=1,ks
    		  var3d(k,i,j) = (tb(k,i,j)-add_offset)/scale_factor
        enddo
          if (fsm(i,j)==0) var3d(:,i,j) = null_value
      enddo
    enddo
    call writenc(ncid, 'temp', var3d, nth_field) 
    print*,'Exported temp-var3d (short INT)',minval(var3d),maxval(var3d)

    do j=1,jm
    	do i=1,im
    		do k=1,ks
    		  var3d(k,i,j) = (sb(k,i,j)-add_offset)/scale_factor
        enddo
          if (fsm(i,j)==0) var3d(:,i,j) = null_value
      enddo
    enddo
    call writenc(ncid, 'salt', var3d, nth_field)
    print*,'Exported salt-var3d(short INT)',minval(var3d),maxval(var3d)

    call close_nc(ncid)
          
    !---------------------------------------------------------------------------------------------
    ! If all the time fields of one data file is completed, reset the export.
    !---------------------------------------------------------------------------------------------    
    if (nth_field==24/interval_export_instant) then
      print*,'Successfully output daily NC file: ', ncfile_instant
      nth_field = 0
      ncfile_instant = ''
      return
    endif

    nth_field = nth_field + 1
    print*,'nth_field = ',nth_field
      
    return
  end subroutine output_instant_data
  
!======================================================================  

  
  subroutine output_instant_ts
  
    use netcdf_mod
    use time_mod
    
    implicit none

    integer ncid
    character(len=200) ncfile_instant
    real(kind=8) :: dtime
    logical existed
    
    dtime = days_now - datenum(1950,1,1,0,0,0)       ! Argo days unit

    !---------------------------------------------------------------------------------------------
    ! Delete the flag file of the previous EAKF assimilation job
    !---------------------------------------------------------------------------------------------    
    if (ens_char=='01') then
    	existed = 0
      inquire(file=trim(adjustl(EAKF_dir))//'/EAKF_job_done.txt',exist=existed)
      if (existed) then
        open(1,file=trim(adjustl(EAKF_dir))//'/EAKF_job_done.txt')
        close(1,status='delete')
        print*,'EAKF_job_done.txt is deleted'
      endif
    endif
                    
    !---------------------------------------------------------------------------------------------
    ! Output the model data for EAKF process
    !---------------------------------------------------------------------------------------------    
    ncfile_instant = trim(adjustl(EAKF_dir))//'/ts_ens'//ens_char//'_masnum.nc'
    print*,'outputing NC instant TS file: ',trim(adjustl(ncfile_instant))
    
    call open_nc(ncid,trim(adjustl(ncfile_instant)),'c')

    call dimension_define(ncid, 'sigma', ks, 'sigma', 5)
    call dimension_define(ncid, 'lon',   im, 'lon',   5)
    call dimension_define(ncid, 'lat',   jm, 'lat',   5)
    call dimension_define(ncid, 'time',   1, 'time',  6)
    
    call set_attribute(ncid, 'units', 'sigma_layers',  'sigma')
    call set_attribute(ncid, 'units', 'degrees_north', 'lat')
    call set_attribute(ncid, 'units', 'degrees_east',  'lon')
    call set_attribute(ncid, 'units', 'days since 1950-1-1 00:00', 'time')
    call set_attribute(ncid, 'date' , date_char, 'time')
    print*,'datestr(dtime)',date_char
      
    call variable_define(ncid,'fsm', 5, ['lon','lat'])
    
    call variable_define(ncid,'depth', 5, ['lon','lat'])
    call variable_define(ncid,'temp' , 5, ['sigma','lon','lat'])
    call variable_define(ncid,'salt' , 5, ['sigma','lon','lat'])
    call end_define(ncid)
    
    call writenc(ncid, 'sigma' , zz   )
    call writenc(ncid, 'lon'   , vlon )
    call writenc(ncid, 'lat'   , vlat )
    call writenc(ncid, 'time'  , dtime)
    call writenc(ncid, 'fsm'   , fsm  )
    call writenc(ncid, 'depth' , H_c  )
    call writenc(ncid, 'temp'  , tb   ) 
    call writenc(ncid, 'salt'  , sb   ) 

    call close_nc(ncid)
          
    print*,'outputing daily instantenous TS at ',date_char
    
    !---------------------------------------------------------------------------------------------
    ! Write flag files implying this ensemble output file is ready
    !---------------------------------------------------------------------------------------------
    open(1,file=trim(adjustl(EAKF_dir))//'/ts_ens'//ens_char//'_masnum.ready')
    write(1,*) date_char
    close(1)
    
    !---------------------------------------------------------------------------------------------
    ! Wait for loading the EAKFed data with a dead-loop inquiring
    !---------------------------------------------------------------------------------------------
    existed = 0
    if ( IF_EAKF_argo ) then
    	print*,'Waiting for EAKF adjusted file to continue......'
    	do while (.not.existed)
    		inquire( file=trim(EAKF_dir)//'/EAKF_job_done.txt', exist=existed )
      enddo
      EAKF_adjusted_file = trim(EAKF_dir)//'/ts_ens'//ens_char//'_EAKFed.nc'
      call open_nc(ncid,trim(adjustl(EAKF_adjusted_file)),'r')
      n_param = get_dimension_len(ncid,'n_param')
      call readnc(ncid, 'temp', tb ) 
      if (n_param>2)  call readnc(ncid, 'salt', sb ) 
      call close_nc(ncid)
      
      ! Delete after finish loading
      open(1,file=trim(adjustl(EAKF_adjusted_file)))
      close(1,status='delete')
      print*,'EAKF data file ',trim(adjustl(EAKF_adjusted_file)),' is loaded and then deleted'
            	
      print*,'----EAKFed TS:'
      print*, 'maxval(tf)', maxloc(tb), maxval(tb)
      print*, 'minval(tf)', minloc(tb), minval(tb)
      print*, 'maxval(sf)', maxloc(sb), maxval(sb)
      print*, 'minval(sf)', minloc(sb), minval(sb)
      print*,'-----------'
      
      stop
      
    endif
    
    return
  end subroutine output_instant_ts
  
!======================================================================  

  subroutine slpmax
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Limits the maximum of:                              *
! *                                                                    *
! *                  <difference of depths>/<sum of depths>            *
! *                                                                    *
! *                for two adjacent cells. The maximum possible value  *
! *                is unity.                                           *
! *                                                                    *
! **********************************************************************
!
    implicit none
    
    real(kind=precision):: mean,del,slmax
    integer loop
    
    slmax = 0.2
    
    do loop=1,10
!     Sweep right:

        do j=2,jm-1
        	
          do i=2,im-1
            if(fsm(i,j)/=0 .and. fsm(i+1,j)/=0) then
              if( abs(H_c(i+1,j)-H_c(i,j))/(H_c(i,j)+H_c(i+1,j))>=slmax ) then
                mean=(H_c(i+1,j)+H_c(i,j))/2.e0
                del=sign(slmax,H_c(i+1,j)-H_c(i,j))
                H_c(i+1,j)=mean*(1.e0+del)
                H_c(i,j)=mean*(1.e0-del)
              endif
            endif
          end do

!    Sweep left:

          do i=im-1,2,-1
            if( fsm(i,j)/=0 .and. fsm(i+1,j)/=0 ) then
              if(abs(H_c(i+1,j)-H_c(i,j))/(H_c(i,j)+H_c(i+1,j))>=slmax) then
                mean=(H_c(i+1,j)+H_c(i,j))/2.e0
                del=sign(slmax,H_c(i+1,j)-H_c(i,j))
                H_c(i+1,j)=mean*(1.e0+del)
                H_c(i,j)=mean*(1.e0-del)
              endif
            endif
          end do

        end do

!   Sweep up:

        do i=2,im-1

          do j=2,jm-1
            if( fsm(i,j)/=0 .and. fsm(i,j+1)/=0 ) then
              if(abs(H_c(i,j+1)-H_c(i,j))/(H_c(i,j)+H_c(i,j+1))>=slmax) then
                mean=(H_c(i,j+1)+H_c(i,j))/2.e0
                del=sign(slmax,H_c(i,j+1)-H_c(i,j))
                H_c(i,j+1)=mean*(1.e0+del)
                H_c(i,j)=mean*(1.e0-del)
              endif
            endif
          enddo

!   Sweep down:

          do j=jm-1,2,-1
            if(fsm(i,j)/=0 .and. fsm(i,j+1)/=0) then
              if(abs(H_c(i,j+1)-H_c(i,j))/(H_c(i,j)+H_c(i,j+1))>=slmax) then
                mean=(H_c(i,j+1)+H_c(i,j))/2.e0
                del=sign(slmax,H_c(i,j+1)-H_c(i,j))
                H_c(i,j+1)=mean*(1.e0+del)
                H_c(i,j)=mean*(1.e0-del)
              endif
            endif
          enddo

        enddo

    enddo

    return

  end subroutine slpmax 
  
!=================================================================================================

  subroutine upd_open_bound(redef_boundpts)
    
    use netcdf_mod
    use time_mod
    
    implicit none
    
    integer, optional :: redef_boundpts
    integer :: if_upd_numbound
    
    integer ncid, ncid2
    character(len=200) ncfile_glb,ncfile_glb2
    
    integer lon_len, lat_len, ks_len
    real(kind=precision), allocatable, target :: lon_t_glb(:), lat_t_glb(:), zz_glb(:)
    real(kind=precision), allocatable, target :: lon_u_glb(:), lat_v_glb(:)
    real(kind=precision), pointer :: lon_glb(:), lat_glb(:)
    real(kind=precision), pointer :: lon_reg(:), lat_reg(:)

    real(kind=precision) dlon_glb, dlat_glb
    real(kind=precision) add_offset, scale_factor, null_value

    real(kind=precision), allocatable, target :: ssh_glb(:,:), ua_glb(:,:), va_glb(:,:)
    real(kind=precision), allocatable, target :: tb_glb(:,:,:), sb_glb(:,:,:)
    real(kind=precision), allocatable, target :: ub_glb(:,:,:), vb_glb(:,:,:)

    real(kind=precision), pointer :: H_glb(:,:), H_reg(:,:)
    real(kind=precision), allocatable, target :: Hc_glb(:,:), Hu_glb(:,:), Hv_glb(:,:)

    real(kind=precision), allocatable, target :: fsm_glb(:,:)
    real(kind=precision), allocatable, target :: dum_glb(:,:), dvm_glb(:,:)
    real(kind=precision), pointer :: mask_glb(:,:)

    integer, allocatable, target :: iloc_t(:), jloc_t(:), iloc_u(:), jloc_v(:)
    integer, pointer :: iloc(:), jloc(:)

    real(kind=precision), allocatable  :: lon_dt(:), lat_dt(:)
    
    real(kind=precision), pointer :: mask(:,:), var2d(:,:), var3d(:,:,:)
    real(kind=precision), pointer :: var2d_obc(:), var3d_obc(:,:)
    
    real(kind=precision), allocatable :: zm(:,:)
    real(kind=precision) :: za
    
    integer irec, dim_num
    real(kind=precision) :: x,y,x1,x3,y1,y3,fin(4),fx
    integer :: i1, i3, j1, j3
    character(len=4) var_char(7)
    real(kind=precision) :: tps_u, tps_v
    
    null_value = 999.
    
    !---------------------------------------------------------------------------------------------
    ! The data file is documented in every 6 hours
    !---------------------------------------------------------------------------------------------
    irec = int(ymdhms(4)/interval_export_instant) + 1
    
    ncfile_glb  = trim(adjustl(glb_dir))//'2012'//date_char(5:8)//'masnum.nc'
    ncfile_glb2 = trim(adjustl(glb_dir))//'common_var.nc'
    print*,'Loading Global Run data ', trim(ncfile_glb), ' at irec=',irec
    
    
    call open_nc(ncid, ncfile_glb, 'r')
    call open_nc(ncid2, ncfile_glb2, 'r')

    lon_len  = get_dimension_len(ncid,'lon')
    lat_len  = get_dimension_len(ncid,'lat')
    ks_len   = get_dimension_len(ncid,'layer')
!    print*,'lon_len',lon_len,lat_len,ks_len
    
    allocate( lon_t_glb(lon_len), lat_t_glb(lat_len), zz_glb(ks_len) )
    allocate( lon_u_glb(lon_len), lat_v_glb(lat_len) )
    
    allocate(  Hc_glb(lon_len,lat_len),  Hu_glb(lon_len,lat_len),  Hv_glb(lon_len,lat_len) )
    allocate( fsm_glb(lon_len,lat_len), dum_glb(lon_len,lat_len), dvm_glb(lon_len,lat_len) )

    allocate( ssh_glb(lon_len,lat_len), ua_glb(lon_len,lat_len), va_glb(lon_len,lat_len) )
    allocate( tb_glb(ks_len,lon_len,lat_len), sb_glb(ks_len,lon_len,lat_len) )
    allocate( ub_glb(ks_len,lon_len,lat_len), vb_glb(ks_len,lon_len,lat_len) )
    
    call readnc(ncid, 'lon',   lon_t_glb)
    call readnc(ncid, 'lat',   lat_t_glb)
    call readnc(ncid, 'layer',    zz_glb)
!    print*,'lon',lon_t_glb
    
    
    dlon_glb = 0.5 * ( lon_t_glb(2) - lon_t_glb(1) )
    dlat_glb = 0.5 * ( lat_t_glb(2) - lat_t_glb(1) )
    
    do i=1,lon_len
      lon_u_glb(i) = lon_t_glb(i) - dlon_glb
    enddo
    do j=1,lat_len
      lat_v_glb(j) = lat_t_glb(j) - dlat_glb
    enddo
    
!    call get_attribute(ncid, 'scale_factor', scale_factor, 'ssh')
!    call get_attribute(ncid, 'add_offset',     add_offset, 'temp')
    scale_factor = 1.0
    add_offset = 0.0
    print*,'add-offset',add_offset,scale_factor

    call readnc(ncid2, 'Hc',   Hc_glb)
!    call readnc(ncid, 'Hu',   Hu_glb)
!    call readnc(ncid, 'Hv',   Hv_glb)
!    call readnc(ncid, 'fsm', fsm_glb)
!    call readnc(ncid, 'dum', dum_glb)
!    call readnc(ncid, 'dvm', dvm_glb)

    !---------------------------------------------------------------------------------------------
    ! Set values of H_u and H_v with H_c
    !---------------------------------------------------------------------------------------------
    do j=1,lat_len
      do i=2,lon_len
        Hu_glb(i,j) = 0.5 * (Hc_glb(i,j)+Hc_glb(i-1,j))
                                            ! Calculate H_u
      enddo
    enddo
    do j=2,lat_len
      do i=1,lon_len
        Hv_glb(i,j) = 0.5 * (Hc_glb(i,j)+Hc_glb(i,j-1))
                                            ! Calculate H_v
      enddo
    enddo
    Hv_glb(:,1) = Hv_glb(:,2)                   ! Define marginal H_v
    Hu_glb(1,:) = Hu_glb(2,:) 

    !---------------------------------------------------------------------------------------------
    ! Generate mask(fsm) with H_c
    !---------------------------------------------------------------------------------------------
    fsm_glb = 1.0
    do j=1,lat_len
      do i=1,lon_len
        if(Hc_glb(i,j)<=Land_Depth) fsm_glb(i,j)=0          ! Land point
      end do
    end do
    
    !---------------------------------------------------------------------------------------------
    ! Set velocity mask: dum&dvm according to fsm
    !---------------------------------------------------------------------------------------------
    do j=1,lat_len
      do i=2,lon_len
        dum_glb(i,j) = fsm_glb(i,j) * fsm_glb(i-1,j)
      enddo
    enddo

    do j=2,lat_len
      do i=1,lon_len
        dvm_glb(i,j) = fsm_glb(i,j) * fsm_glb(i,j-1)
      enddo
    enddo

    print*,'common init'

    call readnc(ncid, 'ssh',  ssh_glb, irec)
    call readnc(ncid, 'ua',   ua_glb,  irec)
    call readnc(ncid, 'va',   va_glb,  irec)
    call readnc(ncid, 'temp', tb_glb,  irec)
    call readnc(ncid, 'salt', sb_glb,  irec)
    call readnc(ncid, 'u',    ub_glb,  irec)
    call readnc(ncid, 'v',    vb_glb,  irec)
    
    call close_nc(ncid)
    
    print*,'lat_len',lon_len,lat_len,ks_len
    

    
    !---------------------------------------------------------------------------------------------
    ! Recover the short INT data to be original values
    !---------------------------------------------------------------------------------------------    
    do j=1,lat_len
    	do i=1,lon_len
    		ssh_glb(i,j) = ssh_glb(i,j)*scale_factor
    		 ua_glb(i,j) =  ua_glb(i,j)*scale_factor
    		 va_glb(i,j) =  va_glb(i,j)*scale_factor
    	  do k=1,ks_len
    	  	tb_glb(k,i,j) = tb_glb(k,i,j)*scale_factor + add_offset
    	  	sb_glb(k,i,j) = sb_glb(k,i,j)*scale_factor + add_offset
    	  	ub_glb(k,i,j) = ub_glb(k,i,j)*scale_factor
    	  	vb_glb(k,i,j) = vb_glb(k,i,j)*scale_factor
    	  enddo
    	enddo
    enddo
    
 !   print*,'i-1'
    
    do j=1,lat_len
    	do i=1,lon_len
        if (fsm_glb(i,j)==0) then
        	ssh_glb(i,j)  = null_value
        	tb_glb(:,i,j) = null_value
        	sb_glb(:,i,j) = null_value
        endif
        if (dum_glb(i,j)==0) then
        	ua_glb(i,j)  = null_value
        	ub_glb(:,i,j) = null_value
        endif
        if (dvm_glb(i,j)==0) then
          va_glb(i,j)  = null_value
        	vb_glb(:,i,j) = null_value
        endif
      enddo
    enddo

  !  print*,'i-2'

    !stop
     
    !---------------------------------------------------------------------------------------------
    ! Interpolate part
    !---------------------------------------------------------------------------------------------    
    allocate( lon_dt(lon_len), lat_dt(lat_len) )
    allocate( iloc_t(numboundary), iloc_u(numboundary) )
    allocate( jloc_t(numboundary), jloc_v(numboundary) )
    
   ! print*,'i-3'
    
    do n=1,numboundary
            
      !---------------------------------------------------------------------------------------------
      ! determine the index of open boundary points in the global model grids
      !---------------------------------------------------------------------------------------------    
      if(reg_xy(1,n) <= 0) then
        lon_dt  = lon_t_glb - reg_xy(1,n) - 360
      else
        lon_dt  =  lon_t_glb - reg_xy(1,n)
      endif
      iloc_t(n) =  maxloc(lon_dt, dim=1, mask=lon_dt<=0)

      if(reg_xy(3,n) <= 0) then
        lon_dt  =  lon_u_glb - reg_xy(3,n) - 360
      else
        lon_dt  =  lon_u_glb - reg_xy(3,n)
      endif
      iloc_u(n) =  maxloc(lon_dt, dim=1, mask=lon_dt<0)  		
                   
      lat_dt    =  lat_t_glb - reg_xy(2,n)
      jloc_t(n) =  maxloc(lat_dt, dim=1, mask=lat_dt<=0)
      lat_dt    =  lat_v_glb - reg_xy(6,n)
      jloc_v(n) =  maxloc(lat_dt, dim=1, mask=lat_dt<=0)
      
    enddo

    !---------------------------------------------------------------------------------------------
    ! If the interploated OBC point of finer model is within N/A area of coarser model,
    ! this originally water point is closed to be land point in this model run.
    ! Below part is provoked in the mask_define subroutine, only once throughout the model run.
    !---------------------------------------------------------------------------------------------    
    if (present(redef_boundpts)) then
    	
      if_upd_numbound = 1
    	
      do n=1,numboundary
        
        ib = ij_obc(1,n)
        jb = ij_obc(2,n)
        i1 = iloc_t(n)
        j1 = jloc_t(n)
        
        if ( fsm_glb(i1+1,j1) + fsm_glb(i1,j1)          &
           + fsm_glb(i1,j1+1) + fsm_glb(i1+1,j1+1) == 0 ) then
          fsm(ib,jb) = 0
          print*,'New land point at boundary is added in the "upd_open_bound" at ',ib,jb
          if_upd_numbound = 0
        endif
               
      enddo
      
      if (if_upd_numbound==1) then
        print*,'No new land boundary point is added in the subroutine "upd_open_bound".'
      endif
      
      deallocate( lon_t_glb, lat_t_glb, zz_glb, lon_u_glb, lat_v_glb )
      deallocate( Hc_glb, Hu_glb, Hv_glb )
      deallocate( fsm_glb, dum_glb, dvm_glb )
      deallocate( ssh_glb, ua_glb, va_glb )
      deallocate( tb_glb, sb_glb, ub_glb, vb_glb )
      deallocate( lon_dt, lat_dt, iloc_t, iloc_u, jloc_t, jloc_v )
      
      
      return
      
    endif

    !print*,'i-5'

    !---------------------------------------------------------------------------------------------
    ! interpolate 2D variables
    !---------------------------------------------------------------------------------------------    
    allocate( zm(ks_len,4) )
    
    var_char(1) = 'ETA'
    var_char(2) = 'TB'
    var_char(3) = 'SB'
    var_char(4) = 'UA'
    var_char(5) = 'Ub'
    var_char(6) = 'VA'
    var_char(7) = 'Vb'
    
!    do i=1,numboundary
!    print*,'obc',n,ij_obc(1,n),ij_obc(2,n)
!    enddo
!    stop
    
    do m=1,7
    !---------------------------------------------------------------------------------------------    
    ! Variables' sequence: [1-3]:ssh,T,S  [4-5]:ua,u  [6-7]:va,v
    !---------------------------------------------------------------------------------------------
    
      !print*,'i-6 m=',m

    
      if (m<=3) then            ! T-cell
        mask => fsm
        mask_glb => fsm_glb
        iloc => iloc_t
        jloc => jloc_t
        lon_glb => lon_t_glb
        lat_glb => lat_t_glb
        lon_reg => vlon
        lat_reg => vlat
        H_glb => Hc_glb
        H_reg => H_c
        
        if (m==1) then
          var2d => ssh_glb
          var2d_obc => et_obc
          dim_num = 2
        elseif (m==2) then
          var3d => tb_glb
          var3d_obc => tb_obc
          dim_num = 3
        elseif (m==3) then
          var3d => sb_glb
          var3d_obc => sb_obc
          dim_num = 3
        endif
        
      elseif (m<=5) then        ! U-cell
        mask => dum01
        mask_glb => dum_glb
        iloc => iloc_u
        jloc => jloc_t
        lon_glb => lon_u_glb
        lat_glb => lat_t_glb
        lon_reg => vlon_u
        lat_reg => vlat
        H_glb => Hu_glb
        H_reg => H_u
        
        if (m==4) then
          var2d => ua_glb
          var2d_obc => ua_obc
          dim_num = 2
        elseif (m==5) then
          var3d => ub_glb
          var3d_obc => ub_obc
          dim_num = 3
        endif
        
      elseif (m<=7) then        ! V-cell
        mask => dvm01
        mask_glb => dvm_glb
        iloc => iloc_t
        jloc => jloc_v
        lon_glb => lon_t_glb
        lat_glb => lat_v_glb
        lon_reg => vlon
        lat_reg => vlat_v
        H_glb => Hv_glb
        H_reg => H_v
        
        if (m==6) then
          var2d => va_glb
          var2d_obc => va_obc
          dim_num = 2
        elseif (m==7) then
          var3d => vb_glb
          var3d_obc => vb_obc
          dim_num = 3
        endif
        
      endif
      
      !---------------------------------------------------------------------------------------------    
      ! Start interpolation for both 2D and 3D variables along open boundary
      ! point by point
      !---------------------------------------------------------------------------------------------
      do n=1,numboundary
      	
       !print*,'i-7 n=',n,m,dim_num
      	
        ib = ij_obc(1,n)
        jb = ij_obc(2,n)

        if (mask(ib,jb)==0) cycle

        !---------------------------------------------------------------------------------------------
        ! Define the "squre box" which the open boundary point is surrounded
        !---------------------------------------------------------------------------------------------    
        i1 = iloc(n)
        j1 = jloc(n)
        i3 = i1 + 1
        j3 = j1 + 1
                
        x  = lon_reg(ib)
        y  = lat_reg(jb)
        x1 = lon_glb(i1)
        y1 = lat_glb(j1)
        x3 = lon_glb(i3)
        y3 = lat_glb(j3)

        if(x <= 0) x = x+360
        if(x<x1 .or. x>x3) then
          print*,'can not find the boundary info!!'
          print*,x,x1,x3
          print*,n,var_char(m)
          print*,i1,i3
          print*,iloc
          stop
        endif
        
        if (dim_num==2) then
        	
          !---------------------------------------------------------------------------------------------
          ! 2D variables : bilinear
          !---------------------------------------------------------------------------------------------    
          fin(1) = var2d(i1,j1)
          fin(2) = var2d(i1,j3)
          fin(3) = var2d(i3,j3)
          fin(4) = var2d(i3,j1)
          
          call bilinear_4points_HL(x,y,x1,x3,y1,y3,fin,fx,null_value)
          var2d_obc(n) = fx
          
          if ( fx==null_value ) then
            print*, '*No value available for interpolation for ',var_char(m)   &
                  , ', simply set to zero at ',ib,jb,i1,j1
!                  print*,'xy',x,y,x1,x3,y1,y3
 !                 print*,'fin',fin
  !                print*,'mask-glb',mask_glb(i1,j1),mask_glb(i1,j3),mask_glb(i3,j3),mask_glb(i3,j1)
  !                print*,H_glb(i1,j1),H_glb(i1,j3),H_glb(i3,j3),H_glb(i3,j1)
          	var2d_obc(n) = ZERO
          endif

        elseif (dim_num==3) then
        	          
          !---------------------------------------------------------------------------------------------
          ! 3D variable interpolation : vertical-interp1d first, then horizontal-bilinear
          !---------------------------------------------------------------------------------------------    
          do k=1,ks_len
            zm(k,1) = H_glb(i1  ,j1  ) * zz_glb(k)       !  Depth < 0
            zm(k,2) = H_glb(i1  ,j1+1) * zz_glb(k)
            zm(k,3) = H_glb(i1+1,j1+1) * zz_glb(k)
            zm(k,4) = H_glb(i1+1,j1  ) * zz_glb(k)
          enddo
          
          var3d(ks_len,:,:) = var3d(ks_len-1,:,:)
          
          do k=ksm1,1,-1
          	za = H_reg(ib,jb) * zz(k)

            call interp1d( ks_len, zm(:,1), var3d(:,i1  ,j1  ), za, fin(1), null_value )
            call interp1d( ks_len, zm(:,2), var3d(:,i1  ,j1+1), za, fin(2), null_value )
            call interp1d( ks_len, zm(:,3), var3d(:,i1+1,j1+1), za, fin(3), null_value )
            call interp1d( ks_len, zm(:,4), var3d(:,i1+1,j1  ), za, fin(4), null_value )
            
            call bilinear_4points_HL(x,y,x1,x3,y1,y3,fin,fx,null_value)
           	var3d_obc(k,n) = fx
           	
           	if (m==-1.and.ib==307.and.jb==54) then
           		print*,'===307 54:',ib,jb,n,k
           		print*,'var3d-1',var3d(:,i1,j1)
           		print*,'var3d-2',var3d(:,i1,j1+1)
           		print*,'var3d-3',var3d(:,i1+1,j1+1)
           		print*,'var3d-4',var3d(:,i1+1,j1)
           		print*,'var3d_obc',tb_obc(:,n)
           		print*,'zm-1',zm(:,1)
           		print*,'zm-2',zm(:,2)
           		print*,'zm-3',zm(:,3)
           		print*,'zm-4',zm(:,4)
           		print*,'H-reg',H_reg(ib,jb) * zz
           		!if (k==1) stop
           	endif
            
            if ( fx==null_value ) then
              if (k==ksm1) then    ! Only Uf,Vf is possible
                print*, '*No value available for interpolation for ',var_char(m)   &
                      , ', simply set to zero at ',ib,jb,k,i1,j1
                var3d_obc(k,n) = ZERO
                if ( m<=3 ) then
                	print*,'Null',var_char(m),' value interpolated at surface. Program stopped!'
                	stop
                endif
              elseif (k/=ksm1) then
            	  var3d_obc(k,n) = var3d_obc(k+1,n)
            	endif
            endif
                        
          enddo
          
           
           !if (m==3.and.ib==307.and.jb==47) then
           	!print*,'tb-obc',tb_obc(:,n)
          !endif

          
        endif
     
      enddo        ! do n=1,numboundary
      
      nullify(mask,mask_glb,iloc,jloc,lon_glb,lat_glb,H_glb,lon_reg,lat_reg,H_reg)
      
      if (dim_num==2) nullify(var2d, var2d_obc)
      if (dim_num==3) nullify(var3d, var3d_obc)
      
       
    enddo          ! do m=1,7
    
    
    !---------------------------------------------------------------------------------------------
    ! Adjust the 3D velocity to be compatible with 2D velocity with the same method used 
    ! in subroutine adjust_uv_sub
    !---------------------------------------------------------------------------------------------    
    do n=1,numboundary
    	
    	ib = ij_obc(1,n)
    	jb = ij_obc(2,n)
    	tps_u = ZERO
    	tps_v = ZERO    	
    	
    	do k=1,ksm1
    	  tps_u = tps_u + ub_obc(k,n)*dz(k)
    	  tps_v = tps_v + vb_obc(k,n)*dz(k)
    	enddo
    	do k=1,ksm1
    		ub_obc(k,n) = ( ub_obc(k,n) - tps_u + ua_obc(n) ) * dum01(ib,jb)
    		vb_obc(k,n) = ( vb_obc(k,n) - tps_v + va_obc(n) ) * dvm01(ib,jb)
      enddo    		
    
    enddo
    
    !---------------------------------------------------------------------------------------------
    ! Recheck 3D variables in case that the 3D field is interploated null.
    ! If so, replace the null_value with initial field (T&S) and motionless condition (U&V)
    !---------------------------------------------------------------------------------------------    
    m=1
    i=1
    do n=1,numboundary
  		ib = ij_obc(1,n)
  		jb = ij_obc(2,n)
  		
    	do k=1,ks
    		if (tb_obc(k,n)==null_value) print*,'*null obc T',n,ib,jb,k,tb_obc(k,n),tclim(k,ib,jb),fsm(ib,jb)
    		if (sb_obc(k,n)==null_value) print*,'*null obc S',n,ib,jb,k,sb_obc(k,n),sclim(k,ib,jb)
    		if (ub_obc(k,n)==null_value) print*,'*null obc U',n,ib,jb,k,ub_obc(k,n),zero
    		if (vb_obc(k,n)==null_value) print*,'*null obc V',n,ib,jb,k,vb_obc(k,n),zero
    		if (tb_obc(k,n)==null_value) tb_obc(k,n) = tclim(k,ib,jb)
    		if (sb_obc(k,n)==null_value) sb_obc(k,n) = sclim(k,ib,jb)
    		if (ub_obc(k,n)==null_value) ub_obc(k,n) = zero
    		if (vb_obc(k,n)==null_value) vb_obc(k,n) = zero
    	enddo
    	    	    	
    enddo    
    


    
    print*,'-----------------------------------------'
    print*,'Open boundary data updated with range:'
    print*,'max(et-obc)',maxval(et_obc),maxloc(et_obc)
    print*,'min(et-obc)',minval(et_obc),minloc(et_obc)
    print*,'max(T-obc) ',maxval(tb_obc),maxloc(tb_obc)
    print*,'min(T-obc) ',minval(tb_obc),minloc(tb_obc)
    print*,'max(S-obc) ',maxval(sb_obc),maxloc(sb_obc)
    print*,'min(S-obc) ',minval(sb_obc),minloc(sb_obc)
    print*,'max(ub-obc)',maxval(ub_obc),maxloc(ub_obc)
    print*,'min(ub-obc)',minval(ub_obc),minloc(ub_obc)
    print*,'max(vb-obc)',maxval(vb_obc),maxloc(vb_obc)
    print*,'min(vb-obc)',minval(vb_obc),minloc(vb_obc)
    print*,'max(ua-obc)',maxval(ua_obc),maxloc(ua_obc)
    print*,'min(ua-obc)',minval(ua_obc),minloc(ua_obc)
    print*,'max(va-obc)',maxval(va_obc),maxloc(va_obc)
    print*,'min(va-obc)',minval(va_obc),minloc(va_obc)
    print*,'-----------------------------------------'
    

  end subroutine upd_open_bound

!=================================================================================================

  subroutine save_data
    use netcdf_mod
    use time_mod
    
    implicit none
    
    integer :: ncid
    real(kind=precision) :: null_value
    character(len=255) :: var_name
    character(len=100) :: file_name_check
    integer :: nfield
    real, allocatable :: var_3d(:,:,:), var_2d(:,:)
    real :: time1
    character(len=100) :: ti

    null_value = 999.0
    time1 = 1.0
    
    allocate(var_3d(ks,im,jm), var_2d(im,jm))
    
    ti = datestr(days_now)

    if(iint == 1) then
      file_name_check = 'history/common_var.nc'
      call open_nc(ncid, file_name_check, 'c')
      
      call dimension_define(ncid, 'lon',   im, 'lon',   5)
      call dimension_define(ncid, 'lat',   jm, 'lat',   5)
      call dimension_define(ncid, 'layer', ks, 'layer', 5)
!      call dimension_define(ncid, 'time',   0, 'time',  5)
      
      call set_attribute(ncid, 'units', 'sigma layers',  'layer')
      call set_attribute(ncid, 'units', 'degrees_north', 'lat')
      call set_attribute(ncid, 'units', 'degrees_east',  'lon')
!      call set_attribute(ncid, 'units', 'days since 1900-1-1 00:00', 'time')
      
      call variable_define(ncid,'fsm', 5, ['lon','lat'])
      call variable_define(ncid,'dum', 5, ['lon','lat'])
      call variable_define(ncid,'dvm', 5, ['lon','lat'])
      call variable_define(ncid,'Hc', 5, ['lon','lat'])
      call variable_define(ncid,'Hu', 5, ['lon','lat'])
      call variable_define(ncid,'Hv', 5, ['lon','lat'])
      
      call variable_define(ncid, 'ssh', 5, ['lon','lat'])
      call variable_define(ncid, 'ua', 5, ['lon','lat'])
      call variable_define(ncid, 'va', 5, ['lon','lat'])
      call variable_define(ncid, 'u', 5, ['layer','lon','lat'])
      call variable_define(ncid, 'v', 5, ['layer','lon','lat'])
      call variable_define(ncid, 'temp', 5, ['layer','lon','lat'])
      call variable_define(ncid, 'salt', 5, ['layer','lon','lat'])
      call variable_define(ncid, 'rho', 5, ['layer','lon','lat'])
      call variable_define(ncid, 'km', 5, ['layer','lon','lat'])
      call variable_define(ncid, 'kh', 5, ['layer','lon','lat'])
      
      call set_attribute(ncid, 'missing_value', null_value, 'ssh')
      call set_attribute(ncid, 'missing_value', null_value, 'ua')
      call set_attribute(ncid, 'missing_value', null_value, 'ua')
      call set_attribute(ncid, 'missing_value', null_value, 'u')
      call set_attribute(ncid, 'missing_value', null_value, 'v')
      call set_attribute(ncid, 'missing_value', null_value, 'temp')
      call set_attribute(ncid, 'missing_value', null_value, 'salt')
      call set_attribute(ncid, 'missing_value', null_value, 'rho')

      call end_define(ncid)

      call writenc(ncid,'fsm', fsm)
      call writenc(ncid,'dum', dum)
      call writenc(ncid,'dvm', dvm)
      call writenc(ncid,'Hc', h_c)
      call writenc(ncid,'Hu', h_u)
      call writenc(ncid,'Hv', h_v)
      
      call writenc(ncid,'lon',vlon)
      call writenc(ncid,'lat',vlat)
      call writenc(ncid,'layer',z)

!      nfield = iint
      
!      call writenc(ncid, 'time', time1, nfield)

      var_2d = etf_c
      do i=1,im
        do j=1,jm
          if(fsm(i,j) == 0) then
            var_2d(i,j) = null_value
          endif
        enddo
      enddo
      call writenc(ncid, 'ssh', var_2d)
      var_2d = ua
      do i=1,im
        do j=1,jm
          if(dum(i,j) == 0) then
            var_2d(i,j) = null_value
          endif
        enddo
      enddo
      call writenc(ncid, 'ua', var_2d)
      var_2d = va
      do i=1,im
        do j=1,jm
          if(dvm(i,j) == 0) then
            var_2d(i,j) = null_value
          endif
        enddo
      enddo
      call writenc(ncid, 'va', var_2d)

      var_3d = uf
      do i=1,im
        do j=1,jm
          if(dum(i,j) == 0) then
            var_3d(:,i,j) = null_value
          endif
        enddo
      enddo
      call writenc(ncid, 'u', var_3d)
      var_3d = vf
      do i=1,im
        do j=1,jm
          if(dvm(i,j) == 0) then
            var_3d(:,i,j) = null_value
          endif
        enddo
      enddo
      call writenc(ncid, 'v', var_3d)
      var_3d = tf
      do i=1,im
        do j=1,jm
          if(fsm(i,j) == 0) then
            var_3d(:,i,j) = null_value
          endif
        enddo
      enddo
      call writenc(ncid, 'temp', var_3d)
      var_3d = sf
      do i=1,im
        do j=1,jm
          if(fsm(i,j) == 0) then
            var_3d(:,i,j) = null_value
          endif
        enddo
      enddo
      call writenc(ncid, 'salt', var_3d)
      var_3d = rho
      do i=1,im
        do j=1,jm
          if(fsm(i,j) == 0) then
            var_3d(:,i,j) = null_value
          endif
        enddo
      enddo
      call writenc(ncid, 'rho', var_3d)
      var_3d = km
      do i=1,im
        do j=1,jm
          if(fsm(i,j) == 0) then
            var_3d(:,i,j) = null_value
          endif
        enddo
      enddo
      call writenc(ncid, 'km', var_3d)
      var_3d = kh
      do i=1,im
        do j=1,jm
          if(fsm(i,j) == 0) then
            var_3d(:,i,j) = null_value
          endif
        enddo
      enddo
      call writenc(ncid, 'kh', var_3d)
      
    call close_nc(ncid)

    else
      file_name_check = 'history/'//ti(1:8)//'masnum.nc'
      call open_nc(ncid, file_name_check, 'c')
      
      call dimension_define(ncid, 'lon',   im, 'lon',   5)
      call dimension_define(ncid, 'lat',   jm, 'lat',   5)
      call dimension_define(ncid, 'layer', ks, 'layer', 5)
!      call dimension_define(ncid, 'time',   0, 'time',  5)
      
      call set_attribute(ncid, 'units', 'sigma layers',  'layer')
      call set_attribute(ncid, 'units', 'degrees_north', 'lat')
      call set_attribute(ncid, 'units', 'degrees_east',  'lon')
!      call set_attribute(ncid, 'units', 'days since 1900-1-1 00:00', 'time')
      
      call variable_define(ncid, 'ssh', 5, ['lon','lat'])
!      call variable_define(ncid, 'ua', 5, ['lon','lat'])
!      call variable_define(ncid, 'va', 5, ['lon','lat'])
      call variable_define(ncid, 'u', 5, ['layer','lon','lat'])
      call variable_define(ncid, 'v', 5, ['layer','lon','lat'])
      call variable_define(ncid, 'temp', 5, ['layer','lon','lat'])
      call variable_define(ncid, 'salt', 5, ['layer','lon','lat'])
      call variable_define(ncid, 'rho', 5, ['layer','lon','lat'])
!      call variable_define(ncid, 'bsmt', 5, ['layer','lon','lat'])
!      call variable_define(ncid, 'bsms', 5, ['layer','lon','lat'])
      
      call set_attribute(ncid, 'missing_value', null_value, 'ssh')
!      call set_attribute(ncid, 'missing_value', null_value, 'ua')
!      call set_attribute(ncid, 'missing_value', null_value, 'ua')
      call set_attribute(ncid, 'missing_value', null_value, 'u')
      call set_attribute(ncid, 'missing_value', null_value, 'v')
      call set_attribute(ncid, 'missing_value', null_value, 'temp')
      call set_attribute(ncid, 'missing_value', null_value, 'salt')
      call set_attribute(ncid, 'missing_value', null_value, 'rho')
!      call set_attribute(ncid, 'missing_value', null_value, 'bsmt')
!      call set_attribute(ncid, 'missing_value', null_value, 'bsms')

      call end_define(ncid)
      
      call writenc(ncid,'lon',vlon)
      call writenc(ncid,'lat',vlat)
      call writenc(ncid,'layer',z)

!      nfield = iint
      
!      call writenc(ncid, 'time', time1, nfield)

      var_2d = etf_c
      do i=1,im
        do j=1,jm
          if(fsm(i,j) == 0) then
            var_2d(i,j) = null_value
          endif
        enddo
      enddo
      call writenc(ncid, 'ssh', var_2d)

      var_3d = uf
      do i=1,im
        do j=1,jm
          if(dum(i,j) == 0) then
            var_3d(:,i,j) = null_value
          endif
        enddo
      enddo
      call writenc(ncid, 'u', var_3d)
      var_3d = vf
      do i=1,im
        do j=1,jm
          if(dvm(i,j) == 0) then
            var_3d(:,i,j) = null_value
          endif
        enddo
      enddo
      call writenc(ncid, 'v', var_3d)
      var_3d = tf
      do i=1,im
        do j=1,jm
          if(fsm(i,j) == 0) then
            var_3d(:,i,j) = null_value
          endif
        enddo
      enddo
      call writenc(ncid, 'temp', var_3d)
      var_3d = sf
      do i=1,im
        do j=1,jm
          if(fsm(i,j) == 0) then
            var_3d(:,i,j) = null_value
          endif
        enddo
      enddo
      call writenc(ncid, 'salt', var_3d)
      var_3d = rho
      do i=1,im
        do j=1,jm
          if(fsm(i,j) == 0) then
            var_3d(:,i,j) = null_value
          endif
        enddo
      enddo
      call writenc(ncid, 'rho', var_3d)

!      var_3d = bsmt
!      do i=1,im
!        do j=1,jm
!          if(fsm(i,j) == 0) then
!            var_3d(:,i,j) = null_value
!          endif
!        enddo
!      enddo
!      call writenc(ncid, 'bsmt', var_3d)
!      var_3d = bsms
!      do i=1,im
!        do j=1,jm
!          if(fsm(i,j) == 0) then
!            var_3d(:,i,j) = null_value
!          endif
!        enddo
!      enddo
!      call writenc(ncid, 'bsms', var_3d)

    call close_nc(ncid)

    endif
    
    deallocate(var_3d, var_2d)
  
    return
  end subroutine save_data
  
  !-----------------------------------------------------------------------------------------------
  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !-----------------------------------------------------------------------------------------------

  subroutine bv_mix

    if (iint==1 .or. ( ymdhms(4)/=ymdhms_prev(4) .and. mod(ymdhms(4),wav_freq)==0 ) ) then
      bv_var = ZERO
      call get_3d_field('wave_mix','bv_wtv2')
    endif
    
    km = km + bv_var
    kh = kh + 2*bv_var
  
  end subroutine bv_mix


  !-----------------------------------------------------------------------------------------------
  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !-----------------------------------------------------------------------------------------------

  subroutine bsm_mix

    if (iint==1 .or. ( ymdhms(4)/=ymdhms_prev(4) .and. mod(ymdhms(4),6)==0 ) ) then
      bsm1 = ZERO
      bsm2 = ZERO
      call get_3d_field('wave_mix','bsm1')
      call get_3d_field('wave_mix','bsm2')
    endif
  
  end subroutine bsm_mix

  !-----------------------------------------------------------------------------------------------
  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !-----------------------------------------------------------------------------------------------

  subroutine bti_mix

    real(kind=4) :: sst

    do i=1,im
      do j=1,jm
        ml_dep(i,j) = 0.0
        sst = tf(ksm1,i,j)
        do k=ks-2,1,-1
          if(sst - tf(k,i,j) >= 1) then
            ml_dep(i,j) = -zz(k) * h_c(i,j)
            exit
          endif
        enddo
        ml_dep(i,j) = min(ml_dep(i,j),500.)
      enddo
    enddo

    if (iint==1 .or. ymdhms(2)/=ymdhms_prev(2)) then
      bti_var = ZERO
      call get_3d_field('bti')

      !!!!!! validate the coefficients!
      bti_var = bti_var * 1e2
      !!!!!!

      do i=1,im
        do j=1,jm
          if(fsm(i,j) == 0) then
            bti_var(:,i,j) = 0.
            continue
          endif
!          bti_var(:,i,j) = bti_var(:,i,j) * 1e-3
        enddo
      enddo

    endif
    
    km = km + bti_var
    kh = kh + 2*bti_var
  
  end subroutine bti_mix

!=================================================================================================

subroutine write_nc(fname,a)

  use netcdf_mod

  character(len=*), intent(in) :: fname
  real(kind=precision), intent(in) :: a(:,:,:)
  integer :: ncid
  real(kind=precision), allocatable :: var(:,:,:)

  allocate(var(ks,im,jm))
    var = ZERO

  call open_nc(ncid,fname,'c')

  call dimension_define(ncid, 'lon',   im, 'lon',   5)
  call dimension_define(ncid, 'lat',   jm, 'lat',   5)
  call dimension_define(ncid, 'layer', ks, 'layer', 5)

  call variable_define(ncid,'tt', 5, ['layer','lon','lat'])

  call end_define(ncid)

  var = a
  do i=1,im
    do j=1,jm
      if(fsm(i,j)==0) var(:,i,j) = 999.
    enddo
  enddo
  call writenc(ncid,'tt',var)

  deallocate(var)

  return

end subroutine write_nc

!=================================================================================================


  subroutine save_monthly
  
    use netcdf_mod

    integer :: ncid
    character(len=100) :: ti
    real(kind=precision), allocatable :: var(:,:,:)
    real(kind=precision) :: null_value

    ti = datestr(ymdhms_prev)
    null_value = 999.0

    um = um + ub
    vm = vm + vb
    tm = tm + tb
    sm = sm + sb
    n_month = n_month + 1

    if(ymdhms(2) /= ymdhms_prev(2)) then

      um = um / n_month
      vm = vm / n_month
      tm = tm / n_month
      sm = sm / n_month
      
      allocate(var(ks,im,jm))
        var = ZERO

      call open_nc(ncid, trim(adjustl('monthly/mas_monthly_'//ti(1:6)//'.nc')),'c')

      call dimension_define(ncid, 'lon',   im, 'lon',   5)
      call dimension_define(ncid, 'lat',   jm, 'lat',   5)
      call dimension_define(ncid, 'layer', ks, 'layer', 5)

      call variable_define(ncid,'fsm', 5, ['lon','lat'])
      call variable_define(ncid,'Hc', 5, ['lon','lat'])
      
      call variable_define(ncid, 'u', 5, ['layer','lon','lat'])
      call variable_define(ncid, 'v', 5, ['layer','lon','lat'])
      call variable_define(ncid, 't', 5, ['layer','lon','lat'])
      call variable_define(ncid, 's', 5, ['layer','lon','lat'])
      
      call set_attribute(ncid, 'missing_value', null_value, 'u')
      call set_attribute(ncid, 'missing_value', null_value, 'v')
      call set_attribute(ncid, 'missing_value', null_value, 't')
      call set_attribute(ncid, 'missing_value', null_value, 's')

      call end_define(ncid)

      call writenc(ncid,'fsm', fsm)
      call writenc(ncid,'Hc', h_c)
      
      call writenc(ncid,'lon',vlon)
      call writenc(ncid,'lat',vlat)
      call writenc(ncid,'layer',z)
      
      var = um
      do i=1,im
        do j=1,jm
          if(fsm(i,j)==0) var(:,i,j) = null_value
        enddo
      enddo
      call writenc(ncid,'u',var)

      var = vm
      do i=1,im
        do j=1,jm
          if(fsm(i,j)==0) var(:,i,j) = null_value
        enddo
      enddo
      call writenc(ncid,'v',var)

      var = tm
      do i=1,im
        do j=1,jm
          if(fsm(i,j)==0) var(:,i,j) = null_value
        enddo
      enddo
      call writenc(ncid,'t',var)

      var = sm
      do i=1,im
        do j=1,jm
          if(fsm(i,j)==0) var(:,i,j) = null_value
        enddo
      enddo
      call writenc(ncid,'s',var)

      call close_nc(ncid)
      
      deallocate(var)

      um = ZERO
      vm = ZERO
      tm = ZERO
      sm = ZERO
      n_month = ZERO

    endif
  
    return

  end subroutine save_monthly

!---------------------------------------------------------------------------------------------

  subroutine write_2d(vname,var)

    use netcdf_mod

    character(len=*), intent(in) :: vname
    real(kind=precision), intent(in) :: var(im,jm)
    real, allocatable :: var_2d(:,:)

    integer :: ncid

    allocate(var_2d(im,jm))
      var_2d = var

    call open_nc(ncid, trim(adjustl(vname)),'c')
    
    call dimension_define(ncid, 'lon',   im, 'lon',   5)
    call dimension_define(ncid, 'lat',   jm, 'lat',   5)

    call variable_define(ncid,'tt', 5, ['lon','lat'])

    call end_define(ncid)

    call writenc(ncid,'lon',vlon)
    call writenc(ncid,'lat',vlat)
    call writenc(ncid,'tt',var_2d)

    call close_nc(ncid)

    deallocate(var_2d)

    return

  end subroutine write_2d


  subroutine write_3d(vname,var)

    use netcdf_mod

    character(len=*), intent(in) :: vname
    real(kind=precision), intent(in) :: var(ks,im,jm)
    real, allocatable :: var_3d(:,:,:)

    integer :: ncid

    allocate(var_3d(ks,im,jm))
      var_3d = var

    call open_nc(ncid, trim(adjustl(vname)),'c')
    
    call dimension_define(ncid, 'lon',   im, 'lon',   5)
    call dimension_define(ncid, 'lat',   jm, 'lat',   5)
    call dimension_define(ncid, 'layer', ks, 'layer', 5)

    call variable_define(ncid,'tt', 5, ['layer','lon','lat'])

    call end_define(ncid)

    call writenc(ncid,'lon',vlon)
    call writenc(ncid,'lat',vlat)
    call writenc(ncid,'layer',z)
    call writenc(ncid,'tt',var)

    call close_nc(ncid)

    deallocate(var_3d)

    return

  end subroutine write_3d

!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
! Contains over
!---------------------------------------------------------------------------------------------
  
end program masnum_circ_model_main
  

  
!=================================================================================================
!=================================================================================================
!======================********************************===========================================
!======================***    MASNUM-3D Code END    ***===========================================
!======================********************************===========================================
!=================================================================================================
