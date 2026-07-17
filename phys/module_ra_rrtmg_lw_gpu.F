!------------------------------------------------------------------
! GPU dispatch wrapper for RRTMG longwave radiation.
! M2 skeleton: mirrors the RRTMG_LWRAD interface and uses currently
! CPU implementation.
!------------------------------------------------------------------
MODULE module_ra_rrtmg_lw_gpu

CONTAINS

!------------------------------------------------------------------
   SUBROUTINE RRTMG_LWRAD_GPU(                                    &
                       rthratenlw,                                &
                       rthratenlwc,                               &
                       lwupt, lwuptc, lwuptcln, lwdnt, lwdntc, lwdntcln, &
                       lwupb, lwupbc, lwupbcln, lwdnb, lwdnbc, lwdnbcln, &
!                      lwupflx, lwupflxc, lwdnflx, lwdnflxc,      &
                       glw, olr, lwcf, emiss,                     &
                       p8w, p3d, pi3d,                            &
                       dz8w, tsk, t3d, t8w, rho3d, r, g,          &
                       icloud, warm_rain, cldfra3d,               &
                       cldovrlp,idcor,xlat,                       & 
                       lradius,iradius,                           & 
                       is_cammgmp_used,                           & 
                       f_ice_phy, f_rain_phy,                     &
                       xland, xice, snow,                         &
                       qv3d, qc3d, qr3d,                          &
                       qi3d, qs3d, qg3d,                          &
                       o3input, o33d,                             &
                       f_qv, f_qc, f_qr, f_qi, f_qs, f_qg,        &
                       re_cloud, re_ice, re_snow,                 &  ! G. Thompson
                       has_reqc, has_reqi, has_reqs,              &  ! G. Thompson
                       tauaerlw1,tauaerlw2,tauaerlw3,tauaerlw4,   & ! czhao 
                       tauaerlw5,tauaerlw6,tauaerlw7,tauaerlw8,   & ! czhao 
                       tauaerlw9,tauaerlw10,tauaerlw11,tauaerlw12,   & ! czhao 
                       tauaerlw13,tauaerlw14,tauaerlw15,tauaerlw16,   & ! czhao 
                       aer_ra_feedback,                           & !czhao
!jdfcz                 progn,prescribe,                           & !czhao
                       progn,calc_clean_atm_diag,                 & !czhao
                       qndrop3d,f_qndrop,                         & !czhao
!ccc added for time varying gases.
                       yr,julian,ghg_input,                       &
!ccc
                       mp_physics,                                &
                       ids,ide, jds,jde, kds,kde,                 & 
                       ims,ime, jms,jme, kms,kme,                 &
                       its,ite, jts,jte, kts,kte,                 &
                       lwupflx, lwupflxc, lwdnflx, lwdnflxc       &
                                                                  )
!------------------------------------------------------------------
!ccc To use clWRF time varying trace gases
   USE MODULE_RA_CLWRF_SUPPORT, ONLY : read_CAMgases

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Imports for the copied RRTMG_LWRAD body. These live at module scope in
! module_ra_rrtmg_lw; here they are pulled in at subroutine scope.
   USE module_wrf_error
#if (HWRF == 1)
   USE module_state_description, ONLY : FER_MP_HIRES, FER_MP_HIRES_ADVECT, ETAMP_HWRF
#else
   USE module_state_description, ONLY : FER_MP_HIRES, FER_MP_HIRES_ADVECT
#endif
   USE parrrtm,             ONLY : nbndlw, ngptlw
   USE rrtmg_lw_rad,        ONLY : rrtmg_lw
   USE mcica_subcol_gen_lw, ONLY : mcica_subcol_lw
! nlayers is set at init by rrtmg_lwinit; import it (do NOT redeclare) so we
! share that value. retab/deltap are the shared table/parameter.
   USE module_ra_rrtmg_lw,  ONLY : nlayers, retab, deltap, &
                                   inirad, relcalc, reicalc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   IMPLICIT NONE
!------------------------------------------------------------------
   LOGICAL, INTENT(IN )      ::        warm_rain
   LOGICAL, INTENT(IN )      ::   is_CAMMGMP_used ! Added for CAM5 RRTMG<->CAMMGMP
!
   INTEGER, INTENT(IN )      ::        ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte

   INTEGER, INTENT(IN )      ::        ICLOUD, GHG_INPUT
   INTEGER, INTENT(IN )      ::        MP_PHYSICS
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme )                 , &
         INTENT(IN   ) ::                                   dz8w, &
                                                             t3d, &
                                                             t8w, &
                                                             p8w, &
                                                             p3d, &
                                                            pi3d, &
                                                           rho3d

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme )                 , &
         INTENT(INOUT)  ::                            RTHRATENLW, &
                                                      RTHRATENLWC

   REAL, DIMENSION( ims:ime, jms:jme )                          , &
         INTENT(INOUT)  ::                                   GLW, &
                                                             OLR, &
                                                            LWCF

   REAL, DIMENSION( ims:ime, jms:jme )                          , &
         INTENT(IN   )  ::                                 EMISS, &
                                                             TSK

   REAL, INTENT(IN  )   ::                                   R,G

   REAL, DIMENSION( ims:ime, jms:jme )                          , &
         INTENT(IN   )  ::                                 XLAND, &
                                                            XICE, &
                                                            SNOW
! Added for time-varying trace gases.
   INTEGER, INTENT(IN    ) ::                                 yr
   REAL, INTENT(IN    ) ::                                julian

!
! Optional
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme )                 , &
         OPTIONAL                                               , &
         INTENT(IN   ) ::                                         &
                                                        CLDFRA3D, &
                                                         LRADIUS, &
                                                         IRADIUS, &

                                                            QV3D, &
                                                            QC3D, &
                                                            QR3D, &
                                                            QI3D, &
                                                            QS3D, &
                                                            QG3D, &
                                                        QNDROP3D

!..Added by G. Thompson to couple cloud physics effective radii.
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN)::         &
                                                        re_cloud, &
                                                          re_ice, &
                                                         re_snow
   INTEGER, INTENT(IN):: has_reqc, has_reqi, has_reqs

   real pi,third,relconst,lwpmin,rhoh2o

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme )                 , &
         OPTIONAL                                               , &
         INTENT(IN   ) ::                                         &
                                                       F_ICE_PHY, &
                                                      F_RAIN_PHY

   LOGICAL, OPTIONAL, INTENT(IN)   ::                             &
                                   F_QV,F_QC,F_QR,F_QI,F_QS,F_QG,F_QNDROP
! Optional
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), OPTIONAL ,       &
         INTENT(IN    ) :: tauaerlw1,tauaerlw2,tauaerlw3,tauaerlw4, & ! czhao 
                           tauaerlw5,tauaerlw6,tauaerlw7,tauaerlw8, & ! czhao 
                           tauaerlw9,tauaerlw10,tauaerlw11,tauaerlw12, & ! czhao 
                           tauaerlw13,tauaerlw14,tauaerlw15,tauaerlw16

   INTEGER,    INTENT(IN  ), OPTIONAL   ::       aer_ra_feedback
!jdfcz   INTEGER,    INTENT(IN  ), OPTIONAL   ::       progn,prescribe
   INTEGER,    INTENT(IN  ), OPTIONAL   ::       progn
   INTEGER,    INTENT(IN  )             ::       calc_clean_atm_diag

!  Ozone
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme )                 , &
         INTENT(INOUT) :: O33D
   INTEGER, INTENT(IN ) :: o3input

      real, parameter :: thresh=1.e-9
      real slope
      character(len=200) :: msg


! Top of atmosphere and surface longwave fluxes (W m-2)
   REAL, DIMENSION( ims:ime, jms:jme ),                           &
         OPTIONAL, INTENT(INOUT) ::                               &
                      LWUPT,LWUPTC,LWUPTCLN,LWDNT,LWDNTC,LWDNTCLN,&
                      LWUPB,LWUPBC,LWUPBCLN,LWDNB,LWDNBC,LWDNBCLN 

! Layer longwave fluxes (including extra layer above model top)
! Vertical ordering is from bottom to top (W m-2)
   REAL, DIMENSION( ims:ime, kms:kme+2, jms:jme ),                &
         OPTIONAL, INTENT(OUT) ::                                 &
                               LWUPFLX,LWUPFLXC,LWDNFLX,LWDNFLXC

!  LOCAL VARS
 
   REAL, DIMENSION( kts:kte+1 ) ::                          Pw1D, &
                                                            Tw1D

   REAL, DIMENSION( kts:kte ) ::                          TTEN1D, &
                                                        CLDFRA1D, &
                                                            DZ1D, &
                                                             P1D, &
                                                             T1D, &
                                                            QV1D, &
                                                            QC1D, &
                                                            QR1D, &
                                                            QI1D, &
                                                           RHO1D, &
                                                            QS1D, &
                                                            QG1D, &
                                                            O31D, &
                                                          qndrop1d 
!BSF: From eq. (5) on p. 2434 in McFarquhar & Heymsfield (1996)
       real, parameter :: re_50C=1250.0/9.917, re_40C=1250.0/9.337, &
                          re_30C=1250.0/9.208, re_20C=1250.0/9.387


! Added local arrays for RRTMG
    integer ::                                              ncol, &
                                                            nlay, &
                                                            icld, &
                                                        cldovrlp, &
                                                           idcor, &
                                                          juldat, &
                                                         inflglw, &
                                                        iceflglw, &
                                                        liqflglw
! Latitude 
    real, dimension( ims:ime,jms:jme ), intent(in) ::       xlat     ! (latitude for cldovrlp=4 or 5)
    real                                           ::        lat

! Dimension with extra layer from model top to TOA
    real, dimension( 1, kts:nlayers+1 )  ::                 plev, &
                                                            tlev
    real, dimension( 1, kts:nlayers )  ::                   play, &
                                                            tlay, &
                                                          h2ovmr, &
                                                           o3vmr, &
                                                          co2vmr, &
                                                           o2vmr, &
                                                          ch4vmr, &
                                                          n2ovmr, &
                                                        cfc11vmr, &
                                                        cfc12vmr, &
                                                        cfc22vmr, &
                                                         ccl4vmr
    real, dimension( kts:nlayers )  ::                     o3mmr
! Add height of each layer for exponential-random cloud overlap
! This will be derived below from the dz in each layer
    real, dimension( 1, kts:nlayers )  ::                   hgt
    real ::                                               dzsum
! For old cloud property specification for rrtm_lw
    real, dimension( kts:kte )  ::                          clwp, &
                                                            ciwp, &
                                                            cswp, &
                                                            plwp, &
                                                            piwp
! Surface emissivity (for 16 LW spectral bands)
    real, dimension( 1, nbndlw )  ::                        emis
! Dimension with extra layer from model top to TOA, 
! though no clouds are allowed in extra layer
    real, dimension( 1, kts:nlayers )  ::                 clwpth, &
                                                          ciwpth, &
                                                          cswpth, &
                                                             rel, &
                                                             rei, &
                                                             res, &
                                                         cldfrac, &
                                                         relqmcl, &
                                                         reicmcl, &
                                                         resnmcl
    real, dimension( nbndlw, 1, kts:nlayers )  ::        taucld
    real, dimension( ngptlw, 1, kts:nlayers )  ::        cldfmcl, &
                                                         clwpmcl, &
                                                         ciwpmcl, &
                                                         cswpmcl, &
                                                         taucmcl
    real, dimension( 1, kts:nlayers, nbndlw )  ::           tauaer

! Output arrays contain extra layer from model top to TOA
    real, dimension( 1, kts:nlayers+1 )  ::                 uflx, &
                                                            dflx, &
                                                           uflxc, &
                                                           dflxc, &
                                                         uflxcln, &
                                                         dflxcln
                                                           
    real, dimension( 1, kts:nlayers )  ::                    hr, &
                                                             hrc

    real, dimension ( 1 ) ::                                tsfc, &
                                                              ps
    real ::                                                   ro, &
                                                              dz
    real:: snow_mass_factor

!..We can use message interface regardless of what options are running,
!.. so let us ask for it here.
      CHARACTER(LEN=256)                           :: message
      LOGICAL, EXTERNAL                            :: wrf_dm_on_monitor

!ccc To add time-varying trace gases (CO2, N2O and CH4). Read the conc.  from file
! then interpolate to date of run.
      REAL(8)                                      :: co2, n2o, ch4, cfc11, cfc12

! Set trace gas volume mixing ratios, 2005 values, IPCC (2007)
! cfc-22 (169 ppt)
    real :: cfc22
    data cfc22 / 0.169e-9 / 
! ccl4 (93 ppt)
    real :: ccl4
    data ccl4 / 0.093e-9 / 
! Set oxygen volume mixing ratio (for o2mmr=0.23143)
    real :: o2
    data o2 / 0.209488 /

    integer :: iplon, irng, permuteseed
    integer :: nb

! For old cloud property specification for rrtm_lw
! Cloud and precipitation absorption coefficients
    real :: abcw,abice,abrn,absn
    data abcw /0.144/
    data abice /0.0735/
    data abrn /0.330e-3/
    data absn /2.34e-3/

! Molecular weights and ratios for converting mmr to vmr units
!    real :: amd       ! Effective molecular weight of dry air (g/mol)  
!    real :: amw       ! Molecular weight of water vapor (g/mol)        
!    real :: amo       ! Molecular weight of ozone (g/mol)              
!    real :: amo2      ! Molecular weight of oxygen (g/mol)              
! Atomic weights for conversion from mass to volume mixing ratios                
!    data amd   /  28.9660   /                                                  
!    data amw   /  18.0160   /                                                  
!    data amo   /  47.9998   /                                                  
!    data amo2  /  31.9999   /
                                                                                 
    real :: amdw     ! Molecular weight of dry air / water vapor  
    real :: amdo     ! Molecular weight of dry air / ozone
    real :: amdo2    ! Molecular weight of dry air / oxygen
    data amdw /  1.607793 /                                                    
    data amdo /  0.603461 /
    data amdo2 / 0.905190 /
    
!!
    real, dimension( 1, 1:kte-kts+1 )  :: pdel         ! Layer pressure thickness (mb)

    real, dimension(1, 1:kte-kts+1) ::   cicewp, &     ! in-cloud cloud ice water path
                                         cliqwp, &     ! in-cloud cloud liquid water path
                                         csnowp, &     ! in-cloud snow water path
                                          reliq, &     ! effective drop radius (microns)
                                          reice        ! effective ice crystal size (microns)
    real, dimension(1, 1:kte-kts+1):: recloud1d, &
                                        reice1d, &
                                       resnow1d

    real :: gliqwp, gicewp, gsnowp, gravmks

!
!    REAL   ::  TSFC,GLW0,OLR0,EMISS0,FP

    real, dimension (1) :: landfrac, landm, snowh, icefrac

    integer :: pcols, pver

!
    INTEGER :: i,j,K, idx_rei
    REAL :: corr
    LOGICAL :: predicate

! Added for top of model adjustment.  Steven Cavallo NCAR/MMM December 2010
    INTEGER, PARAMETER :: nproflevs = 60 ! Constant, from the table
    INTEGER :: L, LL, klev               ! Loop indices      
    REAL, DIMENSION( kts:nlayers+1 ) :: varint
    REAL :: wght,vark,vark1,tem1,tem2,tem3
    REAL :: PPROF(nproflevs), TPROF(nproflevs)            
    ! Weighted mean pressure and temperature profiles from midlatitude 
    ! summer (MLS),midlatitude winter (MLW), sub-Arctic 
    ! winter (SAW),sub-Arctic summer (SAS), and tropical (TROP) 
    ! standard atmospheres.
    DATA PPROF   /1000.00,855.47,731.82,626.05,535.57,458.16,     &
                  391.94,335.29,286.83,245.38,209.91,179.57,      &
                  153.62,131.41,112.42,96.17,82.27,70.38,         &
                  60.21,51.51,44.06,37.69,32.25,27.59,            &
                  23.60,20.19,17.27,14.77,12.64,10.81,            &
                  9.25,7.91,6.77,5.79,4.95,4.24,                  &
                  3.63,3.10,2.65,2.27,1.94,1.66,                  &
                  1.42,1.22,1.04,0.89,0.76,0.65,                  &
                  0.56,0.48,0.41,0.35,0.30,0.26,                  &
                  0.22,0.19,0.16,0.14,0.12,0.10/
    DATA TPROF   /286.96,281.07,275.16,268.11,260.56,253.02,      &
                  245.62,238.41,231.57,225.91,221.72,217.79,      &
                  215.06,212.74,210.25,210.16,210.69,212.14,      &
                  213.74,215.37,216.82,217.94,219.03,220.18,      &
                  221.37,222.64,224.16,225.88,227.63,229.51,      &
                  231.50,233.73,236.18,238.78,241.60,244.44,      &
                  247.35,250.33,253.32,256.30,259.22,262.12,      &
                  264.80,266.50,267.59,268.44,268.69,267.76,      &
                  266.13,263.96,261.54,258.93,256.15,253.23,      &
                  249.89,246.67,243.48,240.25,236.66,233.86/    

! gpu: kernel launcher
   INTERFACE
      SUBROUTINE wrf_gpu_lw(x, n) BIND(C, name="wrf_gpu_lw")
         USE iso_c_binding, ONLY : c_float, c_int
         REAL(c_float)    :: x(*)
         INTEGER(c_int)   :: n
      END SUBROUTINE
   END INTERFACE

! gpu: batching parameters
   INTEGER, PARAMETER :: batch_size = 512  ! max columns per solver batch (= GPU launch width)
   INTEGER :: nx, ny, ncols, batch_start, nbatch, bcol
   INTEGER :: i_of(batch_size), j_of(batch_size)      ! bcol -> (i,j) mapping for batching
! gpu: batched rrtmg_lw inputs
   REAL, ALLOCATABLE, DIMENSION(:,:)   :: play_p, tlay_p,                       & ! (batch, nlay)
                                          plev_p, tlev_p,                       & ! (batch, nlay+1)
                                          h2ovmr_p, o3vmr_p, co2vmr_p,          & ! (batch, nlay)
                                          o2vmr_p, ch4vmr_p, n2ovmr_p,          &
                                          cfc11vmr_p, cfc12vmr_p,               &
                                          cfc22vmr_p, ccl4vmr_p,                &
                                          reicmcl_p, relqmcl_p, resnmcl_p,      & ! (batch, nlay)
                                          emis_p                                  ! (batch, nbndlw)
   REAL, ALLOCATABLE, DIMENSION(:)     :: tsfc_p                                  ! (batch)
   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: cldfmcl_p, taucmcl_p, ciwpmcl_p,      & ! (ngptlw, batch, nlay) - ncol is the MIDDLE dim
                                          clwpmcl_p, cswpmcl_p,                 &
                                          tauaer_p                                ! (batch, nlay, nbndlw)
! gpu: batched rrtmg_lw outputs
   REAL, ALLOCATABLE, DIMENSION(:,:)   :: uflx_p, dflx_p, uflxc_p, dflxc_p,     & ! (batch, nlay+1)
                                          uflxcln_p, dflxcln_p,                 &
                                          hr_p, hrc_p                             ! (batch, nlay)
!------------------------------------------------------------------
#if ( WRF_CHEM == 1 )
      IF ( aer_ra_feedback == 1) then
      IF ( .NOT. &
      ( PRESENT(tauaerlw1) .AND. &
        PRESENT(tauaerlw2) .AND. &
        PRESENT(tauaerlw3) .AND. &
        PRESENT(tauaerlw4) .AND. &
        PRESENT(tauaerlw5) .AND. &
        PRESENT(tauaerlw6) .AND. &
        PRESENT(tauaerlw7) .AND. &
        PRESENT(tauaerlw8) .AND. &
        PRESENT(tauaerlw9) .AND. &
        PRESENT(tauaerlw10) .AND. &
        PRESENT(tauaerlw11) .AND. &
        PRESENT(tauaerlw12) .AND. &
        PRESENT(tauaerlw13) .AND. &
        PRESENT(tauaerlw14) .AND. &
        PRESENT(tauaerlw15) .AND. &
        PRESENT(tauaerlw16) ) ) THEN
      CALL wrf_error_fatal  &
      ('Warning: missing fields required for aerosol radiation' )
      ENDIF
      ENDIF
#endif


!-----CALCULATE LONG WAVE RADIATION
!                                                              
! All fields are ordered vertically from bottom to top
! Pressures are in mb
!

! Read time-varying trace gases concentrations and interpolate them to run date.
!
   IF ( GHG_INPUT .EQ. 1 ) THEN 
      CALL read_CAMgases(yr,julian,.false.,"RRTMG",co2,n2o,ch4,cfc11,cfc12)
      IF ( wrf_dm_on_monitor() ) THEN
        WRITE(message,*)'RRTMG LW CLWRF interpolated GHG values year:',yr,' julian day:',julian
        call wrf_debug( 1, message)
        WRITE(message,*)'  co2vmr: ',co2,' n2ovmr:',n2o,' ch4vmr:',ch4,' cfc11vmr:',cfc11,' cfc12vmr:',cfc12
        call wrf_debug( 1, message)
      ENDIF
   ELSE
! Set trace gas volume mixing ratios, 2005 values, IPCC (2007)
! Annual function for co2 in WRF v4.2
      co2 = (280. + 90.*exp(0.02*(yr-2000)))*1.e-6
!     co2 = 379.e-6
      ch4 = 1774.e-9
      n2o = 319.e-9
      cfc11 = 0.251e-9
      cfc12 = 0.538e-9
   END IF

! gpu: allocate the solver batch (on the heap)
   nx = ite - its + 1
   ny = jte - jts + 1
   ALLOCATE( play_p (batch_size, kts:nlayers),                                 &
             tlay_p (batch_size, kts:nlayers),                                 &
             plev_p (batch_size, kts:nlayers+1),                               &
             tlev_p (batch_size, kts:nlayers+1),                               &
             tsfc_p (batch_size),                                              &
             h2ovmr_p  (batch_size, kts:nlayers),                              &
             o3vmr_p   (batch_size, kts:nlayers),                              &
             co2vmr_p  (batch_size, kts:nlayers),                              &
             o2vmr_p   (batch_size, kts:nlayers),                              &
             ch4vmr_p  (batch_size, kts:nlayers),                              &
             n2ovmr_p  (batch_size, kts:nlayers),                              &
             cfc11vmr_p(batch_size, kts:nlayers),                              &
             cfc12vmr_p(batch_size, kts:nlayers),                              &
             cfc22vmr_p(batch_size, kts:nlayers),                              &
             ccl4vmr_p (batch_size, kts:nlayers),                              &
             emis_p    (batch_size, nbndlw),                                   &
             cldfmcl_p (ngptlw, batch_size, kts:nlayers),                      &
             taucmcl_p (ngptlw, batch_size, kts:nlayers),                      &
             ciwpmcl_p (ngptlw, batch_size, kts:nlayers),                      &
             clwpmcl_p (ngptlw, batch_size, kts:nlayers),                      &
             cswpmcl_p (ngptlw, batch_size, kts:nlayers),                      &
             reicmcl_p (batch_size, kts:nlayers),                              &
             relqmcl_p (batch_size, kts:nlayers),                              &
             resnmcl_p (batch_size, kts:nlayers),                              &
             tauaer_p  (batch_size, kts:nlayers, nbndlw) )
   ALLOCATE( uflx_p   (batch_size, kts:nlayers+1),                             &
             dflx_p   (batch_size, kts:nlayers+1),                             &
             uflxc_p  (batch_size, kts:nlayers+1),                             &
             dflxc_p  (batch_size, kts:nlayers+1),                             &
             uflxcln_p(batch_size, kts:nlayers+1),                             &
             dflxcln_p(batch_size, kts:nlayers+1),                             &
             hr_p     (batch_size, kts:nlayers),                               &
             hrc_p    (batch_size, kts:nlayers) )

! gpu: process tile in batches of columns
   ncols = nx * ny
   batch_loop: do batch_start = 0, ncols-1, batch_size
      nbatch = min( batch_size, ncols - batch_start )
! gpu: pack buffers
      pack_loop: do bcol = 1, nbatch

         i = its + mod( batch_start + bcol - 1, nx )
         j = jts + ( batch_start + bcol - 1 ) / nx
         i_of(bcol) = i
         j_of(bcol) = j

         do k=kts,kte+1
            Pw1D(K) = p8w(I,K,J)/100.
            Tw1D(K) = t8w(I,K,J)
         enddo

         DO K=kts,kte
            QV1D(K)=0.
            QC1D(K)=0.
            QR1D(K)=0.
            QI1D(K)=0.
            QS1D(K)=0.
            CLDFRA1D(k)=0.
         ENDDO

         DO K=kts,kte
            QV1D(K)=QV3D(I,K,J)
            QV1D(K)=max(0.,QV1D(K))
         ENDDO

         IF (o3input.eq.2) THEN
            DO K=kts,kte
               O31D(K)=O33D(I,K,J)
            ENDDO
         ELSE
            DO K=kts,kte
               O31D(K)=0.0
            ENDDO
         ENDIF

         DO K=kts,kte
            TTEN1D(K)=0.
            T1D(K)=T3D(I,K,J)
            P1D(K)=P3D(I,K,J)/100.
            RHO1D(K)=RHO3D(I,K,J)
            DZ1D(K)=dz8w(I,K,J)
         ENDDO

! moist variables

         IF (ICLOUD .ne. 0) THEN
            IF ( PRESENT( CLDFRA3D ) ) THEN
              DO K=kts,kte
                 CLDFRA1D(k)=CLDFRA3D(I,K,J)
              ENDDO
            ENDIF

            IF (PRESENT(F_QC) .AND. PRESENT(QC3D)) THEN
              IF ( F_QC) THEN
                 DO K=kts,kte
                    QC1D(K)=QC3D(I,K,J)
                    QC1D(K)=max(0.,QC1D(K))
                 ENDDO
              ENDIF
            ENDIF

            IF (PRESENT(F_QR) .AND. PRESENT(QR3D)) THEN
              IF ( F_QR) THEN
                 DO K=kts,kte
                    QR1D(K)=QR3D(I,K,J)
                    QR1D(K)=max(0.,QR1D(K))
                 ENDDO
              ENDIF
            ENDIF

            IF ( PRESENT(F_QNDROP).AND.PRESENT(QNDROP3D)) THEN
             IF (F_QNDROP) THEN
              DO K=kts,kte
               qndrop1d(K)=qndrop3d(I,K,J)
              ENDDO
             ENDIF
            ENDIF

! This logic is tortured because cannot test F_QI unless
! it is present, and order of evaluation of expressions
! is not specified in Fortran

            IF ( PRESENT ( F_QI ) ) THEN
              predicate = F_QI
            ELSE
              predicate = .FALSE.
            ENDIF

! For MP option 3
            IF (.NOT. predicate .and. .not. warm_rain) THEN
               DO K=kts,kte
                  IF (T1D(K) .lt. 273.15) THEN
                  QI1D(K)=QC1D(K)
                  QS1D(K)=QR1D(K)
                  QC1D(K)=0.
                  QR1D(K)=0.
                  ENDIF
               ENDDO
            ENDIF

            IF (PRESENT(F_QI) .AND. PRESENT(QI3D)) THEN
               IF (F_QI) THEN
                  DO K=kts,kte
                     QI1D(K)=QI3D(I,K,J)
                     QI1D(K)=max(0.,QI1D(K))
                  ENDDO
               ENDIF
            ENDIF

            IF (PRESENT(F_QS) .AND. PRESENT(QS3D)) THEN
               IF (F_QS) THEN
                  DO K=kts,kte
                     QS1D(K)=QS3D(I,K,J)
                     QS1D(K)=max(0.,QS1D(K))
                  ENDDO
               ENDIF
            ENDIF

            IF (PRESENT(F_QG) .AND. PRESENT(QG3D)) THEN
               IF (F_QG) THEN
                  DO K=kts,kte
                     QG1D(K)=QG3D(I,K,J)
                     QG1D(K)=max(0.,QG1D(K))
                  ENDDO
               ENDIF
            ENDIF

! mji - For MP option 5
            IF ( PRESENT(F_QI) .and. PRESENT(F_QC) .and. PRESENT(F_QS) .and. PRESENT(F_ICE_PHY) ) THEN
               IF ( F_QC .and. .not. F_QI .and. F_QS ) THEN
                  DO K=kts,kte
                     qi1d(k) = 0.1*qs3d(i,k,j)
                     qs1d(k) = 0.9*qs3d(i,k,j)
                     qc1d(k) = qc3d(i,k,j)
                     qi1d(k) = max(0.,qi1d(k))
                     qc1d(k) = max(0.,qc1d(k))
                  ENDDO
               ENDIF
            ENDIF

        ENDIF

!   For mp option=5 or 85  (new Ferrier- Aligo or fer_hires scheme), QI3D saves all
#if (HWRF == 1)
        IF ( mp_physics == FER_MP_HIRES .OR. &
             mp_physics == FER_MP_HIRES_ADVECT .OR. &
             mp_physics == ETAMP_HWRF ) THEN
#else
        IF ( mp_physics == FER_MP_HIRES .OR. &
             mp_physics == FER_MP_HIRES_ADVECT) THEN
#endif
                  DO K=kts,kte
                     qi1d(k) = qi3d(i,k,j)
                     qs1d(k) = 0.0
                     qc1d(k) = qc3d(i,k,j)
                     qi1d(k) = max(0.,qi1d(k))
                     qc1d(k) = max(0.,qc1d(k))
                  ENDDO
        ENDIF

!         EMISS0=EMISS(I,J)
!         GLW0=0. 
!         OLR0=0. 
!         TSFC=TSK(I,J)
         DO K=kts,kte
            QV1D(K)=AMAX1(QV1D(K),1.E-12) 
         ENDDO

! Set up input for longwave
         ncol = 1
! Add extra layer from top of model to top of atmosphere
!         nlay = (kte - kts + 1) + 1
! Edited for top of model adjustment (nlayers = kte + 1).  
! Steven Cavallo, December 2010
          nlay = nlayers ! Keep these indices the same

! Select cloud overlap assumption (1 = random, 2 = maximum-random, 3 = maximum, 4 = exponential, 5 = exponential-random
         icld=cldovrlp ! J. Henderson AER assign namelist variable cldovrlp to existing icld

! Set julian date
         juldat = julian
! Select cloud liquid and ice optics parameterization options
! For passing in cloud optical properties directly:
!         icld = 2
!         inflglw = 0
!         iceflglw = 0
!         liqflglw = 0
! For passing in cloud physical properties; cloud optics parameterized in RRTMG:
         inflglw = 2
         iceflglw = 3
         liqflglw = 1

!Mukul change the flags here with reference to the new effective cloud/ice/snow radius
         IF (ICLOUD .ne. 0) THEN
            IF ( has_reqc .ne. 0) THEN
               inflglw = 3
               DO K=kts,kte
                  recloud1D(ncol,K) = MAX(2.5, re_cloud(I,K,J)*1.E6)
                  if (recloud1D(ncol,K).LE.2.5.AND.cldfra3d(i,k,j).gt.0. &
     &                            .AND. (XLAND(I,J)-1.5).GT.0.) then     !--- Ocean
                     recloud1D(ncol,K) = 10.5
                  elseif(recloud1D(ncol,K).LE.2.5.AND.cldfra3d(i,k,j).gt.0. &
     &                            .AND. (XLAND(I,J)-1.5).LT.0.) then     !--- Land
                     recloud1D(ncol,K) = 7.5
                  endif
               ENDDO
            ELSE
               DO K=kts,kte
                  recloud1D(ncol,K) = 5.0
               ENDDO
            ENDIF

            IF ( has_reqi .ne. 0) THEN
               inflglw  = 4
               iceflglw = 4
               DO K=kts,kte
                  reice1D(ncol,K) = MAX(5., re_ice(I,K,J)*1.E6)
                  if (reice1D(ncol,K).LE.5..AND.cldfra3d(i,k,j).gt.0.) then
                     idx_rei = int(t3d(i,k,j)-179.)
                     idx_rei = min(max(idx_rei,1),75)
                     corr = t3d(i,k,j) - int(t3d(i,k,j))
                     reice1D(ncol,K) = retab(idx_rei)*(1.-corr) +       &
     &                                 retab(idx_rei+1)*corr
                     reice1D(ncol,K) = MAX(reice1D(ncol,K), 5.0)
                  endif
               ENDDO
            ELSE
               DO K=kts,kte
#if (EM_CORE==1) 
                  reice1D(ncol,K) = 10.0
#else
                  tem2 = 25.0  !- was 10.0
                  tem3=1.e3*rho1d(k)*qi1d(k)  !- IWC (g m^-3)
                  if (tem3>thresh) then       !- Only when IWC>1.e-9 gm^-3
                    tem1=t1d(k)-273.15
                    if (tem1 < -50.0) then
                      tem2 = re_50C*tem3**0.109
                    elseif (tem1 < -40.0) then
                      tem2 = re_40C*tem3**0.08
                    elseif (tem1 < -30.0) then
                      tem2 = re_30C*tem3**0.055
                    else
                      tem2 = re_20C*tem3**0.031
                    endif
                    tem2 = max(25.,tem2)
                  endif
                  reice1D(ncol,K) = min(tem2, 135.72)   !- 1.0315*reice<= 140 microns
#endif
               ENDDO
            ENDIF

            IF ( has_reqs .ne. 0) THEN
               inflglw  = 5
               iceflglw = 5
               DO K=kts,kte
                  resnow1D(ncol,K) = MAX(10., re_snow(I,K,J)*1.E6)
               ENDDO
            ELSE
               DO K=kts,kte
                  resnow1D(ncol,K) = 10.0
               ENDDO
            ENDIF

! special case for P3 microphysics
! put ice into snow category for optics, then set ice to zero
            IF (has_reqs .eq. 0 .and. has_reqi .ne. 0 .and. has_reqc .ne. 0) THEN
               inflglw  = 5
               iceflglw = 5
               DO K=kts,kte
                  resnow1D(ncol,K) = MAX(10., re_ice(I,K,J)*1.E6)
                  QS1D(K)=QI3D(I,K,J)
                  QI1D(K)=0.
                  reice1D(ncol,K)=10.
               END DO
            END IF

         ENDIF

! Layer indexing goes bottom to top here for all fields.
! Water vapor and ozone are converted from mmr to vmr. 
! Pressures are in units of mb here. 
         plev(ncol,1) = pw1d(1)
         tlev(ncol,1) = tw1d(1)
         tsfc(ncol) = tsk(i,j)
         do k = kts, kte
            play(ncol,k) = p1d(k)
            plev(ncol,k+1) = pw1d(k+1)
            pdel(ncol,k) = plev(ncol,k) - plev(ncol,k+1)
            tlay(ncol,k) = t1d(k)
            tlev(ncol,k+1) = tw1d(k+1)
            h2ovmr(ncol,k) = qv1d(k) * amdw
            co2vmr(ncol,k) = co2
            o2vmr(ncol,k) = o2
            ch4vmr(ncol,k) = ch4
            n2ovmr(ncol,k) = n2o
            cfc11vmr(ncol,k) = cfc11
            cfc12vmr(ncol,k) = cfc12
            cfc22vmr(ncol,k) = cfc22
            ccl4vmr(ncol,k) = ccl4
         enddo

! Derive height of each layer mid-point from layer thickness.
! Needed for exponential (icld=4) and exponential-random overlap (icld=5) options only.
         dzsum = 0.0
         do k = kts, kte
            dz = dz1d(k)
            hgt(ncol,k) = dzsum + 0.5*dz
            dzsum = dzsum + dz
         enddo

! This section is replaced with a new method to deal with model top
         if ( 1 == 0 ) then

!  Define profile values for extra layer from model top to top of atmosphere. 
!  The top layer temperature for all gridpoints is set to the top layer-1 
!  temperature plus a constant (0 K) that represents an isothermal layer    
!  above ptop.  Top layer interface temperatures are linearly interpolated 
!  from the layer temperatures.  

         play(ncol,kte+1) = 0.5 * plev(ncol,kte+1)
         tlay(ncol,kte+1) = tlev(ncol,kte+1) + 0.0
         plev(ncol,kte+2) = 1.0e-5
         tlev(ncol,kte+2) = tlev(ncol,kte+1) + 0.0
         h2ovmr(ncol,kte+1) = h2ovmr(ncol,kte) 
         co2vmr(ncol,kte+1) = co2vmr(ncol,kte) 
         o2vmr(ncol,kte+1) = o2vmr(ncol,kte) 
         ch4vmr(ncol,kte+1) = ch4vmr(ncol,kte) 
         n2ovmr(ncol,kte+1) = n2ovmr(ncol,kte) 
         cfc11vmr(ncol,kte+1) = cfc11vmr(ncol,kte) 
         cfc12vmr(ncol,kte+1) = cfc12vmr(ncol,kte) 
         cfc22vmr(ncol,kte+1) = cfc22vmr(ncol,kte) 
         ccl4vmr(ncol,kte+1) = ccl4vmr(ncol,kte) 

         endif

!  Set up values for extra layers to the top of the atmosphere.                       
!  Temperature is calculated based on an average temperature profile given
!  here in a table.  The input table data is linearly interpolated to the
!  column pressure.  Mixing ratios are held constant except for ozone.  
!  Caution should be used if model top pressure is less than 5 hPa.
!  Steven Cavallo, NCAR/MMM, December 2010
       ! Calculate the column pressure buffer levels above the 
       ! model top       
       do L=kte+1,nlayers,1
          plev(ncol,L+1) = plev(ncol,L) - deltap
          play(ncol,L) = 0.5*(plev(ncol,L) + plev(ncol,L+1))
! Fill in height array above model top to top of atmosphere using
! dz from model top layer for completeness, though this information is not
! likely to be used by the exponential-random cloud overlap method.
          hgt(ncol,L) = dzsum + 0.5*dz
          dzsum = dzsum + dz
       enddo          
       ! Add zero as top level.  This gets the temperature max at the
       ! stratopause, reducing the downward flux errors in the top 
       ! levels.  If zero happened to be the top level already,
       ! this will add another level with zero, but will not affect
       ! the radiative transfer calculation.
       plev(ncol,nlayers+1) = 0.00
       play(ncol,nlayers) =  0.5*(plev(ncol,nlayers) + plev(ncol,nlayers+1))

       ! Interpolate the table temperatures to column pressure levels    
       do L=1,nlayers+1,1
          if ( PPROF(nproflevs) .lt. plev(ncol,L) ) then
             do LL=2,nproflevs,1       
                if ( PPROF(LL) .lt. plev(ncol,L) ) then           
                   klev = LL - 1
                   exit
                endif
             enddo
          
          else
             klev = nproflevs
          endif  
  
          if (klev .ne. nproflevs ) then
             vark  = TPROF(klev) 
             vark1 = TPROF(klev+1)
             wght=(plev(ncol,L)-PPROF(klev) )/( PPROF(klev+1)-PPROF(klev))
          else
             vark  = TPROF(klev) 
             vark1 = TPROF(klev)
             wght = 0.0
          endif
          varint(L) = wght*(vark1-vark)+vark

       enddo                   
       
       ! Match the interpolated table temperature profile to WRF column                    
       do L=kte+1,nlayers+1,1
          tlev(ncol,L) = varint(L) + (tlev(ncol,kte) - varint(kte))
          !if ( L .le. nlay ) then
          tlay(ncol,L-1) = 0.5*(tlev(ncol,L) + tlev(ncol,L-1))  
          !endif
       enddo 

       ! Now the chemical species (except for ozone)
       do L=kte+1,nlayers,1
          h2ovmr(ncol,L) = h2ovmr(ncol,kte) 
          co2vmr(ncol,L) = co2vmr(ncol,kte) 
          o2vmr(ncol,L) = o2vmr(ncol,kte) 
          ch4vmr(ncol,L) = ch4vmr(ncol,kte) 
          n2ovmr(ncol,L) = n2ovmr(ncol,kte) 
          cfc11vmr(ncol,L) = cfc11vmr(ncol,kte) 
          cfc12vmr(ncol,L) = cfc12vmr(ncol,kte) 
          cfc22vmr(ncol,L) = cfc22vmr(ncol,kte) 
          ccl4vmr(ncol,L) = ccl4vmr(ncol,kte) 
       enddo     
! End top of model buffer 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Get ozone profile including amount in extra layer above model top.
! Steven Cavallo: Must pass nlay-1 into subroutine to get nlayers 
! dimension for o3mmr
         call inirad (o3mmr,plev,kts,nlay-1)

! Steven Cavallo: Changed to nlayers from kte+1
        if(o3input.eq.2) then
         do k = kts, nlayers
            o3vmr(ncol,k) = o3mmr(k) * amdo
            if(k.le.kte)then
               o3vmr(ncol,k) = o31d(k)
            else
! apply shifted climatology profile above model top
               o3vmr(ncol,k) = o31d(kte) - o3mmr(kte)*amdo + o3mmr(k)*amdo
               if(o3vmr(ncol,k) .le. 0.)o3vmr(ncol,k) = o3mmr(k)*amdo
            endif
         enddo
        else
         do k = kts, nlayers
            o3vmr(ncol,k) = o3mmr(k) * amdo
            if (k.le.kte) o31d(k) = o3vmr(ncol,k)
         enddo
        endif

! output o3 for o3input=0
        IF (o3input.ne.2) THEN
            DO K=kts,kte
               O33D(I,K,J)=O31D(K)
            ENDDO
        ENDIF

! Set surface emissivity in each RRTMG longwave band
         do nb = 1, nbndlw
            emis(ncol, nb) = emiss(i,j)
         enddo

! Define cloud optical properties for radiation (inflglw = 0)
! This is approach used with older RRTM_LW;
! Cloud and precipitation paths in g/m2 
! qi=0 if no ice phase
! qs=0 if no ice phase
         if (inflglw .eq. 0) then
            do k = kts,kte
               ro = p1d(k) / (r * t1d(k))*100. 
               dz = dz1d(k)
               clwp(k) = ro*qc1d(k)*dz*1000.         
               ciwp(k) = ro*qi1d(k)*dz*1000.         
               plwp(k) = (ro*qr1d(k))**0.75*dz*1000. 
               piwp(k) = (ro*qs1d(k))**0.75*dz*1000. 
            enddo

! Cloud fraction and cloud optical depth; old approach used with RRTM_LW
            do k = kts, kte
               cldfrac(ncol,k) = cldfra1d(k)
               do nb = 1, nbndlw
                  taucld(nb,ncol,k) = abcw*clwp(k) + abice*ciwp(k) & 
                            +abrn*plwp(k) + absn*piwp(k) 
                  if (taucld(nb,ncol,k) .gt. 0.01) cldfrac(ncol,k) = 1. 
               enddo
            enddo

! Zero out cloud physical property arrays; not used when passing optical properties
! into radiation
            do k = kts, kte
               clwpth(ncol,k) = 0.0
               ciwpth(ncol,k) = 0.0
               rel(ncol,k) = 10.0
               rei(ncol,k) = 10.0
            enddo
         endif

! Define cloud physical properties for radiation (inflglw = 1 or 2)
! Cloud fraction
! Set cloud arrays if passing cloud physical properties into radiation
         if (inflglw .gt. 0) then 
            do k = kts, kte
               cldfrac(ncol,k) = cldfra1d(k)
            enddo

! Compute cloud water/ice paths and particle sizes for input to radiation (CAM method)
            pcols = ncol
            pver = kte - kts + 1
            gravmks = g
            landfrac(ncol) = 2.-XLAND(I,J)
            landm(ncol) = landfrac(ncol)
            snowh(ncol) = 0.001*SNOW(I,J)
            icefrac(ncol) = XICE(I,J)

! From module_ra_cam: Convert liquid and ice mixing ratios to water paths;
! pdel is in mb here; convert back to Pa (*100.)
! Water paths are in units of g/m2
! snow added as ice cloud (JD 091022)
            do k = kts, kte
               gicewp = (qi1d(k)+qs1d(k)) * pdel(ncol,k)*100.0 / gravmks * 1000.0     ! Grid box ice water path.
               gliqwp = qc1d(k) * pdel(ncol,k)*100.0 / gravmks * 1000.0     ! Grid box liquid water path.
               cicewp(ncol,k) = gicewp / max(0.01,cldfrac(ncol,k))               ! In-cloud ice water path.
               cliqwp(ncol,k) = gliqwp / max(0.01,cldfrac(ncol,k))               ! In-cloud liquid water path.
            end do


! Mukul
!..The ice water path is already sum of cloud ice and snow, but when we have explicit
!.. ice effective radius, overwrite the ice path with only the cloud ice variable,
!.. leaving out the snow for its own effect.
           if(iceflglw.ge.4)then
              do k = kts, kte
                     gicewp = qi1d(k) * pdel(ncol,k)*100.0 / gravmks * 1000.0     ! Grid box ice water path.
                     cicewp(ncol,k) = gicewp / max(0.01,cldfrac(ncol,k))               ! In-cloud ice water path.
              end do
           end if

!..Here the snow path is adjusted if (radiation) effective radius of snow is
!.. larger than what we currently have in the lookup tables.  Since mass goes
!.. rather close to diameter squared, adjust the mixing ratio of snow used
!.. to compute its water path in combination with the max diameter.  Not a
!.. perfect fix, but certainly better than using all snow mass when diameter is
!.. far larger than table currently contains and crystal sizes much larger than
!.. about 140 microns have lesser impact than those much smaller sizes.

           if(iceflglw.eq.5)then
              do k = kts, kte
                 snow_mass_factor = 0.99        ! Assume 1% of snow overlaps the cloud ice category
                 gicewp = gicewp + (qs1d(k)*(1.0-snow_mass_factor) * pdel(ncol,k)*100.0 / gravmks * 1000.0)
                 if (resnow1d(ncol,k) .gt. 130.)then
                     snow_mass_factor = MIN(snow_mass_factor,                       &
     &                         (130.0/resnow1d(ncol,k))*(130.0/resnow1d(ncol,k)))
                     resnow1d(ncol,k)   = 130.0
                     IF ( wrf_dm_on_monitor() ) THEN
                       WRITE(message,*)'RRTMG:  reducing snow mass (cloud path) to ', &
                                       nint(snow_mass_factor*100.), ' percent of full value'
                       call wrf_debug(150, message)
                     ENDIF
                 endif
                 gsnowp = qs1d(k) * snow_mass_factor * pdel(ncol,k)*100.0 / gravmks * 1000.0     ! Grid box snow water path.
                 csnowp(ncol,k) = gsnowp / max(0.01,cldfrac(ncol,k))
              end do
           end if


!link the aerosol feedback to cloud  -czhao
  if( PRESENT( progn ) ) then
    if (progn == 1) then
!jdfcz     if(prescribe==0) then

      pi = 4.*atan(1.0)
      third=1./3.
      rhoh2o=1.e3
      relconst=3/(4.*pi*rhoh2o)
!     minimun liquid water path to calculate rel
!     corresponds to optical depth of 1.e-3 for radius 4 microns.
      lwpmin=3.e-5
      do k = kts, kte
         reliq(ncol,k) = 10.
         if( PRESENT( F_QNDROP ) ) then
            if( F_QNDROP ) then
              if ( qc1d(k)*pdel(ncol,k).gt.lwpmin.and. &
                   qndrop1d(k).gt.1000. ) then
               reliq(ncol,k)=(relconst*qc1d(k)/qndrop1d(k))**third ! effective radius in m
!           apply scaling from Martin et al., JAS 51, 1830.
               reliq(ncol,k)=1.1*reliq(ncol,k)
               reliq(ncol,k)=reliq(ncol,k)*1.e6 ! convert from m to microns
               reliq(ncol,k)=max(reliq(ncol,k),4.)
               reliq(ncol,k)=min(reliq(ncol,k),20.)
              end if
            end if
         end if
      end do
!jdfcz     else ! prescribe 
! following Kiehl
!     call relcalc(ncol, pcols, pver, tlay, landfrac, landm, icefrac, reliq, snowh)
!     write(0,*) 'lw prescribe aerosol',maxval(qndrop3d)
!jdfcz     endif
    else  ! progn   
      call relcalc(ncol, pcols, pver, tlay, landfrac, landm, icefrac, reliq, snowh)
    endif
  else   !present(progn) 
      call relcalc(ncol, pcols, pver, tlay, landfrac, landm, icefrac, reliq, snowh)
  endif

! following Kristjansson and Mitchell
            call reicalc(ncol, pcols, pver, tlay, reice)


!..If we already have effective radius of cloud and ice, then just overwrite what
!.. was computed in the relcalc and reicalc subroutines above.

      if (inflglw .ge. 3) then
         do k = kts, kte
            reliq(ncol,k) = recloud1d(ncol,k)
         end do
      endif
#if (EM_CORE==1) 
      if (iceflglw .ge. 4) then
#else
      if (iceflglw .ge. 3) then    !BSF: was .ge. 4
#endif

         do k = kts, kte
            reice(ncol,k) = reice1d(ncol,k)
         end do
      endif

! Limit upper bound of reice for Fu ice parameterization and convert
! from effective radius to generalized effective size (*1.0315; Fu, 1996)
            if (iceflglw .eq. 3) then
               do k = kts, kte
                  reice(ncol,k) = reice(ncol,k) * 1.0315
                  reice(ncol,k) = min(140.0,reice(ncol,k))
               end do
            endif
!if CAMMGMP is used, use output from CAMMGMP
            if(is_CAMMGMP_used) then
               do k = kts, kte
                  if ( qi1d(k) .gt. 1.e-20 .or. qs1d(k) .gt. 1.e-20) then
                     reice(ncol,k) = iradius(i,k,j)
                  else
                     reice(ncol,k) = 25.
                  end if
                  reice(ncol,k) = max(5., min(140.0,reice(ncol,k)))
                  if ( qc1d(k) .gt. 1.e-20) then
                     reliq(ncol,k) = lradius(i,k,j)
                  else
                     reliq(ncol,k) = 10.
                  end if
                  reliq(ncol,k) = max(2.5, min(60.0,reliq(ncol,k)))
               enddo
            endif


! Set cloud physical property arrays
            do k = kts, kte
               clwpth(ncol,k) = cliqwp(ncol,k)
               ciwpth(ncol,k) = cicewp(ncol,k)
               rel(ncol,k) = reliq(ncol,k)
               rei(ncol,k) = reice(ncol,k)
            enddo

!Mukul
            if (inflglw .eq. 5) then
               do k = kts, kte
                  cswpth(ncol,k) = csnowp(ncol,k)
                  res(ncol,k) = resnow1d(ncol,k)
               end do
            else
               do k = kts, kte
                  cswpth(ncol,k) = 0.
                  res(ncol,k) = 10.
               end do
            endif

! Zero out cloud optical properties here; not used when passing physical properties
! to radiation and taucld is calculated in radiation 
            do k = kts, kte
               do nb = 1, nbndlw
                  taucld(nb,ncol,k) = 0.0
               enddo
            enddo
         endif

! No clouds are allowed in the extra layer from model top to TOA
         ! Steven Cavallo: Edited out for buffer adjustment below
         if ( 1 == 0 ) then


         clwpth(ncol,kte+1) = 0.
         ciwpth(ncol,kte+1) = 0.
         cswpth(ncol,kte+1) = 0.
         rel(ncol,kte+1) = 10.
         rei(ncol,kte+1) = 10.
         res(ncol,kte+1) = 10.
         cldfrac(ncol,kte+1) = 0.
         do nb = 1, nbndlw
            taucld(nb,ncol,kte+1) = 0.
         enddo

         endif

         ! Buffer adjustment. Steven Cavallo December 2010
         do k=kte+1,nlayers
            clwpth(ncol,k) = 0.
            ciwpth(ncol,k) = 0.
            cswpth(ncol,k) = 0.
            rel(ncol,k) = 10.
            rei(ncol,k) = 10.
            res(ncol,k) = 10.
            cldfrac(ncol,k) = 0.
            do nb = 1,nbndlw
               taucld(nb,ncol,k) = 0.
            enddo
         enddo   

         iplon = 1
         irng = 0
         permuteseed = 150

! Sub-column generator for McICA
         lat=xlat(i,j) !retrieve scalar latitude for column calculation
         call mcica_subcol_lw(iplon, ncol, nlay, icld, permuteseed, irng, play, &
                       cldfrac, ciwpth, clwpth, cswpth, rei, rel, res, taucld, &
                       hgt, idcor, juldat, lat, &
                       cldfmcl, ciwpmcl, clwpmcl, cswpmcl, reicmcl, relqmcl, resnmcl, taucmcl)

!--------------------------------------------------------------------------
! Aerosol optical depth, single scattering albedo and asymmetry parameter -czhao 03/2010
!--------------------------------------------------------------------------
! Aerosol optical depth by layer for each RRTMG longwave band
! No aerosols in layer above model top (kte+1)
! Steven Cavallo: Upper bound of loop changed to nlayers from kte+1
!        do nb = 1, nbndlw
!           do k = kts, kte+1
!              tauaer(ncol,k,nb) = 0.
!           enddo
!        enddo

! ... Aerosol effects. Added aerosol feedbacks from Chem , 03/2010 -czhao
!
      do nb = 1, nbndlw
      do k = kts,nlayers
         tauaer(ncol,k,nb) = 0.
      end do
      end do

#if ( WRF_CHEM == 1 )
   IF ( AER_RA_FEEDBACK == 1) then
!     do nb = 1, nbndlw 
      do k = kts,kte      !wig
        if(tauaerlw1(i,k,j).gt.thresh .and. tauaerlw16(i,k,j).gt.thresh) then
          tauaer(ncol,k,1)=tauaerlw1(i,k,j)
          tauaer(ncol,k,2)=tauaerlw2(i,k,j)
          tauaer(ncol,k,3)=tauaerlw3(i,k,j)
          tauaer(ncol,k,4)=tauaerlw4(i,k,j)
          tauaer(ncol,k,5)=tauaerlw5(i,k,j)
          tauaer(ncol,k,6)=tauaerlw6(i,k,j)
          tauaer(ncol,k,7)=tauaerlw7(i,k,j)
          tauaer(ncol,k,8)=tauaerlw8(i,k,j)
          tauaer(ncol,k,9)=tauaerlw9(i,k,j)
          tauaer(ncol,k,10)=tauaerlw10(i,k,j)
          tauaer(ncol,k,11)=tauaerlw11(i,k,j)
          tauaer(ncol,k,12)=tauaerlw12(i,k,j)
          tauaer(ncol,k,13)=tauaerlw13(i,k,j)
          tauaer(ncol,k,14)=tauaerlw14(i,k,j)
          tauaer(ncol,k,15)=tauaerlw15(i,k,j)
          tauaer(ncol,k,16)=tauaerlw16(i,k,j)
        endif
      enddo ! k
!     end do ! nb

!wig beg
      do nb = 1, nbndlw
         slope = 0.  !use slope as a sum holder
         do k = kts,kte
            slope = slope + tauaer(ncol,k,nb)
         end do
         if( slope < 0. ) then
            write(msg,'("ERROR: Negative total lw optical depth of ",f8.2," at point i,j,nb=",3i5)') slope,i,j,nb
            call wrf_error_fatal(msg)
         else if( slope > 5. ) then
            call wrf_message("-------------------------")
            write(msg,'("WARNING: Large total lw optical depth of ",f8.2," at point i,j,nb=",3i5)') slope,i,j,nb
            call wrf_message(msg)

            call wrf_message("Diagnostics 1: k, tauaerlw1, tauaerlw16")
            do k=kts,kte
               write(msg,'(i4,2f8.2)') k, tauaerlw1(i,k,j), tauaerlw16(i,k,j)
               call wrf_message(msg)
            end do
            call wrf_message("-------------------------")
         endif
      enddo  ! nb
      endif  ! aer_ra_feedback
#endif

! gpu: pack this column's solver inputs into the batch
         play_p (bcol, :) = play (1, :)
         tlay_p (bcol, :) = tlay (1, :)
         plev_p (bcol, :) = plev (1, :)
         tlev_p (bcol, :) = tlev (1, :)
         tsfc_p (bcol)    = tsfc (1)
         h2ovmr_p  (bcol, :) = h2ovmr  (1, :)
         o3vmr_p   (bcol, :) = o3vmr   (1, :)
         co2vmr_p  (bcol, :) = co2vmr  (1, :)
         o2vmr_p   (bcol, :) = o2vmr   (1, :)
         ch4vmr_p  (bcol, :) = ch4vmr  (1, :)
         n2ovmr_p  (bcol, :) = n2ovmr  (1, :)
         cfc11vmr_p(bcol, :) = cfc11vmr(1, :)
         cfc12vmr_p(bcol, :) = cfc12vmr(1, :)
         cfc22vmr_p(bcol, :) = cfc22vmr(1, :)
         ccl4vmr_p (bcol, :) = ccl4vmr (1, :)
         emis_p    (bcol, :) = emis    (1, :)
         cldfmcl_p (:, bcol, :) = cldfmcl(:, 1, :)   ! ncol is the middle dim here
         taucmcl_p (:, bcol, :) = taucmcl(:, 1, :)
         ciwpmcl_p (:, bcol, :) = ciwpmcl(:, 1, :)
         clwpmcl_p (:, bcol, :) = clwpmcl(:, 1, :)
         cswpmcl_p (:, bcol, :) = cswpmcl(:, 1, :)
         reicmcl_p (bcol, :) = reicmcl(1, :)
         relqmcl_p (bcol, :) = relqmcl(1, :)
         resnmcl_p (bcol, :) = resnmcl(1, :)
         tauaer_p  (bcol, :, :) = tauaer(1, :, :)

      end do pack_loop

! gpu: solve phase - one single GPU launch over the batch (TODO)
      ! call wrf_gpu_lw( play_p, SIZE(play_p) )

      solve_loop: do bcol = 1, nbatch

! Call RRTMG longwave radiation model (into packed arrays)
         call rrtmg_lw &
            (ncol    ,nlay    ,icld    , &
             play_p(bcol:bcol,:)     ,plev_p(bcol:bcol,:)     , &
             tlay_p(bcol:bcol,:)     ,tlev_p(bcol:bcol,:)     , &
             tsfc_p(bcol:bcol)       , &
             h2ovmr_p(bcol:bcol,:)   ,o3vmr_p(bcol:bcol,:)    , &
             co2vmr_p(bcol:bcol,:)   ,ch4vmr_p(bcol:bcol,:)   , &
             n2ovmr_p(bcol:bcol,:)   ,o2vmr_p(bcol:bcol,:)    , &
             cfc11vmr_p(bcol:bcol,:) ,cfc12vmr_p(bcol:bcol,:) , &
             cfc22vmr_p(bcol:bcol,:) ,ccl4vmr_p(bcol:bcol,:)  , &
             emis_p(bcol:bcol,:)     , &
             inflglw ,iceflglw,liqflglw, &
             cldfmcl_p(:,bcol:bcol,:) , &
             taucmcl_p(:,bcol:bcol,:) ,ciwpmcl_p(:,bcol:bcol,:) , &
             clwpmcl_p(:,bcol:bcol,:) ,cswpmcl_p(:,bcol:bcol,:) , &
             reicmcl_p(bcol:bcol,:)  ,relqmcl_p(bcol:bcol,:)  , &
             resnmcl_p(bcol:bcol,:)  , &
             tauaer_p(bcol:bcol,:,:) , &
             uflx_p(bcol:bcol,:)     ,dflx_p(bcol:bcol,:)     , &
             hr_p(bcol:bcol,:)       , &
             uflxc_p(bcol:bcol,:)    ,dflxc_p(bcol:bcol,:)    , &
             hrc_p(bcol:bcol,:)      , &
             uflxcln_p(bcol:bcol,:)  ,dflxcln_p(bcol:bcol,:)  , &
             calc_clean_atm_diag )

      end do solve_loop

! gpu: unpack phase - scatter the batch outputs into the WRF arrays
      unpack_loop: do bcol = 1, nbatch
         i = i_of(bcol)
         j = j_of(bcol)

! Output downard surface flux, and outgoing longwave flux and cloud forcing
! at the top of atmosphere (W/m2)
         glw(i,j) = dflx_p(bcol,1)
!         olr(i,j) = uflx(1,kte+2)
!         lwcf(i,j) = uflxc(1,kte+2) - uflx(1,kte+2)
! Steven Cavallo: Changed OLR to be valid at the top of atmosphere instead 
! of top of model.  Dec 2010.
         olr(i,j) = uflx_p(bcol,nlayers+1)
         lwcf(i,j) = uflxc_p(bcol,nlayers+1) - uflx_p(bcol,nlayers+1)

         if (present(lwupt)) then 
! Output up and down toa fluxes for total and clear sky
! nlayers+1 represents value at 0 mb
            lwupt(i,j)     = uflx_p(bcol,nlayers+1)
            lwuptc(i,j)    = uflxc_p(bcol,nlayers+1)
            lwdnt(i,j)     = dflx_p(bcol,nlayers+1)
            lwdntc(i,j)    = dflxc_p(bcol,nlayers+1)
! Output up and down surface fluxes for total and clear sky
            lwupb(i,j)     = uflx_p(bcol,1)
            lwupbc(i,j)    = uflxc_p(bcol,1)
            lwdnb(i,j)     = dflx_p(bcol,1)
            lwdnbc(i,j)    = dflxc_p(bcol,1)
			if(calc_clean_atm_diag .gt. 0)then
! Output up and down toa fluxes for clean sky
            	lwuptcln(i,j)  = uflxcln_p(bcol,nlayers+1)
            	lwdntcln(i,j)  = dflxcln_p(bcol,nlayers+1)
! Output up and down surface fluxes for clean sky
            	lwupbcln(i,j)  = uflxcln_p(bcol,1)
            	lwdnbcln(i,j)  = dflxcln_p(bcol,1)
			end if
         endif

! Output up and down layer fluxes for total and clear sky.
! Vertical ordering is from bottom to top in units of W m-2. 
         if ( present (lwupflx) ) then
         do k=kts,kte+2
            lwupflx(i,k,j)  = uflx_p(bcol,k)
            lwupflxc(i,k,j) = uflxc_p(bcol,k)
            lwdnflx(i,k,j)  = dflx_p(bcol,k)
            lwdnflxc(i,k,j) = dflxc_p(bcol,k)
         enddo
         endif

! Output heating rate tendency; convert heating rate from K/d to K/s
! Heating rate arrays are ordered vertically from bottom to top here. 
         do k=kts,kte
            tten1d(k) = hr_p(bcol,k)/86400.
            rthratenlw(i,k,j) = tten1d(k)/pi3d(i,k,j)
            tten1d(k) = hrc_p(bcol,k)/86400.
            rthratenlwc(i,k,j) = tten1d(k)/pi3d(i,k,j)
         enddo

!
      end do unpack_loop
   end do batch_loop
   
! gpu: free the solver batch
   DEALLOCATE( play_p, tlay_p, plev_p, tlev_p, tsfc_p,                         &
               h2ovmr_p, o3vmr_p, co2vmr_p, o2vmr_p, ch4vmr_p, n2ovmr_p,       &
               cfc11vmr_p, cfc12vmr_p, cfc22vmr_p, ccl4vmr_p, emis_p,          &
               cldfmcl_p, taucmcl_p, ciwpmcl_p, clwpmcl_p, cswpmcl_p,          &
               reicmcl_p, relqmcl_p, resnmcl_p, tauaer_p,                      &
               uflx_p, dflx_p, uflxc_p, dflxc_p, uflxcln_p, dflxcln_p,         &
               hr_p, hrc_p )

!-------------------------------------------------------------------

   END SUBROUTINE RRTMG_LWRAD_GPU

END MODULE module_ra_rrtmg_lw_gpu
