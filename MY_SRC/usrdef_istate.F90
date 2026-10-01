MODULE usrdef_istate
   !!======================================================================
   !!                   ***  MODULE  usrdef_istate   ***
   !!
   !!                     ===  GYRE configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-03  (S. Flavoni) Original code
   !!                 ! 2020-11  (S. Techene, G. Madec) separate tsuv from ssh
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE splines
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate       ! called in istate.F90
   PUBLIC   usr_def_istate_ssh   ! called by domqco.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here GYRE configuration example : (double gyre with rotated domain)
      !!
      !! ** Method  : - set temprature field
      !!              - set salinity   field
      !!----------------------------------------------------------------------
      REAL(wp),              DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp),              DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp),              DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp),              DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s] 
      REAL(wp),              DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s] 
      REAL(wp), ALLOCATABLE, DIMENSION(:)                               ::   zdep, zsal, ztem
      REAL(wp)                                                          ::   a0, a1, a2, b0, b1, b2
      !
      INTEGER :: ji, jj, jk  ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : analytical definition of initial state '
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   Ocean at rest, with an horizontally uniform T and S profiles'
      !
      pu  (:,:,:) = 0._wp           ! ocean at rest
      pv  (:,:,:) = 0._wp
      !
      ! T/S profiles representative of Weddel Gyre in Feb (from GOSI9, Guiavarc'h et at. 2025)
      ALLOCATE( zdep(75), zsal(75), ztem(75) )

      zdep(:) = (/ &
                   &    0.50,    1.55,    2.66,    3.85,    5.14,    6.54, &
                   &    8.09,    9.82,   11.77,   13.99,   16.52,   19.42, &
                   &   22.75,   26.55,   30.87,   35.74,   41.18,   47.21, &
                   &   53.85,   61.11,   69.02,   77.61,   86.92,   97.04, &
                   &  108.03,  120.00,  133.07,  147.40,  163.16,  180.54, &
                   &  199.78,  221.14,  244.89,  271.35,  300.88,  333.86, &
                   &  370.68,  411.79,  457.62,  508.63,  565.29,  628.02, &
                   &  697.25,  773.36,  856.67,  947.44,  1045.85,1151.99, &
                   & 1265.86, 1387.37, 1516.36, 1652.56, 1795.67, 1945.29, &
                   & 2101.02, 2262.42, 2429.02, 2600.38, 2776.03, 2955.57, &
                   & 3138.56, 3324.64, 3513.44, 3704.65, 3897.98, 4093.15, &
                   & 4289.95, 4488.15, 4687.58, 4888.06, 5089.47, 5291.68, &
                   & 5494.57, 5698.06, 5902.05 &
                   & /)

      ztem(:) = (/ &
                   &  0.0280,  0.0227,  0.0171,  0.0111,  0.0048, -0.0014, &
                   & -0.0083, -0.0159, -0.0240, -0.0331, -0.0470, -0.0654, &
                   & -0.0943, -0.1415, -0.2277, -0.3852, -0.5965, -0.8201, &
                   & -1.0382, -1.2168, -1.3276, -1.3812, -1.3842, -1.3266, &
                   & -1.1830, -1.0074, -0.7929, -0.5426, -0.3335, -0.1374, &
                   &  0.0322,  0.1387,  0.2129,  0.2720,  0.3443,  0.3798, &
                   &  0.4137,  0.4261,  0.4305,  0.4462,  0.4423,  0.4277, &
                   &  0.4053,  0.3749,  0.3372,  0.2916,  0.2520,  0.2050, &
                   &  0.1493,  0.0889,  0.0321, -0.0268, -0.0834, -0.1430, &
                   & -0.2021, -0.2528, -0.3024, -0.3498, -0.3922, -0.4339, &
                   & -0.4736, -0.5114, -0.5480, -0.5851, -0.6255, -0.6662, &
                   & -0.7107, -0.7649, -0.8049, -0.8366, -0.8409, -0.8452, &
                   & -0.7990, -0.7861,  -0.7861 &
                   & /)

      zsal(:) = (/ &
                   & 34.0206, 34.0238, 34.0272, 34.0308, 34.0349, 34.0410, &
                   & 34.0477, 34.0551, 34.0637, 34.0733, 34.0859, 34.1014, &
                   & 34.1274, 34.1614, 34.2048, 34.2588, 34.3211, 34.3761, &
                   & 34.4232, 34.4648, 34.4962, 34.5223, 34.5476, 34.5720, &
                   & 34.5979, 34.6276, 34.6585, 34.6914, 34.7189, 34.7443, &
                   & 34.7653, 34.7798, 34.7920, 34.8024, 34.8138, 34.8221, &
                   & 34.8291, 34.8353, 34.8404, 34.8464, 34.8493, 34.8516, &
                   & 34.8530, 34.8540, 34.8540, 34.8535, 34.8527, 34.8510, &
                   & 34.8487, 34.8466, 34.8452, 34.8435, 34.8418, 34.8399, &
                   & 34.8399, 34.8387, 34.8373, 34.8356, 34.8340, 34.8333, &
                   & 34.8323, 34.8310, 34.8297, 34.8285, 34.8274, 34.8267, &
                   & 34.8257, 34.8245, 34.8228, 34.8227, 34.8228, 34.8238, &
                   & 34.8320, 34.8290, 34.8290 &
                   & /)

      ! Use splines to interpolate
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         pts(ji,jj,:,jp_tem) = spline3(zdep,ztem,pdept(ji,jj,:))
         pts(ji,jj,:,jp_sal) = spline3(zdep,zsal,pdept(ji,jj,:))
      END_2D
      !   
   END SUBROUTINE usr_def_istate

   
   SUBROUTINE usr_def_istate_ssh( ptmask, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate_ssh  ***
      !! 
      !! ** Purpose :   Initialization of ssh
      !!
      !! ** Method  :   Set ssh as null, ptmask is required for test cases
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask   [m]
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height   [m]
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate_ssh : GYRE configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~   Ocean at rest, ssh is zero'
      !
      ! Sea level:
      pssh(:,:) = 0._wp
      !
   END SUBROUTINE usr_def_istate_ssh

   !!======================================================================
END MODULE usrdef_istate
