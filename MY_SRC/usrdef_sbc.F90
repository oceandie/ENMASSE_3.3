MODULE usrdef_sbc
   !!======================================================================
   !!                     ***  MODULE  usrdef_sbc  ***
   !!
   !!                     ===  GYRE configuration  ===
   !!
   !! User defined :   surface forcing of a user configuration
   !!======================================================================
   !! History :  4.0   ! 2016-03  (S. Flavoni, G. Madec)  user defined interface
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usrdef_sbc    : user defined surface bounday conditions in GYRE case
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran    ! for glob_2Dsum

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usrdef_sbc_oce       ! routine called in sbcmod module
   PUBLIC   usrdef_sbc_ice_tau   ! routine called by icestp.F90 for ice dynamics
   PUBLIC   usrdef_sbc_ice_flx   ! routine called by icestp.F90 for ice thermo

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usrdef_sbc_oce( kt, Kbb )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE usrdef_sbc  ***
      !!              
      !! ** Purpose :   provide at each time-step the GYRE surface boundary
      !!              condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   analytical seasonal cycle for GYRE configuration.
      !!                CAUTION : never mask the surface stress field !
      !!
      !! ** Action  : - set the ocean surface boundary condition, i.e.   
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, sfx
      !!
      !! Reference : Hazeleger, W., and S. Drijfhout, JPO, 30, 677-695, 2000.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      INTEGER, INTENT(in) ::   Kbb  ! ocean time index
      !!
      !!----------------------------------------------------------------------

      emp (:,:) = 0.0_wp        ! freshwater flux
      sfx (:,:) = 0.0_wp        ! no salt flux
      qns (:,:) = 0.0_wp        ! non solar heat flux
      qsr (:,:) = 0.0_wp        ! solar heat flux
      utau(:,:) = 0.0_wp
      vtau(:,:) = 0.0_wp
      taum(:,:) = 0.0_wp
      wndm(:,:) = 0.0_wp
      
      !CALL lbc_lnk( 'usrdef_sbc', taum(:,:), 'T', 1. , wndm(:,:), 'T', 1.)
 
      ! ---------------------------------- !
      !  control print at first time-step  !
      ! ---------------------------------- !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'usrdef_sbc_oce : all surface fluxes set to 0'
         WRITE(numout,*)'~~~~~~~~~~~ '
      ENDIF
      !
   END SUBROUTINE usrdef_sbc_oce


   SUBROUTINE usrdef_sbc_ice_tau( kt )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
   END SUBROUTINE usrdef_sbc_ice_tau


   SUBROUTINE usrdef_sbc_ice_flx( kt, phs, phi )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phs    ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phi    ! ice thickness
   END SUBROUTINE usrdef_sbc_ice_flx

   !!======================================================================
END MODULE usrdef_sbc

