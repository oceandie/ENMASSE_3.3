MODULE diahsb
   !!======================================================================
   !!                       ***  MODULE  diahsb  ***
   !! Ocean diagnostics: Heat, salt and volume budgets
   !!======================================================================
   !! History :  3.3  ! 2010-09  (M. Leclair)  Original code
   !!                 ! 2012-10  (C. Rousset)  add iom_put
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_hsb       : Diagnose the conservation of ocean heat and salt contents, and volume
   !!   dia_hsb_rst   : Read or write DIA file in restart file
   !!   dia_hsb_init  : Initialization of the conservation diagnostic
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE sbc_oce        ! surface thermohaline fluxes
   USE isf_oce        ! ice shelf fluxes
   USE sbcrnf         ! river runoff
   USE traqsr         ! penetrative solar radiation
   USE trabbc         ! bottom boundary condition
   USE trabbc         ! bottom boundary condition
   USE restart        ! ocean restart
   USE bdy_oce , ONLY : ln_bdy, nb_bdy, idx_bdy, bdytmask
#if defined key_agrif
   USE agrif_oce, ONLY : ssh_upd_msk, ssh_upd_lost, nn_shift_bar, e1e2t_frac
#endif
   !
   USE iom            ! I/O manager
   USE in_out_manager ! I/O manager
   USE lib_fortran    ! glob_2Dsum
   USE lib_mpp        ! distributed memory computing library
   USE lbclnk         ! lateral boundary conditions (volume budget masks)
   USE timing         ! preformance summary

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_hsb        ! routine called by step.F90
   PUBLIC   dia_hsb_init   ! routine called by nemogcm.F90
   PUBLIC   dia_hsb_post   ! routine called by stprk3.F90 after Agrif_update_all

   LOGICAL, PUBLIC ::   l_diahsb   !: check the heat and salt budgets
   LOGICAL, PUBLIC ::   l_diahsb_init ! =T start budget from kt = nit000

   REAL(wp) ::   surf_tot              ! ocean surface
   REAL(wp) ::   frc_t, frc_s, frc_v   ! global forcing trends
   REAL(wp) ::   frc_wn_t, frc_wn_s    ! global forcing trends
   !
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE ::   surf
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE ::   surf_ini      , ssh_ini          !
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE ::   ssh_hc_loc_ini, ssh_sc_loc_ini   !
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   hc_loc_ini, sc_loc_ini, e3t_ini  !
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   tmask_ini

   ! Barotropic volume budget over the prognostic domain (open boundaries + AGRIF)
   ! All volumes in m3, cumulated from the first time step of the current run segment
   INTEGER  ::   numvbg = -1           ! logical unit of the ASCII budget file
   LOGICAL  ::   ll_vbg_first = .TRUE. ! first call of dia_hsb_post (area of U written to the header)
   REAL(wp) ::   vbg_area              ! A_Omega = sum_Omega( e1e2t ) (m2)
   REAL(wp) ::   vbg_v0                ! sum_Omega( e1e2t * ssh ) at nit000
   REAL(wp) ::   vbg_frc               ! cumulated volume source -sum_Omega( e1e2t * ssh_frc ) * rn_Dt
   REAL(wp) ::   vbg_isf               ! ice-shelf part of vbg_frc (information only)
   REAL(wp) ::   vbg_obc               ! cumulated lateral inflow through the edge of Omega (BDY or AGRIF interface)
   REAL(wp) ::   vbg_upd               ! cumulated increment of sum_Omega( e1e2t * ssh ) due to AGRIF update
   REAL(wp) ::   vbg_frcu              ! cumulated volume source over the AGRIF update region U
   REAL(wp) ::   vbg_obcu              ! cumulated inflow into U through its edge
   REAL(wp) ::   vbg_res               ! residual = V(t) - V(0) - vbg_frc - vbg_obc - vbg_upd
   REAL(wp) ::   vbg_resmax            ! max |per-step pointwise continuity residual| over Omega (m)
   REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   vbg_msk    ! Omega: prognostic ssh points (jpi,jpj)
   REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   ssh_prev   ! ssh at the start of the step (after AGRIF update)
   REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   ssh_pre    ! ssh at the end of the step, before AGRIF update

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_hsb( kt, Kbb, Kmm )
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_hsb  ***
      !!
      !! ** Purpose: Compute the ocean global heat content, salt content and volume conservation
      !!
      !! ** Method : - Compute the deviation of heat content, salt content and volume
      !!	            at the current time step from their values at nit000
      !!	            - Compute the contribution of forcing and remove it from these deviations
      !!
      !!---------------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt         ! ocean time-step index
      INTEGER, INTENT(in) ::   Kbb, Kmm   ! ocean time level indices
      !
      INTEGER    ::   ji, jj, jk                  ! dummy loop indice
      REAL(wp)   ::   zdiff_hc    , zdiff_sc      ! heat and salt content variations
      REAL(wp)   ::   zdiff_hc1   , zdiff_sc1     !  -         -     -        -
      REAL(wp)   ::   zdiff_v1    , zdiff_v2      ! volume variation
      REAL(wp)   ::   zerr_hc1    , zerr_sc1      ! heat and salt content misfit
      REAL(wp)   ::   zvol_tot                    ! volume
      REAL(wp)   ::   z_frc_trd_t , z_frc_trd_s   !    -     -
      REAL(wp)   ::   z_frc_trd_v                 !    -     -
      REAL(wp)   ::   z_wn_trd_t , z_wn_trd_s     !    -     -
      REAL(wp)   ::   z_ssh_hc , z_ssh_sc         !    -     -
      REAL(wp), DIMENSION(A2D(0),17)      ::   ztmp
      REAL(wp), DIMENSION(17)             ::   zbg
      !!---------------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('dia_hsb')
      !
      DO_2D( 0, 0, 0, 0 )
         ztmp (ji,jj,:)   = 0._wp ! should be better coded
      END_2D
      !
      ! ------------------------- !
      ! 1 - Trends due to forcing !
      ! ------------------------- !
      ! prepare trends
      DO_2D( 0, 0, 0, 0 )
                        ztmp(ji,jj,1)  =                 emp(ji,jj)                         ! volume
         IF( ln_rnf )   ztmp(ji,jj,1)  = ztmp(ji,jj,1) - rnf(ji,jj)
         IF( ln_isf )   ztmp(ji,jj,1)  = ztmp(ji,jj,1) - fwfisf_cav(ji,jj) - fwfisf_par(ji,jj)
         ztmp(ji,jj,1)  =  - r1_rho0 * ztmp(ji,jj,1) * surf(ji,jj)
#if defined key_RK3
         ztmp(ji,jj,2)  =   r1_rho0_rcp * qns(ji,jj) * surf(ji,jj)                           ! heat
         ztmp(ji,jj,3)  =   r1_rho0     * sfx(ji,jj) * surf(ji,jj)                           ! salt
         IF( lk_linssh ) THEN
            ztmp(ji,jj,2)  =  r1_rho0 * emp(ji,jj) * ts(ji,jj,1,jp_tem,Kmm) * surf(ji,jj)   ! heat
            ztmp(ji,jj,3)  =  r1_rho0 * emp(ji,jj) * ts(ji,jj,1,jp_sal,Kmm) * surf(ji,jj)   ! salt
         ENDIF
#else
         ztmp(ji,jj,2)  =   sbc_tsc(ji,jj,jp_tem) * surf(ji,jj)                              ! heat
         ztmp(ji,jj,3)  =   sbc_tsc(ji,jj,jp_sal) * surf(ji,jj)                              ! salt
#endif
      END_2D
      IF( ln_rnf     ) THEN
         DO_2D( 0, 0, 0, 0 )
            ztmp(ji,jj,4) = rnf_tsc(ji,jj,jp_tem) * surf(ji,jj)    ! runoff temp
            ztmp(ji,jj,5) = rnf_tsc(ji,jj,jp_sal) * surf(ji,jj)    ! runoff salt
         END_2D
      END IF
      IF( ln_isf     ) THEN
         DO_2D( 0, 0, 0, 0 )
            ztmp(ji,jj,6) = (   risf_cav_tsc(ji,jj,jp_tem) &
                &             + risf_par_tsc(ji,jj,jp_tem) ) * surf(ji,jj) ! isf temp
         END_2D
      END IF
      IF( ln_traqsr  ) THEN
         DO_2D( 0, 0, 0, 0 )
            ztmp(ji,jj,7) = r1_rho0_rcp * qsr(ji,jj) * surf(ji,jj) ! penetrative solar radiation
         END_2D
      END IF
      IF( ln_trabbc  ) THEN
         DO_2D( 0, 0, 0, 0 )
            ztmp(ji,jj,8) = qgh_trd0(ji,jj) * surf(ji,jj)             ! geothermal heat
         END_2D
      END IF
      !
      IF( lk_linssh ) THEN   ! Advection flux through fixed surface (z=0)
         IF( ln_isfcav ) THEN
            DO_2D( 0, 0, 0, 0 )
               ztmp(ji,jj,9 ) = - surf(ji,jj) * ww(ji,jj,mikt(ji,jj)) * ts(ji,jj,mikt(ji,jj),jp_tem,Kbb)
               ztmp(ji,jj,10) = - surf(ji,jj) * ww(ji,jj,mikt(ji,jj)) * ts(ji,jj,mikt(ji,jj),jp_sal,Kbb)
            END_2D
         ELSE
            DO_2D( 0, 0, 0, 0 )
               ztmp(ji,jj,9 ) = - surf(ji,jj) * ww(ji,jj,1) * ts(ji,jj,1,jp_tem,Kbb)
               ztmp(ji,jj,10) = - surf(ji,jj) * ww(ji,jj,1) * ts(ji,jj,1,jp_sal,Kbb)
            END_2D
         END IF
      ENDIF

      ! --------------------------------- !
      ! 2 -  Content variations with ssh  !
      ! --------------------------------- !
      ! glob_2Dsum is needed because you keep only the interior domain to compute the sum (iscpl)
      !
      !                    ! volume variation (calculated with ssh)
      DO_2D( 0, 0, 0, 0 )
         ztmp(ji,jj,11) = surf(ji,jj)*ssh(ji,jj,Kmm) - surf_ini(ji,jj)*ssh_ini(ji,jj)
      END_2D

      !                    ! heat & salt content variation (associated with ssh)
      IF( lk_linssh ) THEN       ! linear free surface case
         IF( ln_isfcav ) THEN          ! ISF case
            DO_2D( 0, 0, 0, 0 )
               ztmp(ji,jj,12) = surf(ji,jj) * ( ts(ji,jj,mikt(ji,jj),jp_tem,Kmm) * ssh(ji,jj,Kmm) - ssh_hc_loc_ini(ji,jj) )
               ztmp(ji,jj,13) = surf(ji,jj) * ( ts(ji,jj,mikt(ji,jj),jp_sal,Kmm) * ssh(ji,jj,Kmm) - ssh_sc_loc_ini(ji,jj) )
            END_2D
         ELSE                          ! no under ice-shelf seas
            DO_2D( 0, 0, 0, 0 )
               ztmp(ji,jj,12) = surf(ji,jj) * ( ts(ji,jj,1,jp_tem,Kmm) * ssh(ji,jj,Kmm) - ssh_hc_loc_ini(ji,jj) )
               ztmp(ji,jj,13) = surf(ji,jj) * ( ts(ji,jj,1,jp_sal,Kmm) * ssh(ji,jj,Kmm) - ssh_sc_loc_ini(ji,jj) )
            END_2D
         END IF
      ENDIF

      ! --------------------------------- !
      ! 3 -  Content variations with e3t  !
      ! --------------------------------- !
      ! glob_2Dsum is needed because you keep only the interior domain to compute the sum (iscpl)
      !
      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
         ! volume
         ztmp(ji,jj,14) = ztmp(ji,jj,14) + surf    (ji,jj) * e3t    (ji,jj,jk,Kmm)*tmask    (ji,jj,jk) &
            &                            - surf_ini(ji,jj) * e3t_ini(ji,jj,jk)    *tmask_ini(ji,jj,jk)
         ! heat
         ztmp(ji,jj,15) = ztmp(ji,jj,15) + ( surf    (ji,jj) * e3t(ji,jj,jk,Kmm)*ts(ji,jj,jk,jp_tem,Kmm) &
            &                            -   surf_ini(ji,jj) * hc_loc_ini(ji,jj,jk) )
         ! salt
         ztmp(ji,jj,16) = ztmp(ji,jj,16) + ( surf    (ji,jj) * e3t(ji,jj,jk,Kmm)*ts(ji,jj,jk,jp_sal,Kmm) &
            &                            -   surf_ini(ji,jj) * sc_loc_ini(ji,jj,jk) )
         ! total ocean volume
         ztmp(ji,jj,17) = ztmp(ji,jj,17) + surf(ji,jj) * e3t(ji,jj,jk,Kmm) * tmask(ji,jj,jk)
      END_3D

      ! ----------
      ! global sum
      ! ----------
      zbg(:) = glob_2Dsum( 'dia_hsb', ztmp, cdelay = 'ocebg' )
      !
      IF ( .NOT.(ln_mppdelay.AND.l_diahsb_init.AND.( kt==nit000 ) ) ) THEN ! skip summation at startup because of delayed global sums

         ! 1)
         z_frc_trd_v = zbg(1)  ! volume fluxes
         z_frc_trd_t = zbg(2)  ! heat fluxes
         z_frc_trd_s = zbg(3)  ! salt fluxes
         IF( ln_rnf    )   z_frc_trd_t = z_frc_trd_t + zbg(4) ! runoff heat
         IF( ln_rnf    )   z_frc_trd_s = z_frc_trd_s + zbg(5) ! runoff salt
         IF( ln_isf    )   z_frc_trd_t = z_frc_trd_t + zbg(6) ! isf heat
         IF( ln_traqsr )   z_frc_trd_t = z_frc_trd_t + zbg(7) ! penetrative solar flux
         IF( ln_trabbc )   z_frc_trd_t = z_frc_trd_t + zbg(8) ! geothermal heat
         !
         frc_v = frc_v + z_frc_trd_v * rn_Dt
         frc_t = frc_t + z_frc_trd_t * rn_Dt
         frc_s = frc_s + z_frc_trd_s * rn_Dt
         !                                          ! Advection flux through fixed surface (z=0)
         IF( lk_linssh ) THEN
            z_wn_trd_t = zbg(9)
            z_wn_trd_s = zbg(10)
            !
            frc_wn_t = frc_wn_t + z_wn_trd_t * rn_Dt
            frc_wn_s = frc_wn_s + z_wn_trd_s * rn_Dt
         ENDIF
      
         ! 2)
         zdiff_v1 = zbg(11)
         !                    ! heat & salt content variation (associated with ssh)
         IF( lk_linssh ) THEN       ! linear free surface case
            z_ssh_hc = zbg(12)
            z_ssh_sc = zbg(13)
         ENDIF
         !
         ! 3)
         zdiff_v2 = zbg(14)     ! glob_2Dsum needed as tmask and tmask_ini could be different
         zdiff_hc = zbg(15)
         zdiff_sc = zbg(16)
         zvol_tot = zbg(17)

         ! ------------------------ !
         ! 4 -  Drifts              !
         ! ------------------------ !
         zdiff_v1 = zdiff_v1 - frc_v
         IF( .NOT.lk_linssh )   zdiff_v2 = zdiff_v2 - frc_v
         zdiff_hc = zdiff_hc - frc_t
         zdiff_sc = zdiff_sc - frc_s
         IF( lk_linssh ) THEN
            zdiff_hc1 = zdiff_hc + z_ssh_hc
            zdiff_sc1 = zdiff_sc + z_ssh_sc
            zerr_hc1  = z_ssh_hc - frc_wn_t
            zerr_sc1  = z_ssh_sc - frc_wn_s
         ENDIF

!!gm to be added ?
!      IF( ln_linssh ) THEN            ! fixed volume, add the ssh contribution
!        zvol_tot = zvol_tot + glob_2Dsum( 'diahsb', surf(:,:) * ssh(:,:,Kmm) )
!      ENDIF
!!gm end
      ELSE
         zvol_tot = -1._wp
         zdiff_v1  = 0._wp
         IF( .NOT. lk_linssh ) THEN
            zdiff_hc  = 0._wp  
            zdiff_sc  = 0._wp  
            zdiff_v2  = 0._wp  
         ELSE
            zdiff_hc1 = 0._wp  
            zdiff_sc1 = 0._wp  
            zerr_hc1  = 0._wp  
            zerr_sc1  = 0._wp  
         ENDIF
      ENDIF

      CALL iom_put(   'bgfrcvol' , frc_v    * 1.e-9    )              ! vol - surface forcing (km3)
      CALL iom_put(   'bgfrctem' , frc_t    * rho0 * rcp * 1.e-20 )   ! hc  - surface forcing (1.e20 J)
      CALL iom_put(   'bgfrchfx' , frc_t    * rho0 * rcp /  &         ! hc  - surface forcing (W/m2)
         &                       ( surf_tot * kt * rn_Dt )        )
      CALL iom_put(   'bgfrcsal' , frc_s    * 1.e-9    )              ! sc  - surface forcing (psu*km3)

      IF( .NOT. lk_linssh ) THEN
         CALL iom_put( 'bgtemper' , zdiff_hc / zvol_tot )              ! Temperature drift     (C)
         CALL iom_put( 'bgsaline' , zdiff_sc / zvol_tot )              ! Salinity    drift     (PSU)
         CALL iom_put( 'bgheatco' , zdiff_hc * 1.e-20 * rho0 * rcp )   ! Heat content drift    (1.e20 J)
         CALL iom_put( 'bgheatfx' , zdiff_hc * rho0 * rcp /  &         ! Heat flux drift       (W/m2)
            &                       ( surf_tot * kt * rn_Dt )        )
         CALL iom_put( 'bgsaltco' , zdiff_sc * 1.e-9    )              ! Salt content drift    (psu*km3)
         CALL iom_put( 'bgvolssh' , zdiff_v1 * 1.e-9    )              ! volume ssh drift      (km3)
         CALL iom_put( 'bgvole3t' , zdiff_v2 * 1.e-9    )              ! volume e3t drift      (km3)
         !
         IF( kt == nitend .AND. lwp ) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'dia_hsb : last time step hsb diagnostics: at it= ', kt-1,' date= ', ndastp
            WRITE(numout,*) '~~~~~~~'
            WRITE(numout,*) '   Temperature drift = ', zdiff_hc / zvol_tot, ' C'
            WRITE(numout,*) '   Salinity    drift = ', zdiff_sc / zvol_tot, ' PSU'
            WRITE(numout,*) '   volume ssh  drift = ', zdiff_v1 * 1.e-9   , ' km^3'
            WRITE(numout,*) '   volume e3t  drift = ', zdiff_v2 * 1.e-9   , ' km^3'
         ENDIF
         !
      ELSE
         CALL iom_put( 'bgtemper' , zdiff_hc1 / zvol_tot)              ! Heat content drift    (C)
         CALL iom_put( 'bgsaline' , zdiff_sc1 / zvol_tot)              ! Salt content drift    (PSU)
         CALL iom_put( 'bgheatco' , zdiff_hc1 * 1.e-20 * rho0 * rcp )  ! Heat content drift    (1.e20 J)
         CALL iom_put( 'bgheatfx' , zdiff_hc1 * rho0 * rcp /  &        ! Heat flux drift       (W/m2)
            &                       ( surf_tot * kt * rn_Dt )         )
         CALL iom_put( 'bgsaltco' , zdiff_sc1 * 1.e-9    )             ! Salt content drift    (psu*km3)
         CALL iom_put( 'bgvolssh' , zdiff_v1 * 1.e-9    )              ! volume ssh drift      (km3)
         CALL iom_put( 'bgmistem' , zerr_hc1 / zvol_tot )              ! hc  - error due to free surface (C)
         CALL iom_put( 'bgmissal' , zerr_sc1 / zvol_tot )              ! sc  - error due to free surface (psu)
      ENDIF
      !
      CALL dia_hsb_vol( kt, Kmm )     ! barotropic volume budget (open boundaries + AGRIF)
      !
      IF( lrst_oce )   CALL dia_hsb_rst( kt, Kmm, 'WRITE' )
      !
      IF( ln_timing )   CALL timing_stop('dia_hsb')
      !
   END SUBROUTINE dia_hsb


   SUBROUTINE dia_hsb_rst( kt, Kmm, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dia_hsb_rst  ***
      !!
      !! ** Purpose :   Read or write DIA file in restart file
      !!
      !! ** Method  :   use of IOM library
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt     ! ocean time-step
      INTEGER         , INTENT(in) ::   Kmm    ! ocean time level index
      CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      !
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      INTEGER ::   id0          ! local integer
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise

         id0 = iom_varid( numror, 'frc_v' , ldstop = .FALSE. ) ! test if this variable exists

#if defined key_RK3
         IF( ln_rstart .AND. id0 > 0 ) THEN      !* Read the restart file
#else
         IF( ln_rstart .AND. id0 > 0 .AND. (.NOT.l_1st_euler) ) THEN      !* Read the restart file
#endif
            !
            l_diahsb_init = .FALSE.
            !
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '   dia_hsb_rst : read hsb restart at it= ', kt,' date= ', ndastp
            IF(lwp) WRITE(numout,*)
            CALL iom_get( numror, 'frc_v', frc_v )
            CALL iom_get( numror, 'frc_t', frc_t )
            CALL iom_get( numror, 'frc_s', frc_s )
            IF( lk_linssh ) THEN
               CALL iom_get( numror, 'frc_wn_t', frc_wn_t )
               CALL iom_get( numror, 'frc_wn_s', frc_wn_s )
            ENDIF
            CALL iom_get( numror, jpdom_auto, 'surf_ini'  , surf_ini   ) ! ice sheet coupling
            CALL iom_get( numror, jpdom_auto, 'ssh_ini'   , ssh_ini    )
            CALL iom_get( numror, jpdom_auto, 'e3t_ini'   , e3t_ini    )
            CALL iom_get( numror, jpdom_auto, 'tmask_ini' , tmask_ini  )
            CALL iom_get( numror, jpdom_auto, 'hc_loc_ini', hc_loc_ini )
            CALL iom_get( numror, jpdom_auto, 'sc_loc_ini', sc_loc_ini )
            IF( lk_linssh ) THEN
               CALL iom_get( numror, jpdom_auto, 'ssh_hc_loc_ini', ssh_hc_loc_ini )
               CALL iom_get( numror, jpdom_auto, 'ssh_sc_loc_ini', ssh_sc_loc_ini )
            ENDIF
         ELSE
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '   dia_hsb_rst : initialise hsb at initial state '
            IF(lwp) WRITE(numout,*)
            !
            l_diahsb_init = .TRUE.
            !
            DO_2D( 0, 0, 0, 0 )
               surf_ini(ji,jj) = e1e2t(ji,jj) * tmask_i(ji,jj)         ! initial ocean surface
               ssh_ini(ji,jj) = ssh(ji,jj,Kmm)                          ! initial ssh
            END_2D
            ! if ice sheet/oceqn coupling, need to mask ini variables here (mask could change at the next NEMO instance).
            DO_3D( 0, 0, 0, 0, 1, jpk )
               e3t_ini   (ji,jj,jk) = e3t(ji,jj,jk,Kmm) * tmask(ji,jj,jk)  ! initial vertical scale factors
               tmask_ini (ji,jj,jk) = tmask(ji,jj,jk)                      ! initial mask
               hc_loc_ini(ji,jj,jk) = ts(ji,jj,jk,jp_tem,Kmm) * e3t(ji,jj,jk,Kmm) * tmask(ji,jj,jk)  ! initial heat content
               sc_loc_ini(ji,jj,jk) = ts(ji,jj,jk,jp_sal,Kmm) * e3t(ji,jj,jk,Kmm) * tmask(ji,jj,jk)  ! initial salt content
            END_3D
            frc_v = 0._wp                                           ! volume       trend due to forcing
            frc_t = 0._wp                                           ! heat content   -    -   -    -
            frc_s = 0._wp                                           ! salt content   -    -   -    -
            IF( lk_linssh ) THEN
               IF( ln_isfcav ) THEN
                  DO_2D( 0, 0, 0, 0 )
                     ssh_hc_loc_ini(ji,jj) = ts(ji,jj,mikt(ji,jj),jp_tem,Kmm) * ssh(ji,jj,Kmm)   ! initial heat content in ssh
                     ssh_sc_loc_ini(ji,jj) = ts(ji,jj,mikt(ji,jj),jp_sal,Kmm) * ssh(ji,jj,Kmm)   ! initial salt content in ssh
                  END_2D
               ELSE
                  DO_2D( 0, 0, 0, 0 )
                     ssh_hc_loc_ini(ji,jj) = ts(ji,jj,1,jp_tem,Kmm) * ssh(ji,jj,Kmm)   ! initial heat content in ssh
                     ssh_sc_loc_ini(ji,jj) = ts(ji,jj,1,jp_sal,Kmm) * ssh(ji,jj,Kmm)   ! initial salt content in ssh
                  END_2D
               END IF
               frc_wn_t = 0._wp    ! initial heat content misfit due to free surface
               frc_wn_s = 0._wp    ! initial salt content misfit due to free surface
            ENDIF
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! -------------------
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   dia_hsb_rst : write restart at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*)
         !
         CALL iom_rstput( kt, nitrst, numrow, 'frc_v', frc_v )
         CALL iom_rstput( kt, nitrst, numrow, 'frc_t', frc_t )
         CALL iom_rstput( kt, nitrst, numrow, 'frc_s', frc_s )
         IF( lk_linssh ) THEN
            CALL iom_rstput( kt, nitrst, numrow, 'frc_wn_t', frc_wn_t )
            CALL iom_rstput( kt, nitrst, numrow, 'frc_wn_s', frc_wn_s )
         ENDIF
         CALL iom_rstput( kt, nitrst, numrow, 'surf_ini'  , surf_ini   )      ! ice sheet coupling
         CALL iom_rstput( kt, nitrst, numrow, 'ssh_ini'   , ssh_ini    )
         CALL iom_rstput( kt, nitrst, numrow, 'e3t_ini'   , e3t_ini    )
         CALL iom_rstput( kt, nitrst, numrow, 'tmask_ini' , tmask_ini  )
         CALL iom_rstput( kt, nitrst, numrow, 'hc_loc_ini', hc_loc_ini )
         CALL iom_rstput( kt, nitrst, numrow, 'sc_loc_ini', sc_loc_ini )
         IF( lk_linssh ) THEN
            CALL iom_rstput( kt, nitrst, numrow, 'ssh_hc_loc_ini', ssh_hc_loc_ini )
            CALL iom_rstput( kt, nitrst, numrow, 'ssh_sc_loc_ini', ssh_sc_loc_ini )
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE dia_hsb_rst


   SUBROUTINE dia_hsb_init( Kmm )
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_hsb  ***
      !!
      !! ** Purpose: Initialization for the heat salt volume budgets
      !!
      !! ** Method : Compute initial heat content, salt content and volume
      !!
      !! ** Action : - Compute initial heat content, salt content and volume
      !!             - Initialize forcing trends
      !!             - Compute coefficients for conversion
      !!---------------------------------------------------------------------------
      INTEGER, INTENT(in) :: Kmm ! time level index
      !
      INTEGER ::   ierror, ios   ! local integer
      INTEGER ::   ji, jj        ! loop index
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dia_hsb_init : heat and salt budgets diagnostics'
         WRITE(numout,*) '~~~~~~~~~~~~ '
      ENDIF

      ! define l_diahsb
      l_diahsb = iom_use('bgfrcvol') .OR. iom_use('bgfrctem') .OR. iom_use('bgfrchfx') .OR. iom_use('bgfrcsal') .OR. &
         &       iom_use('bgtemper') .OR. iom_use('bgsaline') .OR. iom_use('bgheatco') .OR. iom_use('bgheatfx') .OR. &
         &       iom_use('bgsaltco') .OR. iom_use('bgvolssh') .OR. iom_use('bgvole3t') .OR. &
         &       iom_use('bgmistem') .OR. iom_use('bgmissal') .OR. &
         &       iom_use('bgvolres') .OR. iom_use('volres2d')
      !
      IF( .NOT. l_diahsb )   RETURN
      !
      ! ------------------- !
      ! 1 - Allocate memory !
      ! ------------------- !
      ALLOCATE( hc_loc_ini(A2D(0),jpk), sc_loc_ini(A2D(0),jpk), surf_ini(A2D(0)), &
         &      e3t_ini(A2D(0),jpk), surf(A2D(0)),  ssh_ini(A2D(0)), tmask_ini(A2D(0),jpk),STAT=ierror  )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_hsb_init: unable to allocate hc_loc_ini' )   ;   RETURN
      ENDIF

      IF( lk_linssh )   ALLOCATE( ssh_hc_loc_ini(A2D(0)), ssh_sc_loc_ini(A2D(0)),STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_hsb: unable to allocate ssh_hc_loc_ini' )   ;   RETURN
      ENDIF

      ! ----------------------------------------------- !
      ! 2 - Time independant variables and file opening !
      ! ----------------------------------------------- !
      DO_2D( 0, 0, 0, 0 )
         surf(ji,jj) = e1e2t(ji,jj) * smask0_i(ji,jj)               ! masked surface grid cell area
      END_2D
      surf_tot  = glob_2Dsum( 'diahsb', surf(:,:) )         ! total ocean surface area

      IF( ln_bdy ) CALL ctl_warn( 'dia_hsb_init: heat/salt budget does not consider open boundary fluxes' )
      !
      ! ---------------------------------- !
      ! 4 - initial conservation variables !
      ! ---------------------------------- !
      CALL dia_hsb_rst( nit000, Kmm, 'READ' )  !* read or initialize all required files
      !
      CALL dia_hsb_vol_init( Kmm )             !* barotropic volume budget
      !
   END SUBROUTINE dia_hsb_init


   SUBROUTINE dia_hsb_vol_init( Kmm )
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_hsb_vol_init  ***
      !!
      !! ** Purpose : Define the budget domain Omega and initialise the barotropic
      !!              volume budget
      !!
      !! ** Method  : Omega = wet (ssmask) T-points whose ssh is advanced by the
      !!              barotropic continuity equation, i.e. excluding
      !!              - root/BDY : points with bdytmask=0 and the T-points of
      !!                           rims 0 and 1 (overwritten by bdy_ssh)
      !!              - AGRIF child : ghost cells (fluxes/ssh set from parent)
      !!              The budget is reset at the start of each run segment
      !!              (no restart fields).
      !!---------------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kmm   ! time level index
      !
      INTEGER  ::   ji, jj, jb, ib_bdy
#if defined key_agrif
      INTEGER  ::   ii0, ii1, ij0, ij1   ! global index bounds of the child interior
#endif
      REAL(wp), DIMENSION(A2D(0),2) ::   z2d
      REAL(wp), DIMENSION(jpi,jpj)  ::   zint      ! =1 on the grid interior (child: ghost cells excluded), land included
#if defined key_agrif
      REAL(wp), DIMENSION(A2D(0),2) ::   zar       ! child areas: own (a_c) and seen through e1e2t_frac (a_c/e1e2t_frac)
      REAL(wp), DIMENSION(2)        ::   zars
      REAL(wp)                      ::   zfmin, zfmax
#endif
      CHARACTER(len=64)             ::   cl_grid   ! grid description for the file header
      !!---------------------------------------------------------------------------
      !
      ALLOCATE( vbg_msk(jpi,jpj), ssh_prev(A2D(0)), ssh_pre(A2D(0)) )
      !
      vbg_msk(:,:) = ssmask(:,:)
      zint   (:,:) = 1._wp
      !
      IF( ln_bdy ) THEN                        ! remove BDY rims 0 and 1 and the exterior
         vbg_msk(:,:) = vbg_msk(:,:) * bdytmask(:,:)
         DO ib_bdy = 1, nb_bdy
            DO jb = 1, idx_bdy(ib_bdy)%nblenrim(1)
               vbg_msk( idx_bdy(ib_bdy)%nbi(jb,1), idx_bdy(ib_bdy)%nbj(jb,1) ) = 0._wp
            END DO
         END DO
      ENDIF
      !
#if defined key_agrif
      IF( .NOT. Agrif_Root() ) THEN            ! remove AGRIF ghost cells (same bounds as Agrif_dyn_ts_flux)
         ii0 = 1   ;   ii1 = jpiglo   ;   ij0 = 1   ;   ij1 = jpjglo
         IF( lk_west  )   ii0 = nn_hls + nbghostcells + nn_shift_bar * Agrif_irhox() + 1
         IF( lk_east  )   ii1 = jpiglo - ( nn_hls + nbghostcells ) - nn_shift_bar * Agrif_irhox()
         IF( lk_south )   ij0 = nn_hls + nbghostcells + nn_shift_bar * Agrif_irhoy() + 1
         IF( lk_north )   ij1 = jpjglo - ( nn_hls + nbghostcells ) - nn_shift_bar * Agrif_irhoy()
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF(  mig(ji,nn_hls) < ii0 .OR. mig(ji,nn_hls) > ii1 .OR.   &
                  & mjg(jj,nn_hls) < ij0 .OR. mjg(jj,nn_hls) > ij1  ) THEN
                  vbg_msk(ji,jj) = 0._wp
                  zint   (ji,jj) = 0._wp
               ENDIF
            END DO
         END DO
      ENDIF
      !                                        ! parent cells updated from a child (set in updateSSH)
      IF( .NOT. ALLOCATED( ssh_upd_msk ) ) THEN
         ALLOCATE( ssh_upd_msk(jpi,jpj), ssh_upd_lost(jpi,jpj) )
         ssh_upd_msk (:,:) = 0._wp
         ssh_upd_lost(:,:) = 0._wp
      ENDIF
#endif
      CALL lbc_lnk( 'diahsb', vbg_msk, 'T', 1._wp )
      !
      DO_2D( 0, 0, 0, 0 )
         ssh_prev(ji,jj) = ssh(ji,jj,Kmm)
         ssh_pre (ji,jj) = ssh(ji,jj,Kmm)
         z2d     (ji,jj,1) = vbg_msk(ji,jj) * e1e2t(ji,jj) * ssh(ji,jj,Kmm)
         z2d     (ji,jj,2) = vbg_msk(ji,jj) * e1e2t(ji,jj)
      END_2D
      vbg_v0   = glob_2Dsum( 'diahsb', z2d(:,:,1) )
      vbg_area = glob_2Dsum( 'diahsb', z2d(:,:,2) )
#if defined key_agrif
      IF( .NOT. Agrif_Root() ) THEN            ! area weights of the child-to-parent update
         zfmin =  HUGE(1._wp)   ;   zfmax = -HUGE(1._wp)
         DO_2D( 0, 0, 0, 0 )
            zar(ji,jj,1) = zint(ji,jj) * e1e2t(ji,jj)                                        ! a_c
            zar(ji,jj,2) = zint(ji,jj) * e1e2t(ji,jj) / MAX( e1e2t_frac(ji,jj), 1.e-20_wp )  ! T_c/rho^2
            IF( zint(ji,jj) > 0._wp ) THEN
               zfmin = MIN( zfmin, e1e2t_frac(ji,jj) )   ;   zfmax = MAX( zfmax, e1e2t_frac(ji,jj) )
            ENDIF
         END_2D
         ! plain sums: glob_2Dsum ignores land points, and these areas must include them
         zars(1) = SUM( zar(:,:,1) )   ;   CALL mpp_sum( 'diahsb', zars(1) )
         zars(2) = SUM( zar(:,:,2) )   ;   CALL mpp_sum( 'diahsb', zars(2) )
         CALL mpp_min( 'diahsb', zfmin )   ;   CALL mpp_max( 'diahsb', zfmax )
      ENDIF
#endif
      vbg_frc  = 0._wp   ;   vbg_isf  = 0._wp   ;   vbg_obc = 0._wp   ;   vbg_upd = 0._wp
      vbg_frcu = 0._wp   ;   vbg_obcu = 0._wp   ;   vbg_res = 0._wp   ;   vbg_resmax = 0._wp
      !
      IF( lwm ) THEN
         cl_grid = 'single grid (no AGRIF)'
#if defined key_agrif
         IF( Agrif_Root() ) THEN   ;   cl_grid = 'parent (AGRIF root grid)'
         ELSE                      ;   cl_grid = 'AGRIF child grid '//TRIM(Agrif_CFixed())
         ENDIF
#endif
         CALL ctl_opn( numvbg, 'volume_budget.txt', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )
         WRITE(numvbg,'(a)')         '# NEMO barotropic volume budget (diahsb.F90: dia_hsb_vol, dia_hsb_post)'
         WRITE(numvbg,'(a)')         '#'
         WRITE(numvbg,'(a,a)')       '# grid                     : ', TRIM(cl_grid)
         WRITE(numvbg,'(a,ES24.16)') '# time step rn_Dt     [s]  : ', rn_Dt
         WRITE(numvbg,'(a,ES24.16)') '# A_Omega = sum(e1e2t) [m2] : ', vbg_area
         WRITE(numvbg,'(a,ES24.16)') '# V0 = sum(e1e2t*ssh)  [m3] : ', vbg_v0
#if defined key_agrif
         IF( .NOT. Agrif_Root() ) THEN
            WRITE(numvbg,'(a,ES24.16)') '# A_interior = sum(e1e2t), interior incl. land [m2] : ', zars(1)
            WRITE(numvbg,'(a,ES24.16)') '# A_frac = sum(e1e2t/e1e2t_frac), interior incl. land [m2] : ', zars(2)
            WRITE(numvbg,'(a,2ES24.16)') '# e1e2t_frac min max over interior : ', zfmin, zfmax
         ENDIF
#endif
         WRITE(numvbg,'(a)')         '#'
         WRITE(numvbg,'(a)')         '# Omega : wet T-points whose ssh obeys the barotropic continuity equation'
         WRITE(numvbg,'(a)')         '#         (parent: BDY rims 0-1 excluded; child: AGRIF ghost cells excluded)'
         WRITE(numvbg,'(a)')         '# U     : parent T-points overwritten by the child-to-parent update (empty on a child)'
         WRITE(numvbg,'(a)')         '# Budget: dV = c_surf_isf + c_lateral_in + c_agrif_update + residual'
         WRITE(numvbg,'(a)')         '# One row per time step, written after the AGRIF update.'
         WRITE(numvbg,'(a)')         '# c_* terms are cumulated from the first step of the run: sum_n rn_Dt*(...)'
         WRITE(numvbg,'(a)')         '#'
         WRITE(numvbg,'(a)')         '# col name               unit meaning'
         WRITE(numvbg,'(a)')         '#  1  kt                 -    time step index'
         WRITE(numvbg,'(a)')         '#  2  time               s    time since start of run (end of step)'
         WRITE(numvbg,'(a)')         '#  3  V                  m3   sum_Omega e1e2t*ssh'
         WRITE(numvbg,'(a)')         '#  4  dV                 m3   V - V0: volume change since start'
         WRITE(numvbg,'(a)')         '#  5  c_surf_isf         m3   surface freshwater + ISF source:'//  &
            &                        ' -sum_Omega e1e2t*ssh_frc'
         WRITE(numvbg,'(a)')         '#  6  c_isf              m3   ISF melt part of c_surf_isf:'//  &
            &                        ' sum_Omega e1e2t*(fwfisf_cav+fwfisf_par)/rho0'
         WRITE(numvbg,'(a)')         '#  7  c_lateral_in       m3   inflow through the edge of Omega'//  &
            &                        ' (parent: open boundaries; child: AGRIF interface)'
         WRITE(numvbg,'(a)')         '#  8  c_agrif_update     m3   change of V due to the child-to-parent update (0 on a child)'
         WRITE(numvbg,'(a)')         '#  9  residual           m3   dV - c_surf_isf - c_lateral_in - c_agrif_update'
         WRITE(numvbg,'(a)')         '# 10  residual_ssh       m    residual / A_Omega: equivalent sea-level error'
         WRITE(numvbg,'(a)')         '# 11  max_pointwise_res  m    max_Omega |ssh(n+1)-ssh(n)'//  &
            &                        '+rn_Dt*(ssh_frc+div(un_adv,vn_adv))| this step'
         WRITE(numvbg,'(a)')         '# 12  V_zoom             m3   sum_U e1e2t*ssh after update (parent);'//  &
            &                        ' compare with child column 3'
         WRITE(numvbg,'(a)')         '# 13  c_surf_isf_zoom    m3   as c_surf_isf, over U (parent)'
         WRITE(numvbg,'(a)')         '# 14  c_lateral_in_zoom  m3   inflow into U (parent); compare with child column 7'
         WRITE(numvbg,'(a)')         '# 15  V_zoom_land        m3   child volume on parent land points of U,'//  &
            &                        ' discarded by the update (parent); V_zoom + V_zoom_land = child column 3'
         WRITE(numvbg,'(a)')         '#'
         WRITE(numvbg,'(a)')         '#names: kt time V dV c_surf_isf c_isf c_lateral_in c_agrif_update residual'//   &
            &                        ' residual_ssh max_pointwise_res V_zoom c_surf_isf_zoom c_lateral_in_zoom V_zoom_land'
      ENDIF
      !
   END SUBROUTINE dia_hsb_vol_init


   SUBROUTINE dia_hsb_vol( kt, Kmm )
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_hsb_vol  ***
      !!
      !! ** Purpose : Barotropic volume budget over Omega, evaluated at the end of
      !!              the step and before the AGRIF update
      !!
      !! ** Method  : Discrete barotropic continuity over one step (RK3, forward):
      !!                 ssh(n+1) = ssh(n) - rn_Dt * [ ssh_frc + div( un_adv, vn_adv ) ]
      !!              with (un_adv,vn_adv) the step-averaged transports of dyn_spg_ts.
      !!              - pointwise residual (m):
      !!                 r = ssh(n+1) - ssh(n) + rn_Dt * [ ssh_frc + div ]   on Omega
      !!              - lateral inflow into Omega: sum over faces separating Omega
      !!                from its complement (BDY rim or AGRIF interface)
      !!---------------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      INTEGER, INTENT(in) ::   Kmm   ! time level index (end of step)
      !
      INTEGER  ::   ji, jj
      REAL(wp) ::   zdiv, zmax
      REAL(wp), DIMENSION(jpi,jpj)  ::   zres
      REAL(wp), DIMENSION(A2D(0),3) ::   ztmp
      REAL(wp), DIMENSION(3)        ::   zsum
      !!---------------------------------------------------------------------------
      !
      zres(:,:) = 0._wp
      DO_2D( 0, 0, 0, 0 )
         zdiv = (  ( e2u(ji  ,jj) * un_adv(ji  ,jj) - e2u(ji-1,jj  ) * un_adv(ji-1,jj  ) )   &
            &    + ( e1v(ji,jj  ) * vn_adv(ji,jj  ) - e1v(ji  ,jj-1) * vn_adv(ji  ,jj-1) )  ) * r1_e1e2t(ji,jj)
         zres(ji,jj) = vbg_msk(ji,jj) * ( ssh(ji,jj,Kmm) - ssh_prev(ji,jj) + rn_Dt * ( ssh_frc(ji,jj) + zdiv ) )
         !
         ztmp(ji,jj,1) =   vbg_msk(ji,jj) * e1e2t(ji,jj) * ssh(ji,jj,Kmm)          ! volume (ssh part)
         ztmp(ji,jj,2) = - vbg_msk(ji,jj) * e1e2t(ji,jj) * ssh_frc(ji,jj)          ! volume source (m3/s)
         ztmp(ji,jj,3) =   vbg_msk(ji,jj) * (                                    & ! lateral inflow (m3/s)
            &              ( 1._wp - vbg_msk(ji-1,jj) ) * e2u(ji-1,jj) * un_adv(ji-1,jj)   &
            &            - ( 1._wp - vbg_msk(ji+1,jj) ) * e2u(ji  ,jj) * un_adv(ji  ,jj)   &
            &            + ( 1._wp - vbg_msk(ji,jj-1) ) * e1v(ji,jj-1) * vn_adv(ji,jj-1)   &
            &            - ( 1._wp - vbg_msk(ji,jj+1) ) * e1v(ji,jj  ) * vn_adv(ji,jj  )   )
         ssh_pre(ji,jj) = ssh(ji,jj,Kmm)
      END_2D
      !
      zsum(:) = glob_2Dsum( 'dia_hsb_vol', ztmp )     ! no delayed sum: exact at every step
      zmax    = MAXVAL( ABS( zres(A2D(0)) ) )
      CALL mpp_max( 'dia_hsb_vol', zmax )
      !
      IF( ln_isf ) THEN                                ! ice-shelf melt part of the source (information)
         DO_2D( 0, 0, 0, 0 )
            ztmp(ji,jj,1) = vbg_msk(ji,jj) * e1e2t(ji,jj) * r1_rho0 * ( fwfisf_cav(ji,jj) + fwfisf_par(ji,jj) )
         END_2D
         vbg_isf = vbg_isf + rn_Dt * glob_2Dsum( 'dia_hsb_vol', ztmp(:,:,1) )
      ENDIF
      !
      vbg_frc    = vbg_frc + rn_Dt * zsum(2)
      vbg_obc    = vbg_obc + rn_Dt * zsum(3)
      vbg_res    = ( zsum(1) - vbg_v0 ) - vbg_frc - vbg_obc - vbg_upd
      vbg_resmax = zmax
      !
      CALL iom_put( 'bgvolres', vbg_res )              ! cumulated residual (m3)
      CALL iom_put( 'volres2d', zres    )              ! per-step pointwise residual (m)
#if defined key_agrif
      IF( iom_use('sshupdland') ) THEN                 ! child ssh volume per unit area discarded on parent
         zres(:,:) = ssh_upd_msk(:,:) * ssh_upd_lost(:,:)   ! land points by the last update (m)
         CALL iom_put( 'sshupdland', zres )
      ENDIF
#endif
      !
   END SUBROUTINE dia_hsb_vol


   SUBROUTINE dia_hsb_post( kt, Kmm )
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_hsb_post  ***
      !!
      !! ** Purpose : Complete the volume budget after Agrif_update_all:
      !!              - increment of sum_Omega(e1e2t*ssh) due to the child-to-parent update
      !!              - volume, source and inflow over the update region U
      !!              - child volume discarded by the update on parent land points of U
      !!              - store ssh(n+1) for the next pointwise check
      !!              - write one line of volume_budget.txt (1_volume_budget.txt on the child)
      !!---------------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      INTEGER, INTENT(in) ::   Kmm   ! time level index (end of step, after update)
      !
      INTEGER  ::   ji, jj
      REAL(wp), DIMENSION(jpi,jpj)  ::   zupd
      REAL(wp), DIMENSION(A2D(0),6) ::   ztmp
      REAL(wp), DIMENSION(6)        ::   zsum
      REAL(wp)                      ::   zau    ! area of the update region U
      !!---------------------------------------------------------------------------
      !
      zupd(:,:) = 0._wp
      ztmp(:,:,6) = 0._wp
#if defined key_agrif
      zupd(:,:) = ssh_upd_msk(:,:) * vbg_msk(:,:)
      CALL lbc_lnk( 'dia_hsb_post', zupd, 'T', 1._wp )
#endif
      !
      DO_2D( 0, 0, 0, 0 )
         ztmp(ji,jj,1) =   vbg_msk(ji,jj) * e1e2t(ji,jj) * ( ssh(ji,jj,Kmm) - ssh_pre(ji,jj) )   ! update increment
         ztmp(ji,jj,2) =   vbg_msk(ji,jj) * e1e2t(ji,jj) *   ssh(ji,jj,Kmm)                      ! volume after update
         ztmp(ji,jj,3) =      zupd(ji,jj) * e1e2t(ji,jj) *   ssh(ji,jj,Kmm)                      ! volume over U
         ztmp(ji,jj,4) = -    zupd(ji,jj) * e1e2t(ji,jj) *   ssh_frc(ji,jj)                      ! source over U
         ztmp(ji,jj,5) =      zupd(ji,jj) * (                                                  & ! inflow into U
            &              ( 1._wp - zupd(ji-1,jj) ) * e2u(ji-1,jj) * un_adv(ji-1,jj)   &
            &            - ( 1._wp - zupd(ji+1,jj) ) * e2u(ji  ,jj) * un_adv(ji  ,jj)   &
            &            + ( 1._wp - zupd(ji,jj-1) ) * e1v(ji,jj-1) * vn_adv(ji,jj-1)   &
            &            - ( 1._wp - zupd(ji,jj+1) ) * e1v(ji,jj  ) * vn_adv(ji,jj  )   )
         ssh_prev(ji,jj) = ssh(ji,jj,Kmm)
      END_2D
#if defined key_agrif
      DO_2D( 0, 0, 0, 0 )                              ! child volume discarded on parent land points of U
         ztmp(ji,jj,6) = ssh_upd_msk(ji,jj) * e1e2t(ji,jj) * ssh_upd_lost(ji,jj)
      END_2D
#endif
      !
      zsum(:) = glob_2Dsum( 'dia_hsb_post', ztmp )
      ! V_zoom_land lives on parent LAND points, which glob_2Dsum ignores: use a plain sum
      zsum(6) = SUM( ztmp(:,:,6) )   ;   CALL mpp_sum( 'dia_hsb_post', zsum(6) )
#if defined key_agrif
      IF( ll_vbg_first ) THEN                        ! area of the update region U (parent), land included
         DO_2D( 0, 0, 0, 0 )
            ztmp(ji,jj,1) = ssh_upd_msk(ji,jj) * e1e2t(ji,jj)
         END_2D
         zau = SUM( ztmp(:,:,1) )   ;   CALL mpp_sum( 'dia_hsb_post', zau )   ! plain sum: land included
         IF( lwm .AND. Agrif_Root() )   &
            &   WRITE(numvbg,'(a,ES24.16)') '# A_U = sum(e1e2t) over the update region U, incl. land [m2] : ', zau
      ENDIF
#endif
      ll_vbg_first = .FALSE.
      !
      vbg_upd  = vbg_upd  + zsum(1)
      vbg_frcu = vbg_frcu + rn_Dt * zsum(4)
      vbg_obcu = vbg_obcu + rn_Dt * zsum(5)
      vbg_res  = ( zsum(2) - vbg_v0 ) - vbg_frc - vbg_obc - vbg_upd
      !
      IF( lwm ) THEN
         WRITE(numvbg,'(I10,14(1X,ES24.16))') kt, REAL( kt - nit000 + 1, wp ) * rn_Dt,             &
            &     zsum(2), zsum(2) - vbg_v0, vbg_frc, vbg_isf, vbg_obc, vbg_upd, vbg_res,              &
            &     vbg_res / vbg_area, vbg_resmax, zsum(3), vbg_frcu, vbg_obcu, zsum(6)
         FLUSH( numvbg )
      ENDIF
      !
   END SUBROUTINE dia_hsb_post

   !!======================================================================
END MODULE diahsb
