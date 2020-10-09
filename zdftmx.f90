MODULE zdftmx
  USE oce
  USE dom_oce
  USE zdf_oce
  USE lbclnk
  USE eosbn2
  USE phycst
  USE prtctl
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE timing
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: zdf_tmx
  PUBLIC :: zdf_tmx_init
  PUBLIC :: zdf_tmx_alloc
  LOGICAL, PUBLIC, PARAMETER :: lk_zdftmx = .TRUE.
  REAL(KIND = wp) :: rn_htmx
  REAL(KIND = wp) :: rn_n2min
  REAL(KIND = wp) :: rn_tfe
  REAL(KIND = wp) :: rn_me
  LOGICAL :: ln_tmx_itf
  REAL(KIND = wp) :: rn_tfe_itf
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: en_tmx
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: mask_itf
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: az_tmx
  CONTAINS
  INTEGER FUNCTION zdf_tmx_alloc()
    ALLOCATE(en_tmx(jpi, jpj), mask_itf(jpi, jpj), az_tmx(jpi, jpj, jpk), STAT = zdf_tmx_alloc)
    IF (lk_mpp) CALL mpp_sum('zdftmx', zdf_tmx_alloc)
    IF (zdf_tmx_alloc /= 0) CALL ctl_warn('zdf_tmx_alloc: failed to allocate arrays')
  END FUNCTION zdf_tmx_alloc
  SUBROUTINE zdf_tmx(kt, p_avm, p_avt, p_avs)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: p_avm
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: p_avt, p_avs
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: ztpc
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zkz
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zav_tide
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    !$ACC KERNELS
    zav_tide(:, :, :) = MIN(60.E-4, az_tmx(:, :, :) / MAX(rn_n2min, rn2(:, :, :)))
    zkz(:, :) = 0.E0
    DO jk = 2, jpkm1
      zkz(:, :) = zkz(:, :) + e3w_n(:, :, jk) * MAX(0.E0, rn2(:, :, jk)) * rau0 * zav_tide(:, :, jk) * wmask(:, :, jk)
    END DO
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (zkz(ji, jj) /= 0.E0) zkz(ji, jj) = en_tmx(ji, jj) / zkz(ji, jj)
      END DO
    END DO
    DO jk = 2, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          zav_tide(ji, jj, jk) = zav_tide(ji, jj, jk) * MIN(zkz(ji, jj), 30. / 6.) * wmask(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    IF (kt == nit000) THEN
      !$ACC KERNELS
      ztpc = 0.E0
      DO jk = 1, jpk
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            ztpc = ztpc + e3w_n(ji, jj, jk) * e1t(ji, jj) * e2t(ji, jj) * MAX(0.E0, rn2(ji, jj, jk)) * zav_tide(ji, jj, jk) * &
&tmask(ji, jj, jk) * tmask_i(ji, jj)
          END DO
        END DO
      END DO
      ztpc = rau0 / (rn_tfe * rn_me) * ztpc
      !$ACC END KERNELS
      CALL profile_psy_data0 % PreStart('zdf_tmx', 'r0', 0, 0)
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) '          N Total power consumption by av_tide    : ztpc = ', ztpc * 1.E-12, 'TW'
      CALL profile_psy_data0 % PostEnd
    END IF
    IF (ln_tmx_itf) CALL tmx_itf(kt, zav_tide)
    !$ACC KERNELS
    DO jk = 2, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          p_avt(ji, jj, jk) = p_avt(ji, jj, jk) + zav_tide(ji, jj, jk) * wmask(ji, jj, jk)
          p_avs(ji, jj, jk) = p_avs(ji, jj, jk) + zav_tide(ji, jj, jk) * wmask(ji, jj, jk)
          p_avm(ji, jj, jk) = p_avm(ji, jj, jk) + zav_tide(ji, jj, jk) * wmask(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('zdf_tmx', 'r1', 0, 0)
    CALL iom_put("av_tmx", zav_tide)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = zav_tide, clinfo1 = ' tmx - av_tide: ', tab3d_2 = p_avt, clinfo2 = ' p_avt: ', kdim = jpk)
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE zdf_tmx
  SUBROUTINE tmx_itf(kt, pav)
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), INTENT(INOUT), DIMENSION(jpi, jpj, jpk) :: pav
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zcoef, ztpc
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zkz
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zsum1, zsum2, zsum
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zempba_3d_1, zempba_3d_2
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zempba_3d, zdn2dz
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zavt_itf
    !$ACC KERNELS
    zdn2dz(:, :, jpk) = 0.E0
    zempba_3d_1(:, :, jpk) = 0.E0
    zempba_3d_2(:, :, jpk) = 0.E0
    !$ACC END KERNELS
    DO jk = 1, jpkm1
      !$ACC KERNELS
      zdn2dz(:, :, jk) = rn2(:, :, jk) - rn2(:, :, jk + 1)
      zempba_3d_1(:, :, jk) = SQRT(MAX(0.E0, rn2(:, :, jk)))
      zempba_3d_2(:, :, jk) = MAX(0.E0, rn2(:, :, jk))
      !$ACC END KERNELS
    END DO
    !$ACC KERNELS
    zsum(:, :) = 0.E0
    zsum1(:, :) = 0.E0
    zsum2(:, :) = 0.E0
    !$ACC END KERNELS
    DO jk = 2, jpk
      !$ACC KERNELS
      zsum1(:, :) = zsum1(:, :) + zempba_3d_1(:, :, jk) * e3w_n(:, :, jk) * tmask(:, :, jk) * tmask(:, :, jk - 1)
      zsum2(:, :) = zsum2(:, :) + zempba_3d_2(:, :, jk) * e3w_n(:, :, jk) * tmask(:, :, jk) * tmask(:, :, jk - 1)
      !$ACC END KERNELS
    END DO
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (zsum1(ji, jj) /= 0.E0) zsum1(ji, jj) = 1.E0 / zsum1(ji, jj)
        IF (zsum2(ji, jj) /= 0.E0) zsum2(ji, jj) = 1.E0 / zsum2(ji, jj)
      END DO
    END DO
    DO jk = 1, jpk
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          zcoef = 0.5 - SIGN(0.5, zdn2dz(ji, jj, jk))
          ztpc = zempba_3d_1(ji, jj, jk) * zsum1(ji, jj) * zcoef + zempba_3d_2(ji, jj, jk) * zsum2(ji, jj) * (1. - zcoef)
          zempba_3d(ji, jj, jk) = ztpc
          zsum(ji, jj) = zsum(ji, jj) + ztpc * e3w_n(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (zsum(ji, jj) > 0.E0) zsum(ji, jj) = 1.E0 / zsum(ji, jj)
      END DO
    END DO
    zcoef = rn_tfe_itf / (rn_tfe * rau0)
    DO jk = 1, jpk
      zavt_itf(:, :, jk) = MIN(10.E-4, zcoef * en_tmx(:, :) * zsum(:, :) * zempba_3d(:, :, jk) / MAX(rn_n2min, rn2(:, :, jk)) * &
&tmask(:, :, jk))
    END DO
    zkz(:, :) = 0.E0
    DO jk = 2, jpkm1
      zkz(:, :) = zkz(:, :) + e3w_n(:, :, jk) * MAX(0.E0, rn2(:, :, jk)) * rau0 * zavt_itf(:, :, jk) * tmask(:, :, jk) * tmask(:, &
&:, jk - 1)
    END DO
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (zkz(ji, jj) /= 0.E0) zkz(ji, jj) = en_tmx(ji, jj) * rn_tfe_itf / rn_tfe / zkz(ji, jj)
      END DO
    END DO
    DO jk = 2, jpkm1
      zavt_itf(:, :, jk) = zavt_itf(:, :, jk) * MIN(zkz(:, :), 120. / 10.) * tmask(:, :, jk) * tmask(:, :, jk - 1)
    END DO
    !$ACC END KERNELS
    IF (kt == nit000) THEN
      !$ACC KERNELS
      ztpc = 0.E0
      DO jk = 1, jpk
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            ztpc = ztpc + e1t(ji, jj) * e2t(ji, jj) * e3w_n(ji, jj, jk) * MAX(0.E0, rn2(ji, jj, jk)) * zavt_itf(ji, jj, jk) * &
&tmask(ji, jj, jk) * tmask_i(ji, jj)
          END DO
        END DO
      END DO
      ztpc = rau0 * ztpc / (rn_me * rn_tfe_itf)
      !$ACC END KERNELS
      IF (lwp) WRITE(numout, FMT = *) '          N Total power consumption by zavt_itf: ztpc = ', ztpc * 1.E-12, 'TW'
    END IF
    !$ACC KERNELS
    DO jk = 2, jpkm1
      pav(:, :, jk) = pav(:, :, jk) * (1.E0 - mask_itf(:, :)) + zavt_itf(:, :, jk) * mask_itf(:, :)
    END DO
    !$ACC END KERNELS
  END SUBROUTINE tmx_itf
  SUBROUTINE zdf_tmx_init
    INTEGER :: ji, jj, jk
    INTEGER :: inum
    INTEGER :: ios
    REAL(KIND = wp) :: ztpc, ze_z
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zem2, zek1
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zkz
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zfact
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zhdep
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zpc
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zav_tide
    NAMELIST /namzdf_tmx/ rn_htmx, rn_n2min, rn_tfe, rn_me, ln_tmx_itf, rn_tfe_itf
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namzdf_tmx, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namzdf_tmx in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namzdf_tmx, IOSTAT = ios, ERR = 902)
902 IF (ios /= 0) CALL ctl_nam(ios, 'namzdf_tmx in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namzdf_tmx)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'zdf_tmx_init : tidal mixing'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namzdf_tmx : set tidal mixing parameters'
      WRITE(numout, FMT = *) '      Vertical decay scale for turbulence   = ', rn_htmx
      WRITE(numout, FMT = *) '      Brunt-Vaisala frequency threshold     = ', rn_n2min
      WRITE(numout, FMT = *) '      Tidal dissipation efficiency          = ', rn_tfe
      WRITE(numout, FMT = *) '      Mixing efficiency                     = ', rn_me
      WRITE(numout, FMT = *) '      ITF specific parameterisation         = ', ln_tmx_itf
      WRITE(numout, FMT = *) '      ITF tidal dissipation efficiency      = ', rn_tfe_itf
    END IF
    IF (zdf_tmx_alloc() /= 0) CALL ctl_stop('STOP', 'zdf_tmx_init : unable to allocate tmx arrays')
    IF (ln_tmx_itf) THEN
      CALL iom_open('mask_itf', inum)
      CALL iom_get(inum, jpdom_data, 'tmaskitf', mask_itf, 1)
      CALL iom_close(inum)
    END IF
    CALL iom_open('M2rowdrg', inum)
    CALL iom_get(inum, jpdom_data, 'field', zem2, 1)
    CALL iom_close(inum)
    CALL iom_open('K1rowdrg', inum)
    CALL iom_get(inum, jpdom_data, 'field', zek1, 1)
    CALL iom_close(inum)
    !$ACC KERNELS
    en_tmx(:, :) = - rn_tfe * rn_me * (zem2(:, :) * 1.25 + zek1(:, :)) * ssmask(:, :)
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zhdep(ji, jj) = gdepw_0(ji, jj, mbkt(ji, jj) + 1)
        zfact(ji, jj) = rau0 * rn_htmx * (1. - EXP(- zhdep(ji, jj) / rn_htmx))
        IF (zfact(ji, jj) /= 0) zfact(ji, jj) = en_tmx(ji, jj) / zfact(ji, jj)
      END DO
    END DO
    DO jk = 1, jpk
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          az_tmx(ji, jj, jk) = zfact(ji, jj) * EXP(- (zhdep(ji, jj) - gdepw_0(ji, jj, jk)) / rn_htmx) * tmask(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    IF (nprint == 1 .AND. lwp) THEN
      !$ACC KERNELS
      zav_tide(:, :, :) = 0.E0
      DO jk = 2, jpkm1
        zav_tide(:, :, jk) = az_tmx(:, :, jk) / MAX(rn_n2min, rn2(:, :, jk))
      END DO
      ztpc = 0.E0
      zpc(:, :, :) = MAX(rn_n2min, rn2(:, :, :)) * zav_tide(:, :, :)
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            ztpc = ztpc + e3w_0(ji, jj, jk) * e1t(ji, jj) * e2t(ji, jj) * zpc(ji, jj, jk) * wmask(ji, jj, jk) * tmask_i(ji, jj)
          END DO
        END DO
      END DO
      ztpc = rau0 * 1 / (rn_tfe * rn_me) * ztpc
      !$ACC END KERNELS
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) '          Total power consumption of the tidally driven part of Kz : ztpc = ', ztpc * 1.E-12, 'TW'
      !$ACC KERNELS
      zav_tide(:, :, :) = MIN(zav_tide(:, :, :), 60.E-4)
      zkz(:, :) = 0.E0
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            zkz(ji, jj) = zkz(ji, jj) + e3w_0(ji, jj, jk) * MAX(0.E0, rn2(ji, jj, jk)) * rau0 * zav_tide(ji, jj, jk) * wmask(ji, &
&jj, jk)
          END DO
        END DO
      END DO
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (zkz(ji, jj) /= 0.E0) THEN
            zkz(ji, jj) = en_tmx(ji, jj) / zkz(ji, jj)
          END IF
        END DO
      END DO
      ztpc = 1.E50
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (zkz(ji, jj) /= 0.E0) THEN
            ztpc = MIN(zkz(ji, jj), ztpc)
          END IF
        END DO
      END DO
      !$ACC END KERNELS
      WRITE(numout, FMT = *) '          Min de zkz ', ztpc, ' Max = ', MAXVAL(zkz(:, :))
      !$ACC KERNELS
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            zav_tide(ji, jj, jk) = zav_tide(ji, jj, jk) * MIN(zkz(ji, jj), 30. / 6.) * wmask(ji, jj, jk)
          END DO
        END DO
      END DO
      ztpc = 0.E0
      zpc(:, :, :) = MAX(0.E0, rn2(:, :, :)) * zav_tide(:, :, :)
      DO jk = 1, jpk
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            ztpc = ztpc + e3w_0(ji, jj, jk) * e1t(ji, jj) * e2t(ji, jj) * zpc(ji, jj, jk) * wmask(ji, jj, jk) * tmask_i(ji, jj)
          END DO
        END DO
      END DO
      ztpc = rau0 * 1 / (rn_tfe * rn_me) * ztpc
      !$ACC END KERNELS
      WRITE(numout, FMT = *) '          2 Total power consumption of the tidally driven part of Kz : ztpc = ', ztpc * 1.E-12, 'TW'
      DO jk = 1, jpk
        !$ACC KERNELS
        ze_z = SUM(e1t(:, :) * e2t(:, :) * zav_tide(:, :, jk) * tmask_i(:, :)) / MAX(1.E-20, SUM(e1t(:, :) * e2t(:, :) * wmask(:, &
&:, jk) * tmask_i(:, :)))
        ztpc = 1.E50
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            IF (zav_tide(ji, jj, jk) /= 0.E0) ztpc = MIN(ztpc, zav_tide(ji, jj, jk))
          END DO
        END DO
        !$ACC END KERNELS
        WRITE(numout, FMT = *) '            N2 min - jk= ', jk, '   ', ze_z * 1.E4, ' cm2/s min= ', ztpc * 1.E4, 'max= ', &
&MAXVAL(zav_tide(:, :, jk)) * 1.E4, ' cm2/s'
      END DO
      WRITE(numout, FMT = *) '          e_tide : ', SUM(e1t * e2t * en_tmx) / (rn_tfe * rn_me) * 1.E-12, 'TW'
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) '          Initial profile of tidal vertical mixing'
      DO jk = 1, jpk
        !$ACC KERNELS
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            zkz(ji, jj) = az_tmx(ji, jj, jk) / MAX(rn_n2min, rn2(ji, jj, jk))
          END DO
        END DO
        ze_z = SUM(e1t(:, :) * e2t(:, :) * zkz(:, :) * tmask_i(:, :)) / MAX(1.E-20, SUM(e1t(:, :) * e2t(:, :) * wmask(:, :, jk) * &
&tmask_i(:, :)))
        !$ACC END KERNELS
        WRITE(numout, FMT = *) '                jk= ', jk, '   ', ze_z * 1.E4, ' cm2/s'
      END DO
      DO jk = 1, jpk
        !$ACC KERNELS
        zkz(:, :) = az_tmx(:, :, jk) / rn_n2min
        ze_z = SUM(e1t(:, :) * e2t(:, :) * zkz(:, :) * tmask_i(:, :)) / MAX(1.E-20, SUM(e1t(:, :) * e2t(:, :) * wmask(:, :, jk) * &
&tmask_i(:, :)))
        !$ACC END KERNELS
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) '          N2 min - jk= ', jk, '   ', ze_z * 1.E4, ' cm2/s min= ', MINVAL(zkz) * 1.E4, 'max= ', &
&MAXVAL(zkz) * 1.E4, ' cm2/s'
      END DO
    END IF
  END SUBROUTINE zdf_tmx_init
END MODULE zdftmx