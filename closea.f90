MODULE closea
  USE oce
  USE dom_oce
  USE phycst
  USE sbc_oce
  USE iom
  USE in_out_manager
  USE lib_fortran, ONLY: glob_sum
  USE lbclnk
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dom_clo
  PUBLIC :: sbc_clo
  PUBLIC :: clo_rnf
  PUBLIC :: clo_bat
  LOGICAL, PUBLIC :: ln_closea
  LOGICAL, PUBLIC :: l_sbc_clo
  LOGICAL, PUBLIC :: l_clo_rnf
  INTEGER, PUBLIC :: jncs
  INTEGER, PUBLIC :: jncsr
  INTEGER, PUBLIC :: jncse
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: closea_mask
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: closea_mask_rnf
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: closea_mask_empmr
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: surf
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: surfr
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: surfe
  CONTAINS
  SUBROUTINE dom_clo
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: inum
    INTEGER :: ierr
    INTEGER :: id
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdata_in
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data7
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data8
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data9
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data10
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data11
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data12
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data13
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data14
    CALL profile_psy_data0 % PreStart('dom_clo', 'r0', 0, 0)
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) 'dom_clo : read in masks to define closed seas '
    IF (lwp) WRITE(numout, FMT = *) '~~~~~~~'
    CALL profile_psy_data0 % PostEnd
    IF (ln_read_cfg) THEN
      CALL profile_psy_data1 % PreStart('dom_clo', 'r1', 0, 0)
      CALL iom_open(cn_domcfg, inum)
      id = iom_varid(inum, 'closea_mask', ldstop = .FALSE.)
      CALL profile_psy_data1 % PostEnd
      IF (id > 0) THEN
        CALL profile_psy_data2 % PreStart('dom_clo', 'r2', 0, 0)
        l_sbc_clo = .TRUE.
        ALLOCATE(closea_mask(jpi, jpj), STAT = ierr)
        IF (ierr /= 0) CALL ctl_stop('STOP', 'dom_clo: failed to allocate closea_mask array')
        CALL profile_psy_data2 % PostEnd
        !$ACC KERNELS
        zdata_in(:, :) = 0.0
        !$ACC END KERNELS
        CALL iom_get(inum, jpdom_data, 'closea_mask', zdata_in)
        !$ACC KERNELS
        closea_mask(:, :) = NINT(zdata_in(:, :)) * tmask(:, :, 1)
        jncs = MAXVAL(closea_mask(:, :))
        !$ACC END KERNELS
        CALL profile_psy_data3 % PreStart('dom_clo', 'r3', 0, 0)
        CALL mpp_max('closea', jncs)
        IF (jncs > 0) THEN
          IF (lwp) WRITE(numout, FMT = *) 'Number of closed seas : ', jncs
        ELSE
          CALL ctl_stop('Problem with closea_mask field in domain_cfg file. Has no values > 0 so no closed seas defined.')
        END IF
        CALL profile_psy_data3 % PostEnd
      ELSE
        CALL profile_psy_data4 % PreStart('dom_clo', 'r4', 0, 0)
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   closea_mask field not found in domain_cfg file.'
        IF (lwp) WRITE(numout, FMT = *) '           No closed seas defined.'
        IF (lwp) WRITE(numout, FMT = *)
        l_sbc_clo = .FALSE.
        jncs = 0
        CALL profile_psy_data4 % PostEnd
      END IF
      CALL profile_psy_data5 % PreStart('dom_clo', 'r5', 0, 0)
      l_clo_rnf = .FALSE.
      CALL profile_psy_data5 % PostEnd
      IF (l_sbc_clo) THEN
        CALL profile_psy_data6 % PreStart('dom_clo', 'r6', 0, 0)
        id = iom_varid(inum, 'closea_mask_rnf', ldstop = .FALSE.)
        CALL profile_psy_data6 % PostEnd
        IF (id > 0) THEN
          CALL profile_psy_data7 % PreStart('dom_clo', 'r7', 0, 0)
          l_clo_rnf = .TRUE.
          ALLOCATE(closea_mask_rnf(jpi, jpj), STAT = ierr)
          IF (ierr /= 0) CALL ctl_stop('STOP', 'dom_clo: failed to allocate closea_mask_rnf array')
          CALL iom_get(inum, jpdom_data, 'closea_mask_rnf', zdata_in)
          CALL profile_psy_data7 % PostEnd
          !$ACC KERNELS
          closea_mask_rnf(:, :) = NINT(zdata_in(:, :)) * tmask(:, :, 1)
          jncsr = MAXVAL(closea_mask_rnf(:, :))
          !$ACC END KERNELS
          CALL profile_psy_data8 % PreStart('dom_clo', 'r8', 0, 0)
          CALL mpp_max('closea', jncsr)
          IF (jncsr > 0) THEN
            IF (lwp) WRITE(numout, FMT = *) 'Number of closed seas rnf mappings : ', jncsr
          ELSE
            CALL ctl_stop('Problem with closea_mask_rnf field in domain_cfg file. Has no values > 0 so no closed seas rnf mappings &
&defined.')
          END IF
          CALL profile_psy_data8 % PostEnd
        ELSE
          CALL profile_psy_data9 % PreStart('dom_clo', 'r9', 0, 0)
          IF (lwp) WRITE(numout, FMT = *) 'closea_mask_rnf field not found in domain_cfg file. No closed seas rnf mappings defined.'
          jncsr = 0
          CALL profile_psy_data9 % PostEnd
        END IF
        CALL profile_psy_data10 % PreStart('dom_clo', 'r10', 0, 0)
        id = iom_varid(inum, 'closea_mask_empmr', ldstop = .FALSE.)
        CALL profile_psy_data10 % PostEnd
        IF (id > 0) THEN
          CALL profile_psy_data11 % PreStart('dom_clo', 'r11', 0, 0)
          l_clo_rnf = .TRUE.
          ALLOCATE(closea_mask_empmr(jpi, jpj), STAT = ierr)
          IF (ierr /= 0) CALL ctl_stop('STOP', 'dom_clo: failed to allocate closea_mask_empmr array')
          CALL iom_get(inum, jpdom_data, 'closea_mask_empmr', zdata_in)
          CALL profile_psy_data11 % PostEnd
          !$ACC KERNELS
          closea_mask_empmr(:, :) = NINT(zdata_in(:, :)) * tmask(:, :, 1)
          jncse = MAXVAL(closea_mask_empmr(:, :))
          !$ACC END KERNELS
          CALL profile_psy_data12 % PreStart('dom_clo', 'r12', 0, 0)
          CALL mpp_max('closea', jncse)
          IF (jncse > 0) THEN
            IF (lwp) WRITE(numout, FMT = *) 'Number of closed seas empmr mappings : ', jncse
          ELSE
            CALL ctl_stop('Problem with closea_mask_empmr field in domain_cfg file. Has no values > 0 so no closed seas empmr &
&mappings defined.')
          END IF
          CALL profile_psy_data12 % PostEnd
        ELSE
          CALL profile_psy_data13 % PreStart('dom_clo', 'r13', 0, 0)
          IF (lwp) WRITE(numout, FMT = *) 'closea_mask_empmr field not found in domain_cfg file. No closed seas empmr mappings &
&defined.'
          jncse = 0
          CALL profile_psy_data13 % PostEnd
        END IF
      END IF
      CALL iom_close(inum)
    ELSE
      CALL profile_psy_data14 % PreStart('dom_clo', 'r14', 0, 0)
      IF (lwp) WRITE(numout, FMT = *) 'No domain_cfg file so no closed seas defined.'
      l_sbc_clo = .FALSE.
      l_clo_rnf = .FALSE.
      CALL profile_psy_data14 % PostEnd
    END IF
  END SUBROUTINE dom_clo
  SUBROUTINE sbc_clo(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ierr
    INTEGER :: jc, jcr, jce
    REAL(KIND = wp), PARAMETER :: rsmall = 1.E-20_wp
    REAL(KIND = wp) :: zfwf_total, zcoef, zcoef1
    REAL(KIND = wp), DIMENSION(jncs) :: zfwf
    REAL(KIND = wp), DIMENSION(jncsr + 1) :: zfwfr
    REAL(KIND = wp), DIMENSION(jncse + 1) :: zfwfe
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztmp2d
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data7
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data8
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data9
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data10
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data11
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data12
    IF (ln_timing) CALL timing_start('sbc_clo')
    IF (kt == nit000) THEN
      CALL profile_psy_data0 % PreStart('sbc_clo', 'r0', 0, 0)
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'sbc_clo : closed seas '
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~'
      ALLOCATE(surf(jncs + 1), STAT = ierr)
      IF (ierr /= 0) CALL ctl_stop('STOP', 'sbc_clo: failed to allocate surf array')
      CALL profile_psy_data0 % PostEnd
      !$ACC KERNELS
      surf(:) = 0.E0_wp
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('sbc_clo', 'r1', 0, 0)
      ALLOCATE(surfr(jncsr + 1), STAT = ierr)
      IF (ierr /= 0) CALL ctl_stop('STOP', 'sbc_clo: failed to allocate surfr array')
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      surfr(:) = 0.E0_wp
      !$ACC END KERNELS
      CALL profile_psy_data2 % PreStart('sbc_clo', 'r2', 0, 0)
      ALLOCATE(surfe(jncse + 1), STAT = ierr)
      IF (ierr /= 0) CALL ctl_stop('STOP', 'sbc_clo: failed to allocate surfe array')
      CALL profile_psy_data2 % PostEnd
      !$ACC KERNELS
      surfe(:) = 0.E0_wp
      !$ACC END KERNELS
      CALL profile_psy_data3 % PreStart('sbc_clo', 'r3', 0, 0)
      surf(jncs + 1) = glob_sum('closea', e1e2t(:, :))
      CALL profile_psy_data3 % PostEnd
      DO jc = 1, jncs
        !$ACC KERNELS
        ztmp2d(:, :) = 0.E0_wp
        WHERE (closea_mask(:, :) == jc) ztmp2d(:, :) = e1e2t(:, :) * tmask_i(:, :)
        !$ACC END KERNELS
        CALL profile_psy_data4 % PreStart('sbc_clo', 'r4', 0, 0)
        surf(jc) = glob_sum('closea', ztmp2d(:, :))
        CALL profile_psy_data4 % PostEnd
      END DO
      CALL profile_psy_data5 % PreStart('sbc_clo', 'r5', 0, 0)
      surf(jncs + 1) = surf(jncs + 1) - SUM(surf(1 : jncs))
      CALL profile_psy_data5 % PostEnd
      IF (jncsr > 0) THEN
        DO jcr = 1, jncsr
          !$ACC KERNELS
          ztmp2d(:, :) = 0.E0_wp
          WHERE (closea_mask_rnf(:, :) == jcr .AND. closea_mask(:, :) == 0) ztmp2d(:, :) = e1e2t(:, :) * tmask_i(:, :)
          !$ACC END KERNELS
          CALL profile_psy_data6 % PreStart('sbc_clo', 'r6', 0, 0)
          surfr(jcr) = glob_sum('closea', ztmp2d(:, :))
          CALL profile_psy_data6 % PostEnd
        END DO
      END IF
      IF (jncse > 0) THEN
        DO jce = 1, jncse
          !$ACC KERNELS
          ztmp2d(:, :) = 0.E0_wp
          WHERE (closea_mask_empmr(:, :) == jce .AND. closea_mask(:, :) == 0) ztmp2d(:, :) = e1e2t(:, :) * tmask_i(:, :)
          !$ACC END KERNELS
          CALL profile_psy_data7 % PreStart('sbc_clo', 'r7', 0, 0)
          surfe(jce) = glob_sum('closea', ztmp2d(:, :))
          CALL profile_psy_data7 % PostEnd
        END DO
      END IF
      CALL profile_psy_data8 % PreStart('sbc_clo', 'r8', 0, 0)
      IF (lwp) WRITE(numout, FMT = *) '     Closed sea surface areas (km2)'
      DO jc = 1, jncs
        IF (lwp) WRITE(numout, FMT = '(1I3,5X,ES12.2)') jc, surf(jc) * 1.0E-6
      END DO
      IF (lwp) WRITE(numout, FMT = '(A,ES12.2)') 'Global surface area excluding closed seas (km2): ', surf(jncs + 1) * 1.0E-6
      IF (jncsr > 0) THEN
        IF (lwp) WRITE(numout, FMT = *) '     Closed sea target rnf surface areas (km2)'
        DO jcr = 1, jncsr
          IF (lwp) WRITE(numout, FMT = '(1I3,5X,ES12.2)') jcr, surfr(jcr) * 1.0E-6
        END DO
      END IF
      IF (jncse > 0) THEN
        IF (lwp) WRITE(numout, FMT = *) '     Closed sea target empmr surface areas (km2)'
        DO jce = 1, jncse
          IF (lwp) WRITE(numout, FMT = '(1I3,5X,ES12.2)') jce, surfe(jce) * 1.0E-6
        END DO
      END IF
      CALL profile_psy_data8 % PostEnd
    END IF
    !$ACC KERNELS
    zfwf_total = 0._wp
    zfwf(:) = 0.E0_wp
    !$ACC END KERNELS
    DO jc = 1, jncs
      !$ACC KERNELS
      ztmp2d(:, :) = 0.E0_wp
      WHERE (closea_mask(:, :) == jc) ztmp2d(:, :) = e1e2t(:, :) * (emp(:, :) - rnf(:, :)) * tmask_i(:, :)
      !$ACC END KERNELS
      CALL profile_psy_data9 % PreStart('sbc_clo', 'r9', 0, 0)
      zfwf(jc) = glob_sum('closea', ztmp2d(:, :))
      CALL profile_psy_data9 % PostEnd
    END DO
    !$ACC KERNELS
    zfwf_total = SUM(zfwf)
    zfwfr(:) = 0.E0_wp
    !$ACC END KERNELS
    IF (jncsr > 0) THEN
      DO jcr = 1, jncsr
        !$ACC KERNELS
        ztmp2d(:, :) = 0.E0_wp
        WHERE (closea_mask_rnf(:, :) == jcr .AND. closea_mask(:, :) > 0) ztmp2d(:, :) = e1e2t(:, :) * (emp(:, :) - rnf(:, :)) * &
&tmask_i(:, :)
        !$ACC END KERNELS
        CALL profile_psy_data10 % PreStart('sbc_clo', 'r10', 0, 0)
        zfwfr(jcr) = glob_sum('closea', ztmp2d(:, :))
        CALL profile_psy_data10 % PostEnd
        IF (ABS(zfwfr(jcr) / surf(jncs + 1)) > rsmall) THEN
          IF (zfwfr(jcr) < 0.0) THEN
            !$ACC KERNELS
            zfwf_total = zfwf_total - zfwfr(jcr)
            zcoef = zfwfr(jcr) / surfr(jcr)
            zcoef1 = rcp * zcoef
            WHERE (closea_mask_rnf(:, :) == jcr .AND. closea_mask(:, :) == 0.0)
              emp(:, :) = emp(:, :) + zcoef
              qns(:, :) = qns(:, :) - zcoef1 * sst_m(:, :)
            END WHERE
            !$ACC END KERNELS
          END IF
        END IF
      END DO
    END IF
    !$ACC KERNELS
    zfwfe(:) = 0.E0_wp
    !$ACC END KERNELS
    IF (jncse > 0) THEN
      DO jce = 1, jncse
        !$ACC KERNELS
        ztmp2d(:, :) = 0.E0_wp
        WHERE (closea_mask_empmr(:, :) == jce .AND. closea_mask(:, :) > 0) ztmp2d(:, :) = e1e2t(:, :) * (emp(:, :) - rnf(:, :)) * &
&tmask_i(:, :)
        !$ACC END KERNELS
        CALL profile_psy_data11 % PreStart('sbc_clo', 'r11', 0, 0)
        zfwfe(jce) = glob_sum('closea', ztmp2d(:, :))
        CALL profile_psy_data11 % PostEnd
        IF (ABS(zfwfe(jce) / surf(jncs + 1)) > rsmall) THEN
          !$ACC KERNELS
          zfwf_total = zfwf_total - zfwfe(jce)
          zcoef = zfwfe(jce) / surfe(jce)
          zcoef1 = rcp * zcoef
          WHERE (closea_mask_empmr(:, :) == jce .AND. closea_mask(:, :) == 0.0)
            emp(:, :) = emp(:, :) + zcoef
            qns(:, :) = qns(:, :) - zcoef1 * sst_m(:, :)
          END WHERE
          !$ACC END KERNELS
        END IF
      END DO
    END IF
    IF (ABS(zfwf_total / surf(jncs + 1)) > rsmall) THEN
      !$ACC KERNELS
      zcoef = zfwf_total / surf(jncs + 1)
      zcoef1 = rcp * zcoef
      WHERE (closea_mask(:, :) == 0)
        emp(:, :) = emp(:, :) + zcoef
        qns(:, :) = qns(:, :) - zcoef1 * sst_m(:, :)
      END WHERE
      !$ACC END KERNELS
    END IF
    DO jc = 1, jncs
      IF (ABS(zfwf(jc) / surf(jncs + 1)) > rsmall) THEN
        !$ACC KERNELS
        zcoef = zfwf(jc) / surf(jc)
        zcoef1 = rcp * zcoef
        WHERE (closea_mask(:, :) == jc)
          emp(:, :) = emp(:, :) - zcoef
          qns(:, :) = qns(:, :) + zcoef1 * sst_m(:, :)
        END WHERE
        !$ACC END KERNELS
      END IF
    END DO
    !$ACC KERNELS
    emp(:, :) = emp(:, :) * tmask(:, :, 1)
    !$ACC END KERNELS
    CALL profile_psy_data12 % PreStart('sbc_clo', 'r12', 0, 0)
    CALL lbc_lnk('closea', emp, 'T', 1._wp)
    IF (ln_timing) CALL timing_stop('sbc_clo')
    CALL profile_psy_data12 % PostEnd
  END SUBROUTINE sbc_clo
  SUBROUTINE clo_rnf(p_rnfmsk)
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: p_rnfmsk
    !$ACC KERNELS
    IF (jncsr > 0) THEN
      WHERE (closea_mask_rnf(:, :) > 0 .AND. closea_mask(:, :) == 0)
        p_rnfmsk(:, :) = MAX(p_rnfmsk(:, :), 1.0_wp)
      END WHERE
    END IF
    IF (jncse > 0) THEN
      WHERE (closea_mask_empmr(:, :) > 0 .AND. closea_mask(:, :) == 0)
        p_rnfmsk(:, :) = MAX(p_rnfmsk(:, :), 1.0_wp)
      END WHERE
    END IF
    !$ACC END KERNELS
  END SUBROUTINE clo_rnf
  SUBROUTINE clo_bat(k_top, k_bot)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, DIMENSION(:, :), INTENT(INOUT) :: k_top, k_bot
    INTEGER :: inum, id
    INTEGER, DIMENSION(jpi, jpj) :: closea_mask
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdata_in
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('clo_bat', 'r0', 0, 0)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'clo_bat : suppression of closed seas'
      WRITE(numout, FMT = *) '~~~~~~~'
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (ln_read_cfg) THEN
      CALL profile_psy_data1 % PreStart('clo_bat', 'r1', 0, 0)
      CALL iom_open(cn_domcfg, inum)
      id = iom_varid(inum, 'closea_mask', ldstop = .FALSE.)
      CALL profile_psy_data1 % PostEnd
      IF (id > 0) THEN
        CALL profile_psy_data2 % PreStart('clo_bat', 'r2', 0, 0)
        IF (lwp) WRITE(numout, FMT = *) 'Suppressing closed seas in bathymetry based on closea_mask field,'
        CALL iom_get(inum, jpdom_data, 'closea_mask', zdata_in)
        CALL profile_psy_data2 % PostEnd
        !$ACC KERNELS
        closea_mask(:, :) = NINT(zdata_in(:, :))
        WHERE (closea_mask(:, :) > 0)
          k_top(:, :) = 0
          k_bot(:, :) = 0
        END WHERE
        !$ACC END KERNELS
      ELSE
        IF (lwp) WRITE(numout, FMT = *) 'No closea_mask field found in domain_cfg file. No suppression of closed seas.'
      END IF
      CALL iom_close(inum)
    ELSE
      IF (lwp) WRITE(numout, FMT = *) 'No domain_cfg file => no suppression of closed seas.'
    END IF
    CALL profile_psy_data3 % PreStart('clo_bat', 'r3', 0, 0)
    l_sbc_clo = .FALSE.
    l_clo_rnf = .FALSE.
    CALL profile_psy_data3 % PostEnd
  END SUBROUTINE clo_bat
END MODULE closea