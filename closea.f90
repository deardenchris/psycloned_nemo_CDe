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
    INTEGER :: inum
    INTEGER :: ierr
    INTEGER :: id
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdata_in
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) 'dom_clo : read in masks to define closed seas '
    IF (lwp) WRITE(numout, FMT = *) '~~~~~~~'
    IF (ln_read_cfg) THEN
      CALL iom_open(cn_domcfg, inum)
      id = iom_varid(inum, 'closea_mask', ldstop = .FALSE.)
      IF (id > 0) THEN
        l_sbc_clo = .TRUE.
        ALLOCATE(closea_mask(jpi, jpj), STAT = ierr)
        IF (ierr /= 0) CALL ctl_stop('STOP', 'dom_clo: failed to allocate closea_mask array')
        zdata_in(:, :) = 0.0
        CALL iom_get(inum, jpdom_data, 'closea_mask', zdata_in)
        closea_mask(:, :) = NINT(zdata_in(:, :)) * tmask(:, :, 1)
        jncs = MAXVAL(closea_mask(:, :))
        CALL mpp_max('closea', jncs)
        IF (jncs > 0) THEN
          IF (lwp) WRITE(numout, FMT = *) 'Number of closed seas : ', jncs
        ELSE
          CALL ctl_stop('Problem with closea_mask field in domain_cfg file. Has no values > 0 so no closed seas defined.')
        END IF
      ELSE
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   closea_mask field not found in domain_cfg file.'
        IF (lwp) WRITE(numout, FMT = *) '           No closed seas defined.'
        IF (lwp) WRITE(numout, FMT = *)
        l_sbc_clo = .FALSE.
        jncs = 0
      END IF
      l_clo_rnf = .FALSE.
      IF (l_sbc_clo) THEN
        id = iom_varid(inum, 'closea_mask_rnf', ldstop = .FALSE.)
        IF (id > 0) THEN
          l_clo_rnf = .TRUE.
          ALLOCATE(closea_mask_rnf(jpi, jpj), STAT = ierr)
          IF (ierr /= 0) CALL ctl_stop('STOP', 'dom_clo: failed to allocate closea_mask_rnf array')
          CALL iom_get(inum, jpdom_data, 'closea_mask_rnf', zdata_in)
          closea_mask_rnf(:, :) = NINT(zdata_in(:, :)) * tmask(:, :, 1)
          jncsr = MAXVAL(closea_mask_rnf(:, :))
          CALL mpp_max('closea', jncsr)
          IF (jncsr > 0) THEN
            IF (lwp) WRITE(numout, FMT = *) 'Number of closed seas rnf mappings : ', jncsr
          ELSE
            CALL ctl_stop('Problem with closea_mask_rnf field in domain_cfg file. Has no values > 0 so no closed seas rnf mappings &
&defined.')
          END IF
        ELSE
          IF (lwp) WRITE(numout, FMT = *) 'closea_mask_rnf field not found in domain_cfg file. No closed seas rnf mappings defined.'
          jncsr = 0
        END IF
        id = iom_varid(inum, 'closea_mask_empmr', ldstop = .FALSE.)
        IF (id > 0) THEN
          l_clo_rnf = .TRUE.
          ALLOCATE(closea_mask_empmr(jpi, jpj), STAT = ierr)
          IF (ierr /= 0) CALL ctl_stop('STOP', 'dom_clo: failed to allocate closea_mask_empmr array')
          CALL iom_get(inum, jpdom_data, 'closea_mask_empmr', zdata_in)
          closea_mask_empmr(:, :) = NINT(zdata_in(:, :)) * tmask(:, :, 1)
          jncse = MAXVAL(closea_mask_empmr(:, :))
          CALL mpp_max('closea', jncse)
          IF (jncse > 0) THEN
            IF (lwp) WRITE(numout, FMT = *) 'Number of closed seas empmr mappings : ', jncse
          ELSE
            CALL ctl_stop('Problem with closea_mask_empmr field in domain_cfg file. Has no values > 0 so no closed seas empmr &
&mappings defined.')
          END IF
        ELSE
          IF (lwp) WRITE(numout, FMT = *) 'closea_mask_empmr field not found in domain_cfg file. No closed seas empmr mappings &
&defined.'
          jncse = 0
        END IF
      END IF
      CALL iom_close(inum)
    ELSE
      IF (lwp) WRITE(numout, FMT = *) 'No domain_cfg file so no closed seas defined.'
      l_sbc_clo = .FALSE.
      l_clo_rnf = .FALSE.
    END IF
  END SUBROUTINE dom_clo
  SUBROUTINE sbc_clo(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ierr
    INTEGER :: jc, jcr, jce
    REAL(KIND = wp), PARAMETER :: rsmall = 1.E-20_wp
    REAL(KIND = wp) :: zfwf_total, zcoef, zcoef1
    REAL(KIND = wp), DIMENSION(jncs) :: zfwf
    REAL(KIND = wp), DIMENSION(jncsr + 1) :: zfwfr
    REAL(KIND = wp), DIMENSION(jncse + 1) :: zfwfe
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztmp2d
    IF (ln_timing) CALL timing_start('sbc_clo')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'sbc_clo : closed seas '
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~'
      ALLOCATE(surf(jncs + 1), STAT = ierr)
      IF (ierr /= 0) CALL ctl_stop('STOP', 'sbc_clo: failed to allocate surf array')
      surf(:) = 0.E0_wp
      ALLOCATE(surfr(jncsr + 1), STAT = ierr)
      IF (ierr /= 0) CALL ctl_stop('STOP', 'sbc_clo: failed to allocate surfr array')
      surfr(:) = 0.E0_wp
      ALLOCATE(surfe(jncse + 1), STAT = ierr)
      IF (ierr /= 0) CALL ctl_stop('STOP', 'sbc_clo: failed to allocate surfe array')
      surfe(:) = 0.E0_wp
      surf(jncs + 1) = glob_sum('closea', e1e2t(:, :))
      DO jc = 1, jncs
        ztmp2d(:, :) = 0.E0_wp
        WHERE (closea_mask(:, :) == jc) ztmp2d(:, :) = e1e2t(:, :) * tmask_i(:, :)
        surf(jc) = glob_sum('closea', ztmp2d(:, :))
      END DO
      surf(jncs + 1) = surf(jncs + 1) - SUM(surf(1 : jncs))
      IF (jncsr > 0) THEN
        DO jcr = 1, jncsr
          ztmp2d(:, :) = 0.E0_wp
          WHERE (closea_mask_rnf(:, :) == jcr .AND. closea_mask(:, :) == 0) ztmp2d(:, :) = e1e2t(:, :) * tmask_i(:, :)
          surfr(jcr) = glob_sum('closea', ztmp2d(:, :))
        END DO
      END IF
      IF (jncse > 0) THEN
        DO jce = 1, jncse
          ztmp2d(:, :) = 0.E0_wp
          WHERE (closea_mask_empmr(:, :) == jce .AND. closea_mask(:, :) == 0) ztmp2d(:, :) = e1e2t(:, :) * tmask_i(:, :)
          surfe(jce) = glob_sum('closea', ztmp2d(:, :))
        END DO
      END IF
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
    END IF
    zfwf_total = 0._wp
    zfwf(:) = 0.E0_wp
    DO jc = 1, jncs
      ztmp2d(:, :) = 0.E0_wp
      WHERE (closea_mask(:, :) == jc) ztmp2d(:, :) = e1e2t(:, :) * (emp(:, :) - rnf(:, :)) * tmask_i(:, :)
      zfwf(jc) = glob_sum('closea', ztmp2d(:, :))
    END DO
    zfwf_total = SUM(zfwf)
    zfwfr(:) = 0.E0_wp
    IF (jncsr > 0) THEN
      DO jcr = 1, jncsr
        ztmp2d(:, :) = 0.E0_wp
        WHERE (closea_mask_rnf(:, :) == jcr .AND. closea_mask(:, :) > 0) ztmp2d(:, :) = e1e2t(:, :) * (emp(:, :) - rnf(:, :)) * &
&tmask_i(:, :)
        zfwfr(jcr) = glob_sum('closea', ztmp2d(:, :))
        IF (ABS(zfwfr(jcr) / surf(jncs + 1)) > rsmall) THEN
          IF (zfwfr(jcr) < 0.0) THEN
            zfwf_total = zfwf_total - zfwfr(jcr)
            zcoef = zfwfr(jcr) / surfr(jcr)
            zcoef1 = rcp * zcoef
            WHERE (closea_mask_rnf(:, :) == jcr .AND. closea_mask(:, :) == 0.0)
              emp(:, :) = emp(:, :) + zcoef
              qns(:, :) = qns(:, :) - zcoef1 * sst_m(:, :)
            END WHERE
          END IF
        END IF
      END DO
    END IF
    zfwfe(:) = 0.E0_wp
    IF (jncse > 0) THEN
      DO jce = 1, jncse
        ztmp2d(:, :) = 0.E0_wp
        WHERE (closea_mask_empmr(:, :) == jce .AND. closea_mask(:, :) > 0) ztmp2d(:, :) = e1e2t(:, :) * (emp(:, :) - rnf(:, :)) * &
&tmask_i(:, :)
        zfwfe(jce) = glob_sum('closea', ztmp2d(:, :))
        IF (ABS(zfwfe(jce) / surf(jncs + 1)) > rsmall) THEN
          zfwf_total = zfwf_total - zfwfe(jce)
          zcoef = zfwfe(jce) / surfe(jce)
          zcoef1 = rcp * zcoef
          WHERE (closea_mask_empmr(:, :) == jce .AND. closea_mask(:, :) == 0.0)
            emp(:, :) = emp(:, :) + zcoef
            qns(:, :) = qns(:, :) - zcoef1 * sst_m(:, :)
          END WHERE
        END IF
      END DO
    END IF
    IF (ABS(zfwf_total / surf(jncs + 1)) > rsmall) THEN
      zcoef = zfwf_total / surf(jncs + 1)
      zcoef1 = rcp * zcoef
      WHERE (closea_mask(:, :) == 0)
        emp(:, :) = emp(:, :) + zcoef
        qns(:, :) = qns(:, :) - zcoef1 * sst_m(:, :)
      END WHERE
    END IF
    DO jc = 1, jncs
      IF (ABS(zfwf(jc) / surf(jncs + 1)) > rsmall) THEN
        zcoef = zfwf(jc) / surf(jc)
        zcoef1 = rcp * zcoef
        WHERE (closea_mask(:, :) == jc)
          emp(:, :) = emp(:, :) - zcoef
          qns(:, :) = qns(:, :) + zcoef1 * sst_m(:, :)
        END WHERE
      END IF
    END DO
    emp(:, :) = emp(:, :) * tmask(:, :, 1)
    CALL lbc_lnk('closea', emp, 'T', 1._wp)
    IF (ln_timing) CALL timing_stop('sbc_clo')
  END SUBROUTINE sbc_clo
  SUBROUTINE clo_rnf(p_rnfmsk)
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: p_rnfmsk
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
  END SUBROUTINE clo_rnf
  SUBROUTINE clo_bat(k_top, k_bot)
    INTEGER, DIMENSION(:, :), INTENT(INOUT) :: k_top, k_bot
    INTEGER :: inum, id
    INTEGER, DIMENSION(jpi, jpj) :: closea_mask
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdata_in
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'clo_bat : suppression of closed seas'
      WRITE(numout, FMT = *) '~~~~~~~'
    END IF
    IF (ln_read_cfg) THEN
      CALL iom_open(cn_domcfg, inum)
      id = iom_varid(inum, 'closea_mask', ldstop = .FALSE.)
      IF (id > 0) THEN
        IF (lwp) WRITE(numout, FMT = *) 'Suppressing closed seas in bathymetry based on closea_mask field,'
        CALL iom_get(inum, jpdom_data, 'closea_mask', zdata_in)
        closea_mask(:, :) = NINT(zdata_in(:, :))
        WHERE (closea_mask(:, :) > 0)
          k_top(:, :) = 0
          k_bot(:, :) = 0
        END WHERE
      ELSE
        IF (lwp) WRITE(numout, FMT = *) 'No closea_mask field found in domain_cfg file. No suppression of closed seas.'
      END IF
      CALL iom_close(inum)
    ELSE
      IF (lwp) WRITE(numout, FMT = *) 'No domain_cfg file => no suppression of closed seas.'
    END IF
    l_sbc_clo = .FALSE.
    l_clo_rnf = .FALSE.
  END SUBROUTINE clo_bat
END MODULE closea