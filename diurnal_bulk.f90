MODULE diurnal_bulk
  USE par_kind
  USE phycst
  USE dom_oce
  USE lib_mpp
  USE solfrac_mod
  USE in_out_manager
  IMPLICIT NONE
  PRIVATE
  LOGICAL, PUBLIC :: ln_diurnal
  LOGICAL, PUBLIC :: ln_diurnal_only
  REAL(KIND = wp), PRIVATE, PARAMETER :: pp_alpha = 2.0E-4_wp
  REAL(KIND = wp), PRIVATE, PARAMETER :: pp_veltol = 0._wp
  REAL(KIND = wp), PRIVATE, PARAMETER :: pp_min_fvel = 1.E-10_wp
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: x_dsst
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: x_solfrac
  PUBLIC :: diurnal_sst_bulk_init, diurnal_sst_takaya_step
  CONTAINS
  SUBROUTINE diurnal_sst_bulk_init
    INTEGER :: ios
    NAMELIST /namdiu/ ln_diurnal, ln_diurnal_only
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namdiu, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namdiu in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namdiu, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namdiu in configuration namelist', lwp)
    IF (ln_diurnal_only .AND. (.NOT. ln_diurnal)) THEN
      CALL ctl_stop("ln_diurnal_only set, but ln_diurnal = FALSE !")
    END IF
    IF (ln_diurnal) THEN
      ALLOCATE(x_dsst(jpi, jpj), x_solfrac(jpi, jpj))
      x_solfrac = 0._wp
      x_dsst = 0._wp
      IF (ln_diurnal_only) THEN
        CALL ctl_warn("ln_diurnal_only set; only the diurnal component of SST will be calculated")
      END IF
    END IF
  END SUBROUTINE diurnal_sst_bulk_init
  SUBROUTINE diurnal_sst_takaya_step(kt, psolflux, pqflux, ptauflux, prho, p_rdt, pla, pthick, pcoolthick, pmu, p_fvel_bkginc, p_hflux_bkginc)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: psolflux
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pqflux
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: ptauflux
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: prho
    REAL(KIND = wp), INTENT(IN) :: p_rdt
    REAL(KIND = wp), DIMENSION(jpi, jpj), OPTIONAL, INTENT(IN) :: pLa
    REAL(KIND = wp), DIMENSION(jpi, jpj), OPTIONAL, INTENT(IN) :: pthick
    REAL(KIND = wp), DIMENSION(jpi, jpj), OPTIONAL, INTENT(IN) :: pcoolthick
    REAL(KIND = wp), DIMENSION(jpi, jpj), OPTIONAL, INTENT(IN) :: pmu
    REAL(KIND = wp), DIMENSION(jpi, jpj), OPTIONAL, INTENT(IN) :: p_hflux_bkginc
    REAL(KIND = wp), DIMENSION(jpi, jpj), OPTIONAL, INTENT(IN) :: p_fvel_bkginc
    INTEGER :: ji, jj
    LOGICAL :: ll_calcfrac
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z_fvel
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zthick, zcoolthick, zmu, zla
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z_abflux
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z_fla
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    IF (.NOT. PRESENT(pthick)) THEN
      !$ACC KERNELS
      zthick(:, :) = 3._wp
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      zthick(:, :) = pthick(:, :)
      !$ACC END KERNELS
    END IF
    IF (.NOT. PRESENT(pcoolthick)) THEN
      !$ACC KERNELS
      zcoolthick(:, :) = 0._wp
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      zcoolthick(:, :) = pcoolthick(:, :)
      !$ACC END KERNELS
    END IF
    IF (.NOT. PRESENT(pmu)) THEN
      !$ACC KERNELS
      zmu(:, :) = 0.3_wp
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      zmu(:, :) = pmu(:, :)
      !$ACC END KERNELS
    END IF
    IF (.NOT. PRESENT(pla)) THEN
      !$ACC KERNELS
      zla(:, :) = 0.3_wp
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      zla(:, :) = pla(:, :)
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data0 % PreStart('diurnal_sst_takaya_step', 'r0', 0, 0)
    IF (kt == nit000) THEN
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF ((x_solfrac(ji, jj) == 0._wp) .AND. (tmask(ji, jj, 1) == 1._wp)) x_solfrac(ji, jj) = solfrac(zcoolthick(ji, jj), zthick(ji, jj))
        END DO
      END DO
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    WHERE (tmask(:, :, 1) == 1._wp)
      z_abflux(:, :) = (x_solfrac(:, :) * psolflux(:, :)) + pqflux(:, :)
    ELSEWHERE
      z_abflux(:, :) = 0._wp
    END WHERE
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('diurnal_sst_takaya_step', 'r1', 0, 0)
    IF (PRESENT(p_hflux_bkginc)) z_abflux(:, :) = z_abflux(:, :) + p_hflux_bkginc
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    WHERE (ABS(z_abflux(:, :)) < rsmall)
      z_abflux(:, :) = rsmall
    END WHERE
    WHERE ((ptauflux /= 0) .AND. (tmask(:, :, 1) == 1.))
      z_fvel(:, :) = SQRT(ptauflux(:, :) / prho(:, :))
    ELSEWHERE
      z_fvel(:, :) = 0._wp
    END WHERE
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('diurnal_sst_takaya_step', 'r2', 0, 0)
    IF (PRESENT(p_fvel_bkginc)) z_fvel(:, :) = z_fvel(:, :) + p_fvel_bkginc
    CALL profile_psy_data2 % PostEnd
    !$ACC KERNELS
    WHERE (tmask(:, :, 1) == 1.)
      z_fla(:, :) = MAX(1._wp, zla(:, :) ** (- 2._wp / 3._wp))
    ELSEWHERE
      z_fla(:, :) = 0._wp
    END WHERE
    !$ACC END KERNELS
    CALL profile_psy_data3 % PreStart('diurnal_sst_takaya_step', 'r3', 0, 0)
    x_dsst(:, :) = t_imp(x_dsst(:, :), p_rdt, z_abflux(:, :), z_fvel(:, :), z_fla(:, :), zmu(:, :), zthick(:, :), prho(:, :))
    CALL profile_psy_data3 % PostEnd
  END SUBROUTINE diurnal_sst_takaya_step
  FUNCTION t_imp(p_dsst, p_rdt, p_abflux, p_fvel, p_fla, pmu, pthick, prho)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    IMPLICIT NONE
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: t_imp
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: p_dsst
    REAL(KIND = wp), INTENT(IN) :: p_rdt
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: p_abflux
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: p_fvel
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: p_fla
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pmu
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pthick
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: prho
    REAL(KIND = wp) :: z_olength
    REAL(KIND = wp) :: z_sigma, z_sigma2
    REAL(KIND = wp) :: z_term1, z_term2
    REAL(KIND = wp) :: z_stabfunc
    REAL(KIND = wp) :: z_fvel
    CHARACTER(LEN = 200) :: warn_string
    LOGICAL :: lwarn
    INTEGER :: ji, jj
    !TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !CALL profile_psy_data0 % PreStart('t_imp', 'r0', 0, 0)
    lwarn = .FALSE.
    !ARPDBG - not done by PSyclone because of CYCLE?
    !$ACC KERNELS
    !$ACC LOOP COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (tmask(ji, jj, 1) /= 1._wp) THEN
          t_imp(ji, jj) = 0._wp
          CYCLE
        END IF
        IF (p_fvel(ji, jj) < pp_min_fvel) THEN
          z_fvel = pp_min_fvel
          lwarn = .TRUE.
        ELSE
          z_fvel = p_fvel(ji, jj)
        END IF
        IF ((z_fvel < pp_veltol) .AND. (p_dsst(ji, jj) > 0._wp)) THEN
          z_olength = z_fvel / SQRT(p_dsst(ji, jj) * vkarmn * grav * pp_alpha / (5._wp * pthick(ji, jj)))
        ELSE
          z_olength = (prho(ji, jj) * rcp * z_fvel ** 3._wp) / (vkarmn * grav * pp_alpha * p_abflux(ji, jj))
        END IF
        z_sigma = pthick(ji, jj) / z_olength
        z_sigma2 = z_sigma * z_sigma
        IF (z_sigma >= 0.) THEN
          z_stabfunc = 1._wp + ((5._wp * z_sigma + 4._wp * z_sigma2) / (1._wp + 3._wp * z_sigma + 0.25_wp * z_sigma2))
        ELSE
          z_stabfunc = 1._wp / SQRT(1._wp - 16._wp * z_sigma)
        END IF
        z_term1 = (p_abflux(ji, jj) * (pmu(ji, jj) + 1._wp) / (pmu(ji, jj) * pthick(ji, jj) * prho(ji, jj) * rcp))
        z_term2 = - ((pmu(ji, jj) + 1._wp) * (vkarmn * z_fvel * p_fla(ji, jj)) / (pthick(ji, jj) * z_stabfunc))
        t_imp(ji, jj) = (p_dsst(ji, jj) + p_rdt * z_term1) / (1._wp - p_rdt * z_term2)
      END DO
    END DO
    !$ACC END KERNELS
    IF (lwarn) THEN
      WRITE(warn_string, FMT = *) "diurnal_sst_takaya step: " // "friction velocity < minimum\n" // "Setting friction velocity =", pp_min_fvel
      CALL ctl_warn(warn_string)
    END IF
    !CALL profile_psy_data0 % PostEnd
  END FUNCTION t_imp
END MODULE diurnal_bulk
