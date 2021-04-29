MODULE dtatsd
  USE oce
  USE phycst
  USE dom_oce
  USE fldread
  USE in_out_manager
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dta_tsd_init
  PUBLIC :: dta_tsd
  LOGICAL, PUBLIC :: ln_tsd_init
  LOGICAL, PUBLIC :: ln_tsd_dmp
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_tsd
  CONTAINS
  SUBROUTINE dta_tsd_init(ld_tradmp)
    LOGICAL, INTENT(IN), OPTIONAL :: ld_tradmp
    INTEGER :: ios, ierr0, ierr1, ierr2, ierr3
    CHARACTER(LEN = 100) :: cn_dir
    TYPE(FLD_N), DIMENSION(jpts) :: slf_i
    TYPE(FLD_N) :: sn_tem, sn_sal
    NAMELIST /namtsd/ ln_tsd_init, ln_tsd_dmp, cn_dir, sn_tem, sn_sal
    ierr0 = 0
    ierr1 = 0
    ierr2 = 0
    ierr3 = 0
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namtsd, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namtsd in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namtsd, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namtsd in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namtsd)
    IF (PRESENT(ld_tradmp)) ln_tsd_dmp = .TRUE.
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'dta_tsd_init : Temperature & Salinity data '
      WRITE(numout, *) '~~~~~~~~~~~~ '
      WRITE(numout, *) '   Namelist namtsd'
      WRITE(numout, *) '      Initialisation of ocean T & S with T &S input data   ln_tsd_init = ', ln_tsd_init
      WRITE(numout, *) '      damping of ocean T & S toward T &S input data        ln_tsd_dmp  = ', ln_tsd_dmp
      WRITE(numout, *)
      IF (.NOT. ln_tsd_init .AND. .NOT. ln_tsd_dmp) THEN
        WRITE(numout, *)
        WRITE(numout, *) '   ===>>   T & S data not used'
      END IF
    END IF
    IF (ln_rstart .AND. ln_tsd_init) THEN
      CALL ctl_warn('dta_tsd_init: ocean restart and T & S data intialisation, ', 'we keep the restart T & S values and set ln_tsd_init to FALSE')
      ln_tsd_init = .FALSE.
    END IF
    IF (ln_tsd_init .OR. ln_tsd_dmp) THEN
      ALLOCATE(sf_tsd(jpts), STAT = ierr0)
      IF (ierr0 > 0) THEN
        CALL ctl_stop('dta_tsd_init: unable to allocate sf_tsd structure')
        RETURN
      END IF
      ALLOCATE(sf_tsd(jp_tem) % fnow(jpi, jpj, jpk), STAT = ierr0)
      IF (sn_tem % ln_tint) ALLOCATE(sf_tsd(jp_tem) % fdta(jpi, jpj, jpk, 2), STAT = ierr1)
      ALLOCATE(sf_tsd(jp_sal) % fnow(jpi, jpj, jpk), STAT = ierr2)
      IF (sn_sal % ln_tint) ALLOCATE(sf_tsd(jp_sal) % fdta(jpi, jpj, jpk, 2), STAT = ierr3)
      IF (ierr0 + ierr1 + ierr2 + ierr3 > 0) THEN
        CALL ctl_stop('dta_tsd : unable to allocate T & S data arrays')
        RETURN
      END IF
      slf_i(jp_tem) = sn_tem
      slf_i(jp_sal) = sn_sal
      CALL fld_fill(sf_tsd, slf_i, cn_dir, 'dta_tsd', 'Temperature & Salinity data', 'namtsd', no_print)
    END IF
  END SUBROUTINE dta_tsd_init
  SUBROUTINE dta_tsd(kt, ptsd)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, jpts), INTENT(OUT) :: ptsd
    INTEGER :: ji, jj, jk, jl, jkk
    INTEGER :: ik, il0, il1, ii0, ii1, ij0, ij1
    REAL(KIND = wp) :: zl, zi
    REAL(KIND = wp), DIMENSION(jpk) :: ztp, zsp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('dta_tsd', 'r0', 0, 0)
    CALL fld_read(kt, 1, sf_tsd)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    IF (cn_cfg == "orca" .OR. cn_cfg == "ORCA") THEN
      IF (nn_cfg == 2 .AND. ln_tsd_dmp) THEN
        ij0 = 101
        ij1 = 109
        ii0 = 141
        ii1 = 155
        !$ACC loop independent collapse(2)
        DO jj = mj0(ij0), mj1(ij1)
          DO ji = mi0(ii0), mi1(ii1)
            sf_tsd(jp_tem) % fnow(ji, jj, 13 : 13) = sf_tsd(jp_tem) % fnow(ji, jj, 13 : 13) - 0.20_wp
            sf_tsd(jp_tem) % fnow(ji, jj, 14 : 15) = sf_tsd(jp_tem) % fnow(ji, jj, 14 : 15) - 0.35_wp
            sf_tsd(jp_tem) % fnow(ji, jj, 16 : 25) = sf_tsd(jp_tem) % fnow(ji, jj, 16 : 25) - 0.40_wp
            sf_tsd(jp_sal) % fnow(ji, jj, 13 : 13) = sf_tsd(jp_sal) % fnow(ji, jj, 13 : 13) - 0.15_wp
            sf_tsd(jp_sal) % fnow(ji, jj, 14 : 15) = sf_tsd(jp_sal) % fnow(ji, jj, 14 : 15) - 0.25_wp
            sf_tsd(jp_sal) % fnow(ji, jj, 16 : 17) = sf_tsd(jp_sal) % fnow(ji, jj, 16 : 17) - 0.30_wp
            sf_tsd(jp_sal) % fnow(ji, jj, 18 : 25) = sf_tsd(jp_sal) % fnow(ji, jj, 18 : 25) - 0.35_wp
          END DO
        END DO
        ij0 = 87
        ij1 = 96
        ii0 = 148
        ii1 = 160
        sf_tsd(jp_tem) % fnow(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1), 4 : 10) = 7.0_wp
        sf_tsd(jp_tem) % fnow(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1), 11 : 13) = 6.5_wp
        sf_tsd(jp_tem) % fnow(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1), 14 : 20) = 6.0_wp
      END IF
    END IF
    ptsd(:, :, :, jp_tem) = sf_tsd(jp_tem) % fnow(:, :, :)
    ptsd(:, :, :, jp_sal) = sf_tsd(jp_sal) % fnow(:, :, :)
    !$ACC END KERNELS
    IF (ln_sco) THEN
      CALL profile_psy_data1 % PreStart('dta_tsd', 'r1', 0, 0)
      IF (kt == nit000 .AND. lwp) THEN
        WRITE(numout, *)
        WRITE(numout, *) 'dta_tsd: interpolates T & S data onto the s- or mixed s-z-coordinate mesh'
      END IF
      CALL profile_psy_data1 % PostEnd
      DO jj = 1, jpj
        DO ji = 1, jpi
          CALL profile_psy_data2 % PreStart('dta_tsd', 'r2', 0, 0)
          DO jk = 1, jpk
            zl = gdept_0(ji, jj, jk)
            IF (zl < gdept_1d(1)) THEN
              ztp(jk) = ptsd(ji, jj, 1, jp_tem)
              zsp(jk) = ptsd(ji, jj, 1, jp_sal)
            ELSE IF (zl > gdept_1d(jpk)) THEN
              ztp(jk) = ptsd(ji, jj, jpkm1, jp_tem)
              zsp(jk) = ptsd(ji, jj, jpkm1, jp_sal)
            ELSE
              DO jkk = 1, jpkm1
                IF ((zl - gdept_1d(jkk)) * (zl - gdept_1d(jkk + 1)) <= 0._wp) THEN
                  zi = (zl - gdept_1d(jkk)) / (gdept_1d(jkk + 1) - gdept_1d(jkk))
                  ztp(jk) = ptsd(ji, jj, jkk, jp_tem) + (ptsd(ji, jj, jkk + 1, jp_tem) - ptsd(ji, jj, jkk, jp_tem)) * zi
                  zsp(jk) = ptsd(ji, jj, jkk, jp_sal) + (ptsd(ji, jj, jkk + 1, jp_sal) - ptsd(ji, jj, jkk, jp_sal)) * zi
                END IF
              END DO
            END IF
          END DO
          CALL profile_psy_data2 % PostEnd
          !$ACC KERNELS
          DO jk = 1, jpkm1
            ptsd(ji, jj, jk, jp_tem) = ztp(jk) * tmask(ji, jj, jk)
            ptsd(ji, jj, jk, jp_sal) = zsp(jk) * tmask(ji, jj, jk)
          END DO
          ptsd(ji, jj, jpk, jp_tem) = 0._wp
          ptsd(ji, jj, jpk, jp_sal) = 0._wp
          !$ACC END KERNELS
        END DO
      END DO
    ELSE
      !$ACC KERNELS
      ptsd(:, :, :, jp_tem) = ptsd(:, :, :, jp_tem) * tmask(:, :, :)
      ptsd(:, :, :, jp_sal) = ptsd(:, :, :, jp_sal) * tmask(:, :, :)
      !$ACC END KERNELS
      IF (ln_zps) THEN
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            ik = mbkt(ji, jj)
            IF (ik > 1) THEN
              zl = (gdept_1d(ik) - gdept_0(ji, jj, ik)) / (gdept_1d(ik) - gdept_1d(ik - 1))
              ptsd(ji, jj, ik, jp_tem) = (1. - zl) * ptsd(ji, jj, ik, jp_tem) + zl * ptsd(ji, jj, ik - 1, jp_tem)
              ptsd(ji, jj, ik, jp_sal) = (1. - zl) * ptsd(ji, jj, ik, jp_sal) + zl * ptsd(ji, jj, ik - 1, jp_sal)
            END IF
            ik = mikt(ji, jj)
            IF (ik > 1) THEN
              zl = (gdept_0(ji, jj, ik) - gdept_1d(ik)) / (gdept_1d(ik + 1) - gdept_1d(ik))
              ptsd(ji, jj, ik, jp_tem) = (1. - zl) * ptsd(ji, jj, ik, jp_tem) + zl * ptsd(ji, jj, ik + 1, jp_tem)
              ptsd(ji, jj, ik, jp_sal) = (1. - zl) * ptsd(ji, jj, ik, jp_sal) + zl * ptsd(ji, jj, ik + 1, jp_sal)
            END IF
          END DO
        END DO
        !$ACC END KERNELS
      END IF
    END IF
    CALL profile_psy_data3 % PreStart('dta_tsd', 'r3', 0, 0)
    IF (.NOT. ln_tsd_dmp) THEN
      IF (lwp) WRITE(numout, *) 'dta_tsd: deallocte T & S arrays as they are only use to initialize the run'
      DEALLOCATE(sf_tsd(jp_tem) % fnow)
      IF (sf_tsd(jp_tem) % ln_tint) DEALLOCATE(sf_tsd(jp_tem) % fdta)
      DEALLOCATE(sf_tsd(jp_sal) % fnow)
      IF (sf_tsd(jp_sal) % ln_tint) DEALLOCATE(sf_tsd(jp_sal) % fdta)
      DEALLOCATE(sf_tsd)
    END IF
    CALL profile_psy_data3 % PostEnd
  END SUBROUTINE dta_tsd
END MODULE dtatsd