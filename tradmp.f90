MODULE tradmp
  USE oce
  USE dom_oce
  USE c1d
  USE trd_oce
  USE trdtra
  USE zdf_oce
  USE phycst
  USE dtatsd
  USE zdfmxl
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE prtctl
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_dmp
  PUBLIC :: tra_dmp_init
  LOGICAL, PUBLIC :: ln_tradmp
  INTEGER, PUBLIC :: nn_zdmp
  CHARACTER(LEN = 200), PUBLIC :: cn_resto
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: resto
  CONTAINS
  INTEGER FUNCTION tra_dmp_alloc()
    ALLOCATE(resto(jpi, jpj, jpk), STAT = tra_dmp_alloc)
    CALL mpp_sum('tradmp', tra_dmp_alloc)
    IF (tra_dmp_alloc > 0) CALL ctl_warn('tra_dmp_alloc: allocation of arrays failed')
  END FUNCTION tra_dmp_alloc
  SUBROUTINE tra_dmp(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk, jn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, jpts) :: zts_dta
    REAL(KIND = wp), DIMENSION(:, :, :, :), ALLOCATABLE :: ztrdts
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    IF (ln_timing) CALL timing_start('tra_dmp')
    IF (l_trdtra) THEN
      ALLOCATE(ztrdts(jpi, jpj, jpk, jpts))
      !$ACC KERNELS
      ztrdts(:, :, :, :) = tsa(:, :, :, :)
      !$ACC END KERNELS
    END IF
    CALL dta_tsd(kt, zts_dta)
    SELECT CASE (nn_zdmp)
    CASE (0)
      DO jn = 1, jpts
        !$ACC KERNELS
        DO jk = 1, jpkm1
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              tsa(ji, jj, jk, jn) = tsa(ji, jj, jk, jn) + resto(ji, jj, jk) * (zts_dta(ji, jj, jk, jn) - tsb(ji, jj, jk, jn))
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      END DO
    CASE (1)
      !$ACC KERNELS
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            IF (avt(ji, jj, jk) <= avt_c) THEN
              tsa(ji, jj, jk, jp_tem) = tsa(ji, jj, jk, jp_tem) + resto(ji, jj, jk) * (zts_dta(ji, jj, jk, jp_tem) - tsb(ji, jj, &
&jk, jp_tem))
              tsa(ji, jj, jk, jp_sal) = tsa(ji, jj, jk, jp_sal) + resto(ji, jj, jk) * (zts_dta(ji, jj, jk, jp_sal) - tsb(ji, jj, &
&jk, jp_sal))
            END IF
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    CASE (2)
      !$ACC KERNELS
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            IF (gdept_n(ji, jj, jk) >= hmlp(ji, jj)) THEN
              tsa(ji, jj, jk, jp_tem) = tsa(ji, jj, jk, jp_tem) + resto(ji, jj, jk) * (zts_dta(ji, jj, jk, jp_tem) - tsb(ji, jj, &
&jk, jp_tem))
              tsa(ji, jj, jk, jp_sal) = tsa(ji, jj, jk, jp_sal) + resto(ji, jj, jk) * (zts_dta(ji, jj, jk, jp_sal) - tsb(ji, jj, &
&jk, jp_sal))
            END IF
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END SELECT
    IF (l_trdtra) THEN
      !$ACC KERNELS
      ztrdts(:, :, :, :) = tsa(:, :, :, :) - ztrdts(:, :, :, :)
      !$ACC END KERNELS
      CALL profile_psy_data0 % PreStart('tra_dmp', 'r0', 0, 0)
      CALL trd_tra(kt, 'TRA', jp_tem, jptra_dmp, ztrdts(:, :, :, jp_tem))
      CALL trd_tra(kt, 'TRA', jp_sal, jptra_dmp, ztrdts(:, :, :, jp_sal))
      DEALLOCATE(ztrdts)
      CALL profile_psy_data0 % PostEnd
    END IF
    CALL profile_psy_data1 % PreStart('tra_dmp', 'r1', 0, 0)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = tsa(:, :, :, jp_tem), clinfo1 = ' dmp  - Ta: ', mask1 = tmask, tab3d_2 = tsa(:, :, :, &
&jp_sal), clinfo2 = ' Sa: ', mask2 = tmask, clinfo3 = 'tra')
    IF (ln_timing) CALL timing_stop('tra_dmp')
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE tra_dmp
  SUBROUTINE tra_dmp_init
    INTEGER :: ios, imask
    NAMELIST /namtra_dmp/ ln_tradmp, nn_zdmp, cn_resto
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namtra_dmp, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namtra_dmp in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namtra_dmp, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namtra_dmp in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namtra_dmp)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'tra_dmp_init : T and S newtonian relaxation'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namtra_dmp : set relaxation parameters'
      WRITE(numout, FMT = *) '      Apply relaxation   or not       ln_tradmp   = ', ln_tradmp
      WRITE(numout, FMT = *) '         mixed layer damping option      nn_zdmp  = ', nn_zdmp
      WRITE(numout, FMT = *) '         Damping file name               cn_resto = ', cn_resto
      WRITE(numout, FMT = *)
    END IF
    IF (ln_tradmp) THEN
      IF (tra_dmp_alloc() /= 0) CALL ctl_stop('STOP', 'tra_dmp_init: unable to allocate arrays')
      SELECT CASE (nn_zdmp)
      CASE (0)
        IF (lwp) WRITE(numout, FMT = *) '   tracer damping as specified by mask'
      CASE (1)
        IF (lwp) WRITE(numout, FMT = *) '   no tracer damping in the mixing layer (kz > 5 cm2/s)'
      CASE (2)
        IF (lwp) WRITE(numout, FMT = *) '   no tracer damping in the mixed  layer'
      CASE DEFAULT
        CALL ctl_stop('tra_dmp_init : wrong value of nn_zdmp')
      END SELECT
      IF (.NOT. ln_tsd_dmp) THEN
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) '   read T-S data not initialized, we force ln_tsd_dmp=T'
        CALL dta_tsd_init(ld_tradmp = ln_tradmp)
      END IF
      CALL iom_open(cn_resto, imask)
      CALL iom_get(imask, jpdom_autoglo, 'resto', resto)
      CALL iom_close(imask)
    END IF
  END SUBROUTINE tra_dmp_init
END MODULE tradmp