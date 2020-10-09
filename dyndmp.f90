MODULE dyndmp
  USE oce
  USE dom_oce
  USE c1d
  USE tradmp
  USE zdf_oce
  USE phycst
  USE dtauvd
  USE zdfmxl
  USE in_out_manager
  USE lib_mpp
  USE prtctl
  USE timing
  USE iom
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_dmp_init
  PUBLIC :: dyn_dmp
  LOGICAL, PUBLIC :: ln_dyndmp
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: utrdmp
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: vtrdmp
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: resto_uv
  CONTAINS
  INTEGER FUNCTION dyn_dmp_alloc()
    ALLOCATE(utrdmp(jpi, jpj, jpk), vtrdmp(jpi, jpj, jpk), resto_uv(jpi, jpj, jpk), STAT = dyn_dmp_alloc)
    CALL mpp_sum('dyndmp', dyn_dmp_alloc)
    IF (dyn_dmp_alloc > 0) CALL ctl_warn('dyn_dmp_alloc: allocation of arrays failed')
  END FUNCTION dyn_dmp_alloc
  SUBROUTINE dyn_dmp_init
    INTEGER :: ios, imask
    NAMELIST /namc1d_dyndmp/ ln_dyndmp
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namc1d_dyndmp, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namc1d_dyndmp in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namc1d_dyndmp, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namc1d_dyndmp in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namc1d_dyndmp)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dyn_dmp_init : U and V current Newtonian damping'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namc1d_dyndmp : Set damping flag'
      WRITE(numout, FMT = *) '      add a damping term or not       ln_dyndmp = ', ln_dyndmp
      WRITE(numout, FMT = *) '   Namelist namtra_dmp    : Set damping parameters'
      WRITE(numout, FMT = *) '      Apply relaxation   or not       ln_tradmp = ', ln_tradmp
      WRITE(numout, FMT = *) '      mixed layer damping option      nn_zdmp   = ', nn_zdmp
      WRITE(numout, FMT = *) '      Damping file name               cn_resto  = ', cn_resto
      WRITE(numout, FMT = *)
    END IF
    IF (ln_dyndmp) THEN
      IF (dyn_dmp_alloc() /= 0) CALL ctl_stop('STOP', 'dyn_dmp_init: unable to allocate arrays')
      SELECT CASE (nn_zdmp)
      CASE (0)
        IF (lwp) WRITE(numout, FMT = *) '   momentum damping throughout the water column'
      CASE (1)
        IF (lwp) WRITE(numout, FMT = *) '   no momentum damping in the turbocline (avt > 5 cm2/s)'
      CASE (2)
        IF (lwp) WRITE(numout, FMT = *) '   no momentum damping in the mixed layer'
      CASE DEFAULT
        WRITE(ctmp1, FMT = *) '          bad flag value for nn_zdmp = ', nn_zdmp
        CALL ctl_stop(ctmp1)
      END SELECT
      IF (.NOT. ln_uvd_dyndmp) THEN
        CALL ctl_warn('dyn_dmp_init: U & V current read data not initialized, we force ln_uvd_dyndmp=T')
        CALL dta_uvd_init(ld_dyndmp = ln_dyndmp)
      END IF
      !$ACC KERNELS
      utrdmp(:, :, :) = 0._wp
      vtrdmp(:, :, :) = 0._wp
      !$ACC END KERNELS
      CALL iom_open(cn_resto, imask)
      CALL iom_get(imask, jpdom_autoglo, 'resto', resto)
      CALL iom_close(imask)
    END IF
  END SUBROUTINE dyn_dmp_init
  SUBROUTINE dyn_dmp(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zua, zva
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, 2) :: zuv_dta
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('dyn_dmp', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('dyn_dmp')
    CALL dta_uvd(kt, zuv_dta)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    SELECT CASE (nn_zdmp)
    CASE (0)
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zua = resto_uv(ji, jj, jk) * (zuv_dta(ji, jj, jk, 1) - ub(ji, jj, jk))
            zva = resto_uv(ji, jj, jk) * (zuv_dta(ji, jj, jk, 2) - vb(ji, jj, jk))
            ua(ji, jj, jk) = ua(ji, jj, jk) + zua
            va(ji, jj, jk) = va(ji, jj, jk) + zva
            utrdmp(ji, jj, jk) = zua
            vtrdmp(ji, jj, jk) = zva
          END DO
        END DO
      END DO
    CASE (1)
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            IF (avt(ji, jj, jk) <= avt_c) THEN
              zua = resto_uv(ji, jj, jk) * (zuv_dta(ji, jj, jk, 1) - ub(ji, jj, jk))
              zva = resto_uv(ji, jj, jk) * (zuv_dta(ji, jj, jk, 2) - vb(ji, jj, jk))
            ELSE
              zua = 0._wp
              zva = 0._wp
            END IF
            ua(ji, jj, jk) = ua(ji, jj, jk) + zua
            va(ji, jj, jk) = va(ji, jj, jk) + zva
            utrdmp(ji, jj, jk) = zua
            vtrdmp(ji, jj, jk) = zva
          END DO
        END DO
      END DO
    CASE (2)
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            IF (gdept_n(ji, jj, jk) >= hmlp(ji, jj)) THEN
              zua = resto_uv(ji, jj, jk) * (zuv_dta(ji, jj, jk, 1) - ub(ji, jj, jk))
              zva = resto_uv(ji, jj, jk) * (zuv_dta(ji, jj, jk, 2) - vb(ji, jj, jk))
            ELSE
              zua = 0._wp
              zva = 0._wp
            END IF
            ua(ji, jj, jk) = ua(ji, jj, jk) + zua
            va(ji, jj, jk) = va(ji, jj, jk) + zva
            utrdmp(ji, jj, jk) = zua
            vtrdmp(ji, jj, jk) = zva
          END DO
        END DO
      END DO
    END SELECT
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('dyn_dmp', 'r1', 0, 0)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = ua(:, :, :), clinfo1 = ' dmp  - Ua: ', mask1 = umask, tab3d_2 = va(:, :, :), clinfo2 = ' &
&Va: ', mask2 = vmask, clinfo3 = 'dyn')
    IF (ln_timing) CALL timing_stop('dyn_dmp')
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE dyn_dmp
END MODULE dyndmp