MODULE zdfric
  USE oce
  USE dom_oce
  USE zdf_oce
  USE phycst
  USE sbc_oce, ONLY: taum
  USE in_out_manager
  USE iom
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: zdf_ric
  PUBLIC :: ric_rst
  PUBLIC :: zdf_ric_init
  INTEGER :: nn_ric
  REAL(KIND = wp) :: rn_avmri
  REAL(KIND = wp) :: rn_alp
  REAL(KIND = wp) :: rn_ekmfc
  REAL(KIND = wp) :: rn_mldmin
  REAL(KIND = wp) :: rn_mldmax
  REAL(KIND = wp) :: rn_wtmix
  REAL(KIND = wp) :: rn_wvmix
  LOGICAL :: ln_mldw
  CONTAINS
  SUBROUTINE zdf_ric_init
    INTEGER :: ji, jj, jk
    INTEGER :: ios
    NAMELIST /namzdf_ric/ rn_avmri, rn_alp, nn_ric, rn_ekmfc, rn_mldmin, rn_mldmax, rn_wtmix, rn_wvmix, ln_mldw
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namzdf_ric, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namzdf_ric in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namzdf_ric, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namzdf_ric in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namzdf_ric)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'zdf_ric_init : Ri depend vertical mixing scheme'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namzdf_ric : set Kz=F(Ri) parameters'
      WRITE(numout, FMT = *) '      maximum vertical viscosity        rn_avmri  = ', rn_avmri
      WRITE(numout, FMT = *) '      coefficient                       rn_alp    = ', rn_alp
      WRITE(numout, FMT = *) '      exponent                          nn_ric    = ', nn_ric
      WRITE(numout, FMT = *) '      Ekman layer enhanced mixing       ln_mldw   = ', ln_mldw
      WRITE(numout, FMT = *) '         Ekman Factor Coeff             rn_ekmfc  = ', rn_ekmfc
      WRITE(numout, FMT = *) '         minimum mixed layer depth      rn_mldmin = ', rn_mldmin
      WRITE(numout, FMT = *) '         maximum mixed layer depth      rn_mldmax = ', rn_mldmax
      WRITE(numout, FMT = *) '         Vertical eddy Diff. in the ML  rn_wtmix  = ', rn_wtmix
      WRITE(numout, FMT = *) '         Vertical eddy Visc. in the ML  rn_wvmix  = ', rn_wvmix
    END IF
    CALL ric_rst(nit000, 'READ')
    IF (lwxios) THEN
      CALL iom_set_rstw_var_active('avt_k')
      CALL iom_set_rstw_var_active('avm_k')
    END IF
  END SUBROUTINE zdf_ric_init
  SUBROUTINE zdf_ric(kt, pdept, p_sh2, p_avm, p_avt)
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pdept
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: p_sh2
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: p_avm, p_avt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zcfRi, zav, zustar, zhek
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zh_ekm
    !$ACC KERNELS
    DO jk = 2, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zcfRi = 1._wp / (1._wp + rn_alp * MAX(0._wp, avm(ji, jj, jk) * rn2(ji, jj, jk) / (p_sh2(ji, jj, jk) + 1.E-20)))
          zav = rn_avmri * zcfRi ** nn_ric
          p_avm(ji, jj, jk) = MAX(zav, avmb(jk)) * wmask(ji, jj, jk)
          p_avt(ji, jj, jk) = MAX(zav * zcfRi, avtb(jk)) * wmask(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    IF (ln_mldw) THEN
      !$ACC KERNELS
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zustar = SQRT(taum(ji, jj) * r1_rau0)
          zhek = rn_ekmfc * zustar / (ABS(ff_t(ji, jj)) + rsmall)
          zh_ekm(ji, jj) = MAX(rn_mldmin, MIN(zhek, rn_mldmax))
        END DO
      END DO
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            IF (pdept(ji, jj, jk) < zh_ekm(ji, jj)) THEN
              p_avm(ji, jj, jk) = MAX(p_avm(ji, jj, jk), rn_wvmix) * wmask(ji, jj, jk)
              p_avt(ji, jj, jk) = MAX(p_avt(ji, jj, jk), rn_wtmix) * wmask(ji, jj, jk)
            END IF
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END IF
  END SUBROUTINE zdf_ric
  SUBROUTINE ric_rst(kt, cdrw)
    INTEGER, INTENT(IN) :: kt
    CHARACTER(LEN = *), INTENT(IN) :: cdrw
    INTEGER :: jit, jk
    INTEGER :: id1, id2
    IF (TRIM(cdrw) == 'READ') THEN
      IF (ln_rstart) THEN
        id1 = iom_varid(numror, 'avt_k', ldstop = .FALSE.)
        id2 = iom_varid(numror, 'avm_k', ldstop = .FALSE.)
        IF (MIN(id1, id2) > 0) THEN
          CALL iom_get(numror, jpdom_autoglo, 'avt_k', avt_k, ldxios = lrxios)
          CALL iom_get(numror, jpdom_autoglo, 'avm_k', avm_k, ldxios = lrxios)
        END IF
      END IF
    ELSE IF (TRIM(cdrw) == 'WRITE') THEN
      IF (lwp) WRITE(numout, FMT = *) '---- ric-rst ----'
      IF (lwxios) CALL iom_swap(cwxios_context)
      CALL iom_rstput(kt, nitrst, numrow, 'avt_k', avt_k, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'avm_k', avm_k, ldxios = lwxios)
      IF (lwxios) CALL iom_swap(cxios_context)
    END IF
  END SUBROUTINE ric_rst
END MODULE zdfric