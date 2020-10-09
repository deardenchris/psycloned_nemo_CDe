MODULE zdftke
  USE oce
  USE phycst
  USE dom_oce
  USE domvvl
  USE sbc_oce
  USE zdfdrg
  USE zdfmxl
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lbclnk
  USE prtctl
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: zdf_tke
  PUBLIC :: zdf_tke_init
  PUBLIC :: tke_rst
  LOGICAL :: ln_mxl0
  INTEGER :: nn_mxl
  REAL(KIND = wp) :: rn_mxl0
  INTEGER :: nn_pdl
  REAL(KIND = wp) :: rn_ediff
  REAL(KIND = wp) :: rn_ediss
  REAL(KIND = wp) :: rn_ebb
  REAL(KIND = wp) :: rn_emin
  REAL(KIND = wp) :: rn_emin0
  REAL(KIND = wp) :: rn_bshear
  LOGICAL :: ln_drg
  INTEGER :: nn_etau
  INTEGER :: nn_htau
  REAL(KIND = wp) :: rn_efr
  REAL(KIND = wp) :: rn_eice
  LOGICAL :: ln_lc
  REAL(KIND = wp) :: rn_lc
  REAL(KIND = wp) :: ri_cri
  REAL(KIND = wp) :: rmxl_min
  REAL(KIND = wp) :: rhftau_add = 1.E-3_wp
  REAL(KIND = wp) :: rhftau_scl = 1.0_wp
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: htau
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: dissl
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: apdlr
  CONTAINS
  INTEGER FUNCTION zdf_tke_alloc()
    ALLOCATE(htau(jpi, jpj), dissl(jpi, jpj, jpk), apdlr(jpi, jpj, jpk), STAT = zdf_tke_alloc)
    CALL mpp_sum('zdftke', zdf_tke_alloc)
    IF (zdf_tke_alloc /= 0) CALL ctl_stop('STOP', 'zdf_tke_alloc: failed to allocate arrays')
  END FUNCTION zdf_tke_alloc
  SUBROUTINE zdf_tke(kt, p_sh2, p_avm, p_avt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: p_sh2
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: p_avm, p_avt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('zdf_tke', 'r0', 0, 0)
    CALL tke_tke(gdepw_n, e3t_n, e3w_n, p_sh2, p_avm, p_avt)
    CALL tke_avn(gdepw_n, e3t_n, e3w_n, p_avm, p_avt)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE zdf_tke
  SUBROUTINE tke_tke(pdepw, p_e3t, p_e3w, p_sh2, p_avm, p_avt)
    USE zdf_oce, ONLY: en
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pdepw
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: p_e3t, p_e3w
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: p_sh2
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: p_avm, p_avt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zetop, zebot, zmsku, zmskv
    REAL(KIND = wp) :: zrhoa = 1.22
    REAL(KIND = wp) :: zcdrag = 1.5E-3
    REAL(KIND = wp) :: zbbrau, zri
    REAL(KIND = wp) :: zfact1, zfact2, zfact3
    REAL(KIND = wp) :: ztx2, zty2, zcof
    REAL(KIND = wp) :: ztau, zdif
    REAL(KIND = wp) :: zus, zwlc, zind
    REAL(KIND = wp) :: zzd_up, zzd_lw
    INTEGER, DIMENSION(jpi, jpj) :: imlc
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zhlc, zfr_i
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zpelc, zdiag, zd_up, zd_lw
    !$ACC KERNELS
    zbbrau = rn_ebb / rau0
    zfact1 = - .5_wp * rdt
    zfact2 = 1.5_wp * rdt * rn_ediss
    zfact3 = 0.5_wp * rn_ediss
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        en(ji, jj, 1) = MAX(rn_emin0, zbbrau * taum(ji, jj)) * tmask(ji, jj, 1)
      END DO
    END DO
    !$ACC END KERNELS
    IF (ln_isfcav) THEN
      !$ACC KERNELS
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          en(ji, jj, mikt(ji, jj)) = rn_emin * tmask(ji, jj, 1)
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    IF (ln_drg) THEN
      !$ACC KERNELS
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zmsku = (2. - umask(ji - 1, jj, mbkt(ji, jj)) * umask(ji, jj, mbkt(ji, jj)))
          zmskv = (2. - vmask(ji, jj - 1, mbkt(ji, jj)) * vmask(ji, jj, mbkt(ji, jj)))
          zebot = - 0.001875_wp * rCdU_bot(ji, jj) * SQRT((zmsku * (ub(ji, jj, mbkt(ji, jj)) + ub(ji - 1, jj, mbkt(ji, jj)))) ** 2 &
&+ (zmskv * (vb(ji, jj, mbkt(ji, jj)) + vb(ji, jj - 1, mbkt(ji, jj)))) ** 2)
          en(ji, jj, mbkt(ji, jj) + 1) = MAX(zebot, rn_emin) * ssmask(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
      IF (ln_isfcav) THEN
        !$ACC KERNELS
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zmsku = (2. - umask(ji - 1, jj, mikt(ji, jj)) * umask(ji, jj, mikt(ji, jj)))
            zmskv = (2. - vmask(ji, jj - 1, mikt(ji, jj)) * vmask(ji, jj, mikt(ji, jj)))
            zetop = - 0.001875_wp * rCdU_top(ji, jj) * SQRT((zmsku * (ub(ji, jj, mikt(ji, jj)) + ub(ji - 1, jj, mikt(ji, jj)))) ** &
&2 + (zmskv * (vb(ji, jj, mikt(ji, jj)) + vb(ji, jj - 1, mikt(ji, jj)))) ** 2)
            en(ji, jj, mikt(ji, jj)) = MAX(zetop, rn_emin) * (1._wp - tmask(ji, jj, 1))
          END DO
        END DO
        !$ACC END KERNELS
      END IF
    END IF
    IF (ln_lc) THEN
      !$ACC KERNELS
      zpelc(:, :, 1) = MAX(rn2b(:, :, 1), 0._wp) * pdepw(:, :, 1) * p_e3w(:, :, 1)
      DO jk = 2, jpk
        zpelc(:, :, jk) = zpelc(:, :, jk - 1) + MAX(rn2b(:, :, jk), 0._wp) * pdepw(:, :, jk) * p_e3w(:, :, jk)
      END DO
      zcof = 0.5 * 0.016 * 0.016 / (zrhoa * zcdrag)
      imlc(:, :) = mbkt(:, :) + 1
      DO jk = jpkm1, 2, - 1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            zus = zcof * taum(ji, jj)
            IF (zpelc(ji, jj, jk) > zus) imlc(ji, jj) = jk
          END DO
        END DO
      END DO
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          zhlc(ji, jj) = pdepw(ji, jj, imlc(ji, jj))
        END DO
      END DO
      zcof = 0.016 / SQRT(zrhoa * zcdrag)
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zus = zcof * SQRT(taum(ji, jj))
          zfr_i(ji, jj) = (1._wp - 4._wp * fr_i(ji, jj)) * zus * zus * zus * tmask(ji, jj, 1)
          IF (zfr_i(ji, jj) < 0.) zfr_i(ji, jj) = 0.
        END DO
      END DO
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            IF (zfr_i(ji, jj) /= 0.) THEN
              IF (pdepw(ji, jj, jk) - zhlc(ji, jj) < 0 .AND. wmask(ji, jj, jk) /= 0.) THEN
                zwlc = rn_lc * SIN(rpi * pdepw(ji, jj, jk) / zhlc(ji, jj))
                en(ji, jj, jk) = en(ji, jj, jk) + rdt * zfr_i(ji, jj) * (zwlc * zwlc * zwlc) / zhlc(ji, jj)
              END IF
            END IF
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    IF (nn_pdl == 1) THEN
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zri = MAX(rn2b(ji, jj, jk), 0._wp) * p_avm(ji, jj, jk) / (p_sh2(ji, jj, jk) + rn_bshear)
            apdlr(ji, jj, jk) = MAX(0.1_wp, ri_cri / MAX(ri_cri, zri))
          END DO
        END DO
      END DO
    END IF
    DO jk = 2, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zcof = zfact1 * tmask(ji, jj, jk)
          zzd_up = zcof * MAX(p_avm(ji, jj, jk + 1) + p_avm(ji, jj, jk), 2.E-5_wp) / (p_e3t(ji, jj, jk) * p_e3w(ji, jj, jk))
          zzd_lw = zcof * MAX(p_avm(ji, jj, jk) + p_avm(ji, jj, jk - 1), 2.E-5_wp) / (p_e3t(ji, jj, jk - 1) * p_e3w(ji, jj, jk))
          zd_up(ji, jj, jk) = zzd_up
          zd_lw(ji, jj, jk) = zzd_lw
          zdiag(ji, jj, jk) = 1._wp - zzd_lw - zzd_up + zfact2 * dissl(ji, jj, jk) * wmask(ji, jj, jk)
          en(ji, jj, jk) = en(ji, jj, jk) + rdt * (p_sh2(ji, jj, jk) - p_avt(ji, jj, jk) * rn2(ji, jj, jk) + zfact3 * dissl(ji, &
&jj, jk) * en(ji, jj, jk)) * wmask(ji, jj, jk)
        END DO
      END DO
    END DO
    DO jk = 3, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zdiag(ji, jj, jk) = zdiag(ji, jj, jk) - zd_lw(ji, jj, jk) * zd_up(ji, jj, jk - 1) / zdiag(ji, jj, jk - 1)
        END DO
      END DO
    END DO
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zd_lw(ji, jj, 2) = en(ji, jj, 2) - zd_lw(ji, jj, 2) * en(ji, jj, 1)
      END DO
    END DO
    DO jk = 3, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zd_lw(ji, jj, jk) = en(ji, jj, jk) - zd_lw(ji, jj, jk) / zdiag(ji, jj, jk - 1) * zd_lw(ji, jj, jk - 1)
        END DO
      END DO
    END DO
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        en(ji, jj, jpkm1) = zd_lw(ji, jj, jpkm1) / zdiag(ji, jj, jpkm1)
      END DO
    END DO
    DO jk = jpk - 2, 2, - 1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          en(ji, jj, jk) = (zd_lw(ji, jj, jk) - zd_up(ji, jj, jk) * en(ji, jj, jk + 1)) / zdiag(ji, jj, jk)
        END DO
      END DO
    END DO
    DO jk = 2, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          en(ji, jj, jk) = MAX(en(ji, jj, jk), rn_emin) * wmask(ji, jj, jk)
        END DO
      END DO
    END DO
    IF (nn_etau == 1) THEN
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            en(ji, jj, jk) = en(ji, jj, jk) + rn_efr * en(ji, jj, 1) * EXP(- pdepw(ji, jj, jk) / htau(ji, jj)) * MAX(0., 1._wp - &
&rn_eice * fr_i(ji, jj)) * wmask(ji, jj, jk) * tmask(ji, jj, 1)
          END DO
        END DO
      END DO
    ELSE IF (nn_etau == 2) THEN
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          jk = nmln(ji, jj)
          en(ji, jj, jk) = en(ji, jj, jk) + rn_efr * en(ji, jj, 1) * EXP(- pdepw(ji, jj, jk) / htau(ji, jj)) * MAX(0., 1._wp - &
&rn_eice * fr_i(ji, jj)) * wmask(ji, jj, jk) * tmask(ji, jj, 1)
        END DO
      END DO
    ELSE IF (nn_etau == 3) THEN
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            ztx2 = utau(ji - 1, jj) + utau(ji, jj)
            zty2 = vtau(ji, jj - 1) + vtau(ji, jj)
            ztau = 0.5_wp * SQRT(ztx2 * ztx2 + zty2 * zty2) * tmask(ji, jj, 1)
            zdif = taum(ji, jj) - ztau
            zdif = rhftau_scl * MAX(0._wp, zdif + rhftau_add)
            en(ji, jj, jk) = en(ji, jj, jk) + zbbrau * zdif * EXP(- pdepw(ji, jj, jk) / htau(ji, jj)) * MAX(0., 1._wp - rn_eice * &
&fr_i(ji, jj)) * wmask(ji, jj, jk) * tmask(ji, jj, 1)
          END DO
        END DO
      END DO
    END IF
    !$ACC END KERNELS
  END SUBROUTINE tke_tke
  SUBROUTINE tke_avn(pdepw, p_e3t, p_e3w, p_avm, p_avt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE zdf_oce, ONLY: en, avtb, avmb, avtb_2d
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pdepw
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: p_e3t, p_e3w
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: p_avm, p_avt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zrn2, zraug, zcoef, zav
    REAL(KIND = wp) :: zdku, zdkv, zsqen
    REAL(KIND = wp) :: zemxl, zemlm, zemlp
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zmxlm, zmxld
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    zmxlm(:, :, :) = rmxl_min
    zmxld(:, :, :) = rmxl_min
    !$ACC END KERNELS
    IF (ln_mxl0) THEN
      !$ACC KERNELS
      zraug = vkarmn * 2.E5_wp / (rau0 * grav)
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zmxlm(ji, jj, 1) = MAX(rn_mxl0, zraug * taum(ji, jj) * tmask(ji, jj, 1))
        END DO
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      zmxlm(:, :, 1) = rn_mxl0
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    DO jk = 2, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zrn2 = MAX(rn2(ji, jj, jk), rsmall)
          zmxlm(ji, jj, jk) = MAX(rmxl_min, SQRT(2._wp * en(ji, jj, jk) / zrn2))
        END DO
      END DO
    END DO
    zmxld(:, :, 1) = zmxlm(:, :, 1)
    zmxld(:, :, jpk) = rmxl_min
    SELECT CASE (nn_mxl)
    CASE (0)
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zemxl = MIN(pdepw(ji, jj, jk) - pdepw(ji, jj, mikt(ji, jj)), zmxlm(ji, jj, jk), pdepw(ji, jj, mbkt(ji, jj) + 1) - &
&pdepw(ji, jj, jk))
            zmxlm(ji, jj, jk) = zemxl * wmask(ji, jj, jk) + MIN(zmxlm(ji, jj, jk), p_e3w(ji, jj, jk)) * (1 - wmask(ji, jj, jk))
            zmxld(ji, jj, jk) = zemxl * wmask(ji, jj, jk) + MIN(zmxlm(ji, jj, jk), p_e3w(ji, jj, jk)) * (1 - wmask(ji, jj, jk))
          END DO
        END DO
      END DO
    CASE (1)
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zemxl = MIN(p_e3w(ji, jj, jk), zmxlm(ji, jj, jk))
            zmxlm(ji, jj, jk) = zemxl
            zmxld(ji, jj, jk) = zemxl
          END DO
        END DO
      END DO
    CASE (2)
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zmxlm(ji, jj, jk) = MIN(zmxlm(ji, jj, jk - 1) + p_e3t(ji, jj, jk - 1), zmxlm(ji, jj, jk))
          END DO
        END DO
      END DO
      DO jk = jpkm1, 2, - 1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zemxl = MIN(zmxlm(ji, jj, jk + 1) + p_e3t(ji, jj, jk + 1), zmxlm(ji, jj, jk))
            zmxlm(ji, jj, jk) = zemxl
            zmxld(ji, jj, jk) = zemxl
          END DO
        END DO
      END DO
    CASE (3)
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zmxld(ji, jj, jk) = MIN(zmxld(ji, jj, jk - 1) + p_e3t(ji, jj, jk - 1), zmxlm(ji, jj, jk))
          END DO
        END DO
      END DO
      DO jk = jpkm1, 2, - 1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zmxlm(ji, jj, jk) = MIN(zmxlm(ji, jj, jk + 1) + p_e3t(ji, jj, jk + 1), zmxlm(ji, jj, jk))
          END DO
        END DO
      END DO
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zemlm = MIN(zmxld(ji, jj, jk), zmxlm(ji, jj, jk))
            zemlp = SQRT(zmxld(ji, jj, jk) * zmxlm(ji, jj, jk))
            zmxlm(ji, jj, jk) = zemlm
            zmxld(ji, jj, jk) = zemlp
          END DO
        END DO
      END DO
    END SELECT
    DO jk = 1, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zsqen = SQRT(en(ji, jj, jk))
          zav = rn_ediff * zmxlm(ji, jj, jk) * zsqen
          p_avm(ji, jj, jk) = MAX(zav, avmb(jk)) * wmask(ji, jj, jk)
          p_avt(ji, jj, jk) = MAX(zav, avtb_2d(ji, jj) * avtb(jk)) * wmask(ji, jj, jk)
          dissl(ji, jj, jk) = zsqen / zmxld(ji, jj, jk)
        END DO
      END DO
    END DO
    IF (nn_pdl == 1) THEN
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            p_avt(ji, jj, jk) = MAX(apdlr(ji, jj, jk) * p_avt(ji, jj, jk), avtb_2d(ji, jj) * avtb(jk)) * tmask(ji, jj, jk)
          END DO
        END DO
      END DO
    END IF
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('tke_avn', 'r0', 0, 0)
    IF (ln_ctl) THEN
      CALL prt_ctl(tab3d_1 = en, clinfo1 = ' tke  - e: ', tab3d_2 = p_avt, clinfo2 = ' t: ', kdim = jpk)
      CALL prt_ctl(tab3d_1 = p_avm, clinfo1 = ' tke  - m: ', kdim = jpk)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE tke_avn
  SUBROUTINE zdf_tke_init
    USE zdf_oce, ONLY: ln_zdfiwm
    INTEGER :: ji, jj, jk
    INTEGER :: ios
    NAMELIST /namzdf_tke/ rn_ediff, rn_ediss, rn_ebb, rn_emin, rn_emin0, rn_bshear, nn_mxl, ln_mxl0, rn_mxl0, nn_pdl, ln_drg, &
&ln_lc, rn_lc, nn_etau, nn_htau, rn_efr, rn_eice
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namzdf_tke, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namzdf_tke in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namzdf_tke, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namzdf_tke in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namzdf_tke)
    ri_cri = 2._wp / (2._wp + rn_ediss / rn_ediff)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'zdf_tke_init : tke turbulent closure scheme - initialisation'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namzdf_tke : set tke mixing parameters'
      WRITE(numout, FMT = *) '      coef. to compute avt                        rn_ediff  = ', rn_ediff
      WRITE(numout, FMT = *) '      Kolmogoroff dissipation coef.               rn_ediss  = ', rn_ediss
      WRITE(numout, FMT = *) '      tke surface input coef.                     rn_ebb    = ', rn_ebb
      WRITE(numout, FMT = *) '      minimum value of tke                        rn_emin   = ', rn_emin
      WRITE(numout, FMT = *) '      surface minimum value of tke                rn_emin0  = ', rn_emin0
      WRITE(numout, FMT = *) '      prandl number flag                          nn_pdl    = ', nn_pdl
      WRITE(numout, FMT = *) '      background shear (>0)                       rn_bshear = ', rn_bshear
      WRITE(numout, FMT = *) '      mixing length type                          nn_mxl    = ', nn_mxl
      WRITE(numout, FMT = *) '         surface mixing length = F(stress) or not    ln_mxl0   = ', ln_mxl0
      WRITE(numout, FMT = *) '         surface  mixing length minimum value        rn_mxl0   = ', rn_mxl0
      WRITE(numout, FMT = *) '      top/bottom friction forcing flag            ln_drg    = ', ln_drg
      WRITE(numout, FMT = *) '      Langmuir cells parametrization              ln_lc     = ', ln_lc
      WRITE(numout, FMT = *) '         coef to compute vertical velocity of LC     rn_lc  = ', rn_lc
      WRITE(numout, FMT = *) '      test param. to add tke induced by wind      nn_etau   = ', nn_etau
      WRITE(numout, FMT = *) '          type of tke penetration profile            nn_htau   = ', nn_htau
      WRITE(numout, FMT = *) '          fraction of TKE that penetrates            rn_efr    = ', rn_efr
      WRITE(numout, FMT = *) '          below sea-ice:  =0 ON                      rn_eice   = ', rn_eice
      WRITE(numout, FMT = *) '          =4 OFF when ice fraction > 1/4   '
      IF (ln_drg) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) '   Namelist namdrg_top/_bot:   used values:'
        WRITE(numout, FMT = *) '      top    ocean cavity roughness (m)          rn_z0(_top)= ', r_z0_top
        WRITE(numout, FMT = *) '      Bottom seafloor     roughness (m)          rn_z0(_bot)= ', r_z0_bot
      END IF
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) '   ==>>>   critical Richardson nb with your parameters  ri_cri = ', ri_cri
      WRITE(numout, FMT = *)
    END IF
    IF (ln_zdfiwm) THEN
      rn_emin = 1.E-10_wp
      rmxl_min = 1.E-03_wp
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   Internal wave-driven mixing case:   force   rn_emin = 1.e-10 and rmxl_min = 1.e-3'
    ELSE
      rmxl_min = 1.E-6_wp / (rn_ediff * SQRT(rn_emin))
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   minimum mixing length with your parameters rmxl_min = ', rmxl_min
    END IF
    IF (zdf_tke_alloc() /= 0) CALL ctl_stop('STOP', 'zdf_tke_init : unable to allocate arrays')
    IF (nn_mxl < 0 .OR. nn_mxl > 3) CALL ctl_stop('bad flag: nn_mxl is  0, 1 or 2 ')
    IF (nn_pdl < 0 .OR. nn_pdl > 1) CALL ctl_stop('bad flag: nn_pdl is  0 or 1    ')
    IF ((nn_htau < 0 .OR. nn_htau > 1) .AND. nn_htau .NE. 4) CALL ctl_stop('bad flag: nn_htau is 0, 1 or 4 ')
    IF (nn_etau == 3 .AND. .NOT. ln_cpl) CALL ctl_stop('nn_etau == 3 : HF taum only known in coupled mode')
    IF (ln_mxl0) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   use a surface mixing length = F(stress) :   set rn_mxl0 = rmxl_min'
      rn_mxl0 = rmxl_min
    END IF
    IF (nn_etau == 2) CALL zdf_mxl(nit000)
    !$ACC KERNELS
    IF (nn_etau /= 0) THEN
      SELECT CASE (nn_htau)
      CASE (0)
        htau(:, :) = 10._wp
      CASE (1)
        htau(:, :) = MAX(0.5_wp, MIN(30._wp, 45._wp * ABS(SIN(rpi / 180._wp * gphit(:, :)))))
      CASE (4)
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            IF (gphit(ji, jj) <= 0._wp) THEN
              htau(ji, jj) = MAX(0.5_wp, MIN(30._wp, 45._wp * ABS(SIN(rpi / 180._wp * gphit(ji, jj)))))
            ELSE
              htau(ji, jj) = MAX(0.5_wp, MIN(10._wp, 45._wp * ABS(SIN(rpi / 180._wp * gphit(ji, jj)))))
            END IF
          END DO
        END DO
      END SELECT
    END IF
    !$ACC END KERNELS
    CALL tke_rst(nit000, 'READ')
    IF (lwxios) THEN
      CALL iom_set_rstw_var_active('en')
      CALL iom_set_rstw_var_active('avt_k')
      CALL iom_set_rstw_var_active('avm_k')
      CALL iom_set_rstw_var_active('dissl')
    END IF
  END SUBROUTINE zdf_tke_init
  SUBROUTINE tke_rst(kt, cdrw)
    USE zdf_oce, ONLY: en, avt_k, avm_k
    INTEGER, INTENT(IN) :: kt
    CHARACTER(LEN = *), INTENT(IN) :: cdrw
    INTEGER :: jit, jk
    INTEGER :: id1, id2, id3, id4
    IF (TRIM(cdrw) == 'READ') THEN
      IF (ln_rstart) THEN
        id1 = iom_varid(numror, 'en', ldstop = .FALSE.)
        id2 = iom_varid(numror, 'avt_k', ldstop = .FALSE.)
        id3 = iom_varid(numror, 'avm_k', ldstop = .FALSE.)
        id4 = iom_varid(numror, 'dissl', ldstop = .FALSE.)
        IF (MIN(id1, id2, id3, id4) > 0) THEN
          CALL iom_get(numror, jpdom_autoglo, 'en', en, ldxios = lrxios)
          CALL iom_get(numror, jpdom_autoglo, 'avt_k', avt_k, ldxios = lrxios)
          CALL iom_get(numror, jpdom_autoglo, 'avm_k', avm_k, ldxios = lrxios)
          CALL iom_get(numror, jpdom_autoglo, 'dissl', dissl, ldxios = lrxios)
        ELSE
          IF (lwp) WRITE(numout, FMT = *)
          IF (lwp) WRITE(numout, FMT = *) '   ==>>>   previous run without TKE scheme, set en to background values'
          !$ACC KERNELS
          en(:, :, :) = rn_emin * wmask(:, :, :)
          dissl(:, :, :) = 1.E-12_wp
          !$ACC END KERNELS
        END IF
      ELSE
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   start from rest: set en to the background value'
        !$ACC KERNELS
        en(:, :, :) = rn_emin * wmask(:, :, :)
        dissl(:, :, :) = 1.E-12_wp
        !$ACC END KERNELS
      END IF
    ELSE IF (TRIM(cdrw) == 'WRITE') THEN
      IF (lwp) WRITE(numout, FMT = *) '---- tke_rst ----'
      IF (lwxios) CALL iom_swap(cwxios_context)
      CALL iom_rstput(kt, nitrst, numrow, 'en', en, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'avt_k', avt_k, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'avm_k', avm_k, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'dissl', dissl, ldxios = lwxios)
      IF (lwxios) CALL iom_swap(cxios_context)
    END IF
  END SUBROUTINE tke_rst
END MODULE zdftke