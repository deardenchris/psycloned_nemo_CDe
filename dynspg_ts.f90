MODULE dynspg_ts
  USE oce
  USE dom_oce
  USE sbc_oce
  USE zdf_oce
  USE zdfdrg
  USE sbcisf
  USE sbcapr
  USE dynadv, ONLY: ln_dynadv_vec
  USE dynvor
  USE phycst
  USE dynvor
  USE wet_dry
  USE bdy_oce
  USE bdyvol
  USE bdytides
  USE bdydyn2d
  USE sbctide
  USE updtide
  USE sbcwave
  USE diatmb
  USE in_out_manager
  USE lib_mpp
  USE lbclnk
  USE prtctl
  USE iom
  USE restart
  USE diatmb
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_spg_ts
  PUBLIC :: dyn_spg_ts_init
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: un_adv, vn_adv
  INTEGER, SAVE :: icycle
  REAL(KIND = wp), SAVE :: rdtbt
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:) :: wgtbtp1, wgtbtp2
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: zwz
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: ftnw, ftne
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: ftsw, ftse
  REAL(KIND = wp) :: r1_12 = 1._wp / 12._wp
  REAL(KIND = wp) :: r1_8 = 0.125_wp
  REAL(KIND = wp) :: r1_4 = 0.25_wp
  REAL(KIND = wp) :: r1_2 = 0.5_wp
  CONTAINS
  INTEGER FUNCTION dyn_spg_ts_alloc()
    INTEGER :: ierr(3)
    !$ACC KERNELS
    ierr(:) = 0
    !$ACC END KERNELS
    ALLOCATE(wgtbtp1(3 * nn_baro), wgtbtp2(3 * nn_baro), zwz(jpi, jpj), STAT = ierr(1))
    IF (ln_dynvor_een .OR. ln_dynvor_eeT) ALLOCATE(ftnw(jpi, jpj), ftne(jpi, jpj), ftsw(jpi, jpj), ftse(jpi, jpj), STAT = ierr(2))
    ALLOCATE(un_adv(jpi, jpj), vn_adv(jpi, jpj), STAT = ierr(3))
    dyn_spg_ts_alloc = MAXVAL(ierr(:))
    CALL mpp_sum('dynspg_ts', dyn_spg_ts_alloc)
    IF (dyn_spg_ts_alloc /= 0) CALL ctl_stop('STOP', 'dyn_spg_ts_alloc: failed to allocate arrays')
  END FUNCTION dyn_spg_ts_alloc
  SUBROUTINE dyn_spg_ts(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk, jn
    LOGICAL :: ll_fw_start
    LOGICAL :: ll_init
    LOGICAL :: ll_tmp1, ll_tmp2
    INTEGER :: ikbu, iktu, noffset
    INTEGER :: ikbv, iktv
    REAL(KIND = wp) :: r1_2dt_b, z2dt_bf
    REAL(KIND = wp) :: zx1, zx2, zu_spg, zhura, z1_hu
    REAL(KIND = wp) :: zy1, zy2, zv_spg, zhvra, z1_hv
    REAL(KIND = wp) :: za0, za1, za2, za3
    REAL(KIND = wp) :: zmdi, zztmp, z1_ht
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zsshp2_e, zhf
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwx, zu_trd, zu_frc, zssh_frc
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwy, zv_trd, zv_frc, zhdiv
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zsshu_a, zhup2_e, zhust_e, zhtp2_e
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zsshv_a, zhvp2_e, zhvst_e
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zCdU_u, zCdU_v
    REAL(KIND = wp) :: zwdramp
    INTEGER :: iwdg, jwdg, kwdg
    REAL(KIND = wp) :: zepsilon, zgamma
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zcpx, zcpy
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: ztwdmask, zuwdmask, zvwdmask
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zuwdav2, zvwdav2
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
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data15
    CALL profile_psy_data0 % PreStart('dyn_spg_ts', 'r0', 0, 0)
    IF (ln_wd_il) ALLOCATE(zcpx(jpi, jpj), zcpy(jpi, jpj))
    IF (ln_wd_dl) ALLOCATE(ztwdmask(jpi, jpj), zuwdmask(jpi, jpj), zvwdmask(jpi, jpj), zuwdav2(jpi, jpj), zvwdav2(jpi, jpj))
    zmdi = 1.E+20
    zwdramp = r_rn_wdmin1
    IF (kt == nit000 .AND. neuler == 0) THEN
      z2dt_bf = rdt
    ELSE
      z2dt_bf = 2.0_wp * rdt
    END IF
    r1_2dt_b = 1.0_wp / z2dt_bf
    ll_init = ln_bt_av
    ll_fw_start = .FALSE.
    IF (.NOT. ln_bt_fw) THEN
      noffset = - nn_baro
    ELSE
      noffset = 0
    END IF
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'dyn_spg_ts : surface pressure gradient trend'
      IF (lwp) WRITE(numout, *) '~~~~~~~~~~   free surface with time splitting'
      IF (lwp) WRITE(numout, *)
      IF (neuler == 0) ll_init = .TRUE.
      IF (ln_bt_fw .OR. neuler == 0) THEN
        ll_fw_start = .TRUE.
        noffset = 0
      ELSE
        ll_fw_start = .FALSE.
      END IF
      CALL ts_wgt(ln_bt_av, ll_fw_start, icycle, wgtbtp1, wgtbtp2)
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (ln_isfcav) THEN
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zCdU_u(ji, jj) = r1_2 * (rCdU_bot(ji + 1, jj) + rCdU_bot(ji, jj) + rCdU_top(ji + 1, jj) + rCdU_top(ji, jj))
          zCdU_v(ji, jj) = r1_2 * (rCdU_bot(ji, jj + 1) + rCdU_bot(ji, jj) + rCdU_top(ji, jj + 1) + rCdU_top(ji, jj))
        END DO
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zCdU_u(ji, jj) = r1_2 * (rCdU_bot(ji + 1, jj) + rCdU_bot(ji, jj))
          zCdU_v(ji, jj) = r1_2 * (rCdU_bot(ji, jj + 1) + rCdU_bot(ji, jj))
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    IF (kt == nit000 .OR. .NOT. ln_linssh) THEN
      SELECT CASE (nvor_scheme)
      CASE (np_EEN)
        SELECT CASE (nn_een_e3f)
        CASE (0)
          !$ACC KERNELS
          !$ACC loop independent collapse(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              zwz(ji, jj) = (ht_n(ji, jj + 1) + ht_n(ji + 1, jj + 1) + ht_n(ji, jj) + ht_n(ji + 1, jj)) * 0.25_wp
              IF (zwz(ji, jj) /= 0._wp) zwz(ji, jj) = ff_f(ji, jj) / zwz(ji, jj)
            END DO
          END DO
          !$ACC END KERNELS
        CASE (1)
          !$ACC KERNELS
          !$ACC loop independent collapse(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              zwz(ji, jj) = (ht_n(ji, jj + 1) + ht_n(ji + 1, jj + 1) + ht_n(ji, jj) + ht_n(ji + 1, jj)) / (MAX(1._wp, ssmask(ji, jj + 1) + ssmask(ji + 1, jj + 1) + ssmask(ji, jj) + ssmask(ji + 1, jj)))
              IF (zwz(ji, jj) /= 0._wp) zwz(ji, jj) = ff_f(ji, jj) / zwz(ji, jj)
            END DO
          END DO
          !$ACC END KERNELS
        END SELECT
        CALL profile_psy_data1 % PreStart('dyn_spg_ts', 'r1', 0, 0)
        CALL lbc_lnk('dynspg_ts', zwz, 'F', 1._wp)
        CALL profile_psy_data1 % PostEnd
        !$ACC KERNELS
        ftne(1, :) = 0._wp
        ftnw(1, :) = 0._wp
        ftse(1, :) = 0._wp
        ftsw(1, :) = 0._wp
        !$ACC loop independent collapse(2)
        DO jj = 2, jpj
          DO ji = 2, jpi
            ftne(ji, jj) = zwz(ji - 1, jj) + zwz(ji, jj) + zwz(ji, jj - 1)
            ftnw(ji, jj) = zwz(ji - 1, jj - 1) + zwz(ji - 1, jj) + zwz(ji, jj)
            ftse(ji, jj) = zwz(ji, jj) + zwz(ji, jj - 1) + zwz(ji - 1, jj - 1)
            ftsw(ji, jj) = zwz(ji, jj - 1) + zwz(ji - 1, jj - 1) + zwz(ji - 1, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_EET)
        !$ACC KERNELS
        ftne(1, :) = 0._wp
        ftnw(1, :) = 0._wp
        ftse(1, :) = 0._wp
        ftsw(1, :) = 0._wp
        !$ACC loop independent collapse(2)
        DO jj = 2, jpj
          DO ji = 2, jpi
            z1_ht = ssmask(ji, jj) / (ht_n(ji, jj) + 1._wp - ssmask(ji, jj))
            ftne(ji, jj) = (ff_f(ji - 1, jj) + ff_f(ji, jj) + ff_f(ji, jj - 1)) * z1_ht
            ftnw(ji, jj) = (ff_f(ji - 1, jj - 1) + ff_f(ji - 1, jj) + ff_f(ji, jj)) * z1_ht
            ftse(ji, jj) = (ff_f(ji, jj) + ff_f(ji, jj - 1) + ff_f(ji - 1, jj - 1)) * z1_ht
            ftsw(ji, jj) = (ff_f(ji, jj - 1) + ff_f(ji - 1, jj - 1) + ff_f(ji - 1, jj)) * z1_ht
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_ENE, np_ENS, np_MIX)
        !$ACC KERNELS
        zwz(:, :) = 0._wp
        zhf(:, :) = 0._wp
        !$ACC END KERNELS
        IF (.NOT. ln_sco) THEN
        ELSE
          !$ACC KERNELS
          !$ACC loop independent collapse(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              zhf(ji, jj) = (ht_0(ji, jj) + ht_0(ji + 1, jj) + ht_0(ji, jj + 1) + ht_0(ji + 1, jj + 1)) / MAX(ssmask(ji, jj) + ssmask(ji + 1, jj) + ssmask(ji, jj + 1) + ssmask(ji + 1, jj + 1), 1._wp)
            END DO
          END DO
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpi
            zhf(ji, jj) = zhf(ji, jj) * (1._wp - umask(ji, jj, 1) * umask(ji, jj + 1, 1))
          END DO
        END DO
        DO jk = 1, jpkm1
          !$ACC loop independent collapse(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpi
              zhf(ji, jj) = zhf(ji, jj) + e3f_n(ji, jj, jk) * umask(ji, jj, jk) * umask(ji, jj + 1, jk)
            END DO
          END DO
        END DO
        !$ACC END KERNELS
        CALL profile_psy_data2 % PreStart('dyn_spg_ts', 'r2', 0, 0)
        CALL lbc_lnk('dynspg_ts', zhf, 'F', 1._wp)
        CALL profile_psy_data2 % PostEnd
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            IF (zhf(ji, jj) /= 0._wp) zwz(ji, jj) = 1._wp / zhf(ji, jj)
          END DO
        END DO
        zwz(:, :) = ff_f(:, :) * zwz(:, :)
        !$ACC END KERNELS
      END SELECT
    END IF
    CALL profile_psy_data3 % PreStart('dyn_spg_ts', 'r3', 0, 0)
    IF (.NOT. ln_bt_fw .AND. (neuler == 0 .AND. kt == nit000 + 1)) THEN
      ll_fw_start = .FALSE.
      CALL ts_wgt(ln_bt_av, ll_fw_start, icycle, wgtbtp1, wgtbtp2)
    END IF
    CALL profile_psy_data3 % PostEnd
    !$ACC KERNELS
    zu_frc(:, :) = 0._wp
    zv_frc(:, :) = 0._wp
    DO jk = 1, jpkm1
      zu_frc(:, :) = zu_frc(:, :) + e3u_n(:, :, jk) * ua(:, :, jk) * umask(:, :, jk)
      zv_frc(:, :) = zv_frc(:, :) + e3v_n(:, :, jk) * va(:, :, jk) * vmask(:, :, jk)
    END DO
    zu_frc(:, :) = zu_frc(:, :) * r1_hu_n(:, :)
    zv_frc(:, :) = zv_frc(:, :) * r1_hv_n(:, :)
    DO jk = 1, jpkm1
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ua(ji, jj, jk) = ua(ji, jj, jk) - zu_frc(ji, jj) * umask(ji, jj, jk)
          va(ji, jj, jk) = va(ji, jj, jk) - zv_frc(ji, jj) * vmask(ji, jj, jk)
        END DO
      END DO
    END DO
    zwx(:, :) = un_b(:, :) * hu_n(:, :) * e2u(:, :)
    zwy(:, :) = vn_b(:, :) * hv_n(:, :) * e1v(:, :)
    SELECT CASE (nvor_scheme)
    CASE (np_ENT)
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zu_trd(ji, jj) = + r1_4 * r1_e1e2u(ji, jj) * r1_hu_n(ji, jj) * (e1e2t(ji + 1, jj) * ht_n(ji + 1, jj) * ff_t(ji + 1, jj) * (vn_b(ji + 1, jj) + vn_b(ji + 1, jj - 1)) + e1e2t(ji, jj) * ht_n(ji, jj) * ff_t(ji, jj) * (vn_b(ji, jj) + vn_b(ji, jj - 1)))
          zv_trd(ji, jj) = - r1_4 * r1_e1e2v(ji, jj) * r1_hv_n(ji, jj) * (e1e2t(ji, jj + 1) * ht_n(ji, jj + 1) * ff_t(ji, jj + 1) * (un_b(ji, jj + 1) + un_b(ji - 1, jj + 1)) + e1e2t(ji, jj) * ht_n(ji, jj) * ff_t(ji, jj) * (un_b(ji, jj) + un_b(ji - 1, jj)))
        END DO
      END DO
    CASE (np_ENE, np_MIX)
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zy1 = (zwy(ji, jj - 1) + zwy(ji + 1, jj - 1)) * r1_e1u(ji, jj)
          zy2 = (zwy(ji, jj) + zwy(ji + 1, jj)) * r1_e1u(ji, jj)
          zx1 = (zwx(ji - 1, jj) + zwx(ji - 1, jj + 1)) * r1_e2v(ji, jj)
          zx2 = (zwx(ji, jj) + zwx(ji, jj + 1)) * r1_e2v(ji, jj)
          zu_trd(ji, jj) = r1_4 * (zwz(ji, jj - 1) * zy1 + zwz(ji, jj) * zy2)
          zv_trd(ji, jj) = - r1_4 * (zwz(ji - 1, jj) * zx1 + zwz(ji, jj) * zx2)
        END DO
      END DO
    CASE (np_ENS)
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zy1 = r1_8 * (zwy(ji, jj - 1) + zwy(ji + 1, jj - 1) + zwy(ji, jj) + zwy(ji + 1, jj)) * r1_e1u(ji, jj)
          zx1 = - r1_8 * (zwx(ji - 1, jj) + zwx(ji - 1, jj + 1) + zwx(ji, jj) + zwx(ji, jj + 1)) * r1_e2v(ji, jj)
          zu_trd(ji, jj) = zy1 * (zwz(ji, jj - 1) + zwz(ji, jj))
          zv_trd(ji, jj) = zx1 * (zwz(ji - 1, jj) + zwz(ji, jj))
        END DO
      END DO
    CASE (np_EET, np_EEN)
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zu_trd(ji, jj) = + r1_12 * r1_e1u(ji, jj) * (ftne(ji, jj) * zwy(ji, jj) + ftnw(ji + 1, jj) * zwy(ji + 1, jj) + ftse(ji, jj) * zwy(ji, jj - 1) + ftsw(ji + 1, jj) * zwy(ji + 1, jj - 1))
          zv_trd(ji, jj) = - r1_12 * r1_e2v(ji, jj) * (ftsw(ji, jj + 1) * zwx(ji - 1, jj + 1) + ftse(ji, jj + 1) * zwx(ji, jj + 1) + ftnw(ji, jj) * zwx(ji - 1, jj) + ftne(ji, jj) * zwx(ji, jj))
        END DO
      END DO
    END SELECT
    !$ACC END KERNELS
    IF (.NOT. ln_linssh) THEN
      IF (ln_wd_il) THEN
        CALL profile_psy_data4 % PreStart('dyn_spg_ts', 'r4', 0, 0)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            ll_tmp1 = MIN(sshn(ji, jj), sshn(ji + 1, jj)) > MAX(- ht_0(ji, jj), - ht_0(ji + 1, jj)) .AND. MAX(sshn(ji, jj) + ht_0(ji, jj), sshn(ji + 1, jj) + ht_0(ji + 1, jj)) > rn_wdmin1 + rn_wdmin2
            ll_tmp2 = (ABS(sshn(ji + 1, jj) - sshn(ji, jj)) > 1.E-12) .AND. (MAX(sshn(ji, jj), sshn(ji + 1, jj)) > MAX(- ht_0(ji, jj), - ht_0(ji + 1, jj)) + rn_wdmin1 + rn_wdmin2)
            IF (ll_tmp1) THEN
              zcpx(ji, jj) = 1.0_wp
            ELSE IF (ll_tmp2) THEN
              zcpx(ji, jj) = ABS((sshn(ji + 1, jj) + ht_0(ji + 1, jj) - sshn(ji, jj) - ht_0(ji, jj)) / (sshn(ji + 1, jj) - sshn(ji, jj)))
              zcpx(ji, jj) = MAX(MIN(zcpx(ji, jj), 1.0_wp), 0.0_wp)
            ELSE
              zcpx(ji, jj) = 0._wp
            END IF
            ll_tmp1 = MIN(sshn(ji, jj), sshn(ji, jj + 1)) > MAX(- ht_0(ji, jj), - ht_0(ji, jj + 1)) .AND. MAX(sshn(ji, jj) + ht_0(ji, jj), sshn(ji, jj + 1) + ht_0(ji, jj + 1)) > rn_wdmin1 + rn_wdmin2
            ll_tmp2 = (ABS(sshn(ji, jj) - sshn(ji, jj + 1)) > 1.E-12) .AND. (MAX(sshn(ji, jj), sshn(ji, jj + 1)) > MAX(- ht_0(ji, jj), - ht_0(ji, jj + 1)) + rn_wdmin1 + rn_wdmin2)
            IF (ll_tmp1) THEN
              zcpy(ji, jj) = 1.0_wp
            ELSE IF (ll_tmp2) THEN
              zcpy(ji, jj) = ABS((sshn(ji, jj + 1) + ht_0(ji, jj + 1) - sshn(ji, jj) - ht_0(ji, jj)) / (sshn(ji, jj + 1) - sshn(ji, jj)))
              zcpy(ji, jj) = MAX(0._wp, MIN(zcpy(ji, jj), 1.0_wp))
            ELSE
              zcpy(ji, jj) = 0._wp
            END IF
          END DO
        END DO
        CALL profile_psy_data4 % PostEnd
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zu_trd(ji, jj) = zu_trd(ji, jj) - grav * (sshn(ji + 1, jj) - sshn(ji, jj)) * r1_e1u(ji, jj) * zcpx(ji, jj) * wdrampu(ji, jj)
            zv_trd(ji, jj) = zv_trd(ji, jj) - grav * (sshn(ji, jj + 1) - sshn(ji, jj)) * r1_e2v(ji, jj) * zcpy(ji, jj) * wdrampv(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zu_trd(ji, jj) = zu_trd(ji, jj) - grav * (sshn(ji + 1, jj) - sshn(ji, jj)) * r1_e1u(ji, jj)
            zv_trd(ji, jj) = zv_trd(ji, jj) - grav * (sshn(ji, jj + 1) - sshn(ji, jj)) * r1_e2v(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      END IF
    END IF
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zu_frc(ji, jj) = zu_frc(ji, jj) - zu_trd(ji, jj) * ssumask(ji, jj)
        zv_frc(ji, jj) = zv_frc(ji, jj) - zv_trd(ji, jj) * ssvmask(ji, jj)
      END DO
    END DO
    !$ACC END KERNELS
    IF (ln_bt_fw) THEN
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ikbu = mbku(ji, jj)
          ikbv = mbkv(ji, jj)
          zwx(ji, jj) = un(ji, jj, ikbu) - un_b(ji, jj)
          zwy(ji, jj) = vn(ji, jj, ikbv) - vn_b(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ikbu = mbku(ji, jj)
          ikbv = mbkv(ji, jj)
          zwx(ji, jj) = ub(ji, jj, ikbu) - ub_b(ji, jj)
          zwy(ji, jj) = vb(ji, jj, ikbv) - vb_b(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    IF (ln_wd_il) THEN
      !$ACC KERNELS
      zztmp = - 1._wp / rdtbt
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zu_frc(ji, jj) = zu_frc(ji, jj) + MAX(r1_hu_n(ji, jj) * r1_2 * (rCdU_bot(ji + 1, jj) + rCdU_bot(ji, jj)), zztmp) * zwx(ji, jj) * wdrampu(ji, jj)
          zv_frc(ji, jj) = zv_frc(ji, jj) + MAX(r1_hv_n(ji, jj) * r1_2 * (rCdU_bot(ji, jj + 1) + rCdU_bot(ji, jj)), zztmp) * zwy(ji, jj) * wdrampv(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zu_frc(ji, jj) = zu_frc(ji, jj) + r1_hu_n(ji, jj) * r1_2 * (rCdU_bot(ji + 1, jj) + rCdU_bot(ji, jj)) * zwx(ji, jj)
          zv_frc(ji, jj) = zv_frc(ji, jj) + r1_hv_n(ji, jj) * r1_2 * (rCdU_bot(ji, jj + 1) + rCdU_bot(ji, jj)) * zwy(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    IF (ln_isfcav) THEN
      IF (ln_bt_fw) THEN
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            iktu = miku(ji, jj)
            iktv = mikv(ji, jj)
            zwx(ji, jj) = un(ji, jj, iktu) - un_b(ji, jj)
            zwy(ji, jj) = vn(ji, jj, iktv) - vn_b(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            iktu = miku(ji, jj)
            iktv = mikv(ji, jj)
            zwx(ji, jj) = ub(ji, jj, iktu) - ub_b(ji, jj)
            zwy(ji, jj) = vb(ji, jj, iktv) - vb_b(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zu_frc(ji, jj) = zu_frc(ji, jj) + r1_hu_n(ji, jj) * r1_2 * (rCdU_top(ji + 1, jj) + rCdU_top(ji, jj)) * zwx(ji, jj)
          zv_frc(ji, jj) = zv_frc(ji, jj) + r1_hv_n(ji, jj) * r1_2 * (rCdU_top(ji, jj + 1) + rCdU_top(ji, jj)) * zwy(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    IF (ln_bt_fw) THEN
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zu_frc(ji, jj) = zu_frc(ji, jj) + r1_rau0 * utau(ji, jj) * r1_hu_n(ji, jj)
          zv_frc(ji, jj) = zv_frc(ji, jj) + r1_rau0 * vtau(ji, jj) * r1_hv_n(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      zztmp = r1_rau0 * r1_2
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zu_frc(ji, jj) = zu_frc(ji, jj) + zztmp * (utau_b(ji, jj) + utau(ji, jj)) * r1_hu_n(ji, jj)
          zv_frc(ji, jj) = zv_frc(ji, jj) + zztmp * (vtau_b(ji, jj) + vtau(ji, jj)) * r1_hv_n(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    IF (ln_apr_dyn) THEN
      IF (ln_bt_fw) THEN
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zu_spg = grav * (ssh_ib(ji + 1, jj) - ssh_ib(ji, jj)) * r1_e1u(ji, jj)
            zv_spg = grav * (ssh_ib(ji, jj + 1) - ssh_ib(ji, jj)) * r1_e2v(ji, jj)
            zu_frc(ji, jj) = zu_frc(ji, jj) + zu_spg
            zv_frc(ji, jj) = zv_frc(ji, jj) + zv_spg
          END DO
        END DO
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        zztmp = grav * r1_2
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zu_spg = zztmp * (ssh_ib(ji + 1, jj) - ssh_ib(ji, jj) + ssh_ibb(ji + 1, jj) - ssh_ibb(ji, jj)) * r1_e1u(ji, jj)
            zv_spg = zztmp * (ssh_ib(ji, jj + 1) - ssh_ib(ji, jj) + ssh_ibb(ji, jj + 1) - ssh_ibb(ji, jj)) * r1_e2v(ji, jj)
            zu_frc(ji, jj) = zu_frc(ji, jj) + zu_spg
            zv_frc(ji, jj) = zv_frc(ji, jj) + zv_spg
          END DO
        END DO
        !$ACC END KERNELS
      END IF
    END IF
    IF (ln_bt_fw) THEN
      !$ACC KERNELS
      zssh_frc(:, :) = r1_rau0 * (emp(:, :) - rnf(:, :) + fwfisf(:, :))
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      zztmp = r1_rau0 * r1_2
      zssh_frc(:, :) = zztmp * (emp(:, :) + emp_b(:, :) - rnf(:, :) - rnf_b(:, :) + fwfisf(:, :) + fwfisf_b(:, :))
      !$ACC END KERNELS
    END IF
    IF (ln_sdw) THEN
      !$ACC KERNELS
      zssh_frc(:, :) = zssh_frc(:, :) + div_sd(:, :)
      !$ACC END KERNELS
    END IF
    IF (ll_init) THEN
      !$ACC KERNELS
      sshbb_e(:, :) = 0._wp
      ubb_e(:, :) = 0._wp
      vbb_e(:, :) = 0._wp
      sshb_e(:, :) = 0._wp
      ub_e(:, :) = 0._wp
      vb_e(:, :) = 0._wp
      !$ACC END KERNELS
    END IF
    IF (ln_bt_fw) THEN
      !$ACC KERNELS
      sshn_e(:, :) = sshn(:, :)
      un_e(:, :) = un_b(:, :)
      vn_e(:, :) = vn_b(:, :)
      hu_e(:, :) = hu_n(:, :)
      hv_e(:, :) = hv_n(:, :)
      hur_e(:, :) = r1_hu_n(:, :)
      hvr_e(:, :) = r1_hv_n(:, :)
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      sshn_e(:, :) = sshb(:, :)
      un_e(:, :) = ub_b(:, :)
      vn_e(:, :) = vb_b(:, :)
      hu_e(:, :) = hu_b(:, :)
      hv_e(:, :) = hv_b(:, :)
      hur_e(:, :) = r1_hu_b(:, :)
      hvr_e(:, :) = r1_hv_b(:, :)
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    ua_b(:, :) = 0._wp
    va_b(:, :) = 0._wp
    ssha(:, :) = 0._wp
    un_adv(:, :) = 0._wp
    vn_adv(:, :) = 0._wp
    !$ACC END KERNELS
    IF (ln_wd_dl) THEN
      !$ACC KERNELS
      zuwdmask(:, :) = 0._wp
      zvwdmask(:, :) = 0._wp
      zuwdav2(:, :) = 0._wp
      zvwdav2(:, :) = 0._wp
      !$ACC END KERNELS
    END IF
    DO jn = 1, icycle
      CALL profile_psy_data5 % PreStart('dyn_spg_ts', 'r5', 0, 0)
      l_full_nf_update = jn == icycle
      IF (ln_bdy .AND. ln_tide) CALL bdy_dta_tides(kt, kit = jn, time_offset = noffset + 1)
      IF (ln_tide_pot .AND. ln_tide) CALL upd_tide(kt, kit = jn, time_offset = noffset)
      CALL profile_psy_data5 % PostEnd
      !$ACC KERNELS
      IF ((jn < 3) .AND. ll_init) THEN
        za1 = 1._wp
        za2 = 0._wp
        za3 = 0._wp
      ELSE
        za1 = 1.781105_wp
        za2 = - 1.06221_wp
        za3 = 0.281105_wp
      END IF
      ua_e(:, :) = za1 * un_e(:, :) + za2 * ub_e(:, :) + za3 * ubb_e(:, :)
      va_e(:, :) = za1 * vn_e(:, :) + za2 * vb_e(:, :) + za3 * vbb_e(:, :)
      !$ACC END KERNELS
      IF (.NOT. ln_linssh) THEN
        !$ACC KERNELS
        zsshp2_e(:, :) = za1 * sshn_e(:, :) + za2 * sshb_e(:, :) + za3 * sshbb_e(:, :)
        !$ACC END KERNELS
        IF (ln_wd_dl) THEN
          IF (ln_wd_dl_rmp) THEN
            CALL profile_psy_data6 % PreStart('dyn_spg_ts', 'r6', 0, 0)
            DO jj = 1, jpj
              DO ji = 1, jpi
                IF (zsshp2_e(ji, jj) + ht_0(ji, jj) > 2._wp * rn_wdmin1) THEN
                  ztwdmask(ji, jj) = 1._wp
                ELSE IF (zsshp2_e(ji, jj) + ht_0(ji, jj) > rn_wdmin1) THEN
                  ztwdmask(ji, jj) = (TANH(50._wp * ((zsshp2_e(ji, jj) + ht_0(ji, jj) - rn_wdmin1) * r_rn_wdmin1)))
                ELSE
                  ztwdmask(ji, jj) = 0._wp
                END IF
              END DO
            END DO
            CALL profile_psy_data6 % PostEnd
          ELSE
            !$ACC KERNELS
            !$ACC loop independent collapse(2)
            DO jj = 1, jpj
              DO ji = 1, jpi
                IF (zsshp2_e(ji, jj) + ht_0(ji, jj) > rn_wdmin1) THEN
                  ztwdmask(ji, jj) = 1._wp
                ELSE
                  ztwdmask(ji, jj) = 0._wp
                END IF
              END DO
            END DO
            !$ACC END KERNELS
          END IF
        END IF
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zwx(ji, jj) = r1_2 * ssumask(ji, jj) * r1_e1e2u(ji, jj) * (e1e2t(ji, jj) * zsshp2_e(ji, jj) + e1e2t(ji + 1, jj) * zsshp2_e(ji + 1, jj))
            zwy(ji, jj) = r1_2 * ssvmask(ji, jj) * r1_e1e2v(ji, jj) * (e1e2t(ji, jj) * zsshp2_e(ji, jj) + e1e2t(ji, jj + 1) * zsshp2_e(ji, jj + 1))
          END DO
        END DO
        !$ACC END KERNELS
        CALL profile_psy_data7 % PreStart('dyn_spg_ts', 'r7', 0, 0)
        CALL lbc_lnk_multi('dynspg_ts', zwx, 'U', 1._wp, zwy, 'V', 1._wp)
        CALL profile_psy_data7 % PostEnd
        !$ACC KERNELS
        zhup2_e(:, :) = hu_0(:, :) + zwx(:, :)
        zhvp2_e(:, :) = hv_0(:, :) + zwy(:, :)
        zhtp2_e(:, :) = ht_0(:, :) + zsshp2_e(:, :)
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        zhup2_e(:, :) = hu_n(:, :)
        zhvp2_e(:, :) = hv_n(:, :)
        zhtp2_e(:, :) = ht_n(:, :)
        !$ACC END KERNELS
      END IF
      CALL profile_psy_data8 % PreStart('dyn_spg_ts', 'r8', 0, 0)
      IF (ln_bdy .AND. ln_vol) CALL bdy_vol2d(kt, jn, ua_e, va_e, zhup2_e, zhvp2_e)
      CALL profile_psy_data8 % PostEnd
      !$ACC KERNELS
      zwx(:, :) = e2u(:, :) * ua_e(:, :) * zhup2_e(:, :)
      zwy(:, :) = e1v(:, :) * va_e(:, :) * zhvp2_e(:, :)
      !$ACC END KERNELS
      CALL profile_psy_data9 % PreStart('dyn_spg_ts', 'r9', 0, 0)
      IF (ln_wd_il) CALL wad_lmt_bt(zwx, zwy, sshn_e, zssh_frc, rdtbt)
      CALL profile_psy_data9 % PostEnd
      IF (ln_wd_dl) THEN
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            IF (zwx(ji, jj) > 0.0) THEN
              zuwdmask(ji, jj) = ztwdmask(ji, jj)
            ELSE
              zuwdmask(ji, jj) = ztwdmask(ji + 1, jj)
            END IF
            zwx(ji, jj) = zuwdmask(ji, jj) * zwx(ji, jj)
            un_e(ji, jj) = zuwdmask(ji, jj) * un_e(ji, jj)
            IF (zwy(ji, jj) > 0.0) THEN
              zvwdmask(ji, jj) = ztwdmask(ji, jj)
            ELSE
              zvwdmask(ji, jj) = ztwdmask(ji, jj + 1)
            END IF
            zwy(ji, jj) = zvwdmask(ji, jj) * zwy(ji, jj)
            vn_e(ji, jj) = zvwdmask(ji, jj) * vn_e(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      za2 = wgtbtp2(jn)
      un_adv(:, :) = un_adv(:, :) + za2 * zwx(:, :) * r1_e2u(:, :)
      vn_adv(:, :) = vn_adv(:, :) + za2 * zwy(:, :) * r1_e1v(:, :)
      !$ACC END KERNELS
      IF (ln_wd_dl_bc) THEN
        !$ACC KERNELS
        zuwdav2(:, :) = zuwdav2(:, :) + za2 * zuwdmask(:, :)
        zvwdav2(:, :) = zvwdav2(:, :) + za2 * zvwdmask(:, :)
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zhdiv(ji, jj) = (zwx(ji, jj) - zwx(ji - 1, jj) + zwy(ji, jj) - zwy(ji, jj - 1)) * r1_e1e2t(ji, jj)
        END DO
      END DO
      ssha_e(:, :) = (sshn_e(:, :) - rdtbt * (zssh_frc(:, :) + zhdiv(:, :))) * ssmask(:, :)
      !$ACC END KERNELS
      CALL profile_psy_data10 % PreStart('dyn_spg_ts', 'r10', 0, 0)
      CALL lbc_lnk('dynspg_ts', ssha_e, 'T', 1._wp)
      IF (ln_bdy) CALL bdy_ssh(ssha_e)
      CALL profile_psy_data10 % PostEnd
      IF (.NOT. ln_linssh) THEN
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zsshu_a(ji, jj) = r1_2 * ssumask(ji, jj) * r1_e1e2u(ji, jj) * (e1e2t(ji, jj) * ssha_e(ji, jj) + e1e2t(ji + 1, jj) * ssha_e(ji + 1, jj))
            zsshv_a(ji, jj) = r1_2 * ssvmask(ji, jj) * r1_e1e2v(ji, jj) * (e1e2t(ji, jj) * ssha_e(ji, jj) + e1e2t(ji, jj + 1) * ssha_e(ji, jj + 1))
          END DO
        END DO
        !$ACC END KERNELS
        CALL profile_psy_data11 % PreStart('dyn_spg_ts', 'r11', 0, 0)
        CALL lbc_lnk_multi('dynspg_ts', zsshu_a, 'U', 1._wp, zsshv_a, 'V', 1._wp)
        CALL profile_psy_data11 % PostEnd
      END IF
      !$ACC KERNELS
      IF ((jn == 1) .AND. ll_init) THEN
        za0 = 1._wp
        za1 = 0._wp
        za2 = 0._wp
        za3 = 0._wp
      ELSE IF ((jn == 2) .AND. ll_init) THEN
        za0 = 1.0833333333333_wp
        za1 = - 0.1666666666666_wp
        za2 = 0.0833333333333_wp
        za3 = 0._wp
      ELSE
        IF (rn_bt_alpha == 0._wp) THEN
          za0 = 0.614_wp
          za1 = 0.285_wp
          za2 = 0.088_wp
          za3 = 0.013_wp
        ELSE
          zepsilon = 0.00976186_wp - 0.13451357_wp * rn_bt_alpha
          zgamma = 0.08344500_wp - 0.51358400_wp * rn_bt_alpha
          za0 = 0.5_wp + zgamma + 2._wp * rn_bt_alpha + 2._wp * zepsilon
          za1 = 1._wp - za0 - zgamma - zepsilon
          za2 = zgamma
          za3 = zepsilon
        END IF
      END IF
      zsshp2_e(:, :) = za0 * ssha_e(:, :) + za1 * sshn_e(:, :) + za2 * sshb_e(:, :) + za3 * sshbb_e(:, :)
      !$ACC END KERNELS
      CALL profile_psy_data12 % PreStart('dyn_spg_ts', 'r12', 0, 0)
      IF (ln_wd_il) THEN
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            ll_tmp1 = MIN(zsshp2_e(ji, jj), zsshp2_e(ji + 1, jj)) > MAX(- ht_0(ji, jj), - ht_0(ji + 1, jj)) .AND. MAX(zsshp2_e(ji, jj) + ht_0(ji, jj), zsshp2_e(ji + 1, jj) + ht_0(ji + 1, jj)) > rn_wdmin1 + rn_wdmin2
            ll_tmp2 = (ABS(zsshp2_e(ji, jj) - zsshp2_e(ji + 1, jj)) > 1.E-12) .AND. (MAX(zsshp2_e(ji, jj), zsshp2_e(ji + 1, jj)) > MAX(- ht_0(ji, jj), - ht_0(ji + 1, jj)) + rn_wdmin1 + rn_wdmin2)
            IF (ll_tmp1) THEN
              zcpx(ji, jj) = 1.0_wp
            ELSE IF (ll_tmp2) THEN
              zcpx(ji, jj) = ABS((zsshp2_e(ji + 1, jj) + ht_0(ji + 1, jj) - zsshp2_e(ji, jj) - ht_0(ji, jj)) / (zsshp2_e(ji + 1, jj) - zsshp2_e(ji, jj)))
            ELSE
              zcpx(ji, jj) = 0._wp
            END IF
            ll_tmp1 = MIN(zsshp2_e(ji, jj), zsshp2_e(ji, jj + 1)) > MAX(- ht_0(ji, jj), - ht_0(ji, jj + 1)) .AND. MAX(zsshp2_e(ji, jj) + ht_0(ji, jj), zsshp2_e(ji, jj + 1) + ht_0(ji, jj + 1)) > rn_wdmin1 + rn_wdmin2
            ll_tmp2 = (ABS(zsshp2_e(ji, jj) - zsshp2_e(ji, jj + 1)) > 1.E-12) .AND. (MAX(zsshp2_e(ji, jj), zsshp2_e(ji, jj + 1)) > MAX(- ht_0(ji, jj), - ht_0(ji, jj + 1)) + rn_wdmin1 + rn_wdmin2)
            IF (ll_tmp1) THEN
              zcpy(ji, jj) = 1.0_wp
            ELSE IF (ll_tmp2) THEN
              zcpy(ji, jj) = ABS((zsshp2_e(ji, jj + 1) + ht_0(ji, jj + 1) - zsshp2_e(ji, jj) - ht_0(ji, jj)) / (zsshp2_e(ji, jj + 1) - zsshp2_e(ji, jj)))
            ELSE
              zcpy(ji, jj) = 0._wp
            END IF
          END DO
        END DO
      END IF
      CALL profile_psy_data12 % PostEnd
      !$ACC KERNELS
      IF (.NOT. ln_linssh .AND. .NOT. ln_dynadv_vec) THEN
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zx1 = r1_2 * ssumask(ji, jj) * r1_e1e2u(ji, jj) * (e1e2t(ji, jj) * zsshp2_e(ji, jj) + e1e2t(ji + 1, jj) * zsshp2_e(ji + 1, jj))
            zy1 = r1_2 * ssvmask(ji, jj) * r1_e1e2v(ji, jj) * (e1e2t(ji, jj) * zsshp2_e(ji, jj) + e1e2t(ji, jj + 1) * zsshp2_e(ji, jj + 1))
            zhust_e(ji, jj) = hu_0(ji, jj) + zx1
            zhvst_e(ji, jj) = hv_0(ji, jj) + zy1
          END DO
        END DO
      END IF
      SELECT CASE (nvor_scheme)
      CASE (np_ENT)
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            z1_hu = ssumask(ji, jj) / (zhup2_e(ji, jj) + 1._wp - ssumask(ji, jj))
            z1_hv = ssvmask(ji, jj) / (zhvp2_e(ji, jj) + 1._wp - ssvmask(ji, jj))
            zu_trd(ji, jj) = + r1_4 * r1_e1e2u(ji, jj) * z1_hu * (e1e2t(ji + 1, jj) * zhtp2_e(ji + 1, jj) * ff_t(ji + 1, jj) * (va_e(ji + 1, jj) + va_e(ji + 1, jj - 1)) + e1e2t(ji, jj) * zhtp2_e(ji, jj) * ff_t(ji, jj) * (va_e(ji, jj) + va_e(ji, jj - 1)))
            zv_trd(ji, jj) = - r1_4 * r1_e1e2v(ji, jj) * z1_hv * (e1e2t(ji, jj + 1) * zhtp2_e(ji, jj + 1) * ff_t(ji, jj + 1) * (ua_e(ji, jj + 1) + ua_e(ji - 1, jj + 1)) + e1e2t(ji, jj) * zhtp2_e(ji, jj) * ff_t(ji, jj) * (ua_e(ji, jj) + ua_e(ji - 1, jj)))
          END DO
        END DO
      CASE (np_ENE, np_MIX)
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zy1 = (zwy(ji, jj - 1) + zwy(ji + 1, jj - 1)) * r1_e1u(ji, jj)
            zy2 = (zwy(ji, jj) + zwy(ji + 1, jj)) * r1_e1u(ji, jj)
            zx1 = (zwx(ji - 1, jj) + zwx(ji - 1, jj + 1)) * r1_e2v(ji, jj)
            zx2 = (zwx(ji, jj) + zwx(ji, jj + 1)) * r1_e2v(ji, jj)
            zu_trd(ji, jj) = r1_4 * (zwz(ji, jj - 1) * zy1 + zwz(ji, jj) * zy2)
            zv_trd(ji, jj) = - r1_4 * (zwz(ji - 1, jj) * zx1 + zwz(ji, jj) * zx2)
          END DO
        END DO
      CASE (np_ENS)
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zy1 = r1_8 * (zwy(ji, jj - 1) + zwy(ji + 1, jj - 1) + zwy(ji, jj) + zwy(ji + 1, jj)) * r1_e1u(ji, jj)
            zx1 = - r1_8 * (zwx(ji - 1, jj) + zwx(ji - 1, jj + 1) + zwx(ji, jj) + zwx(ji, jj + 1)) * r1_e2v(ji, jj)
            zu_trd(ji, jj) = zy1 * (zwz(ji, jj - 1) + zwz(ji, jj))
            zv_trd(ji, jj) = zx1 * (zwz(ji - 1, jj) + zwz(ji, jj))
          END DO
        END DO
      CASE (np_EET, np_EEN)
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zu_trd(ji, jj) = + r1_12 * r1_e1u(ji, jj) * (ftne(ji, jj) * zwy(ji, jj) + ftnw(ji + 1, jj) * zwy(ji + 1, jj) + ftse(ji, jj) * zwy(ji, jj - 1) + ftsw(ji + 1, jj) * zwy(ji + 1, jj - 1))
            zv_trd(ji, jj) = - r1_12 * r1_e2v(ji, jj) * (ftsw(ji, jj + 1) * zwx(ji - 1, jj + 1) + ftse(ji, jj + 1) * zwx(ji, jj + 1) + ftnw(ji, jj) * zwx(ji - 1, jj) + ftne(ji, jj) * zwx(ji, jj))
          END DO
        END DO
      END SELECT
      IF (ln_tide .AND. ln_tide_pot) THEN
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zu_spg = grav * (pot_astro(ji + 1, jj) - pot_astro(ji, jj)) * r1_e1u(ji, jj)
            zv_spg = grav * (pot_astro(ji, jj + 1) - pot_astro(ji, jj)) * r1_e2v(ji, jj)
            zu_trd(ji, jj) = zu_trd(ji, jj) + zu_spg
            zv_trd(ji, jj) = zv_trd(ji, jj) + zv_spg
          END DO
        END DO
      END IF
      !$ACC END KERNELS
      IF (.NOT. ll_wd) THEN
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zu_trd(ji, jj) = zu_trd(ji, jj) + zCdU_u(ji, jj) * un_e(ji, jj) * hur_e(ji, jj)
            zv_trd(ji, jj) = zv_trd(ji, jj) + zCdU_v(ji, jj) * vn_e(ji, jj) * hvr_e(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      END IF
      IF (ln_wd_il) THEN
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zu_spg = - grav * (zsshp2_e(ji + 1, jj) - zsshp2_e(ji, jj)) * r1_e1u(ji, jj)
            zv_spg = - grav * (zsshp2_e(ji, jj + 1) - zsshp2_e(ji, jj)) * r1_e2v(ji, jj)
            zwx(ji, jj) = (1._wp - rn_scal_load) * zu_spg * zcpx(ji, jj)
            zwy(ji, jj) = (1._wp - rn_scal_load) * zv_spg * zcpy(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zu_spg = - grav * (zsshp2_e(ji + 1, jj) - zsshp2_e(ji, jj)) * r1_e1u(ji, jj)
            zv_spg = - grav * (zsshp2_e(ji, jj + 1) - zsshp2_e(ji, jj)) * r1_e2v(ji, jj)
            zwx(ji, jj) = (1._wp - rn_scal_load) * zu_spg
            zwy(ji, jj) = (1._wp - rn_scal_load) * zv_spg
          END DO
        END DO
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      IF (ln_dynadv_vec .OR. ln_linssh) THEN
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            ua_e(ji, jj) = (un_e(ji, jj) + rdtbt * (zwx(ji, jj) + zu_trd(ji, jj) + zu_frc(ji, jj))) * ssumask(ji, jj)
            va_e(ji, jj) = (vn_e(ji, jj) + rdtbt * (zwy(ji, jj) + zv_trd(ji, jj) + zv_frc(ji, jj))) * ssvmask(ji, jj)
          END DO
        END DO
      ELSE
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zhura = hu_0(ji, jj) + zsshu_a(ji, jj)
            zhvra = hv_0(ji, jj) + zsshv_a(ji, jj)
            zhura = ssumask(ji, jj) / (zhura + 1._wp - ssumask(ji, jj))
            zhvra = ssvmask(ji, jj) / (zhvra + 1._wp - ssvmask(ji, jj))
            ua_e(ji, jj) = (hu_e(ji, jj) * un_e(ji, jj) + rdtbt * (zhust_e(ji, jj) * zwx(ji, jj) + zhup2_e(ji, jj) * zu_trd(ji, jj) + hu_n(ji, jj) * zu_frc(ji, jj))) * zhura
            va_e(ji, jj) = (hv_e(ji, jj) * vn_e(ji, jj) + rdtbt * (zhvst_e(ji, jj) * zwy(ji, jj) + zhvp2_e(ji, jj) * zv_trd(ji, jj) + hv_n(ji, jj) * zv_frc(ji, jj))) * zhvra
          END DO
        END DO
      END IF
      !$ACC END KERNELS
      IF (ll_wd) THEN
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            ua_e(ji, jj) = ua_e(ji, jj) / (1.0 - rdtbt * zCdU_u(ji, jj) * hur_e(ji, jj))
            va_e(ji, jj) = va_e(ji, jj) / (1.0 - rdtbt * zCdU_v(ji, jj) * hvr_e(ji, jj))
          END DO
        END DO
        !$ACC END KERNELS
      END IF
      IF (.NOT. ln_linssh) THEN
        !$ACC KERNELS
        hu_e(:, :) = hu_0(:, :) + zsshu_a(:, :)
        hv_e(:, :) = hv_0(:, :) + zsshv_a(:, :)
        hur_e(:, :) = ssumask(:, :) / (hu_e(:, :) + 1._wp - ssumask(:, :))
        hvr_e(:, :) = ssvmask(:, :) / (hv_e(:, :) + 1._wp - ssvmask(:, :))
        !$ACC END KERNELS
      END IF
      CALL profile_psy_data13 % PreStart('dyn_spg_ts', 'r13', 0, 0)
      CALL lbc_lnk_multi('dynspg_ts', ua_e, 'U', - 1._wp, va_e, 'V', - 1._wp)
      IF (ln_bdy) CALL bdy_dyn2d(jn, ua_e, va_e, un_e, vn_e, hur_e, hvr_e, ssha_e)
      CALL profile_psy_data13 % PostEnd
      !$ACC KERNELS
      ubb_e(:, :) = ub_e(:, :)
      ub_e(:, :) = un_e(:, :)
      un_e(:, :) = ua_e(:, :)
      vbb_e(:, :) = vb_e(:, :)
      vb_e(:, :) = vn_e(:, :)
      vn_e(:, :) = va_e(:, :)
      sshbb_e(:, :) = sshb_e(:, :)
      sshb_e(:, :) = sshn_e(:, :)
      sshn_e(:, :) = ssha_e(:, :)
      za1 = wgtbtp1(jn)
      !$ACC END KERNELS
      IF (ln_dynadv_vec .OR. ln_linssh) THEN
        !$ACC KERNELS
        ua_b(:, :) = ua_b(:, :) + za1 * ua_e(:, :)
        va_b(:, :) = va_b(:, :) + za1 * va_e(:, :)
        !$ACC END KERNELS
      ELSE
        IF (.NOT. ln_wd_dl) THEN
          !$ACC KERNELS
          ua_b(:, :) = ua_b(:, :) + za1 * ua_e(:, :) * hu_e(:, :)
          va_b(:, :) = va_b(:, :) + za1 * va_e(:, :) * hv_e(:, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          ua_b(:, :) = ua_b(:, :) + za1 * ua_e(:, :) * hu_e(:, :) * zuwdmask(:, :)
          va_b(:, :) = va_b(:, :) + za1 * va_e(:, :) * hv_e(:, :) * zvwdmask(:, :)
          !$ACC END KERNELS
        END IF
      END IF
      !$ACC KERNELS
      ssha(:, :) = ssha(:, :) + za1 * ssha_e(:, :)
      !$ACC END KERNELS
    END DO
    IF (ln_bt_fw) THEN
      !$ACC KERNELS
      zwx(:, :) = un_adv(:, :)
      zwy(:, :) = vn_adv(:, :)
      !$ACC END KERNELS
      IF (.NOT. (kt == nit000 .AND. neuler == 0)) THEN
        !$ACC KERNELS
        un_adv(:, :) = r1_2 * (ub2_b(:, :) + zwx(:, :) - atfp * un_bf(:, :))
        vn_adv(:, :) = r1_2 * (vb2_b(:, :) + zwy(:, :) - atfp * vn_bf(:, :))
        un_bf(:, :) = atfp * un_bf(:, :) + (zwx(:, :) - ub2_b(:, :))
        vn_bf(:, :) = atfp * vn_bf(:, :) + (zwy(:, :) - vb2_b(:, :))
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        un_bf(:, :) = 0._wp
        vn_bf(:, :) = 0._wp
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      ub2_b(:, :) = zwx(:, :)
      vb2_b(:, :) = zwy(:, :)
      !$ACC END KERNELS
    END IF
    IF (ln_dynadv_vec .OR. ln_linssh) THEN
      !$ACC KERNELS
      DO jk = 1, jpkm1
        ua(:, :, jk) = ua(:, :, jk) + (ua_b(:, :) - ub_b(:, :)) * r1_2dt_b
        va(:, :, jk) = va(:, :, jk) + (va_b(:, :) - vb_b(:, :)) * r1_2dt_b
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zsshu_a(ji, jj) = r1_2 * ssumask(ji, jj) * r1_e1e2u(ji, jj) * (e1e2t(ji, jj) * ssha(ji, jj) + e1e2t(ji + 1, jj) * ssha(ji + 1, jj))
          zsshv_a(ji, jj) = r1_2 * ssvmask(ji, jj) * r1_e1e2v(ji, jj) * (e1e2t(ji, jj) * ssha(ji, jj) + e1e2t(ji, jj + 1) * ssha(ji, jj + 1))
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data14 % PreStart('dyn_spg_ts', 'r14', 0, 0)
      CALL lbc_lnk_multi('dynspg_ts', zsshu_a, 'U', 1._wp, zsshv_a, 'V', 1._wp)
      CALL profile_psy_data14 % PostEnd
      !$ACC KERNELS
      DO jk = 1, jpkm1
        ua(:, :, jk) = ua(:, :, jk) + r1_hu_n(:, :) * (ua_b(:, :) - ub_b(:, :) * hu_b(:, :)) * r1_2dt_b
        va(:, :, jk) = va(:, :, jk) + r1_hv_n(:, :) * (va_b(:, :) - vb_b(:, :) * hv_b(:, :)) * r1_2dt_b
      END DO
      ua_b(:, :) = ua_b(:, :) / (hu_0(:, :) + zsshu_a(:, :) + 1._wp - ssumask(:, :))
      va_b(:, :) = va_b(:, :) / (hv_0(:, :) + zsshv_a(:, :) + 1._wp - ssvmask(:, :))
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    DO jk = 1, jpkm1
      un(:, :, jk) = (un(:, :, jk) + un_adv(:, :) * r1_hu_n(:, :) - un_b(:, :)) * umask(:, :, jk)
      vn(:, :, jk) = (vn(:, :, jk) + vn_adv(:, :) * r1_hv_n(:, :) - vn_b(:, :)) * vmask(:, :, jk)
    END DO
    IF (ln_wd_dl .AND. ln_wd_dl_bc) THEN
      DO jk = 1, jpkm1
        un(:, :, jk) = (un_adv(:, :) * r1_hu_n(:, :) + zuwdav2(:, :) * (un(:, :, jk) - un_adv(:, :) * r1_hu_n(:, :))) * umask(:, :, jk)
        vn(:, :, jk) = (vn_adv(:, :) * r1_hv_n(:, :) + zvwdav2(:, :) * (vn(:, :, jk) - vn_adv(:, :) * r1_hv_n(:, :))) * vmask(:, :, jk)
      END DO
    END IF
    !$ACC END KERNELS
    CALL profile_psy_data15 % PreStart('dyn_spg_ts', 'r15', 0, 0)
    CALL iom_put("ubar", un_adv(:, :) * r1_hu_n(:, :))
    CALL iom_put("vbar", vn_adv(:, :) * r1_hv_n(:, :))
    IF (lrst_oce .AND. ln_bt_fw) CALL ts_rst(kt, 'WRITE')
    IF (ln_wd_il) DEALLOCATE(zcpx, zcpy)
    IF (ln_wd_dl) DEALLOCATE(ztwdmask, zuwdmask, zvwdmask, zuwdav2, zvwdav2)
    IF (ln_diatmb) THEN
      CALL iom_put("baro_u", un_b * ssumask(:, :) + zmdi * (1. - ssumask(:, :)))
      CALL iom_put("baro_v", vn_b * ssvmask(:, :) + zmdi * (1. - ssvmask(:, :)))
    END IF
    CALL profile_psy_data15 % PostEnd
  END SUBROUTINE dyn_spg_ts
  SUBROUTINE ts_wgt(ll_av, ll_fw, jpit, zwgt1, zwgt2)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    LOGICAL, INTENT(IN) :: ll_av
    LOGICAL, INTENT(IN) :: ll_fw
    INTEGER, INTENT(INOUT) :: jpit
    REAL(KIND = wp), DIMENSION(3 * nn_baro), INTENT(INOUT) :: zwgt1, zwgt2
    INTEGER :: jic, jn, ji
    REAL(KIND = wp) :: za1, za2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    zwgt1(:) = 0._wp
    zwgt2(:) = 0._wp
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('ts_wgt', 'r0', 0, 0)
    IF (ll_fw) THEN
      jic = nn_baro
    ELSE
      jic = 2 * nn_baro
    END IF
    IF (ll_av) THEN
      SELECT CASE (nn_bt_flt)
      CASE (0)
        zwgt1(jic) = 1._wp
        jpit = jic
      CASE (1)
        DO jn = 1, 3 * nn_baro
          za1 = ABS(FLOAT(jn - jic)) / FLOAT(nn_baro)
          IF (za1 < 0.5_wp) THEN
            zwgt1(jn) = 1._wp
            jpit = jn
          END IF
        END DO
      CASE (2)
        DO jn = 1, 3 * nn_baro
          za1 = ABS(FLOAT(jn - jic)) / FLOAT(nn_baro)
          IF (za1 < 1._wp) THEN
            zwgt1(jn) = 1._wp
            jpit = jn
          END IF
        END DO
      CASE DEFAULT
        CALL ctl_stop('unrecognised value for nn_bt_flt')
      END SELECT
    ELSE
      zwgt1(jic) = 1._wp
      jpit = jic
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    DO jn = 1, jpit
      DO ji = jn, jpit
        zwgt2(jn) = zwgt2(jn) + zwgt1(ji)
      END DO
    END DO
    za1 = 1._wp / SUM(zwgt1(1 : jpit))
    za2 = 1._wp / SUM(zwgt2(1 : jpit))
    DO jn = 1, jpit
      zwgt1(jn) = zwgt1(jn) * za1
      zwgt2(jn) = zwgt2(jn) * za2
    END DO
    !$ACC END KERNELS
  END SUBROUTINE ts_wgt
  SUBROUTINE ts_rst(kt, cdrw)
    INTEGER, INTENT(IN) :: kt
    CHARACTER(LEN = *), INTENT(IN) :: cdrw
    IF (TRIM(cdrw) == 'READ') THEN
      IF (ln_rstart .AND. ln_bt_fw .AND. (neuler /= 0)) THEN
        CALL iom_get(numror, jpdom_autoglo, 'ub2_b', ub2_b(:, :), ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'vb2_b', vb2_b(:, :), ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'un_bf', un_bf(:, :), ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'vn_bf', vn_bf(:, :), ldxios = lrxios)
        IF (.NOT. ln_bt_av) THEN
          CALL iom_get(numror, jpdom_autoglo, 'sshbb_e', sshbb_e(:, :), ldxios = lrxios)
          CALL iom_get(numror, jpdom_autoglo, 'ubb_e', ubb_e(:, :), ldxios = lrxios)
          CALL iom_get(numror, jpdom_autoglo, 'vbb_e', vbb_e(:, :), ldxios = lrxios)
          CALL iom_get(numror, jpdom_autoglo, 'sshb_e', sshb_e(:, :), ldxios = lrxios)
          CALL iom_get(numror, jpdom_autoglo, 'ub_e', ub_e(:, :), ldxios = lrxios)
          CALL iom_get(numror, jpdom_autoglo, 'vb_e', vb_e(:, :), ldxios = lrxios)
        END IF
      ELSE
        IF (lwp) WRITE(numout, *)
        IF (lwp) WRITE(numout, *) '   ==>>>   start from rest: set barotropic values to 0'
        !$ACC KERNELS
        ub2_b(:, :) = 0._wp
        vb2_b(:, :) = 0._wp
        un_adv(:, :) = 0._wp
        vn_adv(:, :) = 0._wp
        un_bf(:, :) = 0._wp
        vn_bf(:, :) = 0._wp
        !$ACC END KERNELS
      END IF
    ELSE IF (TRIM(cdrw) == 'WRITE') THEN
      IF (lwp) WRITE(numout, *) '---- ts_rst ----'
      IF (lwxios) CALL iom_swap(cwxios_context)
      CALL iom_rstput(kt, nitrst, numrow, 'ub2_b', ub2_b(:, :), ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'vb2_b', vb2_b(:, :), ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'un_bf', un_bf(:, :), ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'vn_bf', vn_bf(:, :), ldxios = lwxios)
      IF (.NOT. ln_bt_av) THEN
        CALL iom_rstput(kt, nitrst, numrow, 'sshbb_e', sshbb_e(:, :), ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'ubb_e', ubb_e(:, :), ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'vbb_e', vbb_e(:, :), ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'sshb_e', sshb_e(:, :), ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'ub_e', ub_e(:, :), ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'vb_e', vb_e(:, :), ldxios = lwxios)
      END IF
      IF (lwxios) CALL iom_swap(cxios_context)
    END IF
  END SUBROUTINE ts_rst
  SUBROUTINE dyn_spg_ts_init
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zxr2, zyr2, zcmax
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zcu
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zxr2 = r1_e1t(ji, jj) * r1_e1t(ji, jj)
        zyr2 = r1_e2t(ji, jj) * r1_e2t(ji, jj)
        zcu(ji, jj) = SQRT(grav * MAX(ht_0(ji, jj), 0._wp) * (zxr2 + zyr2))
      END DO
    END DO
    zcmax = MAXVAL(zcu(:, :))
    !$ACC END KERNELS
    CALL mpp_max('dynspg_ts', zcmax)
    IF (ln_bt_auto) nn_baro = CEILING(rdt / rn_bt_cmax * zcmax)
    rdtbt = rdt / REAL(nn_baro, wp)
    zcmax = zcmax * rdtbt
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) 'dyn_spg_ts_init : split-explicit free surface'
    IF (lwp) WRITE(numout, *) '~~~~~~~~~~~~~~~'
    IF (ln_bt_auto) THEN
      IF (lwp) WRITE(numout, *) '     ln_ts_auto =.true. Automatically set nn_baro '
      IF (lwp) WRITE(numout, *) '     Max. courant number allowed: ', rn_bt_cmax
    ELSE
      IF (lwp) WRITE(numout, *) '     ln_ts_auto=.false.: Use nn_baro in namelist   nn_baro = ', nn_baro
    END IF
    IF (ln_bt_av) THEN
      IF (lwp) WRITE(numout, *) '     ln_bt_av =.true.  ==> Time averaging over nn_baro time steps is on '
    ELSE
      IF (lwp) WRITE(numout, *) '     ln_bt_av =.false. => No time averaging of barotropic variables '
    END IF
    IF (ln_bt_fw) THEN
      IF (lwp) WRITE(numout, *) '     ln_bt_fw=.true.  => Forward integration of barotropic variables '
    ELSE
      IF (lwp) WRITE(numout, *) '     ln_bt_fw =.false.=> Centred integration of barotropic variables '
    END IF
    IF (lwp) WRITE(numout, *) '     Time filter choice, nn_bt_flt: ', nn_bt_flt
    SELECT CASE (nn_bt_flt)
    CASE (0)
      IF (lwp) WRITE(numout, *) '           Dirac'
    CASE (1)
      IF (lwp) WRITE(numout, *) '           Boxcar: width = nn_baro'
    CASE (2)
      IF (lwp) WRITE(numout, *) '           Boxcar: width = 2*nn_baro'
    CASE DEFAULT
      CALL ctl_stop('unrecognised value for nn_bt_flt: should 0,1, or 2')
    END SELECT
    IF (lwp) WRITE(numout, *) ' '
    IF (lwp) WRITE(numout, *) '     nn_baro = ', nn_baro
    IF (lwp) WRITE(numout, *) '     Barotropic time step [s] is :', rdtbt
    IF (lwp) WRITE(numout, *) '     Maximum Courant number is   :', zcmax
    IF (lwp) WRITE(numout, *) '     Time diffusion parameter rn_bt_alpha: ', rn_bt_alpha
    IF ((ln_bt_av .AND. nn_bt_flt /= 0) .AND. (rn_bt_alpha > 0._wp)) THEN
      CALL ctl_stop('dynspg_ts ERROR: if rn_bt_alpha > 0, remove temporal averaging')
    END IF
    IF (.NOT. ln_bt_av .AND. .NOT. ln_bt_fw) THEN
      CALL ctl_stop('dynspg_ts ERROR: No time averaging => only forward integration is possible')
    END IF
    IF (zcmax > 0.9_wp) THEN
      CALL ctl_stop('dynspg_ts ERROR: Maximum Courant number is greater than 0.9: Inc. nn_baro !')
    END IF
    IF (dyn_spg_ts_alloc() /= 0) CALL ctl_stop('STOP', 'dyn_spg_init: failed to allocate dynspg_ts  arrays')
    CALL ts_rst(nit000, 'READ')
    IF (lwxios) THEN
      CALL iom_set_rstw_var_active('ub2_b')
      CALL iom_set_rstw_var_active('vb2_b')
      CALL iom_set_rstw_var_active('un_bf')
      CALL iom_set_rstw_var_active('vn_bf')
      IF (.NOT. ln_bt_av) THEN
        CALL iom_set_rstw_var_active('sshbb_e')
        CALL iom_set_rstw_var_active('ubb_e')
        CALL iom_set_rstw_var_active('vbb_e')
        CALL iom_set_rstw_var_active('sshb_e')
        CALL iom_set_rstw_var_active('ub_e')
        CALL iom_set_rstw_var_active('vb_e')
      END IF
    END IF
  END SUBROUTINE dyn_spg_ts_init
END MODULE dynspg_ts