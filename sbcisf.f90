MODULE sbcisf
  USE oce
  USE dom_oce
  USE phycst
  USE eosbn2
  USE sbc_oce
  USE zdfdrg
  USE in_out_manager
  USE iom
  USE fldread
  USE lbclnk
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sbc_isf, sbc_isf_init, sbc_isf_div, sbc_isf_alloc
  REAL(KIND = wp), PUBLIC :: rn_hisf_tbl
  INTEGER, PUBLIC :: nn_isf
  INTEGER, PUBLIC :: nn_isfblk
  INTEGER, PUBLIC :: nn_gammablk
  REAL(KIND = wp), PUBLIC :: rn_gammat0
  REAL(KIND = wp), PUBLIC :: rn_gammas0
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: misfkt, misfkb
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: rzisf_tbl
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: rhisf_tbl, rhisf_tbl_0
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: r1_hisf_tbl
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ralpha
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: risfLeff
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ttbl, stbl, utbl, vtbl
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: qisf
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: risf_tsc_b, risf_tsc
  LOGICAL, PUBLIC :: l_isfcpl = .FALSE.
  REAL(KIND = wp), PUBLIC, SAVE :: rcpisf = 2000.0_wp
  REAL(KIND = wp), PUBLIC, SAVE :: rkappa = 1.54E-6_wp
  REAL(KIND = wp), PUBLIC, SAVE :: rhoisf = 920.0_wp
  REAL(KIND = wp), PUBLIC, SAVE :: tsurf = - 20.0_wp
  REAL(KIND = wp), PUBLIC, SAVE :: rLfusisf = 0.334E6_wp
  CHARACTER(LEN = 100), PUBLIC :: cn_dirisf = './'
  TYPE(FLD_N), PUBLIC :: sn_fwfisf
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_fwfisf
  TYPE(FLD_N), PUBLIC :: sn_rnfisf
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_rnfisf
  TYPE(FLD_N), PUBLIC :: sn_depmax_isf
  TYPE(FLD_N), PUBLIC :: sn_depmin_isf
  TYPE(FLD_N), PUBLIC :: sn_Leff_isf
  CONTAINS
  SUBROUTINE sbc_isf(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: ikt, ikb
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zt_frz, zdep
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: zqhcisf2d
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: zfwfisf3d, zqhcisf3d, zqlatisf3d
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
    IF (MOD(kt - 1, nn_fsbc) == 0) THEN
      SELECT CASE (nn_isf)
      CASE (1)
        CALL profile_psy_data0 % PreStart('sbc_isf', 'r0', 0, 0)
        CALL sbc_isf_tbl(tsn(:, :, :, jp_tem), ttbl(:, :), 'T')
        CALL sbc_isf_tbl(tsn(:, :, :, jp_sal), stbl(:, :), 'T')
        CALL sbc_isf_tbl(un(:, :, :), utbl(:, :), 'U')
        CALL sbc_isf_tbl(vn(:, :, :), vtbl(:, :), 'V')
        CALL iom_put('ttbl', ttbl(:, :))
        CALL iom_put('stbl', stbl(:, :))
        CALL iom_put('utbl', utbl(:, :) * (1._wp - tmask(:, :, 1)) * ssmask(:, :))
        CALL iom_put('vtbl', vtbl(:, :) * (1._wp - tmask(:, :, 1)) * ssmask(:, :))
        CALL profile_psy_data0 % PostEnd
        IF (.NOT. l_isfcpl) THEN
          CALL profile_psy_data1 % PreStart('sbc_isf', 'r1', 0, 0)
          CALL sbc_isf_cav(kt)
          CALL profile_psy_data1 % PostEnd
        ELSE
          !$ACC KERNELS
          qisf(:, :) = fwfisf(:, :) * rlfusisf
          !$ACC END KERNELS
        END IF
      CASE (2)
        !$ACC KERNELS
        stbl(:, :) = soce
        !$ACC END KERNELS
        CALL profile_psy_data2 % PreStart('sbc_isf', 'r2', 0, 0)
        CALL sbc_isf_bg03(kt)
        CALL profile_psy_data2 % PostEnd
      CASE (3)
        IF (.NOT. l_isfcpl) THEN
          CALL profile_psy_data3 % PreStart('sbc_isf', 'r3', 0, 0)
          CALL fld_read(kt, nn_fsbc, sf_rnfisf)
          CALL profile_psy_data3 % PostEnd
          !$ACC KERNELS
          fwfisf(:, :) = - sf_rnfisf(1) % fnow(:, :, 1)
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        qisf(:, :) = fwfisf(:, :) * rLfusisf
        stbl(:, :) = soce
        !$ACC END KERNELS
      CASE (4)
        IF (.NOT. l_isfcpl) THEN
          CALL profile_psy_data4 % PreStart('sbc_isf', 'r4', 0, 0)
          CALL fld_read(kt, nn_fsbc, sf_fwfisf)
          CALL profile_psy_data4 % PostEnd
          !$ACC KERNELS
          fwfisf(:, :) = - sf_fwfisf(1) % fnow(:, :, 1)
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        qisf(:, :) = fwfisf(:, :) * rLfusisf
        stbl(:, :) = soce
        !$ACC END KERNELS
      END SELECT
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          zdep(ji, jj) = gdepw_n(ji, jj, misfkt(ji, jj))
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data5 % PreStart('sbc_isf', 'r5', 0, 0)
      CALL eos_fzp(stbl(:, :), zt_frz(:, :), zdep(:, :))
      CALL profile_psy_data5 % PostEnd
      !$ACC KERNELS
      risf_tsc(:, :, jp_tem) = qisf(:, :) * r1_rau0_rcp - fwfisf(:, :) * zt_frz(:, :) * r1_rau0
      risf_tsc(:, :, jp_sal) = 0.0_wp
      !$ACC END KERNELS
      CALL profile_psy_data6 % PreStart('sbc_isf', 'r6', 0, 0)
      CALL lbc_lnk_multi('sbcisf', risf_tsc(:, :, jp_tem), 'T', 1., risf_tsc(:, :, jp_sal), 'T', 1., fwfisf, 'T', 1., qisf, 'T', 1.)
      IF (iom_use('iceshelf_cea')) CALL iom_put('iceshelf_cea', - fwfisf(:, :))
      IF (iom_use('hflx_isf_cea')) CALL iom_put('hflx_isf_cea', risf_tsc(:, :, jp_tem) * rau0 * rcp)
      IF (iom_use('qlatisf')) CALL iom_put('qlatisf', qisf(:, :))
      IF (iom_use('fwfisf')) CALL iom_put('fwfisf', fwfisf(:, :))
      CALL profile_psy_data6 % PostEnd
      IF (iom_use('fwfisf3d') .OR. iom_use('qlatisf3d') .OR. iom_use('qhcisf3d') .OR. iom_use('qhcisf')) THEN
        CALL profile_psy_data7 % PreStart('sbc_isf', 'r7', 0, 0)
        ALLOCATE(zfwfisf3d(jpi, jpj, jpk), zqhcisf3d(jpi, jpj, jpk), zqlatisf3d(jpi, jpj, jpk))
        ALLOCATE(zqhcisf2d(jpi, jpj))
        CALL profile_psy_data7 % PostEnd
        !$ACC KERNELS
        zfwfisf3d(:, :, :) = 0._wp
        zqhcisf3d(:, :, :) = 0._wp
        zqlatisf3d(:, :, :) = 0._wp
        zqhcisf2d(:, :) = fwfisf(:, :) * zt_frz * rcp
        !$ACC loop independent collapse(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            ikt = misfkt(ji, jj)
            ikb = misfkb(ji, jj)
            DO jk = ikt, ikb - 1
              zfwfisf3d(ji, jj, jk) = zfwfisf3d(ji, jj, jk) + fwfisf(ji, jj) * r1_hisf_tbl(ji, jj) * e3t_n(ji, jj, jk)
              zqhcisf3d(ji, jj, jk) = zqhcisf3d(ji, jj, jk) + zqhcisf2d(ji, jj) * r1_hisf_tbl(ji, jj) * e3t_n(ji, jj, jk)
              zqlatisf3d(ji, jj, jk) = zqlatisf3d(ji, jj, jk) + qisf(ji, jj) * r1_hisf_tbl(ji, jj) * e3t_n(ji, jj, jk)
            END DO
            zfwfisf3d(ji, jj, jk) = zfwfisf3d(ji, jj, jk) + fwfisf(ji, jj) * r1_hisf_tbl(ji, jj) * ralpha(ji, jj) * e3t_n(ji, jj, jk)
            zqhcisf3d(ji, jj, jk) = zqhcisf3d(ji, jj, jk) + zqhcisf2d(ji, jj) * r1_hisf_tbl(ji, jj) * ralpha(ji, jj) * e3t_n(ji, jj, jk)
            zqlatisf3d(ji, jj, jk) = zqlatisf3d(ji, jj, jk) + qisf(ji, jj) * r1_hisf_tbl(ji, jj) * ralpha(ji, jj) * e3t_n(ji, jj, jk)
          END DO
        END DO
        !$ACC END KERNELS
        CALL profile_psy_data8 % PreStart('sbc_isf', 'r8', 0, 0)
        CALL iom_put('fwfisf3d', zfwfisf3d(:, :, :))
        CALL iom_put('qlatisf3d', zqlatisf3d(:, :, :))
        CALL iom_put('qhcisf3d', zqhcisf3d(:, :, :))
        CALL iom_put('qhcisf', zqhcisf2d(:, :))
        DEALLOCATE(zfwfisf3d, zqhcisf3d, zqlatisf3d)
        DEALLOCATE(zqhcisf2d)
        CALL profile_psy_data8 % PostEnd
      END IF
    END IF
    IF (kt == nit000) THEN
      IF (ln_rstart .AND. iom_varid(numror, 'fwf_isf_b', ldstop = .FALSE.) > 0) THEN
        CALL profile_psy_data9 % PreStart('sbc_isf', 'r9', 0, 0)
        IF (lwp) WRITE(numout, *) '          nit000-1 isf tracer content forcing fields read in the restart file'
        CALL iom_get(numror, jpdom_autoglo, 'fwf_isf_b', fwfisf_b(:, :), ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'isf_sc_b', risf_tsc_b(:, :, jp_sal), ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'isf_hc_b', risf_tsc_b(:, :, jp_tem), ldxios = lrxios)
        CALL profile_psy_data9 % PostEnd
      ELSE
        !$ACC KERNELS
        fwfisf_b(:, :) = fwfisf(:, :)
        risf_tsc_b(:, :, :) = risf_tsc(:, :, :)
        !$ACC END KERNELS
      END IF
    END IF
    CALL profile_psy_data10 % PreStart('sbc_isf', 'r10', 0, 0)
    IF (lrst_oce) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'sbc : isf surface tracer content forcing fields written in ocean restart file ', 'at it= ', kt, ' date= ', ndastp
      IF (lwp) WRITE(numout, *) '~~~~'
      IF (lwxios) CALL iom_swap(cwxios_context)
      CALL iom_rstput(kt, nitrst, numrow, 'fwf_isf_b', fwfisf(:, :), ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'isf_hc_b', risf_tsc(:, :, jp_tem), ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'isf_sc_b', risf_tsc(:, :, jp_sal), ldxios = lwxios)
      IF (lwxios) CALL iom_swap(cxios_context)
    END IF
    CALL profile_psy_data10 % PostEnd
  END SUBROUTINE sbc_isf
  INTEGER FUNCTION sbc_isf_alloc()
    sbc_isf_alloc = 0
    IF (.NOT. ALLOCATED(qisf)) THEN
      ALLOCATE(risf_tsc(jpi, jpj, jpts), risf_tsc_b(jpi, jpj, jpts), qisf(jpi, jpj), rhisf_tbl(jpi, jpj), r1_hisf_tbl(jpi, jpj), rzisf_tbl(jpi, jpj), ttbl(jpi, jpj), stbl(jpi, jpj), utbl(jpi, jpj), vtbl(jpi, jpj), risfLeff(jpi, jpj), rhisf_tbl_0(jpi, jpj), ralpha(jpi, jpj), misfkt(jpi, jpj), misfkb(jpi, jpj), STAT = sbc_isf_alloc)
      CALL mpp_sum('sbcisf', sbc_isf_alloc)
      IF (sbc_isf_alloc /= 0) CALL ctl_stop('STOP', 'sbc_isf_alloc: failed to allocate arrays.')
    END IF
  END FUNCTION
  SUBROUTINE sbc_isf_init
    INTEGER :: ji, jj, jk
    INTEGER :: ik
    INTEGER :: ikt, ikb
    INTEGER :: inum, ierror
    INTEGER :: ios
    REAL(KIND = wp) :: zhk
    CHARACTER(LEN = 256) :: cvarzisf, cvarhisf
    CHARACTER(LEN = 32) :: cvarLeff
    NAMELIST /namsbc_isf/ nn_isfblk, rn_hisf_tbl, rn_gammat0, rn_gammas0, nn_gammablk, nn_isf, sn_fwfisf, sn_rnfisf, sn_depmax_isf, sn_depmin_isf, sn_Leff_isf
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namsbc_isf, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namsbc_isf in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namsbc_isf, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namsbc_isf in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namsbc_isf)
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) 'sbc_isf_init : heat flux of the ice shelf'
    IF (lwp) WRITE(numout, *) '~~~~~~~~~~~~'
    IF (lwp) WRITE(numout, *) '   Namelist namsbc_isf :'
    IF (lwp) WRITE(numout, *) '      type ice shelf melting/freezing         nn_isf      = ', nn_isf
    IF (lwp) WRITE(numout, *) '      bulk formulation (nn_isf=1 only)        nn_isfblk   = ', nn_isfblk
    IF (lwp) WRITE(numout, *) '      thickness of the top boundary layer     rn_hisf_tbl = ', rn_hisf_tbl
    IF (lwp) WRITE(numout, *) '      gamma formulation                       nn_gammablk = ', nn_gammablk
    IF (lwp) WRITE(numout, *) '      gammat coefficient                      rn_gammat0  = ', rn_gammat0
    IF (lwp) WRITE(numout, *) '      gammas coefficient                      rn_gammas0  = ', rn_gammas0
    IF (lwp) WRITE(numout, *) '      top drag coef. used (from namdrg_top)   rn_Cd0      = ', r_Cdmin_top
    IF (sbc_isf_alloc() /= 0) CALL ctl_stop('STOP', 'sbc_isf : unable to allocate arrays')
    !$ACC KERNELS
    qisf(:, :) = 0._wp
    fwfisf(:, :) = 0._wp
    risf_tsc(:, :, :) = 0._wp
    fwfisf_b(:, :) = 0._wp
    !$ACC END KERNELS
    SELECT CASE (nn_isf)
    CASE (1)
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '      ==>>>   presence of under iceshelf seas (nn_isf = 1)'
      !$ACC KERNELS
      rhisf_tbl(:, :) = rn_hisf_tbl
      misfkt(:, :) = mikt(:, :)
      !$ACC END KERNELS
    CASE (2, 3)
      IF (.NOT. l_isfcpl) THEN
        ALLOCATE(sf_rnfisf(1), STAT = ierror)
        ALLOCATE(sf_rnfisf(1) % fnow(jpi, jpj, 1), sf_rnfisf(1) % fdta(jpi, jpj, 1, 2))
        CALL fld_fill(sf_rnfisf, (/sn_rnfisf/), cn_dirisf, 'sbc_isf_init', 'read fresh water flux isf data', 'namsbc_isf')
      END IF
      IF (nn_isf == 2) THEN
        IF (lwp) WRITE(numout, *)
        IF (lwp) WRITE(numout, *) '      ==>>>   bg03 parametrisation (nn_isf = 2)'
        CALL iom_open(sn_Leff_isf % clname, inum)
        cvarLeff = TRIM(sn_Leff_isf % clvar)
        CALL iom_get(inum, jpdom_data, cvarLeff, risfLeff, 1)
        CALL iom_close(inum)
        risfLeff = risfLeff * 1000.0_wp
      ELSE
        IF (lwp) WRITE(numout, *)
        IF (lwp) WRITE(numout, *) '      ==>>>   rnf file for isf (nn_isf = 3)'
      END IF
      CALL iom_open(sn_depmax_isf % clname, inum)
      cvarhisf = TRIM(sn_depmax_isf % clvar)
      CALL iom_get(inum, jpdom_data, cvarhisf, rhisf_tbl, 1)
      CALL iom_close(inum)
      CALL iom_open(sn_depmin_isf % clname, inum)
      cvarzisf = TRIM(sn_depmin_isf % clvar)
      CALL iom_get(inum, jpdom_data, cvarzisf, rzisf_tbl, 1)
      CALL iom_close(inum)
      !$ACC KERNELS
      rhisf_tbl(:, :) = rhisf_tbl(:, :) - rzisf_tbl(:, :)
      !$ACC END KERNELS
      DO ji = 1, jpi
        DO jj = 1, jpj
          ik = 2
          DO WHILE (ik <= mbkt(ji, jj) .AND. gdepw_n(ji, jj, ik) < rzisf_tbl(ji, jj))
            ik = ik + 1
          END DO
          misfkt(ji, jj) = ik - 1
        END DO
      END DO
    CASE (4)
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '      ==>>>   specified fresh water flux in ISF (nn_isf = 4)'
      !$ACC KERNELS
      rhisf_tbl(:, :) = rn_hisf_tbl
      misfkt(:, :) = mikt(:, :)
      !$ACC END KERNELS
      IF (.NOT. l_isfcpl) THEN
        ALLOCATE(sf_fwfisf(1), STAT = ierror)
        ALLOCATE(sf_fwfisf(1) % fnow(jpi, jpj, 1), sf_fwfisf(1) % fdta(jpi, jpj, 1, 2))
        CALL fld_fill(sf_fwfisf, (/sn_fwfisf/), cn_dirisf, 'sbc_isf_init', 'read fresh water flux isf data', 'namsbc_isf')
      END IF
    CASE DEFAULT
      CALL ctl_stop('sbc_isf_init: wrong value of nn_isf')
    END SELECT
    !$ACC KERNELS
    rhisf_tbl_0(:, :) = rhisf_tbl(:, :)
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        ikt = misfkt(ji, jj)
        ikb = misfkt(ji, jj)
        rhisf_tbl(ji, jj) = MAX(rhisf_tbl_0(ji, jj), e3t_n(ji, jj, ikt))
        DO jk = ikt + 1, mbkt(ji, jj)
          IF ((SUM(e3t_n(ji, jj, ikt : jk - 1)) < rhisf_tbl(ji, jj)) .AND. (tmask(ji, jj, jk) == 1)) ikb = jk
        END DO
        rhisf_tbl(ji, jj) = MIN(rhisf_tbl(ji, jj), SUM(e3t_n(ji, jj, ikt : ikb)))
        misfkb(ji, jj) = ikb
        r1_hisf_tbl(ji, jj) = 1._wp / rhisf_tbl(ji, jj)
        zhk = SUM(e3t_n(ji, jj, ikt : ikb - 1)) * r1_hisf_tbl(ji, jj)
        ralpha(ji, jj) = rhisf_tbl(ji, jj) * (1._wp - zhk) / e3t_n(ji, jj, ikb)
      END DO
    END DO
    !$ACC END KERNELS
    IF (lwxios) THEN
      CALL iom_set_rstw_var_active('fwf_isf_b')
      CALL iom_set_rstw_var_active('isf_hc_b')
      CALL iom_set_rstw_var_active('isf_sc_b')
    END IF
  END SUBROUTINE sbc_isf_init
  SUBROUTINE sbc_isf_bg03(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: ik
    REAL(KIND = wp) :: zt_sum
    REAL(KIND = wp) :: zt_ave
    REAL(KIND = wp) :: zt_frz
    REAL(KIND = wp) :: zpress
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('sbc_isf_bg03', 'r0', 0, 0)
    DO ji = 1, jpi
      DO jj = 1, jpj
        ik = misfkt(ji, jj)
        zt_sum = 0.E0_wp
        IF (ik > 1) THEN
          DO jk = misfkt(ji, jj), misfkb(ji, jj)
            zpress = grav * rau0 * gdept_n(ji, jj, ik) * 1.E-04
            CALL eos_fzp(stbl(ji, jj), zt_frz, zpress)
            zt_sum = zt_sum + (tsn(ji, jj, jk, jp_tem) - zt_frz) * e3t_n(ji, jj, jk) * tmask(ji, jj, jk)
          END DO
          zt_ave = zt_sum / rhisf_tbl(ji, jj)
          qisf(ji, jj) = - rau0 * rcp * rn_gammat0 * risfLeff(ji, jj) * e1t(ji, jj) * zt_ave * r1_e1e2t(ji, jj) * tmask(ji, jj, jk)
          fwfisf(ji, jj) = qisf(ji, jj) / rLfusisf
          fwfisf(ji, jj) = fwfisf(ji, jj) * (soce / stbl(ji, jj))
        ELSE
          qisf(ji, jj) = 0._wp
          fwfisf(ji, jj) = 0._wp
        END IF
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE sbc_isf_bg03
  SUBROUTINE sbc_isf_cav(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj
    INTEGER :: nit
    LOGICAL :: lit
    REAL(KIND = wp) :: zlamb1, zlamb2, zlamb3
    REAL(KIND = wp) :: zeps1, zeps2, zeps3, zeps4, zeps6, zeps7
    REAL(KIND = wp) :: zaqe, zbqe, zcqe, zaqer, zdis, zsfrz, zcfac
    REAL(KIND = wp) :: zeps = 1.E-20_wp
    REAL(KIND = wp) :: zerr
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zfrz
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zgammat, zgammas
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zfwflx, zhtflx, zhtflx_b
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('sbc_isf_cav', 'r0', 0, 0)
    IF (l_useCT) THEN
      zlamb1 = - 0.0564_wp
      zlamb2 = 0.0773_wp
      zlamb3 = - 7.8633E-8 * grav * rau0
    ELSE
      zlamb1 = - 0.0573_wp
      zlamb2 = 0.0832_wp
      zlamb3 = - 7.53E-8 * grav * rau0
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zgammat(:, :) = rn_gammat0
    zgammas(:, :) = rn_gammas0
    zhtflx(:, :) = 0.0_wp
    zhtflx_b(:, :) = 0.0_wp
    zfwflx(:, :) = 0.0_wp
    nit = 1
    lit = .TRUE.
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('sbc_isf_cav', 'r1', 0, 0)
    DO WHILE (lit)
      SELECT CASE (nn_isfblk)
      CASE (1)
        CALL eos_fzp(stbl(:, :), zfrz(:, :), risfdep(:, :))
        CALL sbc_isf_gammats(zgammat, zgammas, zhtflx, zfwflx)
        DO jj = 1, jpj
          DO ji = 1, jpi
            zhtflx(ji, jj) = zgammat(ji, jj) * rcp * rau0 * (ttbl(ji, jj) - zfrz(ji, jj))
            zfwflx(ji, jj) = - zhtflx(ji, jj) / rLfusisf
          END DO
        END DO
        qisf(:, :) = - zhtflx(:, :) * (1._wp - tmask(:, :, 1)) * ssmask(:, :)
        fwfisf(:, :) = zfwflx(:, :) * (1._wp - tmask(:, :, 1)) * ssmask(:, :)
      CASE (2)
        CALL sbc_isf_gammats(zgammat, zgammas, zhtflx, zfwflx)
        DO jj = 1, jpj
          DO ji = 1, jpi
            zeps1 = rcp * rau0 * zgammat(ji, jj)
            zeps2 = rLfusisf * rau0 * zgammas(ji, jj)
            zeps3 = rhoisf * rcpisf * rkappa / MAX(risfdep(ji, jj), zeps)
            zeps4 = zlamb2 + zlamb3 * risfdep(ji, jj)
            zeps6 = zeps4 - ttbl(ji, jj)
            zeps7 = zeps4 - tsurf
            zaqe = zlamb1 * (zeps1 + zeps3)
            zaqer = 0.5_wp / MIN(zaqe, - zeps)
            zbqe = zeps1 * zeps6 + zeps3 * zeps7 - zeps2
            zcqe = zeps2 * stbl(ji, jj)
            zdis = zbqe * zbqe - 4.0_wp * zaqe * zcqe
            zsfrz = (- zbqe - SQRT(zdis)) * zaqer
            IF (zsfrz < 0.0_wp) zsfrz = (- zbqe + SQRT(zdis)) * zaqer
            zfrz(ji, jj) = zeps4 + zlamb1 * zsfrz
            zfwflx(ji, jj) = rau0 * zgammas(ji, jj) * (zsfrz - stbl(ji, jj)) / MAX(zsfrz, zeps)
            zhtflx(ji, jj) = zgammat(ji, jj) * rau0 * rcp * (ttbl(ji, jj) - zfrz(ji, jj))
          END DO
        END DO
        qisf(:, :) = - zhtflx(:, :) * (1._wp - tmask(:, :, 1)) * ssmask(:, :)
        fwfisf(:, :) = zfwflx(:, :) * (1._wp - tmask(:, :, 1)) * ssmask(:, :)
      END SELECT
      IF (nn_gammablk < 2) THEN
        lit = .FALSE.
      ELSE
        IF (nit >= 100) THEN
          CALL ctl_stop('STOP', 'sbc_isf_hol99 : too many iteration ...')
        ELSE
          nit = nit + 1
        END IF
        zerr = MAXVAL(ABS(zhtflx - zhtflx_b))
        IF (zerr <= 0.01_wp) THEN
          lit = .FALSE.
        ELSE
          zhtflx_b(:, :) = zhtflx(:, :)
        END IF
      END IF
    END DO
    CALL iom_put('isfgammat', zgammat)
    CALL iom_put('isfgammas', zgammas)
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE sbc_isf_cav
  SUBROUTINE sbc_isf_gammats(pgt, pgs, pqhisf, pqwisf)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: pgt, pgs
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pqhisf, pqwisf
    INTEGER :: ji, jj
    INTEGER :: ikt
    REAL(KIND = wp) :: zdku, zdkv
    REAL(KIND = wp) :: zPr, zSc, zRc
    REAL(KIND = wp) :: zmob, zmols
    REAL(KIND = wp) :: zbuofdep, zhnu
    REAL(KIND = wp) :: zhmax
    REAL(KIND = wp) :: zetastar
    REAL(KIND = wp) :: zgmolet, zgmoles, zgturb
    REAL(KIND = wp) :: zcoef
    REAL(KIND = wp) :: zdep
    REAL(KIND = wp) :: zeps = 1.0E-20_wp
    REAL(KIND = wp), PARAMETER :: zxsiN = 0.052_wp
    REAL(KIND = wp), PARAMETER :: znu = 1.95E-6_wp
    REAL(KIND = wp), DIMENSION(2) :: zts, zab
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zustar
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    SELECT CASE (nn_gammablk)
    CASE (0)
      !$ACC KERNELS
      pgt(:, :) = rn_gammat0
      pgs(:, :) = rn_gammas0
      !$ACC END KERNELS
    CASE (1)
      !$ACC KERNELS
      zustar(:, :) = SQRT(r_Cdmin_top * (utbl(:, :) * utbl(:, :) + vtbl(:, :) * vtbl(:, :) + r_ke0_top))
      pgt(:, :) = zustar(:, :) * rn_gammat0
      pgs(:, :) = zustar(:, :) * rn_gammas0
      !$ACC END KERNELS
    CASE (2)
      !$ACC KERNELS
      zustar(:, :) = SQRT(r_Cdmin_top * (utbl(:, :) * utbl(:, :) + vtbl(:, :) * vtbl(:, :) + r_ke0_top))
      zPr = 13.8_wp
      zSc = 2432.0_wp
      zgmolet = 12.5_wp * zPr ** (2.0 / 3.0) - 6.0_wp
      zgmoles = 12.5_wp * zSc ** (2.0 / 3.0) - 6.0_wp
      !$ACC END KERNELS
      CALL profile_psy_data0 % PreStart('sbc_isf_gammats', 'r0', 0, 0)
      DO ji = 2, jpi
        DO jj = 2, jpj
          ikt = mikt(ji, jj)
          IF (zustar(ji, jj) == 0._wp) THEN
            pgt = rn_gammat0
            pgs = rn_gammas0
          ELSE
            zcoef = 0.5_wp / e3w_n(ji, jj, ikt + 1)
            zdku = zcoef * (un(ji - 1, jj, ikt) + un(ji, jj, ikt) - un(ji - 1, jj, ikt + 1) - un(ji, jj, ikt + 1))
            zdkv = zcoef * (vn(ji, jj - 1, ikt) + vn(ji, jj, ikt) - vn(ji, jj - 1, ikt + 1) - vn(ji, jj, ikt + 1))
            zRc = rn2(ji, jj, ikt + 1) / MAX(zdku * zdku + zdkv * zdkv, zeps)
            zts(jp_tem) = ttbl(ji, jj)
            zts(jp_sal) = stbl(ji, jj)
            zdep = gdepw_n(ji, jj, ikt)
            CALL eos_rab(zts, zdep, zab)
            zbuofdep = grav * (zab(jp_tem) * pqhisf(ji, jj) - zab(jp_sal) * pqwisf(ji, jj))
            zhmax = gdept_n(ji, jj, mbkt(ji, jj)) - gdepw_n(ji, jj, mikt(ji, jj)) - 0.001_wp
            zmob = zustar(ji, jj) ** 3 / (vkarmn * (zbuofdep + zeps))
            zmols = SIGN(1._wp, zmob) * MIN(ABS(zmob), zhmax) * tmask(ji, jj, ikt)
            zetastar = 1._wp / (SQRT(1._wp + MAX(zxsiN * zustar(ji, jj) / (ABS(ff_f(ji, jj)) * zmols * zRc), 0._wp)))
            zhnu = 5 * znu / zustar(ji, jj)
            zgturb = 1._wp / vkarmn * LOG(zustar(ji, jj) * zxsiN * zetastar * zetastar / (ABS(ff_f(ji, jj)) * zhnu)) + 1._wp / (2 * zxsiN * zetastar) - 1._wp / vkarmn
            pgt(ji, jj) = zustar(ji, jj) / (zgturb + zgmolet)
            pgs(ji, jj) = zustar(ji, jj) / (zgturb + zgmoles)
          END IF
        END DO
      END DO
      CALL lbc_lnk_multi('sbcisf', pgt, 'T', 1., pgs, 'T', 1.)
      CALL profile_psy_data0 % PostEnd
    END SELECT
  END SUBROUTINE sbc_isf_gammats
  SUBROUTINE sbc_isf_tbl(pvarin, pvarout, cd_ptin)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pvarin
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: pvarout
    CHARACTER(LEN = 1), INTENT(IN) :: cd_ptin
    INTEGER :: ji, jj, jk
    INTEGER :: ikt, ikb
    REAL(KIND = wp) :: ze3, zhk
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zhisf_tbl
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    !$ACC KERNELS
    pvarout(:, :) = 0._wp
    !$ACC END KERNELS
    SELECT CASE (cd_ptin)
    CASE ('U')
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          ikt = miku(ji, jj)
          ikb = miku(ji, jj)
          zhisf_tbl(ji, jj) = MAX(rhisf_tbl_0(ji, jj), e3u_n(ji, jj, ikt))
          DO jk = ikt + 1, mbku(ji, jj)
            IF ((SUM(e3u_n(ji, jj, ikt : jk - 1)) < zhisf_tbl(ji, jj)) .AND. (umask(ji, jj, jk) == 1)) ikb = jk
          END DO
          zhisf_tbl(ji, jj) = MIN(zhisf_tbl(ji, jj), SUM(e3u_n(ji, jj, ikt : ikb)))
          DO jk = ikt, ikb - 1
            ze3 = e3u_n(ji, jj, jk)
            pvarout(ji, jj) = pvarout(ji, jj) + pvarin(ji, jj, jk) / zhisf_tbl(ji, jj) * ze3
          END DO
          zhk = SUM(e3u_n(ji, jj, ikt : ikb - 1)) / zhisf_tbl(ji, jj)
          pvarout(ji, jj) = pvarout(ji, jj) + pvarin(ji, jj, ikb) * (1._wp - zhk)
        END DO
      END DO
      !$ACC loop independent collapse(2)
      DO jj = 2, jpj
        DO ji = 2, jpi
          pvarout(ji, jj) = 0.5_wp * (pvarout(ji, jj) + pvarout(ji - 1, jj))
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data0 % PreStart('sbc_isf_tbl', 'r0', 0, 0)
      CALL lbc_lnk('sbcisf', pvarout, 'T', - 1.)
      CALL profile_psy_data0 % PostEnd
    CASE ('V')
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          ikt = mikv(ji, jj)
          ikb = mikv(ji, jj)
          zhisf_tbl(ji, jj) = MAX(rhisf_tbl_0(ji, jj), e3v_n(ji, jj, ikt))
          DO jk = ikt + 1, mbkv(ji, jj)
            IF ((SUM(e3v_n(ji, jj, ikt : jk - 1)) < zhisf_tbl(ji, jj)) .AND. (vmask(ji, jj, jk) == 1)) ikb = jk
          END DO
          zhisf_tbl(ji, jj) = MIN(zhisf_tbl(ji, jj), SUM(e3v_n(ji, jj, ikt : ikb)))
          DO jk = ikt, ikb - 1
            ze3 = e3v_n(ji, jj, jk)
            pvarout(ji, jj) = pvarout(ji, jj) + pvarin(ji, jj, jk) / zhisf_tbl(ji, jj) * ze3
          END DO
          zhk = SUM(e3v_n(ji, jj, ikt : ikb - 1)) / zhisf_tbl(ji, jj)
          pvarout(ji, jj) = pvarout(ji, jj) + pvarin(ji, jj, ikb) * (1._wp - zhk)
        END DO
      END DO
      !$ACC loop independent collapse(2)
      DO jj = 2, jpj
        DO ji = 2, jpi
          pvarout(ji, jj) = 0.5_wp * (pvarout(ji, jj) + pvarout(ji, jj - 1))
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('sbc_isf_tbl', 'r1', 0, 0)
      CALL lbc_lnk('sbcisf', pvarout, 'T', - 1.)
      CALL profile_psy_data1 % PostEnd
    CASE ('T')
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          ikt = misfkt(ji, jj)
          ikb = misfkb(ji, jj)
          DO jk = ikt, ikb - 1
            ze3 = e3t_n(ji, jj, jk)
            pvarout(ji, jj) = pvarout(ji, jj) + pvarin(ji, jj, jk) * r1_hisf_tbl(ji, jj) * ze3
          END DO
          zhk = SUM(e3t_n(ji, jj, ikt : ikb - 1)) * r1_hisf_tbl(ji, jj)
          pvarout(ji, jj) = pvarout(ji, jj) + pvarin(ji, jj, ikb) * (1._wp - zhk)
        END DO
      END DO
      !$ACC END KERNELS
    END SELECT
    !$ACC KERNELS
    pvarout(:, :) = pvarout(:, :) * ssmask(:, :)
    !$ACC END KERNELS
  END SUBROUTINE sbc_isf_tbl
  SUBROUTINE sbc_isf_div(phdivn)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: phdivn
    INTEGER :: ji, jj, jk
    INTEGER :: ikt, ikb
    REAL(KIND = wp) :: zhk
    REAL(KIND = wp) :: zfact
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('sbc_isf_div', 'r0', 0, 0)
    zfact = 0.5_wp
    CALL profile_psy_data0 % PostEnd
    IF (.NOT. ln_linssh) THEN
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          ikt = misfkt(ji, jj)
          ikb = misfkt(ji, jj)
          rhisf_tbl(ji, jj) = MAX(rhisf_tbl_0(ji, jj), e3t_n(ji, jj, ikt))
          DO jk = ikt, mbkt(ji, jj)
            IF ((SUM(e3t_n(ji, jj, ikt : jk - 1)) .LT. rhisf_tbl(ji, jj)) .AND. (tmask(ji, jj, jk) == 1)) ikb = jk
          END DO
          rhisf_tbl(ji, jj) = MIN(rhisf_tbl(ji, jj), SUM(e3t_n(ji, jj, ikt : ikb)))
          misfkb(ji, jj) = ikb
          r1_hisf_tbl(ji, jj) = 1._wp / rhisf_tbl(ji, jj)
          zhk = SUM(e3t_n(ji, jj, ikt : ikb - 1)) * r1_hisf_tbl(ji, jj)
          ralpha(ji, jj) = rhisf_tbl(ji, jj) * (1._wp - zhk) / e3t_n(ji, jj, ikb)
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        ikt = misfkt(ji, jj)
        ikb = misfkb(ji, jj)
        DO jk = ikt, ikb - 1
          phdivn(ji, jj, jk) = phdivn(ji, jj, jk) + (fwfisf(ji, jj) + fwfisf_b(ji, jj)) * r1_hisf_tbl(ji, jj) * r1_rau0 * zfact
        END DO
        phdivn(ji, jj, ikb) = phdivn(ji, jj, ikb) + (fwfisf(ji, jj) + fwfisf_b(ji, jj)) * r1_hisf_tbl(ji, jj) * r1_rau0 * zfact * ralpha(ji, jj)
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE sbc_isf_div
END MODULE sbcisf