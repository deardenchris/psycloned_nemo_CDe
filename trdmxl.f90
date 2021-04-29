MODULE trdmxl
  USE oce
  USE dom_oce
  USE trd_oce
  USE trdmxl_oce
  USE ldftra
  USE zdf_oce
  USE phycst
  USE dianam
  USE ldfslp
  USE zdfmxl
  USE zdfddm
  USE lbclnk
  USE trdmxl_rst
  USE in_out_manager
  USE ioipsl
  USE prtctl
  USE restart
  USE lib_mpp
  USE iom
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: trd_mxl
  PUBLIC :: trd_mxl_init
  PUBLIC :: trd_mxl_zint
  INTEGER :: nkstp
  CHARACTER(LEN = 40) :: clhstnam
  INTEGER :: nh_t, nmoymltrd
  INTEGER :: nidtrd
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: ndextrd1
  INTEGER :: ndimtrd1
  INTEGER :: ionce, icount
  CONTAINS
  INTEGER FUNCTION trd_mxl_alloc()
    ALLOCATE(ndextrd1(jpi * jpj), STAT = trd_mxl_alloc)
    CALL mpp_sum('trdmxl', trd_mxl_alloc)
    IF (trd_mxl_alloc /= 0) CALL ctl_warn('trd_mxl_alloc: failed to allocate array ndextrd1')
  END FUNCTION trd_mxl_alloc
  SUBROUTINE trd_tra_mxl(ptrdx, ptrdy, ktrd, kt, p2dt, kmxln)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: ptrdx
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: ptrdy
    INTEGER, INTENT(IN) :: ktrd
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), INTENT(IN) :: p2dt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: kmxln
    INTEGER :: ji, jj, jk
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    IF (kt /= nkstp) THEN
      nkstp = kt
      tmltrd(:, :, :) = 0._wp
      smltrd(:, :, :) = 0._wp
      wkx(:, :, :) = 0._wp
      DO jk = 1, jpktrd
        !$ACC loop independent collapse(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            IF (jk - kmxln(ji, jj) < 0) wkx(ji, jj, jk) = e3t_n(ji, jj, jk) * tmask(ji, jj, jk)
          END DO
        END DO
      END DO
      hmxl(:, :) = 0._wp
      DO jk = 1, jpktrd
        hmxl(:, :) = hmxl(:, :) + wkx(:, :, jk)
      END DO
      DO jk = 1, jpktrd
        wkx(:, :, jk) = wkx(:, :, jk) / MAX(1.E-20_wp, hmxl(:, :)) * tmask(:, :, 1)
      END DO
      tml(:, :) = 0._wp
      sml(:, :) = 0._wp
      DO jk = 1, jpktrd
        tml(:, :) = tml(:, :) + wkx(:, :, jk) * tsn(:, :, jk, jp_tem)
        sml(:, :) = sml(:, :) + wkx(:, :, jk) * tsn(:, :, jk, jp_sal)
      END DO
    END IF
    tmltrd(:, :, ktrd) = tmltrd(:, :, ktrd) + ptrdx(:, :, jk) * wkx(:, :, jk)
    smltrd(:, :, ktrd) = smltrd(:, :, ktrd) + ptrdy(:, :, jk) * wkx(:, :, jk)
    !$ACC END KERNELS
    SELECT CASE (ktrd)
    CASE (jptra_npc)
    CASE (jptra_zdfp)
      !$ACC KERNELS
      tmltrd(:, :, jpmxl_ldf) = tmltrd(:, :, jpmxl_ldf) + (tmltrd(:, :, jpmxl_zdf) - tmltrd(:, :, jpmxl_zdfp))
      smltrd(:, :, jpmxl_ldf) = smltrd(:, :, jpmxl_ldf) + (smltrd(:, :, jpmxl_zdf) - smltrd(:, :, jpmxl_zdfp))
      tmltrd(:, :, jpmxl_zdf) = tmltrd(:, :, jpmxl_zdfp)
      smltrd(:, :, jpmxl_zdf) = smltrd(:, :, jpmxl_zdfp)
      !$ACC END KERNELS
      IF (ln_zdfnpc) THEN
        !$ACC KERNELS
        tmltrd(:, :, jpmxl_zdf) = tmltrd(:, :, jpmxl_zdf) + tmltrd(:, :, jpmxl_npc)
        smltrd(:, :, jpmxl_zdf) = smltrd(:, :, jpmxl_zdf) + smltrd(:, :, jpmxl_npc)
        !$ACC END KERNELS
      END IF
    CASE (jptra_atf)
      CALL profile_psy_data0 % PreStart('trd_tra_mxl', 'r0', 0, 0)
      CALL profile_psy_data0 % PostEnd
    END SELECT
  END SUBROUTINE trd_tra_mxl
  SUBROUTINE trd_mean(kt, ptrd, ptrdm)
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: ptrd
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: ptrdm
    INTEGER, INTENT(IN) :: kt
    !$ACC KERNELS
    IF (kt == nn_it000) ptrdm(:, :, :) = 0._wp
    ptrdm(:, :, :) = ptrdm(:, :, :) + ptrd(:, :, :)
    IF (MOD(kt - nn_it000 + 1, nn_trd) == 0) THEN
    END IF
    !$ACC END KERNELS
  END SUBROUTINE trd_mean
  SUBROUTINE trd_mxl_zint(pttrdmxl, pstrdmxl, ktrd, ctype)
    INTEGER, INTENT(IN) :: ktrd
    CHARACTER(LEN = 2), INTENT(IN) :: ctype
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pttrdmxl
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pstrdmxl
    INTEGER :: ji, jj, jk, isum
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zvlmsk
    !$ACC KERNELS
    IF (icount == 1) THEN
      IF (nn_ctls == 0) THEN
        nmxl(:, :) = nmln(:, :)
      ELSE IF (nn_ctls == 1) THEN
        nmxl(:, :) = nbol(:, :)
      ELSE IF (nn_ctls >= 2) THEN
        nn_ctls = MIN(nn_ctls, jpktrd - 1)
        nmxl(:, :) = nn_ctls + 1
      END IF
    END IF
    !$ACC END KERNELS
  END SUBROUTINE trd_mxl_zint
  SUBROUTINE trd_mxl(kt, p2dt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), INTENT(IN) :: p2dt
    INTEGER :: ji, jj, jk, jl, ik, it, itmod
    LOGICAL :: lldebug = .TRUE.
    REAL(KIND = wp) :: zavt, zfn, zfn2
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztmltot, zsmltot, ztmlres, zsmlres, ztmlatf, zsmlatf
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztmltot2, zsmltot2, ztmlres2, zsmlres2, ztmlatf2, zsmlatf2, ztmltrdm2, zsmltrdm2
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: ztmltrd2, zsmltrd2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data7
    !$ACC KERNELS
    ztmltrd2(:, :, :) = 0.E0
    zsmltrd2(:, :, :) = 0.E0
    ztmltot2(:, :) = 0.E0
    zsmltot2(:, :) = 0.E0
    ztmlres2(:, :) = 0.E0
    zsmlres2(:, :) = 0.E0
    ztmlatf2(:, :) = 0.E0
    zsmlatf2(:, :) = 0.E0
    IF (kt > nit000) THEN
      tmlb(:, :) = tml(:, :)
      smlb(:, :) = sml(:, :)
      tmlatfn(:, :) = tmltrd(:, :, jpmxl_atf)
      smlatfn(:, :) = smltrd(:, :, jpmxl_atf)
    END IF
    !$ACC END KERNELS
    IF (kt == 2) THEN
      !$ACC KERNELS
      tmlbb(:, :) = tmlb(:, :)
      smlbb(:, :) = smlb(:, :)
      tmlbn(:, :) = tml(:, :)
      smlbn(:, :) = sml(:, :)
      tmlatfb(:, :) = tmlatfn(:, :)
      smlatfb(:, :) = smlatfn(:, :)
      tmltrd_csum_ub(:, :, :) = 0.E0
      smltrd_csum_ub(:, :, :) = 0.E0
      tmltrd_atf_sumb(:, :) = 0.E0
      smltrd_atf_sumb(:, :) = 0.E0
      hmxlbn(:, :) = hmxl(:, :)
      !$ACC END KERNELS
      CALL profile_psy_data0 % PreStart('trd_mxl', 'r0', 0, 0)
      IF (ln_ctl) THEN
        WRITE(numout, *) '             we reach kt == nit000 + 1 = ', nit000 + 1
        CALL prt_ctl(tab2d_1 = tmlbb, clinfo1 = ' tmlbb   -   : ', mask1 = tmask)
        CALL prt_ctl(tab2d_1 = tmlbn, clinfo1 = ' tmlbn   -   : ', mask1 = tmask)
        CALL prt_ctl(tab2d_1 = tmlatfb, clinfo1 = ' tmlatfb -   : ', mask1 = tmask)
      END IF
      CALL profile_psy_data0 % PostEnd
    END IF
    CALL profile_psy_data1 % PreStart('trd_mxl', 'r1', 0, 0)
    IF ((ln_rstart) .AND. (kt == nit000) .AND. (ln_ctl)) THEN
      IF (ln_trdmxl_instant) THEN
        WRITE(numout, *) '             restart from kt == nit000 = ', nit000
        CALL prt_ctl(tab2d_1 = tmlbb, clinfo1 = ' tmlbb   -   : ', mask1 = tmask)
        CALL prt_ctl(tab2d_1 = tmlbn, clinfo1 = ' tmlbn   -   : ', mask1 = tmask)
        CALL prt_ctl(tab2d_1 = tmlatfb, clinfo1 = ' tmlatfb -   : ', mask1 = tmask)
      ELSE
        WRITE(numout, *) '             restart from kt == nit000 = ', nit000
        CALL prt_ctl(tab2d_1 = tmlbn, clinfo1 = ' tmlbn           -  : ', mask1 = tmask)
        CALL prt_ctl(tab2d_1 = hmxlbn, clinfo1 = ' hmxlbn          -  : ', mask1 = tmask)
        CALL prt_ctl(tab2d_1 = tml_sumb, clinfo1 = ' tml_sumb        -  : ', mask1 = tmask)
        CALL prt_ctl(tab2d_1 = tmltrd_atf_sumb, clinfo1 = ' tmltrd_atf_sumb -  : ', mask1 = tmask)
        CALL prt_ctl(tab3d_1 = tmltrd_csum_ub, clinfo1 = ' tmltrd_csum_ub  -  : ', mask1 = tmask, kdim = 1)
      END IF
    END IF
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    IF ((kt >= 2) .OR. (ln_rstart)) THEN
      nmoymltrd = nmoymltrd + 1
      DO jl = 1, jpltrd
        tmltrdm(:, :) = tmltrdm(:, :) + tmltrd(:, :, jl)
        smltrdm(:, :) = smltrdm(:, :) + smltrd(:, :, jl)
      END DO
      tmlatfm(:, :) = tmlatfm(:, :) + tmlatfn(:, :)
      smlatfm(:, :) = smlatfm(:, :) + smlatfn(:, :)
      tmltrd_sum(:, :, :) = tmltrd_sum(:, :, :) + tmltrd(:, :, :)
      tmltrd_csum_ln(:, :, :) = tmltrd_csum_ln(:, :, :) + tmltrd_sum(:, :, :)
      tml_sum(:, :) = tml_sum(:, :) + tml(:, :)
      smltrd_sum(:, :, :) = smltrd_sum(:, :, :) + smltrd(:, :, :)
      smltrd_csum_ln(:, :, :) = smltrd_csum_ln(:, :, :) + smltrd_sum(:, :, :)
      sml_sum(:, :) = sml_sum(:, :) + sml(:, :)
      hmxl_sum(:, :) = hmxl_sum(:, :) + hmxl(:, :)
    END IF
    tmltrd(:, :, :) = tmltrd(:, :, :) * rn_ucf
    smltrd(:, :, :) = smltrd(:, :, :) * rn_ucf
    it = kt
    itmod = kt - nit000 + 1
    !$ACC END KERNELS
    MODULO_NTRD:IF (MOD(itmod, nn_trd) == 0) THEN
      !$ACC KERNELS
      ztmltot(:, :) = 0.E0
      zsmltot(:, :) = 0.E0
      ztmlres(:, :) = 0.E0
      zsmlres(:, :) = 0.E0
      ztmltot2(:, :) = 0.E0
      zsmltot2(:, :) = 0.E0
      ztmlres2(:, :) = 0.E0
      zsmlres2(:, :) = 0.E0
      zfn = REAL(nmoymltrd, wp)
      zfn2 = zfn * zfn
      ztmltot(:, :) = (tml(:, :) - tmlbn(:, :) + tmlb(:, :) - tmlbb(:, :)) / p2dt
      zsmltot(:, :) = (sml(:, :) - smlbn(:, :) + smlb(:, :) - smlbb(:, :)) / p2dt
      ztmlres(:, :) = ztmltot(:, :) - (tmltrdm(:, :) - tmlatfn(:, :) + tmlatfb(:, :))
      zsmlres(:, :) = zsmltot(:, :) - (smltrdm(:, :) - smlatfn(:, :) + smlatfb(:, :))
      ztmlatf(:, :) = tmlatfm(:, :) - tmlatfn(:, :) + tmlatfb(:, :)
      zsmlatf(:, :) = smlatfm(:, :) - smlatfn(:, :) + smlatfb(:, :)
      !$ACC END KERNELS
      CALL profile_psy_data2 % PreStart('trd_mxl', 'r2', 0, 0)
      CALL lbc_lnk_multi('trdmxl', ztmltot, 'T', 1., zsmltot, 'T', 1., ztmlres, 'T', 1., zsmlres, 'T', 1., ztmlatf, 'T', 1., zsmlatf, 'T', 1.)
      CALL profile_psy_data2 % PostEnd
      !$ACC KERNELS
      hmxl_sum(:, :) = hmxlbn(:, :) + 2 * (hmxl_sum(:, :) - hmxl(:, :)) + hmxl(:, :)
      tml_sum(:, :) = tmlbn(:, :) + 2 * (tml_sum(:, :) - tml(:, :)) + tml(:, :)
      ztmltot2(:, :) = (tml_sum(:, :) - tml_sumb(:, :)) / p2dt
      sml_sum(:, :) = smlbn(:, :) + 2 * (sml_sum(:, :) - sml(:, :)) + sml(:, :)
      zsmltot2(:, :) = (sml_sum(:, :) - sml_sumb(:, :)) / p2dt
      DO jl = 1, jpltrd
        ztmltrd2(:, :, jl) = tmltrd_csum_ub(:, :, jl) + tmltrd_csum_ln(:, :, jl)
      END DO
      ztmltrdm2(:, :) = 0.E0
      DO jl = 1, jpltrd
        ztmltrdm2(:, :) = ztmltrdm2(:, :) + ztmltrd2(:, :, jl)
      END DO
      ztmlres2(:, :) = ztmltot2(:, :) - (ztmltrdm2(:, :) - tmltrd_sum(:, :, jpmxl_atf) + tmltrd_atf_sumb(:, :))
      DO jl = 1, jpltrd
        zsmltrd2(:, :, jl) = smltrd_csum_ub(:, :, jl) + smltrd_csum_ln(:, :, jl)
      END DO
      zsmltrdm2(:, :) = 0.
      DO jl = 1, jpltrd
        zsmltrdm2(:, :) = zsmltrdm2(:, :) + zsmltrd2(:, :, jl)
      END DO
      zsmlres2(:, :) = zsmltot2(:, :) - (zsmltrdm2(:, :) - smltrd_sum(:, :, jpmxl_atf) + smltrd_atf_sumb(:, :))
      ztmlatf2(:, :) = ztmltrd2(:, :, jpmxl_atf) - tmltrd_sum(:, :, jpmxl_atf) + tmltrd_atf_sumb(:, :)
      zsmlatf2(:, :) = zsmltrd2(:, :, jpmxl_atf) - smltrd_sum(:, :, jpmxl_atf) + smltrd_atf_sumb(:, :)
      !$ACC END KERNELS
      CALL profile_psy_data3 % PreStart('trd_mxl', 'r3', 0, 0)
      CALL lbc_lnk_multi('trdmxl', ztmltot2, 'T', 1., zsmltot2, 'T', 1., ztmlres2, 'T', 1., zsmlres2, 'T', 1.)
      CALL lbc_lnk_multi('trdmxl', ztmltrd2(:, :, :), 'T', 1., zsmltrd2(:, :, :), 'T', 1.)
      CALL profile_psy_data3 % PostEnd
      !$ACC KERNELS
      tmlbb(:, :) = tmlb(:, :)
      smlbb(:, :) = smlb(:, :)
      tmlbn(:, :) = tml(:, :)
      smlbn(:, :) = sml(:, :)
      tmlatfb(:, :) = tmlatfn(:, :)
      smlatfb(:, :) = smlatfn(:, :)
      tmltrd_csum_ub(:, :, :) = zfn * tmltrd_sum(:, :, :) - tmltrd_csum_ln(:, :, :)
      tml_sumb(:, :) = tml_sum(:, :)
      tmltrd_atf_sumb(:, :) = tmltrd_sum(:, :, jpmxl_atf)
      smltrd_csum_ub(:, :, :) = zfn * smltrd_sum(:, :, :) - smltrd_csum_ln(:, :, :)
      sml_sumb(:, :) = sml_sum(:, :)
      smltrd_atf_sumb(:, :) = smltrd_sum(:, :, jpmxl_atf)
      hmxlbn(:, :) = hmxl(:, :)
      !$ACC END KERNELS
      CALL profile_psy_data4 % PreStart('trd_mxl', 'r4', 0, 0)
      IF (ln_ctl) THEN
        IF (ln_trdmxl_instant) THEN
          CALL prt_ctl(tab2d_1 = tmlbb, clinfo1 = ' tmlbb   -   : ', mask1 = tmask)
          CALL prt_ctl(tab2d_1 = tmlbn, clinfo1 = ' tmlbn   -   : ', mask1 = tmask)
          CALL prt_ctl(tab2d_1 = tmlatfb, clinfo1 = ' tmlatfb -   : ', mask1 = tmask)
        ELSE
          CALL prt_ctl(tab2d_1 = tmlbn, clinfo1 = ' tmlbn           -  : ', mask1 = tmask)
          CALL prt_ctl(tab2d_1 = hmxlbn, clinfo1 = ' hmxlbn          -  : ', mask1 = tmask)
          CALL prt_ctl(tab2d_1 = tml_sumb, clinfo1 = ' tml_sumb        -  : ', mask1 = tmask)
          CALL prt_ctl(tab2d_1 = tmltrd_atf_sumb, clinfo1 = ' tmltrd_atf_sumb -  : ', mask1 = tmask)
          CALL prt_ctl(tab3d_1 = tmltrd_csum_ub, clinfo1 = ' tmltrd_csum_ub  -  : ', mask1 = tmask, kdim = 1)
        END IF
      END IF
      CALL profile_psy_data4 % PostEnd
      !$ACC KERNELS
      ztmltot(:, :) = ztmltot(:, :) * rn_ucf / zfn
      zsmltot(:, :) = zsmltot(:, :) * rn_ucf / zfn
      ztmlres(:, :) = ztmlres(:, :) * rn_ucf / zfn
      zsmlres(:, :) = zsmlres(:, :) * rn_ucf / zfn
      ztmlatf(:, :) = ztmlatf(:, :) * rn_ucf / zfn
      zsmlatf(:, :) = zsmlatf(:, :) * rn_ucf / zfn
      tml_sum(:, :) = tml_sum(:, :) / (2 * zfn)
      sml_sum(:, :) = sml_sum(:, :) / (2 * zfn)
      ztmltot2(:, :) = ztmltot2(:, :) * rn_ucf / zfn2
      zsmltot2(:, :) = zsmltot2(:, :) * rn_ucf / zfn2
      ztmltrd2(:, :, :) = ztmltrd2(:, :, :) * rn_ucf / zfn2
      zsmltrd2(:, :, :) = zsmltrd2(:, :, :) * rn_ucf / zfn2
      ztmlatf2(:, :) = ztmlatf2(:, :) * rn_ucf / zfn2
      zsmlatf2(:, :) = zsmlatf2(:, :) * rn_ucf / zfn2
      ztmlres2(:, :) = ztmlres2(:, :) * rn_ucf / zfn2
      zsmlres2(:, :) = zsmlres2(:, :) * rn_ucf / zfn2
      hmxl_sum(:, :) = hmxl_sum(:, :) / (2 * zfn)
      !$ACC END KERNELS
      CALL profile_psy_data5 % PreStart('trd_mxl', 'r5', 0, 0)
      IF (lldebug) THEN
        WRITE(numout, *)
        WRITE(numout, *) 'trd_mxl : write trends in the Mixed Layer for debugging process:'
        WRITE(numout, *) '~~~~~~~  '
        WRITE(numout, *) '          TRA kt = ', kt, 'nmoymltrd = ', nmoymltrd
        WRITE(numout, *)
        WRITE(numout, *) '          >>>>>>>>>>>>>>>>>>  TRA TEMPERATURE <<<<<<<<<<<<<<<<<<'
        WRITE(numout, *) '          TRA ztmlres    : ', SUM(ztmlres(:, :))
        WRITE(numout, *) '          TRA ztmltot    : ', SUM(ztmltot(:, :))
        WRITE(numout, *) '          TRA tmltrdm    : ', SUM(tmltrdm(:, :))
        WRITE(numout, *) '          TRA tmlatfb    : ', SUM(tmlatfb(:, :))
        WRITE(numout, *) '          TRA tmlatfn    : ', SUM(tmlatfn(:, :))
        DO jl = 1, jpltrd
          WRITE(numout, *) '          * TRA TREND INDEX jpmxl_xxx = jl = ', jl, ' tmltrd : ', SUM(tmltrd(:, :, jl))
        END DO
        WRITE(numout, *) '          TRA ztmlres (jpi/2,jpj/2) : ', ztmlres(jpi / 2, jpj / 2)
        WRITE(numout, *) '          TRA ztmlres2(jpi/2,jpj/2) : ', ztmlres2(jpi / 2, jpj / 2)
        WRITE(numout, *)
        WRITE(numout, *) '          >>>>>>>>>>>>>>>>>>  TRA SALINITY <<<<<<<<<<<<<<<<<<'
        WRITE(numout, *) '          TRA zsmlres    : ', SUM(zsmlres(:, :))
        WRITE(numout, *) '          TRA zsmltot    : ', SUM(zsmltot(:, :))
        WRITE(numout, *) '          TRA smltrdm    : ', SUM(smltrdm(:, :))
        WRITE(numout, *) '          TRA smlatfb    : ', SUM(smlatfb(:, :))
        WRITE(numout, *) '          TRA smlatfn    : ', SUM(smlatfn(:, :))
        DO jl = 1, jpltrd
          WRITE(numout, *) '          * TRA TREND INDEX jpmxl_xxx = jl = ', jl, ' smltrd : ', SUM(smltrd(:, :, jl))
        END DO
        WRITE(numout, *) '          TRA zsmlres (jpi/2,jpj/2) : ', zsmlres(jpi / 2, jpj / 2)
        WRITE(numout, *) '          TRA zsmlres2(jpi/2,jpj/2) : ', zsmlres2(jpi / 2, jpj / 2)
      END IF
      CALL profile_psy_data5 % PostEnd
    END IF MODULO_NTRD
    CALL profile_psy_data6 % PreStart('trd_mxl', 'r6', 0, 0)
    IF (ln_trdmxl_instant) THEN
      CALL iom_put("mxl_depth", hmxl(:, :))
      CALL iom_put("tml", tml(:, :))
      CALL iom_put("tml_tot", ztmltot(:, :))
      CALL iom_put("tml_res", ztmlres(:, :))
      DO jl = 1, jpltrd - 1
        CALL iom_put(TRIM("tml" // ctrd(jl, 2)), tmltrd(:, :, jl))
      END DO
      CALL iom_put(TRIM("tml" // ctrd(jpmxl_atf, 2)), ztmlatf(:, :))
      CALL iom_put("sml", sml(:, :))
      CALL iom_put("sml_tot", zsmltot(:, :))
      CALL iom_put("sml_res", zsmlres(:, :))
      DO jl = 1, jpltrd - 1
        CALL iom_put(TRIM("sml" // ctrd(jl, 2)), smltrd(:, :, jl))
      END DO
      CALL iom_put(TRIM("sml" // ctrd(jpmxl_atf, 2)), zsmlatf(:, :))
    ELSE
      CALL iom_put("mxl_depth", hmxl_sum(:, :))
      CALL iom_put("tml", tml_sum(:, :))
      CALL iom_put("tml_tot", ztmltot2(:, :))
      CALL iom_put("tml_res", ztmlres2(:, :))
      DO jl = 1, jpltrd - 1
        CALL iom_put(TRIM("tml" // ctrd(jl, 2)), ztmltrd2(:, :, jl))
      END DO
      CALL iom_put(TRIM("tml" // ctrd(jpmxl_atf, 2)), ztmlatf2(:, :))
      CALL iom_put("sml", sml_sum(:, :))
      CALL iom_put("sml_tot", zsmltot2(:, :))
      CALL iom_put("sml_res", zsmlres2(:, :))
      DO jl = 1, jpltrd - 1
        CALL iom_put(TRIM("sml" // ctrd(jl, 2)), zsmltrd2(:, :, jl))
      END DO
      CALL iom_put(TRIM("sml" // ctrd(jpmxl_atf, 2)), zsmlatf2(:, :))
    END IF
    CALL profile_psy_data6 % PostEnd
    !$ACC KERNELS
    IF (MOD(itmod, nn_trd) == 0) THEN
      nmoymltrd = 0
      tmltrdm(:, :) = 0.E0
      smltrdm(:, :) = 0.E0
      tmlatfm(:, :) = 0.E0
      smlatfm(:, :) = 0.E0
      tml_sum(:, :) = 0.E0
      sml_sum(:, :) = 0.E0
      tmltrd_csum_ln(:, :, :) = 0.E0
      smltrd_csum_ln(:, :, :) = 0.E0
      tmltrd_sum(:, :, :) = 0.E0
      smltrd_sum(:, :, :) = 0.E0
      hmxl_sum(:, :) = 0.E0
    END IF
    !$ACC END KERNELS
    CALL profile_psy_data7 % PreStart('trd_mxl', 'r7', 0, 0)
    IF (lrst_oce) CALL trd_mxl_rst_write(kt)
    CALL profile_psy_data7 % PostEnd
  END SUBROUTINE trd_mxl
  SUBROUTINE trd_mxl_init
    INTEGER :: jl
    INTEGER :: inum
    INTEGER :: ios
    REAL(KIND = wp) :: zjulian, zsto, zout
    CHARACTER(LEN = 40) :: clop
    CHARACTER(LEN = 12) :: clmxl, cltu, clsu
    NAMELIST /namtrd_mxl/ nn_trd, cn_trdrst_in, ln_trdmxl_restart, nn_ctls, cn_trdrst_out, ln_trdmxl_instant, rn_ucf, rn_rho_c
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namtrd_mxl, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namtrd_mxl in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namtrd_mxl, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namtrd_mxl in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namtrd_mxl)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) ' trd_mxl_init : Mixed-layer trends'
      WRITE(numout, *) ' ~~~~~~~~~~'
      WRITE(numout, *) '   Namelist namtrd : set trends parameters'
      WRITE(numout, *) '      frequency of trends diagnostics (glo)      nn_trd             = ', nn_trd
      WRITE(numout, *) '      density criteria used to defined the MLD   rn_rho_c           = ', rn_rho_c
      WRITE(numout, *) '      control surface type            (mld)      nn_ctls            = ', nn_ctls
      WRITE(numout, *) '      restart for ML diagnostics                 ln_trdmxl_restart  = ', ln_trdmxl_restart
      WRITE(numout, *) '      instantaneous or mean ML T/S               ln_trdmxl_instant  = ', ln_trdmxl_instant
      WRITE(numout, *) '      unit conversion factor                     rn_ucf             = ', rn_ucf
      WRITE(numout, *) '      criteria to compute the MLD                rn_rho_c           = ', rn_rho_c
    END IF
    IF (rn_rho_c /= rho_c) CALL ctl_warn('Unless you have good reason to do so, you should use the value ', 'defined in zdfmxl.F90 module to calculate the mixed layer depth')
    IF (MOD(nitend, nn_trd) /= 0) THEN
      WRITE(numout, cform_err)
      WRITE(numout, *) '                Your nitend parameter, nitend = ', nitend
      WRITE(numout, *) '                is no multiple of the trends diagnostics frequency        '
      WRITE(numout, *) '                          you defined, nn_trd   = ', nn_trd
      WRITE(numout, *) '                This will not allow you to restart from this simulation.  '
      WRITE(numout, *) '                You should reconsider this choice.                        '
      WRITE(numout, *)
      WRITE(numout, *) '                N.B. the nitend parameter is also constrained to be a     '
      WRITE(numout, *) '                     multiple of the nn_fsbc parameter '
      CALL ctl_stop('trd_mxl_init: see comment just above')
    END IF
    IF (trd_mxl_alloc() /= 0) CALL ctl_stop('STOP', 'trd_mxl_init : unable to allocate trdmxl     arrays')
    IF (trdmxl_oce_alloc() /= 0) CALL ctl_stop('STOP', 'trd_mxl_init : unable to allocate trdmxl_oce arrays')
    !$ACC KERNELS
    nkstp = nit000 - 1
    nmoymltrd = 0
    tml(:, :) = 0.E0
    sml(:, :) = 0.E0
    tmltrdm(:, :) = 0.E0
    smltrdm(:, :) = 0.E0
    tmlatfm(:, :) = 0.E0
    smlatfm(:, :) = 0.E0
    tml_sum(:, :) = 0.E0
    sml_sum(:, :) = 0.E0
    tmltrd_sum(:, :, :) = 0.E0
    smltrd_sum(:, :, :) = 0.E0
    tmltrd_csum_ln(:, :, :) = 0.E0
    smltrd_csum_ln(:, :, :) = 0.E0
    hmxl(:, :) = 0.E0
    hmxl_sum(:, :) = 0.E0
    !$ACC END KERNELS
    IF (ln_rstart .AND. ln_trdmxl_restart) THEN
      CALL trd_mxl_rst_read
    ELSE
      !$ACC KERNELS
      tmlb(:, :) = 0.E0
      smlb(:, :) = 0.E0
      tmlbb(:, :) = 0.E0
      smlbb(:, :) = 0.E0
      tmlbn(:, :) = 0.E0
      smlbn(:, :) = 0.E0
      tml_sumb(:, :) = 0.E0
      sml_sumb(:, :) = 0.E0
      tmltrd_csum_ub(:, :, :) = 0.E0
      smltrd_csum_ub(:, :, :) = 0.E0
      tmltrd_atf_sumb(:, :) = 0.E0
      smltrd_atf_sumb(:, :) = 0.E0
      !$ACC END KERNELS
    END IF
    icount = 1
    ionce = 1
    IF (nn_ctls == 1) THEN
      CALL ctl_opn(inum, 'ctlsurf_idx', 'OLD', 'UNFORMATTED', 'SEQUENTIAL', - 1, numout, lwp)
      READ(inum, *) nbol
      CLOSE(UNIT = inum)
    END IF
    IF (nn_ctls == 0) THEN
      clmxl = 'Mixed Layer '
    ELSE IF (nn_ctls == 1) THEN
      clmxl = '      Bowl '
    ELSE IF (nn_ctls >= 2) THEN
      WRITE(clmxl, '(A10,I2,1X)') 'Levels 1 -', nn_ctls
    END IF
    ctrd(jpmxl_xad, 1) = " Zonal advection"
    ctrd(jpmxl_xad, 2) = "_xad"
    ctrd(jpmxl_yad, 1) = " Meridional advection"
    ctrd(jpmxl_yad, 2) = "_yad"
    ctrd(jpmxl_zad, 1) = " Vertical advection"
    ctrd(jpmxl_zad, 2) = "_zad"
    ctrd(jpmxl_ldf, 1) = " Lateral diffusion"
    ctrd(jpmxl_ldf, 2) = "_ldf"
    ctrd(jpmxl_for, 1) = " Forcing"
    ctrd(jpmxl_for, 2) = "_for"
    ctrd(jpmxl_zdf, 1) = " Vertical diff. (Kz)"
    ctrd(jpmxl_zdf, 2) = "_zdf"
    ctrd(jpmxl_bbc, 1) = " Geothermal flux"
    ctrd(jpmxl_bbc, 2) = "_bbc"
    ctrd(jpmxl_bbl, 1) = " Adv/diff. Bottom boundary layer"
    ctrd(jpmxl_bbl, 2) = "_bbl"
    ctrd(jpmxl_dmp, 1) = " Tracer damping"
    ctrd(jpmxl_dmp, 2) = "_dmp"
    ctrd(jpmxl_npc, 1) = " Non penetrative convec. adjust."
    ctrd(jpmxl_npc, 2) = "_npc"
    ctrd(jpmxl_atf, 1) = " Asselin time filter"
    ctrd(jpmxl_atf, 2) = "_atf"
    IF (rn_ucf == 1.) THEN
      cltu = "degC/s"
      clsu = "p.s.u./s"
    ELSE IF (rn_ucf == 3600. * 24.) THEN
      cltu = "degC/day"
      clsu = "p.s.u./day"
    ELSE
      cltu = "unknown?"
      clsu = "unknown?"
    END IF
  END SUBROUTINE trd_mxl_init
END MODULE trdmxl