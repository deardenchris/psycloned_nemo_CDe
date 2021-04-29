MODULE traadv
  USE oce
  USE dom_oce
  USE domvvl
  USE sbcwave
  USE sbc_oce
  USE traadv_cen
  USE traadv_fct
  USE traadv_mus
  USE traadv_ubs
  USE traadv_qck
  USE tramle
  USE ldftra
  USE ldfslp
  USE trd_oce
  USE trdtra
  USE diaptr
  USE in_out_manager
  USE iom
  USE prtctl
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_adv
  PUBLIC :: tra_adv_init
  LOGICAL :: ln_traadv_OFF
  LOGICAL :: ln_traadv_cen
  INTEGER :: nn_cen_h, nn_cen_v
  LOGICAL :: ln_traadv_fct
  INTEGER :: nn_fct_h, nn_fct_v
  LOGICAL :: ln_traadv_mus
  LOGICAL :: ln_mus_ups
  LOGICAL :: ln_traadv_ubs
  INTEGER :: nn_ubs_v
  LOGICAL :: ln_traadv_qck
  INTEGER :: nadv
  INTEGER, PARAMETER :: np_NO_adv = 0
  INTEGER, PARAMETER :: np_CEN = 1
  INTEGER, PARAMETER :: np_FCT = 2
  INTEGER, PARAMETER :: np_MUS = 3
  INTEGER, PARAMETER :: np_UBS = 4
  INTEGER, PARAMETER :: np_QCK = 5
  CONTAINS
  SUBROUTINE tra_adv(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: jk
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zun, zvn, zwn
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ztrdt, ztrds
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    CALL profile_psy_data0 % PreStart('tra_adv', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('tra_adv')
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    IF (neuler == 0 .AND. kt == nit000) THEN
      r2dt = rdt
    ELSE IF (kt <= nit000 + 1) THEN
      r2dt = 2._wp * rdt
    END IF
    zun(:, :, jpk) = 0._wp
    zvn(:, :, jpk) = 0._wp
    zwn(:, :, jpk) = 0._wp
    IF (ln_wave .AND. ln_sdw) THEN
      DO jk = 1, jpkm1
        zun(:, :, jk) = e2u(:, :) * e3u_n(:, :, jk) * (un(:, :, jk) + usd(:, :, jk))
        zvn(:, :, jk) = e1v(:, :) * e3v_n(:, :, jk) * (vn(:, :, jk) + vsd(:, :, jk))
        zwn(:, :, jk) = e1e2t(:, :) * (wn(:, :, jk) + wsd(:, :, jk))
      END DO
    ELSE
      DO jk = 1, jpkm1
        zun(:, :, jk) = e2u(:, :) * e3u_n(:, :, jk) * un(:, :, jk)
        zvn(:, :, jk) = e1v(:, :) * e3v_n(:, :, jk) * vn(:, :, jk)
        zwn(:, :, jk) = e1e2t(:, :) * wn(:, :, jk)
      END DO
    END IF
    IF (ln_vvl_ztilde .OR. ln_vvl_layer) THEN
      zun(:, :, :) = zun(:, :, :) + un_td(:, :, :)
      zvn(:, :, :) = zvn(:, :, :) + vn_td(:, :, :)
    END IF
    zun(:, :, jpk) = 0._wp
    zvn(:, :, jpk) = 0._wp
    zwn(:, :, jpk) = 0._wp
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('tra_adv', 'r1', 0, 0)
    IF (ln_ldfeiv .AND. .NOT. ln_traldf_triad) CALL ldf_eiv_trp(kt, nit000, zun, zvn, zwn, 'TRA')
    IF (ln_mle) CALL tra_mle_trp(kt, nit000, zun, zvn, zwn, 'TRA')
    CALL iom_put("uocetr_eff", zun)
    CALL iom_put("vocetr_eff", zvn)
    CALL iom_put("wocetr_eff", zwn)
    IF (ln_diaptr) CALL dia_ptr(zvn)
    CALL profile_psy_data1 % PostEnd
    IF (l_trdtra) THEN
      ALLOCATE(ztrdt(jpi, jpj, jpk), ztrds(jpi, jpj, jpk))
      !$ACC KERNELS
      ztrdt(:, :, :) = tsa(:, :, :, jp_tem)
      ztrds(:, :, :) = tsa(:, :, :, jp_sal)
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data2 % PreStart('tra_adv', 'r2', 0, 0)
    SELECT CASE (nadv)
    CASE (np_CEN)
      CALL tra_adv_cen(kt, nit000, 'TRA', zun, zvn, zwn, tsn, tsa, jpts, nn_cen_h, nn_cen_v)
    CASE (np_FCT)
      CALL tra_adv_fct(kt, nit000, 'TRA', r2dt, zun, zvn, zwn, tsb, tsn, tsa, jpts, nn_fct_h, nn_fct_v)
    CASE (np_MUS)
      CALL tra_adv_mus(kt, nit000, 'TRA', r2dt, zun, zvn, zwn, tsb, tsa, jpts, ln_mus_ups)
    CASE (np_UBS)
      CALL tra_adv_ubs(kt, nit000, 'TRA', r2dt, zun, zvn, zwn, tsb, tsn, tsa, jpts, nn_ubs_v)
    CASE (np_QCK)
      CALL tra_adv_qck(kt, nit000, 'TRA', r2dt, zun, zvn, zwn, tsb, tsn, tsa, jpts)
    END SELECT
    CALL profile_psy_data2 % PostEnd
    IF (l_trdtra) THEN
      !$ACC KERNELS
      DO jk = 1, jpkm1
        ztrdt(:, :, jk) = tsa(:, :, jk, jp_tem) - ztrdt(:, :, jk)
        ztrds(:, :, jk) = tsa(:, :, jk, jp_sal) - ztrds(:, :, jk)
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data3 % PreStart('tra_adv', 'r3', 0, 0)
      CALL trd_tra(kt, 'TRA', jp_tem, jptra_totad, ztrdt)
      CALL trd_tra(kt, 'TRA', jp_sal, jptra_totad, ztrds)
      DEALLOCATE(ztrdt, ztrds)
      CALL profile_psy_data3 % PostEnd
    END IF
    CALL profile_psy_data4 % PreStart('tra_adv', 'r4', 0, 0)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = tsa(:, :, :, jp_tem), clinfo1 = ' adv  - Ta: ', mask1 = tmask, tab3d_2 = tsa(:, :, :, jp_sal), clinfo2 = ' Sa: ', mask2 = tmask, clinfo3 = 'tra')
    IF (ln_timing) CALL timing_stop('tra_adv')
    CALL profile_psy_data4 % PostEnd
  END SUBROUTINE tra_adv
  SUBROUTINE tra_adv_init
    INTEGER :: ioptio, ios
    NAMELIST /namtra_adv/ ln_traadv_OFF, ln_traadv_cen, nn_cen_h, nn_cen_v, ln_traadv_fct, nn_fct_h, nn_fct_v, ln_traadv_mus, ln_mus_ups, ln_traadv_ubs, nn_ubs_v, ln_traadv_qck
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namtra_adv, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namtra_adv in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namtra_adv, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namtra_adv in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namtra_adv)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'tra_adv_init : choice/control of the tracer advection scheme'
      WRITE(numout, *) '~~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist namtra_adv : chose a advection scheme for tracers'
      WRITE(numout, *) '      No advection on T & S                     ln_traadv_OFF = ', ln_traadv_OFF
      WRITE(numout, *) '      centered scheme                           ln_traadv_cen = ', ln_traadv_cen
      WRITE(numout, *) '            horizontal 2nd/4th order               nn_cen_h   = ', nn_fct_h
      WRITE(numout, *) '            vertical   2nd/4th order               nn_cen_v   = ', nn_fct_v
      WRITE(numout, *) '      Flux Corrected Transport scheme           ln_traadv_fct = ', ln_traadv_fct
      WRITE(numout, *) '            horizontal 2nd/4th order               nn_fct_h   = ', nn_fct_h
      WRITE(numout, *) '            vertical   2nd/4th order               nn_fct_v   = ', nn_fct_v
      WRITE(numout, *) '      MUSCL scheme                              ln_traadv_mus = ', ln_traadv_mus
      WRITE(numout, *) '            + upstream scheme near river mouths    ln_mus_ups = ', ln_mus_ups
      WRITE(numout, *) '      UBS scheme                                ln_traadv_ubs = ', ln_traadv_ubs
      WRITE(numout, *) '            vertical   2nd/4th order               nn_ubs_v   = ', nn_ubs_v
      WRITE(numout, *) '      QUICKEST scheme                           ln_traadv_qck = ', ln_traadv_qck
    END IF
    ioptio = 0
    IF (ln_traadv_off) THEN
      ioptio = ioptio + 1
      nadv = np_no_adv
    END IF
    IF (ln_traadv_cen) THEN
      ioptio = ioptio + 1
      nadv = np_cen
    END IF
    IF (ln_traadv_fct) THEN
      ioptio = ioptio + 1
      nadv = np_fct
    END IF
    IF (ln_traadv_mus) THEN
      ioptio = ioptio + 1
      nadv = np_mus
    END IF
    IF (ln_traadv_ubs) THEN
      ioptio = ioptio + 1
      nadv = np_ubs
    END IF
    IF (ln_traadv_qck) THEN
      ioptio = ioptio + 1
      nadv = np_qck
    END IF
    IF (ioptio /= 1) CALL ctl_stop('tra_adv_init: Choose ONE advection option in namelist namtra_adv')
    IF (ln_traadv_cen .AND. (nn_cen_h /= 2 .AND. nn_cen_h /= 4) .AND. (nn_cen_v /= 2 .AND. nn_cen_v /= 4)) THEN
      CALL ctl_stop('tra_adv_init: CEN scheme, choose 2nd or 4th order')
    END IF
    IF (ln_traadv_fct .AND. (nn_fct_h /= 2 .AND. nn_fct_h /= 4) .AND. (nn_fct_v /= 2 .AND. nn_fct_v /= 4)) THEN
      CALL ctl_stop('tra_adv_init: FCT scheme, choose 2nd or 4th order')
    END IF
    IF (ln_traadv_ubs .AND. (nn_ubs_v /= 2 .AND. nn_ubs_v /= 4)) THEN
      CALL ctl_stop('tra_adv_init: UBS scheme, choose 2nd or 4th order')
    END IF
    IF (ln_traadv_ubs .AND. nn_ubs_v == 4) THEN
      CALL ctl_warn('tra_adv_init: UBS scheme, only 2nd FCT scheme available on the vertical. It will be used')
    END IF
    IF (ln_isfcav) THEN
      IF (ln_traadv_cen .AND. nn_cen_v == 4 .OR. ln_traadv_fct .AND. nn_fct_v == 4) CALL ctl_stop('tra_adv_init: 4th order COMPACT scheme not allowed with ISF')
    END IF
    IF (lwp) THEN
      WRITE(numout, *)
      SELECT CASE (nadv)
      CASE (np_no_adv)
        WRITE(numout, *) '   ==>>>   NO T-S advection'
      CASE (np_cen)
        WRITE(numout, *) '   ==>>>   CEN      scheme is used. Horizontal order: ', nn_cen_h, ' Vertical   order: ', nn_cen_v
      CASE (np_fct)
        WRITE(numout, *) '   ==>>>   FCT      scheme is used. Horizontal order: ', nn_fct_h, ' Vertical   order: ', nn_fct_v
      CASE (np_mus)
        WRITE(numout, *) '   ==>>>   MUSCL    scheme is used'
      CASE (np_ubs)
        WRITE(numout, *) '   ==>>>   UBS      scheme is used'
      CASE (np_qck)
        WRITE(numout, *) '   ==>>>   QUICKEST scheme is used'
      END SELECT
    END IF
    CALL tra_mle_init
  END SUBROUTINE tra_adv_init
END MODULE traadv