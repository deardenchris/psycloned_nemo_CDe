MODULE step
  USE step_oce
  USE iom
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: stp
  CONTAINS
  SUBROUTINE stp(kstp)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kstp
    INTEGER :: ji, jj, jk
    INTEGER :: indic
    INTEGER :: kcall
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('stp', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('stp')
    indic = 0
    IF (kstp == nit000) THEN
      CALL iom_init(cxios_context)
      IF (ln_crs) CALL iom_init(TRIM(cxios_context) // "_crs")
    END IF
    IF (kstp /= nit000) CALL day(kstp)
    CALL iom_setkt(kstp - nit000 + 1, cxios_context)
    IF (ln_crs) CALL iom_setkt(kstp - nit000 + 1, TRIM(cxios_context) // "_crs")
    IF (ln_tide) CALL sbc_tide(kstp)
    IF (ln_apr_dyn) CALL sbc_apr(kstp)
    IF (ln_bdy) CALL bdy_dta(kstp, time_offset = + 1)
    CALL sbc(kstp)
    IF (ln_sto_eos) CALL sto_par(kstp)
    IF (ln_sto_eos) CALL sto_pts(tsn)
    CALL eos_rab(tsb, rab_b)
    CALL eos_rab(tsn, rab_n)
    CALL bn2(tsb, rab_b, rn2b)
    CALL bn2(tsn, rab_n, rn2)
    CALL zdf_phy(kstp)
    IF (l_ldfslp) THEN
      CALL eos(tsb, rhd, gdept_0(:, :, :))
      IF (ln_zps .AND. .NOT. ln_isfcav) CALL zps_hde(kstp, jpts, tsb, gtsu, gtsv, rhd, gru, grv)
      IF (ln_zps .AND. ln_isfcav) CALL zps_hde_isf(kstp, jpts, tsb, gtsu, gtsv, gtui, gtvi, rhd, gru, grv, grui, grvi)
      IF (ln_traldf_triad) THEN
        CALL ldf_slp_triad(kstp)
      ELSE
        CALL ldf_slp(kstp, rhd, rn2b)
      END IF
    END IF
    IF (l_ldftra_time .OR. l_ldfeiv_time) CALL ldf_tra(kstp)
    IF (l_ldfdyn_time) CALL ldf_dyn(kstp)
    CALL ssh_nxt(kstp)
    IF (.NOT. ln_linssh) CALL dom_vvl_sf_nxt(kstp)
    CALL wzv(kstp)
    IF (ln_zad_Aimp) CALL wAimp(kstp)
    CALL eos(tsn, rhd, rhop, gdept_n(:, :, :))
    IF (ln_zps .AND. .NOT. ln_isfcav) CALL zps_hde(kstp, jpts, tsn, gtsu, gtsv, rhd, gru, grv)
    IF (ln_zps .AND. ln_isfcav) CALL zps_hde_isf(kstp, jpts, tsn, gtsu, gtsv, gtui, gtvi, rhd, gru, grv, grui, grvi)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    ua(:, :, :) = 0._wp
    va(:, :, :) = 0._wp
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('stp', 'r1', 0, 0)
    IF (lk_asminc .AND. ln_asmiau .AND. ln_dyninc) CALL dyn_asm_inc(kstp)
    IF (ln_bdy) CALL bdy_dyn3d_dmp(kstp)
    CALL dyn_adv(kstp)
    CALL dyn_vor(kstp)
    CALL dyn_ldf(kstp)
    IF (ln_zdfosm) CALL dyn_osm(kstp)
    CALL dyn_hpg(kstp)
    CALL dyn_spg(kstp)
    IF (ln_dynspg_ts) THEN
      CALL div_hor(kstp)
      IF (.NOT. ln_linssh) CALL dom_vvl_sf_nxt(kstp, kcall = 2)
      CALL wzv(kstp)
      IF (ln_zad_Aimp) CALL wAimp(kstp)
    END IF
    CALL dyn_zdf(kstp)
    IF (ln_diurnal) CALL stp_diurnal(kstp)
    IF (lk_floats) CALL flo_stp(kstp)
    IF (ln_diacfl) CALL dia_cfl(kstp)
    IF (lk_diahth) CALL dia_hth(kstp)
    IF (lk_diadct) CALL dia_dct(kstp)
    CALL dia_ar5(kstp)
    IF (lk_diaharm) CALL dia_harm(kstp)
    CALL dia_prod(kstp)
    CALL dia_wri(kstp)
    IF (ln_crs) CALL crs_fld(kstp)
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    tsa(:, :, :, :) = 0._wp
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('stp', 'r2', 0, 0)
    IF (lk_asminc .AND. ln_asmiau .AND. ln_trainc) CALL tra_asm_inc(kstp)
    CALL tra_sbc(kstp)
    IF (ln_traqsr) CALL tra_qsr(kstp)
    IF (ln_trabbc) CALL tra_bbc(kstp)
    IF (ln_trabbl) CALL tra_bbl(kstp)
    IF (ln_tradmp) CALL tra_dmp(kstp)
    IF (ln_bdy) CALL bdy_tra_dmp(kstp)
    CALL tra_adv(kstp)
    IF (ln_zdfosm) CALL tra_osm(kstp)
    IF (lrst_oce .AND. ln_zdfosm) CALL osm_rst(kstp, 'WRITE')
    CALL tra_ldf(kstp)
    IF (ln_diaptr) CALL dia_ptr
    CALL tra_zdf(kstp)
    IF (ln_zdfnpc) CALL tra_npc(kstp)
    CALL tra_nxt(kstp)
    CALL dyn_nxt(kstp)
    CALL ssh_swp(kstp)
    IF (.NOT. ln_linssh) CALL dom_vvl_sf_swp(kstp)
    IF (ln_diahsb) CALL dia_hsb(kstp)
    IF (lrst_oce) CALL rst_write(kstp)
    IF (ln_sto_eos) CALL sto_rst_write(kstp)
    IF (ln_diaobs) CALL dia_obs(kstp)
    CALL stp_ctl(kstp, indic)
    IF (kstp == nit000) THEN
      CALL iom_close(numror)
      IF (lwm) CALL FLUSH(numond)
      IF (lwm .AND. numoni /= - 1) CALL FLUSH(numoni)
    END IF
    IF (lk_oasis) CALL sbc_cpl_snd(kstp)
    IF (ln_timing) CALL timing_stop('stp')
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE stp
END MODULE step