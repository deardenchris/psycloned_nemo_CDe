MODULE icbdia
  USE par_oce
  USE dom_oce
  USE in_out_manager
  USE lib_mpp
  USE iom
  USE icb_oce
  USE icbutl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: icb_dia_init
  PUBLIC :: icb_dia
  PUBLIC :: icb_dia_step
  PUBLIC :: icb_dia_put
  PUBLIC :: icb_dia_melt
  PUBLIC :: icb_dia_size
  PUBLIC :: icb_dia_speed
  PUBLIC :: icb_dia_calve
  PUBLIC :: icb_dia_income
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE, PUBLIC :: berg_melt
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE, PUBLIC :: berg_melt_hcflx
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE, PUBLIC :: berg_melt_qlat
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE, PUBLIC :: buoy_melt
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE, PUBLIC :: eros_melt
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE, PUBLIC :: conv_melt
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE, PUBLIC :: bits_src
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE, PUBLIC :: bits_melt
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE, PUBLIC :: bits_mass
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE, PUBLIC :: virtual_area
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE, PUBLIC :: berg_mass
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE, PUBLIC :: real_calving
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: tmpc
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: rsumbuf
  INTEGER, DIMENSION(:), ALLOCATABLE :: nsumbuf
  REAL(KIND = wp) :: berg_melt_net
  REAL(KIND = wp) :: bits_src_net
  REAL(KIND = wp) :: bits_melt_net
  REAL(KIND = wp) :: bits_mass_start, bits_mass_end
  REAL(KIND = wp) :: floating_heat_start, floating_heat_end
  REAL(KIND = wp) :: floating_mass_start, floating_mass_end
  REAL(KIND = wp) :: bergs_mass_start, bergs_mass_end
  REAL(KIND = wp) :: stored_start, stored_heat_start
  REAL(KIND = wp) :: stored_end, stored_heat_end
  REAL(KIND = wp) :: calving_src_net, calving_out_net
  REAL(KIND = wp) :: calving_src_heat_net, calving_out_heat_net
  REAL(KIND = wp) :: calving_src_heat_used_net
  REAL(KIND = wp) :: calving_rcv_net, calving_ret_net, calving_used_net
  REAL(KIND = wp) :: heat_to_bergs_net, heat_to_ocean_net, melt_net
  REAL(KIND = wp) :: calving_to_bergs_net
  INTEGER :: nbergs_start, nbergs_end, nbergs_calved
  INTEGER :: nbergs_melted
  INTEGER :: nspeeding_tickets
  INTEGER, DIMENSION(nclasses) :: nbergs_calved_by_class
  CONTAINS
  SUBROUTINE icb_dia_init
    IF (.NOT. ln_bergdia) RETURN
    ALLOCATE(berg_melt(jpi, jpj))
    berg_melt(:, :) = 0._wp
    ALLOCATE(berg_melt_hcflx(jpi, jpj))
    berg_melt_hcflx(:, :) = 0._wp
    ALLOCATE(berg_melt_qlat(jpi, jpj))
    berg_melt_qlat(:, :) = 0._wp
    ALLOCATE(buoy_melt(jpi, jpj))
    buoy_melt(:, :) = 0._wp
    ALLOCATE(eros_melt(jpi, jpj))
    eros_melt(:, :) = 0._wp
    ALLOCATE(conv_melt(jpi, jpj))
    conv_melt(:, :) = 0._wp
    ALLOCATE(bits_src(jpi, jpj))
    bits_src(:, :) = 0._wp
    ALLOCATE(bits_melt(jpi, jpj))
    bits_melt(:, :) = 0._wp
    ALLOCATE(bits_mass(jpi, jpj))
    bits_mass(:, :) = 0._wp
    ALLOCATE(virtual_area(jpi, jpj))
    virtual_area(:, :) = 0._wp
    ALLOCATE(berg_mass(jpi, jpj))
    berg_mass(:, :) = 0._wp
    ALLOCATE(real_calving(jpi, jpj, nclasses))
    real_calving(:, :, :) = 0._wp
    ALLOCATE(tmpc(jpi, jpj))
    tmpc(:, :) = 0._wp
    nbergs_start = 0
    nbergs_end = 0
    stored_end = 0._wp
    nbergs_start = 0._wp
    stored_start = 0._wp
    nbergs_melted = 0
    nbergs_calved = 0
    nbergs_calved_by_class(:) = 0
    nspeeding_tickets = 0
    stored_heat_end = 0._wp
    floating_heat_end = 0._wp
    floating_mass_end = 0._wp
    bergs_mass_end = 0._wp
    bits_mass_end = 0._wp
    stored_heat_start = 0._wp
    floating_heat_start = 0._wp
    floating_mass_start = 0._wp
    bergs_mass_start = 0._wp
    bits_mass_start = 0._wp
    bits_mass_end = 0._wp
    calving_used_net = 0._wp
    calving_to_bergs_net = 0._wp
    heat_to_bergs_net = 0._wp
    heat_to_ocean_net = 0._wp
    calving_rcv_net = 0._wp
    calving_ret_net = 0._wp
    calving_src_net = 0._wp
    calving_out_net = 0._wp
    calving_src_heat_net = 0._wp
    calving_src_heat_used_net = 0._wp
    calving_out_heat_net = 0._wp
    melt_net = 0._wp
    berg_melt_net = 0._wp
    bits_melt_net = 0._wp
    bits_src_net = 0._wp
    floating_mass_start = icb_utl_mass(first_berg)
    bergs_mass_start = icb_utl_mass(first_berg, justbergs = .TRUE.)
    bits_mass_start = icb_utl_mass(first_berg, justbits = .TRUE.)
    IF (lk_mpp) THEN
      ALLOCATE(rsumbuf(23))
      rsumbuf(:) = 0._wp
      ALLOCATE(nsumbuf(4 + nclasses))
      nsumbuf(:) = 0
      rsumbuf(1) = floating_mass_start
      rsumbuf(2) = bergs_mass_start
      rsumbuf(3) = bits_mass_start
      CALL mpp_sum('icbdia', rsumbuf(1 : 3), 3)
      floating_mass_start = rsumbuf(1)
      bergs_mass_start = rsumbuf(2)
      bits_mass_start = rsumbuf(3)
    END IF
  END SUBROUTINE icb_dia_init
  SUBROUTINE icb_dia(ld_budge)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    LOGICAL, INTENT(IN) :: ld_budge
    INTEGER :: ik
    REAL(KIND = wp) :: zunused_calving, ztmpsum, zgrdd_berg_mass, zgrdd_bits_mass
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    IF (.NOT. ln_bergdia) RETURN
    CALL profile_psy_data0 % PreStart('icb_dia', 'r0', 0, 0)
    zunused_calving = SUM(berg_grid % calving(:, :))
    ztmpsum = SUM(berg_grid % floating_melt(:, :) * e1e2t(:, :) * tmask_i(:, :))
    melt_net = melt_net + ztmpsum * berg_dt
    calving_out_net = calving_out_net + (zunused_calving + ztmpsum) * berg_dt
    ztmpsum = SUM(berg_melt(:, :) * e1e2t(:, :) * tmask_i(:, :))
    berg_melt_net = berg_melt_net + ztmpsum * berg_dt
    ztmpsum = SUM(bits_src(:, :) * e1e2t(:, :) * tmask_i(:, :))
    bits_src_net = bits_src_net + ztmpsum * berg_dt
    ztmpsum = SUM(bits_melt(:, :) * e1e2t(:, :) * tmask_i(:, :))
    bits_melt_net = bits_melt_net + ztmpsum * berg_dt
    ztmpsum = SUM(src_calving(:, :) * tmask_i(:, :))
    calving_ret_net = calving_ret_net + ztmpsum * berg_dt
    ztmpsum = SUM(berg_grid % calving_hflx(:, :) * e1e2t(:, :) * tmask_i(:, :))
    calving_out_heat_net = calving_out_heat_net + ztmpsum * berg_dt
    IF (ld_budge) THEN
      stored_end = SUM(berg_grid % stored_ice(:, :, :))
      stored_heat_end = SUM(berg_grid % stored_heat(:, :))
      floating_mass_end = icb_utl_mass(first_berg)
      bergs_mass_end = icb_utl_mass(first_berg, justbergs = .TRUE.)
      bits_mass_end = icb_utl_mass(first_berg, justbits = .TRUE.)
      floating_heat_end = icb_utl_heat(first_berg)
      nbergs_end = icb_utl_count()
      zgrdd_berg_mass = SUM(berg_mass(:, :) * e1e2t(:, :) * tmask_i(:, :))
      zgrdd_bits_mass = SUM(bits_mass(:, :) * e1e2t(:, :) * tmask_i(:, :))
      IF (lk_mpp) THEN
        rsumbuf(1) = stored_end
        rsumbuf(2) = stored_heat_end
        rsumbuf(3) = floating_mass_end
        rsumbuf(4) = bergs_mass_end
        rsumbuf(5) = bits_mass_end
        rsumbuf(6) = floating_heat_end
        rsumbuf(7) = calving_ret_net
        rsumbuf(8) = calving_out_net
        rsumbuf(9) = calving_rcv_net
        rsumbuf(10) = calving_src_net
        rsumbuf(11) = calving_src_heat_net
        rsumbuf(12) = calving_src_heat_used_net
        rsumbuf(13) = calving_out_heat_net
        rsumbuf(14) = calving_used_net
        rsumbuf(15) = calving_to_bergs_net
        rsumbuf(16) = heat_to_bergs_net
        rsumbuf(17) = heat_to_ocean_net
        rsumbuf(18) = melt_net
        rsumbuf(19) = berg_melt_net
        rsumbuf(20) = bits_src_net
        rsumbuf(21) = bits_melt_net
        rsumbuf(22) = zgrdd_berg_mass
        rsumbuf(23) = zgrdd_bits_mass
        CALL mpp_sum('icbdia', rsumbuf(1 : 23), 23)
        stored_end = rsumbuf(1)
        stored_heat_end = rsumbuf(2)
        floating_mass_end = rsumbuf(3)
        bergs_mass_end = rsumbuf(4)
        bits_mass_end = rsumbuf(5)
        floating_heat_end = rsumbuf(6)
        calving_ret_net = rsumbuf(7)
        calving_out_net = rsumbuf(8)
        calving_rcv_net = rsumbuf(9)
        calving_src_net = rsumbuf(10)
        calving_src_heat_net = rsumbuf(11)
        calving_src_heat_used_net = rsumbuf(12)
        calving_out_heat_net = rsumbuf(13)
        calving_used_net = rsumbuf(14)
        calving_to_bergs_net = rsumbuf(15)
        heat_to_bergs_net = rsumbuf(16)
        heat_to_ocean_net = rsumbuf(17)
        melt_net = rsumbuf(18)
        berg_melt_net = rsumbuf(19)
        bits_src_net = rsumbuf(20)
        bits_melt_net = rsumbuf(21)
        zgrdd_berg_mass = rsumbuf(22)
        zgrdd_bits_mass = rsumbuf(23)
        nsumbuf(1) = nbergs_end
        nsumbuf(2) = nbergs_calved
        nsumbuf(3) = nbergs_melted
        nsumbuf(4) = nspeeding_tickets
        DO ik = 1, nclasses
          nsumbuf(4 + ik) = nbergs_calved_by_class(ik)
        END DO
        CALL mpp_sum('icbdia', nsumbuf(1 : nclasses + 4), nclasses + 4)
        nbergs_end = nsumbuf(1)
        nbergs_calved = nsumbuf(2)
        nbergs_melted = nsumbuf(3)
        nspeeding_tickets = nsumbuf(4)
        DO ik = 1, nclasses
          nbergs_calved_by_class(ik) = nsumbuf(4 + ik)
        END DO
      END IF
      CALL report_state('stored ice', 'kg', '', stored_start, '', stored_end, '')
      CALL report_state('floating', 'kg', '', floating_mass_start, '', floating_mass_end, '', nbergs_end)
      CALL report_state('icebergs', 'kg', '', bergs_mass_start, '', bergs_mass_end, '')
      CALL report_state('bits', 'kg', '', bits_mass_start, '', bits_mass_end, '')
      CALL report_istate('berg #', '', nbergs_start, '', nbergs_end, '')
      CALL report_ibudget('berg #', 'calved', nbergs_calved, 'melted', nbergs_melted, '#', nbergs_start, nbergs_end)
      CALL report_budget('stored mass', 'kg', 'calving used', calving_used_net, 'bergs', calving_to_bergs_net, 'stored mass', &
&stored_start, stored_end)
      CALL report_budget('floating mass', 'kg', 'calving used', calving_to_bergs_net, 'bergs', melt_net, 'stored mass', &
&floating_mass_start, floating_mass_end)
      CALL report_budget('berg mass', 'kg', 'calving', calving_to_bergs_net, 'melt+eros', berg_melt_net, 'berg mass', &
&bergs_mass_start, bergs_mass_end)
      CALL report_budget('bits mass', 'kg', 'eros used', bits_src_net, 'bergs', bits_melt_net, 'stored mass', bits_mass_start, &
&bits_mass_end)
      CALL report_budget('net mass', 'kg', 'recvd', calving_rcv_net, 'rtrnd', calving_ret_net, 'net mass', &
&stored_start + floating_mass_start, stored_end + floating_mass_end)
      CALL report_consistant('iceberg mass', 'kg', 'gridded', zgrdd_berg_mass, 'bergs', bergs_mass_end)
      CALL report_consistant('bits mass', 'kg', 'gridded', zgrdd_bits_mass, 'bits', bits_mass_end)
      CALL report_state('net heat', 'J', '', stored_heat_start + floating_heat_start, '', stored_heat_end + floating_heat_end, '')
      CALL report_state('stored heat', 'J', '', stored_heat_start, '', stored_heat_end, '')
      CALL report_state('floating heat', 'J', '', floating_heat_start, '', floating_heat_end, '')
      CALL report_budget('net heat', 'J', 'net heat', calving_src_heat_net, 'net heat', calving_out_heat_net, 'net heat', &
&stored_heat_start + floating_heat_start, stored_heat_end + floating_heat_end)
      CALL report_budget('stored heat', 'J', 'calving used', calving_src_heat_used_net, 'bergs', heat_to_bergs_net, 'net heat', &
&stored_heat_start, stored_heat_end)
      CALL report_budget('flting heat', 'J', 'calved', heat_to_bergs_net, 'melt', heat_to_ocean_net, 'net heat', &
&floating_heat_start, floating_heat_end)
      IF (nn_verbose_level >= 1) THEN
        CALL report_consistant('top interface', 'kg', 'from SIS', calving_src_net, 'received', calving_rcv_net)
        CALL report_consistant('bot interface', 'kg', 'sent', calving_out_net, 'returned', calving_ret_net)
      END IF
      IF (nn_verbose_level > 0) THEN
        WRITE(numicb, '("calved by class = ",i6,20(",",i6))') (nbergs_calved_by_class(ik), ik = 1, nclasses)
        IF (nspeeding_tickets > 0) WRITE(numicb, '("speeding tickets issued = ",i6)') nspeeding_tickets
      END IF
      nbergs_start = nbergs_end
      stored_start = stored_end
      nbergs_melted = 0
      nbergs_calved = 0
      nbergs_calved_by_class(:) = 0
      nspeeding_tickets = 0
      stored_heat_start = stored_heat_end
      floating_heat_start = floating_heat_end
      floating_mass_start = floating_mass_end
      bergs_mass_start = bergs_mass_end
      bits_mass_start = bits_mass_end
      calving_used_net = 0._wp
      calving_to_bergs_net = 0._wp
      heat_to_bergs_net = 0._wp
      heat_to_ocean_net = 0._wp
      calving_rcv_net = 0._wp
      calving_ret_net = 0._wp
      calving_src_net = 0._wp
      calving_out_net = 0._wp
      calving_src_heat_net = 0._wp
      calving_src_heat_used_net = 0._wp
      calving_out_heat_net = 0._wp
      melt_net = 0._wp
      berg_melt_net = 0._wp
      bits_melt_net = 0._wp
      bits_src_net = 0._wp
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_dia
  SUBROUTINE icb_dia_step
    IF (.NOT. ln_bergdia) RETURN
    berg_melt(:, :) = 0._wp
    berg_melt_hcflx(:, :) = 0._wp
    berg_melt_qlat(:, :) = 0._wp
    buoy_melt(:, :) = 0._wp
    eros_melt(:, :) = 0._wp
    conv_melt(:, :) = 0._wp
    bits_src(:, :) = 0._wp
    bits_melt(:, :) = 0._wp
    bits_mass(:, :) = 0._wp
    berg_mass(:, :) = 0._wp
    virtual_area(:, :) = 0._wp
    real_calving(:, :, :) = 0._wp
  END SUBROUTINE icb_dia_step
  SUBROUTINE icb_dia_put
    IF (.NOT. ln_bergdia) RETURN
    CALL iom_put("berg_melt", berg_melt(:, :))
    CALL iom_put("berg_melt_hcflx", berg_melt_hcflx(:, :))
    CALL iom_put("berg_melt_qlat", berg_melt_qlat(:, :))
    CALL iom_put("berg_buoy_melt", buoy_melt(:, :))
    CALL iom_put("berg_eros_melt", eros_melt(:, :))
    CALL iom_put("berg_conv_melt", conv_melt(:, :))
    CALL iom_put("berg_virtual_area", virtual_area(:, :))
    CALL iom_put("bits_src", bits_src(:, :))
    CALL iom_put("bits_melt", bits_melt(:, :))
    CALL iom_put("bits_mass", bits_mass(:, :))
    CALL iom_put("berg_mass", berg_mass(:, :))
    CALL iom_put("berg_real_calving", real_calving(:, :, :))
  END SUBROUTINE icb_dia_put
  SUBROUTINE icb_dia_calve(ki, kj, kn, pcalved, pheated)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: ki, kj, kn
    REAL(KIND = wp), INTENT(IN) :: pcalved
    REAL(KIND = wp), INTENT(IN) :: pheated
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    IF (.NOT. ln_bergdia) RETURN
    CALL profile_psy_data0 % PreStart('icb_dia_calve', 'r0', 0, 0)
    real_calving(ki, kj, kn) = real_calving(ki, kj, kn) + pcalved / berg_dt
    nbergs_calved = nbergs_calved + 1
    nbergs_calved_by_class(kn) = nbergs_calved_by_class(kn) + 1
    calving_to_bergs_net = calving_to_bergs_net + pcalved
    heat_to_bergs_net = heat_to_bergs_net + pheated
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_dia_calve
  SUBROUTINE icb_dia_income(kt, pcalving_used, pheat_used)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), INTENT(IN) :: pcalving_used
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pheat_used
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    IF (.NOT. ln_bergdia) RETURN
    CALL profile_psy_data0 % PreStart('icb_dia_income', 'r0', 0, 0)
    IF (kt == nit000) THEN
      stored_start = SUM(berg_grid % stored_ice(:, :, :))
      CALL mpp_sum('icbdia', stored_start)
      stored_heat_start = SUM(berg_grid % stored_heat(:, :))
      CALL mpp_sum('icbdia', stored_heat_start)
      IF (nn_verbose_level > 0) THEN
        WRITE(numicb, FMT = '(a,es13.6,a)') 'icb_dia_income: initial stored mass=', stored_start, ' kg'
        WRITE(numicb, FMT = '(a,es13.6,a)') 'icb_dia_income: initial stored heat=', stored_heat_start, ' J'
      END IF
    END IF
    calving_rcv_net = calving_rcv_net + SUM(berg_grid % calving(:, :)) * berg_dt
    calving_src_net = calving_rcv_net
    calving_src_heat_net = calving_src_heat_net + SUM(berg_grid % calving_hflx(:, :) * e1e2t(:, :)) * berg_dt
    calving_used_net = calving_used_net + pcalving_used * berg_dt
    calving_src_heat_used_net = calving_src_heat_used_net + SUM(pheat_used(:, :))
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_dia_income
  SUBROUTINE icb_dia_size(ki, kj, pWn, pLn, pAbits, pmass_scale, pMnew, pnMbits, pz1_e1e2)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: ki, kj
    REAL(KIND = wp), INTENT(IN) :: pWn, pLn, pAbits, pmass_scale, pMnew, pnMbits, pz1_e1e2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    IF (.NOT. ln_bergdia) RETURN
    CALL profile_psy_data0 % PreStart('icb_dia_size', 'r0', 0, 0)
    virtual_area(ki, kj) = virtual_area(ki, kj) + (pWn * pLn + pAbits) * pmass_scale
    berg_mass(ki, kj) = berg_mass(ki, kj) + pMnew * pz1_e1e2
    bits_mass(ki, kj) = bits_mass(ki, kj) + pnMbits * pz1_e1e2
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_dia_size
  SUBROUTINE icb_dia_speed
    IF (.NOT. ln_bergdia) RETURN
    nspeeding_tickets = nspeeding_tickets + 1
  END SUBROUTINE icb_dia_speed
  SUBROUTINE icb_dia_melt(ki, kj, pmnew, pheat_hcflux, pheat_latent, pmass_scale, pdM, pdMbitsE, pdMbitsM, pdMb, pdMe, pdMv, &
&pz1_dt_e1e2)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: ki, kj
    REAL(KIND = wp), INTENT(IN) :: pmnew, pheat_hcflux, pheat_latent, pmass_scale
    REAL(KIND = wp), INTENT(IN) :: pdM, pdMbitsE, pdMbitsM, pdMb, pdMe, pdMv, pz1_dt_e1e2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    IF (.NOT. ln_bergdia) RETURN
    CALL profile_psy_data0 % PreStart('icb_dia_melt', 'r0', 0, 0)
    berg_melt(ki, kj) = berg_melt(ki, kj) + pdM * pz1_dt_e1e2
    berg_melt_hcflx(ki, kj) = berg_melt_hcflx(ki, kj) + pheat_hcflux * pz1_dt_e1e2
    berg_melt_qlat(ki, kj) = berg_melt_qlat(ki, kj) + pheat_latent * pz1_dt_e1e2
    bits_src(ki, kj) = bits_src(ki, kj) + pdMbitsE * pz1_dt_e1e2
    bits_melt(ki, kj) = bits_melt(ki, kj) + pdMbitsM * pz1_dt_e1e2
    buoy_melt(ki, kj) = buoy_melt(ki, kj) + pdMb * pz1_dt_e1e2
    eros_melt(ki, kj) = eros_melt(ki, kj) + pdMe * pz1_dt_e1e2
    conv_melt(ki, kj) = conv_melt(ki, kj) + pdMv * pz1_dt_e1e2
    heat_to_ocean_net = heat_to_ocean_net + (pheat_hcflux + pheat_latent) * pmass_scale * berg_dt
    IF (pmnew <= 0._wp) nbergs_melted = nbergs_melted + 1
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_dia_melt
  SUBROUTINE report_state(cd_budgetstr, cd_budgetunits, cd_startstr, pstartval, cd_endstr, pendval, cd_delstr, kbergs)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER*(*), INTENT(IN) :: cd_budgetstr, cd_budgetunits, cd_startstr, cd_endstr, cd_delstr
    REAL(KIND = wp), INTENT(IN) :: pstartval, pendval
    INTEGER, INTENT(IN), OPTIONAL :: kbergs
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    IF (nn_verbose_level == 0) RETURN
    CALL profile_psy_data0 % PreStart('report_state', 'r0', 0, 0)
    IF (PRESENT(kbergs)) THEN
      WRITE(numicb, 100) cd_budgetstr // ' state:', cd_startstr // ' start', pstartval, cd_budgetunits, cd_endstr // ' end', &
&pendval, cd_budgetunits, 'Delta ' // cd_delstr, pendval - pstartval, cd_budgetunits, '# of bergs', kbergs
    ELSE
      WRITE(numicb, 100) cd_budgetstr // ' state:', cd_startstr // ' start', pstartval, cd_budgetunits, cd_endstr // ' end', &
&pendval, cd_budgetunits, cd_delstr // 'Delta', pendval - pstartval, cd_budgetunits
    END IF
100 FORMAT(A19, 3(A18, "=", ES14.7, X, A2, :, ","), A12, I8)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE report_state
  SUBROUTINE report_consistant(cd_budgetstr, cd_budgetunits, cd_startstr, pstartval, cd_endstr, pendval)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER*(*), INTENT(IN) :: cd_budgetstr, cd_budgetunits, cd_startstr, cd_endstr
    REAL(KIND = wp), INTENT(IN) :: pstartval, pendval
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    IF (nn_verbose_level == 0) RETURN
    CALL profile_psy_data0 % PreStart('report_consistant', 'r0', 0, 0)
    WRITE(numicb, 200) cd_budgetstr // ' check:', cd_startstr, pstartval, cd_budgetunits, cd_endstr, pendval, cd_budgetunits, &
&'error', (pendval - pstartval) / ((pendval + pstartval) + 1E-30), 'nd'
200 FORMAT(A19, 10(A18, "=", ES14.7, X, A2, :, ","))
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE report_consistant
  SUBROUTINE report_budget(cd_budgetstr, cd_budgetunits, cd_instr, pinval, cd_outstr, poutval, cd_delstr, pstartval, pendval)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER*(*), INTENT(IN) :: cd_budgetstr, cd_budgetunits, cd_instr, cd_outstr, cd_delstr
    REAL(KIND = wp), INTENT(IN) :: pinval, poutval, pstartval, pendval
    REAL(KIND = wp) :: zval
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    IF (nn_verbose_level == 0) RETURN
    CALL profile_psy_data0 % PreStart('report_budget', 'r0', 0, 0)
    zval = ((pendval - pstartval) - (pinval - poutval)) / MAX(1.E-30, MAX(ABS(pendval - pstartval), ABS(pinval - poutval)))
    WRITE(numicb, 200) cd_budgetstr // ' budget:', cd_instr // ' in', pinval, cd_budgetunits, cd_outstr // ' out', poutval, &
&cd_budgetunits, 'Delta ' // cd_delstr, pinval - poutval, cd_budgetunits, 'error', zval, 'nd'
200 FORMAT(A19, 3(A18, "=", ES14.7, X, A2, :, ","), A8, "=", ES10.3, X, A2)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE report_budget
  SUBROUTINE report_istate(cd_budgetstr, cd_startstr, pstartval, cd_endstr, pendval, cd_delstr)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER*(*), INTENT(IN) :: cd_budgetstr, cd_startstr, cd_endstr, cd_delstr
    INTEGER, INTENT(IN) :: pstartval, pendval
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    IF (nn_verbose_level == 0) RETURN
    CALL profile_psy_data0 % PreStart('report_istate', 'r0', 0, 0)
    WRITE(numicb, 100) cd_budgetstr // ' state:', cd_startstr // ' start', pstartval, cd_endstr // ' end', pendval, cd_delstr // &
&'Delta', pendval - pstartval
100 FORMAT(A19, 3(A18, "=", I14, X, :, ","))
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE report_istate
  SUBROUTINE report_ibudget(cd_budgetstr, cd_instr, pinval, cd_outstr, poutval, cd_delstr, pstartval, pendval)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER*(*), INTENT(IN) :: cd_budgetstr, cd_instr, cd_outstr, cd_delstr
    INTEGER, INTENT(IN) :: pinval, poutval, pstartval, pendval
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    IF (nn_verbose_level == 0) RETURN
    CALL profile_psy_data0 % PreStart('report_ibudget', 'r0', 0, 0)
    WRITE(numicb, 200) cd_budgetstr // ' budget:', cd_instr // ' in', pinval, cd_outstr // ' out', poutval, 'Delta ' // cd_delstr, &
&pinval - poutval, 'error', ((pendval - pstartval) - (pinval - poutval))
200 FORMAT(A19, 10(A18, "=", I14, X, :, ","))
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE report_ibudget
END MODULE icbdia