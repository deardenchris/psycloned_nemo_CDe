MODULE dynspg
  USE oce
  USE dom_oce
  USE c1d
  USE phycst
  USE sbc_oce
  USE sbc_ice, ONLY: snwice_mass, snwice_mass_b
  USE sbcapr
  USE dynspg_exp
  USE dynspg_ts
  USE sbctide
  USE updtide
  USE trd_oce
  USE trddyn
  USE prtctl
  USE in_out_manager
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_spg
  PUBLIC :: dyn_spg_init
  INTEGER :: nspg = 0
  INTEGER, PARAMETER :: np_TS = 1
  INTEGER, PARAMETER :: np_EXP = 0
  INTEGER, PARAMETER :: np_NO = - 1
  CONTAINS
  SUBROUTINE dyn_spg(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: z2dt, zg_2, zintp, zgrau0r, zld
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zpice
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ztrdu, ztrdv
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    IF (ln_timing) CALL timing_start('dyn_spg')
    IF (l_trddyn) THEN
      ALLOCATE(ztrdu(jpi, jpj, jpk), ztrdv(jpi, jpj, jpk))
      !$ACC KERNELS
      ztrdu(:, :, :) = ua(:, :, :)
      ztrdv(:, :, :) = va(:, :, :)
      !$ACC END KERNELS
    END IF
    IF (ln_apr_dyn .OR. (.NOT. ln_dynspg_ts .AND. (ln_tide_pot .AND. ln_tide)) .OR. ln_ice_embd) THEN
      !$ACC KERNELS
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          spgu(ji, jj) = 0._wp
          spgv(ji, jj) = 0._wp
        END DO
      END DO
      IF (ln_apr_dyn .AND. .NOT. ln_dynspg_ts) THEN
        zg_2 = grav * 0.5
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            spgu(ji, jj) = spgu(ji, jj) + zg_2 * (ssh_ib(ji + 1, jj) - ssh_ib(ji, jj) + ssh_ibb(ji + 1, jj) - ssh_ibb(ji, jj)) * &
&r1_e1u(ji, jj)
            spgv(ji, jj) = spgv(ji, jj) + zg_2 * (ssh_ib(ji, jj + 1) - ssh_ib(ji, jj) + ssh_ibb(ji, jj + 1) - ssh_ibb(ji, jj)) * &
&r1_e2v(ji, jj)
          END DO
        END DO
      END IF
      !$ACC END KERNELS
      IF (.NOT. ln_dynspg_ts .AND. (ln_tide_pot .AND. ln_tide)) THEN
        CALL upd_tide(kt)
        !$ACC KERNELS
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            spgu(ji, jj) = spgu(ji, jj) + grav * (pot_astro(ji + 1, jj) - pot_astro(ji, jj)) * r1_e1u(ji, jj)
            spgv(ji, jj) = spgv(ji, jj) + grav * (pot_astro(ji, jj + 1) - pot_astro(ji, jj)) * r1_e2v(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
        IF (ln_scal_load) THEN
          !$ACC KERNELS
          zld = rn_scal_load * grav
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              spgu(ji, jj) = spgu(ji, jj) + zld * (sshn(ji + 1, jj) - sshn(ji, jj)) * r1_e1u(ji, jj)
              spgv(ji, jj) = spgv(ji, jj) + zld * (sshn(ji, jj + 1) - sshn(ji, jj)) * r1_e2v(ji, jj)
            END DO
          END DO
          !$ACC END KERNELS
        END IF
      END IF
      IF (ln_ice_embd) THEN
        ALLOCATE(zpice(jpi, jpj))
        !$ACC KERNELS
        zintp = REAL(MOD(kt - 1, nn_fsbc)) / REAL(nn_fsbc)
        zgrau0r = - grav * r1_rau0
        zpice(:, :) = (zintp * snwice_mass(:, :) + (1. - zintp) * snwice_mass_b(:, :)) * zgrau0r
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            spgu(ji, jj) = spgu(ji, jj) + (zpice(ji + 1, jj) - zpice(ji, jj)) * r1_e1u(ji, jj)
            spgv(ji, jj) = spgv(ji, jj) + (zpice(ji, jj + 1) - zpice(ji, jj)) * r1_e2v(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
        DEALLOCATE(zpice)
      END IF
      !$ACC KERNELS
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            ua(ji, jj, jk) = ua(ji, jj, jk) + spgu(ji, jj)
            va(ji, jj, jk) = va(ji, jj, jk) + spgv(ji, jj)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data0 % PreStart('dyn_spg', 'r0', 0, 0)
    SELECT CASE (nspg)
    CASE (np_exp)
      CALL dyn_spg_exp(kt)
    CASE (np_ts)
      CALL dyn_spg_ts(kt)
    END SELECT
    CALL profile_psy_data0 % PostEnd
    IF (l_trddyn) THEN
      !$ACC KERNELS
      ztrdu(:, :, :) = ua(:, :, :) - ztrdu(:, :, :)
      ztrdv(:, :, :) = va(:, :, :) - ztrdv(:, :, :)
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('dyn_spg', 'r1', 0, 0)
      CALL trd_dyn(ztrdu, ztrdv, jpdyn_spg, kt)
      DEALLOCATE(ztrdu, ztrdv)
      CALL profile_psy_data1 % PostEnd
    END IF
    CALL profile_psy_data2 % PreStart('dyn_spg', 'r2', 0, 0)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = ua, clinfo1 = ' spg  - Ua: ', mask1 = umask, tab3d_2 = va, clinfo2 = ' Va: ', mask2 = &
&vmask, clinfo3 = 'dyn')
    IF (ln_timing) CALL timing_stop('dyn_spg')
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE dyn_spg
  SUBROUTINE dyn_spg_init
    INTEGER :: ioptio, ios
    NAMELIST /namdyn_spg/ ln_dynspg_exp, ln_dynspg_ts, ln_bt_fw, ln_bt_av, ln_bt_auto, nn_baro, rn_bt_cmax, nn_bt_flt, rn_bt_alpha
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dyn_spg_init : choice of the surface pressure gradient scheme'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
    END IF
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namdyn_spg, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namdyn_spg in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namdyn_spg, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namdyn_spg in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namdyn_spg)
    IF (lwp) THEN
      WRITE(numout, FMT = *) '   Namelist : namdyn_spg                    '
      WRITE(numout, FMT = *) '      Explicit free surface                  ln_dynspg_exp = ', ln_dynspg_exp
      WRITE(numout, FMT = *) '      Free surface with time splitting       ln_dynspg_ts  = ', ln_dynspg_ts
    END IF
    nspg = np_no
    ioptio = 0
    IF (ln_dynspg_exp) THEN
      nspg = np_exp
      ioptio = ioptio + 1
    END IF
    IF (ln_dynspg_ts) THEN
      nspg = np_ts
      ioptio = ioptio + 1
    END IF
    IF (ioptio > 1) CALL ctl_stop('Choose only one surface pressure gradient scheme')
    IF (ioptio == 0) CALL ctl_warn('NO surface pressure gradient trend in momentum Eqs.')
    IF (ln_dynspg_exp .AND. ln_isfcav) CALL ctl_stop(' dynspg_exp not tested with ice shelf cavity ')
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      IF (nspg == np_EXP) WRITE(numout, FMT = *) '   ==>>>   explicit free surface'
      IF (nspg == np_TS) WRITE(numout, FMT = *) '   ==>>>   free surface with time splitting scheme'
      IF (nspg == np_NO) WRITE(numout, FMT = *) '   ==>>>   No surface surface pressure gradient trend in momentum Eqs.'
    END IF
    IF (nspg == np_TS) THEN
      CALL dyn_spg_ts_init
    END IF
  END SUBROUTINE dyn_spg_init
END MODULE dynspg