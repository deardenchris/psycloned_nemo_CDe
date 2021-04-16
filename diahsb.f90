MODULE diahsb
  USE oce
  USE dom_oce
  USE phycst
  USE sbc_oce
  USE sbcrnf
  USE sbcisf
  USE domvvl
  USE traqsr
  USE trabbc
  USE trabbc
  USE restart
  USE bdy_oce, ONLY: ln_bdy
  USE iom
  USE in_out_manager
  USE lib_fortran
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dia_hsb
  PUBLIC :: dia_hsb_init
  LOGICAL, PUBLIC :: ln_diahsb
  REAL(KIND = wp) :: surf_tot
  REAL(KIND = wp) :: frc_t, frc_s, frc_v
  REAL(KIND = wp) :: frc_wn_t, frc_wn_s
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: surf
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: surf_ini, ssh_ini
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: ssh_hc_loc_ini, ssh_sc_loc_ini
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: hc_loc_ini, sc_loc_ini, e3t_ini
  CONTAINS
  SUBROUTINE dia_hsb(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zdiff_hc, zdiff_sc
    REAL(KIND = wp) :: zdiff_hc1, zdiff_sc1
    REAL(KIND = wp) :: zdiff_v1, zdiff_v2
    REAL(KIND = wp) :: zerr_hc1, zerr_sc1
    REAL(KIND = wp) :: zvol_tot
    REAL(KIND = wp) :: z_frc_trd_t, z_frc_trd_s
    REAL(KIND = wp) :: z_frc_trd_v
    REAL(KIND = wp) :: z_wn_trd_t, z_wn_trd_s
    REAL(KIND = wp) :: z_ssh_hc, z_ssh_sc
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z2d0, z2d1
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpkm1) :: zwrk
    IF (ln_timing) CALL timing_start('dia_hsb')
    tsn(:, :, :, 1) = tsn(:, :, :, 1) * tmask(:, :, :)
    tsb(:, :, :, 1) = tsb(:, :, :, 1) * tmask(:, :, :)
    tsn(:, :, :, 2) = tsn(:, :, :, 2) * tmask(:, :, :)
    tsb(:, :, :, 2) = tsb(:, :, :, 2) * tmask(:, :, :)
    z_frc_trd_v = r1_rau0 * glob_sum('diahsb', - (emp(:, :) - rnf(:, :) + fwfisf(:, :)) * surf(:, :))
    z_frc_trd_t = glob_sum('diahsb', sbc_tsc(:, :, jp_tem) * surf(:, :))
    z_frc_trd_s = glob_sum('diahsb', sbc_tsc(:, :, jp_sal) * surf(:, :))
    IF (ln_rnf) z_frc_trd_t = z_frc_trd_t + glob_sum('diahsb', rnf_tsc(:, :, jp_tem) * surf(:, :))
    IF (ln_rnf_sal) z_frc_trd_s = z_frc_trd_s + glob_sum('diahsb', rnf_tsc(:, :, jp_sal) * surf(:, :))
    IF (ln_isf) z_frc_trd_t = z_frc_trd_t + glob_sum('diahsb', risf_tsc(:, :, jp_tem) * surf(:, :))
    IF (ln_traqsr) z_frc_trd_t = z_frc_trd_t + r1_rau0_rcp * glob_sum('diahsb', qsr(:, :) * surf(:, :))
    IF (ln_trabbc) z_frc_trd_t = z_frc_trd_t + glob_sum('diahsb', qgh_trd0(:, :) * surf(:, :))
    IF (ln_linssh) THEN
      IF (ln_isfcav) THEN
        DO ji = 1, jpi
          DO jj = 1, jpj
            z2d0(ji, jj) = surf(ji, jj) * wn(ji, jj, mikt(ji, jj)) * tsb(ji, jj, mikt(ji, jj), jp_tem)
            z2d1(ji, jj) = surf(ji, jj) * wn(ji, jj, mikt(ji, jj)) * tsb(ji, jj, mikt(ji, jj), jp_sal)
          END DO
        END DO
      ELSE
        z2d0(:, :) = surf(:, :) * wn(:, :, 1) * tsb(:, :, 1, jp_tem)
        z2d1(:, :) = surf(:, :) * wn(:, :, 1) * tsb(:, :, 1, jp_sal)
      END IF
      z_wn_trd_t = - glob_sum('diahsb', z2d0)
      z_wn_trd_s = - glob_sum('diahsb', z2d1)
    END IF
    frc_v = frc_v + z_frc_trd_v * rdt
    frc_t = frc_t + z_frc_trd_t * rdt
    frc_s = frc_s + z_frc_trd_s * rdt
    IF (ln_linssh) THEN
      frc_wn_t = frc_wn_t + z_wn_trd_t * rdt
      frc_wn_s = frc_wn_s + z_wn_trd_s * rdt
    END IF
    zdiff_v1 = glob_sum_full('diahsb', surf(:, :) * sshn(:, :) - surf_ini(:, :) * ssh_ini(:, :))
    IF (ln_linssh) THEN
      IF (ln_isfcav) THEN
        DO ji = 1, jpi
          DO jj = 1, jpj
            z2d0(ji, jj) = surf(ji, jj) * (tsn(ji, jj, mikt(ji, jj), jp_tem) * sshn(ji, jj) - ssh_hc_loc_ini(ji, jj))
            z2d1(ji, jj) = surf(ji, jj) * (tsn(ji, jj, mikt(ji, jj), jp_sal) * sshn(ji, jj) - ssh_sc_loc_ini(ji, jj))
          END DO
        END DO
      ELSE
        z2d0(:, :) = surf(:, :) * (tsn(:, :, 1, jp_tem) * sshn(:, :) - ssh_hc_loc_ini(:, :))
        z2d1(:, :) = surf(:, :) * (tsn(:, :, 1, jp_sal) * sshn(:, :) - ssh_sc_loc_ini(:, :))
      END IF
      z_ssh_hc = glob_sum_full('diahsb', z2d0)
      z_ssh_sc = glob_sum_full('diahsb', z2d1)
    END IF
    !$OMP parallel default(shared), private(jk)
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      zwrk(:, :, jk) = (surf(:, :) * e3t_n(:, :, jk) - surf_ini(:, :) * e3t_ini(:, :, jk)) * tmask(:, :, jk)
    END DO
    !$OMP end do
    !$OMP end parallel
    zdiff_v2 = glob_sum_full('diahsb', zwrk(:, :, :))
    !$OMP parallel default(shared), private(jk)
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      zwrk(:, :, jk) = (surf(:, :) * e3t_n(:, :, jk) * tsn(:, :, jk, jp_tem) - surf_ini(:, :) * hc_loc_ini(:, :, jk)) * tmask(:, :, jk)
    END DO
    !$OMP end do
    !$OMP end parallel
    zdiff_hc = glob_sum_full('diahsb', zwrk(:, :, :))
    !$OMP parallel default(shared), private(jk)
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      zwrk(:, :, jk) = (surf(:, :) * e3t_n(:, :, jk) * tsn(:, :, jk, jp_sal) - surf_ini(:, :) * sc_loc_ini(:, :, jk)) * tmask(:, :, jk)
    END DO
    !$OMP end do
    !$OMP end parallel
    zdiff_sc = glob_sum_full('diahsb', zwrk(:, :, :))
    zdiff_v1 = zdiff_v1 - frc_v
    IF (.NOT. ln_linssh) zdiff_v2 = zdiff_v2 - frc_v
    zdiff_hc = zdiff_hc - frc_t
    zdiff_sc = zdiff_sc - frc_s
    IF (ln_linssh) THEN
      zdiff_hc1 = zdiff_hc + z_ssh_hc
      zdiff_sc1 = zdiff_sc + z_ssh_sc
      zerr_hc1 = z_ssh_hc - frc_wn_t
      zerr_sc1 = z_ssh_sc - frc_wn_s
    END IF
    !$OMP parallel default(shared), private(jk)
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      zwrk(:, :, jk) = surf(:, :) * e3t_n(:, :, jk) * tmask(:, :, jk)
    END DO
    !$OMP end do
    !$OMP end parallel
    zvol_tot = glob_sum_full('diahsb', zwrk(:, :, :))
    CALL iom_put('bgfrcvol', frc_v * 1.E-9)
    CALL iom_put('bgfrctem', frc_t * rau0 * rcp * 1.E-20)
    CALL iom_put('bgfrchfx', frc_t * rau0 * rcp / (surf_tot * kt * rdt))
    CALL iom_put('bgfrcsal', frc_s * 1.E-9)
    IF (.NOT. ln_linssh) THEN
      CALL iom_put('bgtemper', zdiff_hc / zvol_tot)
      CALL iom_put('bgsaline', zdiff_sc / zvol_tot)
      CALL iom_put('bgheatco', zdiff_hc * 1.E-20 * rau0 * rcp)
      CALL iom_put('bgheatfx', zdiff_hc * rau0 * rcp / (surf_tot * kt * rdt))
      CALL iom_put('bgsaltco', zdiff_sc * 1.E-9)
      CALL iom_put('bgvolssh', zdiff_v1 * 1.E-9)
      CALL iom_put('bgvole3t', zdiff_v2 * 1.E-9)
      IF (kt == nitend .AND. lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) 'dia_hsb : last time step hsb diagnostics: at it= ', kt, ' date= ', ndastp
        WRITE(numout, FMT = *) '~~~~~~~'
        WRITE(numout, FMT = *) '   Temperature drift = ', zdiff_hc / zvol_tot, ' C'
        WRITE(numout, FMT = *) '   Salinity    drift = ', zdiff_sc / zvol_tot, ' PSU'
        WRITE(numout, FMT = *) '   volume ssh  drift = ', zdiff_v1 * 1.E-9, ' km^3'
        WRITE(numout, FMT = *) '   volume e3t  drift = ', zdiff_v2 * 1.E-9, ' km^3'
      END IF
    ELSE
      CALL iom_put('bgtemper', zdiff_hc1 / zvol_tot)
      CALL iom_put('bgsaline', zdiff_sc1 / zvol_tot)
      CALL iom_put('bgheatco', zdiff_hc1 * 1.E-20 * rau0 * rcp)
      CALL iom_put('bgheatfx', zdiff_hc1 * rau0 * rcp / (surf_tot * kt * rdt))
      CALL iom_put('bgsaltco', zdiff_sc1 * 1.E-9)
      CALL iom_put('bgvolssh', zdiff_v1 * 1.E-9)
      CALL iom_put('bgmistem', zerr_hc1 / zvol_tot)
      CALL iom_put('bgmissal', zerr_sc1 / zvol_tot)
    END IF
    IF (lrst_oce) CALL dia_hsb_rst(kt, 'WRITE')
    IF (ln_timing) CALL timing_stop('dia_hsb')
  END SUBROUTINE dia_hsb
  SUBROUTINE dia_hsb_rst(kt, cdrw)
    INTEGER, INTENT(IN) :: kt
    CHARACTER(LEN = *), INTENT(IN) :: cdrw
    INTEGER :: ji, jj, jk
    IF (TRIM(cdrw) == 'READ') THEN
      IF (ln_rstart) THEN
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) '   dia_hsb_rst : read hsb restart at it= ', kt, ' date= ', ndastp
        IF (lwp) WRITE(numout, FMT = *)
        CALL iom_get(numror, 'frc_v', frc_v, ldxios = lrxios)
        CALL iom_get(numror, 'frc_t', frc_t, ldxios = lrxios)
        CALL iom_get(numror, 'frc_s', frc_s, ldxios = lrxios)
        IF (ln_linssh) THEN
          CALL iom_get(numror, 'frc_wn_t', frc_wn_t, ldxios = lrxios)
          CALL iom_get(numror, 'frc_wn_s', frc_wn_s, ldxios = lrxios)
        END IF
        CALL iom_get(numror, jpdom_autoglo, 'surf_ini', surf_ini, ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'ssh_ini', ssh_ini, ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'e3t_ini', e3t_ini, ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'hc_loc_ini', hc_loc_ini, ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'sc_loc_ini', sc_loc_ini, ldxios = lrxios)
        IF (ln_linssh) THEN
          CALL iom_get(numror, jpdom_autoglo, 'ssh_hc_loc_ini', ssh_hc_loc_ini, ldxios = lrxios)
          CALL iom_get(numror, jpdom_autoglo, 'ssh_sc_loc_ini', ssh_sc_loc_ini, ldxios = lrxios)
        END IF
      ELSE
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) '   dia_hsb_rst : initialise hsb at initial state '
        IF (lwp) WRITE(numout, FMT = *)
        surf_ini(:, :) = e1e2t(:, :) * tmask_i(:, :)
        ssh_ini(:, :) = sshn(:, :)
        !$OMP parallel default(shared), private(jk)
        !$OMP do schedule(static)
        DO jk = 1, jpk
          e3t_ini(:, :, jk) = e3t_n(:, :, jk) * tmask(:, :, jk)
          hc_loc_ini(:, :, jk) = tsn(:, :, jk, jp_tem) * e3t_n(:, :, jk) * tmask(:, :, jk)
          sc_loc_ini(:, :, jk) = tsn(:, :, jk, jp_sal) * e3t_n(:, :, jk) * tmask(:, :, jk)
        END DO
        !$OMP end do
        !$OMP end parallel
        frc_v = 0._wp
        frc_t = 0._wp
        frc_s = 0._wp
        IF (ln_linssh) THEN
          IF (ln_isfcav) THEN
            DO ji = 1, jpi
              DO jj = 1, jpj
                ssh_hc_loc_ini(ji, jj) = tsn(ji, jj, mikt(ji, jj), jp_tem) * sshn(ji, jj)
                ssh_sc_loc_ini(ji, jj) = tsn(ji, jj, mikt(ji, jj), jp_sal) * sshn(ji, jj)
              END DO
            END DO
          ELSE
            ssh_hc_loc_ini(:, :) = tsn(:, :, 1, jp_tem) * sshn(:, :)
            ssh_sc_loc_ini(:, :) = tsn(:, :, 1, jp_sal) * sshn(:, :)
          END IF
          frc_wn_t = 0._wp
          frc_wn_s = 0._wp
        END IF
      END IF
    ELSE IF (TRIM(cdrw) == 'WRITE') THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) '   dia_hsb_rst : write restart at it= ', kt, ' date= ', ndastp
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwxios) CALL iom_swap(cwxios_context)
      CALL iom_rstput(kt, nitrst, numrow, 'frc_v', frc_v, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'frc_t', frc_t, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'frc_s', frc_s, ldxios = lwxios)
      IF (ln_linssh) THEN
        CALL iom_rstput(kt, nitrst, numrow, 'frc_wn_t', frc_wn_t, ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'frc_wn_s', frc_wn_s, ldxios = lwxios)
      END IF
      CALL iom_rstput(kt, nitrst, numrow, 'surf_ini', surf_ini, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'ssh_ini', ssh_ini, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'e3t_ini', e3t_ini, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'hc_loc_ini', hc_loc_ini, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'sc_loc_ini', sc_loc_ini, ldxios = lwxios)
      IF (ln_linssh) THEN
        CALL iom_rstput(kt, nitrst, numrow, 'ssh_hc_loc_ini', ssh_hc_loc_ini, ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'ssh_sc_loc_ini', ssh_sc_loc_ini, ldxios = lwxios)
      END IF
      IF (lwxios) CALL iom_swap(cxios_context)
    END IF
  END SUBROUTINE dia_hsb_rst
  SUBROUTINE dia_hsb_init
    INTEGER :: ierror, ios
    NAMELIST /namhsb/ ln_diahsb
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dia_hsb_init : heat and salt budgets diagnostics'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~ '
    END IF
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namhsb, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namhsb in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namhsb, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namhsb in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namhsb)
    IF (lwp) THEN
      WRITE(numout, FMT = *) '   Namelist  namhsb :'
      WRITE(numout, FMT = *) '      check the heat and salt budgets (T) or not (F)       ln_diahsb = ', ln_diahsb
    END IF
    IF (.NOT. ln_diahsb) RETURN
    IF (lwxios) THEN
      CALL iom_set_rstw_var_active('frc_v')
      CALL iom_set_rstw_var_active('frc_t')
      CALL iom_set_rstw_var_active('frc_s')
      CALL iom_set_rstw_var_active('surf_ini')
      CALL iom_set_rstw_var_active('ssh_ini')
      CALL iom_set_rstw_var_active('e3t_ini')
      CALL iom_set_rstw_var_active('hc_loc_ini')
      CALL iom_set_rstw_var_active('sc_loc_ini')
      IF (ln_linssh) THEN
        CALL iom_set_rstw_var_active('ssh_hc_loc_ini')
        CALL iom_set_rstw_var_active('ssh_sc_loc_ini')
        CALL iom_set_rstw_var_active('frc_wn_t')
        CALL iom_set_rstw_var_active('frc_wn_s')
      END IF
    END IF
    ALLOCATE(hc_loc_ini(jpi, jpj, jpk), sc_loc_ini(jpi, jpj, jpk), surf_ini(jpi, jpj), e3t_ini(jpi, jpj, jpk), surf(jpi, jpj), ssh_ini(jpi, jpj), STAT = ierror)
    IF (ierror > 0) THEN
      CALL ctl_stop('dia_hsb_init: unable to allocate hc_loc_ini')
      RETURN
    END IF
    IF (ln_linssh) ALLOCATE(ssh_hc_loc_ini(jpi, jpj), ssh_sc_loc_ini(jpi, jpj), STAT = ierror)
    IF (ierror > 0) THEN
      CALL ctl_stop('dia_hsb: unable to allocate ssh_hc_loc_ini')
      RETURN
    END IF
    surf(:, :) = e1e2t(:, :) * tmask_i(:, :)
    surf_tot = glob_sum('diahsb', surf(:, :))
    IF (ln_bdy) CALL ctl_warn('dia_hsb_init: heat/salt budget does not consider open boundary fluxes')
    CALL dia_hsb_rst(nit000, 'READ')
  END SUBROUTINE dia_hsb_init
END MODULE diahsb