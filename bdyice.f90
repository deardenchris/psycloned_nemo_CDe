MODULE bdyice
  USE oce
  USE ice
  USE icevar
  USE icecor
  USE icectl
  USE phycst
  USE eosbn2
  USE par_oce
  USE dom_oce
  USE sbc_oce
  USE bdy_oce
  USE lbclnk
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: bdy_ice
  PUBLIC :: bdy_ice_dyn
  CONTAINS
  SUBROUTINE bdy_ice(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: jbdy
    IF (ln_timing) CALL timing_start('bdy_ice_thd')
    IF (ln_icediachk) CALL ice_cons_hsm(0, 'bdy_ice_thd', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    CALL ice_var_glo2eqv
    DO jbdy = 1, nb_bdy
      SELECT CASE (cn_ice(jbdy))
      CASE ('none')
        CYCLE
      CASE ('frs')
        CALL bdy_ice_frs(idx_bdy(jbdy), dta_bdy(jbdy), kt, jbdy)
      CASE DEFAULT
        CALL ctl_stop('bdy_ice : unrecognised option for open boundaries for ice fields')
      END SELECT
    END DO
    CALL ice_cor(kt, 0)
    CALL ice_var_agg(1)
    IF (ln_icediachk) CALL ice_cons_hsm(1, 'bdy_ice_thd', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    IF (ln_icectl) CALL ice_prt(kt, iiceprt, jiceprt, 1, ' - ice thermo bdy - ')
    IF (ln_timing) CALL timing_stop('bdy_ice_thd')
  END SUBROUTINE bdy_ice
  SUBROUTINE bdy_ice_frs(idx, dta, kt, jbdy)
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    TYPE(OBC_DATA), INTENT(IN) :: dta
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: jbdy
    INTEGER :: jpbound
    INTEGER :: i_bdy, jgrd
    INTEGER :: ji, jj, jk, jl, ib, jb
    REAL(KIND = wp) :: zwgt, zwgt1
    REAL(KIND = wp) :: ztmelts, zdh
    jgrd = 1
    DO jl = 1, jpl
      DO i_bdy = 1, idx % nblenrim(jgrd)
        ji = idx % nbi(i_bdy, jgrd)
        jj = idx % nbj(i_bdy, jgrd)
        zwgt = idx % nbw(i_bdy, jgrd)
        zwgt1 = 1.E0 - idx % nbw(i_bdy, jgrd)
        a_i(ji, jj, jl) = (a_i(ji, jj, jl) * zwgt1 + dta % a_i(i_bdy, jl) * zwgt) * tmask(ji, jj, 1)
        h_i(ji, jj, jl) = (h_i(ji, jj, jl) * zwgt1 + dta % h_i(i_bdy, jl) * zwgt) * tmask(ji, jj, 1)
        h_s(ji, jj, jl) = (h_s(ji, jj, jl) * zwgt1 + dta % h_s(i_bdy, jl) * zwgt) * tmask(ji, jj, 1)
        zdh = MAX(0._wp, (rhos * h_s(ji, jj, jl) + (rhoi - rau0) * h_i(ji, jj, jl)) * r1_rau0)
        h_i(ji, jj, jl) = MIN(hi_max(jl), h_i(ji, jj, jl) + zdh)
        h_s(ji, jj, jl) = MAX(0._wp, h_s(ji, jj, jl) - zdh * rhoi / rhos)
      END DO
    END DO
    CALL lbc_bdy_lnk('bdyice', a_i(:, :, :), 'T', 1., jbdy)
    CALL lbc_bdy_lnk('bdyice', h_i(:, :, :), 'T', 1., jbdy)
    CALL lbc_bdy_lnk('bdyice', h_s(:, :, :), 'T', 1., jbdy)
    DO jl = 1, jpl
      DO i_bdy = 1, idx % nblenrim(jgrd)
        ji = idx % nbi(i_bdy, jgrd)
        jj = idx % nbj(i_bdy, jgrd)
        jpbound = 0
        ib = ji
        jb = jj
        IF (u_ice(ji, jj) < 0. .AND. umask(ji - 1, jj, 1) == 0.) jpbound = 1
        ib = ji + 1
        IF (u_ice(ji - 1, jj) > 0. .AND. umask(ji, jj, 1) == 0.) jpbound = 1
        ib = ji - 1
        IF (v_ice(ji, jj) < 0. .AND. vmask(ji, jj - 1, 1) == 0.) jpbound = 1
        jb = jj + 1
        IF (v_ice(ji, jj - 1) > 0. .AND. vmask(ji, jj, 1) == 0.) jpbound = 1
        jb = jj - 1
        IF (nn_ice_dta(jbdy) == 0) jpbound = 0
        ib = ji
        jb = jj
        IF (a_i(ib, jb, jl) > 0._wp) THEN
          a_i(ji, jj, jl) = a_i(ib, jb, jl)
          h_i(ji, jj, jl) = h_i(ib, jb, jl)
          h_s(ji, jj, jl) = h_s(ib, jb, jl)
          SELECT CASE (jpbound)
          CASE (0)
            oa_i(ji, jj, jl) = rn_ice_age(jbdy) * a_i(ji, jj, jl)
            a_ip(ji, jj, jl) = 0._wp
            v_ip(ji, jj, jl) = 0._wp
            t_su(ji, jj, jl) = rn_ice_tem(jbdy)
            t_s(ji, jj, :, jl) = rn_ice_tem(jbdy)
            t_i(ji, jj, :, jl) = rn_ice_tem(jbdy)
            s_i(ji, jj, jl) = rn_ice_sal(jbdy)
            sz_i(ji, jj, :, jl) = rn_ice_sal(jbdy)
          CASE (1)
            oa_i(ji, jj, jl) = oa_i(ib, jb, jl)
            a_ip(ji, jj, jl) = a_ip(ib, jb, jl)
            v_ip(ji, jj, jl) = v_ip(ib, jb, jl)
            t_su(ji, jj, jl) = t_su(ib, jb, jl)
            t_s(ji, jj, :, jl) = t_s(ib, jb, :, jl)
            t_i(ji, jj, :, jl) = t_i(ib, jb, :, jl)
            s_i(ji, jj, jl) = s_i(ib, jb, jl)
            sz_i(ji, jj, :, jl) = sz_i(ib, jb, :, jl)
          END SELECT
          IF (nn_icesal == 1) THEN
            s_i(ji, jj, jl) = rn_icesal
            sz_i(ji, jj, :, jl) = rn_icesal
          END IF
          v_i(ji, jj, jl) = h_i(ji, jj, jl) * a_i(ji, jj, jl)
          v_s(ji, jj, jl) = h_s(ji, jj, jl) * a_i(ji, jj, jl)
          sv_i(ji, jj, jl) = MIN(s_i(ji, jj, jl), sss_m(ji, jj)) * v_i(ji, jj, jl)
          DO jk = 1, nlay_s
            e_s(ji, jj, jk, jl) = rhos * (rcpi * (rt0 - t_s(ji, jj, jk, jl)) + rLfus)
            e_s(ji, jj, jk, jl) = e_s(ji, jj, jk, jl) * v_s(ji, jj, jl) * r1_nlay_s
          END DO
          DO jk = 1, nlay_i
            ztmelts = - rTmlt * sz_i(ji, jj, jk, jl)
            t_i(ji, jj, jk, jl) = MIN(t_i(ji, jj, jk, jl), ztmelts + rt0)
            e_i(ji, jj, jk, jl) = rhoi * (rcpi * (ztmelts - (t_i(ji, jj, jk, jl) - rt0)) + rLfus * (1._wp - ztmelts / (t_i(ji, jj, &
&jk, jl) - rt0)) - rcp * ztmelts)
            e_i(ji, jj, jk, jl) = e_i(ji, jj, jk, jl) * v_i(ji, jj, jl) * r1_nlay_i
          END DO
        ELSE
          a_i(ji, jj, jl) = 0._wp
          h_i(ji, jj, jl) = 0._wp
          h_s(ji, jj, jl) = 0._wp
          oa_i(ji, jj, jl) = 0._wp
          a_ip(ji, jj, jl) = 0._wp
          v_ip(ji, jj, jl) = 0._wp
          t_su(ji, jj, jl) = rt0
          t_s(ji, jj, :, jl) = rt0
          t_i(ji, jj, :, jl) = rt0
          IF (nn_icesal == 1) THEN
            s_i(ji, jj, jl) = rn_icesal
            sz_i(ji, jj, :, jl) = rn_icesal
          ELSE
            s_i(ji, jj, jl) = rn_simin
            sz_i(ji, jj, :, jl) = rn_simin
          END IF
          v_i(ji, jj, jl) = 0._wp
          v_s(ji, jj, jl) = 0._wp
          sv_i(ji, jj, jl) = 0._wp
          e_s(ji, jj, :, jl) = 0._wp
          e_i(ji, jj, :, jl) = 0._wp
        END IF
      END DO
    END DO
    CALL lbc_bdy_lnk('bdyice', a_i(:, :, :), 'T', 1., jbdy)
    CALL lbc_bdy_lnk('bdyice', h_i(:, :, :), 'T', 1., jbdy)
    CALL lbc_bdy_lnk('bdyice', h_s(:, :, :), 'T', 1., jbdy)
    CALL lbc_bdy_lnk('bdyice', oa_i(:, :, :), 'T', 1., jbdy)
    CALL lbc_bdy_lnk('bdyice', a_ip(:, :, :), 'T', 1., jbdy)
    CALL lbc_bdy_lnk('bdyice', v_ip(:, :, :), 'T', 1., jbdy)
    CALL lbc_bdy_lnk('bdyice', s_i(:, :, :), 'T', 1., jbdy)
    CALL lbc_bdy_lnk('bdyice', t_su(:, :, :), 'T', 1., jbdy)
    CALL lbc_bdy_lnk('bdyice', v_i(:, :, :), 'T', 1., jbdy)
    CALL lbc_bdy_lnk('bdyice', v_s(:, :, :), 'T', 1., jbdy)
    CALL lbc_bdy_lnk('bdyice', sv_i(:, :, :), 'T', 1., jbdy)
    CALL lbc_bdy_lnk('bdyice', t_s(:, :, :, :), 'T', 1., jbdy)
    CALL lbc_bdy_lnk('bdyice', e_s(:, :, :, :), 'T', 1., jbdy)
    CALL lbc_bdy_lnk('bdyice', t_i(:, :, :, :), 'T', 1., jbdy)
    CALL lbc_bdy_lnk('bdyice', e_i(:, :, :, :), 'T', 1., jbdy)
  END SUBROUTINE bdy_ice_frs
  SUBROUTINE bdy_ice_dyn(cd_type)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    INTEGER :: i_bdy, jgrd
    INTEGER :: ji, jj
    INTEGER :: jbdy
    REAL(KIND = wp) :: zmsk1, zmsk2, zflag
    IF (ln_timing) CALL timing_start('bdy_ice_dyn')
    DO jbdy = 1, nb_bdy
      SELECT CASE (cn_ice(jbdy))
      CASE ('none')
        CYCLE
      CASE ('frs')
        IF (nn_ice_dta(jbdy) == 0) CYCLE
        SELECT CASE (cd_type)
        CASE ('U')
          jgrd = 2
          DO i_bdy = 1, idx_bdy(jbdy) % nblenrim(jgrd)
            ji = idx_bdy(jbdy) % nbi(i_bdy, jgrd)
            jj = idx_bdy(jbdy) % nbj(i_bdy, jgrd)
            zflag = idx_bdy(jbdy) % flagu(i_bdy, jgrd)
            IF (ABS(zflag) == 1.) THEN
              zmsk1 = 1._wp - MAX(0.0_wp, SIGN(1.0_wp, - vt_i(ji + 1, jj)))
              zmsk2 = 1._wp - MAX(0.0_wp, SIGN(1.0_wp, - vt_i(ji, jj)))
              u_ice(ji, jj) = u_ice(ji + 1, jj) * 0.5_wp * ABS(zflag + 1._wp) * zmsk1 + u_ice(ji - 1, jj) * 0.5_wp * ABS(zflag - &
&1._wp) * zmsk2 + u_ice(ji, jj) * (1._wp - MIN(1._wp, zmsk1 + zmsk2))
            ELSE
              u_ice(ji, jj) = 0._wp
            END IF
          END DO
          CALL lbc_bdy_lnk('bdyice', u_ice(:, :), 'U', - 1., jbdy)
        CASE ('V')
          jgrd = 3
          DO i_bdy = 1, idx_bdy(jbdy) % nblenrim(jgrd)
            ji = idx_bdy(jbdy) % nbi(i_bdy, jgrd)
            jj = idx_bdy(jbdy) % nbj(i_bdy, jgrd)
            zflag = idx_bdy(jbdy) % flagv(i_bdy, jgrd)
            IF (ABS(zflag) == 1.) THEN
              zmsk1 = 1._wp - MAX(0.0_wp, SIGN(1.0_wp, - vt_i(ji, jj + 1)))
              zmsk2 = 1._wp - MAX(0.0_wp, SIGN(1.0_wp, - vt_i(ji, jj)))
              v_ice(ji, jj) = v_ice(ji, jj + 1) * 0.5_wp * ABS(zflag + 1._wp) * zmsk1 + v_ice(ji, jj - 1) * 0.5_wp * ABS(zflag - &
&1._wp) * zmsk2 + v_ice(ji, jj) * (1._wp - MIN(1._wp, zmsk1 + zmsk2))
            ELSE
              v_ice(ji, jj) = 0._wp
            END IF
          END DO
          CALL lbc_bdy_lnk('bdyice', v_ice(:, :), 'V', - 1., jbdy)
        END SELECT
      CASE DEFAULT
        CALL ctl_stop('bdy_ice_dyn : unrecognised option for open boundaries for ice fields')
      END SELECT
    END DO
    IF (ln_timing) CALL timing_stop('bdy_ice_dyn')
  END SUBROUTINE bdy_ice_dyn
END MODULE bdyice