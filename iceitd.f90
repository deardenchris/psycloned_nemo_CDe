MODULE iceitd
  USE dom_oce
  USE phycst
  USE ice1D
  USE ice
  USE icevar
  USE icectl
  USE icetab
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  USE prtctl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_itd_init
  PUBLIC :: ice_itd_rem
  PUBLIC :: ice_itd_reb
  INTEGER :: nice_catbnd
  INTEGER, PARAMETER :: np_cathfn = 1
  INTEGER, PARAMETER :: np_catusr = 2
  LOGICAL :: ln_cat_hfn
  REAL(KIND = wp) :: rn_himean
  LOGICAL :: ln_cat_usr
  REAL(KIND = wp), DIMENSION(0 : 100) :: rn_catbnd
  CONTAINS
  SUBROUTINE ice_itd_rem(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jl, jcat
    INTEGER :: ipti
    REAL(KIND = wp) :: zx1, zwk1, zdh0, zetamin, zdamax
    REAL(KIND = wp) :: zx2, zwk2, zda0, zetamax
    REAL(KIND = wp) :: zx3
    REAL(KIND = wp) :: zslope
    INTEGER, DIMENSION(jpij) :: iptidx
    INTEGER, DIMENSION(jpij, jpl - 1) :: jdonor
    REAL(KIND = wp), DIMENSION(jpij, jpl) :: zdhice
    REAL(KIND = wp), DIMENSION(jpij, jpl) :: g0, g1
    REAL(KIND = wp), DIMENSION(jpij, jpl) :: hL, hR
    REAL(KIND = wp), DIMENSION(jpij, jpl - 1) :: zdaice, zdvice
    REAL(KIND = wp), DIMENSION(jpij) :: zhb0, zhb1
    REAL(KIND = wp), DIMENSION(jpij, 0 : jpl) :: zhbnew
    IF (kt == nit000 .AND. lwp) WRITE(numout, FMT = *) '-- ice_itd_rem: remapping ice thickness distribution'
    IF (ln_icediachk) CALL ice_cons_hsm(0, 'iceitd_rem', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    at_i(:, :) = SUM(a_i, dim = 3)
    npti = 0
    nptidx(:) = 0
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (at_i(ji, jj) > epsi10) THEN
          npti = npti + 1
          nptidx(npti) = (jj - 1) * jpi + ji
        END IF
      END DO
    END DO
    IF (npti > 0) THEN
      zdhice(:, :) = 0._wp
      zhbnew(:, :) = 0._wp
      CALL tab_3d_2d(npti, nptidx(1 : npti), h_i_2d(1 : npti, 1 : jpl), h_i)
      CALL tab_3d_2d(npti, nptidx(1 : npti), h_ib_2d(1 : npti, 1 : jpl), h_i_b)
      CALL tab_3d_2d(npti, nptidx(1 : npti), a_i_2d(1 : npti, 1 : jpl), a_i)
      CALL tab_3d_2d(npti, nptidx(1 : npti), a_ib_2d(1 : npti, 1 : jpl), a_i_b)
      DO jl = 1, jpl
        DO ji = 1, npti
          IF (a_i_2d(ji, jl) > epsi10) zdhice(ji, jl) = h_i_2d(ji, jl) - h_ib_2d(ji, jl)
        END DO
      END DO
      DO jl = 1, jpl - 1
        DO ji = 1, npti
          IF (a_ib_2d(ji, jl) > epsi10 .AND. a_ib_2d(ji, jl + 1) > epsi10) THEN
            zslope = (zdhice(ji, jl + 1) - zdhice(ji, jl)) / (h_ib_2d(ji, jl + 1) - h_ib_2d(ji, jl))
            zhbnew(ji, jl) = hi_max(jl) + zdhice(ji, jl) + zslope * (hi_max(jl) - h_ib_2d(ji, jl))
          ELSE IF (a_ib_2d(ji, jl) > epsi10 .AND. a_ib_2d(ji, jl + 1) <= epsi10) THEN
            zhbnew(ji, jl) = hi_max(jl) + zdhice(ji, jl)
          ELSE IF (a_ib_2d(ji, jl) <= epsi10 .AND. a_ib_2d(ji, jl + 1) > epsi10) THEN
            zhbnew(ji, jl) = hi_max(jl) + zdhice(ji, jl + 1)
          ELSE
            zhbnew(ji, jl) = hi_max(jl)
          END IF
          IF (a_i_2d(ji, jl) > epsi10 .AND. h_i_2d(ji, jl) > (zhbnew(ji, jl) - epsi10)) nptidx(ji) = 0
          IF (a_i_2d(ji, jl + 1) > epsi10 .AND. h_i_2d(ji, jl + 1) < (zhbnew(ji, jl) + epsi10)) nptidx(ji) = 0
          IF (zhbnew(ji, jl) < hi_max(jl - 1)) nptidx(ji) = 0
          IF (zhbnew(ji, jl) > hi_max(jl + 1)) nptidx(ji) = 0
        END DO
      END DO
      DO ji = 1, npti
        IF (a_i_2d(ji, jpl) > epsi10) THEN
          zhbnew(ji, jpl) = MAX(hi_max(jpl - 1), 3._wp * h_i_2d(ji, jpl) - 2._wp * zhbnew(ji, jpl - 1))
        ELSE
          zhbnew(ji, jpl) = hi_max(jpl)
        END IF
        IF (h_ib_2d(ji, 1) < (hi_max(0) + epsi10)) nptidx(ji) = 0
        IF (h_ib_2d(ji, 1) > (hi_max(1) - epsi10)) nptidx(ji) = 0
      END DO
      ipti = 0
      iptidx(:) = 0
      DO ji = 1, npti
        IF (nptidx(ji) /= 0) THEN
          ipti = ipti + 1
          iptidx(ipti) = nptidx(ji)
          zhbnew(ipti, :) = zhbnew(ji, :)
        END IF
      END DO
      nptidx(:) = iptidx(:)
      npti = ipti
    END IF
    IF (npti > 0) THEN
      zhb0(:) = hi_max(0)
      zhb1(:) = hi_max(1)
      g0(:, :) = 0._wp
      g1(:, :) = 0._wp
      hl(:, :) = 0._wp
      hr(:, :) = 0._wp
      DO jl = 1, jpl
        CALL tab_2d_1d(npti, nptidx(1 : npti), h_ib_1d(1 : npti), h_i_b(:, :, jl))
        CALL tab_2d_1d(npti, nptidx(1 : npti), h_i_1d(1 : npti), h_i(:, :, jl))
        CALL tab_2d_1d(npti, nptidx(1 : npti), a_i_1d(1 : npti), a_i(:, :, jl))
        CALL tab_2d_1d(npti, nptidx(1 : npti), v_i_1d(1 : npti), v_i(:, :, jl))
        IF (jl == 1) THEN
          CALL itd_glinear(zhb0(1 : npti), zhb1(1 : npti), h_ib_1d(1 : npti), a_i_1d(1 : npti), g0(1 : npti, 1), g1(1 : npti, 1), &
&hL(1 : npti, 1), hR(1 : npti, 1))
          DO ji = 1, npti
            IF (a_i_1d(ji) > epsi10) THEN
              zdh0 = h_i_1d(ji) - h_ib_1d(ji)
              IF (zdh0 < 0.0) THEN
                zdh0 = MIN(- zdh0, hi_max(1))
                zetamax = MIN(zdh0, hR(ji, 1)) - hL(ji, 1)
                IF (zetamax > 0.0) THEN
                  zx1 = zetamax
                  zx2 = 0.5 * zetamax * zetamax
                  zda0 = g1(ji, 1) * zx2 + g0(ji, 1) * zx1
                  zdamax = a_i_1d(ji) * (1.0 - h_i_1d(ji) / h_ib_1d(ji))
                  zda0 = MIN(zda0, zdamax)
                  h_i_1d(ji) = h_i_1d(ji) * a_i_1d(ji) / (a_i_1d(ji) - zda0)
                  a_i_1d(ji) = a_i_1d(ji) - zda0
                  v_i_1d(ji) = a_i_1d(ji) * h_i_1d(ji)
                END IF
              ELSE
                zhbnew(ji, 0) = MIN(zdh0, hi_max(1))
              END IF
            END IF
          END DO
          CALL tab_1d_2d(npti, nptidx(1 : npti), h_i_1d(1 : npti), h_i(:, :, jl))
          CALL tab_1d_2d(npti, nptidx(1 : npti), a_i_1d(1 : npti), a_i(:, :, jl))
          CALL tab_1d_2d(npti, nptidx(1 : npti), v_i_1d(1 : npti), v_i(:, :, jl))
        END IF
        CALL itd_glinear(zhbnew(1 : npti, jl - 1), zhbnew(1 : npti, jl), h_i_1d(1 : npti), a_i_1d(1 : npti), g0(1 : npti, jl), &
&g1(1 : npti, jl), hL(1 : npti, jl), hR(1 : npti, jl))
      END DO
      DO jl = 1, jpl - 1
        DO ji = 1, npti
          IF (zhbnew(ji, jl) > hi_max(jl)) THEN
            zetamin = MAX(hi_max(jl), hL(ji, jl)) - hL(ji, jl)
            zetamax = MIN(zhbnew(ji, jl), hR(ji, jl)) - hL(ji, jl)
            jdonor(ji, jl) = jl
          ELSE
            zetamin = 0.0
            zetamax = MIN(hi_max(jl), hR(ji, jl + 1)) - hL(ji, jl + 1)
            jdonor(ji, jl) = jl + 1
          END IF
          zetamax = MAX(zetamax, zetamin)
          zx1 = zetamax - zetamin
          zwk1 = zetamin * zetamin
          zwk2 = zetamax * zetamax
          zx2 = 0.5 * (zwk2 - zwk1)
          zwk1 = zwk1 * zetamin
          zwk2 = zwk2 * zetamax
          zx3 = 1.0 / 3.0 * (zwk2 - zwk1)
          jcat = jdonor(ji, jl)
          zdaice(ji, jl) = g1(ji, jcat) * zx2 + g0(ji, jcat) * zx1
          zdvice(ji, jl) = g1(ji, jcat) * zx3 + g0(ji, jcat) * zx2 + zdaice(ji, jl) * hL(ji, jcat)
        END DO
      END DO
      CALL itd_shiftice(jdonor(1 : npti, :), zdaice(1 : npti, :), zdvice(1 : npti, :))
      CALL tab_2d_1d(npti, nptidx(1 : npti), h_i_1d(1 : npti), h_i(:, :, 1))
      CALL tab_2d_1d(npti, nptidx(1 : npti), a_i_1d(1 : npti), a_i(:, :, 1))
      CALL tab_2d_1d(npti, nptidx(1 : npti), a_ip_1d(1 : npti), a_ip(:, :, 1))
      DO ji = 1, npti
        IF (a_i_1d(ji) > epsi10 .AND. h_i_1d(ji) < rn_himin) THEN
          a_i_1d(ji) = a_i_1d(ji) * h_i_1d(ji) / rn_himin
          IF (ln_pnd_H12) a_ip_1d(ji) = a_ip_1d(ji) * h_i_1d(ji) / rn_himin
          h_i_1d(ji) = rn_himin
        END IF
      END DO
      CALL tab_1d_2d(npti, nptidx(1 : npti), h_i_1d(1 : npti), h_i(:, :, 1))
      CALL tab_1d_2d(npti, nptidx(1 : npti), a_i_1d(1 : npti), a_i(:, :, 1))
      CALL tab_1d_2d(npti, nptidx(1 : npti), a_ip_1d(1 : npti), a_ip(:, :, 1))
    END IF
    IF (ln_icediachk) CALL ice_cons_hsm(1, 'iceitd_rem', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
  END SUBROUTINE ice_itd_rem
  SUBROUTINE itd_glinear(HbL, Hbr, phice, paice, pg0, pg1, phL, phR)
    REAL(KIND = wp), DIMENSION(:), INTENT(IN) :: HbL, HbR
    REAL(KIND = wp), DIMENSION(:), INTENT(IN) :: phice, paice
    REAL(KIND = wp), DIMENSION(:), INTENT(INOUT) :: pg0, pg1
    REAL(KIND = wp), DIMENSION(:), INTENT(INOUT) :: phL, phR
    INTEGER :: ji
    REAL(KIND = wp) :: z1_3, z2_3
    REAL(KIND = wp) :: zh13
    REAL(KIND = wp) :: zh23
    REAL(KIND = wp) :: zdhr
    REAL(KIND = wp) :: zwk1, zwk2
    z1_3 = 1._wp / 3._wp
    z2_3 = 2._wp / 3._wp
    DO ji = 1, npti
      IF (paice(ji) > epsi10 .AND. phice(ji) > 0._wp) THEN
        phL(ji) = HbL(ji)
        phR(ji) = HbR(ji)
        zh13 = z1_3 * (2._wp * phL(ji) + phR(ji))
        zh23 = z1_3 * (phL(ji) + 2._wp * phR(ji))
        IF (phice(ji) < zh13) THEN
          phr(ji) = 3._wp * phice(ji) - 2._wp * phl(ji)
        ELSE IF (phice(ji) > zh23) THEN
          phl(ji) = 3._wp * phice(ji) - 2._wp * phr(ji)
        END IF
        zdhr = 1._wp / (phR(ji) - phL(ji))
        zwk1 = 6._wp * paice(ji) * zdhr
        zwk2 = (phice(ji) - phL(ji)) * zdhr
        pg0(ji) = zwk1 * (z2_3 - zwk2)
        pg1(ji) = 2._wp * zdhr * zwk1 * (zwk2 - 0.5_wp)
      ELSE
        phL(ji) = 0._wp
        phR(ji) = 0._wp
        pg0(ji) = 0._wp
        pg1(ji) = 0._wp
      END IF
    END DO
  END SUBROUTINE itd_glinear
  SUBROUTINE itd_shiftice(kdonor, pdaice, pdvice)
    INTEGER, DIMENSION(:, :), INTENT(IN) :: kdonor
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pdaice
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pdvice
    INTEGER :: ji, jl, jk
    INTEGER :: jl2, jl1
    REAL(KIND = wp) :: ztrans
    REAL(KIND = wp), DIMENSION(jpij) :: zworka, zworkv
    REAL(KIND = wp), DIMENSION(jpij, jpl) :: zaTsfn
    REAL(KIND = wp), DIMENSION(jpij, nlay_i, jpl) :: ze_i_2d
    REAL(KIND = wp), DIMENSION(jpij, nlay_s, jpl) :: ze_s_2d
    CALL tab_3d_2d(npti, nptidx(1 : npti), h_i_2d(1 : npti, 1 : jpl), h_i)
    CALL tab_3d_2d(npti, nptidx(1 : npti), a_i_2d(1 : npti, 1 : jpl), a_i)
    CALL tab_3d_2d(npti, nptidx(1 : npti), v_i_2d(1 : npti, 1 : jpl), v_i)
    CALL tab_3d_2d(npti, nptidx(1 : npti), v_s_2d(1 : npti, 1 : jpl), v_s)
    CALL tab_3d_2d(npti, nptidx(1 : npti), oa_i_2d(1 : npti, 1 : jpl), oa_i)
    CALL tab_3d_2d(npti, nptidx(1 : npti), sv_i_2d(1 : npti, 1 : jpl), sv_i)
    CALL tab_3d_2d(npti, nptidx(1 : npti), a_ip_2d(1 : npti, 1 : jpl), a_ip)
    CALL tab_3d_2d(npti, nptidx(1 : npti), v_ip_2d(1 : npti, 1 : jpl), v_ip)
    CALL tab_3d_2d(npti, nptidx(1 : npti), t_su_2d(1 : npti, 1 : jpl), t_su)
    DO jl = 1, jpl
      DO jk = 1, nlay_s
        CALL tab_2d_1d(npti, nptidx(1 : npti), ze_s_2d(1 : npti, jk, jl), e_s(:, :, jk, jl))
      END DO
      DO jk = 1, nlay_i
        CALL tab_2d_1d(npti, nptidx(1 : npti), ze_i_2d(1 : npti, jk, jl), e_i(:, :, jk, jl))
      END DO
    END DO
    CALL tab_2d_1d(npti, nptidx(1 : npti), rn_amax_1d(1 : npti), rn_amax_2d)
    DO jl = 1, jpl
      DO ji = 1, npti
        zaTsfn(ji, jl) = a_i_2d(ji, jl) * t_su_2d(ji, jl)
      END DO
    END DO
    DO jl = 1, jpl - 1
      DO ji = 1, npti
        jl1 = kdonor(ji, jl)
        IF (jl1 > 0) THEN
          IF (jl1 == jl) THEN
            jl2 = jl1 + 1
          ELSE
            jl2 = jl
          END IF
          IF (v_i_2d(ji, jl1) >= epsi10) THEN
            zworkv(ji) = pdvice(ji, jl) / v_i_2d(ji, jl1)
          ELSE
            zworkv(ji) = 0._wp
          END IF
          IF (a_i_2d(ji, jl1) >= epsi10) THEN
            zworka(ji) = pdaice(ji, jl) / a_i_2d(ji, jl1)
          ELSE
            zworka(ji) = 0._wp
          END IF
          a_i_2d(ji, jl1) = a_i_2d(ji, jl1) - pdaice(ji, jl)
          a_i_2d(ji, jl2) = a_i_2d(ji, jl2) + pdaice(ji, jl)
          v_i_2d(ji, jl1) = v_i_2d(ji, jl1) - pdvice(ji, jl)
          v_i_2d(ji, jl2) = v_i_2d(ji, jl2) + pdvice(ji, jl)
          ztrans = v_s_2d(ji, jl1) * zworkv(ji)
          v_s_2d(ji, jl1) = v_s_2d(ji, jl1) - ztrans
          v_s_2d(ji, jl2) = v_s_2d(ji, jl2) + ztrans
          ztrans = oa_i_2d(ji, jl1) * zworka(ji)
          oa_i_2d(ji, jl1) = oa_i_2d(ji, jl1) - ztrans
          oa_i_2d(ji, jl2) = oa_i_2d(ji, jl2) + ztrans
          ztrans = sv_i_2d(ji, jl1) * zworkv(ji)
          sv_i_2d(ji, jl1) = sv_i_2d(ji, jl1) - ztrans
          sv_i_2d(ji, jl2) = sv_i_2d(ji, jl2) + ztrans
          ztrans = zaTsfn(ji, jl1) * zworka(ji)
          zaTsfn(ji, jl1) = zaTsfn(ji, jl1) - ztrans
          zaTsfn(ji, jl2) = zaTsfn(ji, jl2) + ztrans
          IF (ln_pnd_H12) THEN
            ztrans = a_ip_2d(ji, jl1) * zworka(ji)
            a_ip_2d(ji, jl1) = a_ip_2d(ji, jl1) - ztrans
            a_ip_2d(ji, jl2) = a_ip_2d(ji, jl2) + ztrans
            ztrans = v_ip_2d(ji, jl1) * zworka(ji)
            v_ip_2d(ji, jl1) = v_ip_2d(ji, jl1) - ztrans
            v_ip_2d(ji, jl2) = v_ip_2d(ji, jl2) + ztrans
          END IF
        END IF
      END DO
      DO jk = 1, nlay_s
        DO ji = 1, npti
          jl1 = kdonor(ji, jl)
          IF (jl1 > 0) THEN
            IF (jl1 == jl) THEN
              jl2 = jl + 1
            ELSE
              jl2 = jl
            END IF
            ztrans = ze_s_2d(ji, jk, jl1) * zworkv(ji)
            ze_s_2d(ji, jk, jl1) = ze_s_2d(ji, jk, jl1) - ztrans
            ze_s_2d(ji, jk, jl2) = ze_s_2d(ji, jk, jl2) + ztrans
          END IF
        END DO
      END DO
      DO jk = 1, nlay_i
        DO ji = 1, npti
          jl1 = kdonor(ji, jl)
          IF (jl1 > 0) THEN
            IF (jl1 == jl) THEN
              jl2 = jl + 1
            ELSE
              jl2 = jl
            END IF
            ztrans = ze_i_2d(ji, jk, jl1) * zworkv(ji)
            ze_i_2d(ji, jk, jl1) = ze_i_2d(ji, jk, jl1) - ztrans
            ze_i_2d(ji, jk, jl2) = ze_i_2d(ji, jk, jl2) + ztrans
          END IF
        END DO
      END DO
    END DO
    CALL ice_var_roundoff(a_i_2d, v_i_2d, v_s_2d, sv_i_2d, oa_i_2d, a_ip_2d, v_ip_2d, ze_s_2d, ze_i_2d)
    zworka(1 : npti) = SUM(a_i_2d(1 : npti, :), dim = 2)
    DO jl = 1, jpl
      WHERE (zworka(1 : npti) > rn_amax_1d(1 : npti)) a_i_2d(1 : npti, jl) = a_i_2d(1 : npti, jl) * rn_amax_1d(1 : npti) / &
&zworka(1 : npti)
    END DO
    WHERE (a_i_2d(1 : npti, :) >= epsi20)
      h_i_2d(1 : npti, :) = v_i_2d(1 : npti, :) / a_i_2d(1 : npti, :)
      t_su_2d(1 : npti, :) = zaTsfn(1 : npti, :) / a_i_2d(1 : npti, :)
    ELSEWHERE
      h_i_2d(1 : npti, :) = 0._wp
      t_su_2d(1 : npti, :) = rt0
    END WHERE
    CALL tab_2d_3d(npti, nptidx(1 : npti), h_i_2d(1 : npti, 1 : jpl), h_i)
    CALL tab_2d_3d(npti, nptidx(1 : npti), a_i_2d(1 : npti, 1 : jpl), a_i)
    CALL tab_2d_3d(npti, nptidx(1 : npti), v_i_2d(1 : npti, 1 : jpl), v_i)
    CALL tab_2d_3d(npti, nptidx(1 : npti), v_s_2d(1 : npti, 1 : jpl), v_s)
    CALL tab_2d_3d(npti, nptidx(1 : npti), oa_i_2d(1 : npti, 1 : jpl), oa_i)
    CALL tab_2d_3d(npti, nptidx(1 : npti), sv_i_2d(1 : npti, 1 : jpl), sv_i)
    CALL tab_2d_3d(npti, nptidx(1 : npti), a_ip_2d(1 : npti, 1 : jpl), a_ip)
    CALL tab_2d_3d(npti, nptidx(1 : npti), v_ip_2d(1 : npti, 1 : jpl), v_ip)
    CALL tab_2d_3d(npti, nptidx(1 : npti), t_su_2d(1 : npti, 1 : jpl), t_su)
    DO jl = 1, jpl
      DO jk = 1, nlay_s
        CALL tab_1d_2d(npti, nptidx(1 : npti), ze_s_2d(1 : npti, jk, jl), e_s(:, :, jk, jl))
      END DO
      DO jk = 1, nlay_i
        CALL tab_1d_2d(npti, nptidx(1 : npti), ze_i_2d(1 : npti, jk, jl), e_i(:, :, jk, jl))
      END DO
    END DO
  END SUBROUTINE itd_shiftice
  SUBROUTINE ice_itd_reb(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jl
    INTEGER, DIMENSION(jpij, jpl - 1) :: jdonor
    REAL(KIND = wp), DIMENSION(jpij, jpl - 1) :: zdaice, zdvice
    IF (kt == nit000 .AND. lwp) WRITE(numout, FMT = *) '-- ice_itd_reb: rebining ice thickness distribution'
    IF (ln_icediachk) CALL ice_cons_hsm(0, 'iceitd_reb', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    jdonor(:, :) = 0
    zdaice(:, :) = 0._wp
    zdvice(:, :) = 0._wp
    DO jl = 1, jpl - 1
      npti = 0
      nptidx(:) = 0
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (a_i(ji, jj, jl) > 0._wp .AND. v_i(ji, jj, jl) > (a_i(ji, jj, jl) * hi_max(jl))) THEN
            npti = npti + 1
            nptidx(npti) = (jj - 1) * jpi + ji
          END IF
        END DO
      END DO
      CALL tab_2d_1d(npti, nptidx(1 : npti), a_i_1d(1 : npti), a_i(:, :, jl))
      CALL tab_2d_1d(npti, nptidx(1 : npti), v_i_1d(1 : npti), v_i(:, :, jl))
      DO ji = 1, npti
        jdonor(ji, jl) = jl
        zdaice(ji, jl) = a_i_1d(ji) * 0.5_wp
        zdvice(ji, jl) = v_i_1d(ji) - zdaice(ji, jl) * (hi_max(jl) + hi_max(jl - 1)) * 0.5_wp
      END DO
      IF (npti > 0) THEN
        CALL itd_shiftice(jdonor(1 : npti, :), zdaice(1 : npti, :), zdvice(1 : npti, :))
        jdonor(1 : npti, jl) = 0
        zdaice(1 : npti, jl) = 0._wp
        zdvice(1 : npti, jl) = 0._wp
      END IF
    END DO
    DO jl = jpl - 1, 1, - 1
      npti = 0
      nptidx(:) = 0
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (a_i(ji, jj, jl + 1) > 0._wp .AND. v_i(ji, jj, jl + 1) <= (a_i(ji, jj, jl + 1) * hi_max(jl))) THEN
            npti = npti + 1
            nptidx(npti) = (jj - 1) * jpi + ji
          END IF
        END DO
      END DO
      CALL tab_2d_1d(npti, nptidx(1 : npti), a_i_1d(1 : npti), a_i(:, :, jl + 1))
      CALL tab_2d_1d(npti, nptidx(1 : npti), v_i_1d(1 : npti), v_i(:, :, jl + 1))
      DO ji = 1, npti
        jdonor(ji, jl) = jl + 1
        zdaice(ji, jl) = a_i_1d(ji)
        zdvice(ji, jl) = v_i_1d(ji)
      END DO
      IF (npti > 0) THEN
        CALL itd_shiftice(jdonor(1 : npti, :), zdaice(1 : npti, :), zdvice(1 : npti, :))
        jdonor(1 : npti, jl) = 0
        zdaice(1 : npti, jl) = 0._wp
        zdvice(1 : npti, jl) = 0._wp
      END IF
    END DO
    IF (ln_icediachk) CALL ice_cons_hsm(1, 'iceitd_reb', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
  END SUBROUTINE ice_itd_reb
  SUBROUTINE ice_itd_init
    INTEGER :: jl
    INTEGER :: ios, ioptio
    REAL(KIND = wp) :: zhmax, znum, zden, zalpha
    NAMELIST /namitd/ ln_cat_hfn, rn_himean, ln_cat_usr, rn_catbnd, rn_himin
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, namitd, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namitd in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, namitd, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namitd in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, namitd)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'ice_itd_init: Initialization of ice cat distribution '
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namitd: '
      WRITE(numout, FMT = *) '      Ice categories are defined by a function of rn_himean**(-0.05)    ln_cat_hfn = ', ln_cat_hfn
      WRITE(numout, FMT = *) '         mean ice thickness in the domain                               rn_himean  = ', rn_himean
      WRITE(numout, FMT = *) '      Ice categories are defined by rn_catbnd                           ln_cat_usr = ', ln_cat_usr
      WRITE(numout, FMT = *) '      minimum ice thickness                                             rn_himin   = ', rn_himin
    END IF
    ioptio = 0
    IF (ln_cat_hfn) THEN
      ioptio = ioptio + 1
      nice_catbnd = np_cathfn
    END IF
    IF (ln_cat_usr) THEN
      ioptio = ioptio + 1
      nice_catbnd = np_catusr
    END IF
    IF (ioptio /= 1) CALL ctl_stop('ice_itd_init: choose one and only one ice categories boundaries')
    SELECT CASE (nice_catbnd)
    CASE (np_cathfn)
      zalpha = 0.05_wp
      zhmax = 3._wp * rn_himean
      hi_max(0) = 0._wp
      DO jl = 1, jpl
        znum = jpl * (zhmax + 1) ** zalpha
        zden = REAL(jpl - jl, wp) * (zhmax + 1._wp) ** zalpha + REAL(jl, wp)
        hi_max(jl) = (znum / zden) ** (1. / zalpha) - 1
      END DO
    CASE (np_catusr)
      DO jl = 0, jpl
        hi_max(jl) = rn_catbnd(jl)
      END DO
    END SELECT
    DO jl = 1, jpl
      hi_mean(jl) = (hi_max(jl) + hi_max(jl - 1)) * 0.5_wp
    END DO
    hi_max(jpl) = 99._wp
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) '   ===>>>   resulting thickness category boundaries :'
    IF (lwp) WRITE(numout, FMT = *) '            hi_max(:)= ', hi_max(0 : jpl)
    IF (hi_max(1) < rn_himin) CALL ctl_stop('ice_itd_init: the upper bound of the 1st category must be bigger than rn_himin')
  END SUBROUTINE ice_itd_init
END MODULE iceitd