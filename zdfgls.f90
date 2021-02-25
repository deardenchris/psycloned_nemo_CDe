MODULE zdfgls
  USE oce
  USE dom_oce
  USE domvvl
  USE zdfdrg, ONLY: r_z0_top, r_z0_bot
  USE zdfdrg, ONLY: rCdU_top, rCdU_bot
  USE sbc_oce
  USE phycst
  USE zdfmxl
  USE sbcwave, ONLY: hsw
  USE lbclnk
  USE lib_mpp
  USE prtctl
  USE in_out_manager
  USE iom
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: zdf_gls
  PUBLIC :: zdf_gls_init
  PUBLIC :: gls_rst
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: hmxl_n
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: zwall
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ustar2_surf
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ustar2_top
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ustar2_bot
  LOGICAL :: ln_length_lim
  LOGICAL :: ln_sigpsi
  INTEGER :: nn_bc_surf
  INTEGER :: nn_bc_bot
  INTEGER :: nn_z0_met
  INTEGER :: nn_stab_func
  INTEGER :: nn_clos
  REAL(KIND = wp) :: rn_clim_galp
  REAL(KIND = wp) :: rn_epsmin
  REAL(KIND = wp) :: rn_emin
  REAL(KIND = wp) :: rn_charn
  REAL(KIND = wp) :: rn_crban
  REAL(KIND = wp) :: rn_hsro
  REAL(KIND = wp) :: rn_frac_hs
  REAL(KIND = wp) :: rcm_sf = 0.73_wp
  REAL(KIND = wp) :: ra_sf = - 2.0_wp
  REAL(KIND = wp) :: rl_sf = 0.2_wp
  REAL(KIND = wp) :: rghmin = - 0.28_wp
  REAL(KIND = wp) :: rgh0 = 0.0329_wp
  REAL(KIND = wp) :: rghcri = 0.03_wp
  REAL(KIND = wp) :: ra1 = 0.92_wp
  REAL(KIND = wp) :: ra2 = 0.74_wp
  REAL(KIND = wp) :: rb1 = 16.60_wp
  REAL(KIND = wp) :: rb2 = 10.10_wp
  REAL(KIND = wp) :: re2 = 1.33_wp
  REAL(KIND = wp) :: rl1 = 0.107_wp
  REAL(KIND = wp) :: rl2 = 0.0032_wp
  REAL(KIND = wp) :: rl3 = 0.0864_wp
  REAL(KIND = wp) :: rl4 = 0.12_wp
  REAL(KIND = wp) :: rl5 = 11.9_wp
  REAL(KIND = wp) :: rl6 = 0.4_wp
  REAL(KIND = wp) :: rl7 = 0.0_wp
  REAL(KIND = wp) :: rl8 = 0.48_wp
  REAL(KIND = wp) :: rm1 = 0.127_wp
  REAL(KIND = wp) :: rm2 = 0.00336_wp
  REAL(KIND = wp) :: rm3 = 0.0906_wp
  REAL(KIND = wp) :: rm4 = 0.101_wp
  REAL(KIND = wp) :: rm5 = 11.2_wp
  REAL(KIND = wp) :: rm6 = 0.4_wp
  REAL(KIND = wp) :: rm7 = 0.0_wp
  REAL(KIND = wp) :: rm8 = 0.318_wp
  REAL(KIND = wp) :: rtrans = 0.1_wp
  REAL(KIND = wp) :: rc02, rc02r, rc03, rc04
  REAL(KIND = wp) :: rsbc_tke1, rsbc_tke2, rfact_tke
  REAL(KIND = wp) :: rsbc_psi1, rsbc_psi2, rfact_psi
  REAL(KIND = wp) :: rsbc_zs1, rsbc_zs2
  REAL(KIND = wp) :: rc0, rc2, rc3, rf6, rcff, rc_diff
  REAL(KIND = wp) :: rs0, rs1, rs2, rs4, rs5, rs6
  REAL(KIND = wp) :: rd0, rd1, rd2, rd3, rd4, rd5
  REAL(KIND = wp) :: rsc_tke, rsc_psi, rpsi1, rpsi2, rpsi3, rsc_psi0
  REAL(KIND = wp) :: rpsi3m, rpsi3p, rpp, rmm, rnn
  REAL(KIND = wp) :: r2_3 = 2._wp / 3._wp
  CONTAINS
  INTEGER FUNCTION zdf_gls_alloc()
    ALLOCATE(hmxl_n(jpi, jpj, jpk), ustar2_surf(jpi, jpj), zwall(jpi, jpj, jpk), ustar2_top(jpi, jpj), ustar2_bot(jpi, jpj), STAT &
&= zdf_gls_alloc)
    CALL mpp_sum('zdfgls', zdf_gls_alloc)
    IF (zdf_gls_alloc /= 0) CALL ctl_stop('STOP', 'zdf_gls_alloc: failed to allocate arrays')
  END FUNCTION zdf_gls_alloc
  SUBROUTINE zdf_gls(kt, p_sh2, p_avm, p_avt)
    USE zdf_oce, ONLY: en, avtb, avmb
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: p_sh2
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: p_avm, p_avt
    INTEGER :: ji, jj, jk
    INTEGER :: ibot, ibotm1
    INTEGER :: itop, itopp1
    REAL(KIND = wp) :: zesh2, zsigpsi, zcoef, zex1, zex2
    REAL(KIND = wp) :: ztx2, zty2, zup, zdown, zcof, zdir
    REAL(KIND = wp) :: zratio, zrn2, zflxb, sh, z_en
    REAL(KIND = wp) :: prod, buoy, diss, zdiss, sm
    REAL(KIND = wp) :: gh, gm, shr, dif, zsqen, zavt, zavm
    REAL(KIND = wp) :: zmsku, zmskv
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdep
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zkar
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zflxs
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zhsro
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: eb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: hmxl_b
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: eps
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwall_psi
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: psi
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zd_lw, zd_up, zdiag
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zstt, zstm
    ustar2_surf(:, :) = 0._wp
    psi(:, :, :) = 0._wp
    ustar2_top(:, :) = 0._wp
    zwall_psi(:, :, :) = 0._wp
    ustar2_bot(:, :) = 0._wp
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        ustar2_surf(ji, jj) = r1_rau0 * taum(ji, jj) * tmask(ji, jj, 1)
        zmsku = (2._wp - umask(ji - 1, jj, mbkt(ji, jj)) * umask(ji, jj, mbkt(ji, jj)))
        zmskv = (2._wp - vmask(ji, jj - 1, mbkt(ji, jj)) * vmask(ji, jj, mbkt(ji, jj)))
        ustar2_bot(ji, jj) = - rCdU_bot(ji, jj) * SQRT((zmsku * (ub(ji, jj, mbkt(ji, jj)) + ub(ji - 1, jj, mbkt(ji, jj)))) ** 2 + &
&(zmskv * (vb(ji, jj, mbkt(ji, jj)) + vb(ji, jj - 1, mbkt(ji, jj)))) ** 2)
      END DO
    END DO
    IF (ln_isfcav) THEN
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zmsku = (2. - umask(ji - 1, jj, mikt(ji, jj)) * umask(ji, jj, mikt(ji, jj)))
          zmskv = (2. - vmask(ji, jj - 1, mikt(ji, jj)) * vmask(ji, jj, mikt(ji, jj)))
          ustar2_top(ji, jj) = - rCdU_top(ji, jj) * SQRT((zmsku * (ub(ji, jj, mikt(ji, jj)) + ub(ji - 1, jj, mikt(ji, jj)))) ** 2 &
&+ (zmskv * (vb(ji, jj, mikt(ji, jj)) + vb(ji, jj - 1, mikt(ji, jj)))) ** 2)
        END DO
      END DO
    END IF
    SELECT CASE (nn_z0_met)
    CASE (0)
      zhsro(:, :) = rn_hsro
    CASE (1)
      zhsro(:, :) = MAX(rsbc_zs1 * ustar2_surf(:, :), rn_hsro)
    CASE (2)
      zdep(:, :) = 30. * TANH(2. * 0.3 / (28. * SQRT(MAX(ustar2_surf(:, :), rsmall))))
      zhsro(:, :) = MAX(rsbc_zs2 * ustar2_surf(:, :) * zdep(:, :) ** 1.5, rn_hsro)
    CASE (3)
      zhsro(:, :) = rn_frac_hs * hsw(:, :)
    END SELECT
    DO jk = 2, jpkm1
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          eps(ji, jj, jk) = rc03 * en(ji, jj, jk) * SQRT(en(ji, jj, jk)) / hmxl_n(ji, jj, jk)
        END DO
      END DO
    END DO
    eb(:, :, :) = en(:, :, :)
    hmxl_b(:, :, :) = hmxl_n(:, :, :)
    IF (nn_clos == 0) THEN
      DO jk = 2, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zup = hmxl_n(ji, jj, jk) * gdepw_n(ji, jj, mbkt(ji, jj) + 1)
            zdown = vkarmn * gdepw_n(ji, jj, jk) * (- gdepw_n(ji, jj, jk) + gdepw_n(ji, jj, mbkt(ji, jj) + 1))
            zcoef = (zup / MAX(zdown, rsmall))
            zwall(ji, jj, jk) = (1._wp + re2 * zcoef * zcoef) * tmask(ji, jj, jk)
          END DO
        END DO
      END DO
    END IF
    DO jk = 2, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          buoy = - p_avt(ji, jj, jk) * rn2(ji, jj, jk)
          diss = eps(ji, jj, jk)
          zdir = 0.5_wp + SIGN(0.5_wp, p_sh2(ji, jj, jk) + buoy)
          zesh2 = zdir * (p_sh2(ji, jj, jk) + buoy) + (1._wp - zdir) * p_sh2(ji, jj, jk)
          zdiss = zdir * (diss / en(ji, jj, jk)) + (1._wp - zdir) * (diss - buoy) / en(ji, jj, jk)
          IF (ln_sigpsi) THEN
            zsigpsi = MIN(1._wp, zesh2 / eps(ji, jj, jk))
            zwall_psi(ji, jj, jk) = rsc_psi / (zsigpsi * rsc_psi + (1._wp - zsigpsi) * rsc_psi0 / MAX(zwall(ji, jj, jk), 1._wp))
          ELSE
            zwall_psi(ji, jj, jk) = 1._wp
          END IF
          zcof = rfact_tke * tmask(ji, jj, jk)
          zd_lw(ji, jj, jk) = zcof * (p_avm(ji, jj, jk) + p_avm(ji, jj, jk - 1)) / (e3t_n(ji, jj, jk - 1) * e3w_n(ji, jj, jk))
          zd_up(ji, jj, jk) = zcof * (p_avm(ji, jj, jk + 1) + p_avm(ji, jj, jk)) / (e3t_n(ji, jj, jk) * e3w_n(ji, jj, jk))
          zdiag(ji, jj, jk) = 1._wp - zd_lw(ji, jj, jk) - zd_up(ji, jj, jk) + rdt * zdiss * wmask(ji, jj, jk)
          en(ji, jj, jk) = en(ji, jj, jk) + rdt * zesh2 * wmask(ji, jj, jk)
        END DO
      END DO
    END DO
    zdiag(:, :, jpk) = 1._wp
    zwall_psi(:, :, 1) = zwall_psi(:, :, 2)
    zwall_psi(:, :, jpk) = 1._wp
    SELECT CASE (nn_bc_surf)
    CASE (0)
      en(:, :, 1) = MAX(rn_emin, rc02r * ustar2_surf(:, :) * (1._wp + rsbc_tke1) ** r2_3)
      zd_lw(:, :, 1) = en(:, :, 1)
      zd_up(:, :, 1) = 0._wp
      zdiag(:, :, 1) = 1._wp
      en(:, :, 2) = MAX(rc02r * ustar2_surf(:, :) * (1._wp + rsbc_tke1 * ((zhsro(:, :) + gdepw_n(:, :, 2)) / zhsro(:, :)) ** &
&(1.5_wp * ra_sf)) ** (2._wp / 3._wp), rn_emin)
      zd_lw(:, :, 2) = 0._wp
      zd_up(:, :, 2) = 0._wp
      zdiag(:, :, 2) = 1._wp
    CASE (1)
      en(:, :, 1) = MAX(rc02r * ustar2_surf(:, :) * (1._wp + rsbc_tke1) ** r2_3, rn_emin)
      zd_lw(:, :, 1) = en(:, :, 1)
      zd_up(:, :, 1) = 0._wp
      zdiag(:, :, 1) = 1._wp
      zdiag(:, :, 2) = zdiag(:, :, 2) + zd_lw(:, :, 2)
      zd_lw(:, :, 2) = 0._wp
      zkar(:, :) = (rl_sf + (vkarmn - rl_sf) * (1. - EXP(- rtrans * gdept_n(:, :, 1) / zhsro(:, :))))
      zflxs(:, :) = rsbc_tke2 * ustar2_surf(:, :) ** 1.5_wp * zkar(:, :) * ((zhsro(:, :) + gdept_n(:, :, 1)) / zhsro(:, :)) ** &
&(1.5_wp * ra_sf)
      en(:, :, 2) = en(:, :, 2) + zflxs(:, :) / e3w_n(:, :, 2)
    END SELECT
    SELECT CASE (nn_bc_bot)
    CASE (0)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ibot = mbkt(ji, jj) + 1
          ibotm1 = mbkt(ji, jj)
          z_en = MAX(rc02r * ustar2_bot(ji, jj), rn_emin)
          zd_lw(ji, jj, ibot) = 0._wp
          zd_lw(ji, jj, ibotm1) = 0._wp
          zd_up(ji, jj, ibot) = 0._wp
          zd_up(ji, jj, ibotm1) = 0._wp
          zdiag(ji, jj, ibot) = 1._wp
          zdiag(ji, jj, ibotm1) = 1._wp
          en(ji, jj, ibot) = z_en
          en(ji, jj, ibotm1) = z_en
        END DO
      END DO
      IF (ln_isfcav) THEN
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            itop = mikt(ji, jj)
            itopp1 = mikt(ji, jj) + 1
            z_en = MAX(rc02r * ustar2_top(ji, jj), rn_emin) * (1._wp - tmask(ji, jj, 1))
            zd_lw(ji, jj, itop) = 0._wp
            zd_lw(ji, jj, itopp1) = 0._wp
            zd_up(ji, jj, itop) = 0._wp
            zd_up(ji, jj, itopp1) = 0._wp
            zdiag(ji, jj, itop) = 1._wp
            zdiag(ji, jj, itopp1) = 1._wp
            en(ji, jj, itop) = z_en
            en(ji, jj, itopp1) = z_en
          END DO
        END DO
      END IF
    CASE (1)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ibot = mbkt(ji, jj) + 1
          ibotm1 = mbkt(ji, jj)
          z_en = MAX(rc02r * ustar2_bot(ji, jj), rn_emin)
          zd_lw(ji, jj, ibot) = 0._wp
          zdiag(ji, jj, ibot) = 1._wp
          zdiag(ji, jj, ibotm1) = zdiag(ji, jj, ibotm1) + zd_up(ji, jj, ibotm1)
          zd_up(ji, jj, ibot) = 0._wp
          zd_up(ji, jj, ibotm1) = 0._wp
        END DO
      END DO
      IF (ln_isfcav) THEN
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            itop = mikt(ji, jj)
            itopp1 = mikt(ji, jj) + 1
            z_en = MAX(rc02r * ustar2_top(ji, jj), rn_emin) * (1._wp - tmask(ji, jj, 1))
            zd_lw(ji, jj, itop) = 0._wp
            zdiag(ji, jj, itop) = 1._wp
            zdiag(ji, jj, itopp1) = zdiag(ji, jj, itopp1) + zd_up(ji, jj, itopp1)
            zd_up(ji, jj, itop) = 0._wp
            zd_up(ji, jj, itopp1) = 0._wp
          END DO
        END DO
      END IF
    END SELECT
    DO jk = 2, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zdiag(ji, jj, jk) = zdiag(ji, jj, jk) - zd_lw(ji, jj, jk) * zd_up(ji, jj, jk - 1) / zdiag(ji, jj, jk - 1)
        END DO
      END DO
    END DO
    DO jk = 2, jpk
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zd_lw(ji, jj, jk) = en(ji, jj, jk) - zd_lw(ji, jj, jk) / zdiag(ji, jj, jk - 1) * zd_lw(ji, jj, jk - 1)
        END DO
      END DO
    END DO
    DO jk = jpk - 1, 2, - 1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          en(ji, jj, jk) = (zd_lw(ji, jj, jk) - zd_up(ji, jj, jk) * en(ji, jj, jk + 1)) / zdiag(ji, jj, jk)
        END DO
      END DO
    END DO
    en(:, :, :) = MAX(en(:, :, :), rn_emin)
    SELECT CASE (nn_clos)
    CASE (0)
      DO jk = 2, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            psi(ji, jj, jk) = eb(ji, jj, jk) * hmxl_b(ji, jj, jk)
          END DO
        END DO
      END DO
    CASE (1)
      DO jk = 2, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            psi(ji, jj, jk) = eps(ji, jj, jk)
          END DO
        END DO
      END DO
    CASE (2)
      DO jk = 2, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            psi(ji, jj, jk) = SQRT(eb(ji, jj, jk)) / (rc0 * hmxl_b(ji, jj, jk))
          END DO
        END DO
      END DO
    CASE (3)
      DO jk = 2, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            psi(ji, jj, jk) = rc02 * eb(ji, jj, jk) * hmxl_b(ji, jj, jk) ** rnn
          END DO
        END DO
      END DO
    END SELECT
    DO jk = 2, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zratio = psi(ji, jj, jk) / eb(ji, jj, jk)
          zdir = 0.5_wp + SIGN(0.5_wp, rn2(ji, jj, jk))
          rpsi3 = zdir * rpsi3m + (1._wp - zdir) * rpsi3p
          prod = rpsi1 * zratio * p_sh2(ji, jj, jk)
          buoy = rpsi3 * zratio * (- p_avt(ji, jj, jk) * rn2(ji, jj, jk))
          diss = rpsi2 * zratio * zwall(ji, jj, jk) * eps(ji, jj, jk)
          zdir = 0.5_wp + SIGN(0.5_wp, prod + buoy)
          zesh2 = zdir * (prod + buoy) + (1._wp - zdir) * prod
          zdiss = zdir * (diss / psi(ji, jj, jk)) + (1._wp - zdir) * (diss - buoy) / psi(ji, jj, jk)
          zcof = rfact_psi * zwall_psi(ji, jj, jk) * tmask(ji, jj, jk)
          zd_lw(ji, jj, jk) = zcof * (p_avm(ji, jj, jk) + p_avm(ji, jj, jk - 1)) / (e3t_n(ji, jj, jk - 1) * e3w_n(ji, jj, jk))
          zd_up(ji, jj, jk) = zcof * (p_avm(ji, jj, jk + 1) + p_avm(ji, jj, jk)) / (e3t_n(ji, jj, jk) * e3w_n(ji, jj, jk))
          zdiag(ji, jj, jk) = 1._wp - zd_lw(ji, jj, jk) - zd_up(ji, jj, jk) + rdt * zdiss * wmask(ji, jj, jk)
          psi(ji, jj, jk) = psi(ji, jj, jk) + rdt * zesh2 * wmask(ji, jj, jk)
        END DO
      END DO
    END DO
    zdiag(:, :, jpk) = 1._wp
    SELECT CASE (nn_bc_surf)
    CASE (0)
      zdep(:, :) = zhsro(:, :) * rl_sf
      psi(:, :, 1) = rc0 ** rpp * en(:, :, 1) ** rmm * zdep(:, :) ** rnn * tmask(:, :, 1)
      zd_lw(:, :, 1) = psi(:, :, 1)
      zd_up(:, :, 1) = 0._wp
      zdiag(:, :, 1) = 1._wp
      zkar(:, :) = (rl_sf + (vkarmn - rl_sf) * (1._wp - EXP(- rtrans * gdepw_n(:, :, 2) / zhsro(:, :))))
      zdep(:, :) = (zhsro(:, :) + gdepw_n(:, :, 2)) * zkar(:, :)
      psi(:, :, 2) = rc0 ** rpp * en(:, :, 2) ** rmm * zdep(:, :) ** rnn * tmask(:, :, 1)
      zd_lw(:, :, 2) = 0._wp
      zd_up(:, :, 2) = 0._wp
      zdiag(:, :, 2) = 1._wp
    CASE (1)
      zdep(:, :) = zhsro(:, :) * rl_sf
      psi(:, :, 1) = rc0 ** rpp * en(:, :, 1) ** rmm * zdep(:, :) ** rnn * tmask(:, :, 1)
      zd_lw(:, :, 1) = psi(:, :, 1)
      zd_up(:, :, 1) = 0._wp
      zdiag(:, :, 1) = 1._wp
      zdiag(:, :, 2) = zdiag(:, :, 2) + zd_lw(:, :, 2)
      zd_lw(:, :, 2) = 0._wp
      zkar(:, :) = rl_sf + (vkarmn - rl_sf) * (1._wp - EXP(- rtrans * gdept_n(:, :, 1) / zhsro(:, :)))
      zdep(:, :) = ((zhsro(:, :) + gdept_n(:, :, 1)) / zhsro(:, :)) ** (rmm * ra_sf)
      zflxs(:, :) = (rnn + rsbc_tke1 * (rnn + rmm * ra_sf) * zdep(:, :)) * (1._wp + rsbc_tke1 * zdep(:, :)) ** (2._wp * rmm / &
&3._wp - 1_wp)
      zdep(:, :) = rsbc_psi1 * (zwall_psi(:, :, 1) * p_avm(:, :, 1) + zwall_psi(:, :, 2) * p_avm(:, :, 2)) * ustar2_surf(:, :) ** &
&rmm * zkar(:, :) ** rnn * (zhsro(:, :) + gdept_n(:, :, 1)) ** (rnn - 1.)
      zflxs(:, :) = zdep(:, :) * zflxs(:, :)
      psi(:, :, 2) = psi(:, :, 2) + zflxs(:, :) / e3w_n(:, :, 2)
    END SELECT
    SELECT CASE (nn_bc_bot)
    CASE (0)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ibot = mbkt(ji, jj) + 1
          ibotm1 = mbkt(ji, jj)
          zdep(ji, jj) = vkarmn * r_z0_bot
          psi(ji, jj, ibot) = rc0 ** rpp * en(ji, jj, ibot) ** rmm * zdep(ji, jj) ** rnn
          zd_lw(ji, jj, ibot) = 0._wp
          zd_up(ji, jj, ibot) = 0._wp
          zdiag(ji, jj, ibot) = 1._wp
          zdep(ji, jj) = vkarmn * (r_z0_bot + e3t_n(ji, jj, ibotm1))
          psi(ji, jj, ibotm1) = rc0 ** rpp * en(ji, jj, ibot) ** rmm * zdep(ji, jj) ** rnn
          zd_lw(ji, jj, ibotm1) = 0._wp
          zd_up(ji, jj, ibotm1) = 0._wp
          zdiag(ji, jj, ibotm1) = 1._wp
        END DO
      END DO
    CASE (1)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ibot = mbkt(ji, jj) + 1
          ibotm1 = mbkt(ji, jj)
          zdep(ji, jj) = vkarmn * r_z0_bot
          psi(ji, jj, ibot) = rc0 ** rpp * en(ji, jj, ibot) ** rmm * zdep(ji, jj) ** rnn
          zd_lw(ji, jj, ibot) = 0._wp
          zd_up(ji, jj, ibot) = 0._wp
          zdiag(ji, jj, ibot) = 1._wp
          zdiag(ji, jj, ibotm1) = zdiag(ji, jj, ibotm1) + zd_up(ji, jj, ibotm1)
          zd_up(ji, jj, ibotm1) = 0.
          zdep(ji, jj) = r_z0_bot + 0.5_wp * e3t_n(ji, jj, ibotm1)
          zflxb = rsbc_psi2 * (p_avm(ji, jj, ibot) + p_avm(ji, jj, ibotm1)) * (0.5_wp * (en(ji, jj, ibot) + en(ji, jj, ibotm1))) &
&** rmm * zdep(ji, jj) ** (rnn - 1._wp)
          psi(ji, jj, ibotm1) = psi(ji, jj, ibotm1) + zflxb / e3w_n(ji, jj, ibotm1)
        END DO
      END DO
    END SELECT
    DO jk = 2, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zdiag(ji, jj, jk) = zdiag(ji, jj, jk) - zd_lw(ji, jj, jk) * zd_up(ji, jj, jk - 1) / zdiag(ji, jj, jk - 1)
        END DO
      END DO
    END DO
    DO jk = 2, jpk
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zd_lw(ji, jj, jk) = psi(ji, jj, jk) - zd_lw(ji, jj, jk) / zdiag(ji, jj, jk - 1) * zd_lw(ji, jj, jk - 1)
        END DO
      END DO
    END DO
    DO jk = jpk - 1, 2, - 1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          psi(ji, jj, jk) = (zd_lw(ji, jj, jk) - zd_up(ji, jj, jk) * psi(ji, jj, jk + 1)) / zdiag(ji, jj, jk)
        END DO
      END DO
    END DO
    SELECT CASE (nn_clos)
    CASE (0)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            eps(ji, jj, jk) = rc03 * en(ji, jj, jk) * en(ji, jj, jk) * SQRT(en(ji, jj, jk)) / MAX(psi(ji, jj, jk), rn_epsmin)
          END DO
        END DO
      END DO
    CASE (1)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            eps(ji, jj, jk) = psi(ji, jj, jk)
          END DO
        END DO
      END DO
    CASE (2)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            eps(ji, jj, jk) = rc04 * en(ji, jj, jk) * psi(ji, jj, jk)
          END DO
        END DO
      END DO
    CASE (3)
      zcoef = rc0 ** (3._wp + rpp / rnn)
      zex1 = (1.5_wp + rmm / rnn)
      zex2 = - 1._wp / rnn
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            eps(ji, jj, jk) = zcoef * en(ji, jj, jk) ** zex1 * psi(ji, jj, jk) ** zex2
          END DO
        END DO
      END DO
    END SELECT
    DO jk = 1, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          eps(ji, jj, jk) = MAX(eps(ji, jj, jk), rn_epsmin)
          hmxl_n(ji, jj, jk) = rc03 * en(ji, jj, jk) * SQRT(en(ji, jj, jk)) / eps(ji, jj, jk)
          zrn2 = MAX(rn2(ji, jj, jk), rsmall)
          IF (ln_length_lim) hmxl_n(ji, jj, jk) = MIN(rn_clim_galp * SQRT(2._wp * en(ji, jj, jk) / zrn2), hmxl_n(ji, jj, jk))
        END DO
      END DO
    END DO
    SELECT CASE (nn_stab_func)
    CASE (0, 1)
      DO jk = 2, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zcof = hmxl_b(ji, jj, jk) * hmxl_b(ji, jj, jk) / (2._wp * eb(ji, jj, jk))
            gh = - rn2(ji, jj, jk) * zcof
            gh = MIN(gh, rgh0)
            gh = MAX(gh, rghmin)
            sh = ra2 * (1._wp - 6._wp * ra1 / rb1) / (1. - 3. * ra2 * gh * (6. * ra1 + rb2 * (1._wp - rc3)))
            sm = (rb1 ** (- 1._wp / 3._wp) + (18._wp * ra1 * ra1 + 9._wp * ra1 * ra2 * (1._wp - rc2)) * sh * gh) / (1._wp - 9._wp &
&* ra1 * ra2 * gh)
            zstt(ji, jj, jk) = rc_diff * sh * tmask(ji, jj, jk)
            zstm(ji, jj, jk) = rc_diff * sm * tmask(ji, jj, jk)
          END DO
        END DO
      END DO
    CASE (2, 3)
      DO jk = 2, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zcof = hmxl_b(ji, jj, jk) * hmxl_b(ji, jj, jk) / (2._wp * eb(ji, jj, jk))
            gh = - rn2(ji, jj, jk) * zcof
            gh = MIN(gh, rgh0)
            gh = MAX(gh, rghmin)
            gh = gh * rf6
            shr = p_sh2(ji, jj, jk) / MAX(p_avm(ji, jj, jk), rsmall)
            gm = MAX(shr * zcof, 1.E-10)
            gm = gm * rf6
            gm = MIN((rd0 - rd1 * gh + rd3 * gh * gh) / (rd2 - rd4 * gh), gm)
            rcff = rd0 - rd1 * gh + rd2 * gm + rd3 * gh * gh - rd4 * gh * gm + rd5 * gm * gm
            sm = (rs0 - rs1 * gh + rs2 * gm) / rcff
            sh = (rs4 - rs5 * gh + rs6 * gm) / rcff
            zstt(ji, jj, jk) = rc_diff * sh * tmask(ji, jj, jk)
            zstm(ji, jj, jk) = rc_diff * sm * tmask(ji, jj, jk)
          END DO
        END DO
      END DO
    END SELECT
    zstm(:, :, 1) = zstm(:, :, 2)
    zstm(:, :, jpk) = 0.
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zstm(ji, jj, mbkt(ji, jj) + 1) = zstm(ji, jj, mbkt(ji, jj))
      END DO
    END DO
    zstt(:, :, 1) = wmask(:, :, 1)
    zstt(:, :, jpk) = wmask(:, :, jpk)
    DO jk = 1, jpk
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zsqen = SQRT(2._wp * en(ji, jj, jk)) * hmxl_n(ji, jj, jk)
          zavt = zsqen * zstt(ji, jj, jk)
          zavm = zsqen * zstm(ji, jj, jk)
          p_avt(ji, jj, jk) = MAX(zavt, avtb(jk)) * wmask(ji, jj, jk)
          p_avm(ji, jj, jk) = MAX(zavm, avmb(jk))
        END DO
      END DO
    END DO
    p_avt(:, :, 1) = 0._wp
    IF (ln_ctl) THEN
      CALL prt_ctl(tab3d_1 = en, clinfo1 = ' gls  - e: ', tab3d_2 = p_avt, clinfo2 = ' t: ', kdim = jpk)
      CALL prt_ctl(tab3d_1 = p_avm, clinfo1 = ' gls  - m: ', kdim = jpk)
    END IF
  END SUBROUTINE zdf_gls
  SUBROUTINE zdf_gls_init
    INTEGER :: jk
    INTEGER :: ios
    REAL(KIND = wp) :: zcr
    NAMELIST /namzdf_gls/ rn_emin, rn_epsmin, ln_length_lim, rn_clim_galp, ln_sigpsi, rn_hsro, rn_crban, rn_charn, rn_frac_hs, &
&nn_bc_surf, nn_bc_bot, nn_z0_met, nn_stab_func, nn_clos
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namzdf_gls, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namzdf_gls in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namzdf_gls, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namzdf_gls in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namzdf_gls)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'zdf_gls_init : GLS turbulent closure scheme'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namzdf_gls : set gls mixing parameters'
      WRITE(numout, FMT = *) '      minimum value of en                           rn_emin        = ', rn_emin
      WRITE(numout, FMT = *) '      minimum value of eps                          rn_epsmin      = ', rn_epsmin
      WRITE(numout, FMT = *) '      Limit dissipation rate under stable stratif.  ln_length_lim  = ', ln_length_lim
      WRITE(numout, FMT = *) '      Galperin limit (Standard: 0.53, Holt: 0.26)   rn_clim_galp   = ', rn_clim_galp
      WRITE(numout, FMT = *) '      TKE Surface boundary condition                nn_bc_surf     = ', nn_bc_surf
      WRITE(numout, FMT = *) '      TKE Bottom boundary condition                 nn_bc_bot      = ', nn_bc_bot
      WRITE(numout, FMT = *) '      Modify psi Schmidt number (wb case)           ln_sigpsi      = ', ln_sigpsi
      WRITE(numout, FMT = *) '      Craig and Banner coefficient                  rn_crban       = ', rn_crban
      WRITE(numout, FMT = *) '      Charnock coefficient                          rn_charn       = ', rn_charn
      WRITE(numout, FMT = *) '      Surface roughness formula                     nn_z0_met      = ', nn_z0_met
      WRITE(numout, FMT = *) '      Wave height frac. (used if nn_z0_met=2)       rn_frac_hs     = ', rn_frac_hs
      WRITE(numout, FMT = *) '      Stability functions                           nn_stab_func   = ', nn_stab_func
      WRITE(numout, FMT = *) '      Type of closure                               nn_clos        = ', nn_clos
      WRITE(numout, FMT = *) '      Surface roughness (m)                         rn_hsro        = ', rn_hsro
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) '   Namelist namdrg_top/_bot:   used values:'
      WRITE(numout, FMT = *) '      top    ocean cavity roughness (m)             rn_z0(_top)   = ', r_z0_top
      WRITE(numout, FMT = *) '      Bottom seafloor     roughness (m)             rn_z0(_bot)   = ', r_z0_bot
      WRITE(numout, FMT = *)
    END IF
    IF (zdf_gls_alloc() /= 0) CALL ctl_stop('STOP', 'zdf_gls_init : unable to allocate arrays')
    IF (nn_bc_surf < 0 .OR. nn_bc_surf > 1) CALL ctl_stop('zdf_gls_init: bad flag: nn_bc_surf is 0 or 1')
    IF (nn_bc_surf < 0 .OR. nn_bc_surf > 1) CALL ctl_stop('zdf_gls_init: bad flag: nn_bc_surf is 0 or 1')
    IF (nn_z0_met < 0 .OR. nn_z0_met > 3) CALL ctl_stop('zdf_gls_init: bad flag: nn_z0_met is 0, 1, 2 or 3')
    IF (nn_z0_met == 3 .AND. .NOT. ln_wave) CALL ctl_stop('zdf_gls_init: nn_z0_met=3 requires ln_wave=T')
    IF (nn_stab_func < 0 .OR. nn_stab_func > 3) CALL ctl_stop('zdf_gls_init: bad flag: nn_stab_func is 0, 1, 2 and 3')
    IF (nn_clos < 0 .OR. nn_clos > 3) CALL ctl_stop('zdf_gls_init: bad flag: nn_clos is 0, 1, 2 or 3')
    SELECT CASE (nn_clos)
    CASE (0)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>   k-kl closure chosen (i.e. closed to the classical Mellor-Yamada)'
      IF (lwp) WRITE(numout, FMT = *)
      rpp = 0._wp
      rmm = 1._wp
      rnn = 1._wp
      rsc_tke = 1.96_wp
      rsc_psi = 1.96_wp
      rpsi1 = 0.9_wp
      rpsi3p = 1._wp
      rpsi2 = 0.5_wp
      SELECT CASE (nn_stab_func)
      CASE (0, 1)
        rpsi3m = 2.53_wp
      CASE (2)
        rpsi3m = 2.62_wp
      CASE (3)
        rpsi3m = 2.38
      END SELECT
    CASE (1)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>   k-eps closure chosen'
      IF (lwp) WRITE(numout, FMT = *)
      rpp = 3._wp
      rmm = 1.5_wp
      rnn = - 1._wp
      rsc_tke = 1._wp
      rsc_psi = 1.2_wp
      rpsi1 = 1.44_wp
      rpsi3p = 1._wp
      rpsi2 = 1.92_wp
      SELECT CASE (nn_stab_func)
      CASE (0, 1)
        rpsi3m = - 0.52_wp
      CASE (2)
        rpsi3m = - 0.629_wp
      CASE (3)
        rpsi3m = - 0.566
      END SELECT
    CASE (2)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>   k-omega closure chosen'
      IF (lwp) WRITE(numout, FMT = *)
      rpp = - 1._wp
      rmm = 0.5_wp
      rnn = - 1._wp
      rsc_tke = 2._wp
      rsc_psi = 2._wp
      rpsi1 = 0.555_wp
      rpsi3p = 1._wp
      rpsi2 = 0.833_wp
      SELECT CASE (nn_stab_func)
      CASE (0, 1)
        rpsi3m = - 0.58_wp
      CASE (2)
        rpsi3m = - 0.64_wp
      CASE (3)
        rpsi3m = - 0.64_wp
      END SELECT
    CASE (3)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>   generic closure chosen'
      IF (lwp) WRITE(numout, FMT = *)
      rpp = 2._wp
      rmm = 1._wp
      rnn = - 0.67_wp
      rsc_tke = 0.8_wp
      rsc_psi = 1.07_wp
      rpsi1 = 1._wp
      rpsi3p = 1._wp
      rpsi2 = 1.22_wp
      SELECT CASE (nn_stab_func)
      CASE (0, 1)
        rpsi3m = 0.1_wp
      CASE (2)
        rpsi3m = 0.05_wp
      CASE (3)
        rpsi3m = 0.05_wp
      END SELECT
    END SELECT
    SELECT CASE (nn_stab_func)
    CASE (0)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>   Stability functions from Galperin'
      rc2 = 0._wp
      rc3 = 0._wp
      rc_diff = 1._wp
      rc0 = 0.5544_wp
      rcm_sf = 0.9884_wp
      rghmin = - 0.28_wp
      rgh0 = 0.0233_wp
      rghcri = 0.02_wp
    CASE (1)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>   Stability functions from Kantha-Clayson'
      rc2 = 0.7_wp
      rc3 = 0.2_wp
      rc_diff = 1._wp
      rc0 = 0.5544_wp
      rcm_sf = 0.9884_wp
      rghmin = - 0.28_wp
      rgh0 = 0.0233_wp
      rghcri = 0.02_wp
    CASE (2)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>   Stability functions from Canuto A'
      rs0 = 1.5_wp * rl1 * rl5 * rl5
      rs1 = - rl4 * (rl6 + rl7) + 2._wp * rl4 * rl5 * (rl1 - (1._wp / 3._wp) * rl2 - rl3) + 1.5_wp * rl1 * rl5 * rl8
      rs2 = - (3._wp / 8._wp) * rl1 * (rl6 * rl6 - rl7 * rl7)
      rs4 = 2._wp * rl5
      rs5 = 2._wp * rl4
      rs6 = (2._wp / 3._wp) * rl5 * (3._wp * rl3 * rl3 - rl2 * rl2) - 0.5_wp * rl5 * rl1 * (3._wp * rl3 - rl2) + 0.75_wp * rl1 * &
&(rl6 - rl7)
      rd0 = 3._wp * rl5 * rl5
      rd1 = rl5 * (7._wp * rl4 + 3._wp * rl8)
      rd2 = rl5 * rl5 * (3._wp * rl3 * rl3 - rl2 * rl2) - 0.75_wp * (rl6 * rl6 - rl7 * rl7)
      rd3 = rl4 * (4._wp * rl4 + 3._wp * rl8)
      rd4 = rl4 * (rl2 * rl6 - 3._wp * rl3 * rl7 - rl5 * (rl2 * rl2 - rl3 * rl3)) + rl5 * rl8 * (3._wp * rl3 * rl3 - rl2 * rl2)
      rd5 = 0.25_wp * (rl2 * rl2 - 3._wp * rl3 * rl3) * (rl6 * rl6 - rl7 * rl7)
      rc0 = 0.5268_wp
      rf6 = 8._wp / (rc0 ** 6._wp)
      rc_diff = SQRT(2._wp) / (rc0 ** 3._wp)
      rcm_sf = 0.7310_wp
      rghmin = - 0.28_wp
      rgh0 = 0.0329_wp
      rghcri = 0.03_wp
    CASE (3)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>   Stability functions from Canuto B'
      rs0 = 1.5_wp * rm1 * rm5 * rm5
      rs1 = - rm4 * (rm6 + rm7) + 2._wp * rm4 * rm5 * (rm1 - (1._wp / 3._wp) * rm2 - rm3) + 1.5_wp * rm1 * rm5 * rm8
      rs2 = - (3._wp / 8._wp) * rm1 * (rm6 * rm6 - rm7 * rm7)
      rs4 = 2._wp * rm5
      rs5 = 2._wp * rm4
      rs6 = (2._wp / 3._wp) * rm5 * (3._wp * rm3 * rm3 - rm2 * rm2) - 0.5_wp * rm5 * rm1 * (3._wp * rm3 - rm2) + 0.75_wp * rm1 * &
&(rm6 - rm7)
      rd0 = 3._wp * rm5 * rm5
      rd1 = rm5 * (7._wp * rm4 + 3._wp * rm8)
      rd2 = rm5 * rm5 * (3._wp * rm3 * rm3 - rm2 * rm2) - 0.75_wp * (rm6 * rm6 - rm7 * rm7)
      rd3 = rm4 * (4._wp * rm4 + 3._wp * rm8)
      rd4 = rm4 * (rm2 * rm6 - 3._wp * rm3 * rm7 - rm5 * (rm2 * rm2 - rm3 * rm3)) + rm5 * rm8 * (3._wp * rm3 * rm3 - rm2 * rm2)
      rd5 = 0.25_wp * (rm2 * rm2 - 3._wp * rm3 * rm3) * (rm6 * rm6 - rm7 * rm7)
      rc0 = 0.5268_wp
      rf6 = 8._wp / (rc0 ** 6._wp)
      rc_diff = SQRT(2._wp) / (rc0 ** 3.)
      rcm_sf = 0.7470_wp
      rghmin = - 0.28_wp
      rgh0 = 0.0444_wp
      rghcri = 0.0414_wp
    END SELECT
    IF (ln_sigpsi) THEN
      ra_sf = - 1.5
      rsc_psi0 = rsc_tke / (24. * rpsi2) * (- 1. + (4. * rnn + ra_sf * (1. + 4. * rmm)) ** 2. / (ra_sf ** 2.))
    ELSE
      rsc_psi0 = rsc_psi
    END IF
    ra_sf = - 4._wp * rnn * SQRT(rsc_tke) / ((1._wp + 4._wp * rmm) * SQRT(rsc_tke) - SQRT(rsc_tke + 24._wp * rsc_psi0 * rpsi2))
    IF (rn_crban == 0._wp) THEN
      rl_sf = vkarmn
    ELSE
      rl_sf = rc0 * SQRT(rc0 / rcm_sf) * SQRT(((1._wp + 4._wp * rmm + 8._wp * rmm ** 2_wp) * rsc_tke + 12._wp * rsc_psi0 * rpsi2 - &
&(1._wp + 4._wp * rmm) * SQRT(rsc_tke * (rsc_tke + 24._wp * rsc_psi0 * rpsi2))) / (12._wp * rnn ** 2.))
    END IF
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) '   Limit values :'
      WRITE(numout, FMT = *) '      Parameter  m = ', rmm
      WRITE(numout, FMT = *) '      Parameter  n = ', rnn
      WRITE(numout, FMT = *) '      Parameter  p = ', rpp
      WRITE(numout, FMT = *) '      rpsi1    = ', rpsi1
      WRITE(numout, FMT = *) '      rpsi2    = ', rpsi2
      WRITE(numout, FMT = *) '      rpsi3m   = ', rpsi3m
      WRITE(numout, FMT = *) '      rpsi3p   = ', rpsi3p
      WRITE(numout, FMT = *) '      rsc_tke  = ', rsc_tke
      WRITE(numout, FMT = *) '      rsc_psi  = ', rsc_psi
      WRITE(numout, FMT = *) '      rsc_psi0 = ', rsc_psi0
      WRITE(numout, FMT = *) '      rc0      = ', rc0
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) '   Shear free turbulence parameters:'
      WRITE(numout, FMT = *) '      rcm_sf   = ', rcm_sf
      WRITE(numout, FMT = *) '      ra_sf    = ', ra_sf
      WRITE(numout, FMT = *) '      rl_sf    = ', rl_sf
    END IF
    rc02 = rc0 * rc0
    rc02r = 1. / rc02
    rc03 = rc02 * rc0
    rc04 = rc03 * rc0
    rsbc_tke1 = - 3._wp / 2._wp * rn_crban * ra_sf * rl_sf
    rsbc_tke2 = rdt * rn_crban / rl_sf
    zcr = MAX(rsmall, rsbc_tke1 ** (1. / (- ra_sf * 3._wp / 2._wp)) - 1._wp)
    rtrans = 0.2_wp / zcr
    rsbc_zs1 = rn_charn / grav
    rsbc_zs2 = rn_frac_hs / 0.85_wp / grav * 665._wp
    rsbc_psi1 = - 0.5_wp * rdt * rc0 ** (rpp - 2._wp * rmm) / rsc_psi
    rsbc_psi2 = - 0.5_wp * rdt * rc0 ** rpp * rnn * vkarmn ** rnn / rsc_psi
    rfact_tke = - 0.5_wp / rsc_tke * rdt
    rfact_psi = - 0.5_wp / rsc_psi * rdt
    zwall(:, :, :) = 1._wp * tmask(:, :, :)
    CALL gls_rst(nit000, 'READ')
    IF (lwxios) THEN
      CALL iom_set_rstw_var_active('en')
      CALL iom_set_rstw_var_active('avt_k')
      CALL iom_set_rstw_var_active('avm_k')
      CALL iom_set_rstw_var_active('hmxl_n')
    END IF
  END SUBROUTINE zdf_gls_init
  SUBROUTINE gls_rst(kt, cdrw)
    USE zdf_oce, ONLY: en, avt_k, avm_k
    INTEGER, INTENT(IN) :: kt
    CHARACTER(LEN = *), INTENT(IN) :: cdrw
    INTEGER :: jit, jk
    INTEGER :: id1, id2, id3, id4
    INTEGER :: ji, jj, ikbu, ikbv
    REAL(KIND = wp) :: cbx, cby
    IF (TRIM(cdrw) == 'READ') THEN
      IF (ln_rstart) THEN
        id1 = iom_varid(numror, 'en', ldstop = .FALSE.)
        id2 = iom_varid(numror, 'avt_k', ldstop = .FALSE.)
        id3 = iom_varid(numror, 'avm_k', ldstop = .FALSE.)
        id4 = iom_varid(numror, 'hmxl_n', ldstop = .FALSE.)
        IF (MIN(id1, id2, id3, id4) > 0) THEN
          CALL iom_get(numror, jpdom_autoglo, 'en', en, ldxios = lrxios)
          CALL iom_get(numror, jpdom_autoglo, 'avt_k', avt_k, ldxios = lrxios)
          CALL iom_get(numror, jpdom_autoglo, 'avm_k', avm_k, ldxios = lrxios)
          CALL iom_get(numror, jpdom_autoglo, 'hmxl_n', hmxl_n, ldxios = lrxios)
        ELSE
          IF (lwp) WRITE(numout, FMT = *)
          IF (lwp) WRITE(numout, FMT = *) '   ==>>   previous run without GLS scheme, set en and hmxl_n to background values'
          en(:, :, :) = rn_emin
          hmxl_n(:, :, :) = 0.05_wp
        END IF
      ELSE
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) '   ==>>   start from rest, set en and hmxl_n by background values'
        en(:, :, :) = rn_emin
        hmxl_n(:, :, :) = 0.05_wp
      END IF
    ELSE IF (TRIM(cdrw) == 'WRITE') THEN
      IF (lwp) WRITE(numout, FMT = *) '---- gls-rst ----'
      IF (lwxios) CALL iom_swap(cwxios_context)
      CALL iom_rstput(kt, nitrst, numrow, 'en', en, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'avt_k', avt_k, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'avm_k', avm_k, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'hmxl_n', hmxl_n, ldxios = lwxios)
      IF (lwxios) CALL iom_swap(cxios_context)
    END IF
  END SUBROUTINE gls_rst
END MODULE zdfgls