MODULE ldfslp
  USE oce
  USE dom_oce
  USE phycst
  USE zdfmxl
  USE eosbn2
  USE in_out_manager
  USE prtctl
  USE lbclnk
  USE lib_mpp
  USE lib_fortran
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ldf_slp
  PUBLIC :: ldf_slp_triad
  PUBLIC :: ldf_slp_init
  LOGICAL, PUBLIC :: l_ldfslp = .FALSE.
  LOGICAL, PUBLIC :: ln_traldf_iso = .TRUE.
  LOGICAL, PUBLIC :: ln_traldf_triad = .FALSE.
  LOGICAL, PUBLIC :: ln_dynldf_iso
  LOGICAL, PUBLIC :: ln_triad_iso = .FALSE.
  LOGICAL, PUBLIC :: ln_botmix_triad = .FALSE.
  REAL(KIND = wp), PUBLIC :: rn_sw_triad = 1._wp
  REAL(KIND = wp), PUBLIC :: rn_slpmax = 0.01_wp
  LOGICAL, PUBLIC :: l_grad_zps = .FALSE.
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: uslp, wslpi
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: vslp, wslpj
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: wslp2
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :, :, :) :: triadi_g, triadj_g
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :, :, :) :: triadi, triadj
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: ah_wslp2
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: akz
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: omlmask
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: uslpml, wslpiml
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: vslpml, wslpjml
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zslpml_hmlpu, zslpml_hmlpv
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zgru, zwz, zdzr
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zgrv, zww
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: z1_mlbw
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zalbet
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :, :) :: zdxrho, zdyrho, zdzrho
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :, :) :: zti_mlb, ztj_mlb
  REAL(KIND = wp) :: repsln = 1.E-25_wp
  CONTAINS
  SUBROUTINE ldf_slp(kt, prd, pn2)
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :) :: prd
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :) :: pn2
    INTEGER :: ji, jj, jk
    INTEGER :: ii0, ii1
    INTEGER :: ij0, ij1
    REAL(KIND = wp) :: zeps, zm1_g, zm1_2g, z1_16, zcofw, z1_slpmax
    REAL(KIND = wp) :: zci, zfi, zau, zbu, zai, zbi
    REAL(KIND = wp) :: zcj, zfj, zav, zbv, zaj, zbj
    REAL(KIND = wp) :: zck, zfk, zbw
    REAL(KIND = wp) :: zdepu, zdepv
    IF (ln_timing) CALL timing_start('ldf_slp')
    zeps = 1.E-20_wp
    z1_16 = 1.0_wp / 16._wp
    zm1_g = - 1.0_wp / grav
    zm1_2g = - 0.5_wp / grav
    z1_slpmax = 1._wp / rn_slpmax
    zww(:, :, :) = 0._wp
    zwz(:, :, :) = 0._wp
    !$OMP parallel default(shared), private(ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = 1, jpk
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zgru(ji, jj, jk) = umask(ji, jj, jk) * (prd(ji + 1, jj, jk) - prd(ji, jj, jk))
          zgrv(ji, jj, jk) = vmask(ji, jj, jk) * (prd(ji, jj + 1, jk) - prd(ji, jj, jk))
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    IF (ln_zps) THEN
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zgru(ji, jj, mbku(ji, jj)) = gru(ji, jj)
          zgrv(ji, jj, mbkv(ji, jj)) = grv(ji, jj)
        END DO
      END DO
    END IF
    IF (ln_zps .AND. ln_isfcav) THEN
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          IF (miku(ji, jj) > 1) zgru(ji, jj, miku(ji, jj)) = grui(ji, jj)
          IF (mikv(ji, jj) > 1) zgrv(ji, jj, mikv(ji, jj)) = grvi(ji, jj)
        END DO
      END DO
    END IF
    zdzr(:, :, 1) = 0._wp
    !$OMP parallel default(shared), private(jk)
    !$OMP do schedule(static)
    DO jk = 2, jpkm1
      zdzr(:, :, jk) = zm1_g * (prd(:, :, jk) + 1._wp) * (pn2(:, :, jk) + pn2(:, :, jk + 1)) * (1._wp - 0.5_wp * tmask(:, :, jk + &
&1))
    END DO
    !$OMP end do
    !$OMP end parallel
    CALL ldf_slp_mxl(prd, pn2, zgru, zgrv, zdzr)
    IF (ln_isfcav) THEN
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zslpml_hmlpu(ji, jj) = uslpml(ji, jj) / (MAX(hmlpt(ji, jj), hmlpt(ji + 1, jj), 5._wp) - MAX(risfdep(ji, jj), risfdep(ji &
&+ 1, jj)))
          zslpml_hmlpv(ji, jj) = vslpml(ji, jj) / (MAX(hmlpt(ji, jj), hmlpt(ji, jj + 1), 5._wp) - MAX(risfdep(ji, jj), risfdep(ji, &
&jj + 1)))
        END DO
      END DO
    ELSE
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zslpml_hmlpu(ji, jj) = uslpml(ji, jj) / MAX(hmlpt(ji, jj), hmlpt(ji + 1, jj), 5._wp)
          zslpml_hmlpv(ji, jj) = vslpml(ji, jj) / MAX(hmlpt(ji, jj), hmlpt(ji, jj + 1), 5._wp)
        END DO
      END DO
    END IF
    !$OMP parallel default(shared), private(ji,jj,jk,zau,zav,zbu,zbv,zdepu,zdepv,zfi,zfj)
    !$OMP do schedule(static)
    DO jk = 2, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zau = zgru(ji, jj, jk) * r1_e1u(ji, jj)
          zav = zgrv(ji, jj, jk) * r1_e2v(ji, jj)
          zbu = 0.5_wp * (zdzr(ji, jj, jk) + zdzr(ji + 1, jj, jk))
          zbv = 0.5_wp * (zdzr(ji, jj, jk) + zdzr(ji, jj + 1, jk))
          zbu = MIN(zbu, - z1_slpmax * ABS(zau), - 7.E+3_wp / e3u_n(ji, jj, jk) * ABS(zau))
          zbv = MIN(zbv, - z1_slpmax * ABS(zav), - 7.E+3_wp / e3v_n(ji, jj, jk) * ABS(zav))
          zfi = MAX(omlmask(ji, jj, jk), omlmask(ji + 1, jj, jk))
          zfj = MAX(omlmask(ji, jj, jk), omlmask(ji, jj + 1, jk))
          zdepu = 0.5_wp * ((gdept_n(ji, jj, jk) + gdept_n(ji + 1, jj, jk)) - 2 * MAX(risfdep(ji, jj), risfdep(ji + 1, jj)) - &
&e3u_n(ji, jj, miku(ji, jj)))
          zdepv = 0.5_wp * ((gdept_n(ji, jj, jk) + gdept_n(ji, jj + 1, jk)) - 2 * MAX(risfdep(ji, jj), risfdep(ji, jj + 1)) - &
&e3v_n(ji, jj, mikv(ji, jj)))
          zwz(ji, jj, jk) = ((1._wp - zfi) * zau / (zbu - zeps) + zfi * zdepu * zslpml_hmlpu(ji, jj)) * umask(ji, jj, jk)
          zww(ji, jj, jk) = ((1._wp - zfj) * zav / (zbv - zeps) + zfj * zdepv * zslpml_hmlpv(ji, jj)) * vmask(ji, jj, jk)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    CALL lbc_lnk_multi('ldfslp', zwz, 'U', - 1., zww, 'V', - 1.)
    !$OMP parallel default(shared), private(ji,jj,jk,zai,zaj,zbi,zbj,zbw,zci,zcj,zck,zfk)
    !$OMP do schedule(static)
    DO jk = 2, jpkm1
      DO jj = 2, jpjm1, MAX(1, jpj - 3)
        DO ji = 2, jpim1
          uslp(ji, jj, jk) = z1_16 * (zwz(ji - 1, jj - 1, jk) + zwz(ji + 1, jj - 1, jk) + zwz(ji - 1, jj + 1, jk) + zwz(ji + 1, jj &
&+ 1, jk) + 2. * (zwz(ji, jj - 1, jk) + zwz(ji - 1, jj, jk) + zwz(ji + 1, jj, jk) + zwz(ji, jj + 1, jk)) + 4. * zwz(ji, jj, jk))
          vslp(ji, jj, jk) = z1_16 * (zww(ji - 1, jj - 1, jk) + zww(ji + 1, jj - 1, jk) + zww(ji - 1, jj + 1, jk) + zww(ji + 1, jj &
&+ 1, jk) + 2. * (zww(ji, jj - 1, jk) + zww(ji - 1, jj, jk) + zww(ji + 1, jj, jk) + zww(ji, jj + 1, jk)) + 4. * zww(ji, jj, jk))
        END DO
      END DO
      DO jj = 3, jpj - 2
        DO ji = 2, jpim1
          uslp(ji, jj, jk) = z1_16 * (zwz(ji - 1, jj - 1, jk) + zwz(ji + 1, jj - 1, jk) + zwz(ji - 1, jj + 1, jk) + zwz(ji + 1, jj &
&+ 1, jk) + 2. * (zwz(ji, jj - 1, jk) + zwz(ji - 1, jj, jk) + zwz(ji + 1, jj, jk) + zwz(ji, jj + 1, jk)) + 4. * zwz(ji, jj, jk))
          vslp(ji, jj, jk) = z1_16 * (zww(ji - 1, jj - 1, jk) + zww(ji + 1, jj - 1, jk) + zww(ji - 1, jj + 1, jk) + zww(ji + 1, jj &
&+ 1, jk) + 2. * (zww(ji, jj - 1, jk) + zww(ji - 1, jj, jk) + zww(ji + 1, jj, jk) + zww(ji, jj + 1, jk)) + 4. * zww(ji, jj, jk))
        END DO
      END DO
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          uslp(ji, jj, jk) = uslp(ji, jj, jk) * (umask(ji, jj + 1, jk) + umask(ji, jj - 1, jk)) * 0.5_wp * (umask(ji, jj, jk) + &
&umask(ji, jj, jk + 1)) * 0.5_wp
          vslp(ji, jj, jk) = vslp(ji, jj, jk) * (vmask(ji + 1, jj, jk) + vmask(ji - 1, jj, jk)) * 0.5_wp * (vmask(ji, jj, jk) + &
&vmask(ji, jj, jk + 1)) * 0.5_wp
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP do schedule(static)
    DO jk = 2, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zbw = zm1_2g * pn2(ji, jj, jk) * (prd(ji, jj, jk) + prd(ji, jj, jk - 1) + 2.)
          zci = MAX(umask(ji - 1, jj, jk) + umask(ji, jj, jk) + umask(ji - 1, jj, jk - 1) + umask(ji, jj, jk - 1), zeps) * e1t(ji, &
&jj)
          zcj = MAX(vmask(ji, jj - 1, jk) + vmask(ji, jj, jk - 1) + vmask(ji, jj - 1, jk - 1) + vmask(ji, jj, jk), zeps) * e2t(ji, &
&jj)
          zai = (zgru(ji - 1, jj, jk) + zgru(ji, jj, jk - 1) + zgru(ji - 1, jj, jk - 1) + zgru(ji, jj, jk)) / zci * wmask(ji, jj, &
&jk)
          zaj = (zgrv(ji, jj - 1, jk) + zgrv(ji, jj, jk - 1) + zgrv(ji, jj - 1, jk - 1) + zgrv(ji, jj, jk)) / zcj * wmask(ji, jj, &
&jk)
          zbi = MIN(zbw, - 100._wp * ABS(zai), - 7.E+3_wp / e3w_n(ji, jj, jk) * ABS(zai))
          zbj = MIN(zbw, - 100._wp * ABS(zaj), - 7.E+3_wp / e3w_n(ji, jj, jk) * ABS(zaj))
          zfk = MAX(omlmask(ji, jj, jk), omlmask(ji, jj, jk - 1))
          zck = (gdepw_n(ji, jj, jk) - gdepw_n(ji, jj, mikt(ji, jj))) / MAX(hmlp(ji, jj) - gdepw_n(ji, jj, mikt(ji, jj)), 10._wp)
          zwz(ji, jj, jk) = (zai / (zbi - zeps) * (1._wp - zfk) + zck * wslpiml(ji, jj) * zfk) * wmask(ji, jj, jk)
          zww(ji, jj, jk) = (zaj / (zbj - zeps) * (1._wp - zfk) + zck * wslpjml(ji, jj) * zfk) * wmask(ji, jj, jk)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    CALL lbc_lnk_multi('ldfslp', zwz, 'T', - 1., zww, 'T', - 1.)
    !$OMP parallel default(shared), private(ji,jj,jk,zck,zcofw)
    !$OMP do schedule(static)
    DO jk = 2, jpkm1
      DO jj = 2, jpjm1, MAX(1, jpj - 3)
        DO ji = 2, jpim1
          zcofw = wmask(ji, jj, jk) * z1_16
          wslpi(ji, jj, jk) = (zwz(ji - 1, jj - 1, jk) + zwz(ji + 1, jj - 1, jk) + zwz(ji - 1, jj + 1, jk) + zwz(ji + 1, jj + 1, &
&jk) + 2. * (zwz(ji, jj - 1, jk) + zwz(ji - 1, jj, jk) + zwz(ji + 1, jj, jk) + zwz(ji, jj + 1, jk)) + 4. * zwz(ji, jj, jk)) * zcofw
          wslpj(ji, jj, jk) = (zww(ji - 1, jj - 1, jk) + zww(ji + 1, jj - 1, jk) + zww(ji - 1, jj + 1, jk) + zww(ji + 1, jj + 1, &
&jk) + 2. * (zww(ji, jj - 1, jk) + zww(ji - 1, jj, jk) + zww(ji + 1, jj, jk) + zww(ji, jj + 1, jk)) + 4. * zww(ji, jj, jk)) * zcofw
        END DO
      END DO
      DO jj = 3, jpj - 2
        DO ji = 2, jpim1
          zcofw = wmask(ji, jj, jk) * z1_16
          wslpi(ji, jj, jk) = (zwz(ji - 1, jj - 1, jk) + zwz(ji + 1, jj - 1, jk) + zwz(ji - 1, jj + 1, jk) + zwz(ji + 1, jj + 1, &
&jk) + 2. * (zwz(ji, jj - 1, jk) + zwz(ji - 1, jj, jk) + zwz(ji + 1, jj, jk) + zwz(ji, jj + 1, jk)) + 4. * zwz(ji, jj, jk)) * zcofw
          wslpj(ji, jj, jk) = (zww(ji - 1, jj - 1, jk) + zww(ji + 1, jj - 1, jk) + zww(ji - 1, jj + 1, jk) + zww(ji + 1, jj + 1, &
&jk) + 2. * (zww(ji, jj - 1, jk) + zww(ji - 1, jj, jk) + zww(ji + 1, jj, jk) + zww(ji, jj + 1, jk)) + 4. * zww(ji, jj, jk)) * zcofw
        END DO
      END DO
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zck = (umask(ji, jj, jk) + umask(ji - 1, jj, jk)) * (vmask(ji, jj, jk) + vmask(ji, jj - 1, jk)) * 0.25
          wslpi(ji, jj, jk) = wslpi(ji, jj, jk) * zck
          wslpj(ji, jj, jk) = wslpj(ji, jj, jk) * zck
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    CALL lbc_lnk_multi('ldfslp', uslp, 'U', - 1., vslp, 'V', - 1., wslpi, 'W', - 1., wslpj, 'W', - 1.)
    IF (ln_ctl) THEN
      CALL prt_ctl(tab3d_1 = uslp, clinfo1 = ' slp  - u : ', tab3d_2 = vslp, clinfo2 = ' v : ', kdim = jpk)
      CALL prt_ctl(tab3d_1 = wslpi, clinfo1 = ' slp  - wi: ', tab3d_2 = wslpj, clinfo2 = ' wj: ', kdim = jpk)
    END IF
    IF (ln_timing) CALL timing_stop('ldf_slp')
  END SUBROUTINE ldf_slp
  SUBROUTINE ldf_slp_triad(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk, jl, ip, jp, kp
    INTEGER :: iku, ikv
    REAL(KIND = wp) :: zfacti, zfactj
    REAL(KIND = wp) :: znot_thru_surface
    REAL(KIND = wp) :: zdit, zdis, zdkt, zbu, zbti, zisw
    REAL(KIND = wp) :: zdjt, zdjs, zdks, zbv, zbtj, zjsw
    REAL(KIND = wp) :: zdxrho_raw, zti_coord, zti_raw, zti_lim, zti_g_raw, zti_g_lim
    REAL(KIND = wp) :: zdyrho_raw, ztj_coord, ztj_raw, ztj_lim, ztj_g_raw, ztj_g_lim
    REAL(KIND = wp) :: zdzrho_raw
    REAL(KIND = wp) :: zbeta0, ze3_e1, ze3_e2
    IF (ln_timing) CALL timing_start('ldf_slp_triad')
    DO jl = 0, 1
      ip = jl
      jp = jl
      DO jk = 1, jpkm1
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zdit = (tsb(ji + 1, jj, jk, jp_tem) - tsb(ji, jj, jk, jp_tem))
            zdis = (tsb(ji + 1, jj, jk, jp_sal) - tsb(ji, jj, jk, jp_sal))
            zdjt = (tsb(ji, jj + 1, jk, jp_tem) - tsb(ji, jj, jk, jp_tem))
            zdjs = (tsb(ji, jj + 1, jk, jp_sal) - tsb(ji, jj, jk, jp_sal))
            zdxrho_raw = (- rab_b(ji + ip, jj, jk, jp_tem) * zdit + rab_b(ji + ip, jj, jk, jp_sal) * zdis) * r1_e1u(ji, jj)
            zdyrho_raw = (- rab_b(ji, jj + jp, jk, jp_tem) * zdjt + rab_b(ji, jj + jp, jk, jp_sal) * zdjs) * r1_e2v(ji, jj)
            zdxrho(ji + ip, jj, jk, 1 - ip) = SIGN(MAX(repsln, ABS(zdxrho_raw)), zdxrho_raw)
            zdyrho(ji, jj + jp, jk, 1 - jp) = SIGN(MAX(repsln, ABS(zdyrho_raw)), zdyrho_raw)
          END DO
        END DO
      END DO
      IF (ln_zps .AND. l_grad_zps) THEN
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            iku = mbku(ji, jj)
            ikv = mbkv(ji, jj)
            zdit = gtsu(ji, jj, jp_tem)
            zdjt = gtsv(ji, jj, jp_tem)
            zdis = gtsu(ji, jj, jp_sal)
            zdjs = gtsv(ji, jj, jp_sal)
            zdxrho_raw = (- rab_b(ji + ip, jj, iku, jp_tem) * zdit + rab_b(ji + ip, jj, iku, jp_sal) * zdis) * r1_e1u(ji, jj)
            zdyrho_raw = (- rab_b(ji, jj + jp, ikv, jp_tem) * zdjt + rab_b(ji, jj + jp, ikv, jp_sal) * zdjs) * r1_e2v(ji, jj)
            zdxrho(ji + ip, jj, iku, 1 - ip) = SIGN(MAX(repsln, ABS(zdxrho_raw)), zdxrho_raw)
            zdyrho(ji, jj + jp, ikv, 1 - jp) = SIGN(MAX(repsln, ABS(zdyrho_raw)), zdyrho_raw)
          END DO
        END DO
      END IF
    END DO
    DO kp = 0, 1
      DO jk = 1, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            IF (jk + kp > 1) THEN
              zdkt = (tsb(ji, jj, jk + kp - 1, jp_tem) - tsb(ji, jj, jk + kp, jp_tem))
              zdks = (tsb(ji, jj, jk + kp - 1, jp_sal) - tsb(ji, jj, jk + kp, jp_sal))
            ELSE
              zdkt = 0._wp
              zdks = 0._wp
            END IF
            zdzrho_raw = (- rab_b(ji, jj, jk + kp, jp_tem) * zdkt + rab_b(ji, jj, jk + kp, jp_sal) * zdks) / e3w_n(ji, jj, jk + kp)
            zdzrho(ji, jj, jk, kp) = - MIN(- repsln, zdzrho_raw)
          END DO
        END DO
      END DO
    END DO
    DO jj = 1, jpj
      DO ji = 1, jpi
        jk = MIN(nmln(ji, jj), mbkt(ji, jj)) + 1
        z1_mlbw(ji, jj) = 1._wp / gdepw_n(ji, jj, jk)
      END DO
    END DO
    wslp2(:, :, :) = 0._wp
    triadi_g(:, :, 1, :, :) = 0._wp
    triadi_g(:, :, jpk, :, :) = 0._wp
    triadj_g(:, :, 1, :, :) = 0._wp
    triadj_g(:, :, jpk, :, :) = 0._wp
    triadi(:, :, 1, :, :) = 0._wp
    triadj(:, :, jpk, :, :) = 0._wp
    triadj(:, :, 1, :, :) = 0._wp
    triadj(:, :, jpk, :, :) = 0._wp
    DO jl = 0, 1
      DO kp = 0, 1
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            ip = jl
            jp = jl
            jk = nmln(ji + ip, jj) + 1
            IF (jk > mbkt(ji + ip, jj)) THEN
              zti_mlb(ji + ip, jj, 1 - ip, kp) = 0.0_wp
            ELSE
              zti_g_raw = (zdxrho(ji + ip, jj, jk - kp, 1 - ip) / zdzrho(ji + ip, jj, jk - kp, kp) - (gdept_n(ji + 1, jj, jk - kp) &
&- gdept_n(ji, jj, jk - kp)) * r1_e1u(ji, jj)) * umask(ji, jj, jk)
              ze3_e1 = e3w_n(ji + ip, jj, jk - kp) * r1_e1u(ji, jj)
              zti_mlb(ji + ip, jj, 1 - ip, kp) = SIGN(MIN(rn_slpmax, 5.0_wp * ze3_e1, ABS(zti_g_raw)), zti_g_raw)
            END IF
            jk = nmln(ji, jj + jp) + 1
            IF (jk > mbkt(ji, jj + jp)) THEN
              ztj_mlb(ji, jj + jp, 1 - jp, kp) = 0.0_wp
            ELSE
              ztj_g_raw = (zdyrho(ji, jj + jp, jk - kp, 1 - jp) / zdzrho(ji, jj + jp, jk - kp, kp) - (gdept_n(ji, jj + 1, jk - kp) &
&- gdept_n(ji, jj, jk - kp)) / e2v(ji, jj)) * vmask(ji, jj, jk)
              ze3_e2 = e3w_n(ji, jj + jp, jk - kp) / e2v(ji, jj)
              ztj_mlb(ji, jj + jp, 1 - jp, kp) = SIGN(MIN(rn_slpmax, 5.0_wp * ze3_e2, ABS(ztj_g_raw)), ztj_g_raw)
            END IF
          END DO
        END DO
      END DO
    END DO
    DO kp = 0, 1
      DO jl = 0, 1
        ip = jl
        jp = jl
        DO jk = 1, jpkm1
          znot_thru_surface = REAL(1 - 1 / (jk + kp), wp)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              zti_raw = zdxrho(ji + ip, jj, jk, 1 - ip) / zdzrho(ji + ip, jj, jk, kp)
              ztj_raw = zdyrho(ji, jj + jp, jk, 1 - jp) / zdzrho(ji, jj + jp, jk, kp)
              zti_coord = znot_thru_surface * (gdept_n(ji + 1, jj, jk) - gdept_n(ji, jj, jk)) * r1_e1u(ji, jj)
              ztj_coord = znot_thru_surface * (gdept_n(ji, jj + 1, jk) - gdept_n(ji, jj, jk)) * r1_e2v(ji, jj)
              zti_g_raw = zti_raw - zti_coord
              ztj_g_raw = ztj_raw - ztj_coord
              ze3_e1 = e3w_n(ji + ip, jj, jk + kp) * r1_e1u(ji, jj)
              ze3_e2 = e3w_n(ji, jj + jp, jk + kp) * r1_e2v(ji, jj)
              zti_g_lim = SIGN(MIN(rn_slpmax, 5.0_wp * ze3_e1, ABS(zti_g_raw)), zti_g_raw)
              ztj_g_lim = SIGN(MIN(rn_slpmax, 5.0_wp * ze3_e2, ABS(ztj_g_raw)), ztj_g_raw)
              zfacti = REAL(1 - 1 / (1 + (jk + kp - 1) / nmln(ji + ip, jj)), wp)
              zfactj = REAL(1 - 1 / (1 + (jk + kp - 1) / nmln(ji, jj + jp)), wp)
              zti_g_lim = (zfacti * zti_g_lim + (1._wp - zfacti) * zti_mlb(ji + ip, jj, 1 - ip, kp) * gdepw_n(ji + ip, jj, jk + &
&kp) * z1_mlbw(ji + ip, jj)) * umask(ji, jj, jk + kp)
              ztj_g_lim = (zfactj * ztj_g_lim + (1._wp - zfactj) * ztj_mlb(ji, jj + jp, 1 - jp, kp) * gdepw_n(ji, jj + jp, jk + &
&kp) * z1_mlbw(ji, jj + jp)) * vmask(ji, jj, jk + kp)
              triadi_g(ji + ip, jj, jk, 1 - ip, kp) = zti_g_lim
              triadj_g(ji, jj + jp, jk, 1 - jp, kp) = ztj_g_lim
              zti_lim = (zti_g_lim + zti_coord) * umask(ji, jj, jk + kp)
              ztj_lim = (ztj_g_lim + ztj_coord) * vmask(ji, jj, jk + kp)
              IF (ln_triad_iso) THEN
                zti_raw = zti_lim * zti_lim / zti_raw
                ztj_raw = ztj_lim * ztj_lim / ztj_raw
                zti_raw = SIGN(MIN(ABS(zti_lim), ABS(zti_raw)), zti_raw)
                ztj_raw = SIGN(MIN(ABS(ztj_lim), ABS(ztj_raw)), ztj_raw)
                zti_lim = zfacti * zti_lim + (1._wp - zfacti) * zti_raw
                ztj_lim = zfactj * ztj_lim + (1._wp - zfactj) * ztj_raw
              END IF
              zisw = (1._wp - rn_sw_triad) + rn_sw_triad * 2._wp * ABS(0.5_wp - kp - (0.5_wp - ip) * SIGN(1._wp, zdxrho(ji + ip, &
&jj, jk, 1 - ip)))
              zjsw = (1._wp - rn_sw_triad) + rn_sw_triad * 2._wp * ABS(0.5_wp - kp - (0.5_wp - jp) * SIGN(1._wp, zdyrho(ji, jj + &
&jp, jk, 1 - jp)))
              triadi(ji + ip, jj, jk, 1 - ip, kp) = zti_lim * zisw
              triadj(ji, jj + jp, jk, 1 - jp, kp) = ztj_lim * zjsw
              zbu = e1e2u(ji, jj) * e3u_n(ji, jj, jk)
              zbv = e1e2v(ji, jj) * e3v_n(ji, jj, jk)
              zbti = e1e2t(ji + ip, jj) * e3w_n(ji + ip, jj, jk + kp)
              zbtj = e1e2t(ji, jj + jp) * e3w_n(ji, jj + jp, jk + kp)
              wslp2(ji + ip, jj, jk + kp) = wslp2(ji + ip, jj, jk + kp) + 0.25_wp * zbu / zbti * zti_g_lim * zti_g_lim
              wslp2(ji, jj + jp, jk + kp) = wslp2(ji, jj + jp, jk + kp) + 0.25_wp * zbv / zbtj * ztj_g_lim * ztj_g_lim
            END DO
          END DO
        END DO
      END DO
    END DO
    wslp2(:, :, 1) = 0._wp
    CALL lbc_lnk('ldfslp', wslp2, 'W', 1.)
    IF (ln_timing) CALL timing_stop('ldf_slp_triad')
  END SUBROUTINE ldf_slp_triad
  SUBROUTINE ldf_slp_mxl(prd, pn2, p_gru, p_grv, p_dzr)
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: prd
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pn2
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: p_gru, p_grv
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: p_dzr
    INTEGER :: ji, jj, jk
    INTEGER :: iku, ikv, ik, ikm1
    REAL(KIND = wp) :: zeps, zm1_g, zm1_2g, z1_slpmax
    REAL(KIND = wp) :: zci, zfi, zau, zbu, zai, zbi
    REAL(KIND = wp) :: zcj, zfj, zav, zbv, zaj, zbj
    REAL(KIND = wp) :: zck, zfk, zbw
    zeps = 1.E-20_wp
    zm1_g = - 1.0_wp / grav
    zm1_2g = - 0.5_wp / grav
    z1_slpmax = 1._wp / rn_slpmax
    uslpml(1, :) = 0._wp
    uslpml(jpi, :) = 0._wp
    vslpml(1, :) = 0._wp
    vslpml(jpi, :) = 0._wp
    wslpiml(1, :) = 0._wp
    wslpiml(jpi, :) = 0._wp
    wslpjml(1, :) = 0._wp
    wslpjml(jpi, :) = 0._wp
    DO jk = 1, jpk
      DO jj = 1, jpj
        DO ji = 1, jpi
          ik = nmln(ji, jj) - 1
          IF (jk <= ik) THEN
            omlmask(ji, jj, jk) = 1._wp
          ELSE
            omlmask(ji, jj, jk) = 0._wp
          END IF
        END DO
      END DO
    END DO
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        iku = MIN(MAX(1, nmln(ji, jj), nmln(ji + 1, jj)), jpkm1)
        ikv = MIN(MAX(1, nmln(ji, jj), nmln(ji, jj + 1)), jpkm1)
        zbu = 0.5_wp * (p_dzr(ji, jj, iku) + p_dzr(ji + 1, jj, iku))
        zbv = 0.5_wp * (p_dzr(ji, jj, ikv) + p_dzr(ji, jj + 1, ikv))
        zau = p_gru(ji, jj, iku) * r1_e1u(ji, jj)
        zav = p_grv(ji, jj, ikv) * r1_e2v(ji, jj)
        zbu = MIN(zbu, - z1_slpmax * ABS(zau), - 7.E+3_wp / e3u_n(ji, jj, iku) * ABS(zau))
        zbv = MIN(zbv, - z1_slpmax * ABS(zav), - 7.E+3_wp / e3v_n(ji, jj, ikv) * ABS(zav))
        uslpml(ji, jj) = zau / (zbu - zeps) * umask(ji, jj, iku)
        vslpml(ji, jj) = zav / (zbv - zeps) * vmask(ji, jj, ikv)
        ik = MIN(nmln(ji, jj) + 1, jpk)
        ikm1 = MAX(1, ik - 1)
        zbw = zm1_2g * pn2(ji, jj, ik) * (prd(ji, jj, ik) + prd(ji, jj, ikm1) + 2.)
        zci = MAX(umask(ji - 1, jj, ik) + umask(ji, jj, ik) + umask(ji - 1, jj, ikm1) + umask(ji, jj, ikm1), zeps) * e1t(ji, jj)
        zcj = MAX(vmask(ji, jj - 1, ik) + vmask(ji, jj, ik) + vmask(ji, jj - 1, ikm1) + vmask(ji, jj, ikm1), zeps) * e2t(ji, jj)
        zai = (p_gru(ji - 1, jj, ik) + p_gru(ji, jj, ik) + p_gru(ji - 1, jj, ikm1) + p_gru(ji, jj, ikm1)) / zci * tmask(ji, jj, ik)
        zaj = (p_grv(ji, jj - 1, ik) + p_grv(ji, jj, ik) + p_grv(ji, jj - 1, ikm1) + p_grv(ji, jj, ikm1)) / zcj * tmask(ji, jj, ik)
        zbi = MIN(zbw, - 100._wp * ABS(zai), - 7.E+3_wp / e3w_n(ji, jj, ik) * ABS(zai))
        zbj = MIN(zbw, - 100._wp * ABS(zaj), - 7.E+3_wp / e3w_n(ji, jj, ik) * ABS(zaj))
        wslpiml(ji, jj) = zai / (zbi - zeps) * tmask(ji, jj, ik)
        wslpjml(ji, jj) = zaj / (zbj - zeps) * tmask(ji, jj, ik)
      END DO
    END DO
    CALL lbc_lnk_multi('ldfslp', uslpml, 'U', - 1., vslpml, 'V', - 1., wslpiml, 'W', - 1., wslpjml, 'W', - 1.)
  END SUBROUTINE ldf_slp_mxl
  SUBROUTINE ldf_slp_init
    INTEGER :: ji, jj, jk
    INTEGER :: ierr
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'ldf_slp_init : direction of lateral mixing'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
    END IF
    ALLOCATE(ah_wslp2(jpi, jpj, jpk), akz(jpi, jpj, jpk), STAT = ierr)
    IF (ierr > 0) CALL ctl_stop('STOP', 'ldf_slp_init : unable to allocate ah_slp2 or akz')
    ALLOCATE(zslpml_hmlpu(jpi, jpj), zslpml_hmlpv(jpi, jpj))
    ALLOCATE(zgru(jpi, jpj, jpk), zwz(jpi, jpj, jpk), zdzr(jpi, jpj, jpk))
    ALLOCATE(zgrv(jpi, jpj, jpk), zww(jpi, jpj, jpk))
    IF (ln_traldf_triad) THEN
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   triad) operator (Griffies)'
      ALLOCATE(triadi_g(jpi, jpj, jpk, 0 : 1, 0 : 1), triadj_g(jpi, jpj, jpk, 0 : 1, 0 : 1), triadi(jpi, jpj, jpk, 0 : 1, 0 : 1), &
&triadj(jpi, jpj, jpk, 0 : 1, 0 : 1), wslp2(jpi, jpj, jpk), STAT = ierr)
      IF (ierr > 0) CALL ctl_stop('STOP', 'ldf_slp_init : unable to allocate Griffies operator slope')
      IF (ln_dynldf_iso) CALL ctl_stop('ldf_slp_init: Griffies operator on momentum not supported')
    ELSE
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   iso operator (Madec)'
      ALLOCATE(omlmask(jpi, jpj, jpk), uslp(jpi, jpj, jpk), uslpml(jpi, jpj), wslpi(jpi, jpj, jpk), wslpiml(jpi, jpj), vslp(jpi, &
&jpj, jpk), vslpml(jpi, jpj), wslpj(jpi, jpj, jpk), wslpjml(jpi, jpj), STAT = ierr)
      IF (ierr > 0) CALL ctl_stop('STOP', 'ldf_slp_init : unable to allocate Madec operator slope ')
      ALLOCATE(z1_mlbw(jpi, jpj))
      ALLOCATE(zalbet(jpi, jpj, jpk))
      ALLOCATE(zdxrho(jpi, jpj, jpk, 0 : 1), zdyrho(jpi, jpj, jpk, 0 : 1), zdzrho(jpi, jpj, jpk, 0 : 1))
      ALLOCATE(zti_mlb(jpi, jpj, 0 : 1, 0 : 1), ztj_mlb(jpi, jpj, 0 : 1, 0 : 1))
      uslp(:, :, :) = 0._wp
      uslpml(:, :) = 0._wp
      vslp(:, :, :) = 0._wp
      vslpml(:, :) = 0._wp
      wslpi(:, :, :) = 0._wp
      wslpiml(:, :) = 0._wp
      wslpj(:, :, :) = 0._wp
      wslpjml(:, :) = 0._wp
    END IF
  END SUBROUTINE ldf_slp_init
END MODULE ldfslp