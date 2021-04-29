MODULE obs_averg_h2d
  USE par_kind, ONLY: wp
  USE par_oce, ONLY: jpi, jpj
  USE phycst, ONLY: rad, ra, rpi
  USE dom_oce, ONLY: e1t, e2t, e1f, e2f, glamt, gphit, nproc
  USE in_out_manager
  USE obs_const, ONLY: obfillflt
  USE obs_utils
  USE lib_mpp, ONLY: ctl_warn, ctl_stop, mpp_min, lk_mpp
  IMPLICIT NONE
  PRIVATE :: obs_avg_h2d_rad, obs_avg_h2d_rec, obs_deg2dist, obs_dist2corners
  PUBLIC :: obs_avg_h2d, obs_avg_h2d_init, obs_max_fpsize
  CONTAINS
  SUBROUTINE obs_avg_h2d_init(kpk, kpk2, kmaxifp, kmaxjfp, k2dint, plam, pphi, pglam, pgphi, pglamf, pgphif, pmask, plamscl, pphiscl, lindegrees, pweig, pobsmask, iminpoints)
    INTEGER, INTENT(IN) :: kpk, kpk2, kmaxifp, kmaxjfp, k2dint
    REAL(KIND = wp), INTENT(INOUT) :: plam, pphi
    REAL(KIND = wp), DIMENSION(kmaxifp, kmaxjfp), INTENT(IN) :: pglam, pgphi
    REAL(KIND = wp), DIMENSION(kmaxifp + 1, kmaxjfp + 1), INTENT(IN) :: pglamf, pgphif
    REAL(KIND = wp), DIMENSION(kmaxifp, kmaxjfp, kpk2), INTENT(IN) :: pmask
    REAL(KIND = wp), INTENT(IN) :: plamscl, pphiscl
    LOGICAL, INTENT(IN) :: lindegrees
    REAL(KIND = wp), DIMENSION(kmaxifp, kmaxjfp, kpk2), INTENT(OUT) :: pweig
    REAL(KIND = wp), DIMENSION(kpk2), INTENT(OUT) :: pobsmask
    INTEGER, INTENT(IN), OPTIONAL :: iminpoints
    INTEGER :: jk
    INTEGER :: ikmax
    IF (kpk2 == 1) THEN
      ikmax = 1
    ELSE IF (kpk2 == kpk) THEN
      ikmax = kpk - 1
    END IF
    SELECT CASE (k2dint)
    CASE (5)
      CALL obs_avg_h2d_rad(kpk2, ikmax, kmaxifp, kmaxjfp, plam, pphi, plamscl, pphiscl, lindegrees, pmask, pglam, pgphi, pglamf, pgphif, pweig)
    CASE (6)
      CALL obs_avg_h2d_rec(kpk2, ikmax, kmaxifp, kmaxjfp, plam, pphi, plamscl, pphiscl, lindegrees, pmask, pglam, pgphi, pglamf, pgphif, pweig)
    END SELECT
  END SUBROUTINE obs_avg_h2d_init
  SUBROUTINE obs_avg_h2d_rad(kpk2, kmax, kmaxifp, kmaxjfp, plam, pphi, plamscl, pphiscl, lindegrees, pmask, pglam, pgphi, pglamf, pgphif, pweig)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE phycst, ONLY: ra, rpi
    INTEGER, INTENT(IN) :: kpk2, kmax
    INTEGER, INTENT(IN) :: kmaxifp, kmaxjfp
    REAL(KIND = wp), INTENT(IN) :: plam, pphi
    REAL(KIND = wp), INTENT(IN) :: plamscl, pphiscl
    LOGICAL, INTENT(IN) :: lindegrees
    REAL(KIND = wp), DIMENSION(kmaxifp, kmaxjfp, kpk2), INTENT(IN) :: pmask
    REAL(KIND = wp), DIMENSION(kmaxifp, kmaxjfp), INTENT(IN) :: pglam, pgphi
    REAL(KIND = wp), DIMENSION(kmaxifp + 1, kmaxjfp + 1), INTENT(IN) :: pglamf, pgphif
    REAL(KIND = wp), DIMENSION(kmaxifp, kmaxjfp, kpk2), INTENT(OUT) :: pweig
    INTEGER :: ji, jj, jk
    INTEGER :: jvert, jis, jjs
    INTEGER :: jnumvert, jnumvertbig
    INTEGER, PARAMETER :: jnumsubgrid = 20
    REAL(KIND = wp), DIMENSION(4) :: zxvert, zyvert, zdist
    REAL(KIND = wp), DIMENSION(4) :: zxgrid, zygrid, zdgrid
    REAL(KIND = wp) :: zdx, zdy, zarea_subbox, zxpos, zypos, zsubdist, zarea_fp, zareabox
    REAL(KIND = wp) :: zphiscl_m, zlamscl_m
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_avg_h2d_rad', 'r0', 0, 0)
    pweig(:, :, :) = 0.0_wp
    IF (pphiscl /= plamscl) THEN
      CALL ctl_warn('obs_avg_h2d_rad:', 'The two components of the obs footprint size are not equal', 'yet the radial option has been selected - using pphiscl here')
    END IF
    DO jk = 1, kmax
      DO ji = 1, kmaxifp
        DO jj = 1, kmaxjfp
          IF (pmask(ji, jj, jk) == 1.0_wp) THEN
            IF (lindegrees) THEN
              CALL obs_deg2dist(1, 1, pglam(ji, jj), pgphi(ji, jj), plamscl, pphiscl, zlamscl_m, zphiscl_m)
            ELSE
              zphiscl_m = pphiscl
            END IF
            CALL obs_dist2corners(pglamf(ji, jj), pglamf(ji + 1, jj), pglamf(ji, jj + 1), pglamf(ji + 1, jj + 1), pgphif(ji, jj), pgphif(ji + 1, jj), pgphif(ji, jj + 1), pgphif(ji + 1, jj + 1), pglam(ji, jj), pgphi(ji, jj), zxgrid, zygrid, zdgrid)
            zareabox = ABS(zxgrid(1) - zxgrid(2)) * ABS(zygrid(1) - zygrid(4))
            CALL obs_dist2corners(pglamf(ji, jj), pglamf(ji + 1, jj), pglamf(ji, jj + 1), pglamf(ji + 1, jj + 1), pgphif(ji, jj), pgphif(ji + 1, jj), pgphif(ji, jj + 1), pgphif(ji + 1, jj + 1), plam, pphi, zxvert, zyvert, zdist)
            jnumvert = 0
            jnumvertbig = 0
            DO jvert = 1, 4
              IF (zdist(jvert) <= (zphiscl_m / 2.0_wp)) jnumvert = jnumvert + 1
              IF (zdist(jvert) - (0.5_wp * zareabox) <= (zphiscl_m / 2.0)) jnumvertbig = jnumvertbig + 1
            END DO
            IF (jnumvertbig == 0) CYCLE
            IF (jnumvert == 4) THEN
              pweig(ji, jj, jk) = 1.0_wp
              CYCLE
            END IF
            zdx = ABS(zxvert(3) - zxvert(4)) / REAL(jnumsubgrid, wp)
            zdy = ABS(zyvert(1) - zyvert(4)) / REAL(jnumsubgrid, wp)
            zarea_subbox = zdx * zdy
            zarea_fp = 0.0_wp
            DO jis = 1, jnumsubgrid
              zxpos = zxvert(4) + (REAL(jis, wp) * zdx) - (0.5_wp * zdx)
              DO jjs = 1, jnumsubgrid
                zypos = zyvert(4) + (REAL(jjs, wp) * zdy) - (0.5_wp * zdy)
                zsubdist = SQRT((zxpos * zxpos) + (zypos * zypos))
                IF (zsubdist < (zphiscl_m / 2.0_wp)) zarea_fp = zarea_fp + zarea_subbox
              END DO
            END DO
            pweig(ji, jj, jk) = MIN(zarea_fp / zareabox, 1.0_wp)
          END IF
        END DO
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_avg_h2d_rad
  SUBROUTINE obs_avg_h2d_rec(kpk2, kmax, kmaxifp, kmaxjfp, plam, pphi, plamscl, pphiscl, lindegrees, pmask, pglam, pgphi, pglamf, pgphif, pweig)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE phycst, ONLY: ra, rpi
    INTEGER, INTENT(IN) :: kpk2, kmax
    INTEGER, INTENT(IN) :: kmaxifp, kmaxjfp
    REAL(KIND = wp), INTENT(IN) :: plam, pphi
    REAL(KIND = wp), INTENT(IN) :: plamscl, pphiscl
    LOGICAL, INTENT(IN) :: lindegrees
    REAL(KIND = wp), DIMENSION(kmaxifp, kmaxjfp, kpk2), INTENT(IN) :: pmask
    REAL(KIND = wp), DIMENSION(kmaxifp, kmaxjfp), INTENT(IN) :: pglam, pgphi
    REAL(KIND = wp), DIMENSION(kmaxifp + 1, kmaxjfp + 1), INTENT(IN) :: pglamf, pgphif
    REAL(KIND = wp), DIMENSION(kmaxifp, kmaxjfp, kpk2), INTENT(OUT) :: pweig
    INTEGER :: ji, jj, jk
    INTEGER :: jvert
    INTEGER, DIMENSION(4) :: jnumvert
    REAL(KIND = wp), DIMENSION(4) :: zxvert, zyvert
    REAL(KIND = wp), DIMENSION(4) :: zdist
    REAL(KIND = wp), DIMENSION(4) :: zxgrid, zygrid, zdgrid
    REAL(KIND = wp) :: zareabox, zarea_fp, zarea_intersect
    REAL(KIND = wp) :: zlamscl_m, zphiscl_m
    REAL(KIND = wp) :: z_awidth, z_aheight, z_cwidth, z_cheight
    REAL(KIND = wp) :: zleft, zright, ztop, zbottom
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_avg_h2d_rec', 'r0', 0, 0)
    pweig(:, :, :) = 0.0_wp
    DO jk = 1, kmax
      DO ji = 1, kmaxifp
        DO jj = 1, kmaxjfp
          IF (pmask(ji, jj, jk) == 1.0_wp) THEN
            IF (lindegrees) THEN
              CALL obs_deg2dist(1, 1, pglam(ji, jj), pgphi(ji, jj), plamscl, pphiscl, zlamscl_m, zphiscl_m)
            ELSE
              zlamscl_m = plamscl
              zphiscl_m = pphiscl
            END IF
            CALL obs_dist2corners(pglamf(ji, jj), pglamf(ji + 1, jj), pglamf(ji, jj + 1), pglamf(ji + 1, jj + 1), pgphif(ji, jj), pgphif(ji + 1, jj), pgphif(ji, jj + 1), pgphif(ji + 1, jj + 1), pglam(ji, jj), pgphi(ji, jj), zxgrid, zygrid, zdgrid)
            z_awidth = ABS(zxgrid(1) - zxgrid(2))
            z_aheight = ABS(zygrid(1) - zygrid(4))
            zareabox = z_awidth * z_aheight
            zarea_fp = zlamscl_m * zphiscl_m
            CALL obs_dist2corners(pglamf(ji, jj), pglamf(ji + 1, jj), pglamf(ji, jj + 1), pglamf(ji + 1, jj + 1), pgphif(ji, jj), pgphif(ji + 1, jj), pgphif(ji, jj + 1), pgphif(ji + 1, jj + 1), plam, pphi, zxvert, zyvert, zdist)
            z_cwidth = MAX(zxvert(1), zxvert(2), - zlamscl_m / 2.0_wp, zlamscl_m / 2.0_wp) - MIN(zxvert(1), zxvert(2), - zlamscl_m / 2.0_wp, zlamscl_m / 2.0_wp)
            z_cheight = MAX(zyvert(1), zyvert(4), zphiscl_m / 2.0_wp, - zphiscl_m / 2.0_wp) - MIN(zyvert(1), zyvert(4), zphiscl_m / 2.0_wp, - zphiscl_m / 2.0_wp)
            IF ((z_cwidth >= z_awidth + zlamscl_m) .OR. (z_cheight >= z_aheight + zphiscl_m)) THEN
              pweig(ji, jj, jk) = 0.0_wp
            ELSE IF ((z_cwidth == zlamscl_m) .AND. (z_cheight == zphiscl_m)) THEN
              pweig(ji, jj, jk) = 1.0_wp
            ELSE IF ((z_cwidth == z_awidth) .AND. (z_cheight == z_aheight)) THEN
              pweig(ji, jj, jk) = zarea_fp / zareabox
            ELSE
              zleft = MAX(zxvert(1), - zlamscl_m / 2.0_wp)
              zright = MIN(zxvert(2), zlamscl_m / 2.0_wp)
              zbottom = MAX(zyvert(4), - zphiscl_m / 2.0_wp)
              ztop = MIN(zyvert(1), zphiscl_m / 2.0_wp)
              IF ((zleft < zright) .AND. (zbottom < ztop)) THEN
                zarea_intersect = (zright - zleft) * (ztop - zbottom)
                pweig(ji, jj, jk) = zarea_intersect / zareabox
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_avg_h2d_rec
  SUBROUTINE obs_avg_h2d(kpk, kpk2, kmaxifp, kmaxjfp, pweig, pmod, pobsk)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kpk, kpk2
    INTEGER, INTENT(IN) :: kmaxifp, kmaxjfp
    REAL(KIND = wp), DIMENSION(kmaxifp, kmaxjfp, kpk2), INTENT(IN) :: pweig
    REAL(KIND = wp), DIMENSION(kmaxifp, kmaxjfp, kpk2), INTENT(IN) :: pmod
    REAL(KIND = wp), DIMENSION(kpk2), INTENT(OUT) :: pobsk
    INTEGER :: jk
    INTEGER :: ikmax
    REAL(KIND = wp) :: zsum
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_avg_h2d', 'r0', 0, 0)
    IF (kpk2 == 1) THEN
      ikmax = 1
    ELSE IF (kpk2 == kpk) THEN
      ikmax = kpk - 1
    END IF
    pobsk = obfillflt
    DO jk = 1, ikmax
      zsum = SUM(pweig(:, :, jk))
      IF (zsum /= 0.0_wp) THEN
        pobsk(jk) = SUM(pweig(:, :, jk) * pmod(:, :, jk), Mask = pweig(:, :, jk) > 0.0_wp)
        pobsk(jk) = pobsk(jk) / zsum
      END IF
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_avg_h2d
  SUBROUTINE obs_max_fpsize(k2dint, plamscl, pphiscl, lindegrees, pmask, kmaxifp, kmaxjfp)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: k2dint
    REAL(KIND = wp), INTENT(IN) :: plamscl, pphiscl
    LOGICAL, INTENT(IN) :: lindegrees
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pmask
    INTEGER, INTENT(OUT) :: kmaxifp, kmaxjfp
    REAL(KIND = wp) :: ze1min, ze2min
    REAL(KIND = wp) :: zphiscl_m, zlamscl_m
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_max_fpsize', 'r0', 0, 0)
    IF (k2dint <= 4) THEN
      kmaxifp = 2
      kmaxjfp = 2
    ELSE
      IF (lindegrees) THEN
        CALL obs_deg2dist(jpi, jpj, glamt, gphit, plamscl, pphiscl, zlamscl_m, zphiscl_m)
      ELSE
        zlamscl_m = plamscl
        zphiscl_m = pphiscl
      END IF
      ze1min = MINVAL(e1t(:, :), mask = pmask(:, :) == 1._wp)
      ze2min = MINVAL(e2t(:, :), mask = pmask(:, :) == 1._wp)
      IF (lk_mpp) THEN
        CALL mpp_min('obs_averg_h2d', ze1min)
        CALL mpp_min('obs_averg_h2d', ze2min)
      END IF
      kmaxifp = CEILING(zlamscl_m / ze1min) + 1
      kmaxjfp = CEILING(zphiscl_m / ze2min) + 1
      kmaxifp = kmaxifp + MOD(kmaxifp, 2)
      kmaxjfp = kmaxjfp + MOD(kmaxjfp, 2)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_max_fpsize
  SUBROUTINE obs_deg2dist(ki, kj, pglam, pgphi, plamscl_deg, pphiscl_deg, plamscl_max, pphiscl_max)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: ki, kj
    REAL(KIND = wp), INTENT(IN), DIMENSION(ki, kj) :: pglam, pgphi
    REAL(KIND = wp), INTENT(IN) :: plamscl_deg, pphiscl_deg
    REAL(KIND = wp), INTENT(OUT) :: plamscl_max, pphiscl_max
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zlon1, zlon2, zlat1, zlat2, zdlat, zdlon
    REAL(KIND = wp) :: za1, za2, za, zc, zd
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_deg2dist', 'r0', 0, 0)
    plamscl_max = - 1.0_wp
    pphiscl_max = - 1.0_wp
    DO ji = 1, ki
      DO jj = 1, kj
        zlon1 = rad * (pglam(ji, jj) + (0.5_wp * plamscl_deg))
        zlon2 = rad * (pglam(ji, jj) - (0.5_wp * plamscl_deg))
        zlat1 = rad * pgphi(ji, jj)
        zlat2 = rad * pgphi(ji, jj)
        zdlon = zlon2 - zlon1
        zdlat = zlat2 - zlat1
        za1 = SIN(zdlat / 2.0_wp)
        za2 = SIN(zdlon / 2.0_wp)
        za = (za1 * za1) + (COS(zlat1) * COS(zlat2) * (za2 * za2))
        zc = 2.0_wp * ATAN2(SQRT(za), SQRT(1.0_wp - za))
        zd = ra * zc
        IF (zd > plamscl_max) plamscl_max = zd
        zlon1 = rad * pglam(ji, jj)
        zlon2 = rad * pglam(ji, jj)
        zlat1 = rad * (pgphi(ji, jj) + (0.5_wp * pphiscl_deg))
        zlat2 = rad * (pgphi(ji, jj) - (0.5_wp * pphiscl_deg))
        zdlon = zlon2 - zlon1
        zdlat = zlat2 - zlat1
        za1 = SIN(zdlat / 2.0_wp)
        za2 = SIN(zdlon / 2.0_wp)
        za = (za1 * za1) + (COS(zlat1) * COS(zlat2) * (za2 * za2))
        zc = 2.0_wp * ATAN2(SQRT(za), SQRT(1.0_wp - za))
        zd = ra * zc
        IF (zd > pphiscl_max) pphiscl_max = zd
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_deg2dist
  SUBROUTINE obs_dist2corners(pglam_bl, pglam_br, pglam_tl, pglam_tr, pgphi_bl, pgphi_br, pgphi_tl, pgphi_tr, plam, pphi, pxvert, pyvert, pdist)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN) :: pglam_bl, pglam_br, pglam_tl, pglam_tr
    REAL(KIND = wp), INTENT(IN) :: pgphi_bl, pgphi_br, pgphi_tl, pgphi_tr
    REAL(KIND = wp), INTENT(IN) :: pphi, plam
    REAL(KIND = wp), DIMENSION(4), INTENT(OUT) :: pxvert, pyvert
    REAL(KIND = wp), DIMENSION(4), INTENT(OUT) :: pdist
    INTEGER :: jvert
    REAL(KIND = wp) :: zphi, zlam
    REAL(KIND = wp) :: za1, za2, zb1, zb2, zc1, zc2
    REAL(KIND = wp) :: zdist_centre_lat, zdist_centre_lon
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_dist2corners', 'r0', 0, 0)
    DO jvert = 1, 4
      SELECT CASE (jvert)
      CASE (1)
        zphi = pgphi_tl
        zlam = pglam_tl
      CASE (2)
        zphi = pgphi_tr
        zlam = pglam_tr
      CASE (3)
        zphi = pgphi_br
        zlam = pglam_br
      CASE (4)
        zphi = pgphi_bl
        zlam = pglam_bl
      END SELECT
      IF (zlam == plam) THEN
        pxvert(jvert) = 0.0_wp
      ELSE
        za1 = SIN(zphi * rad)
        za2 = SIN(zphi * rad)
        zb1 = COS(zphi * rad) * COS(zlam * rad)
        zb2 = COS(zphi * rad) * COS(plam * rad)
        zc1 = COS(zphi * rad) * SIN(zlam * rad)
        zc2 = COS(zphi * rad) * SIN(plam * rad)
        pxvert(jvert) = grt_cir_dis(za1, za2, zb1, zb2, zc1, zc2)
        pxvert(jvert) = ra * pxvert(jvert)
        IF (zlam < plam) pxvert(jvert) = - pxvert(jvert)
      END IF
      IF (zphi == pphi) THEN
        pyvert(jvert) = 0.0_wp
      ELSE
        za1 = SIN(zphi * rad)
        za2 = SIN(pphi * rad)
        zb1 = COS(zphi * rad) * COS(zlam * rad)
        zb2 = COS(pphi * rad) * COS(zlam * rad)
        zc1 = COS(zphi * rad) * SIN(zlam * rad)
        zc2 = COS(pphi * rad) * SIN(zlam * rad)
        pyvert(jvert) = grt_cir_dis(za1, za2, zb1, zb2, zc1, zc2)
        pyvert(jvert) = ra * pyvert(jvert)
        IF (zphi < pphi) pyvert(jvert) = - pyvert(jvert)
      END IF
      pdist(jvert) = SQRT((pxvert(jvert) * pxvert(jvert)) + (pyvert(jvert) * pyvert(jvert)))
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_dist2corners
END MODULE obs_averg_h2d