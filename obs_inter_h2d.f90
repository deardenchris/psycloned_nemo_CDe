MODULE obs_inter_h2d
  USE par_kind, ONLY: wp
  USE phycst, ONLY: rad, rpi
  USE in_out_manager
  USE obs_const, ONLY: obfillflt
  USE obs_utils
  USE lib_mpp, ONLY: ctl_warn, ctl_stop
  IMPLICIT NONE
  PRIVATE :: obs_int_h2d_ds1, obs_int_h2d_ds2, obs_int_h2d_bil, obs_int_h2d_bir, obs_int_h2d_pol, lu_invmat, lu_decomp, lu_backsb, &
&bil_wgt
  PUBLIC :: obs_int_h2d, obs_int_h2d_init
  CONTAINS
  SUBROUTINE obs_int_h2d_init(kpk, kpk2, k2dint, plam, pphi, pglam, pgphi, pmask, pweig, pobsmask, iminpoints)
    INTEGER, INTENT(IN) :: kpk, kpk2, k2dint
    REAL(KIND = wp), INTENT(INOUT) :: plam, pphi
    REAL(KIND = wp), DIMENSION(2, 2), INTENT(IN) :: pglam, pgphi
    REAL(KIND = wp), DIMENSION(2, 2, kpk2), INTENT(IN) :: pmask
    REAL(KIND = wp), DIMENSION(2, 2, kpk2), INTENT(OUT) :: pweig
    REAL(KIND = wp), DIMENSION(kpk2), INTENT(OUT) :: pobsmask
    INTEGER, INTENT(IN), OPTIONAL :: iminpoints
    INTEGER :: jk
    INTEGER :: ikmax, iamb1, iamb2
    REAL(KIND = wp) :: zphimm, zphimp, zphipm, zphipp, zlammm, zlammp, zlampm, zlampp, zphimin, zphimax, zlammin, zlammax
    REAL(KIND = wp), DIMENSION(kpk2) :: z2dmm, z2dmp, z2dpm, z2dpp, z2dmmt, z2dmpt, z2dpmt, z2dppt, zsum
    LOGICAL :: ll_ds1, ll_skip, ll_fail
    iamb1 = 10
    iamb2 = 3
    IF (kpk2 == 1) THEN
      ikmax = 1
    ELSE IF (kpk2 == kpk) THEN
      ikmax = kpk - 1
    END IF
    zphimm = pgphi(1, 1)
    zphimp = pgphi(1, 2)
    zphipm = pgphi(2, 1)
    zphipp = pgphi(2, 2)
    zlammm = pglam(1, 1)
    zlammp = pglam(1, 2)
    zlampm = pglam(2, 1)
    zlampp = pglam(2, 2)
    DO WHILE ((zlammm < 0.0_wp) .OR. (zlammm > 360.0_wp) .OR. (zlampm < 0.0_wp) .OR. (zlampm > 360.0_wp) .OR. (zlampp < 0.0_wp) &
&.OR. (zlampp > 360.0_wp) .OR. (zlammp < 0.0_wp) .OR. (zlammp > 360.0_wp))
      IF (zlammm < 0.0_wp) zlammm = zlammm + 360.0_wp
      IF (zlammm > 360.0_wp) zlammm = zlammm - 360.0_wp
      IF (zlammp < 0.0_wp) zlammp = zlammp + 360.0_wp
      IF (zlammp > 360.0_wp) zlammp = zlammp - 360.0_wp
      IF (zlampm < 0.0_wp) zlampm = zlampm + 360.0_wp
      IF (zlampm > 360.0_wp) zlampm = zlampm - 360.0_wp
      IF (zlampp < 0.0_wp) zlampp = zlampp + 360.0_wp
      IF (zlampp > 360.0_wp) zlampp = zlampp - 360.0_wp
    END DO
    DO WHILE ((plam < 0.0_wp) .OR. (plam > 360.0_wp))
      IF (plam < 0.0_wp) plam = plam + 360.0_wp
      IF (plam > 360.0_wp) plam = plam - 360.0_wp
    END DO
    ll_skip = .FALSE.
    IF ((ABS(zphimm - pphi) < 1.0E-6_wp) .AND. (ABS(zlammm - plam) < 1.0E-6_wp)) THEN
      z2dmm(:) = 1.0_wp
      z2dpm(:) = 0.0_wp
      z2dmp(:) = 0.0_wp
      z2dpp(:) = 0.0_wp
      ll_skip = .TRUE.
    END IF
    IF ((ABS(zphipm - pphi) < 1.0E-6_wp) .AND. (ABS(zlampm - plam) < 1.0E-6_wp)) THEN
      z2dmm(:) = 0.0_wp
      z2dpm(:) = 1.0_wp
      z2dmp(:) = 0.0_wp
      z2dpp(:) = 0.0_wp
      ll_skip = .TRUE.
    END IF
    IF ((ABS(zphimp - pphi) < 1.0E-6_wp) .AND. (ABS(zlammp - plam) < 1.0E-6_wp)) THEN
      z2dmm(:) = 0.0_wp
      z2dpm(:) = 0.0_wp
      z2dmp(:) = 1.0_wp
      z2dpp(:) = 0.0_wp
      ll_skip = .TRUE.
    END IF
    IF ((ABS(zphipp - pphi) < 1.0E-6_wp) .AND. (ABS(zlampp - plam) < 1.0E-6_wp)) THEN
      z2dmm(:) = 0.0_wp
      z2dpm(:) = 0.0_wp
      z2dmp(:) = 0.0_wp
      z2dpp(:) = 1.0_wp
      ll_skip = .TRUE.
    END IF
    IF (.NOT. ll_skip) THEN
      zphimin = MIN(zphimm, zphipm, zphipp, zphimp)
      zphimax = MAX(zphimm, zphipm, zphipp, zphimp)
      zlammin = MIN(zlammm, zlampm, zlampp, zlammp)
      zlammax = MAX(zlammm, zlampm, zlampp, zlammp)
      IF (((zlammax - zlammin) / (zphimax - zphimin)) > iamb1) THEN
        IF (iamb2 * zlammm < zlammax) zlammm = zlammm + 360.0_wp
        IF (iamb2 * zlammp < zlammax) zlammp = zlammp + 360.0_wp
        IF (iamb2 * zlampm < zlammax) zlampm = zlampm + 360.0_wp
        IF (iamb2 * zlampp < zlammax) zlampp = zlampp + 360.0_wp
      END IF
      zlammin = MIN(zlammm, zlampm, zlampp, zlammp)
      IF (zlammm > (zlammin + 180.0_wp)) zlammm = zlammm - 360.0_wp
      IF (zlammp > (zlammin + 180.0_wp)) zlammp = zlammp - 360.0_wp
      IF (zlampm > (zlammin + 180.0_wp)) zlampm = zlampm - 360.0_wp
      IF (zlampp > (zlammin + 180.0_wp)) zlampp = zlampp - 360.0_wp
      IF (plam < zlammin) plam = plam + 360.0_wp
      z2dmm = 0.0_wp
      z2dmp = 0.0_wp
      z2dpm = 0.0_wp
      z2dpp = 0.0_wp
      SELECT CASE (k2dint)
      CASE (0)
        CALL obs_int_h2d_ds1(kpk2, ikmax, pphi, plam, pmask, zphimm, zlammm, zphimp, zlammp, zphipm, zlampm, zphipp, zlampp, &
&z2dmm, z2dmp, z2dpm, z2dpp)
      CASE (1)
        CALL obs_int_h2d_ds2(kpk2, ikmax, pphi, plam, pmask, zphimm, zlammm, zphimp, zlammp, zphipm, zlampm, zphipp, zlampp, &
&z2dmm, z2dmp, z2dpm, z2dpp)
      CASE (2)
        CALL obs_int_h2d_bil(kpk2, ikmax, pphi, plam, pmask, zlammp, zphipm, zphipp, zlampp, z2dmm, z2dmp, z2dpm, z2dpp)
      CASE (3)
        CALL obs_int_h2d_bir(kpk2, ikmax, pphi, plam, pmask, zphimm, zlammm, zphimp, zlammp, zphipm, zlampm, zphipp, zlampp, &
&z2dmm, z2dmp, z2dpm, z2dpp, ll_fail)
        IF (ll_fail) THEN
          IF (lwp) THEN
            WRITE(numout, FMT = *) 'Bilinear weight computation failed'
            WRITE(numout, FMT = *) 'Switching to great circle distance'
            WRITE(numout, FMT = *)
          END IF
          CALL obs_int_h2d_ds1(kpk2, ikmax, pphi, plam, pmask, zphimm, zlammm, zphimp, zlammp, zphipm, zlampm, zphipp, zlampp, &
&z2dmm, z2dmp, z2dpm, z2dpp)
        END IF
      CASE (4)
        CALL obs_int_h2d_pol(kpk2, ikmax, pphi, plam, pmask, zphimm, zlammm, zphimp, zlammp, zphipm, zlampm, zphipp, zlampp, &
&z2dmm, z2dmp, z2dpm, z2dpp)
      END SELECT
    END IF
    pobsmask(:) = 0.0_wp
    pweig(:, :, :) = 0.0_wp
    ll_ds1 = .FALSE.
    DO jk = 1, ikmax
      IF (PRESENT(iminpoints)) THEN
        IF (NINT(SUM(pmask(:, :, jk))) < iminpoints) CYCLE
      END IF
      zsum(jk) = z2dmm(jk) + z2dmp(jk) + z2dpm(jk) + z2dpp(jk)
      IF (zsum(jk) /= 0.0_wp) THEN
        pweig(1, 1, jk) = z2dmm(jk)
        pweig(1, 2, jk) = z2dmp(jk)
        pweig(2, 1, jk) = z2dpm(jk)
        pweig(2, 2, jk) = z2dpp(jk)
        IF (((z2dmm(jk) > 0.0_wp) .AND. (pmask(1, 1, jk) == 1.0_wp)) .OR. ((z2dmp(jk) > 0.0_wp) .AND. (pmask(1, 2, jk) == 1.0_wp)) &
&.OR. ((z2dpm(jk) > 0.0_wp) .AND. (pmask(2, 1, jk) == 1.0_wp)) .OR. ((z2dpp(jk) > 0.0_wp) .AND. (pmask(2, 2, jk) == 1.0_wp))) &
&pobsmask(jk) = 1.0_wp
      ELSE
        IF ((pmask(1, 1, jk) /= 0.0_wp) .OR. (pmask(1, 2, jk) /= 0.0_wp) .OR. (pmask(2, 1, jk) /= 0.0_wp) .OR. (pmask(2, 2, jk) /= &
&0.0_wp)) THEN
          IF (.NOT. ll_ds1) THEN
            CALL obs_int_h2d_ds1(kpk2, ikmax, pphi, plam, pmask, zphimm, zlammm, zphimp, zlammp, zphipm, zlampm, zphipp, zlampp, &
&z2dmmt, z2dmpt, z2dpmt, z2dppt)
            ll_ds1 = .TRUE.
          END IF
          zsum(jk) = z2dmmt(jk) + z2dmpt(jk) + z2dpmt(jk) + z2dppt(jk)
          IF (zsum(jk) /= 0.0_wp) THEN
            pweig(1, 1, jk) = z2dmmt(jk)
            pweig(1, 2, jk) = z2dmpt(jk)
            pweig(2, 1, jk) = z2dpmt(jk)
            pweig(2, 2, jk) = z2dppt(jk)
            IF (((z2dmmt(jk) > 0.0_wp) .AND. (pmask(1, 1, jk) == 1.0_wp)) .OR. ((z2dmpt(jk) > 0.0_wp) .AND. (pmask(1, 2, jk) == &
&1.0_wp)) .OR. ((z2dpmt(jk) > 0.0_wp) .AND. (pmask(2, 1, jk) == 1.0_wp)) .OR. ((z2dppt(jk) > 0.0_wp) .AND. (pmask(2, 2, jk) == &
&1.0_wp))) pobsmask(jk) = 1.0_wp
          END IF
        END IF
      END IF
    END DO
  END SUBROUTINE obs_int_h2d_init
  SUBROUTINE obs_int_h2d(kpk, kpk2, pweig, pmod, pobsk)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kpk, kpk2
    REAL(KIND = wp), DIMENSION(2, 2, kpk2), INTENT(IN) :: pweig
    REAL(KIND = wp), DIMENSION(2, 2, kpk2), INTENT(IN) :: pmod
    REAL(KIND = wp), DIMENSION(kpk2), INTENT(OUT) :: pobsk
    INTEGER :: jk
    INTEGER :: ikmax
    REAL(KIND = wp) :: zsum
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_int_h2d', 'r0', 0, 0)
    IF (kpk2 == 1) THEN
      ikmax = 1
    ELSE IF (kpk2 == kpk) THEN
      ikmax = kpk - 1
    END IF
    pobsk(:) = obfillflt
    DO jk = 1, ikmax
      zsum = pweig(1, 1, jk) + pweig(1, 2, jk) + pweig(2, 1, jk) + pweig(2, 2, jk)
      IF (zsum /= 0.0_wp) THEN
        pobsk(jk) = (pweig(1, 1, jk) * pmod(1, 1, jk) + pweig(1, 2, jk) * pmod(1, 2, jk) + pweig(2, 1, jk) * pmod(2, 1, jk) + &
&pweig(2, 2, jk) * pmod(2, 2, jk)) / zsum
      END IF
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_int_h2d
  SUBROUTINE obs_int_h2d_ds1(kpk2, kmax, pphi, plam, pmask, pphimm, plammm, pphimp, plammp, pphipm, plampm, pphipp, plampp, p2dmm, &
&p2dmp, p2dpm, p2dpp)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kpk2, kmax
    REAL(KIND = wp), INTENT(IN) :: pphi, plam, pphimm, pphimp, pphipm, pphipp, plammm, plammp, plampm, plampp
    REAL(KIND = wp), DIMENSION(2, 2, kpk2), INTENT(IN) :: pmask
    REAL(KIND = wp), DIMENSION(kpk2), INTENT(OUT) :: p2dmm, p2dmp, p2dpm, p2dpp
    INTEGER :: jk
    REAL(KIND = wp) :: zphi2, zlam2, zcola, za2, zb2, zc2, zphimm2, zphimp2, zphipm2, zphipp2, zlammm2, zlammp2, zlampm2, zlampp2, &
&za1mm, za1mp, za1pm, za1pp, zcomm, zcomp, zcopm, zcopp, zb1mm, zb1mp, zb1pm, zb1pp, zc1mm, zc1mp, zc1pm, zc1pp, zsopmpp, zsommmp, &
&zsomm, zsomp, zsopm, zsopp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_int_h2d_ds1', 'r0', 0, 0)
    zphi2 = pphi * rad
    zlam2 = plam * rad
    zcola = COS(zphi2)
    za2 = SIN(zphi2)
    zb2 = zcola * COS(zlam2)
    zc2 = zcola * SIN(zlam2)
    zphimm2 = pphimm * rad
    zphimp2 = pphimp * rad
    zphipm2 = pphipm * rad
    zphipp2 = pphipp * rad
    zlammm2 = plammm * rad
    zlammp2 = plammp * rad
    zlampm2 = plampm * rad
    zlampp2 = plampp * rad
    za1mm = SIN(zphimm2)
    za1mp = SIN(zphimp2)
    za1pm = SIN(zphipm2)
    za1pp = SIN(zphipp2)
    zcomm = COS(zphimm2)
    zcomp = COS(zphimp2)
    zcopm = COS(zphipm2)
    zcopp = COS(zphipp2)
    zb1mm = zcomm * COS(zlammm2)
    zb1mp = zcomp * COS(zlammp2)
    zb1pm = zcopm * COS(zlampm2)
    zb1pp = zcopp * COS(zlampp2)
    zc1mm = zcomm * SIN(zlammm2)
    zc1mp = zcomp * SIN(zlammp2)
    zc1pm = zcopm * SIN(zlampm2)
    zc1pp = zcopp * SIN(zlampp2)
    zsomm = grt_cir_dis(za1mm, za2, zb1mm, zb2, zc1mm, zc2)
    zsomp = grt_cir_dis(za1mp, za2, zb1mp, zb2, zc1mp, zc2)
    zsopm = grt_cir_dis(za1pm, za2, zb1pm, zb2, zc1pm, zc2)
    zsopp = grt_cir_dis(za1pp, za2, zb1pp, zb2, zc1pp, zc2)
    zsopmpp = zsopm * zsopp
    zsommmp = zsomm * zsomp
    DO jk = 1, kmax
      p2dmm(jk) = zsomp * zsopmpp * pmask(1, 1, jk)
      p2dmp(jk) = zsomm * zsopmpp * pmask(1, 2, jk)
      p2dpm(jk) = zsopp * zsommmp * pmask(2, 1, jk)
      p2dpp(jk) = zsopm * zsommmp * pmask(2, 2, jk)
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_int_h2d_ds1
  SUBROUTINE obs_int_h2d_ds2(kpk2, kmax, pphi, plam, pmask, pphimm, plammm, pphimp, plammp, pphipm, plampm, pphipp, plampp, p2dmm, &
&p2dmp, p2dpm, p2dpp)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kpk2, kmax
    REAL(KIND = wp), INTENT(IN) :: pphi, plam, pphimm, pphimp, pphipm, pphipp, plammm, plammp, plampm, plampp
    REAL(KIND = wp), DIMENSION(2, 2, kpk2), INTENT(IN) :: pmask
    REAL(KIND = wp), DIMENSION(kpk2), INTENT(OUT) :: p2dmm, p2dmp, p2dpm, p2dpp
    INTEGER :: jk
    REAL(KIND = wp) :: zcosp, zdlmm, zdlmp, zdlpm, zdlpp, zdpmm, zdpmp, zdppm, zdppp, zsomm, zsomp, zsopm, zsopp, zsopmpp, zsommmp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_int_h2d_ds2', 'r0', 0, 0)
    zcosp = COS(pphi * rad)
    zdlmm = plammm - plam
    zdlmp = plammp - plam
    zdlpm = plampm - plam
    zdlpp = plampp - plam
    zdpmm = pphimm - pphi
    zdpmp = pphimp - pphi
    zdppm = pphipm - pphi
    zdppp = pphipp - pphi
    zsomm = grt_cir_dis_saa(zdlmm, zdpmm, zcosp)
    zsomp = grt_cir_dis_saa(zdlmp, zdpmp, zcosp)
    zsopm = grt_cir_dis_saa(zdlpm, zdppm, zcosp)
    zsopp = grt_cir_dis_saa(zdlpp, zdppp, zcosp)
    zsopmpp = zsopm * zsopp
    zsommmp = zsomm * zsomp
    DO jk = 1, kmax
      p2dmm(jk) = zsomp * zsopmpp * pmask(1, 1, jk)
      p2dmp(jk) = zsomm * zsopmpp * pmask(1, 2, jk)
      p2dpm(jk) = zsopp * zsommmp * pmask(2, 1, jk)
      p2dpp(jk) = zsopm * zsommmp * pmask(2, 2, jk)
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_int_h2d_ds2
  SUBROUTINE obs_int_h2d_bil(kpk2, kmax, pphi, plam, pmask, plammp, pphipm, pphipp, plampp, p2dmm, p2dmp, p2dpm, p2dpp)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kpk2, kmax
    REAL(KIND = wp), INTENT(IN) :: pphi, plam, pphipm, pphipp, plammp, plampp
    REAL(KIND = wp), DIMENSION(2, 2, kpk2), INTENT(IN) :: pmask
    REAL(KIND = wp), DIMENSION(kpk2), INTENT(OUT) :: p2dmm, p2dmp, p2dpm, p2dpp
    INTEGER :: jk
    REAL(KIND = wp) :: zdlmp, zdppm, zdlpp, zdppp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_int_h2d_bil', 'r0', 0, 0)
    zdlmp = ABS(plam - plammp)
    zdppm = ABS(pphi - pphipm)
    zdlpp = ABS(plampp - plam)
    zdppp = ABS(pphipp - pphi)
    DO jk = 1, kmax
      p2dmm(jk) = zdlpp * zdppp * pmask(1, 1, jk)
      p2dmp(jk) = zdlpp * zdppm * pmask(1, 2, jk)
      p2dpm(jk) = zdlmp * zdppp * pmask(2, 1, jk)
      p2dpp(jk) = zdlmp * zdppm * pmask(2, 2, jk)
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_int_h2d_bil
  SUBROUTINE obs_int_h2d_bir(kpk2, kmax, pphi, plam, pmask, pphimm, plammm, pphimp, plammp, pphipm, plampm, pphipp, plampp, p2dmm, &
&p2dmp, p2dpm, p2dpp, ldfail)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kpk2, kmax
    REAL(KIND = wp), INTENT(IN) :: pphi, plam, pphimm, pphimp, pphipm, pphipp, plammm, plammp, plampm, plampp
    REAL(KIND = wp), DIMENSION(2, 2, kpk2), INTENT(IN) :: pmask
    REAL(KIND = wp), DIMENSION(kpk2), INTENT(OUT) :: p2dmm, p2dmp, p2dpm, p2dpp
    LOGICAL, INTENT(OUT) :: ldfail
    INTEGER :: jk
    REAL(KIND = wp) :: zbiwmm, zbiwmp, zbiwpm, zbiwpp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_int_h2d_bir', 'r0', 0, 0)
    CALL bil_wgt(pphimm, pphimp, pphipm, pphipp, plammm, plammp, plampm, plampp, zbiwmm, zbiwmp, zbiwpm, zbiwpp, pphi, plam, ldfail)
    IF (.NOT. ldfail) THEN
      DO jk = 1, kmax
        p2dmm(jk) = zbiwmm * pmask(1, 1, jk)
        p2dmp(jk) = zbiwmp * pmask(1, 2, jk)
        p2dpm(jk) = zbiwpm * pmask(2, 1, jk)
        p2dpp(jk) = zbiwpp * pmask(2, 2, jk)
      END DO
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_int_h2d_bir
  SUBROUTINE obs_int_h2d_pol(kpk2, kmax, pphi, plam, pmask, pphimm, plammm, pphimp, plammp, pphipm, plampm, pphipp, plampp, p2dmm, &
&p2dmp, p2dpm, p2dpp)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kpk2, kmax
    REAL(KIND = wp), INTENT(IN) :: pphi, plam, pphimm, pphimp, pphipm, pphipp, plammm, plammp, plampm, plampp
    REAL(KIND = wp), DIMENSION(2, 2, kpk2), INTENT(IN) :: pmask
    REAL(KIND = wp), DIMENSION(kpk2), INTENT(OUT) :: p2dmm, p2dmp, p2dpm, p2dpp
    INTEGER :: jk
    REAL(KIND = wp) :: zplp
    REAL(KIND = wp), DIMENSION(4, 4) :: zmat, zmati
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_int_h2d_pol', 'r0', 0, 0)
    zmat(1, 1) = 1.0_wp
    zmat(1, 2) = 1.0_wp
    zmat(1, 3) = 1.0_wp
    zmat(1, 4) = 1.0_wp
    zmat(2, 1) = plammm
    zmat(2, 2) = plammp
    zmat(2, 3) = plampm
    zmat(2, 4) = plampp
    zmat(3, 1) = pphimm
    zmat(3, 2) = pphimp
    zmat(3, 3) = pphipm
    zmat(3, 4) = pphipp
    zmat(4, 1) = plammm * pphimm
    zmat(4, 2) = plammp * pphimp
    zmat(4, 3) = plampm * pphipm
    zmat(4, 4) = plampp * pphipp
    CALL lu_invmat(zmat, 4, zmati)
    zplp = plam * pphi
    DO jk = 1, kmax
      p2dmm(jk) = ABS(zmati(1, 1) + zmati(1, 2) * plam + zmati(1, 3) * pphi + zmati(1, 4) * zplp) * pmask(1, 1, jk)
      p2dmp(jk) = ABS(zmati(2, 1) + zmati(2, 2) * plam + zmati(2, 3) * pphi + zmati(2, 4) * zplp) * pmask(1, 2, jk)
      p2dpm(jk) = ABS(zmati(3, 1) + zmati(3, 2) * plam + zmati(3, 3) * pphi + zmati(3, 4) * zplp) * pmask(2, 1, jk)
      p2dpp(jk) = ABS(zmati(4, 1) + zmati(4, 2) * plam + zmati(4, 3) * pphi + zmati(4, 4) * zplp) * pmask(2, 2, jk)
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_int_h2d_pol
  SUBROUTINE bil_wgt(pphimm, pphimp, pphipm, pphipp, plammm, plammp, plampm, plampp, pbiwmm, pbiwmp, pbiwpm, pbiwpp, pphi, plam, &
&ldfail)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN) :: pphi, plam, pphimm, pphimp, pphipm, pphipp, plammm, plammp, plampm, plampp
    REAL(KIND = wp), INTENT(OUT) :: pbiwmm, pbiwmp, pbiwpm, pbiwpp
    LOGICAL, INTENT(OUT) :: ldfail
    INTEGER :: jiter
    INTEGER :: itermax
    REAL(KIND = wp) :: zphi, zlam, zphimm, zphimp, zphipm, zphipp, zlammm, zlammp, zlampm, zlampp, zdth1, zdth2, zdth3, zdthp, &
&zdph1, zdph2, zdph3, zdphp, zmat1, zmat2, zmat3, zmat4, zdeli, zdelj, ziguess, zjguess, zeps, zdeterm, z2pi, zhpi
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bil_wgt', 'r0', 0, 0)
    zphi = pphi * rad
    zlam = plam * rad
    zphimm = pphimm * rad
    zphimp = pphimp * rad
    zphipm = pphipm * rad
    zphipp = pphipp * rad
    zlammm = plammm * rad
    zlammp = plammp * rad
    zlampm = plampm * rad
    zlampp = plampp * rad
    ldfail = .FALSE.
    zdth1 = zphipm - zphimm
    zdth2 = zphimp - zphimm
    zdth3 = zphipp - zphipm - zdth2
    zdph1 = zlampm - zlammm
    zdph2 = zlammp - zlammm
    zdph3 = zlampp - zlampm
    z2pi = 2.0_wp * rpi
    IF (zdph1 > 3.0_wp * rpi) zdph1 = zdph1 - z2pi
    IF (zdph2 > 3.0_wp * rpi) zdph2 = zdph2 - z2pi
    IF (zdph3 > 3.0_wp * rpi) zdph3 = zdph3 - z2pi
    IF (zdph1 < - 3.0_wp * rpi) zdph1 = zdph1 + z2pi
    IF (zdph2 < - 3.0_wp * rpi) zdph2 = zdph2 + z2pi
    IF (zdph3 < - 3.0_wp * rpi) zdph3 = zdph3 + z2pi
    zdph3 = zdph3 - zdph2
    ziguess = 0.5_wp
    zjguess = 0.5_wp
    itermax = 100
    IF (wp == sp) THEN
      zeps = 1.0E-6_wp
    ELSE
      zeps = 1.0E-10_wp
    END IF
    jiter_loop:DO jiter = 1, itermax
      zdthp = zphi - zphimm - zdth1 * ziguess - zdth2 * zjguess - zdth3 * ziguess * zjguess
      zdphp = zlam - zlammm
      zhpi = 0.5_wp * rpi
      IF (zdphp > 3.0_wp * zhpi) zdphp = zdphp - z2pi
      IF (zdphp < - 3.0_wp * zhpi) zdphp = zdphp + z2pi
      zdphp = zdphp - zdph1 * ziguess - zdph2 * zjguess - zdph3 * ziguess * zjguess
      zmat1 = zdth1 + zdth3 * zjguess
      zmat2 = zdth2 + zdth3 * ziguess
      zmat3 = zdph1 + zdph3 * zjguess
      zmat4 = zdph2 + zdph3 * ziguess
      zdeterm = zmat1 * zmat4 - zmat2 * zmat3
      zdeli = (zdthp * zmat4 - zmat2 * zdphp) / zdeterm
      zdelj = (zmat1 * zdphp - zdthp * zmat3) / zdeterm
      IF (ABS(zdeli) < zeps .AND. ABS(zdelj) < zeps) EXIT jiter_loop
      ziguess = ziguess + zdeli
      zjguess = zjguess + zdelj
      IF (ziguess < 0) ziguess = 0.0_wp
      IF (zjguess < 0) zjguess = 0.0_wp
      IF (ziguess > 1) ziguess = 1.0_wp
      IF (zjguess > 1) zjguess = 1.0_wp
    END DO jiter_loop
    IF (jiter <= itermax) THEN
      pbiwmm = (1.0_wp - ziguess) * (1.0_wp - zjguess)
      pbiwmp = (1.0_wp - ziguess) * zjguess
      pbiwpm = ziguess * (1.0_wp - zjguess)
      pbiwpp = ziguess * zjguess
    ELSE IF (jiter > itermax) THEN
      IF (lwp) THEN
        WRITE(numout, FMT = *) 'Obs lat/lon  : ', pphi, plam
        WRITE(numout, FMT = *) 'Grid lats    : ', pphimm, pphimp, pphipm, pphipp
        WRITE(numout, FMT = *) 'Grid lons    : ', plammm, plammp, plampm, plampp
        WRITE(numout, FMT = *) 'Current i,j  : ', ziguess, zjguess
        WRITE(numout, FMT = *) 'jiter        = ', jiter
        WRITE(numout, FMT = *) 'zeps         = ', zeps
        WRITE(numout, FMT = *) 'zdeli, zdelj = ', zdeli, zdelj
        WRITE(numout, FMT = *) ' Iterations for i,j exceed max iteration count!'
        WRITE(numout, FMT = *)
        ldfail = .TRUE.
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE bil_wgt
  SUBROUTINE lu_invmat(pmatin, kdim, pmatou)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kdim
    REAL(KIND = wp), DIMENSION(kdim, kdim), INTENT(IN) :: pmatin
    REAL(KIND = wp), DIMENSION(kdim, kdim), INTENT(OUT) :: pmatou
    INTEGER :: ji, jj
    INTEGER, DIMENSION(kdim) :: indx
    REAL(KIND = wp), DIMENSION(kdim, kdim) :: zmat
    REAL(KIND = wp) :: zd
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('lu_invmat', 'r0', 0, 0)
    DO jj = 1, kdim
      DO ji = 1, kdim
        pmatou(ji, jj) = 0.0_wp
        zmat(ji, jj) = pmatin(ji, jj)
      END DO
      pmatou(jj, jj) = 1.0_wp
    END DO
    CALL lu_decomp(zmat, kdim, kdim, indx, zd)
    DO jj = 1, kdim
      CALL lu_backsb(zmat, kdim, kdim, indx, pmatou(1, jj))
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE lu_invmat
  SUBROUTINE lu_decomp(pmatin, kdim1, kdim2, kindex, pflt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kdim1, kdim2
    INTEGER, DIMENSION(kdim1), INTENT(OUT) :: kindex
    REAL(KIND = wp), INTENT(OUT) :: pflt
    REAL(KIND = wp), DIMENSION(kdim2, kdim2), INTENT(INOUT) :: pmatin
    INTEGER, PARAMETER :: jpmax = 100
    REAL(KIND = wp), PARAMETER :: pptiny = 1.0E-20_wp
    REAL(KIND = wp), DIMENSION(jpmax) :: zvv
    INTEGER :: ji, jj, jk
    INTEGER :: imax
    REAL(KIND = wp) :: zsum, zdum, zaamax
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('lu_decomp', 'r0', 0, 0)
    imax = - 1
    pflt = 1.0_wp
    DO ji = 1, kdim1
      zaamax = 0.0_wp
      DO jj = 1, kdim1
        IF (ABS(pmatin(ji, jj)) > zaamax) zaamax = ABS(pmatin(ji, jj))
      END DO
      IF (zaamax == 0.0_wp) THEN
        CALL ctl_stop('singular matrix')
      END IF
      zvv(ji) = 1.0_wp / zaamax
    END DO
    DO jj = 1, kdim1
      DO ji = 1, jj - 1
        zsum = pmatin(ji, jj)
        DO jk = 1, ji - 1
          zsum = zsum - pmatin(ji, jk) * pmatin(jk, jj)
        END DO
        pmatin(ji, jj) = zsum
      END DO
      zaamax = 0.0_wp
      DO ji = jj, kdim1
        zsum = pmatin(ji, jj)
        DO jk = 1, jj - 1
          zsum = zsum - pmatin(ji, jk) * pmatin(jk, jj)
        END DO
        pmatin(ji, jj) = zsum
        zdum = zvv(ji) * ABS(zsum)
        IF (zdum >= zaamax) THEN
          imax = ji
          zaamax = zdum
        END IF
      END DO
      IF (jj /= imax) THEN
        DO jk = 1, kdim1
          zdum = pmatin(imax, jk)
          pmatin(imax, jk) = pmatin(jj, jk)
          pmatin(jj, jk) = zdum
        END DO
        pflt = - pflt
        zvv(imax) = zvv(jj)
      END IF
      kindex(jj) = imax
      IF (pmatin(jj, jj) == 0.0_wp) pmatin(jj, jj) = pptiny
      IF (jj /= kdim1) THEN
        zdum = 1.0_wp / pmatin(jj, jj)
        DO ji = jj + 1, kdim1
          pmatin(ji, jj) = pmatin(ji, jj) * zdum
        END DO
      END IF
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE lu_decomp
  SUBROUTINE lu_backsb(pmat, kdim1, kdim2, kindex, pvect)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kdim1, kdim2
    INTEGER, DIMENSION(kdim1), INTENT(IN) :: kindex
    REAL(KIND = wp), DIMENSION(kdim1), INTENT(INOUT) :: pvect
    REAL(KIND = wp), DIMENSION(kdim2, kdim2), INTENT(IN) :: pmat
    INTEGER :: ji, jii, jj, jll
    REAL(KIND = wp) :: zsum
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('lu_backsb', 'r0', 0, 0)
    jii = 0
    DO ji = 1, kdim1
      jll = kindex(ji)
      zsum = pvect(jll)
      pvect(jll) = pvect(ji)
      IF (jii /= 0) THEN
        DO jj = jii, ji - 1
          zsum = zsum - pmat(ji, jj) * pvect(jj)
        END DO
      ELSE IF (zsum /= 0.0_wp) THEN
        jii = ji
      END IF
      pvect(ji) = zsum
    END DO
    DO ji = kdim1, 1, - 1
      zsum = pvect(ji)
      DO jj = ji + 1, kdim1
        zsum = zsum - pmat(ji, jj) * pvect(jj)
      END DO
      pvect(ji) = zsum / pmat(ji, ji)
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE lu_backsb
END MODULE obs_inter_h2d