MODULE usrdef_hgr
  USE dom_oce, ONLY: nimpp, njmpp
  USE par_oce
  USE phycst
  USE usrdef_nam
  USE in_out_manager
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: usr_def_hgr
  CONTAINS
  SUBROUTINE usr_def_hgr(plamt, plamu, plamv, plamf, pphit, pphiu, pphiv, pphif, kff, pff_f, pff_t, pe1t, pe1u, pe1v, pe1f, pe2t, pe2u, pe2v, pe2f, ke1e2u_v, pe1e2u, pe1e2v)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: plamt, plamu, plamv, plamf
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: pphit, pphiu, pphiv, pphif
    INTEGER, INTENT(OUT) :: kff
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: pff_f, pff_t
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: pe1t, pe1u, pe1v, pe1f
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: pe2t, pe2u, pe2v, pe2f
    INTEGER, INTENT(OUT) :: ke1e2u_v
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: pe1e2u, pe1e2v
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zlam1, zlam0, zcos_alpha, zim1, zjm1, ze1, ze1deg, zf0
    REAL(KIND = wp) :: zphi1, zphi0, zsin_alpha, zim05, zjm05, zbeta, znorme
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('usr_def_hgr', 'r0', 0, 0)
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) 'usr_def_hgr : GYRE configuration (beta-plane with rotated regular grid-spacing)'
    IF (lwp) WRITE(numout, *) '~~~~~~~~~~~'
    zlam1 = - 85._wp
    zphi1 = 29._wp
    ze1 = 106000._wp / REAL(nn_GYRE, wp)
    zsin_alpha = - SQRT(2._wp) * 0.5_wp
    zcos_alpha = SQRT(2._wp) * 0.5_wp
    ze1deg = ze1 / (ra * rad)
    zlam0 = zlam1 + zcos_alpha * ze1deg * REAL(jpjglo - 2, wp)
    zphi0 = zphi1 + zsin_alpha * ze1deg * REAL(jpjglo - 2, wp)
    IF (ln_bench) THEN
      ze1 = 106000._wp
      CALL ctl_warn(' GYRE used as Benchmark: e1=e2=106km, no need to adjust rdt, ahm,aht ')
    END IF
    IF (nprint == 1 .AND. lwp) THEN
      WRITE(numout, *) 'ze1', ze1, 'cosalpha', zcos_alpha, 'sinalpha', zsin_alpha
      WRITE(numout, *) 'ze1deg', ze1deg, 'zlam0', zlam0, 'zphi0', zphi0
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zim1 = REAL(ji + nimpp - 1) - 1.
        zim05 = REAL(ji + nimpp - 1) - 1.5
        zjm1 = REAL(jj + njmpp - 1) - 1.
        zjm05 = REAL(jj + njmpp - 1) - 1.5
        plamt(ji, jj) = zlam0 + zim05 * ze1deg * zcos_alpha + zjm05 * ze1deg * zsin_alpha
        pphit(ji, jj) = zphi0 - zim05 * ze1deg * zsin_alpha + zjm05 * ze1deg * zcos_alpha
        plamu(ji, jj) = zlam0 + zim1 * ze1deg * zcos_alpha + zjm05 * ze1deg * zsin_alpha
        pphiu(ji, jj) = zphi0 - zim1 * ze1deg * zsin_alpha + zjm05 * ze1deg * zcos_alpha
        plamv(ji, jj) = zlam0 + zim05 * ze1deg * zcos_alpha + zjm1 * ze1deg * zsin_alpha
        pphiv(ji, jj) = zphi0 - zim05 * ze1deg * zsin_alpha + zjm1 * ze1deg * zcos_alpha
        plamf(ji, jj) = zlam0 + zim1 * ze1deg * zcos_alpha + zjm1 * ze1deg * zsin_alpha
        pphif(ji, jj) = zphi0 - zim1 * ze1deg * zsin_alpha + zjm1 * ze1deg * zcos_alpha
      END DO
    END DO
    pe1t(:, :) = ze1
    pe2t(:, :) = ze1
    pe1u(:, :) = ze1
    pe2u(:, :) = ze1
    pe1v(:, :) = ze1
    pe2v(:, :) = ze1
    pe1f(:, :) = ze1
    pe2f(:, :) = ze1
    ke1e2u_v = 0
    pe1e2u(:, :) = 0._wp
    pe1e2v(:, :) = 0._wp
    kff = 1
    zbeta = 2. * omega * COS(rad * zphi1) / ra
    zphi0 = 15._wp
    zf0 = 2. * omega * SIN(rad * zphi0)
    pff_f(:, :) = (zf0 + zbeta * ABS(pphif(:, :) - zphi0) * rad * ra)
    pff_t(:, :) = (zf0 + zbeta * ABS(pphit(:, :) - zphi0) * rad * ra)
    !$ACC END KERNELS
    IF (lwp) WRITE(numout, *) '                           beta-plane used. beta = ', zbeta, ' 1/(s.m)'
  END SUBROUTINE usr_def_hgr
END MODULE usrdef_hgr