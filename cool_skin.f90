MODULE cool_skin
  USE par_kind
  USE phycst
  USE dom_oce
  USE in_out_manager
  USE sbc_oce
  USE lib_mpp
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  REAL(KIND = wp), PRIVATE, PARAMETER :: pp_k = 0.596_wp
  REAL(KIND = wp), PRIVATE, PARAMETER :: pp_v = 1.05E-6_wp
  REAL(KIND = wp), PRIVATE, PARAMETER :: pp_C = 86400
  REAL(KIND = wp), PRIVATE, PARAMETER :: pp_cw = 3993._wp
  REAL(KIND = wp), PRIVATE, PARAMETER :: pp_h = 10._wp
  REAL(KIND = wp), PRIVATE, PARAMETER :: pp_rhoa = 1.20421_wp
  REAL(KIND = wp), PRIVATE, PARAMETER :: pp_cda = 1.45E-3_wp
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: x_csdsst
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: x_csthick
  PUBLIC :: diurnal_sst_coolskin_step, diurnal_sst_coolskin_init
  CONTAINS
  SUBROUTINE diurnal_sst_coolskin_init
    ALLOCATE(x_csdsst(jpi, jpj), x_csthick(jpi, jpj))
    x_csdsst = 0.
    x_csthick = 0.
  END SUBROUTINE diurnal_sst_coolskin_init
  SUBROUTINE diurnal_sst_coolskin_step(psqflux, pstauflux, psrho, rdt)
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj) :: psqflux
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj) :: pstauflux
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj) :: psrho
    REAL(KIND = wp), INTENT(IN) :: rdt
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z_fv
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z_gamma
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z_lamda
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z_wspd
    REAL(KIND = wp) :: z_ztx
    REAL(KIND = wp) :: z_zty
    REAL(KIND = wp) :: z_zmod
    INTEGER :: ji, jj
    IF (.NOT. ln_blk) CALL ctl_stop("cool_skin.f90: diurnal flux processing only implemented for bulk forcing")
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (tmask(ji, jj, 1) == 1. .AND. pstauflux(ji, jj) /= 0 .AND. psrho(ji, jj) /= 0) THEN
          z_fv(ji, jj) = SQRT(pstauflux(ji, jj) / psrho(ji, jj))
          z_wspd(ji, jj) = SQRT(pstauflux(ji, jj) / (pp_cda * pp_rhoa))
        ELSE
          z_fv(ji, jj) = 0.
          z_wspd(ji, jj) = 0.
        END IF
        IF (tmask(ji, jj, 1) == 1.) THEN
          IF ((z_wspd(ji, jj) <= 7.5)) z_gamma(ji, jj) = (0.2 * z_wspd(ji, jj)) + 0.5
          IF ((z_wspd(ji, jj) > 7.5) .AND. (z_wspd(ji, jj) < 10.)) z_gamma(ji, jj) = (1.6 * z_wspd(ji, jj)) - 10.
          IF ((z_wspd(ji, jj) >= 10.)) z_gamma(ji, jj) = 6.
        END IF
        IF (tmask(ji, jj, 1) == 1. .AND. z_fv(ji, jj) /= 0) THEN
          z_lamda(ji, jj) = (z_fv(ji, jj) * pp_k * pp_C) / (z_gamma(ji, jj) * psrho(ji, jj) * pp_cw * pp_h * pp_v)
        ELSE
          z_lamda(ji, jj) = 0.
        END IF
        IF (tmask(ji, jj, 1) == 1. .AND. z_fv(ji, jj) /= 0 .AND. psqflux(ji, jj) < 0) THEN
          x_csthick(ji, jj) = (z_lamda(ji, jj) * pp_v) / z_fv(ji, jj)
        ELSE
          x_csthick(ji, jj) = 0.
        END IF
        IF (tmask(ji, jj, 1) == 1. .AND. x_csthick(ji, jj) /= 0. .AND. psqflux(ji, jj) < 0.) THEN
          x_csdsst(ji, jj) = (psqflux(ji, jj) * x_csthick(ji, jj)) / pp_k
        ELSE
          x_csdsst(ji, jj) = 0.
        END IF
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE diurnal_sst_coolskin_step
END MODULE cool_skin