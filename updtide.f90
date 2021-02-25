MODULE updtide
  USE oce
  USE dom_oce
  USE in_out_manager
  USE phycst
  USE sbctide
  USE tideini, ONLY: ln_tide_ramp, rdttideramp
  IMPLICIT NONE
  PUBLIC
  PUBLIC :: upd_tide
  CONTAINS
  SUBROUTINE upd_tide(kt, kit, time_offset)
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN), OPTIONAL :: kit
    INTEGER, INTENT(IN), OPTIONAL :: time_offset
    INTEGER :: joffset
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zt, zramp
    REAL(KIND = wp), DIMENSION(nb_harmo) :: zwt
    zt = (kt - kt_tide) * rdt
    joffset = 0
    IF (PRESENT(time_offset)) joffset = time_offset
    IF (PRESENT(kit)) THEN
      zt = zt + (kit + joffset - 1) * rdt / REAL(nn_baro, wp)
    ELSE
      zt = zt + joffset * rdt
    END IF
    zwt(:) = omega_tide(:) * zt
    pot_astro(:, :) = 0._wp
    !$OMP parallel default(shared), private(jk)
    !$OMP do schedule(static)
    DO jk = 1, nb_harmo
      pot_astro(:, :) = pot_astro(:, :) + amp_pot(:, :, jk) * COS(zwt(jk) + phi_pot(:, :, jk))
    END DO
    !$OMP end do
    !$OMP end parallel
    IF (ln_tide_ramp) THEN
      zt = (kt - nit000) * rdt
      IF (PRESENT(kit)) zt = zt + (kit + joffset - 1) * rdt / REAL(nn_baro, wp)
      zramp = MIN(MAX(zt / (rdttideramp * rday), 0._wp), 1._wp)
      pot_astro(:, :) = zramp * pot_astro(:, :)
    END IF
  END SUBROUTINE upd_tide
END MODULE updtide