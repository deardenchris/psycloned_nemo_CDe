MODULE icealb
  USE ice, ONLY: jpl
  USE phycst
  USE dom_oce
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_alb_init
  PUBLIC :: ice_alb
  REAL(KIND = wp), PUBLIC, PARAMETER :: rn_alb_oce = 0.066
  REAL(KIND = wp) :: rn_alb_sdry
  REAL(KIND = wp) :: rn_alb_smlt
  REAL(KIND = wp) :: rn_alb_idry
  REAL(KIND = wp) :: rn_alb_imlt
  REAL(KIND = wp) :: rn_alb_dpnd
  CONTAINS
  SUBROUTINE ice_alb(pt_su, ph_ice, ph_snw, ld_pnd_alb, pafrac_pnd, ph_pnd, palb_cs, palb_os)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :) :: pt_su
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :) :: ph_ice
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :) :: ph_snw
    LOGICAL, INTENT(IN) :: ld_pnd_alb
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :) :: pafrac_pnd
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :) :: ph_pnd
    REAL(KIND = wp), INTENT(OUT), DIMENSION(:, :, :) :: palb_cs
    REAL(KIND = wp), INTENT(OUT), DIMENSION(:, :, :) :: palb_os
    INTEGER :: ji, jj, jl
    REAL(KIND = wp) :: z1_c1, z1_c2, z1_c3, z1_c4
    REAL(KIND = wp) :: z1_href_pnd
    REAL(KIND = wp) :: zalb_pnd, zafrac_pnd
    REAL(KIND = wp) :: zalb_ice, zafrac_ice
    REAL(KIND = wp) :: zalb_snw, zafrac_snw
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('ice_alb', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('icealb')
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    z1_href_pnd = 1. / 0.05
    z1_c1 = 1. / (LOG(1.5) - LOG(0.05))
    z1_c2 = 1. / 0.05
    z1_c3 = 1. / 0.02
    z1_c4 = 1. / 0.03
    DO jl = 1, jpl
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (ph_snw(ji, jj, jl) == 0._wp) THEN
            zafrac_snw = 0._wp
            IF (ld_pnd_alb) THEN
              zafrac_pnd = pafrac_pnd(ji, jj, jl)
            ELSE
              zafrac_pnd = 0._wp
            END IF
            zafrac_ice = 1._wp - zafrac_pnd
          ELSE
            zafrac_snw = 1._wp
            zafrac_pnd = 0._wp
            zafrac_ice = 0._wp
          END IF
          IF (ld_pnd_alb) THEN
            zalb_ice = rn_alb_idry
          ELSE
            IF (ph_snw(ji, jj, jl) == 0._wp .AND. pt_su(ji, jj, jl) >= rt0) THEN
              zalb_ice = rn_alb_imlt
            ELSE
              zalb_ice = rn_alb_idry
            END IF
          END IF
          IF (0.05 < ph_ice(ji, jj, jl) .AND. ph_ice(ji, jj, jl) <= 1.5) THEN
            zalb_ice = zalb_ice + (0.18 - zalb_ice) * z1_c1 * (LOG(1.5) - LOG(ph_ice(ji, jj, jl)))
          ELSE IF (ph_ice(ji, jj, jl) <= 0.05) THEN
            zalb_ice = rn_alb_oce + (0.18 - rn_alb_oce) * z1_c2 * ph_ice(ji, jj, jl)
          END IF
          IF (pt_su(ji, jj, jl) < rt0) THEN
            zalb_snw = rn_alb_sdry - (rn_alb_sdry - zalb_ice) * EXP(- ph_snw(ji, jj, jl) * z1_c3)
          ELSE
            zalb_snw = rn_alb_smlt - (rn_alb_smlt - zalb_ice) * EXP(- ph_snw(ji, jj, jl) * z1_c4)
          END IF
          IF (ld_pnd_alb) THEN
            zalb_pnd = rn_alb_dpnd - (rn_alb_dpnd - zalb_ice) * EXP(- ph_pnd(ji, jj, jl) * z1_href_pnd)
          ELSE
            zalb_pnd = rn_alb_dpnd
          END IF
          palb_os(ji, jj, jl) = (zafrac_snw * zalb_snw + zafrac_pnd * zalb_pnd + zafrac_ice * zalb_ice) * tmask(ji, jj, 1)
          palb_cs(ji, jj, jl) = palb_os(ji, jj, jl) - (- 0.1010 * palb_os(ji, jj, jl) * palb_os(ji, jj, jl) + 0.1933 * palb_os(ji, jj, jl) - 0.0148) * tmask(ji, jj, 1)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('ice_alb', 'r1', 0, 0)
    IF (ln_timing) CALL timing_stop('icealb')
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE ice_alb
  SUBROUTINE ice_alb_init
    INTEGER :: ios
    NAMELIST /namalb/ rn_alb_sdry, rn_alb_smlt, rn_alb_idry, rn_alb_imlt, rn_alb_dpnd
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, namalb, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namalb in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, namalb, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namalb in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, namalb)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ice_alb_init: set albedo parameters'
      WRITE(numout, *) '~~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist namalb:'
      WRITE(numout, *) '      albedo of dry snow                   rn_alb_sdry = ', rn_alb_sdry
      WRITE(numout, *) '      albedo of melting snow               rn_alb_smlt = ', rn_alb_smlt
      WRITE(numout, *) '      albedo of dry ice                    rn_alb_idry = ', rn_alb_idry
      WRITE(numout, *) '      albedo of bare puddled ice           rn_alb_imlt = ', rn_alb_imlt
      WRITE(numout, *) '      albedo of ponded ice                 rn_alb_dpnd = ', rn_alb_dpnd
    END IF
  END SUBROUTINE ice_alb_init
END MODULE icealb