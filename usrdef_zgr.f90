MODULE usrdef_zgr
  USE oce
  USE dom_oce
  USE depth_e3
  USE in_out_manager
  USE lbclnk
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: usr_def_zgr
  CONTAINS
  SUBROUTINE usr_def_zgr(ld_zco, ld_zps, ld_sco, ld_isfcav, pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d, pdept, pdepw, pe3t, pe3u, pe3v, &
&pe3f, pe3w, pe3uw, pe3vw, k_top, k_bot)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    LOGICAL, INTENT(OUT) :: ld_zco, ld_zps, ld_sco
    LOGICAL, INTENT(OUT) :: ld_isfcav
    REAL(KIND = wp), DIMENSION(:), INTENT(OUT) :: pdept_1d, pdepw_1d
    REAL(KIND = wp), DIMENSION(:), INTENT(OUT) :: pe3t_1d, pe3w_1d
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: pdept, pdepw
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: pe3t, pe3u, pe3v, pe3f
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: pe3w, pe3uw, pe3vw
    INTEGER, DIMENSION(:, :), INTENT(OUT) :: k_top, k_bot
    INTEGER :: inum
    REAL(KIND = WP) :: z_zco, z_zps, z_sco, z_cav
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z2d
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('usr_def_zgr', 'r0', 0, 0)
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) 'usr_def_zgr : GYRE configuration (z-coordinate closed flat box ocean without cavities)'
    IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~'
    ld_zco = .TRUE.
    ld_zps = .FALSE.
    ld_sco = .FALSE.
    ld_isfcav = .FALSE.
    CALL zgr_z(pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d)
    CALL zgr_msk_top_bot(k_top, k_bot)
    CALL zgr_zco(pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d, pdept, pdepw, pe3t, pe3u, pe3v, pe3f, pe3w, pe3uw, pe3vw)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE usr_def_zgr
  SUBROUTINE zgr_z(pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:), INTENT(OUT) :: pdept_1d, pdepw_1d
    REAL(KIND = wp), DIMENSION(:), INTENT(OUT) :: pe3t_1d, pe3w_1d
    INTEGER :: jk
    REAL(KIND = wp) :: zt, zw
    REAL(KIND = wp) :: zsur, za0, za1, zkth, zacr
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('zgr_z', 'r0', 0, 0)
    zsur = - 2033.194295283385_wp
    za0 = 155.8325369664153_wp
    za1 = 146.3615918601890_wp
    zkth = 17.28520372419791_wp
    zacr = 5.0_wp
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) '    zgr_z   : Reference vertical z-coordinates '
      WRITE(numout, FMT = *) '    ~~~~~~~'
      WRITE(numout, FMT = *) '       GYRE case : MI96 function with the following coefficients :'
      WRITE(numout, FMT = *) '                 zsur = ', zsur
      WRITE(numout, FMT = *) '                 za0  = ', za0
      WRITE(numout, FMT = *) '                 za1  = ', za1
      WRITE(numout, FMT = *) '                 zkth = ', zkth
      WRITE(numout, FMT = *) '                 zacr = ', zacr
    END IF
    DO jk = 1, jpk
      zw = REAL(jk, wp)
      zt = REAL(jk, wp) + 0.5_wp
      pdepw_1d(jk) = (zsur + za0 * zw + za1 * zacr * LOG(COSH((zw - zkth) / zacr)))
      pdept_1d(jk) = (zsur + za0 * zt + za1 * zacr * LOG(COSH((zt - zkth) / zacr)))
    END DO
    CALL depth_to_e3(pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d)
    CALL e3_to_depth(pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) '              Reference 1D z-coordinate depth and scale factors:'
      WRITE(numout, FMT = "(9x,' level  gdept_1d  gdepw_1d  e3t_1d   e3w_1d  ')")
      WRITE(numout, FMT = "(10x, i4, 4f9.2)") (jk, pdept_1d(jk), pdepw_1d(jk), pe3t_1d(jk), pe3w_1d(jk), jk = 1, jpk)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE zgr_z
  SUBROUTINE zgr_msk_top_bot(k_top, k_bot)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, DIMENSION(:, :), INTENT(OUT) :: k_top, k_bot
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z2d
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('zgr_msk_top_bot', 'r0', 0, 0)
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) '    zgr_top_bot : defines the top and bottom wet ocean levels.'
    IF (lwp) WRITE(numout, FMT = *) '    ~~~~~~~~~~~'
    IF (lwp) WRITE(numout, FMT = *) '       GYRE case : closed flat box ocean without ocean cavities'
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    z2d(:, :) = REAL(jpkm1, wp)
    !$ACC END KERNELS
    CALL lbc_lnk('usrdef_zgr', z2d, 'T', 1.)
    !$ACC KERNELS
    k_bot(:, :) = NINT(z2d(:, :))
    k_top(:, :) = MIN(1, k_bot(:, :))
    !$ACC END KERNELS
  END SUBROUTINE zgr_msk_top_bot
  SUBROUTINE zgr_zco(pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d, pdept, pdepw, pe3t, pe3u, pe3v, pe3f, pe3w, pe3uw, pe3vw)
    REAL(KIND = wp), DIMENSION(:), INTENT(IN) :: pdept_1d, pdepw_1d
    REAL(KIND = wp), DIMENSION(:), INTENT(IN) :: pe3t_1d, pe3w_1d
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: pdept, pdepw
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: pe3t, pe3u, pe3v, pe3f
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: pe3w, pe3uw, pe3vw
    INTEGER :: jk
    DO jk = 1, jpk
      !$ACC KERNELS
      pdept(:, :, jk) = pdept_1d(jk)
      pdepw(:, :, jk) = pdepw_1d(jk)
      pe3t(:, :, jk) = pe3t_1d(jk)
      pe3u(:, :, jk) = pe3t_1d(jk)
      pe3v(:, :, jk) = pe3t_1d(jk)
      pe3f(:, :, jk) = pe3t_1d(jk)
      pe3w(:, :, jk) = pe3w_1d(jk)
      pe3uw(:, :, jk) = pe3w_1d(jk)
      pe3vw(:, :, jk) = pe3w_1d(jk)
      !$ACC END KERNELS
    END DO
  END SUBROUTINE zgr_zco
END MODULE usrdef_zgr