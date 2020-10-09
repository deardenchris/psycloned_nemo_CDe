MODULE icbthm
  USE par_oce
  USE dom_oce
  USE in_out_manager
  USE lib_mpp
  USE phycst
  USE sbc_oce
  USE icb_oce
  USE icbutl
  USE icbdia
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: icb_thm
  CONTAINS
  SUBROUTINE icb_thm(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ii, ij
    REAL(KIND = wp) :: zM, zT, zW, zL, zSST, zVol, zLn, zWn, zTn, znVol, zIC, zDn
    REAL(KIND = wp) :: zMv, zMe, zMb, zmelt, zdvo, zdva, zdM, zSs, zdMe, zdMb, zdMv
    REAL(KIND = wp) :: zMnew, zMnew1, zMnew2, zheat_hcflux, zheat_latent, z1_12
    REAL(KIND = wp) :: zMbits, znMbits, zdMbitsE, zdMbitsM, zLbits, zAbits, zMbb
    REAL(KIND = wp) :: zxi, zyj, zff, z1_rday, z1_e1e2, zdt, z1_dt, z1_dt_e1e2
    TYPE(iceberg), POINTER :: this, next
    TYPE(point), POINTER :: pt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_thm', 'r0', 0, 0)
    z1_rday = 1._wp / rday
    z1_12 = 1._wp / 12._wp
    zdt = berg_dt
    z1_dt = 1._wp / zdt
    berg_grid % floating_melt(:, :) = 0._wp
    berg_grid % calving_hflx(:, :) = 0._wp
    this => first_berg
    DO WHILE (ASSOCIATED(this))
      pt => this % current_point
      nknberg = this % number(1)
      CALL icb_utl_interp(pt % xi, pt % e1, pt % uo, pt % ui, pt % ua, pt % ssh_x, pt % yj, pt % e2, pt % vo, pt % vi, pt % va, &
&pt % ssh_y, pt % sst, pt % cn, pt % hi, zff)
      zSST = pt % sst
      zIC = MIN(1._wp, pt % cn + rn_sicn_shift)
      zM = pt % mass
      zT = pt % thickness
      zW = pt % width
      zL = pt % length
      zxi = pt % xi
      zyj = pt % yj
      ii = INT(zxi + 0.5)
      ii = mi1(ii)
      ij = INT(zyj + 0.5)
      ij = mj1(ij)
      zVol = zT * zW * zL
      zdvo = SQRT((pt % uvel - pt % uo) ** 2 + (pt % vvel - pt % vo) ** 2)
      zdva = SQRT((pt % ua - pt % uo) ** 2 + (pt % va - pt % vo) ** 2)
      zSs = 1.5_wp * SQRT(zdva) + 0.1_wp * zdva
      zMv = MAX(7.62D-3 * zSST + 1.29D-3 * (zSST ** 2), 0._wp) * z1_rday
      zMb = MAX(0.58_wp * (zdvo ** 0.8_wp) * (zSST + 4.0_wp) / (zL ** 0.2_wp), 0._wp) * z1_rday
      zMe = MAX(z1_12 * (zSST + 2.) * zSs * (1._wp + COS(rpi * (zIC ** 3))), 0._wp) * z1_rday
      IF (ln_operator_splitting) THEN
        zTn = MAX(zT - zMb * zdt, 0._wp)
        znVol = zTn * zW * zL
        zMnew1 = (znVol / zVol) * zM
        zdMb = zM - zMnew1
        zLn = MAX(zL - zMv * zdt, 0._wp)
        zWn = MAX(zW - zMv * zdt, 0._wp)
        znVol = zTn * zWn * zLn
        zMnew2 = (znVol / zVol) * zM
        zdMv = zMnew1 - zMnew2
        zLn = MAX(zLn - zMe * zdt, 0._wp)
        zWn = MAX(zWn - zMe * zdt, 0._wp)
        znVol = zTn * zWn * zLn
        zMnew = (znVol / zVol) * zM
        zdMe = zMnew2 - zMnew
        zdM = zM - zMnew
      ELSE
        zLn = MAX(zL - (zMv + zMe) * zdt, 0._wp)
        zWn = MAX(zW - (zMv + zMe) * zdt, 0._wp)
        zTn = MAX(zT - zMb * zdt, 0._wp)
        znVol = zTn * zWn * zLn
        zMnew = (znVol / zVol) * zM
        zdM = zM - zMnew
        zdMb = (zM / zVol) * (zW * zL) * zMb * zdt
        zdMe = (zM / zVol) * (zT * (zW + zL)) * zMe * zdt
        zdMv = (zM / zVol) * (zT * (zW + zL)) * zMv * zdt
      END IF
      IF (rn_bits_erosion_fraction > 0._wp) THEN
        zMbits = pt % mass_of_bits
        zdMbitsE = rn_bits_erosion_fraction * zdMe
        znMbits = zMbits + zdMbitsE
        zLbits = MIN(zL, zW, zT, 40._wp)
        zAbits = (zMbits / rn_rho_bergs) / zLbits
        zMbb = MAX(0.58_wp * (zdvo ** 0.8_wp) * (zSST + 2._wp) / (zLbits ** 0.2_wp), 0._wp) * z1_rday
        zMbb = rn_rho_bergs * zAbits * zMbb
        zdMbitsM = MIN(zMbb * zdt, znMbits)
        znMbits = znMbits - zdMbitsM
        IF (zMnew == 0._wp) THEN
          zdMbitsM = zdMbitsM + znMbits
          znMbits = 0._wp
        END IF
      ELSE
        zAbits = 0._wp
        zdMbitsE = 0._wp
        zdMbitsM = 0._wp
        znMbits = pt % mass_of_bits
      END IF
      IF (tmask(ii, ij, 1) /= 0._wp) THEN
        z1_e1e2 = r1_e1e2t(ii, ij) * this % mass_scaling
        z1_dt_e1e2 = z1_dt * z1_e1e2
        zmelt = (zdM - (zdMbitsE - zdMbitsM)) * z1_dt
        berg_grid % floating_melt(ii, ij) = berg_grid % floating_melt(ii, ij) + zmelt * z1_e1e2
        zheat_hcflux = zmelt * pt % heat_density
        zheat_latent = - zmelt * rLfus
        berg_grid % calving_hflx(ii, ij) = berg_grid % calving_hflx(ii, ij) + (zheat_hcflux + zheat_latent) * z1_e1e2
        CALL icb_dia_melt(ii, ij, zMnew, zheat_hcflux, zheat_latent, this % mass_scaling, zdM, zdMbitsE, zdMbitsM, zdMb, zdMe, &
&zdMv, z1_dt_e1e2)
      ELSE
        WRITE(numout, FMT = *) 'icb_thm: berg ', this % number(:), ' appears to have grounded  at ', narea, ii, ij
        CALL icb_utl_print_berg(this, kt)
        WRITE(numout, FMT = *) 'msk=', tmask(ii, ij, 1), e1e2t(ii, ij)
        CALL ctl_stop('icb_thm', 'berg appears to have grounded!')
      END IF
      zDn = (rn_rho_bergs / pp_rho_seawater) * zTn
      IF (zDn > 0._wp .AND. MAX(zWn, zLn) < SQRT(0.92 * (zDn ** 2) + 58.32 * zDn)) THEN
        zT = zTn
        zTn = zWn
        zWn = zT
      END IF
      pt % mass = zMnew
      pt % mass_of_bits = znMbits
      pt % thickness = zTn
      pt % width = MIN(zWn, zLn)
      pt % length = MAX(zWn, zLn)
      next => this % next
      IF (zMnew <= 0._wp) THEN
        CALL icb_utl_delete(first_berg, this)
      ELSE
        z1_e1e2 = r1_e1e2t(ii, ij) * this % mass_scaling
        CALL icb_dia_size(ii, ij, zWn, zLn, zAbits, this % mass_scaling, zMnew, znMbits, z1_e1e2)
      END IF
      this => next
    END DO
    IF (.NOT. ln_passive_mode) THEN
      emp(:, :) = emp(:, :) - berg_grid % floating_melt(:, :)
      qns(:, :) = qns(:, :) + berg_grid % calving_hflx(:, :)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_thm
END MODULE icbthm