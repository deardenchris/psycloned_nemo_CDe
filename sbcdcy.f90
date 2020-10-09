MODULE sbcdcy
  USE oce
  USE phycst
  USE dom_oce
  USE sbc_oce
  USE in_out_manager
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  INTEGER, PUBLIC :: nday_qsr
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: raa, rbb, rcc, rab
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: rtmd, rdawn, rdusk, rscal
  PUBLIC :: sbc_dcy
  CONTAINS
  INTEGER FUNCTION sbc_dcy_alloc()
    ALLOCATE(raa(jpi, jpj), rbb(jpi, jpj), rcc(jpi, jpj), rab(jpi, jpj), rtmd(jpi, jpj), rdawn(jpi, jpj), rdusk(jpi, jpj), &
&rscal(jpi, jpj), STAT = sbc_dcy_alloc)
    CALL mpp_sum('sbcdcy', sbc_dcy_alloc)
    IF (sbc_dcy_alloc /= 0) CALL ctl_stop('STOP', 'sbc_dcy_alloc: failed to allocate arrays')
  END FUNCTION sbc_dcy_alloc
  FUNCTION sbc_dcy(pqsrin, l_mask) RESULT(zqsrout)
    LOGICAL, OPTIONAL, INTENT(IN) :: l_mask
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pqsrin
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zqsrout
    INTEGER :: ji, jj
    INTEGER, DIMENSION(jpi, jpj) :: imask_night
    REAL(KIND = wp) :: ztwopi, zinvtwopi, zconvrad
    REAL(KIND = wp) :: zlo, zup, zlousd, zupusd
    REAL(KIND = wp) :: zdsws, zdecrad, ztx, zsin, zcos
    REAL(KIND = wp) :: ztmp, ztmp1, ztmp2, ztest
    REAL(KIND = wp) :: ztmpm, ztmpm1, ztmpm2
    REAL(KIND = wp) :: fintegral, pt1, pt2, paaa, pbbb, pccc
    fintegral(pt1, pt2, paaa, pbbb, pccc) = paaa * pt2 + zinvtwopi * pbbb * SIN(pccc + ztwopi * pt2) - paaa * pt1 - zinvtwopi * &
&pbbb * SIN(pccc + ztwopi * pt1)
    ztwopi = 2._wp * rpi
    zinvtwopi = 1._wp / ztwopi
    zconvrad = ztwopi / 360._wp
    zlo = (REAL(nsec_day, wp) - 0.5_wp * rdt) / rday
    zup = zlo + (REAL(nn_fsbc, wp) * rdt) / rday
    IF (nday_qsr == - 1) THEN
      IF (lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) 'sbc_dcy : introduce diurnal cycle from daily mean qsr'
        WRITE(numout, FMT = *) '~~~~~~~'
        WRITE(numout, FMT = *)
      END IF
      IF (sbc_dcy_alloc() /= 0) CALL ctl_stop('STOP', 'sbc_dcy_alloc : unable to allocate arrays')
      !$ACC KERNELS
      rcc(:, :) = zconvrad * glamt(:, :) - rpi
      rtmd(:, :) = 0.5_wp - glamt(:, :) / 360._wp
      rtmd(:, :) = MOD((rtmd(:, :) + 1._wp), 1._wp)
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    IF (nday_qsr /= nday) THEN
      nday_qsr = nday
      zdsws = REAL(11 + nday_year, wp)
      zdecrad = (- 23.5_wp * zconvrad) * COS(zdsws * ztwopi / REAL(nyear_len(1), wp))
      zsin = SIN(zdecrad)
      zcos = COS(zdecrad)
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          ztmp = zconvrad * gphit(ji, jj)
          raa(ji, jj) = SIN(ztmp) * zsin
          rbb(ji, jj) = COS(ztmp) * zcos
        END DO
      END DO
      rab(:, :) = - raa(:, :) / rbb(:, :)
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (ABS(rab(ji, jj)) < 1._wp) THEN
            ztx = zinvtwopi * (ACOS(rab(ji, jj)) - rcc(ji, jj))
            ztest = - rbb(ji, jj) * SIN(rcc(ji, jj) + ztwopi * ztx)
            IF (ztest > 0._wp) THEN
              rdawn(ji, jj) = ztx
              rdusk(ji, jj) = rtmd(ji, jj) + (rtmd(ji, jj) - rdawn(ji, jj))
            ELSE
              rdusk(ji, jj) = ztx
              rdawn(ji, jj) = rtmd(ji, jj) - (rdusk(ji, jj) - rtmd(ji, jj))
            END IF
          ELSE
            rdawn(ji, jj) = rtmd(ji, jj) + 0.5_wp
            rdusk(ji, jj) = rdawn(ji, jj)
          END IF
        END DO
      END DO
      rdawn(:, :) = MOD((rdawn(:, :) + 1._wp), 1._wp)
      rdusk(:, :) = MOD((rdusk(:, :) + 1._wp), 1._wp)
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (ABS(rab(ji, jj)) < 1._wp) THEN
            rscal(ji, jj) = 0.0_wp
            IF (rdawn(ji, jj) < rdusk(ji, jj)) THEN
              IF ((rdusk(ji, jj) - rdawn(ji, jj)) .GE. 0.001_wp) THEN
                rscal(ji, jj) = fintegral(rdawn(ji, jj), rdusk(ji, jj), raa(ji, jj), rbb(ji, jj), rcc(ji, jj))
                rscal(ji, jj) = 1._wp / rscal(ji, jj)
              END IF
            ELSE
              IF ((rdusk(ji, jj) + (1._wp - rdawn(ji, jj))) .GE. 0.001_wp) THEN
                rscal(ji, jj) = fintegral(0._wp, rdusk(ji, jj), raa(ji, jj), rbb(ji, jj), rcc(ji, jj)) + fintegral(rdawn(ji, jj), &
&1._wp, raa(ji, jj), rbb(ji, jj), rcc(ji, jj))
                rscal(ji, jj) = 1. / rscal(ji, jj)
              END IF
            END IF
          ELSE
            IF (raa(ji, jj) > rbb(ji, jj)) THEN
              rscal(ji, jj) = fintegral(0._wp, 1._wp, raa(ji, jj), rbb(ji, jj), rcc(ji, jj))
              rscal(ji, jj) = 1._wp / rscal(ji, jj)
            ELSE
              rscal(ji, jj) = 0.0_wp
            END IF
          END IF
        END DO
      END DO
      ztmp = rday / (rdt * REAL(nn_fsbc, wp))
      rscal(:, :) = rscal(:, :) * ztmp
    END IF
    imask_night(:, :) = 0
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        ztmpm = 0._wp
        IF (ABS(rab(ji, jj)) < 1.) THEN
          IF (rdawn(ji, jj) < rdusk(ji, jj)) THEN
            zlousd = MAX(zlo, rdawn(ji, jj))
            zlousd = MIN(zlousd, zup)
            zupusd = MIN(zup, rdusk(ji, jj))
            zupusd = MAX(zupusd, zlo)
            ztmp = fintegral(zlousd, zupusd, raa(ji, jj), rbb(ji, jj), rcc(ji, jj))
            zqsrout(ji, jj) = pqsrin(ji, jj) * ztmp * rscal(ji, jj)
            ztmpm = zupusd - zlousd
            IF (ztmpm .EQ. 0) imask_night(ji, jj) = 1
          ELSE
            zlousd = MIN(zlo, rdusk(ji, jj))
            zupusd = MIN(zup, rdusk(ji, jj))
            ztmp1 = fintegral(zlousd, zupusd, raa(ji, jj), rbb(ji, jj), rcc(ji, jj))
            ztmpm1 = zupusd - zlousd
            zlousd = MAX(zlo, rdawn(ji, jj))
            zupusd = MAX(zup, rdawn(ji, jj))
            ztmp2 = fintegral(zlousd, zupusd, raa(ji, jj), rbb(ji, jj), rcc(ji, jj))
            ztmpm2 = zupusd - zlousd
            ztmp = ztmp1 + ztmp2
            ztmpm = ztmpm1 + ztmpm2
            zqsrout(ji, jj) = pqsrin(ji, jj) * ztmp * rscal(ji, jj)
            IF (ztmpm .EQ. 0.) imask_night(ji, jj) = 1
          END IF
        ELSE
          IF (raa(ji, jj) > rbb(ji, jj)) THEN
            ztmp = fintegral(zlo, zup, raa(ji, jj), rbb(ji, jj), rcc(ji, jj))
            zqsrout(ji, jj) = pqsrin(ji, jj) * ztmp * rscal(ji, jj)
            imask_night(ji, jj) = 0
          ELSE
            zqsrout(ji, jj) = 0.0_wp
            imask_night(ji, jj) = 1
          END IF
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    IF (PRESENT(l_mask)) THEN
      IF (l_mask) THEN
        !$ACC KERNELS
        zqsrout(:, :) = FLOAT(imask_night(:, :))
        !$ACC END KERNELS
      END IF
    END IF
  END FUNCTION sbc_dcy
END MODULE sbcdcy