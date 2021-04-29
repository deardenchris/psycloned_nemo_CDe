MODULE sbctide
  USE oce
  USE dom_oce
  USE phycst
  USE daymod
  USE tideini
  USE in_out_manager
  USE iom
  USE ioipsl
  USE lbclnk
  IMPLICIT NONE
  PUBLIC
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: pot_astro
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: amp_pot, phi_pot
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: amp_load, phi_load
  CONTAINS
  SUBROUTINE sbc_tide(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: jk
    INTEGER :: nsec_day_orig
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    IF (nsec_day == NINT(0.5_wp * rdt) .OR. kt == nit000) THEN
      CALL profile_psy_data0 % PreStart('sbc_tide', 'r0', 0, 0)
      IF (kt == nit000) THEN
        ALLOCATE(amp_pot(jpi, jpj, nb_harmo), phi_pot(jpi, jpj, nb_harmo), pot_astro(jpi, jpj))
        IF (ln_read_load) THEN
          ALLOCATE(amp_load(jpi, jpj, nb_harmo), phi_load(jpi, jpj, nb_harmo))
          CALL tide_init_load
        END IF
      END IF
      CALL profile_psy_data0 % PostEnd
      IF (ln_read_load) THEN
        !$ACC KERNELS
        amp_pot(:, :, :) = amp_load(:, :, :)
        phi_pot(:, :, :) = phi_load(:, :, :)
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        amp_pot(:, :, :) = 0._wp
        phi_pot(:, :, :) = 0._wp
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      pot_astro(:, :) = 0._wp
      nsec_day_orig = nsec_day
      IF (nsec_day /= NINT(0.5_wp * rdt)) THEN
        kt_tide = kt - (nsec_day - 0.5_wp * rdt) / rdt
        nsec_day = NINT(0.5_wp * rdt)
      ELSE
        kt_tide = kt
      END IF
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('sbc_tide', 'r1', 0, 0)
      CALL tide_harmo(omega_tide, v0tide, utide, ftide, ntide, nb_harmo)
      IF (lwp) THEN
        WRITE(numout, *)
        WRITE(numout, *) 'sbc_tide : Update of the components and (re)Init. the potential at kt=', kt
        WRITE(numout, *) '~~~~~~~~ '
        DO jk = 1, nb_harmo
          WRITE(numout, *) Wave(ntide(jk)) % cname_tide, utide(jk), ftide(jk), v0tide(jk), omega_tide(jk)
        END DO
      END IF
      IF (ln_tide_pot) CALL tide_init_potential
      nsec_day = nsec_day_orig
      CALL profile_psy_data1 % PostEnd
    END IF
  END SUBROUTINE sbc_tide
  SUBROUTINE tide_init_potential
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zcons, ztmp1, ztmp2, zlat, zlon, ztmp, zamp, zcs
    DO jk = 1, nb_harmo
      zcons = 0.7_wp * Wave(ntide(jk)) % equitide * ftide(jk)
      DO ji = 1, jpi
        DO jj = 1, jpj
          ztmp1 = ftide(jk) * amp_pot(ji, jj, jk) * COS(phi_pot(ji, jj, jk) + v0tide(jk) + utide(jk))
          ztmp2 = - ftide(jk) * amp_pot(ji, jj, jk) * SIN(phi_pot(ji, jj, jk) + v0tide(jk) + utide(jk))
          zlat = gphit(ji, jj) * rad
          zlon = glamt(ji, jj) * rad
          ztmp = v0tide(jk) + utide(jk) + Wave(ntide(jk)) % nutide * zlon
          IF (wave(ntide(jk)) % nutide == 1) THEN
            zcs = zcons * SIN(2._wp * zlat)
          ELSE IF (wave(ntide(jk)) % nutide == 2) THEN
            zcs = zcons * COS(zlat) ** 2
          ELSE
            zcs = 0._wp
          END IF
          ztmp1 = ztmp1 + zcs * COS(ztmp)
          ztmp2 = ztmp2 - zcs * SIN(ztmp)
          zamp = SQRT(ztmp1 * ztmp1 + ztmp2 * ztmp2)
          amp_pot(ji, jj, jk) = zamp
          phi_pot(ji, jj, jk) = ATAN2(- ztmp2 / MAX(1.E-10_wp, zamp), ztmp1 / MAX(1.E-10_wp, zamp))
        END DO
      END DO
    END DO
  END SUBROUTINE tide_init_potential
  SUBROUTINE tide_init_load
    INTEGER :: inum
    INTEGER :: ji, jj, itide
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztr, zti
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'tide_init_load : Initialization of load potential from file'
      WRITE(numout, *) '~~~~~~~~~~~~~~ '
    END IF
    CALL iom_open(cn_tide_load, inum)
    DO itide = 1, nb_harmo
      CALL iom_get(inum, jpdom_data, TRIM(Wave(ntide(itide)) % cname_tide) // '_z1', ztr(:, :))
      CALL iom_get(inum, jpdom_data, TRIM(Wave(ntide(itide)) % cname_tide) // '_z2', zti(:, :))
      DO ji = 1, jpi
        DO jj = 1, jpj
          amp_load(ji, jj, itide) = SQRT(ztr(ji, jj) ** 2. + zti(ji, jj) ** 2.)
          phi_load(ji, jj, itide) = ATAN2(- zti(ji, jj), ztr(ji, jj))
        END DO
      END DO
    END DO
    CALL iom_close(inum)
  END SUBROUTINE tide_init_load
END MODULE sbctide