MODULE bdytides
  USE oce
  USE dom_oce
  USE phycst
  USE bdy_oce
  USE tideini
  USE daymod
  USE in_out_manager
  USE iom
  USE fldread
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: bdytide_init
  PUBLIC :: bdytide_update
  PUBLIC :: bdy_dta_tides
  TYPE, PUBLIC :: TIDES_DATA
    REAL(KIND = wp), POINTER, DIMENSION(:, :, :) :: ssh0
    REAL(KIND = wp), POINTER, DIMENSION(:, :, :) :: u0, v0
    REAL(KIND = wp), POINTER, DIMENSION(:, :, :) :: ssh
    REAL(KIND = wp), POINTER, DIMENSION(:, :, :) :: u, v
  END TYPE TIDES_DATA
  TYPE(TIDES_DATA), PUBLIC, DIMENSION(jp_bdy), TARGET :: tides
  TYPE(OBC_DATA), PUBLIC, DIMENSION(jp_bdy) :: dta_bdy_s
  CONTAINS
  SUBROUTINE bdytide_init
    CHARACTER(LEN = 80) :: filtide
    LOGICAL :: ln_bdytide_2ddta
    LOGICAL :: ln_bdytide_conj
    INTEGER :: ib_bdy, itide, ib
    INTEGER :: ii, ij
    INTEGER :: inum, igrd
    INTEGER, DIMENSION(3) :: ilen0
    INTEGER, POINTER, DIMENSION(:) :: nblen, nblenrim
    INTEGER :: ios
    CHARACTER(LEN = 80) :: clfile
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: dta_read
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: ztr, zti
    TYPE(TIDES_DATA), POINTER :: td
    TYPE(MAP_POINTER), DIMENSION(jpbgrd) :: ibmap_ptr
    NAMELIST /nambdy_tide/ filtide, ln_bdytide_2ddta, ln_bdytide_conj
    IF (nb_bdy > 0) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'bdytide_init : initialization of tidal harmonic forcing at open boundaries'
      IF (lwp) WRITE(numout, *) '~~~~~~~~~~~~'
    END IF
    REWIND(UNIT = numnam_cfg)
    DO ib_bdy = 1, nb_bdy
      IF (nn_dyn2d_dta(ib_bdy) >= 2) THEN
        td => tides(ib_bdy)
        nblen => idx_bdy(ib_bdy) % nblen
        nblenrim => idx_bdy(ib_bdy) % nblenrim
        filtide(:) = ''
        READ(numnam_ref, nambdy_tide, IOSTAT = ios, ERR = 901)
901     IF (ios /= 0) CALL ctl_nam(ios, 'nambdy_tide in reference namelist', lwp)
        READ(numnam_cfg, nambdy_tide, IOSTAT = ios, ERR = 902)
902     IF (ios > 0) CALL ctl_nam(ios, 'nambdy_tide in configuration namelist', lwp)
        IF (lwm) WRITE(numond, nambdy_tide)
        IF (lwp) WRITE(numout, *) '  '
        IF (lwp) WRITE(numout, *) '          Namelist nambdy_tide : tidal harmonic forcing at open boundaries'
        IF (lwp) WRITE(numout, *) '             read tidal data in 2d files: ', ln_bdytide_2ddta
        IF (lwp) WRITE(numout, *) '             assume complex conjugate   : ', ln_bdytide_conj
        IF (lwp) WRITE(numout, *) '             Number of tidal components to read: ', nb_harmo
        IF (lwp) THEN
          WRITE(numout, *) '             Tidal components: '
          DO itide = 1, nb_harmo
            WRITE(numout, *) '                 ', Wave(ntide(itide)) % cname_tide
          END DO
        END IF
        IF (lwp) WRITE(numout, *) ' '
        IF (cn_dyn2d(ib_bdy) == 'frs') THEN
          ilen0(:) = nblen(:)
        ELSE
          ilen0(:) = nblenrim(:)
        END IF
        ALLOCATE(td % ssh0(ilen0(1), nb_harmo, 2))
        ALLOCATE(td % ssh(ilen0(1), nb_harmo, 2))
        ALLOCATE(td % u0(ilen0(2), nb_harmo, 2))
        ALLOCATE(td % u(ilen0(2), nb_harmo, 2))
        ALLOCATE(td % v0(ilen0(3), nb_harmo, 2))
        ALLOCATE(td % v(ilen0(3), nb_harmo, 2))
        td % ssh0(:, :, :) = 0._wp
        td % ssh(:, :, :) = 0._wp
        td % u0(:, :, :) = 0._wp
        td % u(:, :, :) = 0._wp
        td % v0(:, :, :) = 0._wp
        td % v(:, :, :) = 0._wp
        IF (ln_bdytide_2ddta) THEN
          ALLOCATE(zti(jpi, jpj), ztr(jpi, jpj))
          clfile = TRIM(filtide) // '_grid_T.nc'
          CALL iom_open(clfile, inum)
          igrd = 1
          DO itide = 1, nb_harmo
            CALL iom_get(inum, jpdom_autoglo, TRIM(Wave(ntide(itide)) % cname_tide) // '_z1', ztr(:, :))
            CALL iom_get(inum, jpdom_autoglo, TRIM(Wave(ntide(itide)) % cname_tide) // '_z2', zti(:, :))
            DO ib = 1, ilen0(igrd)
              ii = idx_bdy(ib_bdy) % nbi(ib, igrd)
              ij = idx_bdy(ib_bdy) % nbj(ib, igrd)
              td % ssh0(ib, itide, 1) = ztr(ii, ij)
              td % ssh0(ib, itide, 2) = zti(ii, ij)
            END DO
          END DO
          CALL iom_close(inum)
          clfile = TRIM(filtide) // '_grid_U.nc'
          CALL iom_open(clfile, inum)
          igrd = 2
          DO itide = 1, nb_harmo
            CALL iom_get(inum, jpdom_autoglo, TRIM(Wave(ntide(itide)) % cname_tide) // '_u1', ztr(:, :))
            CALL iom_get(inum, jpdom_autoglo, TRIM(Wave(ntide(itide)) % cname_tide) // '_u2', zti(:, :))
            DO ib = 1, ilen0(igrd)
              ii = idx_bdy(ib_bdy) % nbi(ib, igrd)
              ij = idx_bdy(ib_bdy) % nbj(ib, igrd)
              td % u0(ib, itide, 1) = ztr(ii, ij)
              td % u0(ib, itide, 2) = zti(ii, ij)
            END DO
          END DO
          CALL iom_close(inum)
          clfile = TRIM(filtide) // '_grid_V.nc'
          CALL iom_open(clfile, inum)
          igrd = 3
          DO itide = 1, nb_harmo
            CALL iom_get(inum, jpdom_autoglo, TRIM(Wave(ntide(itide)) % cname_tide) // '_v1', ztr(:, :))
            CALL iom_get(inum, jpdom_autoglo, TRIM(Wave(ntide(itide)) % cname_tide) // '_v2', zti(:, :))
            DO ib = 1, ilen0(igrd)
              ii = idx_bdy(ib_bdy) % nbi(ib, igrd)
              ij = idx_bdy(ib_bdy) % nbj(ib, igrd)
              td % v0(ib, itide, 1) = ztr(ii, ij)
              td % v0(ib, itide, 2) = zti(ii, ij)
            END DO
          END DO
          CALL iom_close(inum)
          DEALLOCATE(ztr, zti)
        ELSE
          ALLOCATE(dta_read(MAXVAL(ilen0(1 : 3)), 1, 1))
          ibmap_ptr(1) % ptr => idx_bdy(ib_bdy) % nbmap(:, 1)
          ibmap_ptr(1) % ll_unstruc = ln_coords_file(ib_bdy)
          ibmap_ptr(2) % ptr => idx_bdy(ib_bdy) % nbmap(:, 2)
          ibmap_ptr(2) % ll_unstruc = ln_coords_file(ib_bdy)
          ibmap_ptr(3) % ptr => idx_bdy(ib_bdy) % nbmap(:, 3)
          ibmap_ptr(3) % ll_unstruc = ln_coords_file(ib_bdy)
          DO itide = 1, nb_harmo
            clfile = TRIM(filtide) // TRIM(Wave(ntide(itide)) % cname_tide) // '_grid_T.nc'
            CALL iom_open(clfile, inum)
            CALL fld_map(inum, 'z1', dta_read(1 : ilen0(1), 1 : 1, 1 : 1), 1, ibmap_ptr(1))
            td % ssh0(:, itide, 1) = dta_read(1 : ilen0(1), 1, 1)
            CALL fld_map(inum, 'z2', dta_read(1 : ilen0(1), 1 : 1, 1 : 1), 1, ibmap_ptr(1))
            td % ssh0(:, itide, 2) = dta_read(1 : ilen0(1), 1, 1)
            CALL iom_close(inum)
            clfile = TRIM(filtide) // TRIM(Wave(ntide(itide)) % cname_tide) // '_grid_U.nc'
            CALL iom_open(clfile, inum)
            CALL fld_map(inum, 'u1', dta_read(1 : ilen0(2), 1 : 1, 1 : 1), 1, ibmap_ptr(2))
            td % u0(:, itide, 1) = dta_read(1 : ilen0(2), 1, 1)
            CALL fld_map(inum, 'u2', dta_read(1 : ilen0(2), 1 : 1, 1 : 1), 1, ibmap_ptr(2))
            td % u0(:, itide, 2) = dta_read(1 : ilen0(2), 1, 1)
            CALL iom_close(inum)
            clfile = TRIM(filtide) // TRIM(Wave(ntide(itide)) % cname_tide) // '_grid_V.nc'
            CALL iom_open(clfile, inum)
            CALL fld_map(inum, 'v1', dta_read(1 : ilen0(3), 1 : 1, 1 : 1), 1, ibmap_ptr(3))
            td % v0(:, itide, 1) = dta_read(1 : ilen0(3), 1, 1)
            CALL fld_map(inum, 'v2', dta_read(1 : ilen0(3), 1 : 1, 1 : 1), 1, ibmap_ptr(3))
            td % v0(:, itide, 2) = dta_read(1 : ilen0(3), 1, 1)
            CALL iom_close(inum)
          END DO
          DEALLOCATE(dta_read)
        END IF
        IF (ln_bdytide_conj) THEN
          td % ssh0(:, :, 2) = - td % ssh0(:, :, 2)
          td % u0(:, :, 2) = - td % u0(:, :, 2)
          td % v0(:, :, 2) = - td % v0(:, :, 2)
        END IF
        ALLOCATE(dta_bdy_s(ib_bdy) % ssh(ilen0(1)))
        ALLOCATE(dta_bdy_s(ib_bdy) % u2d(ilen0(2)))
        ALLOCATE(dta_bdy_s(ib_bdy) % v2d(ilen0(3)))
        dta_bdy_s(ib_bdy) % ssh(:) = 0._wp
        dta_bdy_s(ib_bdy) % u2d(:) = 0._wp
        dta_bdy_s(ib_bdy) % v2d(:) = 0._wp
      END IF
    END DO
  END SUBROUTINE bdytide_init
  SUBROUTINE bdytide_update(kt, idx, dta, td, jit, time_offset)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    TYPE(OBC_DATA), INTENT(INOUT) :: dta
    TYPE(TIDES_DATA), INTENT(INOUT) :: td
    INTEGER, OPTIONAL, INTENT(IN) :: jit
    INTEGER, OPTIONAL, INTENT(IN) :: time_offset
    INTEGER :: itide, igrd, ib
    INTEGER :: time_add
    INTEGER, DIMENSION(3) :: ilen0
    REAL(KIND = wp) :: z_arg, z_sarg, zflag, zramp
    REAL(KIND = wp), DIMENSION(jpmax_harmo) :: z_sist, z_cost
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bdytide_update', 'r0', 0, 0)
    ilen0(1) = SIZE(td % ssh(:, 1, 1))
    ilen0(2) = SIZE(td % u(:, 1, 1))
    ilen0(3) = SIZE(td % v(:, 1, 1))
    zflag = 1
    IF (PRESENT(jit)) THEN
      IF (jit /= 1) zflag = 0
    END IF
    IF ((nsec_day == NINT(0.5_wp * rdt) .OR. kt == nit000) .AND. zflag == 1) THEN
      kt_tide = kt - (nsec_day - 0.5_wp * rdt) / rdt
      IF (lwp) THEN
        WRITE(numout, *)
        WRITE(numout, *) 'bdytide_update : (re)Initialization of the tidal bdy forcing at kt=', kt
        WRITE(numout, *) '~~~~~~~~~~~~~~ '
      END IF
      CALL tide_init_elevation(idx, td)
      CALL tide_init_velocities(idx, td)
    END IF
    time_add = 0
    IF (PRESENT(time_offset)) THEN
      time_add = time_offset
    END IF
    IF (PRESENT(jit)) THEN
      z_arg = ((kt - kt_tide) * rdt + (jit + 0.5_wp * (time_add - 1)) * rdt / REAL(nn_baro, wp))
    ELSE
      z_arg = ((kt - kt_tide) + time_add) * rdt
    END IF
    zramp = 1._wp
    IF (ln_tide_ramp) zramp = MIN(MAX((z_arg + (kt_tide - nit000) * rdt) / (rdttideramp * rday), 0._wp), 1._wp)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    DO itide = 1, nb_harmo
      z_sarg = z_arg * omega_tide(itide)
      z_cost(itide) = COS(z_sarg)
      z_sist(itide) = SIN(z_sarg)
    END DO
    DO itide = 1, nb_harmo
      igrd = 1
      DO ib = 1, ilen0(igrd)
        dta % ssh(ib) = dta % ssh(ib) + zramp * (td % ssh(ib, itide, 1) * z_cost(itide) + td % ssh(ib, itide, 2) * z_sist(itide))
      END DO
      igrd = 2
      DO ib = 1, ilen0(igrd)
        dta % u2d(ib) = dta % u2d(ib) + zramp * (td % u(ib, itide, 1) * z_cost(itide) + td % u(ib, itide, 2) * z_sist(itide))
      END DO
      igrd = 3
      DO ib = 1, ilen0(igrd)
        dta % v2d(ib) = dta % v2d(ib) + zramp * (td % v(ib, itide, 1) * z_cost(itide) + td % v(ib, itide, 2) * z_sist(itide))
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE bdytide_update
  SUBROUTINE bdy_dta_tides(kt, kit, time_offset)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, OPTIONAL, INTENT(IN) :: kit
    INTEGER, OPTIONAL, INTENT(IN) :: time_offset
    LOGICAL :: lk_first_btstp
    INTEGER :: itide, ib_bdy, ib, igrd
    INTEGER :: time_add
    INTEGER, DIMENSION(jpbgrd) :: ilen0
    INTEGER, DIMENSION(1 : jpbgrd) :: nblen, nblenrim
    REAL(KIND = wp) :: z_arg, z_sarg, zramp, zoff, z_cost, z_sist
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('bdy_dta_tides', 'r0', 0, 0)
    lk_first_btstp = .TRUE.
    IF (PRESENT(kit) .AND. (kit /= 1)) THEN
      lk_first_btstp = .FALSE.
    END IF
    time_add = 0
    IF (PRESENT(time_offset)) THEN
      time_add = time_offset
    END IF
    IF (PRESENT(kit)) THEN
      z_arg = (kt + (kit + time_add - 1) / REAL(nn_baro, wp)) * rdt
    ELSE
      z_arg = (kt + time_add) * rdt
    END IF
    zramp = 1.
    IF (ln_tide_ramp) zramp = MIN(MAX((z_arg - nit000 * rdt) / (rdttideramp * rday), 0.), 1.)
    CALL profile_psy_data0 % PostEnd
    DO ib_bdy = 1, nb_bdy
      IF (nn_dyn2d_dta(ib_bdy) >= 2) THEN
        !$ACC KERNELS
        nblen(1 : jpbgrd) = idx_bdy(ib_bdy) % nblen(1 : jpbgrd)
        nblenrim(1 : jpbgrd) = idx_bdy(ib_bdy) % nblenrim(1 : jpbgrd)
        !$ACC END KERNELS
        IF (cn_dyn2d(ib_bdy) == 'frs') THEN
          !$ACC KERNELS
          ilen0(:) = nblen(:)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          ilen0(:) = nblenrim(:)
          !$ACC END KERNELS
        END IF
        CALL profile_psy_data1 % PreStart('bdy_dta_tides', 'r1', 0, 0)
        IF ((nsec_day == NINT(0.5_wp * rdt) .OR. kt == nit000) .AND. lk_first_btstp) THEN
          kt_tide = kt - (nsec_day - 0.5_wp * rdt) / rdt
          IF (lwp) THEN
            WRITE(numout, *)
            WRITE(numout, *) 'bdy_tide_dta : Refresh nodal factors for tidal open bdy data at kt=', kt
            WRITE(numout, *) '~~~~~~~~~~~~~~ '
          END IF
          CALL tide_init_elevation(idx = idx_bdy(ib_bdy), td = tides(ib_bdy))
          CALL tide_init_velocities(idx = idx_bdy(ib_bdy), td = tides(ib_bdy))
        END IF
        zoff = - kt_tide * rdt
        IF (PRESENT(kit)) THEN
          IF (dta_bdy(ib_bdy) % ll_ssh) dta_bdy(ib_bdy) % ssh(1 : ilen0(1)) = dta_bdy_s(ib_bdy) % ssh(1 : ilen0(1))
          IF (dta_bdy(ib_bdy) % ll_u2d) dta_bdy(ib_bdy) % u2d(1 : ilen0(2)) = dta_bdy_s(ib_bdy) % u2d(1 : ilen0(2))
          IF (dta_bdy(ib_bdy) % ll_v2d) dta_bdy(ib_bdy) % v2d(1 : ilen0(3)) = dta_bdy_s(ib_bdy) % v2d(1 : ilen0(3))
        END IF
        CALL profile_psy_data1 % PostEnd
        DO itide = 1, nb_harmo
          CALL profile_psy_data2 % PreStart('bdy_dta_tides', 'r2', 0, 0)
          z_sarg = (z_arg + zoff) * omega_tide(itide)
          z_cost = zramp * COS(z_sarg)
          z_sist = zramp * SIN(z_sarg)
          CALL profile_psy_data2 % PostEnd
          IF (dta_bdy(ib_bdy) % ll_ssh) THEN
            !$ACC KERNELS
            igrd = 1
            DO ib = 1, ilen0(igrd)
              dta_bdy(ib_bdy) % ssh(ib) = dta_bdy(ib_bdy) % ssh(ib) + (tides(ib_bdy) % ssh(ib, itide, 1) * z_cost + tides(ib_bdy) % ssh(ib, itide, 2) * z_sist)
            END DO
            !$ACC END KERNELS
          END IF
          IF (dta_bdy(ib_bdy) % ll_u2d) THEN
            !$ACC KERNELS
            igrd = 2
            DO ib = 1, ilen0(igrd)
              dta_bdy(ib_bdy) % u2d(ib) = dta_bdy(ib_bdy) % u2d(ib) + (tides(ib_bdy) % u(ib, itide, 1) * z_cost + tides(ib_bdy) % u(ib, itide, 2) * z_sist)
            END DO
            !$ACC END KERNELS
          END IF
          IF (dta_bdy(ib_bdy) % ll_v2d) THEN
            !$ACC KERNELS
            igrd = 3
            DO ib = 1, ilen0(igrd)
              dta_bdy(ib_bdy) % v2d(ib) = dta_bdy(ib_bdy) % v2d(ib) + (tides(ib_bdy) % v(ib, itide, 1) * z_cost + tides(ib_bdy) % v(ib, itide, 2) * z_sist)
            END DO
            !$ACC END KERNELS
          END IF
        END DO
      END IF
    END DO
  END SUBROUTINE bdy_dta_tides
  SUBROUTINE tide_init_elevation(idx, td)
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    TYPE(TIDES_DATA), INTENT(INOUT) :: td
    INTEGER :: itide, igrd, ib
    INTEGER, DIMENSION(1) :: ilen0
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:) :: mod_tide, phi_tide
    igrd = 1
    ilen0(1) = SIZE(td % ssh0(:, 1, 1))
    ALLOCATE(mod_tide(ilen0(igrd)), phi_tide(ilen0(igrd)))
    DO itide = 1, nb_harmo
      DO ib = 1, ilen0(igrd)
        mod_tide(ib) = SQRT(td % ssh0(ib, itide, 1) ** 2. + td % ssh0(ib, itide, 2) ** 2.)
        phi_tide(ib) = ATAN2(- td % ssh0(ib, itide, 2), td % ssh0(ib, itide, 1))
      END DO
      !$ACC KERNELS
      DO ib = 1, ilen0(igrd)
        mod_tide(ib) = mod_tide(ib) * ftide(itide)
        phi_tide(ib) = phi_tide(ib) + v0tide(itide) + utide(itide)
      END DO
      DO ib = 1, ilen0(igrd)
        td % ssh(ib, itide, 1) = mod_tide(ib) * COS(phi_tide(ib))
        td % ssh(ib, itide, 2) = - mod_tide(ib) * SIN(phi_tide(ib))
      END DO
      !$ACC END KERNELS
    END DO
    DEALLOCATE(mod_tide, phi_tide)
  END SUBROUTINE tide_init_elevation
  SUBROUTINE tide_init_velocities(idx, td)
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    TYPE(TIDES_DATA), INTENT(INOUT) :: td
    INTEGER :: itide, igrd, ib
    INTEGER, DIMENSION(3) :: ilen0
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:) :: mod_tide, phi_tide
    ilen0(2) = SIZE(td % u0(:, 1, 1))
    ilen0(3) = SIZE(td % v0(:, 1, 1))
    igrd = 2
    ALLOCATE(mod_tide(ilen0(igrd)), phi_tide(ilen0(igrd)))
    DO itide = 1, nb_harmo
      DO ib = 1, ilen0(igrd)
        mod_tide(ib) = SQRT(td % u0(ib, itide, 1) ** 2. + td % u0(ib, itide, 2) ** 2.)
        phi_tide(ib) = ATAN2(- td % u0(ib, itide, 2), td % u0(ib, itide, 1))
      END DO
      !$ACC KERNELS
      DO ib = 1, ilen0(igrd)
        mod_tide(ib) = mod_tide(ib) * ftide(itide)
        phi_tide(ib) = phi_tide(ib) + v0tide(itide) + utide(itide)
      END DO
      DO ib = 1, ilen0(igrd)
        td % u(ib, itide, 1) = mod_tide(ib) * COS(phi_tide(ib))
        td % u(ib, itide, 2) = - mod_tide(ib) * SIN(phi_tide(ib))
      END DO
      !$ACC END KERNELS
    END DO
    DEALLOCATE(mod_tide, phi_tide)
    igrd = 3
    ALLOCATE(mod_tide(ilen0(igrd)), phi_tide(ilen0(igrd)))
    DO itide = 1, nb_harmo
      DO ib = 1, ilen0(igrd)
        mod_tide(ib) = SQRT(td % v0(ib, itide, 1) ** 2. + td % v0(ib, itide, 2) ** 2.)
        phi_tide(ib) = ATAN2(- td % v0(ib, itide, 2), td % v0(ib, itide, 1))
      END DO
      !$ACC KERNELS
      DO ib = 1, ilen0(igrd)
        mod_tide(ib) = mod_tide(ib) * ftide(itide)
        phi_tide(ib) = phi_tide(ib) + v0tide(itide) + utide(itide)
      END DO
      DO ib = 1, ilen0(igrd)
        td % v(ib, itide, 1) = mod_tide(ib) * COS(phi_tide(ib))
        td % v(ib, itide, 2) = - mod_tide(ib) * SIN(phi_tide(ib))
      END DO
      !$ACC END KERNELS
    END DO
    DEALLOCATE(mod_tide, phi_tide)
  END SUBROUTINE tide_init_velocities
END MODULE bdytides