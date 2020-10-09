MODULE obs_grid
  USE par_kind, ONLY: wp
  USE par_oce, ONLY: jpk, jpni, jpnj, jpnij
  USE dom_oce
  USE obs_mpp, ONLY: obs_mpp_find_obs_proc, mpp_global_max, obs_mpp_max_integer
  USE phycst, ONLY: rad
  USE obs_utils, ONLY: grt_cir_dis, chkerr
  USE in_out_manager
  USE netcdf
  USE obs_const, ONLY: obfillflt
  USE lib_mpp, ONLY: ctl_warn, ctl_stop
  IMPLICIT NONE
  PUBLIC :: obs_grid_setup, obs_grid_search, obs_grid_deallocate, obs_level_search
  PRIVATE :: linquad, maxdist, obs_grd_bruteforce, obs_grd_lookup
  REAL, PUBLIC :: rn_gridsearchres = 0.5
  INTEGER, PRIVATE :: gsearch_nlons_def
  INTEGER, PRIVATE :: gsearch_nlats_def
  REAL(KIND = wp), PRIVATE :: gsearch_lonmin_def
  REAL(KIND = wp), PRIVATE :: gsearch_latmin_def
  REAL(KIND = wp), PRIVATE :: gsearch_dlon_def
  REAL(KIND = wp), PRIVATE :: gsearch_dlat_def
  INTEGER, PRIVATE :: nlons
  INTEGER, PRIVATE :: nlats
  REAL(KIND = wp), PRIVATE :: lonmin
  REAL(KIND = wp), PRIVATE :: latmin
  REAL(KIND = wp), PRIVATE :: dlon
  REAL(KIND = wp), PRIVATE :: dlat
  INTEGER, PRIVATE :: maxxdiff, maxydiff
  INTEGER, PRIVATE :: limxdiff, limydiff
  REAL(KIND = wp), PRIVATE, DIMENSION(:, :), ALLOCATABLE :: lons, lats
  INTEGER, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: ixpos, iypos, iprocn
  LOGICAL, PUBLIC :: ln_grid_search_lookup
  LOGICAL, PUBLIC :: ln_grid_global
  CHARACTER(LEN = 44), PUBLIC :: cn_gridsearchfile
  CONTAINS
  SUBROUTINE obs_grid_search(kobsin, plam, pphi, kobsi, kobsj, kproc, cdgrid)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: kobsin
    REAL(KIND = wp), DIMENSION(kobsin), INTENT(IN) :: plam, pphi
    INTEGER, DIMENSION(kobsin), INTENT(OUT) :: kobsi, kobsj, kproc
    CHARACTER(LEN = 1) :: cdgrid
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_grid_search', 'r0', 0, 0)
    IF (kobsin > 0) THEN
      IF (ln_grid_search_lookup .AND. (cdgrid == 'T')) THEN
        CALL obs_grd_lookup(kobsin, plam, pphi, kobsi, kobsj, kproc)
      ELSE
        IF (cdgrid == 'T') THEN
          CALL obs_grd_bruteforce(jpi, jpj, jpiglo, jpjglo, 1, nlci, 1, nlcj, nproc, jpnij, glamt, gphit, tmask, kobsin, plam, &
&pphi, kobsi, kobsj, kproc)
        ELSE IF (cdgrid == 'U') THEN
          CALL obs_grd_bruteforce(jpi, jpj, jpiglo, jpjglo, 1, nlci, 1, nlcj, nproc, jpnij, glamu, gphiu, umask, kobsin, plam, &
&pphi, kobsi, kobsj, kproc)
        ELSE IF (cdgrid == 'V') THEN
          CALL obs_grd_bruteforce(jpi, jpj, jpiglo, jpjglo, 1, nlci, 1, nlcj, nproc, jpnij, glamv, gphiv, vmask, kobsin, plam, &
&pphi, kobsi, kobsj, kproc)
        ELSE IF (cdgrid == 'F') THEN
          CALL obs_grd_bruteforce(jpi, jpj, jpiglo, jpjglo, 1, nlci, 1, nlcj, nproc, jpnij, glamf, gphif, fmask, kobsin, plam, &
&pphi, kobsi, kobsj, kproc)
        ELSE
          CALL ctl_stop('Grid not supported')
        END IF
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_grid_search
  SUBROUTINE obs_grd_bruteforce(kpi, kpj, kpiglo, kpjglo, kldi, klei, kldj, klej, kmyproc, ktotproc, pglam, pgphi, pmask, kobs, &
&plam, pphi, kobsi, kobsj, kproc)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kpi
    INTEGER, INTENT(IN) :: kpj
    INTEGER, INTENT(IN) :: kpiglo
    INTEGER, INTENT(IN) :: kpjglo
    INTEGER, INTENT(IN) :: kldi
    INTEGER, INTENT(IN) :: klei
    INTEGER, INTENT(IN) :: kldj
    INTEGER, INTENT(IN) :: klej
    INTEGER, INTENT(IN) :: kmyproc
    INTEGER, INTENT(IN) :: ktotproc
    REAL(KIND = wp), DIMENSION(kpi, kpj), INTENT(IN) :: pglam, pgphi, pmask
    INTEGER, INTENT(IN) :: kobs
    REAL(KIND = wp), DIMENSION(kobs), INTENT(IN) :: plam, pphi
    INTEGER, DIMENSION(kobs), INTENT(OUT) :: kobsi, kobsj, kproc
    REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: zplam, zpphi
    REAL(KIND = wp) :: zlammax
    REAL(KIND = wp) :: zlam
    INTEGER :: ji
    INTEGER :: jj
    INTEGER :: jk
    INTEGER :: jo
    INTEGER :: jlon
    INTEGER :: jlat
    INTEGER :: joffset
    INTEGER :: jostride
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: zlamg, zphig, zmskg, zphitmax, zphitmin, zlamtmax, zlamtmin
    LOGICAL, DIMENSION(:, :), ALLOCATABLE :: llinvalidcell
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: zlamtm, zphitm
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_grd_bruteforce', 'r0', 0, 0)
    IF (ln_grid_global) THEN
      jlon = kpiglo
      jlat = kpjglo
      joffset = kmyproc
      jostride = ktotproc
    ELSE
      jlon = kpi
      jlat = kpj
      joffset = 0
      jostride = 1
    END IF
    ALLOCATE(zlamg(jlon, jlat), zphig(jlon, jlat), zmskg(jlon, jlat), zphitmax(jlon - 1, jlat - 1), zphitmin(jlon - 1, jlat - 1), &
&zlamtmax(jlon - 1, jlat - 1), zlamtmin(jlon - 1, jlat - 1), llinvalidcell(jlon - 1, jlat - 1), zlamtm(4, jlon - 1, jlat - 1), &
&zphitm(4, jlon - 1, jlat - 1))
    IF (ln_grid_global) THEN
      zlamg(:, :) = - 1.E+10
      zphig(:, :) = - 1.E+10
      zmskg(:, :) = - 1.E+10
      DO jj = kldj, klej
        DO ji = kldi, klei
          zlamg(mig(ji), mjg(jj)) = pglam(ji, jj)
          zphig(mig(ji), mjg(jj)) = pgphi(ji, jj)
          zmskg(mig(ji), mjg(jj)) = pmask(ji, jj)
        END DO
      END DO
      CALL mpp_global_max(zlamg)
      CALL mpp_global_max(zphig)
      CALL mpp_global_max(zmskg)
    ELSE
      DO jj = 1, jlat
        DO ji = 1, jlon
          zlamg(ji, jj) = pglam(ji, jj)
          zphig(ji, jj) = pgphi(ji, jj)
          zmskg(ji, jj) = pmask(ji, jj)
        END DO
      END DO
    END IF
    ALLOCATE(zplam(kobs), zpphi(kobs))
    DO jo = 1, kobs
      zplam(jo) = plam(jo)
      zpphi(jo) = pphi(jo)
    END DO
    kproc(:) = - 1
    kobsi(:) = - 1
    kobsj(:) = - 1
    DO jj = 1, jlat - 1
      DO ji = 1, jlon - 1
        zlamtm(1, ji, jj) = zlamg(ji, jj)
        zphitm(1, ji, jj) = zphig(ji, jj)
        zlamtm(2, ji, jj) = zlamg(ji + 1, jj)
        zphitm(2, ji, jj) = zphig(ji + 1, jj)
        zlamtm(3, ji, jj) = zlamg(ji + 1, jj + 1)
        zphitm(3, ji, jj) = zphig(ji + 1, jj + 1)
        zlamtm(4, ji, jj) = zlamg(ji, jj + 1)
        zphitm(4, ji, jj) = zphig(ji, jj + 1)
      END DO
    END DO
    WHERE (zlamtm(:, :, :) < 0.0_wp)
      zlamtm(:, :, :) = zlamtm(:, :, :) + 360.0_wp
    END WHERE
    WHERE (zlamtm(:, :, :) > 360.0_wp)
      zlamtm(:, :, :) = zlamtm(:, :, :) - 360.0_wp
    END WHERE
    DO jj = 1, jlat - 1
      DO ji = 1, jlon - 1
        zlammax = MAXVAL(zlamtm(:, ji, jj))
        WHERE (zlammax - zlamtm(:, ji, jj) > 180) zlamtm(:, ji, jj) = zlamtm(:, ji, jj) + 360._wp
        zphitmax(ji, jj) = MAXVAL(zphitm(:, ji, jj))
        zphitmin(ji, jj) = MINVAL(zphitm(:, ji, jj))
        zlamtmax(ji, jj) = MAXVAL(zlamtm(:, ji, jj))
        zlamtmin(ji, jj) = MINVAL(zlamtm(:, ji, jj))
      END DO
    END DO
    llinvalidcell(:, :) = .FALSE.
    DO jj = 1, jlat - 1
      DO ji = 1, jlon - 1
        llinvalidcell(ji, jj) = zmskg(ji, jj) == 0.0_wp .AND. zmskg(ji + 1, jj) == 0.0_wp .AND. zmskg(ji + 1, jj + 1) == 0.0_wp &
&.AND. zmskg(ji, jj + 1) == 0.0_wp
      END DO
    END DO
    DO jo = 1 + joffset, kobs, jostride
      IF (zplam(jo) < 0.0_wp) zplam(jo) = zplam(jo) + 360.0_wp
      IF (zplam(jo) > 360.0_wp) zplam(jo) = zplam(jo) - 360.0_wp
      gridloop:DO jj = 1, jlat - 1
        DO ji = 1, jlon - 1
          IF (ABS(zphig(ji, jj) - zpphi(jo)) < 1E-6) THEN
            zlam = zlamg(ji, jj)
            IF (zlam < 0.0_wp) zlam = zlam + 360.0_wp
            IF (zlam > 360.0_wp) zlam = zlam - 360.0_wp
            IF (ABS(zlam - zplam(jo)) < 1E-6) THEN
              IF (llinvalidcell(ji, jj)) THEN
                kproc(jo) = kmyproc + 1000000
                kobsi(jo) = ji + 1
                kobsj(jo) = jj + 1
                CYCLE
              ELSE
                kproc(jo) = kmyproc
                kobsi(jo) = ji + 1
                kobsj(jo) = jj + 1
                EXIT gridloop
              END IF
            END IF
          END IF
        END DO
      END DO gridloop
      IF (zplam(jo) > 180) zplam(jo) = zplam(jo) - 360.0_wp
      IF (kproc(jo) == - 1) THEN
        gridpoints:DO jj = 1, jlat - 1
          DO ji = 1, jlon - 1
            IF ((zplam(jo) > zlamtmax(ji, jj)) .OR. (zplam(jo) < zlamtmin(ji, jj))) CYCLE
            IF (ABS(zpphi(jo)) < 85) THEN
              IF ((zpphi(jo) > zphitmax(ji, jj)) .OR. (zpphi(jo) < zphitmin(ji, jj))) CYCLE
            END IF
            IF (linquad(zplam(jo), zpphi(jo), zlamtm(:, ji, jj), zphitm(:, ji, jj))) THEN
              IF (llinvalidcell(ji, jj)) THEN
                kproc(jo) = kmyproc + 1000000
                kobsi(jo) = ji + 1
                kobsj(jo) = jj + 1
                CYCLE
              ELSE
                kproc(jo) = kmyproc
                kobsi(jo) = ji + 1
                kobsj(jo) = jj + 1
                EXIT gridpoints
              END IF
            END IF
          END DO
        END DO gridpoints
      END IF
      IF (kproc(jo) == - 1) THEN
        gridpoints_greenwich:DO jj = 1, jlat - 1
          DO ji = 1, jlon - 1
            IF ((zplam(jo) + 360.0_wp > zlamtmax(ji, jj)) .OR. (zplam(jo) + 360.0_wp < zlamtmin(ji, jj))) CYCLE
            IF (ABS(zpphi(jo)) < 85) THEN
              IF ((zpphi(jo) > zphitmax(ji, jj)) .OR. (zpphi(jo) < zphitmin(ji, jj))) CYCLE
            END IF
            IF (linquad(zplam(jo) + 360.0_wp, zpphi(jo), zlamtm(:, ji, jj), zphitm(:, ji, jj))) THEN
              IF (llinvalidcell(ji, jj)) THEN
                kproc(jo) = kmyproc + 1000000
                kobsi(jo) = ji + 1
                kobsj(jo) = jj + 1
                CYCLE
              ELSE
                kproc(jo) = kmyproc
                kobsi(jo) = ji + 1
                kobsj(jo) = jj + 1
                EXIT gridpoints_greenwich
              END IF
            END IF
          END DO
        END DO gridpoints_greenwich
      END IF
    END DO
    IF (ln_grid_global) THEN
      CALL obs_mpp_max_integer(kproc, kobs)
      CALL obs_mpp_max_integer(kobsi, kobs)
      CALL obs_mpp_max_integer(kobsj, kobs)
    ELSE
      CALL obs_mpp_find_obs_proc(kproc, kobs)
    END IF
    WHERE (kproc(:) >= 1000000)
      kproc(:) = kproc(:) - 1000000
    END WHERE
    DEALLOCATE(zlamg, zphig, zmskg, zphitmax, zphitmin, zlamtmax, zlamtmin, llinvalidcell, zlamtm, zphitm, zplam, zpphi)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_grd_bruteforce
  SUBROUTINE obs_grd_lookup(kobs, plam, pphi, kobsi, kobsj, kproc)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: kobs
    REAL(KIND = wp), DIMENSION(kobs), INTENT(IN) :: plam, pphi
    INTEGER, DIMENSION(kobs), INTENT(OUT) :: kobsi, kobsj, kproc
    REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: zplam
    REAL(KIND = wp) :: zlammax
    REAL(KIND = wp) :: zlam
    INTEGER :: ji
    INTEGER :: jj
    INTEGER :: jk
    INTEGER :: jo
    INTEGER :: isx
    INTEGER :: isy
    INTEGER :: jimin
    INTEGER :: jimax
    INTEGER :: jjmin
    INTEGER :: jjmax
    INTEGER :: jojimin
    INTEGER :: jojimax
    INTEGER :: jojjmin
    INTEGER :: jojjmax
    INTEGER :: ipx1
    INTEGER :: ipy1
    INTEGER :: ip
    INTEGER :: jp
    INTEGER :: ipx
    INTEGER :: ipy
    INTEGER :: ipmx
    INTEGER :: jlon
    INTEGER :: jlat
    INTEGER :: joffset
    INTEGER :: jostride
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: zlamg, zphig, zmskg, zphitmax, zphitmin, zlamtmax, zlamtmin
    LOGICAL, DIMENSION(:, :), ALLOCATABLE :: llinvalidcell
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: zlamtm, zphitm
    LOGICAL :: llfourflag
    INTEGER :: ifourflagcountt
    INTEGER :: ifourflagcountf
    INTEGER, DIMENSION(5) :: ifourflagcountr
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_grd_lookup', 'r0', 0, 0)
    IF (ln_grid_global) THEN
      jlon = jpiglo
      jlat = jpjglo
      joffset = nproc
      jostride = jpnij
    ELSE
      jlon = jpi
      jlat = jpj
      joffset = 0
      jostride = 1
    END IF
    ALLOCATE(zlamg(jlon, jlat), zphig(jlon, jlat), zmskg(jlon, jlat), zphitmax(jlon - 1, jlat - 1), zphitmin(jlon - 1, jlat - 1), &
&zlamtmax(jlon - 1, jlat - 1), zlamtmin(jlon - 1, jlat - 1), llinvalidcell(jlon - 1, jlat - 1), zlamtm(4, jlon - 1, jlat - 1), &
&zphitm(4, jlon - 1, jlat - 1))
    IF (ln_grid_global) THEN
      zlamg(:, :) = - 1.E+10
      zphig(:, :) = - 1.E+10
      zmskg(:, :) = - 1.E+10
      DO jj = 1, nlcj
        DO ji = 1, nlci
          zlamg(mig(ji), mjg(jj)) = glamt(ji, jj)
          zphig(mig(ji), mjg(jj)) = gphit(ji, jj)
          zmskg(mig(ji), mjg(jj)) = tmask(ji, jj, 1)
        END DO
      END DO
      CALL mpp_global_max(zlamg)
      CALL mpp_global_max(zphig)
      CALL mpp_global_max(zmskg)
    ELSE
      DO jj = 1, jlat
        DO ji = 1, jlon
          zlamg(ji, jj) = glamt(ji, jj)
          zphig(ji, jj) = gphit(ji, jj)
          zmskg(ji, jj) = tmask(ji, jj, 1)
        END DO
      END DO
    END IF
    ALLOCATE(zplam(kobs))
    DO jo = 1, kobs
      zplam(jo) = plam(jo)
    END DO
    kproc(:) = - 1
    kobsi(:) = - 1
    kobsj(:) = - 1
    DO jj = 1, jlat - 1
      DO ji = 1, jlon - 1
        zlamtm(1, ji, jj) = zlamg(ji, jj)
        zphitm(1, ji, jj) = zphig(ji, jj)
        zlamtm(2, ji, jj) = zlamg(ji + 1, jj)
        zphitm(2, ji, jj) = zphig(ji + 1, jj)
        zlamtm(3, ji, jj) = zlamg(ji + 1, jj + 1)
        zphitm(3, ji, jj) = zphig(ji + 1, jj + 1)
        zlamtm(4, ji, jj) = zlamg(ji, jj + 1)
        zphitm(4, ji, jj) = zphig(ji, jj + 1)
      END DO
    END DO
    WHERE (zlamtm(:, :, :) < 0.0_wp)
      zlamtm(:, :, :) = zlamtm(:, :, :) + 360.0_wp
    END WHERE
    WHERE (zlamtm(:, :, :) > 360.0_wp)
      zlamtm(:, :, :) = zlamtm(:, :, :) - 360.0_wp
    END WHERE
    DO jj = 1, jlat - 1
      DO ji = 1, jlon - 1
        zlammax = MAXVAL(zlamtm(:, ji, jj))
        WHERE (zlammax - zlamtm(:, ji, jj) > 180) zlamtm(:, ji, jj) = zlamtm(:, ji, jj) + 360._wp
        zphitmax(ji, jj) = MAXVAL(zphitm(:, ji, jj))
        zphitmin(ji, jj) = MINVAL(zphitm(:, ji, jj))
        zlamtmax(ji, jj) = MAXVAL(zlamtm(:, ji, jj))
        zlamtmin(ji, jj) = MINVAL(zlamtm(:, ji, jj))
      END DO
    END DO
    llinvalidcell(:, :) = .FALSE.
    DO jj = 1, jlat - 1
      DO ji = 1, jlon - 1
        llinvalidcell(ji, jj) = zmskg(ji, jj) == 0.0_wp .AND. zmskg(ji + 1, jj) == 0.0_wp .AND. zmskg(ji + 1, jj + 1) == 0.0_wp &
&.AND. zmskg(ji, jj + 1) == 0.0_wp
      END DO
    END DO
    IF (lwp) WRITE(numout, FMT = *) 'obs_grid_lookup do coordinate search using lookup table'
    ifourflagcountt = 0
    ifourflagcountf = 0
    ifourflagcountr(:) = 0
    gpkobs:DO jo = 1 + joffset, kobs, jostride
      ipx1 = INT((zplam(jo) - lonmin) / dlon + 1.0)
      ipy1 = INT((pphi(jo) - latmin) / dlat + 1.0)
      ipx = ipx1 + 1
      ipy = ipy1 + 1
      llfourflag = .FALSE.
      IF ((ipx1 > nlons) .OR. (ipy1 > nlats) .OR. (ipx < 1) .OR. (ipy < 1)) THEN
        CYCLE
      END IF
      IF ((ipx > nlons) .OR. (ipy > nlats) .OR. (ipx1 < 1) .OR. (ipy1 < 1)) THEN
        llfourflag = .TRUE.
        ifourflagcountr(1) = ifourflagcountr(1) + 1
      END IF
      IF (.NOT. llfourflag) THEN
        IF (MAXVAL(ixpos(ipx1 : ipx, ipy1 : ipy)) == - 1) CYCLE
      END IF
      jimin = 0
      jimax = 0
      jjmin = 0
      jjmax = 0
      IF (.NOT. llfourflag) THEN
        jojimin = MINVAL(ixpos(ipx1 : ipx, ipy1 : ipy)) - 1
        jojimax = MAXVAL(ixpos(ipx1 : ipx, ipy1 : ipy)) + 1
        jojjmin = MINVAL(iypos(ipx1 : ipx, ipy1 : ipy)) - 1
        jojjmax = MAXVAL(iypos(ipx1 : ipx, ipy1 : ipy)) + 1
        jimin = jojimin - 1
        jimax = jojimax + 1
        jjmin = jojjmin - 1
        jjmax = jojjmax + 1
        IF (jojimin < 0 .OR. jojjmin < 0) THEN
          llfourflag = .TRUE.
          ifourflagcountr(2) = ifourflagcountr(2) + 1
        END IF
        IF (jojimax - jojimin > maxxdiff) THEN
          llfourflag = .TRUE.
          ifourflagcountr(3) = ifourflagcountr(3) + 1
        END IF
        IF (jojjmax - jojjmin > maxydiff) THEN
          llfourflag = .TRUE.
          ifourflagcountr(4) = ifourflagcountr(4) + 1
        END IF
      END IF
      ipmx = 0
      IF (llfourflag) ipmx = 1
      IF (llfourflag) THEN
        ifourflagcountt = ifourflagcountt + 1
      ELSE
        ifourflagcountf = ifourflagcountf + 1
      END IF
      gridpointsn:DO ip = 0, ipmx
        DO jp = 0, ipmx
          IF (kproc(jo) /= - 1) EXIT gridpointsn
          ipx = ipx1 + ip
          ipy = ipy1 + jp
          IF (llfourflag) THEN
            IF (ipx > nlons) ipx = 1
            IF (ipy > nlats) ipy = 1
            IF (ipx < 1) ipx = nlons
            IF (ipy < 1) ipy = nlats
            isx = ixpos(ipx, ipy)
            isy = iypos(ipx, ipy)
            jimin = isx - maxxdiff - 1
            jimax = isx + maxxdiff + 1
            jjmin = isy - maxydiff - 1
            jjmax = isy + maxydiff + 1
          END IF
          IF (jimin < 1) jimin = 1
          IF (jimax > jlon - 1) jimax = jlon - 1
          IF (jjmin < 1) jjmin = 1
          IF (jjmax > jlat - 1) jjmax = jlat - 1
          IF (zplam(jo) < 0.0_wp) zplam(jo) = zplam(jo) + 360.0_wp
          IF (zplam(jo) > 360.0_wp) zplam(jo) = zplam(jo) - 360.0_wp
          gridloop:DO jj = jjmin, jjmax
            DO ji = jimin, jimax
              IF (ABS(zphig(ji, jj) - pphi(jo)) < 1E-6) THEN
                zlam = zlamg(ji, jj)
                IF (zlam < 0.0_wp) zlam = zlam + 360.0_wp
                IF (zlam > 360.0_wp) zlam = zlam - 360.0_wp
                IF (ABS(zlam - zplam(jo)) < 1E-6) THEN
                  IF (llinvalidcell(ji, jj)) THEN
                    kproc(jo) = nproc + 1000000
                    kobsi(jo) = ji + 1
                    kobsj(jo) = jj + 1
                    CYCLE
                  ELSE
                    kproc(jo) = nproc
                    kobsi(jo) = ji + 1
                    kobsj(jo) = jj + 1
                    EXIT gridloop
                  END IF
                END IF
              END IF
            END DO
          END DO gridloop
          IF (zplam(jo) > 180) zplam(jo) = zplam(jo) - 360.0_wp
          IF (kproc(jo) == - 1) THEN
            gridpoints:DO jj = jjmin, jjmax
              DO ji = jimin, jimax
                IF ((zplam(jo) > zlamtmax(ji, jj)) .OR. (zplam(jo) < zlamtmin(ji, jj))) CYCLE
                IF (ABS(pphi(jo)) < 85) THEN
                  IF ((pphi(jo) > zphitmax(ji, jj)) .OR. (pphi(jo) < zphitmin(ji, jj))) CYCLE
                END IF
                IF (linquad(zplam(jo), pphi(jo), zlamtm(:, ji, jj), zphitm(:, ji, jj))) THEN
                  IF (llinvalidcell(ji, jj)) THEN
                    kproc(jo) = nproc + 1000000
                    kobsi(jo) = ji + 1
                    kobsj(jo) = jj + 1
                    CYCLE
                  ELSE
                    kproc(jo) = nproc
                    kobsi(jo) = ji + 1
                    kobsj(jo) = jj + 1
                    EXIT gridpoints
                  END IF
                END IF
              END DO
            END DO gridpoints
          END IF
          IF (kproc(jo) == - 1) THEN
            gridpoints_greenwich:DO jj = jjmin, jjmax
              DO ji = jimin, jimax
                IF ((zplam(jo) + 360.0_wp > zlamtmax(ji, jj)) .OR. (zplam(jo) + 360.0_wp < zlamtmin(ji, jj))) CYCLE
                IF (ABS(pphi(jo)) < 85) THEN
                  IF ((pphi(jo) > zphitmax(ji, jj)) .OR. (pphi(jo) < zphitmin(ji, jj))) CYCLE
                END IF
                IF (linquad(zplam(jo) + 360.0_wp, pphi(jo), zlamtm(:, ji, jj), zphitm(:, ji, jj))) THEN
                  IF (llinvalidcell(ji, jj)) THEN
                    kproc(jo) = nproc + 1000000
                    kobsi(jo) = ji + 1
                    kobsj(jo) = jj + 1
                    CYCLE
                  ELSE
                    kproc(jo) = nproc
                    kobsi(jo) = ji + 1
                    kobsj(jo) = jj + 1
                    EXIT gridpoints_greenwich
                  END IF
                END IF
              END DO
            END DO gridpoints_greenwich
          END IF
        END DO
      END DO gridpointsn
    END DO gpkobs
    IF (ln_grid_global) THEN
      CALL obs_mpp_max_integer(kproc, kobs)
      CALL obs_mpp_max_integer(kobsi, kobs)
      CALL obs_mpp_max_integer(kobsj, kobs)
    ELSE
      CALL obs_mpp_find_obs_proc(kproc, kobs)
    END IF
    WHERE (kproc(:) >= 1000000)
      kproc(:) = kproc(:) - 1000000
    END WHERE
    DEALLOCATE(zlamg, zphig, zmskg, zphitmax, zphitmin, zlamtmax, zlamtmin, llinvalidcell, zlamtm, zphitm, zplam)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_grd_lookup
  SUBROUTINE obs_grid_setup
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = 15), PARAMETER :: cpname = 'obs_grid_setup'
    CHARACTER(LEN = 40) :: cfname
    INTEGER :: ji
    INTEGER :: jj
    INTEGER :: jk
    INTEGER :: jo
    INTEGER :: idfile, idny, idnx, idxpos, idypos
    INTEGER :: idlat, idlon, fileexist
    INTEGER, DIMENSION(2) :: incdim
    CHARACTER(LEN = 20) :: datestr = " ", timestr = " "
    REAL(KIND = wp) :: tmpx1, tmpx2, tmpy1, tmpy2
    REAL(KIND = wp) :: meanxdiff, meanydiff
    REAL(KIND = wp) :: meanxdiff1, meanydiff1
    REAL(KIND = wp) :: meanxdiff2, meanydiff2
    INTEGER :: numx1, numx2, numy1, numy2, df
    INTEGER :: jimin, jimax, jjmin, jjmax
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: lonsi, latsi
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: ixposi, iyposi, iproci
    INTEGER, PARAMETER :: histsize = 90
    INTEGER, DIMENSION(histsize) :: histx1, histx2, histy1, histy2
    REAL, DIMENSION(histsize) :: fhistx1, fhistx2, fhisty1, fhisty2
    REAL(KIND = wp) :: histtol
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_grid_setup', 'r0', 0, 0)
    IF (ln_grid_search_lookup) THEN
      WRITE(numout, FMT = *) 'Calling obs_grid_setup'
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'Grid search resolution : ', rn_gridsearchres
      gsearch_nlons_def = NINT(360.0_wp / rn_gridsearchres)
      gsearch_nlats_def = NINT(180.0_wp / rn_gridsearchres)
      gsearch_lonmin_def = - 180.0_wp + 0.5_wp * rn_gridsearchres
      gsearch_latmin_def = - 90.0_wp + 0.5_wp * rn_gridsearchres
      gsearch_dlon_def = rn_gridsearchres
      gsearch_dlat_def = rn_gridsearchres
      IF (lwp) THEN
        WRITE(numout, FMT = *) 'Grid search gsearch_nlons_def  = ', gsearch_nlons_def
        WRITE(numout, FMT = *) 'Grid search gsearch_nlats_def  = ', gsearch_nlats_def
        WRITE(numout, FMT = *) 'Grid search gsearch_lonmin_def = ', gsearch_lonmin_def
        WRITE(numout, FMT = *) 'Grid search gsearch_latmin_def = ', gsearch_latmin_def
        WRITE(numout, FMT = *) 'Grid search gsearch_dlon_def   = ', gsearch_dlon_def
        WRITE(numout, FMT = *) 'Grid search gsearch_dlat_def   = ', gsearch_dlat_def
      END IF
      IF (ln_grid_global) THEN
        WRITE(cfname, FMT = "(A,'_',A)") TRIM(cn_gridsearchfile), 'global.nc'
      ELSE
        WRITE(cfname, FMT = "(A,'_',I4.4,'of',I4.4,'by',I4.4,'.nc')") TRIM(cn_gridsearchfile), nproc, jpni, jpnj
      END IF
      fileexist = nf90_open(TRIM(cfname), nf90_nowrite, idfile)
      IF (fileexist == nf90_noerr) THEN
        WRITE(numout, FMT = *) 'Reading: ', cfname
        CALL chkerr(nf90_open(TRIM(cfname), nf90_nowrite, idfile), cpname, 729)
        CALL chkerr(nf90_get_att(idfile, nf90_global, 'maxxdiff', maxxdiff), cpname, 731)
        CALL chkerr(nf90_get_att(idfile, nf90_global, 'maxydiff', maxydiff), cpname, 733)
        CALL chkerr(nf90_get_att(idfile, nf90_global, 'dlon', dlon), cpname, 735)
        CALL chkerr(nf90_get_att(idfile, nf90_global, 'dlat', dlat), cpname, 737)
        CALL chkerr(nf90_get_att(idfile, nf90_global, 'lonmin', lonmin), cpname, 739)
        CALL chkerr(nf90_get_att(idfile, nf90_global, 'latmin', latmin), cpname, 741)
        CALL chkerr(nf90_inq_dimid(idfile, 'nx', idnx), cpname, 744)
        CALL chkerr(nf90_inquire_dimension(idfile, idnx, len = nlons), cpname, 746)
        CALL chkerr(nf90_inq_dimid(idfile, 'ny', idny), cpname, 748)
        CALL chkerr(nf90_inquire_dimension(idfile, idny, len = nlats), cpname, 750)
        ALLOCATE(lons(nlons, nlats), lats(nlons, nlats), ixpos(nlons, nlats), iypos(nlons, nlats), iprocn(nlons, nlats))
        CALL chkerr(nf90_inq_varid(idfile, 'XPOS', idxpos), cpname, 761)
        CALL chkerr(nf90_get_var(idfile, idxpos, ixpos), cpname, 763)
        CALL chkerr(nf90_inq_varid(idfile, 'YPOS', idypos), cpname, 765)
        CALL chkerr(nf90_get_var(idfile, idypos, iypos), cpname, 767)
        CALL chkerr(nf90_close(idfile), cpname, 769)
        DO ji = 1, nlons
          DO jj = 1, nlats
            lons(ji, jj) = lonmin + (ji - 1) * dlon
            lats(ji, jj) = latmin + (jj - 1) * dlat
          END DO
        END DO
      ELSE
        IF (lwp) THEN
          WRITE(numout, FMT = *) 'creating: ', cfname
          WRITE(numout, FMT = *) 'calling obs_grid_search: ', nlons * nlats
        END IF
        nlons = gsearch_nlons_def
        nlats = gsearch_nlats_def
        lonmin = gsearch_lonmin_def
        latmin = gsearch_latmin_def
        dlon = gsearch_dlon_def
        dlat = gsearch_dlat_def
        ALLOCATE(lonsi(nlons, nlats), latsi(nlons, nlats), ixposi(nlons, nlats), iyposi(nlons, nlats), iproci(nlons, nlats))
        DO ji = 1, nlons
          DO jj = 1, nlats
            lonsi(ji, jj) = lonmin + (ji - 1) * dlon
            latsi(ji, jj) = latmin + (jj - 1) * dlat
          END DO
        END DO
        CALL obs_grd_bruteforce(jpi, jpj, jpiglo, jpjglo, 1, nlci, 1, nlcj, nproc, jpnij, glamt, gphit, tmask, nlons * nlats, &
&lonsi, latsi, ixposi, iyposi, iproci)
        jimin = 1
        jimax = nlons
        jjmin = 1
        jjmax = nlats
        minlon_xpos:DO ji = 1, nlons
          IF (COUNT(ixposi(ji, :) >= 0) > 0) THEN
            jimin = ji
            EXIT minlon_xpos
          END IF
        END DO minlon_xpos
        maxlon_xpos:DO ji = nlons, 1, - 1
          IF (COUNT(ixposi(ji, :) >= 0) > 0) THEN
            jimax = ji
            EXIT maxlon_xpos
          END IF
        END DO maxlon_xpos
        minlat_xpos:DO jj = 1, nlats
          IF (COUNT(ixposi(:, jj) >= 0) > 0) THEN
            jjmin = jj
            EXIT minlat_xpos
          END IF
        END DO minlat_xpos
        maxlat_xpos:DO jj = nlats, 1, - 1
          IF (COUNT(ixposi(:, jj) >= 0) > 0) THEN
            jjmax = jj
            EXIT maxlat_xpos
          END IF
        END DO maxlat_xpos
        lonmin = lonsi(jimin, jjmin)
        latmin = latsi(jimin, jjmin)
        nlons = jimax - jimin + 1
        nlats = jjmax - jjmin + 1
        ALLOCATE(lons(nlons, nlats), lats(nlons, nlats), ixpos(nlons, nlats), iypos(nlons, nlats), iprocn(nlons, nlats))
        lons(:, :) = lonsi(jimin : jimax, jjmin : jjmax)
        lats(:, :) = latsi(jimin : jimax, jjmin : jjmax)
        ixpos(:, :) = ixposi(jimin : jimax, jjmin : jjmax)
        iypos(:, :) = iyposi(jimin : jimax, jjmin : jjmax)
        iprocn(:, :) = iproci(jimin : jimax, jjmin : jjmax)
        DEALLOCATE(lonsi, latsi, ixposi, iyposi, iproci)
        maxxdiff = 1
        maxydiff = 1
        tmpx1 = 0
        tmpx2 = 0
        tmpy1 = 0
        tmpy2 = 0
        numx1 = 0
        numx2 = 0
        numy1 = 0
        numy2 = 0
        DO ji = 1, nlons - 1
          DO jj = 1, nlats - 1
            IF (ixpos(ji, jj) > 0 .AND. iypos(ji, jj) > 0) THEN
              IF (ixpos(ji + 1, jj) > 0) THEN
                df = ABS(ixpos(ji + 1, jj) - ixpos(ji, jj))
                tmpx1 = tmpx1 + df
                numx1 = numx1 + 1
                IF (df < histsize) histx1(df + 1) = histx1(df + 1) + 1
              END IF
              IF (ixpos(ji, jj + 1) > 0) THEN
                df = ABS(ixpos(ji, jj + 1) - ixpos(ji, jj))
                tmpx2 = tmpx2 + df
                numx2 = numx2 + 1
                IF (df < histsize) histx2(df + 1) = histx2(df + 1) + 1
              END IF
              IF (iypos(ji + 1, jj) > 0) THEN
                df = ABS(iypos(ji + 1, jj) - iypos(ji, jj))
                tmpy1 = tmpy1 + df
                numy1 = numy1 + 1
                IF (df < histsize) histy1(df + 1) = histy1(df + 1) + 1
              END IF
              IF (iypos(ji, jj + 1) > 0) THEN
                df = ABS(iypos(ji, jj + 1) - iypos(ji, jj))
                tmpy2 = tmpy2 + df
                numy2 = numy2 + 1
                IF (df < histsize) histy2(df + 1) = histy2(df + 1) + 1
              END IF
            END IF
          END DO
        END DO
        IF (lwp) THEN
          WRITE(numout, FMT = *) 'histograms'
          WRITE(numout, FMT = *) '0   1   2   3   4   5   6   7   8   9   10 ...'
          WRITE(numout, FMT = *) 'histx1'
          WRITE(numout, FMT = *) histx1
          WRITE(numout, FMT = *) 'histx2'
          WRITE(numout, FMT = *) histx2
          WRITE(numout, FMT = *) 'histy1'
          WRITE(numout, FMT = *) histy1
          WRITE(numout, FMT = *) 'histy2'
          WRITE(numout, FMT = *) histy2
        END IF
        meanxdiff1 = tmpx1 / numx1
        meanydiff1 = tmpy1 / numy1
        meanxdiff2 = tmpx2 / numx2
        meanydiff2 = tmpy2 / numy2
        meanxdiff = MAXVAL((/meanxdiff1, meanxdiff2/))
        meanydiff = MAXVAL((/meanydiff1, meanydiff2/))
        IF (lwp) THEN
          WRITE(numout, FMT = *) tmpx1, tmpx2, tmpy1, tmpy2
          WRITE(numout, FMT = *) numx1, numx2, numy1, numy2
          WRITE(numout, FMT = *) 'meanxdiff: ', meanxdiff, meanxdiff1, meanxdiff2
          WRITE(numout, FMT = *) 'meanydiff: ', meanydiff, meanydiff1, meanydiff2
        END IF
        tmpx1 = 0
        tmpx2 = 0
        tmpy1 = 0
        tmpy2 = 0
        numx1 = 0
        numx2 = 0
        numy1 = 0
        numy2 = 0
        histx1(:) = 0
        histx2(:) = 0
        histy1(:) = 0
        histy2(:) = 0
        limxdiff = meanxdiff * 4
        limydiff = meanydiff * 4
        DO ji = 1, nlons - 1
          DO jj = 1, nlats - 1
            IF (ixpos(ji, jj) > 0 .AND. iypos(ji, jj) > 0) THEN
              IF (ixpos(ji + 1, jj) > 0) THEN
                df = ABS(ixpos(ji + 1, jj) - ixpos(ji, jj))
                tmpx1 = df
                IF (df < limxdiff) numx1 = numx1 + 1
                IF (df < histsize) histx1(df + 1) = histx1(df + 1) + 1
              END IF
              IF (ixpos(ji, jj + 1) > 0) THEN
                df = ABS(ixpos(ji, jj + 1) - ixpos(ji, jj))
                tmpx2 = df
                IF (df < limxdiff) numx2 = numx2 + 1
                IF (df < histsize) histx2(df + 1) = histx2(df + 1) + 1
              END IF
              IF (iypos(ji + 1, jj) > 0) THEN
                df = ABS(iypos(ji + 1, jj) - iypos(ji, jj))
                tmpy1 = df
                IF (df < limydiff) numy1 = numy1 + 1
                IF (df < histsize) histy1(df + 1) = histy1(df + 1) + 1
              END IF
              IF (iypos(ji, jj + 1) > 0) THEN
                df = ABS(iypos(ji, jj + 1) - iypos(ji, jj))
                tmpy2 = df
                IF (df < limydiff) numy2 = numy2 + 1
                IF (df < histsize) histy2(df + 1) = histy2(df + 1) + 1
              END IF
              IF (maxxdiff < tmpx1 .AND. tmpx1 < limxdiff) maxxdiff = tmpx1
              IF (maxxdiff < tmpx2 .AND. tmpx2 < limxdiff) maxxdiff = tmpx2
              IF (maxydiff < tmpy1 .AND. tmpy1 < limydiff) maxydiff = tmpy1
              IF (maxydiff < tmpy2 .AND. tmpy2 < limydiff) maxydiff = tmpy2
            END IF
          END DO
        END DO
        DO ji = 1, histsize - 1
          histx1(ji + 1) = histx1(ji + 1) + histx1(ji)
          histx2(ji + 1) = histx2(ji + 1) + histx2(ji)
          histy1(ji + 1) = histy1(ji + 1) + histy1(ji)
          histy2(ji + 1) = histy2(ji + 1) + histy2(ji)
        END DO
        fhistx1(:) = histx1(:) * 1.0 / numx1
        fhistx2(:) = histx2(:) * 1.0 / numx2
        fhisty1(:) = histy1(:) * 1.0 / numy1
        fhisty2(:) = histy2(:) * 1.0 / numy2
        IF (lwp) THEN
          WRITE(numout, FMT = *) 'cumulative histograms'
          WRITE(numout, FMT = *) '0   1   2   3   4   5   6   7   8   9   10 ...'
          WRITE(numout, FMT = *) 'fhistx1'
          WRITE(numout, FMT = *) fhistx1
          WRITE(numout, FMT = *) 'fhistx2'
          WRITE(numout, FMT = *) fhistx2
          WRITE(numout, FMT = *) 'fhisty1'
          WRITE(numout, FMT = *) fhisty1
          WRITE(numout, FMT = *) 'fhisty2'
          WRITE(numout, FMT = *) fhisty2
        END IF
        histtol = 0.999
        tmpx1 = MAXVAL(MAXLOC(fhistx1(:), mask = (fhistx1(:) <= histtol)))
        tmpx2 = MAXVAL(MAXLOC(fhistx2(:), mask = (fhistx2(:) <= histtol)))
        tmpy1 = MAXVAL(MAXLOC(fhisty1(:), mask = (fhisty1(:) <= histtol)))
        tmpy2 = MAXVAL(MAXLOC(fhisty2(:), mask = (fhisty2(:) <= histtol)))
        maxxdiff = MAXVAL((/tmpx1, tmpx2/)) + 1
        maxydiff = MAXVAL((/tmpy1, tmpy2/)) + 1
        IF ((.NOT. ln_grid_global) .OR. ((ln_grid_global) .AND. (nproc == 0))) THEN
          CALL chkerr(nf90_create(TRIM(cfname), nf90_clobber, idfile), cpname, 1072)
          CALL chkerr(nf90_put_att(idfile, nf90_global, 'title', 'Mapping file from lon/lat to model grid point'), cpname, 1075)
          CALL chkerr(nf90_put_att(idfile, nf90_global, 'maxxdiff', maxxdiff), cpname, 1078)
          CALL chkerr(nf90_put_att(idfile, nf90_global, 'maxydiff', maxydiff), cpname, 1081)
          CALL chkerr(nf90_put_att(idfile, nf90_global, 'dlon', dlon), cpname, 1083)
          CALL chkerr(nf90_put_att(idfile, nf90_global, 'dlat', dlat), cpname, 1085)
          CALL chkerr(nf90_put_att(idfile, nf90_global, 'lonmin', lonmin), cpname, 1088)
          CALL chkerr(nf90_put_att(idfile, nf90_global, 'latmin', latmin), cpname, 1091)
          CALL chkerr(nf90_def_dim(idfile, 'nx', nlons, idnx), cpname, 1094)
          CALL chkerr(nf90_def_dim(idfile, 'ny', nlats, idny), cpname, 1096)
          incdim(1) = idnx
          incdim(2) = idny
          CALL chkerr(nf90_def_var(idfile, 'LON', nf90_float, incdim, idlon), cpname, 1103)
          CALL chkerr(nf90_put_att(idfile, idlon, 'long_name', 'longitude'), cpname, 1106)
          CALL chkerr(nf90_def_var(idfile, 'LAT', nf90_float, incdim, idlat), cpname, 1110)
          CALL chkerr(nf90_put_att(idfile, idlat, 'long_name', 'latitude'), cpname, 1113)
          CALL chkerr(nf90_def_var(idfile, 'XPOS', nf90_int, incdim, idxpos), cpname, 1117)
          CALL chkerr(nf90_put_att(idfile, idxpos, 'long_name', 'x position'), cpname, 1120)
          CALL chkerr(nf90_put_att(idfile, idxpos, '_FillValue', - 1), cpname, 1122)
          CALL chkerr(nf90_def_var(idfile, 'YPOS', nf90_int, incdim, idypos), cpname, 1126)
          CALL chkerr(nf90_put_att(idfile, idypos, 'long_name', 'y position'), cpname, 1129)
          CALL chkerr(nf90_put_att(idfile, idypos, '_FillValue', - 1), cpname, 1131)
          CALL chkerr(nf90_enddef(idfile), cpname, 1133)
          CALL chkerr(nf90_put_var(idfile, idlon, lons), cpname, 1136)
          CALL chkerr(nf90_put_var(idfile, idlat, lats), cpname, 1138)
          CALL chkerr(nf90_put_var(idfile, idxpos, ixpos), cpname, 1140)
          CALL chkerr(nf90_put_var(idfile, idypos, iypos), cpname, 1142)
          CALL chkerr(nf90_close(idfile), cpname, 1144)
        END IF
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_grid_setup
  SUBROUTINE obs_grid_deallocate
    IF (ln_grid_search_lookup) THEN
      DEALLOCATE(lons, lats, ixpos, iypos, iprocn)
    END IF
  END SUBROUTINE obs_grid_deallocate
  SUBROUTINE obs_level_search(kgrd, pgrddep, kobs, pobsdep, kobsk)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kgrd
    REAL(KIND = wp), DIMENSION(kgrd), INTENT(INOUT) :: pgrddep
    INTEGER, INTENT(IN) :: kobs
    REAL(KIND = wp), DIMENSION(kobs), INTENT(INOUT) :: pobsdep
    INTEGER, DIMENSION(kobs), INTENT(OUT) :: kobsk
    INTEGER :: ji
    INTEGER :: jk
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_level_search', 'r0', 0, 0)
    DO ji = 1, kobs
      kobsk(ji) = 1
      depk:DO jk = 2, kgrd
        IF (pgrddep(jk) >= pobsdep(ji)) EXIT depk
      END DO depk
      kobsk(ji) = jk
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_level_search
  LOGICAL FUNCTION linquad(px, py, pxv, pyv)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN) :: px
    REAL(KIND = wp), INTENT(IN) :: py
    REAL(KIND = wp), DIMENSION(4), INTENT(IN) :: pxv, pyv
    REAL(KIND = wp) :: zst1
    REAL(KIND = wp) :: zst2
    REAL(KIND = wp) :: zst3
    REAL(KIND = wp) :: zst4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('linquad', 'r0', 0, 0)
    linquad = .FALSE.
    zst1 = (px - pxv(1)) * (py - pyv(4)) - (py - pyv(1)) * (px - pxv(4))
    IF (zst1 <= 0.0) THEN
      zst2 = (px - pxv(4)) * (py - pyv(3)) - (py - pyv(4)) * (px - pxv(3))
      IF (zst2 <= 0.0) THEN
        zst3 = (px - pxv(3)) * (py - pyv(2)) - (py - pyv(3)) * (px - pxv(2))
        IF (zst3 <= 0.0) THEN
          zst4 = (px - pxv(2)) * (py - pyv(1)) - (py - pyv(2)) * (px - pxv(1))
          IF (zst4 <= 0.0) linquad = .TRUE.
        END IF
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END FUNCTION linquad
  REAL FUNCTION maxdist(pxv, pyv)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(4), INTENT(IN) :: pxv, pyv
    REAL(KIND = wp), DIMENSION(4) :: zxv, zyv, za, zb, zc
    REAL(KIND = wp) :: zdist
    INTEGER :: ji
    INTEGER :: jj
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('maxdist', 'r0', 0, 0)
    DO ji = 1, 4
      zxv(ji) = pxv(ji) * rad
      zyv(ji) = pyv(ji) * rad
    END DO
    DO ji = 1, 4
      za(ji) = SIN(zyv(ji))
      zb(ji) = COS(zyv(ji)) * COS(zxv(ji))
      zc(ji) = COS(zyv(ji)) * SIN(zxv(ji))
    END DO
    maxdist = 0.0
    DO jj = 1, 4
      DO ji = jj + 1, 4
        zdist = grt_cir_dis(za(jj), za(ji), zb(jj), zb(ji), zc(jj), zc(ji))
        IF (zdist > maxdist) THEN
          maxdist = zdist
        END IF
      END DO
    END DO
    maxdist = maxdist / rad
    CALL profile_psy_data0 % PostEnd
  END FUNCTION maxdist
  SUBROUTINE find_obs_proc(kldi, klei, kldj, klej, kmyproc, kobsp, kobsi, kobsj, kno)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kldi
    INTEGER, INTENT(IN) :: klei
    INTEGER, INTENT(IN) :: kldj
    INTEGER, INTENT(IN) :: klej
    INTEGER, INTENT(IN) :: kmyproc
    INTEGER, INTENT(IN) :: kno
    INTEGER, DIMENSION(kno), INTENT(IN) :: kobsi
    INTEGER, DIMENSION(kno), INTENT(IN) :: kobsj
    INTEGER, DIMENSION(kno), INTENT(INOUT) :: kobsp
    INTEGER :: ji
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('find_obs_proc', 'r0', 0, 0)
    DO ji = 1, kno
      IF (kobsi(ji) < kldi .OR. kobsj(ji) < kldj .OR. kobsi(ji) > klei .OR. kobsj(ji) > klej) THEN
        IF (lwp .AND. kobsp(ji) /= - 1) WRITE(numout, FMT = *) 'kobs: ', kobsi(ji), kobsj(ji), kobsp(ji)
        kobsp(ji) = 1000000
      END IF
    END DO
    WHERE (kobsp(:) /= kmyproc) kobsp(:) = 1000000
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE find_obs_proc
END MODULE obs_grid