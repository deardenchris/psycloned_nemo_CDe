MODULE obs_write
  USE par_kind, ONLY: wp
  USE in_out_manager
  USE dom_oce
  USE obs_types
  USE julian, ONLY: greg2jul
  USE obs_utils, ONLY: chkerr
  USE obs_profiles_def
  USE obs_surf_def
  USE obs_fbm
  USE obs_grid
  USE obs_conv
  USE obs_const
  USE obs_mpp
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: obs_wri_prof, obs_wri_surf, obswriinfo
  TYPE :: obswriinfo
    INTEGER :: inum
    INTEGER, POINTER, DIMENSION(:) :: ipoint
    CHARACTER(LEN = ilenname), POINTER, DIMENSION(:) :: cdname
    CHARACTER(LEN = ilenlong), POINTER, DIMENSION(:, :) :: cdlong
    CHARACTER(LEN = ilenunit), POINTER, DIMENSION(:, :) :: cdunit
  END TYPE obswriinfo
  CONTAINS
  SUBROUTINE obs_wri_prof(profdata, padd, pext)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(obs_prof), INTENT(INOUT) :: profdata
    TYPE(obswriinfo), OPTIONAL :: padd
    TYPE(obswriinfo), OPTIONAL :: pext
    TYPE(obfbdata) :: fbdata
    CHARACTER(LEN = 40) :: clfname
    CHARACTER(LEN = 10) :: clfiletype
    INTEGER :: ilevel
    INTEGER :: jvar
    INTEGER :: jo
    INTEGER :: jk
    INTEGER :: ik
    INTEGER :: ja
    INTEGER :: je
    INTEGER :: iadd
    INTEGER :: iext
    REAL(KIND = wp) :: zpres
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_wri_prof', 'r0', 0, 0)
    IF (PRESENT(padd)) THEN
      iadd = padd % inum
    ELSE
      iadd = 0
    END IF
    IF (PRESENT(pext)) THEN
      iext = pext % inum
    ELSE
      iext = 0
    END IF
    CALL init_obfbdata(fbdata)
    ilevel = 0
    DO jvar = 1, 2
      ilevel = MAX(ilevel, MAXVAL(profdata % var(jvar) % nvlidx(:)))
    END DO
    SELECT CASE (TRIM(profdata % cvars(1)))
    CASE ('POTM')
      clfiletype = 'profb'
      CALL alloc_obfbdata(fbdata, 2, profdata % nprof, ilevel, 1 + iadd, 1 + iext, .TRUE.)
      fbdata % cname(1) = profdata % cvars(1)
      fbdata % cname(2) = profdata % cvars(2)
      fbdata % coblong(1) = 'Potential temperature'
      fbdata % coblong(2) = 'Practical salinity'
      fbdata % cobunit(1) = 'Degrees centigrade'
      fbdata % cobunit(2) = 'PSU'
      fbdata % cextname(1) = 'TEMP'
      fbdata % cextlong(1) = 'Insitu temperature'
      fbdata % cextunit(1) = 'Degrees centigrade'
      fbdata % caddlong(1, 1) = 'Model interpolated potential temperature'
      fbdata % caddlong(1, 2) = 'Model interpolated practical salinity'
      fbdata % caddunit(1, 1) = 'Degrees centigrade'
      fbdata % caddunit(1, 2) = 'PSU'
      fbdata % cgrid(:) = 'T'
      DO je = 1, iext
        fbdata % cextname(1 + je) = pext % cdname(je)
        fbdata % cextlong(1 + je) = pext % cdlong(je, 1)
        fbdata % cextunit(1 + je) = pext % cdunit(je, 1)
      END DO
      DO ja = 1, iadd
        fbdata % caddname(1 + ja) = padd % cdname(ja)
        DO jvar = 1, 2
          fbdata % caddlong(1 + ja, jvar) = padd % cdlong(ja, jvar)
          fbdata % caddunit(1 + ja, jvar) = padd % cdunit(ja, jvar)
        END DO
      END DO
    CASE ('UVEL')
      clfiletype = 'velfb'
      CALL alloc_obfbdata(fbdata, 2, profdata % nprof, ilevel, 1, 0, .TRUE.)
      fbdata % cname(1) = profdata % cvars(1)
      fbdata % cname(2) = profdata % cvars(2)
      fbdata % coblong(1) = 'Zonal velocity'
      fbdata % coblong(2) = 'Meridional velocity'
      fbdata % cobunit(1) = 'm/s'
      fbdata % cobunit(2) = 'm/s'
      DO je = 1, iext
        fbdata % cextname(je) = pext % cdname(je)
        fbdata % cextlong(je) = pext % cdlong(je, 1)
        fbdata % cextunit(je) = pext % cdunit(je, 1)
      END DO
      fbdata % caddlong(1, 1) = 'Model interpolated zonal velocity'
      fbdata % caddlong(1, 2) = 'Model interpolated meridional velocity'
      fbdata % caddunit(1, 1) = 'm/s'
      fbdata % caddunit(1, 2) = 'm/s'
      fbdata % cgrid(1) = 'U'
      fbdata % cgrid(2) = 'V'
      DO ja = 1, iadd
        fbdata % caddname(1 + ja) = padd % cdname(ja)
        fbdata % caddlong(1 + ja, 1) = padd % cdlong(ja, 1)
        fbdata % caddunit(1 + ja, 1) = padd % cdunit(ja, 1)
      END DO
    END SELECT
    fbdata % caddname(1) = 'Hx'
    WRITE(clfname, FMT = "(A,'_fdbk_',I4.4,'.nc')") TRIM(clfiletype), nproc
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'obs_wri_prof :'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~~'
      WRITE(numout, FMT = *) 'Writing ' // TRIM(clfiletype) // ' feedback file : ', TRIM(clfname)
    END IF
    fbdata % cdjuldref = '19500101000000'
    DO jo = 1, profdata % nprof
      fbdata % plam(jo) = profdata % rlam(jo)
      fbdata % pphi(jo) = profdata % rphi(jo)
      WRITE(fbdata % cdtyp(jo), FMT = '(I4)') profdata % ntyp(jo)
      fbdata % ivqc(jo, :) = profdata % ivqc(jo, :)
      fbdata % ivqcf(:, jo, :) = profdata % ivqcf(:, jo, :)
      IF (profdata % nqc(jo) > 255) THEN
        fbdata % ioqc(jo) = IBSET(profdata % nqc(jo), 2)
        fbdata % ioqcf(1, jo) = profdata % nqcf(1, jo)
        fbdata % ioqcf(2, jo) = profdata % nqc(jo)
      ELSE
        fbdata % ioqc(jo) = profdata % nqc(jo)
        fbdata % ioqcf(:, jo) = profdata % nqcf(:, jo)
      END IF
      fbdata % ipqc(jo) = profdata % ipqc(jo)
      fbdata % ipqcf(:, jo) = profdata % ipqcf(:, jo)
      fbdata % itqc(jo) = profdata % itqc(jo)
      fbdata % itqcf(:, jo) = profdata % itqcf(:, jo)
      fbdata % cdwmo(jo) = profdata % cwmo(jo)
      fbdata % kindex(jo) = profdata % npfil(jo)
      DO jvar = 1, profdata % nvar
        IF (ln_grid_global) THEN
          fbdata % iobsi(jo, jvar) = profdata % mi(jo, jvar)
          fbdata % iobsj(jo, jvar) = profdata % mj(jo, jvar)
        ELSE
          fbdata % iobsi(jo, jvar) = mig(profdata % mi(jo, jvar))
          fbdata % iobsj(jo, jvar) = mjg(profdata % mj(jo, jvar))
        END IF
      END DO
      CALL greg2jul(0, profdata % nmin(jo), profdata % nhou(jo), profdata % nday(jo), profdata % nmon(jo), profdata % nyea(jo), &
&fbdata % ptim(jo), krefdate = 19500101)
      DO jvar = 1, 2
        DO jk = profdata % npvsta(jo, jvar), profdata % npvend(jo, jvar)
          ik = profdata % var(jvar) % nvlidx(jk)
          fbdata % padd(ik, jo, 1, jvar) = profdata % var(jvar) % vmod(jk)
          fbdata % pob(ik, jo, jvar) = profdata % var(jvar) % vobs(jk)
          fbdata % pdep(ik, jo) = profdata % var(jvar) % vdep(jk)
          fbdata % idqc(ik, jo) = profdata % var(jvar) % idqc(jk)
          fbdata % idqcf(:, ik, jo) = profdata % var(jvar) % idqcf(:, jk)
          IF (profdata % var(jvar) % nvqc(jk) > 255) THEN
            fbdata % ivlqc(ik, jo, jvar) = IBSET(profdata % var(jvar) % nvqc(jk), 2)
            fbdata % ivlqcf(1, ik, jo, jvar) = profdata % var(jvar) % nvqcf(1, jk)
            fbdata % ivlqcf(2, ik, jo, jvar) = IAND(profdata % var(jvar) % nvqc(jk), B'0000000011111111')
          ELSE
            fbdata % ivlqc(ik, jo, jvar) = profdata % var(jvar) % nvqc(jk)
            fbdata % ivlqcf(:, ik, jo, jvar) = profdata % var(jvar) % nvqcf(:, jk)
          END IF
          fbdata % iobsk(ik, jo, jvar) = profdata % var(jvar) % mvk(jk)
          DO ja = 1, iadd
            fbdata % padd(ik, jo, 1 + ja, jvar) = profdata % var(jvar) % vext(jk, padd % ipoint(ja))
          END DO
          DO je = 1, iext
            fbdata % pext(ik, jo, 1 + je) = profdata % var(jvar) % vext(jk, pext % ipoint(je))
          END DO
          IF ((jvar == 1) .AND. (TRIM(profdata % cvars(1)) == 'POTM')) THEN
            fbdata % pext(ik, jo, 1) = profdata % var(jvar) % vext(jk, 1)
          END IF
        END DO
      END DO
    END DO
    IF (TRIM(profdata % cvars(1)) == 'POTM') THEN
      DO jo = 1, fbdata % nobs
        IF (fbdata % pphi(jo) < 9999.0) THEN
          DO jk = 1, fbdata % nlev
            IF ((fbdata % pob(jk, jo, 1) >= 9999.0) .AND. (fbdata % pdep(jk, jo) < 9999.0) .AND. (fbdata % padd(jk, jo, 1, 2) < &
&9999.0) .AND. (fbdata % pext(jk, jo, 1) < 9999.0)) THEN
              zpres = dep_to_p(REAL(fbdata % pdep(jk, jo), wp), REAL(fbdata % pphi(jo), wp))
              fbdata % pob(jk, jo, 1) = potemp(REAL(fbdata % padd(jk, jo, 1, 2), wp), REAL(fbdata % pext(jk, jo, 1), wp), zpres, &
&0.0_wp)
            END IF
          END DO
        END IF
      END DO
    END IF
    CALL write_obfbdata(clfname, fbdata)
    CALL obs_wri_stats(fbdata)
    CALL dealloc_obfbdata(fbdata)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_wri_prof
  SUBROUTINE obs_wri_surf(surfdata, padd, pext)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    IMPLICIT NONE
    TYPE(obs_surf), INTENT(INOUT) :: surfdata
    TYPE(obswriinfo), OPTIONAL :: padd
    TYPE(obswriinfo), OPTIONAL :: pext
    TYPE(obfbdata) :: fbdata
    CHARACTER(LEN = 40) :: clfname
    CHARACTER(LEN = 10) :: clfiletype
    CHARACTER(LEN = 12), PARAMETER :: cpname = 'obs_wri_surf'
    INTEGER :: jo
    INTEGER :: ja
    INTEGER :: je
    INTEGER :: iadd
    INTEGER :: iext
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_wri_surf', 'r0', 0, 0)
    IF (PRESENT(padd)) THEN
      iadd = padd % inum
    ELSE
      iadd = 0
    END IF
    IF (PRESENT(pext)) THEN
      iext = pext % inum
    ELSE
      iext = 0
    END IF
    CALL init_obfbdata(fbdata)
    SELECT CASE (TRIM(surfdata % cvars(1)))
    CASE ('SLA')
      CALL alloc_obfbdata(fbdata, 1, surfdata % nsurf, 1, 2 + iadd, 1 + iext, .TRUE.)
      clfiletype = 'slafb'
      fbdata % cname(1) = surfdata % cvars(1)
      fbdata % coblong(1) = 'Sea level anomaly'
      fbdata % cobunit(1) = 'Metres'
      fbdata % cextname(1) = 'MDT'
      fbdata % cextlong(1) = 'Mean dynamic topography'
      fbdata % cextunit(1) = 'Metres'
      DO je = 1, iext
        fbdata % cextname(je) = pext % cdname(je)
        fbdata % cextlong(je) = pext % cdlong(je, 1)
        fbdata % cextunit(je) = pext % cdunit(je, 1)
      END DO
      fbdata % caddlong(1, 1) = 'Model interpolated SSH - MDT'
      fbdata % caddunit(1, 1) = 'Metres'
      fbdata % caddname(2) = 'SSH'
      fbdata % caddlong(2, 1) = 'Model Sea surface height'
      fbdata % caddunit(2, 1) = 'Metres'
      fbdata % cgrid(1) = 'T'
      DO ja = 1, iadd
        fbdata % caddname(2 + ja) = padd % cdname(ja)
        fbdata % caddlong(2 + ja, 1) = padd % cdlong(ja, 1)
        fbdata % caddunit(2 + ja, 1) = padd % cdunit(ja, 1)
      END DO
    CASE ('SST')
      CALL alloc_obfbdata(fbdata, 1, surfdata % nsurf, 1, 1 + iadd, iext, .TRUE.)
      clfiletype = 'sstfb'
      fbdata % cname(1) = surfdata % cvars(1)
      fbdata % coblong(1) = 'Sea surface temperature'
      fbdata % cobunit(1) = 'Degree centigrade'
      DO je = 1, iext
        fbdata % cextname(je) = pext % cdname(je)
        fbdata % cextlong(je) = pext % cdlong(je, 1)
        fbdata % cextunit(je) = pext % cdunit(je, 1)
      END DO
      fbdata % caddlong(1, 1) = 'Model interpolated SST'
      fbdata % caddunit(1, 1) = 'Degree centigrade'
      fbdata % cgrid(1) = 'T'
      DO ja = 1, iadd
        fbdata % caddname(1 + ja) = padd % cdname(ja)
        fbdata % caddlong(1 + ja, 1) = padd % cdlong(ja, 1)
        fbdata % caddunit(1 + ja, 1) = padd % cdunit(ja, 1)
      END DO
    CASE ('ICECONC')
      CALL alloc_obfbdata(fbdata, 1, surfdata % nsurf, 1, 1 + iadd, iext, .TRUE.)
      clfiletype = 'sicfb'
      fbdata % cname(1) = surfdata % cvars(1)
      fbdata % coblong(1) = 'Sea ice'
      fbdata % cobunit(1) = 'Fraction'
      DO je = 1, iext
        fbdata % cextname(je) = pext % cdname(je)
        fbdata % cextlong(je) = pext % cdlong(je, 1)
        fbdata % cextunit(je) = pext % cdunit(je, 1)
      END DO
      fbdata % caddlong(1, 1) = 'Model interpolated ICE'
      fbdata % caddunit(1, 1) = 'Fraction'
      fbdata % cgrid(1) = 'T'
      DO ja = 1, iadd
        fbdata % caddname(1 + ja) = padd % cdname(ja)
        fbdata % caddlong(1 + ja, 1) = padd % cdlong(ja, 1)
        fbdata % caddunit(1 + ja, 1) = padd % cdunit(ja, 1)
      END DO
    CASE ('SSS')
      CALL alloc_obfbdata(fbdata, 1, surfdata % nsurf, 1, 1 + iadd, iext, .TRUE.)
      clfiletype = 'sssfb'
      fbdata % cname(1) = surfdata % cvars(1)
      fbdata % coblong(1) = 'Sea surface salinity'
      fbdata % cobunit(1) = 'psu'
      DO je = 1, iext
        fbdata % cextname(je) = pext % cdname(je)
        fbdata % cextlong(je) = pext % cdlong(je, 1)
        fbdata % cextunit(je) = pext % cdunit(je, 1)
      END DO
      fbdata % caddlong(1, 1) = 'Model interpolated SSS'
      fbdata % caddunit(1, 1) = 'psu'
      fbdata % cgrid(1) = 'T'
      DO ja = 1, iadd
        fbdata % caddname(1 + ja) = padd % cdname(ja)
        fbdata % caddlong(1 + ja, 1) = padd % cdlong(ja, 1)
        fbdata % caddunit(1 + ja, 1) = padd % cdunit(ja, 1)
      END DO
    CASE DEFAULT
      CALL ctl_stop('Unknown observation type ' // TRIM(surfdata % cvars(1)) // ' in obs_wri_surf')
    END SELECT
    fbdata % caddname(1) = 'Hx'
    WRITE(clfname, FMT = "(A,'_fdbk_',I4.4,'.nc')") TRIM(clfiletype), nproc
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'obs_wri_surf :'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~~'
      WRITE(numout, FMT = *) 'Writing ' // TRIM(surfdata % cvars(1)) // ' feedback file : ', TRIM(clfname)
    END IF
    fbdata % cdjuldref = '19500101000000'
    DO jo = 1, surfdata % nsurf
      fbdata % plam(jo) = surfdata % rlam(jo)
      fbdata % pphi(jo) = surfdata % rphi(jo)
      WRITE(fbdata % cdtyp(jo), FMT = '(I4)') surfdata % ntyp(jo)
      fbdata % ivqc(jo, :) = 0
      fbdata % ivqcf(:, jo, :) = 0
      IF (surfdata % nqc(jo) > 255) THEN
        fbdata % ioqc(jo) = 4
        fbdata % ioqcf(1, jo) = 0
        fbdata % ioqcf(2, jo) = IAND(surfdata % nqc(jo), B'0000000011111111')
      ELSE
        fbdata % ioqc(jo) = surfdata % nqc(jo)
        fbdata % ioqcf(:, jo) = 0
      END IF
      fbdata % ipqc(jo) = 0
      fbdata % ipqcf(:, jo) = 0
      fbdata % itqc(jo) = 0
      fbdata % itqcf(:, jo) = 0
      fbdata % cdwmo(jo) = surfdata % cwmo(jo)
      fbdata % kindex(jo) = surfdata % nsfil(jo)
      IF (ln_grid_global) THEN
        fbdata % iobsi(jo, 1) = surfdata % mi(jo)
        fbdata % iobsj(jo, 1) = surfdata % mj(jo)
      ELSE
        fbdata % iobsi(jo, 1) = mig(surfdata % mi(jo))
        fbdata % iobsj(jo, 1) = mjg(surfdata % mj(jo))
      END IF
      CALL greg2jul(0, surfdata % nmin(jo), surfdata % nhou(jo), surfdata % nday(jo), surfdata % nmon(jo), surfdata % nyea(jo), &
&fbdata % ptim(jo), krefdate = 19500101)
      fbdata % padd(1, jo, 1, 1) = surfdata % rmod(jo, 1)
      IF (TRIM(surfdata % cvars(1)) == 'SLA') fbdata % padd(1, jo, 2, 1) = surfdata % rext(jo, 1)
      fbdata % pob(1, jo, 1) = surfdata % robs(jo, 1)
      fbdata % pdep(1, jo) = 0.0
      fbdata % idqc(1, jo) = 0
      fbdata % idqcf(:, 1, jo) = 0
      IF (surfdata % nqc(jo) > 255) THEN
        fbdata % ivqc(jo, 1) = 4
        fbdata % ivlqc(1, jo, 1) = 4
        fbdata % ivlqcf(1, 1, jo, 1) = 0
        fbdata % ivlqcf(2, 1, jo, 1) = IAND(surfdata % nqc(jo), B'0000000011111111')
      ELSE
        fbdata % ivqc(jo, 1) = surfdata % nqc(jo)
        fbdata % ivlqc(1, jo, 1) = surfdata % nqc(jo)
        fbdata % ivlqcf(:, 1, jo, 1) = 0
      END IF
      fbdata % iobsk(1, jo, 1) = 0
      IF (TRIM(surfdata % cvars(1)) == 'SLA') fbdata % pext(1, jo, 1) = surfdata % rext(jo, 2)
      DO ja = 1, iadd
        fbdata % padd(1, jo, 2 + ja, 1) = surfdata % rext(jo, padd % ipoint(ja))
      END DO
      DO je = 1, iext
        fbdata % pext(1, jo, 1 + je) = surfdata % rext(jo, pext % ipoint(je))
      END DO
    END DO
    CALL write_obfbdata(clfname, fbdata)
    CALL obs_wri_stats(fbdata)
    CALL dealloc_obfbdata(fbdata)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_wri_surf
  SUBROUTINE obs_wri_stats(fbdata)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(obfbdata) :: fbdata
    INTEGER :: jvar
    INTEGER :: jo
    INTEGER :: jk
    INTEGER :: inumgoodobs
    INTEGER :: inumgoodobsmpp
    REAL(KIND = wp) :: zsumx
    REAL(KIND = wp) :: zsumx2
    REAL(KIND = wp) :: zomb
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_wri_stats', 'r0', 0, 0)
    IF (lwp) THEN
      WRITE(numout, FMT = *) ''
      WRITE(numout, FMT = *) 'obs_wri_stats :'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~~~~'
    END IF
    DO jvar = 1, fbdata % nvar
      zsumx = 0.0_wp
      zsumx2 = 0.0_wp
      inumgoodobs = 0
      DO jo = 1, fbdata % nobs
        DO jk = 1, fbdata % nlev
          IF ((fbdata % pob(jk, jo, jvar) < 9999.0) .AND. (fbdata % pdep(jk, jo) < 9999.0) .AND. (fbdata % padd(jk, jo, 1, jvar) < &
&9999.0)) THEN
            zomb = fbdata % pob(jk, jo, jvar) - fbdata % padd(jk, jo, 1, jvar)
            zsumx = zsumx + zomb
            zsumx2 = zsumx2 + zomb ** 2
            inumgoodobs = inumgoodobs + 1
          END IF
        END DO
      END DO
      CALL obs_mpp_sum_integer(inumgoodobs, inumgoodobsmpp)
      CALL mpp_sum('obs_write', zsumx)
      CALL mpp_sum('obs_write', zsumx2)
      IF (lwp) THEN
        WRITE(numout, FMT = *) 'Type: ', fbdata % cname(jvar), '  Total number of good observations: ', inumgoodobsmpp
        WRITE(numout, FMT = *) 'Overall mean obs minus model of the good observations: ', zsumx / inumgoodobsmpp
        WRITE(numout, FMT = *) 'Overall RMS obs minus model of the good observations: ', SQRT(zsumx2 / inumgoodobsmpp)
        WRITE(numout, FMT = *) ''
      END IF
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_wri_stats
END MODULE obs_write