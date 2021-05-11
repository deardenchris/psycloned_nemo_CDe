MODULE fldread
  USE oce
  USE dom_oce
  USE phycst
  USE sbc_oce
  USE geo2ocean
  USE in_out_manager
  USE iom
  USE ioipsl, ONLY: ymds2ju, ju2ymds
  USE lib_mpp
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: fld_map
  PUBLIC :: fld_read, fld_fill
  PUBLIC :: fld_clopn
  TYPE, PUBLIC :: FLD_N
    CHARACTER(LEN = 256) :: clname
    REAL(KIND = wp) :: nfreqh
    CHARACTER(LEN = 34) :: clvar
    LOGICAL :: ln_tint
    LOGICAL :: ln_clim
    CHARACTER(LEN = 8) :: cltype
    CHARACTER(LEN = 256) :: wname
    CHARACTER(LEN = 34) :: vcomp
    CHARACTER(LEN = 34) :: lname
  END TYPE FLD_N
  TYPE, PUBLIC :: FLD
    CHARACTER(LEN = 256) :: clrootname
    CHARACTER(LEN = 256) :: clname
    REAL(KIND = wp) :: nfreqh
    CHARACTER(LEN = 34) :: clvar
    LOGICAL :: ln_tint
    LOGICAL :: ln_clim
    CHARACTER(LEN = 8) :: cltype
    INTEGER :: num
    INTEGER, DIMENSION(2) :: nrec_b
    INTEGER, DIMENSION(2) :: nrec_a
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: fnow
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :, :) :: fdta
    CHARACTER(LEN = 256) :: wgtname
    CHARACTER(LEN = 34) :: vcomp
    LOGICAL, DIMENSION(2) :: rotn
    INTEGER :: nreclast
    CHARACTER(LEN = 256) :: lsmname
    INTEGER :: igrd
    INTEGER :: ibdy
  END TYPE FLD
  TYPE, PUBLIC :: MAP_POINTER
    INTEGER, POINTER, DIMENSION(:) :: ptr
    LOGICAL :: ll_unstruc
  END TYPE MAP_POINTER
  TYPE :: WGT
    CHARACTER(LEN = 256) :: wgtname
    INTEGER, DIMENSION(2) :: ddims
    INTEGER, DIMENSION(2) :: botleft
    INTEGER, DIMENSION(2) :: topright
    INTEGER :: jpiwgt
    INTEGER :: jpjwgt
    INTEGER :: numwgt
    INTEGER :: nestid
    INTEGER :: overlap
    LOGICAL :: cyclic
    INTEGER, DIMENSION(:, :, :), POINTER :: data_jpi
    INTEGER, DIMENSION(:, :, :), POINTER :: data_jpj
    REAL(KIND = wp), DIMENSION(:, :, :), POINTER :: data_wgt
    REAL(KIND = wp), DIMENSION(:, :, :), POINTER :: fly_dta
    REAL(KIND = wp), DIMENSION(:, :, :), POINTER :: col
  END TYPE WGT
  INTEGER, PARAMETER :: tot_wgts = 20
  TYPE(WGT), DIMENSION(tot_wgts) :: ref_wgts
  INTEGER :: nxt_wgt = 1
  REAL(KIND = wp), PARAMETER :: undeff_lsm = - 999.00_wp
  CONTAINS
  SUBROUTINE fld_read(kt, kn_fsbc, sd, map, kit, kt_offset, jpk_bdy, fvl)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kn_fsbc
    TYPE(FLD), INTENT(INOUT), DIMENSION(:) :: sd
    TYPE(MAP_POINTER), INTENT(IN), OPTIONAL, DIMENSION(:) :: map
    INTEGER, INTENT(IN), OPTIONAL :: kit
    INTEGER, INTENT(IN), OPTIONAL :: kt_offset
    INTEGER, INTENT(IN), OPTIONAL :: jpk_bdy
    LOGICAL, INTENT(IN), OPTIONAL :: fvl
    INTEGER :: itmp
    INTEGER :: imf
    INTEGER :: jf
    INTEGER :: isecend
    INTEGER :: isecsbc
    INTEGER :: it_offset
    LOGICAL :: llnxtyr
    LOGICAL :: llnxtmth
    LOGICAL :: llstop
    LOGICAL :: ll_firstcall
    REAL(KIND = wp) :: ztinta
    REAL(KIND = wp) :: ztintb
    CHARACTER(LEN = 1000) :: clfmt
    TYPE(MAP_POINTER) :: imap
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('fld_read', 'r0', 0, 0)
    ll_firstcall = kt == nit000
    IF (PRESENT(kit)) ll_firstcall = ll_firstcall .AND. kit == 1
    IF (nn_components == jp_iam_sas) THEN
      it_offset = nn_fsbc
    ELSE
      it_offset = 0
    END IF
    IF (PRESENT(kt_offset)) it_offset = kt_offset
    imap % ptr => NULL()
    IF (PRESENT(kit)) THEN
      isecsbc = nsec_year + nsec1jan000 + (kit + it_offset) * NINT(rdt / REAL(nn_baro, wp))
    ELSE
      isecsbc = nsec_year + nsec1jan000 + NINT(0.5 * REAL(kn_fsbc - 1, wp) * rdt) + it_offset * NINT(rdt)
    END IF
    imf = SIZE(sd)
    IF (ll_firstcall) THEN
      DO jf = 1, imf
        IF (PRESENT(map)) imap = map(jf)
        IF (PRESENT(jpk_bdy)) THEN
          CALL fld_init(kn_fsbc, sd(jf), imap, jpk_bdy, fvl)
        ELSE
          CALL fld_init(kn_fsbc, sd(jf), imap)
        END IF
      END DO
      IF (lwp) CALL wgt_print
    END IF
    IF (MOD(kt - 1, kn_fsbc) == 0) THEN
      DO jf = 1, imf
        IF (isecsbc > sd(jf) % nrec_a(2) .OR. ll_firstcall) THEN
          IF (PRESENT(map)) imap = map(jf)
          sd(jf) % nrec_b(:) = sd(jf) % nrec_a(:)
          sd(jf) % rotn(1) = sd(jf) % rotn(2)
          IF (sd(jf) % ln_tint) sd(jf) % fdta(:, :, :, 1) = sd(jf) % fdta(:, :, :, 2)
          CALL fld_rec(kn_fsbc, sd(jf), kt_offset = it_offset, kit = kit)
          IF (.NOT. ll_firstcall .AND. sd(jf) % ln_tint .AND. sd(jf) % nrec_b(1) /= sd(jf) % nreclast .AND. MOD(sd(jf) % nrec_a(1), sd(jf) % nreclast) == 1) THEN
            itmp = sd(jf) % nrec_a(1)
            sd(jf) % nrec_a(1) = sd(jf) % nreclast
            CALL fld_get(sd(jf), imap)
            sd(jf) % fdta(:, :, :, 1) = sd(jf) % fdta(:, :, :, 2)
            sd(jf) % nrec_b(1) = sd(jf) % nrec_a(1)
            sd(jf) % nrec_b(2) = sd(jf) % nrec_a(2) - NINT(sd(jf) % nfreqh * 3600)
            sd(jf) % rotn(1) = sd(jf) % rotn(2)
            sd(jf) % nrec_a(1) = itmp
          END IF
          CALL fld_clopn(sd(jf))
          IF (sd(jf) % ln_tint) THEN
            IF (.NOT. ll_firstcall .AND. MOD(sd(jf) % nrec_a(1), sd(jf) % nreclast) /= 1 .AND. sd(jf) % nrec_b(1) /= sd(jf) % nrec_a(1) - 1) THEN
              sd(jf) % nrec_a(1) = sd(jf) % nrec_a(1) - 1
              CALL fld_get(sd(jf), imap)
              sd(jf) % fdta(:, :, :, 1) = sd(jf) % fdta(:, :, :, 2)
              sd(jf) % nrec_b(1) = sd(jf) % nrec_a(1)
              sd(jf) % nrec_b(2) = sd(jf) % nrec_a(2) - NINT(sd(jf) % nfreqh * 3600)
              sd(jf) % rotn(1) = sd(jf) % rotn(2)
              sd(jf) % nrec_a(1) = sd(jf) % nrec_a(1) + 1
            END IF
          END IF
          IF (sd(jf) % nrec_a(1) > sd(jf) % nreclast) THEN
            sd(jf) % nrec_a(1) = sd(jf) % nrec_a(1) - sd(jf) % nreclast
            IF (.NOT. (sd(jf) % ln_clim .AND. sd(jf) % cltype == 'yearly')) THEN
              llnxtmth = sd(jf) % cltype == 'monthly' .OR. nday == nmonth_len(nmonth)
              llnxtyr = sd(jf) % cltype == 'yearly' .OR. (nmonth == 12 .AND. llnxtmth)
              isecend = nsec_year + nsec1jan000 + (nitend - kt) * NINT(rdt)
              llstop = isecend > sd(jf) % nrec_a(2)
              CALL fld_clopn(sd(jf), nyear + COUNT((/llnxtyr/)), nmonth + COUNT((/llnxtmth/)) - 12 * COUNT((/llnxtyr/)), nday + 1 - nmonth_len(nmonth) * COUNT((/llnxtmth/)), llstop)
              IF (sd(jf) % num <= 0 .AND. .NOT. llstop) THEN
                CALL ctl_warn('next year/month/week/day file: ' // TRIM(sd(jf) % clname) // ' not present -> back to current year/month/day')
                CALL fld_clopn(sd(jf))
                sd(jf) % nrec_a(1) = sd(jf) % nreclast
              END IF
            END IF
          END IF
          IF (PRESENT(jpk_bdy)) THEN
            CALL fld_get(sd(jf), imap, jpk_bdy, fvl)
          ELSE
            CALL fld_get(sd(jf), imap)
          END IF
        END IF
      END DO
      CALL fld_rot(kt, sd)
      DO jf = 1, imf
        IF (sd(jf) % ln_tint) THEN
          IF (lwp .AND. kt - nit000 <= 100) THEN
            clfmt = "('   fld_read: var ', a, ' kt = ', i8, ' (', f9.4,' days), Y/M/D = ', i4.4,'/', i2.2,'/', i2.2," // "', records b/a: ', i6.4, '/', i6.4, ' (days ', f9.4,'/', f9.4, ')')"
            WRITE(numout, clfmt) TRIM(sd(jf) % clvar), kt, REAL(isecsbc, wp) / rday, nyear, nmonth, nday, sd(jf) % nrec_b(1), sd(jf) % nrec_a(1), REAL(sd(jf) % nrec_b(2), wp) / rday, REAL(sd(jf) % nrec_a(2), wp) / rday
            WRITE(numout, *) '      it_offset is : ', it_offset
          END IF
          ztinta = REAL(isecsbc - sd(jf) % nrec_b(2), wp) / REAL(sd(jf) % nrec_a(2) - sd(jf) % nrec_b(2), wp)
          ztintb = 1. - ztinta
          sd(jf) % fnow(:, :, :) = ztintb * sd(jf) % fdta(:, :, :, 1) + ztinta * sd(jf) % fdta(:, :, :, 2)
        ELSE
          IF (lwp .AND. kt - nit000 <= 100) THEN
            clfmt = "('   fld_read: var ', a, ' kt = ', i8,' (', f9.4,' days), Y/M/D = ', i4.4,'/', i2.2,'/', i2.2," // "', record: ', i6.4, ' (days ', f9.4, ' <-> ', f9.4, ')')"
            WRITE(numout, clfmt) TRIM(sd(jf) % clvar), kt, REAL(isecsbc, wp) / rday, nyear, nmonth, nday, sd(jf) % nrec_a(1), REAL(sd(jf) % nrec_b(2), wp) / rday, REAL(sd(jf) % nrec_a(2), wp) / rday
          END IF
        END IF
        IF (kt == nitend - kn_fsbc + 1) CALL iom_close(sd(jf) % num)
      END DO
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE fld_read
  SUBROUTINE fld_init(kn_fsbc, sdjf, map, jpk_bdy, fvl)
    INTEGER, INTENT(IN) :: kn_fsbc
    TYPE(FLD), INTENT(INOUT) :: sdjf
    TYPE(MAP_POINTER), INTENT(IN) :: map
    INTEGER, INTENT(IN), OPTIONAL :: jpk_bdy
    LOGICAL, INTENT(IN), OPTIONAL :: fvl
    LOGICAL :: llprevyr
    LOGICAL :: llprevmth
    LOGICAL :: llprevweek
    LOGICAL :: llprevday
    LOGICAL :: llprev
    INTEGER :: idvar
    INTEGER :: inrec
    INTEGER :: iyear, imonth, iday
    INTEGER :: isec_week
    CHARACTER(LEN = 1000) :: clfmt
    llprevyr = .FALSE.
    llprevmth = .FALSE.
    llprevweek = .FALSE.
    llprevday = .FALSE.
    isec_week = 0
    CALL fld_rec(kn_fsbc, sdjf, ldbefore = .TRUE.)
    IF (sdjf % ln_tint) THEN
      IF (sdjf % nrec_a(1) == 0) THEN
        IF (sdjf % nfreqh == - 12) THEN
          IF (sdjf % cltype == 'yearly') THEN
            sdjf % nrec_a(1) = 1
            llprevyr = .NOT. sdjf % ln_clim
          ELSE
            CALL ctl_stop("fld_init: yearly mean file must be in a yearly type of file: " // TRIM(sdjf % clrootname))
          END IF
        ELSE IF (sdjf % nfreqh == - 1) THEN
          IF (sdjf % cltype == 'monthly') THEN
            sdjf % nrec_a(1) = 1
            llprevmth = .TRUE.
            llprevyr = llprevmth .AND. nmonth == 1
          ELSE
            sdjf % nrec_a(1) = 12
            llprevyr = .NOT. sdjf % ln_clim
          END IF
        ELSE
          IF (sdjf % cltype == 'monthly') THEN
            sdjf % nrec_a(1) = NINT(24 * nmonth_len(nmonth - 1) / sdjf % nfreqh)
            llprevmth = .TRUE.
            llprevyr = llprevmth .AND. nmonth == 1
          ELSE IF (sdjf % cltype(1 : 4) == 'week') THEN
            llprevweek = .TRUE.
            sdjf % nrec_a(1) = NINT(24 * 7 / sdjf % nfreqh)
            isec_week = NINT(rday) * 7
          ELSE IF (sdjf % cltype == 'daily') THEN
            sdjf % nrec_a(1) = NINT(24 / sdjf % nfreqh)
            llprevday = .TRUE.
            llprevmth = llprevday .AND. nday == 1
            llprevyr = llprevmth .AND. nmonth == 1
          ELSE
            sdjf % nrec_a(1) = NINT(24 * nyear_len(0) / sdjf % nfreqh)
            llprevyr = .NOT. sdjf % ln_clim
          END IF
        END IF
      END IF
      IF (sdjf % cltype(1 : 4) == 'week') THEN
        isec_week = isec_week + ksec_week(sdjf % cltype(6 : 8))
        llprevmth = isec_week > nsec_month
        llprevyr = llprevmth .AND. nmonth == 1
      END IF
      llprev = llprevyr .OR. llprevmth .OR. llprevweek .OR. llprevday
      iyear = nyear - COUNT((/llprevyr/))
      imonth = nmonth - COUNT((/llprevmth/)) + 12 * COUNT((/llprevyr/))
      iday = nday - COUNT((/llprevday/)) + nmonth_len(nmonth - 1) * COUNT((/llprevmth/)) - isec_week / NINT(rday)
      CALL fld_clopn(sdjf, iyear, imonth, iday, .NOT. llprev)
      IF (llprev .AND. sdjf % num <= 0) THEN
        CALL ctl_warn('previous year/month/week/day file: ' // TRIM(sdjf % clrootname) // ' not present -> back to current year/month/week/day')
        llprev = .FALSE.
        sdjf % nrec_a(1) = 1
        CALL fld_clopn(sdjf)
      END IF
      IF (llprev) THEN
        idvar = iom_varid(sdjf % num, sdjf % clvar)
        IF (idvar <= 0) RETURN
        inrec = iom_file(sdjf % num) % dimsz(iom_file(sdjf % num) % ndims(idvar), idvar)
        sdjf % nrec_a(1) = MIN(sdjf % nrec_a(1), inrec)
      END IF
      IF (PRESENT(jpk_bdy)) THEN
        CALL fld_get(sdjf, map, jpk_bdy, fvl)
      ELSE
        CALL fld_get(sdjf, map)
      END IF
      clfmt = "('   fld_init : time-interpolation for ', a, ' read previous record = ', i6, ' at time = ', f7.2, ' days')"
      IF (lwp) WRITE(numout, clfmt) TRIM(sdjf % clvar), sdjf % nrec_a(1), REAL(sdjf % nrec_a(2), wp) / rday
    END IF
  END SUBROUTINE fld_init
  SUBROUTINE fld_rec(kn_fsbc, sdjf, ldbefore, kit, kt_offset)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kn_fsbc
    TYPE(FLD), INTENT(INOUT) :: sdjf
    LOGICAL, INTENT(IN), OPTIONAL :: ldbefore
    INTEGER, INTENT(IN), OPTIONAL :: kit
    INTEGER, INTENT(IN), OPTIONAL :: kt_offset
    LOGICAL :: llbefore
    INTEGER :: iendrec
    INTEGER :: imth
    INTEGER :: ifreq_sec
    INTEGER :: isec_week
    INTEGER :: it_offset
    REAL(KIND = wp) :: ztmp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('fld_rec', 'r0', 0, 0)
    IF (PRESENT(ldbefore)) THEN
      llbefore = ldbefore .AND. sdjf % ln_tint
    ELSE
      llbefore = .FALSE.
    END IF
    IF (nn_components == jp_iam_sas) THEN
      it_offset = nn_fsbc
    ELSE
      it_offset = 0
    END IF
    IF (PRESENT(kt_offset)) it_offset = kt_offset
    IF (PRESENT(kit)) THEN
      it_offset = (kit + it_offset) * NINT(rdt / REAL(nn_baro, wp))
    ELSE
      it_offset = it_offset * NINT(rdt)
    END IF
    IF (sdjf % nfreqh == - 12) THEN
      IF (sdjf % ln_tint) THEN
        ztmp = REAL(nsec_year, wp) / (REAL(nyear_len(1), wp) * rday) + 0.5 + REAL(it_offset, wp) / (REAL(nyear_len(1), wp) * rday)
        sdjf % nrec_a(1) = 1 + INT(ztmp) - COUNT((/llbefore/))
        IF (llbefore) THEN
          sdjf % nrec_a(2) = nsec1jan000 - (1 - INT(ztmp)) * NINT(0.5 * rday) * nyear_len(0) + INT(ztmp) * NINT(0.5 * rday) * nyear_len(1)
        ELSE
          sdjf % nrec_a(2) = nsec1jan000 + (1 - INT(ztmp)) * NINT(0.5 * rday) * nyear_len(1) + INT(ztmp) * INT(rday) * nyear_len(1) + INT(ztmp) * NINT(0.5 * rday) * nyear_len(2)
        END IF
      ELSE
        sdjf % nrec_a(1) = 1
        sdjf % nrec_a(2) = NINT(rday) * nyear_len(1) + nsec1jan000
        sdjf % nrec_b(2) = nsec1jan000
      END IF
    ELSE IF (sdjf % nfreqh == - 1) THEN
      IF (sdjf % ln_tint) THEN
        ztmp = REAL(nsec_month, wp) / (REAL(nmonth_len(nmonth), wp) * rday) + 0.5 + REAL(it_offset, wp) / (REAL(nmonth_len(nmonth), wp) * rday)
        imth = nmonth + INT(ztmp) - COUNT((/llbefore/))
        IF (sdjf % cltype == 'monthly') THEN
          sdjf % nrec_a(1) = 1 + INT(ztmp) - COUNT((/llbefore/))
        ELSE
          sdjf % nrec_a(1) = imth
        END IF
        sdjf % nrec_a(2) = nmonth_half(imth) + nsec1jan000
      ELSE
        IF (sdjf % cltype == 'monthly') THEN
          sdjf % nrec_a(1) = 1
        ELSE
          sdjf % nrec_a(1) = nmonth
        END IF
        sdjf % nrec_a(2) = nmonth_end(nmonth) + nsec1jan000
        sdjf % nrec_b(2) = nmonth_end(nmonth - 1) + nsec1jan000
      END IF
    ELSE
      ifreq_sec = NINT(sdjf % nfreqh * 3600)
      IF (sdjf % cltype(1 : 4) == 'week') isec_week = ksec_week(sdjf % cltype(6 : 8))
      IF (sdjf % cltype == 'monthly') THEN
        ztmp = REAL(nsec_month, wp)
      ELSE IF (sdjf % cltype(1 : 4) == 'week') THEN
        ztmp = REAL(isec_week, wp)
      ELSE IF (sdjf % cltype == 'daily') THEN
        ztmp = REAL(nsec_day, wp)
      ELSE
        ztmp = REAL(nsec_year, wp)
      END IF
      ztmp = ztmp + 0.5 * REAL(kn_fsbc - 1, wp) * rdt + REAL(it_offset, wp)
      ztmp = ztmp + 0.01 * rdt
      IF (sdjf % ln_tint) THEN
        ztmp = ztmp / REAL(ifreq_sec, wp) + 0.5
      ELSE
        ztmp = ztmp / REAL(ifreq_sec, wp)
      END IF
      sdjf % nrec_a(1) = 1 + INT(ztmp) - COUNT((/llbefore/))
      iendrec = ifreq_sec * sdjf % nrec_a(1) + nsec1jan000
      IF (sdjf % cltype == 'monthly') iendrec = iendrec + NINT(rday) * SUM(nmonth_len(1 : nmonth - 1))
      IF (sdjf % cltype(1 : 4) == 'week') iendrec = iendrec + (nsec_year - isec_week)
      IF (sdjf % cltype == 'daily') iendrec = iendrec + NINT(rday) * (nday_year - 1)
      IF (sdjf % ln_tint) THEN
        sdjf % nrec_a(2) = iendrec - ifreq_sec / 2
      ELSE
        sdjf % nrec_a(2) = iendrec
        sdjf % nrec_b(2) = iendrec - ifreq_sec
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE fld_rec
  SUBROUTINE fld_get(sdjf, map, jpk_bdy, fvl)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(FLD), INTENT(INOUT) :: sdjf
    TYPE(MAP_POINTER), INTENT(IN) :: map
    INTEGER, INTENT(IN), OPTIONAL :: jpk_bdy
    LOGICAL, INTENT(IN), OPTIONAL :: fvl
    INTEGER :: ipk
    INTEGER :: iw
    INTEGER :: ipdom
    INTEGER :: idvar
    INTEGER :: idmspc
    LOGICAL :: lmoor
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('fld_get', 'r0', 0, 0)
    ipk = SIZE(sdjf % fnow, 3)
    IF (ASSOCIATED(map % ptr)) THEN
      IF (PRESENT(jpk_bdy)) THEN
        IF (sdjf % ln_tint) THEN
          CALL fld_map(sdjf % num, sdjf % clvar, sdjf % fdta(:, :, :, 2), sdjf % nrec_a(1), map, sdjf % igrd, sdjf % ibdy, jpk_bdy, fvl)
        ELSE
          CALL fld_map(sdjf % num, sdjf % clvar, sdjf % fnow(:, :, :), sdjf % nrec_a(1), map, sdjf % igrd, sdjf % ibdy, jpk_bdy, fvl)
        END IF
      ELSE
        IF (sdjf % ln_tint) THEN
          CALL fld_map(sdjf % num, sdjf % clvar, sdjf % fdta(:, :, :, 2), sdjf % nrec_a(1), map)
        ELSE
          CALL fld_map(sdjf % num, sdjf % clvar, sdjf % fnow(:, :, :), sdjf % nrec_a(1), map)
        END IF
      END IF
    ELSE IF (LEN(TRIM(sdjf % wgtname)) > 0) THEN
      CALL wgt_list(sdjf, iw)
      IF (sdjf % ln_tint) THEN
        CALL fld_interp(sdjf % num, sdjf % clvar, iw, ipk, sdjf % fdta(:, :, :, 2), sdjf % nrec_a(1), sdjf % lsmname)
      ELSE
        CALL fld_interp(sdjf % num, sdjf % clvar, iw, ipk, sdjf % fnow(:, :, :), sdjf % nrec_a(1), sdjf % lsmname)
      END IF
    ELSE
      IF (SIZE(sdjf % fnow, 1) == jpi) THEN
        ipdom = jpdom_data
      ELSE
        ipdom = jpdom_unknown
      END IF
      idvar = iom_varid(sdjf % num, sdjf % clvar)
      idmspc = iom_file(sdjf % num) % ndims(idvar)
      IF (iom_file(sdjf % num) % luld(idvar)) idmspc = idmspc - 1
      lmoor = (idmspc == 0 .OR. PRODUCT(iom_file(sdjf % num) % dimsz(1 : MAX(idmspc, 1), idvar)) == ipk)
      SELECT CASE (ipk)
      CASE (1)
        IF (lk_c1d .AND. lmoor) THEN
          IF (sdjf % ln_tint) THEN
            CALL iom_get(sdjf % num, sdjf % clvar, sdjf % fdta(2, 2, 1, 2), sdjf % nrec_a(1))
            CALL lbc_lnk('fldread', sdjf % fdta(:, :, 1, 2), 'Z', 1.)
          ELSE
            CALL iom_get(sdjf % num, sdjf % clvar, sdjf % fnow(2, 2, 1), sdjf % nrec_a(1))
            CALL lbc_lnk('fldread', sdjf % fnow(:, :, 1), 'Z', 1.)
          END IF
        ELSE
          IF (sdjf % ln_tint) THEN
            CALL iom_get(sdjf % num, ipdom, sdjf % clvar, sdjf % fdta(:, :, 1, 2), sdjf % nrec_a(1))
          ELSE
            CALL iom_get(sdjf % num, ipdom, sdjf % clvar, sdjf % fnow(:, :, 1), sdjf % nrec_a(1))
          END IF
        END IF
      CASE DEFAULT
        IF (lk_c1d .AND. lmoor) THEN
          IF (sdjf % ln_tint) THEN
            CALL iom_get(sdjf % num, jpdom_unknown, sdjf % clvar, sdjf % fdta(2, 2, :, 2), sdjf % nrec_a(1))
            CALL lbc_lnk('fldread', sdjf % fdta(:, :, :, 2), 'Z', 1.)
          ELSE
            CALL iom_get(sdjf % num, jpdom_unknown, sdjf % clvar, sdjf % fnow(2, 2, :), sdjf % nrec_a(1))
            CALL lbc_lnk('fldread', sdjf % fnow(:, :, :), 'Z', 1.)
          END IF
        ELSE
          IF (sdjf % ln_tint) THEN
            CALL iom_get(sdjf % num, ipdom, sdjf % clvar, sdjf % fdta(:, :, :, 2), sdjf % nrec_a(1))
          ELSE
            CALL iom_get(sdjf % num, ipdom, sdjf % clvar, sdjf % fnow(:, :, :), sdjf % nrec_a(1))
          END IF
        END IF
      END SELECT
    END IF
    sdjf % rotn(2) = .FALSE.
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE fld_get
  SUBROUTINE fld_map(num, clvar, dta, nrec, map, igrd, ibdy, jpk_bdy, fvl)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE bdy_oce, ONLY: ln_bdy, idx_bdy, dta_global, dta_global_z, dta_global_dz, dta_global2, dta_global2_z, dta_global2_dz
    INTEGER, INTENT(IN) :: num
    CHARACTER(LEN = *), INTENT(IN) :: clvar
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: dta
    INTEGER, INTENT(IN) :: nrec
    TYPE(MAP_POINTER), INTENT(IN) :: map
    INTEGER, INTENT(IN), OPTIONAL :: igrd, ibdy, jpk_bdy
    LOGICAL, INTENT(IN), OPTIONAL :: fvl
    INTEGER :: jpkm1_bdy
    INTEGER :: ipi
    INTEGER :: ipj
    INTEGER :: ipk
    INTEGER :: ilendta
    INTEGER :: idvar
    INTEGER :: ib, ik, ji, jj
    INTEGER :: ierr
    REAL(KIND = wp) :: fv
    REAL(KIND = wp), POINTER, DIMENSION(:, :, :) :: dta_read
    REAL(KIND = wp), POINTER, DIMENSION(:, :, :) :: dta_read_z
    REAL(KIND = wp), POINTER, DIMENSION(:, :, :) :: dta_read_dz
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('fld_map', 'r0', 0, 0)
    ipi = SIZE(dta, 1)
    ipj = 1
    ipk = SIZE(dta, 3)
    idvar = iom_varid(num, clvar)
    ilendta = iom_file(num) % dimsz(1, idvar)
    IF (ln_bdy) THEN
      ipj = iom_file(num) % dimsz(2, idvar)
      IF (map % ll_unstruc) THEN
        dta_read => dta_global
        IF (PRESENT(jpk_bdy)) THEN
          IF (jpk_bdy > 0) THEN
            dta_read_z => dta_global_z
            dta_read_dz => dta_global_dz
            jpkm1_bdy = jpk_bdy - 1
          END IF
        END IF
      ELSE
        dta_read => dta_global2
        IF (PRESENT(jpk_bdy)) THEN
          IF (jpk_bdy > 0) THEN
            dta_read_z => dta_global2_z
            dta_read_dz => dta_global2_dz
            jpkm1_bdy = jpk_bdy - 1
          END IF
        END IF
      END IF
    END IF
    IF (lwp) WRITE(numout, *) 'Dim size for ', TRIM(clvar), ' is ', ilendta
    IF (lwp) WRITE(numout, *) 'Number of levels for ', TRIM(clvar), ' is ', ipk
    SELECT CASE (ipk)
    CASE (1)
      CALL iom_get(num, jpdom_unknown, clvar, dta_read(1 : ilendta, 1 : ipj, 1), nrec)
      IF (map % ll_unstruc) THEN
        DO ib = 1, ipi
          DO ik = 1, ipk
            dta(ib, 1, ik) = dta_read(map % ptr(ib), 1, ik)
          END DO
        END DO
      ELSE
        DO ib = 1, ipi
          jj = 1 + FLOOR(REAL(map % ptr(ib) - 1) / REAL(ilendta))
          ji = map % ptr(ib) - (jj - 1) * ilendta
          DO ik = 1, ipk
            dta(ib, 1, ik) = dta_read(ji, jj, ik)
          END DO
        END DO
      END IF
    CASE DEFAULT
      IF (PRESENT(jpk_bdy) .AND. jpk_bdy > 0) THEN
        CALL iom_getatt(num, '_FillValue', fv, cdvar = clvar)
        dta_read(:, :, :) = - ABS(fv)
        dta_read_z(:, :, :) = 0._wp
        dta_read_dz(:, :, :) = 0._wp
        CALL iom_get(num, jpdom_unknown, clvar, dta_read(1 : ilendta, 1 : ipj, 1 : jpk_bdy), nrec)
        SELECT CASE (igrd)
        CASE (1)
          CALL iom_get(num, jpdom_unknown, 'gdept', dta_read_z(1 : ilendta, 1 : ipj, 1 : jpk_bdy))
          CALL iom_get(num, jpdom_unknown, 'e3t', dta_read_dz(1 : ilendta, 1 : ipj, 1 : jpk_bdy))
        CASE (2)
          CALL iom_get(num, jpdom_unknown, 'gdepu', dta_read_z(1 : ilendta, 1 : ipj, 1 : jpk_bdy))
          CALL iom_get(num, jpdom_unknown, 'e3u', dta_read_dz(1 : ilendta, 1 : ipj, 1 : jpk_bdy))
        CASE (3)
          CALL iom_get(num, jpdom_unknown, 'gdepv', dta_read_z(1 : ilendta, 1 : ipj, 1 : jpk_bdy))
          CALL iom_get(num, jpdom_unknown, 'e3v', dta_read_dz(1 : ilendta, 1 : ipj, 1 : jpk_bdy))
        END SELECT
        IF (ln_bdy) CALL fld_bdy_interp(dta_read, dta_read_z, dta_read_dz, map, jpk_bdy, igrd, ibdy, fv, dta, fvl, ilendta)
      ELSE
        CALL iom_get(num, jpdom_unknown, clvar, dta_read(1 : ilendta, 1 : ipj, 1 : ipk), nrec)
        IF (map % ll_unstruc) THEN
          DO ib = 1, ipi
            DO ik = 1, ipk
              dta(ib, 1, ik) = dta_read(map % ptr(ib), 1, ik)
            END DO
          END DO
        ELSE
          DO ib = 1, ipi
            jj = 1 + FLOOR(REAL(map % ptr(ib) - 1) / REAL(ilendta))
            ji = map % ptr(ib) - (jj - 1) * ilendta
            DO ik = 1, ipk
              dta(ib, 1, ik) = dta_read(ji, jj, ik)
            END DO
          END DO
        END IF
      END IF
    END SELECT
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE fld_map
  SUBROUTINE fld_bdy_interp(dta_read, dta_read_z, dta_read_dz, map, jpk_bdy, igrd, ibdy, fv, dta, fvl, ilendta)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE bdy_oce, ONLY: idx_bdy
    REAL(KIND = wp), POINTER, DIMENSION(:, :, :), INTENT(IN) :: dta_read
    REAL(KIND = wp), POINTER, DIMENSION(:, :, :), INTENT(IN) :: dta_read_z
    REAL(KIND = wp), POINTER, DIMENSION(:, :, :), INTENT(IN) :: dta_read_dz
    REAL(KIND = wp), INTENT(IN) :: fv
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: dta
    TYPE(MAP_POINTER), INTENT(IN) :: map
    LOGICAL, INTENT(IN), OPTIONAL :: fvl
    INTEGER, INTENT(IN) :: igrd, ibdy, jpk_bdy
    INTEGER, INTENT(IN) :: ilendta
    INTEGER :: ipi
    INTEGER :: ipj
    INTEGER :: ipk
    INTEGER :: jpkm1_bdy
    INTEGER :: ib, ik, ikk
    INTEGER :: ji, jj, zij, zjj
    REAL(KIND = wp) :: zl, zi, zh
    REAL(KIND = wp) :: fv_alt, ztrans, ztrans_new
    CHARACTER(LEN = 10) :: ibstr
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('fld_bdy_interp', 'r0', 0, 0)
    ipi = SIZE(dta, 1)
    ipj = SIZE(dta_read, 2)
    ipk = SIZE(dta, 3)
    jpkm1_bdy = jpk_bdy - 1
    fv_alt = - ABS(fv)
    DO ib = 1, ipi
      zij = idx_bdy(ibdy) % nbi(ib, igrd)
      zjj = idx_bdy(ibdy) % nbj(ib, igrd)
      IF (narea == 2) WRITE(*, *) 'MAPI', ib, igrd, map % ptr(ib), narea - 1, zij, zjj
    END DO
    IF (map % ll_unstruc) THEN
      DO ib = 1, ipi
        DO ik = 1, jpk_bdy
          IF ((dta_read(map % ptr(ib), 1, ik) == fv)) THEN
            dta_read_z(map % ptr(ib), 1, ik) = fv_alt
            dta_read_dz(map % ptr(ib), 1, ik) = 0._wp
          END IF
        END DO
      END DO
      DO ib = 1, ipi
        zij = idx_bdy(ibdy) % nbi(ib, igrd)
        zjj = idx_bdy(ibdy) % nbj(ib, igrd)
        zh = SUM(dta_read_dz(map % ptr(ib), 1, :))
        SELECT CASE (igrd)
        CASE (1)
          IF (ABS((zh - ht_n(zij, zjj)) / ht_n(zij, zjj)) * tmask(zij, zjj, 1) > 0.01_wp) THEN
            WRITE(ibstr, "(I10.10)") map % ptr(ib)
            CALL ctl_warn('fld_bdy_interp: T depths differ between grids at BDY point ' // TRIM(ibstr) // ' by more than 1%')
          END IF
        CASE (2)
          IF (ABS((zh - hu_n(zij, zjj)) * r1_hu_n(zij, zjj)) * umask(zij, zjj, 1) > 0.01_wp) THEN
            WRITE(ibstr, "(I10.10)") map % ptr(ib)
            CALL ctl_warn('fld_bdy_interp: U depths differ between grids at BDY point ' // TRIM(ibstr) // ' by more than 1%')
            IF (lwp) WRITE(*, *) 'DEPTHU', zh, SUM(e3u_n(zij, zjj, :), mask = umask(zij, zjj, :) == 1), SUM(umask(zij, zjj, :)), hu_n(zij, zjj), map % ptr(ib), ib, zij, zjj, narea - 1, dta_read(map % ptr(ib), 1, :)
          END IF
        CASE (3)
          IF (ABS((zh - hv_n(zij, zjj)) * r1_hv_n(zij, zjj)) * vmask(zij, zjj, 1) > 0.01_wp) THEN
            WRITE(ibstr, "(I10.10)") map % ptr(ib)
            CALL ctl_warn('fld_bdy_interp: V depths differ between grids at BDY point ' // TRIM(ibstr) // ' by more than 1%')
          END IF
        END SELECT
        DO ik = 1, ipk
          SELECT CASE (igrd)
          CASE (1)
            zl = gdept_n(zij, zjj, ik)
          CASE (2)
            IF (ln_sco) THEN
              zl = (gdept_n(zij, zjj, ik) + gdept_n(zij + 1, zjj, ik)) * 0.5_wp
            ELSE
              zl = MIN(gdept_n(zij, zjj, ik), gdept_n(zij + 1, zjj, ik))
            END IF
          CASE (3)
            IF (ln_sco) THEN
              zl = (gdept_n(zij, zjj, ik) + gdept_n(zij, zjj + 1, ik)) * 0.5_wp
            ELSE
              zl = MIN(gdept_n(zij, zjj, ik), gdept_n(zij, zjj + 1, ik))
            END IF
          END SELECT
          IF (zl < dta_read_z(map % ptr(ib), 1, 1)) THEN
            dta(ib, 1, ik) = dta_read(map % ptr(ib), 1, 1)
          ELSE IF (zl > MAXVAL(dta_read_z(map % ptr(ib), 1, :), 1)) THEN
            dta(ib, 1, ik) = dta_read(map % ptr(ib), 1, MAXLOC(dta_read_z(map % ptr(ib), 1, :), 1))
          ELSE
            DO ikk = 1, jpkm1_bdy
              IF (((zl - dta_read_z(map % ptr(ib), 1, ikk)) * (zl - dta_read_z(map % ptr(ib), 1, ikk + 1)) <= 0._wp) .AND. (dta_read_z(map % ptr(ib), 1, ikk + 1) /= fv_alt)) THEN
                zi = (zl - dta_read_z(map % ptr(ib), 1, ikk)) / (dta_read_z(map % ptr(ib), 1, ikk + 1) - dta_read_z(map % ptr(ib), 1, ikk))
                dta(ib, 1, ik) = dta_read(map % ptr(ib), 1, ikk) + (dta_read(map % ptr(ib), 1, ikk + 1) - dta_read(map % ptr(ib), 1, ikk)) * zi
              END IF
            END DO
          END IF
        END DO
      END DO
      IF (igrd == 2) THEN
        DO ib = 1, ipi
          zij = idx_bdy(ibdy) % nbi(ib, igrd)
          zjj = idx_bdy(ibdy) % nbj(ib, igrd)
          zh = SUM(dta_read_dz(map % ptr(ib), 1, :))
          ztrans = 0._wp
          ztrans_new = 0._wp
          DO ik = 1, jpk_bdy
            ztrans = ztrans + dta_read(map % ptr(ib), 1, ik) * dta_read_dz(map % ptr(ib), 1, ik)
          END DO
          DO ik = 1, ipk
            ztrans_new = ztrans_new + dta(ib, 1, ik) * e3u_n(zij, zjj, ik) * umask(zij, zjj, ik)
          END DO
          DO ik = 1, ipk
            IF (fvl) THEN
              dta(ib, 1, ik) = (dta(ib, 1, ik) + (ztrans - ztrans_new) * r1_hu_n(zij, zjj)) * umask(zij, zjj, ik)
            ELSE
              IF (ABS(ztrans * r1_hu_n(zij, zjj)) > 0.01_wp) CALL ctl_warn('fld_bdy_interp: barotropic component of > 0.01 ms-1 found in baroclinic velocities at')
              dta(ib, 1, ik) = dta(ib, 1, ik) + (0._wp - ztrans_new) * r1_hu_n(zij, zjj) * umask(zij, zjj, ik)
            END IF
          END DO
        END DO
      END IF
      IF (igrd == 3) THEN
        DO ib = 1, ipi
          zij = idx_bdy(ibdy) % nbi(ib, igrd)
          zjj = idx_bdy(ibdy) % nbj(ib, igrd)
          zh = SUM(dta_read_dz(map % ptr(ib), 1, :))
          ztrans = 0._wp
          ztrans_new = 0._wp
          DO ik = 1, jpk_bdy
            ztrans = ztrans + dta_read(map % ptr(ib), 1, ik) * dta_read_dz(map % ptr(ib), 1, ik)
          END DO
          DO ik = 1, ipk
            ztrans_new = ztrans_new + dta(ib, 1, ik) * e3v_n(zij, zjj, ik) * vmask(zij, zjj, ik)
          END DO
          DO ik = 1, ipk
            IF (fvl) THEN
              dta(ib, 1, ik) = (dta(ib, 1, ik) + (ztrans - ztrans_new) * r1_hv_n(zij, zjj)) * vmask(zij, zjj, ik)
            ELSE
              dta(ib, 1, ik) = dta(ib, 1, ik) + (0._wp - ztrans_new) * r1_hv_n(zij, zjj) * vmask(zij, zjj, ik)
            END IF
          END DO
        END DO
      END IF
    ELSE
      DO ib = 1, ipi
        jj = 1 + FLOOR(REAL(map % ptr(ib) - 1) / REAL(ilendta))
        ji = map % ptr(ib) - (jj - 1) * ilendta
        DO ik = 1, jpk_bdy
          IF ((dta_read(ji, jj, ik) == fv)) THEN
            dta_read_z(ji, jj, ik) = fv_alt
            dta_read_dz(ji, jj, ik) = 0._wp
          END IF
        END DO
      END DO
      DO ib = 1, ipi
        jj = 1 + FLOOR(REAL(map % ptr(ib) - 1) / REAL(ilendta))
        ji = map % ptr(ib) - (jj - 1) * ilendta
        zij = idx_bdy(ibdy) % nbi(ib, igrd)
        zjj = idx_bdy(ibdy) % nbj(ib, igrd)
        zh = SUM(dta_read_dz(ji, jj, :))
        SELECT CASE (igrd)
        CASE (1)
          IF (ABS((zh - ht_n(zij, zjj)) / ht_n(zij, zjj)) * tmask(zij, zjj, 1) > 0.01_wp) THEN
            WRITE(ibstr, "(I10.10)") map % ptr(ib)
            CALL ctl_warn('fld_bdy_interp: T depths differ between grids at BDY point ' // TRIM(ibstr) // ' by more than 1%')
          END IF
        CASE (2)
          IF (ABS((zh - hu_n(zij, zjj)) * r1_hu_n(zij, zjj)) * umask(zij, zjj, 1) > 0.01_wp) THEN
            WRITE(ibstr, "(I10.10)") map % ptr(ib)
            CALL ctl_warn('fld_bdy_interp: U depths differ between grids at BDY point ' // TRIM(ibstr) // ' by more than 1%')
          END IF
        CASE (3)
          IF (ABS((zh - hv_n(zij, zjj)) * r1_hv_n(zij, zjj)) * vmask(zij, zjj, 1) > 0.01_wp) THEN
            WRITE(ibstr, "(I10.10)") map % ptr(ib)
            CALL ctl_warn('fld_bdy_interp: V depths differ between grids at BDY point ' // TRIM(ibstr) // ' by more than 1%')
          END IF
        END SELECT
        DO ik = 1, ipk
          SELECT CASE (igrd)
          CASE (1)
            zl = gdept_n(zij, zjj, ik)
          CASE (2)
            IF (ln_sco) THEN
              zl = (gdept_n(zij, zjj, ik) + gdept_n(zij + 1, zjj, ik)) * 0.5_wp
            ELSE
              zl = MIN(gdept_n(zij, zjj, ik), gdept_n(zij + 1, zjj, ik))
            END IF
          CASE (3)
            IF (ln_sco) THEN
              zl = (gdept_n(zij, zjj, ik) + gdept_n(zij, zjj + 1, ik)) * 0.5_wp
            ELSE
              zl = MIN(gdept_n(zij, zjj, ik), gdept_n(zij, zjj + 1, ik))
            END IF
          END SELECT
          IF (zl < dta_read_z(ji, jj, 1)) THEN
            dta(ib, 1, ik) = dta_read(ji, jj, 1)
          ELSE IF (zl > MAXVAL(dta_read_z(ji, jj, :), 1)) THEN
            dta(ib, 1, ik) = dta_read(ji, jj, MAXLOC(dta_read_z(ji, jj, :), 1))
          ELSE
            DO ikk = 1, jpkm1_bdy
              IF (((zl - dta_read_z(ji, jj, ikk)) * (zl - dta_read_z(ji, jj, ikk + 1)) <= 0._wp) .AND. (dta_read_z(ji, jj, ikk + 1) /= fv_alt)) THEN
                zi = (zl - dta_read_z(ji, jj, ikk)) / (dta_read_z(ji, jj, ikk + 1) - dta_read_z(ji, jj, ikk))
                dta(ib, 1, ik) = dta_read(ji, jj, ikk) + (dta_read(ji, jj, ikk + 1) - dta_read(ji, jj, ikk)) * zi
              END IF
            END DO
          END IF
        END DO
      END DO
      IF (igrd == 2) THEN
        DO ib = 1, ipi
          jj = 1 + FLOOR(REAL(map % ptr(ib) - 1) / REAL(ilendta))
          ji = map % ptr(ib) - (jj - 1) * ilendta
          zij = idx_bdy(ibdy) % nbi(ib, igrd)
          zjj = idx_bdy(ibdy) % nbj(ib, igrd)
          zh = SUM(dta_read_dz(ji, jj, :))
          ztrans = 0._wp
          ztrans_new = 0._wp
          DO ik = 1, jpk_bdy
            ztrans = ztrans + dta_read(ji, jj, ik) * dta_read_dz(ji, jj, ik)
          END DO
          DO ik = 1, ipk
            ztrans_new = ztrans_new + dta(ib, 1, ik) * e3u_n(zij, zjj, ik) * umask(zij, zjj, ik)
          END DO
          DO ik = 1, ipk
            IF (fvl) THEN
              dta(ib, 1, ik) = (dta(ib, 1, ik) + (ztrans - ztrans_new) * r1_hu_n(zij, zjj)) * umask(zij, zjj, ik)
            ELSE
              dta(ib, 1, ik) = (dta(ib, 1, ik) + (0._wp - ztrans_new) * r1_hu_n(zij, zjj)) * umask(zij, zjj, ik)
            END IF
          END DO
        END DO
      END IF
      IF (igrd == 3) THEN
        DO ib = 1, ipi
          jj = 1 + FLOOR(REAL(map % ptr(ib) - 1) / REAL(ilendta))
          ji = map % ptr(ib) - (jj - 1) * ilendta
          zij = idx_bdy(ibdy) % nbi(ib, igrd)
          zjj = idx_bdy(ibdy) % nbj(ib, igrd)
          zh = SUM(dta_read_dz(ji, jj, :))
          ztrans = 0._wp
          ztrans_new = 0._wp
          DO ik = 1, jpk_bdy
            ztrans = ztrans + dta_read(ji, jj, ik) * dta_read_dz(ji, jj, ik)
          END DO
          DO ik = 1, ipk
            ztrans_new = ztrans_new + dta(ib, 1, ik) * e3v_n(zij, zjj, ik) * vmask(zij, zjj, ik)
          END DO
          DO ik = 1, ipk
            IF (fvl) THEN
              dta(ib, 1, ik) = (dta(ib, 1, ik) + (ztrans - ztrans_new) * r1_hv_n(zij, zjj)) * vmask(zij, zjj, ik)
            ELSE
              dta(ib, 1, ik) = (dta(ib, 1, ik) + (0._wp - ztrans_new) * r1_hv_n(zij, zjj)) * vmask(zij, zjj, ik)
            END IF
          END DO
        END DO
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE fld_bdy_interp
  SUBROUTINE fld_rot(kt, sd)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    TYPE(FLD), DIMENSION(:), INTENT(INOUT) :: sd
    INTEGER :: ju, jv, jk, jn
    INTEGER :: imf
    INTEGER :: ill
    INTEGER :: iv
    CHARACTER(LEN = 100) :: clcomp
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: utmp, vtmp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('fld_rot', 'r0', 0, 0)
    imf = SIZE(sd)
    DO ju = 1, imf
      ill = LEN_TRIM(sd(ju) % vcomp)
      DO jn = 2 - COUNT((/sd(ju) % ln_tint/)), 2
        IF (ill > 0 .AND. .NOT. sd(ju) % rotn(jn)) THEN
          IF (sd(ju) % vcomp(1 : 1) == 'U') THEN
            clcomp = 'V' // sd(ju) % vcomp(2 : ill)
            iv = - 1
            DO jv = 1, imf
              IF (TRIM(sd(jv) % vcomp) == TRIM(clcomp)) iv = jv
            END DO
            IF (iv > 0) THEN
              DO jk = 1, SIZE(sd(ju) % fnow, 3)
                IF (sd(ju) % ln_tint) THEN
                  CALL rot_rep(sd(ju) % fdta(:, :, jk, jn), sd(iv) % fdta(:, :, jk, jn), 'T', 'en->i', utmp(:, :))
                  CALL rot_rep(sd(ju) % fdta(:, :, jk, jn), sd(iv) % fdta(:, :, jk, jn), 'T', 'en->j', vtmp(:, :))
                  sd(ju) % fdta(:, :, jk, jn) = utmp(:, :)
                  sd(iv) % fdta(:, :, jk, jn) = vtmp(:, :)
                ELSE
                  CALL rot_rep(sd(ju) % fnow(:, :, jk), sd(iv) % fnow(:, :, jk), 'T', 'en->i', utmp(:, :))
                  CALL rot_rep(sd(ju) % fnow(:, :, jk), sd(iv) % fnow(:, :, jk), 'T', 'en->j', vtmp(:, :))
                  sd(ju) % fnow(:, :, jk) = utmp(:, :)
                  sd(iv) % fnow(:, :, jk) = vtmp(:, :)
                END IF
              END DO
              sd(ju) % rotn(jn) = .TRUE.
              IF (lwp .AND. kt == nit000) WRITE(numout, *) 'fld_read: vector pair (' // TRIM(sd(ju) % clvar) // ', ' // TRIM(sd(iv) % clvar) // ') rotated on to model grid'
            END IF
          END IF
        END IF
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE fld_rot
  SUBROUTINE fld_clopn(sdjf, kyear, kmonth, kday, ldstop)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(FLD), INTENT(INOUT) :: sdjf
    INTEGER, OPTIONAL, INTENT(IN) :: kyear
    INTEGER, OPTIONAL, INTENT(IN) :: kmonth
    INTEGER, OPTIONAL, INTENT(IN) :: kday
    LOGICAL, OPTIONAL, INTENT(IN) :: ldstop
    LOGICAL :: llprevyr
    LOGICAL :: llprevmth
    INTEGER :: iyear, imonth, iday
    INTEGER :: isec_week
    INTEGER :: indexyr
    INTEGER :: iyear_len, imonth_len
    CHARACTER(LEN = 256) :: clname
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('fld_clopn', 'r0', 0, 0)
    IF (PRESENT(kyear)) THEN
      iyear = kyear
      imonth = kmonth
      iday = kday
      IF (sdjf % cltype(1 : 4) == 'week') THEN
        isec_week = ksec_week(sdjf % cltype(6 : 8)) - (86400 * 8)
        llprevmth = isec_week > nsec_month
        llprevyr = llprevmth .AND. nmonth == 1
        iyear = nyear - COUNT((/llprevyr/))
        imonth = nmonth - COUNT((/llprevmth/)) + 12 * COUNT((/llprevyr/))
        iday = nday + nmonth_len(nmonth - 1) * COUNT((/llprevmth/)) - isec_week / NINT(rday)
      END IF
    ELSE
      IF (sdjf % cltype(1 : 4) == 'week') THEN
        isec_week = ksec_week(sdjf % cltype(6 : 8))
        llprevmth = isec_week > nsec_month
        llprevyr = llprevmth .AND. nmonth == 1
      ELSE
        isec_week = 0
        llprevmth = .FALSE.
        llprevyr = .FALSE.
      END IF
      iyear = nyear - COUNT((/llprevyr/))
      imonth = nmonth - COUNT((/llprevmth/)) + 12 * COUNT((/llprevyr/))
      iday = nday + nmonth_len(nmonth - 1) * COUNT((/llprevmth/)) - isec_week / NINT(rday)
    END IF
    clname = TRIM(sdjf % clrootname)
    IF (.NOT. sdjf % ln_clim) THEN
      WRITE(clname, '(a,"_y",i4.4)') TRIM(sdjf % clrootname), iyear
      IF (sdjf % cltype /= 'yearly') WRITE(clname, '(a,"m" ,i2.2)') TRIM(clname), imonth
    ELSE
      IF (sdjf % cltype /= 'yearly') WRITE(clname, '(a,"_m",i2.2)') TRIM(sdjf % clrootname), imonth
    END IF
    IF (sdjf % cltype == 'daily' .OR. sdjf % cltype(1 : 4) == 'week') WRITE(clname, '(a,"d" ,i2.2)') TRIM(clname), iday
    IF (TRIM(clname) /= TRIM(sdjf % clname) .OR. sdjf % num == 0) THEN
      sdjf % clname = TRIM(clname)
      IF (sdjf % num /= 0) CALL iom_close(sdjf % num)
      CALL iom_open(sdjf % clname, sdjf % num, ldstop = ldstop, ldiof = LEN(TRIM(sdjf % wgtname)) > 0)
      indexyr = iyear - nyear + 1
      iyear_len = nyear_len(indexyr)
      SELECT CASE (indexyr)
      CASE (0)
        imonth_len = 31
      CASE (1)
        imonth_len = nmonth_len(imonth)
      CASE (2)
        imonth_len = 31
      END SELECT
      IF (sdjf % nfreqh == - 12) THEN
        sdjf % nreclast = 1
      ELSE IF (sdjf % nfreqh == - 1) THEN
        IF (sdjf % cltype == 'monthly') THEN
          sdjf % nreclast = 1
        ELSE
          sdjf % nreclast = 12
        END IF
      ELSE
        IF (sdjf % cltype == 'monthly') THEN
          sdjf % nreclast = NINT(24 * imonth_len / sdjf % nfreqh)
        ELSE IF (sdjf % cltype(1 : 4) == 'week') THEN
          sdjf % nreclast = NINT(24 * 7 / sdjf % nfreqh)
        ELSE IF (sdjf % cltype == 'daily') THEN
          sdjf % nreclast = NINT(24 / sdjf % nfreqh)
        ELSE
          sdjf % nreclast = NINT(24 * iyear_len / sdjf % nfreqh)
        END IF
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE fld_clopn
  SUBROUTINE fld_fill(sdf, sdf_n, cdir, cdcaller, cdtitle, cdnam, knoprint)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(FLD), DIMENSION(:), INTENT(INOUT) :: sdf
    TYPE(FLD_N), DIMENSION(:), INTENT(IN) :: sdf_n
    CHARACTER(LEN = *), INTENT(IN) :: cdir
    CHARACTER(LEN = *), INTENT(IN) :: cdcaller
    CHARACTER(LEN = *), INTENT(IN) :: cdtitle
    CHARACTER(LEN = *), INTENT(IN) :: cdnam
    INTEGER, OPTIONAL, INTENT(IN) :: knoprint
    INTEGER :: jf
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('fld_fill', 'r0', 0, 0)
    DO jf = 1, SIZE(sdf)
      sdf(jf) % clrootname = TRIM(cdir) // TRIM(sdf_n(jf) % clname)
      sdf(jf) % clname = "not yet defined"
      sdf(jf) % nfreqh = sdf_n(jf) % nfreqh
      sdf(jf) % clvar = sdf_n(jf) % clvar
      sdf(jf) % ln_tint = sdf_n(jf) % ln_tint
      sdf(jf) % ln_clim = sdf_n(jf) % ln_clim
      sdf(jf) % cltype = sdf_n(jf) % cltype
      sdf(jf) % num = - 1
      sdf(jf) % wgtname = " "
      IF (LEN(TRIM(sdf_n(jf) % wname)) > 0) sdf(jf) % wgtname = TRIM(cdir) // TRIM(sdf_n(jf) % wname)
      sdf(jf) % lsmname = " "
      IF (LEN(TRIM(sdf_n(jf) % lname)) > 0) sdf(jf) % lsmname = TRIM(cdir) // TRIM(sdf_n(jf) % lname)
      sdf(jf) % vcomp = sdf_n(jf) % vcomp
      sdf(jf) % rotn(:) = .TRUE.
      IF (sdf(jf) % cltype(1 : 4) == 'week' .AND. nn_leapy == 0) CALL ctl_stop('fld_clopn: weekly file (' // TRIM(sdf(jf) % clrootname) // ') needs nn_leapy = 1')
      IF (sdf(jf) % cltype(1 : 4) == 'week' .AND. sdf(jf) % ln_clim) CALL ctl_stop('fld_clopn: weekly file (' // TRIM(sdf(jf) % clrootname) // ') needs ln_clim = .FALSE.')
      sdf(jf) % nreclast = - 1
    END DO
    IF (lwp) THEN
      WRITE(numout, *)
      IF (.NOT. PRESENT(knoprint)) THEN
        WRITE(numout, *) TRIM(cdcaller) // ' : ' // TRIM(cdtitle)
        WRITE(numout, *) (/('~', jf = 1, LEN_TRIM(cdcaller))/)
      END IF
      WRITE(numout, *) '   fld_fill : fill data structure with information from namelist ' // TRIM(cdnam)
      WRITE(numout, *) '   ~~~~~~~~'
      WRITE(numout, *) '      list of files and frequency (>0: in hours ; <0 in months)'
      DO jf = 1, SIZE(sdf)
        WRITE(numout, *) '      root filename: ', TRIM(sdf(jf) % clrootname), '   variable name: ', TRIM(sdf(jf) % clvar)
        WRITE(numout, *) '         frequency: ', sdf(jf) % nfreqh, '   time interp: ', sdf(jf) % ln_tint, '   climatology: ', sdf(jf) % ln_clim
        WRITE(numout, *) '         weights: ', TRIM(sdf(jf) % wgtname), '   pairing: ', TRIM(sdf(jf) % vcomp), '   data type: ', sdf(jf) % cltype, '   land/sea mask:', TRIM(sdf(jf) % lsmname)
        CALL flush(numout)
      END DO
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE fld_fill
  SUBROUTINE wgt_list(sd, kwgt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(FLD), INTENT(IN) :: sd
    INTEGER, INTENT(INOUT) :: kwgt
    INTEGER :: kw, nestid
    LOGICAL :: found
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wgt_list', 'r0', 0, 0)
    found = .FALSE.
    nestid = 0
    DO kw = 1, nxt_wgt - 1
      IF (TRIM(ref_wgts(kw) % wgtname) == TRIM(sd % wgtname) .AND. ref_wgts(kw) % nestid == nestid) THEN
        kwgt = kw
        found = .TRUE.
        EXIT
      END IF
    END DO
    IF (.NOT. found) THEN
      kwgt = nxt_wgt
      CALL fld_weight(sd)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE wgt_list
  SUBROUTINE wgt_print
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: kw
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wgt_print', 'r0', 0, 0)
    DO kw = 1, nxt_wgt - 1
      WRITE(numout, *) 'weight file:  ', TRIM(ref_wgts(kw) % wgtname)
      WRITE(numout, *) '      ddims:  ', ref_wgts(kw) % ddims(1), ref_wgts(kw) % ddims(2)
      WRITE(numout, *) '     numwgt:  ', ref_wgts(kw) % numwgt
      WRITE(numout, *) '     jpiwgt:  ', ref_wgts(kw) % jpiwgt
      WRITE(numout, *) '     jpjwgt:  ', ref_wgts(kw) % jpjwgt
      WRITE(numout, *) '    botleft:  ', ref_wgts(kw) % botleft
      WRITE(numout, *) '   topright:  ', ref_wgts(kw) % topright
      IF (ref_wgts(kw) % cyclic) THEN
        WRITE(numout, *) '       cyclical'
        IF (ref_wgts(kw) % overlap > 0) WRITE(numout, *) '              with overlap of ', ref_wgts(kw) % overlap
      ELSE
        WRITE(numout, *) '       not cyclical'
      END IF
      IF (ASSOCIATED(ref_wgts(kw) % data_wgt)) WRITE(numout, *) '       allocated'
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE wgt_print
  SUBROUTINE fld_weight(sd)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(FLD), INTENT(IN) :: sd
    INTEGER :: jn
    INTEGER :: inum
    INTEGER :: id
    INTEGER :: ipk
    INTEGER :: zwrap
    LOGICAL :: cyclical
    CHARACTER(LEN = 5) :: aname
    INTEGER, DIMENSION(:), ALLOCATABLE :: ddims
    INTEGER, DIMENSION(jpi, jpj) :: data_src
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: data_tmp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('fld_weight', 'r0', 0, 0)
    IF (nxt_wgt > tot_wgts) THEN
      CALL ctl_stop("fld_weight: weights array size exceeded, increase tot_wgts")
    END IF
    CALL iom_open(sd % clname, inum, ldiof = LEN(TRIM(sd % wgtname)) > 0)
    IF (SIZE(sd % fnow, 3) > 1) THEN
      ALLOCATE(ddims(4))
    ELSE
      ALLOCATE(ddims(3))
    END IF
    id = iom_varid(inum, sd % clvar, ddims)
    CALL iom_close(inum)
    CALL iom_open(sd % wgtname, inum)
    IF (inum > 0) THEN
      CALL iom_getatt(inum, 'ew_wrap', zwrap)
      IF (zwrap >= 0) THEN
        cyclical = .TRUE.
      ELSE IF (zwrap == - 999) THEN
        cyclical = .TRUE.
        zwrap = 0
      ELSE
        cyclical = .FALSE.
      END IF
      ref_wgts(nxt_wgt) % ddims(1) = ddims(1)
      ref_wgts(nxt_wgt) % ddims(2) = ddims(2)
      ref_wgts(nxt_wgt) % wgtname = sd % wgtname
      ref_wgts(nxt_wgt) % overlap = zwrap
      ref_wgts(nxt_wgt) % cyclic = cyclical
      ref_wgts(nxt_wgt) % nestid = 0
      id = iom_varid(inum, 'src05', ldstop = .FALSE.)
      IF (id <= 0) THEN
        ref_wgts(nxt_wgt) % numwgt = 4
      ELSE
        ref_wgts(nxt_wgt) % numwgt = 16
      END IF
      ALLOCATE(ref_wgts(nxt_wgt) % data_jpi(jpi, jpj, 4))
      ALLOCATE(ref_wgts(nxt_wgt) % data_jpj(jpi, jpj, 4))
      ALLOCATE(ref_wgts(nxt_wgt) % data_wgt(jpi, jpj, ref_wgts(nxt_wgt) % numwgt))
      DO jn = 1, 4
        aname = ' '
        WRITE(aname, '(a3,i2.2)') 'src', jn
        data_tmp(:, :) = 0
        CALL iom_get(inum, jpdom_data, aname, data_tmp(:, :))
        data_src(:, :) = INT(data_tmp(:, :))
        ref_wgts(nxt_wgt) % data_jpj(:, :, jn) = 1 + (data_src(:, :) - 1) / ref_wgts(nxt_wgt) % ddims(1)
        ref_wgts(nxt_wgt) % data_jpi(:, :, jn) = data_src(:, :) - ref_wgts(nxt_wgt) % ddims(1) * (ref_wgts(nxt_wgt) % data_jpj(:, :, jn) - 1)
      END DO
      DO jn = 1, ref_wgts(nxt_wgt) % numwgt
        aname = ' '
        WRITE(aname, '(a3,i2.2)') 'wgt', jn
        ref_wgts(nxt_wgt) % data_wgt(:, :, jn) = 0.0
        CALL iom_get(inum, jpdom_data, aname, ref_wgts(nxt_wgt) % data_wgt(:, :, jn))
      END DO
      CALL iom_close(inum)
      ref_wgts(nxt_wgt) % botleft(1) = MINVAL(ref_wgts(nxt_wgt) % data_jpi(:, :, :))
      ref_wgts(nxt_wgt) % botleft(2) = MINVAL(ref_wgts(nxt_wgt) % data_jpj(:, :, :))
      ref_wgts(nxt_wgt) % topright(1) = MAXVAL(ref_wgts(nxt_wgt) % data_jpi(:, :, :))
      ref_wgts(nxt_wgt) % topright(2) = MAXVAL(ref_wgts(nxt_wgt) % data_jpj(:, :, :))
      ref_wgts(nxt_wgt) % jpiwgt = ref_wgts(nxt_wgt) % topright(1) - ref_wgts(nxt_wgt) % botleft(1) + 1
      ref_wgts(nxt_wgt) % jpjwgt = ref_wgts(nxt_wgt) % topright(2) - ref_wgts(nxt_wgt) % botleft(2) + 1
      ref_wgts(nxt_wgt) % data_jpi(:, :, :) = ref_wgts(nxt_wgt) % data_jpi(:, :, :) - ref_wgts(nxt_wgt) % botleft(1) + 1
      ref_wgts(nxt_wgt) % data_jpj(:, :, :) = ref_wgts(nxt_wgt) % data_jpj(:, :, :) - ref_wgts(nxt_wgt) % botleft(2) + 1
      ipk = SIZE(sd % fnow, 3)
      ALLOCATE(ref_wgts(nxt_wgt) % fly_dta(ref_wgts(nxt_wgt) % jpiwgt + 3, ref_wgts(nxt_wgt) % jpjwgt + 3, ipk))
      IF (ref_wgts(nxt_wgt) % cyclic) ALLOCATE(ref_wgts(nxt_wgt) % col(1, ref_wgts(nxt_wgt) % jpjwgt + 3, ipk))
      nxt_wgt = nxt_wgt + 1
    ELSE
      CALL ctl_stop('    fld_weight : unable to read the file ')
    END IF
    DEALLOCATE(ddims)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE fld_weight
  SUBROUTINE apply_seaoverland(clmaskfile, zfieldo, jpi1_lsm, jpi2_lsm, jpj1_lsm, jpj2_lsm, itmpi, itmpj, itmpz, rec1_lsm, recn_lsm)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: itmpi, itmpj, itmpz
    INTEGER, INTENT(IN) :: jpi1_lsm, jpi2_lsm, jpj1_lsm, jpj2_lsm
    INTEGER, DIMENSION(3), INTENT(IN) :: rec1_lsm, recn_lsm
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: zfieldo
    CHARACTER(LEN = 100), INTENT(IN) :: clmaskfile
    INTEGER :: inum, jni, jnj, jnz, jc
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: zslmec1
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: zfieldn
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: zfield
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('apply_seaoverland', 'r0', 0, 0)
    ALLOCATE(zslmec1(itmpi, itmpj, itmpz), zfieldn(itmpi, itmpj), zfield(itmpi, itmpj))
    CALL iom_open(clmaskfile, inum)
    SELECT CASE (SIZE(zfieldo(jpi1_lsm : jpi2_lsm, jpj1_lsm : jpj2_lsm, :), 3))
    CASE (1)
      CALL iom_get(inum, jpdom_unknown, 'LSM', zslmec1(jpi1_lsm : jpi2_lsm, jpj1_lsm : jpj2_lsm, 1), 1, rec1_lsm, recn_lsm)
    CASE DEFAULT
      CALL iom_get(inum, jpdom_unknown, 'LSM', zslmec1(jpi1_lsm : jpi2_lsm, jpj1_lsm : jpj2_lsm, :), 1, rec1_lsm, recn_lsm)
    END SELECT
    CALL iom_close(inum)
    DO jnz = 1, rec1_lsm(3)
      DO jni = 1, itmpi
        DO jnj = 1, itmpj
          zfieldn(jni, jnj) = zfieldo(jni, jnj, jnz)
          IF (zslmec1(jni, jnj, jnz) == 1.) zfieldn(jni, jnj) = undeff_lsm
        END DO
      END DO
      CALL seaoverland(zfieldn, itmpi, itmpj, zfield)
      DO jc = 1, nn_lsm
        CALL seaoverland(zfield, itmpi, itmpj, zfield)
      END DO
      IF (ANY(zfield == undeff_lsm)) THEN
        DO jni = 1, itmpi
          DO jnj = 1, itmpj
            IF (zfield(jni, jnj) == undeff_lsm) zfield(jni, jnj) = zfieldo(jni, jnj, jnz)
          END DO
        END DO
      END IF
      zfieldo(:, :, jnz) = zfield(:, :)
    END DO
    DEALLOCATE(zslmec1, zfieldn, zfield)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE apply_seaoverland
  SUBROUTINE seaoverland(zfieldn, ileni, ilenj, zfield)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: ileni, ilenj
    REAL, DIMENSION(ileni, ilenj), INTENT(IN) :: zfieldn
    REAL, DIMENSION(ileni, ilenj), INTENT(OUT) :: zfield
    REAL, DIMENSION(ileni, ilenj) :: zmat1, zmat2, zmat3, zmat4
    REAL, DIMENSION(ileni, ilenj) :: zmat5, zmat6, zmat7, zmat8
    REAL, DIMENSION(ileni, ilenj) :: zlsm2d
    REAL, DIMENSION(ileni, ilenj, 8) :: zlsm3d
    LOGICAL, DIMENSION(ileni, ilenj, 8) :: ll_msknan3d
    LOGICAL, DIMENSION(ileni, ilenj) :: ll_msknan2d
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('seaoverland', 'r0', 0, 0)
    zmat8 = EOSHIFT(zfieldn, SHIFT = - 1, BOUNDARY = (/zfieldn(:, 1)/), DIM = 2)
    zmat1 = EOSHIFT(zmat8, SHIFT = - 1, BOUNDARY = (/zmat8(1, :)/), DIM = 1)
    zmat2 = EOSHIFT(zfieldn, SHIFT = - 1, BOUNDARY = (/zfieldn(1, :)/), DIM = 1)
    zmat4 = EOSHIFT(zfieldn, SHIFT = 1, BOUNDARY = (/zfieldn(:, ilenj)/), DIM = 2)
    zmat3 = EOSHIFT(zmat4, SHIFT = - 1, BOUNDARY = (/zmat4(1, :)/), DIM = 1)
    zmat5 = EOSHIFT(zmat4, SHIFT = 1, BOUNDARY = (/zmat4(ileni, :)/), DIM = 1)
    zmat6 = EOSHIFT(zfieldn, SHIFT = 1, BOUNDARY = (/zfieldn(ileni, :)/), DIM = 1)
    zmat7 = EOSHIFT(zmat8, SHIFT = 1, BOUNDARY = (/zmat8(ileni, :)/), DIM = 1)
    zlsm3d = RESHAPE((/zmat1, zmat2, zmat3, zmat4, zmat5, zmat6, zmat7, zmat8/), (/ileni, ilenj, 8/))
    ll_msknan3d = .NOT. (zlsm3d == undeff_lsm)
    ll_msknan2d = .NOT. (zfieldn == undeff_lsm)
    zlsm2d = SUM(zlsm3d, 3, ll_msknan3d) / MAX(1, COUNT(ll_msknan3d, 3))
    WHERE (COUNT(ll_msknan3d, 3) == 0._wp) zlsm2d = undeff_lsm
    zfield = MERGE(zfieldn, zlsm2d, ll_msknan2d)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE seaoverland
  SUBROUTINE fld_interp(num, clvar, kw, kk, dta, nrec, lsmfile)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: num
    CHARACTER(LEN = *), INTENT(IN) :: clvar
    INTEGER, INTENT(IN) :: kw
    INTEGER, INTENT(IN) :: kk
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: dta
    INTEGER, INTENT(IN) :: nrec
    CHARACTER(LEN = *), INTENT(IN) :: lsmfile
    INTEGER, DIMENSION(3) :: rec1, recn
    INTEGER, DIMENSION(3) :: rec1_lsm, recn_lsm
    INTEGER :: ii_lsm1, ii_lsm2, ij_lsm1, ij_lsm2
    INTEGER :: jk, jn, jm, jir, jjr
    INTEGER :: ni, nj
    INTEGER :: jpimin, jpiwid
    INTEGER :: jpimin_lsm, jpiwid_lsm
    INTEGER :: jpjmin, jpjwid
    INTEGER :: jpjmin_lsm, jpjwid_lsm
    INTEGER :: jpi1, jpi2, jpj1, jpj2
    INTEGER :: jpi1_lsm, jpi2_lsm, jpj1_lsm, jpj2_lsm
    INTEGER :: itmpi, itmpj, itmpz
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ztmp_fly_dta
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('fld_interp', 'r0', 0, 0)
    !$ACC KERNELS ! CDe
    jpimin = ref_wgts(kw) % botleft(1)
    jpjmin = ref_wgts(kw) % botleft(2)
    jpiwid = ref_wgts(kw) % jpiwgt
    jpjwid = ref_wgts(kw) % jpjwgt
    rec1(1) = MAX(jpimin - 1, 1)
    rec1(2) = MAX(jpjmin - 1, 1)
    rec1(3) = 1
    recn(1) = MIN(jpiwid + 2, ref_wgts(kw) % ddims(1) - rec1(1) + 1)
    recn(2) = MIN(jpjwid + 2, ref_wgts(kw) % ddims(2) - rec1(2) + 1)
    recn(3) = kk
    jpi1 = 2 + rec1(1) - jpimin
    jpj1 = 2 + rec1(2) - jpjmin
    jpi2 = jpi1 + recn(1) - 1
    jpj2 = jpj1 + recn(2) - 1
    !$ACC END KERNELS
    IF (LEN(TRIM(lsmfile)) > 0) THEN
      !$ACC KERNELS ! CDe
      rec1_lsm(1) = MAX(rec1(1) - nn_lsm, 1)
      rec1_lsm(2) = MAX(rec1(2) - nn_lsm, 1)
      rec1_lsm(3) = 1
      recn_lsm(1) = MIN(rec1(1) - rec1_lsm(1) + recn(1) + nn_lsm, ref_wgts(kw) % ddims(1) - rec1_lsm(1))
      recn_lsm(2) = MIN(rec1(2) - rec1_lsm(2) + recn(2) + nn_lsm, ref_wgts(kw) % ddims(2) - rec1_lsm(2))
      recn_lsm(3) = kk
      jpimin_lsm = MAX(rec1_lsm(1) + 1, 1)
      jpjmin_lsm = MAX(rec1_lsm(2) + 1, 1)
      jpiwid_lsm = MIN(recn_lsm(1) - 2, ref_wgts(kw) % ddims(1) - rec1(1) + 1)
      jpjwid_lsm = MIN(recn_lsm(2) - 2, ref_wgts(kw) % ddims(2) - rec1(2) + 1)
      jpi1_lsm = 2 + rec1_lsm(1) - jpimin_lsm
      jpj1_lsm = 2 + rec1_lsm(2) - jpjmin_lsm
      jpi2_lsm = jpi1_lsm + recn_lsm(1) - 1
      jpj2_lsm = jpj1_lsm + recn_lsm(2) - 1
      itmpi = jpi2_lsm - jpi1_lsm + 1
      itmpj = jpj2_lsm - jpj1_lsm + 1
      itmpz = kk
      !$ACC END KERNELS
      ALLOCATE(ztmp_fly_dta(itmpi, itmpj, itmpz))
      !$ACC KERNELS ! CDe
      ztmp_fly_dta(:, :, :) = 0.0
      !$ACC END KERNELS
      SELECT CASE (SIZE(ztmp_fly_dta(jpi1_lsm : jpi2_lsm, jpj1_lsm : jpj2_lsm, :), 3))
      CASE (1)
        CALL iom_get(num, jpdom_unknown, clvar, ztmp_fly_dta(jpi1_lsm : jpi2_lsm, jpj1_lsm : jpj2_lsm, 1), nrec, rec1_lsm, recn_lsm)
      CASE DEFAULT
        CALL iom_get(num, jpdom_unknown, clvar, ztmp_fly_dta(jpi1_lsm : jpi2_lsm, jpj1_lsm : jpj2_lsm, :), nrec, rec1_lsm, recn_lsm)
      END SELECT
      CALL apply_seaoverland(lsmfile, ztmp_fly_dta(jpi1_lsm : jpi2_lsm, jpj1_lsm : jpj2_lsm, :), jpi1_lsm, jpi2_lsm, jpj1_lsm, jpj2_lsm, itmpi, itmpj, itmpz, rec1_lsm, recn_lsm)
      !$ACC KERNELS ! CDe
      ii_lsm1 = (rec1(1) - rec1_lsm(1)) + 1
      ii_lsm2 = (ii_lsm1 + recn(1)) - 1
      ij_lsm1 = (rec1(2) - rec1_lsm(2)) + 1
      ij_lsm2 = (ij_lsm1 + recn(2)) - 1
      ref_wgts(kw) % fly_dta(:, :, :) = 0.0
      ref_wgts(kw) % fly_dta(jpi1 : jpi2, jpj1 : jpj2, :) = ztmp_fly_dta(ii_lsm1 : ii_lsm2, ij_lsm1 : ij_lsm2, :)
      !$ACC END KERNELS
      DEALLOCATE(ztmp_fly_dta)
    ELSE
      !$ACC KERNELS ! CDe      
      ref_wgts(kw) % fly_dta(:, :, :) = 0.0
      !$ACC END KERNELS
      SELECT CASE (SIZE(ref_wgts(kw) % fly_dta(jpi1 : jpi2, jpj1 : jpj2, :), 3))
      CASE (1)
        CALL iom_get(num, jpdom_unknown, clvar, ref_wgts(kw) % fly_dta(jpi1 : jpi2, jpj1 : jpj2, 1), nrec, rec1, recn)
      CASE DEFAULT
        CALL iom_get(num, jpdom_unknown, clvar, ref_wgts(kw) % fly_dta(jpi1 : jpi2, jpj1 : jpj2, :), nrec, rec1, recn)
      END SELECT
    END IF
    !$ACC KERNELS ! CDe
    dta(:, :, :) = 0.0
    DO jk = 1, 4
      DO jn = 1, jpj
        DO jm = 1, jpi
          ni = ref_wgts(kw) % data_jpi(jm, jn, jk)
          nj = ref_wgts(kw) % data_jpj(jm, jn, jk)
          dta(jm, jn, :) = dta(jm, jn, :) + ref_wgts(kw) % data_wgt(jm, jn, jk) * ref_wgts(kw) % fly_dta(ni + 1, nj + 1, :)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    IF (ref_wgts(kw) % numwgt .EQ. 16) THEN
      IF (jpi1 == 2) THEN
        ref_wgts(kw) % fly_dta(jpi1 - 1, :, :) = ref_wgts(kw) % fly_dta(jpi1, :, :)
      END IF
      IF (jpi2 + jpimin - 1 == ref_wgts(kw) % ddims(1) + 1) THEN
        ref_wgts(kw) % fly_dta(jpi2 + 1, :, :) = ref_wgts(kw) % fly_dta(jpi2, :, :)
      END IF
      IF (jpj1 == 2) THEN
        ref_wgts(kw) % fly_dta(:, jpj1 - 1, :) = ref_wgts(kw) % fly_dta(:, jpj1, :)
      END IF
      IF (jpj2 + jpjmin - 1 == ref_wgts(kw) % ddims(2) + 1 .AND. jpj2 .LT. jpjwid + 2) THEN
        ref_wgts(kw) % fly_dta(:, jpj2 + 1, :) = 2.0 * ref_wgts(kw) % fly_dta(:, jpj2, :) - ref_wgts(kw) % fly_dta(:, jpj2 - 1, :)
      END IF
      IF (ref_wgts(kw) % cyclic) THEN
        rec1(2) = MAX(jpjmin - 1, 1)
        recn(1) = 1
        recn(2) = MIN(jpjwid + 2, ref_wgts(kw) % ddims(2) - rec1(2) + 1)
        jpj1 = 2 + rec1(2) - jpjmin
        jpj2 = jpj1 + recn(2) - 1
        IF (jpi1 == 2) THEN
          rec1(1) = ref_wgts(kw) % ddims(1) - ref_wgts(kw) % overlap
          SELECT CASE (SIZE(ref_wgts(kw) % col(:, jpj1 : jpj2, :), 3))
          CASE (1)
            CALL iom_get(num, jpdom_unknown, clvar, ref_wgts(kw) % col(:, jpj1 : jpj2, 1), nrec, rec1, recn)
          CASE DEFAULT
            CALL iom_get(num, jpdom_unknown, clvar, ref_wgts(kw) % col(:, jpj1 : jpj2, :), nrec, rec1, recn)
          END SELECT
          ref_wgts(kw) % fly_dta(jpi1 - 1, jpj1 : jpj2, :) = ref_wgts(kw) % col(1, jpj1 : jpj2, :)
        END IF
        IF (jpi2 + jpimin - 1 == ref_wgts(kw) % ddims(1) + 1) THEN
          rec1(1) = 1 + ref_wgts(kw) % overlap
          SELECT CASE (SIZE(ref_wgts(kw) % col(:, jpj1 : jpj2, :), 3))
          CASE (1)
            CALL iom_get(num, jpdom_unknown, clvar, ref_wgts(kw) % col(:, jpj1 : jpj2, 1), nrec, rec1, recn)
          CASE DEFAULT
            CALL iom_get(num, jpdom_unknown, clvar, ref_wgts(kw) % col(:, jpj1 : jpj2, :), nrec, rec1, recn)
          END SELECT
          ref_wgts(kw) % fly_dta(jpi2 + 1, jpj1 : jpj2, :) = ref_wgts(kw) % col(1, jpj1 : jpj2, :)
        END IF
      END IF
      !$ACC KERNELS ! CDe
      DO jk = 1, 4
        DO jn = 1, jpj
          DO jm = 1, jpi
            ni = ref_wgts(kw) % data_jpi(jm, jn, jk)
            nj = ref_wgts(kw) % data_jpj(jm, jn, jk)
            dta(jm, jn, :) = dta(jm, jn, :) + ref_wgts(kw) % data_wgt(jm, jn, jk + 4) * 0.5 * (ref_wgts(kw) % fly_dta(ni + 2, nj + 1, :) - ref_wgts(kw) % fly_dta(ni, nj + 1, :))
          END DO
        END DO
      END DO
      DO jk = 1, 4
        DO jn = 1, jpj
          DO jm = 1, jpi
            ni = ref_wgts(kw) % data_jpi(jm, jn, jk)
            nj = ref_wgts(kw) % data_jpj(jm, jn, jk)
            dta(jm, jn, :) = dta(jm, jn, :) + ref_wgts(kw) % data_wgt(jm, jn, jk + 8) * 0.5 * (ref_wgts(kw) % fly_dta(ni + 1, nj + 2, :) - ref_wgts(kw) % fly_dta(ni + 1, nj, :))
          END DO
        END DO
      END DO
      DO jk = 1, 4
        DO jn = 1, jpj
          DO jm = 1, jpi
            ni = ref_wgts(kw) % data_jpi(jm, jn, jk)
            nj = ref_wgts(kw) % data_jpj(jm, jn, jk)
            dta(jm, jn, :) = dta(jm, jn, :) + ref_wgts(kw) % data_wgt(jm, jn, jk + 12) * 0.25 * ((ref_wgts(kw) % fly_dta(ni + 2, nj + 2, :) - ref_wgts(kw) % fly_dta(ni, nj + 2, :)) - (ref_wgts(kw) % fly_dta(ni + 2, nj, :) - ref_wgts(kw) % fly_dta(ni, nj, :)))
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE fld_interp
  FUNCTION ksec_week(cdday)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdday
    INTEGER :: ksec_week
    INTEGER :: ijul, ishift
    CHARACTER(LEN = 3), DIMENSION(7) :: cl_week
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ksec_week', 'r0', 0, 0)
    cl_week = (/"sun", "sat", "fri", "thu", "wed", "tue", "mon"/)
    DO ijul = 1, 7
      IF (cl_week(ijul) == TRIM(cdday)) EXIT
    END DO
    IF (ijul .GT. 7) CALL ctl_stop('ksec_week: wrong day for sdjf%cltype(6:8): ' // TRIM(cdday))
    ishift = ijul * NINT(rday)
    ksec_week = nsec_week + ishift
    ksec_week = MOD(ksec_week, 7 * NINT(rday))
    CALL profile_psy_data0 % PostEnd
  END FUNCTION ksec_week
END MODULE fldread
