MODULE dianam
  USE dom_oce
  USE phycst
  USE in_out_manager
  USE ioipsl, ONLY: ju2ymds
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dia_nam
  CONTAINS
  SUBROUTINE dia_nam(cdfnam, kfreq, cdsuff, ldfsec)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(OUT) :: cdfnam
    CHARACTER(LEN = *), INTENT(IN) :: cdsuff
    INTEGER, INTENT(IN) :: kfreq
    LOGICAL, INTENT(IN), OPTIONAL :: ldfsec
    CHARACTER(LEN = 20) :: clfmt, clfmt0
    CHARACTER(LEN = 20) :: clave
    CHARACTER(LEN = 20) :: cldate1
    CHARACTER(LEN = 20) :: cldate2
    LOGICAL :: llfsec
    INTEGER :: iyear1, imonth1, iday1
    INTEGER :: iyear2, imonth2, iday2
    INTEGER :: indg
    INTEGER :: inbsec, inbmn, inbhr
    INTEGER :: inbday, inbmo, inbyr
    INTEGER :: iyyss, iddss, ihhss, immss
    INTEGER :: iyymo
    REAL(KIND = wp) :: zsec1, zsec2
    REAL(KIND = wp) :: zdrun, zjul
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('dia_nam', 'r0', 0, 0)
    IF (PRESENT(ldfsec)) THEN
      llfsec = ldfsec
    ELSE
      llfsec = .FALSE.
    END IF
    IF (llfsec .OR. kfreq < 0) THEN
      inbsec = kfreq
    ELSE
      inbsec = kfreq * NINT(rdt)
    END IF
    iddss = NINT(rday)
    ihhss = NINT(rmmss * rhhmm)
    immss = NINT(rmmss)
    iyymo = NINT(raamo)
    iyyss = iddss * nyear_len(1)
    clfmt0 = "('(a,i',i1,',a)')"
    IF (inbsec == 0) THEN
      clave = ''
    ELSE IF (inbsec < 0) THEN
      inbmo = - inbsec
      IF (MOD(inbmo, iyymo) == 0) THEN
        inbyr = inbmo / iyymo
        indg = INT(LOG10(REAL(inbyr, wp))) + 1
        WRITE(clfmt, clfmt0) indg
        WRITE(clave, clfmt) '_', inbyr, 'y'
      ELSE
        indg = INT(LOG10(REAL(inbmo, wp))) + 1
        WRITE(clfmt, clfmt0) indg
        WRITE(clave, clfmt) '_', inbmo, 'm'
      END IF
    ELSE IF (MOD(inbsec, iyyss) == 0) THEN
      inbyr = inbsec / iyyss
      indg = INT(LOG10(REAL(inbyr, wp))) + 1
      WRITE(clfmt, clfmt0) indg
      WRITE(clave, clfmt) '_', inbyr, 'y'
    ELSE IF (MOD(inbsec, iddss) == 0) THEN
      inbday = inbsec / iddss
      indg = INT(LOG10(REAL(inbday, wp))) + 1
      WRITE(clfmt, clfmt0) indg
      WRITE(clave, clfmt) '_', inbday, 'd'
      IF (inbday == nmonth_len(nmonth)) clave = '_1m'
    ELSE IF (MOD(inbsec, ihhss) == 0) THEN
      inbhr = inbsec / ihhss
      indg = INT(LOG10(REAL(inbhr, wp))) + 1
      WRITE(clfmt, clfmt0) indg
      WRITE(clave, clfmt) '_', inbhr, 'h'
    ELSE IF (MOD(inbsec, immss) == 0) THEN
      inbmn = inbsec / immss
      indg = INT(LOG10(REAL(inbmn, wp))) + 1
      WRITE(clfmt, clfmt0) indg
      WRITE(clave, clfmt) '_', inbmn, 'mn'
    ELSE
      indg = INT(LOG10(REAL(inbsec, wp))) + 1
      WRITE(clfmt, clfmt0) indg
      WRITE(clave, clfmt) '_', inbsec, 's'
    END IF
    zdrun = rdt / rday * REAL(nitend - nit000, wp)
    zjul = fjulday - rdt / rday
    CALL ju2ymds(zjul, iyear1, imonth1, iday1, zsec1)
    CALL ju2ymds(zjul + zdrun, iyear2, imonth2, iday2, zsec2)
    IF (iyear2 < 10000) THEN
      clfmt = "(i4.4,2i2.2)"
    ELSE
      WRITE(clfmt, FMT = "('(i',i1,',2i2.2)')") INT(LOG10(REAL(iyear2, wp))) + 1
    END IF
    WRITE(cldate1, clfmt) iyear1, imonth1, iday1
    WRITE(cldate2, clfmt) iyear2, imonth2, iday2
    cdfnam = TRIM(cexper) // TRIM(clave) // "_" // TRIM(cldate1) // "_" // TRIM(cldate2) // "_" // TRIM(cdsuff)
    IF (.NOT. Agrif_Root()) cdfnam = TRIM(Agrif_CFixed()) // '_' // TRIM(cdfnam)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE dia_nam
END MODULE dianam