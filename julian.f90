MODULE julian
  USE par_kind, ONLY: wp, dp
  USE lib_mpp, ONLY: ctl_warn, ctl_stop
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: jul2greg, greg2jul
  CONTAINS
  RECURSIVE SUBROUTINE jul2greg(ksec, kminut, khour, kday, kmonth, kyear, prelday, krefdate)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN), OPTIONAL :: krefdate
    INTEGER, INTENT(OUT) :: ksec, kminut, khour, kday, kmonth, kyear
    REAL(KIND = dp), INTENT(IN) :: prelday
    INTEGER, PARAMETER :: jpgreg = 2299161, jporef = 2433283, jparef = 2415021
    INTEGER :: ijulian, ij1, ija, ijb, ijc, ijd, ije, isec, imin, ihou, iday, imon, iyea, iref
    REAL(KIND = wp) :: zday, zref
    CHARACTER(LEN = 200) :: cerr
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('jul2greg', 'r0', 0, 0)
    IF (PRESENT(krefdate)) THEN
      SELECT CASE (krefdate)
      CASE (0)
        iref = jpgreg
      CASE (19500101)
        iref = jporef
      CASE (19000101)
        iref = jparef
      CASE DEFAULT
        WRITE(cerr, FMT = '(A,I8.8)') 'jul2greg: Unknown krefdate:', krefdate
        CALL ctl_stop(cerr)
      END SELECT
    ELSE
      iref = jporef
    END IF
    zday = prelday
    ksec = FLOOR(86400. * MOD(zday, 1.))
    IF (ksec < 0.) ksec = 86400. + ksec
    khour = ksec / 3600
    kminut = (ksec - 3600 * khour) / 60
    ksec = MOD(ksec, 60)
    ijulian = iref + INT(zday)
    IF (zday < 0.) ijulian = ijulian - 1
    IF (ijulian >= jpgreg) THEN
      ij1 = INT((DBLE(ijulian - 1867216) - 0.25) / 36524.25)
      ija = ijulian + 1 + ij1 - INT((0.25 * ij1))
    ELSE
      ija = ijulian
    END IF
    ijb = ija + 1524
    ijc = INT(6680. + (DBLE(ijb - 2439870) - 122.1) / 365.25)
    ijd = 365 * ijc + INT(0.25 * ijc)
    ije = INT((ijb - ijd) / 30.6001)
    kday = ijb - ijd - INT(30.6001 * ije)
    kmonth = ije - 1
    IF (kmonth > 12) kmonth = kmonth - 12
    kyear = ijc - 4715
    IF (kmonth > 2) kyear = kyear - 1
    IF (kyear <= 0) kyear = kyear - 1
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE jul2greg
  SUBROUTINE greg2jul(ksec, kmin, khour, kday, kmonth, kyear, pjulian, krefdate)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: ksec, kmin, khour, kday, kmonth, kyear
    REAL(KIND = dp), INTENT(OUT) :: pjulian
    INTEGER, INTENT(IN), OPTIONAL :: krefdate
    INTEGER, PARAMETER :: jpgreg = 15 + 31 * (10 + 12 * 1582), jporef = 2433283, jparef = 2415021, jpgref = 2299161
    INTEGER :: ija, ijy, ijm, ijultmp, ijyear, iref
    CHARACTER(LEN = 200) :: cerr
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('greg2jul', 'r0', 0, 0)
    IF (PRESENT(krefdate)) THEN
      SELECT CASE (krefdate)
      CASE (0)
        iref = jpgref
      CASE (19500101)
        iref = jporef
      CASE (19000101)
        iref = jparef
      CASE DEFAULT
        WRITE(cerr, FMT = '(A,I8.8)') 'greg2jul: Unknown krefdate:', krefdate
        CALL ctl_stop(cerr)
      END SELECT
    ELSE
      iref = jporef
    END IF
    ijyear = kyear
    IF (ijyear < 0) ijyear = ijyear + 1
    IF (kmonth > 2) THEN
      ijy = ijyear
      ijm = kmonth + 1
    ELSE
      ijy = ijyear - 1
      ijm = kmonth + 13
    END IF
    ijultmp = INT(365.25 * ijy) + INT(30.6001 * ijm) + kday + 1720995
    IF (kday + 31 * (kmonth + 12 * ijyear) >= jpgreg) THEN
      ija = INT(0.01 * ijy)
      ijultmp = ijultmp + 2 - ija + INT(0.25 * ija)
    END IF
    pjulian = (ijultmp - iref) + ((60 * khour + kmin) * 60 + ksec) / 86400.
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE greg2jul
END MODULE julian