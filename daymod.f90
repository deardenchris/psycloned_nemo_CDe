MODULE daymod
  USE dom_oce
  USE phycst
  USE ioipsl, ONLY: ymds2ju
  USE trc_oce, ONLY: l_offline
  USE in_out_manager
  USE prtctl
  USE iom
  USE timing
  USE restart
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: day
  PUBLIC :: day_init
  PUBLIC :: day_mth
  INTEGER, PUBLIC :: nsecd, nsecd05, ndt, ndt05
  CONTAINS
  SUBROUTINE day_init
    INTEGER :: inbday, idweek
    REAL(KIND = wp) :: zjul
    IF (REAL(nitend - nit000 + 1) * rdt > REAL(HUGE(nsec1jan000))) THEN
      CALL ctl_stop('The number of seconds between each restart exceeds the integer 4 max value: 2^31-1. ', &
&'You must do a restart at higher frequency (or remove this stop and recompile the code in I8)')
    END IF
    nsecd = NINT(rday)
    nsecd05 = NINT(0.5 * rday)
    ndt = NINT(rdt)
    ndt05 = NINT(0.5 * rdt)
    IF (.NOT. l_offline) CALL day_rst(nit000, 'READ')
    nyear = ndastp / 10000
    nmonth = (ndastp - (nyear * 10000)) / 100
    nday = ndastp - (nyear * 10000) - (nmonth * 100)
    nhour = nn_time0 / 100
    nminute = (nn_time0 - nhour * 100)
    CALL ymds2ju(nyear, nmonth, nday, nhour * 3600._wp + nminute * 60._wp, fjulday)
    IF (ABS(fjulday - REAL(NINT(fjulday), wp)) < 0.1 / rday) fjulday = REAL(NINT(fjulday), wp)
    IF (nn_time0 * 3600 - ndt05 .LT. 0) fjulday = fjulday + 1.
    nsec1jan000 = 0
    CALL day_mth
    IF (nday == 0) THEN
      nmonth = nmonth - 1
      nday = nmonth_len(nmonth)
    END IF
    IF (nmonth == 0) THEN
      nmonth = 12
      nyear = nyear - 1
      nsec1jan000 = nsec1jan000 - nsecd * nyear_len(0)
      IF (nleapy == 1) CALL day_mth
    END IF
    nday_year = nday + SUM(nmonth_len(1 : nmonth - 1))
    CALL ymds2ju(1900, 01, 01, 0.0, zjul)
    inbday = FLOOR(fjulday - zjul)
    idweek = MOD(inbday, 7)
    IF (idweek .LT. 0) idweek = idweek + 7
    IF (nhour * 3600 + nminute * 60 - ndt05 .GT. 0) THEN
      nsec_year = (nday_year - 1) * nsecd + nhour * 3600 + nminute * 60 - ndt05
      nsec_month = (nday - 1) * nsecd + nhour * 3600 + nminute * 60 - ndt05
    ELSE
      nsec_year = nday_year * nsecd + nhour * 3600 + nminute * 60 - ndt05
      nsec_month = nday * nsecd + nhour * 3600 + nminute * 60 - ndt05
    END IF
    nsec_week = idweek * nsecd + nhour * 3600 + nminute * 60 - ndt05
    nsec_day = nhour * 3600 + nminute * 60 - ndt05
    IF (nsec_day .LT. 0) nsec_day = nsec_day + nsecd
    IF (nsec_week .LT. 0) nsec_week = nsec_week + nsecd * 7
    IF (lwp) WRITE(numout, FMT = '(a,i6,a,i2,a,i2,a,i8,a,i8,a,i8,a,i8)') ' =======>> 1/2 time step before the start of the run &
&DATE Y/M/D = ', nyear, '/', nmonth, '/', nday, '  nsec_day:', nsec_day, '  nsec_week:', nsec_week, '                     &
&nsec_month:', nsec_month, '  nsec_year:', nsec_year
    CALL day(nit000)
    IF (lwxios) THEN
      CALL iom_set_rstw_var_active('kt')
      CALL iom_set_rstw_var_active('ndastp')
      CALL iom_set_rstw_var_active('adatrj')
      CALL iom_set_rstw_var_active('ntime')
    END IF
  END SUBROUTINE day_init
  SUBROUTINE day_mth
    INTEGER :: jm
    IF (nleapy < 2) THEN
      nmonth_len(:) = (/31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31/)
      nyear_len(:) = 365
      IF (nleapy == 1) THEN
        IF (MOD(nyear - 1, 4) == 0 .AND. (MOD(nyear - 1, 400) == 0 .OR. MOD(nyear - 1, 100) /= 0)) THEN
          nyear_len(0) = 366
        END IF
        IF (MOD(nyear, 4) == 0 .AND. (MOD(nyear, 400) == 0 .OR. MOD(nyear, 100) /= 0)) THEN
          nmonth_len(2) = 29
          nyear_len(1) = 366
        END IF
        IF (MOD(nyear + 1, 4) == 0 .AND. (MOD(nyear + 1, 400) == 0 .OR. MOD(nyear + 1, 100) /= 0)) THEN
          nyear_len(2) = 366
        END IF
      END IF
    ELSE
      nmonth_len(:) = nleapy
      nyear_len(:) = 12 * nleapy
    END IF
    nmonth_half(0) = - nsecd05 * nmonth_len(0)
    DO jm = 1, 13
      nmonth_half(jm) = nmonth_half(jm - 1) + nsecd05 * (nmonth_len(jm - 1) + nmonth_len(jm))
    END DO
    nmonth_end(0) = 0
    DO jm = 1, 13
      nmonth_end(jm) = nmonth_end(jm - 1) + nsecd * nmonth_len(jm)
    END DO
  END SUBROUTINE
  SUBROUTINE day(kt)
    INTEGER, INTENT(IN) :: kt
    CHARACTER(LEN = 25) :: charout
    REAL(KIND = wp) :: zprec
    IF (ln_timing) CALL timing_start('day')
    zprec = 0.1 / rday
    nsec_year = nsec_year + ndt
    nsec_month = nsec_month + ndt
    nsec_week = nsec_week + ndt
    nsec_day = nsec_day + ndt
    adatrj = adatrj + rdt / rday
    fjulday = fjulday + rdt / rday
    IF (ABS(fjulday - REAL(NINT(fjulday), wp)) < zprec) fjulday = REAL(NINT(fjulday), wp)
    IF (ABS(adatrj - REAL(NINT(adatrj), wp)) < zprec) adatrj = REAL(NINT(adatrj), wp)
    IF (nsec_day > nsecd) THEN
      nday = nday + 1
      nday_year = nday_year + 1
      nsec_day = ndt05
      IF (nday == nmonth_len(nmonth) + 1) THEN
        nday = 1
        nmonth = nmonth + 1
        nsec_month = ndt05
        IF (nmonth == 13) THEN
          nyear = nyear + 1
          nmonth = 1
          nday_year = 1
          nsec_year = ndt05
          nsec1jan000 = nsec1jan000 + nsecd * nyear_len(1)
          IF (nleapy == 1) CALL day_mth
        END IF
      END IF
      ndastp = nyear * 10000 + nmonth * 100 + nday
      CALL ymds2ju(nyear, 01, 01, 0.0, fjulstartyear)
      IF (lwp) WRITE(numout, FMT = '(a,i8,a,i4.4,a,i2.2,a,i2.2,a,i3.3)') '======>> time-step =', kt, '      New day, DATE Y/M/D = &
&', nyear, '/', nmonth, '/', nday, '      nday_year = ', nday_year
      IF (lwp) WRITE(numout, FMT = '(a,i8,a,i7,a,i5)') '         nsec_year = ', nsec_year, '   nsec_month = ', nsec_month, '   &
&nsec_day = ', nsec_day, '   nsec_week = ', nsec_week
    END IF
    IF (nsec_week > 7 * nsecd) nsec_week = ndt05
    IF (ln_ctl) THEN
      WRITE(charout, FMT = "('kt =', I4,'  d/m/y =',I2,I2,I4)") kt, nday, nmonth, nyear
      CALL prt_ctl_info(charout)
    END IF
    IF (.NOT. l_offline) CALL rst_opn(kt)
    IF (lrst_oce) CALL day_rst(kt, 'WRITE')
    IF (ln_timing) CALL timing_stop('day')
  END SUBROUTINE day
  SUBROUTINE day_rst(kt, cdrw)
    INTEGER, INTENT(IN) :: kt
    CHARACTER(LEN = *), INTENT(IN) :: cdrw
    REAL(KIND = wp) :: zkt, zndastp, zdayfrac, ksecs, ktime
    INTEGER :: ihour, iminute
    IF (TRIM(cdrw) == 'READ') THEN
      IF (iom_varid(numror, 'kt', ldstop = .FALSE.) > 0) THEN
        CALL iom_get(numror, 'kt', zkt, ldxios = lrxios)
        IF (lwp) THEN
          WRITE(numout, FMT = *) ' *** Info read in restart : '
          WRITE(numout, FMT = *) '   previous time-step                               : ', NINT(zkt)
          WRITE(numout, FMT = *) ' *** restart option'
          SELECT CASE (nrstdt)
          CASE (0)
            WRITE(numout, FMT = *) ' nrstdt = 0 : no control of nit000'
          CASE (1)
            WRITE(numout, FMT = *) ' nrstdt = 1 : no control the date at nit000 (use ndate0 read in the namelist)'
          CASE (2)
            WRITE(numout, FMT = *) ' nrstdt = 2 : calendar parameters read in restart'
          END SELECT
          WRITE(numout, FMT = *)
        END IF
        IF (nit000 - NINT(zkt) /= 1 .AND. nrstdt /= 0) CALL ctl_stop(' ===>>>> : problem with nit000 for the restart', ' verify &
&the restart file or rerun with nrstdt = 0 (namelist)')
        IF (nrstdt == 2) THEN
          CALL iom_get(numror, 'ndastp', zndastp, ldxios = lrxios)
          ndastp = NINT(zndastp)
          CALL iom_get(numror, 'adatrj', adatrj, ldxios = lrxios)
          CALL iom_get(numror, 'ntime', ktime, ldxios = lrxios)
          nn_time0 = INT(ktime)
          zdayfrac = adatrj - INT(adatrj)
          ksecs = NINT(zdayfrac * 86400)
          ihour = INT(ksecs / 3600)
          iminute = ksecs / 60 - ihour * 60
          nhour = nn_time0 / 100
          nminute = (nn_time0 - nhour * 100)
          nminute = nminute + iminute
          IF (nminute >= 60) THEN
            nminute = nminute - 60
            nhour = nhour + 1
          END IF
          nhour = nhour + ihour
          IF (nhour >= 24) THEN
            nhour = nhour - 24
            adatrj = adatrj + 1
          END IF
          nn_time0 = nhour * 100 + nminute
          adatrj = INT(adatrj)
        ELSE
          ndastp = ndate0
          nhour = nn_time0 / 100
          nminute = (nn_time0 - nhour * 100)
          IF (nhour * 3600 + nminute * 60 - ndt05 .LT. 0) ndastp = ndastp - 1
          adatrj = (REAL(nit000 - 1, wp) * rdt) / rday
        END IF
      ELSE
        ndastp = ndate0
        nhour = nn_time0 / 100
        nminute = (nn_time0 - nhour * 100)
        IF (nhour * 3600 + nminute * 60 - ndt05 .LT. 0) ndastp = ndastp - 1
        adatrj = (REAL(nit000 - 1, wp) * rdt) / rday
      END IF
      IF (ABS(adatrj - REAL(NINT(adatrj), wp)) < 0.1 / rday) adatrj = REAL(NINT(adatrj), wp)
      IF (lwp) THEN
        WRITE(numout, FMT = *) ' *** Info used values : '
        WRITE(numout, FMT = *) '   date ndastp                                      : ', ndastp
        WRITE(numout, FMT = *) '   number of elapsed days since the begining of run : ', adatrj
        WRITE(numout, FMT = *) '   nn_time0                                         : ', nn_time0
        WRITE(numout, FMT = *)
      END IF
    ELSE IF (TRIM(cdrw) == 'WRITE') THEN
      IF (kt == nitrst) THEN
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) 'rst_write : write oce restart file  kt =', kt
        IF (lwp) WRITE(numout, FMT = *) '~~~~~~~'
      END IF
      IF (lwxios) CALL iom_swap(cwxios_context)
      CALL iom_rstput(kt, nitrst, numrow, 'kt', REAL(kt, wp), ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'ndastp', REAL(ndastp, wp), ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'adatrj', adatrj, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'ntime', REAL(nn_time0, wp), ldxios = lwxios)
      IF (lwxios) CALL iom_swap(cxios_context)
    END IF
  END SUBROUTINE day_rst
END MODULE daymod