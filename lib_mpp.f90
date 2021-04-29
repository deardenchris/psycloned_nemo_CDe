MODULE lib_mpp
  USE in_out_manager
  INTERFACE mpp_sum
    MODULE PROCEDURE mppsum_int, mppsum_a_int, mppsum_real, mppsum_a_real, mppsum_realdd, mppsum_a_realdd
  END INTERFACE
  INTERFACE mpp_max
    MODULE PROCEDURE mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real
  END INTERFACE
  INTERFACE mpp_min
    MODULE PROCEDURE mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real
  END INTERFACE
  INTERFACE mpp_minloc
    MODULE PROCEDURE mpp_minloc2d, mpp_minloc3d
  END INTERFACE
  INTERFACE mpp_maxloc
    MODULE PROCEDURE mpp_maxloc2d, mpp_maxloc3d
  END INTERFACE
  LOGICAL, PUBLIC, PARAMETER :: lk_mpp = .FALSE.
  LOGICAL, PUBLIC :: ln_nnogather
  INTEGER, PUBLIC :: mpi_comm_oce
  INTEGER, PARAMETER, PUBLIC :: nbdelay = 0
  CHARACTER(LEN = 32), DIMENSION(1), PUBLIC :: c_delaylist = 'empty'
  CHARACTER(LEN = 32), DIMENSION(1), PUBLIC :: c_delaycpnt = 'empty'
  LOGICAL, PUBLIC :: l_full_nf_update = .TRUE.
  TYPE :: DELAYARR
    REAL(KIND = wp), POINTER, DIMENSION(:) :: z1d => NULL()
    COMPLEX(KIND = wp), POINTER, DIMENSION(:) :: y1d => NULL()
  END TYPE DELAYARR
  TYPE(DELAYARR), DIMENSION(1), PUBLIC :: todelay
  INTEGER, PUBLIC, DIMENSION(1) :: ndelayid = - 1
  CONTAINS
  INTEGER FUNCTION lib_mpp_alloc(kumout)
    INTEGER, INTENT(IN) :: kumout
    lib_mpp_alloc = 0
  END FUNCTION lib_mpp_alloc
  FUNCTION mynode(ldtxt, ldname, kumnam_ref, knumnam_cfg, kumond, kstop, localComm) RESULT(function_value)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, OPTIONAL, INTENT(IN) :: localComm
    CHARACTER(LEN = *), DIMENSION(:) :: ldtxt
    CHARACTER(LEN = *) :: ldname
    INTEGER :: kumnam_ref, knumnam_cfg, kumond, kstop
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mynode', 'r0', 0, 0)
    IF (PRESENT(localComm)) mpi_comm_oce = localComm
    function_value = 0
    IF (.FALSE.) ldtxt(:) = 'never done'
    CALL ctl_opn(kumond, TRIM(ldname), 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', - 1, 6, .FALSE., 1)
    CALL profile_psy_data0 % PostEnd
  END FUNCTION mynode
  SUBROUTINE mppsync
  END SUBROUTINE mppsync
  SUBROUTINE mppmax_int(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    INTEGER, INTENT(INOUT) :: ptab
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    RETURN
  END SUBROUTINE mppmax_int
  SUBROUTINE mppmax_a_int(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    INTEGER, INTENT(INOUT) :: ptab(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    RETURN
  END SUBROUTINE mppmax_a_int
  SUBROUTINE mppmax_real(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    RETURN
  END SUBROUTINE mppmax_real
  SUBROUTINE mppmax_a_real(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    RETURN
  END SUBROUTINE mppmax_a_real
  SUBROUTINE mppmin_int(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    INTEGER, INTENT(INOUT) :: ptab
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    RETURN
  END SUBROUTINE mppmin_int
  SUBROUTINE mppmin_a_int(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    INTEGER, INTENT(INOUT) :: ptab(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    RETURN
  END SUBROUTINE mppmin_a_int
  SUBROUTINE mppmin_real(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    RETURN
  END SUBROUTINE mppmin_real
  SUBROUTINE mppmin_a_real(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    RETURN
  END SUBROUTINE mppmin_a_real
  SUBROUTINE mppsum_int(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    INTEGER, INTENT(INOUT) :: ptab
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    RETURN
  END SUBROUTINE mppsum_int
  SUBROUTINE mppsum_a_int(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    INTEGER, INTENT(INOUT) :: ptab(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    RETURN
  END SUBROUTINE mppsum_a_int
  SUBROUTINE mppsum_real(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    RETURN
  END SUBROUTINE mppsum_real
  SUBROUTINE mppsum_a_real(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    RETURN
  END SUBROUTINE mppsum_a_real
  SUBROUTINE mppsum_realdd(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    COMPLEX, INTENT(INOUT) :: ptab
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    RETURN
  END SUBROUTINE mppsum_realdd
  SUBROUTINE mppsum_a_realdd(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    COMPLEX, INTENT(INOUT) :: ptab(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    RETURN
  END SUBROUTINE mppsum_a_realdd
  SUBROUTINE mpp_minloc2d(cdname, ptab, pmask, pmin, kindex)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :)
    REAL(KIND = wp), INTENT(IN) :: pmask(:, :)
    REAL(KIND = wp), INTENT(OUT) :: pmin
    INTEGER, INTENT(OUT) :: kindex(2)
    kindex = 0
    pmin = 0.
    WRITE(*, *) 'ROUTINE_LOC: You should not have seen this print! error?'
  END SUBROUTINE mpp_minloc2d
  SUBROUTINE mpp_minloc3d(cdname, ptab, pmask, pmin, kindex)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :, :)
    REAL(KIND = wp), INTENT(IN) :: pmask(:, :, :)
    REAL(KIND = wp), INTENT(OUT) :: pmin
    INTEGER, INTENT(OUT) :: kindex(3)
    kindex = 0
    pmin = 0.
    WRITE(*, *) 'ROUTINE_LOC: You should not have seen this print! error?'
  END SUBROUTINE mpp_minloc3d
  SUBROUTINE mpp_maxloc2d(cdname, ptab, pmask, pmin, kindex)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :)
    REAL(KIND = wp), INTENT(IN) :: pmask(:, :)
    REAL(KIND = wp), INTENT(OUT) :: pmin
    INTEGER, INTENT(OUT) :: kindex(2)
    kindex = 0
    pmin = 0.
    WRITE(*, *) 'ROUTINE_LOC: You should not have seen this print! error?'
  END SUBROUTINE mpp_maxloc2d
  SUBROUTINE mpp_maxloc3d(cdname, ptab, pmask, pmin, kindex)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :, :)
    REAL(KIND = wp), INTENT(IN) :: pmask(:, :, :)
    REAL(KIND = wp), INTENT(OUT) :: pmin
    INTEGER, INTENT(OUT) :: kindex(3)
    kindex = 0
    pmin = 0.
    WRITE(*, *) 'ROUTINE_LOC: You should not have seen this print! error?'
  END SUBROUTINE mpp_maxloc3d
  SUBROUTINE mpp_delay_sum(cdname, cdelay, y_in, pout, ldlast, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    CHARACTER(LEN = *), INTENT(IN) :: cdelay
    COMPLEX(KIND = wp), INTENT(IN), DIMENSION(:) :: y_in
    REAL(KIND = wp), INTENT(OUT), DIMENSION(:) :: pout
    LOGICAL, INTENT(IN) :: ldlast
    INTEGER, INTENT(IN), OPTIONAL :: kcom
    !$ACC KERNELS
    pout(:) = REAL(y_in(:), wp)
    !$ACC END KERNELS
  END SUBROUTINE mpp_delay_sum
  SUBROUTINE mpp_delay_max(cdname, cdelay, p_in, pout, ldlast, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    CHARACTER(LEN = *), INTENT(IN) :: cdelay
    REAL(KIND = wp), INTENT(IN), DIMENSION(:) :: p_in
    REAL(KIND = wp), INTENT(OUT), DIMENSION(:) :: pout
    LOGICAL, INTENT(IN) :: ldlast
    INTEGER, INTENT(IN), OPTIONAL :: kcom
    !$ACC KERNELS
    pout(:) = p_in(:)
    !$ACC END KERNELS
  END SUBROUTINE mpp_delay_max
  SUBROUTINE mpp_delay_rcv(kid)
    INTEGER, INTENT(IN) :: kid
    WRITE(*, *) 'mpp_delay_rcv: You should not have seen this print! error?', kid
  END SUBROUTINE mpp_delay_rcv
  SUBROUTINE mppstop(ldfinal, ld_force_abort)
    LOGICAL, OPTIONAL, INTENT(IN) :: ldfinal
    LOGICAL, OPTIONAL, INTENT(IN) :: ld_force_abort
    STOP
  END SUBROUTINE mppstop
  SUBROUTINE mpp_ini_znl(knum)
    INTEGER :: knum
    WRITE(*, *) 'mpp_ini_znl: You should not have seen this print! error?', knum
  END SUBROUTINE mpp_ini_znl
  SUBROUTINE mpp_comm_free(kcom)
    INTEGER :: kcom
    WRITE(*, *) 'mpp_comm_free: You should not have seen this print! error?', kcom
  END SUBROUTINE mpp_comm_free
  SUBROUTINE ctl_stop(cd1, cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cd1, cd2, cd3, cd4, cd5
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cd6, cd7, cd8, cd9, cd10
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ctl_stop', 'r0', 0, 0)
    nstop = nstop + 1
    IF (numout == 6) CALL ctl_opn(numout, 'ocean.output', 'APPEND', 'FORMATTED', 'SEQUENTIAL', - 1, 6, .FALSE.)
    WRITE(numout, cform_err)
    IF (PRESENT(cd1)) WRITE(numout, *) TRIM(cd1)
    IF (PRESENT(cd2)) WRITE(numout, *) TRIM(cd2)
    IF (PRESENT(cd3)) WRITE(numout, *) TRIM(cd3)
    IF (PRESENT(cd4)) WRITE(numout, *) TRIM(cd4)
    IF (PRESENT(cd5)) WRITE(numout, *) TRIM(cd5)
    IF (PRESENT(cd6)) WRITE(numout, *) TRIM(cd6)
    IF (PRESENT(cd7)) WRITE(numout, *) TRIM(cd7)
    IF (PRESENT(cd8)) WRITE(numout, *) TRIM(cd8)
    IF (PRESENT(cd9)) WRITE(numout, *) TRIM(cd9)
    IF (PRESENT(cd10)) WRITE(numout, *) TRIM(cd10)
    CALL FLUSH(numout)
    IF (numstp /= - 1) CALL FLUSH(numstp)
    IF (numrun /= - 1) CALL FLUSH(numrun)
    IF (numevo_ice /= - 1) CALL FLUSH(numevo_ice)
    IF (cd1 == 'STOP') THEN
      WRITE(numout, *) 'huge E-R-R-O-R : immediate stop'
      CALL mppstop(ld_force_abort = .TRUE.)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ctl_stop
  SUBROUTINE ctl_warn(cd1, cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cd1, cd2, cd3, cd4, cd5
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cd6, cd7, cd8, cd9, cd10
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ctl_warn', 'r0', 0, 0)
    nwarn = nwarn + 1
    IF (lwp) THEN
      WRITE(numout, cform_war)
      IF (PRESENT(cd1)) WRITE(numout, *) TRIM(cd1)
      IF (PRESENT(cd2)) WRITE(numout, *) TRIM(cd2)
      IF (PRESENT(cd3)) WRITE(numout, *) TRIM(cd3)
      IF (PRESENT(cd4)) WRITE(numout, *) TRIM(cd4)
      IF (PRESENT(cd5)) WRITE(numout, *) TRIM(cd5)
      IF (PRESENT(cd6)) WRITE(numout, *) TRIM(cd6)
      IF (PRESENT(cd7)) WRITE(numout, *) TRIM(cd7)
      IF (PRESENT(cd8)) WRITE(numout, *) TRIM(cd8)
      IF (PRESENT(cd9)) WRITE(numout, *) TRIM(cd9)
      IF (PRESENT(cd10)) WRITE(numout, *) TRIM(cd10)
    END IF
    CALL FLUSH(numout)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ctl_warn
  SUBROUTINE ctl_opn(knum, cdfile, cdstat, cdform, cdacce, klengh, kout, ldwp, karea)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(OUT) :: knum
    CHARACTER(LEN = *), INTENT(IN) :: cdfile
    CHARACTER(LEN = *), INTENT(IN) :: cdstat
    CHARACTER(LEN = *), INTENT(IN) :: cdform
    CHARACTER(LEN = *), INTENT(IN) :: cdacce
    INTEGER, INTENT(IN) :: klengh
    INTEGER, INTENT(IN) :: kout
    LOGICAL, INTENT(IN) :: ldwp
    INTEGER, OPTIONAL, INTENT(IN) :: karea
    CHARACTER(LEN = 80) :: clfile
    INTEGER :: iost
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ctl_opn', 'r0', 0, 0)
    clfile = TRIM(cdfile)
    IF (PRESENT(karea)) THEN
      IF (karea > 1) WRITE(clfile, "(a,'_',i6.6)") TRIM(clfile), karea - 1
    END IF
    knum = get_unit()
    IF (TRIM(cdfile) == '/dev/null') clfile = TRIM(cdfile)
    iost = 0
    IF (cdacce(1 : 6) == 'DIRECT') THEN
      OPEN(UNIT = knum, FILE = clfile, FORM = cdform, ACCESS = cdacce, STATUS = cdstat, RECL = klengh, ERR = 100, IOSTAT = iost)
    ELSE IF (TRIM(cdstat) == 'APPEND') THEN
      OPEN(UNIT = knum, FILE = clfile, FORM = cdform, ACCESS = cdacce, STATUS = 'UNKNOWN', POSITION = 'APPEND', ERR = 100, IOSTAT = iost)
    ELSE
      OPEN(UNIT = knum, FILE = clfile, FORM = cdform, ACCESS = cdacce, STATUS = cdstat, ERR = 100, IOSTAT = iost)
    END IF
    IF (iost /= 0 .AND. TRIM(clfile) == '/dev/null') OPEN(UNIT = knum, FILE = 'NUL', FORM = cdform, ACCESS = cdacce, STATUS = cdstat, ERR = 100, IOSTAT = iost)
    IF (iost == 0) THEN
      IF (ldwp) THEN
        WRITE(kout, *) '     file   : ', TRIM(clfile), ' open ok'
        WRITE(kout, *) '     unit   = ', knum
        WRITE(kout, *) '     status = ', cdstat
        WRITE(kout, *) '     form   = ', cdform
        WRITE(kout, *) '     access = ', cdacce
        WRITE(kout, *)
      END IF
    END IF
100 CONTINUE
    IF (iost /= 0) THEN
      IF (ldwp) THEN
        WRITE(kout, *)
        WRITE(kout, *) ' ===>>>> : bad opening file: ', TRIM(clfile)
        WRITE(kout, *) ' =======   ===  '
        WRITE(kout, *) '           unit   = ', knum
        WRITE(kout, *) '           status = ', cdstat
        WRITE(kout, *) '           form   = ', cdform
        WRITE(kout, *) '           access = ', cdacce
        WRITE(kout, *) '           iostat = ', iost
        WRITE(kout, *) '           we stop. verify the file '
        WRITE(kout, *)
      ELSE
        WRITE(*, *)
        WRITE(*, *) ' ===>>>> : bad opening file: ', TRIM(clfile)
        WRITE(*, *) ' =======   ===  '
        WRITE(*, *) '           unit   = ', knum
        WRITE(*, *) '           status = ', cdstat
        WRITE(*, *) '           form   = ', cdform
        WRITE(*, *) '           access = ', cdacce
        WRITE(*, *) '           iostat = ', iost
        WRITE(*, *) '           we stop. verify the file '
        WRITE(*, *)
      END IF
      CALL FLUSH(kout)
      STOP 'ctl_opn bad opening'
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ctl_opn
  SUBROUTINE ctl_nam(kios, cdnam, ldwp)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(INOUT) :: kios
    CHARACTER(LEN = *), INTENT(IN) :: cdnam
    CHARACTER(LEN = 5) :: clios
    LOGICAL, INTENT(IN) :: ldwp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ctl_nam', 'r0', 0, 0)
    WRITE(clios, '(I5.0)') kios
    IF (kios < 0) THEN
      CALL ctl_warn('end of record or file while reading namelist ' // TRIM(cdnam) // ' iostat = ' // TRIM(clios))
    END IF
    IF (kios > 0) THEN
      CALL ctl_stop('misspelled variable in namelist ' // TRIM(cdnam) // ' iostat = ' // TRIM(clios))
    END IF
    kios = 0
    CALL profile_psy_data0 % PostEnd
    RETURN
  END SUBROUTINE ctl_nam
  INTEGER FUNCTION get_unit()
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    LOGICAL :: llopn
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('get_unit', 'r0', 0, 0)
    get_unit = 15
    llopn = .TRUE.
    DO WHILE ((get_unit < 998) .AND. llopn)
      get_unit = get_unit + 1
      INQUIRE(UNIT = get_unit, OPENED = llopn)
    END DO
    IF ((get_unit == 999) .AND. llopn) THEN
      CALL ctl_stop('get_unit: All logical units until 999 are used...')
      get_unit = - 1
    END IF
    CALL profile_psy_data0 % PostEnd
  END FUNCTION get_unit
END MODULE lib_mpp