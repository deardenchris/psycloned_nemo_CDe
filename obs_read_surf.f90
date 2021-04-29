MODULE obs_read_surf
  USE par_kind
  USE in_out_manager
  USE dom_oce
  USE obs_mpp
  USE julian
  USE obs_utils
  USE obs_grid
  USE obs_sort
  USE obs_surf_def
  USE obs_types
  USE obs_fbm
  USE netcdf
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: obs_rea_surf
  CONTAINS
  SUBROUTINE obs_rea_surf(surfdata, knumfiles, cdfilenames, kvars, kextr, kstp, ddobsini, ddobsend, ldignmis, ldmod, ldnightav)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(obs_surf), INTENT(INOUT) :: surfdata
    INTEGER, INTENT(IN) :: knumfiles
    CHARACTER(LEN = 128), INTENT(IN) :: cdfilenames(knumfiles)
    INTEGER, INTENT(IN) :: kvars
    INTEGER, INTENT(IN) :: kextr
    INTEGER, INTENT(IN) :: kstp
    LOGICAL, INTENT(IN) :: ldignmis
    LOGICAL, INTENT(IN) :: ldmod
    LOGICAL, INTENT(IN) :: ldnightav
    REAL(KIND = dp), INTENT(IN) :: ddobsini
    REAL(KIND = dp), INTENT(IN) :: ddobsend
    CHARACTER(LEN = 11), PARAMETER :: cpname = 'obs_rea_surf'
    CHARACTER(LEN = 8) :: clrefdate
    CHARACTER(LEN = 8), DIMENSION(:), ALLOCATABLE :: clvars
    INTEGER :: ji
    INTEGER :: jj
    INTEGER :: jk
    INTEGER :: iflag
    INTEGER :: inobf
    INTEGER :: i_file_id
    INTEGER :: inowin
    INTEGER :: iyea
    INTEGER :: imon
    INTEGER :: iday
    INTEGER :: ihou
    INTEGER :: imin
    INTEGER :: isec
    INTEGER :: itype
    INTEGER :: iobsmpp
    INTEGER :: iobs
    INTEGER :: iobstot
    INTEGER :: ios
    INTEGER :: ioserrcount
    INTEGER, PARAMETER :: jpsurfmaxtype = 1024
    INTEGER, DIMENSION(knumfiles) :: irefdate
    INTEGER, DIMENSION(jpsurfmaxtype + 1) :: ityp, itypmpp
    INTEGER, DIMENSION(:), ALLOCATABLE :: iobsi, iobsj, iproc, iindx, ifileidx, isurfidx
    REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: zphi, zlam
    REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: zdat
    REAL(KIND = wp), DIMENSION(knumfiles) :: djulini, djulend
    LOGICAL :: llvalprof
    TYPE(obfbdata), POINTER, DIMENSION(:) :: inpfiles
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    CALL profile_psy_data0 % PreStart('obs_rea_surf', 'r0', 0, 0)
    iobs = 0
    inobf = knumfiles
    ALLOCATE(inpfiles(inobf))
    CALL profile_psy_data0 % PostEnd
    surf_files:DO jj = 1, inobf
      CALL profile_psy_data1 % PreStart('obs_rea_surf', 'r1', 0, 0)
      IF (lwp) THEN
        WRITE(numout, *)
        WRITE(numout, *) ' obs_rea_surf : Reading from file = ', TRIM(TRIM(cdfilenames(jj)))
        WRITE(numout, *) ' ~~~~~~~~~~~'
        WRITE(numout, *)
      END IF
      iflag = nf90_open(TRIM(TRIM(cdfilenames(jj))), nf90_nowrite, i_file_id)
      CALL profile_psy_data1 % PostEnd
      IF (iflag /= nf90_noerr) THEN
        CALL profile_psy_data2 % PreStart('obs_rea_surf', 'r2', 0, 0)
        IF (ldignmis) THEN
          inpfiles(jj) % nobs = 0
          CALL ctl_warn('File ' // TRIM(TRIM(cdfilenames(jj))) // ' not found')
        ELSE
          CALL ctl_stop('File ' // TRIM(TRIM(cdfilenames(jj))) // ' not found')
        END IF
        CALL profile_psy_data2 % PostEnd
      ELSE
        CALL profile_psy_data3 % PreStart('obs_rea_surf', 'r3', 0, 0)
        iflag = nf90_close(i_file_id)
        CALL init_obfbdata(inpfiles(jj))
        CALL read_obfbdata(TRIM(cdfilenames(jj)), inpfiles(jj), ldgrid = .TRUE.)
        CALL profile_psy_data3 % PostEnd
        IF (ldmod .AND. (inpfiles(jj) % nadd == 0)) THEN
          CALL profile_psy_data4 % PreStart('obs_rea_surf', 'r4', 0, 0)
          CALL ctl_stop('Model not in input data')
          CALL profile_psy_data4 % PostEnd
          RETURN
        END IF
        CALL profile_psy_data5 % PreStart('obs_rea_surf', 'r5', 0, 0)
        IF (jj == 1) THEN
          ALLOCATE(clvars(inpfiles(jj) % nvar))
          DO ji = 1, inpfiles(jj) % nvar
            clvars(ji) = inpfiles(jj) % cname(ji)
          END DO
        ELSE
          DO ji = 1, inpfiles(jj) % nvar
            IF (inpfiles(jj) % cname(ji) /= clvars(ji)) THEN
              CALL ctl_stop('Feedback file variables not consistent', ' with previous files for this type')
            END IF
          END DO
        END IF
        IF (lwp) WRITE(numout, *) 'Observation file contains ', inpfiles(jj) % nobs, ' observations'
        DO ji = 1, inpfiles(jj) % nobs
          IF (inpfiles(jj) % plam(ji) < - 180.) inpfiles(jj) % plam(ji) = inpfiles(jj) % plam(ji) + 360.
          IF (inpfiles(jj) % plam(ji) > 180.) inpfiles(jj) % plam(ji) = inpfiles(jj) % plam(ji) - 360.
        END DO
        clrefdate = inpfiles(jj) % cdjuldref(1 : 8)
        READ(clrefdate, '(I8)') irefdate(jj)
        CALL ddatetoymdhms(ddobsini, iyea, imon, iday, ihou, imin, isec)
        CALL greg2jul(isec, imin, ihou, iday, imon, iyea, djulini(jj), krefdate = irefdate(jj))
        CALL ddatetoymdhms(ddobsend, iyea, imon, iday, ihou, imin, isec)
        CALL greg2jul(isec, imin, ihou, iday, imon, iyea, djulend(jj), krefdate = irefdate(jj))
        IF (ldnightav) THEN
          IF (lwp) THEN
            WRITE(numout, *) 'Resetting time of night-time averaged observations', ' to the end of the day'
          END IF
          DO ji = 1, inpfiles(jj) % nobs
            IF (inpfiles(jj) % ptim(ji) >= 0.) THEN
              inpfiles(jj) % ptim(ji) = INT(inpfiles(jj) % ptim(ji)) + 0.9999
            ELSE
              inpfiles(jj) % ptim(ji) = INT(inpfiles(jj) % ptim(ji)) - 0.0001
            END IF
          END DO
        END IF
        IF (inpfiles(jj) % nobs > 0) THEN
          inpfiles(jj) % iproc = - 1
          inpfiles(jj) % iobsi = - 1
          inpfiles(jj) % iobsj = - 1
        END IF
        inowin = 0
        DO ji = 1, inpfiles(jj) % nobs
          IF ((inpfiles(jj) % ptim(ji) > djulini(jj)) .AND. (inpfiles(jj) % ptim(ji) <= djulend(jj))) THEN
            inowin = inowin + 1
          END IF
        END DO
        ALLOCATE(zlam(inowin))
        ALLOCATE(zphi(inowin))
        ALLOCATE(iobsi(inowin))
        ALLOCATE(iobsj(inowin))
        ALLOCATE(iproc(inowin))
        inowin = 0
        DO ji = 1, inpfiles(jj) % nobs
          IF ((inpfiles(jj) % ptim(ji) > djulini(jj)) .AND. (inpfiles(jj) % ptim(ji) <= djulend(jj))) THEN
            inowin = inowin + 1
            zlam(inowin) = inpfiles(jj) % plam(ji)
            zphi(inowin) = inpfiles(jj) % pphi(ji)
          END IF
        END DO
        CALL obs_grid_search(inowin, zlam, zphi, iobsi, iobsj, iproc, 'T')
        inowin = 0
        DO ji = 1, inpfiles(jj) % nobs
          IF ((inpfiles(jj) % ptim(ji) > djulini(jj)) .AND. (inpfiles(jj) % ptim(ji) <= djulend(jj))) THEN
            inowin = inowin + 1
            inpfiles(jj) % iproc(ji, 1) = iproc(inowin)
            inpfiles(jj) % iobsi(ji, 1) = iobsi(inowin)
            inpfiles(jj) % iobsj(ji, 1) = iobsj(inowin)
          END IF
        END DO
        DEALLOCATE(zlam, zphi, iobsi, iobsj, iproc)
        DO ji = 1, inpfiles(jj) % nobs
          IF ((inpfiles(jj) % ptim(ji) > djulini(jj)) .AND. (inpfiles(jj) % ptim(ji) <= djulend(jj))) THEN
            IF (nproc == 0) THEN
              IF (inpfiles(jj) % iproc(ji, 1) > nproc) CYCLE
            ELSE
              IF (inpfiles(jj) % iproc(ji, 1) /= nproc) CYCLE
            END IF
            llvalprof = .FALSE.
            IF (.NOT. BTEST(inpfiles(jj) % ivlqc(1, ji, 1), 2)) THEN
              iobs = iobs + 1
            END IF
          END IF
        END DO
        CALL profile_psy_data5 % PostEnd
      END IF
    END DO surf_files
    CALL profile_psy_data6 % PreStart('obs_rea_surf', 'r6', 0, 0)
    iobstot = 0
    DO jj = 1, inobf
      DO ji = 1, inpfiles(jj) % nobs
        IF ((inpfiles(jj) % ptim(ji) > djulini(jj)) .AND. (inpfiles(jj) % ptim(ji) <= djulend(jj))) THEN
          iobstot = iobstot + 1
        END IF
      END DO
    END DO
    ALLOCATE(iindx(iobstot), ifileidx(iobstot), isurfidx(iobstot), zdat(iobstot))
    jk = 0
    DO jj = 1, inobf
      DO ji = 1, inpfiles(jj) % nobs
        IF ((inpfiles(jj) % ptim(ji) > djulini(jj)) .AND. (inpfiles(jj) % ptim(ji) <= djulend(jj))) THEN
          jk = jk + 1
          ifileidx(jk) = jj
          isurfidx(jk) = ji
          zdat(jk) = inpfiles(jj) % ptim(ji)
        END IF
      END DO
    END DO
    CALL sort_dp_indx(iobstot, zdat, iindx)
    CALL obs_surf_alloc(surfdata, iobs, kvars, kextr, kstp, jpi, jpj)
    iobs = 0
    surfdata % cvars(:) = clvars(:)
    ityp(:) = 0
    itypmpp(:) = 0
    ioserrcount = 0
    DO jk = 1, iobstot
      jj = ifileidx(iindx(jk))
      ji = isurfidx(iindx(jk))
      IF ((inpfiles(jj) % ptim(ji) > djulini(jj)) .AND. (inpfiles(jj) % ptim(ji) <= djulend(jj))) THEN
        IF (nproc == 0) THEN
          IF (inpfiles(jj) % iproc(ji, 1) > nproc) CYCLE
        ELSE
          IF (inpfiles(jj) % iproc(ji, 1) /= nproc) CYCLE
        END IF
        IF (.NOT. BTEST(inpfiles(jj) % ivlqc(1, ji, 1), 2)) THEN
          iobs = iobs + 1
          CALL jul2greg(isec, imin, ihou, iday, imon, iyea, inpfiles(jj) % ptim(ji), irefdate(jj))
          surfdata % nyea(iobs) = iyea
          surfdata % nmon(iobs) = imon
          surfdata % nday(iobs) = iday
          surfdata % nhou(iobs) = ihou
          surfdata % nmin(iobs) = imin
          surfdata % rlam(iobs) = inpfiles(jj) % plam(ji)
          surfdata % rphi(iobs) = inpfiles(jj) % pphi(ji)
          surfdata % mi(iobs) = inpfiles(jj) % iobsi(ji, 1)
          surfdata % mj(iobs) = inpfiles(jj) % iobsj(ji, 1)
          surfdata % cwmo(iobs) = inpfiles(jj) % cdwmo(ji)
          READ(inpfiles(jj) % cdtyp(ji), '(I4)', IOSTAT = ios, ERR = 901) itype
901       IF (ios /= 0) THEN
            IF (ioserrcount == 0) THEN
              CALL ctl_warn('Problem converting an instrument type ', 'to integer. Setting type to zero')
            END IF
            ioserrcount = ioserrcount + 1
            itype = 0
          END IF
          surfdata % ntyp(iobs) = itype
          IF (itype < jpsurfmaxtype + 1) THEN
            ityp(itype + 1) = ityp(itype + 1) + 1
          ELSE
            IF (lwp) WRITE(numout, *) 'WARNING:Increase jpsurfmaxtype in ', cpname
          END IF
          surfdata % nsidx(iobs) = iobs
          surfdata % nsfil(iobs) = iindx(jk)
          surfdata % nqc(iobs) = inpfiles(jj) % ivqc(ji, 1)
          surfdata % robs(iobs, 1) = inpfiles(jj) % pob(1, ji, 1)
          IF (ldmod) THEN
            surfdata % rmod(iobs, 1) = inpfiles(jj) % padd(1, ji, 1, 1)
            IF (TRIM(surfdata % cvars(1)) == 'SLA') THEN
              surfdata % rext(iobs, 1) = inpfiles(jj) % padd(1, ji, 2, 1)
              surfdata % rext(iobs, 2) = inpfiles(jj) % pext(1, ji, 1)
            END IF
          ELSE
            surfdata % rmod(iobs, 1) = fbrmdi
            IF (TRIM(surfdata % cvars(1)) == 'SLA') surfdata % rext(iobs, :) = fbrmdi
          END IF
        END IF
      END IF
    END DO
    CALL obs_mpp_sum_integer(iobs, iobsmpp)
    CALL obs_mpp_sum_integers(ityp, itypmpp, jpsurfmaxtype + 1)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, '(1X,A)') TRIM(surfdata % cvars(1)) // ' data'
      WRITE(numout, '(1X,A)') '--------------'
      DO jj = 1, 8
        IF (itypmpp(jj) > 0) THEN
          WRITE(numout, '(1X,A4,I4,A3,I10)') 'Type ', jj, ' = ', itypmpp(jj)
        END IF
      END DO
      WRITE(numout, '(1X,A)') '---------------------------------------------------------------'
      WRITE(numout, '(1X,A,I8)') 'Total data for variable ' // TRIM(surfdata % cvars(1)) // '           = ', iobsmpp
      WRITE(numout, '(1X,A)') '---------------------------------------------------------------'
      WRITE(numout, *)
    END IF
    DEALLOCATE(ifileidx, isurfidx, zdat, clvars)
    DO jj = 1, inobf
      IF (inpfiles(jj) % lalloc) THEN
        CALL dealloc_obfbdata(inpfiles(jj))
      END IF
    END DO
    DEALLOCATE(inpfiles)
    CALL profile_psy_data6 % PostEnd
  END SUBROUTINE obs_rea_surf
END MODULE obs_read_surf