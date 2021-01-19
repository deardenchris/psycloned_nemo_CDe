MODULE timing
  USE in_out_manager
  USE dom_oce
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: timing_init, timing_finalize
  PUBLIC :: timing_reset
  PUBLIC :: timing_start, timing_stop
  INCLUDE 'mpif.h'
  TYPE :: timer
    CHARACTER(LEN = 20) :: cname
    CHARACTER(LEN = 20) :: surname
    INTEGER :: rank
    REAL(KIND = wp) :: t_cpu, t_clock, tsum_cpu, tsum_clock, tmax_cpu, tmax_clock, tmin_cpu, tmin_clock, tsub_cpu, tsub_clock
    INTEGER :: ncount, ncount_max, ncount_rate
    INTEGER :: niter
    LOGICAL :: l_tdone
    TYPE(timer), POINTER :: next => NULL()
    TYPE(timer), POINTER :: prev => NULL()
    TYPE(timer), POINTER :: parent_section => NULL()
  END TYPE timer
  TYPE :: alltimer
    CHARACTER(LEN = 20), DIMENSION(:), POINTER :: cname => NULL()
    REAL(KIND = wp), DIMENSION(:), POINTER :: tsum_cpu => NULL()
    REAL(KIND = wp), DIMENSION(:), POINTER :: tsum_clock => NULL()
    INTEGER, DIMENSION(:), POINTER :: niter => NULL()
    TYPE(alltimer), POINTER :: next => NULL()
    TYPE(alltimer), POINTER :: prev => NULL()
  END TYPE alltimer
  TYPE(timer), POINTER :: s_timer_root => NULL()
  TYPE(timer), POINTER :: s_timer => NULL()
  TYPE(timer), POINTER :: s_timer_old => NULL()
  TYPE(timer), POINTER :: s_wrk => NULL()
  REAL(KIND = wp) :: t_overclock, t_overcpu
  LOGICAL :: l_initdone = .FALSE.
  INTEGER :: nsize
  REAL(KIND = wp) :: tot_etime, tot_ctime
  REAL(KIND = wp), DIMENSION(2) :: t_elaps, t_cpu
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:) :: all_etime, all_ctime
  INTEGER :: nfinal_count, ncount, ncount_rate, ncount_max
  INTEGER, DIMENSION(8) :: nvalues
  CHARACTER(LEN = 8), DIMENSION(2) :: cdate
  CHARACTER(LEN = 10), DIMENSION(2) :: ctime
  CHARACTER(LEN = 5) :: czone
  LOGICAL :: ln_onefile = .TRUE.
  LOGICAL :: lwriter
  CONTAINS
  SUBROUTINE timing_start(cdinfo)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdinfo
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('timing_start', 'r0', 0, 0)
    IF (ASSOCIATED(s_timer)) s_timer_old => s_timer
    CALL timing_ini_var(cdinfo)
    IF (.NOT. s_timer_old % l_tdone) THEN
      s_timer % parent_section => s_timer_old
    ELSE
      s_timer % parent_section => NULL()
    END IF
    s_timer % l_tdone = .FALSE.
    s_timer % niter = s_timer % niter + 1
    s_timer % t_cpu = 0.
    s_timer % t_clock = 0.
    CALL CPU_TIME(s_timer % t_cpu)
    s_timer % t_clock = MPI_Wtime()
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE timing_start
  SUBROUTINE timing_stop(cdinfo, csection)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdinfo
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: csection
    INTEGER :: ifinal_count, iperiods
    REAL(KIND = wp) :: zcpu_end, zmpitime, zcpu_raw, zclock_raw
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('timing_stop', 'r0', 0, 0)
    s_wrk => NULL()
    zmpitime = MPI_Wtime()
    CALL CPU_TIME(zcpu_end)
    zcpu_raw = zcpu_end - s_timer % t_cpu - t_overcpu
    s_timer % t_cpu = zcpu_raw - s_timer % tsub_cpu
    zclock_raw = zmpitime - s_timer % t_clock - t_overclock
    s_timer % t_clock = zclock_raw - t_overclock - s_timer % tsub_clock
    IF (.NOT. PRESENT(csection)) THEN
      IF (ASSOCIATED(s_timer % parent_section)) THEN
        s_timer % parent_section % tsub_cpu = zcpu_raw + s_timer % parent_section % tsub_cpu
        s_timer % parent_section % tsub_clock = zclock_raw + s_timer % parent_section % tsub_clock
      END IF
    END IF
    s_timer % tsum_clock = s_timer % tsum_clock + s_timer % t_clock
    s_timer % tsum_cpu = s_timer % tsum_cpu + s_timer % t_cpu
    s_timer % tsub_clock = 0.
    s_timer % tsub_cpu = 0.
    s_timer % l_tdone = .TRUE.
    IF (ASSOCIATED(s_timer % parent_section)) s_timer => s_timer % parent_section
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE timing_stop
  SUBROUTINE timing_init
    INTEGER :: iperiods, istart_count, ifinal_count
    REAL(KIND = wp) :: zdum
    LOGICAL :: ll_f
    IF (ln_onefile) THEN
      IF (lwp) CALL ctl_opn(numtime, 'timing.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, numout, .TRUE., narea)
      lwriter = lwp
    ELSE
      CALL ctl_opn(numtime, 'timing.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, numout, .FALSE., narea)
      lwriter = .TRUE.
    END IF
    IF (lwriter) THEN
      WRITE(numtime, FMT = *)
      WRITE(numtime, FMT = *) '      CNRS - NERC - Met OFFICE - MERCATOR-ocean - CMCC - INGV'
      WRITE(numtime, FMT = *) '                             NEMO team'
      WRITE(numtime, FMT = *) '                  Ocean General Circulation Model'
      WRITE(numtime, FMT = *) '                        version 4.0  (2019) '
      WRITE(numtime, FMT = *)
      WRITE(numtime, FMT = *) '                        Timing Informations '
      WRITE(numtime, FMT = *)
      WRITE(numtime, FMT = *)
    END IF
    t_overclock = MPI_WTIME()
    t_overclock = MPI_WTIME() - t_overclock
    CALL CPU_TIME(zdum)
    CALL CPU_TIME(t_overcpu)
    t_overcpu = t_overcpu - zdum
    t_overclock = t_overcpu + t_overclock
    CALL DATE_AND_TIME(cdate(1), ctime(1), czone, nvalues)
    CALL CPU_TIME(t_cpu(1))
    t_elaps(1) = MPI_WTIME()
  END SUBROUTINE timing_init
  SUBROUTINE timing_finalize
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(timer), POINTER :: s_temp
    INTEGER :: idum, iperiods, icode
    INTEGER :: ji
    LOGICAL :: ll_ord, ll_averep
    CHARACTER(LEN = 120) :: clfmt
    REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: timing_glob
    REAL(KIND = wp) :: zsypd
    REAL(KIND = wp) :: zperc, ztot
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('timing_finalize', 'r0', 0, 0)
    ll_averep = .TRUE.
    CALL CPU_TIME(t_cpu(2))
    t_cpu(2) = t_cpu(2) - t_cpu(1) - t_overcpu
    t_elaps(2) = MPI_WTIME() - t_elaps(1) - t_overclock
    CALL DATE_AND_TIME(cdate(2), ctime(2), czone, nvalues)
    nsize = 0
    s_timer => s_timer_root
    DO WHILE (ASSOCIATED(s_timer))
      nsize = nsize + 1
      s_timer => s_timer % next
    END DO
    idum = nsize
    CALL mpp_sum('timing', idum)
    IF (idum / jpnij /= nsize) THEN
      IF (lwriter) WRITE(numtime, FMT = *) '        ===> W A R N I N G: '
      IF (lwriter) WRITE(numtime, FMT = *) ' Some CPU have different number of routines instrumented for timing'
      IF (lwriter) WRITE(numtime, FMT = *) ' No detailed report on averaged timing can be provided'
      IF (lwriter) WRITE(numtime, FMT = *) ' The following detailed report only deals with the current processor'
      IF (lwriter) WRITE(numtime, FMT = *)
      ll_averep = .FALSE.
    END IF
    ALLOCATE(all_etime(jpnij), all_ctime(jpnij))
    CALL MPI_ALLGATHER(t_elaps(2), 1, MPI_DOUBLE_PRECISION, all_etime, 1, MPI_DOUBLE_PRECISION, MPI_COMM_OCE, icode)
    CALL MPI_ALLGATHER(t_cpu(2), 1, MPI_DOUBLE_PRECISION, all_ctime, 1, MPI_DOUBLE_PRECISION, MPI_COMM_OCE, icode)
    tot_etime = SUM(all_etime(:))
    tot_ctime = SUM(all_ctime(:))
    IF (lwriter) WRITE(numtime, FMT = *) 'Total timing (sum) :'
    IF (lwriter) WRITE(numtime, FMT = *) '--------------------'
    IF (lwriter) WRITE(numtime, FMT = "('Elapsed Time (s)  CPU Time (s)')")
    IF (lwriter) WRITE(numtime, FMT = '(5x,f12.3,1x,f12.3)') tot_etime, tot_ctime
    IF (lwriter) WRITE(numtime, FMT = *)
    IF (ll_averep) CALL waver_info
    CALL wmpi_info
    IF (lwriter) CALL wcurrent_info
    clfmt = '(1X,"Timing started on ",2(A2,"/"),A4," at ",2(A2,":"),A2," MET ",A3,":",A2," from GMT")'
    IF (lwriter) WRITE(numtime, TRIM(clfmt)) cdate(1)(7 : 8), cdate(1)(5 : 6), cdate(1)(1 : 4), ctime(1)(1 : 2), ctime(1)(3 : 4), &
&ctime(1)(5 : 6), czone(1 : 3), czone(4 : 5)
    clfmt = '(1X,  "Timing   ended on ",2(A2,"/"),A4," at ",2(A2,":"),A2," MET ",A3,":",A2," from GMT")'
    IF (lwriter) WRITE(numtime, TRIM(clfmt)) cdate(2)(7 : 8), cdate(2)(5 : 6), cdate(2)(1 : 4), ctime(2)(1 : 2), ctime(2)(3 : 4), &
&ctime(2)(5 : 6), czone(1 : 3), czone(4 : 5)
    ALLOCATE(timing_glob(4 * jpnij), STAT = icode)
    CALL MPI_GATHER((/compute_time, waiting_time(1), waiting_time(2), elapsed_time/), 4, MPI_DOUBLE_PRECISION, timing_glob, 4, &
&MPI_DOUBLE_PRECISION, 0, MPI_COMM_OCE, icode)
    IF (narea == 1) THEN
      WRITE(numtime, FMT = *) ' '
      WRITE(numtime, FMT = *) ' Report on time spent on waiting MPI messages '
      WRITE(numtime, FMT = *) '    total timing measured between nit000+1 and nitend-1 '
      WRITE(numtime, FMT = *) '    warning: includes restarts writing time if output before nitend... '
      WRITE(numtime, FMT = *) ' '
      DO ji = 1, jpnij
        ztot = SUM(timing_glob(4 * ji - 3 : 4 * ji - 1))
        WRITE(numtime, FMT = '(A28,F11.6,            A34,I8)') 'Computing       time : ', timing_glob(4 * ji - 3), ' on MPI rank : &
&', ji
        IF (ztot /= 0.) zperc = timing_glob(4 * ji - 2) / ztot * 100.
        WRITE(numtime, FMT = '(A28,F11.6,A2, F4.1,A3,A25,I8)') 'Waiting lbc_lnk time : ', timing_glob(4 * ji - 2), ' (', zperc, ' &
&%)', ' on MPI rank : ', ji
        IF (ztot /= 0.) zperc = timing_glob(4 * ji - 1) / ztot * 100.
        WRITE(numtime, FMT = '(A28,F11.6,A2, F4.1,A3,A25,I8)') 'Waiting  global time : ', timing_glob(4 * ji - 1), ' (', zperc, ' &
&%)', ' on MPI rank : ', ji
        zsypd = rn_rdt * REAL(nitend - nit000 - 1, wp) / (timing_glob(4 * ji) * 365.)
        WRITE(numtime, FMT = '(A28,F11.6,A7,F10.3,A2,A15,I8)') 'Total           time : ', timing_glob(4 * ji), ' (SYPD: ', zsypd, &
&')', ' on MPI rank : ', ji
      END DO
    END IF
    DEALLOCATE(timing_glob)
    IF (lwriter) CLOSE(UNIT = numtime)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE timing_finalize
  SUBROUTINE wcurrent_info
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    LOGICAL :: ll_ord
    CHARACTER(LEN = 2048) :: clfmt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wcurrent_info', 'r0', 0, 0)
    s_wrk => NULL()
    s_timer => s_timer_root
    DO
      ll_ord = .TRUE.
      s_timer => s_timer_root
      DO WHILE (ASSOCIATED(s_timer % next))
        IF (.NOT. ASSOCIATED(s_timer % next)) EXIT
        IF (s_timer % tsum_clock < s_timer % next % tsum_clock) THEN
          ALLOCATE(s_wrk)
          s_wrk = s_timer % next
          CALL insert(s_timer, s_timer_root, s_wrk)
          CALL suppress(s_timer % next)
          ll_ord = .FALSE.
          CYCLE
        END IF
        IF (ASSOCIATED(s_timer % next)) s_timer => s_timer % next
      END DO
      IF (ll_ord) EXIT
    END DO
    WRITE(numtime, FMT = *) 'Detailed timing for proc :', narea - 1
    WRITE(numtime, FMT = *) '--------------------------'
    WRITE(numtime, FMT = *) 'Section             ', 'Elapsed Time (s)  ', 'Elapsed Time (%)  ', 'CPU Time(s)  ', 'CPU Time (%)  ', &
&'CPU/Elapsed  ', 'Frequency'
    s_timer => s_timer_root
    clfmt = '(1x,a,4x,f12.3,6x,f12.3,x,f12.3,2x,f12.3,6x,f7.3,2x,i9)'
    DO WHILE (ASSOCIATED(s_timer))
      WRITE(numtime, TRIM(clfmt)) s_timer % cname, s_timer % tsum_clock, s_timer % tsum_clock * 100. / t_elaps(2), s_timer % &
&tsum_cpu, s_timer % tsum_cpu * 100. / t_cpu(2), s_timer % tsum_cpu / s_timer % tsum_clock, s_timer % niter
      s_timer => s_timer % next
    END DO
    WRITE(numtime, FMT = *)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE wcurrent_info
  SUBROUTINE waver_info
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(alltimer), POINTER :: sl_timer_glob_root => NULL()
    TYPE(alltimer), POINTER :: sl_timer_glob => NULL()
    TYPE(timer), POINTER :: sl_timer_ave_root => NULL()
    TYPE(timer), POINTER :: sl_timer_ave => NULL()
    INTEGER :: icode
    INTEGER :: ierr
    LOGICAL :: ll_ord
    CHARACTER(LEN = 200) :: clfmt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    ALLOCATE(sl_timer_glob_root, STAT = ierr)
    IF (ierr /= 0) THEN
      WRITE(numtime, FMT = *) 'Failed to allocate global timing structure in waver_info'
      RETURN
    END IF
    ALLOCATE(sl_timer_glob_root % cname(jpnij), sl_timer_glob_root % tsum_cpu(jpnij), sl_timer_glob_root % tsum_clock(jpnij), &
&sl_timer_glob_root % niter(jpnij), STAT = ierr)
    IF (ierr /= 0) THEN
      WRITE(numtime, FMT = *) 'Failed to allocate global timing structure in waver_info'
      RETURN
    END IF
    CALL profile_psy_data0 % PreStart('waver_info', 'r0', 0, 0)
    sl_timer_glob_root % cname(:) = ''
    sl_timer_glob_root % tsum_cpu(:) = 0._wp
    sl_timer_glob_root % tsum_clock(:) = 0._wp
    sl_timer_glob_root % niter(:) = 0
    sl_timer_glob_root % next => NULL()
    sl_timer_glob_root % prev => NULL()
    sl_timer_glob => sl_timer_glob_root
    IF (narea .EQ. 1) THEN
      ALLOCATE(sl_timer_ave_root)
      sl_timer_ave_root % cname = ''
      sl_timer_ave_root % t_cpu = 0._wp
      sl_timer_ave_root % t_clock = 0._wp
      sl_timer_ave_root % tsum_cpu = 0._wp
      sl_timer_ave_root % tsum_clock = 0._wp
      sl_timer_ave_root % tmax_cpu = 0._wp
      sl_timer_ave_root % tmax_clock = 0._wp
      sl_timer_ave_root % tmin_cpu = 0._wp
      sl_timer_ave_root % tmin_clock = 0._wp
      sl_timer_ave_root % tsub_cpu = 0._wp
      sl_timer_ave_root % tsub_clock = 0._wp
      sl_timer_ave_root % ncount = 0
      sl_timer_ave_root % ncount_rate = 0
      sl_timer_ave_root % ncount_max = 0
      sl_timer_ave_root % niter = 0
      sl_timer_ave_root % l_tdone = .FALSE.
      sl_timer_ave_root % next => NULL()
      sl_timer_ave_root % prev => NULL()
      ALLOCATE(sl_timer_ave)
      sl_timer_ave => sl_timer_ave_root
    END IF
    s_timer => s_timer_root
    DO WHILE (ASSOCIATED(s_timer))
      CALL MPI_GATHER(s_timer % cname, 20, MPI_CHARACTER, sl_timer_glob % cname, 20, MPI_CHARACTER, 0, MPI_COMM_OCE, icode)
      CALL MPI_GATHER(s_timer % tsum_clock, 1, MPI_DOUBLE_PRECISION, sl_timer_glob % tsum_clock, 1, MPI_DOUBLE_PRECISION, 0, &
&MPI_COMM_OCE, icode)
      CALL MPI_GATHER(s_timer % tsum_cpu, 1, MPI_DOUBLE_PRECISION, sl_timer_glob % tsum_cpu, 1, MPI_DOUBLE_PRECISION, 0, &
&MPI_COMM_OCE, icode)
      CALL MPI_GATHER(s_timer % niter, 1, MPI_INTEGER, sl_timer_glob % niter, 1, MPI_INTEGER, 0, MPI_COMM_OCE, icode)
      IF (narea == 1 .AND. ASSOCIATED(s_timer % next)) THEN
        ALLOCATE(sl_timer_glob % next)
        ALLOCATE(sl_timer_glob % next % cname(jpnij))
        ALLOCATE(sl_timer_glob % next % tsum_cpu(jpnij))
        ALLOCATE(sl_timer_glob % next % tsum_clock(jpnij))
        ALLOCATE(sl_timer_glob % next % niter(jpnij))
        sl_timer_glob % next % prev => sl_timer_glob
        sl_timer_glob % next % next => NULL()
        sl_timer_glob => sl_timer_glob % next
      END IF
      s_timer => s_timer % next
    END DO
    IF (narea == 1) THEN
      sl_timer_glob => sl_timer_glob_root
      DO WHILE (ASSOCIATED(sl_timer_glob))
        sl_timer_ave % cname = sl_timer_glob % cname(1)
        sl_timer_ave % tsum_cpu = SUM(sl_timer_glob % tsum_cpu(:)) / jpnij
        sl_timer_ave % tsum_clock = SUM(sl_timer_glob % tsum_clock(:)) / jpnij
        sl_timer_ave % tmax_cpu = MAXVAL(sl_timer_glob % tsum_cpu(:))
        sl_timer_ave % tmax_clock = MAXVAL(sl_timer_glob % tsum_clock(:))
        sl_timer_ave % tmin_cpu = MINVAL(sl_timer_glob % tsum_cpu(:))
        sl_timer_ave % tmin_clock = MINVAL(sl_timer_glob % tsum_clock(:))
        sl_timer_ave % niter = SUM(sl_timer_glob % niter(:))
        IF (ASSOCIATED(sl_timer_glob % next)) THEN
          ALLOCATE(sl_timer_ave % next)
          sl_timer_ave % next % prev => sl_timer_ave
          sl_timer_ave % next % next => NULL()
          sl_timer_ave => sl_timer_ave % next
        END IF
        sl_timer_glob => sl_timer_glob % next
      END DO
      s_wrk => NULL()
      sl_timer_ave => sl_timer_ave_root
      DO
        ll_ord = .TRUE.
        sl_timer_ave => sl_timer_ave_root
        DO WHILE (ASSOCIATED(sl_timer_ave % next))
          IF (.NOT. ASSOCIATED(sl_timer_ave % next)) EXIT
          IF (sl_timer_ave % tsum_clock < sl_timer_ave % next % tsum_clock) THEN
            ALLOCATE(s_wrk)
            s_wrk = sl_timer_ave % next
            CALL insert(sl_timer_ave, sl_timer_ave_root, s_wrk)
            CALL suppress(sl_timer_ave % next)
            ll_ord = .FALSE.
            CYCLE
          END IF
          IF (ASSOCIATED(sl_timer_ave % next)) sl_timer_ave => sl_timer_ave % next
        END DO
        IF (ll_ord) EXIT
      END DO
      WRITE(numtime, FMT = "('Averaged timing on all processors :')")
      WRITE(numtime, FMT = "('-----------------------------------')")
      WRITE(numtime, FMT = "('Section',13x,'Elap. Time(s)',2x,'Elap. Time(%)',2x,    'CPU Time(s)',2x,'CPU &
&Time(%)',2x,'CPU/Elap',1x,      'Max elap(%)',2x,'Min elap(%)',2x,               'Freq')")
      sl_timer_ave => sl_timer_ave_root
      clfmt = '((A),E15.7,2x,f6.2,5x,f12.2,5x,f6.2,5x,f7.2,2x,f12.2,4x,f6.2,2x,f9.2)'
      DO WHILE (ASSOCIATED(sl_timer_ave))
        WRITE(numtime, TRIM(clfmt)) sl_timer_ave % cname(1 : 18), sl_timer_ave % tsum_clock, sl_timer_ave % tsum_clock * 100. * &
&jpnij / tot_etime, sl_timer_ave % tsum_cpu, sl_timer_ave % tsum_cpu * 100. * jpnij / tot_ctime, sl_timer_ave % tsum_cpu / &
&sl_timer_ave % tsum_clock, sl_timer_ave % tmax_clock * 100. * jpnij / tot_etime, sl_timer_ave % tmin_clock * 100. * jpnij / &
&tot_etime, sl_timer_ave % niter / REAL(jpnij)
        sl_timer_ave => sl_timer_ave % next
      END DO
      WRITE(numtime, FMT = *)
      DEALLOCATE(sl_timer_ave_root)
    END IF
    DEALLOCATE(sl_timer_glob_root)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE waver_info
  SUBROUTINE wmpi_info
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: idum, icode
    INTEGER, ALLOCATABLE, DIMENSION(:) :: iall_rank
    REAL(KIND = wp) :: ztot_ratio
    REAL(KIND = wp) :: zmax_etime, zmax_ctime, zmax_ratio, zmin_etime, zmin_ctime, zmin_ratio
    REAL(KIND = wp) :: zavg_etime, zavg_ctime, zavg_ratio
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:) :: zall_ratio
    CHARACTER(LEN = 128), DIMENSION(8) :: cllignes
    CHARACTER(LEN = 128) :: clhline, clstart_date, clfinal_date
    CHARACTER(LEN = 2048) :: clfmt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wmpi_info', 'r0', 0, 0)
    ALLOCATE(zall_ratio(jpnij), iall_rank(jpnij))
    IF (narea == 1) THEN
      iall_rank(:) = (/(idum, idum = 0, jpnij - 1)/)
      zavg_etime = tot_etime / REAL(jpnij, wp)
      zmax_etime = MAXVAL(all_etime(:))
      zmin_etime = MINVAL(all_etime(:))
      zavg_ctime = tot_ctime / REAL(jpnij, wp)
      zmax_ctime = MAXVAL(all_ctime(:))
      zmin_ctime = MINVAL(all_ctime(:))
      zall_ratio(:) = all_ctime(:) / all_etime(:)
      ztot_ratio = SUM(zall_ratio(:))
      zavg_ratio = ztot_ratio / REAL(jpnij, wp)
      zmax_ratio = MAXVAL(zall_ratio(:))
      zmin_ratio = MINVAL(zall_ratio(:))
      clhline = '1x,13("-"),"|",18("-"),"|",14("-"),"|",18("-"),/,'
      cllignes(1) = '(1x,"MPI summary report :",/,'
      cllignes(2) = '1x,"--------------------",//,'
      cllignes(3) = '1x,"Process Rank |"," Elapsed Time (s) |"," CPU Time (s) |"," Ratio CPU/Elapsed",/,'
      cllignes(4) = '      (1x,i6,7x,"|",f12.3,6x,"|",f12.3,2x,"|",4x,f7.3,/),'
      WRITE(cllignes(4)(1 : 6), FMT = '(I6)') jpnij
      cllignes(5) = '1x,"Total        |",f12.3,6x,"|",F12.3,2x,"|",4x,f7.3,/,'
      cllignes(6) = '1x,"Minimum      |",f12.3,6x,"|",F12.3,2x,"|",4x,f7.3,/,'
      cllignes(7) = '1x,"Maximum      |",f12.3,6x,"|",F12.3,2x,"|",4x,f7.3,/,'
      cllignes(8) = '1x,"Average      |",f12.3,6x,"|",F12.3,2x,"|",4x,f7.3)'
      clfmt = TRIM(cllignes(1)) // TRIM(cllignes(2)) // TRIM(cllignes(3)) // TRIM(clhline) // TRIM(cllignes(4)) // TRIM(clhline) &
&// TRIM(cllignes(5)) // TRIM(clhline) // TRIM(cllignes(6)) // TRIM(clhline) // TRIM(cllignes(7)) // TRIM(clhline) // &
&TRIM(cllignes(8))
      WRITE(numtime, TRIM(clfmt)) (iall_rank(idum), all_etime(idum), all_ctime(idum), zall_ratio(idum), idum = 1, jpnij), &
&tot_etime, tot_ctime, ztot_ratio, zmin_etime, zmin_ctime, zmin_ratio, zmax_etime, zmax_ctime, zmax_ratio, zavg_etime, zavg_ctime, &
&zavg_ratio
      WRITE(numtime, FMT = *)
    END IF
    DEALLOCATE(zall_ratio, iall_rank)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE wmpi_info
  SUBROUTINE timing_ini_var(cdinfo)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdinfo
    LOGICAL :: ll_section
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('timing_ini_var', 'r0', 0, 0)
    IF (.NOT. ASSOCIATED(s_timer_root)) THEN
      ALLOCATE(s_timer_root)
      s_timer_root % cname = cdinfo
      s_timer_root % t_cpu = 0._wp
      s_timer_root % t_clock = 0._wp
      s_timer_root % tsum_cpu = 0._wp
      s_timer_root % tsum_clock = 0._wp
      s_timer_root % tmax_cpu = 0._wp
      s_timer_root % tmax_clock = 0._wp
      s_timer_root % tmin_cpu = 0._wp
      s_timer_root % tmin_clock = 0._wp
      s_timer_root % tsub_cpu = 0._wp
      s_timer_root % tsub_clock = 0._wp
      s_timer_root % ncount = 0
      s_timer_root % ncount_rate = 0
      s_timer_root % ncount_max = 0
      s_timer_root % niter = 0
      s_timer_root % l_tdone = .FALSE.
      s_timer_root % next => NULL()
      s_timer_root % prev => NULL()
      s_timer => s_timer_root
      ALLOCATE(s_wrk)
      s_wrk => NULL()
      ALLOCATE(s_timer_old)
      s_timer_old % cname = cdinfo
      s_timer_old % t_cpu = 0._wp
      s_timer_old % t_clock = 0._wp
      s_timer_old % tsum_cpu = 0._wp
      s_timer_old % tsum_clock = 0._wp
      s_timer_old % tmax_cpu = 0._wp
      s_timer_old % tmax_clock = 0._wp
      s_timer_old % tmin_cpu = 0._wp
      s_timer_old % tmin_clock = 0._wp
      s_timer_old % tsub_cpu = 0._wp
      s_timer_old % tsub_clock = 0._wp
      s_timer_old % ncount = 0
      s_timer_old % ncount_rate = 0
      s_timer_old % ncount_max = 0
      s_timer_old % niter = 0
      s_timer_old % l_tdone = .TRUE.
      s_timer_old % next => NULL()
      s_timer_old % prev => NULL()
    ELSE
      s_timer => s_timer_root
      DO WHILE (ASSOCIATED(s_timer))
        IF (TRIM(s_timer % cname) .EQ. TRIM(cdinfo)) THEN
          RETURN
        END IF
        s_timer => s_timer % next
      END DO
      s_timer => s_timer_root
      DO WHILE (ASSOCIATED(s_timer % next))
        s_timer => s_timer % next
      END DO
      ALLOCATE(s_timer % next)
      s_timer % next % cname = cdinfo
      s_timer % next % t_cpu = 0._wp
      s_timer % next % t_clock = 0._wp
      s_timer % next % tsum_cpu = 0._wp
      s_timer % next % tsum_clock = 0._wp
      s_timer % next % tmax_cpu = 0._wp
      s_timer % next % tmax_clock = 0._wp
      s_timer % next % tmin_cpu = 0._wp
      s_timer % next % tmin_clock = 0._wp
      s_timer % next % tsub_cpu = 0._wp
      s_timer % next % tsub_clock = 0._wp
      s_timer % next % ncount = 0
      s_timer % next % ncount_rate = 0
      s_timer % next % ncount_max = 0
      s_timer % next % niter = 0
      s_timer % next % l_tdone = .FALSE.
      s_timer % next % parent_section => NULL()
      s_timer % next % prev => s_timer
      s_timer % next % next => NULL()
      s_timer => s_timer % next
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE timing_ini_var
  SUBROUTINE timing_reset
    l_initdone = .TRUE.
    CALL timing_list(s_timer_root)
  END SUBROUTINE timing_reset
  RECURSIVE SUBROUTINE timing_list(ptr)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(timer), POINTER, INTENT(INOUT) :: ptr
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('timing_list', 'r0', 0, 0)
    IF (ASSOCIATED(ptr % next)) CALL timing_list(ptr % next)
    IF (lwp) WRITE(numout, FMT = *) '   ', ptr % cname
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE timing_list
  SUBROUTINE insert(sd_current, sd_root, sd_ptr)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(timer), POINTER, INTENT(INOUT) :: sd_current, sd_root, sd_ptr
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('insert', 'r0', 0, 0)
    IF (ASSOCIATED(sd_current, sd_root)) THEN
      sd_root => sd_ptr
    ELSE
      sd_current % prev % next => sd_ptr
    END IF
    sd_ptr % next => sd_current
    sd_ptr % prev => sd_current % prev
    sd_current % prev => sd_ptr
    sd_ptr => NULL()
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE insert
  SUBROUTINE suppress(sd_ptr)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(timer), POINTER, INTENT(INOUT) :: sd_ptr
    TYPE(timer), POINTER :: sl_temp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('suppress', 'r0', 0, 0)
    sl_temp => sd_ptr
    sd_ptr => sd_ptr % next
    IF (ASSOCIATED(sl_temp % next)) sl_temp % next % prev => sl_temp % prev
    DEALLOCATE(sl_temp)
    sl_temp => NULL()
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE suppress
END MODULE timing