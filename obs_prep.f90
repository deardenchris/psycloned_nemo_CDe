MODULE obs_prep
  USE par_kind, ONLY: wp
  USE in_out_manager
  USE obs_profiles_def
  USE obs_surf_def
  USE obs_mpp, ONLY: obs_mpp_sum_integer, obs_mpp_sum_integers
  USE obs_inter_sup
  USE obs_oper
  USE lib_mpp, ONLY: ctl_warn, ctl_stop
  USE bdy_oce, ONLY: idx_bdy, nb_bdy, ln_bdy
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: obs_pre_prof
  PUBLIC :: obs_pre_surf
  PUBLIC :: calc_month_len
  CONTAINS
  SUBROUTINE obs_pre_surf(surfdata, surfdataqc, ld_nea, ld_bound_reject, kqc_cutoff)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE par_oce
    USE dom_oce, ONLY: glamt, gphit, tmask, nproc
    TYPE(obs_surf), INTENT(INOUT) :: surfdata
    TYPE(obs_surf), INTENT(INOUT) :: surfdataqc
    LOGICAL, INTENT(IN) :: ld_nea
    LOGICAL, INTENT(IN) :: ld_bound_reject
    INTEGER, INTENT(IN), OPTIONAL :: kqc_cutoff
    INTEGER :: iqc_cutoff = 255
    INTEGER :: iyea0
    INTEGER :: imon0
    INTEGER :: iday0
    INTEGER :: ihou0
    INTEGER :: imin0
    INTEGER :: icycle
    INTEGER :: iotdobs
    INTEGER :: iosdsobs
    INTEGER :: ilansobs
    INTEGER :: inlasobs
    INTEGER :: igrdobs
    INTEGER :: ibdysobs
    INTEGER :: iotdobsmpp
    INTEGER :: iosdsobsmpp
    INTEGER :: ilansobsmpp
    INTEGER :: inlasobsmpp
    INTEGER :: igrdobsmpp
    INTEGER :: ibdysobsmpp
    LOGICAL, DIMENSION(:), ALLOCATABLE :: llvalid
    INTEGER :: jobs
    INTEGER :: jstp
    INTEGER :: inrc
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_pre_surf', 'r0', 0, 0)
    IF (lwp) WRITE(numout, FMT = *) 'obs_pre_surf : Preparing the surface observations...'
    IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~~'
    iyea0 = ndate0 / 10000
    imon0 = (ndate0 - iyea0 * 10000) / 100
    iday0 = ndate0 - iyea0 * 10000 - imon0 * 100
    ihou0 = nn_time0 / 100
    imin0 = (nn_time0 - ihou0 * 100)
    icycle = nn_no
    iotdobs = 0
    igrdobs = 0
    iosdsobs = 0
    ilansobs = 0
    inlasobs = 0
    ibdysobs = 0
    IF (PRESENT(kqc_cutoff)) iqc_cutoff = kqc_cutoff
    CALL obs_coo_tim(icycle, iyea0, imon0, iday0, ihou0, imin0, surfdata % nsurf, surfdata % nyea, surfdata % nmon, &
&surfdata % nday, surfdata % nhou, surfdata % nmin, surfdata % nqc, surfdata % mstp, iotdobs)
    CALL obs_mpp_sum_integer(iotdobs, iotdobsmpp)
    CALL obs_coo_grd(surfdata % nsurf, surfdata % mi, surfdata % mj, surfdata % nqc, igrdobs)
    CALL obs_mpp_sum_integer(igrdobs, igrdobsmpp)
    CALL obs_coo_spc_2d(surfdata % nsurf, jpi, jpj, surfdata % mi, surfdata % mj, surfdata % rlam, surfdata % rphi, glamt, gphit, &
&tmask(:, :, 1), surfdata % nqc, iosdsobs, ilansobs, inlasobs, ld_nea, ibdysobs, ld_bound_reject, iqc_cutoff)
    CALL obs_mpp_sum_integer(iosdsobs, iosdsobsmpp)
    CALL obs_mpp_sum_integer(ilansobs, ilansobsmpp)
    CALL obs_mpp_sum_integer(inlasobs, inlasobsmpp)
    CALL obs_mpp_sum_integer(ibdysobs, ibdysobsmpp)
    ALLOCATE(llvalid(surfdata % nsurf))
    llvalid(:) = (surfdata % nqc(:) <= iqc_cutoff)
    CALL obs_surf_compress(surfdata, surfdataqc, .TRUE., numout, lvalid = llvalid)
    DEALLOCATE(llvalid)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) ' ' // surfdataqc % cvars(1) // ' data outside time domain                  = ', iotdobsmpp
      WRITE(numout, FMT = *) ' Remaining ' // surfdataqc % cvars(1) // ' data that failed grid search    = ', igrdobsmpp
      WRITE(numout, FMT = *) ' Remaining ' // surfdataqc % cvars(1) // ' data outside space domain       = ', iosdsobsmpp
      WRITE(numout, FMT = *) ' Remaining ' // surfdataqc % cvars(1) // ' data at land points             = ', ilansobsmpp
      IF (ld_nea) THEN
        WRITE(numout, FMT = *) ' Remaining ' // surfdataqc % cvars(1) // ' data near land points (removed) = ', inlasobsmpp
      ELSE
        WRITE(numout, FMT = *) ' Remaining ' // surfdataqc % cvars(1) // ' data near land points (kept)    = ', inlasobsmpp
      END IF
      WRITE(numout, FMT = *) ' Remaining ' // surfdataqc % cvars(1) // ' data near open boundary (removed) = ', ibdysobsmpp
      WRITE(numout, FMT = *) ' ' // surfdataqc % cvars(1) // ' data accepted                             = ', surfdataqc % nsurfmpp
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) ' Number of observations per time step :'
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = '(10X,A,10X,A)') 'Time step', surfdataqc % cvars(1)
      WRITE(numout, FMT = '(10X,A,5X,A)') '---------', '-----------------'
      CALL FLUSH(numout)
    END IF
    DO jobs = 1, surfdataqc % nsurf
      inrc = surfdataqc % mstp(jobs) + 2 - nit000
      surfdataqc % nsstp(inrc) = surfdataqc % nsstp(inrc) + 1
    END DO
    CALL obs_mpp_sum_integers(surfdataqc % nsstp, surfdataqc % nsstpmpp, nitend - nit000 + 2)
    IF (lwp) THEN
      DO jstp = nit000 - 1, nitend
        inrc = jstp - nit000 + 2
        WRITE(numout, 1999) jstp, surfdataqc % nsstpmpp(inrc)
        CALL FLUSH(numout)
      END DO
    END IF
1999 FORMAT(10X, I9, 5X, I17)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_pre_surf
  SUBROUTINE obs_pre_prof(profdata, prodatqc, ld_var1, ld_var2, kpi, kpj, kpk, zmask1, pglam1, pgphi1, zmask2, pglam2, pgphi2, &
&ld_nea, ld_bound_reject, kdailyavtypes, kqc_cutoff)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE par_oce
    USE dom_oce, ONLY: gdept_1d, nproc
    TYPE(obs_prof), INTENT(INOUT) :: profdata
    TYPE(obs_prof), INTENT(INOUT) :: prodatqc
    LOGICAL, INTENT(IN) :: ld_var1
    LOGICAL, INTENT(IN) :: ld_var2
    LOGICAL, INTENT(IN) :: ld_nea
    LOGICAL, INTENT(IN) :: ld_bound_reject
    INTEGER, INTENT(IN) :: kpi, kpj, kpk
    INTEGER, DIMENSION(imaxavtypes), OPTIONAL :: kdailyavtypes
    REAL(KIND = wp), INTENT(IN), DIMENSION(kpi, kpj, kpk) :: zmask1, zmask2
    REAL(KIND = wp), INTENT(IN), DIMENSION(kpi, kpj) :: pglam1, pglam2, pgphi1, pgphi2
    INTEGER, INTENT(IN), OPTIONAL :: kqc_cutoff
    INTEGER :: iqc_cutoff = 255
    INTEGER :: iyea0
    INTEGER :: imon0
    INTEGER :: iday0
    INTEGER :: ihou0
    INTEGER :: imin0
    INTEGER :: icycle
    INTEGER :: iotdobs
    INTEGER :: iosdv1obs
    INTEGER :: iosdv2obs
    INTEGER :: ilanv1obs
    INTEGER :: ilanv2obs
    INTEGER :: inlav1obs
    INTEGER :: inlav2obs
    INTEGER :: ibdyv1obs
    INTEGER :: ibdyv2obs
    INTEGER :: igrdobs
    INTEGER :: iuvchku
    INTEGER :: iuvchkv
    INTEGER :: iotdobsmpp
    INTEGER :: iosdv1obsmpp
    INTEGER :: iosdv2obsmpp
    INTEGER :: ilanv1obsmpp
    INTEGER :: ilanv2obsmpp
    INTEGER :: inlav1obsmpp
    INTEGER :: inlav2obsmpp
    INTEGER :: ibdyv1obsmpp
    INTEGER :: ibdyv2obsmpp
    INTEGER :: igrdobsmpp
    INTEGER :: iuvchkumpp
    INTEGER :: iuvchkvmpp
    TYPE(obs_prof_valid) :: llvalid
    TYPE(obs_prof_valid), DIMENSION(profdata % nvar) :: llvvalid
    INTEGER :: jvar
    INTEGER :: jobs
    INTEGER :: jstp
    INTEGER :: inrc
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_pre_prof', 'r0', 0, 0)
    IF (lwp) WRITE(numout, FMT = *) 'obs_pre_prof: Preparing the profile data...'
    IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~'
    iyea0 = ndate0 / 10000
    imon0 = (ndate0 - iyea0 * 10000) / 100
    iday0 = ndate0 - iyea0 * 10000 - imon0 * 100
    ihou0 = nn_time0 / 100
    imin0 = (nn_time0 - ihou0 * 100)
    icycle = nn_no
    iotdobs = 0
    igrdobs = 0
    iosdv1obs = 0
    iosdv2obs = 0
    ilanv1obs = 0
    ilanv2obs = 0
    inlav1obs = 0
    inlav2obs = 0
    ibdyv1obs = 0
    ibdyv2obs = 0
    iuvchku = 0
    iuvchkv = 0
    IF (PRESENT(kqc_cutoff)) iqc_cutoff = kqc_cutoff
    IF (PRESENT(kdailyavtypes)) THEN
      CALL obs_coo_tim_prof(icycle, iyea0, imon0, iday0, ihou0, imin0, profdata % nprof, profdata % nyea, profdata % nmon, &
&profdata % nday, profdata % nhou, profdata % nmin, profdata % ntyp, profdata % nqc, profdata % mstp, iotdobs, &
&kdailyavtypes = kdailyavtypes, kqc_cutoff = iqc_cutoff)
    ELSE
      CALL obs_coo_tim_prof(icycle, iyea0, imon0, iday0, ihou0, imin0, profdata % nprof, profdata % nyea, profdata % nmon, &
&profdata % nday, profdata % nhou, profdata % nmin, profdata % ntyp, profdata % nqc, profdata % mstp, iotdobs, &
&kqc_cutoff = iqc_cutoff)
    END IF
    CALL obs_mpp_sum_integer(iotdobs, iotdobsmpp)
    CALL obs_coo_grd(profdata % nprof, profdata % mi(:, 1), profdata % mj(:, 1), profdata % nqc, igrdobs)
    CALL obs_coo_grd(profdata % nprof, profdata % mi(:, 2), profdata % mj(:, 2), profdata % nqc, igrdobs)
    CALL obs_mpp_sum_integer(igrdobs, igrdobsmpp)
    CALL obs_pro_rej(profdata, kqc_cutoff = iqc_cutoff)
    CALL obs_coo_spc_3d(profdata % nprof, profdata % nvprot(1), profdata % npvsta(:, 1), profdata % npvend(:, 1), jpi, jpj, jpk, &
&profdata % mi, profdata % mj, profdata % var(1) % mvk, profdata % rlam, profdata % rphi, profdata % var(1) % vdep, pglam1, &
&pgphi1, gdept_1d, zmask1, profdata % nqc, profdata % var(1) % nvqc, iosdv1obs, ilanv1obs, inlav1obs, ld_nea, ibdyv1obs, &
&ld_bound_reject, iqc_cutoff)
    CALL obs_mpp_sum_integer(iosdv1obs, iosdv1obsmpp)
    CALL obs_mpp_sum_integer(ilanv1obs, ilanv1obsmpp)
    CALL obs_mpp_sum_integer(inlav1obs, inlav1obsmpp)
    CALL obs_mpp_sum_integer(ibdyv1obs, ibdyv1obsmpp)
    CALL obs_coo_spc_3d(profdata % nprof, profdata % nvprot(2), profdata % npvsta(:, 2), profdata % npvend(:, 2), jpi, jpj, jpk, &
&profdata % mi, profdata % mj, profdata % var(2) % mvk, profdata % rlam, profdata % rphi, profdata % var(2) % vdep, pglam2, &
&pgphi2, gdept_1d, zmask2, profdata % nqc, profdata % var(2) % nvqc, iosdv2obs, ilanv2obs, inlav2obs, ld_nea, ibdyv2obs, &
&ld_bound_reject, iqc_cutoff)
    CALL obs_mpp_sum_integer(iosdv2obs, iosdv2obsmpp)
    CALL obs_mpp_sum_integer(ilanv2obs, ilanv2obsmpp)
    CALL obs_mpp_sum_integer(inlav2obs, inlav2obsmpp)
    CALL obs_mpp_sum_integer(ibdyv2obs, ibdyv2obsmpp)
    IF (TRIM(profdata % cvars(1)) == 'UVEL') THEN
      CALL obs_uv_rej(profdata, iuvchku, iuvchkv, iqc_cutoff)
      CALL obs_mpp_sum_integer(iuvchku, iuvchkumpp)
      CALL obs_mpp_sum_integer(iuvchkv, iuvchkvmpp)
    END IF
    ALLOCATE(llvalid % luse(profdata % nprof))
    DO jvar = 1, profdata % nvar
      ALLOCATE(llvvalid(jvar) % luse(profdata % nvprot(jvar)))
    END DO
    llvalid % luse(:) = (profdata % nqc(:) <= iqc_cutoff)
    DO jvar = 1, profdata % nvar
      llvvalid(jvar) % luse(:) = (profdata % var(jvar) % nvqc(:) <= iqc_cutoff)
    END DO
    CALL obs_prof_compress(profdata, prodatqc, .TRUE., numout, lvalid = llvalid, lvvalid = llvvalid)
    DEALLOCATE(llvalid % luse)
    DO jvar = 1, profdata % nvar
      DEALLOCATE(llvvalid(jvar) % luse)
    END DO
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) ' Profiles outside time domain                     = ', iotdobsmpp
      WRITE(numout, FMT = *) ' Remaining profiles that failed grid search       = ', igrdobsmpp
      WRITE(numout, FMT = *) ' Remaining ' // prodatqc % cvars(1) // ' data outside space domain       = ', iosdv1obsmpp
      WRITE(numout, FMT = *) ' Remaining ' // prodatqc % cvars(1) // ' data at land points             = ', ilanv1obsmpp
      IF (ld_nea) THEN
        WRITE(numout, FMT = *) ' Remaining ' // prodatqc % cvars(1) // ' data near land points (removed) = ', inlav1obsmpp
      ELSE
        WRITE(numout, FMT = *) ' Remaining ' // prodatqc % cvars(1) // ' data near land points (kept)    = ', inlav1obsmpp
      END IF
      IF (TRIM(profdata % cvars(1)) == 'UVEL') THEN
        WRITE(numout, FMT = *) ' U observation rejected since V rejected     = ', iuvchku
      END IF
      WRITE(numout, FMT = *) ' Remaining ' // prodatqc % cvars(1) // ' data near open boundary (removed) = ', ibdyv1obsmpp
      WRITE(numout, FMT = *) ' ' // prodatqc % cvars(1) // ' data accepted                             = ', prodatqc % nvprotmpp(1)
      WRITE(numout, FMT = *) ' Remaining ' // prodatqc % cvars(2) // ' data outside space domain       = ', iosdv2obsmpp
      WRITE(numout, FMT = *) ' Remaining ' // prodatqc % cvars(2) // ' data at land points             = ', ilanv2obsmpp
      IF (ld_nea) THEN
        WRITE(numout, FMT = *) ' Remaining ' // prodatqc % cvars(2) // ' data near land points (removed) = ', inlav2obsmpp
      ELSE
        WRITE(numout, FMT = *) ' Remaining ' // prodatqc % cvars(2) // ' data near land points (kept)    = ', inlav2obsmpp
      END IF
      IF (TRIM(profdata % cvars(1)) == 'UVEL') THEN
        WRITE(numout, FMT = *) ' V observation rejected since U rejected     = ', iuvchkv
      END IF
      WRITE(numout, FMT = *) ' Remaining ' // prodatqc % cvars(2) // ' data near open boundary (removed) = ', ibdyv2obsmpp
      WRITE(numout, FMT = *) ' ' // prodatqc % cvars(2) // ' data accepted                             = ', prodatqc % nvprotmpp(2)
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) ' Number of observations per time step :'
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = '(10X,A,5X,A,5X,A,A)') 'Time step', 'Profiles', '     ' // prodatqc % cvars(1) // '     ', '     ' // &
&prodatqc % cvars(2) // '     '
      WRITE(numout, 998)
    END IF
    DO jobs = 1, prodatqc % nprof
      inrc = prodatqc % mstp(jobs) + 2 - nit000
      prodatqc % npstp(inrc) = prodatqc % npstp(inrc) + 1
      DO jvar = 1, prodatqc % nvar
        IF (prodatqc % npvend(jobs, jvar) > 0) THEN
          prodatqc % nvstp(inrc, jvar) = prodatqc % nvstp(inrc, jvar) + (prodatqc % npvend(jobs, jvar) - prodatqc % npvsta(jobs, &
&jvar) + 1)
        END IF
      END DO
    END DO
    CALL obs_mpp_sum_integers(prodatqc % npstp, prodatqc % npstpmpp, nitend - nit000 + 2)
    DO jvar = 1, prodatqc % nvar
      CALL obs_mpp_sum_integers(prodatqc % nvstp(:, jvar), prodatqc % nvstpmpp(:, jvar), nitend - nit000 + 2)
    END DO
    IF (lwp) THEN
      DO jstp = nit000 - 1, nitend
        inrc = jstp - nit000 + 2
        WRITE(numout, 999) jstp, prodatqc % npstpmpp(inrc), prodatqc % nvstpmpp(inrc, 1), prodatqc % nvstpmpp(inrc, 2)
      END DO
    END IF
998 FORMAT(10X, '---------', 5X, '--------', 5X, '-----------', 5X, '----------------')
999 FORMAT(10X, I9, 5X, I8, 5X, I11, 5X, I8)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_pre_prof
  SUBROUTINE obs_coo_tim(kcycle, kyea0, kmon0, kday0, khou0, kmin0, kobsno, kobsyea, kobsmon, kobsday, kobshou, kobsmin, kobsqc, &
&kobsstp, kotdobs)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE dom_oce, ONLY: rdt
    USE phycst, ONLY: rday, rmmss, rhhmm
    INTEGER, INTENT(IN) :: kcycle
    INTEGER, INTENT(IN) :: kyea0
    INTEGER, INTENT(IN) :: kmon0
    INTEGER, INTENT(IN) :: kday0
    INTEGER, INTENT(IN) :: khou0
    INTEGER, INTENT(IN) :: kmin0
    INTEGER, INTENT(IN) :: kobsno
    INTEGER, INTENT(INOUT) :: kotdobs
    INTEGER, DIMENSION(kobsno), INTENT(IN) :: kobsyea, kobsmon, kobsday, kobshou, kobsmin
    INTEGER, DIMENSION(kobsno), INTENT(INOUT) :: kobsqc
    INTEGER, DIMENSION(kobsno), INTENT(OUT) :: kobsstp
    INTEGER :: jyea
    INTEGER :: jmon
    INTEGER :: jday
    INTEGER :: jobs
    INTEGER :: iyeastr
    INTEGER :: iyeaend
    INTEGER :: imonstr
    INTEGER :: imonend
    INTEGER :: idaystr
    INTEGER :: idayend
    INTEGER :: iskip
    INTEGER :: idaystp
    REAL(KIND = wp) :: zminstp
    REAL(KIND = wp) :: zhoustp
    REAL(KIND = wp) :: zobsstp
    INTEGER, DIMENSION(12) :: imonth_len
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    CALL profile_psy_data0 % PreStart('obs_coo_tim', 'r0', 0, 0)
    idaystp = NINT(rday / rdt)
    CALL profile_psy_data0 % PostEnd
    DO jobs = 1, kobsno
      CALL profile_psy_data1 % PreStart('obs_coo_tim', 'r1', 0, 0)
      kobsstp(jobs) = nit000 - 1
      IF ((kobsyea(jobs) < kyea0) .OR. ((kobsyea(jobs) == kyea0) .AND. (kobsmon(jobs) < kmon0)) .OR. ((kobsyea(jobs) == kyea0) &
&.AND. (kobsmon(jobs) == kmon0) .AND. (kobsday(jobs) < kday0)) .OR. ((kobsyea(jobs) == kyea0) .AND. (kobsmon(jobs) == kmon0) .AND. &
&(kobsday(jobs) == kday0) .AND. (kobshou(jobs) < khou0)) .OR. ((kobsyea(jobs) == kyea0) .AND. (kobsmon(jobs) == kmon0) .AND. &
&(kobsday(jobs) == kday0) .AND. (kobshou(jobs) == khou0) .AND. (kobsmin(jobs) <= kmin0))) THEN
        kobsstp(jobs) = - 1
        kobsqc(jobs) = IBSET(kobsqc(jobs), 13)
        kotdobs = kotdobs + 1
        CYCLE
      END IF
      iyeastr = kyea0
      iyeaend = kobsyea(jobs)
      CALL profile_psy_data1 % PostEnd
      DO jyea = iyeastr, iyeaend
        CALL profile_psy_data2 % PreStart('obs_coo_tim', 'r2', 0, 0)
        CALL calc_month_len(jyea, imonth_len)
        imonstr = 1
        IF (jyea == kyea0) imonstr = kmon0
        imonend = 12
        IF (jyea == kobsyea(jobs)) imonend = kobsmon(jobs)
        CALL profile_psy_data2 % PostEnd
        DO jmon = imonstr, imonend
          CALL profile_psy_data3 % PreStart('obs_coo_tim', 'r3', 0, 0)
          idaystr = 1
          IF ((jmon == kmon0) .AND. (jyea == kyea0)) idaystr = kday0
          idayend = imonth_len(jmon)
          IF ((jmon == kobsmon(jobs)) .AND. (jyea == kobsyea(jobs))) idayend = kobsday(jobs) - 1
          CALL profile_psy_data3 % PostEnd
          !$ACC KERNELS
          DO jday = idaystr, idayend
            kobsstp(jobs) = kobsstp(jobs) + idaystp
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      CALL profile_psy_data4 % PreStart('obs_coo_tim', 'r4', 0, 0)
      zminstp = rmmss / rdt
      zhoustp = rhhmm * zminstp
      zobsstp = REAL(kobsmin(jobs) - kmin0, KIND = wp) * zminstp + REAL(kobshou(jobs) - khou0, KIND = wp) * zhoustp
      kobsstp(jobs) = kobsstp(jobs) + NINT(zobsstp)
      IF ((kobsstp(jobs) < (nit000 - 1)) .OR. (kobsstp(jobs) > nitend)) THEN
        kobsqc(jobs) = IBSET(kobsqc(jobs), 13)
        kotdobs = kotdobs + 1
        CYCLE
      END IF
      CALL profile_psy_data4 % PostEnd
    END DO
  END SUBROUTINE obs_coo_tim
  SUBROUTINE calc_month_len(iyear, imonth_len)
    INTEGER, DIMENSION(12) :: imonth_len
    INTEGER :: iyear
    !$ACC KERNELS
    IF (nleapy < 2) THEN
      imonth_len(:) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
      IF (nleapy == 1) THEN
        IF (MOD(iyear, 4) == 0 .AND. (MOD(iyear, 400) == 0 .OR. MOD(iyear, 100) /= 0)) THEN
          imonth_len(2) = 29
        END IF
      END IF
    ELSE
      imonth_len(:) = nleapy
    END IF
    !$ACC END KERNELS
  END SUBROUTINE
  SUBROUTINE obs_coo_tim_prof(kcycle, kyea0, kmon0, kday0, khou0, kmin0, kobsno, kobsyea, kobsmon, kobsday, kobshou, kobsmin, &
&ktyp, kobsqc, kobsstp, kotdobs, kdailyavtypes, kqc_cutoff)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kcycle
    INTEGER, INTENT(IN) :: kyea0
    INTEGER, INTENT(IN) :: kmon0
    INTEGER, INTENT(IN) :: kday0
    INTEGER, INTENT(IN) :: khou0
    INTEGER, INTENT(IN) :: kmin0
    INTEGER, INTENT(IN) :: kobsno
    INTEGER, INTENT(INOUT) :: kotdobs
    INTEGER, DIMENSION(kobsno), INTENT(IN) :: kobsyea, kobsmon, kobsday, kobshou, kobsmin, ktyp
    INTEGER, DIMENSION(kobsno), INTENT(INOUT) :: kobsqc
    INTEGER, DIMENSION(kobsno), INTENT(OUT) :: kobsstp
    INTEGER, DIMENSION(imaxavtypes), OPTIONAL :: kdailyavtypes
    INTEGER, OPTIONAL, INTENT(IN) :: kqc_cutoff
    INTEGER :: jobs
    INTEGER :: iqc_cutoff = 255
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_coo_tim_prof', 'r0', 0, 0)
    CALL obs_coo_tim(kcycle, kyea0, kmon0, kday0, khou0, kmin0, kobsno, kobsyea, kobsmon, kobsday, kobshou, kobsmin, kobsqc, &
&kobsstp, kotdobs)
    IF (PRESENT(kdailyavtypes)) THEN
      DO jobs = 1, kobsno
        IF (kobsqc(jobs) <= iqc_cutoff) THEN
          IF ((kobsstp(jobs) == (nit000 - 1)) .AND. (ANY(kdailyavtypes(:) == ktyp(jobs)))) THEN
            kobsqc(jobs) = IBSET(kobsqc(jobs), 13)
            kotdobs = kotdobs + 1
            CYCLE
          END IF
        END IF
      END DO
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_coo_tim_prof
  SUBROUTINE obs_coo_grd(kobsno, kobsi, kobsj, kobsqc, kgrdobs)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kobsno
    INTEGER, DIMENSION(kobsno), INTENT(IN) :: kobsi, kobsj
    INTEGER, INTENT(INOUT) :: kgrdobs
    INTEGER, DIMENSION(kobsno), INTENT(INOUT) :: kobsqc
    INTEGER :: jobs
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_coo_grd', 'r0', 0, 0)
    DO jobs = 1, kobsno
      IF ((kobsi(jobs) <= 0) .AND. (kobsj(jobs) <= 0)) THEN
        kobsqc(jobs) = IBSET(kobsqc(jobs), 12)
        kgrdobs = kgrdobs + 1
      END IF
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_coo_grd
  SUBROUTINE obs_coo_spc_2d(kobsno, kpi, kpj, kobsi, kobsj, pobslam, pobsphi, plam, pphi, pmask, kobsqc, kosdobs, klanobs, &
&knlaobs, ld_nea, kbdyobs, ld_bound_reject, kqc_cutoff)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kobsno
    INTEGER, INTENT(IN) :: kpi, kpj
    INTEGER, INTENT(IN), DIMENSION(kobsno) :: kobsi, kobsj
    REAL(KIND = wp), INTENT(IN), DIMENSION(kobsno) :: pobslam, pobsphi
    REAL(KIND = wp), INTENT(IN), DIMENSION(kpi, kpj) :: plam, pphi
    REAL(KIND = wp), INTENT(IN), DIMENSION(kpi, kpj) :: pmask
    INTEGER, INTENT(INOUT), DIMENSION(kobsno) :: kobsqc
    INTEGER, INTENT(INOUT) :: kosdobs
    INTEGER, INTENT(INOUT) :: klanobs
    INTEGER, INTENT(INOUT) :: knlaobs
    INTEGER, INTENT(INOUT) :: kbdyobs
    LOGICAL, INTENT(IN) :: ld_nea
    LOGICAL, INTENT(IN) :: ld_bound_reject
    INTEGER, INTENT(IN) :: kqc_cutoff
    REAL(KIND = wp), DIMENSION(2, 2, kobsno) :: zgmsk
    REAL(KIND = wp), DIMENSION(2, 2, kobsno) :: zbmsk
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zbdymask
    REAL(KIND = wp), DIMENSION(2, 2, kobsno) :: zglam, zgphi
    INTEGER, DIMENSION(2, 2, kobsno) :: igrdi, igrdj
    LOGICAL :: lgridobs
    INTEGER :: iig, ijg
    INTEGER :: jobs, ji, jj
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('obs_coo_spc_2d', 'r0', 0, 0)
    DO jobs = 1, kobsno
      IF (kobsqc(jobs) >= kqc_cutoff) THEN
        igrdi(1, 1, jobs) = 1
        igrdj(1, 1, jobs) = 1
        igrdi(1, 2, jobs) = 1
        igrdj(1, 2, jobs) = 2
        igrdi(2, 1, jobs) = 2
        igrdj(2, 1, jobs) = 1
        igrdi(2, 2, jobs) = 2
        igrdj(2, 2, jobs) = 2
      ELSE
        igrdi(1, 1, jobs) = kobsi(jobs) - 1
        igrdj(1, 1, jobs) = kobsj(jobs) - 1
        igrdi(1, 2, jobs) = kobsi(jobs) - 1
        igrdj(1, 2, jobs) = kobsj(jobs)
        igrdi(2, 1, jobs) = kobsi(jobs)
        igrdj(2, 1, jobs) = kobsj(jobs) - 1
        igrdi(2, 2, jobs) = kobsi(jobs)
        igrdj(2, 2, jobs) = kobsj(jobs)
      END IF
    END DO
    CALL profile_psy_data0 % PostEnd
    IF (ln_bdy) THEN
      IF (ld_bound_reject) THEN
        !$ACC KERNELS
        zbdymask(:, :) = 1.0_wp
        !$ACC END KERNELS
        CALL profile_psy_data1 % PreStart('obs_coo_spc_2d', 'r1', 0, 0)
        DO ji = 1, nb_bdy
          DO jj = 1, idx_bdy(ji) % nblen(1)
            zbdymask(idx_bdy(ji) % nbi(jj, 1), idx_bdy(ji) % nbj(jj, 1)) = 0.0_wp
          END DO
        END DO
        CALL obs_int_comm_2d(2, 2, kobsno, kpi, kpj, igrdi, igrdj, zbdymask, zbmsk)
        CALL profile_psy_data1 % PostEnd
      END IF
    END IF
    CALL profile_psy_data2 % PreStart('obs_coo_spc_2d', 'r2', 0, 0)
    CALL obs_int_comm_2d(2, 2, kobsno, kpi, kpj, igrdi, igrdj, pmask, zgmsk)
    CALL obs_int_comm_2d(2, 2, kobsno, kpi, kpj, igrdi, igrdj, plam, zglam)
    CALL obs_int_comm_2d(2, 2, kobsno, kpi, kpj, igrdi, igrdj, pphi, zgphi)
    DO jobs = 1, kobsno
      IF (kobsqc(jobs) >= kqc_cutoff) CYCLE
      IF ((pobslam(jobs) < - 180.) .OR. (pobslam(jobs) > 180.) .OR. (pobsphi(jobs) < - 90.) .OR. (pobsphi(jobs) > 90.)) THEN
        kobsqc(jobs) = IBSET(kobsqc(jobs), 11)
        kosdobs = kosdobs + 1
        CYCLE
      END IF
      IF (SUM(zgmsk(1 : 2, 1 : 2, jobs)) == 0.0_wp) THEN
        kobsqc(jobs) = IBSET(kobsqc(jobs), 10)
        klanobs = klanobs + 1
        CYCLE
      END IF
      lgridobs = .FALSE.
      iig = - 1
      ijg = - 1
      DO jj = 1, 2
        DO ji = 1, 2
          IF ((ABS(zgphi(ji, jj, jobs) - pobsphi(jobs)) < 1.0E-6_wp) .AND. (ABS(MOD(zglam(ji, jj, jobs) - pobslam(jobs), 360.0)) < &
&1.0E-6_wp)) THEN
            lgridobs = .TRUE.
            iig = ji
            ijg = jj
          END IF
        END DO
      END DO
      IF (lgridobs) THEN
        IF (zgmsk(iig, ijg, jobs) == 0.0_wp) THEN
          kobsqc(jobs) = IBSET(kobsqc(jobs), 10)
          klanobs = klanobs + 1
          CYCLE
        END IF
      END IF
      IF (MINVAL(zgmsk(1 : 2, 1 : 2, jobs)) == 0.0_wp) THEN
        knlaobs = knlaobs + 1
        IF (ld_nea) THEN
          kobsqc(jobs) = IBSET(kobsqc(jobs), 9)
          CYCLE
        END IF
      END IF
      IF (ln_bdy) THEN
        IF (ld_bound_reject) THEN
          IF (MINVAL(zbmsk(1 : 2, 1 : 2, jobs)) == 0.0_wp) THEN
            kobsqc(jobs) = IBSET(kobsqc(jobs), 8)
            kbdyobs = kbdyobs + 1
            CYCLE
          END IF
          IF (lgridobs) THEN
            IF (zbmsk(iig, ijg, jobs) == 0.0_wp) THEN
              kobsqc(jobs) = IBSET(kobsqc(jobs), 8)
              kbdyobs = kbdyobs + 1
              CYCLE
            END IF
          END IF
        END IF
      END IF
    END DO
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE obs_coo_spc_2d
  SUBROUTINE obs_coo_spc_3d(kprofno, kobsno, kpstart, kpend, kpi, kpj, kpk, kobsi, kobsj, kobsk, pobslam, pobsphi, pobsdep, plam, &
&pphi, pdep, pmask, kpobsqc, kobsqc, kosdobs, klanobs, knlaobs, ld_nea, kbdyobs, ld_bound_reject, kqc_cutoff)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE dom_oce, ONLY: gdepw_1d, gdepw_0, gdepw_n, gdept_n, ln_zco, ln_zps
    INTEGER, INTENT(IN) :: kprofno
    INTEGER, INTENT(IN) :: kobsno
    INTEGER, INTENT(IN) :: kpi
    INTEGER, INTENT(IN) :: kpj
    INTEGER, INTENT(IN) :: kpk
    INTEGER, DIMENSION(kprofno), INTENT(IN) :: kpstart, kpend
    INTEGER, DIMENSION(kprofno), INTENT(IN) :: kobsi, kobsj
    INTEGER, DIMENSION(kobsno), INTENT(IN) :: kobsk
    REAL(KIND = wp), DIMENSION(kprofno), INTENT(IN) :: pobslam, pobsphi
    REAL(KIND = wp), DIMENSION(kobsno), INTENT(INOUT) :: pobsdep
    REAL(KIND = wp), DIMENSION(kpi, kpj), INTENT(IN) :: plam, pphi
    REAL(KIND = wp), DIMENSION(kpk), INTENT(IN) :: pdep
    REAL(KIND = wp), DIMENSION(kpi, kpj, kpk), INTENT(IN) :: pmask
    INTEGER, DIMENSION(kprofno), INTENT(INOUT) :: kpobsqc
    INTEGER, DIMENSION(kobsno), INTENT(INOUT) :: kobsqc
    INTEGER, INTENT(INOUT) :: kosdobs
    INTEGER, INTENT(INOUT) :: klanobs
    INTEGER, INTENT(INOUT) :: knlaobs
    INTEGER, INTENT(INOUT) :: kbdyobs
    LOGICAL, INTENT(IN) :: ld_nea
    LOGICAL, INTENT(IN) :: ld_bound_reject
    INTEGER, INTENT(IN) :: kqc_cutoff
    REAL(KIND = wp), DIMENSION(2, 2, kpk, kprofno) :: zgmsk
    REAL(KIND = wp), DIMENSION(2, 2, kprofno) :: zbmsk
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zbdymask
    REAL(KIND = wp), DIMENSION(2, 2, kpk, kprofno) :: zgdepw
    REAL(KIND = wp), DIMENSION(2, 2, kprofno) :: zglam, zgphi
    INTEGER, DIMENSION(2, 2, kprofno) :: igrdi, igrdj
    LOGICAL :: lgridobs
    LOGICAL :: ll_next_to_land
    INTEGER :: iig, ijg
    INTEGER :: jobs, jobsp, jk, ji, jj
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('obs_coo_spc_3d', 'r0', 0, 0)
    DO jobs = 1, kprofno
      IF (kpobsqc(jobs) >= kqc_cutoff) THEN
        igrdi(1, 1, jobs) = 1
        igrdj(1, 1, jobs) = 1
        igrdi(1, 2, jobs) = 1
        igrdj(1, 2, jobs) = 2
        igrdi(2, 1, jobs) = 2
        igrdj(2, 1, jobs) = 1
        igrdi(2, 2, jobs) = 2
        igrdj(2, 2, jobs) = 2
      ELSE
        igrdi(1, 1, jobs) = kobsi(jobs) - 1
        igrdj(1, 1, jobs) = kobsj(jobs) - 1
        igrdi(1, 2, jobs) = kobsi(jobs) - 1
        igrdj(1, 2, jobs) = kobsj(jobs)
        igrdi(2, 1, jobs) = kobsi(jobs)
        igrdj(2, 1, jobs) = kobsj(jobs) - 1
        igrdi(2, 2, jobs) = kobsi(jobs)
        igrdj(2, 2, jobs) = kobsj(jobs)
      END IF
    END DO
    CALL profile_psy_data0 % PostEnd
    IF (ln_bdy) THEN
      IF (ld_bound_reject) THEN
        !$ACC KERNELS
        zbdymask(:, :) = 1.0_wp
        !$ACC END KERNELS
        CALL profile_psy_data1 % PreStart('obs_coo_spc_3d', 'r1', 0, 0)
        DO ji = 1, nb_bdy
          DO jj = 1, idx_bdy(ji) % nblen(1)
            zbdymask(idx_bdy(ji) % nbi(jj, 1), idx_bdy(ji) % nbj(jj, 1)) = 0.0_wp
          END DO
        END DO
        CALL profile_psy_data1 % PostEnd
      END IF
      CALL obs_int_comm_2d(2, 2, kprofno, kpi, kpj, igrdi, igrdj, zbdymask, zbmsk)
    END IF
    CALL profile_psy_data2 % PreStart('obs_coo_spc_3d', 'r2', 0, 0)
    CALL obs_int_comm_3d(2, 2, kprofno, kpi, kpj, kpk, igrdi, igrdj, pmask, zgmsk)
    CALL obs_int_comm_2d(2, 2, kprofno, kpi, kpj, igrdi, igrdj, plam, zglam)
    CALL obs_int_comm_2d(2, 2, kprofno, kpi, kpj, igrdi, igrdj, pphi, zgphi)
    CALL obs_int_comm_3d(2, 2, kprofno, kpi, kpj, kpk, igrdi, igrdj, gdepw_n(:, :, :), zgdepw)
    DO jobs = 1, kprofno
      IF (kpobsqc(jobs) >= kqc_cutoff) CYCLE
      lgridobs = .FALSE.
      iig = - 1
      ijg = - 1
      DO jj = 1, 2
        DO ji = 1, 2
          IF ((ABS(zgphi(ji, jj, jobs) - pobsphi(jobs)) < 1.0E-6_wp) .AND. (ABS(MOD(zglam(ji, jj, jobs) - pobslam(jobs), 360.0)) < &
&1.0E-6_wp)) THEN
            lgridobs = .TRUE.
            iig = ji
            ijg = jj
          END IF
        END DO
      END DO
      IF (ANY(zgmsk(1 : 2, 1 : 2, 1, jobs) == 0.0_wp)) THEN
        ll_next_to_land = .TRUE.
      ELSE
        ll_next_to_land = .FALSE.
      END IF
      DO jobsp = kpstart(jobs), kpend(jobs)
        IF ((pobslam(jobs) < - 180.) .OR. (pobslam(jobs) > 180.) .OR. (pobsphi(jobs) < - 90.) .OR. (pobsphi(jobs) > 90.) .OR. &
&(pobsdep(jobsp) < 0.0) .OR. (pobsdep(jobsp) > gdepw_1d(kpk))) THEN
          kobsqc(jobsp) = IBSET(kobsqc(jobsp), 11)
          kosdobs = kosdobs + 1
          CYCLE
        END IF
        IF (ALL(zgdepw(1 : 2, 1 : 2, kpk, jobs) < pobsdep(jobsp)) .OR. (SUM(zgmsk(1 : 2, 1 : 2, kobsk(jobsp) - 1 : kobsk(jobsp), &
&jobs)) == 0.0_wp)) THEN
          kobsqc(jobsp) = IBSET(kobsqc(jobsp), 10)
          klanobs = klanobs + 1
          CYCLE
        END IF
        IF (ll_next_to_land) THEN
          knlaobs = knlaobs + 1
          IF (ld_nea) THEN
            kobsqc(jobsp) = IBSET(kobsqc(jobsp), 10)
          END IF
        END IF
        IF (lgridobs) THEN
          IF (zgmsk(iig, ijg, kobsk(jobsp) - 1, jobs) == 0.0_wp) THEN
            kobsqc(jobsp) = IBSET(kobsqc(jobs), 10)
            klanobs = klanobs + 1
            CYCLE
          END IF
        END IF
        IF (MINVAL(zgmsk(1 : 2, 1 : 2, kobsk(jobsp) - 1 : kobsk(jobsp), jobs)) == 0.0_wp) THEN
          IF (ld_nea) kobsqc(jobsp) = kobsqc(jobsp) + 14
          knlaobs = knlaobs + 1
        END IF
        IF (pobsdep(jobsp) <= pdep(1)) THEN
          pobsdep(jobsp) = pdep(1)
        END IF
        IF (ln_bdy) THEN
          IF (ld_bound_reject) THEN
            IF (MINVAL(zbmsk(1 : 2, 1 : 2, jobs)) == 0.0_wp) THEN
              kobsqc(jobsp) = IBSET(kobsqc(jobs), 8)
              kbdyobs = kbdyobs + 1
              CYCLE
            END IF
            IF (lgridobs) THEN
              IF (zbmsk(iig, ijg, jobs) == 0.0_wp) THEN
                kobsqc(jobsp) = IBSET(kobsqc(jobs), 8)
                kbdyobs = kbdyobs + 1
                CYCLE
              END IF
            END IF
          END IF
        END IF
      END DO
    END DO
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE obs_coo_spc_3d
  SUBROUTINE obs_pro_rej(profdata, kqc_cutoff)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(obs_prof), INTENT(INOUT) :: profdata
    INTEGER, INTENT(IN) :: kqc_cutoff
    INTEGER :: jprof
    INTEGER :: jvar
    INTEGER :: jobs
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_pro_rej', 'r0', 0, 0)
    DO jprof = 1, profdata % nprof
      IF (profdata % nqc(jprof) > kqc_cutoff) THEN
        DO jvar = 1, profdata % nvar
          DO jobs = profdata % npvsta(jprof, jvar), profdata % npvend(jprof, jvar)
            profdata % var(jvar) % nvqc(jobs) = IBSET(profdata % var(jvar) % nvqc(jobs), 14)
          END DO
        END DO
      END IF
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_pro_rej
  SUBROUTINE obs_uv_rej(profdata, knumu, knumv, kqc_cutoff)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(obs_prof), INTENT(INOUT) :: profdata
    INTEGER, INTENT(INOUT) :: knumu
    INTEGER, INTENT(INOUT) :: knumv
    INTEGER, INTENT(IN) :: kqc_cutoff
    INTEGER :: jprof
    INTEGER :: jvar
    INTEGER :: jobs
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    DO jprof = 1, profdata % nprof
      IF ((profdata % npvsta(jprof, 1) /= profdata % npvsta(jprof, 2)) .OR. (profdata % npvend(jprof, 1) /= profdata % &
&npvend(jprof, 2))) THEN
        CALL ctl_stop('U,V profiles inconsistent in obs_uv_rej')
        RETURN
      END IF
      CALL profile_psy_data0 % PreStart('obs_uv_rej', 'r0', 0, 0)
      DO jobs = profdata % npvsta(jprof, 1), profdata % npvend(jprof, 1)
        IF ((profdata % var(1) % nvqc(jobs) > kqc_cutoff) .AND. (profdata % var(2) % nvqc(jobs) <= kqc_cutoff)) THEN
          profdata % var(2) % nvqc(jobs) = IBSET(profdata % var(1) % nvqc(jobs), 15)
          knumv = knumv + 1
        END IF
        IF ((profdata % var(2) % nvqc(jobs) > kqc_cutoff) .AND. (profdata % var(1) % nvqc(jobs) <= kqc_cutoff)) THEN
          profdata % var(1) % nvqc(jobs) = IBSET(profdata % var(1) % nvqc(jobs), 15)
          knumu = knumu + 1
        END IF
      END DO
      CALL profile_psy_data0 % PostEnd
    END DO
  END SUBROUTINE obs_uv_rej
END MODULE obs_prep