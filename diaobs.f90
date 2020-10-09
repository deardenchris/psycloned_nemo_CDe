MODULE diaobs
  USE par_kind
  USE in_out_manager
  USE par_oce
  USE dom_oce
  USE sbc_oce
  USE obs_read_prof
  USE obs_read_surf
  USE obs_sstbias
  USE obs_readmdt
  USE obs_prep
  USE obs_oper
  USE obs_write
  USE obs_grid
  USE obs_read_altbias
  USE obs_profiles_def
  USE obs_surf_def
  USE obs_types
  USE mpp_map
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dia_obs_init
  PUBLIC :: dia_obs
  PUBLIC :: dia_obs_wri
  PUBLIC :: dia_obs_dealloc
  PUBLIC :: calc_date
  LOGICAL, PUBLIC :: ln_diaobs
  LOGICAL :: ln_sstnight
  LOGICAL :: ln_sla_fp_indegs
  LOGICAL :: ln_sst_fp_indegs
  LOGICAL :: ln_sss_fp_indegs
  LOGICAL :: ln_sic_fp_indegs
  REAL(KIND = wp) :: rn_sla_avglamscl
  REAL(KIND = wp) :: rn_sla_avgphiscl
  REAL(KIND = wp) :: rn_sst_avglamscl
  REAL(KIND = wp) :: rn_sst_avgphiscl
  REAL(KIND = wp) :: rn_sss_avglamscl
  REAL(KIND = wp) :: rn_sss_avgphiscl
  REAL(KIND = wp) :: rn_sic_avglamscl
  REAL(KIND = wp) :: rn_sic_avgphiscl
  INTEGER :: nn_1dint
  INTEGER :: nn_2dint
  INTEGER :: nn_2dint_sla
  INTEGER :: nn_2dint_sst
  INTEGER :: nn_2dint_sss
  INTEGER :: nn_2dint_sic
  INTEGER, DIMENSION(imaxavtypes) :: nn_profdavtypes
  INTEGER :: nproftypes
  INTEGER :: nsurftypes
  INTEGER, DIMENSION(:), ALLOCATABLE :: nvarsprof, nvarssurf
  INTEGER, DIMENSION(:), ALLOCATABLE :: nextrprof, nextrsurf
  INTEGER, DIMENSION(:), ALLOCATABLE :: n2dintsurf
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: zavglamscl, zavgphiscl
  LOGICAL, DIMENSION(:), ALLOCATABLE :: lfpindegs
  LOGICAL, DIMENSION(:), ALLOCATABLE :: llnightav
  TYPE(obs_surf), PUBLIC, POINTER, DIMENSION(:) :: surfdata
  TYPE(obs_surf), PUBLIC, POINTER, DIMENSION(:) :: surfdataqc
  TYPE(obs_prof), PUBLIC, POINTER, DIMENSION(:) :: profdata
  TYPE(obs_prof), PUBLIC, POINTER, DIMENSION(:) :: profdataqc
  CHARACTER(LEN = 6), PUBLIC, DIMENSION(:), ALLOCATABLE :: cobstypesprof, cobstypessurf
  CONTAINS
  SUBROUTINE dia_obs_init
    INTEGER, PARAMETER :: jpmaxnfiles = 1000
    INTEGER, DIMENSION(:), ALLOCATABLE :: ifilesprof, ifilessurf
    INTEGER :: ios
    INTEGER :: jtype
    INTEGER :: jvar
    INTEGER :: jfile
    INTEGER :: jnumsstbias
    CHARACTER(LEN = 128), DIMENSION(jpmaxnfiles) :: cn_profbfiles, cn_sstfbfiles, cn_sssfbfiles, cn_slafbfiles, cn_sicfbfiles, &
&cn_velfbfiles, cn_sstbiasfiles
    CHARACTER(LEN = 128) :: cn_altbiasfile
    CHARACTER(LEN = 128), DIMENSION(:, :), ALLOCATABLE :: clproffiles, clsurffiles
    LOGICAL :: ln_t3d
    LOGICAL :: ln_s3d
    LOGICAL :: ln_sla
    LOGICAL :: ln_sst
    LOGICAL :: ln_sss
    LOGICAL :: ln_sic
    LOGICAL :: ln_vel3d
    LOGICAL :: ln_nea
    LOGICAL :: ln_altbias
    LOGICAL :: ln_sstbias
    LOGICAL :: ln_ignmis
    LOGICAL :: ln_s_at_t
    LOGICAL :: ln_bound_reject
    LOGICAL :: llvar1
    LOGICAL :: llvar2
    LOGICAL, DIMENSION(jpmaxnfiles) :: lmask
    REAL(KIND = dp) :: rn_dobsini
    REAL(KIND = dp) :: rn_dobsend
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zglam1, zglam2
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zgphi1, zgphi2
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zmask1, zmask2
    NAMELIST /namobs/ ln_diaobs, ln_t3d, ln_s3d, ln_sla, ln_sst, ln_sic, ln_sss, ln_vel3d, ln_altbias, ln_sstbias, ln_nea, &
&ln_grid_global, ln_grid_search_lookup, ln_ignmis, ln_s_at_t, ln_bound_reject, ln_sstnight, ln_sla_fp_indegs, ln_sst_fp_indegs, &
&ln_sss_fp_indegs, ln_sic_fp_indegs, cn_profbfiles, cn_slafbfiles, cn_sstfbfiles, cn_sicfbfiles, cn_velfbfiles, cn_sssfbfiles, &
&cn_sstbiasfiles, cn_altbiasfile, cn_gridsearchfile, rn_gridsearchres, rn_dobsini, rn_dobsend, rn_sla_avglamscl, rn_sla_avgphiscl, &
&rn_sst_avglamscl, rn_sst_avgphiscl, rn_sss_avglamscl, rn_sss_avgphiscl, rn_sic_avglamscl, rn_sic_avgphiscl, nn_1dint, nn_2dint, &
&nn_2dint_sla, nn_2dint_sst, nn_2dint_sss, nn_2dint_sic, nn_msshc, rn_mdtcorr, rn_mdtcutoff, nn_profdavtypes
    cn_profbfiles(:) = ''
    cn_slafbfiles(:) = ''
    cn_sstfbfiles(:) = ''
    cn_sicfbfiles(:) = ''
    cn_velfbfiles(:) = ''
    cn_sssfbfiles(:) = ''
    cn_sstbiasfiles(:) = ''
    nn_profdavtypes(:) = - 1
    CALL ini_date(rn_dobsini)
    CALL fin_date(rn_dobsend)
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namobs, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namobs in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namobs, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namobs in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namobs)
    IF (.NOT. ln_diaobs) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'dia_obs_init : NO Observation diagnostic used'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      RETURN
    END IF
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dia_obs_init : Observation diagnostic initialization'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namobs : set observation diagnostic parameters'
      WRITE(numout, FMT = *) '      Logical switch for T profile observations                ln_t3d = ', ln_t3d
      WRITE(numout, FMT = *) '      Logical switch for S profile observations                ln_s3d = ', ln_s3d
      WRITE(numout, FMT = *) '      Logical switch for SLA observations                      ln_sla = ', ln_sla
      WRITE(numout, FMT = *) '      Logical switch for SST observations                      ln_sst = ', ln_sst
      WRITE(numout, FMT = *) '      Logical switch for Sea Ice observations                  ln_sic = ', ln_sic
      WRITE(numout, FMT = *) '      Logical switch for velocity observations               ln_vel3d = ', ln_vel3d
      WRITE(numout, FMT = *) '      Logical switch for SSS observations                      ln_sss = ', ln_sss
      WRITE(numout, FMT = *) '      Global distribution of observations              ln_grid_global = ', ln_grid_global
      WRITE(numout, FMT = *) '      Logical switch for obs grid search lookup ln_grid_search_lookup = ', ln_grid_search_lookup
      IF (ln_grid_search_lookup) WRITE(numout, FMT = *) '      Grid search lookup file header                cn_gridsearchfile = &
&', cn_gridsearchfile
      WRITE(numout, FMT = *) '      Initial date in window YYYYMMDD.HHMMSS               rn_dobsini = ', rn_dobsini
      WRITE(numout, FMT = *) '      Final date in window YYYYMMDD.HHMMSS                 rn_dobsend = ', rn_dobsend
      WRITE(numout, FMT = *) '      Type of vertical interpolation method                  nn_1dint = ', nn_1dint
      WRITE(numout, FMT = *) '      Type of horizontal interpolation method                nn_2dint = ', nn_2dint
      WRITE(numout, FMT = *) '      Rejection of observations near land switch               ln_nea = ', ln_nea
      WRITE(numout, FMT = *) '      Rejection of obs near open bdys                 ln_bound_reject = ', ln_bound_reject
      WRITE(numout, FMT = *) '      MSSH correction scheme                                 nn_msshc = ', nn_msshc
      WRITE(numout, FMT = *) '      MDT  correction                                      rn_mdtcorr = ', rn_mdtcorr
      WRITE(numout, FMT = *) '      MDT cutoff for computed correction                 rn_mdtcutoff = ', rn_mdtcutoff
      WRITE(numout, FMT = *) '      Logical switch for alt bias                          ln_altbias = ', ln_altbias
      WRITE(numout, FMT = *) '      Logical switch for sst bias                          ln_sstbias = ', ln_sstbias
      WRITE(numout, FMT = *) '      Logical switch for ignoring missing files             ln_ignmis = ', ln_ignmis
      WRITE(numout, FMT = *) '      Daily average types                             nn_profdavtypes = ', nn_profdavtypes
      WRITE(numout, FMT = *) '      Logical switch for night-time SST obs               ln_sstnight = ', ln_sstnight
    END IF
    nproftypes = COUNT((/ln_t3d .OR. ln_s3d, ln_vel3d/))
    nsurftypes = COUNT((/ln_sla, ln_sst, ln_sic, ln_sss/))
    IF (ln_sstbias) THEN
      lmask(:) = .FALSE.
      WHERE (cn_sstbiasfiles(:) /= '') lmask(:) = .TRUE.
      jnumsstbias = COUNT(lmask)
      lmask(:) = .FALSE.
    END IF
    IF (nproftypes == 0 .AND. nsurftypes == 0) THEN
      CALL ctl_warn('dia_obs_init: ln_diaobs is set to true, but all obs operator logical flags', ' (ln_t3d, ln_s3d, ln_sla, &
&ln_sst, ln_sic, ln_vel3d)', ' are set to .FALSE. so turning off calls to dia_obs')
      ln_diaobs = .FALSE.
      RETURN
    END IF
    IF (nproftypes > 0) THEN
      ALLOCATE(cobstypesprof(nproftypes))
      ALLOCATE(ifilesprof(nproftypes))
      ALLOCATE(clproffiles(nproftypes, jpmaxnfiles))
      jtype = 0
      IF (ln_t3d .OR. ln_s3d) THEN
        jtype = jtype + 1
        CALL obs_settypefiles(nproftypes, jpmaxnfiles, jtype, 'prof  ', cn_profbfiles, ifilesprof, cobstypesprof, clproffiles)
      END IF
      IF (ln_vel3d) THEN
        jtype = jtype + 1
        CALL obs_settypefiles(nproftypes, jpmaxnfiles, jtype, 'vel   ', cn_velfbfiles, ifilesprof, cobstypesprof, clproffiles)
      END IF
    END IF
    IF (nsurftypes > 0) THEN
      ALLOCATE(cobstypessurf(nsurftypes))
      ALLOCATE(ifilessurf(nsurftypes))
      ALLOCATE(clsurffiles(nsurftypes, jpmaxnfiles))
      ALLOCATE(n2dintsurf(nsurftypes))
      ALLOCATE(zavglamscl(nsurftypes))
      ALLOCATE(zavgphiscl(nsurftypes))
      ALLOCATE(lfpindegs(nsurftypes))
      ALLOCATE(llnightav(nsurftypes))
      jtype = 0
      IF (ln_sla) THEN
        jtype = jtype + 1
        CALL obs_settypefiles(nsurftypes, jpmaxnfiles, jtype, 'sla   ', cn_slafbfiles, ifilessurf, cobstypessurf, clsurffiles)
        CALL obs_setinterpopts(nsurftypes, jtype, 'sla   ', nn_2dint, nn_2dint_sla, rn_sla_avglamscl, rn_sla_avgphiscl, &
&ln_sla_fp_indegs, .FALSE., n2dintsurf, zavglamscl, zavgphiscl, lfpindegs, llnightav)
      END IF
      IF (ln_sst) THEN
        jtype = jtype + 1
        CALL obs_settypefiles(nsurftypes, jpmaxnfiles, jtype, 'sst   ', cn_sstfbfiles, ifilessurf, cobstypessurf, clsurffiles)
        CALL obs_setinterpopts(nsurftypes, jtype, 'sst   ', nn_2dint, nn_2dint_sst, rn_sst_avglamscl, rn_sst_avgphiscl, &
&ln_sst_fp_indegs, ln_sstnight, n2dintsurf, zavglamscl, zavgphiscl, lfpindegs, llnightav)
      END IF
      IF (ln_sss) THEN
        jtype = jtype + 1
        CALL obs_settypefiles(nsurftypes, jpmaxnfiles, jtype, 'sss   ', cn_sssfbfiles, ifilessurf, cobstypessurf, clsurffiles)
        CALL obs_setinterpopts(nsurftypes, jtype, 'sss   ', nn_2dint, nn_2dint_sss, rn_sss_avglamscl, rn_sss_avgphiscl, &
&ln_sss_fp_indegs, .FALSE., n2dintsurf, zavglamscl, zavgphiscl, lfpindegs, llnightav)
      END IF
    END IF
    IF (ln_vel3d .AND. .NOT. ln_grid_global) THEN
      CALL ctl_stop('Velocity data only works with ln_grid_global=.true.')
      RETURN
    END IF
    IF (ln_grid_global) THEN
      CALL ctl_warn('dia_obs_init: ln_grid_global=T may cause memory issues when used with a large number of processors')
    END IF
    IF (nn_1dint < 0 .OR. nn_1dint > 1) THEN
      CALL ctl_stop('dia_obs_init: Choice of vertical (1D) interpolation method is not available')
    END IF
    IF (nn_2dint < 0 .OR. nn_2dint > 6) THEN
      CALL ctl_stop('dia_obs_init: Choice of horizontal (2D) interpolation method is not available')
    END IF
    CALL obs_typ_init
    IF (ln_grid_global) CALL mppmap_init
    CALL obs_grid_setup
    IF (nproftypes > 0) THEN
      ALLOCATE(profdata(nproftypes), nvarsprof(nproftypes))
      ALLOCATE(profdataqc(nproftypes), nextrprof(nproftypes))
      DO jtype = 1, nproftypes
        nvarsprof(jtype) = 2
        IF (TRIM(cobstypesprof(jtype)) == 'prof') THEN
          nextrprof(jtype) = 1
          llvar1 = ln_t3d
          llvar2 = ln_s3d
          zglam1 = glamt
          zgphi1 = gphit
          zmask1 = tmask
          zglam2 = glamt
          zgphi2 = gphit
          zmask2 = tmask
        END IF
        IF (TRIM(cobstypesprof(jtype)) == 'vel') THEN
          nextrprof(jtype) = 2
          llvar1 = ln_vel3d
          llvar2 = ln_vel3d
          zglam1 = glamu
          zgphi1 = gphiu
          zmask1 = umask
          zglam2 = glamv
          zgphi2 = gphiv
          zmask2 = vmask
        END IF
        CALL obs_rea_prof(profdata(jtype), ifilesprof(jtype), clproffiles(jtype, 1 : ifilesprof(jtype)), nvarsprof(jtype), &
&nextrprof(jtype), nitend - nit000 + 2, rn_dobsini, rn_dobsend, llvar1, llvar2, ln_ignmis, ln_s_at_t, .FALSE., &
&kdailyavtypes = nn_profdavtypes)
        DO jvar = 1, nvarsprof(jtype)
          CALL obs_prof_staend(profdata(jtype), jvar)
        END DO
        CALL obs_pre_prof(profdata(jtype), profdataqc(jtype), llvar1, llvar2, jpi, jpj, jpk, zmask1, zglam1, zgphi1, zmask2, &
&zglam2, zgphi2, ln_nea, ln_bound_reject, kdailyavtypes = nn_profdavtypes)
      END DO
      DEALLOCATE(ifilesprof, clproffiles)
    END IF
    IF (nsurftypes > 0) THEN
      ALLOCATE(surfdata(nsurftypes), nvarssurf(nsurftypes))
      ALLOCATE(surfdataqc(nsurftypes), nextrsurf(nsurftypes))
      DO jtype = 1, nsurftypes
        nvarssurf(jtype) = 1
        nextrsurf(jtype) = 0
        llnightav(jtype) = .FALSE.
        IF (TRIM(cobstypessurf(jtype)) == 'sla') nextrsurf(jtype) = 2
        IF (TRIM(cobstypessurf(jtype)) == 'sst') llnightav(jtype) = ln_sstnight
        CALL obs_rea_surf(surfdata(jtype), ifilessurf(jtype), clsurffiles(jtype, 1 : ifilessurf(jtype)), nvarssurf(jtype), &
&nextrsurf(jtype), nitend - nit000 + 2, rn_dobsini, rn_dobsend, ln_ignmis, .FALSE., llnightav(jtype))
        CALL obs_pre_surf(surfdata(jtype), surfdataqc(jtype), ln_nea, ln_bound_reject)
        IF (TRIM(cobstypessurf(jtype)) == 'sla') THEN
          CALL obs_rea_mdt(surfdataqc(jtype), n2dintsurf(jtype))
          IF (ln_altbias) CALL obs_rea_altbias(surfdataqc(jtype), n2dintsurf(jtype), cn_altbiasfile)
        END IF
        IF (TRIM(cobstypessurf(jtype)) == 'sst' .AND. ln_sstbias) THEN
          jnumsstbias = 0
          DO jfile = 1, jpmaxnfiles
            IF (TRIM(cn_sstbiasfiles(jfile)) /= '') jnumsstbias = jnumsstbias + 1
          END DO
          IF (jnumsstbias == 0) CALL ctl_stop("ln_sstbias set but no bias files to read in")
          CALL obs_app_sstbias(surfdataqc(jtype), n2dintsurf(jtype), jnumsstbias, cn_sstbiasfiles(1 : jnumsstbias))
        END IF
      END DO
      DEALLOCATE(ifilessurf, clsurffiles)
    END IF
  END SUBROUTINE dia_obs_init
  SUBROUTINE dia_obs(kstp)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE dom_oce, ONLY: gdept_n, gdept_1d
    USE phycst, ONLY: rday
    USE oce, ONLY: tsn, un, vn, sshn
    USE phycst, ONLY: rday
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kstp
    INTEGER :: idaystp
    INTEGER :: jtype
    INTEGER :: jvar
    INTEGER :: ji, jj
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zprofvar1, zprofvar2
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zprofmask1, zprofmask2
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zsurfvar, zsurfmask
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zglam1, zglam2, zgphi1, zgphi2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('dia_obs', 'r0', 0, 0)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dia_obs : Call the observation operators', kstp
      WRITE(numout, FMT = *) '~~~~~~~'
    END IF
    idaystp = NINT(rday / rdt)
    IF (nproftypes > 0) THEN
      DO jtype = 1, nproftypes
        SELECT CASE (TRIM(cobstypesprof(jtype)))
        CASE ('prof')
          zprofvar1(:, :, :) = tsn(:, :, :, jp_tem)
          zprofvar2(:, :, :) = tsn(:, :, :, jp_sal)
          zprofmask1(:, :, :) = tmask(:, :, :)
          zprofmask2(:, :, :) = tmask(:, :, :)
          zglam1(:, :) = glamt(:, :)
          zglam2(:, :) = glamt(:, :)
          zgphi1(:, :) = gphit(:, :)
          zgphi2(:, :) = gphit(:, :)
        CASE ('vel')
          zprofvar1(:, :, :) = un(:, :, :)
          zprofvar2(:, :, :) = vn(:, :, :)
          zprofmask1(:, :, :) = umask(:, :, :)
          zprofmask2(:, :, :) = vmask(:, :, :)
          zglam1(:, :) = glamu(:, :)
          zglam2(:, :) = glamv(:, :)
          zgphi1(:, :) = gphiu(:, :)
          zgphi2(:, :) = gphiv(:, :)
        CASE DEFAULT
          CALL ctl_stop('Unknown profile observation type ' // TRIM(cobstypesprof(jtype)) // ' in dia_obs')
        END SELECT
        CALL obs_prof_opt(profdataqc(jtype), kstp, jpi, jpj, jpk, nit000, idaystp, zprofvar1, zprofvar2, gdept_n(:, :, :), &
&gdepw_n(:, :, :), zprofmask1, zprofmask2, zglam1, zglam2, zgphi1, zgphi2, nn_1dint, nn_2dint, kdailyavtypes = nn_profdavtypes)
      END DO
    END IF
    IF (nsurftypes > 0) THEN
      DO jtype = 1, nsurftypes
        zsurfmask(:, :) = tmask(:, :, 1)
        SELECT CASE (TRIM(cobstypessurf(jtype)))
        CASE ('sst')
          zsurfvar(:, :) = tsn(:, :, 1, jp_tem)
        CASE ('sla')
          zsurfvar(:, :) = sshn(:, :)
        CASE ('sss')
          zsurfvar(:, :) = tsn(:, :, 1, jp_sal)
        CASE ('sic')
          IF (kstp == 0) THEN
            IF (lwp .AND. surfdataqc(jtype) % nsstpmpp(1) > 0) THEN
              CALL ctl_warn('Sea-ice not initialised on zeroth ' // 'time-step but some obs are valid then.')
              WRITE(numout, FMT = *) surfdataqc(jtype) % nsstpmpp(1), ' sea-ice obs will be missed'
            END IF
            surfdataqc(jtype) % nsurfup = surfdataqc(jtype) % nsurfup + surfdataqc(jtype) % nsstp(1)
            CYCLE
          ELSE
            CALL ctl_stop(' Trying to run sea-ice observation operator', ' but no sea-ice model appears to have been defined')
          END IF
        END SELECT
        CALL obs_surf_opt(surfdataqc(jtype), kstp, jpi, jpj, nit000, idaystp, zsurfvar, zsurfmask, n2dintsurf(jtype), &
&llnightav(jtype), zavglamscl(jtype), zavgphiscl(jtype), lfpindegs(jtype))
      END DO
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE dia_obs
  SUBROUTINE dia_obs_wri
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE obs_rot_vel
    IMPLICIT NONE
    INTEGER :: jtype
    INTEGER :: jo, jvar, jk
    REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: zu, zv
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('dia_obs_wri', 'r0', 0, 0)
    IF (nproftypes > 0) THEN
      DO jtype = 1, nproftypes
        IF (TRIM(cobstypesprof(jtype)) == 'vel') THEN
          ALLOCATE(zu(profdataqc(jtype) % nvprot(1)), zv(profdataqc(jtype) % nvprot(2)))
          CALL obs_rotvel(profdataqc(jtype), nn_2dint, zu, zv)
          DO jo = 1, profdataqc(jtype) % nprof
            DO jvar = 1, 2
              DO jk = profdataqc(jtype) % npvsta(jo, jvar), profdataqc(jtype) % npvend(jo, jvar)
                IF (jvar == 1) THEN
                  profdataqc(jtype) % var(jvar) % vmod(jk) = zu(jk)
                ELSE
                  profdataqc(jtype) % var(jvar) % vmod(jk) = zv(jk)
                END IF
              END DO
            END DO
          END DO
          DEALLOCATE(zu)
          DEALLOCATE(zv)
        END IF
        CALL obs_prof_decompress(profdataqc(jtype), profdata(jtype), .TRUE., numout)
        CALL obs_wri_prof(profdata(jtype))
      END DO
    END IF
    IF (nsurftypes > 0) THEN
      DO jtype = 1, nsurftypes
        CALL obs_surf_decompress(surfdataqc(jtype), surfdata(jtype), .TRUE., numout)
        CALL obs_wri_surf(surfdata(jtype))
      END DO
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE dia_obs_wri
  SUBROUTINE dia_obs_dealloc
    IMPLICIT NONE
    CALL obs_grid_deallocate
    IF (nproftypes > 0) DEALLOCATE(cobstypesprof, profdata, profdataqc, nvarsprof, nextrprof)
    IF (nsurftypes > 0) DEALLOCATE(cobstypessurf, surfdata, surfdataqc, nvarssurf, nextrsurf, n2dintsurf, zavglamscl, zavgphiscl, &
&lfpindegs, llnightav)
  END SUBROUTINE dia_obs_dealloc
  SUBROUTINE calc_date(kstp, ddobs)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE phycst, ONLY: rday
    USE dom_oce, ONLY: rdt
    IMPLICIT NONE
    REAL(KIND = dp), INTENT(OUT) :: ddobs
    INTEGER :: kstp
    INTEGER :: iyea
    INTEGER :: imon
    INTEGER :: iday
    INTEGER :: ihou
    INTEGER :: imin
    INTEGER :: imday
    REAL(KIND = wp) :: zdayfrc
    INTEGER, DIMENSION(12) :: imonth_len
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('calc_date', 'r0', 0, 0)
    iyea = ndate0 / 10000
    imon = (ndate0 - iyea * 10000) / 100
    iday = ndate0 - iyea * 10000 - imon * 100
    ihou = nn_time0 / 100
    imin = (nn_time0 - ihou * 100)
    zdayfrc = kstp * rdt / rday
    zdayfrc = zdayfrc - AINT(zdayfrc)
    imin = imin + INT(zdayfrc * 24 * 60)
    DO WHILE (imin >= 60)
      imin = imin - 60
      ihou = ihou + 1
    END DO
    DO WHILE (ihou >= 24)
      ihou = ihou - 24
      iday = iday + 1
    END DO
    iday = iday + kstp * rdt / rday
    CALL calc_month_len(iyea, imonth_len)
    DO WHILE (iday > imonth_len(imon))
      iday = iday - imonth_len(imon)
      imon = imon + 1
      IF (imon > 12) THEN
        imon = 1
        iyea = iyea + 1
        CALL calc_month_len(iyea, imonth_len)
      END IF
    END DO
    ddobs = iyea * 10000_dp + imon * 100_dp + iday + ihou * 0.01_dp + imin * 0.0001_dp
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE calc_date
  SUBROUTINE ini_date(ddobsini)
    IMPLICIT NONE
    REAL(KIND = dp), INTENT(OUT) :: ddobsini
    CALL calc_date(nit000 - 1, ddobsini)
  END SUBROUTINE ini_date
  SUBROUTINE fin_date(ddobsfin)
    IMPLICIT NONE
    REAL(KIND = dp), INTENT(OUT) :: ddobsfin
    CALL calc_date(nitend, ddobsfin)
  END SUBROUTINE fin_date
  SUBROUTINE obs_settypefiles(ntypes, jpmaxnfiles, jtype, ctypein, cfilestype, ifiles, cobstypes, cfiles)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: ntypes
    INTEGER, INTENT(IN) :: jpmaxnfiles
    INTEGER, INTENT(IN) :: jtype
    INTEGER, DIMENSION(ntypes), INTENT(INOUT) :: ifiles
    CHARACTER(LEN = 6), INTENT(IN) :: ctypein
    CHARACTER(LEN = 128), DIMENSION(jpmaxnfiles), INTENT(IN) :: cfilestype
    CHARACTER(LEN = 6), DIMENSION(ntypes), INTENT(INOUT) :: cobstypes
    CHARACTER(LEN = 128), DIMENSION(ntypes, jpmaxnfiles), INTENT(INOUT) :: cfiles
    INTEGER :: jfile
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_settypefiles', 'r0', 0, 0)
    cfiles(jtype, :) = cfilestype(:)
    cobstypes(jtype) = ctypein
    ifiles(jtype) = 0
    DO jfile = 1, jpmaxnfiles
      IF (TRIM(cfiles(jtype, jfile)) /= '') ifiles(jtype) = ifiles(jtype) + 1
    END DO
    IF (ifiles(jtype) == 0) THEN
      CALL ctl_stop('Logical for observation type ' // TRIM(ctypein) // ' set to true but no files available to read')
    END IF
    IF (lwp) THEN
      WRITE(numout, FMT = *) '             ' // cobstypes(jtype) // ' input observation file names:'
      DO jfile = 1, ifiles(jtype)
        WRITE(numout, FMT = *) '                ' // TRIM(cfiles(jtype, jfile))
      END DO
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_settypefiles
  SUBROUTINE obs_setinterpopts(ntypes, jtype, ctypein, n2dint_default, n2dint_type, zavglamscl_type, zavgphiscl_type, &
&lfp_indegs_type, lavnight_type, n2dint, zavglamscl, zavgphiscl, lfpindegs, lavnight)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: ntypes
    INTEGER, INTENT(IN) :: jtype
    INTEGER, INTENT(IN) :: n2dint_default
    INTEGER, INTENT(IN) :: n2dint_type
    REAL(KIND = wp), INTENT(IN) :: zavglamscl_type, zavgphiscl_type
    LOGICAL, INTENT(IN) :: lfp_indegs_type
    LOGICAL, INTENT(IN) :: lavnight_type
    CHARACTER(LEN = 6), INTENT(IN) :: ctypein
    INTEGER, DIMENSION(ntypes), INTENT(INOUT) :: n2dint
    REAL(KIND = wp), DIMENSION(ntypes), INTENT(INOUT) :: zavglamscl, zavgphiscl
    LOGICAL, DIMENSION(ntypes), INTENT(INOUT) :: lfpindegs, lavnight
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_setinterpopts', 'r0', 0, 0)
    lavnight(jtype) = lavnight_type
    IF ((n2dint_type >= 1) .AND. (n2dint_type <= 6)) THEN
      n2dint(jtype) = n2dint_type
    ELSE
      n2dint(jtype) = n2dint_default
    END IF
    IF ((n2dint(jtype) > 4) .AND. (n2dint(jtype) <= 6)) THEN
      IF (zavglamscl_type > 0._wp) THEN
        zavglamscl(jtype) = zavglamscl_type
      ELSE
        CALL ctl_stop('Incorrect value set for averaging footprint ' // 'scale (zavglamscl) for observation type ' // TRIM(ctypein))
      END IF
      IF (zavgphiscl_type > 0._wp) THEN
        zavgphiscl(jtype) = zavgphiscl_type
      ELSE
        CALL ctl_stop('Incorrect value set for averaging footprint ' // 'scale (zavgphiscl) for observation type ' // TRIM(ctypein))
      END IF
      lfpindegs(jtype) = lfp_indegs_type
    END IF
    IF (lwp) THEN
      IF (n2dint(jtype) <= 4) THEN
        WRITE(numout, FMT = *) '             ' // TRIM(ctypein) // ' model counterparts will be interpolated horizontally'
      ELSE IF (n2dint(jtype) <= 6) THEN
        WRITE(numout, FMT = *) '             ' // TRIM(ctypein) // ' model counterparts will be averaged horizontally'
        WRITE(numout, FMT = *) '             ' // '    with E/W scale: ', zavglamscl(jtype)
        WRITE(numout, FMT = *) '             ' // '    with N/S scale: ', zavgphiscl(jtype)
        IF (lfpindegs(jtype)) THEN
          WRITE(numout, FMT = *) '             ' // '    (in degrees)'
        ELSE
          WRITE(numout, FMT = *) '             ' // '    (in metres)'
        END IF
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_setinterpopts
END MODULE diaobs