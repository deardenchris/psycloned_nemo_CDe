MODULE icbini
  USE dom_oce
  USE in_out_manager
  USE lib_mpp
  USE sbc_oce
  USE sbc_ice
  USE iom
  USE fldread
  USE lbclnk
  USE icb_oce
  USE icbutl
  USE icbrst
  USE icbtrj
  USE icbdia
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: icb_init
  CHARACTER(LEN = 100) :: cn_dir = './'
  TYPE(FLD_N) :: sn_icb
  TYPE(FLD), PUBLIC, ALLOCATABLE, DIMENSION(:) :: sf_icb
  CONTAINS
  SUBROUTINE icb_init(pdt, kt)
    REAL(KIND = wp), INTENT(IN) :: pdt
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jn
    INTEGER :: i1, i2, i3
    INTEGER :: ii, inum, ivar
    INTEGER :: istat1, istat2, istat3
    CHARACTER(LEN = 300) :: cl_sdist
    CALL icb_nam
    IF (.NOT. ln_icebergs) RETURN
    IF (icb_alloc() /= 0) CALL ctl_stop('STOP', 'icb_alloc : unable to allocate arrays')
    uo_e(:, :) = 0._wp
    vo_e(:, :) = 0._wp
    ua_e(:, :) = 0._wp
    va_e(:, :) = 0._wp
    ff_e(:, :) = 0._wp
    tt_e(:, :) = 0._wp
    fr_e(:, :) = 0._wp
    hi_e(:,:) = 0._wp
    ui_e(:,:) = 0._wp ; vi_e(:,:) = 0._wp
    ssh_e(:, :) = 0._wp
    IF (nn_verbose_level > 0) THEN
      CALL ctl_opn(numicb, 'icebergs.stat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, numout, lwp, narea)
    END IF
    berg_dt = pdt
    first_width(:) = SQRT(rn_initial_mass(:) / (rn_LoW_ratio * rn_rho_bergs * rn_initial_thickness(:)))
    first_length(:) = rn_LoW_ratio * first_width(:)
    berg_grid(1) % calving(:, :) = 0._wp
    berg_grid(1) % calving_hflx(:, :) = 0._wp
    berg_grid(1) % stored_heat(:, :) = 0._wp
    berg_grid(1) % floating_melt(:, :) = 0._wp
    berg_grid(1) % maxclass(:, :) = nclasses
    berg_grid(1) % stored_ice(:, :, :) = 0._wp
    berg_grid(1) % tmp(:, :) = 0._wp
    src_calving(:, :) = 0._wp
    src_calving_hflx(:, :) = 0._wp
    IF (lk_mpp .AND. jpni == 1) CALL ctl_stop('icbinit: having ONE processor in x currently does not work')
    nicbpack = 10000
    IF (jpiglo >= nicbpack) CALL ctl_stop('icbini: processor index packing failure')
    nicbfldproc(:) = - 1
    DO jj = 1, jpj
      DO ji = 1, jpi
        src_calving_hflx(ji, jj) = narea
        src_calving(ji, jj) = nicbpack * mjg(jj) + mig(ji)
      END DO
    END DO
    CALL lbc_lnk('icbini', src_calving_hflx, 'T', 1._wp)
    CALL lbc_lnk('icbini', src_calving, 'T', 1._wp)
    jj = nlcj / 2
    nicbdi = - 1
    nicbei = - 1
    DO ji = 1, jpi
      i3 = INT(src_calving(ji, jj))
      i2 = INT(i3 / nicbpack)
      i1 = i3 - i2 * nicbpack
      i3 = INT(src_calving_hflx(ji, jj))
      IF (i1 == mig(ji) .AND. i3 == narea) THEN
        IF (nicbdi < 0) THEN
          nicbdi = ji
        ELSE
          nicbei = ji
        END IF
      END IF
    END DO
    ji = nlci / 2
    nicbdj = - 1
    nicbej = - 1
    DO jj = 1, jpj
      i3 = INT(src_calving(ji, jj))
      i2 = INT(i3 / nicbpack)
      i1 = i3 - i2 * nicbpack
      i3 = INT(src_calving_hflx(ji, jj))
      IF (i2 == mjg(jj) .AND. i3 == narea) THEN
        IF (nicbdj < 0) THEN
          nicbdj = jj
        ELSE
          nicbej = jj
        END IF
      END IF
    END DO
    i1 = MAX(nicbdi - 1, 1)
    i3 = INT(src_calving(i1, nlcj / 2))
    jj = INT(i3 / nicbpack)
    ricb_left = REAL(i3 - nicbpack * jj, wp)
    i1 = MIN(nicbei + 1, jpi)
    i3 = INT(src_calving(i1, nlcj / 2))
    jj = INT(i3 / nicbpack)
    ricb_right = REAL(i3 - nicbpack * jj, wp)
    IF (npolj > 0) THEN
      nicbfldpts(:) = INT(src_calving(:, nicbej + 1))
      nicbflddest(:) = INT(src_calving_hflx(:, nicbej + 1))
      DO ji = nicbdi, nicbei
        ii = nicbflddest(ji)
        IF (ii .GT. 0) THEN
          DO jn = 1, jpni
            IF (nicbfldproc(jn) == - 1) THEN
              nicbfldproc(jn) = ii
              EXIT
            END IF
            IF (nicbfldproc(jn) == ii) EXIT
          END DO
        END IF
      END DO
    END IF
    IF (nn_verbose_level > 0) THEN
      WRITE(numicb, FMT = *) 'processor ', narea
      WRITE(numicb, FMT = *) 'jpi, jpj   ', jpi, jpj
      WRITE(numicb, FMT = *) 'nldi, nlei ', nldi, nlei
      WRITE(numicb, FMT = *) 'nldj, nlej ', nldj, nlej
      WRITE(numicb, FMT = *) 'berg i interior ', nicbdi, nicbei
      WRITE(numicb, FMT = *) 'berg j interior ', nicbdj, nicbej
      WRITE(numicb, FMT = *) 'berg left       ', ricb_left
      WRITE(numicb, FMT = *) 'berg right      ', ricb_right
      jj = nlcj / 2
      WRITE(numicb, FMT = *) "central j line:"
      WRITE(numicb, FMT = *) "i processor"
      WRITE(numicb, FMT = *) (INT(src_calving_hflx(ji, jj)), ji = 1, jpi)
      WRITE(numicb, FMT = *) "i point"
      WRITE(numicb, FMT = *) (INT(src_calving(ji, jj)), ji = 1, jpi)
      ji = nlci / 2
      WRITE(numicb, FMT = *) "central i line:"
      WRITE(numicb, FMT = *) "j processor"
      WRITE(numicb, FMT = *) (INT(src_calving_hflx(ji, jj)), jj = 1, jpj)
      WRITE(numicb, FMT = *) "j point"
      WRITE(numicb, FMT = *) (INT(src_calving(ji, jj)), jj = 1, jpj)
      IF (npolj > 0) THEN
        WRITE(numicb, FMT = *) 'north fold destination points '
        WRITE(numicb, FMT = *) nicbfldpts
        WRITE(numicb, FMT = *) 'north fold destination procs  '
        WRITE(numicb, FMT = *) nicbflddest
        WRITE(numicb, FMT = *) 'north fold destination proclist  '
        WRITE(numicb, FMT = *) nicbfldproc
      END IF
      CALL flush(numicb)
    END IF
    src_calving(:, :) = 0._wp
    src_calving_hflx(:, :) = 0._wp
    tmask_e(:, :) = 0._wp
    tmask_e(1 : jpi, 1 : jpj) = tmask(:, :, 1)
    umask_e(:, :) = 0._wp
    umask_e(1 : jpi, 1 : jpj) = umask(:, :, 1)
    vmask_e(:, :) = 0._wp
    vmask_e(1 : jpi, 1 : jpj) = vmask(:, :, 1)
    CALL lbc_lnk_icb('icbini', tmask_e, 'T', + 1._wp, 1, 1)
    CALL lbc_lnk_icb('icbini', umask_e, 'T', + 1._wp, 1, 1)
    CALL lbc_lnk_icb('icbini', vmask_e, 'T', + 1._wp, 1, 1)
    num_bergs(:) = 0
    num_bergs(1) = narea - jpnij
    IF (nn_test_icebergs < 0 .OR. ln_use_calving) THEN
      cl_sdist = TRIM(cn_dir) // TRIM(sn_icb % clname)
      CALL iom_open(cl_sdist, inum)
      ivar = iom_varid(inum, 'maxclass', ldstop = .FALSE.)
      IF (ivar > 0) THEN
        CALL iom_get(inum, jpdom_data, 'maxclass', src_calving)
        berg_grid(1) % maxclass(:, :) = INT(src_calving)
        src_calving(:, :) = 0._wp
      END IF
      CALL iom_close(inum)
      IF (nn_verbose_level > 0) THEN
        WRITE(numicb, FMT = *)
        WRITE(numicb, FMT = *) '          calving read in a file'
      END IF
      ALLOCATE(sf_icb(1), STAT = istat1)
      ALLOCATE(sf_icb(1) % fnow(jpi, jpj, 1), STAT = istat2)
      ALLOCATE(sf_icb(1) % fdta(jpi, jpj, 1, 2), STAT = istat3)
      IF (istat1 + istat2 + istat3 > 0) THEN
        CALL ctl_stop('sbc_icb: unable to allocate sf_icb structure')
        RETURN
      END IF
      CALL fld_fill(sf_icb, (/sn_icb/), cn_dir, 'icb_init', 'read calving data', 'namicb')
    END IF
    IF (.NOT. ln_rstart) THEN
      IF (nn_test_icebergs > 0) CALL icb_ini_gen
    ELSE
      IF (nn_test_icebergs > 0) THEN
        CALL icb_ini_gen
      ELSE
        CALL icb_rst_read
        l_restarted_bergs = .TRUE.
      END IF
    END IF
    IF (nn_sample_rate .GT. 0) CALL icb_trj_init(nitend)
    CALL icb_dia_init
    IF (nn_verbose_level >= 2) CALL icb_utl_print('icb_init, initial status', nit000 - 1)
  END SUBROUTINE icb_init
  SUBROUTINE icb_ini_gen
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jj, ibergs
    TYPE(iceberg) :: localberg
    TYPE(point) :: localpt
    INTEGER :: iyr, imon, iday, ihr, imin, isec
    INTEGER :: iberg
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_ini_gen', 'r0', 0, 0)
    iberg = nn_test_icebergs
    iyr = nyear
    imon = nmonth
    iday = nday
    ihr = INT(nsec_day / 3600)
    imin = INT((nsec_day - ihr * 3600) / 60)
    isec = nsec_day - ihr * 3600 - imin * 60
    DO jj = nicbdj, nicbej
      DO ji = nicbdi, nicbei
        IF (tmask(ji, jj, 1) > 0._wp .AND. rn_test_box(1) < glamt(ji, jj) .AND. glamt(ji, jj) < rn_test_box(2) .AND. &
&rn_test_box(3) < gphit(ji, jj) .AND. gphit(ji, jj) < rn_test_box(4)) THEN
          localberg % mass_scaling = rn_mass_scaling(iberg)
          localpt % xi = REAL(mig(ji), wp)
          localpt % yj = REAL(mjg(jj), wp)
          localpt % lon = icb_utl_bilin(glamt, localpt % xi, localpt % yj, 'T')
          localpt % lat = icb_utl_bilin(gphit, localpt % xi, localpt % yj, 'T')
          localpt % mass = rn_initial_mass(iberg)
          localpt % thickness = rn_initial_thickness(iberg)
          localpt % width = first_width(iberg)
          localpt % length = first_length(iberg)
          localpt % year = iyr
          localpt % day = REAL(iday, wp) + (REAL(ihr, wp) + REAL(imin, wp) / 60._wp) / 24._wp
          localpt % mass_of_bits = 0._wp
          localpt % heat_density = 0._wp
          localpt % uvel = 0._wp
          localpt % vvel = 0._wp
          CALL icb_utl_incr
          localberg % number(:) = num_bergs(:)
          CALL icb_utl_add(localberg, localpt)
        END IF
      END DO
    END DO
    ibergs = icb_utl_count()
    CALL mpp_sum('icbini', ibergs)
    IF (nn_verbose_level > 0) THEN
      WRITE(numicb, FMT = '(a,i6,a)') 'diamonds, icb_ini_gen: ', ibergs, ' were generated'
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_ini_gen
  SUBROUTINE icb_nam
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: jn
    INTEGER :: ios
    REAL(KIND = wp) :: zfact
    NAMELIST /namberg/ ln_icebergs, ln_bergdia, nn_sample_rate, rn_initial_mass, rn_distribution, rn_mass_scaling, &
&rn_initial_thickness, nn_verbose_write, rn_rho_bergs, rn_LoW_ratio, nn_verbose_level, ln_operator_splitting, &
&rn_bits_erosion_fraction, rn_sicn_shift, ln_passive_mode, ln_time_average_weight, nn_test_icebergs, rn_test_box, ln_use_calving, &
&rn_speed_limit, cn_dir, sn_icb
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('icb_nam', 'r0', 0, 0)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'icb_nam : iceberg initialization through namberg namelist read'
      WRITE(numout, FMT = *) '~~~~~~~~ '
    END IF
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namberg, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namberg in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namberg, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namberg in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namberg)
    IF (lwp) WRITE(numout, FMT = *)
    CALL profile_psy_data0 % PostEnd
    IF (ln_icebergs) THEN
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   icebergs are used'
    ELSE
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   No icebergs used'
      RETURN
    END IF
    CALL profile_psy_data1 % PreStart('icb_nam', 'r1', 0, 0)
    IF (nn_test_icebergs > nclasses) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   Resetting of nn_test_icebergs to ', nclasses
      nn_test_icebergs = nclasses
    END IF
    IF (nn_test_icebergs < 0 .AND. .NOT. ln_use_calving) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   Resetting ln_use_calving to .true. since we are not using test icebergs'
      ln_use_calving = .TRUE.
    END IF
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'icb_nam : iceberg initialization through namberg namelist read'
      WRITE(numout, FMT = *) '~~~~~~~~ '
      WRITE(numout, FMT = *) '   Calculate budgets                                            ln_bergdia       = ', ln_bergdia
      WRITE(numout, FMT = *) '   Period between sampling of position for trajectory storage   nn_sample_rate = ', nn_sample_rate
      WRITE(numout, FMT = *) '   Mass thresholds between iceberg classes (kg)                 rn_initial_mass     ='
      DO jn = 1, nclasses
        WRITE(numout, FMT = '(a,f15.2)') '                                                                ', rn_initial_mass(jn)
      END DO
      WRITE(numout, FMT = *) '   Fraction of calving to apply to this class (non-dim)         rn_distribution     ='
      DO jn = 1, nclasses
        WRITE(numout, FMT = '(a,f10.4)') '                                                                ', rn_distribution(jn)
      END DO
      WRITE(numout, FMT = *) '   Ratio between effective and real iceberg mass (non-dim)      rn_mass_scaling     = '
      DO jn = 1, nclasses
        WRITE(numout, FMT = '(a,f10.2)') '                                                                ', rn_mass_scaling(jn)
      END DO
      WRITE(numout, FMT = *) '   Total thickness of newly calved bergs (m)                    rn_initial_thickness = '
      DO jn = 1, nclasses
        WRITE(numout, FMT = '(a,f10.2)') '                                                                ', &
&rn_initial_thickness(jn)
      END DO
      WRITE(numout, FMT = *) '   Timesteps between verbose messages                           nn_verbose_write    = ', &
&nn_verbose_write
      WRITE(numout, FMT = *) '   Density of icebergs                           rn_rho_bergs  = ', rn_rho_bergs
      WRITE(numout, FMT = *) '   Initial ratio L/W for newly calved icebergs   rn_LoW_ratio  = ', rn_LoW_ratio
      WRITE(numout, FMT = *) '   Turn on more verbose output                          level  = ', nn_verbose_level
      WRITE(numout, FMT = *) '   Use first order operator splitting for thermodynamics    ', 'use_operator_splitting = ', &
&ln_operator_splitting
      WRITE(numout, FMT = *) '   Fraction of erosion melt flux to divert to bergy bits    ', 'bits_erosion_fraction = ', &
&rn_bits_erosion_fraction
      WRITE(numout, FMT = *) '   Shift of sea-ice concentration in erosion flux modulation ', '(0<sicn_shift<1)    rn_sicn_shift  &
&= ', rn_sicn_shift
      WRITE(numout, FMT = *) '   Do not add freshwater flux from icebergs to ocean                ', '                  &
&passive_mode            = ', ln_passive_mode
      WRITE(numout, FMT = *) '   Time average the weight on the ocean   time_average_weight       = ', ln_time_average_weight
      WRITE(numout, FMT = *) '   Create icebergs in absence of a restart file   nn_test_icebergs  = ', nn_test_icebergs
      WRITE(numout, FMT = *) '                   in lon/lat box                                   = ', rn_test_box
      WRITE(numout, FMT = *) '   Use calving data even if nn_test_icebergs > 0    ln_use_calving  = ', ln_use_calving
      WRITE(numout, FMT = *) '   CFL speed limit for a berg            speed_limit                = ', rn_speed_limit
      WRITE(numout, FMT = *) '   Writing Iceberg status information to icebergs.stat file        '
    END IF
    zfact = SUM(rn_distribution)
    IF (zfact /= 1._wp .AND. 0_wp /= zfact) THEN
      rn_distribution(:) = rn_distribution(:) / zfact
      IF (lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) '      ==>>> CAUTION:    sum of berg input distribution = ', zfact
        WRITE(numout, FMT = *) '            *******     redistribution has been rescaled'
        WRITE(numout, FMT = *) '                        updated berg distribution is :'
        DO jn = 1, nclasses
          WRITE(numout, FMT = '(a,f10.4)') '                                   ', rn_distribution(jn)
        END DO
      END IF
    END IF
    IF (MINVAL(rn_distribution(:)) < 0._wp) THEN
      CALL ctl_stop('icb_nam: a negative rn_distribution value encountered ==>> change your namelist namberg')
    END IF
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE icb_nam
END MODULE icbini
