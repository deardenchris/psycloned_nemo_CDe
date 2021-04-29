MODULE bdyini
  USE oce
  USE dom_oce
  USE bdy_oce
  USE bdydta
  USE bdytides
  USE sbctide
  USE phycst, ONLY: rday
  USE in_out_manager
  USE lbclnk
  USE lib_mpp
  USE iom
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: bdy_init
  INTEGER, PARAMETER :: jp_nseg = 100
  INTEGER, PARAMETER :: nrimmax = 20
  INTEGER :: nbdysege, nbdysegw, nbdysegn, nbdysegs
  INTEGER, DIMENSION(jp_nseg) :: jpieob, jpjedt, jpjeft, npckge
  INTEGER, DIMENSION(jp_nseg) :: jpiwob, jpjwdt, jpjwft, npckgw
  INTEGER, DIMENSION(jp_nseg) :: jpjnob, jpindt, jpinft, npckgn
  INTEGER, DIMENSION(jp_nseg) :: jpjsob, jpisdt, jpisft, npckgs
  CONTAINS
  SUBROUTINE bdy_init
    NAMELIST /nambdy/ ln_bdy, nb_bdy, ln_coords_file, cn_coords_file, ln_mask_file, cn_mask_file, cn_dyn2d, nn_dyn2d_dta, cn_dyn3d, nn_dyn3d_dta, cn_tra, nn_tra_dta, ln_tra_dmp, ln_dyn3d_dmp, rn_time_dmp, rn_time_dmp_out, cn_ice, nn_ice_dta, rn_ice_tem, rn_ice_sal, rn_ice_age, ln_vol, nn_volctl, nn_rimwidth, nb_jpk_bdy
    INTEGER :: ios
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, nambdy, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'nambdy in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, nambdy, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'nambdy in configuration namelist', lwp)
    IF (lwm) WRITE(numond, nambdy)
    IF (.NOT. Agrif_Root()) ln_bdy = .FALSE.
    IF (ln_bdy) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'bdy_init : initialization of open boundaries'
      IF (lwp) WRITE(numout, *) '~~~~~~~~'
      CALL bdy_segs
      CALL bdy_dta_init
      IF (ln_tide) CALL bdytide_init
    ELSE
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'bdy_init : open boundaries not used (ln_bdy = F)'
      IF (lwp) WRITE(numout, *) '~~~~~~~~'
    END IF
  END SUBROUTINE bdy_init
  SUBROUTINE bdy_segs
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ib_bdy, ii, ij, ik, igrd, ib, ir, iseg
    INTEGER :: icount, icountr, ibr_max, ilen1, ibm1
    INTEGER :: iwe, ies, iso, ino, inum, id_dummy
    INTEGER :: igrd_start, igrd_end, jpbdta
    INTEGER :: jpbdtau, jpbdtas
    INTEGER :: ib_bdy1, ib_bdy2, ib1, ib2
    INTEGER :: i_offset, j_offset
    INTEGER, POINTER :: nbi, nbj, nbr
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: pmask
    REAL(KIND = wp) :: zefl, zwfl, znfl, zsfl
    INTEGER, DIMENSION(2) :: kdimsz
    INTEGER, DIMENSION(jpbgrd, jp_bdy) :: nblendta
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: nbidta, nbjdta
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: nbrdta
    CHARACTER(LEN = 1), DIMENSION(jpbgrd) :: cgrid
    INTEGER :: com_east, com_west, com_south, com_north, jpk_max
    INTEGER :: com_east_b, com_west_b, com_south_b, com_north_b
    INTEGER :: iw_b(4), ie_b(4), is_b(4), in_b(4)
    REAL(KIND = wp), TARGET, DIMENSION(jpi, jpj) :: zfmask
    CHARACTER(LEN = 1) :: ctypebdy
    INTEGER :: nbdyind, nbdybeg, nbdyend
    NAMELIST /nambdy_index/ ctypebdy, nbdyind, nbdybeg, nbdyend
    INTEGER :: ios
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bdy_segs', 'r0', 0, 0)
    cgrid = (/'t', 'u', 'v'/)
    IF (jperio /= 0) CALL ctl_stop('bdy_segs: Cyclic or symmetric,', ' and general open boundary condition are not compatible')
    IF (nb_bdy == 0) THEN
      IF (lwp) WRITE(numout, *) 'nb_bdy = 0, NO OPEN BOUNDARIES APPLIED.'
    ELSE
      IF (lwp) WRITE(numout, *) 'Number of open boundary sets : ', nb_bdy
    END IF
    DO ib_bdy = 1, nb_bdy
      IF (lwp) WRITE(numout, *) ' '
      IF (lwp) WRITE(numout, *) '------ Open boundary data set ', ib_bdy, '------'
      IF (ln_coords_file(ib_bdy)) THEN
        IF (lwp) WRITE(numout, *) 'Boundary definition read from file ' // TRIM(cn_coords_file(ib_bdy))
      ELSE
        IF (lwp) WRITE(numout, *) 'Boundary defined in namelist.'
      END IF
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'Boundary conditions for barotropic solution:  '
      SELECT CASE (cn_dyn2d(ib_bdy))
      CASE ('none')
        IF (lwp) WRITE(numout, *) '      no open boundary condition'
        dta_bdy(ib_bdy) % ll_ssh = .FALSE.
        dta_bdy(ib_bdy) % ll_u2d = .FALSE.
        dta_bdy(ib_bdy) % ll_v2d = .FALSE.
      CASE ('frs')
        IF (lwp) WRITE(numout, *) '      Flow Relaxation Scheme'
        dta_bdy(ib_bdy) % ll_ssh = .FALSE.
        dta_bdy(ib_bdy) % ll_u2d = .TRUE.
        dta_bdy(ib_bdy) % ll_v2d = .TRUE.
      CASE ('flather')
        IF (lwp) WRITE(numout, *) '      Flather radiation condition'
        dta_bdy(ib_bdy) % ll_ssh = .TRUE.
        dta_bdy(ib_bdy) % ll_u2d = .TRUE.
        dta_bdy(ib_bdy) % ll_v2d = .TRUE.
      CASE ('orlanski')
        IF (lwp) WRITE(numout, *) '      Orlanski (fully oblique) radiation condition with adaptive nudging'
        dta_bdy(ib_bdy) % ll_ssh = .FALSE.
        dta_bdy(ib_bdy) % ll_u2d = .TRUE.
        dta_bdy(ib_bdy) % ll_v2d = .TRUE.
      CASE ('orlanski_npo')
        IF (lwp) WRITE(numout, *) '      Orlanski (NPO) radiation condition with adaptive nudging'
        dta_bdy(ib_bdy) % ll_ssh = .FALSE.
        dta_bdy(ib_bdy) % ll_u2d = .TRUE.
        dta_bdy(ib_bdy) % ll_v2d = .TRUE.
      CASE DEFAULT
        CALL ctl_stop('unrecognised value for cn_dyn2d')
      END SELECT
      IF (cn_dyn2d(ib_bdy) /= 'none') THEN
        SELECT CASE (nn_dyn2d_dta(ib_bdy))
        CASE (0)
          IF (lwp) WRITE(numout, *) '      initial state used for bdy data'
        CASE (1)
          IF (lwp) WRITE(numout, *) '      boundary data taken from file'
        CASE (2)
          IF (lwp) WRITE(numout, *) '      tidal harmonic forcing taken from file'
        CASE (3)
          IF (lwp) WRITE(numout, *) '      boundary data AND tidal harmonic forcing taken from files'
        CASE DEFAULT
          CALL ctl_stop('nn_dyn2d_dta must be between 0 and 3')
        END SELECT
        IF ((nn_dyn2d_dta(ib_bdy) .GE. 2) .AND. (.NOT. ln_tide)) THEN
          CALL ctl_stop('You must activate with ln_tide to add tidal forcing at open boundaries')
        END IF
      END IF
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'Boundary conditions for baroclinic velocities:  '
      SELECT CASE (cn_dyn3d(ib_bdy))
      CASE ('none')
        IF (lwp) WRITE(numout, *) '      no open boundary condition'
        dta_bdy(ib_bdy) % ll_u3d = .FALSE.
        dta_bdy(ib_bdy) % ll_v3d = .FALSE.
      CASE ('frs')
        IF (lwp) WRITE(numout, *) '      Flow Relaxation Scheme'
        dta_bdy(ib_bdy) % ll_u3d = .TRUE.
        dta_bdy(ib_bdy) % ll_v3d = .TRUE.
      CASE ('specified')
        IF (lwp) WRITE(numout, *) '      Specified value'
        dta_bdy(ib_bdy) % ll_u3d = .TRUE.
        dta_bdy(ib_bdy) % ll_v3d = .TRUE.
      CASE ('neumann')
        IF (lwp) WRITE(numout, *) '      Neumann conditions'
        dta_bdy(ib_bdy) % ll_u3d = .FALSE.
        dta_bdy(ib_bdy) % ll_v3d = .FALSE.
      CASE ('zerograd')
        IF (lwp) WRITE(numout, *) '      Zero gradient for baroclinic velocities'
        dta_bdy(ib_bdy) % ll_u3d = .FALSE.
        dta_bdy(ib_bdy) % ll_v3d = .FALSE.
      CASE ('zero')
        IF (lwp) WRITE(numout, *) '      Zero baroclinic velocities (runoff case)'
        dta_bdy(ib_bdy) % ll_u3d = .FALSE.
        dta_bdy(ib_bdy) % ll_v3d = .FALSE.
      CASE ('orlanski')
        IF (lwp) WRITE(numout, *) '      Orlanski (fully oblique) radiation condition with adaptive nudging'
        dta_bdy(ib_bdy) % ll_u3d = .TRUE.
        dta_bdy(ib_bdy) % ll_v3d = .TRUE.
      CASE ('orlanski_npo')
        IF (lwp) WRITE(numout, *) '      Orlanski (NPO) radiation condition with adaptive nudging'
        dta_bdy(ib_bdy) % ll_u3d = .TRUE.
        dta_bdy(ib_bdy) % ll_v3d = .TRUE.
      CASE DEFAULT
        CALL ctl_stop('unrecognised value for cn_dyn3d')
      END SELECT
      IF (cn_dyn3d(ib_bdy) /= 'none') THEN
        SELECT CASE (nn_dyn3d_dta(ib_bdy))
        CASE (0)
          IF (lwp) WRITE(numout, *) '      initial state used for bdy data'
        CASE (1)
          IF (lwp) WRITE(numout, *) '      boundary data taken from file'
        CASE DEFAULT
          CALL ctl_stop('nn_dyn3d_dta must be 0 or 1')
        END SELECT
      END IF
      IF (ln_dyn3d_dmp(ib_bdy)) THEN
        IF (cn_dyn3d(ib_bdy) == 'none') THEN
          IF (lwp) WRITE(numout, *) 'No open boundary condition for baroclinic velocities: ln_dyn3d_dmp is set to .false.'
          ln_dyn3d_dmp(ib_bdy) = .FALSE.
        ELSE IF (cn_dyn3d(ib_bdy) == 'frs') THEN
          CALL ctl_stop('Use FRS OR relaxation')
        ELSE
          IF (lwp) WRITE(numout, *) '      + baroclinic velocities relaxation zone'
          IF (lwp) WRITE(numout, *) '      Damping time scale: ', rn_time_dmp(ib_bdy), ' days'
          IF ((lwp) .AND. rn_time_dmp(ib_bdy) < 0) CALL ctl_stop('Time scale must be positive')
          dta_bdy(ib_bdy) % ll_u3d = .TRUE.
          dta_bdy(ib_bdy) % ll_v3d = .TRUE.
        END IF
      ELSE
        IF (lwp) WRITE(numout, *) '      NO relaxation on baroclinic velocities'
      END IF
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'Boundary conditions for temperature and salinity:  '
      SELECT CASE (cn_tra(ib_bdy))
      CASE ('none')
        IF (lwp) WRITE(numout, *) '      no open boundary condition'
        dta_bdy(ib_bdy) % ll_tem = .FALSE.
        dta_bdy(ib_bdy) % ll_sal = .FALSE.
      CASE ('frs')
        IF (lwp) WRITE(numout, *) '      Flow Relaxation Scheme'
        dta_bdy(ib_bdy) % ll_tem = .TRUE.
        dta_bdy(ib_bdy) % ll_sal = .TRUE.
      CASE ('specified')
        IF (lwp) WRITE(numout, *) '      Specified value'
        dta_bdy(ib_bdy) % ll_tem = .TRUE.
        dta_bdy(ib_bdy) % ll_sal = .TRUE.
      CASE ('neumann')
        IF (lwp) WRITE(numout, *) '      Neumann conditions'
        dta_bdy(ib_bdy) % ll_tem = .FALSE.
        dta_bdy(ib_bdy) % ll_sal = .FALSE.
      CASE ('runoff')
        IF (lwp) WRITE(numout, *) '      Runoff conditions : Neumann for T and specified to 0.1 for salinity'
        dta_bdy(ib_bdy) % ll_tem = .FALSE.
        dta_bdy(ib_bdy) % ll_sal = .FALSE.
      CASE ('orlanski')
        IF (lwp) WRITE(numout, *) '      Orlanski (fully oblique) radiation condition with adaptive nudging'
        dta_bdy(ib_bdy) % ll_tem = .TRUE.
        dta_bdy(ib_bdy) % ll_sal = .TRUE.
      CASE ('orlanski_npo')
        IF (lwp) WRITE(numout, *) '      Orlanski (NPO) radiation condition with adaptive nudging'
        dta_bdy(ib_bdy) % ll_tem = .TRUE.
        dta_bdy(ib_bdy) % ll_sal = .TRUE.
      CASE DEFAULT
        CALL ctl_stop('unrecognised value for cn_tra')
      END SELECT
      IF (cn_tra(ib_bdy) /= 'none') THEN
        SELECT CASE (nn_tra_dta(ib_bdy))
        CASE (0)
          IF (lwp) WRITE(numout, *) '      initial state used for bdy data'
        CASE (1)
          IF (lwp) WRITE(numout, *) '      boundary data taken from file'
        CASE DEFAULT
          CALL ctl_stop('nn_tra_dta must be 0 or 1')
        END SELECT
      END IF
      IF (ln_tra_dmp(ib_bdy)) THEN
        IF (cn_tra(ib_bdy) == 'none') THEN
          IF (lwp) WRITE(numout, *) 'No open boundary condition for tracers: ln_tra_dmp is set to .false.'
          ln_tra_dmp(ib_bdy) = .FALSE.
        ELSE IF (cn_tra(ib_bdy) == 'frs') THEN
          CALL ctl_stop('Use FRS OR relaxation')
        ELSE
          IF (lwp) WRITE(numout, *) '      + T/S relaxation zone'
          IF (lwp) WRITE(numout, *) '      Damping time scale: ', rn_time_dmp(ib_bdy), ' days'
          IF (lwp) WRITE(numout, *) '      Outflow damping time scale: ', rn_time_dmp_out(ib_bdy), ' days'
          IF ((lwp) .AND. rn_time_dmp(ib_bdy) < 0) CALL ctl_stop('Time scale must be positive')
          dta_bdy(ib_bdy) % ll_tem = .TRUE.
          dta_bdy(ib_bdy) % ll_sal = .TRUE.
        END IF
      ELSE
        IF (lwp) WRITE(numout, *) '      NO T/S relaxation'
      END IF
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '      Width of relaxation zone = ', nn_rimwidth(ib_bdy)
      IF (lwp) WRITE(numout, *)
    END DO
    IF (nb_bdy > 0) THEN
      IF (ln_vol) THEN
        IF (lwp) WRITE(numout, *) 'Volume correction applied at open boundaries'
        IF (lwp) WRITE(numout, *)
        SELECT CASE (nn_volctl)
        CASE (1)
          IF (lwp) WRITE(numout, *) '      The total volume will be constant'
        CASE (0)
          IF (lwp) WRITE(numout, *) '      The total volume will vary according to the surface E-P flux'
        CASE DEFAULT
          CALL ctl_stop('nn_volctl must be 0 or 1')
        END SELECT
        IF (lwp) WRITE(numout, *)
        IF (ln_tide) THEN
          IF (lwp) WRITE(numout, *) ' The total volume correction is not working with tides. '
          IF (lwp) WRITE(numout, *) ' Set ln_vol to .FALSE. '
          IF (lwp) WRITE(numout, *) ' or '
          IF (lwp) WRITE(numout, *) ' equilibriate your bdy input files '
          CALL ctl_stop('The total volume correction is not working with tides.')
        END IF
      ELSE
        IF (lwp) WRITE(numout, *) 'No volume correction applied at open boundaries'
        IF (lwp) WRITE(numout, *)
      END IF
      IF (nb_jpk_bdy(ib_bdy) > 0) THEN
        IF (lwp) WRITE(numout, *) '*** open boundary will be interpolate in the vertical onto the native grid ***'
      ELSE
        IF (lwp) WRITE(numout, *) '*** open boundary will be read straight onto the native grid without vertical interpolation ***'
      END IF
    END IF
    REWIND(UNIT = numnam_cfg)
    nblendta(:, :) = 0
    nbdysege = 0
    nbdysegw = 0
    nbdysegn = 0
    nbdysegs = 0
    icount = 0
    jpbdtas = 1
    jpbdtau = 1
    DO ib_bdy = 1, nb_bdy
      IF (.NOT. ln_coords_file(ib_bdy)) THEN
        icount = icount + 1
        READ(numnam_cfg, nambdy_index, IOSTAT = ios, ERR = 904)
904     IF (ios /= 0) CALL ctl_nam(ios, 'nambdy_index in configuration namelist', lwp)
        IF (lwm) WRITE(numond, nambdy_index)
        SELECT CASE (TRIM(ctypebdy))
        CASE ('N')
          IF (nbdyind == - 1) THEN
            nbdyind = jpjglo - 2
            nbdybeg = 2
            nbdyend = jpiglo - 1
          END IF
          nbdysegn = nbdysegn + 1
          npckgn(nbdysegn) = ib_bdy
          jpjnob(nbdysegn) = nbdyind
          jpindt(nbdysegn) = nbdybeg
          jpinft(nbdysegn) = nbdyend
        CASE ('S')
          IF (nbdyind == - 1) THEN
            nbdyind = 2
            nbdybeg = 2
            nbdyend = jpiglo - 1
          END IF
          nbdysegs = nbdysegs + 1
          npckgs(nbdysegs) = ib_bdy
          jpjsob(nbdysegs) = nbdyind
          jpisdt(nbdysegs) = nbdybeg
          jpisft(nbdysegs) = nbdyend
        CASE ('E')
          IF (nbdyind == - 1) THEN
            nbdyind = jpiglo - 2
            nbdybeg = 2
            nbdyend = jpjglo - 1
          END IF
          nbdysege = nbdysege + 1
          npckge(nbdysege) = ib_bdy
          jpieob(nbdysege) = nbdyind
          jpjedt(nbdysege) = nbdybeg
          jpjeft(nbdysege) = nbdyend
        CASE ('W')
          IF (nbdyind == - 1) THEN
            nbdyind = 2
            nbdybeg = 2
            nbdyend = jpjglo - 1
          END IF
          nbdysegw = nbdysegw + 1
          npckgw(nbdysegw) = ib_bdy
          jpiwob(nbdysegw) = nbdyind
          jpjwdt(nbdysegw) = nbdybeg
          jpjwft(nbdysegw) = nbdyend
        CASE DEFAULT
          CALL ctl_stop('ctypebdy must be N, S, E or W')
        END SELECT
        nblendta(1 : jpbgrd, ib_bdy) = (nbdyend - nbdybeg + 1) * nn_rimwidth(ib_bdy)
        jpbdtas = MAX(jpbdtas, (nbdyend - nbdybeg + 1))
        IF (lwp .AND. (nn_rimwidth(ib_bdy) > nrimmax)) CALL ctl_stop('rimwidth must be lower than nrimmax')
      ELSE
        CALL iom_open(cn_coords_file(ib_bdy), inum)
        DO igrd = 1, jpbgrd
          id_dummy = iom_varid(inum, 'nbi' // cgrid(igrd), kdimsz = kdimsz)
          nblendta(igrd, ib_bdy) = MAXVAL(kdimsz)
          jpbdtau = MAX(jpbdtau, MAXVAL(kdimsz))
        END DO
        CALL iom_close(inum)
      END IF
    END DO
    IF (nb_bdy > 0) THEN
      jpbdta = MAXVAL(nblendta(1 : jpbgrd, 1 : nb_bdy))
      ALLOCATE(nbidta(jpbdta, jpbgrd, nb_bdy), nbjdta(jpbdta, jpbgrd, nb_bdy), nbrdta(jpbdta, jpbgrd, nb_bdy))
      jpk_max = MAXVAL(nb_jpk_bdy)
      jpk_max = MAX(jpk_max, jpk)
      ALLOCATE(dta_global(jpbdtau, 1, jpk_max))
      ALLOCATE(dta_global_z(jpbdtau, 1, jpk_max))
      ALLOCATE(dta_global_dz(jpbdtau, 1, jpk_max))
      IF (icount > 0) THEN
        ALLOCATE(dta_global2(jpbdtas, nrimmax, jpk_max))
        ALLOCATE(dta_global2_z(jpbdtas, nrimmax, jpk_max))
        ALLOCATE(dta_global2_dz(jpbdtas, nrimmax, jpk_max))
      END IF
    END IF
    IF (icount > 0) CALL bdy_ctl_seg
    DO ib_bdy = 1, nb_bdy
      IF (ln_coords_file(ib_bdy)) THEN
        CALL iom_open(cn_coords_file(ib_bdy), inum)
        DO igrd = 1, jpbgrd
          CALL iom_get(inum, jpdom_unknown, 'nbi' // cgrid(igrd), dta_global(1 : nblendta(igrd, ib_bdy), :, 1))
          DO ii = 1, nblendta(igrd, ib_bdy)
            nbidta(ii, igrd, ib_bdy) = INT(dta_global(ii, 1, 1))
          END DO
          CALL iom_get(inum, jpdom_unknown, 'nbj' // cgrid(igrd), dta_global(1 : nblendta(igrd, ib_bdy), :, 1))
          DO ii = 1, nblendta(igrd, ib_bdy)
            nbjdta(ii, igrd, ib_bdy) = INT(dta_global(ii, 1, 1))
          END DO
          CALL iom_get(inum, jpdom_unknown, 'nbr' // cgrid(igrd), dta_global(1 : nblendta(igrd, ib_bdy), :, 1))
          DO ii = 1, nblendta(igrd, ib_bdy)
            nbrdta(ii, igrd, ib_bdy) = INT(dta_global(ii, 1, 1))
          END DO
          ibr_max = MAXVAL(nbrdta(:, igrd, ib_bdy))
          IF (lwp) WRITE(numout, *)
          IF (lwp) WRITE(numout, *) ' Maximum rimwidth in file is ', ibr_max
          IF (lwp) WRITE(numout, *) ' nn_rimwidth from namelist is ', nn_rimwidth(ib_bdy)
          IF (ibr_max < nn_rimwidth(ib_bdy)) CALL ctl_stop('nn_rimwidth is larger than maximum rimwidth in file', cn_coords_file(ib_bdy))
        END DO
        CALL iom_close(inum)
      END IF
    END DO
    DO iseg = 1, nbdysege
      ib_bdy = npckge(iseg)
      igrd = 1
      icount = 0
      DO ir = 1, nn_rimwidth(ib_bdy)
        DO ij = jpjedt(iseg), jpjeft(iseg)
          icount = icount + 1
          nbidta(icount, igrd, ib_bdy) = jpieob(iseg) + 2 - ir
          nbjdta(icount, igrd, ib_bdy) = ij
          nbrdta(icount, igrd, ib_bdy) = ir
        END DO
      END DO
      igrd = 2
      icount = 0
      DO ir = 1, nn_rimwidth(ib_bdy)
        DO ij = jpjedt(iseg), jpjeft(iseg)
          icount = icount + 1
          nbidta(icount, igrd, ib_bdy) = jpieob(iseg) + 1 - ir
          nbjdta(icount, igrd, ib_bdy) = ij
          nbrdta(icount, igrd, ib_bdy) = ir
        END DO
      END DO
      igrd = 3
      icount = 0
      DO ir = 1, nn_rimwidth(ib_bdy)
        DO ij = jpjedt(iseg), jpjeft(iseg)
          icount = icount + 1
          nbidta(icount, igrd, ib_bdy) = jpieob(iseg) + 2 - ir
          nbjdta(icount, igrd, ib_bdy) = ij
          nbrdta(icount, igrd, ib_bdy) = ir
        END DO
        nbidta(icount, igrd, ib_bdy) = - ib_bdy
        nbjdta(icount, igrd, ib_bdy) = - ib_bdy
      END DO
    END DO
    DO iseg = 1, nbdysegw
      ib_bdy = npckgw(iseg)
      igrd = 1
      icount = 0
      DO ir = 1, nn_rimwidth(ib_bdy)
        DO ij = jpjwdt(iseg), jpjwft(iseg)
          icount = icount + 1
          nbidta(icount, igrd, ib_bdy) = jpiwob(iseg) + ir - 1
          nbjdta(icount, igrd, ib_bdy) = ij
          nbrdta(icount, igrd, ib_bdy) = ir
        END DO
      END DO
      igrd = 2
      icount = 0
      DO ir = 1, nn_rimwidth(ib_bdy)
        DO ij = jpjwdt(iseg), jpjwft(iseg)
          icount = icount + 1
          nbidta(icount, igrd, ib_bdy) = jpiwob(iseg) + ir - 1
          nbjdta(icount, igrd, ib_bdy) = ij
          nbrdta(icount, igrd, ib_bdy) = ir
        END DO
      END DO
      igrd = 3
      icount = 0
      DO ir = 1, nn_rimwidth(ib_bdy)
        DO ij = jpjwdt(iseg), jpjwft(iseg)
          icount = icount + 1
          nbidta(icount, igrd, ib_bdy) = jpiwob(iseg) + ir - 1
          nbjdta(icount, igrd, ib_bdy) = ij
          nbrdta(icount, igrd, ib_bdy) = ir
        END DO
        nbidta(icount, igrd, ib_bdy) = - ib_bdy
        nbjdta(icount, igrd, ib_bdy) = - ib_bdy
      END DO
    END DO
    DO iseg = 1, nbdysegn
      ib_bdy = npckgn(iseg)
      igrd = 1
      icount = 0
      DO ir = 1, nn_rimwidth(ib_bdy)
        DO ii = jpindt(iseg), jpinft(iseg)
          icount = icount + 1
          nbidta(icount, igrd, ib_bdy) = ii
          nbjdta(icount, igrd, ib_bdy) = jpjnob(iseg) + 2 - ir
          nbrdta(icount, igrd, ib_bdy) = ir
        END DO
      END DO
      igrd = 2
      icount = 0
      DO ir = 1, nn_rimwidth(ib_bdy)
        DO ii = jpindt(iseg), jpinft(iseg)
          icount = icount + 1
          nbidta(icount, igrd, ib_bdy) = ii
          nbjdta(icount, igrd, ib_bdy) = jpjnob(iseg) + 2 - ir
          nbrdta(icount, igrd, ib_bdy) = ir
        END DO
        nbidta(icount, igrd, ib_bdy) = - ib_bdy
        nbjdta(icount, igrd, ib_bdy) = - ib_bdy
      END DO
      igrd = 3
      icount = 0
      DO ir = 1, nn_rimwidth(ib_bdy)
        DO ii = jpindt(iseg), jpinft(iseg)
          icount = icount + 1
          nbidta(icount, igrd, ib_bdy) = ii
          nbjdta(icount, igrd, ib_bdy) = jpjnob(iseg) + 1 - ir
          nbrdta(icount, igrd, ib_bdy) = ir
        END DO
      END DO
    END DO
    DO iseg = 1, nbdysegs
      ib_bdy = npckgs(iseg)
      igrd = 1
      icount = 0
      DO ir = 1, nn_rimwidth(ib_bdy)
        DO ii = jpisdt(iseg), jpisft(iseg)
          icount = icount + 1
          nbidta(icount, igrd, ib_bdy) = ii
          nbjdta(icount, igrd, ib_bdy) = jpjsob(iseg) + ir - 1
          nbrdta(icount, igrd, ib_bdy) = ir
        END DO
      END DO
      igrd = 2
      icount = 0
      DO ir = 1, nn_rimwidth(ib_bdy)
        DO ii = jpisdt(iseg), jpisft(iseg)
          icount = icount + 1
          nbidta(icount, igrd, ib_bdy) = ii
          nbjdta(icount, igrd, ib_bdy) = jpjsob(iseg) + ir - 1
          nbrdta(icount, igrd, ib_bdy) = ir
        END DO
        nbidta(icount, igrd, ib_bdy) = - ib_bdy
        nbjdta(icount, igrd, ib_bdy) = - ib_bdy
      END DO
      igrd = 3
      icount = 0
      DO ir = 1, nn_rimwidth(ib_bdy)
        DO ii = jpisdt(iseg), jpisft(iseg)
          icount = icount + 1
          nbidta(icount, igrd, ib_bdy) = ii
          nbjdta(icount, igrd, ib_bdy) = jpjsob(iseg) + ir - 1
          nbrdta(icount, igrd, ib_bdy) = ir
        END DO
      END DO
    END DO
    DO igrd = 1, jpbgrd
      DO ib_bdy1 = 1, nb_bdy
        DO ib_bdy2 = 1, nb_bdy
          IF (ib_bdy1 /= ib_bdy2) THEN
            DO ib1 = 1, nblendta(igrd, ib_bdy1)
              DO ib2 = 1, nblendta(igrd, ib_bdy2)
                IF ((nbidta(ib1, igrd, ib_bdy1) == nbidta(ib2, igrd, ib_bdy2)) .AND. (nbjdta(ib1, igrd, ib_bdy1) == nbjdta(ib2, igrd, ib_bdy2))) THEN
                  IF (nbrdta(ib1, igrd, ib_bdy1) < nbrdta(ib2, igrd, ib_bdy2)) THEN
                    nbidta(ib2, igrd, ib_bdy2) = - ib_bdy2
                    nbjdta(ib2, igrd, ib_bdy2) = - ib_bdy2
                  ELSE IF (nbrdta(ib1, igrd, ib_bdy1) > nbrdta(ib2, igrd, ib_bdy2)) THEN
                    nbidta(ib1, igrd, ib_bdy1) = - ib_bdy1
                    nbjdta(ib1, igrd, ib_bdy1) = - ib_bdy1
                  ELSE
                    nbidta(ib1, igrd, ib_bdy1) = - ib_bdy1
                    nbjdta(ib1, igrd, ib_bdy1) = - ib_bdy1
                  END IF
                END IF
              END DO
            END DO
          END IF
        END DO
      END DO
    END DO
    iwe = mig(1) - 1 + 2
    ies = mig(1) + nlci - 1 - 1
    iso = mjg(1) - 1 + 2
    ino = mjg(1) + nlcj - 1 - 1
    ALLOCATE(nbondi_bdy(nb_bdy))
    ALLOCATE(nbondj_bdy(nb_bdy))
    nbondi_bdy(:) = 2
    nbondj_bdy(:) = 2
    ALLOCATE(nbondi_bdy_b(nb_bdy))
    ALLOCATE(nbondj_bdy_b(nb_bdy))
    nbondi_bdy_b(:) = 2
    nbondj_bdy_b(:) = 2
    IF (nbondi == 0) THEN
      iw_b(1) = 1 + nimppt(nowe + 1)
      ie_b(1) = 1 + nimppt(nowe + 1) + nlcit(nowe + 1) - 3
      is_b(1) = 1 + njmppt(nowe + 1)
      in_b(1) = 1 + njmppt(nowe + 1) + nlcjt(nowe + 1) - 3
      iw_b(2) = 1 + nimppt(noea + 1)
      ie_b(2) = 1 + nimppt(noea + 1) + nlcit(noea + 1) - 3
      is_b(2) = 1 + njmppt(noea + 1)
      in_b(2) = 1 + njmppt(noea + 1) + nlcjt(noea + 1) - 3
    ELSE IF (nbondi == 1) THEN
      iw_b(1) = 1 + nimppt(nowe + 1)
      ie_b(1) = 1 + nimppt(nowe + 1) + nlcit(nowe + 1) - 3
      is_b(1) = 1 + njmppt(nowe + 1)
      in_b(1) = 1 + njmppt(nowe + 1) + nlcjt(nowe + 1) - 3
    ELSE IF (nbondi == - 1) THEN
      iw_b(2) = 1 + nimppt(noea + 1)
      ie_b(2) = 1 + nimppt(noea + 1) + nlcit(noea + 1) - 3
      is_b(2) = 1 + njmppt(noea + 1)
      in_b(2) = 1 + njmppt(noea + 1) + nlcjt(noea + 1) - 3
    END IF
    IF (nbondj == 0) THEN
      iw_b(3) = 1 + nimppt(noso + 1)
      ie_b(3) = 1 + nimppt(noso + 1) + nlcit(noso + 1) - 3
      is_b(3) = 1 + njmppt(noso + 1)
      in_b(3) = 1 + njmppt(noso + 1) + nlcjt(noso + 1) - 3
      iw_b(4) = 1 + nimppt(nono + 1)
      ie_b(4) = 1 + nimppt(nono + 1) + nlcit(nono + 1) - 3
      is_b(4) = 1 + njmppt(nono + 1)
      in_b(4) = 1 + njmppt(nono + 1) + nlcjt(nono + 1) - 3
    ELSE IF (nbondj == 1) THEN
      iw_b(3) = 1 + nimppt(noso + 1)
      ie_b(3) = 1 + nimppt(noso + 1) + nlcit(noso + 1) - 3
      is_b(3) = 1 + njmppt(noso + 1)
      in_b(3) = 1 + njmppt(noso + 1) + nlcjt(noso + 1) - 3
    ELSE IF (nbondj == - 1) THEN
      iw_b(4) = 1 + nimppt(nono + 1)
      ie_b(4) = 1 + nimppt(nono + 1) + nlcit(nono + 1) - 3
      is_b(4) = 1 + njmppt(nono + 1)
      in_b(4) = 1 + njmppt(nono + 1) + nlcjt(nono + 1) - 3
    END IF
    DO ib_bdy = 1, nb_bdy
      DO igrd = 1, jpbgrd
        icount = 0
        icountr = 0
        idx_bdy(ib_bdy) % nblen(igrd) = 0
        idx_bdy(ib_bdy) % nblenrim(igrd) = 0
        DO ib = 1, nblendta(igrd, ib_bdy)
          ibm1 = MAX(1, ib - 1)
          IF (lwp) THEN
            IF (nbrdta(ib, igrd, ib_bdy) < nbrdta(ibm1, igrd, ib_bdy)) THEN
              CALL ctl_stop('bdy_segs : ERROR : boundary data in file must be defined ', ' in order of distance from edge nbr A utility for re-ordering ', ' boundary coordinates and data files exists in the TOOLS/OBC directory')
            END IF
          END IF
          IF (nbidta(ib, igrd, ib_bdy) >= iwe .AND. nbidta(ib, igrd, ib_bdy) <= ies .AND. nbjdta(ib, igrd, ib_bdy) >= iso .AND. nbjdta(ib, igrd, ib_bdy) <= ino) THEN
            icount = icount + 1
            IF (nbrdta(ib, igrd, ib_bdy) == 1) icountr = icountr + 1
          END IF
        END DO
        idx_bdy(ib_bdy) % nblenrim(igrd) = icountr
        idx_bdy(ib_bdy) % nblen(igrd) = icount
      END DO
      ilen1 = MAXVAL(idx_bdy(ib_bdy) % nblen(:))
      ALLOCATE(idx_bdy(ib_bdy) % nbi(ilen1, jpbgrd), idx_bdy(ib_bdy) % nbj(ilen1, jpbgrd), idx_bdy(ib_bdy) % nbr(ilen1, jpbgrd), idx_bdy(ib_bdy) % nbd(ilen1, jpbgrd), idx_bdy(ib_bdy) % nbdout(ilen1, jpbgrd), idx_bdy(ib_bdy) % nbmap(ilen1, jpbgrd), idx_bdy(ib_bdy) % nbw(ilen1, jpbgrd), idx_bdy(ib_bdy) % flagu(ilen1, jpbgrd), idx_bdy(ib_bdy) % flagv(ilen1, jpbgrd))
      com_east = 0
      com_west = 0
      com_south = 0
      com_north = 0
      com_east_b = 0
      com_west_b = 0
      com_south_b = 0
      com_north_b = 0
      DO igrd = 1, jpbgrd
        icount = 0
        DO ir = 1, nn_rimwidth(ib_bdy)
          DO ib = 1, nblendta(igrd, ib_bdy)
            IF (nbidta(ib, igrd, ib_bdy) >= iwe .AND. nbidta(ib, igrd, ib_bdy) <= ies .AND. nbjdta(ib, igrd, ib_bdy) >= iso .AND. nbjdta(ib, igrd, ib_bdy) <= ino .AND. nbrdta(ib, igrd, ib_bdy) == ir) THEN
              icount = icount + 1
              idx_bdy(ib_bdy) % nbi(icount, igrd) = nbidta(ib, igrd, ib_bdy) - mig(1) + 1
              idx_bdy(ib_bdy) % nbj(icount, igrd) = nbjdta(ib, igrd, ib_bdy) - mjg(1) + 1
              ii = idx_bdy(ib_bdy) % nbi(icount, igrd)
              ij = idx_bdy(ib_bdy) % nbj(icount, igrd)
              IF ((com_east .NE. 1) .AND. (ii == (nlci - 1)) .AND. (nbondi .LE. 0)) THEN
                com_east = 1
              ELSE IF ((com_west .NE. 1) .AND. (ii == 2) .AND. (nbondi .GE. 0) .AND. (nbondi .NE. 2)) THEN
                com_west = 1
              END IF
              IF ((com_south .NE. 1) .AND. (ij == 2) .AND. (nbondj .GE. 0) .AND. (nbondj .NE. 2)) THEN
                com_south = 1
              ELSE IF ((com_north .NE. 1) .AND. (ij == (nlcj - 1)) .AND. (nbondj .LE. 0)) THEN
                com_north = 1
              END IF
              idx_bdy(ib_bdy) % nbr(icount, igrd) = nbrdta(ib, igrd, ib_bdy)
              idx_bdy(ib_bdy) % nbmap(icount, igrd) = ib
            END IF
            IF (nbondi == 0) THEN
              IF (nbidta(ib, igrd, ib_bdy) >= iw_b(1) .AND. nbidta(ib, igrd, ib_bdy) <= ie_b(1) .AND. nbjdta(ib, igrd, ib_bdy) >= is_b(1) .AND. nbjdta(ib, igrd, ib_bdy) <= in_b(1) .AND. nbrdta(ib, igrd, ib_bdy) == ir) THEN
                ii = nbidta(ib, igrd, ib_bdy) - iw_b(1) + 2
                IF (ii == (nlcit(nowe + 1) - 1)) THEN
                  ij = nbjdta(ib, igrd, ib_bdy) - is_b(1) + 2
                  IF ((ij == 2) .AND. (nbondj == 0 .OR. nbondj == 1)) THEN
                    com_south = 1
                  ELSE IF ((ij == nlcjt(nowe + 1) - 1) .AND. (nbondj == 0 .OR. nbondj == - 1)) THEN
                    com_north = 1
                  END IF
                  com_west_b = 1
                END IF
              END IF
              IF (nbidta(ib, igrd, ib_bdy) >= iw_b(2) .AND. nbidta(ib, igrd, ib_bdy) <= ie_b(2) .AND. nbjdta(ib, igrd, ib_bdy) >= is_b(2) .AND. nbjdta(ib, igrd, ib_bdy) <= in_b(2) .AND. nbrdta(ib, igrd, ib_bdy) == ir) THEN
                ii = nbidta(ib, igrd, ib_bdy) - iw_b(2) + 2
                IF (ii == 2) THEN
                  ij = nbjdta(ib, igrd, ib_bdy) - is_b(2) + 2
                  IF ((ij == 2) .AND. (nbondj == 0 .OR. nbondj == 1)) THEN
                    com_south = 1
                  ELSE IF ((ij == nlcjt(noea + 1) - 1) .AND. (nbondj == 0 .OR. nbondj == - 1)) THEN
                    com_north = 1
                  END IF
                  com_east_b = 1
                END IF
              END IF
            ELSE IF (nbondi == 1) THEN
              IF (nbidta(ib, igrd, ib_bdy) >= iw_b(1) .AND. nbidta(ib, igrd, ib_bdy) <= ie_b(1) .AND. nbjdta(ib, igrd, ib_bdy) >= is_b(1) .AND. nbjdta(ib, igrd, ib_bdy) <= in_b(1) .AND. nbrdta(ib, igrd, ib_bdy) == ir) THEN
                ii = nbidta(ib, igrd, ib_bdy) - iw_b(1) + 2
                IF (ii == (nlcit(nowe + 1) - 1)) THEN
                  ij = nbjdta(ib, igrd, ib_bdy) - is_b(1) + 2
                  IF ((ij == 2) .AND. (nbondj == 0 .OR. nbondj == 1)) THEN
                    com_south = 1
                  ELSE IF ((ij == nlcjt(nowe + 1) - 1) .AND. (nbondj == 0 .OR. nbondj == - 1)) THEN
                    com_north = 1
                  END IF
                  com_west_b = 1
                END IF
              END IF
            ELSE IF (nbondi == - 1) THEN
              IF (nbidta(ib, igrd, ib_bdy) >= iw_b(2) .AND. nbidta(ib, igrd, ib_bdy) <= ie_b(2) .AND. nbjdta(ib, igrd, ib_bdy) >= is_b(2) .AND. nbjdta(ib, igrd, ib_bdy) <= in_b(2) .AND. nbrdta(ib, igrd, ib_bdy) == ir) THEN
                ii = nbidta(ib, igrd, ib_bdy) - iw_b(2) + 2
                IF (ii == 2) THEN
                  ij = nbjdta(ib, igrd, ib_bdy) - is_b(2) + 2
                  IF ((ij == 2) .AND. (nbondj == 0 .OR. nbondj == 1)) THEN
                    com_south = 1
                  ELSE IF ((ij == nlcjt(noea + 1) - 1) .AND. (nbondj == 0 .OR. nbondj == - 1)) THEN
                    com_north = 1
                  END IF
                  com_east_b = 1
                END IF
              END IF
            END IF
            IF (nbondj == 0) THEN
              IF (com_north_b .NE. 1 .AND. (nbidta(ib, igrd, ib_bdy) == iw_b(4) - 1 .OR. nbidta(ib, igrd, ib_bdy) == ie_b(4) + 1) .AND. nbjdta(ib, igrd, ib_bdy) == is_b(4) .AND. nbrdta(ib, igrd, ib_bdy) == ir) THEN
                com_north_b = 1
              END IF
              IF (com_south_b .NE. 1 .AND. (nbidta(ib, igrd, ib_bdy) == iw_b(3) - 1 .OR. nbidta(ib, igrd, ib_bdy) == ie_b(3) + 1) .AND. nbjdta(ib, igrd, ib_bdy) == in_b(3) .AND. nbrdta(ib, igrd, ib_bdy) == ir) THEN
                com_south_b = 1
              END IF
              IF (nbidta(ib, igrd, ib_bdy) >= iw_b(3) .AND. nbidta(ib, igrd, ib_bdy) <= ie_b(3) .AND. nbjdta(ib, igrd, ib_bdy) >= is_b(3) .AND. nbjdta(ib, igrd, ib_bdy) <= in_b(3) .AND. nbrdta(ib, igrd, ib_bdy) == ir) THEN
                ij = nbjdta(ib, igrd, ib_bdy) - is_b(3) + 2
                IF ((com_south_b .NE. 1) .AND. (ij == (nlcjt(noso + 1) - 1))) THEN
                  com_south_b = 1
                END IF
              END IF
              IF (nbidta(ib, igrd, ib_bdy) >= iw_b(4) .AND. nbidta(ib, igrd, ib_bdy) <= ie_b(4) .AND. nbjdta(ib, igrd, ib_bdy) >= is_b(4) .AND. nbjdta(ib, igrd, ib_bdy) <= in_b(4) .AND. nbrdta(ib, igrd, ib_bdy) == ir) THEN
                ij = nbjdta(ib, igrd, ib_bdy) - is_b(4) + 2
                IF ((com_north_b .NE. 1) .AND. (ij == 2)) THEN
                  com_north_b = 1
                END IF
              END IF
            ELSE IF (nbondj == 1) THEN
              IF (com_south_b .NE. 1 .AND. (nbidta(ib, igrd, ib_bdy) == iw_b(3) - 1 .OR. nbidta(ib, igrd, ib_bdy) == ie_b(3) + 1) .AND. nbjdta(ib, igrd, ib_bdy) == in_b(3) .AND. nbrdta(ib, igrd, ib_bdy) == ir) THEN
                com_south_b = 1
              END IF
              IF (nbidta(ib, igrd, ib_bdy) >= iw_b(3) .AND. nbidta(ib, igrd, ib_bdy) <= ie_b(3) .AND. nbjdta(ib, igrd, ib_bdy) >= is_b(3) .AND. nbjdta(ib, igrd, ib_bdy) <= in_b(3) .AND. nbrdta(ib, igrd, ib_bdy) == ir) THEN
                ij = nbjdta(ib, igrd, ib_bdy) - is_b(3) + 2
                IF ((com_south_b .NE. 1) .AND. (ij == (nlcjt(noso + 1) - 1))) THEN
                  com_south_b = 1
                END IF
              END IF
            ELSE IF (nbondj == - 1) THEN
              IF (com_north_b .NE. 1 .AND. (nbidta(ib, igrd, ib_bdy) == iw_b(4) - 1 .OR. nbidta(ib, igrd, ib_bdy) == ie_b(4) + 1) .AND. nbjdta(ib, igrd, ib_bdy) == is_b(4) .AND. nbrdta(ib, igrd, ib_bdy) == ir) THEN
                com_north_b = 1
              END IF
              IF (nbidta(ib, igrd, ib_bdy) >= iw_b(4) .AND. nbidta(ib, igrd, ib_bdy) <= ie_b(4) .AND. nbjdta(ib, igrd, ib_bdy) >= is_b(4) .AND. nbjdta(ib, igrd, ib_bdy) <= in_b(4) .AND. nbrdta(ib, igrd, ib_bdy) == ir) THEN
                ij = nbjdta(ib, igrd, ib_bdy) - is_b(4) + 2
                IF ((com_north_b .NE. 1) .AND. (ij == 2)) THEN
                  com_north_b = 1
                END IF
              END IF
            END IF
          END DO
        END DO
      END DO
      IF ((com_east == 1) .AND. (com_west == 1)) THEN
        nbondi_bdy(ib_bdy) = 0
      ELSE IF ((com_east == 1) .AND. (com_west == 0)) THEN
        nbondi_bdy(ib_bdy) = - 1
      ELSE IF ((com_east == 0) .AND. (com_west == 1)) THEN
        nbondi_bdy(ib_bdy) = 1
      END IF
      IF ((com_north == 1) .AND. (com_south == 1)) THEN
        nbondj_bdy(ib_bdy) = 0
      ELSE IF ((com_north == 1) .AND. (com_south == 0)) THEN
        nbondj_bdy(ib_bdy) = - 1
      ELSE IF ((com_north == 0) .AND. (com_south == 1)) THEN
        nbondj_bdy(ib_bdy) = 1
      END IF
      IF ((com_east_b == 1) .AND. (com_west_b == 1)) THEN
        nbondi_bdy_b(ib_bdy) = 0
      ELSE IF ((com_east_b == 1) .AND. (com_west_b == 0)) THEN
        nbondi_bdy_b(ib_bdy) = - 1
      ELSE IF ((com_east_b == 0) .AND. (com_west_b == 1)) THEN
        nbondi_bdy_b(ib_bdy) = 1
      END IF
      IF ((com_north_b == 1) .AND. (com_south_b == 1)) THEN
        nbondj_bdy_b(ib_bdy) = 0
      ELSE IF ((com_north_b == 1) .AND. (com_south_b == 0)) THEN
        nbondj_bdy_b(ib_bdy) = - 1
      ELSE IF ((com_north_b == 0) .AND. (com_south_b == 1)) THEN
        nbondj_bdy_b(ib_bdy) = 1
      END IF
      DO igrd = 1, jpbgrd
        DO ib = 1, idx_bdy(ib_bdy) % nblen(igrd)
          nbr => idx_bdy(ib_bdy) % nbr(ib, igrd)
          idx_bdy(ib_bdy) % nbw(ib, igrd) = 1. - TANH(REAL(nbr - 1) * 0.5)
        END DO
      END DO
      DO igrd = 1, jpbgrd
        DO ib = 1, idx_bdy(ib_bdy) % nblen(igrd)
          nbr => idx_bdy(ib_bdy) % nbr(ib, igrd)
          idx_bdy(ib_bdy) % nbd(ib, igrd) = 1. / (rn_time_dmp(ib_bdy) * rday) * (REAL(nn_rimwidth(ib_bdy) + 1 - nbr) / REAL(nn_rimwidth(ib_bdy))) ** 2.
          idx_bdy(ib_bdy) % nbdout(ib, igrd) = 1. / (rn_time_dmp_out(ib_bdy) * rday) * (REAL(nn_rimwidth(ib_bdy) + 1 - nbr) / REAL(nn_rimwidth(ib_bdy))) ** 2.
        END DO
      END DO
    END DO
    bdytmask(:, :) = ssmask(:, :)
    bdyumask(:, :) = 0._wp
    bdyvmask(:, :) = 0._wp
    DO ij = 1, jpjm1
      DO ii = 1, jpim1
        bdyumask(ii, ij) = bdytmask(ii, ij) * bdytmask(ii + 1, ij)
        bdyvmask(ii, ij) = bdytmask(ii, ij) * bdytmask(ii, ij + 1)
      END DO
    END DO
    CALL lbc_lnk_multi('bdyini', bdyumask, 'U', 1., bdyvmask, 'V', 1.)
    igrd = 1
    DO ib_bdy = 1, nb_bdy
      DO ib = 1, idx_bdy(ib_bdy) % nblenrim(igrd)
        bdytmask(idx_bdy(ib_bdy) % nbi(ib, igrd), idx_bdy(ib_bdy) % nbj(ib, igrd)) = 0._wp
      END DO
    END DO
    igrd = 2
    DO ib_bdy = 1, nb_bdy
      DO ib = 1, idx_bdy(ib_bdy) % nblenrim(igrd)
        bdyumask(idx_bdy(ib_bdy) % nbi(ib, igrd), idx_bdy(ib_bdy) % nbj(ib, igrd)) = 0._wp
      END DO
    END DO
    igrd = 3
    DO ib_bdy = 1, nb_bdy
      DO ib = 1, idx_bdy(ib_bdy) % nblenrim(igrd)
        bdyvmask(idx_bdy(ib_bdy) % nbi(ib, igrd), idx_bdy(ib_bdy) % nbj(ib, igrd)) = 0._wp
      END DO
    END DO
    zfmask(:, :) = 0
    DO ij = 2, jpjm1
      DO ii = 2, jpim1
        zfmask(ii, ij) = tmask(ii, ij, 1) * tmask(ii + 1, ij, 1) * tmask(ii, ij + 1, 1) * tmask(ii + 1, ij + 1, 1)
      END DO
    END DO
    CALL lbc_lnk('bdyini', zfmask, 'F', 1.)
    CALL lbc_lnk_multi('bdyini', bdyumask, 'U', 1., bdyvmask, 'V', 1., bdytmask, 'T', 1.)
    DO ib_bdy = 1, nb_bdy
      idx_bdy(ib_bdy) % flagu(:, :) = 0._wp
      idx_bdy(ib_bdy) % flagv(:, :) = 0._wp
      icount = 0
      DO igrd = 1, jpbgrd
        SELECT CASE (igrd)
        CASE (1)
          pmask => umask(:, :, 1)
          i_offset = 0
        CASE (2)
          pmask => bdytmask(:, :)
          i_offset = 1
        CASE (3)
          pmask => zfmask(:, :)
          i_offset = 0
        END SELECT
        icount = 0
        DO ib = 1, idx_bdy(ib_bdy) % nblenrim(igrd)
          nbi => idx_bdy(ib_bdy) % nbi(ib, igrd)
          nbj => idx_bdy(ib_bdy) % nbj(ib, igrd)
          zefl = pmask(nbi + i_offset - 1, nbj)
          zwfl = pmask(nbi + i_offset, nbj)
          IF (i_offset == 1 .AND. zefl + zwfl == 2) THEN
            icount = icount + 1
            IF (lwp) WRITE(numout, *) 'Problem with igrd = ', igrd, ' at (global) nbi, nbj : ', mig(nbi), mjg(nbj)
          ELSE
            idx_bdy(ib_bdy) % flagu(ib, igrd) = - zefl + zwfl
          END IF
        END DO
        IF (icount /= 0) THEN
          WRITE(ctmp1, *) ' E R R O R : Some ', cgrid(igrd), ' grid points,', ' are not boundary points (flagu calculation). Check nbi, nbj, indices for boundary set ', ib_bdy
          WRITE(ctmp2, *) ' ========== '
          CALL ctl_stop(' ', ctmp1, ctmp2, ' ')
        END IF
      END DO
      DO igrd = 1, jpbgrd
        SELECT CASE (igrd)
        CASE (1)
          pmask => vmask(:, :, 1)
          j_offset = 0
        CASE (2)
          pmask => zfmask(:, :)
          j_offset = 0
        CASE (3)
          pmask => bdytmask
          j_offset = 1
        END SELECT
        icount = 0
        DO ib = 1, idx_bdy(ib_bdy) % nblenrim(igrd)
          nbi => idx_bdy(ib_bdy) % nbi(ib, igrd)
          nbj => idx_bdy(ib_bdy) % nbj(ib, igrd)
          znfl = pmask(nbi, nbj + j_offset - 1)
          zsfl = pmask(nbi, nbj + j_offset)
          IF (j_offset == 1 .AND. znfl + zsfl == 2) THEN
            IF (lwp) WRITE(numout, *) 'Problem with igrd = ', igrd, ' at (global) nbi, nbj : ', mig(nbi), mjg(nbj)
            icount = icount + 1
          ELSE
            idx_bdy(ib_bdy) % flagv(ib, igrd) = - znfl + zsfl
          END IF
        END DO
        IF (icount /= 0) THEN
          WRITE(ctmp1, *) ' E R R O R : Some ', cgrid(igrd), ' grid points,', ' are not boundary points (flagv calculation). Check nbi, nbj, indices for boundary set ', ib_bdy
          WRITE(ctmp2, *) ' ========== '
          CALL ctl_stop(' ', ctmp1, ctmp2, ' ')
        END IF
      END DO
    END DO
    IF (nb_bdy > 0) DEALLOCATE(nbidta, nbjdta, nbrdta)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE bdy_segs
  SUBROUTINE bdy_ctl_seg
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ib, ib1, ib2, ji, jj, itest
    INTEGER, DIMENSION(jp_nseg, 2) :: icorne, icornw, icornn, icorns
    REAL(KIND = wp), DIMENSION(2) :: ztestmask
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bdy_ctl_seg', 'r0', 0, 0)
    IF (lwp) WRITE(numout, *) ' '
    IF (lwp) WRITE(numout, *) 'bdy_ctl_seg: Check analytical segments'
    IF (lwp) WRITE(numout, *) '~~~~~~~~~~~~'
    IF (lwp) WRITE(numout, *) 'Number of east  segments     : ', nbdysege
    IF (lwp) WRITE(numout, *) 'Number of west  segments     : ', nbdysegw
    IF (lwp) WRITE(numout, *) 'Number of north segments     : ', nbdysegn
    IF (lwp) WRITE(numout, *) 'Number of south segments     : ', nbdysegs
    DO ib = 1, nbdysegn
      IF (lwp) WRITE(numout, *) '**check north seg bounds pckg: ', npckgn(ib)
      IF ((jpjnob(ib) .GE. jpjglo - 1) .OR. (jpjnob(ib) .LE. 1)) CALL ctl_stop('nbdyind out of domain')
      IF (jpindt(ib) .GE. jpinft(ib)) CALL ctl_stop('Bdy start index is greater than end index')
      IF (jpindt(ib) .LE. 1) CALL ctl_stop('Start index out of domain')
      IF (jpinft(ib) .GE. jpiglo) CALL ctl_stop('End index out of domain')
    END DO
    DO ib = 1, nbdysegs
      IF (lwp) WRITE(numout, *) '**check south seg bounds pckg: ', npckgs(ib)
      IF ((jpjsob(ib) .GE. jpjglo - 1) .OR. (jpjsob(ib) .LE. 1)) CALL ctl_stop('nbdyind out of domain')
      IF (jpisdt(ib) .GE. jpisft(ib)) CALL ctl_stop('Bdy start index is greater than end index')
      IF (jpisdt(ib) .LE. 1) CALL ctl_stop('Start index out of domain')
      IF (jpisft(ib) .GE. jpiglo) CALL ctl_stop('End index out of domain')
    END DO
    DO ib = 1, nbdysege
      IF (lwp) WRITE(numout, *) '**check east  seg bounds pckg: ', npckge(ib)
      IF ((jpieob(ib) .GE. jpiglo - 1) .OR. (jpieob(ib) .LE. 1)) CALL ctl_stop('nbdyind out of domain')
      IF (jpjedt(ib) .GE. jpjeft(ib)) CALL ctl_stop('Bdy start index is greater than end index')
      IF (jpjedt(ib) .LE. 1) CALL ctl_stop('Start index out of domain')
      IF (jpjeft(ib) .GE. jpjglo) CALL ctl_stop('End index out of domain')
    END DO
    DO ib = 1, nbdysegw
      IF (lwp) WRITE(numout, *) '**check west  seg bounds pckg: ', npckgw(ib)
      IF ((jpiwob(ib) .GE. jpiglo - 1) .OR. (jpiwob(ib) .LE. 1)) CALL ctl_stop('nbdyind out of domain')
      IF (jpjwdt(ib) .GE. jpjwft(ib)) CALL ctl_stop('Bdy start index is greater than end index')
      IF (jpjwdt(ib) .LE. 1) CALL ctl_stop('Start index out of domain')
      IF (jpjwft(ib) .GE. jpjglo) CALL ctl_stop('End index out of domain')
    END DO
    IF (lwp) WRITE(numout, *) '**Look for segments corners  :'
    itest = 0
    icorne(:, :) = 0.
    icornw(:, :) = 0.
    icornn(:, :) = 0.
    icorns(:, :) = 0.
    IF ((nbdysegw > 0) .AND. (nbdysegs > 0)) THEN
      DO ib1 = 1, nbdysegw
        DO ib2 = 1, nbdysegs
          IF ((jpisdt(ib2) <= jpiwob(ib1)) .AND. (jpisft(ib2) >= jpiwob(ib1)) .AND. (jpjwdt(ib1) <= jpjsob(ib2)) .AND. (jpjwft(ib1) >= jpjsob(ib2))) THEN
            IF ((jpjwdt(ib1) == jpjsob(ib2)) .AND. (jpisdt(ib2) == jpiwob(ib1))) THEN
              icornw(ib1, 1) = npckgs(ib2)
              icorns(ib2, 1) = npckgw(ib1)
            ELSE IF ((jpisft(ib2) == jpiwob(ib1)) .AND. (jpjwft(ib1) == jpjsob(ib2))) THEN
              WRITE(ctmp1, *) ' E R R O R : Found an acute open boundary corner at point (i,j)= ', jpisft(ib2), jpjwft(ib1)
              WRITE(ctmp2, *) ' ==========  Not allowed yet'
              WRITE(ctmp3, *) '             Crossing problem with West segment: ', npckgw(ib1), ' and South segment: ', npckgs(ib2)
              CALL ctl_stop(' ', ctmp1, ctmp2, ctmp3, ' ')
            ELSE
              WRITE(ctmp1, *) ' E R R O R : Check South and West Open boundary indices'
              WRITE(ctmp2, *) ' ==========  Crossing problem with West segment: ', npckgw(ib1), ' and South segment: ', npckgs(ib2)
              CALL ctl_stop(' ', ctmp1, ctmp2, ' ')
            END IF
          END IF
        END DO
      END DO
    END IF
    IF ((nbdysege > 0) .AND. (nbdysegs > 0)) THEN
      DO ib1 = 1, nbdysege
        DO ib2 = 1, nbdysegs
          IF ((jpisdt(ib2) <= jpieob(ib1) + 1) .AND. (jpisft(ib2) >= jpieob(ib1) + 1) .AND. (jpjedt(ib1) <= jpjsob(ib2)) .AND. (jpjeft(ib1) >= jpjsob(ib2))) THEN
            IF ((jpjedt(ib1) == jpjsob(ib2)) .AND. (jpisft(ib2) == jpieob(ib1) + 1)) THEN
              icorne(ib1, 1) = npckgs(ib2)
              icorns(ib2, 2) = npckge(ib1)
            ELSE IF ((jpjeft(ib1) == jpjsob(ib2)) .AND. (jpisdt(ib2) == jpieob(ib1) + 1)) THEN
              WRITE(ctmp1, *) ' E R R O R : Found an acute open boundary corner at point (i,j)= ', jpisdt(ib2), jpjeft(ib1)
              WRITE(ctmp2, *) ' ==========  Not allowed yet'
              WRITE(ctmp3, *) '             Crossing problem with East segment: ', npckge(ib1), ' and South segment: ', npckgs(ib2)
              CALL ctl_stop(' ', ctmp1, ctmp2, ctmp3, ' ')
            ELSE
              WRITE(ctmp1, *) ' E R R O R : Check South and East Open boundary indices'
              WRITE(ctmp2, *) ' ==========  Crossing problem with East segment: ', npckge(ib1), ' and South segment: ', npckgs(ib2)
              CALL ctl_stop(' ', ctmp1, ctmp2, ' ')
            END IF
          END IF
        END DO
      END DO
    END IF
    IF ((nbdysegn > 0) .AND. (nbdysegw > 0)) THEN
      DO ib1 = 1, nbdysegw
        DO ib2 = 1, nbdysegn
          IF ((jpindt(ib2) <= jpiwob(ib1)) .AND. (jpinft(ib2) >= jpiwob(ib1)) .AND. (jpjwdt(ib1) <= jpjnob(ib2) + 1) .AND. (jpjwft(ib1) >= jpjnob(ib2) + 1)) THEN
            IF ((jpjwft(ib1) == jpjnob(ib2) + 1) .AND. (jpindt(ib2) == jpiwob(ib1))) THEN
              icornw(ib1, 2) = npckgn(ib2)
              icornn(ib2, 1) = npckgw(ib1)
            ELSE IF ((jpjwdt(ib1) == jpjnob(ib2) + 1) .AND. (jpinft(ib2) == jpiwob(ib1))) THEN
              WRITE(ctmp1, *) ' E R R O R : Found an acute open boundary corner at point (i,j)= ', jpinft(ib2), jpjwdt(ib1)
              WRITE(ctmp2, *) ' ==========  Not allowed yet'
              WRITE(ctmp3, *) '             Crossing problem with West segment: ', npckgw(ib1), ' and North segment: ', npckgn(ib2)
              CALL ctl_stop(' ', ctmp1, ctmp2, ctmp3, ' ')
            ELSE
              WRITE(ctmp1, *) ' E R R O R : Check North and West Open boundary indices'
              WRITE(ctmp2, *) ' ==========  Crossing problem with West segment: ', npckgw(ib1), ' and North segment: ', npckgn(ib2)
              CALL ctl_stop(' ', ctmp1, ctmp2, ' ')
            END IF
          END IF
        END DO
      END DO
    END IF
    IF ((nbdysegn > 0) .AND. (nbdysege > 0)) THEN
      DO ib1 = 1, nbdysege
        DO ib2 = 1, nbdysegn
          IF ((jpindt(ib2) <= jpieob(ib1) + 1) .AND. (jpinft(ib2) >= jpieob(ib1) + 1) .AND. (jpjedt(ib1) <= jpjnob(ib2) + 1) .AND. (jpjeft(ib1) >= jpjnob(ib2) + 1)) THEN
            IF ((jpjeft(ib1) == jpjnob(ib2) + 1) .AND. (jpinft(ib2) == jpieob(ib1) + 1)) THEN
              icorne(ib1, 2) = npckgn(ib2)
              icornn(ib2, 2) = npckge(ib1)
            ELSE IF ((jpjedt(ib1) == jpjnob(ib2) + 1) .AND. (jpindt(ib2) == jpieob(ib1) + 1)) THEN
              WRITE(ctmp1, *) ' E R R O R : Found an acute open boundary corner at point (i,j)= ', jpindt(ib2), jpjedt(ib1)
              WRITE(ctmp2, *) ' ==========  Not allowed yet'
              WRITE(ctmp3, *) '             Crossing problem with East segment: ', npckge(ib1), ' and North segment: ', npckgn(ib2)
              CALL ctl_stop(' ', ctmp1, ctmp2, ctmp3, ' ')
            ELSE
              WRITE(ctmp1, *) ' E R R O R : Check North and East Open boundary indices'
              WRITE(ctmp2, *) ' ==========  Crossing problem with East segment: ', npckge(ib1), ' and North segment: ', npckgn(ib2)
              CALL ctl_stop(' ', ctmp1, ctmp2, ' ')
            END IF
          END IF
        END DO
      END DO
    END IF
    DO ib = 1, nbdysegw
      ztestmask(1 : 2) = 0.
      DO ji = 1, jpi
        DO jj = 1, jpj
          IF (((ji + nimpp - 1) == jpiwob(ib)) .AND. ((jj + njmpp - 1) == jpjwdt(ib))) ztestmask(1) = tmask(ji, jj, 1)
          IF (((ji + nimpp - 1) == jpiwob(ib)) .AND. ((jj + njmpp - 1) == jpjwft(ib))) ztestmask(2) = tmask(ji, jj, 1)
        END DO
      END DO
      CALL mpp_sum('bdyini', ztestmask, 2)
      IF (ztestmask(1) == 1) THEN
        IF (icornw(ib, 1) == 0) THEN
          WRITE(ctmp1, *) ' E R R O R : Open boundary segment ', npckgw(ib)
          WRITE(ctmp2, *) ' ==========  does not start on land or on a corner'
          CALL ctl_stop(' ', ctmp1, ctmp2, ' ')
        ELSE
          IF (lwp) WRITE(numout, *) 'Found a South-West corner at (i,j): ', jpiwob(ib), jpjwdt(ib)
          CALL bdy_ctl_corn(npckgw(ib), icornw(ib, 1))
          itest = itest + 1
        END IF
      END IF
      IF (ztestmask(2) == 1) THEN
        IF (icornw(ib, 2) == 0) THEN
          WRITE(ctmp1, *) ' E R R O R : Open boundary segment ', npckgw(ib)
          WRITE(ctmp2, *) ' ==========  does not end on land or on a corner'
          CALL ctl_stop(' ', ctmp1, ctmp2, ' ')
        ELSE
          IF (lwp) WRITE(numout, *) 'Found a North-West corner at (i,j): ', jpiwob(ib), jpjwft(ib)
          CALL bdy_ctl_corn(npckgw(ib), icornw(ib, 2))
          itest = itest + 1
        END IF
      END IF
    END DO
    DO ib = 1, nbdysege
      ztestmask(1 : 2) = 0.
      DO ji = 1, jpi
        DO jj = 1, jpj
          IF (((ji + nimpp - 1) == jpieob(ib) + 1) .AND. ((jj + njmpp - 1) == jpjedt(ib))) ztestmask(1) = tmask(ji, jj, 1)
          IF (((ji + nimpp - 1) == jpieob(ib) + 1) .AND. ((jj + njmpp - 1) == jpjeft(ib))) ztestmask(2) = tmask(ji, jj, 1)
        END DO
      END DO
      CALL mpp_sum('bdyini', ztestmask, 2)
      IF (ztestmask(1) == 1) THEN
        IF (icorne(ib, 1) == 0) THEN
          WRITE(ctmp1, *) ' E R R O R : Open boundary segment ', npckge(ib)
          WRITE(ctmp2, *) ' ==========  does not start on land or on a corner'
          CALL ctl_stop(' ', ctmp1, ctmp2, ' ')
        ELSE
          IF (lwp) WRITE(numout, *) 'Found a South-East corner at (i,j): ', jpieob(ib) + 1, jpjedt(ib)
          CALL bdy_ctl_corn(npckge(ib), icorne(ib, 1))
          itest = itest + 1
        END IF
      END IF
      IF (ztestmask(2) == 1) THEN
        IF (icorne(ib, 2) == 0) THEN
          WRITE(ctmp1, *) ' E R R O R : Open boundary segment ', npckge(ib)
          WRITE(ctmp2, *) ' ==========  does not end on land or on a corner'
          CALL ctl_stop(' ', ctmp1, ctmp2, ' ')
        ELSE
          IF (lwp) WRITE(numout, *) 'Found a North-East corner at (i,j): ', jpieob(ib) + 1, jpjeft(ib)
          CALL bdy_ctl_corn(npckge(ib), icorne(ib, 2))
          itest = itest + 1
        END IF
      END IF
    END DO
    DO ib = 1, nbdysegs
      ztestmask(1 : 2) = 0.
      DO ji = 1, jpi
        DO jj = 1, jpj
          IF (((jj + njmpp - 1) == jpjsob(ib)) .AND. ((ji + nimpp - 1) == jpisdt(ib))) ztestmask(1) = tmask(ji, jj, 1)
          IF (((jj + njmpp - 1) == jpjsob(ib)) .AND. ((ji + nimpp - 1) == jpisft(ib))) ztestmask(2) = tmask(ji, jj, 1)
        END DO
      END DO
      CALL mpp_sum('bdyini', ztestmask, 2)
      IF ((ztestmask(1) == 1) .AND. (icorns(ib, 1) == 0)) THEN
        WRITE(ctmp1, *) ' E R R O R : Open boundary segment ', npckgs(ib)
        WRITE(ctmp2, *) ' ==========  does not start on land or on a corner'
        CALL ctl_stop(' ', ctmp1, ctmp2, ' ')
      END IF
      IF ((ztestmask(2) == 1) .AND. (icorns(ib, 2) == 0)) THEN
        WRITE(ctmp1, *) ' E R R O R : Open boundary segment ', npckgs(ib)
        WRITE(ctmp2, *) ' ==========  does not end on land or on a corner'
        CALL ctl_stop(' ', ctmp1, ctmp2, ' ')
      END IF
    END DO
    DO ib = 1, nbdysegn
      ztestmask(1 : 2) = 0.
      DO ji = 1, jpi
        DO jj = 1, jpj
          IF (((jj + njmpp - 1) == jpjnob(ib) + 1) .AND. ((ji + nimpp - 1) == jpindt(ib))) ztestmask(1) = tmask(ji, jj, 1)
          IF (((jj + njmpp - 1) == jpjnob(ib) + 1) .AND. ((ji + nimpp - 1) == jpinft(ib))) ztestmask(2) = tmask(ji, jj, 1)
        END DO
      END DO
      CALL mpp_sum('bdyini', ztestmask, 2)
      IF ((ztestmask(1) == 1) .AND. (icornn(ib, 1) == 0)) THEN
        WRITE(ctmp1, *) ' E R R O R : Open boundary segment ', npckgn(ib)
        WRITE(ctmp2, *) ' ==========  does not start on land'
        CALL ctl_stop(' ', ctmp1, ctmp2, ' ')
      END IF
      IF ((ztestmask(2) == 1) .AND. (icornn(ib, 2) == 0)) THEN
        WRITE(ctmp1, *) ' E R R O R : Open boundary segment ', npckgn(ib)
        WRITE(ctmp2, *) ' ==========  does not end on land'
        CALL ctl_stop(' ', ctmp1, ctmp2, ' ')
      END IF
    END DO
    IF ((itest == 0) .AND. (lwp)) WRITE(numout, *) 'NO open boundary corner found'
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE bdy_ctl_seg
  SUBROUTINE bdy_ctl_corn(ib1, ib2)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: ib1, ib2
    INTEGER :: itest
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bdy_ctl_corn', 'r0', 0, 0)
    itest = 0
    IF (cn_dyn2d(ib1) /= cn_dyn2d(ib2)) itest = itest + 1
    IF (cn_dyn3d(ib1) /= cn_dyn3d(ib2)) itest = itest + 1
    IF (cn_tra(ib1) /= cn_tra(ib2)) itest = itest + 1
    IF (nn_dyn2d_dta(ib1) /= nn_dyn2d_dta(ib2)) itest = itest + 1
    IF (nn_dyn3d_dta(ib1) /= nn_dyn3d_dta(ib2)) itest = itest + 1
    IF (nn_tra_dta(ib1) /= nn_tra_dta(ib2)) itest = itest + 1
    IF (nn_rimwidth(ib1) /= nn_rimwidth(ib2)) itest = itest + 1
    IF (itest > 0) THEN
      WRITE(ctmp1, *) ' E R R O R : Segments ', ib1, 'and ', ib2
      WRITE(ctmp2, *) ' ==========  have different open bdy schemes'
      CALL ctl_stop(' ', ctmp1, ctmp2, ' ')
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE bdy_ctl_corn
END MODULE bdyini