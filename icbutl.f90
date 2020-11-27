MODULE icbutl
  USE par_oce
  USE dom_oce
  USE in_out_manager
  USE lbclnk
  USE lib_mpp
  USE icb_oce
  USE sbc_oce
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: icb_utl_copy
  PUBLIC :: icb_utl_interp
  PUBLIC :: icb_utl_bilin
  PUBLIC :: icb_utl_bilin_x
  PUBLIC :: icb_utl_add
  PUBLIC :: icb_utl_delete
  PUBLIC :: icb_utl_destroy
  PUBLIC :: icb_utl_track
  PUBLIC :: icb_utl_print_berg
  PUBLIC :: icb_utl_print
  PUBLIC :: icb_utl_count
  PUBLIC :: icb_utl_incr
  PUBLIC :: icb_utl_yearday
  PUBLIC :: icb_utl_mass
  PUBLIC :: icb_utl_heat
  CONTAINS
  SUBROUTINE icb_utl_copy
    !$ACC KERNELS ! CDe added      
    uo_e(1 : jpi, 1 : jpj) = ssu_m(:, :) * umask(:, :, 1)
    vo_e(1 : jpi, 1 : jpj) = ssv_m(:, :) * vmask(:, :, 1)
    ff_e(1 : jpi, 1 : jpj) = ff_f(:, :)
    tt_e(1 : jpi, 1 : jpj) = sst_m(:, :)
    fr_e(1 : jpi, 1 : jpj) = fr_i(:, :)
    ua_e(1 : jpi, 1 : jpj) = utau(:, :) * umask(:, :, 1)
    va_e(1 : jpi, 1 : jpj) = vtau(:, :) * vmask(:, :, 1)
    !$ACC END KERNELS
    CALL lbc_lnk_icb('icbutl', uo_e, 'U', - 1._wp, 1, 1)
    CALL lbc_lnk_icb('icbutl', vo_e, 'V', - 1._wp, 1, 1)
    CALL lbc_lnk_icb('icbutl', ff_e, 'F', + 1._wp, 1, 1)
    CALL lbc_lnk_icb('icbutl', ua_e, 'U', - 1._wp, 1, 1)
    CALL lbc_lnk_icb('icbutl', va_e, 'V', - 1._wp, 1, 1)
    CALL lbc_lnk_icb('icbutl', fr_e, 'T', + 1._wp, 1, 1)
    CALL lbc_lnk_icb('icbutl', tt_e, 'T', + 1._wp, 1, 1)
    !$ACC KERNELS
    ssh_e(1 : jpi, 1 : jpj) = ssh_m(:, :) * tmask(:, :, 1)
    !$ACC END KERNELS
    CALL lbc_lnk_icb('icbutl', ssh_e, 'T', + 1._wp, 1, 1)
  END SUBROUTINE icb_utl_copy
  SUBROUTINE icb_utl_interp(pi, pe1, puo, pui, pua, pssh_i, pj, pe2, pvo, pvi, pva, pssh_j, psst, pcn, phi, pff)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN) :: pi, pj
    REAL(KIND = wp), INTENT(OUT) :: pe1, pe2
    REAL(KIND = wp), INTENT(OUT) :: puo, pvo, pui, pvi, pua, pva
    REAL(KIND = wp), INTENT(OUT) :: pssh_i, pssh_j
    REAL(KIND = wp), INTENT(OUT) :: psst, pcn, phi, pff
    REAL(KIND = wp) :: zcd, zmod
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_utl_interp', 'r0', 0, 0)
    pe1 = icb_utl_bilin_e(e1t, e1u, e1v, e1f, pi, pj)
    pe2 = icb_utl_bilin_e(e2t, e2u, e2v, e2f, pi, pj)
    puo = icb_utl_bilin_h(uo_e, pi, pj, 'U', .FALSE.)
    pvo = icb_utl_bilin_h(vo_e, pi, pj, 'V', .FALSE.)
    psst = icb_utl_bilin_h(tt_e, pi, pj, 'T', .TRUE.)
    pcn = icb_utl_bilin_h(fr_e, pi, pj, 'T', .TRUE.)
    pff = icb_utl_bilin_h(ff_e, pi, pj, 'F', .FALSE.)
    pua = icb_utl_bilin_h(ua_e, pi, pj, 'U', .TRUE.)
    pva = icb_utl_bilin_h(va_e, pi, pj, 'V', .TRUE.)
    zcd = 1.22_wp * 1.5E-3_wp
    zmod = 1._wp / MAX(1.E-20, SQRT(zcd * SQRT(pua * pua + pva * pva)))
    pua = pua * zmod
    pva = pva * zmod
    pui = 0._wp
    pvi = 0._wp
    phi = 0._wp
    pssh_i = (icb_utl_bilin_h(ssh_e, pi + 0.1_wp, pj, 'T', .TRUE.) - icb_utl_bilin_h(ssh_e, pi - 0.1_wp, pj, 'T', .TRUE.)) / &
&(0.2_wp * pe1)
    pssh_j = (icb_utl_bilin_h(ssh_e, pi, pj + 0.1_wp, 'T', .TRUE.) - icb_utl_bilin_h(ssh_e, pi, pj - 0.1_wp, 'T', .TRUE.)) / &
&(0.2_wp * pe2)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_utl_interp
  REAL(KIND = wp) FUNCTION icb_utl_bilin_h(pfld, pi, pj, cd_type, plmask)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(0 : jpi + 1, 0 : jpj + 1), INTENT(IN) :: pfld
    REAL(KIND = wp), INTENT(IN) :: pi, pj
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    LOGICAL, INTENT(IN) :: plmask
    INTEGER :: ii, ij
    REAL(KIND = wp) :: zi, zj
    REAL(KIND = wp) :: zw1, zw2, zw3, zw4
    REAL(KIND = wp), DIMENSION(4) :: zmask
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wp', 'r0', 0, 0)
    SELECT CASE (cd_type)
    CASE ('T')
      ii = MAX(0, INT(pi))
      ij = MAX(0, INT(pj))
      zi = pi - REAL(ii, wp)
      zj = pj - REAL(ij, wp)
    CASE ('U')
      ii = MAX(0, INT(pi - 0.5_wp))
      ij = MAX(0, INT(pj))
      zi = pi - 0.5_wp - REAL(ii, wp)
      zj = pj - REAL(ij, wp)
    CASE ('V')
      ii = MAX(0, INT(pi))
      ij = MAX(0, INT(pj - 0.5_wp))
      zi = pi - REAL(ii, wp)
      zj = pj - 0.5_wp - REAL(ij, wp)
    CASE ('F')
      ii = MAX(0, INT(pi - 0.5_wp))
      ij = MAX(0, INT(pj - 0.5_wp))
      zi = pi - 0.5_wp - REAL(ii, wp)
      zj = pj - 0.5_wp - REAL(ij, wp)
    END SELECT
    IF (ii <= mig(1) - 1) THEN
      ii = 0
    ELSE IF (ii > mig(jpi)) THEN
      ii = jpi
    ELSE
      ii = mi1(ii)
    END IF
    IF (ij <= mjg(1) - 1) THEN
      ij = 0
    ELSE IF (ij > mjg(jpj)) THEN
      ij = jpj
    ELSE
      ij = mj1(ij)
    END IF
    IF (plmask) THEN
      SELECT CASE (cd_type)
      CASE ('T')
        zmask = (/tmask_e(ii, ij), tmask_e(ii + 1, ij), tmask_e(ii, ij + 1), tmask_e(ii + 1, ij + 1)/)
      CASE ('U')
        zmask = (/umask_e(ii, ij), umask_e(ii + 1, ij), umask_e(ii, ij + 1), umask_e(ii + 1, ij + 1)/)
      CASE ('V')
        zmask = (/vmask_e(ii, ij), vmask_e(ii + 1, ij), vmask_e(ii, ij + 1), vmask_e(ii + 1, ij + 1)/)
      CASE ('F')
        zmask = 1.
      END SELECT
    ELSE
      zmask = 1.
    END IF
    zw1 = zmask(1) * (1._wp - zi) * (1._wp - zj)
    zw2 = zmask(2) * zi * (1._wp - zj)
    zw3 = zmask(3) * (1._wp - zi) * zj
    zw4 = zmask(4) * zi * zj
    icb_utl_bilin_h = (pfld(ii, ij) * zw1 + pfld(ii + 1, ij) * zw2 + pfld(ii, ij + 1) * zw3 + pfld(ii + 1, ij + 1) * zw4) / &
&MAX(1.E-20, zw1 + zw2 + zw3 + zw4)
    CALL profile_psy_data0 % PostEnd
  END FUNCTION icb_utl_bilin_h
  REAL(KIND = wp) FUNCTION icb_utl_bilin(pfld, pi, pj, cd_type)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pfld
    REAL(KIND = wp), INTENT(IN) :: pi, pj
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    INTEGER :: ii, ij
    REAL(KIND = wp) :: zi, zj
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wp', 'r0', 0, 0)
    SELECT CASE (cd_type)
    CASE ('T')
      ii = MAX(1, INT(pi))
      ij = MAX(1, INT(pj))
      zi = pi - REAL(ii, wp)
      zj = pj - REAL(ij, wp)
    CASE ('U')
      ii = MAX(1, INT(pi - 0.5))
      ij = MAX(1, INT(pj))
      zi = pi - 0.5 - REAL(ii, wp)
      zj = pj - REAL(ij, wp)
    CASE ('V')
      ii = MAX(1, INT(pi))
      ij = MAX(1, INT(pj - 0.5))
      zi = pi - REAL(ii, wp)
      zj = pj - 0.5 - REAL(ij, wp)
    CASE ('F')
      ii = MAX(1, INT(pi - 0.5))
      ij = MAX(1, INT(pj - 0.5))
      zi = pi - 0.5 - REAL(ii, wp)
      zj = pj - 0.5 - REAL(ij, wp)
    END SELECT
    IF (ii < mig(1)) THEN
      ii = 1
    ELSE IF (ii > mig(jpi)) THEN
      ii = jpi
    ELSE
      ii = mi1(ii)
    END IF
    IF (ij < mjg(1)) THEN
      ij = 1
    ELSE IF (ij > mjg(jpj)) THEN
      ij = jpj
    ELSE
      ij = mj1(ij)
    END IF
    IF (ii == jpi) ii = ii - 1
    IF (ij == jpj) ij = ij - 1
    icb_utl_bilin = (pfld(ii, ij) * (1. - zi) + pfld(ii + 1, ij) * zi) * (1. - zj) + (pfld(ii, ij + 1) * (1. - zi) + pfld(ii + 1, &
&ij + 1) * zi) * zj
    CALL profile_psy_data0 % PostEnd
  END FUNCTION icb_utl_bilin
  REAL(KIND = wp) FUNCTION icb_utl_bilin_x(pfld, pi, pj)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pfld
    REAL(KIND = wp), INTENT(IN) :: pi, pj
    INTEGER :: ii, ij
    REAL(KIND = wp) :: zi, zj
    REAL(KIND = wp) :: zret
    REAL(KIND = wp), DIMENSION(4) :: z4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wp', 'r0', 0, 0)
    ii = MAX(1, INT(pi))
    ij = MAX(1, INT(pj))
    zi = pi - REAL(ii, wp)
    zj = pj - REAL(ij, wp)
    IF (ii < mig(1)) THEN
      ii = 1
    ELSE IF (ii > mig(jpi)) THEN
      ii = jpi
    ELSE
      ii = mi1(ii)
    END IF
    IF (ij < mjg(1)) THEN
      ij = 1
    ELSE IF (ij > mjg(jpj)) THEN
      ij = jpj
    ELSE
      ij = mj1(ij)
    END IF
    IF (ii == jpi) ii = ii - 1
    IF (ij == jpj) ij = ij - 1
    z4(1) = pfld(ii, ij)
    z4(2) = pfld(ii + 1, ij)
    z4(3) = pfld(ii, ij + 1)
    z4(4) = pfld(ii + 1, ij + 1)
    IF (MAXVAL(z4) - MINVAL(z4) > 90._wp) THEN
      WHERE (z4 < 0._wp) z4 = z4 + 360._wp
    END IF
    zret = (z4(1) * (1. - zi) + z4(2) * zi) * (1. - zj) + (z4(3) * (1. - zi) + z4(4) * zi) * zj
    IF (zret > 180._wp) zret = zret - 360._wp
    icb_utl_bilin_x = zret
    CALL profile_psy_data0 % PostEnd
  END FUNCTION icb_utl_bilin_x
  REAL(KIND = wp) FUNCTION icb_utl_bilin_e(pet, peu, pev, pef, pi, pj)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pet, peu, pev, pef
    REAL(KIND = wp), INTENT(IN) :: pi, pj
    INTEGER :: ii, ij, icase, ierr
    REAL(KIND = wp) :: zi, zj
    REAL(KIND = wp) :: ze00, ze10, ze01, ze11
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wp', 'r0', 0, 0)
    ii = MAX(1, INT(pi))
    ij = MAX(1, INT(pj))
    zi = pi - REAL(ii, wp)
    zj = pj - REAL(ij, wp)
    ierr = 0
    IF (ii < mig(1)) THEN
      ii = 1
      ierr = ierr + 1
    ELSE IF (ii > mig(jpi)) THEN
      ii = jpi
      ierr = ierr + 1
    ELSE
      ii = mi1(ii)
    END IF
    IF (ij < mjg(1)) THEN
      ij = 1
      ierr = ierr + 1
    ELSE IF (ij > mjg(jpj)) THEN
      ij = jpj
      ierr = ierr + 1
    ELSE
      ij = mj1(ij)
    END IF
    IF (ii == jpi) THEN
      ii = ii - 1
      ierr = ierr + 1
    END IF
    IF (ij == jpj) THEN
      ij = ij - 1
      ierr = ierr + 1
    END IF
    IF (ierr > 0) CALL ctl_stop('STOP', 'icb_utl_bilin_e: an icebergs coordinates is out of valid range (out of bound error)')
    IF (0.0_wp <= zi .AND. zi < 0.5_wp) THEN
      IF (0.0_wp <= zj .AND. zj < 0.5_wp) THEN
        ze01 = pev(ii, ij)
        ze11 = pef(ii, ij)
        ze00 = pet(ii, ij)
        ze10 = peu(ii, ij)
        zi = 2._wp * zi
        zj = 2._wp * zj
      ELSE
        ze01 = pet(ii, ij + 1)
        ze11 = peu(ii, ij + 1)
        ze00 = pev(ii, ij)
        ze10 = pef(ii, ij)
        zi = 2._wp * zi
        zj = 2._wp * (zj - 0.5_wp)
      END IF
    ELSE
      IF (0.0_wp <= zj .AND. zj < 0.5_wp) THEN
        ze01 = pef(ii, ij)
        ze11 = pev(ii + 1, ij)
        ze00 = peu(ii, ij)
        ze10 = pet(ii + 1, ij)
        zi = 2._wp * (zi - 0.5_wp)
        zj = 2._wp * zj
      ELSE
        ze01 = peu(ii, ij + 1)
        ze11 = pet(ii + 1, ij + 1)
        ze00 = pef(ii, ij)
        ze10 = pev(ii + 1, ij)
        zi = 2._wp * (zi - 0.5_wp)
        zj = 2._wp * (zj - 0.5_wp)
      END IF
    END IF
    icb_utl_bilin_e = (ze01 * (1._wp - zi) + ze11 * zi) * zj + (ze00 * (1._wp - zi) + ze10 * zi) * (1._wp - zj)
    CALL profile_psy_data0 % PostEnd
  END FUNCTION icb_utl_bilin_e
  SUBROUTINE icb_utl_add(bergvals, ptvals)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(iceberg), INTENT(IN) :: bergvals
    TYPE(point), INTENT(IN) :: ptvals
    TYPE(iceberg), POINTER :: new => NULL()
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_utl_add', 'r0', 0, 0)
    new => NULL()
    CALL icb_utl_create(new, bergvals, ptvals)
    CALL icb_utl_insert(new)
    new => NULL()
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_utl_add
  SUBROUTINE icb_utl_create(berg, bergvals, ptvals)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(iceberg), INTENT(IN) :: bergvals
    TYPE(point), INTENT(IN) :: ptvals
    TYPE(iceberg), POINTER :: berg
    TYPE(point), POINTER :: pt
    INTEGER :: istat
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_utl_create', 'r0', 0, 0)
    IF (ASSOCIATED(berg)) CALL ctl_stop('icebergs, icb_utl_create: berg already associated')
    ALLOCATE(berg, STAT = istat)
    IF (istat /= 0) CALL ctl_stop('failed to allocate iceberg')
    berg % number(:) = bergvals % number(:)
    berg % mass_scaling = bergvals % mass_scaling
    berg % prev => NULL()
    berg % next => NULL()
    ALLOCATE(pt, STAT = istat)
    IF (istat /= 0) CALL ctl_stop('failed to allocate first iceberg point')
    pt = ptvals
    berg % current_point => pt
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_utl_create
  SUBROUTINE icb_utl_insert(newberg)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(iceberg), POINTER :: newberg
    TYPE(iceberg), POINTER :: this, prev, last
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_utl_insert', 'r0', 0, 0)
    IF (ASSOCIATED(first_berg)) THEN
      last => first_berg
      DO WHILE (ASSOCIATED(last % next))
        last => last % next
      END DO
      newberg % prev => last
      last % next => newberg
      last => newberg
    ELSE
      first_berg => newberg
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_utl_insert
  REAL(KIND = wp) FUNCTION icb_utl_yearday(kmon, kday, khr, kmin, ksec)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kmon, kday, khr, kmin, ksec
    INTEGER, DIMENSION(12) :: imonths = (/0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30/)
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wp', 'r0', 0, 0)
    icb_utl_yearday = REAL(SUM(imonths(1 : kmon)), wp)
    icb_utl_yearday = icb_utl_yearday + REAL(kday - 1, wp) + (REAL(khr, wp) + (REAL(kmin, wp) + REAL(ksec, wp) / 60.) / 60.) / 24.
    CALL profile_psy_data0 % PostEnd
  END FUNCTION icb_utl_yearday
  SUBROUTINE icb_utl_delete(first, berg)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(iceberg), POINTER :: first, berg
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_utl_delete', 'r0', 0, 0)
    IF (ASSOCIATED(berg % prev)) THEN
      berg % prev % next => berg % next
    ELSE
      first => berg % next
    END IF
    IF (ASSOCIATED(berg % next)) berg % next % prev => berg % prev
    CALL icb_utl_destroy(berg)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_utl_delete
  SUBROUTINE icb_utl_destroy(berg)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(iceberg), POINTER :: berg
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_utl_destroy', 'r0', 0, 0)
    IF (ASSOCIATED(berg % current_point)) DEALLOCATE(berg % current_point)
    DEALLOCATE(berg)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_utl_destroy
  SUBROUTINE icb_utl_track(knum, cd_label, kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, DIMENSION(nkounts) :: knum
    CHARACTER(LEN = *) :: cd_label
    INTEGER :: kt
    TYPE(iceberg), POINTER :: this
    LOGICAL :: match
    INTEGER :: k
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_utl_track', 'r0', 0, 0)
    this => first_berg
    DO WHILE (ASSOCIATED(this))
      match = .TRUE.
      DO k = 1, nkounts
        IF (this % number(k) /= knum(k)) match = .FALSE.
      END DO
      IF (match) CALL icb_utl_print_berg(this, kt)
      this => this % next
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_utl_track
  SUBROUTINE icb_utl_print_berg(berg, kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(iceberg), POINTER :: berg
    TYPE(point), POINTER :: pt
    INTEGER :: kt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    IF (nn_verbose_level == 0) RETURN
    CALL profile_psy_data0 % PreStart('icb_utl_print_berg', 'r0', 0, 0)
    pt => berg % current_point
    WRITE(numicb, 9200) kt, berg % number(1), pt % xi, pt % yj, pt % lon, pt % lat, pt % uvel, pt % vvel, pt % uo, pt % vo, pt % &
&ua, pt % va, pt % ui, pt % vi
    CALL flush(numicb)
9200 FORMAT(5X, I5, 2X, I10, 6(2X, 2F10.4))
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_utl_print_berg
  SUBROUTINE icb_utl_print(cd_label, kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *) :: cd_label
    INTEGER :: kt
    INTEGER :: ibergs, inbergs
    TYPE(iceberg), POINTER :: this
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    IF (nn_verbose_level == 0) RETURN
    CALL profile_psy_data0 % PreStart('icb_utl_print', 'r0', 0, 0)
    this => first_berg
    IF (ASSOCIATED(this)) THEN
      WRITE(numicb, '(a," pe=(",i3,")")') cd_label, narea
      WRITE(numicb, FMT = '(a8,4x,a6,12x,a5,15x,a7,19x,a3,17x,a5,17x,a5,17x,a5)') 'timestep', 'number', 'xi,yj', 'lon,lat', 'u,v', &
&'uo,vo', 'ua,va', 'ui,vi'
    END IF
    DO WHILE (ASSOCIATED(this))
      CALL icb_utl_print_berg(this, kt)
      this => this % next
    END DO
    ibergs = icb_utl_count()
    inbergs = ibergs
    CALL mpp_sum('icbutl', inbergs)
    IF (ibergs > 0) WRITE(numicb, FMT = '(a," there are",i5," bergs out of",i6," on PE ",i4)') cd_label, ibergs, inbergs, narea
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_utl_print
  SUBROUTINE icb_utl_incr
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ii, ibig
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_utl_incr', 'r0', 0, 0)
    ibig = HUGE(num_bergs(1))
    IF (ibig - jpnij < num_bergs(1)) THEN
      num_bergs(1) = narea
      DO ii = 2, nkounts
        IF (num_bergs(ii) == ibig) THEN
          num_bergs(ii) = 0
          IF (ii == nkounts) CALL ctl_stop('Sorry, run out of iceberg number space')
        ELSE
          num_bergs(ii) = num_bergs(ii) + 1
          EXIT
        END IF
      END DO
    ELSE
      num_bergs(1) = num_bergs(1) + jpnij
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_utl_incr
  INTEGER FUNCTION icb_utl_count()
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(iceberg), POINTER :: this
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_utl_count', 'r0', 0, 0)
    icb_utl_count = 0
    this => first_berg
    DO WHILE (ASSOCIATED(this))
      icb_utl_count = icb_utl_count + 1
      this => this % next
    END DO
    CALL profile_psy_data0 % PostEnd
  END FUNCTION icb_utl_count
  REAL(KIND = wp) FUNCTION icb_utl_mass(first, justbits, justbergs)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(iceberg), POINTER :: first
    TYPE(point), POINTER :: pt
    LOGICAL, INTENT(IN), OPTIONAL :: justbits, justbergs
    TYPE(iceberg), POINTER :: this
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wp', 'r0', 0, 0)
    icb_utl_mass = 0._wp
    this => first
    IF (PRESENT(justbergs)) THEN
      DO WHILE (ASSOCIATED(this))
        pt => this % current_point
        icb_utl_mass = icb_utl_mass + pt % mass * this % mass_scaling
        this => this % next
      END DO
    ELSE IF (PRESENT(justbits)) THEN
      DO WHILE (ASSOCIATED(this))
        pt => this % current_point
        icb_utl_mass = icb_utl_mass + pt % mass_of_bits * this % mass_scaling
        this => this % next
      END DO
    ELSE
      DO WHILE (ASSOCIATED(this))
        pt => this % current_point
        icb_utl_mass = icb_utl_mass + (pt % mass + pt % mass_of_bits) * this % mass_scaling
        this => this % next
      END DO
    END IF
    CALL profile_psy_data0 % PostEnd
  END FUNCTION icb_utl_mass
  REAL(KIND = wp) FUNCTION icb_utl_heat(first, justbits, justbergs)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(iceberg), POINTER :: first
    LOGICAL, INTENT(IN), OPTIONAL :: justbits, justbergs
    TYPE(iceberg), POINTER :: this
    TYPE(point), POINTER :: pt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wp', 'r0', 0, 0)
    icb_utl_heat = 0._wp
    this => first
    IF (PRESENT(justbergs)) THEN
      DO WHILE (ASSOCIATED(this))
        pt => this % current_point
        icb_utl_heat = icb_utl_heat + pt % mass * this % mass_scaling * pt % heat_density
        this => this % next
      END DO
    ELSE IF (PRESENT(justbits)) THEN
      DO WHILE (ASSOCIATED(this))
        pt => this % current_point
        icb_utl_heat = icb_utl_heat + pt % mass_of_bits * this % mass_scaling * pt % heat_density
        this => this % next
      END DO
    ELSE
      DO WHILE (ASSOCIATED(this))
        pt => this % current_point
        icb_utl_heat = icb_utl_heat + (pt % mass + pt % mass_of_bits) * this % mass_scaling * pt % heat_density
        this => this % next
      END DO
    END IF
    CALL profile_psy_data0 % PostEnd
  END FUNCTION icb_utl_heat
END MODULE icbutl
