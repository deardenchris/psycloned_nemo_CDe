MODULE trdvor
  USE oce
  USE dom_oce
  USE trd_oce
  USE zdf_oce
  USE sbc_oce
  USE phycst
  USE ldfdyn
  USE dianam
  USE zdfmxl
  USE in_out_manager
  USE ioipsl
  USE lbclnk
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  INTERFACE trd_vor_zint
    MODULE PROCEDURE trd_vor_zint_2d, trd_vor_zint_3d
  END INTERFACE
  PUBLIC :: trd_vor
  PUBLIC :: trd_vor_init
  PUBLIC :: trd_vor_alloc
  INTEGER :: nh_t, nmoydpvor, nidvor, nhoridvor, ndimvor1, icount
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndexvor1
  INTEGER :: ndebug
  REAL(KIND = wp), SAVE, ALLOCATABLE, DIMENSION(:, :) :: vor_avr
  REAL(KIND = wp), SAVE, ALLOCATABLE, DIMENSION(:, :) :: vor_avrb
  REAL(KIND = wp), SAVE, ALLOCATABLE, DIMENSION(:, :) :: vor_avrbb
  REAL(KIND = wp), SAVE, ALLOCATABLE, DIMENSION(:, :) :: vor_avrbn
  REAL(KIND = wp), SAVE, ALLOCATABLE, DIMENSION(:, :) :: rotot
  REAL(KIND = wp), SAVE, ALLOCATABLE, DIMENSION(:, :) :: vor_avrtot
  REAL(KIND = wp), SAVE, ALLOCATABLE, DIMENSION(:, :) :: vor_avrres
  REAL(KIND = wp), SAVE, ALLOCATABLE, DIMENSION(:, :, :) :: vortrd
  CHARACTER(LEN = 12) :: cvort
  CONTAINS
  INTEGER FUNCTION trd_vor_alloc()
    ALLOCATE(vor_avr(jpi, jpj), vor_avrb(jpi, jpj), vor_avrbb(jpi, jpj), vor_avrbn(jpi, jpj), rotot(jpi, jpj), vor_avrtot(jpi, &
&jpj), vor_avrres(jpi, jpj), vortrd(jpi, jpj, jpltot_vor), ndexvor1(jpi * jpj), STAT = trd_vor_alloc)
    CALL mpp_sum('trdvor', trd_vor_alloc)
    IF (trd_vor_alloc /= 0) CALL ctl_stop('STOP', 'trd_vor_alloc: failed to allocate arrays')
  END FUNCTION trd_vor_alloc
  SUBROUTINE trd_vor(putrd, pvtrd, ktrd, kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: putrd, pvtrd
    INTEGER, INTENT(IN) :: ktrd
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztswu, ztswv
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    SELECT CASE (ktrd)
    CASE (jpdyn_hpg)
      CALL trd_vor_zint(putrd, pvtrd, jpvor_prg)
    CASE (jpdyn_keg)
      CALL trd_vor_zint(putrd, pvtrd, jpvor_keg)
    CASE (jpdyn_rvo)
      CALL trd_vor_zint(putrd, pvtrd, jpvor_rvo)
    CASE (jpdyn_pvo)
      CALL trd_vor_zint(putrd, pvtrd, jpvor_pvo)
    CASE (jpdyn_ldf)
      CALL trd_vor_zint(putrd, pvtrd, jpvor_ldf)
    CASE (jpdyn_zad)
      CALL trd_vor_zint(putrd, pvtrd, jpvor_zad)
    CASE (jpdyn_spg)
      CALL trd_vor_zint(putrd, pvtrd, jpvor_spg)
    CASE (jpdyn_zdf)
      !$ACC KERNELS
      ztswu(:, :) = 0.E0
      ztswv(:, :) = 0.E0
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ztswu(ji, jj) = 0.5 * (utau_b(ji, jj) + utau(ji, jj)) / (e3u_n(ji, jj, 1) * rau0)
          ztswv(ji, jj) = 0.5 * (vtau_b(ji, jj) + vtau(ji, jj)) / (e3v_n(ji, jj, 1) * rau0)
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data0 % PreStart('trd_vor', 'r0', 0, 0)
      CALL trd_vor_zint(putrd, pvtrd, jpvor_zdf)
      CALL trd_vor_zint(ztswu, ztswv, jpvor_swf)
      CALL profile_psy_data0 % PostEnd
    CASE (jpdyn_bfr)
      CALL profile_psy_data1 % PreStart('trd_vor', 'r1', 0, 0)
      CALL trd_vor_zint(putrd, pvtrd, jpvor_bfr)
      CALL profile_psy_data1 % PostEnd
    CASE (jpdyn_atf)
      CALL trd_vor_iom(kt)
    END SELECT
  END SUBROUTINE trd_vor
  SUBROUTINE trd_vor_zint_2d(putrdvor, pvtrdvor, ktrd)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: ktrd
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: putrdvor
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: pvtrdvor
    INTEGER :: ji, jj
    INTEGER :: ikbu, ikbv
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zudpvor, zvdpvor
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    zudpvor(:, :) = 0._wp
    zvdpvor(:, :) = 0._wp
    !$ACC END KERNELS
    CALL lbc_lnk_multi('trdvor', putrdvor, 'U', - 1., pvtrdvor, 'V', - 1.)
    !$ACC KERNELS
    SELECT CASE (ktrd)
    CASE (jpvor_bfr)
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ikbu = mbkv(ji, jj)
          ikbv = mbkv(ji, jj)
          zudpvor(ji, jj) = putrdvor(ji, jj) * e3u_n(ji, jj, ikbu) * e1u(ji, jj) * umask(ji, jj, ikbu)
          zvdpvor(ji, jj) = pvtrdvor(ji, jj) * e3v_n(ji, jj, ikbv) * e2v(ji, jj) * vmask(ji, jj, ikbv)
        END DO
      END DO
    CASE (jpvor_swf)
      zudpvor(:, :) = putrdvor(:, :) * e3u_n(:, :, 1) * e1u(:, :) * umask(:, :, 1)
      zvdpvor(:, :) = pvtrdvor(:, :) * e3v_n(:, :, 1) * e2v(:, :) * vmask(:, :, 1)
    END SELECT
    zudpvor(:, :) = zudpvor(:, :) * r1_hu_n(:, :)
    zvdpvor(:, :) = zvdpvor(:, :) * r1_hv_n(:, :)
    DO ji = 1, jpim1
      DO jj = 1, jpjm1
        vortrd(ji, jj, ktrd) = (zvdpvor(ji + 1, jj) - zvdpvor(ji, jj) - (zudpvor(ji, jj + 1) - zudpvor(ji, jj))) / (e1f(ji, jj) * &
&e2f(ji, jj))
      END DO
    END DO
    vortrd(:, :, ktrd) = vortrd(:, :, ktrd) * fmask(:, :, 1)
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('trd_vor_zint_2d', 'r0', 0, 0)
    IF (ndebug /= 0) THEN
      IF (lwp) WRITE(numout, FMT = *) ' debuging trd_vor_zint: I done'
      CALL FLUSH(numout)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE trd_vor_zint_2d
  SUBROUTINE trd_vor_zint_3d(putrdvor, pvtrdvor, ktrd)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: ktrd
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: putrdvor
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pvtrdvor
    INTEGER :: ji, jj, jk
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zubet, zvbet
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zudpvor, zvdpvor
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    zubet(:, :) = 0._wp
    zvbet(:, :) = 0._wp
    zudpvor(:, :) = 0._wp
    zvdpvor(:, :) = 0._wp
    !$ACC END KERNELS
    CALL lbc_lnk_multi('trdvor', putrdvor, 'U', - 1., pvtrdvor, 'V', - 1.)
    DO jk = 1, jpk
      !$ACC KERNELS
      zudpvor(:, :) = zudpvor(:, :) + putrdvor(:, :, jk) * e3u_n(:, :, jk) * e1u(:, :) * umask(:, :, jk)
      zvdpvor(:, :) = zvdpvor(:, :) + pvtrdvor(:, :, jk) * e3v_n(:, :, jk) * e2v(:, :) * vmask(:, :, jk)
      !$ACC END KERNELS
    END DO
    !$ACC KERNELS
    IF (ktrd == jpvor_pvo) THEN
      zubet(:, :) = zudpvor(:, :)
      zvbet(:, :) = zvdpvor(:, :)
      DO ji = 1, jpim1
        DO jj = 1, jpjm1
          vortrd(ji, jj, jpvor_bev) = (zvbet(ji + 1, jj) - zvbet(ji, jj) - (zubet(ji, jj + 1) - zubet(ji, jj))) / (e1f(ji, jj) * &
&e2f(ji, jj))
        END DO
      END DO
      vortrd(:, :, jpvor_bev) = vortrd(:, :, jpvor_bev) * r1_hu_n(:, :) * fmask(:, :, 1)
    END IF
    zudpvor(:, :) = zudpvor(:, :) * r1_hu_n(:, :)
    zvdpvor(:, :) = zvdpvor(:, :) * r1_hv_n(:, :)
    DO ji = 1, jpim1
      DO jj = 1, jpjm1
        vortrd(ji, jj, ktrd) = (zvdpvor(ji + 1, jj) - zvdpvor(ji, jj) - (zudpvor(ji, jj + 1) - zudpvor(ji, jj))) / (e1f(ji, jj) * &
&e2f(ji, jj))
      END DO
    END DO
    vortrd(:, :, ktrd) = vortrd(:, :, ktrd) * fmask(:, :, 1)
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('trd_vor_zint_3d', 'r0', 0, 0)
    IF (ndebug /= 0) THEN
      IF (lwp) WRITE(numout, FMT = *) ' debuging trd_vor_zint: I done'
      CALL FLUSH(numout)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE trd_vor_zint_3d
  SUBROUTINE trd_vor_iom(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk, jl
    INTEGER :: it, itmod
    REAL(KIND = wp) :: zmean
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zun, zvn
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    IF (kt > nit000) vor_avrb(:, :) = vor_avr(:, :)
    vor_avr(:, :) = 0._wp
    zun(:, :) = 0._wp
    zvn(:, :) = 0._wp
    vor_avrtot(:, :) = 0._wp
    vor_avrres(:, :) = 0._wp
    !$ACC END KERNELS
    DO jk = 1, jpk - 1
      !$ACC KERNELS
      zun(:, :) = zun(:, :) + e1u(:, :) * un(:, :, jk) * e3u_n(:, :, jk)
      zvn(:, :) = zvn(:, :) + e2v(:, :) * vn(:, :, jk) * e3v_n(:, :, jk)
      !$ACC END KERNELS
    END DO
    !$ACC KERNELS
    zun(:, :) = zun(:, :) * r1_hu_n(:, :)
    zvn(:, :) = zvn(:, :) * r1_hv_n(:, :)
    DO ji = 1, jpim1
      DO jj = 1, jpjm1
        vor_avr(ji, jj) = ((zvn(ji + 1, jj) - zvn(ji, jj)) - (zun(ji, jj + 1) - zun(ji, jj))) / (e1f(ji, jj) * e2f(ji, jj)) * &
&fmask(ji, jj, 1)
      END DO
    END DO
    IF (kt == nit000 + 1) THEN
      vor_avrbb(:, :) = vor_avrb(:, :)
      vor_avrbn(:, :) = vor_avr(:, :)
    END IF
    IF (kt >= nit000 + 2) THEN
      nmoydpvor = nmoydpvor + 1
      DO jl = 1, jpltot_vor
        IF (jl /= 9) THEN
          rotot(:, :) = rotot(:, :) + vortrd(:, :, jl)
        END IF
      END DO
    END IF
    it = kt
    itmod = kt - nit000 + 1
    !$ACC END KERNELS
    IF (MOD(it, nn_trd) == 0) THEN
      !$ACC KERNELS
      zmean = 1._wp / (REAL(nmoydpvor, wp) * 2._wp * rdt)
      vor_avrtot(:, :) = (vor_avr(:, :) - vor_avrbn(:, :) + vor_avrb(:, :) - vor_avrbb(:, :)) * zmean
      zmean = 1._wp / REAL(nmoydpvor, wp)
      vor_avrres(:, :) = vor_avrtot(:, :) - rotot(:, :) / zmean
      !$ACC END KERNELS
      CALL lbc_lnk_multi('trdvor', vor_avrtot, 'F', 1., vor_avrres, 'F', 1.)
      !$ACC KERNELS
      vor_avrbb(:, :) = vor_avrb(:, :)
      vor_avrbn(:, :) = vor_avr(:, :)
      nmoydpvor = 0
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data0 % PreStart('trd_vor_iom', 'r0', 0, 0)
    IF (kt >= nit000 + 1) THEN
      IF (lwp .AND. MOD(itmod, nn_trd) == 0) THEN
        WRITE(numout, FMT = *) ''
        WRITE(numout, FMT = *) 'trd_vor : write trends in the NetCDF file at kt = ', kt
        WRITE(numout, FMT = *) '~~~~~~~  '
      END IF
      CALL histwrite(nidvor, "sovortPh", it, vortrd(:, :, jpvor_prg), ndimvor1, ndexvor1)
      CALL histwrite(nidvor, "sovortEk", it, vortrd(:, :, jpvor_keg), ndimvor1, ndexvor1)
      CALL histwrite(nidvor, "sovozeta", it, vortrd(:, :, jpvor_rvo), ndimvor1, ndexvor1)
      CALL histwrite(nidvor, "sovortif", it, vortrd(:, :, jpvor_pvo), ndimvor1, ndexvor1)
      CALL histwrite(nidvor, "sovodifl", it, vortrd(:, :, jpvor_ldf), ndimvor1, ndexvor1)
      CALL histwrite(nidvor, "sovoadvv", it, vortrd(:, :, jpvor_zad), ndimvor1, ndexvor1)
      CALL histwrite(nidvor, "sovodifv", it, vortrd(:, :, jpvor_zdf), ndimvor1, ndexvor1)
      CALL histwrite(nidvor, "sovortPs", it, vortrd(:, :, jpvor_spg), ndimvor1, ndexvor1)
      CALL histwrite(nidvor, "sovortbv", it, vortrd(:, :, jpvor_bev), ndimvor1, ndexvor1)
      CALL histwrite(nidvor, "sovowind", it, vortrd(:, :, jpvor_swf), ndimvor1, ndexvor1)
      CALL histwrite(nidvor, "sovobfri", it, vortrd(:, :, jpvor_bfr), ndimvor1, ndexvor1)
      CALL histwrite(nidvor, "1st_mbre", it, vor_avrtot, ndimvor1, ndexvor1)
      CALL histwrite(nidvor, "sovorgap", it, vor_avrres, ndimvor1, ndexvor1)
      IF (ndebug /= 0) THEN
        WRITE(numout, FMT = *) ' debuging trd_vor: III.4 done'
        CALL FLUSH(numout)
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    IF (MOD(it, nn_trd) == 0) rotot(:, :) = 0
    !$ACC END KERNELS
    IF (kt == nitend) CALL histclo(nidvor)
  END SUBROUTINE trd_vor_iom
  SUBROUTINE trd_vor_init
    REAL(KIND = wp) :: zjulian, zsto, zout
    CHARACTER(LEN = 40) :: clhstnam
    CHARACTER(LEN = 40) :: clop
    cvort = 'averaged-vor'
    ndebug = 0
    IF (lwp) THEN
      WRITE(numout, FMT = *) ' '
      WRITE(numout, FMT = *) ' trd_vor_init: vorticity trends'
      WRITE(numout, FMT = *) ' ~~~~~~~~~~~~'
      WRITE(numout, FMT = *) ' '
      WRITE(numout, FMT = *) '               ##########################################################################'
      WRITE(numout, FMT = *) '                CAUTION: The interpretation of the vorticity trends is'
      WRITE(numout, FMT = *) '                not obvious, please contact Anne-Marie TREGUIER at: treguier@ifremer.fr '
      WRITE(numout, FMT = *) '               ##########################################################################'
      WRITE(numout, FMT = *) ' '
    END IF
    IF (trd_vor_alloc() /= 0) CALL ctl_stop('STOP', 'trd_vor_init : unable to allocate trdvor arrays')
    !$ACC KERNELS
    nmoydpvor = 0
    rotot(:, :) = 0
    vor_avrtot(:, :) = 0
    vor_avrres(:, :) = 0
    !$ACC END KERNELS
    IF (ndebug /= 0) THEN
      WRITE(numout, FMT = *) ' debuging trd_vor_init: I. done'
      CALL FLUSH(numout)
    END IF
    IF (ln_mskland) THEN
      clop = "only(x)"
    ELSE
      clop = "x"
    END IF
    zsto = rdt
    clop = "ave(" // TRIM(clop) // ")"
    zout = nn_trd * rdt
    IF (lwp) WRITE(numout, FMT = *) '               netCDF initialization'
    CALL ymds2ju(nyear, nmonth, nday, rdt, zjulian)
    zjulian = zjulian - adatrj
    IF (lwp) WRITE(numout, FMT = *) ' '
    IF (lwp) WRITE(numout, FMT = *) '               Date 0 used :', nit000, ' YEAR ', nyear, ' MONTH ', nmonth, ' DAY ', nday, &
&'Julian day : ', zjulian
    CALL dia_nam(clhstnam, nn_trd, 'vort')
    IF (lwp) WRITE(numout, FMT = *) ' Name of NETCDF file ', clhstnam
    CALL histbeg(clhstnam, jpi, glamf, jpj, gphif, 1, jpi, 1, jpj, nit000 - 1, zjulian, rdt, nh_t, nidvor, domain_id = nidom, &
&snc4chunks = snc4set)
    CALL wheneq(jpi * jpj, fmask, 1, 1., ndexvor1, ndimvor1)
    CALL histdef(nidvor, "sovortPh", cvort // "grad Ph", "s-2", jpi, jpj, nh_t, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(nidvor, "sovortEk", cvort // "Energy", "s-2", jpi, jpj, nh_t, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(nidvor, "sovozeta", cvort // "rel vorticity", "s-2", jpi, jpj, nh_t, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(nidvor, "sovortif", cvort // "coriolis", "s-2", jpi, jpj, nh_t, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(nidvor, "sovodifl", cvort // "lat diff ", "s-2", jpi, jpj, nh_t, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(nidvor, "sovoadvv", cvort // "vert adv", "s-2", jpi, jpj, nh_t, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(nidvor, "sovodifv", cvort // "vert diff", "s-2", jpi, jpj, nh_t, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(nidvor, "sovortPs", cvort // "grad Ps", "s-2", jpi, jpj, nh_t, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(nidvor, "sovortbv", cvort // "Beta V", "s-2", jpi, jpj, nh_t, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(nidvor, "sovowind", cvort // "wind stress", "s-2", jpi, jpj, nh_t, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(nidvor, "sovobfri", cvort // "bottom friction", "s-2", jpi, jpj, nh_t, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(nidvor, "1st_mbre", cvort // "1st mbre", "s-2", jpi, jpj, nh_t, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(nidvor, "sovorgap", cvort // "gap", "s-2", jpi, jpj, nh_t, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histend(nidvor, snc4set)
    IF (ndebug /= 0) THEN
      WRITE(numout, FMT = *) ' debuging trd_vor_init: II. done'
      CALL FLUSH(numout)
    END IF
  END SUBROUTINE trd_vor_init
END MODULE trdvor