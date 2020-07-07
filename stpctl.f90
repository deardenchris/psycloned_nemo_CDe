MODULE stpctl
  USE oce
  USE dom_oce
  USE c1d
  USE diawri
  USE in_out_manager
  USE lbclnk
  USE lib_mpp
  USE wet_dry, ONLY: ll_wd, ssh_ref
  USE ieee_arithmetic
  USE netcdf
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: stp_ctl
  INTEGER :: idrun, idtime, idssh, idu, ids1, ids2, istatus
  CONTAINS
  SUBROUTINE stp_ctl(kt, kindic)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(INOUT) :: kindic
    INTEGER :: ji, jj, jk
    INTEGER :: iih, ijh
    INTEGER :: iiu, iju, iku
    INTEGER :: iis1, ijs1, iks1
    INTEGER :: iis2, ijs2, iks2
    REAL(KIND = wp) :: zzz
    INTEGER, DIMENSION(3) :: ilocu, ilocs1, ilocs2
    INTEGER, DIMENSION(2) :: iloch
    REAL(KIND = wp), DIMENSION(5) :: zmax
    CHARACTER(LEN = 20) :: clname
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('stp_ctl', 'r0', 0, 0)
    IF (kt == nit000 .AND. lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'stp_ctl : time-stepping control'
      WRITE(numout, FMT = *) '~~~~~~~'
      CALL ctl_opn(numstp, 'time.step', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, numout, lwp, narea)
      CALL ctl_opn(numrun, 'run.stat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, numout, lwp, narea)
      IF (lwm) THEN
        clname = 'run.stat.nc'
        IF (.NOT. Agrif_Root()) clname = TRIM(Agrif_CFixed()) // "_" // TRIM(clname)
        istatus = NF90_CREATE(TRIM(clname), NF90_CLOBBER, idrun)
        istatus = NF90_DEF_DIM(idrun, 'time', NF90_UNLIMITED, idtime)
        istatus = NF90_DEF_VAR(idrun, 'abs_ssh_max', NF90_DOUBLE, (/idtime/), idssh)
        istatus = NF90_DEF_VAR(idrun, 'abs_u_max', NF90_DOUBLE, (/idtime/), idu)
        istatus = NF90_DEF_VAR(idrun, 's_min', NF90_DOUBLE, (/idtime/), ids1)
        istatus = NF90_DEF_VAR(idrun, 's_max', NF90_DOUBLE, (/idtime/), ids2)
        istatus = NF90_ENDDEF(idrun)
      END IF
    END IF
    IF (lwp) THEN
      WRITE(numstp, FMT = '(1x, i8)') kt
      REWIND(UNIT = numstp)
    END IF
    IF (ll_wd) THEN
      !$ACC KERNELS ! CDe added
      zmax(1) = MAXVAL(ABS(sshn(:, :) + ssh_ref * tmask(:, :, 1)))
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS ! CDe added      
      zmax(1) = MAXVAL(ABS(sshn(:, :)))
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS !CDe added
    zmax(2) = MAXVAL(ABS(un(:, :, :)))
    zmax(3) = MAXVAL(- tsn(:, :, :, jp_sal), mask = tmask(:, :, :) == 1._wp)
    zmax(4) = MAXVAL(tsn(:, :, :, jp_sal), mask = tmask(:, :, :) == 1._wp)
    zmax(5) = REAL(nstop, wp)
    !$ACC END KERNELS
    IF (lk_mpp) THEN
      CALL mpp_max_multiple(zmax(:), 5)
      nstop = NINT(zmax(5))
    END IF
    IF (MOD(kt, nwrite) == 1 .AND. lwp) THEN
      WRITE(numout, FMT = *) ' ==>> time-step= ', kt, ' |ssh| max: ', zmax(1), ' |U| max: ', zmax(2), ' S min: ', - zmax(3), ' S max: ', zmax(4)
    END IF
    IF (zmax(1) > 15._wp .OR. zmax(2) > 10._wp .OR. zmax(3) >= 0._wp .OR. zmax(4) >= 100._wp .OR. zmax(4) < 0._wp .OR. ieee_is_nan(zmax(1) + zmax(2) + zmax(3))) THEN
      IF (lk_mpp) THEN
        CALL mpp_maxloc(ABS(sshn), ssmask(:, :), zzz, iih, ijh)
        CALL mpp_maxloc(ABS(un), umask(:, :, :), zzz, iiu, iju, iku)
        CALL mpp_minloc(tsn(:, :, :, jp_sal), tmask(:, :, :), zzz, iis1, ijs1, iks1)
        CALL mpp_maxloc(tsn(:, :, :, jp_sal), tmask(:, :, :), zzz, iis2, ijs2, iks2)
      ELSE
        iloch = MINLOC(ABS(sshn(:, :)))
        ilocu = MAXLOC(ABS(un(:, :, :)))
        ilocs1 = MINLOC(tsn(:, :, :, jp_sal), mask = tmask(:, :, :) == 1._wp)
        ilocs2 = MAXLOC(tsn(:, :, :, jp_sal), mask = tmask(:, :, :) == 1._wp)
        iih = iloch(1) + nimpp - 1
        ijh = iloch(2) + njmpp - 1
        iiu = ilocu(1) + nimpp - 1
        iju = ilocu(2) + njmpp - 1
        iku = ilocu(3)
        iis1 = ilocs1(1) + nimpp - 1
        ijs1 = ilocs1(2) + njmpp - 1
        iks1 = ilocs1(3)
        iis2 = ilocs2(1) + nimpp - 1
        ijs2 = ilocs2(2) + njmpp - 1
        iks2 = ilocs2(3)
      END IF
      IF (lwp) THEN
        WRITE(numout, cform_err)
        WRITE(numout, FMT = *) ' stp_ctl: |ssh| > 10 m  or  |U| > 10 m/s  or  S <= 0  or  S >= 100  or  NaN encounter in the tests'
        WRITE(numout, FMT = *) ' ======= '
        WRITE(numout, 9100) kt, zmax(1), iih, ijh
        WRITE(numout, 9200) kt, zmax(2), iiu, iju, iku
        WRITE(numout, 9300) kt, - zmax(3), iis1, ijs1, iks1
        WRITE(numout, 9400) kt, zmax(4), iis2, ijs2, iks2
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) '          output of last computed fields in output.abort.nc file'
      END IF
      kindic = - 3
      nstop = nstop + 1
      CALL dia_wri_state('output.abort', kt)
    END IF
9100 FORMAT(' kt=', I8, '   |ssh| max: ', 1P, G11.4, ', at  i j  : ', 2I5)
9200 FORMAT(' kt=', I8, '   |U|   max: ', 1P, G11.4, ', at  i j k: ', 3I5)
9300 FORMAT(' kt=', I8, '   S     min: ', 1P, G11.4, ', at  i j k: ', 3I5)
9400 FORMAT(' kt=', I8, '   S     max: ', 1P, G11.4, ', at  i j k: ', 3I5)
    IF (lwp) WRITE(numrun, 9500) kt, zmax(1), zmax(2), - zmax(3), zmax(4)
    IF (lwm) THEN
      istatus = NF90_PUT_VAR(idrun, idssh, (/zmax(1)/), (/kt/), (/1/))
      istatus = NF90_PUT_VAR(idrun, idu, (/zmax(2)/), (/kt/), (/1/))
      istatus = NF90_PUT_VAR(idrun, ids1, (/- zmax(3)/), (/kt/), (/1/))
      istatus = NF90_PUT_VAR(idrun, ids2, (/zmax(4)/), (/kt/), (/1/))
      IF (MOD(kt, 100) == 0) istatus = NF90_SYNC(idrun)
      IF (kt == nitend) istatus = NF90_CLOSE(idrun)
    END IF
9500 FORMAT(' it :', I8, '    |ssh|_max: ', D23.16, ' |U|_max: ', D23.16, ' S_min: ', D23.16, ' S_max: ', D23.16)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE stp_ctl
END MODULE stpctl
