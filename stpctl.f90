MODULE stpctl
  USE oce
  USE dom_oce
  USE c1d
  USE diawri
  USE in_out_manager
  USE lbclnk
  USE lib_mpp
  USE zdf_oce, ONLY: ln_zad_Aimp
  USE wet_dry, ONLY: ll_wd, ssh_ref
  USE ieee_arithmetic
  USE netcdf
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: stp_ctl
  INTEGER :: idrun, idtime, idssh, idu, ids1, ids2, idt1, idt2, idc1, idw1, istatus
  LOGICAL :: lsomeoce
  CONTAINS
  SUBROUTINE stp_ctl(kt, kindic)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(INOUT) :: kindic
    INTEGER :: ji, jj, jk
    INTEGER, DIMENSION(2) :: ih
    INTEGER, DIMENSION(3) :: iu, is1, is2
    REAL(KIND = wp) :: zzz
    !REAL(KIND = wp), DIMENSION(9) :: zmax
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:) :: zmax ! CDe
    LOGICAL :: ll_wrtstp, ll_colruns, ll_wrtruns
    CHARACTER(LEN = 20) :: clname
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    ALLOCATE(zmax(9)) ! CDe
    CALL profile_psy_data0 % PreStart('stp_ctl', 'r0', 0, 0)
    ll_wrtstp = (MOD(kt, sn_cfctl % ptimincr) == 0) .OR. (kt == nitend)
    ll_colruns = ll_wrtstp .AND. (ln_ctl .OR. sn_cfctl % l_runstat)
    ll_wrtruns = ll_colruns .AND. lwm
    CALL profile_psy_data0 % PostEnd
    IF (kt == nit000 .AND. lwp) THEN
      CALL profile_psy_data1 % PreStart('stp_ctl', 'r1', 0, 0)
      WRITE(numout, *)
      WRITE(numout, *) 'stp_ctl : time-stepping control'
      WRITE(numout, *) '~~~~~~~'
      IF (lwm) CALL ctl_opn(numstp, 'time.step', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, numout, lwp, narea)
      CALL profile_psy_data1 % PostEnd
      IF (lwm .AND. (ln_ctl .OR. sn_cfctl % l_runstat)) THEN
        CALL profile_psy_data2 % PreStart('stp_ctl', 'r2', 0, 0)
        CALL ctl_opn(numrun, 'run.stat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, numout, lwp, narea)
        clname = 'run.stat.nc'
        IF (.NOT. Agrif_Root()) clname = TRIM(Agrif_CFixed()) // "_" // TRIM(clname)
        istatus = NF90_CREATE(TRIM(clname), NF90_CLOBBER, idrun)
        istatus = NF90_DEF_DIM(idrun, 'time', NF90_UNLIMITED, idtime)
        istatus = NF90_DEF_VAR(idrun, 'abs_ssh_max', NF90_DOUBLE, (/idtime/), idssh)
        istatus = NF90_DEF_VAR(idrun, 'abs_u_max', NF90_DOUBLE, (/idtime/), idu)
        istatus = NF90_DEF_VAR(idrun, 's_min', NF90_DOUBLE, (/idtime/), ids1)
        istatus = NF90_DEF_VAR(idrun, 's_max', NF90_DOUBLE, (/idtime/), ids2)
        istatus = NF90_DEF_VAR(idrun, 't_min', NF90_DOUBLE, (/idtime/), idt1)
        istatus = NF90_DEF_VAR(idrun, 't_max', NF90_DOUBLE, (/idtime/), idt2)
        IF (ln_zad_Aimp) THEN
          istatus = NF90_DEF_VAR(idrun, 'abs_wi_max', NF90_DOUBLE, (/idtime/), idw1)
          istatus = NF90_DEF_VAR(idrun, 'Cu_max', NF90_DOUBLE, (/idtime/), idc1)
        END IF
        CALL profile_psy_data2 % PostEnd
        istatus = NF90_ENDDEF(idrun)
        !$ACC KERNELS
        zmax(8 : 9) = 0._wp
        !$ACC END KERNELS
      END IF
    END IF
    CALL profile_psy_data3 % PreStart('stp_ctl', 'r3', 0, 0)
    IF (kt == nit000) lsomeoce = COUNT(ssmask(:, :) == 1._wp) > 0
    IF (lwm .AND. ll_wrtstp) THEN
      WRITE(numstp, '(1x, i8)') kt
      REWIND(UNIT = numstp)
    END IF
    CALL profile_psy_data3 % PostEnd ! CDe added
    IF (ll_wd) THEN
      !$ACC KERNELS ! CDe added      
      zmax(1) = MAXVAL(ABS(sshn(:, :) + ssh_ref * tmask(:, :, 1)))
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS ! CDe added      
      zmax(1) = MAXVAL(ABS(sshn(:, :)))
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS ! CDe added
    zmax(2) = MAXVAL(ABS(un(:, :, :)))
    zmax(3) = MAXVAL(- tsn(:, :, :, jp_sal), mask = tmask(:, :, :) == 1._wp)
    zmax(4) = MAXVAL(tsn(:, :, :, jp_sal), mask = tmask(:, :, :) == 1._wp)
    zmax(5) = MAXVAL(- tsn(:, :, :, jp_tem), mask = tmask(:, :, :) == 1._wp)
    zmax(6) = MAXVAL(tsn(:, :, :, jp_tem), mask = tmask(:, :, :) == 1._wp)
    zmax(7) = REAL(nstop, wp)
    !$ACC END KERNELS
    IF (ln_zad_Aimp) THEN
      !$ACC KERNELS ! CDe added      
      zmax(8) = MAXVAL(ABS(wi(:, :, :)), mask = wmask(:, :, :) == 1._wp)
      zmax(9) = MAXVAL(Cu_adv(:, :, :), mask = tmask(:, :, :) == 1._wp)
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data4 % PreStart('stp_ctl', 'r4', 0, 0)
    IF (ll_colruns) THEN
      CALL mpp_max("stpctl", zmax)
      nstop = NINT(zmax(7))
    END IF
    IF (ll_wrtruns) THEN
      WRITE(numrun, 9500) kt, zmax(1), zmax(2), - zmax(3), zmax(4)
      istatus = NF90_PUT_VAR(idrun, idssh, (/zmax(1)/), (/kt/), (/1/))
      istatus = NF90_PUT_VAR(idrun, idu, (/zmax(2)/), (/kt/), (/1/))
      istatus = NF90_PUT_VAR(idrun, ids1, (/- zmax(3)/), (/kt/), (/1/))
      istatus = NF90_PUT_VAR(idrun, ids2, (/zmax(4)/), (/kt/), (/1/))
      istatus = NF90_PUT_VAR(idrun, idt1, (/- zmax(5)/), (/kt/), (/1/))
      istatus = NF90_PUT_VAR(idrun, idt2, (/zmax(6)/), (/kt/), (/1/))
      IF (ln_zad_Aimp) THEN
        istatus = NF90_PUT_VAR(idrun, idw1, (/zmax(8)/), (/kt/), (/1/))
        istatus = NF90_PUT_VAR(idrun, idc1, (/zmax(9)/), (/kt/), (/1/))
      END IF
      IF (MOD(kt, 100) == 0) istatus = NF90_SYNC(idrun)
      IF (kt == nitend) istatus = NF90_CLOSE(idrun)
    END IF
    IF ((ln_ctl .OR. lsomeoce) .AND. (zmax(1) > 20._wp .OR. zmax(2) > 10._wp .OR. zmax(3) >= 0._wp .OR. zmax(4) >= 100._wp .OR. zmax(4) < 0._wp .OR. IEEE_IS_NAN(zmax(1) + zmax(2) + zmax(3)))) THEN
      IF (lk_mpp .AND. ln_ctl) THEN
        CALL mpp_maxloc('stpctl', ABS(sshn), ssmask(:, :), zzz, ih)
        CALL mpp_maxloc('stpctl', ABS(un), umask(:, :, :), zzz, iu)
        CALL mpp_minloc('stpctl', tsn(:, :, :, jp_sal), tmask(:, :, :), zzz, is1)
        CALL mpp_maxloc('stpctl', tsn(:, :, :, jp_sal), tmask(:, :, :), zzz, is2)
      ELSE
        ih(:) = MAXLOC(ABS(sshn(:, :))) + (/nimpp - 1, njmpp - 1/)
        iu(:) = MAXLOC(ABS(un(:, :, :))) + (/nimpp - 1, njmpp - 1, 0/)
        is1(:) = MINLOC(tsn(:, :, :, jp_sal), mask = tmask(:, :, :) == 1._wp) + (/nimpp - 1, njmpp - 1, 0/)
        is2(:) = MAXLOC(tsn(:, :, :, jp_sal), mask = tmask(:, :, :) == 1._wp) + (/nimpp - 1, njmpp - 1, 0/)
      END IF
      WRITE(ctmp1, *) ' stp_ctl: |ssh| > 20 m  or  |U| > 10 m/s  or  S <= 0  or  S >= 100  or  NaN encounter in the tests'
      WRITE(ctmp2, 9100) kt, zmax(1), ih(1), ih(2)
      WRITE(ctmp3, 9200) kt, zmax(2), iu(1), iu(2), iu(3)
      WRITE(ctmp4, 9300) kt, - zmax(3), is1(1), is1(2), is1(3)
      WRITE(ctmp5, 9400) kt, zmax(4), is2(1), is2(2), is2(3)
      WRITE(ctmp6, *) '      ===> output of last computed fields in output.abort.nc file'
      CALL dia_wri_state('output.abort')
      IF (.NOT. ln_ctl) THEN
        WRITE(ctmp8, *) 'E R R O R message from sub-domain: ', narea
        CALL ctl_stop('STOP', ctmp1, ' ', ctmp8, ' ', ctmp2, ctmp3, ctmp4, ctmp5, ctmp6)
      ELSE
        CALL ctl_stop(ctmp1, ' ', ctmp2, ctmp3, ctmp4, ctmp5, ' ', ctmp6, ' ')
      END IF
      kindic = - 3
    END IF
9100 FORMAT(' kt=', I8, '   |ssh| max: ', 1P, G11.4, ', at  i j  : ', 2I5)
9200 FORMAT(' kt=', I8, '   |U|   max: ', 1P, G11.4, ', at  i j k: ', 3I5)
9300 FORMAT(' kt=', I8, '   S     min: ', 1P, G11.4, ', at  i j k: ', 3I5)
9400 FORMAT(' kt=', I8, '   S     max: ', 1P, G11.4, ', at  i j k: ', 3I5)
9500 FORMAT(' it :', I8, '    |ssh|_max: ', D23.16, ' |U|_max: ', D23.16, ' S_min: ', D23.16, ' S_max: ', D23.16)
    CALL profile_psy_data3 % PostEnd
  END SUBROUTINE stp_ctl
END MODULE stpctl
