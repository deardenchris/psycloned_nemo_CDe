MODULE dia25h
  USE oce
  USE dom_oce
  USE zdf_oce
  USE zdfgls, ONLY: hmxl_n
  USE in_out_manager
  USE iom
  USE wet_dry
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dia_25h_init
  PUBLIC :: dia_25h
  LOGICAL, PUBLIC :: ln_dia25h
  INTEGER, SAVE :: cnt_25h
  REAL(KIND = wp), SAVE :: r1_25 = 0.04_wp
  REAL(KIND = wp), SAVE, ALLOCATABLE, DIMENSION(:, :, :) :: tn_25h, sn_25h
  REAL(KIND = wp), SAVE, ALLOCATABLE, DIMENSION(:, :) :: sshn_25h
  REAL(KIND = wp), SAVE, ALLOCATABLE, DIMENSION(:, :, :) :: un_25h, vn_25h, wn_25h
  REAL(KIND = wp), SAVE, ALLOCATABLE, DIMENSION(:, :, :) :: avt_25h, avm_25h
  REAL(KIND = wp), SAVE, ALLOCATABLE, DIMENSION(:, :, :) :: en_25h, rmxln_25h
  CONTAINS
  SUBROUTINE dia_25h_init
    INTEGER :: ios
    INTEGER :: ierror
    NAMELIST /nam_dia25h/ ln_dia25h
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, nam_dia25h, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'nam_dia25h in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, nam_dia25h, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'nam_dia25h in configuration namelist', lwp)
    IF (lwm) WRITE(numond, nam_dia25h)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dia_25h_init : Output 25 hour mean diagnostics'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist nam_dia25h : set 25h outputs '
      WRITE(numout, FMT = *) '      Switch for 25h diagnostics (T) or not (F)  ln_dia25h  = ', ln_dia25h
    END IF
    IF (.NOT. ln_dia25h) RETURN
    ALLOCATE(tn_25h(jpi, jpj, jpk), sn_25h(jpi, jpj, jpk), sshn_25h(jpi, jpj), un_25h(jpi, jpj, jpk), vn_25h(jpi, jpj, jpk), &
&wn_25h(jpi, jpj, jpk), avt_25h(jpi, jpj, jpk), avm_25h(jpi, jpj, jpk), STAT = ierror)
    IF (ierror > 0) THEN
      CALL ctl_stop('dia_25h: unable to allocate ocean arrays')
      RETURN
    END IF
    IF (ln_zdftke) THEN
      ALLOCATE(en_25h(jpi, jpj, jpk), STAT = ierror)
      IF (ierror > 0) THEN
        CALL ctl_stop('dia_25h: unable to allocate en_25h')
        RETURN
      END IF
    END IF
    IF (ln_zdfgls) THEN
      ALLOCATE(en_25h(jpi, jpj, jpk), rmxln_25h(jpi, jpj, jpk), STAT = ierror)
      IF (ierror > 0) THEN
        CALL ctl_stop('dia_25h: unable to allocate en_25h and rmxln_25h')
        RETURN
      END IF
    END IF
    cnt_25h = 1
    tn_25h(:, :, :) = tsb(:, :, :, jp_tem)
    sn_25h(:, :, :) = tsb(:, :, :, jp_sal)
    sshn_25h(:, :) = sshb(:, :)
    un_25h(:, :, :) = ub(:, :, :)
    vn_25h(:, :, :) = vb(:, :, :)
    avt_25h(:, :, :) = avt(:, :, :)
    avm_25h(:, :, :) = avm(:, :, :)
    IF (ln_zdftke) THEN
      en_25h(:, :, :) = en(:, :, :)
    END IF
    IF (ln_zdfgls) THEN
      en_25h(:, :, :) = en(:, :, :)
      rmxln_25h(:, :, :) = hmxl_n(:, :, :)
    END IF
    CALL ctl_stop('STOP', 'dia_25h not setup yet to do tidemean ice')
  END SUBROUTINE dia_25h_init
  SUBROUTINE dia_25h(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: iyear0, nimonth0, iday0
    LOGICAL :: ll_print = .FALSE.
    REAL(KIND = wp) :: zsto, zout, zmax, zjulian, zmdi
    INTEGER :: i_steps
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zw2d, un_dm, vn_dm
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zw3d
    REAL(KIND = wp), DIMENSION(jpi, jpj, 3) :: zwtmb
    IF (MOD(3600, NINT(rdt)) == 0) THEN
      i_steps = 3600 / NINT(rdt)
    ELSE
      CALL ctl_stop('STOP', 'dia_wri_tide: timestep must give MOD(3600,rdt) = 0 otherwise no hourly values are possible')
    END IF
    ll_print = ll_print .AND. lwp
    IF (kt == nn_it000) THEN
      wn_25h(:, :, :) = wn(:, :, :)
    END IF
    IF (MOD(kt, i_steps) == 0 .AND. kt /= nn_it000) THEN
      IF (lwp) THEN
        WRITE(numout, FMT = *) 'dia_wri_tide : Summing instantaneous hourly diagnostics at timestep ', kt
        WRITE(numout, FMT = *) '~~~~~~~~~~~~ '
      END IF
      tn_25h(:, :, :) = tn_25h(:, :, :) + tsn(:, :, :, jp_tem)
      sn_25h(:, :, :) = sn_25h(:, :, :) + tsn(:, :, :, jp_sal)
      sshn_25h(:, :) = sshn_25h(:, :) + sshn(:, :)
      un_25h(:, :, :) = un_25h(:, :, :) + un(:, :, :)
      vn_25h(:, :, :) = vn_25h(:, :, :) + vn(:, :, :)
      wn_25h(:, :, :) = wn_25h(:, :, :) + wn(:, :, :)
      avt_25h(:, :, :) = avt_25h(:, :, :) + avt(:, :, :)
      avm_25h(:, :, :) = avm_25h(:, :, :) + avm(:, :, :)
      IF (ln_zdftke) THEN
        en_25h(:, :, :) = en_25h(:, :, :) + en(:, :, :)
      END IF
      IF (ln_zdfgls) THEN
        en_25h(:, :, :) = en_25h(:, :, :) + en(:, :, :)
        rmxln_25h(:, :, :) = rmxln_25h(:, :, :) + hmxl_n(:, :, :)
      END IF
      cnt_25h = cnt_25h + 1
      IF (lwp) THEN
        WRITE(numout, FMT = *) 'dia_tide : Summed the following number of hourly values so far', cnt_25h
      END IF
    END IF
    IF (cnt_25h == 25 .AND. MOD(kt, i_steps * 24) == 0 .AND. kt /= nn_it000) THEN
      IF (lwp) THEN
        WRITE(numout, FMT = *) 'dia_wri_tide : Writing 25 hour mean tide diagnostics at timestep', kt
        WRITE(numout, FMT = *) '~~~~~~~~~~~~ '
      END IF
      tn_25h(:, :, :) = tn_25h(:, :, :) * r1_25
      sn_25h(:, :, :) = sn_25h(:, :, :) * r1_25
      sshn_25h(:, :) = sshn_25h(:, :) * r1_25
      un_25h(:, :, :) = un_25h(:, :, :) * r1_25
      vn_25h(:, :, :) = vn_25h(:, :, :) * r1_25
      wn_25h(:, :, :) = wn_25h(:, :, :) * r1_25
      avt_25h(:, :, :) = avt_25h(:, :, :) * r1_25
      avm_25h(:, :, :) = avm_25h(:, :, :) * r1_25
      IF (ln_zdftke) THEN
        en_25h(:, :, :) = en_25h(:, :, :) * r1_25
      END IF
      IF (ln_zdfgls) THEN
        en_25h(:, :, :) = en_25h(:, :, :) * r1_25
        rmxln_25h(:, :, :) = rmxln_25h(:, :, :) * r1_25
      END IF
      IF (lwp) WRITE(numout, FMT = *) 'dia_wri_tide : Mean calculated by dividing 25 hour sums and writing output'
      zmdi = 1.E+20
      zw3d(:, :, :) = tn_25h(:, :, :) * tmask(:, :, :) + zmdi * (1.0 - tmask(:, :, :))
      CALL iom_put("temper25h", zw3d)
      zw3d(:, :, :) = sn_25h(:, :, :) * tmask(:, :, :) + zmdi * (1.0 - tmask(:, :, :))
      CALL iom_put("salin25h", zw3d)
      zw2d(:, :) = sshn_25h(:, :) * tmask(:, :, 1) + zmdi * (1.0 - tmask(:, :, 1))
      IF (ll_wd) THEN
        CALL iom_put("ssh25h", zw2d + ssh_ref)
      ELSE
        CALL iom_put("ssh25h", zw2d)
      END IF
      zw3d(:, :, :) = un_25h(:, :, :) * umask(:, :, :) + zmdi * (1.0 - umask(:, :, :))
      CALL iom_put("vozocrtx25h", zw3d)
      zw3d(:, :, :) = vn_25h(:, :, :) * vmask(:, :, :) + zmdi * (1.0 - vmask(:, :, :))
      CALL iom_put("vomecrty25h", zw3d)
      zw3d(:, :, :) = wn_25h(:, :, :) * wmask(:, :, :) + zmdi * (1.0 - tmask(:, :, :))
      CALL iom_put("vovecrtz25h", zw3d)
      zw3d(:, :, :) = avt_25h(:, :, :) * wmask(:, :, :) + zmdi * (1.0 - tmask(:, :, :))
      CALL iom_put("avt25h", zw3d)
      zw3d(:, :, :) = avm_25h(:, :, :) * wmask(:, :, :) + zmdi * (1.0 - tmask(:, :, :))
      CALL iom_put("avm25h", zw3d)
      IF (ln_zdftke) THEN
        zw3d(:, :, :) = en_25h(:, :, :) * wmask(:, :, :) + zmdi * (1.0 - tmask(:, :, :))
        CALL iom_put("tke25h", zw3d)
      END IF
      IF (ln_zdfgls) THEN
        zw3d(:, :, :) = en_25h(:, :, :) * wmask(:, :, :) + zmdi * (1.0 - tmask(:, :, :))
        CALL iom_put("tke25h", zw3d)
        zw3d(:, :, :) = rmxln_25h(:, :, :) * wmask(:, :, :) + zmdi * (1.0 - tmask(:, :, :))
        CALL iom_put("mxln25h", zw3d)
      END IF
      tn_25h(:, :, :) = tsn(:, :, :, jp_tem)
      sn_25h(:, :, :) = tsn(:, :, :, jp_sal)
      sshn_25h(:, :) = sshn(:, :)
      un_25h(:, :, :) = un(:, :, :)
      vn_25h(:, :, :) = vn(:, :, :)
      wn_25h(:, :, :) = wn(:, :, :)
      avt_25h(:, :, :) = avt(:, :, :)
      avm_25h(:, :, :) = avm(:, :, :)
      IF (ln_zdftke) THEN
        en_25h(:, :, :) = en(:, :, :)
      END IF
      IF (ln_zdfgls) THEN
        en_25h(:, :, :) = en(:, :, :)
        rmxln_25h(:, :, :) = hmxl_n(:, :, :)
      END IF
      cnt_25h = 1
      IF (lwp) WRITE(numout, FMT = *) 'dia_wri_tide :       After 25hr mean write, reset sum to current value and cnt_25h to one &
&for overlapping average', cnt_25h
    END IF
  END SUBROUTINE dia_25h
END MODULE dia25h