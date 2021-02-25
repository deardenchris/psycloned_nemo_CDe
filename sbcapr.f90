MODULE sbcapr
  USE dom_oce
  USE sbc_oce
  USE phycst
  USE fldread
  USE in_out_manager
  USE lib_fortran
  USE iom
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sbc_apr
  PUBLIC :: sbc_apr_init
  LOGICAL, PUBLIC :: ln_apr_obc
  LOGICAL, PUBLIC :: ln_ref_apr
  REAL(KIND = wp) :: rn_pref
  REAL(KIND = wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:, :) :: ssh_ib
  REAL(KIND = wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:, :) :: ssh_ibb
  REAL(KIND = wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:, :) :: apr
  REAL(KIND = wp) :: tarea
  REAL(KIND = wp) :: r1_grau
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_apr
  CONTAINS
  SUBROUTINE sbc_apr_init
    INTEGER :: ierror
    INTEGER :: ios
    CHARACTER(LEN = 100) :: cn_dir
    TYPE(FLD_N) :: sn_apr
    LOGICAL :: lrxios
    NAMELIST /namsbc_apr/ cn_dir, sn_apr, ln_ref_apr, rn_pref, ln_apr_obc
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namsbc_apr, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namsbc_apr in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namsbc_apr, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namsbc_apr in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namsbc_apr)
    ALLOCATE(sf_apr(1), STAT = ierror)
    IF (ierror > 0) CALL ctl_stop('STOP', 'sbc_apr: unable to allocate sf_apr structure')
    CALL fld_fill(sf_apr, (/sn_apr/), cn_dir, 'sbc_apr', 'Atmospheric pressure ', 'namsbc_apr')
    ALLOCATE(sf_apr(1) % fnow(jpi, jpj, 1))
    IF (sn_apr % ln_tint) ALLOCATE(sf_apr(1) % fdta(jpi, jpj, 1, 2))
    ALLOCATE(ssh_ib(jpi, jpj), ssh_ibb(jpi, jpj))
    ALLOCATE(apr(jpi, jpj))
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) '   Namelist namsbc_apr : Atmospheric PRessure as extrenal forcing'
      WRITE(numout, FMT = *) '      ref. pressure: global mean Patm (T) or a constant (F)  ln_ref_apr = ', ln_ref_apr
    END IF
    IF (ln_ref_apr) THEN
      tarea = glob_sum('sbcapr', e1e2t(:, :))
      IF (lwp) WRITE(numout, FMT = *) '         Variable ref. Patm computed over a ocean surface of ', tarea * 1E-6, 'km2'
    ELSE
      IF (lwp) WRITE(numout, FMT = *) '         Reference Patm used : ', rn_pref, ' N/m2'
    END IF
    r1_grau = 1.E0 / (grav * rau0)
    IF (ln_apr_obc) THEN
      IF (lwp) WRITE(numout, FMT = *) '         Inverse barometer added to OBC ssh data'
    END IF
    IF (ln_apr_obc .AND. .NOT. ln_apr_dyn) CALL ctl_warn('sbc_apr: use inverse barometer ssh at open boundary ONLY requires &
&ln_apr_dyn=T')
    IF (lwxios) THEN
      CALL iom_set_rstw_var_active('ssh_ibb')
    END IF
  END SUBROUTINE sbc_apr_init
  SUBROUTINE sbc_apr(kt)
    INTEGER, INTENT(IN) :: kt
    IF (MOD(kt - 1, nn_fsbc) == 0) THEN
      IF (kt /= nit000) ssh_ibb(:, :) = ssh_ib(:, :)
      CALL fld_read(kt, nn_fsbc, sf_apr)
      IF (ln_ref_apr) rn_pref = glob_sum('sbcapr', sf_apr(1) % fnow(:, :, 1) * e1e2t(:, :)) / tarea
      ssh_ib(:, :) = - (sf_apr(1) % fnow(:, :, 1) - rn_pref) * r1_grau
      apr(:, :) = sf_apr(1) % fnow(:, :, 1)
      CALL iom_put("ssh_ib", ssh_ib)
    END IF
    IF (kt == nit000) THEN
      IF (ln_rstart .AND. iom_varid(numror, 'ssh_ibb', ldstop = .FALSE.) > 0) THEN
        IF (lwp) WRITE(numout, FMT = *) 'sbc_apr:   ssh_ibb read in the restart file'
        CALL iom_get(numror, jpdom_autoglo, 'ssh_ibb', ssh_ibb, ldxios = lrxios)
      ELSE
        IF (lwp) WRITE(numout, FMT = *) 'sbc_apr:   ssh_ibb set to nit000 values'
        ssh_ibb(:, :) = ssh_ib(:, :)
      END IF
    END IF
    IF (lrst_oce) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'sbc_apr : ssh_ib written in ocean restart file at it= ', kt, ' date= ', ndastp
      IF (lwp) WRITE(numout, FMT = *) '~~~~'
      IF (lwxios) CALL iom_swap(cwxios_context)
      CALL iom_rstput(kt, nitrst, numrow, 'ssh_ibb', ssh_ib, ldxios = lwxios)
      IF (lwxios) CALL iom_swap(cxios_context)
    END IF
  END SUBROUTINE sbc_apr
END MODULE sbcapr