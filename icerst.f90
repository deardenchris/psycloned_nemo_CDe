MODULE icerst
  USE ice
  USE dom_oce
  USE sbc_oce, ONLY: nn_fsbc, ln_cpl
  USE icectl
  USE in_out_manager
  USE iom
  USE ioipsl, ONLY: ju2ymds
  USE lib_mpp
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_rst_opn
  PUBLIC :: ice_rst_write
  PUBLIC :: ice_rst_read
  CONTAINS
  SUBROUTINE ice_rst_opn(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: iyear, imonth, iday
    REAL(KIND = wp) :: zsec
    REAL(KIND = wp) :: zfjulday
    CHARACTER(LEN = 20) :: clkt
    CHARACTER(LEN = 50) :: clname
    CHARACTER(LEN = 256) :: clpath
    IF (kt == nit000) lrst_ice = .FALSE.
    IF (kt == nitrst - 2 * nn_fsbc + 1 .OR. nstock == nn_fsbc .OR. (kt == nitend - nn_fsbc + 1 .AND. .NOT. lrst_ice)) THEN
      IF (nitrst <= nitend .AND. nitrst > 0) THEN
        IF (ln_rstdate) THEN
          zfjulday = fjulday + (2 * nn_fsbc + 1) * rdt / rday
          IF (ABS(zfjulday - REAL(NINT(zfjulday), wp)) < 0.1 / rday) zfjulday = REAL(NINT(zfjulday), wp)
          CALL ju2ymds(zfjulday, iyear, imonth, iday, zsec)
          WRITE(clkt, '(i4.4,2i2.2)') iyear, imonth, iday
        ELSE
          IF (nitrst > 99999999) THEN
            WRITE(clkt, *) nitrst
          ELSE
            WRITE(clkt, '(i8.8)') nitrst
          END IF
        END IF
        clname = TRIM(cexper) // "_" // TRIM(ADJUSTL(clkt)) // "_" // TRIM(cn_icerst_out)
        clpath = TRIM(cn_icerst_outdir)
        IF (clpath(LEN_TRIM(clpath) :) /= '/') clpath = TRIM(clpath) // '/'
        IF (lwp) THEN
          WRITE(numout, *)
          WRITE(numout, *) '             open ice restart NetCDF file: ', TRIM(clpath) // clname
          IF (kt == nitrst - 2 * nn_fsbc + 1) THEN
            WRITE(numout, *) '             kt = nitrst - 2*nn_fsbc + 1 = ', kt, ' date= ', ndastp
          ELSE
            WRITE(numout, *) '             kt = ', kt, ' date= ', ndastp
          END IF
        END IF
        CALL iom_open(TRIM(clpath) // TRIM(clname), numriw, ldwrt = .TRUE., kdlev = jpl)
        lrst_ice = .TRUE.
      END IF
    END IF
    IF (ln_icectl) CALL ice_prt(kt, iiceprt, jiceprt, 1, ' - Beginning the time step - ')
  END SUBROUTINE ice_rst_opn
  SUBROUTINE ice_rst_write(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: jk
    INTEGER :: iter
    CHARACTER(LEN = 25) :: znam
    CHARACTER(LEN = 2) :: zchar, zchar1
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: z3d
    iter = kt + nn_fsbc - 1
    IF (iter == nitrst) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'ice_rst_write : write ice restart file  kt =', kt
      IF (lwp) WRITE(numout, *) '~~~~~~~~~~~~~~'
    END IF
    CALL iom_rstput(iter, nitrst, numriw, 'nn_fsbc', REAL(nn_fsbc, wp))
    CALL iom_rstput(iter, nitrst, numriw, 'kt_ice', REAL(iter, wp))
    CALL iom_delay_rst('WRITE', 'ICE', numriw)
    CALL iom_rstput(iter, nitrst, numriw, 'v_i', v_i)
    CALL iom_rstput(iter, nitrst, numriw, 'v_s', v_s)
    CALL iom_rstput(iter, nitrst, numriw, 'sv_i', sv_i)
    CALL iom_rstput(iter, nitrst, numriw, 'oa_i', oa_i)
    CALL iom_rstput(iter, nitrst, numriw, 'a_i', a_i)
    CALL iom_rstput(iter, nitrst, numriw, 't_su', t_su)
    CALL iom_rstput(iter, nitrst, numriw, 'a_ip', a_ip)
    CALL iom_rstput(iter, nitrst, numriw, 'v_ip', v_ip)
    DO jk = 1, nlay_s
      WRITE(zchar1, '(I2.2)') jk
      znam = 'e_s' // '_l' // zchar1
      !$ACC KERNELS
      z3d(:, :, :) = e_s(:, :, jk, :)
      !$ACC END KERNELS
      CALL iom_rstput(iter, nitrst, numriw, znam, z3d)
    END DO
    DO jk = 1, nlay_i
      WRITE(zchar1, '(I2.2)') jk
      znam = 'e_i' // '_l' // zchar1
      !$ACC KERNELS
      z3d(:, :, :) = e_i(:, :, jk, :)
      !$ACC END KERNELS
      CALL iom_rstput(iter, nitrst, numriw, znam, z3d)
    END DO
    CALL iom_rstput(iter, nitrst, numriw, 'u_ice', u_ice)
    CALL iom_rstput(iter, nitrst, numriw, 'v_ice', v_ice)
    IF (ln_cpl) THEN
      CALL iom_rstput(iter, nitrst, numriw, 'cnd_ice', cnd_ice)
      CALL iom_rstput(iter, nitrst, numriw, 't1_ice', t1_ice)
    END IF
    IF (iter == nitrst) THEN
      CALL iom_close(numriw)
      lrst_ice = .FALSE.
    END IF
  END SUBROUTINE ice_rst_write
  SUBROUTINE ice_rst_read
    INTEGER :: jk
    LOGICAL :: llok
    INTEGER :: id1
    CHARACTER(LEN = 25) :: znam
    CHARACTER(LEN = 2) :: zchar, zchar1
    REAL(KIND = wp) :: zfice, ziter
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: z3d
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ice_rst_read: read ice NetCDF restart file'
      WRITE(numout, *) '~~~~~~~~~~~~'
    END IF
    CALL iom_open(TRIM(cn_icerst_indir) // '/' // cn_icerst_in, numrir, kdlev = jpl)
    CALL iom_get(numrir, 'nn_fsbc', zfice)
    CALL iom_get(numrir, 'kt_ice', ziter)
    IF (lwp) WRITE(numout, *) '   read ice restart file at time step    : ', ziter
    IF (lwp) WRITE(numout, *) '   in any case we force it to nit000 - 1 : ', nit000 - 1
    IF ((nit000 - NINT(ziter)) /= 1 .AND. ABS(nrstdt) == 1) CALL ctl_stop('ice_rst_read ===>>>> : problem with nit000 in ice restart', '   verify the file or rerun with the value 0 for the', '   control of time parameter  nrstdt')
    IF (NINT(zfice) /= nn_fsbc .AND. ABS(nrstdt) == 1) CALL ctl_stop('ice_rst_read ===>>>> : problem with nn_fsbc in ice restart', '   verify the file or rerun with the value 0 for the', '   control of time parameter  nrstdt')
    CALL iom_get(numrir, jpdom_autoglo, 'v_i', v_i)
    CALL iom_get(numrir, jpdom_autoglo, 'v_s', v_s)
    CALL iom_get(numrir, jpdom_autoglo, 'sv_i', sv_i)
    CALL iom_get(numrir, jpdom_autoglo, 'oa_i', oa_i)
    CALL iom_get(numrir, jpdom_autoglo, 'a_i', a_i)
    CALL iom_get(numrir, jpdom_autoglo, 't_su', t_su)
    id1 = iom_varid(numrir, 'a_ip', ldstop = .FALSE.)
    IF (id1 > 0) THEN
      CALL iom_get(numrir, jpdom_autoglo, 'a_ip', a_ip)
      CALL iom_get(numrir, jpdom_autoglo, 'v_ip', v_ip)
    ELSE
      IF (lwp) WRITE(numout, *) '   ==>>   previous run without melt ponds output then set it to zero'
      !$ACC KERNELS
      a_ip(:, :, :) = 0._wp
      v_ip(:, :, :) = 0._wp
      !$ACC END KERNELS
    END IF
    DO jk = 1, nlay_s
      WRITE(zchar1, '(I2.2)') jk
      znam = 'e_s' // '_l' // zchar1
      CALL iom_get(numrir, jpdom_autoglo, znam, z3d)
      !$ACC KERNELS
      e_s(:, :, jk, :) = z3d(:, :, :)
      !$ACC END KERNELS
    END DO
    DO jk = 1, nlay_i
      WRITE(zchar1, '(I2.2)') jk
      znam = 'e_i' // '_l' // zchar1
      CALL iom_get(numrir, jpdom_autoglo, znam, z3d)
      !$ACC KERNELS
      e_i(:, :, jk, :) = z3d(:, :, :)
      !$ACC END KERNELS
    END DO
    CALL iom_get(numrir, jpdom_autoglo, 'u_ice', u_ice)
    CALL iom_get(numrir, jpdom_autoglo, 'v_ice', v_ice)
    CALL iom_delay_rst('READ', 'ICE', numrir)
    IF (ln_cpl) THEN
      CALL iom_get(numrir, jpdom_autoglo, 'cnd_ice', cnd_ice)
      CALL iom_get(numrir, jpdom_autoglo, 't1_ice', t1_ice)
    END IF
  END SUBROUTINE ice_rst_read
END MODULE icerst