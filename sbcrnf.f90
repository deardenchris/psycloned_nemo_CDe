MODULE sbcrnf
  USE dom_oce
  USE phycst
  USE sbc_oce
  USE sbcisf
  USE eosbn2
  USE closea
  USE in_out_manager
  USE fldread
  USE iom
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sbc_rnf
  PUBLIC :: sbc_rnf_div
  PUBLIC :: sbc_rnf_alloc
  PUBLIC :: sbc_rnf_init
  CHARACTER(LEN = 100) :: cn_dir
  LOGICAL, PUBLIC :: ln_rnf_depth
  LOGICAL :: ln_rnf_depth_ini
  REAL(KIND = wp) :: rn_rnf_max
  REAL(KIND = wp) :: rn_dep_max
  INTEGER :: nn_rnf_depth_file
  LOGICAL :: ln_rnf_tem
  LOGICAL, PUBLIC :: ln_rnf_sal
  TYPE(FLD_N), PUBLIC :: sn_rnf
  TYPE(FLD_N) :: sn_cnf
  TYPE(FLD_N) :: sn_s_rnf
  TYPE(FLD_N) :: sn_t_rnf
  TYPE(FLD_N) :: sn_dep_rnf
  LOGICAL, PUBLIC :: ln_rnf_mouth
  REAL(KIND = wp) :: rn_hrnf
  REAL(KIND = wp), PUBLIC :: rn_avt_rnf
  REAL(KIND = wp), PUBLIC :: rn_rfact
  LOGICAL, PUBLIC :: l_rnfcpl = .FALSE.
  INTEGER, PUBLIC :: nkrnf = 0
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: rnfmsk
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: rnfmsk_z
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: h_rnf
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: nk_rnf
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: rnf_tsc_b, rnf_tsc
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_rnf
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_s_rnf
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_t_rnf
  CONTAINS
  INTEGER FUNCTION sbc_rnf_alloc()
    ALLOCATE(rnfmsk(jpi, jpj), rnfmsk_z(jpk), h_rnf(jpi, jpj), nk_rnf(jpi, jpj), rnf_tsc_b(jpi, jpj, jpts), rnf_tsc(jpi, jpj, jpts), STAT = sbc_rnf_alloc)
    CALL mpp_sum('sbcrnf', sbc_rnf_alloc)
    IF (sbc_rnf_alloc > 0) CALL ctl_warn('sbc_rnf_alloc: allocation of arrays failed')
  END FUNCTION sbc_rnf_alloc
  SUBROUTINE sbc_rnf(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj
    INTEGER :: z_err = 0
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztfrz
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    CALL profile_psy_data0 % PreStart('sbc_rnf', 'r0', 0, 0)
    IF (.NOT. l_rnfcpl) CALL fld_read(kt, nn_fsbc, sf_rnf)
    IF (ln_rnf_tem) CALL fld_read(kt, nn_fsbc, sf_t_rnf)
    IF (ln_rnf_sal) CALL fld_read(kt, nn_fsbc, sf_s_rnf)
    CALL profile_psy_data0 % PostEnd
    IF (MOD(kt - 1, nn_fsbc) == 0) THEN
      CALL profile_psy_data1 % PreStart('sbc_rnf', 'r1', 0, 0)
      !$ACC KERNELS ! CDe
      IF (.NOT. l_rnfcpl) rnf(:, :) = rn_rfact * (sf_rnf(1) % fnow(:, :, 1)) * tmask(:, :, 1)
      !$ACC END KERNELS
      CALL profile_psy_data1 % PostEnd
      IF (ln_rnf_tem) THEN
        !$ACC KERNELS
        rnf_tsc(:, :, jp_tem) = (sf_t_rnf(1) % fnow(:, :, 1)) * rnf(:, :) * r1_rau0
        !$ACC END KERNELS
        CALL profile_psy_data2 % PreStart('sbc_rnf', 'r2', 0, 0)
        CALL eos_fzp(sss_m(:, :), ztfrz(:, :))
        WHERE (sf_t_rnf(1) % fnow(:, :, 1) == - 999._wp)
          rnf_tsc(:, :, jp_tem) = sst_m(:, :) * rnf(:, :) * r1_rau0
        END WHERE
        WHERE (sf_t_rnf(1) % fnow(:, :, 1) == - 222._wp)
          rnf_tsc(:, :, jp_tem) = ztfrz(:, :) * rnf(:, :) * r1_rau0 - rnf(:, :) * rLfusisf * r1_rau0_rcp
        END WHERE
        CALL profile_psy_data2 % PostEnd
      ELSE
        !$ACC KERNELS
        rnf_tsc(:, :, jp_tem) = MAX(sst_m(:, :), 0.0_wp) * rnf(:, :) * r1_rau0
        !$ACC END KERNELS
      END IF
      CALL profile_psy_data3 % PreStart('sbc_rnf', 'r3', 0, 0)
      IF (ln_rnf_sal) rnf_tsc(:, :, jp_sal) = (sf_s_rnf(1) % fnow(:, :, 1)) * rnf(:, :) * r1_rau0
      IF (iom_use('runoffs')) CALL iom_put('runoffs', rnf(:, :))
      IF (iom_use('hflx_rnf_cea')) CALL iom_put('hflx_rnf_cea', rnf_tsc(:, :, jp_tem) * rau0 * rcp)
      CALL profile_psy_data3 % PostEnd
    END IF
    IF (kt == nit000) THEN
      IF (ln_rstart .AND. iom_varid(numror, 'rnf_b', ldstop = .FALSE.) > 0) THEN
        CALL profile_psy_data4 % PreStart('sbc_rnf', 'r4', 0, 0)
        IF (lwp) WRITE(numout, *) '          nit000-1 runoff forcing fields red in the restart file', lrxios
        CALL iom_get(numror, jpdom_autoglo, 'rnf_b', rnf_b, ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'rnf_hc_b', rnf_tsc_b(:, :, jp_tem), ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'rnf_sc_b', rnf_tsc_b(:, :, jp_sal), ldxios = lrxios)
        CALL profile_psy_data4 % PostEnd
      ELSE
        IF (lwp) WRITE(numout, *) '          nit000-1 runoff forcing fields set to nit000'
        !$ACC KERNELS
        rnf_b(:, :) = rnf(:, :)
        rnf_tsc_b(:, :, :) = rnf_tsc(:, :, :)
        !$ACC END KERNELS
      END IF
    END IF
    CALL profile_psy_data5 % PreStart('sbc_rnf', 'r5', 0, 0)
    IF (lrst_oce) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'sbcrnf : runoff forcing fields written in ocean restart file ', 'at it= ', kt, ' date= ', ndastp
      IF (lwp) WRITE(numout, *) '~~~~'
      IF (lwxios) CALL iom_swap(cwxios_context)
      CALL iom_rstput(kt, nitrst, numrow, 'rnf_b', rnf, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'rnf_hc_b', rnf_tsc(:, :, jp_tem), ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'rnf_sc_b', rnf_tsc(:, :, jp_sal), ldxios = lwxios)
      IF (lwxios) CALL iom_swap(cxios_context)
    END IF
    CALL profile_psy_data5 % PostEnd
  END SUBROUTINE sbc_rnf
  SUBROUTINE sbc_rnf_div(phdivn)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: phdivn
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zfact
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('sbc_rnf_div', 'r0', 0, 0)
    !$ACC KERNELS ! CDe
    zfact = 0.5_wp
    !$ACC END KERNELS
    CALL profile_psy_data0 % PostEnd
    IF (ln_rnf_depth .OR. ln_rnf_depth_ini) THEN
      IF (ln_linssh) THEN              
        DO jj = 1, jpj
          DO ji = 1, jpi
            !$ACC KERNELS
            DO jk = 1, nk_rnf(ji, jj)
              phdivn(ji, jj, jk) = phdivn(ji, jj, jk) - (rnf(ji, jj) + rnf_b(ji, jj)) * zfact * r1_rau0 / h_rnf(ji, jj)
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      ELSE
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            h_rnf(ji, jj) = 0._wp
            DO jk = 1, nk_rnf(ji, jj)
              h_rnf(ji, jj) = h_rnf(ji, jj) + e3t_n(ji, jj, jk)
            END DO
            DO jk = 1, nk_rnf(ji, jj)
              phdivn(ji, jj, jk) = phdivn(ji, jj, jk) - (rnf(ji, jj) + rnf_b(ji, jj)) * zfact * r1_rau0 / h_rnf(ji, jj)
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      END IF
    ELSE
      !$ACC KERNELS
      h_rnf(:, :) = e3t_n(:, :, 1)
      phdivn(:, :, 1) = phdivn(:, :, 1) - (rnf(:, :) + rnf_b(:, :)) * zfact * r1_rau0 / e3t_n(:, :, 1)
      !$ACC END KERNELS
    END IF
  END SUBROUTINE sbc_rnf_div
  SUBROUTINE sbc_rnf_init
    CHARACTER(LEN = 32) :: rn_dep_file
    INTEGER :: ji, jj, jk, jm
    INTEGER :: ierror, inum
    INTEGER :: ios
    INTEGER :: nbrec
    REAL(KIND = wp) :: zacoef
    REAL(KIND = wp), DIMENSION(jpi, jpj, 2) :: zrnfcl
    NAMELIST /namsbc_rnf/ cn_dir, ln_rnf_depth, ln_rnf_tem, ln_rnf_sal, sn_rnf, sn_cnf, sn_s_rnf, sn_t_rnf, sn_dep_rnf, ln_rnf_mouth, rn_hrnf, rn_avt_rnf, rn_rfact, ln_rnf_depth_ini, rn_dep_max, rn_rnf_max, nn_rnf_depth_file
    IF (sbc_rnf_alloc() /= 0) CALL ctl_stop('STOP', 'sbc_rnf_alloc : unable to allocate arrays')
    IF (.NOT. ln_rnf) THEN
      !$ACC KERNELS
      ln_rnf_mouth = .FALSE.
      nkrnf = 0
      rnf(:, :) = 0.0_wp
      rnf_b(:, :) = 0.0_wp
      rnfmsk(:, :) = 0.0_wp
      rnfmsk_z(:) = 0.0_wp
      !$ACC END KERNELS
      RETURN
    END IF
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namsbc_rnf, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namsbc_rnf in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namsbc_rnf, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namsbc_rnf in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namsbc_rnf)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'sbc_rnf_init : runoff '
      WRITE(numout, *) '~~~~~~~~~~~~ '
      WRITE(numout, *) '   Namelist namsbc_rnf'
      WRITE(numout, *) '      specific river mouths treatment            ln_rnf_mouth = ', ln_rnf_mouth
      WRITE(numout, *) '      river mouth additional Kz                  rn_avt_rnf   = ', rn_avt_rnf
      WRITE(numout, *) '      depth of river mouth additional mixing     rn_hrnf      = ', rn_hrnf
      WRITE(numout, *) '      multiplicative factor for runoff           rn_rfact     = ', rn_rfact
    END IF
    IF (.NOT. l_rnfcpl) THEN
      ALLOCATE(sf_rnf(1), STAT = ierror)
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   ==>>>   runoffs inflow read in a file'
      IF (ierror > 0) THEN
        CALL ctl_stop('sbc_rnf_init: unable to allocate sf_rnf structure')
        RETURN
      END IF
      ALLOCATE(sf_rnf(1) % fnow(jpi, jpj, 1))
      IF (sn_rnf % ln_tint) ALLOCATE(sf_rnf(1) % fdta(jpi, jpj, 1, 2))
      CALL fld_fill(sf_rnf, (/sn_rnf/), cn_dir, 'sbc_rnf_init', 'read runoffs data', 'namsbc_rnf', no_print)
    END IF
    IF (ln_rnf_tem) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   ==>>>   runoffs temperatures read in a file'
      ALLOCATE(sf_t_rnf(1), STAT = ierror)
      IF (ierror > 0) THEN
        CALL ctl_stop('sbc_rnf_init: unable to allocate sf_t_rnf structure')
        RETURN
      END IF
      ALLOCATE(sf_t_rnf(1) % fnow(jpi, jpj, 1))
      IF (sn_t_rnf % ln_tint) ALLOCATE(sf_t_rnf(1) % fdta(jpi, jpj, 1, 2))
      CALL fld_fill(sf_t_rnf, (/sn_t_rnf/), cn_dir, 'sbc_rnf_init', 'read runoff temperature data', 'namsbc_rnf', no_print)
    END IF
    IF (ln_rnf_sal) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   ==>>>   runoffs salinities read in a file'
      ALLOCATE(sf_s_rnf(1), STAT = ierror)
      IF (ierror > 0) THEN
        CALL ctl_stop('sbc_rnf_init: unable to allocate sf_s_rnf structure')
        RETURN
      END IF
      ALLOCATE(sf_s_rnf(1) % fnow(jpi, jpj, 1))
      IF (sn_s_rnf % ln_tint) ALLOCATE(sf_s_rnf(1) % fdta(jpi, jpj, 1, 2))
      CALL fld_fill(sf_s_rnf, (/sn_s_rnf/), cn_dir, 'sbc_rnf_init', 'read runoff salinity data', 'namsbc_rnf', no_print)
    END IF
    IF (ln_rnf_depth) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   ==>>>   runoffs depth read in a file'
      rn_dep_file = TRIM(cn_dir) // TRIM(sn_dep_rnf % clname)
      IF (.NOT. sn_dep_rnf % ln_clim) THEN
        WRITE(rn_dep_file, '(a,"_y",i4)') TRIM(rn_dep_file), nyear
        IF (sn_dep_rnf % cltype == 'monthly') WRITE(rn_dep_file, '(a,"m",i2)') TRIM(rn_dep_file), nmonth
      END IF
      CALL iom_open(rn_dep_file, inum)
      CALL iom_get(inum, jpdom_data, sn_dep_rnf % clvar, h_rnf)
      CALL iom_close(inum)
      !$ACC KERNELS
      nk_rnf(:, :) = 0
      !$ACC END KERNELS
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (h_rnf(ji, jj) > 0._wp) THEN
            jk = 2
            DO WHILE (jk /= mbkt(ji, jj) .AND. gdept_0(ji, jj, jk) < h_rnf(ji, jj))
              jk = jk + 1
            END DO
            nk_rnf(ji, jj) = jk
          ELSE IF (h_rnf(ji, jj) == - 1._wp) THEN
            nk_rnf(ji, jj) = 1
          ELSE IF (h_rnf(ji, jj) == - 999._wp) THEN
            nk_rnf(ji, jj) = mbkt(ji, jj)
          ELSE
            CALL ctl_stop('sbc_rnf_init: runoff depth not positive, and not -999 or -1, rnf value in file fort.999')
            WRITE(999, *) 'ji, jj, h_rnf(ji,jj) :', ji, jj, h_rnf(ji, jj)
          END IF
        END DO
      END DO
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          h_rnf(ji, jj) = 0._wp
          DO jk = 1, nk_rnf(ji, jj)
            h_rnf(ji, jj) = h_rnf(ji, jj) + e3t_n(ji, jj, jk)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    ELSE IF (ln_rnf_depth_ini) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   ==>>>   depth of runoff computed once from max value of runoff'
      IF (lwp) WRITE(numout, *) '        max value of the runoff climatologie (over global domain) rn_rnf_max = ', rn_rnf_max
      IF (lwp) WRITE(numout, *) '        depth over which runoffs is spread                        rn_dep_max = ', rn_dep_max
      IF (lwp) WRITE(numout, *) '        create (=1) a runoff depth file or not (=0)      nn_rnf_depth_file  = ', nn_rnf_depth_file
      CALL iom_open(TRIM(sn_rnf % clname), inum)
      !$ACC KERNELS
      nbrec = iom_getszuld(inum)
      zrnfcl(:, :, 1) = 0._wp
      !$ACC END KERNELS
      DO jm = 1, nbrec
        CALL iom_get(inum, jpdom_data, TRIM(sn_rnf % clvar), zrnfcl(:, :, 2), jm)
        zrnfcl(:, :, 1) = MAXVAL(zrnfcl(:, :, :), DIM = 3)
      END DO
      CALL iom_close(inum)
      !$ACC KERNELS
      h_rnf(:, :) = 1.
      zacoef = rn_dep_max / rn_rnf_max
      WHERE (zrnfcl(:, :, 1) > 0._wp) h_rnf(:, :) = zacoef * zrnfcl(:, :, 1)
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (zrnfcl(ji, jj, 1) > 0._wp) THEN
            jk = mbkt(ji, jj)
            h_rnf(ji, jj) = MIN(h_rnf(ji, jj), gdept_0(ji, jj, jk))
          END IF
        END DO
      END DO
      nk_rnf(:, :) = 0
      !$ACC END KERNELS
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (zrnfcl(ji, jj, 1) > 0._wp) THEN
            jk = 2
            DO WHILE (jk /= mbkt(ji, jj) .AND. gdept_0(ji, jj, jk) < h_rnf(ji, jj))
              jk = jk + 1
            END DO
            nk_rnf(ji, jj) = jk
          ELSE
            nk_rnf(ji, jj) = 1
          END IF
        END DO
      END DO
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          h_rnf(ji, jj) = 0._wp
          DO jk = 1, nk_rnf(ji, jj)
            h_rnf(ji, jj) = h_rnf(ji, jj) + e3t_n(ji, jj, jk)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      IF (nn_rnf_depth_file == 1) THEN
        IF (lwp) WRITE(numout, *) '   ==>>>   create runoff depht file'
        CALL iom_open(TRIM(sn_dep_rnf % clname), inum, ldwrt = .TRUE.)
        CALL iom_rstput(0, 0, inum, 'rodepth', h_rnf)
        CALL iom_close(inum)
      END IF
    ELSE
      !$ACC KERNELS
      nk_rnf(:, :) = 1
      h_rnf(:, :) = e3t_n(:, :, 1)
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    rnf(:, :) = 0._wp
    rnf_tsc(:, :, :) = 0._wp
    !$ACC END KERNELS
    IF (ln_rnf_mouth) THEN
      IF (ln_rnf_depth) CALL ctl_warn('sbc_rnf_init: increased mixing turned on but effects may already', 'be spread through depth by ln_rnf_depth')
      nkrnf = 0
      IF (rn_hrnf > 0._wp) THEN
        nkrnf = 2
        DO WHILE (nkrnf /= jpkm1 .AND. gdepw_1d(nkrnf + 1) < rn_hrnf)
          nkrnf = nkrnf + 1
        END DO
        IF (ln_sco) CALL ctl_warn('sbc_rnf_init: number of levels over which Kz is increased is computed for zco...')
      END IF
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   ==>>>   Specific treatment used in vicinity of river mouths :'
      IF (lwp) WRITE(numout, *) '             - Increase Kz in surface layers (if rn_hrnf > 0 )'
      IF (lwp) WRITE(numout, *) '               by ', rn_avt_rnf, ' m2/s  over ', nkrnf, ' w-levels'
      IF (lwp) WRITE(numout, *) '             - set to zero SSS damping       (if ln_ssr=T)'
      IF (lwp) WRITE(numout, *) '             - mixed upstream-centered       (if ln_traadv_cen2=T)'
      CALL rnf_mouth
    ELSE
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   ==>>>   No specific treatment at river mouths'
      !$ACC KERNELS
      rnfmsk(:, :) = 0._wp
      rnfmsk_z(:) = 0._wp
      nkrnf = 0
      !$ACC END KERNELS
    END IF
    IF (lwxios) THEN
      CALL iom_set_rstw_var_active('rnf_b')
      CALL iom_set_rstw_var_active('rnf_hc_b')
      CALL iom_set_rstw_var_active('rnf_sc_b')
    END IF
  END SUBROUTINE sbc_rnf_init
  SUBROUTINE rnf_mouth
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: inum
    CHARACTER(LEN = 140) :: cl_rnfile
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('rnf_mouth', 'r0', 0, 0)
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) '   rnf_mouth : river mouth mask'
    IF (lwp) WRITE(numout, *) '   ~~~~~~~~~ '
    cl_rnfile = TRIM(cn_dir) // TRIM(sn_cnf % clname)
    IF (.NOT. sn_cnf % ln_clim) THEN
      WRITE(cl_rnfile, '(a,"_y",i4)') TRIM(cl_rnfile), nyear
      IF (sn_cnf % cltype == 'monthly') WRITE(cl_rnfile, '(a,"m",i2)') TRIM(cl_rnfile), nmonth
    END IF
    CALL iom_open(cl_rnfile, inum)
    CALL iom_get(inum, jpdom_data, sn_cnf % clvar, rnfmsk)
    CALL iom_close(inum)
    IF (l_clo_rnf) CALL clo_rnf(rnfmsk)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    rnfmsk_z(:) = 0._wp
    rnfmsk_z(1) = 1.0
    rnfmsk_z(2) = 1.0
    rnfmsk_z(3) = 0.5
    rnfmsk_z(4) = 0.25
    rnfmsk_z(5) = 0.125
    !$ACC END KERNELS
  END SUBROUTINE rnf_mouth
END MODULE sbcrnf
