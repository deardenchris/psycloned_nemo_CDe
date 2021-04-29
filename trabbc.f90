MODULE trabbc
  USE oce
  USE dom_oce
  USE phycst
  USE trd_oce
  USE trdtra
  USE in_out_manager
  USE iom
  USE fldread
  USE lbclnk
  USE lib_mpp
  USE prtctl
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_bbc
  PUBLIC :: tra_bbc_init
  LOGICAL, PUBLIC :: ln_trabbc
  INTEGER :: nn_geoflx
  REAL(KIND = wp) :: rn_geoflx_cst
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: qgh_trd0
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_qgh
  CONTAINS
  SUBROUTINE tra_bbc(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ztrdt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('tra_bbc', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('tra_bbc')
    CALL profile_psy_data0 % PostEnd
    IF (l_trdtra) THEN
      ALLOCATE(ztrdt(jpi, jpj, jpk))
      !$ACC KERNELS
      ztrdt(:, :, :) = tsa(:, :, :, jp_tem)
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        tsa(ji, jj, mbkt(ji, jj), jp_tem) = tsa(ji, jj, mbkt(ji, jj), jp_tem) + qgh_trd0(ji, jj) / e3t_n(ji, jj, mbkt(ji, jj))
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('tra_bbc', 'r1', 0, 0)
    CALL lbc_lnk('trabbc', tsa(:, :, :, jp_tem), 'T', 1.)
    CALL profile_psy_data1 % PostEnd
    IF (l_trdtra) THEN
      !$ACC KERNELS
      ztrdt(:, :, :) = tsa(:, :, :, jp_tem) - ztrdt(:, :, :)
      !$ACC END KERNELS
      CALL profile_psy_data2 % PreStart('tra_bbc', 'r2', 0, 0)
      CALL trd_tra(kt, 'TRA', jp_tem, jptra_bbc, ztrdt)
      DEALLOCATE(ztrdt)
      CALL profile_psy_data2 % PostEnd
    END IF
    CALL profile_psy_data3 % PreStart('tra_bbc', 'r3', 0, 0)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = tsa(:, :, :, jp_tem), clinfo1 = ' bbc  - Ta: ', mask1 = tmask, clinfo3 = 'tra-ta')
    IF (ln_timing) CALL timing_stop('tra_bbc')
    CALL profile_psy_data3 % PostEnd
  END SUBROUTINE tra_bbc
  SUBROUTINE tra_bbc_init
    INTEGER :: ji, jj
    INTEGER :: inum
    INTEGER :: ios
    INTEGER :: ierror
    TYPE(FLD_N) :: sn_qgh
    CHARACTER(LEN = 256) :: cn_dir
    NAMELIST /nambbc/ ln_trabbc, nn_geoflx, rn_geoflx_cst, sn_qgh, cn_dir
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, nambbc, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'nambbc in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, nambbc, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'nambbc in configuration namelist', lwp)
    IF (lwm) WRITE(numond, nambbc)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'tra_bbc : Bottom Boundary Condition (bbc), apply a Geothermal heating'
      WRITE(numout, *) '~~~~~~~   '
      WRITE(numout, *) '   Namelist nambbc : set bbc parameters'
      WRITE(numout, *) '      Apply a geothermal heating at ocean bottom   ln_trabbc     = ', ln_trabbc
      WRITE(numout, *) '      type of geothermal flux                      nn_geoflx     = ', nn_geoflx
      WRITE(numout, *) '      Constant geothermal flux value               rn_geoflx_cst = ', rn_geoflx_cst
      WRITE(numout, *)
    END IF
    IF (ln_trabbc) THEN
      ALLOCATE(qgh_trd0(jpi, jpj))
      SELECT CASE (nn_geoflx)
      CASE (1)
        IF (lwp) WRITE(numout, *) '   ==>>>   constant heat flux  =   ', rn_geoflx_cst
        !$ACC KERNELS
        qgh_trd0(:, :) = r1_rau0_rcp * rn_geoflx_cst
        !$ACC END KERNELS
      CASE (2)
        IF (lwp) WRITE(numout, *) '   ==>>>   variable geothermal heat flux'
        ALLOCATE(sf_qgh(1), STAT = ierror)
        IF (ierror > 0) THEN
          CALL ctl_stop('tra_bbc_init: unable to allocate sf_qgh structure')
          RETURN
        END IF
        ALLOCATE(sf_qgh(1) % fnow(jpi, jpj, 1))
        IF (sn_qgh % ln_tint) ALLOCATE(sf_qgh(1) % fdta(jpi, jpj, 1, 2))
        CALL fld_fill(sf_qgh, (/sn_qgh/), cn_dir, 'tra_bbc_init', 'bottom temperature boundary condition', 'nambbc', no_print)
        CALL fld_read(nit000, 1, sf_qgh)
        !$ACC KERNELS
        qgh_trd0(:, :) = r1_rau0_rcp * sf_qgh(1) % fnow(:, :, 1) * 1.E-3
        !$ACC END KERNELS
      CASE DEFAULT
        WRITE(ctmp1, *) '     bad flag value for nn_geoflx = ', nn_geoflx
        CALL ctl_stop(ctmp1)
      END SELECT
    ELSE
      IF (lwp) WRITE(numout, *) '   ==>>>   no geothermal heat flux'
    END IF
  END SUBROUTINE tra_bbc_init
END MODULE trabbc