MODULE usrdef_fmask
  USE oce
  USE dom_oce
  USE in_out_manager
  USE lbclnk
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: usr_def_fmask
  CONTAINS
  SUBROUTINE usr_def_fmask(cd_cfg, kcfg, pfmsk)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cd_cfg
    INTEGER, INTENT(IN) :: kcfg
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pfmsk
    INTEGER :: iif, iil, ii0, ii1, ii
    INTEGER :: ijf, ijl, ij0, ij1
    INTEGER :: isrow
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    IF (TRIM(cd_cfg) == "orca") THEN
      SELECT CASE (kcfg)
      CASE (2)
        CALL profile_psy_data0 % PreStart('usr_def_fmask', 'r0', 0, 0)
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) 'usr_def_fmask : ORCA_R2: increase lateral friction near the following straits:'
        IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~~~'
        IF (lwp) WRITE(numout, FMT = *) '      Gibraltar '
        CALL profile_psy_data0 % PostEnd
        !$ACC KERNELS
        ij0 = 101
        ij1 = 101
        ii0 = 139
        ii1 = 140
        pfmsk(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1), 1 : jpk) = 0.5_wp
        ij0 = 102
        ij1 = 102
        ii0 = 139
        ii1 = 140
        pfmsk(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1), 1 : jpk) = 0.5_wp
        !$ACC END KERNELS
        IF (lwp) WRITE(numout, FMT = *) '      Bab el Mandeb '
        !$ACC KERNELS
        ij0 = 87
        ij1 = 88
        ii0 = 160
        ii1 = 160
        pfmsk(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1), 1 : jpk) = 1._wp
        ij0 = 88
        ij1 = 88
        ii0 = 159
        ii1 = 159
        pfmsk(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1), 1 : jpk) = 1._wp
        !$ACC END KERNELS
      CASE (1)
        CALL profile_psy_data1 % PreStart('usr_def_fmask', 'r1', 0, 0)
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) 'usr_def_fmask : ORCA_R1: increase lateral friction near the following straits:'
        IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~~~'
        isrow = 332 - jpjglo
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) '   orca_r1: increase friction near the following straits : '
        IF (lwp) WRITE(numout, FMT = *) '      Gibraltar '
        CALL profile_psy_data1 % PostEnd
        !$ACC KERNELS
        ii0 = 282
        ii1 = 283
        ij0 = 241 - isrow
        ij1 = 241 - isrow
        pfmsk(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1), 1 : jpk) = 2._wp
        !$ACC END KERNELS
        IF (lwp) WRITE(numout, FMT = *) '      Bhosporus '
        !$ACC KERNELS
        ii0 = 314
        ii1 = 315
        ij0 = 248 - isrow
        ij1 = 248 - isrow
        pfmsk(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1), 1 : jpk) = 2._wp
        !$ACC END KERNELS
        IF (lwp) WRITE(numout, FMT = *) '      Makassar (Top) '
        !$ACC KERNELS
        ii0 = 48
        ii1 = 48
        ij0 = 189 - isrow
        ij1 = 190 - isrow
        pfmsk(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1), 1 : jpk) = 3._wp
        !$ACC END KERNELS
        IF (lwp) WRITE(numout, FMT = *) '      Lombok '
        !$ACC KERNELS
        ii0 = 44
        ii1 = 44
        ij0 = 164 - isrow
        ij1 = 165 - isrow
        pfmsk(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1), 1 : jpk) = 2._wp
        !$ACC END KERNELS
        IF (lwp) WRITE(numout, FMT = *) '      Ombai '
        !$ACC KERNELS
        ii0 = 53
        ii1 = 53
        ij0 = 164 - isrow
        ij1 = 165 - isrow
        pfmsk(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1), 1 : jpk) = 2._wp
        !$ACC END KERNELS
        IF (lwp) WRITE(numout, FMT = *) '      Timor Passage '
        !$ACC KERNELS
        ii0 = 56
        ii1 = 56
        ij0 = 164 - isrow
        ij1 = 165 - isrow
        pfmsk(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1), 1 : jpk) = 2._wp
        !$ACC END KERNELS
        IF (lwp) WRITE(numout, FMT = *) '      West Halmahera '
        !$ACC KERNELS
        ii0 = 58
        ii1 = 58
        ij0 = 181 - isrow
        ij1 = 182 - isrow
        pfmsk(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1), 1 : jpk) = 3._wp
        !$ACC END KERNELS
        IF (lwp) WRITE(numout, FMT = *) '      East Halmahera '
        !$ACC KERNELS
        ii0 = 55
        ii1 = 55
        ij0 = 181 - isrow
        ij1 = 182 - isrow
        pfmsk(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1), 1 : jpk) = 3._wp
        !$ACC END KERNELS
      CASE DEFAULT
        CALL profile_psy_data2 % PreStart('usr_def_fmask', 'r2', 0, 0)
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) 'usr_def_fmask : ORCA_R', kcfg, ' : NO alteration of fmask in specific straits '
        IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~~~'
        CALL profile_psy_data2 % PostEnd
      END SELECT
    ELSE
      CALL profile_psy_data3 % PreStart('usr_def_fmask', 'r3', 0, 0)
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'usr_def_fmask : NO alteration of fmask in specific straits '
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~~~'
      CALL profile_psy_data3 % PostEnd
    END IF
    CALL lbc_lnk('usrdef_fmask', pfmsk, 'F', 1._wp)
  END SUBROUTINE usr_def_fmask
END MODULE usrdef_fmask