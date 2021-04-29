MODULE dynldf
  USE oce
  USE dom_oce
  USE phycst
  USE ldfdyn
  USE dynldf_lap_blp
  USE dynldf_iso
  USE trd_oce
  USE trddyn
  USE prtctl
  USE in_out_manager
  USE lib_mpp
  USE lbclnk
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_ldf
  PUBLIC :: dyn_ldf_init
  CONTAINS
  SUBROUTINE dyn_ldf(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ztrdu, ztrdv
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('dyn_ldf', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('dyn_ldf')
    CALL profile_psy_data0 % PostEnd
    IF (l_trddyn) THEN
      ALLOCATE(ztrdu(jpi, jpj, jpk), ztrdv(jpi, jpj, jpk))
      !$ACC KERNELS
      ztrdu(:, :, :) = ua(:, :, :)
      ztrdv(:, :, :) = va(:, :, :)
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data1 % PreStart('dyn_ldf', 'r1', 0, 0)
    SELECT CASE (nldf_dyn)
    CASE (np_lap)
      CALL dyn_ldf_lap(kt, ub, vb, ua, va, 1)
    CASE (np_lap_i)
      CALL dyn_ldf_iso(kt)
    CASE (np_blp)
      CALL dyn_ldf_blp(kt, ub, vb, ua, va)
    END SELECT
    CALL profile_psy_data1 % PostEnd
    IF (l_trddyn) THEN
      !$ACC KERNELS
      ztrdu(:, :, :) = ua(:, :, :) - ztrdu(:, :, :)
      ztrdv(:, :, :) = va(:, :, :) - ztrdv(:, :, :)
      !$ACC END KERNELS
      CALL profile_psy_data2 % PreStart('dyn_ldf', 'r2', 0, 0)
      CALL trd_dyn(ztrdu, ztrdv, jpdyn_ldf, kt)
      DEALLOCATE(ztrdu, ztrdv)
      CALL profile_psy_data2 % PostEnd
    END IF
    CALL profile_psy_data3 % PreStart('dyn_ldf', 'r3', 0, 0)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = ua, clinfo1 = ' ldf  - Ua: ', mask1 = umask, tab3d_2 = va, clinfo2 = ' Va: ', mask2 = vmask, clinfo3 = 'dyn')
    IF (ln_timing) CALL timing_stop('dyn_ldf')
    CALL profile_psy_data3 % PostEnd
  END SUBROUTINE dyn_ldf
  SUBROUTINE dyn_ldf_init
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'dyn_ldf_init : Choice of the lateral diffusive operator on dynamics'
      WRITE(numout, *) '~~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist namdyn_ldf: already read in ldfdyn module'
      WRITE(numout, *) '      see ldf_dyn_init report for lateral mixing parameters'
      WRITE(numout, *)
      SELECT CASE (nldf_dyn)
      CASE (np_no_ldf)
        WRITE(numout, *) '   ==>>>   NO lateral viscosity'
      CASE (np_lap)
        WRITE(numout, *) '   ==>>>   iso-level laplacian operator'
      CASE (np_lap_i)
        WRITE(numout, *) '   ==>>>   rotated laplacian operator with iso-level background'
      CASE (np_blp)
        WRITE(numout, *) '   ==>>>   iso-level bi-laplacian operator'
      END SELECT
    END IF
  END SUBROUTINE dyn_ldf_init
END MODULE dynldf