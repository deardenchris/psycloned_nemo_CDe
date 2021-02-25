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
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ztrdu, ztrdv
    IF (ln_timing) CALL timing_start('dyn_ldf')
    IF (l_trddyn) THEN
      ALLOCATE(ztrdu(jpi, jpj, jpk), ztrdv(jpi, jpj, jpk))
      ztrdu(:, :, :) = ua(:, :, :)
      ztrdv(:, :, :) = va(:, :, :)
    END IF
    SELECT CASE (nldf_dyn)
    CASE (np_lap)
      CALL dyn_ldf_lap(kt, ub, vb, ua, va, 1)
    CASE (np_lap_i)
      CALL dyn_ldf_iso(kt)
    CASE (np_blp)
      CALL dyn_ldf_blp(kt, ub, vb, ua, va)
    END SELECT
    IF (l_trddyn) THEN
      ztrdu(:, :, :) = ua(:, :, :) - ztrdu(:, :, :)
      ztrdv(:, :, :) = va(:, :, :) - ztrdv(:, :, :)
      CALL trd_dyn(ztrdu, ztrdv, jpdyn_ldf, kt)
      DEALLOCATE(ztrdu, ztrdv)
    END IF
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = ua, clinfo1 = ' ldf  - Ua: ', mask1 = umask, tab3d_2 = va, clinfo2 = ' Va: ', mask2 = &
&vmask, clinfo3 = 'dyn')
    IF (ln_timing) CALL timing_stop('dyn_ldf')
  END SUBROUTINE dyn_ldf
  SUBROUTINE dyn_ldf_init
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dyn_ldf_init : Choice of the lateral diffusive operator on dynamics'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namdyn_ldf: already read in ldfdyn module'
      WRITE(numout, FMT = *) '      see ldf_dyn_init report for lateral mixing parameters'
      WRITE(numout, FMT = *)
      SELECT CASE (nldf_dyn)
      CASE (np_no_ldf)
        WRITE(numout, FMT = *) '   ==>>>   NO lateral viscosity'
      CASE (np_lap)
        WRITE(numout, FMT = *) '   ==>>>   iso-level laplacian operator'
      CASE (np_lap_i)
        WRITE(numout, FMT = *) '   ==>>>   rotated laplacian operator with iso-level background'
      CASE (np_blp)
        WRITE(numout, FMT = *) '   ==>>>   iso-level bi-laplacian operator'
      END SELECT
    END IF
  END SUBROUTINE dyn_ldf_init
END MODULE dynldf