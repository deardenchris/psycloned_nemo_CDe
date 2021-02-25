MODULE traldf
  USE oce
  USE dom_oce
  USE phycst
  USE ldftra
  USE ldfslp
  USE traldf_lap_blp
  USE traldf_iso
  USE traldf_triad
  USE trd_oce
  USE trdtra
  USE prtctl
  USE in_out_manager
  USE lib_mpp
  USE lbclnk
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_ldf
  PUBLIC :: tra_ldf_init
  CONTAINS
  SUBROUTINE tra_ldf(kt)
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ztrdt, ztrds
    IF (ln_timing) CALL timing_start('tra_ldf')
    IF (l_trdtra) THEN
      ALLOCATE(ztrdt(jpi, jpj, jpk), ztrds(jpi, jpj, jpk))
      ztrdt(:, :, :) = tsa(:, :, :, jp_tem)
      ztrds(:, :, :) = tsa(:, :, :, jp_sal)
    END IF
    SELECT CASE (nldf_tra)
    CASE (np_lap)
      CALL tra_ldf_lap(kt, nit000, 'TRA', ahtu, ahtv, gtsu, gtsv, gtui, gtvi, tsb, tsa, jpts, 1)
    CASE (np_lap_i)
      CALL tra_ldf_iso(kt, nit000, 'TRA', ahtu, ahtv, gtsu, gtsv, gtui, gtvi, tsb, tsb, tsa, jpts, 1)
    CASE (np_lap_it)
      CALL tra_ldf_triad(kt, nit000, 'TRA', ahtu, ahtv, gtsu, gtsv, gtui, gtvi, tsb, tsb, tsa, jpts, 1)
    CASE (np_blp, np_blp_i, np_blp_it)
      CALL tra_ldf_blp(kt, nit000, 'TRA', ahtu, ahtv, gtsu, gtsv, gtui, gtvi, tsb, tsa, jpts, nldf_tra)
    END SELECT
    IF (l_trdtra) THEN
      ztrdt(:, :, :) = tsa(:, :, :, jp_tem) - ztrdt(:, :, :)
      ztrds(:, :, :) = tsa(:, :, :, jp_sal) - ztrds(:, :, :)
      CALL trd_tra(kt, 'TRA', jp_tem, jptra_ldf, ztrdt)
      CALL trd_tra(kt, 'TRA', jp_sal, jptra_ldf, ztrds)
      DEALLOCATE(ztrdt, ztrds)
    END IF
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = tsa(:, :, :, jp_tem), clinfo1 = ' ldf  - Ta: ', mask1 = tmask, tab3d_2 = tsa(:, :, :, &
&jp_sal), clinfo2 = ' Sa: ', mask2 = tmask, clinfo3 = 'tra')
    IF (ln_timing) CALL timing_stop('tra_ldf')
  END SUBROUTINE tra_ldf
  SUBROUTINE tra_ldf_init
    INTEGER :: ioptio, ierr
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'tra_ldf_init : lateral tracer diffusive operator'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namtra_ldf: already read in ldftra module'
      WRITE(numout, FMT = *) '      see ldf_tra_init report for lateral mixing parameters'
      WRITE(numout, FMT = *)
      SELECT CASE (nldf_tra)
      CASE (np_no_ldf)
        WRITE(numout, FMT = *) '   ==>>>   NO lateral diffusion'
      CASE (np_lap)
        WRITE(numout, FMT = *) '   ==>>>   laplacian iso-level operator'
      CASE (np_lap_i)
        WRITE(numout, FMT = *) '   ==>>>   Rotated laplacian operator (standard)'
      CASE (np_lap_it)
        WRITE(numout, FMT = *) '   ==>>>   Rotated laplacian operator (triad)'
      CASE (np_blp)
        WRITE(numout, FMT = *) '   ==>>>   bilaplacian iso-level operator'
      CASE (np_blp_i)
        WRITE(numout, FMT = *) '   ==>>>   Rotated bilaplacian operator (standard)'
      CASE (np_blp_it)
        WRITE(numout, FMT = *) '   ==>>>   Rotated bilaplacian operator (triad)'
      END SELECT
    END IF
  END SUBROUTINE tra_ldf_init
END MODULE traldf