MODULE bdytra
  USE oce
  USE dom_oce
  USE bdy_oce
  USE bdylib
  USE in_out_manager
  USE lbclnk
  USE lib_mpp, ONLY: ctl_stop
  USE timing
  IMPLICIT NONE
  PRIVATE
  TYPE, PUBLIC :: ztrabdy
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: tra
  END TYPE
  PUBLIC :: bdy_tra
  PUBLIC :: bdy_tra_dmp
  CONTAINS
  SUBROUTINE bdy_tra(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ib_bdy, jn, igrd
    TYPE(ztrabdy), DIMENSION(jpts) :: zdta
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bdy_tra', 'r0', 0, 0)
    igrd = 1
    DO ib_bdy = 1, nb_bdy
      zdta(1) % tra => dta_bdy(ib_bdy) % tem
      zdta(2) % tra => dta_bdy(ib_bdy) % sal
      DO jn = 1, jpts
        SELECT CASE (TRIM(cn_tra(ib_bdy)))
        CASE ('none')
          CYCLE
        CASE ('frs')
          CALL bdy_frs(idx_bdy(ib_bdy), tsa(:, :, :, jn), zdta(jn) % tra)
        CASE ('specified')
          CALL bdy_spe(idx_bdy(ib_bdy), tsa(:, :, :, jn), zdta(jn) % tra)
        CASE ('neumann')
          CALL bdy_nmn(idx_bdy(ib_bdy), igrd, tsa(:, :, :, jn))
        CASE ('orlanski')
          CALL bdy_orl(idx_bdy(ib_bdy), tsb(:, :, :, jn), tsa(:, :, :, jn), zdta(jn) % tra, ll_npo = .FALSE.)
        CASE ('orlanski_npo')
          CALL bdy_orl(idx_bdy(ib_bdy), tsb(:, :, :, jn), tsa(:, :, :, jn), zdta(jn) % tra, ll_npo = .TRUE.)
        CASE ('runoff')
          CALL bdy_rnf(idx_bdy(ib_bdy), tsa(:, :, :, jn), jn)
        CASE DEFAULT
          CALL ctl_stop('bdy_tra : unrecognised option for open boundaries for T and S')
        END SELECT
        CALL lbc_bdy_lnk('bdytra', tsa(:, :, :, jn), 'T', 1., ib_bdy)
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE bdy_tra
  SUBROUTINE bdy_rnf(idx, pta, jpa)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pta
    INTEGER, INTENT(IN) :: jpa
    REAL(KIND = wp) :: zwgt
    INTEGER :: ib, ik, igrd
    INTEGER :: ii, ij, ip, jp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('bdy_rnf', 'r0', 0, 0)
    igrd = 1
    CALL profile_psy_data0 % PostEnd
    DO ib = 1, idx % nblenrim(igrd)
      CALL profile_psy_data1 % PreStart('bdy_rnf', 'r1', 0, 0)
      ii = idx % nbi(ib, igrd)
      ij = idx % nbj(ib, igrd)
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      DO ik = 1, jpkm1
        ip = bdytmask(ii + 1, ij) - bdytmask(ii - 1, ij)
        jp = bdytmask(ii, ij + 1) - bdytmask(ii, ij - 1)
        IF (jpa == jp_tem) pta(ii, ij, ik) = pta(ii + ip, ij + jp, ik) * tmask(ii, ij, ik)
        IF (jpa == jp_sal) pta(ii, ij, ik) = 0.1 * tmask(ii, ij, ik)
      END DO
      !$ACC END KERNELS
    END DO
  END SUBROUTINE bdy_rnf
  SUBROUTINE bdy_tra_dmp(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp) :: zwgt
    REAL(KIND = wp) :: zta, zsa, ztime
    INTEGER :: ib, ik, igrd
    INTEGER :: ii, ij
    INTEGER :: ib_bdy
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bdy_tra_dmp', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('bdy_tra_dmp')
    DO ib_bdy = 1, nb_bdy
      IF (ln_tra_dmp(ib_bdy)) THEN
        igrd = 1
        DO ib = 1, idx_bdy(ib_bdy) % nblen(igrd)
          ii = idx_bdy(ib_bdy) % nbi(ib, igrd)
          ij = idx_bdy(ib_bdy) % nbj(ib, igrd)
          zwgt = idx_bdy(ib_bdy) % nbd(ib, igrd)
          DO ik = 1, jpkm1
            zta = zwgt * (dta_bdy(ib_bdy) % tem(ib, ik) - tsb(ii, ij, ik, jp_tem)) * tmask(ii, ij, ik)
            zsa = zwgt * (dta_bdy(ib_bdy) % sal(ib, ik) - tsb(ii, ij, ik, jp_sal)) * tmask(ii, ij, ik)
            tsa(ii, ij, ik, jp_tem) = tsa(ii, ij, ik, jp_tem) + zta
            tsa(ii, ij, ik, jp_sal) = tsa(ii, ij, ik, jp_sal) + zsa
          END DO
        END DO
      END IF
    END DO
    IF (ln_timing) CALL timing_stop('bdy_tra_dmp')
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE bdy_tra_dmp
END MODULE bdytra