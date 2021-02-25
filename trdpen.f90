MODULE trdpen
  USE oce
  USE dom_oce
  USE sbc_oce
  USE zdf_oce
  USE trd_oce
  USE eosbn2
  USE ldftra
  USE zdfddm
  USE phycst
  USE in_out_manager
  USE iom
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: trd_pen
  PUBLIC :: trd_pen_init
  INTEGER :: nkstp
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :, :) :: rab_pe
  CONTAINS
  INTEGER FUNCTION trd_pen_alloc()
    ALLOCATE(rab_pe(jpi, jpj, jpk, jpts), STAT = trd_pen_alloc)
    CALL mpp_sum('trdpen', trd_pen_alloc)
    IF (trd_pen_alloc /= 0) CALL ctl_stop('STOP', 'trd_pen_alloc: failed to allocate arrays')
  END FUNCTION trd_pen_alloc
  SUBROUTINE trd_pen(ptrdx, ptrdy, ktrd, kt, pdt)
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: ptrdx, ptrdy
    INTEGER, INTENT(IN) :: ktrd
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), INTENT(IN) :: pdt
    INTEGER :: jk
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: z2d
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zpe
    zpe(:, :, :) = 0._wp
    IF (kt /= nkstp) THEN
      nkstp = kt
      CALL eos_pen(tsn, rab_PE, zpe)
      CALL iom_put("alphaPE", rab_pe(:, :, :, jp_tem))
      CALL iom_put("betaPE", rab_pe(:, :, :, jp_sal))
      CALL iom_put("PEanom", zpe)
    END IF
    zpe(:, :, jpk) = 0._wp
    DO jk = 1, jpkm1
      zpe(:, :, jk) = (- (rab_n(:, :, jk, jp_tem) + rab_pe(:, :, jk, jp_tem)) * ptrdx(:, :, jk) + (rab_n(:, :, jk, jp_sal) + &
&rab_pe(:, :, jk, jp_sal)) * ptrdy(:, :, jk))
    END DO
    SELECT CASE (ktrd)
    CASE (jptra_xad)
      CALL iom_put("petrd_xad", zpe)
    CASE (jptra_yad)
      CALL iom_put("petrd_yad", zpe)
    CASE (jptra_zad)
      CALL iom_put("petrd_zad", zpe)
      IF (ln_linssh) THEN
        ALLOCATE(z2d(jpi, jpj))
        z2d(:, :) = wn(:, :, 1) * (- (rab_n(:, :, 1, jp_tem) + rab_pe(:, :, 1, jp_tem)) * tsn(:, :, 1, jp_tem) + (rab_n(:, :, 1, &
&jp_sal) + rab_pe(:, :, 1, jp_sal)) * tsn(:, :, 1, jp_sal)) / e3t_n(:, :, 1)
        CALL iom_put("petrd_sad", z2d)
        DEALLOCATE(z2d)
      END IF
    CASE (jptra_ldf)
      CALL iom_put("petrd_ldf", zpe)
    CASE (jptra_zdf)
      CALL iom_put("petrd_zdf", zpe)
    CASE (jptra_zdfp)
      CALL iom_put("petrd_zdfp", zpe)
    CASE (jptra_dmp)
      CALL iom_put("petrd_dmp", zpe)
    CASE (jptra_bbl)
      CALL iom_put("petrd_bbl", zpe)
    CASE (jptra_npc)
      CALL iom_put("petrd_npc", zpe)
    CASE (jptra_nsr)
      CALL iom_put("petrd_nsr", zpe)
    CASE (jptra_qsr)
      CALL iom_put("petrd_qsr", zpe)
    CASE (jptra_bbc)
      CALL iom_put("petrd_bbc", zpe)
    CASE (jptra_atf)
      CALL iom_put("petrd_atf", zpe)
    END SELECT
  END SUBROUTINE trd_pen
  SUBROUTINE trd_pen_init
    INTEGER :: ji, jj, jk
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'trd_pen_init : 3D Potential ENergy trends'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~~'
    END IF
    IF (trd_pen_alloc() /= 0) CALL ctl_stop('trd_pen_alloc: failed to allocate arrays')
    rab_pe(:, :, :, :) = 0._wp
    IF (.NOT. ln_linssh) CALL ctl_stop('trd_pen_init : PE trends not coded for variable volume')
    nkstp = nit000 - 1
  END SUBROUTINE trd_pen_init
END MODULE trdpen