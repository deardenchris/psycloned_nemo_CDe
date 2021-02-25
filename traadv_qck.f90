MODULE traadv_qck
  USE oce
  USE dom_oce
  USE trc_oce
  USE trd_oce
  USE trdtra
  USE diaptr
  USE in_out_manager
  USE lib_mpp
  USE lbclnk
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_adv_qck
  REAL(KIND = wp) :: r1_6 = 1. / 6.
  LOGICAL :: l_trd
  LOGICAL :: l_ptr
  CONTAINS
  SUBROUTINE tra_adv_qck(kt, kit000, cdtype, p2dt, pun, pvn, pwn, ptb, ptn, pta, kjpt)
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kit000
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    INTEGER, INTENT(IN) :: kjpt
    REAL(KIND = wp), INTENT(IN) :: p2dt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pun, pvn, pwn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: ptb, ptn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    IF (kt == kit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'tra_adv_qck : 3rd order quickest advection scheme on ', cdtype
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      IF (lwp) WRITE(numout, FMT = *)
    END IF
    l_trd = .FALSE.
    l_ptr = .FALSE.
    IF ((cdtype == 'TRA' .AND. l_trdtra) .OR. (cdtype == 'TRC' .AND. l_trdtrc)) l_trd = .TRUE.
    IF (cdtype == 'TRA' .AND. ln_diaptr) l_ptr = .TRUE.
    CALL tra_adv_qck_i(kt, cdtype, p2dt, pun, ptb, ptn, pta, kjpt)
    CALL tra_adv_qck_j(kt, cdtype, p2dt, pvn, ptb, ptn, pta, kjpt)
    CALL tra_adv_cen2_k(kt, cdtype, pwn, ptn, pta, kjpt)
  END SUBROUTINE tra_adv_qck
  SUBROUTINE tra_adv_qck_i(kt, cdtype, p2dt, pun, ptb, ptn, pta, kjpt)
    INTEGER, INTENT(IN) :: kt
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    INTEGER, INTENT(IN) :: kjpt
    REAL(KIND = wp), INTENT(IN) :: p2dt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pun
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: ptb, ptn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    INTEGER :: ji, jj, jk, jn
    REAL(KIND = wp) :: ztra, zbtr, zdir, zdx, zmsk
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwx, zfu, zfc, zfd
    DO jn = 1, kjpt
      zfu(:, :, :) = 0._wp
      zfc(:, :, :) = 0._wp
      zfd(:, :, :) = 0._wp
      zwx(:, :, :) = 0._wp
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zfc(ji, jj, jk) = ptb(ji - 1, jj, jk, jn)
            zfd(ji, jj, jk) = ptb(ji + 1, jj, jk, jn)
          END DO
        END DO
      END DO
      CALL lbc_lnk_multi('traadv_qck', zfc(:, :, :), 'T', 1., zfd(:, :, :), 'T', 1.)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zdir = 0.5 + SIGN(0.5, pun(ji, jj, jk))
            zfu(ji, jj, jk) = zdir * zfc(ji, jj, jk) + (1. - zdir) * zfd(ji + 1, jj, jk)
          END DO
        END DO
      END DO
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zdir = 0.5 + SIGN(0.5, pun(ji, jj, jk))
            zdx = (zdir * e1t(ji, jj) + (1. - zdir) * e1t(ji + 1, jj)) * e2u(ji, jj) * e3u_n(ji, jj, jk)
            zwx(ji, jj, jk) = ABS(pun(ji, jj, jk)) * p2dt / zdx
            zfc(ji, jj, jk) = zdir * ptb(ji, jj, jk, jn) + (1. - zdir) * ptb(ji + 1, jj, jk, jn)
            zfd(ji, jj, jk) = zdir * ptb(ji + 1, jj, jk, jn) + (1. - zdir) * ptb(ji, jj, jk, jn)
          END DO
        END DO
      END DO
      CALL lbc_lnk_multi('traadv_qck', zfu(:, :, :), 'T', 1., zfd(:, :, :), 'T', 1., zfc(:, :, :), 'T', 1., zwx(:, :, :), 'T', 1.)
      CALL quickest(zfu, zfd, zfc, zwx)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zfu(ji, jj, jk) = tmask(ji - 1, jj, jk) + tmask(ji, jj, jk) + tmask(ji + 1, jj, jk) - 2.
          END DO
        END DO
      END DO
      CALL lbc_lnk('traadv_qck', zfu(:, :, :), 'T', 1.)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zdir = 0.5 + SIGN(0.5, pun(ji, jj, jk))
            zmsk = zdir * zfu(ji, jj, jk) + (1. - zdir) * zfu(ji + 1, jj, jk)
            zwx(ji, jj, jk) = zmsk * zwx(ji, jj, jk) + (1. - zmsk) * zfc(ji, jj, jk)
            zwx(ji, jj, jk) = zwx(ji, jj, jk) * pun(ji, jj, jk)
          END DO
        END DO
      END DO
      CALL lbc_lnk('traadv_qck', zwx(:, :, :), 'T', 1.)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zbtr = r1_e1e2t(ji, jj) / e3t_n(ji, jj, jk)
            ztra = - zbtr * (zwx(ji, jj, jk) - zwx(ji - 1, jj, jk))
            pta(ji, jj, jk, jn) = pta(ji, jj, jk, jn) + ztra
          END DO
        END DO
      END DO
      IF (l_trd) CALL trd_tra(kt, cdtype, jn, jptra_xad, zwx, pun, ptn(:, :, :, jn))
    END DO
  END SUBROUTINE tra_adv_qck_i
  SUBROUTINE tra_adv_qck_j(kt, cdtype, p2dt, pvn, ptb, ptn, pta, kjpt)
    INTEGER, INTENT(IN) :: kt
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    INTEGER, INTENT(IN) :: kjpt
    REAL(KIND = wp), INTENT(IN) :: p2dt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pvn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: ptb, ptn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    INTEGER :: ji, jj, jk, jn
    REAL(KIND = wp) :: ztra, zbtr, zdir, zdx, zmsk
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwy, zfu, zfc, zfd
    DO jn = 1, kjpt
      zfu(:, :, :) = 0.0
      zfc(:, :, :) = 0.0
      zfd(:, :, :) = 0.0
      zwy(:, :, :) = 0.0
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zfc(ji, jj, jk) = ptb(ji, jj - 1, jk, jn)
            zfd(ji, jj, jk) = ptb(ji, jj + 1, jk, jn)
          END DO
        END DO
      END DO
      CALL lbc_lnk_multi('traadv_qck', zfc(:, :, :), 'T', 1., zfd(:, :, :), 'T', 1.)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zdir = 0.5 + SIGN(0.5, pvn(ji, jj, jk))
            zfu(ji, jj, jk) = zdir * zfc(ji, jj, jk) + (1. - zdir) * zfd(ji, jj + 1, jk)
          END DO
        END DO
      END DO
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zdir = 0.5 + SIGN(0.5, pvn(ji, jj, jk))
            zdx = (zdir * e2t(ji, jj) + (1. - zdir) * e2t(ji, jj + 1)) * e1v(ji, jj) * e3v_n(ji, jj, jk)
            zwy(ji, jj, jk) = ABS(pvn(ji, jj, jk)) * p2dt / zdx
            zfc(ji, jj, jk) = zdir * ptb(ji, jj, jk, jn) + (1. - zdir) * ptb(ji, jj + 1, jk, jn)
            zfd(ji, jj, jk) = zdir * ptb(ji, jj + 1, jk, jn) + (1. - zdir) * ptb(ji, jj, jk, jn)
          END DO
        END DO
      END DO
      CALL lbc_lnk_multi('traadv_qck', zfu(:, :, :), 'T', 1., zfd(:, :, :), 'T', 1., zfc(:, :, :), 'T', 1., zwy(:, :, :), 'T', 1.)
      CALL quickest(zfu, zfd, zfc, zwy)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zfu(ji, jj, jk) = tmask(ji, jj - 1, jk) + tmask(ji, jj, jk) + tmask(ji, jj + 1, jk) - 2.
          END DO
        END DO
      END DO
      CALL lbc_lnk('traadv_qck', zfu(:, :, :), 'T', 1.)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zdir = 0.5 + SIGN(0.5, pvn(ji, jj, jk))
            zmsk = zdir * zfu(ji, jj, jk) + (1. - zdir) * zfu(ji, jj + 1, jk)
            zwy(ji, jj, jk) = zmsk * zwy(ji, jj, jk) + (1. - zmsk) * zfc(ji, jj, jk)
            zwy(ji, jj, jk) = zwy(ji, jj, jk) * pvn(ji, jj, jk)
          END DO
        END DO
      END DO
      CALL lbc_lnk('traadv_qck', zwy(:, :, :), 'T', 1.)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zbtr = r1_e1e2t(ji, jj) / e3t_n(ji, jj, jk)
            ztra = - zbtr * (zwy(ji, jj, jk) - zwy(ji, jj - 1, jk))
            pta(ji, jj, jk, jn) = pta(ji, jj, jk, jn) + ztra
          END DO
        END DO
      END DO
      IF (l_trd) CALL trd_tra(kt, cdtype, jn, jptra_yad, zwy, pvn, ptn(:, :, :, jn))
      IF (l_ptr) CALL dia_ptr_hst(jn, 'adv', zwy(:, :, :))
    END DO
  END SUBROUTINE tra_adv_qck_j
  SUBROUTINE tra_adv_cen2_k(kt, cdtype, pwn, ptn, pta, kjpt)
    INTEGER, INTENT(IN) :: kt
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    INTEGER, INTENT(IN) :: kjpt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pwn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: ptn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    INTEGER :: ji, jj, jk, jn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwz
    zwz(:, :, 1) = 0._wp
    zwz(:, :, jpk) = 0._wp
    DO jn = 1, kjpt
      DO jk = 2, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zwz(ji, jj, jk) = 0.5 * pwn(ji, jj, jk) * (ptn(ji, jj, jk - 1, jn) + ptn(ji, jj, jk, jn)) * wmask(ji, jj, jk)
          END DO
        END DO
      END DO
      IF (ln_linssh) THEN
        IF (ln_isfcav) THEN
          DO jj = 1, jpj
            DO ji = 1, jpi
              zwz(ji, jj, mikt(ji, jj)) = pwn(ji, jj, mikt(ji, jj)) * ptn(ji, jj, mikt(ji, jj), jn)
            END DO
          END DO
        ELSE
          zwz(:, :, 1) = pwn(:, :, 1) * ptn(:, :, 1, jn)
        END IF
      END IF
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            pta(ji, jj, jk, jn) = pta(ji, jj, jk, jn) - (zwz(ji, jj, jk) - zwz(ji, jj, jk + 1)) * r1_e1e2t(ji, jj) / e3t_n(ji, jj, &
&jk)
          END DO
        END DO
      END DO
      IF (l_trd) CALL trd_tra(kt, cdtype, jn, jptra_zad, zwz, pwn, ptn(:, :, :, jn))
    END DO
  END SUBROUTINE tra_adv_cen2_k
  SUBROUTINE quickest(pfu, pfd, pfc, puc)
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pfu
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pfd
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pfc
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: puc
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zcoef1, zcoef2, zcoef3
    REAL(KIND = wp) :: zc, zcurv, zfho
    DO jk = 1, jpkm1
      DO jj = 1, jpj
        DO ji = 1, jpi
          zc = puc(ji, jj, jk)
          zcurv = pfd(ji, jj, jk) + pfu(ji, jj, jk) - 2. * pfc(ji, jj, jk)
          zcoef1 = 0.5 * (pfc(ji, jj, jk) + pfd(ji, jj, jk))
          zcoef2 = 0.5 * zc * (pfd(ji, jj, jk) - pfc(ji, jj, jk))
          zcoef3 = (1. - (zc * zc)) * r1_6 * zcurv
          zfho = zcoef1 - zcoef2 - zcoef3
          zcoef1 = pfd(ji, jj, jk) - pfu(ji, jj, jk)
          zcoef2 = ABS(zcoef1)
          zcoef3 = ABS(zcurv)
          IF (zcoef3 >= zcoef2) THEN
            zfho = pfc(ji, jj, jk)
          ELSE
            zcoef3 = pfu(ji, jj, jk) + ((pfc(ji, jj, jk) - pfu(ji, jj, jk)) / MAX(zc, 1.E-9))
            IF (zcoef1 >= 0.) THEN
              zfho = MAX(pfc(ji, jj, jk), zfho)
              zfho = MIN(zfho, MIN(zcoef3, pfd(ji, jj, jk)))
            ELSE
              zfho = MIN(pfc(ji, jj, jk), zfho)
              zfho = MAX(zfho, MAX(zcoef3, pfd(ji, jj, jk)))
            END IF
          END IF
          puc(ji, jj, jk) = zfho
        END DO
      END DO
    END DO
  END SUBROUTINE quickest
END MODULE traadv_qck