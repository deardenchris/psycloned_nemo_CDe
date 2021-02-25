MODULE traadv_cen
  USE dom_oce
  USE eosbn2
  USE traadv_fct
  USE trd_oce
  USE trdtra
  USE diaptr
  USE diaar5
  USE in_out_manager
  USE iom
  USE trc_oce
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_adv_cen
  REAL(KIND = wp) :: r1_6 = 1._wp / 6._wp
  LOGICAL :: l_trd
  LOGICAL :: l_ptr
  LOGICAL :: l_hst
  CONTAINS
  SUBROUTINE tra_adv_cen(kt, kit000, cdtype, pun, pvn, pwn, ptn, pta, kjpt, kn_cen_h, kn_cen_v)
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kit000
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    INTEGER, INTENT(IN) :: kjpt
    INTEGER, INTENT(IN) :: kn_cen_h
    INTEGER, INTENT(IN) :: kn_cen_v
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pun, pvn, pwn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: ptn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    INTEGER :: ji, jj, jk, jn
    INTEGER :: ierr
    REAL(KIND = wp) :: zC2t_u, zC4t_u
    REAL(KIND = wp) :: zC2t_v, zC4t_v
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwx, zwy, zwz, ztu, ztv, ztw
    IF (kt == kit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'tra_adv_cen : centered advection scheme on ', cdtype, ' order h/v =', kn_cen_h, '/', kn_cen_v
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~~ '
    END IF
    l_trd = .FALSE.
    l_hst = .FALSE.
    l_ptr = .FALSE.
    IF ((cdtype == 'TRA' .AND. l_trdtra) .OR. (cdtype == 'TRC' .AND. l_trdtrc)) l_trd = .TRUE.
    IF (cdtype == 'TRA' .AND. ln_diaptr) l_ptr = .TRUE.
    IF (cdtype == 'TRA' .AND. (iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR. iom_use("uadv_salttr") .OR. &
&iom_use("vadv_salttr"))) l_hst = .TRUE.
    zwz(:, :, 1) = 0._wp
    zwz(:, :, jpk) = 0._wp
    DO jn = 1, kjpt
      SELECT CASE (kn_cen_h)
      CASE (2)
        !$OMP parallel default(shared), private(ji,jj,jk)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              zwx(ji, jj, jk) = 0.5_wp * pun(ji, jj, jk) * (ptn(ji, jj, jk, jn) + ptn(ji + 1, jj, jk, jn))
              zwy(ji, jj, jk) = 0.5_wp * pvn(ji, jj, jk) * (ptn(ji, jj, jk, jn) + ptn(ji, jj + 1, jk, jn))
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
      CASE (4)
        ztu(:, :, jpk) = 0._wp
        ztv(:, :, jpk) = 0._wp
        !$OMP parallel default(shared), private(ji,jj,jk)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              ztu(ji, jj, jk) = (ptn(ji + 1, jj, jk, jn) - ptn(ji, jj, jk, jn)) * umask(ji, jj, jk)
              ztv(ji, jj, jk) = (ptn(ji, jj + 1, jk, jn) - ptn(ji, jj, jk, jn)) * vmask(ji, jj, jk)
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
        CALL lbc_lnk_multi('traadv_cen', ztu, 'U', - 1., ztv, 'V', - 1.)
        !$OMP parallel default(shared), private(ji,jj,jk,zc2t_u,zc2t_v,zc4t_u,zc4t_v)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 2, jpjm1
            DO ji = 1, jpim1
              zC2t_u = ptn(ji, jj, jk, jn) + ptn(ji + 1, jj, jk, jn)
              zC2t_v = ptn(ji, jj, jk, jn) + ptn(ji, jj + 1, jk, jn)
              zC4t_u = zC2t_u + r1_6 * (ztu(ji - 1, jj, jk) - ztu(ji + 1, jj, jk))
              zC4t_v = zC2t_v + r1_6 * (ztv(ji, jj - 1, jk) - ztv(ji, jj + 1, jk))
              zwx(ji, jj, jk) = 0.5_wp * pun(ji, jj, jk) * zC4t_u
              zwy(ji, jj, jk) = 0.5_wp * pvn(ji, jj, jk) * zC4t_v
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
      CASE DEFAULT
        CALL ctl_stop('traadv_fct: wrong value for nn_fct')
      END SELECT
      SELECT CASE (kn_cen_v)
      CASE (2)
        !$OMP parallel default(shared), private(ji,jj,jk)
        !$OMP do schedule(static)
        DO jk = 2, jpk
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zwz(ji, jj, jk) = 0.5 * pwn(ji, jj, jk) * (ptn(ji, jj, jk, jn) + ptn(ji, jj, jk - 1, jn)) * wmask(ji, jj, jk)
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
      CASE (4)
        CALL interp_4th_cpt(ptn(:, :, :, jn), ztw)
        !$OMP parallel default(shared), private(ji,jj,jk)
        !$OMP do schedule(static)
        DO jk = 2, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zwz(ji, jj, jk) = pwn(ji, jj, jk) * ztw(ji, jj, jk) * wmask(ji, jj, jk)
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
      END SELECT
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
      !$OMP parallel default(shared), private(ji,jj,jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            pta(ji, jj, jk, jn) = pta(ji, jj, jk, jn) - (zwx(ji, jj, jk) - zwx(ji - 1, jj, jk) + zwy(ji, jj, jk) - zwy(ji, jj - 1, &
&jk) + zwz(ji, jj, jk) - zwz(ji, jj, jk + 1)) * r1_e1e2t(ji, jj) / e3t_n(ji, jj, jk)
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
      IF (l_trd) THEN
        CALL trd_tra(kt, cdtype, jn, jptra_xad, zwx, pun, ptn(:, :, :, jn))
        CALL trd_tra(kt, cdtype, jn, jptra_yad, zwy, pvn, ptn(:, :, :, jn))
        CALL trd_tra(kt, cdtype, jn, jptra_zad, zwz, pwn, ptn(:, :, :, jn))
      END IF
      IF (l_ptr) CALL dia_ptr_hst(jn, 'adv', zwy(:, :, :))
      IF (l_hst) CALL dia_ar5_hst(jn, 'adv', zwx(:, :, :), zwy(:, :, :))
    END DO
  END SUBROUTINE tra_adv_cen
END MODULE traadv_cen