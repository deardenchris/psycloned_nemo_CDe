MODULE traadv_ubs
  USE oce
  USE dom_oce
  USE trc_oce
  USE trd_oce
  USE traadv_fct
  USE trdtra
  USE diaptr
  USE diaar5
  USE iom
  USE in_out_manager
  USE lib_mpp
  USE lbclnk
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_adv_ubs
  LOGICAL :: l_trd
  LOGICAL :: l_ptr
  LOGICAL :: l_hst
  CONTAINS
  SUBROUTINE tra_adv_ubs(kt, kit000, cdtype, p2dt, pun, pvn, pwn, ptb, ptn, pta, kjpt, kn_ubs_v)
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kit000
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    INTEGER, INTENT(IN) :: kjpt
    INTEGER, INTENT(IN) :: kn_ubs_v
    REAL(KIND = wp), INTENT(IN) :: p2dt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pun, pvn, pwn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: ptb, ptn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    INTEGER :: ji, jj, jk, jn
    REAL(KIND = wp) :: ztra, zbtr, zcoef
    REAL(KIND = wp) :: zfp_ui, zfm_ui, zcenut, ztak, zfp_wk, zfm_wk
    REAL(KIND = wp) :: zfp_vj, zfm_vj, zcenvt, zeeu, zeev, z_hdivn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: ztu, ztv, zltu, zltv, zti, ztw
    IF (kt == kit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'tra_adv_ubs :  horizontal UBS advection scheme on ', cdtype
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~~'
    END IF
    l_trd = .FALSE.
    l_hst = .FALSE.
    l_ptr = .FALSE.
    IF ((cdtype == 'TRA' .AND. l_trdtra) .OR. (cdtype == 'TRC' .AND. l_trdtrc)) l_trd = .TRUE.
    IF (cdtype == 'TRA' .AND. ln_diaptr) l_ptr = .TRUE.
    IF (cdtype == 'TRA' .AND. (iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR. iom_use("uadv_salttr") .OR. iom_use("vadv_salttr"))) l_hst = .TRUE.
    ztw(:, :, 1) = 0._wp
    zltu(:, :, jpk) = 0._wp
    zltv(:, :, jpk) = 0._wp
    ztw(:, :, jpk) = 0._wp
    zti(:, :, jpk) = 0._wp
    DO jn = 1, kjpt
      !$OMP parallel default(shared), private(ji,jj,jk,zcoef,zeeu,zeev)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zeeu = e2_e1u(ji, jj) * e3u_n(ji, jj, jk) * umask(ji, jj, jk)
            zeev = e1_e2v(ji, jj) * e3v_n(ji, jj, jk) * vmask(ji, jj, jk)
            ztu(ji, jj, jk) = zeeu * (ptb(ji + 1, jj, jk, jn) - ptb(ji, jj, jk, jn))
            ztv(ji, jj, jk) = zeev * (ptb(ji, jj + 1, jk, jn) - ptb(ji, jj, jk, jn))
          END DO
        END DO
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zcoef = 1._wp / (6._wp * e3t_n(ji, jj, jk))
            zltu(ji, jj, jk) = (ztu(ji, jj, jk) - ztu(ji - 1, jj, jk)) * zcoef
            zltv(ji, jj, jk) = (ztv(ji, jj, jk) - ztv(ji, jj - 1, jk)) * zcoef
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
      CALL lbc_lnk('traadv_ubs', zltu, 'T', 1.)
      CALL lbc_lnk('traadv_ubs', zltv, 'T', 1.)
      !$OMP parallel default(shared), private(ji,jj,jk,zcenut,zcenvt,zfm_ui,zfm_vj,zfp_ui,zfp_vj)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zfp_ui = pun(ji, jj, jk) + ABS(pun(ji, jj, jk))
            zfm_ui = pun(ji, jj, jk) - ABS(pun(ji, jj, jk))
            zfp_vj = pvn(ji, jj, jk) + ABS(pvn(ji, jj, jk))
            zfm_vj = pvn(ji, jj, jk) - ABS(pvn(ji, jj, jk))
            zcenut = pun(ji, jj, jk) * (ptn(ji, jj, jk, jn) + ptn(ji + 1, jj, jk, jn))
            zcenvt = pvn(ji, jj, jk) * (ptn(ji, jj, jk, jn) + ptn(ji, jj + 1, jk, jn))
            ztu(ji, jj, jk) = 0.5 * (zcenut - zfp_ui * zltu(ji, jj, jk) - zfm_ui * zltu(ji + 1, jj, jk))
            ztv(ji, jj, jk) = 0.5 * (zcenvt - zfp_vj * zltv(ji, jj, jk) - zfm_vj * zltv(ji, jj + 1, jk))
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
      zltu(:, :, :) = pta(:, :, :, jn)
      !$OMP parallel default(shared), private(ji,jj,jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            pta(ji, jj, jk, jn) = pta(ji, jj, jk, jn) - (ztu(ji, jj, jk) - ztu(ji - 1, jj, jk) + ztv(ji, jj, jk) - ztv(ji, jj - 1, jk)) * r1_e1e2t(ji, jj) / e3t_n(ji, jj, jk)
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
      zltu(:, :, :) = pta(:, :, :, jn) - zltu(:, :, :)
      IF (l_trd) THEN
        CALL trd_tra(kt, cdtype, jn, jptra_xad, ztu, pun, ptn(:, :, :, jn))
        CALL trd_tra(kt, cdtype, jn, jptra_yad, ztv, pvn, ptn(:, :, :, jn))
      END IF
      IF (l_ptr) CALL dia_ptr_hst(jn, 'adv', ztv(:, :, :))
      IF (l_hst) CALL dia_ar5_hst(jn, 'adv', ztu(:, :, :), ztv(:, :, :))
      SELECT CASE (kn_ubs_v)
      CASE (2)
        IF (l_trd) zltv(:, :, :) = pta(:, :, :, jn)
        !$OMP parallel default(shared), private(ji,jj,jk,zfm_wk,zfp_wk)
        !$OMP do schedule(static)
        DO jk = 2, jpkm1
          DO jj = 1, jpj
            DO ji = 1, jpi
              zfp_wk = pwn(ji, jj, jk) + ABS(pwn(ji, jj, jk))
              zfm_wk = pwn(ji, jj, jk) - ABS(pwn(ji, jj, jk))
              ztw(ji, jj, jk) = 0.5_wp * (zfp_wk * ptb(ji, jj, jk, jn) + zfm_wk * ptb(ji, jj, jk - 1, jn)) * wmask(ji, jj, jk)
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
        IF (ln_linssh) THEN
          IF (ln_isfcav) THEN
            DO jj = 1, jpj
              DO ji = 1, jpi
                ztw(ji, jj, mikt(ji, jj)) = pwn(ji, jj, mikt(ji, jj)) * ptb(ji, jj, mikt(ji, jj), jn)
              END DO
            END DO
          ELSE
            ztw(:, :, 1) = pwn(:, :, 1) * ptb(:, :, 1, jn)
          END IF
        END IF
        !$OMP parallel default(shared), private(ji,jj,jk,ztak)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              ztak = - (ztw(ji, jj, jk) - ztw(ji, jj, jk + 1)) * r1_e1e2t(ji, jj) / e3t_n(ji, jj, jk)
              pta(ji, jj, jk, jn) = pta(ji, jj, jk, jn) + ztak
              zti(ji, jj, jk) = (ptb(ji, jj, jk, jn) + p2dt * (ztak + zltu(ji, jj, jk))) * tmask(ji, jj, jk)
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
        CALL lbc_lnk('traadv_ubs', zti, 'T', 1.)
        !$OMP parallel default(shared), private(ji,jj,jk)
        !$OMP do schedule(static)
        DO jk = 2, jpkm1
          DO jj = 1, jpj
            DO ji = 1, jpi
              ztw(ji, jj, jk) = (0.5_wp * pwn(ji, jj, jk) * (ptn(ji, jj, jk, jn) + ptn(ji, jj, jk - 1, jn)) - ztw(ji, jj, jk)) * wmask(ji, jj, jk)
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
        IF (ln_linssh) ztw(:, :, 1) = 0._wp
        CALL nonosc_z(ptb(:, :, :, jn), ztw, zti, p2dt)
      CASE (4)
        CALL interp_4th_cpt(ptn(:, :, :, jn), ztw)
        !$OMP parallel default(shared), private(ji,jj,jk)
        !$OMP do schedule(static)
        DO jk = 2, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              ztw(ji, jj, jk) = pwn(ji, jj, jk) * ztw(ji, jj, jk) * wmask(ji, jj, jk)
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
        IF (ln_linssh) ztw(:, :, 1) = pwn(:, :, 1) * ptn(:, :, 1, jn)
      END SELECT
      !$OMP parallel default(shared), private(ji,jj,jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            pta(ji, jj, jk, jn) = pta(ji, jj, jk, jn) - (ztw(ji, jj, jk) - ztw(ji, jj, jk + 1)) * r1_e1e2t(ji, jj) / e3t_n(ji, jj, jk)
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
      IF (l_trd) THEN
        !$OMP parallel default(shared), private(ji,jj,jk)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zltv(ji, jj, jk) = pta(ji, jj, jk, jn) - zltv(ji, jj, jk) + ptn(ji, jj, jk, jn) * (pwn(ji, jj, jk) - pwn(ji, jj, jk + 1)) * r1_e1e2t(ji, jj) / e3t_n(ji, jj, jk)
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
        CALL trd_tra(kt, cdtype, jn, jptra_zad, zltv)
      END IF
    END DO
  END SUBROUTINE tra_adv_ubs
  SUBROUTINE nonosc_z(pbef, pcc, paft, p2dt)
    REAL(KIND = wp), INTENT(IN) :: p2dt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: pbef
    REAL(KIND = wp), INTENT(INOUT), DIMENSION(jpi, jpj, jpk) :: paft
    REAL(KIND = wp), INTENT(INOUT), DIMENSION(jpi, jpj, jpk) :: pcc
    INTEGER :: ji, jj, jk
    INTEGER :: ikm1
    REAL(KIND = wp) :: zpos, zneg, zbt, za, zb, zc, zbig, zrtrn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zbetup, zbetdo
    zbig = 1.E+40_wp
    zrtrn = 1.E-15_wp
    zbetup(:, :, :) = 0._wp
    zbetdo(:, :, :) = 0._wp
    pbef(:, :, :) = pbef(:, :, :) * tmask(:, :, :) - zbig * (1.E0 - tmask(:, :, :))
    paft(:, :, :) = paft(:, :, :) * tmask(:, :, :) - zbig * (1.E0 - tmask(:, :, :))
    !$OMP parallel default(shared), private(ikm1,ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      ikm1 = MAX(jk - 1, 1)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zbetup(ji, jj, jk) = MAX(pbef(ji, jj, jk), paft(ji, jj, jk), pbef(ji, jj, ikm1), pbef(ji, jj, jk + 1), paft(ji, jj, ikm1), paft(ji, jj, jk + 1))
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    pbef(:, :, :) = pbef(:, :, :) * tmask(:, :, :) + zbig * (1.E0 - tmask(:, :, :))
    paft(:, :, :) = paft(:, :, :) * tmask(:, :, :) + zbig * (1.E0 - tmask(:, :, :))
    !$OMP parallel default(shared), private(ikm1,ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      ikm1 = MAX(jk - 1, 1)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zbetdo(ji, jj, jk) = MIN(pbef(ji, jj, jk), paft(ji, jj, jk), pbef(ji, jj, ikm1), pbef(ji, jj, jk + 1), paft(ji, jj, ikm1), paft(ji, jj, jk + 1))
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    pbef(:, :, :) = pbef(:, :, :) * tmask(:, :, :)
    paft(:, :, :) = paft(:, :, :) * tmask(:, :, :)
    !$OMP parallel default(shared), private(ji,jj,jk,za,zb,zbt,zc,zneg,zpos)
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zpos = MAX(0., pcc(ji, jj, jk + 1)) - MIN(0., pcc(ji, jj, jk))
          zneg = MAX(0., pcc(ji, jj, jk)) - MIN(0., pcc(ji, jj, jk + 1))
          zbt = e1e2t(ji, jj) * e3t_n(ji, jj, jk) / p2dt
          zbetup(ji, jj, jk) = (zbetup(ji, jj, jk) - paft(ji, jj, jk)) / (zpos + zrtrn) * zbt
          zbetdo(ji, jj, jk) = (paft(ji, jj, jk) - zbetdo(ji, jj, jk)) / (zneg + zrtrn) * zbt
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP do schedule(static)
    DO jk = 2, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          za = MIN(1., zbetdo(ji, jj, jk), zbetup(ji, jj, jk - 1))
          zb = MIN(1., zbetup(ji, jj, jk), zbetdo(ji, jj, jk - 1))
          zc = 0.5 * (1.E0 + SIGN(1.E0, pcc(ji, jj, jk)))
          pcc(ji, jj, jk) = pcc(ji, jj, jk) * (zc * za + (1.E0 - zc) * zb)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
  END SUBROUTINE nonosc_z
END MODULE traadv_ubs