MODULE traadv_fct
  USE oce
  USE dom_oce
  USE trc_oce
  USE trd_oce
  USE trdtra
  USE diaptr
  USE diaar5
  USE phycst, ONLY: rau0_rcp
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lbclnk
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_adv_fct
  PUBLIC :: interp_4th_cpt
  LOGICAL :: l_trd
  LOGICAL :: l_ptr
  LOGICAL :: l_hst
  REAL(KIND = wp) :: r1_6 = 1._wp / 6._wp
  INTEGER, PARAMETER :: np_NH = 0
  INTEGER, PARAMETER :: np_CEN2 = 1
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: zbetup, zbetdo, zbup, zbdo
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: zwi, zwx, zwy, zwz, ztu, ztv, zltu, zltv, ztw
  CONTAINS
  SUBROUTINE tra_adv_fct(kt, kit000, cdtype, p2dt, pun, pvn, pwn, ptb, ptn, pta, kjpt, kn_fct_h, kn_fct_v)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kit000
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    INTEGER, INTENT(IN) :: kjpt
    INTEGER, INTENT(IN) :: kn_fct_h
    INTEGER, INTENT(IN) :: kn_fct_v
    REAL(KIND = wp), INTENT(IN) :: p2dt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pun, pvn, pwn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: ptb, ptn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    INTEGER :: ji, jj, jk, jn
    REAL(KIND = wp) :: ztra
    REAL(KIND = wp) :: zfp_ui, zfp_vj, zfp_wk, zC2t_u, zC4t_u
    REAL(KIND = wp) :: zfm_ui, zfm_vj, zfm_wk, zC2t_v, zC4t_v
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ztrdx, ztrdy, ztrdz, zptry
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('tra_adv_fct', 'r0', 0, 0)
    IF (kt == kit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'tra_adv_fct : FCT advection scheme on ', cdtype
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~'
      ALLOCATE(zbetup(jpi, jpj, jpk), zbetdo(jpi, jpj, jpk), zbup(jpi, jpj, jpk), zbdo(jpi, jpj, jpk))
      ALLOCATE(zwi(jpi, jpj, jpk), zwx(jpi, jpj, jpk), zwy(jpi, jpj, jpk), zwz(jpi, jpj, jpk), ztu(jpi, jpj, jpk), ztv(jpi, jpj, &
&jpk), zltu(jpi, jpj, jpk), zltv(jpi, jpj, jpk), ztw(jpi, jpj, jpk))
    END IF
    l_trd = .FALSE.
    l_hst = .FALSE.
    l_ptr = .FALSE.
    IF ((cdtype == 'TRA' .AND. l_trdtra) .OR. (cdtype == 'TRC' .AND. l_trdtrc)) l_trd = .TRUE.
    IF (cdtype == 'TRA' .AND. ln_diaptr) l_ptr = .TRUE.
    IF (cdtype == 'TRA' .AND. (iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR. iom_use("uadv_salttr") .OR. &
&iom_use("vadv_salttr"))) l_hst = .TRUE.
    CALL profile_psy_data0 % PostEnd
    IF (l_trd .OR. l_hst) THEN
      ALLOCATE(ztrdx(jpi, jpj, jpk), ztrdy(jpi, jpj, jpk), ztrdz(jpi, jpj, jpk))
      !$ACC KERNELS
      ztrdx(:, :, :) = 0._wp
      ztrdy(:, :, :) = 0._wp
      ztrdz(:, :, :) = 0._wp
      !$ACC END KERNELS
    END IF
    IF (l_ptr) THEN
      ALLOCATE(zptry(jpi, jpj, jpk))
      !$ACC KERNELS
      zptry(:, :, :) = 0._wp
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    zwz(:, :, 1) = 0._wp
    zwx(:, :, jpk) = 0._wp
    zwy(:, :, jpk) = 0._wp
    zwz(:, :, jpk) = 0._wp
    zwi(:, :, :) = 0._wp
    !$ACC END KERNELS
    DO jn = 1, kjpt
      !$ACC KERNELS
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zfp_ui = pun(ji, jj, jk) + ABS(pun(ji, jj, jk))
            zfm_ui = pun(ji, jj, jk) - ABS(pun(ji, jj, jk))
            zfp_vj = pvn(ji, jj, jk) + ABS(pvn(ji, jj, jk))
            zfm_vj = pvn(ji, jj, jk) - ABS(pvn(ji, jj, jk))
            zwx(ji, jj, jk) = 0.5 * (zfp_ui * ptb(ji, jj, jk, jn) + zfm_ui * ptb(ji + 1, jj, jk, jn))
            zwy(ji, jj, jk) = 0.5 * (zfp_vj * ptb(ji, jj, jk, jn) + zfm_vj * ptb(ji, jj + 1, jk, jn))
          END DO
        END DO
      END DO
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            zfp_wk = pwn(ji, jj, jk) + ABS(pwn(ji, jj, jk))
            zfm_wk = pwn(ji, jj, jk) - ABS(pwn(ji, jj, jk))
            zwz(ji, jj, jk) = 0.5 * (zfp_wk * ptb(ji, jj, jk, jn) + zfm_wk * ptb(ji, jj, jk - 1, jn)) * wmask(ji, jj, jk)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      IF (ln_linssh) THEN
        IF (ln_isfcav) THEN
          !$ACC KERNELS
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 1, jpj
            DO ji = 1, jpi
              zwz(ji, jj, mikt(ji, jj)) = pwn(ji, jj, mikt(ji, jj)) * ptb(ji, jj, mikt(ji, jj), jn)
            END DO
          END DO
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          zwz(:, :, 1) = pwn(:, :, 1) * ptb(:, :, 1, jn)
          !$ACC END KERNELS
        END IF
      END IF
      !$ACC KERNELS
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            ztra = - (zwx(ji, jj, jk) - zwx(ji - 1, jj, jk) + zwy(ji, jj, jk) - zwy(ji, jj - 1, jk) + zwz(ji, jj, jk) - zwz(ji, &
&jj, jk + 1)) * r1_e1e2t(ji, jj)
            pta(ji, jj, jk, jn) = pta(ji, jj, jk, jn) + ztra / e3t_n(ji, jj, jk) * tmask(ji, jj, jk)
            zwi(ji, jj, jk) = (e3t_b(ji, jj, jk) * ptb(ji, jj, jk, jn) + p2dt * ztra) / e3t_a(ji, jj, jk) * tmask(ji, jj, jk)
          END DO
        END DO
      END DO
      IF (l_trd .OR. l_hst) THEN
        ztrdx(:, :, :) = zwx(:, :, :)
        ztrdy(:, :, :) = zwy(:, :, :)
        ztrdz(:, :, :) = zwz(:, :, :)
      END IF
      IF (l_ptr) zptry(:, :, :) = zwy(:, :, :)
      !$ACC END KERNELS
      SELECT CASE (kn_fct_h)
      CASE (2)
        !$ACC KERNELS
        DO jk = 1, jpkm1
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              zwx(ji, jj, jk) = 0.5_wp * pun(ji, jj, jk) * (ptn(ji, jj, jk, jn) + ptn(ji + 1, jj, jk, jn)) - zwx(ji, jj, jk)
              zwy(ji, jj, jk) = 0.5_wp * pvn(ji, jj, jk) * (ptn(ji, jj, jk, jn) + ptn(ji, jj + 1, jk, jn)) - zwy(ji, jj, jk)
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      CASE (4)
        !$ACC KERNELS
        zltu(:, :, jpk) = 0._wp
        zltv(:, :, jpk) = 0._wp
        !$ACC END KERNELS
        DO jk = 1, jpkm1
          !$ACC KERNELS
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              ztu(ji, jj, jk) = (ptn(ji + 1, jj, jk, jn) - ptn(ji, jj, jk, jn)) * umask(ji, jj, jk)
              ztv(ji, jj, jk) = (ptn(ji, jj + 1, jk, jn) - ptn(ji, jj, jk, jn)) * vmask(ji, jj, jk)
            END DO
          END DO
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zltu(ji, jj, jk) = (ztu(ji, jj, jk) + ztu(ji - 1, jj, jk)) * r1_6
              zltv(ji, jj, jk) = (ztv(ji, jj, jk) + ztv(ji, jj - 1, jk)) * r1_6
            END DO
          END DO
          !$ACC END KERNELS
        END DO
        CALL lbc_lnk_multi('traadv_fct', zltu, 'T', 1., zltv, 'T', 1.)
        !$ACC KERNELS
        DO jk = 1, jpkm1
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              zC2t_u = ptn(ji, jj, jk, jn) + ptn(ji + 1, jj, jk, jn)
              zC2t_v = ptn(ji, jj, jk, jn) + ptn(ji, jj + 1, jk, jn)
              zwx(ji, jj, jk) = 0.5_wp * pun(ji, jj, jk) * (zC2t_u + zltu(ji, jj, jk) - zltu(ji + 1, jj, jk)) - zwx(ji, jj, jk)
              zwy(ji, jj, jk) = 0.5_wp * pvn(ji, jj, jk) * (zC2t_v + zltv(ji, jj, jk) - zltv(ji, jj + 1, jk)) - zwy(ji, jj, jk)
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      CASE (41)
        !$ACC KERNELS
        ztu(:, :, jpk) = 0._wp
        ztv(:, :, jpk) = 0._wp
        DO jk = 1, jpkm1
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              ztu(ji, jj, jk) = (ptn(ji + 1, jj, jk, jn) - ptn(ji, jj, jk, jn)) * umask(ji, jj, jk)
              ztv(ji, jj, jk) = (ptn(ji, jj + 1, jk, jn) - ptn(ji, jj, jk, jn)) * vmask(ji, jj, jk)
            END DO
          END DO
        END DO
        !$ACC END KERNELS
        CALL lbc_lnk_multi('traadv_fct', ztu, 'U', - 1., ztv, 'V', - 1.)
        !$ACC KERNELS
        DO jk = 1, jpkm1
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zC2t_u = ptn(ji, jj, jk, jn) + ptn(ji + 1, jj, jk, jn)
              zC2t_v = ptn(ji, jj, jk, jn) + ptn(ji, jj + 1, jk, jn)
              zC4t_u = zC2t_u + r1_6 * (ztu(ji - 1, jj, jk) - ztu(ji + 1, jj, jk))
              zC4t_v = zC2t_v + r1_6 * (ztv(ji, jj - 1, jk) - ztv(ji, jj + 1, jk))
              zwx(ji, jj, jk) = 0.5_wp * pun(ji, jj, jk) * zC4t_u - zwx(ji, jj, jk)
              zwy(ji, jj, jk) = 0.5_wp * pvn(ji, jj, jk) * zC4t_v - zwy(ji, jj, jk)
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      END SELECT
      SELECT CASE (kn_fct_v)
      CASE (2)
        !$ACC KERNELS
        DO jk = 2, jpkm1
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zwz(ji, jj, jk) = (pwn(ji, jj, jk) * 0.5_wp * (ptn(ji, jj, jk, jn) + ptn(ji, jj, jk - 1, jn)) - zwz(ji, jj, jk)) * &
&wmask(ji, jj, jk)
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      CASE (4)
        CALL interp_4th_cpt(ptn(:, :, :, jn), ztw)
        !$ACC KERNELS
        DO jk = 2, jpkm1
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zwz(ji, jj, jk) = (pwn(ji, jj, jk) * ztw(ji, jj, jk) - zwz(ji, jj, jk)) * wmask(ji, jj, jk)
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      END SELECT
      IF (ln_linssh) THEN
        !$ACC KERNELS
        zwz(:, :, 1) = 0._wp
        !$ACC END KERNELS
      END IF
      CALL profile_psy_data1 % PreStart('tra_adv_fct', 'r1', 0, 0)
      CALL lbc_lnk_multi('traadv_fct', zwi, 'T', 1., zwx, 'U', - 1., zwy, 'V', - 1., zwz, 'W', 1.)
      CALL nonosc(ptb(:, :, :, jn), zwx, zwy, zwz, zwi, p2dt)
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            pta(ji, jj, jk, jn) = pta(ji, jj, jk, jn) - (zwx(ji, jj, jk) - zwx(ji - 1, jj, jk) + zwy(ji, jj, jk) - zwy(ji, jj - 1, &
&jk) + zwz(ji, jj, jk) - zwz(ji, jj, jk + 1)) * r1_e1e2t(ji, jj) / e3t_n(ji, jj, jk)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      IF (l_trd .OR. l_hst) THEN
        !$ACC KERNELS
        ztrdx(:, :, :) = ztrdx(:, :, :) + zwx(:, :, :)
        ztrdy(:, :, :) = ztrdy(:, :, :) + zwy(:, :, :)
        ztrdz(:, :, :) = ztrdz(:, :, :) + zwz(:, :, :)
        !$ACC END KERNELS
        CALL profile_psy_data2 % PreStart('tra_adv_fct', 'r2', 0, 0)
        IF (l_trd) THEN
          CALL trd_tra(kt, cdtype, jn, jptra_xad, ztrdx, pun, ptn(:, :, :, jn))
          CALL trd_tra(kt, cdtype, jn, jptra_yad, ztrdy, pvn, ptn(:, :, :, jn))
          CALL trd_tra(kt, cdtype, jn, jptra_zad, ztrdz, pwn, ptn(:, :, :, jn))
        END IF
        IF (l_hst) CALL dia_ar5_hst(jn, 'adv', ztrdx(:, :, :), ztrdy(:, :, :))
        CALL profile_psy_data2 % PostEnd
      END IF
      IF (l_ptr) THEN
        !$ACC KERNELS
        zptry(:, :, :) = zptry(:, :, :) + zwy(:, :, :)
        !$ACC END KERNELS
        CALL dia_ptr_hst(jn, 'adv', zptry(:, :, :))
      END IF
    END DO
    CALL profile_psy_data3 % PreStart('tra_adv_fct', 'r3', 0, 0)
    IF (l_trd .OR. l_hst) THEN
      DEALLOCATE(ztrdx, ztrdy, ztrdz)
    END IF
    IF (l_ptr) THEN
      DEALLOCATE(zptry)
    END IF
    CALL profile_psy_data3 % PostEnd
  END SUBROUTINE tra_adv_fct
  SUBROUTINE nonosc(pbef, paa, pbb, pcc, paft, p2dt)
    REAL(KIND = wp), INTENT(IN) :: p2dt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pbef, paft
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: paa, pbb, pcc
    INTEGER :: ji, jj, jk
    INTEGER :: ikm1
    REAL(KIND = wp) :: zpos, zneg, zbt, za, zb, zc, zbig, zrtrn
    REAL(KIND = wp) :: zau, zbu, zcu, zav, zbv, zcv, zup, zdo
    !$ACC KERNELS
    zbig = 1.E+40_wp
    zrtrn = 1.E-15_wp
    zbetup(:, :, :) = 0._wp
    zbetdo(:, :, :) = 0._wp
    zbup = MAX(pbef * tmask - zbig * (1._wp - tmask), paft * tmask - zbig * (1._wp - tmask))
    zbdo = MIN(pbef * tmask + zbig * (1._wp - tmask), paft * tmask + zbig * (1._wp - tmask))
    DO jk = 1, jpkm1
      ikm1 = MAX(jk - 1, 1)
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zup = MAX(zbup(ji, jj, jk), zbup(ji - 1, jj, jk), zbup(ji + 1, jj, jk), zbup(ji, jj - 1, jk), zbup(ji, jj + 1, jk), &
&zbup(ji, jj, ikm1), zbup(ji, jj, jk + 1))
          zdo = MIN(zbdo(ji, jj, jk), zbdo(ji - 1, jj, jk), zbdo(ji + 1, jj, jk), zbdo(ji, jj - 1, jk), zbdo(ji, jj + 1, jk), &
&zbdo(ji, jj, ikm1), zbdo(ji, jj, jk + 1))
          zpos = MAX(0., paa(ji - 1, jj, jk)) - MIN(0., paa(ji, jj, jk)) + MAX(0., pbb(ji, jj - 1, jk)) - MIN(0., pbb(ji, jj, jk)) &
&+ MAX(0., pcc(ji, jj, jk + 1)) - MIN(0., pcc(ji, jj, jk))
          zneg = MAX(0., paa(ji, jj, jk)) - MIN(0., paa(ji - 1, jj, jk)) + MAX(0., pbb(ji, jj, jk)) - MIN(0., pbb(ji, jj - 1, jk)) &
&+ MAX(0., pcc(ji, jj, jk)) - MIN(0., pcc(ji, jj, jk + 1))
          zbt = e1e2t(ji, jj) * e3t_n(ji, jj, jk) / p2dt
          zbetup(ji, jj, jk) = (zup - paft(ji, jj, jk)) / (zpos + zrtrn) * zbt
          zbetdo(ji, jj, jk) = (paft(ji, jj, jk) - zdo) / (zneg + zrtrn) * zbt
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk_multi('traadv_fct', zbetup, 'T', 1., zbetdo, 'T', 1.)
    !$ACC KERNELS
    DO jk = 1, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zau = MIN(1._wp, zbetdo(ji, jj, jk), zbetup(ji + 1, jj, jk))
          zbu = MIN(1._wp, zbetup(ji, jj, jk), zbetdo(ji + 1, jj, jk))
          zcu = (0.5 + SIGN(0.5, paa(ji, jj, jk)))
          paa(ji, jj, jk) = paa(ji, jj, jk) * (zcu * zau + (1._wp - zcu) * zbu)
          zav = MIN(1._wp, zbetdo(ji, jj, jk), zbetup(ji, jj + 1, jk))
          zbv = MIN(1._wp, zbetup(ji, jj, jk), zbetdo(ji, jj + 1, jk))
          zcv = (0.5 + SIGN(0.5, pbb(ji, jj, jk)))
          pbb(ji, jj, jk) = pbb(ji, jj, jk) * (zcv * zav + (1._wp - zcv) * zbv)
          za = MIN(1., zbetdo(ji, jj, jk + 1), zbetup(ji, jj, jk))
          zb = MIN(1., zbetup(ji, jj, jk + 1), zbetdo(ji, jj, jk))
          zc = (0.5 + SIGN(0.5, pcc(ji, jj, jk + 1)))
          pcc(ji, jj, jk + 1) = pcc(ji, jj, jk + 1) * (zc * za + (1._wp - zc) * zb)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk_multi('traadv_fct', paa, 'U', - 1., pbb, 'V', - 1.)
  END SUBROUTINE nonosc
  SUBROUTINE interp_4th_cpt_org(pt_in, pt_out)
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pt_in
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(OUT) :: pt_out
    INTEGER :: ji, jj, jk
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwd, zwi, zws, zwrm, zwt
    !$ACC KERNELS
    DO jk = 3, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          zwd(ji, jj, jk) = 4._wp
          zwi(ji, jj, jk) = 1._wp
          zws(ji, jj, jk) = 1._wp
          zwrm(ji, jj, jk) = 3._wp * (pt_in(ji, jj, jk - 1) + pt_in(ji, jj, jk))
          IF (tmask(ji, jj, jk + 1) == 0._wp) THEN
            zwd(ji, jj, jk) = 1._wp
            zwi(ji, jj, jk) = 0._wp
            zws(ji, jj, jk) = 0._wp
            zwrm(ji, jj, jk) = 0.5 * (pt_in(ji, jj, jk - 1) + pt_in(ji, jj, jk))
          END IF
        END DO
      END DO
    END DO
    jk = 2
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zwd(ji, jj, jk) = 1._wp
        zwi(ji, jj, jk) = 0._wp
        zws(ji, jj, jk) = 0._wp
        zwrm(ji, jj, jk) = 0.5 * (pt_in(ji, jj, jk - 1) + pt_in(ji, jj, jk))
      END DO
    END DO
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zwt(ji, jj, 2) = zwd(ji, jj, 2)
      END DO
    END DO
    DO jk = 3, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          zwt(ji, jj, jk) = zwd(ji, jj, jk) - zwi(ji, jj, jk) * zws(ji, jj, jk - 1) / zwt(ji, jj, jk - 1)
        END DO
      END DO
    END DO
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        pt_out(ji, jj, 2) = zwrm(ji, jj, 2)
      END DO
    END DO
    DO jk = 3, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          pt_out(ji, jj, jk) = zwrm(ji, jj, jk) - zwi(ji, jj, jk) / zwt(ji, jj, jk - 1) * pt_out(ji, jj, jk - 1)
        END DO
      END DO
    END DO
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        pt_out(ji, jj, jpkm1) = pt_out(ji, jj, jpkm1) / zwt(ji, jj, jpkm1)
      END DO
    END DO
    DO jk = jpk - 2, 2, - 1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          pt_out(ji, jj, jk) = (pt_out(ji, jj, jk) - zws(ji, jj, jk) * pt_out(ji, jj, jk + 1)) / zwt(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE interp_4th_cpt_org
  SUBROUTINE interp_4th_cpt(pt_in, pt_out)
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pt_in
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(OUT) :: pt_out
    INTEGER :: ji, jj, jk
    INTEGER :: ikt, ikb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwd, zwi, zws, zwrm, zwt
    !$ACC KERNELS
    DO jk = 3, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zwd(ji, jj, jk) = 3._wp * wmask(ji, jj, jk) + 1._wp
          zwi(ji, jj, jk) = wmask(ji, jj, jk)
          zws(ji, jj, jk) = wmask(ji, jj, jk)
          zwrm(ji, jj, jk) = 3._wp * wmask(ji, jj, jk) * (pt_in(ji, jj, jk) + pt_in(ji, jj, jk - 1))
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    IF (ln_isfcav) THEN
      !$ACC KERNELS
      zwd(:, :, 2) = 1._wp
      zwi(:, :, 2) = 0._wp
      zws(:, :, 2) = 0._wp
      zwrm(:, :, 2) = 0._wp
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        ikt = mikt(ji, jj) + 1
        ikb = mbkt(ji, jj)
        zwd(ji, jj, ikt) = 1._wp
        zwi(ji, jj, ikt) = 0._wp
        zws(ji, jj, ikt) = 0._wp
        zwrm(ji, jj, ikt) = 0.5_wp * (pt_in(ji, jj, ikt - 1) + pt_in(ji, jj, ikt))
        zwd(ji, jj, ikb) = 1._wp
        zwi(ji, jj, ikb) = 0._wp
        zws(ji, jj, ikb) = 0._wp
        zwrm(ji, jj, ikb) = 0.5_wp * (pt_in(ji, jj, ikb - 1) + pt_in(ji, jj, ikb))
      END DO
    END DO
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zwt(ji, jj, 2) = zwd(ji, jj, 2)
      END DO
    END DO
    DO jk = 3, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zwt(ji, jj, jk) = zwd(ji, jj, jk) - zwi(ji, jj, jk) * zws(ji, jj, jk - 1) / zwt(ji, jj, jk - 1)
        END DO
      END DO
    END DO
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        pt_out(ji, jj, 2) = zwrm(ji, jj, 2)
      END DO
    END DO
    DO jk = 3, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          pt_out(ji, jj, jk) = zwrm(ji, jj, jk) - zwi(ji, jj, jk) / zwt(ji, jj, jk - 1) * pt_out(ji, jj, jk - 1)
        END DO
      END DO
    END DO
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        pt_out(ji, jj, jpkm1) = pt_out(ji, jj, jpkm1) / zwt(ji, jj, jpkm1)
      END DO
    END DO
    DO jk = jpk - 2, 2, - 1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          pt_out(ji, jj, jk) = (pt_out(ji, jj, jk) - zws(ji, jj, jk) * pt_out(ji, jj, jk + 1)) / zwt(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE interp_4th_cpt
  SUBROUTINE tridia_solver(pD, pU, pL, pRHS, pt_out, klev)
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pD, pU, PL
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pRHS
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: pt_out
    INTEGER, INTENT(IN) :: klev
    INTEGER :: ji, jj, jk
    INTEGER :: kstart
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwt
    !$ACC KERNELS
    kstart = 1 + klev
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zwt(ji, jj, kstart) = pD(ji, jj, kstart)
      END DO
    END DO
    DO jk = kstart + 1, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zwt(ji, jj, jk) = pD(ji, jj, jk) - pL(ji, jj, jk) * pU(ji, jj, jk - 1) / zwt(ji, jj, jk - 1)
        END DO
      END DO
    END DO
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        pt_out(ji, jj, kstart) = pRHS(ji, jj, kstart)
      END DO
    END DO
    DO jk = kstart + 1, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          pt_out(ji, jj, jk) = pRHS(ji, jj, jk) - pL(ji, jj, jk) / zwt(ji, jj, jk - 1) * pt_out(ji, jj, jk - 1)
        END DO
      END DO
    END DO
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        pt_out(ji, jj, jpkm1) = pt_out(ji, jj, jpkm1) / zwt(ji, jj, jpkm1)
      END DO
    END DO
    DO jk = jpk - 2, kstart, - 1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          pt_out(ji, jj, jk) = (pt_out(ji, jj, jk) - pU(ji, jj, jk) * pt_out(ji, jj, jk + 1)) / zwt(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE tridia_solver
END MODULE traadv_fct