MODULE tranxt
  USE oce
  USE dom_oce
  USE sbc_oce
  USE sbcrnf
  USE sbcisf
  USE zdf_oce
  USE domvvl
  USE trd_oce
  USE trdtra
  USE traqsr
  USE phycst
  USE ldftra
  USE ldfslp
  USE bdy_oce, ONLY: ln_bdy
  USE bdytra
  USE in_out_manager
  USE lbclnk
  USE prtctl
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_nxt
  PUBLIC :: tra_nxt_fix
  PUBLIC :: tra_nxt_vvl
  CONTAINS
  SUBROUTINE tra_nxt(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk, jn
    REAL(KIND = wp) :: zfact
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ztrdt, ztrds
    IF (ln_timing) CALL timing_start('tra_nxt')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'tra_nxt : achieve the time stepping by Asselin filter and array swap'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~'
    END IF
    CALL lbc_lnk_multi('tranxt', tsa(:, :, :, jp_tem), 'T', 1., tsa(:, :, :, jp_sal), 'T', 1.)
    IF (ln_bdy) CALL bdy_tra(kt)
    IF (neuler == 0 .AND. kt == nit000) THEN
      r2dt = rdt
    ELSE IF (kt <= nit000 + 1) THEN
      r2dt = 2._wp * rdt
    END IF
    IF (l_trdtra) THEN
      ALLOCATE(ztrdt(jpi, jpj, jpk), ztrds(jpi, jpj, jpk))
      ztrdt(:, :, jpk) = 0._wp
      ztrds(:, :, jpk) = 0._wp
      IF (ln_traldf_iso) THEN
        CALL trd_tra(kt, 'TRA', jp_tem, jptra_zdfp, ztrdt)
        CALL trd_tra(kt, 'TRA', jp_sal, jptra_zdfp, ztrds)
      END IF
      zfact = 1.0 / rdt
      !$OMP parallel default(shared), private(jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        ztrdt(:, :, jk) = (tsa(:, :, jk, jp_tem) * e3t_a(:, :, jk) / e3t_n(:, :, jk) - tsn(:, :, jk, jp_tem)) * zfact
        ztrds(:, :, jk) = (tsa(:, :, jk, jp_sal) * e3t_a(:, :, jk) / e3t_n(:, :, jk) - tsn(:, :, jk, jp_sal)) * zfact
      END DO
      !$OMP end do
      !$OMP end parallel
      CALL trd_tra(kt, 'TRA', jp_tem, jptra_tot, ztrdt)
      CALL trd_tra(kt, 'TRA', jp_sal, jptra_tot, ztrds)
      IF (ln_linssh) THEN
        ztrdt(:, :, :) = tsn(:, :, :, jp_tem)
        ztrds(:, :, :) = tsn(:, :, :, jp_sal)
      END IF
    END IF
    IF (neuler == 0 .AND. kt == nit000) THEN
      DO jn = 1, jpts
        !$OMP parallel default(shared), private(jk)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          tsn(:, :, jk, jn) = tsa(:, :, jk, jn)
        END DO
        !$OMP end do
        !$OMP end parallel
      END DO
      IF (l_trdtra .AND. .NOT. ln_linssh) THEN
        ztrdt(:, :, :) = 0._wp
        ztrds(:, :, :) = 0._wp
        CALL trd_tra(kt, 'TRA', jp_tem, jptra_atf, ztrdt)
        CALL trd_tra(kt, 'TRA', jp_sal, jptra_atf, ztrds)
      END IF
    ELSE
      IF (ln_linssh) THEN
        CALL tra_nxt_fix(kt, nit000, 'TRA', tsb, tsn, tsa, jpts)
      ELSE
        CALL tra_nxt_vvl(kt, nit000, rdt, 'TRA', tsb, tsn, tsa, sbc_tsc, sbc_tsc_b, jpts)
      END IF
      CALL lbc_lnk_multi('tranxt', tsb(:, :, :, jp_tem), 'T', 1., tsb(:, :, :, jp_sal), 'T', 1., tsn(:, :, :, jp_tem), 'T', 1., &
&tsn(:, :, :, jp_sal), 'T', 1., tsa(:, :, :, jp_tem), 'T', 1., tsa(:, :, :, jp_sal), 'T', 1.)
    END IF
    IF (l_trdtra .AND. ln_linssh) THEN
      zfact = 1._wp / r2dt
      !$OMP parallel default(shared), private(jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        ztrdt(:, :, jk) = (tsb(:, :, jk, jp_tem) - ztrdt(:, :, jk)) * zfact
        ztrds(:, :, jk) = (tsb(:, :, jk, jp_sal) - ztrds(:, :, jk)) * zfact
      END DO
      !$OMP end do
      !$OMP end parallel
      CALL trd_tra(kt, 'TRA', jp_tem, jptra_atf, ztrdt)
      CALL trd_tra(kt, 'TRA', jp_sal, jptra_atf, ztrds)
    END IF
    IF (l_trdtra) DEALLOCATE(ztrdt, ztrds)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = tsn(:, :, :, jp_tem), clinfo1 = ' nxt  - Tn: ', mask1 = tmask, tab3d_2 = tsn(:, :, :, &
&jp_sal), clinfo2 = ' Sn: ', mask2 = tmask)
    IF (ln_timing) CALL timing_stop('tra_nxt')
  END SUBROUTINE tra_nxt
  SUBROUTINE tra_nxt_fix(kt, kit000, cdtype, ptb, ptn, pta, kjpt)
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kit000
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    INTEGER, INTENT(IN) :: kjpt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: ptb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: ptn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    INTEGER :: ji, jj, jk, jn
    REAL(KIND = wp) :: ztn, ztd
    IF (kt == kit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'tra_nxt_fix : time stepping', cdtype
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~'
    END IF
    DO jn = 1, kjpt
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            ztn = ptn(ji, jj, jk, jn)
            ztd = pta(ji, jj, jk, jn) - 2._wp * ztn + ptb(ji, jj, jk, jn)
            ptb(ji, jj, jk, jn) = ztn + atfp * ztd
            ptn(ji, jj, jk, jn) = pta(ji, jj, jk, jn)
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE tra_nxt_fix
  SUBROUTINE tra_nxt_vvl(kt, kit000, p2dt, cdtype, ptb, ptn, pta, psbc_tc, psbc_tc_b, kjpt)
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kit000
    REAL(KIND = wp), INTENT(IN) :: p2dt
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    INTEGER, INTENT(IN) :: kjpt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: ptb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: ptn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt), INTENT(IN) :: psbc_tc
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt), INTENT(IN) :: psbc_tc_b
    LOGICAL :: ll_traqsr, ll_rnf, ll_isf
    INTEGER :: ji, jj, jk, jn
    REAL(KIND = wp) :: zfact, zfact1, ztc_a, ztc_n, ztc_b, ztc_f, ztc_d
    REAL(KIND = wp) :: zfact2, ze3t_b, ze3t_n, ze3t_a, ze3t_f, ze3t_d
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :, :) :: ztrd_atf
    IF (kt == kit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'tra_nxt_vvl : time stepping', cdtype
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~'
    END IF
    IF (cdtype == 'TRA') THEN
      ll_traqsr = ln_traqsr
      ll_rnf = ln_rnf
      ll_isf = ln_isf
    ELSE
      ll_traqsr = .FALSE.
      ll_rnf = .FALSE.
      ll_isf = .FALSE.
    END IF
    IF ((l_trdtra .AND. cdtype == 'TRA') .OR. (l_trdtrc .AND. cdtype == 'TRC')) THEN
      ALLOCATE(ztrd_atf(jpi, jpj, jpk, kjpt))
      ztrd_atf(:, :, :, :) = 0.0_wp
    END IF
    zfact = 1._wp / p2dt
    zfact1 = atfp * p2dt
    zfact2 = zfact1 * r1_rau0
    DO jn = 1, kjpt
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            ze3t_b = e3t_b(ji, jj, jk)
            ze3t_n = e3t_n(ji, jj, jk)
            ze3t_a = e3t_a(ji, jj, jk)
            ztc_b = ptb(ji, jj, jk, jn) * ze3t_b
            ztc_n = ptn(ji, jj, jk, jn) * ze3t_n
            ztc_a = pta(ji, jj, jk, jn) * ze3t_a
            ze3t_d = ze3t_a - 2. * ze3t_n + ze3t_b
            ztc_d = ztc_a - 2. * ztc_n + ztc_b
            ze3t_f = ze3t_n + atfp * ze3t_d
            ztc_f = ztc_n + atfp * ztc_d
            IF (jk == mikt(ji, jj)) THEN
              ze3t_f = ze3t_f - zfact2 * ((emp_b(ji, jj) - emp(ji, jj)) + (fwfisf_b(ji, jj) - fwfisf(ji, jj)))
              ztc_f = ztc_f - zfact1 * (psbc_tc(ji, jj, jn) - psbc_tc_b(ji, jj, jn))
            END IF
            IF (ln_rnf_depth) THEN
              IF (mikt(ji, jj) <= jk .AND. jk <= nk_rnf(ji, jj)) THEN
                ze3t_f = ze3t_f - zfact2 * (- (rnf_b(ji, jj) - rnf(ji, jj))) * (e3t_n(ji, jj, jk) / h_rnf(ji, jj))
              END IF
            ELSE
              IF (jk == mikt(ji, jj)) THEN
                ze3t_f = ze3t_f - zfact2 * (- (rnf_b(ji, jj) - rnf(ji, jj)))
              END IF
            END IF
            IF (ll_traqsr .AND. jn == jp_tem .AND. jk <= nksr) ztc_f = ztc_f - zfact1 * (qsr_hc(ji, jj, jk) - qsr_hc_b(ji, jj, jk))
            IF (ll_rnf .AND. jk <= nk_rnf(ji, jj)) ztc_f = ztc_f - zfact1 * (rnf_tsc(ji, jj, jn) - rnf_tsc_b(ji, jj, jn)) * &
&e3t_n(ji, jj, jk) / h_rnf(ji, jj)
            IF (ll_isf) THEN
              IF (jk >= misfkt(ji, jj) .AND. jk < misfkb(ji, jj)) ztc_f = ztc_f - zfact1 * (risf_tsc(ji, jj, jn) - risf_tsc_b(ji, &
&jj, jn)) * e3t_n(ji, jj, jk) * r1_hisf_tbl(ji, jj)
              IF (jk == misfkb(ji, jj)) ztc_f = ztc_f - zfact1 * (risf_tsc(ji, jj, jn) - risf_tsc_b(ji, jj, jn)) * e3t_n(ji, jj, &
&jk) * r1_hisf_tbl(ji, jj) * ralpha(ji, jj)
            END IF
            ze3t_f = 1.E0 / ze3t_f
            ptb(ji, jj, jk, jn) = ztc_f * ze3t_f
            ptn(ji, jj, jk, jn) = pta(ji, jj, jk, jn)
            IF ((l_trdtra .AND. cdtype == 'TRA') .OR. (l_trdtrc .AND. cdtype == 'TRC')) THEN
              ztrd_atf(ji, jj, jk, jn) = (ztc_f - ztc_n) * zfact / ze3t_n
            END IF
          END DO
        END DO
      END DO
    END DO
    IF ((l_trdtra .AND. cdtype == 'TRA') .OR. (l_trdtrc .AND. cdtype == 'TRC')) THEN
      IF (l_trdtra .AND. cdtype == 'TRA') THEN
        CALL trd_tra(kt, cdtype, jp_tem, jptra_atf, ztrd_atf(:, :, :, jp_tem))
        CALL trd_tra(kt, cdtype, jp_sal, jptra_atf, ztrd_atf(:, :, :, jp_sal))
      END IF
      IF (l_trdtrc .AND. cdtype == 'TRC') THEN
        DO jn = 1, kjpt
          CALL trd_tra(kt, cdtype, jn, jptra_atf, ztrd_atf(:, :, :, jn))
        END DO
      END IF
      DEALLOCATE(ztrd_atf)
    END IF
  END SUBROUTINE tra_nxt_vvl
END MODULE tranxt