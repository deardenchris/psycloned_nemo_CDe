MODULE trabbl
  USE oce
  USE dom_oce
  USE phycst
  USE eosbn2
  USE trd_oce
  USE trdtra
  USE iom
  USE in_out_manager
  USE lbclnk
  USE prtctl
  USE timing
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_bbl
  PUBLIC :: tra_bbl_init
  PUBLIC :: tra_bbl_dif
  PUBLIC :: tra_bbl_adv
  PUBLIC :: bbl
  LOGICAL, PUBLIC :: ln_trabbl
  INTEGER, PUBLIC :: nn_bbl_ldf
  INTEGER, PUBLIC :: nn_bbl_adv
  REAL(KIND = wp), PUBLIC :: rn_ahtbbl
  REAL(KIND = wp), PUBLIC :: rn_gambbl
  LOGICAL, PUBLIC :: l_bbl
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :), PUBLIC :: utr_bbl, vtr_bbl
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :), PUBLIC :: ahu_bbl, ahv_bbl
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:, :), PUBLIC :: mbku_d, mbkv_d
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:, :), PUBLIC :: mgrhu, mgrhv
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: ahu_bbl_0, ahv_bbl_0
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :), PUBLIC :: e3u_bbl_0, e3v_bbl_0
  CONTAINS
  INTEGER FUNCTION tra_bbl_alloc()
    ALLOCATE(utr_bbl(jpi, jpj), ahu_bbl(jpi, jpj), mbku_d(jpi, jpj), mgrhu(jpi, jpj), vtr_bbl(jpi, jpj), ahv_bbl(jpi, jpj), &
&mbkv_d(jpi, jpj), mgrhv(jpi, jpj), ahu_bbl_0(jpi, jpj), ahv_bbl_0(jpi, jpj), e3u_bbl_0(jpi, jpj), e3v_bbl_0(jpi, jpj), STAT = &
&tra_bbl_alloc)
    CALL mpp_sum('trabbl', tra_bbl_alloc)
    IF (tra_bbl_alloc > 0) CALL ctl_warn('tra_bbl_alloc: allocation of arrays failed.')
  END FUNCTION tra_bbl_alloc
  SUBROUTINE tra_bbl(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ztrdt, ztrds
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    IF (ln_timing) CALL timing_start('tra_bbl')
    IF (l_trdtra) THEN
      ALLOCATE(ztrdt(jpi, jpj, jpk), ztrds(jpi, jpj, jpk))
      !$ACC KERNELS
      ztrdt(:, :, :) = tsa(:, :, :, jp_tem)
      ztrds(:, :, :) = tsa(:, :, :, jp_sal)
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data0 % PreStart('tra_bbl', 'r0', 0, 0)
    IF (l_bbl) CALL bbl(kt, nit000, 'TRA')
    IF (nn_bbl_ldf == 1) THEN
      CALL tra_bbl_dif(tsb, tsa, jpts)
      IF (ln_ctl) CALL prt_ctl(tab3d_1 = tsa(:, :, :, jp_tem), clinfo1 = ' bbl_ldf  - Ta: ', mask1 = tmask, tab3d_2 = tsa(:, :, :, &
&jp_sal), clinfo2 = ' Sa: ', mask2 = tmask, clinfo3 = 'tra')
      CALL lbc_lnk_multi('trabbl', ahu_bbl, 'U', 1., ahv_bbl, 'V', 1.)
      CALL iom_put("ahu_bbl", ahu_bbl)
      CALL iom_put("ahv_bbl", ahv_bbl)
    END IF
    IF (nn_bbl_adv /= 0) THEN
      CALL tra_bbl_adv(tsb, tsa, jpts)
      IF (ln_ctl) CALL prt_ctl(tab3d_1 = tsa(:, :, :, jp_tem), clinfo1 = ' bbl_adv  - Ta: ', mask1 = tmask, tab3d_2 = tsa(:, :, :, &
&jp_sal), clinfo2 = ' Sa: ', mask2 = tmask, clinfo3 = 'tra')
      CALL lbc_lnk_multi('trabbl', utr_bbl, 'U', 1., vtr_bbl, 'V', 1.)
      CALL iom_put("uoce_bbl", utr_bbl)
      CALL iom_put("voce_bbl", vtr_bbl)
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (l_trdtra) THEN
      !$ACC KERNELS
      ztrdt(:, :, :) = tsa(:, :, :, jp_tem) - ztrdt(:, :, :)
      ztrds(:, :, :) = tsa(:, :, :, jp_sal) - ztrds(:, :, :)
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('tra_bbl', 'r1', 0, 0)
      CALL trd_tra(kt, 'TRA', jp_tem, jptra_bbl, ztrdt)
      CALL trd_tra(kt, 'TRA', jp_sal, jptra_bbl, ztrds)
      DEALLOCATE(ztrdt, ztrds)
      CALL profile_psy_data1 % PostEnd
    END IF
    IF (ln_timing) CALL timing_stop('tra_bbl')
  END SUBROUTINE tra_bbl
  SUBROUTINE tra_bbl_dif(ptb, pta, kjpt)
    INTEGER, INTENT(IN) :: kjpt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: ptb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    INTEGER :: ji, jj, jn
    INTEGER :: ik
    REAL(KIND = wp) :: zbtr
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zptb
    !$ACC KERNELS
    DO jn = 1, kjpt
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          ik = mbkt(ji, jj)
          zptb(ji, jj) = ptb(ji, jj, ik, jn)
        END DO
      END DO
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ik = mbkt(ji, jj)
          pta(ji, jj, ik, jn) = pta(ji, jj, ik, jn) + (ahu_bbl(ji, jj) * (zptb(ji + 1, jj) - zptb(ji, jj)) - ahu_bbl(ji - 1, jj) * &
&(zptb(ji, jj) - zptb(ji - 1, jj)) + ahv_bbl(ji, jj) * (zptb(ji, jj + 1) - zptb(ji, jj)) - ahv_bbl(ji, jj - 1) * (zptb(ji, jj) - &
&zptb(ji, jj - 1))) * r1_e1e2t(ji, jj) / e3t_n(ji, jj, ik)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE tra_bbl_dif
  SUBROUTINE tra_bbl_adv(ptb, pta, kjpt)
    INTEGER, INTENT(IN) :: kjpt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: ptb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    INTEGER :: ji, jj, jk, jn
    INTEGER :: iis, iid, ijs, ijd
    INTEGER :: ikus, ikud, ikvs, ikvd
    REAL(KIND = wp) :: zbtr, ztra
    REAL(KIND = wp) :: zu_bbl, zv_bbl
    !$ACC KERNELS
    DO jn = 1, kjpt
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          IF (utr_bbl(ji, jj) /= 0.E0) THEN
            iid = ji + MAX(0, mgrhu(ji, jj))
            iis = ji + 1 - MAX(0, mgrhu(ji, jj))
            ikud = mbku_d(ji, jj)
            ikus = mbku(ji, jj)
            zu_bbl = ABS(utr_bbl(ji, jj))
            zbtr = r1_e1e2t(iis, jj) / e3t_n(iis, jj, ikus)
            ztra = zu_bbl * (ptb(iid, jj, ikus, jn) - ptb(iis, jj, ikus, jn)) * zbtr
            pta(iis, jj, ikus, jn) = pta(iis, jj, ikus, jn) + ztra
            DO jk = ikus, ikud - 1
              zbtr = r1_e1e2t(iid, jj) / e3t_n(iid, jj, jk)
              ztra = zu_bbl * (ptb(iid, jj, jk + 1, jn) - ptb(iid, jj, jk, jn)) * zbtr
              pta(iid, jj, jk, jn) = pta(iid, jj, jk, jn) + ztra
            END DO
            zbtr = r1_e1e2t(iid, jj) / e3t_n(iid, jj, ikud)
            ztra = zu_bbl * (ptb(iis, jj, ikus, jn) - ptb(iid, jj, ikud, jn)) * zbtr
            pta(iid, jj, ikud, jn) = pta(iid, jj, ikud, jn) + ztra
          END IF
          IF (vtr_bbl(ji, jj) /= 0.E0) THEN
            ijd = jj + MAX(0, mgrhv(ji, jj))
            ijs = jj + 1 - MAX(0, mgrhv(ji, jj))
            ikvd = mbkv_d(ji, jj)
            ikvs = mbkv(ji, jj)
            zv_bbl = ABS(vtr_bbl(ji, jj))
            zbtr = r1_e1e2t(ji, ijs) / e3t_n(ji, ijs, ikvs)
            ztra = zv_bbl * (ptb(ji, ijd, ikvs, jn) - ptb(ji, ijs, ikvs, jn)) * zbtr
            pta(ji, ijs, ikvs, jn) = pta(ji, ijs, ikvs, jn) + ztra
            DO jk = ikvs, ikvd - 1
              zbtr = r1_e1e2t(ji, ijd) / e3t_n(ji, ijd, jk)
              ztra = zv_bbl * (ptb(ji, ijd, jk + 1, jn) - ptb(ji, ijd, jk, jn)) * zbtr
              pta(ji, ijd, jk, jn) = pta(ji, ijd, jk, jn) + ztra
            END DO
            zbtr = r1_e1e2t(ji, ijd) / e3t_n(ji, ijd, ikvd)
            ztra = zv_bbl * (ptb(ji, ijs, ikvs, jn) - ptb(ji, ijd, ikvd, jn)) * zbtr
            pta(ji, ijd, ikvd, jn) = pta(ji, ijd, ikvd, jn) + ztra
          END IF
        END DO
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE tra_bbl_adv
  SUBROUTINE bbl(kt, kit000, cdtype)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kit000
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    INTEGER :: ji, jj
    INTEGER :: ik
    INTEGER :: iis, iid, ikus, ikud
    INTEGER :: ijs, ijd, ikvs, ikvd
    REAL(KIND = wp) :: za, zb, zgdrho
    REAL(KIND = wp) :: zsign, zsigna, zgbbl
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpts) :: zts, zab
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zub, zvb, zdep
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bbl', 'r0', 0, 0)
    IF (kt == kit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'trabbl:bbl : Compute bbl velocities and diffusive coefficients in ', cdtype
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~'
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        ik = mbkt(ji, jj)
        zts(ji, jj, jp_tem) = tsb(ji, jj, ik, jp_tem)
        zts(ji, jj, jp_sal) = tsb(ji, jj, ik, jp_sal)
        zdep(ji, jj) = gdept_n(ji, jj, ik)
        zub(ji, jj) = un(ji, jj, mbku(ji, jj))
        zvb(ji, jj) = vn(ji, jj, mbkv(ji, jj))
      END DO
    END DO
    !$ACC END KERNELS
    CALL eos_rab(zts, zdep, zab)
    !$ACC KERNELS
    IF (nn_bbl_ldf == 1) THEN
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          za = zab(ji + 1, jj, jp_tem) + zab(ji, jj, jp_tem)
          zb = zab(ji + 1, jj, jp_sal) + zab(ji, jj, jp_sal)
          zgdrho = (za * (zts(ji + 1, jj, jp_tem) - zts(ji, jj, jp_tem)) - zb * (zts(ji + 1, jj, jp_sal) - zts(ji, jj, jp_sal))) * &
&umask(ji, jj, 1)
          zsign = SIGN(0.5, - zgdrho * REAL(mgrhu(ji, jj)))
          ahu_bbl(ji, jj) = (0.5 - zsign) * ahu_bbl_0(ji, jj)
          za = zab(ji, jj + 1, jp_tem) + zab(ji, jj, jp_tem)
          zb = zab(ji, jj + 1, jp_sal) + zab(ji, jj, jp_sal)
          zgdrho = (za * (zts(ji, jj + 1, jp_tem) - zts(ji, jj, jp_tem)) - zb * (zts(ji, jj + 1, jp_sal) - zts(ji, jj, jp_sal))) * &
&vmask(ji, jj, 1)
          zsign = SIGN(0.5, - zgdrho * REAL(mgrhv(ji, jj)))
          ahv_bbl(ji, jj) = (0.5 - zsign) * ahv_bbl_0(ji, jj)
        END DO
      END DO
    END IF
    IF (nn_bbl_adv /= 0) THEN
      SELECT CASE (nn_bbl_adv)
      CASE (1)
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            za = zab(ji + 1, jj, jp_tem) + zab(ji, jj, jp_tem)
            zb = zab(ji + 1, jj, jp_sal) + zab(ji, jj, jp_sal)
            zgdrho = (za * (zts(ji + 1, jj, jp_tem) - zts(ji, jj, jp_tem)) - zb * (zts(ji + 1, jj, jp_sal) - zts(ji, jj, jp_sal))) &
&* umask(ji, jj, 1)
            zsign = SIGN(0.5, - zgdrho * REAL(mgrhu(ji, jj)))
            zsigna = SIGN(0.5, zub(ji, jj) * REAL(mgrhu(ji, jj)))
            utr_bbl(ji, jj) = (0.5 + zsigna) * (0.5 - zsign) * e2u(ji, jj) * e3u_bbl_0(ji, jj) * zub(ji, jj)
            za = zab(ji, jj + 1, jp_tem) + zab(ji, jj, jp_tem)
            zb = zab(ji, jj + 1, jp_sal) + zab(ji, jj, jp_sal)
            zgdrho = (za * (zts(ji, jj + 1, jp_tem) - zts(ji, jj, jp_tem)) - zb * (zts(ji, jj + 1, jp_sal) - zts(ji, jj, jp_sal))) &
&* vmask(ji, jj, 1)
            zsign = SIGN(0.5, - zgdrho * REAL(mgrhv(ji, jj)))
            zsigna = SIGN(0.5, zvb(ji, jj) * REAL(mgrhv(ji, jj)))
            vtr_bbl(ji, jj) = (0.5 + zsigna) * (0.5 - zsign) * e1v(ji, jj) * e3v_bbl_0(ji, jj) * zvb(ji, jj)
          END DO
        END DO
      CASE (2)
        zgbbl = grav * rn_gambbl
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            iid = ji + MAX(0, mgrhu(ji, jj))
            iis = ji + 1 - MAX(0, mgrhu(ji, jj))
            ikud = mbku_d(ji, jj)
            ikus = mbku(ji, jj)
            za = zab(ji + 1, jj, jp_tem) + zab(ji, jj, jp_tem)
            zb = zab(ji + 1, jj, jp_sal) + zab(ji, jj, jp_sal)
            zgdrho = 0.5 * (za * (zts(iid, jj, jp_tem) - zts(iis, jj, jp_tem)) - zb * (zts(iid, jj, jp_sal) - zts(iis, jj, &
&jp_sal))) * umask(ji, jj, 1)
            zgdrho = MAX(0.E0, zgdrho)
            utr_bbl(ji, jj) = e2u(ji, jj) * e3u_bbl_0(ji, jj) * zgbbl * zgdrho * REAL(mgrhu(ji, jj))
            ijd = jj + MAX(0, mgrhv(ji, jj))
            ijs = jj + 1 - MAX(0, mgrhv(ji, jj))
            ikvd = mbkv_d(ji, jj)
            ikvs = mbkv(ji, jj)
            za = zab(ji, jj + 1, jp_tem) + zab(ji, jj, jp_tem)
            zb = zab(ji, jj + 1, jp_sal) + zab(ji, jj, jp_sal)
            zgdrho = 0.5 * (za * (zts(ji, ijd, jp_tem) - zts(ji, ijs, jp_tem)) - zb * (zts(ji, ijd, jp_sal) - zts(ji, ijs, &
&jp_sal))) * vmask(ji, jj, 1)
            zgdrho = MAX(0.E0, zgdrho)
            vtr_bbl(ji, jj) = e1v(ji, jj) * e3v_bbl_0(ji, jj) * zgbbl * zgdrho * REAL(mgrhv(ji, jj))
          END DO
        END DO
      END SELECT
    END IF
    !$ACC END KERNELS
  END SUBROUTINE bbl
  SUBROUTINE tra_bbl_init
    INTEGER :: ji, jj
    INTEGER :: ii0, ii1, ij0, ij1, ios
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zmbku, zmbkv
    NAMELIST /nambbl/ ln_trabbl, nn_bbl_ldf, nn_bbl_adv, rn_ahtbbl, rn_gambbl
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, nambbl, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'nambbl in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, nambbl, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'nambbl in configuration namelist', lwp)
    IF (lwm) WRITE(numond, nambbl)
    l_bbl = .TRUE.
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'tra_bbl_init : bottom boundary layer initialisation'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '       Namelist nambbl : set bbl parameters'
      WRITE(numout, FMT = *) '          bottom boundary layer flag          ln_trabbl  = ', ln_trabbl
    END IF
    IF (.NOT. ln_trabbl) RETURN
    IF (lwp) THEN
      WRITE(numout, FMT = *) '          diffusive bbl (=1)   or not (=0)    nn_bbl_ldf = ', nn_bbl_ldf
      WRITE(numout, FMT = *) '          advective bbl (=1/2) or not (=0)    nn_bbl_adv = ', nn_bbl_adv
      WRITE(numout, FMT = *) '          diffusive bbl coefficient           rn_ahtbbl  = ', rn_ahtbbl, ' m2/s'
      WRITE(numout, FMT = *) '          advective bbl coefficient           rn_gambbl  = ', rn_gambbl, ' s'
    END IF
    IF (tra_bbl_alloc() /= 0) CALL ctl_stop('STOP', 'tra_bbl_init : unable to allocate arrays')
    IF (lwp) THEN
      IF (nn_bbl_adv == 1) WRITE(numout, FMT = *) '       * Advective BBL using upper velocity'
      IF (nn_bbl_adv == 2) WRITE(numout, FMT = *) '       * Advective BBL using velocity = F( delta rho)'
    END IF
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpjm1
      DO ji = 1, jpim1
        mbku_d(ji, jj) = MAX(mbkt(ji + 1, jj), mbkt(ji, jj))
        mbkv_d(ji, jj) = MAX(mbkt(ji, jj + 1), mbkt(ji, jj))
      END DO
    END DO
    zmbku(:, :) = REAL(mbku_d(:, :), wp)
    zmbkv(:, :) = REAL(mbkv_d(:, :), wp)
    !$ACC END KERNELS
    CALL lbc_lnk_multi('trabbl', zmbku, 'U', 1., zmbkv, 'V', 1.)
    !$ACC KERNELS
    mbku_d(:, :) = MAX(INT(zmbku(:, :)), 1)
    mbkv_d(:, :) = MAX(NINT(zmbkv(:, :)), 1)
    mgrhu(:, :) = 0
    mgrhv(:, :) = 0
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpjm1
      DO ji = 1, jpim1
        IF (gdept_0(ji + 1, jj, mbkt(ji + 1, jj)) - gdept_0(ji, jj, mbkt(ji, jj)) /= 0._wp) THEN
          mgrhu(ji, jj) = INT(SIGN(1.E0, gdept_0(ji + 1, jj, mbkt(ji + 1, jj)) - gdept_0(ji, jj, mbkt(ji, jj))))
        END IF
        IF (gdept_0(ji, jj + 1, mbkt(ji, jj + 1)) - gdept_0(ji, jj, mbkt(ji, jj)) /= 0._wp) THEN
          mgrhv(ji, jj) = INT(SIGN(1.E0, gdept_0(ji, jj + 1, mbkt(ji, jj + 1)) - gdept_0(ji, jj, mbkt(ji, jj))))
        END IF
      END DO
    END DO
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpjm1
      DO ji = 1, jpim1
        e3u_bbl_0(ji, jj) = MIN(e3u_0(ji, jj, mbkt(ji + 1, jj)), e3u_0(ji, jj, mbkt(ji, jj)))
        e3v_bbl_0(ji, jj) = MIN(e3v_0(ji, jj, mbkt(ji, jj + 1)), e3v_0(ji, jj, mbkt(ji, jj)))
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk_multi('trabbl', e3u_bbl_0, 'U', 1., e3v_bbl_0, 'V', 1.)
    !$ACC KERNELS
    ahu_bbl_0(:, :) = rn_ahtbbl * e2_e1u(:, :) * e3u_bbl_0(:, :) * umask(:, :, 1)
    ahv_bbl_0(:, :) = rn_ahtbbl * e1_e2v(:, :) * e3v_bbl_0(:, :) * vmask(:, :, 1)
    !$ACC END KERNELS
  END SUBROUTINE tra_bbl_init
END MODULE trabbl