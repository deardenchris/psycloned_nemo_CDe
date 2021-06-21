MODULE icethd_do
  USE dom_oce
  USE phycst
  USE sbc_oce, ONLY: sss_m
  USE sbc_ice, ONLY: utau_ice, vtau_ice
  USE ice1D
  USE ice
  USE icetab
  USE icectl
  USE icethd_ent
  USE icevar
  USE icethd_sal
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_thd_do
  PUBLIC :: ice_thd_do_init
  REAL(KIND = wp) :: rn_hinew
  LOGICAL :: ln_frazil
  REAL(KIND = wp) :: rn_maxfraz
  REAL(KIND = wp) :: rn_vfraz
  REAL(KIND = wp) :: rn_Cfraz
  CONTAINS
  SUBROUTINE ice_thd_do
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jj, jk, jl
    INTEGER :: iter
    REAL(KIND = wp) :: ztmelts, zfrazb, zweight, zde
    REAL(KIND = wp) :: zgamafr, zvfrx, zvgx, ztaux, ztwogp, zf
    REAL(KIND = wp) :: ztenagm, zvfry, zvgy, ztauy, zvrel2, zfp, zsqcd, zhicrit
    REAL(KIND = wp) :: zQm
    REAL(KIND = wp) :: zEi
    REAL(KIND = wp) :: zEw
    REAL(KIND = wp) :: zfmdt
    REAL(KIND = wp) :: zv_newfra
    INTEGER, DIMENSION(jpij) :: jcat
    REAL(KIND = wp), DIMENSION(jpij) :: zswinew
    REAL(KIND = wp), DIMENSION(jpij) :: zv_newice
    REAL(KIND = wp), DIMENSION(jpij) :: za_newice
    REAL(KIND = wp), DIMENSION(jpij) :: zh_newice
    REAL(KIND = wp), DIMENSION(jpij) :: ze_newice
    REAL(KIND = wp), DIMENSION(jpij) :: zs_newice
    REAL(KIND = wp), DIMENSION(jpij) :: zo_newice
    REAL(KIND = wp), DIMENSION(jpij) :: zdv_res
    REAL(KIND = wp), DIMENSION(jpij) :: zda_res
    REAL(KIND = wp), DIMENSION(jpij) :: zv_frazb
    REAL(KIND = wp), DIMENSION(jpij) :: zvrel_1d
    REAL(KIND = wp), DIMENSION(jpij, jpl) :: zv_b
    REAL(KIND = wp), DIMENSION(jpij, jpl) :: za_b
    REAL(KIND = wp), DIMENSION(jpij, nlay_i, jpl) :: ze_i_2d
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zvrel
    REAL(KIND = wp) :: zcai = 1.4E-3_wp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data7
    CALL profile_psy_data0 % PreStart('ice_thd_do', 'r0', 0, 0)
    IF (ln_icediachk) CALL ice_cons_hsm(0, 'icethd_do', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    at_i(:, :) = SUM(a_i, dim = 3)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zvrel(:, :) = 0._wp
    WHERE (qlead(:, :) < 0._wp)
      ht_i_new(:, :) = rn_hinew
    ELSEWHERE
      ht_i_new(:, :) = 0._wp
    END WHERE
    !$ACC END KERNELS
    IF (ln_frazil) THEN
      !$ACC KERNELS
      ht_i_new(:, :) = 0._wp
      zhicrit = 0.04
      ztwogp = 2. * rau0 / (grav * 0.3 * (rau0 - rhoi))
      zsqcd = 1.0 / SQRT(1.3 * zcai)
      zgamafr = 0.03
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('ice_thd_do', 'r1', 0, 0)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          IF (qlead(ji, jj) < 0._wp .AND. tau_icebfr(ji, jj) == 0._wp) THEN
            ztaux = (utau_ice(ji - 1, jj) * umask(ji - 1, jj, 1) + utau_ice(ji, jj) * umask(ji, jj, 1)) * 0.5_wp
            ztauy = (vtau_ice(ji, jj - 1) * vmask(ji, jj - 1, 1) + vtau_ice(ji, jj) * vmask(ji, jj, 1)) * 0.5_wp
            ztenagm = SQRT(SQRT(ztaux * ztaux + ztauy * ztauy))
            rswitch = MAX(0._wp, SIGN(1._wp, ztenagm - epsi10))
            zvfrx = rswitch * zgamafr * zsqcd * ztaux / MAX(ztenagm, epsi10)
            zvfry = rswitch * zgamafr * zsqcd * ztauy / MAX(ztenagm, epsi10)
            zvgx = (u_ice(ji - 1, jj) * umask(ji - 1, jj, 1) + u_ice(ji, jj) * umask(ji, jj, 1)) * 0.5_wp
            zvgy = (v_ice(ji, jj - 1) * vmask(ji, jj - 1, 1) + v_ice(ji, jj) * vmask(ji, jj, 1)) * 0.5_wp
            rswitch = MAX(0._wp, SIGN(1._wp, at_i(ji, jj) - epsi10))
            zvrel2 = MAX((zvfrx - zvgx) * (zvfrx - zvgx) + (zvfry - zvgy) * (zvfry - zvgy), 0.15 * 0.15) * rswitch
            zvrel(ji, jj) = SQRT(zvrel2)
            ht_i_new(ji, jj) = zhicrit + (zhicrit + 0.1) / ((zhicrit + 0.1) * (zhicrit + 0.1) - zhicrit * zhicrit) * ztwogp * zvrel2
            iter = 1
            DO WHILE (iter < 20)
              zf = (ht_i_new(ji, jj) - zhicrit) * (ht_i_new(ji, jj) * ht_i_new(ji, jj) - zhicrit * zhicrit) - ht_i_new(ji, jj) * zhicrit * ztwogp * zvrel2
              zfp = (ht_i_new(ji, jj) - zhicrit) * (3.0 * ht_i_new(ji, jj) + zhicrit) - zhicrit * ztwogp * zvrel2
              ht_i_new(ji, jj) = ht_i_new(ji, jj) - zf / MAX(zfp, epsi20)
              iter = iter + 1
            END DO
          END IF
        END DO
      END DO
      CALL lbc_lnk_multi('icethd_do', zvrel, 'T', 1., ht_i_new, 'T', 1.)
      CALL profile_psy_data1 % PostEnd
    END IF
    !$ACC KERNELS
    npti = 0
    nptidx(:) = 0
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (qlead(ji, jj) < 0._wp .AND. tau_icebfr(ji, jj) == 0._wp) THEN
          npti = npti + 1
          nptidx(npti) = (jj - 1) * jpi + ji
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    IF (npti > 0) THEN
      CALL profile_psy_data2 % PreStart('ice_thd_do', 'r2', 0, 0)
      CALL tab_2d_1d(npti, nptidx(1 : npti), at_i_1d(1 : npti), at_i)
      CALL tab_3d_2d(npti, nptidx(1 : npti), a_i_2d(1 : npti, 1 : jpl), a_i(:, :, :))
      CALL tab_3d_2d(npti, nptidx(1 : npti), v_i_2d(1 : npti, 1 : jpl), v_i(:, :, :))
      CALL tab_3d_2d(npti, nptidx(1 : npti), sv_i_2d(1 : npti, 1 : jpl), sv_i(:, :, :))
      DO jl = 1, jpl
        DO jk = 1, nlay_i
          CALL tab_2d_1d(npti, nptidx(1 : npti), ze_i_2d(1 : npti, jk, jl), e_i(:, :, jk, jl))
        END DO
      END DO
      CALL tab_2d_1d(npti, nptidx(1 : npti), qlead_1d(1 : npti), qlead)
      CALL tab_2d_1d(npti, nptidx(1 : npti), t_bo_1d(1 : npti), t_bo)
      CALL tab_2d_1d(npti, nptidx(1 : npti), sfx_opw_1d(1 : npti), sfx_opw)
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_opw_1d(1 : npti), wfx_opw)
      CALL tab_2d_1d(npti, nptidx(1 : npti), zh_newice(1 : npti), ht_i_new)
      CALL tab_2d_1d(npti, nptidx(1 : npti), zvrel_1d(1 : npti), zvrel)
      CALL tab_2d_1d(npti, nptidx(1 : npti), hfx_thd_1d(1 : npti), hfx_thd)
      CALL tab_2d_1d(npti, nptidx(1 : npti), hfx_opw_1d(1 : npti), hfx_opw)
      CALL tab_2d_1d(npti, nptidx(1 : npti), rn_amax_1d(1 : npti), rn_amax_2d)
      CALL tab_2d_1d(npti, nptidx(1 : npti), sss_1d(1 : npti), sss_m)
      DO jl = 1, jpl
        DO jk = 1, nlay_i
          WHERE (v_i_2d(1 : npti, jl) > 0._wp)
            ze_i_2d(1 : npti, jk, jl) = ze_i_2d(1 : npti, jk, jl) / v_i_2d(1 : npti, jl) * REAL(nlay_i)
          ELSEWHERE
            ze_i_2d(1 : npti, jk, jl) = 0._wp
          END WHERE
        END DO
      END DO
      CALL profile_psy_data2 % PostEnd
      !$ACC KERNELS
      zv_b(1 : npti, :) = v_i_2d(1 : npti, :)
      za_b(1 : npti, :) = a_i_2d(1 : npti, :)
      SELECT CASE (nn_icesal)
      CASE (1)
        zs_newice(1 : npti) = rn_icesal
      CASE (2)
        DO ji = 1, npti
          zs_newice(ji) = MIN(4.606 + 0.91 / zh_newice(ji), rn_simax, 0.5 * sss_1d(ji))
        END DO
      CASE (3)
        zs_newice(1 : npti) = 2.3
      END SELECT
      DO ji = 1, npti
        ztmelts = - rTmlt * zs_newice(ji)
        ze_newice(ji) = rhoi * (rcpi * (ztmelts - (t_bo_1d(ji) - rt0)) + rLfus * (1.0 - ztmelts / MIN(t_bo_1d(ji) - rt0, - epsi10)) - rcp * ztmelts)
      END DO
      zo_newice(1 : npti) = 0._wp
      DO ji = 1, npti
        zEi = - ze_newice(ji) * r1_rhoi
        zEw = rcp * (t_bo_1d(ji) - rt0)
        zdE = zEi - zEw
        zfmdt = - qlead_1d(ji) / zdE
        zv_newice(ji) = - zfmdt * r1_rhoi
        zQm = zfmdt * zEw
        hfx_thd_1d(ji) = hfx_thd_1d(ji) + zfmdt * zEw * r1_rdtice
        hfx_opw_1d(ji) = hfx_opw_1d(ji) - zfmdt * zdE * r1_rdtice
        wfx_opw_1d(ji) = wfx_opw_1d(ji) - zv_newice(ji) * rhoi * r1_rdtice
        sfx_opw_1d(ji) = sfx_opw_1d(ji) - zv_newice(ji) * rhoi * zs_newice(ji) * r1_rdtice
      END DO
      zv_frazb(1 : npti) = 0._wp
      !$ACC END KERNELS
      CALL profile_psy_data3 % PreStart('ice_thd_do', 'r3', 0, 0)
      IF (ln_frazil) THEN
        DO ji = 1, npti
          rswitch = 1._wp - MAX(0._wp, SIGN(1._wp, - at_i_1d(ji)))
          zfrazb = rswitch * (TANH(rn_Cfraz * (zvrel_1d(ji) - rn_vfraz)) + 1.0) * 0.5 * rn_maxfraz
          zv_frazb(ji) = zfrazb * zv_newice(ji)
          zv_newice(ji) = (1.0 - zfrazb) * zv_newice(ji)
        END DO
      END IF
      CALL profile_psy_data3 % PostEnd
      !$ACC KERNELS
      DO ji = 1, npti
        za_newice(ji) = zv_newice(ji) / zh_newice(ji)
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data4 % PreStart('ice_thd_do', 'r4', 0, 0)
      DO ji = 1, npti
        IF (za_newice(ji) > MAX(0._wp, rn_amax_1d(ji) - at_i_1d(ji))) THEN
          zda_res(ji) = za_newice(ji) - MAX(0._wp, rn_amax_1d(ji) - at_i_1d(ji))
          zdv_res(ji) = zda_res(ji) * zh_newice(ji)
          za_newice(ji) = MAX(0._wp, za_newice(ji) - zda_res(ji))
          zv_newice(ji) = MAX(0._wp, zv_newice(ji) - zdv_res(ji))
        ELSE
          zda_res(ji) = 0._wp
          zdv_res(ji) = 0._wp
        END IF
      END DO
      DO jl = 1, jpl
        DO ji = 1, npti
          IF (zh_newice(ji) > hi_max(jl - 1) .AND. zh_newice(ji) <= hi_max(jl)) THEN
            a_i_2d(ji, jl) = a_i_2d(ji, jl) + za_newice(ji)
            v_i_2d(ji, jl) = v_i_2d(ji, jl) + zv_newice(ji)
            jcat(ji) = jl
          END IF
        END DO
      END DO
      at_i_1d(1 : npti) = SUM(a_i_2d(1 : npti, :), dim = 2)
      CALL profile_psy_data4 % PostEnd
      !$ACC KERNELS
      DO ji = 1, npti
        jl = jcat(ji)
        zswinew(ji) = MAX(0._wp, SIGN(1._wp, - za_b(ji, jl)))
      END DO
      DO jk = 1, nlay_i
        DO ji = 1, npti
          jl = jcat(ji)
          rswitch = MAX(0._wp, SIGN(1._wp, v_i_2d(ji, jl) - epsi20))
          ze_i_2d(ji, jk, jl) = zswinew(ji) * ze_newice(ji) + (1.0 - zswinew(ji)) * (ze_newice(ji) * zv_newice(ji) + ze_i_2d(ji, jk, jl) * zv_b(ji, jl)) * rswitch / MAX(v_i_2d(ji, jl), epsi20)
        END DO
      END DO
      !$ACC END KERNELS
      DO jl = 1, jpl
        !$ACC KERNELS
        h_i_old(1 : npti, 0 : nlay_i + 1) = 0._wp
        eh_i_old(1 : npti, 0 : nlay_i + 1) = 0._wp
        DO jk = 1, nlay_i
          DO ji = 1, npti
            h_i_old(ji, jk) = v_i_2d(ji, jl) * r1_nlay_i
            eh_i_old(ji, jk) = ze_i_2d(ji, jk, jl) * h_i_old(ji, jk)
          END DO
        END DO
        DO ji = 1, npti
          rswitch = MAX(0._wp, SIGN(1._wp, at_i_1d(ji) - epsi20))
          zv_newfra = rswitch * (zdv_res(ji) + zv_frazb(ji)) * a_i_2d(ji, jl) / MAX(at_i_1d(ji), epsi20)
          a_i_2d(ji, jl) = rswitch * a_i_2d(ji, jl)
          v_i_2d(ji, jl) = v_i_2d(ji, jl) + zv_newfra
          h_i_old(ji, nlay_i + 1) = zv_newfra
          eh_i_old(ji, nlay_i + 1) = ze_newice(ji) * zv_newfra
        END DO
        !$ACC END KERNELS
        CALL profile_psy_data5 % PreStart('ice_thd_do', 'r5', 0, 0)
        CALL ice_thd_ent(ze_i_2d(1 : npti, :, jl))
        CALL profile_psy_data5 % PostEnd
      END DO
      !$ACC KERNELS
      DO jl = 1, jpl
        DO ji = 1, npti
          sv_i_2d(ji, jl) = sv_i_2d(ji, jl) + zs_newice(ji) * (v_i_2d(ji, jl) - zv_b(ji, jl))
        END DO
      END DO
      !$ACC END KERNELS
      DO jl = 1, jpl
        !$ACC KERNELS
        DO jk = 1, nlay_i
          ze_i_2d(1 : npti, jk, jl) = ze_i_2d(1 : npti, jk, jl) * v_i_2d(1 : npti, jl) * r1_nlay_i
        END DO
        !$ACC END KERNELS
      END DO
      CALL profile_psy_data6 % PreStart('ice_thd_do', 'r6', 0, 0)
      CALL tab_2d_3d(npti, nptidx(1 : npti), a_i_2d(1 : npti, 1 : jpl), a_i(:, :, :))
      CALL tab_2d_3d(npti, nptidx(1 : npti), v_i_2d(1 : npti, 1 : jpl), v_i(:, :, :))
      CALL tab_2d_3d(npti, nptidx(1 : npti), sv_i_2d(1 : npti, 1 : jpl), sv_i(:, :, :))
      DO jl = 1, jpl
        DO jk = 1, nlay_i
          CALL tab_1d_2d(npti, nptidx(1 : npti), ze_i_2d(1 : npti, jk, jl), e_i(:, :, jk, jl))
        END DO
      END DO
      CALL tab_1d_2d(npti, nptidx(1 : npti), sfx_opw_1d(1 : npti), sfx_opw)
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_opw_1d(1 : npti), wfx_opw)
      CALL tab_1d_2d(npti, nptidx(1 : npti), hfx_thd_1d(1 : npti), hfx_thd)
      CALL tab_1d_2d(npti, nptidx(1 : npti), hfx_opw_1d(1 : npti), hfx_opw)
      CALL profile_psy_data6 % PostEnd
    END IF
    CALL profile_psy_data7 % PreStart('ice_thd_do', 'r7', 0, 0)
    IF (ln_icediachk) CALL ice_cons_hsm(1, 'icethd_do', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    CALL profile_psy_data7 % PostEnd
  END SUBROUTINE ice_thd_do
  SUBROUTINE ice_thd_do_init
    INTEGER :: ios
    NAMELIST /namthd_do/ rn_hinew, ln_frazil, rn_maxfraz, rn_vfraz, rn_Cfraz
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, namthd_do, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namthd_do in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, namthd_do, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namthd_do in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, namthd_do)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ice_thd_do_init: Ice growth in open water'
      WRITE(numout, *) '~~~~~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist namthd_do:'
      WRITE(numout, *) '      ice thickness for lateral accretion                       rn_hinew   = ', rn_hinew
      WRITE(numout, *) '      Frazil ice thickness as a function of wind or not         ln_frazil  = ', ln_frazil
      WRITE(numout, *) '      Maximum proportion of frazil ice collecting at bottom     rn_maxfraz = ', rn_maxfraz
      WRITE(numout, *) '      Threshold relative drift speed for collection of frazil   rn_vfraz   = ', rn_vfraz
      WRITE(numout, *) '      Squeezing coefficient for collection of frazil            rn_Cfraz   = ', rn_Cfraz
    END IF
    IF (rn_hinew < rn_himin) CALL ctl_stop('ice_thd_do_init : rn_hinew should be >= rn_himin')
  END SUBROUTINE ice_thd_do_init
END MODULE icethd_do