MODULE icedyn_rhg_evp
  USE phycst
  USE dom_oce
  USE sbc_oce, ONLY: ln_ice_embd, nn_fsbc, ssh_m
  USE sbc_ice, ONLY: utau_ice, vtau_ice, snwice_mass, snwice_mass_b
  USE ice
  USE icevar
  USE icedyn_rdgrft
  USE bdy_oce, ONLY: ln_bdy
  USE bdyice
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  USE lbclnk
  USE prtctl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_dyn_rhg_evp
  PUBLIC :: rhg_evp_rst
  CONTAINS
  SUBROUTINE ice_dyn_rhg_evp(kt, pstress1_i, pstress2_i, pstress12_i, pshear_i, pdivu_i, pdelta_i)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: pstress1_i, pstress2_i, pstress12_i
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: pshear_i, pdivu_i, pdelta_i
    LOGICAL, PARAMETER :: ll_bdy_substep = .TRUE.
    INTEGER :: ji, jj
    INTEGER :: jter
    REAL(KIND = wp) :: zrhoco
    REAL(KIND = wp) :: zdtevp, z1_dtevp
    REAL(KIND = wp) :: ecc2, z1_ecc2
    REAL(KIND = wp) :: zalph1, z1_alph1, zalph2, z1_alph2
    REAL(KIND = wp) :: zm1, zm2, zm3, zmassU, zmassV, zvU, zvV
    REAL(KIND = wp) :: zdelta, zp_delf, zds2, zdt, zdt2, zdiv, zdiv2
    REAL(KIND = wp) :: zTauO, zTauB, zTauE, zvel
    REAL(KIND = wp) :: zkt
    REAL(KIND = wp) :: zvCr
    REAL(KIND = wp) :: zresm
    REAL(KIND = wp) :: zintb, zintn
    REAL(KIND = wp) :: zfac_x, zfac_y
    REAL(KIND = wp) :: zshear, zdum1, zdum2
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z1_e1t0, z1_e2t0
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zp_delt
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zbeta
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdt_m
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zaU, zaV
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zmU_t, zmV_t
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zmf
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zTauU_ia, ztauV_ia
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zTauU_ib, ztauV_ib
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zspgU, zspgV
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: v_oceU, u_oceV, v_iceU, u_iceV
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zfU, zfV
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zds
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zs1, zs2, zs12
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zsshdyn
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zCorx, zCory
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztaux_oi, ztauy_oi
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zswitchU, zswitchV
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zmaskU, zmaskV
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zfmask, zwf
    REAL(KIND = wp), PARAMETER :: zepsi = 1.0E-20_wp
    REAL(KIND = wp), PARAMETER :: zmmin = 1._wp
    REAL(KIND = wp), PARAMETER :: zamin = 0.001_wp
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zswi
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zsig1, zsig2, zsig3
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_sig1
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_sig2
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_dssh_dx
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_dssh_dy
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_corstrx
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_corstry
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_intstrx
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_intstry
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_utau_oi
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_vtau_oi
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_xmtrp_ice
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_ymtrp_ice
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_xmtrp_snw
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_ymtrp_snw
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_xatrp
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_yatrp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data7
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data8
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data9
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data10
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data11
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data12
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data13
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data14
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data15
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data16
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data17
    IF (kt == nit000 .AND. lwp) WRITE(numout, *) '-- ice_dyn_rhg_evp: EVP sea-ice rheology'
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 1, jpjm1
      DO ji = 1, jpim1
        zfmask(ji, jj) = tmask(ji, jj, 1) * tmask(ji + 1, jj, 1) * tmask(ji, jj + 1, 1) * tmask(ji + 1, jj + 1, 1)
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('ice_dyn_rhg_evp', 'r0', 0, 0)
    CALL lbc_lnk('icedyn_rhg_evp', zfmask, 'F', 1._wp)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zwf(:, :) = zfmask(:, :)
    !$ACC loop independent collapse(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        IF (zfmask(ji, jj) == 0._wp) THEN
          zfmask(ji, jj) = rn_ishlat * MIN(1._wp, MAX(zwf(ji + 1, jj), zwf(ji, jj + 1), zwf(ji - 1, jj), zwf(ji, jj - 1)))
        END IF
      END DO
    END DO
    !$ACC loop independent
    DO jj = 2, jpjm1
      IF (zfmask(1, jj) == 0._wp) THEN
        zfmask(1, jj) = rn_ishlat * MIN(1._wp, MAX(zwf(2, jj), zwf(1, jj + 1), zwf(1, jj - 1)))
      END IF
      IF (zfmask(jpi, jj) == 0._wp) THEN
        zfmask(jpi, jj) = rn_ishlat * MIN(1._wp, MAX(zwf(jpi, jj + 1), zwf(jpim1, jj), zwf(jpi, jj - 1)))
      END IF
    END DO
    !$ACC loop independent
    DO ji = 2, jpim1
      IF (zfmask(ji, 1) == 0._wp) THEN
        zfmask(ji, 1) = rn_ishlat * MIN(1._wp, MAX(zwf(ji + 1, 1), zwf(ji, 2), zwf(ji - 1, 1)))
      END IF
      IF (zfmask(ji, jpj) == 0._wp) THEN
        zfmask(ji, jpj) = rn_ishlat * MIN(1._wp, MAX(zwf(ji + 1, jpj), zwf(ji - 1, jpj), zwf(ji, jpjm1)))
      END IF
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('ice_dyn_rhg_evp', 'r1', 0, 0)
    CALL lbc_lnk('icedyn_rhg_evp', zfmask, 'F', 1._wp)
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    zrhoco = rau0 * rn_cio
    ecc2 = rn_ecc * rn_ecc
    z1_ecc2 = 1._wp / ecc2
    zdtevp = rdt_ice / REAL(nn_nevp)
    z1_dtevp = 1._wp / zdtevp
    IF (.NOT. ln_aEVP) THEN
      zalph1 = (2._wp * rn_relast * rdt_ice) * z1_dtevp
      zalph2 = zalph1 * z1_ecc2
      z1_alph1 = 1._wp / (zalph1 + 1._wp)
      z1_alph2 = 1._wp / (zalph2 + 1._wp)
    END IF
    zs1(:, :) = pstress1_i(:, :)
    zs2(:, :) = pstress2_i(:, :)
    zs12(:, :) = pstress12_i(:, :)
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('ice_dyn_rhg_evp', 'r2', 0, 0)
    CALL ice_strength
    CALL profile_psy_data2 % PostEnd
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        z1_e1t0(ji, jj) = 1._wp / (e1t(ji + 1, jj) + e1t(ji, jj))
        z1_e2t0(ji, jj) = 1._wp / (e2t(ji, jj + 1) + e2t(ji, jj))
      END DO
    END DO
    IF (ln_landfast_l16 .OR. ln_landfast_home) THEN
      zkt = rn_tensile
    ELSE
      zkt = 0._wp
    END IF
    !$ACC END KERNELS
    CALL profile_psy_data3 % PreStart('ice_dyn_rhg_evp', 'r3', 0, 0)
    zsshdyn(:, :) = ice_var_sshdyn(ssh_m, snwice_mass, snwice_mass_b)
    CALL profile_psy_data3 % PostEnd
    !$ACC KERNELS
    !$ACC loop independent collapse(2) private(zmassU, zmassV)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zaU(ji, jj) = 0.5_wp * (at_i(ji, jj) * e1e2t(ji, jj) + at_i(ji + 1, jj) * e1e2t(ji + 1, jj)) * r1_e1e2u(ji, jj) * umask(ji, jj, 1)
        zaV(ji, jj) = 0.5_wp * (at_i(ji, jj) * e1e2t(ji, jj) + at_i(ji, jj + 1) * e1e2t(ji, jj + 1)) * r1_e1e2v(ji, jj) * vmask(ji, jj, 1)
        zm1 = (rhos * vt_s(ji, jj) + rhoi * vt_i(ji, jj))
        zm2 = (rhos * vt_s(ji + 1, jj) + rhoi * vt_i(ji + 1, jj))
        zm3 = (rhos * vt_s(ji, jj + 1) + rhoi * vt_i(ji, jj + 1))
        zmassU = 0.5_wp * (zm1 * e1e2t(ji, jj) + zm2 * e1e2t(ji + 1, jj)) * r1_e1e2u(ji, jj) * umask(ji, jj, 1)
        zmassV = 0.5_wp * (zm1 * e1e2t(ji, jj) + zm3 * e1e2t(ji, jj + 1)) * r1_e1e2v(ji, jj) * vmask(ji, jj, 1)
        v_oceU(ji, jj) = 0.5_wp * ((v_oce(ji, jj) + v_oce(ji, jj - 1)) * e1t(ji + 1, jj) + (v_oce(ji + 1, jj) + v_oce(ji + 1, jj - 1)) * e1t(ji, jj)) * z1_e1t0(ji, jj) * umask(ji, jj, 1)
        u_oceV(ji, jj) = 0.5_wp * ((u_oce(ji, jj) + u_oce(ji - 1, jj)) * e2t(ji, jj + 1) + (u_oce(ji, jj + 1) + u_oce(ji - 1, jj + 1)) * e2t(ji, jj)) * z1_e2t0(ji, jj) * vmask(ji, jj, 1)
        zmf(ji, jj) = zm1 * ff_t(ji, jj)
        zdt_m(ji, jj) = zdtevp / MAX(zm1, zmmin)
        zmU_t(ji, jj) = zmassU * z1_dtevp
        zmV_t(ji, jj) = zmassV * z1_dtevp
        zTauU_ia(ji, jj) = zaU(ji, jj) * utau_ice(ji, jj)
        zTauV_ia(ji, jj) = zaV(ji, jj) * vtau_ice(ji, jj)
        zspgU(ji, jj) = - zmassU * grav * (zsshdyn(ji + 1, jj) - zsshdyn(ji, jj)) * r1_e1u(ji, jj)
        zspgV(ji, jj) = - zmassV * grav * (zsshdyn(ji, jj + 1) - zsshdyn(ji, jj)) * r1_e2v(ji, jj)
        zmaskU(ji, jj) = 1._wp - MAX(0._wp, SIGN(1._wp, - zmassU))
        zmaskV(ji, jj) = 1._wp - MAX(0._wp, SIGN(1._wp, - zmassV))
        IF (zmassu <= zmmin .AND. zau(ji, jj) <= zamin) THEN
          zswitchu(ji, jj) = 0._wp
        ELSE
          zswitchu(ji, jj) = 1._wp
        END IF
        IF (zmassv <= zmmin .AND. zav(ji, jj) <= zamin) THEN
          zswitchv(ji, jj) = 0._wp
        ELSE
          zswitchv(ji, jj) = 1._wp
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data4 % PreStart('ice_dyn_rhg_evp', 'r4', 0, 0)
    CALL lbc_lnk_multi('icedyn_rhg_evp', zmf, 'T', 1., zdt_m, 'T', 1.)
    CALL profile_psy_data4 % PostEnd
    IF (ln_landfast_L16) THEN
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zvU = 0.5_wp * (vt_i(ji, jj) * e1e2t(ji, jj) + vt_i(ji + 1, jj) * e1e2t(ji + 1, jj)) * r1_e1e2u(ji, jj) * umask(ji, jj, 1)
          zvV = 0.5_wp * (vt_i(ji, jj) * e1e2t(ji, jj) + vt_i(ji, jj + 1) * e1e2t(ji, jj + 1)) * r1_e1e2v(ji, jj) * vmask(ji, jj, 1)
          zvCr = zaU(ji, jj) * rn_depfra * hu_n(ji, jj)
          zTauU_ib(ji, jj) = rn_icebfr * MAX(0._wp, zvU - zvCr) * EXP(- rn_crhg * (1._wp - zaU(ji, jj)))
          zvCr = zaV(ji, jj) * rn_depfra * hv_n(ji, jj)
          zTauV_ib(ji, jj) = rn_icebfr * MAX(0._wp, zvV - zvCr) * EXP(- rn_crhg * (1._wp - zaV(ji, jj)))
          zvCr = at_i(ji, jj) * rn_depfra * ht_n(ji, jj)
          tau_icebfr(ji, jj) = rn_icebfr * MAX(0._wp, vt_i(ji, jj) - zvCr) * EXP(- rn_crhg * (1._wp - at_i(ji, jj)))
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data5 % PreStart('ice_dyn_rhg_evp', 'r5', 0, 0)
      CALL lbc_lnk('icedyn_rhg_evp', tau_icebfr(:, :), 'T', 1.)
      CALL profile_psy_data5 % PostEnd
    ELSE IF (ln_landfast_home) THEN
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zTauU_ib(ji, jj) = tau_icebfr(ji, jj)
          zTauV_ib(ji, jj) = tau_icebfr(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zTauU_ib(ji, jj) = 0._wp
          zTauV_ib(ji, jj) = 0._wp
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data6 % PreStart('ice_dyn_rhg_evp', 'r6', 0, 0)
    IF (iom_use('tau_icebfr')) CALL iom_put('tau_icebfr', tau_icebfr(:, :))
    CALL profile_psy_data6 % PostEnd
    DO jter = 1, nn_nevp
      !$ACC KERNELS
      l_full_nf_update = jter == nn_nevp
      !$ACC loop independent collapse(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zds(ji, jj) = ((u_ice(ji, jj + 1) * r1_e1u(ji, jj + 1) - u_ice(ji, jj) * r1_e1u(ji, jj)) * e1f(ji, jj) * e1f(ji, jj) + (v_ice(ji + 1, jj) * r1_e2v(ji + 1, jj) - v_ice(ji, jj) * r1_e2v(ji, jj)) * e2f(ji, jj) * e2f(ji, jj)) * r1_e1e2f(ji, jj) * zfmask(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data7 % PreStart('ice_dyn_rhg_evp', 'r7', 0, 0)
      CALL lbc_lnk('icedyn_rhg_evp', zds, 'F', 1.)
      CALL profile_psy_data7 % PostEnd
      !$ACC KERNELS
      ! CDe - private clause needed to parallelise the nested loop
      !$ACC loop independent collapse(2) private(zalph2, z1_alph2)
      DO jj = 2, jpj
        DO ji = 2, jpi
          zds2 = (zds(ji, jj) * zds(ji, jj) * e1e2f(ji, jj) + zds(ji - 1, jj) * zds(ji - 1, jj) * e1e2f(ji - 1, jj) + zds(ji, jj - 1) * zds(ji, jj - 1) * e1e2f(ji, jj - 1) + zds(ji - 1, jj - 1) * zds(ji - 1, jj - 1) * e1e2f(ji - 1, jj - 1)) * 0.25_wp * r1_e1e2t(ji, jj)
          zdiv = (e2u(ji, jj) * u_ice(ji, jj) - e2u(ji - 1, jj) * u_ice(ji - 1, jj) + e1v(ji, jj) * v_ice(ji, jj) - e1v(ji, jj - 1) * v_ice(ji, jj - 1)) * r1_e1e2t(ji, jj)
          zdiv2 = zdiv * zdiv
          zdt = ((u_ice(ji, jj) * r1_e2u(ji, jj) - u_ice(ji - 1, jj) * r1_e2u(ji - 1, jj)) * e2t(ji, jj) * e2t(ji, jj) - (v_ice(ji, jj) * r1_e1v(ji, jj) - v_ice(ji, jj - 1) * r1_e1v(ji, jj - 1)) * e1t(ji, jj) * e1t(ji, jj)) * r1_e1e2t(ji, jj)
          zdt2 = zdt * zdt
          zdelta = SQRT(zdiv2 + (zdt2 + zds2) * z1_ecc2)
          zp_delt(ji, jj) = strength(ji, jj) / (zdelta + rn_creepl)
          IF (ln_aEVP) THEN
            zalph1 = MAX(50._wp, rpi * SQRT(0.5_wp * zp_delt(ji, jj) * r1_e1e2t(ji, jj) * zdt_m(ji, jj)))
            z1_alph1 = 1._wp / (zalph1 + 1._wp)
            zalph2 = zalph1
            z1_alph2 = z1_alph1
          END IF
          zs1(ji, jj) = (zs1(ji, jj) * zalph1 + zp_delt(ji, jj) * (zdiv * (1._wp + zkt) - zdelta * (1._wp - zkt))) * z1_alph1
          zs2(ji, jj) = (zs2(ji, jj) * zalph2 + zp_delt(ji, jj) * (zdt * z1_ecc2 * (1._wp + zkt))) * z1_alph2
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data8 % PreStart('ice_dyn_rhg_evp', 'r8', 0, 0)
      CALL lbc_lnk('icedyn_rhg_evp', zp_delt, 'T', 1.)
      CALL profile_psy_data8 % PostEnd
      !$ACC KERNELS
      !$ACC loop independent collapse(2) private(zalph2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          IF (ln_aEVP) THEN
            zalph2 = MAX(50._wp, rpi * SQRT(0.5_wp * zp_delt(ji, jj) * r1_e1e2t(ji, jj) * zdt_m(ji, jj)))
            z1_alph2 = 1._wp / (zalph2 + 1._wp)
            zbeta(ji, jj) = zalph2
          END IF
          zp_delf = 0.25_wp * (zp_delt(ji, jj) + zp_delt(ji + 1, jj) + zp_delt(ji, jj + 1) + zp_delt(ji + 1, jj + 1))
          zs12(ji, jj) = (zs12(ji, jj) * zalph2 + zp_delf * (zds(ji, jj) * z1_ecc2 * (1._wp + zkt)) * 0.5_wp) * z1_alph2
        END DO
      END DO
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zfU(ji, jj) = 0.5_wp * ((zs1(ji + 1, jj) - zs1(ji, jj)) * e2u(ji, jj) + (zs2(ji + 1, jj) * e2t(ji + 1, jj) * e2t(ji + 1, jj) - zs2(ji, jj) * e2t(ji, jj) * e2t(ji, jj)) * r1_e2u(ji, jj) + (zs12(ji, jj) * e1f(ji, jj) * e1f(ji, jj) - zs12(ji, jj - 1) * e1f(ji, jj - 1) * e1f(ji, jj - 1)) * 2._wp * r1_e1u(ji, jj)) * r1_e1e2u(ji, jj)
          zfV(ji, jj) = 0.5_wp * ((zs1(ji, jj + 1) - zs1(ji, jj)) * e1v(ji, jj) - (zs2(ji, jj + 1) * e1t(ji, jj + 1) * e1t(ji, jj + 1) - zs2(ji, jj) * e1t(ji, jj) * e1t(ji, jj)) * r1_e1v(ji, jj) + (zs12(ji, jj) * e2f(ji, jj) * e2f(ji, jj) - zs12(ji - 1, jj) * e2f(ji - 1, jj) * e2f(ji - 1, jj)) * 2._wp * r1_e2v(ji, jj)) * r1_e1e2v(ji, jj)
          u_iceV(ji, jj) = 0.5_wp * ((u_ice(ji, jj) + u_ice(ji - 1, jj)) * e2t(ji, jj + 1) + (u_ice(ji, jj + 1) + u_ice(ji - 1, jj + 1)) * e2t(ji, jj)) * z1_e2t0(ji, jj) * vmask(ji, jj, 1)
          v_iceU(ji, jj) = 0.5_wp * ((v_ice(ji, jj) + v_ice(ji, jj - 1)) * e1t(ji + 1, jj) + (v_ice(ji + 1, jj) + v_ice(ji + 1, jj - 1)) * e1t(ji, jj)) * z1_e1t0(ji, jj) * umask(ji, jj, 1)
        END DO
      END DO
      !$ACC END KERNELS
      IF (MOD(jter, 2) == 0) THEN
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zTauO = zaV(ji, jj) * zrhoco * SQRT((v_ice(ji, jj) - v_oce(ji, jj)) * (v_ice(ji, jj) - v_oce(ji, jj)) + (u_iceV(ji, jj) - u_oceV(ji, jj)) * (u_iceV(ji, jj) - u_oceV(ji, jj)))
            ztauy_oi(ji, jj) = zTauO * (v_oce(ji, jj) - v_ice(ji, jj))
            zvel = 5.E-05_wp + SQRT(v_ice(ji, jj) * v_ice(ji, jj) + u_iceV(ji, jj) * u_iceV(ji, jj))
            zTauB = - zTauV_ib(ji, jj) / zvel
            zCory(ji, jj) = - 0.25_wp * r1_e2v(ji, jj) * (zmf(ji, jj) * (e2u(ji, jj) * u_ice(ji, jj) + e2u(ji - 1, jj) * u_ice(ji - 1, jj)) + zmf(ji, jj + 1) * (e2u(ji, jj + 1) * u_ice(ji, jj + 1) + e2u(ji - 1, jj + 1) * u_ice(ji - 1, jj + 1)))
            zTauE = zfV(ji, jj) + zTauV_ia(ji, jj) + zCory(ji, jj) + zspgV(ji, jj) + ztauy_oi(ji, jj)
            rswitch = 1._wp - MIN(1._wp, ABS(SIGN(1._wp, ztauE - zTauV_ib(ji, jj)) - SIGN(1._wp, zTauE)))
            IF (ln_aEVP) THEN
              v_ice(ji, jj) = ((rswitch * (zmV_t(ji, jj) * (zbeta(ji, jj) * v_ice(ji, jj) + v_ice_b(ji, jj)) + zTauE + zTauO * v_ice(ji, jj)) / MAX(zepsi, zmV_t(ji, jj) * (zbeta(ji, jj) + 1._wp) + zTauO - zTauB) + (1._wp - rswitch) * v_ice(ji, jj) * MAX(0._wp, 1._wp - zdtevp * rn_lfrelax)) * zswitchV(ji, jj) + v_oce(ji, jj) * 0.01_wp * (1._wp - zswitchV(ji, jj))) * zmaskV(ji, jj)
            ELSE
              v_ice(ji, jj) = ((rswitch * (zmV_t(ji, jj) * v_ice(ji, jj) + zTauE + zTauO * v_ice(ji, jj)) / MAX(zepsi, zmV_t(ji, jj) + zTauO - zTauB) + (1._wp - rswitch) * v_ice(ji, jj) * MAX(0._wp, 1._wp - zdtevp * rn_lfrelax)) * zswitchV(ji, jj) + v_oce(ji, jj) * 0.01_wp * (1._wp - zswitchV(ji, jj))) * zmaskV(ji, jj)
            END IF
          END DO
        END DO
        !$ACC END KERNELS
        CALL profile_psy_data9 % PreStart('ice_dyn_rhg_evp', 'r9', 0, 0)
        CALL lbc_lnk('icedyn_rhg_evp', v_ice, 'V', - 1.)
        IF (ln_bdy .AND. ll_bdy_substep) CALL bdy_ice_dyn('V')
        CALL profile_psy_data9 % PostEnd
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zTauO = zaU(ji, jj) * zrhoco * SQRT((u_ice(ji, jj) - u_oce(ji, jj)) * (u_ice(ji, jj) - u_oce(ji, jj)) + (v_iceU(ji, jj) - v_oceU(ji, jj)) * (v_iceU(ji, jj) - v_oceU(ji, jj)))
            ztaux_oi(ji, jj) = zTauO * (u_oce(ji, jj) - u_ice(ji, jj))
            zvel = 5.E-05_wp + SQRT(v_iceU(ji, jj) * v_iceU(ji, jj) + u_ice(ji, jj) * u_ice(ji, jj))
            zTauB = - zTauU_ib(ji, jj) / zvel
            zCorx(ji, jj) = 0.25_wp * r1_e1u(ji, jj) * (zmf(ji, jj) * (e1v(ji, jj) * v_ice(ji, jj) + e1v(ji, jj - 1) * v_ice(ji, jj - 1)) + zmf(ji + 1, jj) * (e1v(ji + 1, jj) * v_ice(ji + 1, jj) + e1v(ji + 1, jj - 1) * v_ice(ji + 1, jj - 1)))
            zTauE = zfU(ji, jj) + zTauU_ia(ji, jj) + zCorx(ji, jj) + zspgU(ji, jj) + ztaux_oi(ji, jj)
            rswitch = 1._wp - MIN(1._wp, ABS(SIGN(1._wp, ztauE - zTauU_ib(ji, jj)) - SIGN(1._wp, zTauE)))
            IF (ln_aEVP) THEN
              u_ice(ji, jj) = ((rswitch * (zmU_t(ji, jj) * (zbeta(ji, jj) * u_ice(ji, jj) + u_ice_b(ji, jj)) + zTauE + zTauO * u_ice(ji, jj)) / MAX(zepsi, zmU_t(ji, jj) * (zbeta(ji, jj) + 1._wp) + zTauO - zTauB) + (1._wp - rswitch) * u_ice(ji, jj) * MAX(0._wp, 1._wp - zdtevp * rn_lfrelax)) * zswitchU(ji, jj) + u_oce(ji, jj) * 0.01_wp * (1._wp - zswitchU(ji, jj))) * zmaskU(ji, jj)
            ELSE
              u_ice(ji, jj) = ((rswitch * (zmU_t(ji, jj) * u_ice(ji, jj) + zTauE + zTauO * u_ice(ji, jj)) / MAX(zepsi, zmU_t(ji, jj) + zTauO - zTauB) + (1._wp - rswitch) * u_ice(ji, jj) * MAX(0._wp, 1._wp - zdtevp * rn_lfrelax)) * zswitchU(ji, jj) + u_oce(ji, jj) * 0.01_wp * (1._wp - zswitchU(ji, jj))) * zmaskU(ji, jj)
            END IF
          END DO
        END DO
        !$ACC END KERNELS
        CALL profile_psy_data10 % PreStart('ice_dyn_rhg_evp', 'r10', 0, 0)
        CALL lbc_lnk('icedyn_rhg_evp', u_ice, 'U', - 1.)
        IF (ln_bdy .AND. ll_bdy_substep) CALL bdy_ice_dyn('U')
        CALL profile_psy_data10 % PostEnd
      ELSE
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zTauO = zaU(ji, jj) * zrhoco * SQRT((u_ice(ji, jj) - u_oce(ji, jj)) * (u_ice(ji, jj) - u_oce(ji, jj)) + (v_iceU(ji, jj) - v_oceU(ji, jj)) * (v_iceU(ji, jj) - v_oceU(ji, jj)))
            ztaux_oi(ji, jj) = zTauO * (u_oce(ji, jj) - u_ice(ji, jj))
            zvel = 5.E-05_wp + SQRT(v_iceU(ji, jj) * v_iceU(ji, jj) + u_ice(ji, jj) * u_ice(ji, jj))
            zTauB = - zTauU_ib(ji, jj) / zvel
            zCorx(ji, jj) = 0.25_wp * r1_e1u(ji, jj) * (zmf(ji, jj) * (e1v(ji, jj) * v_ice(ji, jj) + e1v(ji, jj - 1) * v_ice(ji, jj - 1)) + zmf(ji + 1, jj) * (e1v(ji + 1, jj) * v_ice(ji + 1, jj) + e1v(ji + 1, jj - 1) * v_ice(ji + 1, jj - 1)))
            zTauE = zfU(ji, jj) + zTauU_ia(ji, jj) + zCorx(ji, jj) + zspgU(ji, jj) + ztaux_oi(ji, jj)
            rswitch = 1._wp - MIN(1._wp, ABS(SIGN(1._wp, ztauE - zTauU_ib(ji, jj)) - SIGN(1._wp, zTauE)))
            IF (ln_aEVP) THEN
              u_ice(ji, jj) = ((rswitch * (zmU_t(ji, jj) * (zbeta(ji, jj) * u_ice(ji, jj) + u_ice_b(ji, jj)) + zTauE + zTauO * u_ice(ji, jj)) / MAX(zepsi, zmU_t(ji, jj) * (zbeta(ji, jj) + 1._wp) + zTauO - zTauB) + (1._wp - rswitch) * u_ice(ji, jj) * MAX(0._wp, 1._wp - zdtevp * rn_lfrelax)) * zswitchU(ji, jj) + u_oce(ji, jj) * 0.01_wp * (1._wp - zswitchU(ji, jj))) * zmaskU(ji, jj)
            ELSE
              u_ice(ji, jj) = ((rswitch * (zmU_t(ji, jj) * u_ice(ji, jj) + zTauE + zTauO * u_ice(ji, jj)) / MAX(zepsi, zmU_t(ji, jj) + zTauO - zTauB) + (1._wp - rswitch) * u_ice(ji, jj) * MAX(0._wp, 1._wp - zdtevp * rn_lfrelax)) * zswitchU(ji, jj) + u_oce(ji, jj) * 0.01_wp * (1._wp - zswitchU(ji, jj))) * zmaskU(ji, jj)
            END IF
          END DO
        END DO
        !$ACC END KERNELS
        CALL profile_psy_data11 % PreStart('ice_dyn_rhg_evp', 'r11', 0, 0)
        CALL lbc_lnk('icedyn_rhg_evp', u_ice, 'U', - 1.)
        IF (ln_bdy .AND. ll_bdy_substep) CALL bdy_ice_dyn('U')
        CALL profile_psy_data11 % PostEnd
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zTauO = zaV(ji, jj) * zrhoco * SQRT((v_ice(ji, jj) - v_oce(ji, jj)) * (v_ice(ji, jj) - v_oce(ji, jj)) + (u_iceV(ji, jj) - u_oceV(ji, jj)) * (u_iceV(ji, jj) - u_oceV(ji, jj)))
            ztauy_oi(ji, jj) = zTauO * (v_oce(ji, jj) - v_ice(ji, jj))
            zvel = 5.E-05_wp + SQRT(v_ice(ji, jj) * v_ice(ji, jj) + u_iceV(ji, jj) * u_iceV(ji, jj))
            zTauB = - zTauV_ib(ji, jj) / zvel
            zCory(ji, jj) = - 0.25_wp * r1_e2v(ji, jj) * (zmf(ji, jj) * (e2u(ji, jj) * u_ice(ji, jj) + e2u(ji - 1, jj) * u_ice(ji - 1, jj)) + zmf(ji, jj + 1) * (e2u(ji, jj + 1) * u_ice(ji, jj + 1) + e2u(ji - 1, jj + 1) * u_ice(ji - 1, jj + 1)))
            zTauE = zfV(ji, jj) + zTauV_ia(ji, jj) + zCory(ji, jj) + zspgV(ji, jj) + ztauy_oi(ji, jj)
            rswitch = 1._wp - MIN(1._wp, ABS(SIGN(1._wp, zTauE - zTauV_ib(ji, jj)) - SIGN(1._wp, zTauE)))
            IF (ln_aEVP) THEN
              v_ice(ji, jj) = ((rswitch * (zmV_t(ji, jj) * (zbeta(ji, jj) * v_ice(ji, jj) + v_ice_b(ji, jj)) + zTauE + zTauO * v_ice(ji, jj)) / MAX(zepsi, zmV_t(ji, jj) * (zbeta(ji, jj) + 1._wp) + zTauO - zTauB) + (1._wp - rswitch) * v_ice(ji, jj) * MAX(0._wp, 1._wp - zdtevp * rn_lfrelax)) * zswitchV(ji, jj) + v_oce(ji, jj) * 0.01_wp * (1._wp - zswitchV(ji, jj))) * zmaskV(ji, jj)
            ELSE
              v_ice(ji, jj) = ((rswitch * (zmV_t(ji, jj) * v_ice(ji, jj) + zTauE + zTauO * v_ice(ji, jj)) / MAX(zepsi, zmV_t(ji, jj) + zTauO - zTauB) + (1._wp - rswitch) * v_ice(ji, jj) * MAX(0._wp, 1._wp - zdtevp * rn_lfrelax)) * zswitchV(ji, jj) + v_oce(ji, jj) * 0.01_wp * (1._wp - zswitchV(ji, jj))) * zmaskV(ji, jj)
            END IF
          END DO
        END DO
        !$ACC END KERNELS
        CALL profile_psy_data12 % PreStart('ice_dyn_rhg_evp', 'r12', 0, 0)
        CALL lbc_lnk('icedyn_rhg_evp', v_ice, 'V', - 1.)
        IF (ln_bdy .AND. ll_bdy_substep) CALL bdy_ice_dyn('V')
        CALL profile_psy_data12 % PostEnd
      END IF
    END DO
    CALL profile_psy_data13 % PreStart('ice_dyn_rhg_evp', 'r13', 0, 0)
    IF (ln_bdy .AND. .NOT. ll_bdy_substep) THEN
      CALL bdy_ice_dyn('U')
      CALL bdy_ice_dyn('V')
    END IF
    CALL profile_psy_data13 % PostEnd
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 1, jpjm1
      DO ji = 1, jpim1
        zds(ji, jj) = ((u_ice(ji, jj + 1) * r1_e1u(ji, jj + 1) - u_ice(ji, jj) * r1_e1u(ji, jj)) * e1f(ji, jj) * e1f(ji, jj) + (v_ice(ji + 1, jj) * r1_e2v(ji + 1, jj) - v_ice(ji, jj) * r1_e2v(ji, jj)) * e2f(ji, jj) * e2f(ji, jj)) * r1_e1e2f(ji, jj) * zfmask(ji, jj)
      END DO
    END DO
    !$ACC loop independent collapse(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zdt = ((u_ice(ji, jj) * r1_e2u(ji, jj) - u_ice(ji - 1, jj) * r1_e2u(ji - 1, jj)) * e2t(ji, jj) * e2t(ji, jj) - (v_ice(ji, jj) * r1_e1v(ji, jj) - v_ice(ji, jj - 1) * r1_e1v(ji, jj - 1)) * e1t(ji, jj) * e1t(ji, jj)) * r1_e1e2t(ji, jj)
        zdt2 = zdt * zdt
        zds2 = (zds(ji, jj) * zds(ji, jj) * e1e2f(ji, jj) + zds(ji - 1, jj) * zds(ji - 1, jj) * e1e2f(ji - 1, jj) + zds(ji, jj - 1) * zds(ji, jj - 1) * e1e2f(ji, jj - 1) + zds(ji - 1, jj - 1) * zds(ji - 1, jj - 1) * e1e2f(ji - 1, jj - 1)) * 0.25_wp * r1_e1e2t(ji, jj)
        pshear_i(ji, jj) = SQRT(zdt2 + zds2)
        pdivu_i(ji, jj) = (e2u(ji, jj) * u_ice(ji, jj) - e2u(ji - 1, jj) * u_ice(ji - 1, jj) + e1v(ji, jj) * v_ice(ji, jj) - e1v(ji, jj - 1) * v_ice(ji, jj - 1)) * r1_e1e2t(ji, jj)
        zdelta = SQRT(pdivu_i(ji, jj) * pdivu_i(ji, jj) + (zdt2 + zds2) * z1_ecc2)
        rswitch = 1._wp - MAX(0._wp, SIGN(1._wp, - zdelta))
        pdelta_i(ji, jj) = zdelta + rn_creepl * rswitch
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data14 % PreStart('ice_dyn_rhg_evp', 'r14', 0, 0)
    CALL lbc_lnk_multi('icedyn_rhg_evp', pshear_i, 'T', 1., pdivu_i, 'T', 1., pdelta_i, 'T', 1.)
    CALL lbc_lnk_multi('icedyn_rhg_evp', zs1, 'T', 1., zs2, 'T', 1., zs12, 'F', 1.)
    CALL profile_psy_data14 % PostEnd
    !$ACC KERNELS
    pstress1_i(:, :) = zs1(:, :)
    pstress2_i(:, :) = zs2(:, :)
    pstress12_i(:, :) = zs12(:, :)
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zswi(ji, jj) = MAX(0._wp, SIGN(1._wp, at_i(ji, jj) - epsi06))
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data15 % PreStart('ice_dyn_rhg_evp', 'r15', 0, 0)
    IF (iom_use('icediv')) CALL iom_put("icediv", pdivu_i(:, :) * zswi(:, :))
    IF (iom_use('iceshe')) CALL iom_put("iceshe", pshear_i(:, :) * zswi(:, :))
    IF (iom_use('icestr')) CALL iom_put("icestr", strength(:, :) * zswi(:, :))
    CALL profile_psy_data15 % PostEnd
    IF (iom_use('isig1') .OR. iom_use('isig2') .OR. iom_use('isig3')) THEN
      ALLOCATE(zsig1(jpi, jpj), zsig2(jpi, jpj), zsig3(jpi, jpj))
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zdum1 = (zswi(ji - 1, jj) * pstress12_i(ji - 1, jj) + zswi(ji, jj - 1) * pstress12_i(ji, jj - 1) + zswi(ji, jj) * pstress12_i(ji, jj) + zswi(ji - 1, jj - 1) * pstress12_i(ji - 1, jj - 1)) / MAX(1._wp, zswi(ji - 1, jj) + zswi(ji, jj - 1) + zswi(ji, jj) + zswi(ji - 1, jj - 1))
          zshear = SQRT(pstress2_i(ji, jj) * pstress2_i(ji, jj) + 4._wp * zdum1 * zdum1)
          zdum2 = zswi(ji, jj) / MAX(1._wp, strength(ji, jj))
          zsig1(ji, jj) = 0.5_wp * zdum2 * (pstress1_i(ji, jj))
          zsig2(ji, jj) = 0.5_wp * zdum2 * (zshear)
          zsig3(ji, jj) = zdum2 ** 2 * ((pstress1_i(ji, jj) + strength(ji, jj)) ** 2 + (rn_ecc * zshear) ** 2)
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data16 % PreStart('ice_dyn_rhg_evp', 'r16', 0, 0)
      CALL lbc_lnk_multi('icedyn_rhg_evp', zsig1, 'T', 1., zsig2, 'T', 1., zsig3, 'T', 1.)
      IF (iom_use('isig1')) CALL iom_put("isig1", zsig1)
      IF (iom_use('isig2')) CALL iom_put("isig2", zsig2)
      IF (iom_use('isig3')) CALL iom_put("isig3", zsig3)
      DEALLOCATE(zsig1, zsig2, zsig3)
      CALL profile_psy_data16 % PostEnd
    END IF
    IF (iom_use('normstr') .OR. iom_use('sheastr') .OR. iom_use('dssh_dx') .OR. iom_use('dssh_dy') .OR. iom_use('corstrx') .OR. iom_use('corstry') .OR. iom_use('intstrx') .OR. iom_use('intstry') .OR. iom_use('utau_oi') .OR. iom_use('vtau_oi') .OR. iom_use('xmtrpice') .OR. iom_use('ymtrpice') .OR. iom_use('xmtrpsnw') .OR. iom_use('ymtrpsnw') .OR. iom_use('xatrp') .OR. iom_use('yatrp')) THEN
      ALLOCATE(zdiag_sig1(jpi, jpj), zdiag_sig2(jpi, jpj), zdiag_dssh_dx(jpi, jpj), zdiag_dssh_dy(jpi, jpj), zdiag_corstrx(jpi, jpj), zdiag_corstry(jpi, jpj), zdiag_intstrx(jpi, jpj), zdiag_intstry(jpi, jpj), zdiag_utau_oi(jpi, jpj), zdiag_vtau_oi(jpi, jpj), zdiag_xmtrp_ice(jpi, jpj), zdiag_ymtrp_ice(jpi, jpj), zdiag_xmtrp_snw(jpi, jpj), zdiag_ymtrp_snw(jpi, jpj), zdiag_xatrp(jpi, jpj), zdiag_yatrp(jpi, jpj))
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          rswitch = MAX(0._wp, SIGN(1._wp, at_i(ji, jj) - epsi06))
          zdiag_sig1(ji, jj) = (zs1(ji, jj) + zs2(ji, jj)) * rswitch
          zdiag_sig2(ji, jj) = SQRT((zs1(ji, jj) - zs2(ji, jj)) ** 2 + 4 * zs12(ji, jj) ** 2) * rswitch
          zdiag_dssh_dx(ji, jj) = zspgU(ji, jj) * rswitch
          zdiag_dssh_dy(ji, jj) = zspgV(ji, jj) * rswitch
          zdiag_corstrx(ji, jj) = zCorx(ji, jj) * rswitch
          zdiag_corstry(ji, jj) = zCory(ji, jj) * rswitch
          zdiag_intstrx(ji, jj) = zfU(ji, jj) * rswitch
          zdiag_intstry(ji, jj) = zfV(ji, jj) * rswitch
          zdiag_utau_oi(ji, jj) = ztaux_oi(ji, jj) * rswitch
          zdiag_vtau_oi(ji, jj) = ztauy_oi(ji, jj) * rswitch
          zfac_x = 0.5 * u_ice(ji, jj) * e2u(ji, jj) * rswitch
          zfac_y = 0.5 * v_ice(ji, jj) * e1v(ji, jj) * rswitch
          zdiag_xmtrp_ice(ji, jj) = rhoi * zfac_x * (vt_i(ji + 1, jj) + vt_i(ji, jj))
          zdiag_ymtrp_ice(ji, jj) = rhoi * zfac_y * (vt_i(ji, jj + 1) + vt_i(ji, jj))
          zdiag_xmtrp_snw(ji, jj) = rhos * zfac_x * (vt_s(ji + 1, jj) + vt_s(ji, jj))
          zdiag_ymtrp_snw(ji, jj) = rhos * zfac_y * (vt_s(ji, jj + 1) + vt_s(ji, jj))
          zdiag_xatrp(ji, jj) = zfac_x * (at_i(ji + 1, jj) + at_i(ji, jj))
          zdiag_yatrp(ji, jj) = zfac_y * (at_i(ji, jj + 1) + at_i(ji, jj))
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data17 % PreStart('ice_dyn_rhg_evp', 'r17', 0, 0)
      CALL lbc_lnk_multi('icedyn_rhg_evp', zdiag_sig1, 'T', 1., zdiag_sig2, 'T', 1., zdiag_dssh_dx, 'U', - 1., zdiag_dssh_dy, 'V', - 1., zdiag_corstrx, 'U', - 1., zdiag_corstry, 'V', - 1., zdiag_intstrx, 'U', - 1., zdiag_intstry, 'V', - 1.)
      CALL lbc_lnk_multi('icedyn_rhg_evp', zdiag_utau_oi, 'U', - 1., zdiag_vtau_oi, 'V', - 1., zdiag_xmtrp_ice, 'U', - 1., zdiag_xmtrp_snw, 'U', - 1., zdiag_xatrp, 'U', - 1., zdiag_ymtrp_ice, 'V', - 1., zdiag_ymtrp_snw, 'V', - 1., zdiag_yatrp, 'V', - 1.)
      IF (iom_use('normstr')) CALL iom_put('normstr', zdiag_sig1(:, :))
      IF (iom_use('sheastr')) CALL iom_put('sheastr', zdiag_sig2(:, :))
      IF (iom_use('dssh_dx')) CALL iom_put('dssh_dx', zdiag_dssh_dx(:, :))
      IF (iom_use('dssh_dy')) CALL iom_put('dssh_dy', zdiag_dssh_dy(:, :))
      IF (iom_use('corstrx')) CALL iom_put('corstrx', zdiag_corstrx(:, :))
      IF (iom_use('corstry')) CALL iom_put('corstry', zdiag_corstry(:, :))
      IF (iom_use('intstrx')) CALL iom_put('intstrx', zdiag_intstrx(:, :))
      IF (iom_use('intstry')) CALL iom_put('intstry', zdiag_intstry(:, :))
      IF (iom_use('utau_oi')) CALL iom_put('utau_oi', zdiag_utau_oi(:, :))
      IF (iom_use('vtau_oi')) CALL iom_put('vtau_oi', zdiag_vtau_oi(:, :))
      IF (iom_use('xmtrpice')) CALL iom_put('xmtrpice', zdiag_xmtrp_ice(:, :))
      IF (iom_use('ymtrpice')) CALL iom_put('ymtrpice', zdiag_ymtrp_ice(:, :))
      IF (iom_use('xmtrpsnw')) CALL iom_put('xmtrpsnw', zdiag_xmtrp_snw(:, :))
      IF (iom_use('ymtrpsnw')) CALL iom_put('ymtrpsnw', zdiag_ymtrp_snw(:, :))
      IF (iom_use('xatrp')) CALL iom_put('xatrp', zdiag_xatrp(:, :))
      IF (iom_use('yatrp')) CALL iom_put('yatrp', zdiag_yatrp(:, :))
      DEALLOCATE(zdiag_sig1, zdiag_sig2, zdiag_dssh_dx, zdiag_dssh_dy, zdiag_corstrx, zdiag_corstry, zdiag_intstrx, zdiag_intstry, zdiag_utau_oi, zdiag_vtau_oi, zdiag_xmtrp_ice, zdiag_ymtrp_ice, zdiag_xmtrp_snw, zdiag_ymtrp_snw, zdiag_xatrp, zdiag_yatrp)
      CALL profile_psy_data17 % PostEnd
    END IF
  END SUBROUTINE ice_dyn_rhg_evp
  SUBROUTINE rhg_evp_rst(cdrw, kt)
    CHARACTER(LEN = *), INTENT(IN) :: cdrw
    INTEGER, OPTIONAL, INTENT(IN) :: kt
    INTEGER :: iter
    INTEGER :: id1, id2, id3
    IF (TRIM(cdrw) == 'READ') THEN
      IF (ln_rstart) THEN
        id1 = iom_varid(numrir, 'stress1_i', ldstop = .FALSE.)
        id2 = iom_varid(numrir, 'stress2_i', ldstop = .FALSE.)
        id3 = iom_varid(numrir, 'stress12_i', ldstop = .FALSE.)
        IF (MIN(id1, id2, id3) > 0) THEN
          CALL iom_get(numrir, jpdom_autoglo, 'stress1_i', stress1_i)
          CALL iom_get(numrir, jpdom_autoglo, 'stress2_i', stress2_i)
          CALL iom_get(numrir, jpdom_autoglo, 'stress12_i', stress12_i)
        ELSE
          IF (lwp) WRITE(numout, *)
          IF (lwp) WRITE(numout, *) '   ==>>>   previous run without rheology, set stresses to 0'
          !$ACC KERNELS
          stress1_i(:, :) = 0._wp
          stress2_i(:, :) = 0._wp
          stress12_i(:, :) = 0._wp
          !$ACC END KERNELS
        END IF
      ELSE
        IF (lwp) WRITE(numout, *)
        IF (lwp) WRITE(numout, *) '   ==>>>   start from rest: set stresses to 0'
        !$ACC KERNELS
        stress1_i(:, :) = 0._wp
        stress2_i(:, :) = 0._wp
        stress12_i(:, :) = 0._wp
        !$ACC END KERNELS
      END IF
    ELSE IF (TRIM(cdrw) == 'WRITE') THEN
      IF (lwp) WRITE(numout, *) '---- rhg-rst ----'
      iter = kt + nn_fsbc - 1
      CALL iom_rstput(iter, nitrst, numriw, 'stress1_i', stress1_i)
      CALL iom_rstput(iter, nitrst, numriw, 'stress2_i', stress2_i)
      CALL iom_rstput(iter, nitrst, numriw, 'stress12_i', stress12_i)
    END IF
  END SUBROUTINE rhg_evp_rst
END MODULE icedyn_rhg_evp
