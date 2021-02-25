MODULE tramle
  USE oce
  USE dom_oce
  USE phycst
  USE zdfmxl
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_mle_trp
  PUBLIC :: tra_mle_init
  LOGICAL, PUBLIC :: ln_mle
  INTEGER :: nn_mle
  INTEGER :: nn_mld_uv
  INTEGER :: nn_conv
  REAL(KIND = wp) :: rn_ce
  REAL(KIND = wp) :: rn_lf
  REAL(KIND = wp) :: rn_time
  REAL(KIND = wp) :: rn_lat
  REAL(KIND = wp) :: rn_rho_c_mle
  REAL(KIND = wp) :: r5_21 = 5.E0 / 21.E0
  REAL(KIND = wp) :: rb_c
  REAL(KIND = wp) :: rc_f
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: rfu, rfv
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: r1_ft
  CONTAINS
  SUBROUTINE tra_mle_trp(kt, kit000, pu, pv, pw, cdtype)
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kit000
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pu
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pv
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pw
    INTEGER :: ji, jj, jk
    INTEGER :: ii, ij, ik, ikmax
    REAL(KIND = wp) :: zcuw, zmuw, zc
    REAL(KIND = wp) :: zcvw, zmvw
    INTEGER, DIMENSION(jpi, jpj) :: inml_mle
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zpsim_u, zpsim_v, zmld, zbm, zhu, zhv, zn2, zLf_NH, zLf_MH
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zpsi_uw, zpsi_vw
    inml_mle(:, :) = mbkt(:, :) + 1
    IF (nla10 > 0) THEN
      !$OMP parallel default(shared), private(ji,jj,jk)
      !$OMP do schedule(static)
      DO jk = jpkm1, nlb10, - 1
        DO jj = 1, jpj
          DO ji = 1, jpi
            IF (rhop(ji, jj, jk) > rhop(ji, jj, nla10) + rn_rho_c_mle) inml_mle(ji, jj) = jk
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
    END IF
    ikmax = MIN(MAXVAL(inml_mle(:, :)), jpkm1)
    zmld(:, :) = 0._wp
    zbm(:, :) = 0._wp
    zn2(:, :) = 0._wp
    !$OMP parallel default(shared), private(ji,jj,jk,zc)
    !$OMP do schedule(static)
    DO jk = 1, ikmax
      DO jj = 1, jpj
        DO ji = 1, jpi
          zc = e3t_n(ji, jj, jk) * REAL(MIN(MAX(0, inml_mle(ji, jj) - jk), 1))
          zmld(ji, jj) = zmld(ji, jj) + zc
          zbm(ji, jj) = zbm(ji, jj) + zc * (rau0 - rhop(ji, jj, jk)) * r1_rau0
          zn2(ji, jj) = zn2(ji, jj) + zc * (rn2(ji, jj, jk) + rn2(ji, jj, jk + 1)) * 0.5_wp
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    SELECT CASE (nn_mld_uv)
    CASE (0)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zhu(ji, jj) = MIN(zmld(ji + 1, jj), zmld(ji, jj))
          zhv(ji, jj) = MIN(zmld(ji, jj + 1), zmld(ji, jj))
        END DO
      END DO
    CASE (1)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zhu(ji, jj) = (zmld(ji + 1, jj) + zmld(ji, jj)) * 0.5_wp
          zhv(ji, jj) = (zmld(ji, jj + 1) + zmld(ji, jj)) * 0.5_wp
        END DO
      END DO
    CASE (2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zhu(ji, jj) = MAX(zmld(ji + 1, jj), zmld(ji, jj))
          zhv(ji, jj) = MAX(zmld(ji, jj + 1), zmld(ji, jj))
        END DO
      END DO
    END SELECT
    zbm(:, :) = + grav * zbm(:, :) / MAX(e3t_n(:, :, 1), zmld(:, :))
    IF (nn_mle == 0) THEN
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zpsim_u(ji, jj) = rn_ce * zhu(ji, jj) * zhu(ji, jj) * e2_e1u(ji, jj) * (zbm(ji + 1, jj) - zbm(ji, jj)) * MIN(111.E3_wp, &
&e1u(ji, jj)) / (MAX(rn_lf * rfu(ji, jj), SQRT(rb_c * zhu(ji, jj))))
          zpsim_v(ji, jj) = rn_ce * zhv(ji, jj) * zhv(ji, jj) * e1_e2v(ji, jj) * (zbm(ji, jj + 1) - zbm(ji, jj)) * MIN(111.E3_wp, &
&e2v(ji, jj)) / (MAX(rn_lf * rfv(ji, jj), SQRT(rb_c * zhv(ji, jj))))
        END DO
      END DO
    ELSE IF (nn_mle == 1) THEN
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zpsim_u(ji, jj) = rc_f * zhu(ji, jj) * zhu(ji, jj) * e2_e1u(ji, jj) * (zbm(ji + 1, jj) - zbm(ji, jj)) * MIN(111.E3_wp, &
&e1u(ji, jj))
          zpsim_v(ji, jj) = rc_f * zhv(ji, jj) * zhv(ji, jj) * e1_e2v(ji, jj) * (zbm(ji, jj + 1) - zbm(ji, jj)) * MIN(111.E3_wp, &
&e2v(ji, jj))
        END DO
      END DO
    END IF
    IF (nn_conv == 1) THEN
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          IF (MIN(zn2(ji, jj), zn2(ji + 1, jj)) < 0._wp) zpsim_u(ji, jj) = 0._wp
          IF (MIN(zn2(ji, jj), zn2(ji, jj + 1)) < 0._wp) zpsim_v(ji, jj) = 0._wp
        END DO
      END DO
    END IF
    DO jj = 1, jpjm1
      DO ji = 1, jpim1
        zhu(ji, jj) = 1._wp / zhu(ji, jj)
        zhv(ji, jj) = 1._wp / zhv(ji, jj)
      END DO
    END DO
    zpsi_uw(:, :, :) = 0._wp
    zpsi_vw(:, :, :) = 0._wp
    !$OMP parallel default(shared), private(ji,jj,jk,zcuw,zcvw,zmuw,zmvw)
    !$OMP do schedule(static)
    DO jk = 2, ikmax
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zcuw = 1._wp - (gdepw_n(ji + 1, jj, jk) + gdepw_n(ji, jj, jk)) * zhu(ji, jj)
          zcvw = 1._wp - (gdepw_n(ji, jj + 1, jk) + gdepw_n(ji, jj, jk)) * zhv(ji, jj)
          zcuw = zcuw * zcuw
          zcvw = zcvw * zcvw
          zmuw = MAX(0._wp, (1._wp - zcuw) * (1._wp + r5_21 * zcuw))
          zmvw = MAX(0._wp, (1._wp - zcvw) * (1._wp + r5_21 * zcvw))
          zpsi_uw(ji, jj, jk) = zpsim_u(ji, jj) * zmuw * umask(ji, jj, jk)
          zpsi_vw(ji, jj, jk) = zpsim_v(ji, jj) * zmvw * vmask(ji, jj, jk)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP do schedule(static)
    DO jk = 1, ikmax
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          pu(ji, jj, jk) = pu(ji, jj, jk) + (zpsi_uw(ji, jj, jk) - zpsi_uw(ji, jj, jk + 1))
          pv(ji, jj, jk) = pv(ji, jj, jk) + (zpsi_vw(ji, jj, jk) - zpsi_vw(ji, jj, jk + 1))
        END DO
      END DO
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          pw(ji, jj, jk) = pw(ji, jj, jk) - (zpsi_uw(ji, jj, jk) - zpsi_uw(ji - 1, jj, jk) + zpsi_vw(ji, jj, jk) - zpsi_vw(ji, jj &
&- 1, jk))
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    IF (cdtype == 'TRA') THEN
      zLf_NH(:, :) = SQRT(rb_c * zmld(:, :)) * r1_ft(:, :)
      CALL iom_put("Lf_NHpf", zLf_NH)
      !$OMP parallel default(shared), private(jk)
      !$OMP do schedule(static)
      DO jk = 1, ikmax + 1
        zpsi_uw(:, :, jk) = zpsi_uw(:, :, jk) * r1_e2u(:, :)
        zpsi_vw(:, :, jk) = zpsi_vw(:, :, jk) * r1_e1v(:, :)
      END DO
      !$OMP end do
      !$OMP end parallel
      CALL iom_put("psiu_mle", zpsi_uw)
      CALL iom_put("psiv_mle", zpsi_vw)
    END IF
  END SUBROUTINE tra_mle_trp
  SUBROUTINE tra_mle_init
    INTEGER :: ji, jj, jk
    INTEGER :: ierr
    INTEGER :: ios
    REAL(KIND = wp) :: z1_t2, zfu, zfv
    NAMELIST /namtra_mle/ ln_mle, nn_mle, rn_ce, rn_lf, rn_time, rn_lat, nn_mld_uv, nn_conv, rn_rho_c_mle
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namtra_mle, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namtra_mle in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namtra_mle, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namtra_mle in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namtra_mle)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'tra_mle_init : mixed layer eddy (MLE) advection acting on tracers'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namtra_mle : mixed layer eddy advection applied on tracers'
      WRITE(numout, FMT = *) '      use mixed layer eddy (MLE, i.e. Fox-Kemper param) (T/F)      ln_mle       = ', ln_mle
      WRITE(numout, FMT = *) '         MLE type: =0 standard Fox-Kemper ; =1 new formulation        nn_mle    = ', nn_mle
      WRITE(numout, FMT = *) '         magnitude of the MLE (typical value: 0.06 to 0.08)           rn_ce     = ', rn_ce
      WRITE(numout, FMT = *) '         scale of ML front (ML radius of deformation) (rn_mle=0)      rn_lf     = ', rn_lf, 'm'
      WRITE(numout, FMT = *) '         maximum time scale of MLE                    (rn_mle=0)      rn_time   = ', rn_time, 's'
      WRITE(numout, FMT = *) '         reference latitude (degrees) of MLE coef.    (rn_mle=1)      rn_lat    = ', rn_lat, 'deg'
      WRITE(numout, FMT = *) '         space interp. of MLD at u-(v-)pts (0=min,1=averaged,2=max)   nn_mld_uv = ', nn_mld_uv
      WRITE(numout, FMT = *) '         =1 no MLE in case of convection ; =0 always MLE              nn_conv   = ', nn_conv
      WRITE(numout, FMT = *) '         Density difference used to define ML for FK              rn_rho_c_mle  = ', rn_rho_c_mle
    END IF
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      IF (ln_mle) THEN
        WRITE(numout, FMT = *) '   ==>>>   Mixed Layer Eddy induced transport added to tracer advection'
        IF (nn_mle == 0) WRITE(numout, FMT = *) '              Fox-Kemper et al 2010 formulation'
        IF (nn_mle == 1) WRITE(numout, FMT = *) '              New formulation'
      ELSE
        WRITE(numout, FMT = *) '   ==>>>   Mixed Layer Eddy parametrisation NOT used'
      END IF
    END IF
    IF (ln_mle) THEN
      rb_c = grav * rn_rho_c_mle / rau0
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) '      ML buoyancy criteria = ', rb_c, ' m/s2 '
      IF (lwp) WRITE(numout, FMT = *) '      associated ML density criteria defined in zdfmxl = ', rho_c, 'kg/m3'
      IF (nn_mle == 0) THEN
        ALLOCATE(rfu(jpi, jpj), rfv(jpi, jpj), STAT = ierr)
        IF (ierr /= 0) CALL ctl_stop('tra_adv_mle_init: failed to allocate arrays')
        z1_t2 = 1._wp / (rn_time * rn_time)
        DO jj = 2, jpj
          DO ji = 2, jpi
            zfu = (ff_f(ji, jj) + ff_f(ji, jj - 1)) * 0.5_wp
            zfv = (ff_f(ji, jj) + ff_f(ji - 1, jj)) * 0.5_wp
            rfu(ji, jj) = SQRT(zfu * zfu + z1_t2)
            rfv(ji, jj) = SQRT(zfv * zfv + z1_t2)
          END DO
        END DO
        CALL lbc_lnk_multi('tramle', rfu, 'U', 1., rfv, 'V', 1.)
      ELSE IF (nn_mle == 1) THEN
        rc_f = rn_ce / (5.E3_wp * 2._wp * omega * SIN(rad * rn_lat))
      END IF
      ALLOCATE(r1_ft(jpi, jpj), STAT = ierr)
      IF (ierr /= 0) CALL ctl_stop('tra_adv_mle_init: failed to allocate r1_ft array')
      z1_t2 = 1._wp / (rn_time * rn_time)
      r1_ft(:, :) = 1._wp / SQRT(ff_t(:, :) * ff_t(:, :) + z1_t2)
    END IF
  END SUBROUTINE tra_mle_init
END MODULE tramle