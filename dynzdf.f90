MODULE dynzdf
  USE oce
  USE phycst
  USE dom_oce
  USE sbc_oce
  USE zdf_oce
  USE zdfdrg
  USE dynadv, ONLY: ln_dynadv_vec
  USE dynldf_iso, ONLY: akzu, akzv
  USE ldfdyn
  USE trd_oce
  USE trddyn
  USE in_out_manager
  USE lib_mpp
  USE prtctl
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_zdf
  REAL(KIND = wp) :: r_vvl
  CONTAINS
  SUBROUTINE dyn_zdf(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: iku, ikv
    REAL(KIND = wp) :: zzwi, ze3ua, zdt
    REAL(KIND = wp) :: zzws, ze3va
    REAL(KIND = wp) :: z1_e3ua, z1_e3va
    REAL(KIND = wp) :: zWu, zWv
    REAL(KIND = wp) :: zWui, zWvi
    REAL(KIND = wp) :: zWus, zWvs
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwi, zwd, zws
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ztrdu, ztrdv
    IF (ln_timing) CALL timing_start('dyn_zdf')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'dyn_zdf_imp : vertical momentum diffusion implicit operator'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~ '
      IF (ln_linssh) THEN
        r_vvl = 0._wp
      ELSE
        r_vvl = 1._wp
      END IF
    END IF
    IF (neuler == 0 .AND. kt == nit000) THEN
      r2dt = rdt
    ELSE IF (kt <= nit000 + 1) THEN
      r2dt = 2. * rdt
    END IF
    IF (.NOT. ln_drgimp) CALL zdf_drg_exp(kt, ub, vb, ua, va)
    IF (l_trddyn) THEN
      ALLOCATE(ztrdu(jpi, jpj, jpk), ztrdv(jpi, jpj, jpk))
      ztrdu(:, :, :) = ua(:, :, :)
      ztrdv(:, :, :) = va(:, :, :)
    END IF
    IF (ln_dynadv_vec .OR. ln_linssh) THEN
      !$OMP parallel default(shared), private(jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        ua(:, :, jk) = (ub(:, :, jk) + r2dt * ua(:, :, jk)) * umask(:, :, jk)
        va(:, :, jk) = (vb(:, :, jk) + r2dt * va(:, :, jk)) * vmask(:, :, jk)
      END DO
      !$OMP end do
      !$OMP end parallel
    ELSE
      !$OMP parallel default(shared), private(jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        ua(:, :, jk) = (e3u_b(:, :, jk) * ub(:, :, jk) + r2dt * e3u_n(:, :, jk) * ua(:, :, jk)) / e3u_a(:, :, jk) * umask(:, :, jk)
        va(:, :, jk) = (e3v_b(:, :, jk) * vb(:, :, jk) + r2dt * e3v_n(:, :, jk) * va(:, :, jk)) / e3v_a(:, :, jk) * vmask(:, :, jk)
      END DO
      !$OMP end do
      !$OMP end parallel
    END IF
    IF (ln_drgimp .AND. ln_dynspg_ts) THEN
      !$OMP parallel default(shared), private(jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        ua(:, :, jk) = (ua(:, :, jk) - ua_b(:, :)) * umask(:, :, jk)
        va(:, :, jk) = (va(:, :, jk) - va_b(:, :)) * vmask(:, :, jk)
      END DO
      !$OMP end do
      !$OMP end parallel
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          iku = mbku(ji, jj)
          ikv = mbkv(ji, jj)
          ze3ua = (1._wp - r_vvl) * e3u_n(ji, jj, iku) + r_vvl * e3u_a(ji, jj, iku)
          ze3va = (1._wp - r_vvl) * e3v_n(ji, jj, ikv) + r_vvl * e3v_a(ji, jj, ikv)
          ua(ji, jj, iku) = ua(ji, jj, iku) + r2dt * 0.5 * (rCdU_bot(ji + 1, jj) + rCdU_bot(ji, jj)) * ua_b(ji, jj) / ze3ua
          va(ji, jj, ikv) = va(ji, jj, ikv) + r2dt * 0.5 * (rCdU_bot(ji, jj + 1) + rCdU_bot(ji, jj)) * va_b(ji, jj) / ze3va
        END DO
      END DO
      IF (ln_isfcav) THEN
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            iku = miku(ji, jj)
            ikv = mikv(ji, jj)
            ze3ua = (1._wp - r_vvl) * e3u_n(ji, jj, iku) + r_vvl * e3u_a(ji, jj, iku)
            ze3va = (1._wp - r_vvl) * e3v_n(ji, jj, ikv) + r_vvl * e3v_a(ji, jj, ikv)
            ua(ji, jj, iku) = ua(ji, jj, iku) + r2dt * 0.5 * (rCdU_top(ji + 1, jj) + rCdU_top(ji, jj)) * ua_b(ji, jj) / ze3ua
            va(ji, jj, ikv) = va(ji, jj, ikv) + r2dt * 0.5 * (rCdU_top(ji + 1, jj) + rCdU_top(ji, jj)) * va_b(ji, jj) / ze3va
          END DO
        END DO
      END IF
    END IF
    zdt = r2dt * 0.5
    IF (ln_zad_Aimp) THEN
      SELECT CASE (nldf_dyn)
      CASE (np_lap_i)
        !$OMP parallel default(shared), private(ji,jj,jk,ze3ua,zwui,zwus,zzwi,zzws)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              ze3ua = (1._wp - r_vvl) * e3u_n(ji, jj, jk) + r_vvl * e3u_a(ji, jj, jk)
              zzwi = - zdt * (avm(ji + 1, jj, jk) + avm(ji, jj, jk) + akzu(ji, jj, jk)) / (ze3ua * e3uw_n(ji, jj, jk)) * &
&wumask(ji, jj, jk)
              zzws = - zdt * (avm(ji + 1, jj, jk + 1) + avm(ji, jj, jk + 1) + akzu(ji, jj, jk + 1)) / (ze3ua * e3uw_n(ji, jj, jk + &
&1)) * wumask(ji, jj, jk + 1)
              zWui = 0.5_wp * (wi(ji, jj, jk) + wi(ji + 1, jj, jk))
              zWus = 0.5_wp * (wi(ji, jj, jk + 1) + wi(ji + 1, jj, jk + 1))
              zwi(ji, jj, jk) = zzwi + zdt * MIN(zWui, 0._wp)
              zws(ji, jj, jk) = zzws - zdt * MAX(zWus, 0._wp)
              zwd(ji, jj, jk) = 1._wp - zzwi - zzws + zdt * (MAX(zWui, 0._wp) - MIN(zWus, 0._wp))
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
      CASE DEFAULT
        !$OMP parallel default(shared), private(ji,jj,jk,ze3ua,zwui,zwus,zzwi,zzws)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              ze3ua = (1._wp - r_vvl) * e3u_n(ji, jj, jk) + r_vvl * e3u_a(ji, jj, jk)
              zzwi = - zdt * (avm(ji + 1, jj, jk) + avm(ji, jj, jk)) / (ze3ua * e3uw_n(ji, jj, jk)) * wumask(ji, jj, jk)
              zzws = - zdt * (avm(ji + 1, jj, jk + 1) + avm(ji, jj, jk + 1)) / (ze3ua * e3uw_n(ji, jj, jk + 1)) * wumask(ji, jj, &
&jk + 1)
              zWui = 0.5_wp * (wi(ji, jj, jk) + wi(ji + 1, jj, jk))
              zWus = 0.5_wp * (wi(ji, jj, jk + 1) + wi(ji + 1, jj, jk + 1))
              zwi(ji, jj, jk) = zzwi + zdt * MIN(zWui, 0._wp)
              zws(ji, jj, jk) = zzws - zdt * MAX(zWus, 0._wp)
              zwd(ji, jj, jk) = 1._wp - zzwi - zzws + zdt * (MAX(zWui, 0._wp) - MIN(zWus, 0._wp))
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
      END SELECT
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zwi(ji, jj, 1) = 0._wp
          ze3ua = (1._wp - r_vvl) * e3u_n(ji, jj, 1) + r_vvl * e3u_a(ji, jj, 1)
          zzws = - zdt * (avm(ji + 1, jj, 2) + avm(ji, jj, 2)) / (ze3ua * e3uw_n(ji, jj, 2)) * wumask(ji, jj, 2)
          zWus = 0.5_wp * (wi(ji, jj, 2) + wi(ji + 1, jj, 2))
          zws(ji, jj, 1) = zzws - zdt * MAX(zWus, 0._wp)
          zwd(ji, jj, 1) = 1._wp - zzws - zdt * (MIN(zWus, 0._wp))
        END DO
      END DO
    ELSE
      SELECT CASE (nldf_dyn)
      CASE (np_lap_i)
        !$OMP parallel default(shared), private(ji,jj,jk,ze3ua,zzwi,zzws)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              ze3ua = (1._wp - r_vvl) * e3u_n(ji, jj, jk) + r_vvl * e3u_a(ji, jj, jk)
              zzwi = - zdt * (avm(ji + 1, jj, jk) + avm(ji, jj, jk) + akzu(ji, jj, jk)) / (ze3ua * e3uw_n(ji, jj, jk)) * &
&wumask(ji, jj, jk)
              zzws = - zdt * (avm(ji + 1, jj, jk + 1) + avm(ji, jj, jk + 1) + akzu(ji, jj, jk + 1)) / (ze3ua * e3uw_n(ji, jj, jk + &
&1)) * wumask(ji, jj, jk + 1)
              zwi(ji, jj, jk) = zzwi
              zws(ji, jj, jk) = zzws
              zwd(ji, jj, jk) = 1._wp - zzwi - zzws
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
      CASE DEFAULT
        !$OMP parallel default(shared), private(ji,jj,jk,ze3ua,zzwi,zzws)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              ze3ua = (1._wp - r_vvl) * e3u_n(ji, jj, jk) + r_vvl * e3u_a(ji, jj, jk)
              zzwi = - zdt * (avm(ji + 1, jj, jk) + avm(ji, jj, jk)) / (ze3ua * e3uw_n(ji, jj, jk)) * wumask(ji, jj, jk)
              zzws = - zdt * (avm(ji + 1, jj, jk + 1) + avm(ji, jj, jk + 1)) / (ze3ua * e3uw_n(ji, jj, jk + 1)) * wumask(ji, jj, &
&jk + 1)
              zwi(ji, jj, jk) = zzwi
              zws(ji, jj, jk) = zzws
              zwd(ji, jj, jk) = 1._wp - zzwi - zzws
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
      END SELECT
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zwi(ji, jj, 1) = 0._wp
          zwd(ji, jj, 1) = 1._wp - zws(ji, jj, 1)
        END DO
      END DO
    END IF
    IF (ln_drgimp) THEN
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          iku = mbku(ji, jj)
          ze3ua = (1._wp - r_vvl) * e3u_n(ji, jj, iku) + r_vvl * e3u_a(ji, jj, iku)
          zwd(ji, jj, iku) = zwd(ji, jj, iku) - r2dt * 0.5 * (rCdU_bot(ji + 1, jj) + rCdU_bot(ji, jj)) / ze3ua
        END DO
      END DO
      IF (ln_isfcav) THEN
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            iku = miku(ji, jj)
            ze3ua = (1._wp - r_vvl) * e3u_n(ji, jj, iku) + r_vvl * e3u_a(ji, jj, iku)
            zwd(ji, jj, iku) = zwd(ji, jj, iku) - r2dt * 0.5 * (rCdU_top(ji + 1, jj) + rCdU_top(ji, jj)) / ze3ua
          END DO
        END DO
      END IF
    END IF
    !$OMP parallel default(shared), private(ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = 2, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zwd(ji, jj, jk) = zwd(ji, jj, jk) - zwi(ji, jj, jk) * zws(ji, jj, jk - 1) / zwd(ji, jj, jk - 1)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        ze3ua = (1._wp - r_vvl) * e3u_n(ji, jj, 1) + r_vvl * e3u_a(ji, jj, 1)
        ua(ji, jj, 1) = ua(ji, jj, 1) + r2dt * 0.5_wp * (utau_b(ji, jj) + utau(ji, jj)) / (ze3ua * rau0) * umask(ji, jj, 1)
      END DO
    END DO
    !$OMP parallel default(shared), private(ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = 2, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ua(ji, jj, jk) = ua(ji, jj, jk) - zwi(ji, jj, jk) / zwd(ji, jj, jk - 1) * ua(ji, jj, jk - 1)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        ua(ji, jj, jpkm1) = ua(ji, jj, jpkm1) / zwd(ji, jj, jpkm1)
      END DO
    END DO
    !$OMP parallel default(shared), private(ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = jpk - 2, 1, - 1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ua(ji, jj, jk) = (ua(ji, jj, jk) - zws(ji, jj, jk) * ua(ji, jj, jk + 1)) / zwd(ji, jj, jk)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    zdt = r2dt * 0.5
    IF (ln_zad_Aimp) THEN
      SELECT CASE (nldf_dyn)
      CASE (np_lap_i)
        !$OMP parallel default(shared), private(ji,jj,jk,ze3va,zwvi,zwvs,zzwi,zzws)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              ze3va = (1._wp - r_vvl) * e3v_n(ji, jj, jk) + r_vvl * e3v_a(ji, jj, jk)
              zzwi = - zdt * (avm(ji, jj + 1, jk) + avm(ji, jj, jk) + akzv(ji, jj, jk)) / (ze3va * e3vw_n(ji, jj, jk)) * &
&wvmask(ji, jj, jk)
              zzws = - zdt * (avm(ji, jj + 1, jk + 1) + avm(ji, jj, jk + 1) + akzv(ji, jj, jk + 1)) / (ze3va * e3vw_n(ji, jj, jk + &
&1)) * wvmask(ji, jj, jk + 1)
              zWvi = 0.5_wp * (wi(ji, jj, jk) + wi(ji, jj + 1, jk)) * wvmask(ji, jj, jk)
              zWvs = 0.5_wp * (wi(ji, jj, jk + 1) + wi(ji, jj + 1, jk + 1)) * wvmask(ji, jj, jk + 1)
              zwi(ji, jj, jk) = zzwi + zdt * MIN(zWvi, 0._wp)
              zws(ji, jj, jk) = zzws - zdt * MAX(zWvs, 0._wp)
              zwd(ji, jj, jk) = 1._wp - zzwi - zzws - zdt * (- MAX(zWvi, 0._wp) + MIN(zWvs, 0._wp))
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
      CASE DEFAULT
        !$OMP parallel default(shared), private(ji,jj,jk,ze3va,zwvi,zwvs,zzwi,zzws)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              ze3va = (1._wp - r_vvl) * e3v_n(ji, jj, jk) + r_vvl * e3v_a(ji, jj, jk)
              zzwi = - zdt * (avm(ji, jj + 1, jk) + avm(ji, jj, jk)) / (ze3va * e3vw_n(ji, jj, jk)) * wvmask(ji, jj, jk)
              zzws = - zdt * (avm(ji, jj + 1, jk + 1) + avm(ji, jj, jk + 1)) / (ze3va * e3vw_n(ji, jj, jk + 1)) * wvmask(ji, jj, &
&jk + 1)
              zWvi = 0.5_wp * (wi(ji, jj, jk) + wi(ji, jj + 1, jk)) * wvmask(ji, jj, jk)
              zWvs = 0.5_wp * (wi(ji, jj, jk + 1) + wi(ji, jj + 1, jk + 1)) * wvmask(ji, jj, jk + 1)
              zwi(ji, jj, jk) = zzwi + zdt * MIN(zWvi, 0._wp)
              zws(ji, jj, jk) = zzws - zdt * MAX(zWvs, 0._wp)
              zwd(ji, jj, jk) = 1._wp - zzwi - zzws - zdt * (- MAX(zWvi, 0._wp) + MIN(zWvs, 0._wp))
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
      END SELECT
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zwi(ji, jj, 1) = 0._wp
          ze3va = (1._wp - r_vvl) * e3v_n(ji, jj, 1) + r_vvl * e3v_a(ji, jj, 1)
          zzws = - zdt * (avm(ji, jj + 1, 2) + avm(ji, jj, 2)) / (ze3va * e3vw_n(ji, jj, 2)) * wvmask(ji, jj, 2)
          zWvs = 0.5_wp * (wi(ji, jj, 2) + wi(ji, jj + 1, 2))
          zws(ji, jj, 1) = zzws - zdt * MAX(zWvs, 0._wp)
          zwd(ji, jj, 1) = 1._wp - zzws - zdt * (MIN(zWvs, 0._wp))
        END DO
      END DO
    ELSE
      SELECT CASE (nldf_dyn)
      CASE (np_lap_i)
        !$OMP parallel default(shared), private(ji,jj,jk,ze3va,zzwi,zzws)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              ze3va = (1._wp - r_vvl) * e3v_n(ji, jj, jk) + r_vvl * e3v_a(ji, jj, jk)
              zzwi = - zdt * (avm(ji, jj + 1, jk) + avm(ji, jj, jk) + akzv(ji, jj, jk)) / (ze3va * e3vw_n(ji, jj, jk)) * &
&wvmask(ji, jj, jk)
              zzws = - zdt * (avm(ji, jj + 1, jk + 1) + avm(ji, jj, jk + 1) + akzv(ji, jj, jk + 1)) / (ze3va * e3vw_n(ji, jj, jk + &
&1)) * wvmask(ji, jj, jk + 1)
              zwi(ji, jj, jk) = zzwi
              zws(ji, jj, jk) = zzws
              zwd(ji, jj, jk) = 1._wp - zzwi - zzws
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
      CASE DEFAULT
        !$OMP parallel default(shared), private(ji,jj,jk,ze3va,zzwi,zzws)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              ze3va = (1._wp - r_vvl) * e3v_n(ji, jj, jk) + r_vvl * e3v_a(ji, jj, jk)
              zzwi = - zdt * (avm(ji, jj + 1, jk) + avm(ji, jj, jk)) / (ze3va * e3vw_n(ji, jj, jk)) * wvmask(ji, jj, jk)
              zzws = - zdt * (avm(ji, jj + 1, jk + 1) + avm(ji, jj, jk + 1)) / (ze3va * e3vw_n(ji, jj, jk + 1)) * wvmask(ji, jj, &
&jk + 1)
              zwi(ji, jj, jk) = zzwi
              zws(ji, jj, jk) = zzws
              zwd(ji, jj, jk) = 1._wp - zzwi - zzws
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
      END SELECT
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zwi(ji, jj, 1) = 0._wp
          zwd(ji, jj, 1) = 1._wp - zws(ji, jj, 1)
        END DO
      END DO
    END IF
    IF (ln_drgimp) THEN
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ikv = mbkv(ji, jj)
          ze3va = (1._wp - r_vvl) * e3v_n(ji, jj, ikv) + r_vvl * e3v_a(ji, jj, ikv)
          zwd(ji, jj, ikv) = zwd(ji, jj, ikv) - r2dt * 0.5 * (rCdU_bot(ji, jj + 1) + rCdU_bot(ji, jj)) / ze3va
        END DO
      END DO
      IF (ln_isfcav) THEN
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            ikv = mikv(ji, jj)
            ze3va = (1._wp - r_vvl) * e3v_n(ji, jj, ikv) + r_vvl * e3v_a(ji, jj, ikv)
            zwd(ji, jj, iku) = zwd(ji, jj, iku) - r2dt * 0.5 * (rCdU_top(ji + 1, jj) + rCdU_top(ji, jj)) / ze3va
          END DO
        END DO
      END IF
    END IF
    !$OMP parallel default(shared), private(ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = 2, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zwd(ji, jj, jk) = zwd(ji, jj, jk) - zwi(ji, jj, jk) * zws(ji, jj, jk - 1) / zwd(ji, jj, jk - 1)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        ze3va = (1._wp - r_vvl) * e3v_n(ji, jj, 1) + r_vvl * e3v_a(ji, jj, 1)
        va(ji, jj, 1) = va(ji, jj, 1) + r2dt * 0.5_wp * (vtau_b(ji, jj) + vtau(ji, jj)) / (ze3va * rau0) * vmask(ji, jj, 1)
      END DO
    END DO
    !$OMP parallel default(shared), private(ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = 2, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          va(ji, jj, jk) = va(ji, jj, jk) - zwi(ji, jj, jk) / zwd(ji, jj, jk - 1) * va(ji, jj, jk - 1)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        va(ji, jj, jpkm1) = va(ji, jj, jpkm1) / zwd(ji, jj, jpkm1)
      END DO
    END DO
    !$OMP parallel default(shared), private(ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = jpk - 2, 1, - 1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          va(ji, jj, jk) = (va(ji, jj, jk) - zws(ji, jj, jk) * va(ji, jj, jk + 1)) / zwd(ji, jj, jk)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    IF (l_trddyn) THEN
      ztrdu(:, :, :) = (ua(:, :, :) - ub(:, :, :)) / r2dt - ztrdu(:, :, :)
      ztrdv(:, :, :) = (va(:, :, :) - vb(:, :, :)) / r2dt - ztrdv(:, :, :)
      CALL trd_dyn(ztrdu, ztrdv, jpdyn_zdf, kt)
      DEALLOCATE(ztrdu, ztrdv)
    END IF
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = ua, clinfo1 = ' zdf  - Ua: ', mask1 = umask, tab3d_2 = va, clinfo2 = ' Va: ', mask2 = &
&vmask, clinfo3 = 'dyn')
    IF (ln_timing) CALL timing_stop('dyn_zdf')
  END SUBROUTINE dyn_zdf
END MODULE dynzdf