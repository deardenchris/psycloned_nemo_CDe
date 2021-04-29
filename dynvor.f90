MODULE dynvor
  USE oce
  USE dom_oce
  USE dommsk
  USE dynadv
  USE trd_oce
  USE trddyn
  USE sbcwave
  USE sbc_oce, ONLY: ln_stcor
  USE lbclnk
  USE prtctl
  USE in_out_manager
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_vor
  PUBLIC :: dyn_vor_init
  LOGICAL, PUBLIC :: ln_dynvor_ens
  LOGICAL, PUBLIC :: ln_dynvor_ene
  LOGICAL, PUBLIC :: ln_dynvor_enT
  LOGICAL, PUBLIC :: ln_dynvor_eeT
  LOGICAL, PUBLIC :: ln_dynvor_een
  INTEGER, PUBLIC :: nn_een_e3f
  LOGICAL, PUBLIC :: ln_dynvor_mix
  LOGICAL, PUBLIC :: ln_dynvor_msk
  INTEGER, PUBLIC :: nvor_scheme
  INTEGER, PUBLIC, PARAMETER :: np_ENS = 0
  INTEGER, PUBLIC, PARAMETER :: np_ENE = 1
  INTEGER, PUBLIC, PARAMETER :: np_ENT = 2
  INTEGER, PUBLIC, PARAMETER :: np_EET = 3
  INTEGER, PUBLIC, PARAMETER :: np_EEN = 4
  INTEGER, PUBLIC, PARAMETER :: np_MIX = 5
  INTEGER :: ncor, nrvm, ntot
  INTEGER, PUBLIC, PARAMETER :: np_COR = 1
  INTEGER, PUBLIC, PARAMETER :: np_RVO = 2
  INTEGER, PUBLIC, PARAMETER :: np_MET = 3
  INTEGER, PUBLIC, PARAMETER :: np_CRV = 4
  INTEGER, PUBLIC, PARAMETER :: np_CME = 5
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: di_e2u_2
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: dj_e1v_2
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: di_e2v_2e1e2f
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: dj_e1u_2e1e2f
  REAL(KIND = wp) :: r1_4 = 0.250_wp
  REAL(KIND = wp) :: r1_8 = 0.125_wp
  REAL(KIND = wp) :: r1_12 = 1._wp / 12._wp
  CONTAINS
  SUBROUTINE dyn_vor(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ztrdu, ztrdv
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    CALL profile_psy_data0 % PreStart('dyn_vor', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('dyn_vor')
    CALL profile_psy_data0 % PostEnd
    IF (l_trddyn) THEN
      ALLOCATE(ztrdu(jpi, jpj, jpk), ztrdv(jpi, jpj, jpk))
      !$ACC KERNELS
      ztrdu(:, :, :) = ua(:, :, :)
      ztrdv(:, :, :) = va(:, :, :)
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('dyn_vor', 'r1', 0, 0)
      SELECT CASE (nvor_scheme)
      CASE (np_ens)
        CALL vor_ens(kt, ncor, un, vn, ua, va)
        IF (ln_stcor) CALL vor_ens(kt, ncor, usd, vsd, ua, va)
      CASE (np_ene, np_mix)
        CALL vor_ene(kt, ncor, un, vn, ua, va)
        IF (ln_stcor) CALL vor_ene(kt, ncor, usd, vsd, ua, va)
      CASE (np_ent)
        CALL vor_ent(kt, ncor, un, vn, ua, va)
        IF (ln_stcor) CALL vor_enT(kt, ncor, usd, vsd, ua, va)
      CASE (np_eet)
        CALL vor_eet(kt, ncor, un, vn, ua, va)
        IF (ln_stcor) CALL vor_eeT(kt, ncor, usd, vsd, ua, va)
      CASE (np_een)
        CALL vor_een(kt, ncor, un, vn, ua, va)
        IF (ln_stcor) CALL vor_een(kt, ncor, usd, vsd, ua, va)
      END SELECT
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      ztrdu(:, :, :) = ua(:, :, :) - ztrdu(:, :, :)
      ztrdv(:, :, :) = va(:, :, :) - ztrdv(:, :, :)
      !$ACC END KERNELS
      CALL profile_psy_data2 % PreStart('dyn_vor', 'r2', 0, 0)
      CALL trd_dyn(ztrdu, ztrdv, jpdyn_pvo, kt)
      CALL profile_psy_data2 % PostEnd
      IF (n_dynadv /= np_LIN_dyn) THEN
        !$ACC KERNELS
        ztrdu(:, :, :) = ua(:, :, :)
        ztrdv(:, :, :) = va(:, :, :)
        !$ACC END KERNELS
        CALL profile_psy_data3 % PreStart('dyn_vor', 'r3', 0, 0)
        SELECT CASE (nvor_scheme)
        CASE (np_ent)
          CALL vor_ent(kt, nrvm, un, vn, ua, va)
        CASE (np_eet)
          CALL vor_eet(kt, nrvm, un, vn, ua, va)
        CASE (np_ene)
          CALL vor_ene(kt, nrvm, un, vn, ua, va)
        CASE (np_ens, np_mix)
          CALL vor_ens(kt, nrvm, un, vn, ua, va)
        CASE (np_een)
          CALL vor_een(kt, nrvm, un, vn, ua, va)
        END SELECT
        CALL profile_psy_data3 % PostEnd
        !$ACC KERNELS
        ztrdu(:, :, :) = ua(:, :, :) - ztrdu(:, :, :)
        ztrdv(:, :, :) = va(:, :, :) - ztrdv(:, :, :)
        !$ACC END KERNELS
        CALL profile_psy_data4 % PreStart('dyn_vor', 'r4', 0, 0)
        CALL trd_dyn(ztrdu, ztrdv, jpdyn_rvo, kt)
        CALL profile_psy_data4 % PostEnd
      END IF
      DEALLOCATE(ztrdu, ztrdv)
    ELSE
      CALL profile_psy_data5 % PreStart('dyn_vor', 'r5', 0, 0)
      SELECT CASE (nvor_scheme)
      CASE (np_ENT)
        CALL vor_enT(kt, ntot, un, vn, ua, va)
        IF (ln_stcor) CALL vor_enT(kt, ncor, usd, vsd, ua, va)
      CASE (np_EET)
        CALL vor_eeT(kt, ntot, un, vn, ua, va)
        IF (ln_stcor) CALL vor_eeT(kt, ncor, usd, vsd, ua, va)
      CASE (np_ENE)
        CALL vor_ene(kt, ntot, un, vn, ua, va)
        IF (ln_stcor) CALL vor_ene(kt, ncor, usd, vsd, ua, va)
      CASE (np_ENS)
        CALL vor_ens(kt, ntot, un, vn, ua, va)
        IF (ln_stcor) CALL vor_ens(kt, ncor, usd, vsd, ua, va)
      CASE (np_MIX)
        CALL vor_ens(kt, nrvm, un, vn, ua, va)
        CALL vor_ene(kt, ncor, un, vn, ua, va)
        IF (ln_stcor) CALL vor_ene(kt, ncor, usd, vsd, ua, va)
      CASE (np_EEN)
        CALL vor_een(kt, ntot, un, vn, ua, va)
        IF (ln_stcor) CALL vor_een(kt, ncor, usd, vsd, ua, va)
      END SELECT
      CALL profile_psy_data5 % PostEnd
    END IF
    CALL profile_psy_data6 % PreStart('dyn_vor', 'r6', 0, 0)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = ua, clinfo1 = ' vor  - Ua: ', mask1 = umask, tab3d_2 = va, clinfo2 = ' Va: ', mask2 = vmask, clinfo3 = 'dyn')
    IF (ln_timing) CALL timing_stop('dyn_vor')
    CALL profile_psy_data6 % PostEnd
  END SUBROUTINE dyn_vor
  SUBROUTINE vor_enT(kt, kvor, pu, pv, pu_rhs, pv_rhs)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kvor
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pu, pv
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pu_rhs, pv_rhs
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zx1, zy1, zx2, zy2
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwx, zwy, zwt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwz
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('vor_ent', 'r0', 0, 0)
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'dyn:vor_enT : vorticity term: t-point energy conserving scheme'
      IF (lwp) WRITE(numout, *) '~~~~~~~~~~~'
    END IF
    CALL profile_psy_data0 % PostEnd
    SELECT CASE (kvor)
    CASE (np_RVO)
      DO jk = 1, jpkm1
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj, jk) = (e2v(ji + 1, jj) * pv(ji + 1, jj, jk) - e2v(ji, jj) * pv(ji, jj, jk) - e1u(ji, jj + 1) * pu(ji, jj + 1, jk) + e1u(ji, jj) * pu(ji, jj, jk)) * r1_e1e2f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
        IF (ln_dynvor_msk) THEN
          !$ACC KERNELS
          !$ACC loop independent collapse(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              zwz(ji, jj, jk) = zwz(ji, jj, jk) * fmask(ji, jj, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END IF
      END DO
      CALL profile_psy_data1 % PreStart('vor_ent', 'r1', 0, 0)
      CALL lbc_lnk('dynvor', zwz, 'F', 1.)
      CALL profile_psy_data1 % PostEnd
    CASE (np_CRV)
      DO jk = 1, jpkm1
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj, jk) = (e2v(ji + 1, jj) * pv(ji + 1, jj, jk) - e2v(ji, jj) * pv(ji, jj, jk) - e1u(ji, jj + 1) * pu(ji, jj + 1, jk) + e1u(ji, jj) * pu(ji, jj, jk)) * r1_e1e2f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
        IF (ln_dynvor_msk) THEN
          !$ACC KERNELS
          !$ACC loop independent collapse(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              zwz(ji, jj, jk) = zwz(ji, jj, jk) * fmask(ji, jj, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END IF
      END DO
      CALL profile_psy_data2 % PreStart('vor_ent', 'r2', 0, 0)
      CALL lbc_lnk('dynvor', zwz, 'F', 1.)
      CALL profile_psy_data2 % PostEnd
    END SELECT
    DO jk = 1, jpkm1
      SELECT CASE (kvor)
      CASE (np_COR)
        !$ACC KERNELS
        zwt(:, :) = ff_t(:, :) * e1e2t(:, :) * e3t_n(:, :, jk)
        !$ACC END KERNELS
      CASE (np_RVO)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpj
          DO ji = 2, jpi
            zwt(ji, jj) = r1_4 * (zwz(ji - 1, jj, jk) + zwz(ji, jj, jk) + zwz(ji - 1, jj - 1, jk) + zwz(ji, jj - 1, jk)) * e1e2t(ji, jj) * e3t_n(ji, jj, jk)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_MET)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpj
          DO ji = 2, jpi
            zwt(ji, jj) = ((pv(ji, jj, jk) + pv(ji, jj - 1, jk)) * di_e2u_2(ji, jj) - (pu(ji, jj, jk) + pu(ji - 1, jj, jk)) * dj_e1v_2(ji, jj)) * e3t_n(ji, jj, jk)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_CRV)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpj
          DO ji = 2, jpi
            zwt(ji, jj) = (ff_t(ji, jj) + r1_4 * (zwz(ji - 1, jj, jk) + zwz(ji, jj, jk) + zwz(ji - 1, jj - 1, jk) + zwz(ji, jj - 1, jk))) * e1e2t(ji, jj) * e3t_n(ji, jj, jk)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_CME)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpj
          DO ji = 2, jpi
            zwt(ji, jj) = (ff_t(ji, jj) * e1e2t(ji, jj) + (pv(ji, jj, jk) + pv(ji, jj - 1, jk)) * di_e2u_2(ji, jj) - (pu(ji, jj, jk) + pu(ji - 1, jj, jk)) * dj_e1v_2(ji, jj)) * e3t_n(ji, jj, jk)
          END DO
        END DO
        !$ACC END KERNELS
      CASE DEFAULT
        CALL profile_psy_data3 % PreStart('vor_ent', 'r3', 0, 0)
        CALL ctl_stop('STOP', 'dyn_vor: wrong value for kvor')
        CALL profile_psy_data3 % PostEnd
      END SELECT
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          pu_rhs(ji, jj, jk) = pu_rhs(ji, jj, jk) + r1_4 * r1_e1e2u(ji, jj) / e3u_n(ji, jj, jk) * (zwt(ji + 1, jj) * (pv(ji + 1, jj, jk) + pv(ji + 1, jj - 1, jk)) + zwt(ji, jj) * (pv(ji, jj, jk) + pv(ji, jj - 1, jk)))
          pv_rhs(ji, jj, jk) = pv_rhs(ji, jj, jk) - r1_4 * r1_e1e2v(ji, jj) / e3v_n(ji, jj, jk) * (zwt(ji, jj + 1) * (pu(ji, jj + 1, jk) + pu(ji - 1, jj + 1, jk)) + zwt(ji, jj) * (pu(ji, jj, jk) + pu(ji - 1, jj, jk)))
        END DO
      END DO
      !$ACC END KERNELS
    END DO
  END SUBROUTINE vor_enT
  SUBROUTINE vor_ene(kt, kvor, pun, pvn, pua, pva)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kvor
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pun, pvn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pua, pva
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zx1, zy1, zx2, zy2
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwx, zwy, zwz
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('vor_ene', 'r0', 0, 0)
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'dyn:vor_ene : vorticity term: energy conserving scheme'
      IF (lwp) WRITE(numout, *) '~~~~~~~~~~~'
    END IF
    CALL profile_psy_data0 % PostEnd
    DO jk = 1, jpkm1
      SELECT CASE (kvor)
      CASE (np_COR)
        !$ACC KERNELS
        zwz(:, :) = ff_f(:, :)
        !$ACC END KERNELS
      CASE (np_RVO)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj) = (e2v(ji + 1, jj) * pvn(ji + 1, jj, jk) - e2v(ji, jj) * pvn(ji, jj, jk) - e1u(ji, jj + 1) * pun(ji, jj + 1, jk) + e1u(ji, jj) * pun(ji, jj, jk)) * r1_e1e2f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_MET)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj) = (pvn(ji + 1, jj, jk) + pvn(ji, jj, jk)) * di_e2v_2e1e2f(ji, jj) - (pun(ji, jj + 1, jk) + pun(ji, jj, jk)) * dj_e1u_2e1e2f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_CRV)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj) = ff_f(ji, jj) + (e2v(ji + 1, jj) * pvn(ji + 1, jj, jk) - e2v(ji, jj) * pvn(ji, jj, jk) - e1u(ji, jj + 1) * pun(ji, jj + 1, jk) + e1u(ji, jj) * pun(ji, jj, jk)) * r1_e1e2f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_CME)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj) = ff_f(ji, jj) + (pvn(ji + 1, jj, jk) + pvn(ji, jj, jk)) * di_e2v_2e1e2f(ji, jj) - (pun(ji, jj + 1, jk) + pun(ji, jj, jk)) * dj_e1u_2e1e2f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE DEFAULT
        CALL profile_psy_data1 % PreStart('vor_ene', 'r1', 0, 0)
        CALL ctl_stop('STOP', 'dyn_vor: wrong value for kvor')
        CALL profile_psy_data1 % PostEnd
      END SELECT
      IF (ln_dynvor_msk) THEN
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj) = zwz(ji, jj) * fmask(ji, jj, jk)
          END DO
        END DO
        !$ACC END KERNELS
      END IF
      IF (ln_sco) THEN
        !$ACC KERNELS
        zwz(:, :) = zwz(:, :) / e3f_n(:, :, jk)
        zwx(:, :) = e2u(:, :) * e3u_n(:, :, jk) * pun(:, :, jk)
        zwy(:, :) = e1v(:, :) * e3v_n(:, :, jk) * pvn(:, :, jk)
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        zwx(:, :) = e2u(:, :) * pun(:, :, jk)
        zwy(:, :) = e1v(:, :) * pvn(:, :, jk)
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zy1 = zwy(ji, jj - 1) + zwy(ji + 1, jj - 1)
          zy2 = zwy(ji, jj) + zwy(ji + 1, jj)
          zx1 = zwx(ji - 1, jj) + zwx(ji - 1, jj + 1)
          zx2 = zwx(ji, jj) + zwx(ji, jj + 1)
          pua(ji, jj, jk) = pua(ji, jj, jk) + r1_4 * r1_e1u(ji, jj) * (zwz(ji, jj - 1) * zy1 + zwz(ji, jj) * zy2)
          pva(ji, jj, jk) = pva(ji, jj, jk) - r1_4 * r1_e2v(ji, jj) * (zwz(ji - 1, jj) * zx1 + zwz(ji, jj) * zx2)
        END DO
      END DO
      !$ACC END KERNELS
    END DO
  END SUBROUTINE vor_ene
  SUBROUTINE vor_ens(kt, kvor, pun, pvn, pua, pva)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kvor
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pun, pvn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pua, pva
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zuav, zvau
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwx, zwy, zwz, zww
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('vor_ens', 'r0', 0, 0)
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'dyn:vor_ens : vorticity term: enstrophy conserving scheme'
      IF (lwp) WRITE(numout, *) '~~~~~~~~~~~'
    END IF
    CALL profile_psy_data0 % PostEnd
    DO jk = 1, jpkm1
      SELECT CASE (kvor)
      CASE (np_COR)
        !$ACC KERNELS
        zwz(:, :) = ff_f(:, :)
        !$ACC END KERNELS
      CASE (np_RVO)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj) = (e2v(ji + 1, jj) * pvn(ji + 1, jj, jk) - e2v(ji, jj) * pvn(ji, jj, jk) - e1u(ji, jj + 1) * pun(ji, jj + 1, jk) + e1u(ji, jj) * pun(ji, jj, jk)) * r1_e1e2f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_MET)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj) = (pvn(ji + 1, jj, jk) + pvn(ji, jj, jk)) * di_e2v_2e1e2f(ji, jj) - (pun(ji, jj + 1, jk) + pun(ji, jj, jk)) * dj_e1u_2e1e2f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_CRV)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj) = ff_f(ji, jj) + (e2v(ji + 1, jj) * pvn(ji + 1, jj, jk) - e2v(ji, jj) * pvn(ji, jj, jk) - e1u(ji, jj + 1) * pun(ji, jj + 1, jk) + e1u(ji, jj) * pun(ji, jj, jk)) * r1_e1e2f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_CME)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj) = ff_f(ji, jj) + (pvn(ji + 1, jj, jk) + pvn(ji, jj, jk)) * di_e2v_2e1e2f(ji, jj) - (pun(ji, jj + 1, jk) + pun(ji, jj, jk)) * dj_e1u_2e1e2f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE DEFAULT
        CALL profile_psy_data1 % PreStart('vor_ens', 'r1', 0, 0)
        CALL ctl_stop('STOP', 'dyn_vor: wrong value for kvor')
        CALL profile_psy_data1 % PostEnd
      END SELECT
      IF (ln_dynvor_msk) THEN
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj) = zwz(ji, jj) * fmask(ji, jj, jk)
          END DO
        END DO
        !$ACC END KERNELS
      END IF
      IF (ln_sco) THEN
        !$ACC KERNELS
        zwz(:, :) = zwz(:, :) / e3f_n(:, :, jk)
        zwx(:, :) = e2u(:, :) * e3u_n(:, :, jk) * pun(:, :, jk)
        zwy(:, :) = e1v(:, :) * e3v_n(:, :, jk) * pvn(:, :, jk)
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        zwx(:, :) = e2u(:, :) * pun(:, :, jk)
        zwy(:, :) = e1v(:, :) * pvn(:, :, jk)
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zuav = r1_8 * r1_e1u(ji, jj) * (zwy(ji, jj - 1) + zwy(ji + 1, jj - 1) + zwy(ji, jj) + zwy(ji + 1, jj))
          zvau = - r1_8 * r1_e2v(ji, jj) * (zwx(ji - 1, jj) + zwx(ji - 1, jj + 1) + zwx(ji, jj) + zwx(ji, jj + 1))
          pua(ji, jj, jk) = pua(ji, jj, jk) + zuav * (zwz(ji, jj - 1) + zwz(ji, jj))
          pva(ji, jj, jk) = pva(ji, jj, jk) + zvau * (zwz(ji - 1, jj) + zwz(ji, jj))
        END DO
      END DO
      !$ACC END KERNELS
    END DO
  END SUBROUTINE vor_ens
  SUBROUTINE vor_een(kt, kvor, pun, pvn, pua, pva)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kvor
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pun, pvn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pua, pva
    INTEGER :: ji, jj, jk
    INTEGER :: ierr
    REAL(KIND = wp) :: zua, zva
    REAL(KIND = wp) :: zmsk, ze3f
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwx, zwy, z1_e3f
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztnw, ztne, ztsw, ztse
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwz
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('vor_een', 'r0', 0, 0)
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'dyn:vor_een : vorticity term: energy and enstrophy conserving scheme'
      IF (lwp) WRITE(numout, *) '~~~~~~~~~~~'
    END IF
    CALL profile_psy_data0 % PostEnd
    DO jk = 1, jpkm1
      !$ACC KERNELS
      SELECT CASE (nn_een_e3f)
      CASE (0)
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            ze3f = (e3t_n(ji, jj + 1, jk) * tmask(ji, jj + 1, jk) + e3t_n(ji + 1, jj + 1, jk) * tmask(ji + 1, jj + 1, jk) + e3t_n(ji, jj, jk) * tmask(ji, jj, jk) + e3t_n(ji + 1, jj, jk) * tmask(ji + 1, jj, jk))
            IF (ze3f /= 0._wp) THEN
              z1_e3f(ji, jj) = 4._wp / ze3f
            ELSE
              z1_e3f(ji, jj) = 0._wp
            END IF
          END DO
        END DO
      CASE (1)
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            ze3f = (e3t_n(ji, jj + 1, jk) * tmask(ji, jj + 1, jk) + e3t_n(ji + 1, jj + 1, jk) * tmask(ji + 1, jj + 1, jk) + e3t_n(ji, jj, jk) * tmask(ji, jj, jk) + e3t_n(ji + 1, jj, jk) * tmask(ji + 1, jj, jk))
            zmsk = (tmask(ji, jj + 1, jk) + tmask(ji + 1, jj + 1, jk) + tmask(ji, jj, jk) + tmask(ji + 1, jj, jk))
            IF (ze3f /= 0._wp) THEN
              z1_e3f(ji, jj) = zmsk / ze3f
            ELSE
              z1_e3f(ji, jj) = 0._wp
            END IF
          END DO
        END DO
      END SELECT
      !$ACC END KERNELS
      SELECT CASE (kvor)
      CASE (np_COR)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj, jk) = ff_f(ji, jj) * z1_e3f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_RVO)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj, jk) = (e2v(ji + 1, jj) * pvn(ji + 1, jj, jk) - e2v(ji, jj) * pvn(ji, jj, jk) - e1u(ji, jj + 1) * pun(ji, jj + 1, jk) + e1u(ji, jj) * pun(ji, jj, jk)) * r1_e1e2f(ji, jj) * z1_e3f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_MET)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj, jk) = ((pvn(ji + 1, jj, jk) + pvn(ji, jj, jk)) * di_e2v_2e1e2f(ji, jj) - (pun(ji, jj + 1, jk) + pun(ji, jj, jk)) * dj_e1u_2e1e2f(ji, jj)) * z1_e3f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_CRV)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj, jk) = (ff_f(ji, jj) + (e2v(ji + 1, jj) * pvn(ji + 1, jj, jk) - e2v(ji, jj) * pvn(ji, jj, jk) - e1u(ji, jj + 1) * pun(ji, jj + 1, jk) + e1u(ji, jj) * pun(ji, jj, jk)) * r1_e1e2f(ji, jj)) * z1_e3f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_CME)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj, jk) = (ff_f(ji, jj) + (pvn(ji + 1, jj, jk) + pvn(ji, jj, jk)) * di_e2v_2e1e2f(ji, jj) - (pun(ji, jj + 1, jk) + pun(ji, jj, jk)) * dj_e1u_2e1e2f(ji, jj)) * z1_e3f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE DEFAULT
        CALL profile_psy_data1 % PreStart('vor_een', 'r1', 0, 0)
        CALL ctl_stop('STOP', 'dyn_vor: wrong value for kvor')
        CALL profile_psy_data1 % PostEnd
      END SELECT
      IF (ln_dynvor_msk) THEN
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj, jk) = zwz(ji, jj, jk) * fmask(ji, jj, jk)
          END DO
        END DO
        !$ACC END KERNELS
      END IF
    END DO
    CALL profile_psy_data2 % PreStart('vor_een', 'r2', 0, 0)
    CALL lbc_lnk('dynvor', zwz, 'F', 1.)
    CALL profile_psy_data2 % PostEnd
    DO jk = 1, jpkm1
      !$ACC KERNELS
      zwx(:, :) = e2u(:, :) * e3u_n(:, :, jk) * pun(:, :, jk)
      zwy(:, :) = e1v(:, :) * e3v_n(:, :, jk) * pvn(:, :, jk)
      jj = 2
      ztne(1, :) = 0
      ztnw(1, :) = 0
      ztse(1, :) = 0
      ztsw(1, :) = 0
      DO ji = 2, jpi
        ztne(ji, jj) = zwz(ji - 1, jj, jk) + zwz(ji, jj, jk) + zwz(ji, jj - 1, jk)
        ztnw(ji, jj) = zwz(ji - 1, jj - 1, jk) + zwz(ji - 1, jj, jk) + zwz(ji, jj, jk)
        ztse(ji, jj) = zwz(ji, jj, jk) + zwz(ji, jj - 1, jk) + zwz(ji - 1, jj - 1, jk)
        ztsw(ji, jj) = zwz(ji, jj - 1, jk) + zwz(ji - 1, jj - 1, jk) + zwz(ji - 1, jj, jk)
      END DO
      !$ACC loop independent collapse(2)
      DO jj = 3, jpj
        DO ji = 2, jpi
          ztne(ji, jj) = zwz(ji - 1, jj, jk) + zwz(ji, jj, jk) + zwz(ji, jj - 1, jk)
          ztnw(ji, jj) = zwz(ji - 1, jj - 1, jk) + zwz(ji - 1, jj, jk) + zwz(ji, jj, jk)
          ztse(ji, jj) = zwz(ji, jj, jk) + zwz(ji, jj - 1, jk) + zwz(ji - 1, jj - 1, jk)
          ztsw(ji, jj) = zwz(ji, jj - 1, jk) + zwz(ji - 1, jj - 1, jk) + zwz(ji - 1, jj, jk)
        END DO
      END DO
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zua = + r1_12 * r1_e1u(ji, jj) * (ztne(ji, jj) * zwy(ji, jj) + ztnw(ji + 1, jj) * zwy(ji + 1, jj) + ztse(ji, jj) * zwy(ji, jj - 1) + ztsw(ji + 1, jj) * zwy(ji + 1, jj - 1))
          zva = - r1_12 * r1_e2v(ji, jj) * (ztsw(ji, jj + 1) * zwx(ji - 1, jj + 1) + ztse(ji, jj + 1) * zwx(ji, jj + 1) + ztnw(ji, jj) * zwx(ji - 1, jj) + ztne(ji, jj) * zwx(ji, jj))
          pua(ji, jj, jk) = pua(ji, jj, jk) + zua
          pva(ji, jj, jk) = pva(ji, jj, jk) + zva
        END DO
      END DO
      !$ACC END KERNELS
    END DO
  END SUBROUTINE vor_een
  SUBROUTINE vor_eeT(kt, kvor, pun, pvn, pua, pva)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kvor
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pun, pvn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pua, pva
    INTEGER :: ji, jj, jk
    INTEGER :: ierr
    REAL(KIND = wp) :: zua, zva
    REAL(KIND = wp) :: zmsk, z1_e3t
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwx, zwy
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztnw, ztne, ztsw, ztse
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwz
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('vor_eet', 'r0', 0, 0)
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'dyn:vor_een : vorticity term: energy and enstrophy conserving scheme'
      IF (lwp) WRITE(numout, *) '~~~~~~~~~~~'
    END IF
    CALL profile_psy_data0 % PostEnd
    DO jk = 1, jpkm1
      SELECT CASE (kvor)
      CASE (np_COR)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj, jk) = ff_f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_RVO)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj, jk) = (e2v(ji + 1, jj) * pvn(ji + 1, jj, jk) - e2v(ji, jj) * pvn(ji, jj, jk) - e1u(ji, jj + 1) * pun(ji, jj + 1, jk) + e1u(ji, jj) * pun(ji, jj, jk)) * r1_e1e2f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_MET)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj, jk) = (pvn(ji + 1, jj, jk) + pvn(ji, jj, jk)) * di_e2v_2e1e2f(ji, jj) - (pun(ji, jj + 1, jk) + pun(ji, jj, jk)) * dj_e1u_2e1e2f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_CRV)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj, jk) = (ff_f(ji, jj) + (e2v(ji + 1, jj) * pvn(ji + 1, jj, jk) - e2v(ji, jj) * pvn(ji, jj, jk) - e1u(ji, jj + 1) * pun(ji, jj + 1, jk) + e1u(ji, jj) * pun(ji, jj, jk)) * r1_e1e2f(ji, jj))
          END DO
        END DO
        !$ACC END KERNELS
      CASE (np_CME)
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj, jk) = ff_f(ji, jj) + (pvn(ji + 1, jj, jk) + pvn(ji, jj, jk)) * di_e2v_2e1e2f(ji, jj) - (pun(ji, jj + 1, jk) + pun(ji, jj, jk)) * dj_e1u_2e1e2f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
      CASE DEFAULT
        CALL profile_psy_data1 % PreStart('vor_eet', 'r1', 0, 0)
        CALL ctl_stop('STOP', 'dyn_vor: wrong value for kvor')
        CALL profile_psy_data1 % PostEnd
      END SELECT
      IF (ln_dynvor_msk) THEN
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwz(ji, jj, jk) = zwz(ji, jj, jk) * fmask(ji, jj, jk)
          END DO
        END DO
        !$ACC END KERNELS
      END IF
    END DO
    CALL profile_psy_data2 % PreStart('vor_eet', 'r2', 0, 0)
    CALL lbc_lnk('dynvor', zwz, 'F', 1.)
    CALL profile_psy_data2 % PostEnd
    DO jk = 1, jpkm1
      !$ACC KERNELS
      zwx(:, :) = e2u(:, :) * e3u_n(:, :, jk) * pun(:, :, jk)
      zwy(:, :) = e1v(:, :) * e3v_n(:, :, jk) * pvn(:, :, jk)
      jj = 2
      ztne(1, :) = 0
      ztnw(1, :) = 0
      ztse(1, :) = 0
      ztsw(1, :) = 0
      DO ji = 2, jpi
        z1_e3t = 1._wp / e3t_n(ji, jj, jk)
        ztne(ji, jj) = (zwz(ji - 1, jj, jk) + zwz(ji, jj, jk) + zwz(ji, jj - 1, jk)) * z1_e3t
        ztnw(ji, jj) = (zwz(ji - 1, jj - 1, jk) + zwz(ji - 1, jj, jk) + zwz(ji, jj, jk)) * z1_e3t
        ztse(ji, jj) = (zwz(ji, jj, jk) + zwz(ji, jj - 1, jk) + zwz(ji - 1, jj - 1, jk)) * z1_e3t
        ztsw(ji, jj) = (zwz(ji, jj - 1, jk) + zwz(ji - 1, jj - 1, jk) + zwz(ji - 1, jj, jk)) * z1_e3t
      END DO
      !$ACC loop independent collapse(2)
      DO jj = 3, jpj
        DO ji = 2, jpi
          z1_e3t = 1._wp / e3t_n(ji, jj, jk)
          ztne(ji, jj) = (zwz(ji - 1, jj, jk) + zwz(ji, jj, jk) + zwz(ji, jj - 1, jk)) * z1_e3t
          ztnw(ji, jj) = (zwz(ji - 1, jj - 1, jk) + zwz(ji - 1, jj, jk) + zwz(ji, jj, jk)) * z1_e3t
          ztse(ji, jj) = (zwz(ji, jj, jk) + zwz(ji, jj - 1, jk) + zwz(ji - 1, jj - 1, jk)) * z1_e3t
          ztsw(ji, jj) = (zwz(ji, jj - 1, jk) + zwz(ji - 1, jj - 1, jk) + zwz(ji - 1, jj, jk)) * z1_e3t
        END DO
      END DO
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zua = + r1_12 * r1_e1u(ji, jj) * (ztne(ji, jj) * zwy(ji, jj) + ztnw(ji + 1, jj) * zwy(ji + 1, jj) + ztse(ji, jj) * zwy(ji, jj - 1) + ztsw(ji + 1, jj) * zwy(ji + 1, jj - 1))
          zva = - r1_12 * r1_e2v(ji, jj) * (ztsw(ji, jj + 1) * zwx(ji - 1, jj + 1) + ztse(ji, jj + 1) * zwx(ji, jj + 1) + ztnw(ji, jj) * zwx(ji - 1, jj) + ztne(ji, jj) * zwx(ji, jj))
          pua(ji, jj, jk) = pua(ji, jj, jk) + zua
          pva(ji, jj, jk) = pva(ji, jj, jk) + zva
        END DO
      END DO
      !$ACC END KERNELS
    END DO
  END SUBROUTINE vor_eeT
  SUBROUTINE dyn_vor_init
    INTEGER :: ji, jj, jk
    INTEGER :: ioptio, ios
    NAMELIST /namdyn_vor/ ln_dynvor_ens, ln_dynvor_ene, ln_dynvor_enT, ln_dynvor_eeT, ln_dynvor_een, nn_een_e3f, ln_dynvor_mix, ln_dynvor_msk
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'dyn_vor_init : vorticity term : read namelist and control the consistency'
      WRITE(numout, *) '~~~~~~~~~~~~'
    END IF
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namdyn_vor, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namdyn_vor in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namdyn_vor, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namdyn_vor in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namdyn_vor)
    IF (lwp) THEN
      WRITE(numout, *) '   Namelist namdyn_vor : choice of the vorticity term scheme'
      WRITE(numout, *) '      enstrophy conserving scheme                    ln_dynvor_ens = ', ln_dynvor_ens
      WRITE(numout, *) '      f-point energy conserving scheme               ln_dynvor_ene = ', ln_dynvor_ene
      WRITE(numout, *) '      t-point energy conserving scheme               ln_dynvor_enT = ', ln_dynvor_enT
      WRITE(numout, *) '      energy conserving scheme  (een using e3t)      ln_dynvor_eeT = ', ln_dynvor_eeT
      WRITE(numout, *) '      enstrophy and energy conserving scheme         ln_dynvor_een = ', ln_dynvor_een
      WRITE(numout, *) '         e3f = averaging /4 (=0) or /sum(tmask) (=1)    nn_een_e3f = ', nn_een_e3f
      WRITE(numout, *) '      mixed enstrophy/energy conserving scheme       ln_dynvor_mix = ', ln_dynvor_mix
      WRITE(numout, *) '      masked (=T) or unmasked(=F) vorticity          ln_dynvor_msk = ', ln_dynvor_msk
    END IF
    IF (ln_dynvor_msk) CALL ctl_stop('dyn_vor_init:   masked vorticity is not currently not available')
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) '      change fmask value in the angles (T)           ln_vorlat = ', ln_vorlat
    IF (ln_vorlat .AND. (ln_dynvor_ene .OR. ln_dynvor_ens .OR. ln_dynvor_mix)) THEN
      !$ACC KERNELS
      DO jk = 1, jpk
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            IF (tmask(ji, jj + 1, jk) + tmask(ji + 1, jj + 1, jk) + tmask(ji, jj, jk) + tmask(ji + 1, jj + 1, jk) == 3._wp) fmask(ji, jj, jk) = 1._wp
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk('dynvor', fmask, 'F', 1._wp)
    END IF
    ioptio = 0
    IF (ln_dynvor_ens) THEN
      ioptio = ioptio + 1
      nvor_scheme = np_ens
    END IF
    IF (ln_dynvor_ene) THEN
      ioptio = ioptio + 1
      nvor_scheme = np_ene
    END IF
    IF (ln_dynvor_ent) THEN
      ioptio = ioptio + 1
      nvor_scheme = np_ent
    END IF
    IF (ln_dynvor_eet) THEN
      ioptio = ioptio + 1
      nvor_scheme = np_eet
    END IF
    IF (ln_dynvor_een) THEN
      ioptio = ioptio + 1
      nvor_scheme = np_een
    END IF
    IF (ln_dynvor_mix) THEN
      ioptio = ioptio + 1
      nvor_scheme = np_mix
    END IF
    IF (ioptio /= 1) CALL ctl_stop(' use ONE and ONLY one vorticity scheme')
    IF (lwp) WRITE(numout, *)
    ncor = np_COR
    SELECT CASE (n_dynadv)
    CASE (np_LIN_dyn)
      IF (lwp) WRITE(numout, *) '   ==>>>   linear dynamics : total vorticity = Coriolis'
      nrvm = np_COR
      ntot = np_COR
    CASE (np_VEC_c2)
      IF (lwp) WRITE(numout, *) '   ==>>>   vector form dynamics : total vorticity = Coriolis + relative vorticity'
      nrvm = np_RVO
      ntot = np_CRV
    CASE (np_FLX_c2, np_FLX_ubs)
      IF (lwp) WRITE(numout, *) '   ==>>>   flux form dynamics : total vorticity = Coriolis + metric term'
      nrvm = np_MET
      ntot = np_CME
      SELECT CASE (nvor_scheme)
      CASE (np_ENT)
        ALLOCATE(di_e2u_2(jpi, jpj), dj_e1v_2(jpi, jpj))
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            di_e2u_2(ji, jj) = (e2u(ji, jj) - e2u(ji - 1, jj)) * 0.5_wp
            dj_e1v_2(ji, jj) = (e1v(ji, jj) - e1v(ji, jj - 1)) * 0.5_wp
          END DO
        END DO
        !$ACC END KERNELS
        CALL lbc_lnk_multi('dynvor', di_e2u_2, 'T', - 1., dj_e1v_2, 'T', - 1.)
      CASE DEFAULT
        ALLOCATE(di_e2v_2e1e2f(jpi, jpj), dj_e1u_2e1e2f(jpi, jpj))
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            di_e2v_2e1e2f(ji, jj) = (e2v(ji + 1, jj) - e2v(ji, jj)) * 0.5 * r1_e1e2f(ji, jj)
            dj_e1u_2e1e2f(ji, jj) = (e1u(ji, jj + 1) - e1u(ji, jj)) * 0.5 * r1_e1e2f(ji, jj)
          END DO
        END DO
        !$ACC END KERNELS
        CALL lbc_lnk_multi('dynvor', di_e2v_2e1e2f, 'F', - 1., dj_e1u_2e1e2f, 'F', - 1.)
      END SELECT
    END SELECT
    IF (lwp) THEN
      WRITE(numout, *)
      SELECT CASE (nvor_scheme)
      CASE (np_ens)
        WRITE(numout, *) '   ==>>>   enstrophy conserving scheme (ENS)'
      CASE (np_ene)
        WRITE(numout, *) '   ==>>>   energy conserving scheme (Coriolis at F-points) (ENE)'
      CASE (np_ent)
        WRITE(numout, *) '   ==>>>   energy conserving scheme (Coriolis at T-points) (ENT)'
      CASE (np_eet)
        WRITE(numout, *) '   ==>>>   energy conserving scheme (EEN scheme using e3t) (EET)'
      CASE (np_een)
        WRITE(numout, *) '   ==>>>   energy and enstrophy conserving scheme (EEN)'
      CASE (np_mix)
        WRITE(numout, *) '   ==>>>   mixed enstrophy/energy conserving scheme (MIX)'
      END SELECT
    END IF
  END SUBROUTINE dyn_vor_init
END MODULE dynvor