MODULE stopar
  USE storng
  USE par_oce
  USE dom_oce
  USE lbclnk
  USE in_out_manager
  USE iom
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sto_par_init
  PUBLIC :: sto_par
  PUBLIC :: sto_rst_write
  LOGICAL :: ln_rststo = .FALSE.
  LOGICAL :: ln_rstseed = .FALSE.
  CHARACTER(LEN = 32) :: cn_storst_in = "restart_sto"
  CHARACTER(LEN = 32) :: cn_storst_out = "restart_sto"
  INTEGER :: numstor, numstow
  INTEGER :: jpsto2d = 0
  INTEGER :: jpsto3d = 0
  REAL(KIND = wp), PUBLIC, DIMENSION(:, :, :), ALLOCATABLE :: sto2d
  REAL(KIND = wp), PUBLIC, DIMENSION(:, :, :, :), ALLOCATABLE :: sto3d
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: sto_tmp
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: sto2d_abc
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: sto3d_abc
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: sto2d_ave
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: sto3d_ave
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: sto2d_std
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: sto3d_std
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: sto2d_lim
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: sto3d_lim
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: sto2d_tcor
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: sto3d_tcor
  INTEGER, DIMENSION(:), ALLOCATABLE :: sto2d_ord
  INTEGER, DIMENSION(:), ALLOCATABLE :: sto3d_ord
  CHARACTER(LEN = 1), DIMENSION(:), ALLOCATABLE :: sto2d_typ
  CHARACTER(LEN = 1), DIMENSION(:), ALLOCATABLE :: sto3d_typ
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: sto2d_sgn
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: sto3d_sgn
  INTEGER, DIMENSION(:), ALLOCATABLE :: sto2d_flt
  INTEGER, DIMENSION(:), ALLOCATABLE :: sto3d_flt
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: sto2d_fac
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: sto3d_fac
  LOGICAL, PUBLIC :: ln_sto_ldf = .FALSE.
  INTEGER, PUBLIC :: jsto_ldf
  REAL(KIND = wp) :: rn_ldf_std
  REAL(KIND = wp) :: rn_ldf_tcor
  LOGICAL, PUBLIC :: ln_sto_hpg = .FALSE.
  INTEGER, PUBLIC :: jsto_hpgi
  INTEGER, PUBLIC :: jsto_hpgj
  REAL(KIND = wp) :: rn_hpg_std
  REAL(KIND = wp) :: rn_hpg_tcor
  LOGICAL, PUBLIC :: ln_sto_pstar = .FALSE.
  INTEGER, PUBLIC :: jsto_pstar
  REAL(KIND = wp), PUBLIC :: rn_pstar_std
  REAL(KIND = wp) :: rn_pstar_tcor
  INTEGER :: nn_pstar_flt = 0
  INTEGER :: nn_pstar_ord = 1
  LOGICAL, PUBLIC :: ln_sto_trd = .FALSE.
  INTEGER, PUBLIC :: jsto_trd
  REAL(KIND = wp) :: rn_trd_std
  REAL(KIND = wp) :: rn_trd_tcor
  LOGICAL, PUBLIC :: ln_sto_eos = .FALSE.
  INTEGER, PUBLIC :: nn_sto_eos = 1
  INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: jsto_eosi
  INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: jsto_eosj
  INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: jsto_eosk
  REAL(KIND = wp) :: rn_eos_stdxy
  REAL(KIND = wp) :: rn_eos_stdz
  REAL(KIND = wp) :: rn_eos_tcor
  REAL(KIND = wp) :: rn_eos_lim = 3.0_wp
  INTEGER :: nn_eos_flt = 0
  INTEGER :: nn_eos_ord = 1
  LOGICAL, PUBLIC :: ln_sto_trc = .FALSE.
  INTEGER, PUBLIC :: nn_sto_trc = 1
  INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: jsto_trci
  INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: jsto_trcj
  INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: jsto_trck
  REAL(KIND = wp) :: rn_trc_stdxy
  REAL(KIND = wp) :: rn_trc_stdz
  REAL(KIND = wp) :: rn_trc_tcor
  REAL(KIND = wp) :: rn_trc_lim = 3.0_wp
  INTEGER :: nn_trc_flt = 0
  INTEGER :: nn_trc_ord = 1
  CONTAINS
  SUBROUTINE sto_par(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk, jsto, jflt
    REAL(KIND = wp) :: stomax
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('sto_par', 'r0', 0, 0)
    DO jsto = 1, jpsto2d
      sto_tmp(:, :) = sto2d(:, :, jsto)
      IF (sto2d_ord(jsto) == 1) THEN
        CALL sto_par_white(sto2d(:, :, jsto))
        DO jflt = 1, sto2d_flt(jsto)
          CALL lbc_lnk('stopar', sto2d(:, :, jsto), sto2d_typ(jsto), sto2d_sgn(jsto))
          CALL sto_par_flt(sto2d(:, :, jsto))
        END DO
        sto2d(:, :, jsto) = sto2d(:, :, jsto) * sto2d_fac(jsto)
      ELSE
        sto2d(:, :, jsto) = sto2d(:, :, jsto - 1)
      END IF
      sto2d(:, :, jsto) = sto2d(:, :, jsto) * sto2d_abc(jsto, 2)
      sto2d(:, :, jsto) = sto2d(:, :, jsto) + sto_tmp(:, :) * sto2d_abc(jsto, 1)
      sto2d(:, :, jsto) = sto2d(:, :, jsto) + sto2d_abc(jsto, 3)
      stomax = sto2d_std(jsto) * sto2d_lim(jsto)
      sto2d(:, :, jsto) = sto2d(:, :, jsto) - sto2d_ave(jsto)
      sto2d(:, :, jsto) = SIGN(MIN(stomax, ABS(sto2d(:, :, jsto))), sto2d(:, :, jsto))
      sto2d(:, :, jsto) = sto2d(:, :, jsto) + sto2d_ave(jsto)
      CALL lbc_lnk('stopar', sto2d(:, :, jsto), sto2d_typ(jsto), sto2d_sgn(jsto))
    END DO
    DO jsto = 1, jpsto3d
      DO jk = 1, jpk
        sto_tmp(:, :) = sto3d(:, :, jk, jsto)
        IF (sto3d_ord(jsto) == 1) THEN
          CALL sto_par_white(sto3d(:, :, jk, jsto))
          DO jflt = 1, sto3d_flt(jsto)
            CALL lbc_lnk('stopar', sto3d(:, :, jk, jsto), sto3d_typ(jsto), sto3d_sgn(jsto))
            CALL sto_par_flt(sto3d(:, :, jk, jsto))
          END DO
          sto3d(:, :, jk, jsto) = sto3d(:, :, jk, jsto) * sto3d_fac(jsto)
        ELSE
          sto3d(:, :, jk, jsto) = sto3d(:, :, jk, jsto - 1)
        END IF
        sto3d(:, :, jk, jsto) = sto3d(:, :, jk, jsto) * sto3d_abc(jsto, 2)
        sto3d(:, :, jk, jsto) = sto3d(:, :, jk, jsto) + sto_tmp(:, :) * sto3d_abc(jsto, 1)
        sto3d(:, :, jk, jsto) = sto3d(:, :, jk, jsto) + sto3d_abc(jsto, 3)
        stomax = sto3d_std(jsto) * sto3d_lim(jsto)
        sto3d(:, :, jk, jsto) = sto3d(:, :, jk, jsto) - sto3d_ave(jsto)
        sto3d(:, :, jk, jsto) = SIGN(MIN(stomax, ABS(sto3d(:, :, jk, jsto))), sto3d(:, :, jk, jsto))
        sto3d(:, :, jk, jsto) = sto3d(:, :, jk, jsto) + sto3d_ave(jsto)
      END DO
      CALL lbc_lnk('stopar', sto3d(:, :, :, jsto), sto3d_typ(jsto), sto3d_sgn(jsto))
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE sto_par
  SUBROUTINE sto_par_init
    NAMELIST /namsto/ ln_sto_ldf, rn_ldf_std, rn_ldf_tcor, ln_sto_hpg, rn_hpg_std, rn_hpg_tcor, ln_sto_pstar, rn_pstar_std, rn_pstar_tcor, nn_pstar_flt, nn_pstar_ord, ln_sto_trd, rn_trd_std, rn_trd_tcor, ln_sto_eos, nn_sto_eos, rn_eos_stdxy, rn_eos_stdz, rn_eos_tcor, nn_eos_ord, nn_eos_flt, rn_eos_lim, ln_sto_trc, nn_sto_trc, rn_trc_stdxy, rn_trc_stdz, rn_trc_tcor, nn_trc_ord, nn_trc_flt, rn_trc_lim, ln_rststo, ln_rstseed, cn_storst_in, cn_storst_out
    INTEGER :: jsto, jmem, jarea, jdof, jord, jordm1, jk, jflt
    INTEGER(KIND = 8) :: zseed1, zseed2, zseed3, zseed4
    REAL(KIND = wp) :: rinflate
    INTEGER :: ios
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namsto, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namsto in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namsto, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namsto in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namsto)
    IF (.NOT. ln_rststo) THEN
      IF (lwp) THEN
        WRITE(numout, *)
        WRITE(numout, *) 'sto_par_init : NO use of stochastic parameterization'
        WRITE(numout, *) '~~~~~~~~~~~~'
      END IF
      RETURN
    END IF
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'sto_par_init : stochastic parameterization'
      WRITE(numout, *) '~~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist namsto : stochastic parameterization'
      WRITE(numout, *) '      restart stochastic parameters           ln_rststo     = ', ln_rststo
      WRITE(numout, *) '      read seed of RNG from restart file      ln_rstseed    = ', ln_rstseed
      WRITE(numout, *) '      suffix of sto restart name (input)      cn_storst_in  = ', cn_storst_in
      WRITE(numout, *) '      suffix of sto restart name (output)     cn_storst_out = ', cn_storst_out
      WRITE(numout, *) '      stochastic equation of state            ln_sto_eos    = ', ln_sto_eos
      WRITE(numout, *) '      number of degrees of freedom            nn_sto_eos    = ', nn_sto_eos
      WRITE(numout, *) '      random walk horz. std (in grid points)  rn_eos_stdxy  = ', rn_eos_stdxy
      WRITE(numout, *) '      random walk vert. std (in grid points)  rn_eos_stdz   = ', rn_eos_stdz
      WRITE(numout, *) '      random walk tcor (in timesteps)         rn_eos_tcor   = ', rn_eos_tcor
      WRITE(numout, *) '      order of autoregressive  processes      nn_eos_ord    = ', nn_eos_ord
      WRITE(numout, *) '      passes of Laplacian filter              nn_eos_flt    = ', nn_eos_flt
      WRITE(numout, *) '      limitation factor                       rn_eos_lim    = ', rn_eos_lim
    END IF
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) '   stochastic parameterization :'
    jpsto2d = 0
    IF (ln_sto_ldf) THEN
      IF (lwp) WRITE(numout, *) '       - stochastic lateral diffusion'
      jpsto2d = jpsto2d + 1
      jsto_ldf = jpsto2d
    END IF
    IF (ln_sto_pstar) THEN
      IF (lwp) WRITE(numout, *) '       - stochastic ice strength'
      jpsto2d = jpsto2d + 1 * nn_pstar_ord
      jsto_pstar = jpsto2d
    END IF
    IF (ln_sto_eos) THEN
      IF (lk_agrif) CALL ctl_stop('EOS stochastic parametrization is not compatible with AGRIF')
      IF (lwp) WRITE(numout, *) '       - stochastic equation of state'
      ALLOCATE(jsto_eosi(nn_sto_eos))
      ALLOCATE(jsto_eosj(nn_sto_eos))
      ALLOCATE(jsto_eosk(nn_sto_eos))
      DO jdof = 1, nn_sto_eos
        jpsto2d = jpsto2d + 3 * nn_eos_ord
        jsto_eosi(jdof) = jpsto2d - 2 * nn_eos_ord
        jsto_eosj(jdof) = jpsto2d - 1 * nn_eos_ord
        jsto_eosk(jdof) = jpsto2d
      END DO
    ELSE
      nn_sto_eos = 0
    END IF
    IF (ln_sto_trc) THEN
      IF (lwp) WRITE(numout, *) '       - stochastic tracers dynamics'
      ALLOCATE(jsto_trci(nn_sto_trc))
      ALLOCATE(jsto_trcj(nn_sto_trc))
      ALLOCATE(jsto_trck(nn_sto_trc))
      DO jdof = 1, nn_sto_trc
        jpsto2d = jpsto2d + 3 * nn_trc_ord
        jsto_trci(jdof) = jpsto2d - 2 * nn_trc_ord
        jsto_trcj(jdof) = jpsto2d - 1 * nn_trc_ord
        jsto_trck(jdof) = jpsto2d
      END DO
    ELSE
      nn_sto_trc = 0
    END IF
    jpsto3d = 0
    IF (ln_sto_hpg) THEN
      IF (lwp) WRITE(numout, *) '       - stochastic horizontal pressure gradient'
      jpsto3d = jpsto3d + 2
      jsto_hpgi = jpsto3d - 1
      jsto_hpgj = jpsto3d
    END IF
    IF (ln_sto_trd) THEN
      IF (lwp) WRITE(numout, *) '       - stochastic trend'
      jpsto3d = jpsto3d + 1
      jsto_trd = jpsto3d
    END IF
    IF (jpsto2d > 0) THEN
      ALLOCATE(sto2d(jpi, jpj, jpsto2d))
      ALLOCATE(sto2d_abc(jpsto2d, 3))
      ALLOCATE(sto2d_ave(jpsto2d))
      ALLOCATE(sto2d_std(jpsto2d))
      ALLOCATE(sto2d_lim(jpsto2d))
      ALLOCATE(sto2d_tcor(jpsto2d))
      ALLOCATE(sto2d_ord(jpsto2d))
      ALLOCATE(sto2d_typ(jpsto2d))
      ALLOCATE(sto2d_sgn(jpsto2d))
      ALLOCATE(sto2d_flt(jpsto2d))
      ALLOCATE(sto2d_fac(jpsto2d))
    END IF
    IF (jpsto3d > 0) THEN
      ALLOCATE(sto3d(jpi, jpj, jpk, jpsto3d))
      ALLOCATE(sto3d_abc(jpsto3d, 3))
      ALLOCATE(sto3d_ave(jpsto3d))
      ALLOCATE(sto3d_std(jpsto3d))
      ALLOCATE(sto3d_lim(jpsto3d))
      ALLOCATE(sto3d_tcor(jpsto3d))
      ALLOCATE(sto3d_ord(jpsto3d))
      ALLOCATE(sto3d_typ(jpsto3d))
      ALLOCATE(sto3d_sgn(jpsto3d))
      ALLOCATE(sto3d_flt(jpsto3d))
      ALLOCATE(sto3d_fac(jpsto3d))
    END IF
    IF (jpsto2d > 0 .OR. jpsto3d > 0) THEN
      ALLOCATE(sto_tmp(jpi, jpj))
      sto_tmp(:, :) = 0._wp
    END IF
    DO jsto = 1, jpsto2d
      sto2d_typ(jsto) = 'T'
      sto2d_sgn(jsto) = 1._wp
      sto2d_flt(jsto) = 0
      sto2d_ord(jsto) = 1
      DO jord = 0, nn_pstar_ord - 1
        IF (jsto + jord == jsto_pstar) THEN
          sto2d_ord(jsto) = nn_pstar_ord - jord
          sto2d_flt(jsto) = nn_pstar_flt
        END IF
      END DO
      DO jdof = 1, nn_sto_eos
        DO jord = 0, nn_eos_ord - 1
          IF (jsto + jord == jsto_eosi(jdof)) THEN
            sto2d_ord(jsto) = nn_eos_ord - jord
            sto2d_sgn(jsto) = - 1._wp
            sto2d_flt(jsto) = nn_eos_flt
          END IF
          IF (jsto + jord == jsto_eosj(jdof)) THEN
            sto2d_ord(jsto) = nn_eos_ord - jord
            sto2d_sgn(jsto) = - 1._wp
            sto2d_flt(jsto) = nn_eos_flt
          END IF
          IF (jsto + jord == jsto_eosk(jdof)) THEN
            sto2d_ord(jsto) = nn_eos_ord - jord
            sto2d_flt(jsto) = nn_eos_flt
          END IF
        END DO
      END DO
      DO jdof = 1, nn_sto_trc
        DO jord = 0, nn_trc_ord - 1
          IF (jsto + jord == jsto_trci(jdof)) THEN
            sto2d_ord(jsto) = nn_trc_ord - jord
            sto2d_sgn(jsto) = - 1._wp
            sto2d_flt(jsto) = nn_trc_flt
          END IF
          IF (jsto + jord == jsto_trcj(jdof)) THEN
            sto2d_ord(jsto) = nn_trc_ord - jord
            sto2d_sgn(jsto) = - 1._wp
            sto2d_flt(jsto) = nn_trc_flt
          END IF
          IF (jsto + jord == jsto_trck(jdof)) THEN
            sto2d_ord(jsto) = nn_trc_ord - jord
            sto2d_flt(jsto) = nn_trc_flt
          END IF
        END DO
      END DO
      sto2d_fac(jsto) = sto_par_flt_fac(sto2d_flt(jsto))
    END DO
    DO jsto = 1, jpsto3d
      sto3d_typ(jsto) = 'T'
      sto3d_sgn(jsto) = 1._wp
      sto3d_flt(jsto) = 0
      sto3d_ord(jsto) = 1
      IF (jsto == jsto_hpgi) THEN
        sto3d_typ(jsto) = 'U'
      END IF
      IF (jsto == jsto_hpgj) THEN
        sto3d_typ(jsto) = 'V'
      END IF
      sto3d_fac(jsto) = sto_par_flt_fac(sto3d_flt(jsto))
    END DO
    DO jsto = 1, jpsto2d
      sto2d_ave(jsto) = 0._wp
      sto2d_std(jsto) = 1._wp
      sto2d_tcor(jsto) = 1._wp
      sto2d_lim(jsto) = 3._wp
      IF (jsto == jsto_ldf) THEN
        sto2d_ave(jsto) = 1._wp
        sto2d_std(jsto) = rn_ldf_std
        sto2d_tcor(jsto) = rn_ldf_tcor
      END IF
      DO jord = 0, nn_pstar_ord - 1
        IF (jsto + jord == jsto_pstar) THEN
          sto2d_std(jsto) = 1._wp
          sto2d_tcor(jsto) = rn_pstar_tcor
        END IF
      END DO
      DO jdof = 1, nn_sto_eos
        DO jord = 0, nn_eos_ord - 1
          IF (jsto + jord == jsto_eosi(jdof)) THEN
            sto2d_std(jsto) = rn_eos_stdxy
            sto2d_tcor(jsto) = rn_eos_tcor
            sto2d_lim(jsto) = rn_eos_lim
          END IF
          IF (jsto + jord == jsto_eosj(jdof)) THEN
            sto2d_std(jsto) = rn_eos_stdxy
            sto2d_tcor(jsto) = rn_eos_tcor
            sto2d_lim(jsto) = rn_eos_lim
          END IF
          IF (jsto + jord == jsto_eosk(jdof)) THEN
            sto2d_std(jsto) = rn_eos_stdz
            sto2d_tcor(jsto) = rn_eos_tcor
            sto2d_lim(jsto) = rn_eos_lim
          END IF
        END DO
      END DO
      DO jdof = 1, nn_sto_trc
        DO jord = 0, nn_trc_ord - 1
          IF (jsto + jord == jsto_trci(jdof)) THEN
            sto2d_std(jsto) = rn_trc_stdxy
            sto2d_tcor(jsto) = rn_trc_tcor
            sto2d_lim(jsto) = rn_trc_lim
          END IF
          IF (jsto + jord == jsto_trcj(jdof)) THEN
            sto2d_std(jsto) = rn_trc_stdxy
            sto2d_tcor(jsto) = rn_trc_tcor
            sto2d_lim(jsto) = rn_trc_lim
          END IF
          IF (jsto + jord == jsto_trck(jdof)) THEN
            sto2d_std(jsto) = rn_trc_stdz
            sto2d_tcor(jsto) = rn_trc_tcor
            sto2d_lim(jsto) = rn_trc_lim
          END IF
        END DO
      END DO
    END DO
    DO jsto = 1, jpsto3d
      sto3d_ave(jsto) = 0._wp
      sto3d_std(jsto) = 1._wp
      sto3d_tcor(jsto) = 1._wp
      sto3d_lim(jsto) = 3._wp
      IF (jsto == jsto_hpgi) THEN
        sto3d_ave(jsto) = 1._wp
        sto3d_std(jsto) = rn_hpg_std
        sto3d_tcor(jsto) = rn_hpg_tcor
      END IF
      IF (jsto == jsto_hpgj) THEN
        sto3d_ave(jsto) = 1._wp
        sto3d_std(jsto) = rn_hpg_std
        sto3d_tcor(jsto) = rn_hpg_tcor
      END IF
      IF (jsto == jsto_trd) THEN
        sto3d_ave(jsto) = 1._wp
        sto3d_std(jsto) = rn_trd_std
        sto3d_tcor(jsto) = rn_trd_tcor
      END IF
    END DO
    DO jsto = 1, jpsto2d
      IF (sto2d_tcor(jsto) == 0._wp) THEN
        sto2d_abc(jsto, 1) = 0._wp
      ELSE
        sto2d_abc(jsto, 1) = EXP(- 1._wp / sto2d_tcor(jsto))
      END IF
      IF (sto2d_ord(jsto) == 1) THEN
        rinflate = sto2d_std(jsto)
      ELSE
        jordm1 = sto2d_ord(jsto) - 1
        rinflate = SQRT(REAL(jordm1, wp) / REAL(2 * (2 * jordm1 - 1), wp))
      END IF
      sto2d_abc(jsto, 2) = rinflate * SQRT(1._wp - sto2d_abc(jsto, 1) * sto2d_abc(jsto, 1))
      sto2d_abc(jsto, 3) = sto2d_ave(jsto) * (1._wp - sto2d_abc(jsto, 1))
    END DO
    DO jsto = 1, jpsto3d
      IF (sto3d_tcor(jsto) == 0._wp) THEN
        sto3d_abc(jsto, 1) = 0._wp
      ELSE
        sto3d_abc(jsto, 1) = EXP(- 1._wp / sto3d_tcor(jsto))
      END IF
      IF (sto3d_ord(jsto) == 1) THEN
        rinflate = sto3d_std(jsto)
      ELSE
        jordm1 = sto3d_ord(jsto) - 1
        rinflate = SQRT(REAL(jordm1, wp) / REAL(2 * (2 * jordm1 - 1), wp))
      END IF
      sto3d_abc(jsto, 2) = rinflate * SQRT(1._wp - sto3d_abc(jsto, 1) * sto3d_abc(jsto, 1))
      sto3d_abc(jsto, 3) = sto3d_ave(jsto) * (1._wp - sto3d_abc(jsto, 1))
    END DO
    CALL kiss_reset
    DO jarea = 1, narea
      zseed1 = kiss()
      zseed2 = kiss()
      zseed3 = kiss()
      zseed4 = kiss()
    END DO
    CALL kiss_seed(zseed1, zseed2, zseed3, zseed4)
    DO jsto = 1, jpsto2d
      CALL sto_par_white(sto2d(:, :, jsto))
      DO jflt = 1, sto2d_flt(jsto)
        CALL lbc_lnk('stopar', sto2d(:, :, jsto), sto2d_typ(jsto), sto2d_sgn(jsto))
        CALL sto_par_flt(sto2d(:, :, jsto))
      END DO
      sto2d(:, :, jsto) = sto2d(:, :, jsto) * sto2d_fac(jsto)
      sto2d(:, :, jsto) = SIGN(MIN(sto2d_lim(jsto), ABS(sto2d(:, :, jsto))), sto2d(:, :, jsto))
      sto2d(:, :, jsto) = sto2d(:, :, jsto) * sto2d_std(jsto) + sto2d_ave(jsto)
    END DO
    DO jsto = 1, jpsto3d
      DO jk = 1, jpk
        CALL sto_par_white(sto3d(:, :, jk, jsto))
        DO jflt = 1, sto3d_flt(jsto)
          CALL lbc_lnk('stopar', sto3d(:, :, jk, jsto), sto3d_typ(jsto), sto3d_sgn(jsto))
          CALL sto_par_flt(sto3d(:, :, jk, jsto))
        END DO
        sto3d(:, :, jk, jsto) = sto3d(:, :, jk, jsto) * sto3d_fac(jsto)
        sto3d(:, :, jk, jsto) = SIGN(MIN(sto3d_lim(jsto), ABS(sto3d(:, :, jk, jsto))), sto3d(:, :, jk, jsto))
        sto3d(:, :, jk, jsto) = sto3d(:, :, jk, jsto) * sto3d_std(jsto) + sto3d_ave(jsto)
      END DO
    END DO
    IF (ln_rststo) CALL sto_rst_read
  END SUBROUTINE sto_par_init
  SUBROUTINE sto_rst_read
    INTEGER :: jsto, jseed
    INTEGER(KIND = 8) :: ziseed(4)
    REAL(KIND = 8) :: zrseed(4)
    CHARACTER(LEN = 9) :: clsto2d = 'sto2d_000'
    CHARACTER(LEN = 9) :: clsto3d = 'sto3d_000'
    CHARACTER(LEN = 10) :: clseed = 'seed0_0000'
    IF (jpsto2d > 0 .OR. jpsto3d > 0) THEN
      IF (lwp) THEN
        WRITE(numout, *)
        WRITE(numout, *) 'sto_rst_read : read stochastic parameters from restart file'
        WRITE(numout, *) '~~~~~~~~~~~~'
      END IF
      CALL iom_open(cn_storst_in, numstor)
      DO jsto = 1, jpsto2d
        WRITE(clsto2d(7 : 9), '(i3.3)') jsto
        CALL iom_get(numstor, jpdom_autoglo, clsto2d, sto2d(:, :, jsto))
      END DO
      DO jsto = 1, jpsto3d
        WRITE(clsto3d(7 : 9), '(i3.3)') jsto
        CALL iom_get(numstor, jpdom_autoglo, clsto3d, sto3d(:, :, :, jsto))
      END DO
      IF (ln_rstseed) THEN
        DO jseed = 1, 4
          WRITE(clseed(5 : 5), '(i1.1)') jseed
          WRITE(clseed(7 : 10), '(i4.4)') narea
          CALL iom_get(numstor, clseed, zrseed(jseed))
        END DO
        ziseed = TRANSFER(zrseed, ziseed)
        CALL kiss_seed(ziseed(1), ziseed(2), ziseed(3), ziseed(4))
      END IF
      CALL iom_close(numstor)
    END IF
  END SUBROUTINE sto_rst_read
  SUBROUTINE sto_rst_write(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: jsto, jseed
    INTEGER(KIND = 8) :: ziseed(4)
    REAL(KIND = 8) :: zrseed(4)
    CHARACTER(LEN = 20) :: clkt
    CHARACTER(LEN = 50) :: clname
    CHARACTER(LEN = 9) :: clsto2d = 'sto2d_000'
    CHARACTER(LEN = 9) :: clsto3d = 'sto3d_000'
    CHARACTER(LEN = 10) :: clseed = 'seed0_0000'
    IF (jpsto2d > 0 .OR. jpsto3d > 0) THEN
      IF (kt == nitrst .OR. kt == nitend) THEN
        IF (lwp) THEN
          WRITE(numout, *)
          WRITE(numout, *) 'sto_rst_write : write stochastic parameters in restart file'
          WRITE(numout, *) '~~~~~~~~~~~~~'
        END IF
      END IF
      IF (kt > nit000) THEN
        IF (kt == nitrst .OR. kt == nitend) THEN
          CALL kiss_state(ziseed(1), ziseed(2), ziseed(3), ziseed(4))
          zrseed = TRANSFER(ziseed, zrseed)
          DO jseed = 1, 4
            WRITE(clseed(5 : 5), '(i1.1)') jseed
            WRITE(clseed(7 : 10), '(i4.4)') narea
            CALL iom_rstput(kt, nitrst, numstow, clseed, zrseed(jseed))
          END DO
          DO jsto = 1, jpsto2d
            WRITE(clsto2d(7 : 9), '(i3.3)') jsto
            CALL iom_rstput(kt, nitrst, numstow, clsto2d, sto2d(:, :, jsto))
          END DO
          DO jsto = 1, jpsto3d
            WRITE(clsto3d(7 : 9), '(i3.3)') jsto
            CALL iom_rstput(kt, nitrst, numstow, clsto3d, sto3d(:, :, :, jsto))
          END DO
          CALL iom_close(numstow)
        END IF
      END IF
      IF (kt < nitend) THEN
        IF (kt == nitrst - 1 .OR. nstock == 1 .OR. kt == nitend - 1) THEN
          IF (nitrst > 999999999) THEN
            WRITE(clkt, *) nitrst
          ELSE
            WRITE(clkt, '(i8.8)') nitrst
          END IF
          clname = TRIM(cexper) // "_" // TRIM(ADJUSTL(clkt)) // "_" // TRIM(cn_storst_out)
          IF (lwp) THEN
            WRITE(numout, *) '             open stochastic parameters restart file: ' // clname
            IF (kt == nitrst - 1) THEN
              WRITE(numout, *) '             kt = nitrst - 1 = ', kt
            ELSE
              WRITE(numout, *) '             kt = ', kt
            END IF
          END IF
          CALL iom_open(clname, numstow, ldwrt = .TRUE.)
        END IF
      END IF
    END IF
  END SUBROUTINE sto_rst_write
  SUBROUTINE sto_par_white(psto)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(OUT) :: psto
    INTEGER :: ji, jj
    REAL(KIND = 8) :: gran
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('sto_par_white', 'r0', 0, 0)
    DO jj = 1, jpj
      DO ji = 1, jpi
        CALL kiss_gaussian(gran)
        psto(ji, jj) = gran
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE sto_par_white
  SUBROUTINE sto_par_flt(psto)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(OUT) :: psto
    INTEGER :: ji, jj
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('sto_par_flt', 'r0', 0, 0)
    DO jj = 2, jpj - 1
      DO ji = 2, jpi - 1
        psto(ji, jj) = 0.5_wp * psto(ji, jj) + 0.125_wp * (psto(ji - 1, jj) + psto(ji + 1, jj) + psto(ji, jj - 1) + psto(ji, jj + 1))
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE sto_par_flt
  FUNCTION sto_par_flt_fac(kpasses)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kpasses
    REAL(KIND = wp) :: sto_par_flt_fac
    INTEGER :: jpasses, ji, jj, jflti, jfltj
    INTEGER, DIMENSION(- 1 : 1, - 1 : 1) :: pflt0
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: pfltb
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: pflta
    REAL(KIND = wp) :: ratio
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('sto_par_flt_fac', 'r0', 0, 0)
    pflt0(- 1, - 1) = 0
    pflt0(- 1, 0) = 1
    pflt0(- 1, 1) = 0
    pflt0(0, - 1) = 1
    pflt0(0, 0) = 4
    pflt0(0, 1) = 1
    pflt0(1, - 1) = 0
    pflt0(1, 0) = 1
    pflt0(1, 1) = 0
    ALLOCATE(pfltb(- kpasses - 1 : kpasses + 1, - kpasses - 1 : kpasses + 1))
    ALLOCATE(pflta(- kpasses - 1 : kpasses + 1, - kpasses - 1 : kpasses + 1))
    pfltb(:, :) = 0
    pfltb(0, 0) = 1
    DO jpasses = 1, kpasses
      pflta(:, :) = 0
      DO jflti = - 1, 1
        DO jfltj = - 1, 1
          DO ji = - kpasses, kpasses
            DO jj = - kpasses, kpasses
              pflta(ji, jj) = pflta(ji, jj) + pfltb(ji + jflti, jj + jfltj) * pflt0(jflti, jfltj)
            END DO
          END DO
        END DO
      END DO
      pfltb(:, :) = pflta(:, :)
    END DO
    ratio = SUM(pfltb(:, :))
    ratio = ratio * ratio / SUM(pfltb(:, :) * pfltb(:, :))
    ratio = SQRT(ratio)
    DEALLOCATE(pfltb, pflta)
    sto_par_flt_fac = ratio
    CALL profile_psy_data0 % PostEnd
  END FUNCTION sto_par_flt_fac
END MODULE stopar