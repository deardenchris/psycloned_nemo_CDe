MODULE ldfdyn
  USE oce
  USE dom_oce
  USE phycst
  USE ldfslp
  USE ldfc1d_c2d
  USE in_out_manager
  USE iom
  USE timing
  USE lib_mpp
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ldf_dyn_init
  PUBLIC :: ldf_dyn
  LOGICAL, PUBLIC :: ln_dynldf_OFF
  LOGICAL, PUBLIC :: ln_dynldf_lap
  LOGICAL, PUBLIC :: ln_dynldf_blp
  LOGICAL, PUBLIC :: ln_dynldf_lev
  LOGICAL, PUBLIC :: ln_dynldf_hor
  INTEGER, PUBLIC :: nn_ahm_ijk_t
  REAL(KIND = wp), PUBLIC :: rn_Uv
  REAL(KIND = wp), PUBLIC :: rn_Lv
  REAL(KIND = wp), PUBLIC :: rn_csmc
  REAL(KIND = wp), PUBLIC :: rn_minfac
  REAL(KIND = wp), PUBLIC :: rn_maxfac
  REAL(KIND = wp), PUBLIC :: rn_ahm_b
  INTEGER, PARAMETER, PUBLIC :: np_ERROR = - 10
  INTEGER, PARAMETER, PUBLIC :: np_no_ldf = 00
  INTEGER, PARAMETER, PUBLIC :: np_lap = 10, np_blp = 20
  INTEGER, PARAMETER, PUBLIC :: np_lap_i = 11
  INTEGER, PUBLIC :: nldf_dyn
  LOGICAL, PUBLIC :: l_ldfdyn_time
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: ahmt, ahmf
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: dtensq
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: dshesq
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: esqt, esqf
  REAL(KIND = wp) :: r1_2 = 0.5_wp
  REAL(KIND = wp) :: r1_4 = 0.25_wp
  REAL(KIND = wp) :: r1_8 = 0.125_wp
  REAL(KIND = wp) :: r1_12 = 1._wp / 12._wp
  REAL(KIND = wp) :: r1_288 = 1._wp / 288._wp
  CONTAINS
  SUBROUTINE ldf_dyn_init
    INTEGER :: ji, jj, jk
    INTEGER :: ioptio, ierr, inum, ios, inn
    REAL(KIND = wp) :: zah0, zah_max, zUfac
    CHARACTER(LEN = 5) :: cl_Units
    NAMELIST /namdyn_ldf/ ln_dynldf_OFF, ln_dynldf_lap, ln_dynldf_blp, ln_dynldf_lev, ln_dynldf_hor, ln_dynldf_iso, nn_ahm_ijk_t, &
&rn_Uv, rn_Lv, rn_ahm_b, rn_csmc, rn_minfac, rn_maxfac
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namdyn_ldf, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namdyn_ldf in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namdyn_ldf, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namdyn_ldf in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namdyn_ldf)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'ldf_dyn : lateral momentum physics'
      WRITE(numout, FMT = *) '~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namdyn_ldf : set lateral mixing parameters'
      WRITE(numout, FMT = *) '      type :'
      WRITE(numout, FMT = *) '         no explicit diffusion                ln_dynldf_OFF = ', ln_dynldf_OFF
      WRITE(numout, FMT = *) '         laplacian operator                   ln_dynldf_lap = ', ln_dynldf_lap
      WRITE(numout, FMT = *) '         bilaplacian operator                 ln_dynldf_blp = ', ln_dynldf_blp
      WRITE(numout, FMT = *) '      direction of action :'
      WRITE(numout, FMT = *) '         iso-level                            ln_dynldf_lev = ', ln_dynldf_lev
      WRITE(numout, FMT = *) '         horizontal (geopotential)            ln_dynldf_hor = ', ln_dynldf_hor
      WRITE(numout, FMT = *) '         iso-neutral                          ln_dynldf_iso = ', ln_dynldf_iso
      WRITE(numout, FMT = *) '      coefficients :'
      WRITE(numout, FMT = *) '         type of time-space variation         nn_ahm_ijk_t  = ', nn_ahm_ijk_t
      WRITE(numout, FMT = *) '         lateral viscous velocity  (if cst)      rn_Uv      = ', rn_Uv, ' m/s'
      WRITE(numout, FMT = *) '         lateral viscous length    (if cst)      rn_Lv      = ', rn_Lv, ' m'
      WRITE(numout, FMT = *) '         background viscosity (iso-lap case)     rn_ahm_b   = ', rn_ahm_b, ' m2/s'
      WRITE(numout, FMT = *) '      Smagorinsky settings (nn_ahm_ijk_t  = 32) :'
      WRITE(numout, FMT = *) '         Smagorinsky coefficient              rn_csmc       = ', rn_csmc
      WRITE(numout, FMT = *) '         factor multiplier for eddy visc.'
      WRITE(numout, FMT = *) '            lower limit (default 1.0)         rn_minfac    = ', rn_minfac
      WRITE(numout, FMT = *) '            upper limit (default 1.0)         rn_maxfac    = ', rn_maxfac
    END IF
    nldf_dyn = np_ERROR
    ioptio = 0
    IF (ln_dynldf_off) THEN
      nldf_dyn = np_no_ldf
      ioptio = ioptio + 1
    END IF
    IF (ln_dynldf_lap) THEN
      ioptio = ioptio + 1
    END IF
    IF (ln_dynldf_blp) THEN
      ioptio = ioptio + 1
    END IF
    IF (ioptio /= 1) CALL ctl_stop('dyn_ldf_init: use ONE of the 3 operator options (NONE/lap/blp)')
    IF (.NOT. ln_dynldf_OFF) THEN
      ioptio = 0
      IF (ln_dynldf_lev) ioptio = ioptio + 1
      IF (ln_dynldf_hor) ioptio = ioptio + 1
      IF (ln_dynldf_iso) ioptio = ioptio + 1
      IF (ioptio /= 1) CALL ctl_stop('dyn_ldf_init: use ONE of the 3 direction options (level/hor/iso)')
      ierr = 0
      IF (ln_dynldf_lap) THEN
        IF (ln_zco) THEN
          IF (ln_dynldf_lev) nldf_dyn = np_lap
          IF (ln_dynldf_hor) nldf_dyn = np_lap
          IF (ln_dynldf_iso) nldf_dyn = np_lap_i
        END IF
        IF (ln_zps) THEN
          IF (ln_dynldf_lev) nldf_dyn = np_lap
          IF (ln_dynldf_hor) nldf_dyn = np_lap
          IF (ln_dynldf_iso) nldf_dyn = np_lap_i
        END IF
        IF (ln_sco) THEN
          IF (ln_dynldf_lev) nldf_dyn = np_lap
          IF (ln_dynldf_hor) nldf_dyn = np_lap_i
          IF (ln_dynldf_iso) nldf_dyn = np_lap_i
        END IF
      END IF
      IF (ln_dynldf_blp) THEN
        IF (ln_zco) THEN
          IF (ln_dynldf_lev) nldf_dyn = np_blp
          IF (ln_dynldf_hor) nldf_dyn = np_blp
          IF (ln_dynldf_iso) ierr = 2
        END IF
        IF (ln_zps) THEN
          IF (ln_dynldf_lev) nldf_dyn = np_blp
          IF (ln_dynldf_hor) nldf_dyn = np_blp
          IF (ln_dynldf_iso) ierr = 2
        END IF
        IF (ln_sco) THEN
          IF (ln_dynldf_lev) nldf_dyn = np_blp
          IF (ln_dynldf_hor) ierr = 2
          IF (ln_dynldf_iso) ierr = 2
        END IF
      END IF
      IF (ierr == 2) CALL ctl_stop('rotated bi-laplacian operator does not exist')
      IF (nldf_dyn == np_lap_i) l_ldfslp = .TRUE.
    END IF
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      SELECT CASE (nldf_dyn)
      CASE (np_no_ldf)
        WRITE(numout, FMT = *) '   ==>>>   NO lateral viscosity'
      CASE (np_lap)
        WRITE(numout, FMT = *) '   ==>>>   iso-level laplacian operator'
      CASE (np_lap_i)
        WRITE(numout, FMT = *) '   ==>>>   rotated laplacian operator with iso-level background'
      CASE (np_blp)
        WRITE(numout, FMT = *) '   ==>>>   iso-level bi-laplacian operator'
      END SELECT
      WRITE(numout, FMT = *)
    END IF
    l_ldfdyn_time = .FALSE.
    IF (ln_dynldf_OFF) THEN
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   No viscous operator selected. ahmt and ahmf are not allocated'
      RETURN
    ELSE
      ALLOCATE(ahmt(jpi, jpj, jpk), ahmf(jpi, jpj, jpk), STAT = ierr)
      IF (ierr /= 0) CALL ctl_stop('STOP', 'ldf_dyn_init: failed to allocate arrays')
      !$ACC KERNELS
      ahmt(:, :, :) = 0._wp
      ahmf(:, :, :) = 0._wp
      !$ACC END KERNELS
      IF (ln_dynldf_lap) THEN
        zufac = r1_2 * rn_uv
        inn = 1
        cl_units = ' m2/s'
      ELSE IF (ln_dynldf_blp) THEN
        zufac = r1_12 * rn_uv
        inn = 3
        cl_units = ' m4/s'
      END IF
      zah0 = zUfac * rn_Lv ** inn
      zah_max = zUfac * (ra * rad) ** inn
      SELECT CASE (nn_ahm_ijk_t)
      CASE (0)
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   eddy viscosity. = constant = ', zah0, cl_Units
        !$ACC KERNELS
        ahmt(:, :, 1 : jpkm1) = zah0
        ahmf(:, :, 1 : jpkm1) = zah0
        !$ACC END KERNELS
      CASE (10)
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   eddy viscosity = F( depth )'
        IF (lwp) WRITE(numout, FMT = *) '           surface viscous coef. = constant = ', zah0, cl_Units
        !$ACC KERNELS
        ahmt(:, :, 1) = zah0
        ahmf(:, :, 1) = zah0
        !$ACC END KERNELS
        CALL ldf_c1d('DYN', ahmt(:, :, 1), ahmf(:, :, 1), ahmt, ahmf)
      CASE (- 20)
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   eddy viscosity = F(i,j) read in eddy_viscosity.nc file'
        CALL iom_open('eddy_viscosity_2D.nc', inum)
        CALL iom_get(inum, jpdom_data, 'ahmt_2d', ahmt(:, :, 1))
        CALL iom_get(inum, jpdom_data, 'ahmf_2d', ahmf(:, :, 1))
        CALL iom_close(inum)
        DO jk = 2, jpkm1
          !$ACC KERNELS
          ahmt(:, :, jk) = ahmt(:, :, 1)
          ahmf(:, :, jk) = ahmf(:, :, 1)
          !$ACC END KERNELS
        END DO
      CASE (20)
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   eddy viscosity = F( e1, e2 ) or F( e1^3, e2^3 ) (lap. or blp. case)'
        IF (lwp) WRITE(numout, FMT = *) '           using a fixed viscous velocity = ', rn_Uv, ' m/s   and   Lv = Max(e1,e2)'
        IF (lwp) WRITE(numout, FMT = *) '           maximum reachable coefficient (at the Equator) = ', zah_max, cl_Units, '  for &
&e1=1°)'
        CALL ldf_c2d('DYN', zUfac, inn, ahmt, ahmf)
      CASE (- 30)
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   eddy viscosity = F(i,j,k) read in eddy_viscosity_3D.nc file'
        CALL iom_open('eddy_viscosity_3D.nc', inum)
        CALL iom_get(inum, jpdom_data, 'ahmt_3d', ahmt)
        CALL iom_get(inum, jpdom_data, 'ahmf_3d', ahmf)
        CALL iom_close(inum)
      CASE (30)
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   eddy viscosity = F( latitude, longitude, depth )'
        IF (lwp) WRITE(numout, FMT = *) '           using a fixed viscous velocity = ', rn_Uv, ' m/s   and   Ld = Max(e1,e2)'
        IF (lwp) WRITE(numout, FMT = *) '           maximum reachable coefficient (at the Equator) = ', zah_max, cl_Units, '  for &
&e1=1°)'
        CALL ldf_c2d('DYN', zUfac, inn, ahmt, ahmf)
        CALL ldf_c1d('DYN', ahmt(:, :, 1), ahmf(:, :, 1), ahmt, ahmf)
      CASE (31)
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   eddy viscosity = F( latitude, longitude, depth , time )'
        IF (lwp) WRITE(numout, FMT = *) '           proportional to the local velocity : 1/2 |u|e (lap) or 1/12 |u|e^3 (blp)'
        l_ldfdyn_time = .TRUE.
      CASE (32)
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   eddy viscosity = F( latitude, longitude, depth , time )'
        IF (lwp) WRITE(numout, FMT = *) '           proportional to the local deformation rate and gridscale (Smagorinsky)'
        l_ldfdyn_time = .TRUE.
        ALLOCATE(dtensq(jpi, jpj, jpk), dshesq(jpi, jpj, jpk), esqt(jpi, jpj), esqf(jpi, jpj), STAT = ierr)
        IF (ierr /= 0) CALL ctl_stop('STOP', 'ldf_dyn_init: failed to allocate Smagorinsky arrays')
        !$ACC KERNELS
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            esqt(ji, jj) = (e1e2t(ji, jj) / (e1t(ji, jj) + e2t(ji, jj))) ** 2
            esqf(ji, jj) = (e1e2f(ji, jj) / (e1f(ji, jj) + e2f(ji, jj))) ** 2
          END DO
        END DO
        !$ACC END KERNELS
      CASE DEFAULT
        CALL ctl_stop('ldf_dyn_init: wrong choice for nn_ahm_ijk_t, the type of space-time variation of ahm')
      END SELECT
      IF (.NOT. l_ldfdyn_time) THEN
        IF (ln_dynldf_lap) THEN
          !$ACC KERNELS
          ahmt(:, :, 1 : jpkm1) = ahmt(:, :, 1 : jpkm1) * tmask(:, :, 1 : jpkm1)
          ahmf(:, :, 1 : jpkm1) = ahmf(:, :, 1 : jpkm1) * fmask(:, :, 1 : jpkm1)
          !$ACC END KERNELS
        ELSE IF (ln_dynldf_blp) THEN
          !$ACC KERNELS
          ahmt(:, :, 1 : jpkm1) = SQRT(ahmt(:, :, 1 : jpkm1)) * tmask(:, :, 1 : jpkm1)
          ahmf(:, :, 1 : jpkm1) = SQRT(ahmf(:, :, 1 : jpkm1)) * fmask(:, :, 1 : jpkm1)
          !$ACC END KERNELS
        END IF
      END IF
    END IF
  END SUBROUTINE ldf_dyn_init
  SUBROUTINE ldf_dyn(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zu2pv2_ij_p1, zu2pv2_ij, zu2pv2_ij_m1, zemax
    REAL(KIND = wp) :: zcmsmag, zstabf_lo, zstabf_up, zdelta, zdb
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    IF (ln_timing) CALL timing_start('ldf_dyn')
    SELECT CASE (nn_ahm_ijk_t)
    CASE (31)
      IF (ln_dynldf_lap) THEN
        DO jk = 1, jpkm1
          !$ACC KERNELS
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zu2pv2_ij = ub(ji, jj, jk) * ub(ji, jj, jk) + vb(ji, jj, jk) * vb(ji, jj, jk)
              zu2pv2_ij_m1 = ub(ji - 1, jj, jk) * ub(ji - 1, jj, jk) + vb(ji, jj - 1, jk) * vb(ji, jj - 1, jk)
              zemax = MAX(e1t(ji, jj), e2t(ji, jj))
              ahmt(ji, jj, jk) = SQRT((zu2pv2_ij + zu2pv2_ij_m1) * r1_288) * zemax * tmask(ji, jj, jk)
            END DO
          END DO
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              zu2pv2_ij_p1 = ub(ji, jj + 1, jk) * ub(ji, jj + 1, jk) + vb(ji + 1, jj, jk) * vb(ji + 1, jj, jk)
              zu2pv2_ij = ub(ji, jj, jk) * ub(ji, jj, jk) + vb(ji, jj, jk) * vb(ji, jj, jk)
              zemax = MAX(e1f(ji, jj), e2f(ji, jj))
              ahmf(ji, jj, jk) = SQRT((zu2pv2_ij + zu2pv2_ij_p1) * r1_288) * zemax * fmask(ji, jj, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      ELSE IF (ln_dynldf_blp) THEN
        DO jk = 1, jpkm1
          !$ACC KERNELS
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zu2pv2_ij = ub(ji, jj, jk) * ub(ji, jj, jk) + vb(ji, jj, jk) * vb(ji, jj, jk)
              zu2pv2_ij_m1 = ub(ji - 1, jj, jk) * ub(ji - 1, jj, jk) + vb(ji, jj - 1, jk) * vb(ji, jj - 1, jk)
              zemax = MAX(e1t(ji, jj), e2t(ji, jj))
              ahmt(ji, jj, jk) = SQRT(SQRT((zu2pv2_ij + zu2pv2_ij_m1) * r1_288) * zemax) * zemax * tmask(ji, jj, jk)
            END DO
          END DO
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              zu2pv2_ij_p1 = ub(ji, jj + 1, jk) * ub(ji, jj + 1, jk) + vb(ji + 1, jj, jk) * vb(ji + 1, jj, jk)
              zu2pv2_ij = ub(ji, jj, jk) * ub(ji, jj, jk) + vb(ji, jj, jk) * vb(ji, jj, jk)
              zemax = MAX(e1f(ji, jj), e2f(ji, jj))
              ahmf(ji, jj, jk) = SQRT(SQRT((zu2pv2_ij + zu2pv2_ij_p1) * r1_288) * zemax) * zemax * fmask(ji, jj, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END IF
      CALL lbc_lnk_multi('ldfdyn', ahmt, 'T', 1., ahmf, 'F', 1.)
    CASE (32)
      IF (ln_dynldf_lap .OR. ln_dynldf_blp) THEN
        CALL profile_psy_data0 % PreStart('ldf_dyn', 'r0', 0, 0)
        zcmsmag = (rn_csmc / rpi) ** 2
        zstabf_lo = rn_minfac * rn_minfac / (2._wp * 4._wp * zcmsmag)
        zstabf_up = rn_maxfac / (4._wp * zcmsmag * 2._wp * rdt)
        IF (ln_dynldf_blp) zstabf_lo = (16._wp / 9._wp) * zstabf_lo
        CALL profile_psy_data0 % PostEnd
        DO jk = 1, jpkm1
          !$ACC KERNELS
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zdb = (ub(ji, jj, jk) * r1_e2u(ji, jj) - ub(ji - 1, jj, jk) * r1_e2u(ji - 1, jj)) * r1_e1t(ji, jj) * e2t(ji, jj) - &
&(vb(ji, jj, jk) * r1_e1v(ji, jj) - vb(ji, jj - 1, jk) * r1_e1v(ji, jj - 1)) * r1_e2t(ji, jj) * e1t(ji, jj)
              dtensq(ji, jj, jk) = zdb * zdb * tmask(ji, jj, jk)
            END DO
          END DO
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              zdb = (ub(ji, jj + 1, jk) * r1_e1u(ji, jj + 1) - ub(ji, jj, jk) * r1_e1u(ji, jj)) * r1_e2f(ji, jj) * e1f(ji, jj) + &
&(vb(ji + 1, jj, jk) * r1_e2v(ji + 1, jj) - vb(ji, jj, jk) * r1_e2v(ji, jj)) * r1_e1f(ji, jj) * e2f(ji, jj)
              dshesq(ji, jj, jk) = zdb * zdb * fmask(ji, jj, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
        CALL lbc_lnk_multi('ldfdyn', dtensq, 'T', 1.)
        DO jk = 1, jpkm1
          !$ACC KERNELS
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zu2pv2_ij = ub(ji, jj, jk) * ub(ji, jj, jk) + vb(ji, jj, jk) * vb(ji, jj, jk)
              zu2pv2_ij_m1 = ub(ji - 1, jj, jk) * ub(ji - 1, jj, jk) + vb(ji, jj - 1, jk) * vb(ji, jj - 1, jk)
              zdelta = zcmsmag * esqt(ji, jj)
              ahmt(ji, jj, jk) = zdelta * SQRT(dtensq(ji, jj, jk) + r1_4 * (dshesq(ji, jj, jk) + dshesq(ji, jj - 1, jk) + &
&dshesq(ji - 1, jj, jk) + dshesq(ji - 1, jj - 1, jk)))
              ahmt(ji, jj, jk) = MAX(ahmt(ji, jj, jk), SQRT((zu2pv2_ij + zu2pv2_ij_m1) * zdelta * zstabf_lo))
              ahmt(ji, jj, jk) = MIN(ahmt(ji, jj, jk), zdelta * zstabf_up)
            END DO
          END DO
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              zu2pv2_ij_p1 = ub(ji, jj + 1, jk) * ub(ji, jj + 1, jk) + vb(ji + 1, jj, jk) * vb(ji + 1, jj, jk)
              zu2pv2_ij = ub(ji, jj, jk) * ub(ji, jj, jk) + vb(ji, jj, jk) * vb(ji, jj, jk)
              zdelta = zcmsmag * esqf(ji, jj)
              ahmf(ji, jj, jk) = zdelta * SQRT(dshesq(ji, jj, jk) + r1_4 * (dtensq(ji, jj, jk) + dtensq(ji, jj + 1, jk) + &
&dtensq(ji + 1, jj, jk) + dtensq(ji + 1, jj + 1, jk)))
              ahmf(ji, jj, jk) = MAX(ahmf(ji, jj, jk), SQRT((zu2pv2_ij + zu2pv2_ij_p1) * zdelta * zstabf_lo))
              ahmf(ji, jj, jk) = MIN(ahmf(ji, jj, jk), zdelta * zstabf_up)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END IF
      IF (ln_dynldf_blp) THEN
        DO jk = 1, jpkm1
          !$ACC KERNELS
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              ahmt(ji, jj, jk) = SQRT(r1_8 * esqt(ji, jj) * ahmt(ji, jj, jk))
            END DO
          END DO
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              ahmf(ji, jj, jk) = SQRT(r1_8 * esqf(ji, jj) * ahmf(ji, jj, jk))
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END IF
      CALL lbc_lnk_multi('ldfdyn', ahmt, 'T', 1., ahmf, 'F', 1.)
    END SELECT
    CALL profile_psy_data1 % PreStart('ldf_dyn', 'r1', 0, 0)
    CALL iom_put("ahmt_2d", ahmt(:, :, 1))
    CALL iom_put("ahmf_2d", ahmf(:, :, 1))
    CALL iom_put("ahmt_3d", ahmt(:, :, :))
    CALL iom_put("ahmf_3d", ahmf(:, :, :))
    IF (ln_timing) CALL timing_stop('ldf_dyn')
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE ldf_dyn
END MODULE ldfdyn