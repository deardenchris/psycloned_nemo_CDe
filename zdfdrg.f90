MODULE zdfdrg
  USE oce
  USE phycst, ONLY: vkarmn
  USE dom_oce
  USE zdf_oce
  USE trd_oce
  USE trddyn
  USE in_out_manager
  USE iom
  USE lbclnk
  USE lib_mpp
  USE prtctl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: zdf_drg
  PUBLIC :: zdf_drg_exp
  PUBLIC :: zdf_drg_init
  LOGICAL :: ln_OFF
  LOGICAL :: ln_lin
  LOGICAL :: ln_non_lin
  LOGICAL :: ln_loglayer
  LOGICAL, PUBLIC :: ln_drgimp
  REAL(KIND = wp) :: rn_Cd0
  REAL(KIND = wp) :: rn_Uc0
  REAL(KIND = wp) :: rn_Cdmax
  REAL(KIND = wp) :: rn_z0
  REAL(KIND = wp) :: rn_ke0
  LOGICAL :: ln_boost
  REAL(KIND = wp) :: rn_boost
  REAL(KIND = wp), PUBLIC :: r_Cdmin_top, r_Cdmax_top, r_z0_top, r_ke0_top
  REAL(KIND = wp), PUBLIC :: r_Cdmin_bot, r_Cdmax_bot, r_z0_bot, r_ke0_bot
  INTEGER :: ndrg
  INTEGER, PARAMETER :: np_OFF = 0
  INTEGER, PARAMETER :: np_lin = 1
  INTEGER, PARAMETER :: np_non_lin = 2
  INTEGER, PARAMETER :: np_loglayer = 3
  LOGICAL, PUBLIC :: l_zdfdrg
  LOGICAL :: l_log_not_linssh
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :), PUBLIC :: rCd0_top, rCd0_bot
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :), PUBLIC :: rCdU_top, rCdU_bot
  CONTAINS
  SUBROUTINE zdf_drg(kt, k_mk, pCdmin, pCdmax, pz0, pke0, pCd0, pCdU)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, DIMENSION(:, :), INTENT(IN) :: k_mk
    REAL(KIND = wp), INTENT(IN) :: pCdmin
    REAL(KIND = wp), INTENT(IN) :: pCdmax
    REAL(KIND = wp), INTENT(IN) :: pz0
    REAL(KIND = wp), INTENT(IN) :: pke0
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pCd0
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: pCdU
    INTEGER :: ji, jj
    INTEGER :: imk
    REAL(KIND = wp) :: zzz, zut, zvt, zcd
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    IF (l_log_not_linssh) THEN
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          imk = k_mk(ji, jj)
          zut = un(ji, jj, imk) + un(ji - 1, jj, imk)
          zvt = vn(ji, jj, imk) + vn(ji, jj - 1, imk)
          zzz = 0.5_wp * e3t_n(ji, jj, imk)
          zcd = (vkarmn / LOG(zzz / pz0)) ** 2
          zcd = pCd0(ji, jj) * MIN(MAX(pCdmin, zcd), pCdmax)
          pCdU(ji, jj) = - zcd * SQRT(0.25 * (zut * zut + zvt * zvt) + pke0)
        END DO
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          imk = k_mk(ji, jj)
          zut = un(ji, jj, imk) + un(ji - 1, jj, imk)
          zvt = vn(ji, jj, imk) + vn(ji, jj - 1, imk)
          pCdU(ji, jj) = - pCd0(ji, jj) * SQRT(0.25 * (zut * zut + zvt * zvt) + pke0)
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data0 % PreStart('zdf_drg', 'r0', 0, 0)
    IF (ln_ctl) CALL prt_ctl(tab2d_1 = pCdU, clinfo1 = ' Cd*U ')
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE zdf_drg
  SUBROUTINE zdf_drg_exp(kt, pub, pvb, pua, pva)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pub, pvb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pua, pva
    INTEGER :: ji, jj
    INTEGER :: ikbu, ikbv
    REAL(KIND = wp) :: zm1_2dt
    REAL(KIND = wp) :: zCdu, zCdv
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ztrdu, ztrdv
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('zdf_drg_exp', 'r0', 0, 0)
    zm1_2dt = - 1._wp / (2._wp * rdt)
    CALL profile_psy_data0 % PostEnd
    IF (l_trddyn) THEN
      ALLOCATE(ztrdu(jpi, jpj, jpk), ztrdv(jpi, jpj, jpk))
      !$ACC KERNELS
      ztrdu(:, :, :) = pua(:, :, :)
      ztrdv(:, :, :) = pva(:, :, :)
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        ikbu = mbku(ji, jj)
        ikbv = mbkv(ji, jj)
        zCdu = 0.5 * (rCdU_bot(ji + 1, jj) + rCdU_bot(ji, jj)) / e3u_n(ji, jj, ikbu)
        zCdv = 0.5 * (rCdU_bot(ji, jj + 1) + rCdU_bot(ji, jj)) / e3v_n(ji, jj, ikbv)
        pua(ji, jj, ikbu) = pua(ji, jj, ikbu) + MAX(zCdu, zm1_2dt) * pub(ji, jj, ikbu)
        pva(ji, jj, ikbv) = pva(ji, jj, ikbv) + MAX(zCdv, zm1_2dt) * pvb(ji, jj, ikbv)
      END DO
    END DO
    !$ACC END KERNELS
    IF (ln_isfcav) THEN
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ikbu = miku(ji, jj)
          ikbv = mikv(ji, jj)
          zCdu = 0.5 * (rCdU_top(ji + 1, jj) + rCdU_top(ji, jj)) / e3u_n(ji, jj, ikbu)
          zCdv = 0.5 * (rCdU_top(ji, jj + 1) + rCdU_top(ji, jj)) / e3v_n(ji, jj, ikbv)
          pua(ji, jj, ikbu) = pua(ji, jj, ikbu) + MAX(zCdu, zm1_2dt) * pub(ji, jj, ikbu)
          pva(ji, jj, ikbv) = pva(ji, jj, ikbv) + MAX(zCdv, zm1_2dt) * pvb(ji, jj, ikbv)
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    IF (l_trddyn) THEN
      !$ACC KERNELS
      ztrdu(:, :, :) = pua(:, :, :) - ztrdu(:, :, :)
      ztrdv(:, :, :) = pva(:, :, :) - ztrdv(:, :, :)
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('zdf_drg_exp', 'r1', 0, 0)
      CALL trd_dyn(ztrdu(:, :, :), ztrdv(:, :, :), jpdyn_bfr, kt)
      DEALLOCATE(ztrdu, ztrdv)
      CALL profile_psy_data1 % PostEnd
    END IF
    CALL profile_psy_data2 % PreStart('zdf_drg_exp', 'r2', 0, 0)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = pua, clinfo1 = ' bfr  - Ua: ', mask1 = umask, tab3d_2 = pva, clinfo2 = ' Va: ', mask2 = vmask, clinfo3 = 'dyn')
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE zdf_drg_exp
  SUBROUTINE zdf_drg_init
    INTEGER :: ji, jj
    INTEGER :: ios, ioptio
    NAMELIST /namdrg/ ln_OFF, ln_lin, ln_non_lin, ln_loglayer, ln_drgimp
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namdrg, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namdrg in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namdrg, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namdrg in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namdrg)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'zdf_drg_init : top and/or bottom drag setting'
      WRITE(numout, *) '~~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist namdrg : top/bottom friction choices'
      WRITE(numout, *) '      free-slip       : Cd = 0                  ln_OFF      = ', ln_OFF
      WRITE(numout, *) '      linear  drag    : Cd = Cd0                ln_lin      = ', ln_lin
      WRITE(numout, *) '      non-linear  drag: Cd = Cd0_nl |U|         ln_non_lin  = ', ln_non_lin
      WRITE(numout, *) '      logarithmic drag: Cd = vkarmn/log(z/z0)   ln_loglayer = ', ln_loglayer
      WRITE(numout, *) '      implicit friction                         ln_drgimp   = ', ln_drgimp
    END IF
    ioptio = 0
    IF (ln_off) THEN
      ndrg = np_off
      ioptio = ioptio + 1
    END IF
    IF (ln_lin) THEN
      ndrg = np_lin
      ioptio = ioptio + 1
    END IF
    IF (ln_non_lin) THEN
      ndrg = np_non_lin
      ioptio = ioptio + 1
    END IF
    IF (ln_loglayer) THEN
      ndrg = np_loglayer
      ioptio = ioptio + 1
    END IF
    IF (ioptio /= 1) CALL ctl_stop('zdf_drg_init: Choose ONE type of drag coef in namdrg')
    ALLOCATE(rCd0_bot(jpi, jpj), rCdU_bot(jpi, jpj))
    CALL drg_init('BOTTOM', mbkt, r_Cdmin_bot, r_Cdmax_bot, r_z0_bot, r_ke0_bot, rCd0_bot, rCdU_bot)
    IF (ln_isfcav) THEN
      ALLOCATE(rCd0_top(jpi, jpj), rCdU_top(jpi, jpj))
      CALL drg_init('TOP   ', mikt, r_Cdmin_top, r_Cdmax_top, r_z0_top, r_ke0_top, rCd0_top, rCdU_top)
    END IF
  END SUBROUTINE zdf_drg_init
  SUBROUTINE drg_init(cd_topbot, k_mk, pCdmin, pCdmax, pz0, pke0, pCd0, pCdU)
    CHARACTER(LEN = 6), INTENT(IN) :: cd_topbot
    INTEGER, DIMENSION(:, :), INTENT(IN) :: k_mk
    REAL(KIND = wp), INTENT(OUT) :: pCdmin, pCdmax
    REAL(KIND = wp), INTENT(OUT) :: pz0
    REAL(KIND = wp), INTENT(OUT) :: pke0
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: pCd0
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: pCdU
    CHARACTER(LEN = 40) :: cl_namdrg, cl_file, cl_varname, cl_namref, cl_namcfg
    INTEGER :: ji, jj
    LOGICAL :: ll_top, ll_bot
    INTEGER :: ios, inum, imk
    REAL(KIND = wp) :: zmsk, zzz, zcd
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zmsk_boost
    NAMELIST /namdrg_top/ rn_Cd0, rn_Uc0, rn_Cdmax, rn_ke0, rn_z0, ln_boost, rn_boost
    NAMELIST /namdrg_bot/ rn_Cd0, rn_Uc0, rn_Cdmax, rn_ke0, rn_z0, ln_boost, rn_boost
    ll_top = .FALSE.
    ll_bot = .FALSE.
    SELECT CASE (cd_topbot)
    CASE ('TOP   ')
      ll_top = .TRUE.
      cl_namdrg = 'namdrg_top'
      cl_namref = 'namdrg_top in reference     namelist'
      cl_namcfg = 'namdrg_top in configuration namelist'
      cl_file = 'tfr_coef.nc'
      cl_varname = 'tfr_coef'
    CASE ('BOTTOM')
      ll_bot = .TRUE.
      cl_namdrg = 'namdrg_bot'
      cl_namref = 'namdrg_bot  in reference     namelist'
      cl_namcfg = 'namdrg_bot  in configuration namelist'
      cl_file = 'bfr_coef.nc'
      cl_varname = 'bfr_coef'
    CASE DEFAULT
      CALL ctl_stop('drg_init: bad value for cd_topbot ')
    END SELECT
    REWIND(UNIT = numnam_ref)
    IF (ll_top) READ(numnam_ref, namdrg_top, IOSTAT = ios, ERR = 901)
    IF (ll_bot) READ(numnam_ref, namdrg_bot, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, TRIM(cl_namref), lwp)
    REWIND(UNIT = numnam_cfg)
    IF (ll_top) READ(numnam_cfg, namdrg_top, IOSTAT = ios, ERR = 902)
    IF (ll_bot) READ(numnam_cfg, namdrg_bot, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, TRIM(cl_namcfg), lwp)
    IF (lwm .AND. ll_top) WRITE(numond, namdrg_top)
    IF (lwm .AND. ll_bot) WRITE(numond, namdrg_bot)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) '   Namelist ', TRIM(cl_namdrg), ' : set ', TRIM(cd_topbot), ' friction parameters'
      WRITE(numout, *) '      drag coefficient                        rn_Cd0   = ', rn_Cd0
      WRITE(numout, *) '      characteristic velocity (linear case)   rn_Uc0   = ', rn_Uc0, ' m/s'
      WRITE(numout, *) '      non-linear drag maximum                 rn_Cdmax = ', rn_Cdmax
      WRITE(numout, *) '      background kinetic energy  (n-l case)   rn_ke0   = ', rn_ke0
      WRITE(numout, *) '      bottom roughness           (n-l case)   rn_z0    = ', rn_z0
      WRITE(numout, *) '      set a regional boost of Cd0             ln_boost = ', ln_boost
      WRITE(numout, *) '         associated boost factor              rn_boost = ', rn_boost
    END IF
    pCdmin = rn_Cd0
    pCdmax = rn_Cdmax
    pz0 = rn_z0
    pke0 = rn_ke0
    IF (ln_boost) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   ==>>>   use a regional boost read in ', TRIM(cl_file), ' file'
      IF (lwp) WRITE(numout, *) '           using enhancement factor of ', rn_boost
      CALL iom_open(TRIM(cl_file), inum)
      CALL iom_get(inum, jpdom_data, TRIM(cl_varname), zmsk_boost, 1)
      CALL iom_close(inum)
      !$ACC KERNELS
      zmsk_boost(:, :) = 1._wp + rn_boost * zmsk_boost(:, :)
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      zmsk_boost(:, :) = 1._wp
      !$ACC END KERNELS
    END IF
    IF (ll_top) zmsk_boost(:, :) = zmsk_boost(:, :) * ssmask(:, :) * (1. - tmask(:, :, 1))
    IF (ll_bot) zmsk_boost(:, :) = zmsk_boost(:, :) * ssmask(:, :)
    SELECT CASE (ndrg)
    CASE (np_OFF)
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   ==>>>   ', TRIM(cd_topbot), ' free-slip, friction set to zero'
      !$ACC KERNELS
      l_zdfdrg = .FALSE.
      pCdU(:, :) = 0._wp
      pCd0(:, :) = 0._wp
      !$ACC END KERNELS
    CASE (np_lin)
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   ==>>>   linear ', TRIM(cd_topbot), ' friction (constant coef = Cd0*Uc0 = ', rn_Cd0 * rn_Uc0, ')'
      !$ACC KERNELS
      l_zdfdrg = .FALSE.
      pCd0(:, :) = rn_Cd0 * zmsk_boost(:, :)
      pCdU(:, :) = - pCd0(:, :) * rn_Uc0
      !$ACC END KERNELS
    CASE (np_non_lin)
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   ==>>>   quadratic ', TRIM(cd_topbot), ' friction (propotional to module of the velocity)'
      IF (lwp) WRITE(numout, *) '   with    a drag coefficient Cd0 = ', rn_Cd0, ', and'
      IF (lwp) WRITE(numout, *) '           a background velocity module of (rn_ke0)^1/2 = ', SQRT(rn_ke0), 'm/s)'
      !$ACC KERNELS
      l_zdfdrg = .TRUE.
      pCd0(:, :) = rn_Cd0 * zmsk_boost(:, :)
      pCdU(:, :) = 0._wp
      !$ACC END KERNELS
    CASE (np_loglayer)
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   ==>>>   quadratic ', TRIM(cd_topbot), ' drag (propotional to module of the velocity)'
      IF (lwp) WRITE(numout, *) '   with   a logarithmic Cd0 formulation Cd0 = ( vkarman log(z/z0) )^2 ,'
      IF (lwp) WRITE(numout, *) '          a background velocity module of (rn_ke0)^1/2 = ', SQRT(pke0), 'm/s), '
      IF (lwp) WRITE(numout, *) '          a logarithmic formulation: a roughness of ', pz0, ' meters,   and '
      IF (lwp) WRITE(numout, *) '          a proportionality factor bounded by min/max values of ', pCdmin, pCdmax
      l_zdfdrg = .TRUE.
      IF (ln_linssh) THEN
        IF (lwp) WRITE(numout, *)
        IF (lwp) WRITE(numout, *) '   N.B.   linear free surface case, Cd0 computed one for all'
        !$ACC KERNELS
        l_log_not_linssh = .FALSE.
        !$ACC loop independent collapse(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            zzz = 0.5_wp * e3t_0(ji, jj, k_mk(ji, jj))
            zcd = (vkarmn / LOG(zzz / rn_z0)) ** 2
            pCd0(ji, jj) = zmsk_boost(ji, jj) * MIN(MAX(rn_Cd0, zcd), rn_Cdmax)
          END DO
        END DO
        !$ACC END KERNELS
      ELSE
        IF (lwp) WRITE(numout, *)
        IF (lwp) WRITE(numout, *) '   N.B.   non-linear free surface case, Cd0 updated at each time-step '
        !$ACC KERNELS
        l_log_not_linssh = .TRUE.
        pCd0(:, :) = zmsk_boost(:, :)
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      pCdU(:, :) = 0._wp
      !$ACC END KERNELS
    CASE DEFAULT
      CALL ctl_stop('drg_init: bad flag value for ndrg ')
    END SELECT
  END SUBROUTINE drg_init
END MODULE zdfdrg