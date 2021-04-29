MODULE ldftra
  USE oce
  USE dom_oce
  USE phycst
  USE ldfslp
  USE ldfc1d_c2d
  USE diaptr
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ldf_tra_init
  PUBLIC :: ldf_tra
  PUBLIC :: ldf_eiv_init
  PUBLIC :: ldf_eiv
  PUBLIC :: ldf_eiv_trp
  PUBLIC :: ldf_eiv_dia
  LOGICAL, PUBLIC :: ln_traldf_OFF
  LOGICAL, PUBLIC :: ln_traldf_lap
  LOGICAL, PUBLIC :: ln_traldf_blp
  LOGICAL, PUBLIC :: ln_traldf_lev
  LOGICAL, PUBLIC :: ln_traldf_hor
  LOGICAL, PUBLIC :: ln_traldf_msc
  INTEGER, PUBLIC :: nn_aht_ijk_t
  REAL(KIND = wp), PUBLIC :: rn_Ud
  REAL(KIND = wp), PUBLIC :: rn_Ld
  LOGICAL, PUBLIC :: ln_ldfeiv
  LOGICAL, PUBLIC :: ln_ldfeiv_dia
  INTEGER, PUBLIC :: nn_aei_ijk_t
  REAL(KIND = wp), PUBLIC :: rn_Ue
  REAL(KIND = wp), PUBLIC :: rn_Le
  INTEGER, PARAMETER, PUBLIC :: np_ERROR = - 10
  INTEGER, PARAMETER, PUBLIC :: np_no_ldf = 00
  INTEGER, PARAMETER, PUBLIC :: np_lap = 10, np_blp = 20
  INTEGER, PARAMETER, PUBLIC :: np_lap_i = 11, np_blp_i = 21
  INTEGER, PARAMETER, PUBLIC :: np_lap_it = 12, np_blp_it = 22
  INTEGER, PUBLIC :: nldf_tra = 0
  LOGICAL, PUBLIC :: l_ldftra_time = .FALSE.
  LOGICAL, PUBLIC :: l_ldfeiv_time = .FALSE.
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: ahtu, ahtv
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: aeiu, aeiv
  REAL(KIND = wp) :: aht0, aei0
  REAL(KIND = wp) :: r1_2 = 0.5_wp
  REAL(KIND = wp) :: r1_4 = 0.25_wp
  REAL(KIND = wp) :: r1_12 = 1._wp / 12._wp
  CONTAINS
  SUBROUTINE ldf_tra_init
    INTEGER :: jk
    INTEGER :: ioptio, ierr, inum, ios, inn
    REAL(KIND = wp) :: zah_max, zUfac
    CHARACTER(LEN = 5) :: cl_Units
    NAMELIST /namtra_ldf/ ln_traldf_OFF, ln_traldf_lap, ln_traldf_blp, ln_traldf_lev, ln_traldf_hor, ln_traldf_triad, ln_traldf_iso, ln_traldf_msc, rn_slpmax, ln_triad_iso, ln_botmix_triad, rn_sw_triad, nn_aht_ijk_t, rn_Ud, rn_Ld
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ldf_tra_init : lateral tracer diffusion'
      WRITE(numout, *) '~~~~~~~~~~~~ '
    END IF
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namtra_ldf, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namtra_ldf in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namtra_ldf, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namtra_ldf in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namtra_ldf)
    IF (lwp) THEN
      WRITE(numout, *) '   Namelist : namtra_ldf --- lateral mixing parameters (type, direction, coefficients)'
      WRITE(numout, *) '      type :'
      WRITE(numout, *) '         no explicit diffusion                   ln_traldf_OFF   = ', ln_traldf_OFF
      WRITE(numout, *) '         laplacian operator                      ln_traldf_lap   = ', ln_traldf_lap
      WRITE(numout, *) '         bilaplacian operator                    ln_traldf_blp   = ', ln_traldf_blp
      WRITE(numout, *) '      direction of action :'
      WRITE(numout, *) '         iso-level                               ln_traldf_lev   = ', ln_traldf_lev
      WRITE(numout, *) '         horizontal (geopotential)               ln_traldf_hor   = ', ln_traldf_hor
      WRITE(numout, *) '         iso-neutral Madec operator              ln_traldf_iso   = ', ln_traldf_iso
      WRITE(numout, *) '         iso-neutral triad operator              ln_traldf_triad = ', ln_traldf_triad
      WRITE(numout, *) '            use the Method of Stab. Correction   ln_traldf_msc   = ', ln_traldf_msc
      WRITE(numout, *) '            maximum isoppycnal slope             rn_slpmax       = ', rn_slpmax
      WRITE(numout, *) '            pure lateral mixing in ML            ln_triad_iso    = ', ln_triad_iso
      WRITE(numout, *) '            switching triad or not               rn_sw_triad     = ', rn_sw_triad
      WRITE(numout, *) '            lateral mixing on bottom             ln_botmix_triad = ', ln_botmix_triad
      WRITE(numout, *) '      coefficients :'
      WRITE(numout, *) '         type of time-space variation            nn_aht_ijk_t    = ', nn_aht_ijk_t
      WRITE(numout, *) '            lateral diffusive velocity (if cst)  rn_Ud           = ', rn_Ud, ' m/s'
      WRITE(numout, *) '            lateral diffusive length   (if cst)  rn_Ld           = ', rn_Ld, ' m'
    END IF
    nldf_tra = np_ERROR
    ioptio = 0
    IF (ln_traldf_off) THEN
      nldf_tra = np_no_ldf
      ioptio = ioptio + 1
    END IF
    IF (ln_traldf_lap) THEN
      ioptio = ioptio + 1
    END IF
    IF (ln_traldf_blp) THEN
      ioptio = ioptio + 1
    END IF
    IF (ioptio /= 1) CALL ctl_stop('tra_ldf_init: use ONE of the 3 operator options (NONE/lap/blp)')
    IF (.NOT. ln_traldf_OFF) THEN
      ioptio = 0
      IF (ln_traldf_lev) ioptio = ioptio + 1
      IF (ln_traldf_hor) ioptio = ioptio + 1
      IF (ln_traldf_iso) ioptio = ioptio + 1
      IF (ln_traldf_triad) ioptio = ioptio + 1
      IF (ioptio /= 1) CALL ctl_stop('tra_ldf_init: use ONE direction (level/hor/iso/triad)')
      ierr = 0
      IF (ln_traldf_lap) THEN
        IF (ln_zco) THEN
          IF (ln_traldf_lev) nldf_tra = np_lap
          IF (ln_traldf_hor) nldf_tra = np_lap
          IF (ln_traldf_iso) nldf_tra = np_lap_i
          IF (ln_traldf_triad) nldf_tra = np_lap_it
        END IF
        IF (ln_zps) THEN
          IF (ln_traldf_lev) ierr = 1
          IF (ln_traldf_hor) nldf_tra = np_lap
          IF (ln_traldf_iso) nldf_tra = np_lap_i
          IF (ln_traldf_triad) nldf_tra = np_lap_it
        END IF
        IF (ln_sco) THEN
          IF (ln_traldf_lev) nldf_tra = np_lap
          IF (ln_traldf_hor) nldf_tra = np_lap_i
          IF (ln_traldf_iso) nldf_tra = np_lap_i
          IF (ln_traldf_triad) nldf_tra = np_lap_it
        END IF
      END IF
      IF (ln_traldf_blp) THEN
        IF (ln_zco) THEN
          IF (ln_traldf_lev) nldf_tra = np_blp
          IF (ln_traldf_hor) nldf_tra = np_blp
          IF (ln_traldf_iso) nldf_tra = np_blp_i
          IF (ln_traldf_triad) nldf_tra = np_blp_it
        END IF
        IF (ln_zps) THEN
          IF (ln_traldf_lev) ierr = 1
          IF (ln_traldf_hor) nldf_tra = np_blp
          IF (ln_traldf_iso) nldf_tra = np_blp_i
          IF (ln_traldf_triad) nldf_tra = np_blp_it
        END IF
        IF (ln_sco) THEN
          IF (ln_traldf_lev) nldf_tra = np_blp
          IF (ln_traldf_hor) nldf_tra = np_blp_it
          IF (ln_traldf_iso) nldf_tra = np_blp_i
          IF (ln_traldf_triad) nldf_tra = np_blp_it
        END IF
      END IF
      IF (ierr == 1) CALL ctl_stop('iso-level in z-partial step, not allowed')
    END IF
    IF (ln_ldfeiv .AND. .NOT. (ln_traldf_iso .OR. ln_traldf_triad)) CALL ctl_stop('ln_ldfeiv=T requires iso-neutral laplacian diffusion')
    IF (ln_isfcav .AND. ln_traldf_triad) CALL ctl_stop(' ice shelf cavity and traldf_triad not tested')
    IF (nldf_tra == np_lap_i .OR. nldf_tra == np_lap_it .OR. nldf_tra == np_blp_i .OR. nldf_tra == np_blp_it) l_ldfslp = .TRUE.
    IF (ln_traldf_blp .AND. (ln_traldf_iso .OR. ln_traldf_triad)) THEN
      IF (.NOT. ln_traldf_msc) CALL ctl_stop('tra_ldf_init: iso-neutral bilaplacian requires ln_traldf_msc=.true.')
    END IF
    IF (lwp) THEN
      WRITE(numout, *)
      SELECT CASE (nldf_tra)
      CASE (np_no_ldf)
        WRITE(numout, *) '   ==>>>   NO lateral diffusion'
      CASE (np_lap)
        WRITE(numout, *) '   ==>>>   laplacian iso-level operator'
      CASE (np_lap_i)
        WRITE(numout, *) '   ==>>>   Rotated laplacian operator (standard)'
      CASE (np_lap_it)
        WRITE(numout, *) '   ==>>>   Rotated laplacian operator (triad)'
      CASE (np_blp)
        WRITE(numout, *) '   ==>>>   bilaplacian iso-level operator'
      CASE (np_blp_i)
        WRITE(numout, *) '   ==>>>   Rotated bilaplacian operator (standard)'
      CASE (np_blp_it)
        WRITE(numout, *) '   ==>>>   Rotated bilaplacian operator (triad)'
      END SELECT
      WRITE(numout, *)
    END IF
    l_ldftra_time = .FALSE.
    IF (ln_traldf_OFF) THEN
      IF (lwp) WRITE(numout, *) '   ==>>>   No diffusive operator selected. ahtu and ahtv are not allocated'
      RETURN
    ELSE
      ALLOCATE(ahtu(jpi, jpj, jpk), ahtv(jpi, jpj, jpk), STAT = ierr)
      IF (ierr /= 0) CALL ctl_stop('STOP', 'ldf_tra_init: failed to allocate arrays')
      !$ACC KERNELS
      ahtu(:, :, jpk) = 0._wp
      ahtv(:, :, jpk) = 0._wp
      !$ACC END KERNELS
      IF (ln_traldf_lap) THEN
        zufac = r1_2 * rn_ud
        inn = 1
        cl_units = ' m2/s'
      ELSE IF (ln_traldf_blp) THEN
        zufac = r1_12 * rn_ud
        inn = 3
        cl_units = ' m4/s'
      END IF
      aht0 = zUfac * rn_Ld ** inn
      zah_max = zUfac * (ra * rad) ** inn
      SELECT CASE (nn_aht_ijk_t)
      CASE (0)
        IF (lwp) WRITE(numout, *) '   ==>>>   eddy diffusivity = constant = ', aht0, cl_Units
        !$ACC KERNELS
        ahtu(:, :, 1 : jpkm1) = aht0
        ahtv(:, :, 1 : jpkm1) = aht0
        !$ACC END KERNELS
      CASE (10)
        IF (lwp) WRITE(numout, *) '   ==>>>   eddy diffusivity = F( depth )'
        IF (lwp) WRITE(numout, *) '           surface eddy diffusivity = constant = ', aht0, cl_Units
        !$ACC KERNELS
        ahtu(:, :, 1) = aht0
        ahtv(:, :, 1) = aht0
        !$ACC END KERNELS
        CALL ldf_c1d('TRA', ahtu(:, :, 1), ahtv(:, :, 1), ahtu, ahtv)
      CASE (- 20)
        IF (lwp) WRITE(numout, *) '   ==>>>   eddy diffusivity = F(i,j) read in eddy_diffusivity.nc file'
        CALL iom_open('eddy_diffusivity_2D.nc', inum)
        CALL iom_get(inum, jpdom_data, 'ahtu_2D', ahtu(:, :, 1))
        CALL iom_get(inum, jpdom_data, 'ahtv_2D', ahtv(:, :, 1))
        CALL iom_close(inum)
        !$ACC KERNELS
        DO jk = 2, jpkm1
          ahtu(:, :, jk) = ahtu(:, :, 1)
          ahtv(:, :, jk) = ahtv(:, :, 1)
        END DO
        !$ACC END KERNELS
      CASE (20)
        IF (lwp) WRITE(numout, *) '   ==>>>   eddy diffusivity = F( e1, e2 ) or F( e1^3, e2^3 ) (lap or blp case)'
        IF (lwp) WRITE(numout, *) '           using a fixed diffusive velocity = ', rn_Ud, ' m/s   and   Ld = Max(e1,e2)'
        IF (lwp) WRITE(numout, *) '           maximum reachable coefficient (at the Equator) = ', zah_max, cl_Units, '  for e1=1°)'
        CALL ldf_c2d('TRA', zUfac, inn, ahtu, ahtv)
      CASE (21)
        IF (lwp) WRITE(numout, *) '   ==>>>   eddy diffusivity = F( latitude, longitude, time )'
        IF (lwp) WRITE(numout, *) '                            = F( growth rate of baroclinic instability )'
        IF (lwp) WRITE(numout, *) '            min value = 0.2 * aht0 (with aht0= 1/2 rn_Ud*rn_Ld)'
        IF (lwp) WRITE(numout, *) '            max value =       aei0 (with aei0=1/2 rn_Ue*Le  increased to aht0 within 20N-20S'
        l_ldftra_time = .TRUE.
        IF (ln_traldf_blp) CALL ctl_stop('ldf_tra_init: aht=F( growth rate of baroc. insta .)', '              incompatible with bilaplacian operator')
      CASE (- 30)
        IF (lwp) WRITE(numout, *) '   ==>>>   eddy diffusivity = F(i,j,k) read in eddy_diffusivity.nc file'
        CALL iom_open('eddy_diffusivity_3D.nc', inum)
        CALL iom_get(inum, jpdom_data, 'ahtu_3D', ahtu)
        CALL iom_get(inum, jpdom_data, 'ahtv_3D', ahtv)
        CALL iom_close(inum)
      CASE (30)
        IF (lwp) WRITE(numout, *) '   ==>>>   eddy diffusivity = F( latitude, longitude, depth )'
        IF (lwp) WRITE(numout, *) '           using a fixed diffusive velocity = ', rn_Ud, ' m/s   and   Ld = Max(e1,e2)'
        IF (lwp) WRITE(numout, *) '           maximum reachable coefficient (at the Equator) = ', zah_max, cl_Units, '  for e1=1°)'
        CALL ldf_c2d('TRA', zUfac, inn, ahtu, ahtv)
        CALL ldf_c1d('TRA', ahtu(:, :, 1), ahtv(:, :, 1), ahtu, ahtv)
      CASE (31)
        IF (lwp) WRITE(numout, *) '   ==>>>   eddy diffusivity = F( latitude, longitude, depth , time )'
        IF (lwp) WRITE(numout, *) '           proportional to the velocity : 1/2 |u|e or 1/12 |u|e^3'
        l_ldftra_time = .TRUE.
      CASE DEFAULT
        CALL ctl_stop('ldf_tra_init: wrong choice for nn_aht_ijk_t, the type of space-time variation of aht')
      END SELECT
      IF (.NOT. l_ldftra_time) THEN
        IF (ln_traldf_lap) THEN
          !$ACC KERNELS
          ahtu(:, :, 1 : jpkm1) = ahtu(:, :, 1 : jpkm1) * umask(:, :, 1 : jpkm1)
          ahtv(:, :, 1 : jpkm1) = ahtv(:, :, 1 : jpkm1) * vmask(:, :, 1 : jpkm1)
          !$ACC END KERNELS
        ELSE IF (ln_traldf_blp) THEN
          !$ACC KERNELS
          ahtu(:, :, 1 : jpkm1) = SQRT(ahtu(:, :, 1 : jpkm1)) * umask(:, :, 1 : jpkm1)
          ahtv(:, :, 1 : jpkm1) = SQRT(ahtv(:, :, 1 : jpkm1)) * vmask(:, :, 1 : jpkm1)
          !$ACC END KERNELS
        END IF
      END IF
    END IF
  END SUBROUTINE ldf_tra_init
  SUBROUTINE ldf_tra(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zaht, zahf, zaht_min, zDaht, z1_f20
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('ldf_tra', 'r0', 0, 0)
    IF (ln_ldfeiv .AND. nn_aei_ijk_t == 21) THEN
      CALL ldf_eiv(kt, aei0, aeiu, aeiv)
    END IF
    CALL profile_psy_data0 % PostEnd
    SELECT CASE (nn_aht_ijk_t)
    CASE (21)
      IF (ln_ldfeiv .AND. nn_aei_ijk_t == 21) THEN
        !$ACC KERNELS
        ahtu(:, :, 1) = aeiu(:, :, 1)
        ahtv(:, :, 1) = aeiv(:, :, 1)
        !$ACC END KERNELS
      ELSE
        CALL profile_psy_data1 % PreStart('ldf_tra', 'r1', 0, 0)
        CALL ldf_eiv(kt, aht0, ahtu, ahtv)
        CALL profile_psy_data1 % PostEnd
      END IF
      !$ACC KERNELS
      z1_f20 = 1._wp / (2._wp * omega * SIN(rad * 20._wp))
      zaht_min = 0.2_wp * aht0
      zDaht = aht0 - zaht_min
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          zaht = (1._wp - MIN(1._wp, ABS(ff_t(ji, jj) * z1_f20))) * zDaht
          zahf = (1._wp - MIN(1._wp, ABS(ff_f(ji, jj) * z1_f20))) * zDaht
          ahtu(ji, jj, 1) = (MAX(zaht_min, ahtu(ji, jj, 1)) + zaht)
          ahtv(ji, jj, 1) = (MAX(zaht_min, ahtv(ji, jj, 1)) + zahf)
        END DO
      END DO
      DO jk = 1, jpkm1
        ahtu(:, :, jk) = ahtu(:, :, 1) * umask(:, :, jk)
        ahtv(:, :, jk) = ahtv(:, :, 1) * vmask(:, :, jk)
      END DO
      !$ACC END KERNELS
    CASE (31)
      IF (ln_traldf_lap) THEN
        !$ACC KERNELS
        DO jk = 1, jpkm1
          ahtu(:, :, jk) = ABS(ub(:, :, jk)) * e1u(:, :) * r1_12
          ahtv(:, :, jk) = ABS(vb(:, :, jk)) * e2v(:, :) * r1_12
        END DO
        !$ACC END KERNELS
      ELSE IF (ln_traldf_blp) THEN
        !$ACC KERNELS
        DO jk = 1, jpkm1
          ahtu(:, :, jk) = SQRT(ABS(ub(:, :, jk)) * e1u(:, :) * r1_12) * e1u(:, :)
          ahtv(:, :, jk) = SQRT(ABS(vb(:, :, jk)) * e2v(:, :) * r1_12) * e2v(:, :)
        END DO
        !$ACC END KERNELS
      END IF
    END SELECT
    CALL profile_psy_data2 % PreStart('ldf_tra', 'r2', 0, 0)
    CALL iom_put("ahtu_2d", ahtu(:, :, 1))
    CALL iom_put("ahtv_2d", ahtv(:, :, 1))
    CALL iom_put("ahtu_3d", ahtu(:, :, :))
    CALL iom_put("ahtv_3d", ahtv(:, :, :))
    IF (ln_ldfeiv) THEN
      CALL iom_put("aeiu_2d", aeiu(:, :, 1))
      CALL iom_put("aeiv_2d", aeiv(:, :, 1))
      CALL iom_put("aeiu_3d", aeiu(:, :, :))
      CALL iom_put("aeiv_3d", aeiv(:, :, :))
    END IF
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE ldf_tra
  SUBROUTINE ldf_eiv_init
    INTEGER :: jk
    INTEGER :: ierr, inum, ios, inn
    REAL(KIND = wp) :: zah_max, zUfac
    NAMELIST /namtra_eiv/ ln_ldfeiv, ln_ldfeiv_dia, nn_aei_ijk_t, rn_Ue, rn_Le
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ldf_eiv_init : eddy induced velocity parametrization'
      WRITE(numout, *) '~~~~~~~~~~~~ '
    END IF
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namtra_eiv, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namtra_eiv in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namtra_eiv, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namtra_eiv in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namtra_eiv)
    IF (lwp) THEN
      WRITE(numout, *) '   Namelist namtra_eiv : '
      WRITE(numout, *) '      Eddy Induced Velocity (eiv) param.         ln_ldfeiv     = ', ln_ldfeiv
      WRITE(numout, *) '      eiv streamfunction & velocity diag.        ln_ldfeiv_dia = ', ln_ldfeiv_dia
      WRITE(numout, *) '      coefficients :'
      WRITE(numout, *) '         type of time-space variation            nn_aei_ijk_t  = ', nn_aht_ijk_t
      WRITE(numout, *) '         lateral diffusive velocity (if cst)     rn_Ue         = ', rn_Ue, ' m/s'
      WRITE(numout, *) '         lateral diffusive length   (if cst)     rn_Le         = ', rn_Le, ' m'
      WRITE(numout, *)
    END IF
    l_ldfeiv_time = .FALSE.
    IF (.NOT. ln_ldfeiv) THEN
      IF (lwp) WRITE(numout, *) '   ==>>>   eddy induced velocity param is NOT used'
      ln_ldfeiv_dia = .FALSE.
    ELSE
      IF (lwp) WRITE(numout, *) '   ==>>>   use eddy induced velocity parametrization'
      IF (lwp) WRITE(numout, *)
      IF (ln_traldf_blp) CALL ctl_stop('ldf_eiv_init: eddy induced velocity ONLY with laplacian diffusivity')
      ALLOCATE(aeiu(jpi, jpj, jpk), aeiv(jpi, jpj, jpk), STAT = ierr)
      IF (ierr /= 0) CALL ctl_stop('STOP', 'ldf_eiv: failed to allocate arrays')
      !$ACC KERNELS
      aeiu(:, :, jpk) = 0._wp
      aeiv(:, :, jpk) = 0._wp
      zUfac = r1_2 * rn_Ue
      inn = 1
      aei0 = zUfac * rn_Le ** inn
      zah_max = zUfac * (ra * rad) ** inn
      !$ACC END KERNELS
      SELECT CASE (nn_aei_ijk_t)
      CASE (0)
        IF (lwp) WRITE(numout, *) '   ==>>>   eddy induced velocity coef. = constant = ', aei0, ' m2/s'
        !$ACC KERNELS
        aeiu(:, :, 1 : jpkm1) = aei0
        aeiv(:, :, 1 : jpkm1) = aei0
        !$ACC END KERNELS
      CASE (10)
        IF (lwp) WRITE(numout, *) '   ==>>>   eddy induced velocity coef. = F( depth )'
        IF (lwp) WRITE(numout, *) '           surface eddy diffusivity = constant = ', aht0, ' m2/s'
        !$ACC KERNELS
        aeiu(:, :, 1) = aei0
        aeiv(:, :, 1) = aei0
        !$ACC END KERNELS
        CALL ldf_c1d('TRA', aeiu(:, :, 1), aeiv(:, :, 1), aeiu, aeiv)
      CASE (- 20)
        IF (lwp) WRITE(numout, *) '   ==>>>   eddy induced velocity coef. = F(i,j) read in eddy_diffusivity_2D.nc file'
        CALL iom_open('eddy_induced_velocity_2D.nc', inum)
        CALL iom_get(inum, jpdom_data, 'aeiu', aeiu(:, :, 1))
        CALL iom_get(inum, jpdom_data, 'aeiv', aeiv(:, :, 1))
        CALL iom_close(inum)
        !$ACC KERNELS
        DO jk = 2, jpkm1
          aeiu(:, :, jk) = aeiu(:, :, 1)
          aeiv(:, :, jk) = aeiv(:, :, 1)
        END DO
        !$ACC END KERNELS
      CASE (20)
        IF (lwp) WRITE(numout, *) '   ==>>>   eddy induced velocity coef. = F( e1, e2 )'
        IF (lwp) WRITE(numout, *) '           using a fixed diffusive velocity = ', rn_Ue, ' m/s   and   Le = Max(e1,e2)'
        IF (lwp) WRITE(numout, *) '           maximum reachable coefficient (at the Equator) = ', zah_max, ' m2/s   for e1=1°)'
        CALL ldf_c2d('TRA', zUfac, inn, aeiu, aeiv)
      CASE (21)
        IF (lwp) WRITE(numout, *) '   ==>>>   eddy induced velocity coef. = F( latitude, longitude, time )'
        IF (lwp) WRITE(numout, *) '                                       = F( growth rate of baroclinic instability )'
        IF (lwp) WRITE(numout, *) '           maximum allowed value: aei0 = ', aei0, ' m2/s'
        l_ldfeiv_time = .TRUE.
      CASE (- 30)
        IF (lwp) WRITE(numout, *) '   ==>>>   eddy induced velocity coef. = F(i,j,k) read in eddy_diffusivity_3D.nc file'
        CALL iom_open('eddy_induced_velocity_3D.nc', inum)
        CALL iom_get(inum, jpdom_data, 'aeiu', aeiu)
        CALL iom_get(inum, jpdom_data, 'aeiv', aeiv)
        CALL iom_close(inum)
      CASE (30)
        IF (lwp) WRITE(numout, *) '   ==>>>   eddy induced velocity coef. = F( latitude, longitude, depth )'
        CALL ldf_c2d('TRA', zUfac, inn, aeiu, aeiv)
        CALL ldf_c1d('TRA', aeiu(:, :, 1), aeiv(:, :, 1), aeiu, aeiv)
      CASE DEFAULT
        CALL ctl_stop('ldf_tra_init: wrong choice for nn_aei_ijk_t, the type of space-time variation of aei')
      END SELECT
      IF (.NOT. l_ldfeiv_time) THEN
        !$ACC KERNELS
        DO jk = 1, jpkm1
          aeiu(:, :, jk) = aeiu(:, :, jk) * umask(:, :, jk)
          ahtv(:, :, jk) = ahtv(:, :, jk) * vmask(:, :, jk)
        END DO
        !$ACC END KERNELS
      END IF
    END IF
  END SUBROUTINE ldf_eiv_init
  SUBROUTINE ldf_eiv(kt, paei0, paeiu, paeiv)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), INTENT(INOUT) :: paei0
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: paeiu, paeiv
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zfw, ze3w, zn2, z1_f20, zaht, zaht_min, zzaei
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zn, zah, zhw, zRo, zaeiw
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    !$ACC KERNELS
    zn(:, :) = 0._wp
    zhw(:, :) = 5._wp
    zah(:, :) = 0._wp
    zRo(:, :) = 0._wp
    !$ACC END KERNELS
    IF (ln_traldf_triad) THEN
      !$ACC KERNELS
      DO jk = 1, jpk
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zn2 = MAX(rn2b(ji, jj, jk), 0._wp)
            zn(ji, jj) = zn(ji, jj) + SQRT(zn2) * e3w_n(ji, jj, jk)
            ze3w = e3w_n(ji, jj, jk) * tmask(ji, jj, jk)
            zah(ji, jj) = zah(ji, jj) + zn2 * wslp2(ji, jj, jk) * ze3w
            zhw(ji, jj) = zhw(ji, jj) + ze3w
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      DO jk = 1, jpk
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zn2 = MAX(rn2b(ji, jj, jk), 0._wp)
            zn(ji, jj) = zn(ji, jj) + SQRT(zn2) * e3w_n(ji, jj, jk)
            ze3w = e3w_n(ji, jj, jk) * tmask(ji, jj, jk)
            zah(ji, jj) = zah(ji, jj) + zn2 * (wslpi(ji, jj, jk) * wslpi(ji, jj, jk) + wslpj(ji, jj, jk) * wslpj(ji, jj, jk)) * ze3w
            zhw(ji, jj) = zhw(ji, jj) + ze3w
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zfw = MAX(ABS(2. * omega * SIN(rad * gphit(ji, jj))), 1.E-10)
        zRo(ji, jj) = MAX(2.E3, MIN(.4 * zn(ji, jj) / zfw, 40.E3))
        zaeiw(ji, jj) = zRo(ji, jj) * zRo(ji, jj) * SQRT(zah(ji, jj) / zhw(ji, jj)) * tmask(ji, jj, 1)
      END DO
    END DO
    z1_f20 = 1._wp / (2._wp * omega * SIN(rad * 20._wp))
    !$ACC loop independent collapse(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zzaei = MIN(1._wp, ABS(ff_t(ji, jj) * z1_f20)) * zaeiw(ji, jj)
        zaeiw(ji, jj) = MIN(zzaei, paei0)
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('ldf_eiv', 'r0', 0, 0)
    CALL lbc_lnk('ldftra', zaeiw(:, :), 'W', 1.)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        paeiu(ji, jj, 1) = 0.5_wp * (zaeiw(ji, jj) + zaeiw(ji + 1, jj)) * umask(ji, jj, 1)
        paeiv(ji, jj, 1) = 0.5_wp * (zaeiw(ji, jj) + zaeiw(ji, jj + 1)) * vmask(ji, jj, 1)
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('ldf_eiv', 'r1', 0, 0)
    CALL lbc_lnk_multi('ldftra', paeiu(:, :, 1), 'U', 1., paeiv(:, :, 1), 'V', 1.)
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    DO jk = 2, jpkm1
      paeiu(:, :, jk) = paeiu(:, :, 1) * umask(:, :, jk)
      paeiv(:, :, jk) = paeiv(:, :, 1) * vmask(:, :, jk)
    END DO
    !$ACC END KERNELS
  END SUBROUTINE ldf_eiv
  SUBROUTINE ldf_eiv_trp(kt, kit000, pun, pvn, pwn, cdtype)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kit000
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pun
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pvn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pwn
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zuwk, zuwk1, zuwi, zuwi1
    REAL(KIND = wp) :: zvwk, zvwk1, zvwj, zvwj1
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zpsi_uw, zpsi_vw
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('ldf_eiv_trp', 'r0', 0, 0)
    IF (kt == kit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'ldf_eiv_trp : eddy induced advection on ', cdtype, ' :'
      IF (lwp) WRITE(numout, *) '~~~~~~~~~~~   add to velocity fields the eiv component'
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zpsi_uw(:, :, 1) = 0._wp
    zpsi_vw(:, :, 1) = 0._wp
    zpsi_uw(:, :, jpk) = 0._wp
    zpsi_vw(:, :, jpk) = 0._wp
    DO jk = 2, jpkm1
      !$ACC loop independent collapse(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zpsi_uw(ji, jj, jk) = - r1_4 * e2u(ji, jj) * (wslpi(ji, jj, jk) + wslpi(ji + 1, jj, jk)) * (aeiu(ji, jj, jk - 1) + aeiu(ji, jj, jk)) * umask(ji, jj, jk)
          zpsi_vw(ji, jj, jk) = - r1_4 * e1v(ji, jj) * (wslpj(ji, jj, jk) + wslpj(ji, jj + 1, jk)) * (aeiv(ji, jj, jk - 1) + aeiv(ji, jj, jk)) * vmask(ji, jj, jk)
        END DO
      END DO
    END DO
    DO jk = 1, jpkm1
      !$ACC loop independent collapse(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          pun(ji, jj, jk) = pun(ji, jj, jk) - (zpsi_uw(ji, jj, jk) - zpsi_uw(ji, jj, jk + 1))
          pvn(ji, jj, jk) = pvn(ji, jj, jk) - (zpsi_vw(ji, jj, jk) - zpsi_vw(ji, jj, jk + 1))
        END DO
      END DO
    END DO
    DO jk = 1, jpkm1
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          pwn(ji, jj, jk) = pwn(ji, jj, jk) + (zpsi_uw(ji, jj, jk) - zpsi_uw(ji - 1, jj, jk) + zpsi_vw(ji, jj, jk) - zpsi_vw(ji, jj - 1, jk))
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('ldf_eiv_trp', 'r1', 0, 0)
    IF (ln_ldfeiv_dia .AND. cdtype == 'TRA') CALL ldf_eiv_dia(zpsi_uw, zpsi_vw)
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE ldf_eiv_trp
  SUBROUTINE ldf_eiv_dia(psi_uw, psi_vw)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: psi_uw, psi_vw
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zztmp
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zw2d
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zw3d
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
    CALL profile_psy_data0 % PreStart('ldf_eiv_dia', 'r0', 0, 0)
    CALL lbc_lnk_multi('ldftra', psi_uw, 'U', - 1., psi_vw, 'V', - 1.)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zw3d(:, :, jpk) = 0._wp
    !$ACC END KERNELS
    IF (iom_use("uoce_eiv")) THEN
      !$ACC KERNELS
      DO jk = 1, jpkm1
        zw3d(:, :, jk) = (psi_uw(:, :, jk + 1) - psi_uw(:, :, jk)) / (e2u(:, :) * e3u_n(:, :, jk))
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('ldf_eiv_dia', 'r1', 0, 0)
      CALL iom_put("uoce_eiv", zw3d)
      CALL profile_psy_data1 % PostEnd
    END IF
    IF (iom_use("voce_eiv")) THEN
      !$ACC KERNELS
      DO jk = 1, jpkm1
        zw3d(:, :, jk) = (psi_vw(:, :, jk + 1) - psi_vw(:, :, jk)) / (e1v(:, :) * e3v_n(:, :, jk))
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data2 % PreStart('ldf_eiv_dia', 'r2', 0, 0)
      CALL iom_put("voce_eiv", zw3d)
      CALL profile_psy_data2 % PostEnd
    END IF
    IF (iom_use("woce_eiv")) THEN
      !$ACC KERNELS
      DO jk = 1, jpkm1
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zw3d(ji, jj, jk) = (psi_vw(ji, jj, jk) - psi_vw(ji, jj - 1, jk) + psi_uw(ji, jj, jk) - psi_uw(ji - 1, jj, jk)) / e1e2t(ji, jj)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data3 % PreStart('ldf_eiv_dia', 'r3', 0, 0)
      CALL lbc_lnk('ldftra', zw3d, 'T', 1.)
      CALL iom_put("woce_eiv", zw3d)
      CALL profile_psy_data3 % PostEnd
    END IF
    CALL profile_psy_data4 % PreStart('ldf_eiv_dia', 'r4', 0, 0)
    zztmp = 0.5_wp * rau0 * rcp
    CALL profile_psy_data4 % PostEnd
    IF (iom_use('ueiv_heattr') .OR. iom_use('ueiv_heattr3d')) THEN
      !$ACC KERNELS
      zw2d(:, :) = 0._wp
      zw3d(:, :, :) = 0._wp
      DO jk = 1, jpkm1
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zw3d(ji, jj, jk) = zw3d(ji, jj, jk) + (psi_uw(ji, jj, jk + 1) - psi_uw(ji, jj, jk)) * (tsn(ji, jj, jk, jp_tem) + tsn(ji + 1, jj, jk, jp_tem))
            zw2d(ji, jj) = zw2d(ji, jj) + zw3d(ji, jj, jk)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data5 % PreStart('ldf_eiv_dia', 'r5', 0, 0)
      CALL lbc_lnk('ldftra', zw2d, 'U', - 1.)
      CALL lbc_lnk('ldftra', zw3d, 'U', - 1.)
      CALL iom_put("ueiv_heattr", zztmp * zw2d)
      CALL iom_put("ueiv_heattr3d", zztmp * zw3d)
      CALL profile_psy_data5 % PostEnd
    END IF
    IF (iom_use('veiv_heattr')) THEN
      !$ACC KERNELS
      zw2d(:, :) = 0._wp
      zw3d(:, :, :) = 0._wp
      DO jk = 1, jpkm1
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zw3d(ji, jj, jk) = zw3d(ji, jj, jk) + (psi_vw(ji, jj, jk + 1) - psi_vw(ji, jj, jk)) * (tsn(ji, jj, jk, jp_tem) + tsn(ji, jj + 1, jk, jp_tem))
            zw2d(ji, jj) = zw2d(ji, jj) + zw3d(ji, jj, jk)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data6 % PreStart('ldf_eiv_dia', 'r6', 0, 0)
      CALL lbc_lnk('ldftra', zw2d, 'V', - 1.)
      CALL iom_put("veiv_heattr", zztmp * zw2d)
      CALL iom_put("veiv_heattr", zztmp * zw3d)
      CALL profile_psy_data6 % PostEnd
    END IF
    CALL profile_psy_data7 % PreStart('ldf_eiv_dia', 'r7', 0, 0)
    IF (ln_diaptr) CALL dia_ptr_hst(jp_tem, 'eiv', 0.5 * zw3d)
    zztmp = 0.5_wp * 0.5
    CALL profile_psy_data7 % PostEnd
    IF (iom_use('ueiv_salttr') .OR. iom_use('ueiv_salttr3d')) THEN
      !$ACC KERNELS
      zw2d(:, :) = 0._wp
      zw3d(:, :, :) = 0._wp
      DO jk = 1, jpkm1
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zw3d(ji, jj, jk) = zw3d(ji, jj, jk) * (psi_uw(ji, jj, jk + 1) - psi_uw(ji, jj, jk)) * (tsn(ji, jj, jk, jp_sal) + tsn(ji + 1, jj, jk, jp_sal))
            zw2d(ji, jj) = zw2d(ji, jj) + zw3d(ji, jj, jk)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data8 % PreStart('ldf_eiv_dia', 'r8', 0, 0)
      CALL lbc_lnk('ldftra', zw2d, 'U', - 1.)
      CALL lbc_lnk('ldftra', zw3d, 'U', - 1.)
      CALL iom_put("ueiv_salttr", zztmp * zw2d)
      CALL iom_put("ueiv_salttr3d", zztmp * zw3d)
      CALL profile_psy_data8 % PostEnd
    END IF
    IF (iom_use('veiv_salttr')) THEN
      !$ACC KERNELS
      zw2d(:, :) = 0._wp
      zw3d(:, :, :) = 0._wp
      DO jk = 1, jpkm1
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zw3d(ji, jj, jk) = zw3d(ji, jj, jk) + (psi_vw(ji, jj, jk + 1) - psi_vw(ji, jj, jk)) * (tsn(ji, jj, jk, jp_sal) + tsn(ji, jj + 1, jk, jp_sal))
            zw2d(ji, jj) = zw2d(ji, jj) + zw3d(ji, jj, jk)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data9 % PreStart('ldf_eiv_dia', 'r9', 0, 0)
      CALL lbc_lnk('ldftra', zw2d, 'V', - 1.)
      CALL iom_put("veiv_salttr", zztmp * zw2d)
      CALL iom_put("veiv_salttr", zztmp * zw3d)
      CALL profile_psy_data9 % PostEnd
    END IF
    CALL profile_psy_data10 % PreStart('ldf_eiv_dia', 'r10', 0, 0)
    IF (ln_diaptr) CALL dia_ptr_hst(jp_sal, 'eiv', 0.5 * zw3d)
    CALL profile_psy_data10 % PostEnd
  END SUBROUTINE ldf_eiv_dia
END MODULE ldftra