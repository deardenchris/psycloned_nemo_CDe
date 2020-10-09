MODULE dynadv
  USE dom_oce
  USE dynadv_cen2
  USE dynadv_ubs
  USE dynkeg
  USE dynzad
  USE in_out_manager
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_adv
  PUBLIC :: dyn_adv_init
  LOGICAL, PUBLIC :: ln_dynadv_OFF
  LOGICAL, PUBLIC :: ln_dynadv_vec
  INTEGER, PUBLIC :: nn_dynkeg
  LOGICAL, PUBLIC :: ln_dynadv_cen2
  LOGICAL, PUBLIC :: ln_dynadv_ubs
  INTEGER, PUBLIC :: n_dynadv
  INTEGER, PUBLIC, PARAMETER :: np_LIN_dyn = 0
  INTEGER, PUBLIC, PARAMETER :: np_VEC_c2 = 1
  INTEGER, PUBLIC, PARAMETER :: np_FLX_c2 = 2
  INTEGER, PUBLIC, PARAMETER :: np_FLX_ubs = 3
  CONTAINS
  SUBROUTINE dyn_adv(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('dyn_adv', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('dyn_adv')
    SELECT CASE (n_dynadv)
    CASE (np_VEC_c2)
      CALL dyn_keg(kt, nn_dynkeg)
      CALL dyn_zad(kt)
    CASE (np_FLX_c2)
      CALL dyn_adv_cen2(kt)
    CASE (np_FLX_ubs)
      CALL dyn_adv_ubs(kt)
    END SELECT
    IF (ln_timing) CALL timing_stop('dyn_adv')
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE dyn_adv
  SUBROUTINE dyn_adv_init
    INTEGER :: ioptio, ios
    NAMELIST /namdyn_adv/ ln_dynadv_OFF, ln_dynadv_vec, nn_dynkeg, ln_dynadv_cen2, ln_dynadv_ubs
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dyn_adv_init : choice/control of the momentum advection scheme'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
    END IF
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namdyn_adv, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namdyn_adv in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namdyn_adv, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namdyn_adv in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namdyn_adv)
    IF (lwp) THEN
      WRITE(numout, FMT = *) '   Namelist namdyn_adv : chose a advection formulation & scheme for momentum'
      WRITE(numout, FMT = *) '      linear dynamics : no momentum advection          ln_dynadv_OFF  = ', ln_dynadv_OFF
      WRITE(numout, FMT = *) '      Vector form: 2nd order centered scheme           ln_dynadv_vec  = ', ln_dynadv_vec
      WRITE(numout, FMT = *) '         with Hollingsworth scheme (=1) or not (=0)       nn_dynkeg   = ', nn_dynkeg
      WRITE(numout, FMT = *) '      flux form: 2nd order centred scheme              ln_dynadv_cen2 = ', ln_dynadv_cen2
      WRITE(numout, FMT = *) '                 3rd order UBS scheme                  ln_dynadv_ubs  = ', ln_dynadv_ubs
    END IF
    ioptio = 0
    IF (ln_dynadv_off) THEN
      ioptio = ioptio + 1
      n_dynadv = np_lin_dyn
    END IF
    IF (ln_dynadv_vec) THEN
      ioptio = ioptio + 1
      n_dynadv = np_vec_c2
    END IF
    IF (ln_dynadv_cen2) THEN
      ioptio = ioptio + 1
      n_dynadv = np_flx_c2
    END IF
    IF (ln_dynadv_ubs) THEN
      ioptio = ioptio + 1
      n_dynadv = np_flx_ubs
    END IF
    IF (ioptio /= 1) CALL ctl_stop('choose ONE and only ONE advection scheme')
    IF (nn_dynkeg /= nkeg_C2 .AND. nn_dynkeg /= nkeg_HW) CALL ctl_stop('KEG scheme wrong value of nn_dynkeg')
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      SELECT CASE (n_dynadv)
      CASE (np_lin_dyn)
        WRITE(numout, FMT = *) '   ==>>>   linear dynamics : no momentum advection used'
      CASE (np_vec_c2)
        WRITE(numout, FMT = *) '   ==>>>   vector form : keg + zad + vor is used'
        IF (nn_dynkeg == nkeg_C2) WRITE(numout, FMT = *) '              with Centered standard keg scheme'
        IF (nn_dynkeg == nkeg_HW) WRITE(numout, FMT = *) '              with Hollingsworth keg scheme'
      CASE (np_flx_c2)
        WRITE(numout, FMT = *) '   ==>>>   flux form   : 2nd order scheme is used'
      CASE (np_flx_ubs)
        WRITE(numout, FMT = *) '   ==>>>   flux form   : UBS       scheme is used'
      END SELECT
    END IF
  END SUBROUTINE dyn_adv_init
END MODULE dynadv