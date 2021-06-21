MODULE icethd_zdf
  USE dom_oce
  USE phycst
  USE ice
  USE icethd_zdf_BL99
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_thd_zdf
  PUBLIC :: ice_thd_zdf_init
  INTEGER :: nice_zdf
  INTEGER, PARAMETER :: np_BL99 = 1
  LOGICAL :: ln_zdf_BL99
  CONTAINS
  SUBROUTINE ice_thd_zdf
    SELECT CASE (nice_zdf)
    CASE (np_BL99)
      IF (.NOT. ln_cndflx) THEN
        CALL ice_thd_zdf_BL99(np_cnd_OFF)
      ELSE IF (ln_cndflx .AND. .NOT. ln_cndemulate) THEN
        CALL ice_thd_zdf_BL99(np_cnd_ON)
      ELSE IF (ln_cndflx .AND. ln_cndemulate) THEN
        CALL ice_thd_zdf_BL99(np_cnd_EMU)
        CALL ice_thd_zdf_BL99(np_cnd_ON)
      END IF
    END SELECT
  END SUBROUTINE ice_thd_zdf
  SUBROUTINE ice_thd_zdf_init
    INTEGER :: ios, ioptio
    NAMELIST /namthd_zdf/ ln_zdf_BL99, ln_cndi_U64, ln_cndi_P07, rn_cnd_s, rn_kappa_i
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, namthd_zdf, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namthd_zdf in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, namthd_zdf, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namthd_zdf in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, namthd_zdf)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ice_thd_zdf_init: Ice vertical heat diffusion'
      WRITE(numout, *) '~~~~~~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist namthd_zdf:'
      WRITE(numout, *) '      Bitz and Lipscomb (1999) formulation                    ln_zdf_BL99  = ', ln_zdf_BL99
      WRITE(numout, *) '      thermal conductivity in the ice (Untersteiner 1964)     ln_cndi_U64  = ', ln_cndi_U64
      WRITE(numout, *) '      thermal conductivity in the ice (Pringle et al 2007)    ln_cndi_P07  = ', ln_cndi_P07
      WRITE(numout, *) '      thermal conductivity in the snow                        rn_cnd_s     = ', rn_cnd_s
      WRITE(numout, *) '      extinction radiation parameter in sea ice               rn_kappa_i   = ', rn_kappa_i
    END IF
    IF ((ln_cndi_U64 .AND. ln_cndi_P07) .OR. (.NOT. ln_cndi_U64 .AND. .NOT. ln_cndi_P07)) THEN
      CALL ctl_stop('ice_thd_zdf_init: choose 1 and only 1 formulation for thermal conduction (ln_cndi_U64 or ln_cndi_P07)')
    END IF
    ioptio = 0
    IF (ln_zdf_bl99) THEN
      ioptio = ioptio + 1
      nice_zdf = np_bl99
    END IF
    IF (ioptio /= 1) CALL ctl_stop('ice_thd_init: one and only one ice thermo option has to be defined ')
  END SUBROUTINE ice_thd_zdf_init
END MODULE icethd_zdf