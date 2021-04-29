MODULE trdini
  USE dom_oce
  USE trd_oce
  USE trdken
  USE trdpen
  USE trdglo
  USE trdmxl
  USE trdvor
  USE in_out_manager
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: trd_init
  CONTAINS
  SUBROUTINE trd_init
    INTEGER :: ios
    NAMELIST /namtrd/ ln_dyn_trd, ln_KE_trd, ln_vor_trd, ln_dyn_mxl, ln_tra_trd, ln_PE_trd, ln_glo_trd, ln_tra_mxl, nn_trd
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namtrd, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namtrd in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namtrd, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namtrd in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namtrd)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'trd_init : Momentum/Tracers trends'
      WRITE(numout, *) '~~~~~~~~'
      WRITE(numout, *) '   Namelist namtrd : set trends parameters'
      WRITE(numout, *) '      global domain averaged dyn & tra trends   ln_glo_trd  = ', ln_glo_trd
      WRITE(numout, *) '      U & V trends: 3D output                   ln_dyn_trd  = ', ln_dyn_trd
      WRITE(numout, *) '      U & V trends: Mixed Layer averaged        ln_dyn_mxl  = ', ln_dyn_mxl
      WRITE(numout, *) '      T & S trends: 3D output                   ln_tra_trd  = ', ln_tra_trd
      WRITE(numout, *) '      T & S trends: Mixed Layer averaged        ln_tra_mxl  = ', ln_tra_mxl
      WRITE(numout, *) '      Kinetic   Energy trends                   ln_KE_trd   = ', ln_KE_trd
      WRITE(numout, *) '      Potential Energy trends                   ln_PE_trd   = ', ln_PE_trd
      WRITE(numout, *) '      Barotropic vorticity trends               ln_vor_trd  = ', ln_vor_trd
      WRITE(numout, *) '      frequency of trends diagnostics (glo)     nn_trd      = ', nn_trd
    END IF
    l_trdtra = .FALSE.
    IF (ln_tra_trd .OR. ln_PE_trd .OR. ln_tra_mxl .OR. ln_glo_trd) l_trdtra = .TRUE.
    l_trddyn = .FALSE.
    IF (ln_dyn_trd .OR. ln_KE_trd .OR. ln_dyn_mxl .OR. ln_vor_trd .OR. ln_glo_trd) l_trddyn = .TRUE.
    IF (ln_dyn_mxl) CALL ctl_stop('ML diag on momentum are not yet coded we stop')
    IF (ln_tra_mxl .OR. ln_vor_trd) CALL ctl_stop('ML tracer and Barotropic vorticity diags are still using old IOIPSL')
    IF (ln_glo_trd) CALL trd_glo_init
    IF (ln_tra_mxl) CALL trd_mxl_init
    IF (ln_vor_trd) CALL trd_vor_init
    IF (ln_KE_trd) CALL trd_ken_init
    IF (ln_PE_trd) CALL trd_pen_init
  END SUBROUTINE trd_init
END MODULE trdini