MODULE usrdef_nam
  USE dom_oce, ONLY: nimpp, njmpp
  USE par_oce
  USE phycst
  USE in_out_manager
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: usr_def_nam
  LOGICAL, PUBLIC :: ln_bench
  INTEGER, PUBLIC :: nn_GYRE
  CONTAINS
  SUBROUTINE usr_def_nam(ldtxt, ldnam, cd_cfg, kk_cfg, kpi, kpj, kpk, kperio)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), DIMENSION(:), INTENT(OUT) :: ldtxt, ldnam
    CHARACTER(LEN = *), INTENT(OUT) :: cd_cfg
    INTEGER, INTENT(OUT) :: kk_cfg
    INTEGER, INTENT(OUT) :: kpi, kpj, kpk
    INTEGER, INTENT(OUT) :: kperio
    INTEGER :: ios, ii
    NAMELIST /namusr_def/ nn_GYRE, ln_bench, jpkglo
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('usr_def_nam', 'r0', 0, 0)
    ii = 1
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namusr_def, IOSTAT = ios, ERR = 902)
902 IF (ios /= 0) CALL ctl_nam(ios, 'namusr_def in configuration namelist', .TRUE.)
    WRITE(ldnam(:), namusr_def)
    cd_cfg = 'GYRE'
    kk_cfg = nn_GYRE
    kpi = 30 * nn_GYRE + 2
    kpj = 20 * nn_GYRE + 2
    kpk = jpkglo
    WRITE(ldtxt(ii), *) '   '
    ii = ii + 1
    WRITE(ldtxt(ii), *) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'
    ii = ii + 1
    WRITE(ldtxt(ii), *) '~~~~~~~~~~~ '
    ii = ii + 1
    WRITE(ldtxt(ii), *) '   Namelist namusr_def : GYRE case'
    ii = ii + 1
    WRITE(ldtxt(ii), *) '      GYRE used as Benchmark (=T)                      ln_bench  = ', ln_bench
    ii = ii + 1
    WRITE(ldtxt(ii), *) '      inverse resolution & implied domain size         nn_GYRE   = ', nn_gyre
    ii = ii + 1
    WRITE(ldtxt(ii), *) '         jpiglo = 30*nn_GYRE+2                            jpiglo = ', kpi
    ii = ii + 1
    WRITE(ldtxt(ii), *) '         jpjglo = 20*nn_GYRE+2                            jpjglo = ', kpj
    ii = ii + 1
    WRITE(ldtxt(ii), *) '      number of model levels                              jpkglo = ', kpk
    ii = ii + 1
    kperio = 0
    WRITE(ldtxt(ii), *) '   '
    ii = ii + 1
    WRITE(ldtxt(ii), *) '   Lateral b.c. of the global domain set to closed     jperio = ', kperio
    ii = ii + 1
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE usr_def_nam
END MODULE usrdef_nam