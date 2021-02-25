MODULE solfrac_mod
  USE par_kind
  IMPLICIT NONE
  REAL(KIND = wp), PRIVATE, PARAMETER, DIMENSION(9) :: pp_wgt = (/0.2370, 0.36, 0.1790, 0.087, 0.08, 0.025, 0.025, 0.007, 0.0004/)
  REAL(KIND = wp), PRIVATE, PARAMETER, DIMENSION(9) :: pp_len = (/34.84, 2.266, 0.0315, 0.0055, 8.32E-4, 1.26E-4, 3.13E-4, &
&7.82E-4, 1.44E-5/)
  PUBLIC :: solfrac
  CONTAINS
  REAL(KIND = dp) FUNCTION solfrac(ptop, pbottom)
    REAL(KIND = wp), INTENT(IN) :: ptop, pbottom
    INTEGER :: jt
    solfrac = 0._wp
    DO jt = 1, 9
      solfrac = solfrac + pp_wgt(jt) * (EXP(- ptop / pp_len(jt)) - EXP(- pbottom / pp_len(jt)))
    END DO
  END FUNCTION
END MODULE solfrac_mod