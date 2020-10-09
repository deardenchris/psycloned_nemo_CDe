MODULE bdyice
  CONTAINS
  SUBROUTINE bdy_ice(kt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kt
    WRITE(*, FMT = *) 'bdy_ice: You should not have seen this print! error?', kt
  END SUBROUTINE bdy_ice
END MODULE bdyice