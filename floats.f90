MODULE floats
  CONTAINS
  SUBROUTINE flo_stp(kt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kt
    WRITE(*, FMT = *) 'flo_stp: You should not have seen this print! error?', kt
  END SUBROUTINE flo_stp
  SUBROUTINE flo_init
    IMPLICIT NONE
  END SUBROUTINE flo_init
END MODULE floats