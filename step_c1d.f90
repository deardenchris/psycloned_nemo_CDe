MODULE step_c1d
  CONTAINS
  SUBROUTINE stp_c1d(kt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kt
    WRITE(*, FMT = *) 'stp_c1d: You should not have seen this print! error?', kt
  END SUBROUTINE stp_c1d
END MODULE step_c1d