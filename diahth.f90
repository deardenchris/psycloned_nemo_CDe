MODULE diahth
  LOGICAL, PUBLIC, PARAMETER :: lk_diahth = .FALSE.
  CONTAINS
  SUBROUTINE dia_hth(kt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kt
    WRITE(*, *) 'dia_hth: You should not have seen this print! error?', kt
  END SUBROUTINE dia_hth
END MODULE diahth