MODULE diaharm
  LOGICAL, PUBLIC, PARAMETER :: lk_diaharm = .FALSE.
  CONTAINS
  SUBROUTINE dia_harm(kt)
    INTEGER, INTENT(IN) :: kt
    WRITE(*, FMT = *) 'dia_harm: you should not have seen this print'
  END SUBROUTINE dia_harm
END MODULE diaharm