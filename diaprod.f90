MODULE diaprod
  LOGICAL, PUBLIC, PARAMETER :: lk_diaprod = .FALSE.
  CONTAINS
  SUBROUTINE dia_prod(kt)
    INTEGER :: kt
    IF (lk_diaprod) WRITE(*, FMT = *) 'dia_prod: You should not have seen this print! error?', kt
  END SUBROUTINE dia_prod
END MODULE diaprod