MODULE diadct
  LOGICAL, PUBLIC, PARAMETER :: lk_diadct = .FALSE.
  PUBLIC
  CONTAINS
  SUBROUTINE dia_dct_init
    IMPLICIT NONE
    WRITE(*, FMT = *) 'dia_dct_init: You should not have seen this print! error?'
  END SUBROUTINE dia_dct_init
  SUBROUTINE dia_dct(kt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kt
    WRITE(*, FMT = *) 'dia_dct: You should not have seen this print! error?', kt
  END SUBROUTINE dia_dct
END MODULE diadct