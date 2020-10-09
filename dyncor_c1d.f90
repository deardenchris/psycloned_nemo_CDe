MODULE dyncor_c1d
  CONTAINS
  SUBROUTINE cor_c1d
    IMPLICIT NONE
  END SUBROUTINE cor_c1d
  SUBROUTINE dyn_cor_c1d(kt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kt
    WRITE(*, FMT = *) 'dyn_cor_c1d: You should not have seen this print! error?', kt
  END SUBROUTINE dyn_cor_c1d
END MODULE dyncor_c1d