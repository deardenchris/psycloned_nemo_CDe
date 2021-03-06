MODULE par_kind
  IMPLICIT NONE
  PRIVATE
  INTEGER, PUBLIC, PARAMETER :: jpbyt = 8
  INTEGER, PUBLIC, PARAMETER :: jpbytda = 4
  INTEGER, PUBLIC, PARAMETER :: sp = SELECTED_REAL_KIND(6, 37)
  INTEGER, PUBLIC, PARAMETER :: dp = SELECTED_REAL_KIND(12, 307)
  INTEGER, PUBLIC, PARAMETER :: wp = dp
  INTEGER, PUBLIC, PARAMETER :: i4 = SELECTED_INT_KIND(9)
  INTEGER, PUBLIC, PARAMETER :: i8 = SELECTED_INT_KIND(14)
  INTEGER, PUBLIC, PARAMETER :: lc = 256
END MODULE par_kind