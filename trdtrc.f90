MODULE trdtrc
  CONTAINS
  SUBROUTINE trd_trc(ptrtrd, kjn, ktrd, kt)
    INTEGER :: kt, kjn, ktrd
    REAL :: ptrtrd(:, :, :)
    WRITE(*, FMT = *) 'trd_trc : You should not have seen this print! error?', ptrtrd(1, 1, 1)
    WRITE(*, FMT = *) '  "      "      : You should not have seen this print! error?', kjn, ktrd, kt
  END SUBROUTINE trd_trc
END MODULE trdtrc