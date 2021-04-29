MODULE trdtrc
  CONTAINS
  SUBROUTINE trd_trc(ptrtrd, kjn, ktrd, kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: kt, kjn, ktrd
    REAL :: ptrtrd(:, :, :)
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('trd_trc', 'r0', 0, 0)
    WRITE(*, *) 'trd_trc : You should not have seen this print! error?', ptrtrd(1, 1, 1)
    WRITE(*, *) '  "      "      : You should not have seen this print! error?', kjn, ktrd, kt
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE trd_trc
END MODULE trdtrc