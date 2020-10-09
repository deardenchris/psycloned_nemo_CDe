MODULE icestp
  CONTAINS
  SUBROUTINE ice_stp(kt, ksbc)
    INTEGER, INTENT(IN) :: kt, ksbc
    WRITE(*, FMT = *) 'ice_stp: You should not have seen this print! error?', kt
  END SUBROUTINE ice_stp
  SUBROUTINE ice_init
  END SUBROUTINE ice_init
END MODULE icestp