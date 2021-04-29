MODULE sbcice_cice
  CONTAINS
  SUBROUTINE sbc_ice_cice(kt, ksbc)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kt, ksbc
    WRITE(*, *) 'sbc_ice_cice: You should not have seen this print! error?', kt
  END SUBROUTINE sbc_ice_cice
  SUBROUTINE cice_sbc_init(ksbc)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ksbc
    WRITE(*, *) 'cice_sbc_init: You should not have seen this print! error?', ksbc
  END SUBROUTINE cice_sbc_init
  SUBROUTINE cice_sbc_final
    IMPLICIT NONE
    WRITE(*, *) 'cice_sbc_final: You should not have seen this print! error?'
  END SUBROUTINE cice_sbc_final
END MODULE sbcice_cice