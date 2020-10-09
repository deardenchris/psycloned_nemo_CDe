MODULE c1d
  USE par_kind
  LOGICAL, PUBLIC, PARAMETER :: lk_c1d = .FALSE.
  REAL(KIND = wp) :: rn_lat1d, rn_lon1d
  LOGICAL, PUBLIC :: ln_c1d_locpt = .FALSE.
  CONTAINS
  SUBROUTINE c1d_init
  END SUBROUTINE c1d_init
END MODULE c1d