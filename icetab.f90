MODULE icetab
  USE par_oce
  USE ice, ONLY: jpl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tab_3d_2d
  PUBLIC :: tab_2d_1d
  PUBLIC :: tab_2d_3d
  PUBLIC :: tab_1d_2d
  CONTAINS
  SUBROUTINE tab_3d_2d(ndim1d, tab_ind, tab1d, tab2d)
    INTEGER, INTENT(IN) :: ndim1d
    INTEGER, DIMENSION(ndim1d), INTENT(IN) :: tab_ind
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl), INTENT(IN) :: tab2d
    REAL(KIND = wp), DIMENSION(ndim1d, jpl), INTENT(OUT) :: tab1d
    INTEGER :: jl, jn, jid, jjd
    !$ACC KERNELS
    !$ACC loop independent
    DO jl = 1, jpl
      DO jn = 1, ndim1d
        jid = MOD(tab_ind(jn) - 1, jpi) + 1
        jjd = (tab_ind(jn) - 1) / jpi + 1
        tab1d(jn, jl) = tab2d(jid, jjd, jl)
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE tab_3d_2d
  SUBROUTINE tab_2d_1d(ndim1d, tab_ind, tab1d, tab2d)
    INTEGER, INTENT(IN) :: ndim1d
    INTEGER, DIMENSION(ndim1d), INTENT(IN) :: tab_ind
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: tab2d
    REAL(KIND = wp), DIMENSION(ndim1d), INTENT(OUT) :: tab1d
    INTEGER :: jn, jid, jjd
    !$ACC KERNELS
    !$ACC loop independent
    DO jn = 1, ndim1d
      jid = MOD(tab_ind(jn) - 1, jpi) + 1
      jjd = (tab_ind(jn) - 1) / jpi + 1
      tab1d(jn) = tab2d(jid, jjd)
    END DO
    !$ACC END KERNELS
  END SUBROUTINE tab_2d_1d
  SUBROUTINE tab_2d_3d(ndim1d, tab_ind, tab1d, tab2d)
    INTEGER, INTENT(IN) :: ndim1d
    INTEGER, DIMENSION(ndim1d), INTENT(IN) :: tab_ind
    REAL(KIND = wp), DIMENSION(ndim1d, jpl), INTENT(IN) :: tab1d
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl), INTENT(OUT) :: tab2d
    INTEGER :: jl, jn, jid, jjd
    !$ACC KERNELS
    !$ACC loop independent
    DO jl = 1, jpl
      DO jn = 1, ndim1d
        jid = MOD(tab_ind(jn) - 1, jpi) + 1
        jjd = (tab_ind(jn) - 1) / jpi + 1
        tab2d(jid, jjd, jl) = tab1d(jn, jl)
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE tab_2d_3d
  SUBROUTINE tab_1d_2d(ndim1d, tab_ind, tab1d, tab2d)
    INTEGER, INTENT(IN) :: ndim1d
    INTEGER, DIMENSION(ndim1d), INTENT(IN) :: tab_ind
    REAL(KIND = wp), DIMENSION(ndim1d), INTENT(IN) :: tab1d
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(OUT) :: tab2d
    INTEGER :: jn, jid, jjd
    !$ACC KERNELS
    !$ACC loop independent
    DO jn = 1, ndim1d
      jid = MOD(tab_ind(jn) - 1, jpi) + 1
      jjd = (tab_ind(jn) - 1) / jpi + 1
      tab2d(jid, jjd) = tab1d(jn)
    END DO
    !$ACC END KERNELS
  END SUBROUTINE tab_1d_2d
END MODULE icetab