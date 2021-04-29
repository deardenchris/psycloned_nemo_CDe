MODULE obs_mpp
  USE dom_oce, ONLY: nproc, mig, mjg
  USE mpp_map, ONLY: mppmap
  USE in_out_manager
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: obs_mpp_bcast_integer, obs_mpp_max_integer, obs_mpp_find_obs_proc, obs_mpp_sum_integers, obs_mpp_sum_integer, mpp_alltoall_int, mpp_alltoallv_int, mpp_alltoallv_real, mpp_global_max
  CONTAINS
  SUBROUTINE obs_mpp_bcast_integer(kvals, kno, kroot)
    INTEGER, INTENT(IN) :: kno
    INTEGER, INTENT(IN) :: kroot
    INTEGER, DIMENSION(kno), INTENT(INOUT) :: kvals
  END SUBROUTINE obs_mpp_bcast_integer
  SUBROUTINE obs_mpp_max_integer(kvals, kno)
    INTEGER, INTENT(IN) :: kno
    INTEGER, DIMENSION(kno), INTENT(INOUT) :: kvals
  END SUBROUTINE obs_mpp_max_integer
  SUBROUTINE obs_mpp_find_obs_proc(kobsp, kno)
    INTEGER, INTENT(IN) :: kno
    INTEGER, DIMENSION(kno), INTENT(INOUT) :: kobsp
  END SUBROUTINE obs_mpp_find_obs_proc
  SUBROUTINE obs_mpp_sum_integers(kvalsin, kvalsout, kno)
    INTEGER, INTENT(IN) :: kno
    INTEGER, DIMENSION(kno), INTENT(IN) :: kvalsin
    INTEGER, DIMENSION(kno), INTENT(OUT) :: kvalsout
    !$ACC KERNELS
    kvalsout(:) = kvalsin(:)
    !$ACC END KERNELS
  END SUBROUTINE obs_mpp_sum_integers
  SUBROUTINE obs_mpp_sum_integer(kvalin, kvalout)
    INTEGER, INTENT(IN) :: kvalin
    INTEGER, INTENT(OUT) :: kvalout
    kvalout = kvalin
  END SUBROUTINE obs_mpp_sum_integer
  SUBROUTINE mpp_global_max(pval)
    REAL(KIND = wp), DIMENSION(jpiglo, jpjglo), INTENT(INOUT) :: pval
    INTEGER :: ierr
  END SUBROUTINE mpp_global_max
  SUBROUTINE mpp_alltoall_int(kno, kvalsin, kvalsout)
    INTEGER, INTENT(IN) :: kno
    INTEGER, DIMENSION(kno * jpnij), INTENT(IN) :: kvalsin
    INTEGER, DIMENSION(kno * jpnij), INTENT(OUT) :: kvalsout
    INTEGER :: ierr
    kvalsout = kvalsin
  END SUBROUTINE mpp_alltoall_int
  SUBROUTINE mpp_alltoallv_int(kvalsin, knoin, kinv, kvalsout, knoout, koutv)
    INTEGER, INTENT(IN) :: knoin
    INTEGER, INTENT(IN) :: knoout
    INTEGER, DIMENSION(jpnij) :: kinv, koutv
    INTEGER, DIMENSION(knoin), INTENT(IN) :: kvalsin
    INTEGER, DIMENSION(knoout), INTENT(OUT) :: kvalsout
    INTEGER :: ierr
    INTEGER :: jproc
    kvalsout = kvalsin
  END SUBROUTINE mpp_alltoallv_int
  SUBROUTINE mpp_alltoallv_real(pvalsin, knoin, kinv, pvalsout, knoout, koutv)
    INTEGER, INTENT(IN) :: knoin
    INTEGER, INTENT(IN) :: knoout
    INTEGER, DIMENSION(jpnij) :: kinv, koutv
    REAL(KIND = wp), DIMENSION(knoin), INTENT(IN) :: pvalsin
    REAL(KIND = wp), DIMENSION(knoout), INTENT(OUT) :: pvalsout
    INTEGER :: ierr
    INTEGER :: jproc
    pvalsout = pvalsin
  END SUBROUTINE mpp_alltoallv_real
END MODULE obs_mpp