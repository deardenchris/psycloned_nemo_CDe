MODULE obs_mpp
  USE dom_oce, ONLY: nproc, mig, mjg
  USE mpp_map, ONLY: mppmap
  USE in_out_manager
  USE lib_mpp, ONLY: mpi_comm_oce
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: obs_mpp_bcast_integer, obs_mpp_max_integer, obs_mpp_find_obs_proc, obs_mpp_sum_integers, obs_mpp_sum_integer, &
&mpp_alltoall_int, mpp_alltoallv_int, mpp_alltoallv_real, mpp_global_max
  CONTAINS
  SUBROUTINE obs_mpp_bcast_integer(kvals, kno, kroot)
    INTEGER, INTENT(IN) :: kno
    INTEGER, INTENT(IN) :: kroot
    INTEGER, DIMENSION(kno), INTENT(INOUT) :: kvals
    INTEGER :: ierr
    INCLUDE 'mpif.h'
    CALL mpi_bcast(kvals, kno, mpi_integer, kroot, mpi_comm_oce, ierr)
  END SUBROUTINE obs_mpp_bcast_integer
  SUBROUTINE obs_mpp_max_integer(kvals, kno)
    INTEGER, INTENT(IN) :: kno
    INTEGER, DIMENSION(kno), INTENT(INOUT) :: kvals
    INTEGER :: ierr
    INTEGER, DIMENSION(kno) :: ivals
    INCLUDE 'mpif.h'
    CALL mpi_allreduce(kvals, ivals, kno, mpi_integer, mpi_max, mpi_comm_oce, ierr)
    !$ACC KERNELS
    kvals(:) = ivals(:)
    !$ACC END KERNELS
  END SUBROUTINE obs_mpp_max_integer
  SUBROUTINE obs_mpp_find_obs_proc(kobsp, kno)
    INTEGER, INTENT(IN) :: kno
    INTEGER, DIMENSION(kno), INTENT(INOUT) :: kobsp
    INTEGER :: ji, isum
    INTEGER, DIMENSION(kno) :: iobsp
    !$ACC KERNELS
    iobsp(:) = kobsp(:)
    WHERE (iobsp(:) == - 1)
      iobsp(:) = 9999999
    END WHERE
    iobsp(:) = - 1 * iobsp(:)
    !$ACC END KERNELS
    CALL obs_mpp_max_integer(iobsp, kno)
    !$ACC KERNELS
    kobsp(:) = - 1 * iobsp(:)
    isum = 0
    !$ACC END KERNELS
    DO ji = 1, kno
      IF (kobsp(ji) == 9999999) THEN
        isum = isum + 1
        kobsp(ji) = - 1
      END IF
    END DO
    IF (isum > 0) THEN
      IF (lwp) WRITE(numout, FMT = *) isum, ' observations failed the grid search.'
      IF (lwp) WRITE(numout, FMT = *) 'If ln_grid_search_lookup=.TRUE., try reducing grid_search_res'
    END IF
  END SUBROUTINE obs_mpp_find_obs_proc
  SUBROUTINE obs_mpp_sum_integers(kvalsin, kvalsout, kno)
    INTEGER, INTENT(IN) :: kno
    INTEGER, DIMENSION(kno), INTENT(IN) :: kvalsin
    INTEGER, DIMENSION(kno), INTENT(OUT) :: kvalsout
    INTEGER :: ierr
    INCLUDE 'mpif.h'
    CALL mpi_allreduce(kvalsin, kvalsout, kno, mpi_integer, mpi_sum, mpi_comm_oce, ierr)
  END SUBROUTINE obs_mpp_sum_integers
  SUBROUTINE obs_mpp_sum_integer(kvalin, kvalout)
    INTEGER, INTENT(IN) :: kvalin
    INTEGER, INTENT(OUT) :: kvalout
    INTEGER :: ierr
    INCLUDE 'mpif.h'
    CALL mpi_allreduce(kvalin, kvalout, 1, mpi_integer, mpi_sum, mpi_comm_oce, ierr)
  END SUBROUTINE obs_mpp_sum_integer
  SUBROUTINE mpp_global_max(pval)
    REAL(KIND = wp), DIMENSION(jpiglo, jpjglo), INTENT(INOUT) :: pval
    INTEGER :: ierr
    INCLUDE 'mpif.h'
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: zcp
    ALLOCATE(zcp(jpiglo, jpjglo))
    !$ACC KERNELS
    zcp(:, :) = pval(:, :)
    !$ACC END KERNELS
    CALL mpi_allreduce(zcp, pval, jpiglo * jpjglo, mpi_double_precision, mpi_max, mpi_comm_oce, ierr)
    DEALLOCATE(zcp)
  END SUBROUTINE mpp_global_max
  SUBROUTINE mpp_alltoall_int(kno, kvalsin, kvalsout)
    INTEGER, INTENT(IN) :: kno
    INTEGER, DIMENSION(kno * jpnij), INTENT(IN) :: kvalsin
    INTEGER, DIMENSION(kno * jpnij), INTENT(OUT) :: kvalsout
    INTEGER :: ierr
    INCLUDE 'mpif.h'
    CALL mpi_alltoall(kvalsin, kno, mpi_integer, kvalsout, kno, mpi_integer, mpi_comm_oce, ierr)
  END SUBROUTINE mpp_alltoall_int
  SUBROUTINE mpp_alltoallv_int(kvalsin, knoin, kinv, kvalsout, knoout, koutv)
    INTEGER, INTENT(IN) :: knoin
    INTEGER, INTENT(IN) :: knoout
    INTEGER, DIMENSION(jpnij) :: kinv, koutv
    INTEGER, DIMENSION(knoin), INTENT(IN) :: kvalsin
    INTEGER, DIMENSION(knoout), INTENT(OUT) :: kvalsout
    INTEGER :: ierr
    INTEGER :: jproc
    INCLUDE 'mpif.h'
    INTEGER, DIMENSION(jpnij) :: irdsp, isdsp
    !$ACC KERNELS
    irdsp(1) = 0
    isdsp(1) = 0
    DO jproc = 2, jpnij
      isdsp(jproc) = isdsp(jproc - 1) + kinv(jproc - 1)
      irdsp(jproc) = irdsp(jproc - 1) + koutv(jproc - 1)
    END DO
    !$ACC END KERNELS
    CALL mpi_alltoallv(kvalsin, kinv, isdsp, mpi_integer, kvalsout, koutv, irdsp, mpi_integer, mpi_comm_oce, ierr)
  END SUBROUTINE mpp_alltoallv_int
  SUBROUTINE mpp_alltoallv_real(pvalsin, knoin, kinv, pvalsout, knoout, koutv)
    INTEGER, INTENT(IN) :: knoin
    INTEGER, INTENT(IN) :: knoout
    INTEGER, DIMENSION(jpnij) :: kinv, koutv
    REAL(KIND = wp), DIMENSION(knoin), INTENT(IN) :: pvalsin
    REAL(KIND = wp), DIMENSION(knoout), INTENT(OUT) :: pvalsout
    INTEGER :: ierr
    INTEGER :: jproc
    INCLUDE 'mpif.h'
    INTEGER, DIMENSION(jpnij) :: irdsp, isdsp
    !$ACC KERNELS
    irdsp(1) = 0
    isdsp(1) = 0
    DO jproc = 2, jpnij
      isdsp(jproc) = isdsp(jproc - 1) + kinv(jproc - 1)
      irdsp(jproc) = irdsp(jproc - 1) + koutv(jproc - 1)
    END DO
    !$ACC END KERNELS
    CALL mpi_alltoallv(pvalsin, kinv, isdsp, mpi_double_precision, pvalsout, koutv, irdsp, mpi_double_precision, mpi_comm_oce, ierr)
  END SUBROUTINE mpp_alltoallv_real
END MODULE obs_mpp