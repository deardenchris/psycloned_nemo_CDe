MODULE mpp_map
  USE par_kind, ONLY: wp
  USE par_oce, ONLY: jpi, jpj
  USE dom_oce, ONLY: mig, mjg, nldi, nlei, nldj, nlej, nlci, nlcj, narea
  USE lib_mpp, ONLY: mpi_comm_oce
  USE in_out_manager
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: mppmap_init, mppmap
  INTEGER, DIMENSION(:, :), ALLOCATABLE :: mppmap
  CONTAINS
  SUBROUTINE mppmap_init
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: imppmap
    INTEGER :: ierr
    INCLUDE 'mpif.h'
    IF (.NOT. ALLOCATED(mppmap)) THEN
      ALLOCATE(mppmap(jpiglo, jpjglo))
    END IF
    ALLOCATE(imppmap(jpiglo, jpjglo))
    !$ACC KERNELS
    imppmap(:, :) = 0
    imppmap(mig(1) : mig(nlci), mjg(1) : mjg(nlcj)) = narea
    !$ACC END KERNELS
    CALL mpi_allreduce(imppmap, mppmap, jpiglo * jpjglo, mpi_integer, mpi_max, mpi_comm_oce, ierr)
  END SUBROUTINE mppmap_init
END MODULE mpp_map