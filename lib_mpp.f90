MODULE lib_mpp
  USE dom_oce
  USE lbcnfd
  USE in_out_manager
  IMPLICIT NONE
  PRIVATE
  INTERFACE mpp_nfd
    MODULE PROCEDURE mpp_nfd_2d, mpp_nfd_3d, mpp_nfd_4d
    MODULE PROCEDURE mpp_nfd_2d_ptr, mpp_nfd_3d_ptr, mpp_nfd_4d_ptr
  END INTERFACE
  PUBLIC :: mpp_lnk_2d, mpp_lnk_3d, mpp_lnk_4d
  PUBLIC :: mpp_lnk_2d_ptr, mpp_lnk_3d_ptr, mpp_lnk_4d_ptr
  PUBLIC :: mpp_nfd_2d, mpp_nfd_3d, mpp_nfd_4d
  PUBLIC :: mpp_nfd_2d_ptr, mpp_nfd_3d_ptr, mpp_nfd_4d_ptr
  PUBLIC :: ctl_stop, ctl_warn, get_unit, ctl_opn, ctl_nam
  PUBLIC :: mynode, mppstop, mppsync, mpp_comm_free
  PUBLIC :: mpp_ini_north
  PUBLIC :: mpp_lnk_2d_icb
  PUBLIC :: mpp_lbc_north_icb
  PUBLIC :: mpp_min, mpp_max, mpp_sum, mpp_minloc, mpp_maxloc
  PUBLIC :: mpp_delay_max, mpp_delay_sum, mpp_delay_rcv
  PUBLIC :: mppscatter, mppgather
  PUBLIC :: mpp_ini_znl
  PUBLIC :: mppsend, mpprecv
  PUBLIC :: mpp_lnk_bdy_2d, mpp_lnk_bdy_3d, mpp_lnk_bdy_4d
  INTERFACE mpp_min
    MODULE PROCEDURE mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real
  END INTERFACE
  INTERFACE mpp_max
    MODULE PROCEDURE mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real
  END INTERFACE
  INTERFACE mpp_sum
    MODULE PROCEDURE mppsum_a_int, mppsum_int, mppsum_a_real, mppsum_real, mppsum_realdd, mppsum_a_realdd
  END INTERFACE
  INTERFACE mpp_minloc
    MODULE PROCEDURE mpp_minloc2d, mpp_minloc3d
  END INTERFACE
  INTERFACE mpp_maxloc
    MODULE PROCEDURE mpp_maxloc2d, mpp_maxloc3d
  END INTERFACE
  INCLUDE 'mpif.h'
  LOGICAL, PUBLIC, PARAMETER :: lk_mpp = .TRUE.
  INTEGER, PARAMETER :: nprocmax = 2 ** 10
  INTEGER, PUBLIC :: mppsize
  INTEGER, PUBLIC :: mpprank
  INTEGER, PUBLIC :: mpi_comm_oce
  INTEGER :: MPI_SUMDD
  INTEGER, PUBLIC :: ncomm_znl
  LOGICAL, PUBLIC :: l_znl_root
  INTEGER :: ngrp_znl
  INTEGER :: ndim_rank_znl
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: nrank_znl
  INTEGER, PUBLIC :: ngrp_world
  INTEGER, PUBLIC :: ngrp_opa
  INTEGER, PUBLIC :: ngrp_north
  INTEGER, PUBLIC :: ncomm_north
  INTEGER, PUBLIC :: ndim_rank_north
  INTEGER, PUBLIC :: njmppmax
  INTEGER, PUBLIC :: north_root
  INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE, SAVE :: nrank_north
  CHARACTER(LEN = 1), PUBLIC :: cn_mpi_send
  LOGICAL, PUBLIC :: l_isend = .FALSE.
  INTEGER, PUBLIC :: nn_buffer
  CHARACTER(LEN = 128), DIMENSION(:), ALLOCATABLE :: crname_lbc
  CHARACTER(LEN = 128), DIMENSION(:), ALLOCATABLE :: crname_glb
  CHARACTER(LEN = 128), DIMENSION(:), ALLOCATABLE :: crname_dlg
  INTEGER, PUBLIC :: ncom_stp = 0
  INTEGER, PUBLIC :: ncom_fsbc = 1
  INTEGER, PUBLIC :: ncom_dttrc = 1
  INTEGER, PUBLIC :: ncom_freq
  INTEGER, PUBLIC, DIMENSION(:, :), ALLOCATABLE :: ncomm_sequence
  INTEGER, PARAMETER, PUBLIC :: ncom_rec_max = 5000
  INTEGER, PUBLIC :: n_sequence_lbc = 0
  INTEGER, PUBLIC :: n_sequence_glb = 0
  INTEGER, PUBLIC :: n_sequence_dlg = 0
  INTEGER, PUBLIC :: numcom = - 1
  LOGICAL, PUBLIC :: l_full_nf_update = .TRUE.
  INTEGER, PARAMETER, PUBLIC :: nbdelay = 2
  CHARACTER(LEN = 32), DIMENSION(nbdelay), PUBLIC :: c_delaylist = (/'cflice', 'fwb   '/)
  CHARACTER(LEN = 3), DIMENSION(nbdelay), PUBLIC :: c_delaycpnt = (/'ICE', 'OCE'/)
  TYPE, PUBLIC :: DELAYARR
    REAL(KIND = wp), POINTER, DIMENSION(:) :: z1d => NULL()
    COMPLEX(KIND = wp), POINTER, DIMENSION(:) :: y1d => NULL()
  END TYPE DELAYARR
  TYPE(DELAYARR), DIMENSION(nbdelay), PUBLIC, SAVE :: todelay
  INTEGER, DIMENSION(nbdelay), PUBLIC :: ndelayid = - 1
  REAL(KIND = wp), DIMENSION(2), PUBLIC :: waiting_time = 0._wp
  REAL(KIND = wp), PUBLIC :: compute_time = 0._wp, elapsed_time = 0._wp
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE, SAVE :: tampon
  LOGICAL, PUBLIC :: ln_nnogather
  LOGICAL, PUBLIC :: l_north_nogather = .FALSE.
  CONTAINS
  FUNCTION mynode(ldtxt, ldname, kumnam_ref, kumnam_cfg, kumond, kstop, localComm)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), DIMENSION(:), INTENT(OUT) :: ldtxt
    CHARACTER(LEN = *), INTENT(IN) :: ldname
    INTEGER, INTENT(IN) :: kumnam_ref
    INTEGER, INTENT(IN) :: kumnam_cfg
    INTEGER, INTENT(INOUT) :: kumond
    INTEGER, INTENT(INOUT) :: kstop
    INTEGER, OPTIONAL, INTENT(IN) :: localComm
    INTEGER :: mynode, ierr, code, ji, ii, ios
    LOGICAL :: mpi_was_called
    NAMELIST /nammpp/ cn_mpi_send, nn_buffer, jpni, jpnj, ln_nnogather
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mynode', 'r0', 0, 0)
    ii = 1
    WRITE(ldtxt(ii), FMT = *)
    ii = ii + 1
    WRITE(ldtxt(ii), FMT = *) 'mynode : mpi initialisation'
    ii = ii + 1
    WRITE(ldtxt(ii), FMT = *) '~~~~~~ '
    ii = ii + 1
    REWIND(UNIT = kumnam_ref)
    READ(kumnam_ref, nammpp, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'nammpp in reference namelist', lwp)
    REWIND(UNIT = kumnam_cfg)
    READ(kumnam_cfg, nammpp, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'nammpp in configuration namelist', lwp)
    WRITE(ldtxt(ii), FMT = *) '   Namelist nammpp'
    ii = ii + 1
    WRITE(ldtxt(ii), FMT = *) '      mpi send type          cn_mpi_send = ', cn_mpi_send
    ii = ii + 1
    WRITE(ldtxt(ii), FMT = *) '      size exported buffer   nn_buffer   = ', nn_buffer, ' bytes'
    ii = ii + 1
    IF (jpni < 1 .OR. jpnj < 1) THEN
      WRITE(ldtxt(ii), FMT = *) '      jpni and jpnj will be calculated automatically'
      ii = ii + 1
    ELSE
      WRITE(ldtxt(ii), FMT = *) '      processor grid extent in i         jpni = ', jpni
      ii = ii + 1
      WRITE(ldtxt(ii), FMT = *) '      processor grid extent in j         jpnj = ', jpnj
      ii = ii + 1
    END IF
    WRITE(ldtxt(ii), FMT = *) '      avoid use of mpi_allgather at the north fold  ln_nnogather = ', ln_nnogather
    ii = ii + 1
    CALL mpi_initialized(mpi_was_called, code)
    IF (code /= MPI_SUCCESS) THEN
      DO ji = 1, SIZE(ldtxt)
        IF (TRIM(ldtxt(ji)) /= '') WRITE(*, FMT = *) ldtxt(ji)
      END DO
      WRITE(*, cform_err)
      WRITE(*, FMT = *) 'lib_mpp: Error in routine mpi_initialized'
      CALL mpi_abort(mpi_comm_world, code, ierr)
    END IF
    IF (mpi_was_called) THEN
      SELECT CASE (cn_mpi_send)
      CASE ('S')
        WRITE(ldtxt(ii), FMT = *) '           Standard blocking mpi send (send)'
        ii = ii + 1
      CASE ('B')
        WRITE(ldtxt(ii), FMT = *) '           Buffer blocking mpi send (bsend)'
        ii = ii + 1
        IF (Agrif_Root()) CALL mpi_init_oce(ldtxt, ii, ierr)
      CASE ('I')
        WRITE(ldtxt(ii), FMT = *) '           Immediate non-blocking send (isend)'
        ii = ii + 1
        l_isend = .TRUE.
      CASE DEFAULT
        WRITE(ldtxt(ii), cform_err)
        ii = ii + 1
        WRITE(ldtxt(ii), FMT = *) '           bad value for cn_mpi_send = ', cn_mpi_send
        ii = ii + 1
        kstop = kstop + 1
      END SELECT
    ELSE IF (PRESENT(localComm) .AND. .NOT. mpi_was_called) THEN
      WRITE(ldtxt(ii), cform_err)
      ii = ii + 1
      WRITE(ldtxt(ii), FMT = *) ' lib_mpp: You cannot provide a local communicator '
      ii = ii + 1
      WRITE(ldtxt(ii), FMT = *) '          without calling MPI_Init before ! '
      ii = ii + 1
      kstop = kstop + 1
    ELSE
      SELECT CASE (cn_mpi_send)
      CASE ('S')
        WRITE(ldtxt(ii), FMT = *) '           Standard blocking mpi send (send)'
        ii = ii + 1
        CALL mpi_init(ierr)
      CASE ('B')
        WRITE(ldtxt(ii), FMT = *) '           Buffer blocking mpi send (bsend)'
        ii = ii + 1
        IF (Agrif_Root()) CALL mpi_init_oce(ldtxt, ii, ierr)
      CASE ('I')
        WRITE(ldtxt(ii), FMT = *) '           Immediate non-blocking send (isend)'
        ii = ii + 1
        l_isend = .TRUE.
        CALL mpi_init(ierr)
      CASE DEFAULT
        WRITE(ldtxt(ii), cform_err)
        ii = ii + 1
        WRITE(ldtxt(ii), FMT = *) '           bad value for cn_mpi_send = ', cn_mpi_send
        ii = ii + 1
        kstop = kstop + 1
      END SELECT
    END IF
    IF (PRESENT(localComm)) THEN
      IF (Agrif_Root()) THEN
        mpi_comm_oce = localComm
      END IF
    ELSE
      CALL mpi_comm_dup(mpi_comm_world, mpi_comm_oce, code)
      IF (code /= MPI_SUCCESS) THEN
        DO ji = 1, SIZE(ldtxt)
          IF (TRIM(ldtxt(ji)) /= '') WRITE(*, FMT = *) ldtxt(ji)
        END DO
        WRITE(*, cform_err)
        WRITE(*, FMT = *) ' lib_mpp: Error in routine mpi_comm_dup'
        CALL mpi_abort(mpi_comm_world, code, ierr)
      END IF
    END IF
    CALL mpi_comm_rank(mpi_comm_oce, mpprank, ierr)
    CALL mpi_comm_size(mpi_comm_oce, mppsize, ierr)
    mynode = mpprank
    IF (mynode == 0) THEN
      CALL ctl_opn(kumond, TRIM(ldname), 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', - 1, 6, .FALSE., 1)
      WRITE(kumond, nammpp)
    END IF
    CALL MPI_OP_CREATE(DDPDD_MPI, .TRUE., MPI_SUMDD, ierr)
    CALL profile_psy_data0 % PostEnd
  END FUNCTION mynode
  SUBROUTINE mpp_lnk_2d(cdname, ptab, cd_nat, psgn, cd_mpp, pval)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    INTEGER :: ierr
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = 1
    ipl = 1
    ipf = 1
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ipk, ipl, ipf, ld_lbc = .TRUE.)
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    IF (.NOT. PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(1, :) = ptab(jpim1, :)
          ptab(jpi, :) = ptab(2, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :) = zland
          ptab(nlci - nn_hls + 1 : jpi, :) = zland
          !$ACC END KERNELS
        END IF
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(:, 1) = ptab(:, jpjm1)
          ptab(:, jpj) = ptab(:, 2)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls) = zland
          ptab(:, nlcj - nn_hls + 1 : jpj) = zland
          !$ACC END KERNELS
        END IF
      END DO
    END IF
    IF (ABS(nbondi) == 1) ALLOCATE(zt3ew(jpj, nn_hls, ipk, ipl, ipf, 1), zt3we(jpj, nn_hls, ipk, ipl, ipf, 1))
    IF (nbondi == 0) ALLOCATE(zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    SELECT CASE (nbondi)
    CASE (- 1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3we(:, jh, jk, jl, jf, 1) = ptab(iihom + jh, :)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(nn_hls + jh, :)
              zt3we(:, jh, jk, jl, jf, 1) = ptab(iihom + jh, :)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(nn_hls + jh, :)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpj * ipk * ipl * ipf
    IF (ln_timing) CALL tic_tac(.TRUE.)
    SELECT CASE (nbondi)
    CASE (- 1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, noea)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 1), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    IF (ln_timing) CALL tic_tac(.FALSE.)
    iihom = nlci - nn_hls
    SELECT CASE (nbondi)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(iihom + jh, :) = zt3ew(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :) = zt3we(:, jh, jk, jl, jf, 2)
              ptab(iihom + jh, :) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :) = zt3we(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    IF (nbondi /= 2) DEALLOCATE(zt3ew, zt3we)
    IF (npolj /= 0 .AND. .NOT. PRESENT(cd_mpp)) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat, psgn)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat, psgn)
      END SELECT
    END IF
    IF (ABS(nbondj) == 1) ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 1), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 1))
    IF (nbondj == 0) ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2))
    SELECT CASE (nbondj)
    CASE (- 1)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(:, ijhom + jh)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(:, ijhom + jh)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(:, nn_hls + jh)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(:, nn_hls + jh)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpi * ipk * ipl * ipf
    IF (ln_timing) CALL tic_tac(.TRUE.)
    SELECT CASE (nbondj)
    CASE (- 1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, nono)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    IF (ln_timing) CALL tic_tac(.FALSE.)
    ijhom = nlcj - nn_hls
    SELECT CASE (nbondj)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, ijhom + jh) = zt3ns(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(:, ijhom + jh) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh) = zt3sn(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    IF (nbondj /= 2) DEALLOCATE(zt3ns, zt3sn)
  END SUBROUTINE mpp_lnk_2d
  SUBROUTINE mpp_lnk_2d_ptr(cdname, ptab, cd_nat, psgn, kfld, cd_mpp, pval)
    INTEGER, INTENT(IN) :: kfld
    TYPE(PTR_2D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    INTEGER :: ierr
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = 1
    ipl = 1
    ipf = kfld
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ipk, ipl, ipf, ld_lbc = .TRUE.)
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    IF (.NOT. PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          !$ACC KERNELS     
          ptab(jf) % pt2d(1, :) = ptab(jf) % pt2d(jpim1, :)
          ptab(jf) % pt2d(jpi, :) = ptab(jf) % pt2d(2, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS      
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt2d(1 : nn_hls, :) = zland
          ptab(jf) % pt2d(nlci - nn_hls + 1 : jpi, :) = zland
          !$ACC END KERNELS
        END IF
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(jf) % pt2d(:, 1) = ptab(jf) % pt2d(:, jpjm1)
          ptab(jf) % pt2d(:, jpj) = ptab(jf) % pt2d(:, 2)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS      
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt2d(:, 1 : nn_hls) = zland
          !$ACC END KERNELS
          !$ACC KERNELS
          ptab(jf) % pt2d(:, nlcj - nn_hls + 1 : jpj) = zland
          !$ACC END KERNELS
        END IF
      END DO
    END IF
    IF (ABS(nbondi) == 1) ALLOCATE(zt3ew(jpj, nn_hls, ipk, ipl, ipf, 1), zt3we(jpj, nn_hls, ipk, ipl, ipf, 1))
    IF (nbondi == 0) ALLOCATE(zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    SELECT CASE (nbondi)
    CASE (- 1)
      !$ACC KERNELS
      iihom = nlci - nreci
      !$ACC LOOP GANG VECTOR COLLAPSE(5)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              DO jj = 1, jpj
                zt3we(jj, jh, jk, jl, jf, 1) = ptab(jf) % pt2d(iihom + jh, jj)
              END DO
            END DO
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    CASE (0)
      !$ACC KERNELS
      iihom = nlci - nreci
      !$ACC LOOP SEQ
      DO jf = 1, ipf ! 2-8
        !$ACC LOOP SEQ
        DO jl = 1, ipl ! 1
          !$ACC LOOP SEQ
          DO jk = 1, ipk ! 1
            !$ACC LOOP SEQ
            DO jh = 1, nn_hls ! 1
              !$ACC LOOP GANG VECTOR
              DO jj = 1, jpj
                zt3ew(jj, jh, jk, jl, jf, 1) = ptab(jf) % pt2d(nn_hls + jh, jj)
                zt3we(jj, jh, jk, jl, jf, 1) = ptab(jf) % pt2d(iihom + jh, jj)
              END DO
            END DO
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    CASE (1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(jf) % pt2d(nn_hls + jh, :)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpj * ipk * ipl * ipf
    IF (ln_timing) CALL tic_tac(.TRUE.)
    SELECT CASE (nbondi)
    CASE (- 1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, noea)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 1), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    IF (ln_timing) CALL tic_tac(.FALSE.)
    iihom = nlci - nn_hls
    SELECT CASE (nbondi)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt2d(iihom + jh, :) = zt3ew(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      !$ACC KERNELS
      !$ACC LOOP SEQ
      DO jf = 1, ipf
        !$ACC LOOP SEQ
        DO jl = 1, ipl
          !$ACC LOOP SEQ
          DO jk = 1, ipk
            !$ACC LOOP SEQ
            DO jh = 1, nn_hls
              !$ACC LOOP GANG VECTOR INDEPENDENT
              DO jj = 1, jpj
                ptab(jf) % pt2d(jh, jj) = zt3we(jj, jh, jk, jl, jf, 2)
                ptab(jf) % pt2d(iihom + jh, jj) = zt3ew(jj, jh, jk, jl, jf, 2)
              END DO
            END DO
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt2d(jh, :) = zt3we(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    IF (nbondi /= 2) DEALLOCATE(zt3ew, zt3we)
    IF (npolj /= 0 .AND. .NOT. PRESENT(cd_mpp)) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat(:), psgn(:), ipf)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat(:), psgn(:), ipf)
      END SELECT
    END IF
    IF (ABS(nbondj) == 1) ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 1), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 1))
    IF (nbondj == 0) ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2))
    SELECT CASE (nbondj)
    CASE (- 1)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(jf) % pt2d(:, ijhom + jh)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(jf) % pt2d(:, ijhom + jh)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(jf) % pt2d(:, nn_hls + jh)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(jf) % pt2d(:, nn_hls + jh)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpi * ipk * ipl * ipf
    IF (ln_timing) CALL tic_tac(.TRUE.)
    SELECT CASE (nbondj)
    CASE (- 1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, nono)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    IF (ln_timing) CALL tic_tac(.FALSE.)
    ijhom = nlcj - nn_hls
    SELECT CASE (nbondj)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt2d(:, ijhom + jh) = zt3ns(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt2d(:, jh) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(jf) % pt2d(:, ijhom + jh) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt2d(:, jh) = zt3sn(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    IF (nbondj /= 2) DEALLOCATE(zt3ns, zt3sn)
  END SUBROUTINE mpp_lnk_2d_ptr
  SUBROUTINE mpp_lnk_3d(cdname, ptab, cd_nat, psgn, cd_mpp, pval)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    INTEGER :: ierr
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = SIZE(ptab, 3)
    ipl = 1
    ipf = 1
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ipk, ipl, ipf, ld_lbc = .TRUE.)
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    IF (.NOT. PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(1, :, :) = ptab(jpim1, :, :)
          ptab(jpi, :, :) = ptab(2, :, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :, :) = zland
          ptab(nlci - nn_hls + 1 : jpi, :, :) = zland
          !$ACC END KERNELS
        END IF
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(:, 1, :) = ptab(:, jpjm1, :)
          ptab(:, jpj, :) = ptab(:, 2, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls, :) = zland
          ptab(:, nlcj - nn_hls + 1 : jpj, :) = zland
          !$ACC END KERNELS
        END IF
      END DO
    END IF
    IF (ABS(nbondi) == 1) ALLOCATE(zt3ew(jpj, nn_hls, ipk, ipl, ipf, 1), zt3we(jpj, nn_hls, ipk, ipl, ipf, 1))
    IF (nbondi == 0) ALLOCATE(zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    SELECT CASE (nbondi)
    CASE (- 1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3we(:, jh, jk, jl, jf, 1) = ptab(iihom + jh, :, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(nn_hls + jh, :, jk)
              zt3we(:, jh, jk, jl, jf, 1) = ptab(iihom + jh, :, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(nn_hls + jh, :, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpj * ipk * ipl * ipf
    IF (ln_timing) CALL tic_tac(.TRUE.)
    SELECT CASE (nbondi)
    CASE (- 1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, noea)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 1), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    IF (ln_timing) CALL tic_tac(.FALSE.)
    iihom = nlci - nn_hls
    SELECT CASE (nbondi)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(iihom + jh, :, jk) = zt3ew(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :, jk) = zt3we(:, jh, jk, jl, jf, 2)
              ptab(iihom + jh, :, jk) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :, jk) = zt3we(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    IF (nbondi /= 2) DEALLOCATE(zt3ew, zt3we)
    IF (npolj /= 0 .AND. .NOT. PRESENT(cd_mpp)) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat, psgn)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat, psgn)
      END SELECT
    END IF
    IF (ABS(nbondj) == 1) ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 1), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 1))
    IF (nbondj == 0) ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2))
    SELECT CASE (nbondj)
    CASE (- 1)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(:, ijhom + jh, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(:, ijhom + jh, jk)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(:, nn_hls + jh, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(:, nn_hls + jh, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpi * ipk * ipl * ipf
    IF (ln_timing) CALL tic_tac(.TRUE.)
    SELECT CASE (nbondj)
    CASE (- 1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, nono)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    IF (ln_timing) CALL tic_tac(.FALSE.)
    ijhom = nlcj - nn_hls
    SELECT CASE (nbondj)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, ijhom + jh, jk) = zt3ns(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh, jk) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(:, ijhom + jh, jk) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh, jk) = zt3sn(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    IF (nbondj /= 2) DEALLOCATE(zt3ns, zt3sn)
  END SUBROUTINE mpp_lnk_3d
  SUBROUTINE mpp_lnk_3d_ptr(cdname, ptab, cd_nat, psgn, kfld, cd_mpp, pval)
    INTEGER, INTENT(IN) :: kfld
    TYPE(PTR_3D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    INTEGER :: ierr
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = SIZE(ptab(1) % pt3d, 3)
    ipl = 1
    ipf = kfld
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ipk, ipl, ipf, ld_lbc = .TRUE.)
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    IF (.NOT. PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(jf) % pt3d(1, :, :) = ptab(jf) % pt3d(jpim1, :, :)
          ptab(jf) % pt3d(jpi, :, :) = ptab(jf) % pt3d(2, :, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS      
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt3d(1 : nn_hls, :, :) = zland
          ptab(jf) % pt3d(nlci - nn_hls + 1 : jpi, :, :) = zland
          !$ACC END KERNELS
        END IF
        IF (l_Jperio) THEN
          !$ACC KERNELS      
          ptab(jf) % pt3d(:, 1, :) = ptab(jf) % pt3d(:, jpjm1, :)
          ptab(jf) % pt3d(:, jpj, :) = ptab(jf) % pt3d(:, 2, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS      
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt3d(:, 1 : nn_hls, :) = zland
          ptab(jf) % pt3d(:, nlcj - nn_hls + 1 : jpj, :) = zland
          !$ACC END KERNELS
        END IF
      END DO
    END IF
    IF (ABS(nbondi) == 1) ALLOCATE(zt3ew(jpj, nn_hls, ipk, ipl, ipf, 1), zt3we(jpj, nn_hls, ipk, ipl, ipf, 1))
    IF (nbondi == 0) ALLOCATE(zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    SELECT CASE (nbondi)
    CASE (- 1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3we(:, jh, jk, jl, jf, 1) = ptab(jf) % pt3d(iihom + jh, :, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              !$ACC LOOP VECTOR
              DO jj = 1, jpj
                zt3ew(jj, jh, jk, jl, jf, 1) = ptab(jf) % pt3d(nn_hls + jh, jj, jk) ! Invoked ? Gang vector? Or collapse all loops with explicit syntax
                zt3we(jj, jh, jk, jl, jf, 1) = ptab(jf) % pt3d(iihom + jh, jj, jk)  ! Invoked ! But bad read stride
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(jf) % pt3d(nn_hls + jh, :, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpj * ipk * ipl * ipf
    IF (ln_timing) CALL tic_tac(.TRUE.)
    SELECT CASE (nbondi)
    CASE (- 1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, noea)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 1), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    IF (ln_timing) CALL tic_tac(.FALSE.)
    iihom = nlci - nn_hls
    SELECT CASE (nbondi)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt3d(iihom + jh, :, jk) = zt3ew(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          !$ACC LOOP GANG INDEPENDENT
          DO jk = 1, ipk
          !$ACC LOOP SEQ
            DO jh = 1, nn_hls
              !$ACC LOOP VECTOR INDEPENDENT
              DO jj = 1, jpj
                ptab(jf) % pt3d(jh, jj, jk) = zt3we(jj, jh, jk, jl, jf, 2)
              END DO
              !$ACC LOOP VECTOR INDEPENDENT
              DO jj = 1, jpj
                ptab(jf) % pt3d(iihom + jh, jj, jk) = zt3ew(jj, jh, jk, jl, jf, 2)
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt3d(jh, :, jk) = zt3we(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    IF (nbondi /= 2) DEALLOCATE(zt3ew, zt3we)
    IF (npolj /= 0 .AND. .NOT. PRESENT(cd_mpp)) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat(:), psgn(:), ipf)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat(:), psgn(:), ipf)
      END SELECT
    END IF
    IF (ABS(nbondj) == 1) ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 1), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 1))
    IF (nbondj == 0) ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2))
    SELECT CASE (nbondj)
    CASE (- 1)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(jf) % pt3d(:, ijhom + jh, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(jf) % pt3d(:, ijhom + jh, jk)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(jf) % pt3d(:, nn_hls + jh, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(jf) % pt3d(:, nn_hls + jh, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpi * ipk * ipl * ipf
    IF (ln_timing) CALL tic_tac(.TRUE.)
    SELECT CASE (nbondj)
    CASE (- 1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, nono)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    IF (ln_timing) CALL tic_tac(.FALSE.)
    ijhom = nlcj - nn_hls
    SELECT CASE (nbondj)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt3d(:, ijhom + jh, jk) = zt3ns(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt3d(:, jh, jk) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(jf) % pt3d(:, ijhom + jh, jk) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt3d(:, jh, jk) = zt3sn(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    IF (nbondj /= 2) DEALLOCATE(zt3ns, zt3sn)
  END SUBROUTINE mpp_lnk_3d_ptr
  SUBROUTINE mpp_lnk_4d(cdname, ptab, cd_nat, psgn, cd_mpp, pval)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :, :)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    INTEGER :: ierr
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = SIZE(ptab, 3)
    ipl = SIZE(ptab, 4)
    ipf = 1
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ipk, ipl, ipf, ld_lbc = .TRUE.)
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    IF (.NOT. PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(1, :, :, :) = ptab(jpim1, :, :, :)
          ptab(jpi, :, :, :) = ptab(2, :, :, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :, :, :) = zland
          ptab(nlci - nn_hls + 1 : jpi, :, :, :) = zland
          !$ACC END KERNELS
        END IF
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(:, 1, :, :) = ptab(:, jpjm1, :, :)
          ptab(:, jpj, :, :) = ptab(:, 2, :, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls, :, :) = zland
          ptab(:, nlcj - nn_hls + 1 : jpj, :, :) = zland
          !$ACC END KERNELS
        END IF
      END DO
    END IF
    IF (ABS(nbondi) == 1) ALLOCATE(zt3ew(jpj, nn_hls, ipk, ipl, ipf, 1), zt3we(jpj, nn_hls, ipk, ipl, ipf, 1))
    IF (nbondi == 0) ALLOCATE(zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    SELECT CASE (nbondi)
    CASE (- 1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3we(:, jh, jk, jl, jf, 1) = ptab(iihom + jh, :, jk, jl)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(nn_hls + jh, :, jk, jl)
              zt3we(:, jh, jk, jl, jf, 1) = ptab(iihom + jh, :, jk, jl)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(nn_hls + jh, :, jk, jl)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpj * ipk * ipl * ipf
    IF (ln_timing) CALL tic_tac(.TRUE.)
    SELECT CASE (nbondi)
    CASE (- 1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, noea)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 1), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    IF (ln_timing) CALL tic_tac(.FALSE.)
    iihom = nlci - nn_hls
    SELECT CASE (nbondi)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(iihom + jh, :, jk, jl) = zt3ew(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :, jk, jl) = zt3we(:, jh, jk, jl, jf, 2)
              ptab(iihom + jh, :, jk, jl) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :, jk, jl) = zt3we(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    IF (nbondi /= 2) DEALLOCATE(zt3ew, zt3we)
    IF (npolj /= 0 .AND. .NOT. PRESENT(cd_mpp)) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat, psgn)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat, psgn)
      END SELECT
    END IF
    IF (ABS(nbondj) == 1) ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 1), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 1))
    IF (nbondj == 0) ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2))
    SELECT CASE (nbondj)
    CASE (- 1)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(:, ijhom + jh, jk, jl)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(:, ijhom + jh, jk, jl)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(:, nn_hls + jh, jk, jl)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(:, nn_hls + jh, jk, jl)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpi * ipk * ipl * ipf
    IF (ln_timing) CALL tic_tac(.TRUE.)
    SELECT CASE (nbondj)
    CASE (- 1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, nono)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    IF (ln_timing) CALL tic_tac(.FALSE.)
    ijhom = nlcj - nn_hls
    SELECT CASE (nbondj)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, ijhom + jh, jk, jl) = zt3ns(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh, jk, jl) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(:, ijhom + jh, jk, jl) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh, jk, jl) = zt3sn(:, jh, jk, jl, jf, 1)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    IF (nbondj /= 2) DEALLOCATE(zt3ns, zt3sn)
  END SUBROUTINE mpp_lnk_4d
  SUBROUTINE mpp_lnk_4d_ptr(cdname, ptab, cd_nat, psgn, kfld, cd_mpp, pval)
    INTEGER, INTENT(IN) :: kfld
    TYPE(PTR_4D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    INTEGER :: ierr
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = SIZE(ptab(1) % pt4d, 3)
    ipl = SIZE(ptab(1) % pt4d, 4)
    ipf = kfld
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ipk, ipl, ipf, ld_lbc = .TRUE.)
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    IF (.NOT. PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          ptab(jf) % pt4d(1, :, :, :) = ptab(jf) % pt4d(jpim1, :, :, :)
          ptab(jf) % pt4d(jpi, :, :, :) = ptab(jf) % pt4d(2, :, :, :)
        ELSE
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt4d(1 : nn_hls, :, :, :) = zland
          ptab(jf) % pt4d(nlci - nn_hls + 1 : jpi, :, :, :) = zland
        END IF
        IF (l_Jperio) THEN
          ptab(jf) % pt4d(:, 1, :, :) = ptab(jf) % pt4d(:, jpjm1, :, :)
          ptab(jf) % pt4d(:, jpj, :, :) = ptab(jf) % pt4d(:, 2, :, :)
        ELSE
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt4d(:, 1 : nn_hls, :, :) = zland
          ptab(jf) % pt4d(:, nlcj - nn_hls + 1 : jpj, :, :) = zland
        END IF
      END DO
    END IF
    IF (ABS(nbondi) == 1) ALLOCATE(zt3ew(jpj, nn_hls, ipk, ipl, ipf, 1), zt3we(jpj, nn_hls, ipk, ipl, ipf, 1))
    IF (nbondi == 0) ALLOCATE(zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    SELECT CASE (nbondi)
    CASE (- 1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3we(:, jh, jk, jl, jf, 1) = ptab(jf) % pt4d(iihom + jh, :, jk, jl)
            END DO
          END DO
        END DO
      END DO
    CASE (0)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(jf) % pt4d(nn_hls + jh, :, jk, jl)
              zt3we(:, jh, jk, jl, jf, 1) = ptab(jf) % pt4d(iihom + jh, :, jk, jl)
            END DO
          END DO
        END DO
      END DO
    CASE (1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(jf) % pt4d(nn_hls + jh, :, jk, jl)
            END DO
          END DO
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpj * ipk * ipl * ipf
    IF (ln_timing) CALL tic_tac(.TRUE.)
    SELECT CASE (nbondi)
    CASE (- 1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, noea)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 1), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    IF (ln_timing) CALL tic_tac(.FALSE.)
    iihom = nlci - nn_hls
    SELECT CASE (nbondi)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt4d(iihom + jh, :, jk, jl) = zt3ew(:, jh, jk, jl, jf, 1)
            END DO
          END DO
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt4d(jh, :, jk, jl) = zt3we(:, jh, jk, jl, jf, 2)
              ptab(jf) % pt4d(iihom + jh, :, jk, jl) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt4d(jh, :, jk, jl) = zt3we(:, jh, jk, jl, jf, 1)
            END DO
          END DO
        END DO
      END DO
    END SELECT
    IF (nbondi /= 2) DEALLOCATE(zt3ew, zt3we)
    IF (npolj /= 0 .AND. .NOT. PRESENT(cd_mpp)) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat(:), psgn(:), ipf)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat(:), psgn(:), ipf)
      END SELECT
    END IF
    IF (ABS(nbondj) == 1) ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 1), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 1))
    IF (nbondj == 0) ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2))
    SELECT CASE (nbondj)
    CASE (- 1)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(jf) % pt4d(:, ijhom + jh, jk, jl)
            END DO
          END DO
        END DO
      END DO
    CASE (0)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(jf) % pt4d(:, ijhom + jh, jk, jl)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(jf) % pt4d(:, nn_hls + jh, jk, jl)
            END DO
          END DO
        END DO
      END DO
    CASE (1)
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(jf) % pt4d(:, nn_hls + jh, jk, jl)
            END DO
          END DO
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpi * ipk * ipl * ipf
    IF (ln_timing) CALL tic_tac(.TRUE.)
    SELECT CASE (nbondj)
    CASE (- 1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, nono)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    IF (ln_timing) CALL tic_tac(.FALSE.)
    ijhom = nlcj - nn_hls
    SELECT CASE (nbondj)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt4d(:, ijhom + jh, jk, jl) = zt3ns(:, jh, jk, jl, jf, 1)
            END DO
          END DO
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt4d(:, jh, jk, jl) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(jf) % pt4d(:, ijhom + jh, jk, jl) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt4d(:, jh, jk, jl) = zt3sn(:, jh, jk, jl, jf, 1)
            END DO
          END DO
        END DO
      END DO
    END SELECT
    IF (nbondj /= 2) DEALLOCATE(zt3ns, zt3sn)
  END SUBROUTINE mpp_lnk_4d_ptr
  SUBROUTINE mpp_nfd_2d(ptab, cd_nat, psgn, kfld)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf, jr
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ierr, ibuffsize, ilci, ildi, ilei, iilb
    INTEGER :: ij, iproc
    INTEGER, DIMENSION(jpmaxngh) :: ml_req_nf
    INTEGER :: ml_err
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    INTEGER :: ipf_j
    INTEGER :: js
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: jj_s
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipj_s
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ztabl
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: ztab, ztabr
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: znorthloc, zfoldwk
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: znorthgloio
    ipk = 1
    ipl = 1
    ipf = 1
    IF (l_north_nogather) THEN
      ALLOCATE(ipj_s(ipf))
      !$ACC KERNELS
      ipj = 2
      ipj_s(:) = 1
      !$ACC END KERNELS
      ALLOCATE(jj_s(ipf, 2))
      !$ACC KERNELS
      l_full_nf_update = .TRUE.
      IF (l_full_nf_update .OR. (ncom_stp == nit000 .AND. .NOT. ln_rstart)) ipj_s(:) = 2
      DO jf = 1, ipf
        SELECT CASE (npolj)
        CASE (3, 4)
          SELECT CASE (cd_nat)
          CASE ('T', 'W', 'U')
            jj_s(jf, 1) = nlcj - 2
            jj_s(jf, 2) = nlcj - 1
          CASE ('V', 'F')
            jj_s(jf, 1) = nlcj - 3
            jj_s(jf, 2) = nlcj - 2
          END SELECT
        CASE (5, 6)
          SELECT CASE (cd_nat)
          CASE ('T', 'W', 'U')
            jj_s(jf, 1) = nlcj - 1
            ipj_s(jf) = 1
          CASE ('V', 'F')
            jj_s(jf, 1) = nlcj - 2
            jj_s(jf, 2) = nlcj - 1
          END SELECT
        END SELECT
      END DO
      !$ACC END KERNELS
      ipf_j = SUM(ipj_s(:))
      ALLOCATE(znorthloc(jpimax, ipf_j, ipk, ipl, 1))
      js = 0
      DO jf = 1, ipf
        DO jj = 1, ipj_s(jf)
          js = js + 1
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              znorthloc(1 : jpi, js, jk, jl, 1) = ptab(1 : jpi, jj_s(jf, jj))
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      ibuffsize = jpimax * ipf_j * ipk * ipl
      ALLOCATE(zfoldwk(jpimax, ipf_j, ipk, ipl, 1))
      ALLOCATE(ztabr(jpimax * jpmaxngh, ipj, ipk, ipl, ipf))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztabr(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      IF (ln_timing) CALL tic_tac(.TRUE.)
      DO jr = 1, nsndto
        IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
          CALL mppsend(5, znorthloc, ibuffsize, nfipproc(isendto(jr), jpnj), ml_req_nf(jr))
        END IF
      END DO
      DO jr = 1, nsndto
        iproc = nfipproc(isendto(jr), jpnj)
        IF (iproc /= - 1) THEN
          iilb = nimppt(iproc + 1)
          ilci = nlcit(iproc + 1)
          ildi = nldit(iproc + 1)
          ilei = nleit(iproc + 1)
          IF (iilb == 1) ildi = 1
          IF (iilb + ilci - 1 == jpiglo) ilei = ilci
          iilb = nfiimpp(isendto(jr), jpnj) - nfiimpp(isendto(1), jpnj)
        END IF
        IF (iproc /= narea - 1 .AND. iproc /= - 1) THEN
          CALL mpprecv(5, zfoldwk, ibuffsize, iproc)
          js = 0
          DO jf = 1, ipf
            DO jj = 1, ipj_s(jf)
              js = js + 1
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = zfoldwk(ji, js, jk, jl, 1)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END DO
          END DO
        ELSE IF (iproc == narea - 1) THEN
          DO jf = 1, ipf
            DO jj = 1, ipj_s(jf)
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = ptab(ji, jj_s(jf, jj))
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END DO
          END DO
        END IF
      END DO
      IF (l_isend) THEN
        DO jr = 1, nsndto
          IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
            CALL mpi_wait(ml_req_nf(jr), ml_stat, ml_err)
          END IF
        END DO
      END IF
      IF (ln_timing) CALL tic_tac(.FALSE.)
      DO jf = 1, ipf
        CALL lbc_nfd_nogather(ptab(:, :), ztabr(:, 1 : ipj_s(jf), :, :, jf), cd_nat, psgn)
      END DO
      DEALLOCATE(zfoldwk)
      DEALLOCATE(ztabr)
      DEALLOCATE(jj_s)
      DEALLOCATE(ipj_s)
    ELSE
      ipj = 4
      ALLOCATE(znorthloc(jpimax, ipj, ipk, ipl, ipf))
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              znorthloc(1 : jpi, ij, jk, jl, jf) = ptab(1 : jpi, jj)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      ibuffsize = jpimax * ipj * ipk * ipl * ipf
      ALLOCATE(ztab(jpiglo, ipj, ipk, ipl, ipf))
      ALLOCATE(znorthgloio(jpimax, ipj, ipk, ipl, ipf, jpni))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztab(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      IF (ln_timing) CALL tic_tac(.TRUE.)
      CALL MPI_ALLGATHER(znorthloc, ibuffsize, MPI_DOUBLE_PRECISION, znorthgloio, ibuffsize, MPI_DOUBLE_PRECISION, ncomm_north, &
&ierr)
      IF (ln_timing) CALL tic_tac(.FALSE.)
      DO jr = 1, ndim_rank_north
        iproc = nrank_north(jr) + 1
        iilb = nimppt(iproc)
        ilci = nlcit(iproc)
        ildi = nldit(iproc)
        ilei = nleit(iproc)
        IF (iilb == 1) ildi = 1
        IF (iilb + ilci - 1 == jpiglo) ilei = ilci
        DO jf = 1, ipf
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              !$ACC LOOP INDEPENDENT COLLAPSE(2)
              DO jj = 1, ipj
                DO ji = ildi, ilei
                  ztab(ji + iilb - 1, jj, jk, jl, jf) = znorthgloio(ji, jj, jk, jl, jf, jr)
                END DO
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      DO jf = 1, ipf
        CALL lbc_nfd(ztab(:, :, :, :, jf), cd_nat, psgn)
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(ji, jj) = ztab(ji + nimpp - 1, ij, jk, jl, jf)
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      DEALLOCATE(ztab)
      DEALLOCATE(znorthgloio)
    END IF
    DEALLOCATE(znorthloc)
  END SUBROUTINE mpp_nfd_2d
  SUBROUTINE mpp_nfd_2d_ptr(ptab, cd_nat, psgn, kfld)
    TYPE(PTR_2D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf, jr
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ierr, ibuffsize, ilci, ildi, ilei, iilb
    INTEGER :: ij, iproc
    INTEGER, DIMENSION(jpmaxngh) :: ml_req_nf
    INTEGER :: ml_err
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    INTEGER :: ipf_j
    INTEGER :: js
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: jj_s
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipj_s
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ztabl
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: ztab, ztabr
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: znorthloc, zfoldwk
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: znorthgloio
    ipk = 1
    ipl = 1
    ipf = kfld
    IF (l_north_nogather) THEN
      ALLOCATE(ipj_s(ipf))
      !$ACC KERNELS
      ipj = 2
      ipj_s(:) = 1
      !$ACC END KERNELS
      ALLOCATE(jj_s(ipf, 2))
      !$ACC KERNELS
      l_full_nf_update = .TRUE.
      IF (l_full_nf_update .OR. (ncom_stp == nit000 .AND. .NOT. ln_rstart)) ipj_s(:) = 2
      !$ACC END KERNELS
      DO jf = 1, ipf
        SELECT CASE (npolj)
        CASE (3, 4)
          SELECT CASE (cd_nat(jf))
          CASE ('T', 'W', 'U')
            jj_s(jf, 1) = nlcj - 2
            jj_s(jf, 2) = nlcj - 1
          CASE ('V', 'F')
            jj_s(jf, 1) = nlcj - 3
            jj_s(jf, 2) = nlcj - 2
          END SELECT
        CASE (5, 6)
          SELECT CASE (cd_nat(jf))
          CASE ('T', 'W', 'U')
            jj_s(jf, 1) = nlcj - 1
            ipj_s(jf) = 1
          CASE ('V', 'F')
            jj_s(jf, 1) = nlcj - 2
            jj_s(jf, 2) = nlcj - 1
          END SELECT
        END SELECT
      END DO
      ipf_j = SUM(ipj_s(:))
      ALLOCATE(znorthloc(jpimax, ipf_j, ipk, ipl, 1))
      js = 0
      DO jf = 1, ipf
        DO jj = 1, ipj_s(jf)
          js = js + 1
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              znorthloc(1 : jpi, js, jk, jl, 1) = ptab(jf) % pt2d(1 : jpi, jj_s(jf, jj))
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      ibuffsize = jpimax * ipf_j * ipk * ipl
      ALLOCATE(zfoldwk(jpimax, ipf_j, ipk, ipl, 1))
      ALLOCATE(ztabr(jpimax * jpmaxngh, ipj, ipk, ipl, ipf))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztabr(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      IF (ln_timing) CALL tic_tac(.TRUE.)
      DO jr = 1, nsndto
        IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
          CALL mppsend(5, znorthloc, ibuffsize, nfipproc(isendto(jr), jpnj), ml_req_nf(jr))
        END IF
      END DO
      DO jr = 1, nsndto
        iproc = nfipproc(isendto(jr), jpnj)
        IF (iproc /= - 1) THEN
          iilb = nimppt(iproc + 1)
          ilci = nlcit(iproc + 1)
          ildi = nldit(iproc + 1)
          ilei = nleit(iproc + 1)
          IF (iilb == 1) ildi = 1
          IF (iilb + ilci - 1 == jpiglo) ilei = ilci
          iilb = nfiimpp(isendto(jr), jpnj) - nfiimpp(isendto(1), jpnj)
        END IF
        IF (iproc /= narea - 1 .AND. iproc /= - 1) THEN
          CALL mpprecv(5, zfoldwk, ibuffsize, iproc)
          js = 0
          DO jf = 1, ipf
            DO jj = 1, ipj_s(jf)
              js = js + 1
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = zfoldwk(ji, js, jk, jl, 1)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END DO
          END DO
        ELSE IF (iproc == narea - 1) THEN
          DO jf = 1, ipf
            DO jj = 1, ipj_s(jf)
              DO jl = 1, ipl
                 !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = ptab(jf) % pt2d(ji, jj_s(jf, jj))
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END DO
          END DO
        END IF
      END DO
      IF (l_isend) THEN
        DO jr = 1, nsndto
          IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
            CALL mpi_wait(ml_req_nf(jr), ml_stat, ml_err)
          END IF
        END DO
      END IF
      IF (ln_timing) CALL tic_tac(.FALSE.)
      DO jf = 1, ipf
        CALL lbc_nfd_nogather(ptab(jf) % pt2d(:, :), ztabr(:, 1 : ipj_s(jf), :, :, jf), cd_nat(jf), psgn(jf))
      END DO
      DEALLOCATE(zfoldwk)
      DEALLOCATE(ztabr)
      DEALLOCATE(jj_s)
      DEALLOCATE(ipj_s)
    ELSE
      ipj = 4
      ALLOCATE(znorthloc(jpimax, ipj, ipk, ipl, ipf))
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              znorthloc(1 : jpi, ij, jk, jl, jf) = ptab(jf) % pt2d(1 : jpi, jj)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      ibuffsize = jpimax * ipj * ipk * ipl * ipf
      ALLOCATE(ztab(jpiglo, ipj, ipk, ipl, ipf))
      ALLOCATE(znorthgloio(jpimax, ipj, ipk, ipl, ipf, jpni))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztab(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      IF (ln_timing) CALL tic_tac(.TRUE.)
      CALL MPI_ALLGATHER(znorthloc, ibuffsize, MPI_DOUBLE_PRECISION, znorthgloio, ibuffsize, MPI_DOUBLE_PRECISION, ncomm_north, &
&ierr)
      IF (ln_timing) CALL tic_tac(.FALSE.)
      DO jr = 1, ndim_rank_north
        iproc = nrank_north(jr) + 1
        iilb = nimppt(iproc)
        ilci = nlcit(iproc)
        ildi = nldit(iproc)
        ilei = nleit(iproc)
        IF (iilb == 1) ildi = 1
        IF (iilb + ilci - 1 == jpiglo) ilei = ilci
        DO jf = 1, ipf
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              !$ACC LOOP INDEPENDENT COLLAPSE(2)
              DO jj = 1, ipj
                DO ji = ildi, ilei
                  ztab(ji + iilb - 1, jj, jk, jl, jf) = znorthgloio(ji, jj, jk, jl, jf, jr)
                END DO
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      DO jf = 1, ipf
        CALL lbc_nfd(ztab(:, :, :, :, jf), cd_nat(jf), psgn(jf))
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(jf) % pt2d(ji, jj) = ztab(ji + nimpp - 1, ij, jk, jl, jf)
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      DEALLOCATE(ztab)
      DEALLOCATE(znorthgloio)
    END IF
    DEALLOCATE(znorthloc)
  END SUBROUTINE mpp_nfd_2d_ptr
  SUBROUTINE mpp_nfd_3d(ptab, cd_nat, psgn, kfld)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf, jr
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ierr, ibuffsize, ilci, ildi, ilei, iilb
    INTEGER :: ij, iproc
    INTEGER, DIMENSION(jpmaxngh) :: ml_req_nf
    INTEGER :: ml_err
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    INTEGER :: ipf_j
    INTEGER :: js
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: jj_s
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipj_s
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ztabl
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: ztab, ztabr
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: znorthloc, zfoldwk
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: znorthgloio
    ipk = SIZE(ptab, 3)
    ipl = 1
    ipf = 1
    IF (l_north_nogather) THEN
      ALLOCATE(ipj_s(ipf))
      !$ACC KERNELS
      ipj = 2
      ipj_s(:) = 1
      !$ACC END KERNELS
      ALLOCATE(jj_s(ipf, 2))
      !$ACC KERNELS
      l_full_nf_update = .TRUE.
      IF (l_full_nf_update .OR. (ncom_stp == nit000 .AND. .NOT. ln_rstart)) ipj_s(:) = 2
      DO jf = 1, ipf
        SELECT CASE (npolj)
        CASE (3, 4)
          SELECT CASE (cd_nat)
          CASE ('T', 'W', 'U')
            jj_s(jf, 1) = nlcj - 2
            jj_s(jf, 2) = nlcj - 1
          CASE ('V', 'F')
            jj_s(jf, 1) = nlcj - 3
            jj_s(jf, 2) = nlcj - 2
          END SELECT
        CASE (5, 6)
          SELECT CASE (cd_nat)
          CASE ('T', 'W', 'U')
            jj_s(jf, 1) = nlcj - 1
            ipj_s(jf) = 1
          CASE ('V', 'F')
            jj_s(jf, 1) = nlcj - 2
            jj_s(jf, 2) = nlcj - 1
          END SELECT
        END SELECT
      END DO
      !$ACC END KERNELS
      ipf_j = SUM(ipj_s(:))
      ALLOCATE(znorthloc(jpimax, ipf_j, ipk, ipl, 1))
      js = 0
      DO jf = 1, ipf
        DO jj = 1, ipj_s(jf)
          js = js + 1
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              znorthloc(1 : jpi, js, jk, jl, 1) = ptab(1 : jpi, jj_s(jf, jj), jk)
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      ibuffsize = jpimax * ipf_j * ipk * ipl
      ALLOCATE(zfoldwk(jpimax, ipf_j, ipk, ipl, 1))
      ALLOCATE(ztabr(jpimax * jpmaxngh, ipj, ipk, ipl, ipf))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztabr(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      IF (ln_timing) CALL tic_tac(.TRUE.)
      DO jr = 1, nsndto
        IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
          CALL mppsend(5, znorthloc, ibuffsize, nfipproc(isendto(jr), jpnj), ml_req_nf(jr))
        END IF
      END DO
      DO jr = 1, nsndto
        iproc = nfipproc(isendto(jr), jpnj)
        IF (iproc /= - 1) THEN
          iilb = nimppt(iproc + 1)
          ilci = nlcit(iproc + 1)
          ildi = nldit(iproc + 1)
          ilei = nleit(iproc + 1)
          IF (iilb == 1) ildi = 1
          IF (iilb + ilci - 1 == jpiglo) ilei = ilci
          iilb = nfiimpp(isendto(jr), jpnj) - nfiimpp(isendto(1), jpnj)
        END IF
        IF (iproc /= narea - 1 .AND. iproc /= - 1) THEN
          CALL mpprecv(5, zfoldwk, ibuffsize, iproc)
          js = 0
          DO jf = 1, ipf
            DO jj = 1, ipj_s(jf)
              js = js + 1
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = zfoldwk(ji, js, jk, jl, 1)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END DO
          END DO
        ELSE IF (iproc == narea - 1) THEN
          DO jf = 1, ipf
            DO jj = 1, ipj_s(jf)
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = ptab(ji, jj_s(jf, jj), jk)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END DO
          END DO
        END IF
      END DO
      IF (l_isend) THEN
        DO jr = 1, nsndto
          IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
            CALL mpi_wait(ml_req_nf(jr), ml_stat, ml_err)
          END IF
        END DO
      END IF
      IF (ln_timing) CALL tic_tac(.FALSE.)
      DO jf = 1, ipf
        CALL lbc_nfd_nogather(ptab(:, :, :), ztabr(:, 1 : ipj_s(jf), :, :, jf), cd_nat, psgn)
      END DO
      DEALLOCATE(zfoldwk)
      DEALLOCATE(ztabr)
      DEALLOCATE(jj_s)
      DEALLOCATE(ipj_s)
    ELSE
      ipj = 4
      ALLOCATE(znorthloc(jpimax, ipj, ipk, ipl, ipf))
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              znorthloc(1 : jpi, ij, jk, jl, jf) = ptab(1 : jpi, jj, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      ibuffsize = jpimax * ipj * ipk * ipl * ipf
      ALLOCATE(ztab(jpiglo, ipj, ipk, ipl, ipf))
      ALLOCATE(znorthgloio(jpimax, ipj, ipk, ipl, ipf, jpni))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztab(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      IF (ln_timing) CALL tic_tac(.TRUE.)
      CALL MPI_ALLGATHER(znorthloc, ibuffsize, MPI_DOUBLE_PRECISION, znorthgloio, ibuffsize, MPI_DOUBLE_PRECISION, ncomm_north, &
&ierr)
      IF (ln_timing) CALL tic_tac(.FALSE.)
      DO jr = 1, ndim_rank_north
        iproc = nrank_north(jr) + 1
        iilb = nimppt(iproc)
        ilci = nlcit(iproc)
        ildi = nldit(iproc)
        ilei = nleit(iproc)
        IF (iilb == 1) ildi = 1
        IF (iilb + ilci - 1 == jpiglo) ilei = ilci
        DO jf = 1, ipf
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              !$ACC LOOP INDEPENDENT COLLAPSE(2)
              DO jj = 1, ipj
                DO ji = ildi, ilei
                  ztab(ji + iilb - 1, jj, jk, jl, jf) = znorthgloio(ji, jj, jk, jl, jf, jr)
                END DO
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      DO jf = 1, ipf
        CALL lbc_nfd(ztab(:, :, :, :, jf), cd_nat, psgn)
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(ji, jj, jk) = ztab(ji + nimpp - 1, ij, jk, jl, jf)
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      DEALLOCATE(ztab)
      DEALLOCATE(znorthgloio)
    END IF
    DEALLOCATE(znorthloc)
  END SUBROUTINE mpp_nfd_3d
  SUBROUTINE mpp_nfd_3d_ptr(ptab, cd_nat, psgn, kfld)
    TYPE(PTR_3D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf, jr
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ierr, ibuffsize, ilci, ildi, ilei, iilb
    INTEGER :: ij, iproc
    INTEGER, DIMENSION(jpmaxngh) :: ml_req_nf
    INTEGER :: ml_err
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    INTEGER :: ipf_j
    INTEGER :: js
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: jj_s
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipj_s
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ztabl
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: ztab, ztabr
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: znorthloc, zfoldwk
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: znorthgloio
    ipk = SIZE(ptab(1) % pt3d, 3)
    ipl = 1
    ipf = kfld
    IF (l_north_nogather) THEN
      ALLOCATE(ipj_s(ipf))
      !$ACC KERNELS
      ipj = 2
      ipj_s(:) = 1
      !$ACC END KERNELS
      ALLOCATE(jj_s(ipf, 2))
      !$ACC KERNELS
      l_full_nf_update = .TRUE.
      IF (l_full_nf_update .OR. (ncom_stp == nit000 .AND. .NOT. ln_rstart)) ipj_s(:) = 2
      !$ACC END KERNELS
      DO jf = 1, ipf
        SELECT CASE (npolj)
        CASE (3, 4)
          SELECT CASE (cd_nat(jf))
          CASE ('T', 'W', 'U')
            jj_s(jf, 1) = nlcj - 2
            jj_s(jf, 2) = nlcj - 1
          CASE ('V', 'F')
            jj_s(jf, 1) = nlcj - 3
            jj_s(jf, 2) = nlcj - 2
          END SELECT
        CASE (5, 6)
          SELECT CASE (cd_nat(jf))
          CASE ('T', 'W', 'U')
            jj_s(jf, 1) = nlcj - 1
            ipj_s(jf) = 1
          CASE ('V', 'F')
            jj_s(jf, 1) = nlcj - 2
            jj_s(jf, 2) = nlcj - 1
          END SELECT
        END SELECT
      END DO
      ipf_j = SUM(ipj_s(:))
      ALLOCATE(znorthloc(jpimax, ipf_j, ipk, ipl, 1))
      js = 0
      DO jf = 1, ipf
        DO jj = 1, ipj_s(jf)
          js = js + 1
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              znorthloc(1 : jpi, js, jk, jl, 1) = ptab(jf) % pt3d(1 : jpi, jj_s(jf, jj), jk)
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      ibuffsize = jpimax * ipf_j * ipk * ipl
      ALLOCATE(zfoldwk(jpimax, ipf_j, ipk, ipl, 1))
      ALLOCATE(ztabr(jpimax * jpmaxngh, ipj, ipk, ipl, ipf))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztabr(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      IF (ln_timing) CALL tic_tac(.TRUE.)
      DO jr = 1, nsndto
        IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
          CALL mppsend(5, znorthloc, ibuffsize, nfipproc(isendto(jr), jpnj), ml_req_nf(jr))
        END IF
      END DO
      DO jr = 1, nsndto
        iproc = nfipproc(isendto(jr), jpnj)
        IF (iproc /= - 1) THEN
          iilb = nimppt(iproc + 1)
          ilci = nlcit(iproc + 1)
          ildi = nldit(iproc + 1)
          ilei = nleit(iproc + 1)
          IF (iilb == 1) ildi = 1
          IF (iilb + ilci - 1 == jpiglo) ilei = ilci
          iilb = nfiimpp(isendto(jr), jpnj) - nfiimpp(isendto(1), jpnj)
        END IF
        IF (iproc /= narea - 1 .AND. iproc /= - 1) THEN
          CALL mpprecv(5, zfoldwk, ibuffsize, iproc)
          js = 0
          DO jf = 1, ipf
            DO jj = 1, ipj_s(jf)
              js = js + 1
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = zfoldwk(ji, js, jk, jl, 1)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END DO
          END DO
        ELSE IF (iproc == narea - 1) THEN
          DO jf = 1, ipf
            DO jj = 1, ipj_s(jf)
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = ptab(jf) % pt3d(ji, jj_s(jf, jj), jk)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END DO
          END DO
        END IF
      END DO
      IF (l_isend) THEN
        DO jr = 1, nsndto
          IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
            CALL mpi_wait(ml_req_nf(jr), ml_stat, ml_err)
          END IF
        END DO
      END IF
      IF (ln_timing) CALL tic_tac(.FALSE.)
      DO jf = 1, ipf
        CALL lbc_nfd_nogather(ptab(jf) % pt3d(:, :, :), ztabr(:, 1 : ipj_s(jf), :, :, jf), cd_nat(jf), psgn(jf))
      END DO
      DEALLOCATE(zfoldwk)
      DEALLOCATE(ztabr)
      DEALLOCATE(jj_s)
      DEALLOCATE(ipj_s)
    ELSE
      ipj = 4
      ALLOCATE(znorthloc(jpimax, ipj, ipk, ipl, ipf))
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              znorthloc(1 : jpi, ij, jk, jl, jf) = ptab(jf) % pt3d(1 : jpi, jj, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      ibuffsize = jpimax * ipj * ipk * ipl * ipf
      ALLOCATE(ztab(jpiglo, ipj, ipk, ipl, ipf))
      ALLOCATE(znorthgloio(jpimax, ipj, ipk, ipl, ipf, jpni))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztab(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      IF (ln_timing) CALL tic_tac(.TRUE.)
      CALL MPI_ALLGATHER(znorthloc, ibuffsize, MPI_DOUBLE_PRECISION, znorthgloio, ibuffsize, MPI_DOUBLE_PRECISION, ncomm_north, &
&ierr)
      IF (ln_timing) CALL tic_tac(.FALSE.)
      DO jr = 1, ndim_rank_north
        iproc = nrank_north(jr) + 1
        iilb = nimppt(iproc)
        ilci = nlcit(iproc)
        ildi = nldit(iproc)
        ilei = nleit(iproc)
        IF (iilb == 1) ildi = 1
        IF (iilb + ilci - 1 == jpiglo) ilei = ilci
        DO jf = 1, ipf
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              !$ACC LOOP INDEPENDENT COLLAPSE(2)
              DO jj = 1, ipj
                DO ji = ildi, ilei
                  ztab(ji + iilb - 1, jj, jk, jl, jf) = znorthgloio(ji, jj, jk, jl, jf, jr)
                END DO
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      DO jf = 1, ipf
        CALL lbc_nfd(ztab(:, :, :, :, jf), cd_nat(jf), psgn(jf))
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(jf) % pt3d(ji, jj, jk) = ztab(ji + nimpp - 1, ij, jk, jl, jf)
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      DEALLOCATE(ztab)
      DEALLOCATE(znorthgloio)
    END IF
    DEALLOCATE(znorthloc)
  END SUBROUTINE mpp_nfd_3d_ptr
  SUBROUTINE mpp_nfd_4d(ptab, cd_nat, psgn, kfld)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf, jr
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ierr, ibuffsize, ilci, ildi, ilei, iilb
    INTEGER :: ij, iproc
    INTEGER, DIMENSION(jpmaxngh) :: ml_req_nf
    INTEGER :: ml_err
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    INTEGER :: ipf_j
    INTEGER :: js
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: jj_s
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipj_s
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ztabl
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: ztab, ztabr
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: znorthloc, zfoldwk
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: znorthgloio
    ipk = SIZE(ptab, 3)
    ipl = SIZE(ptab, 4)
    ipf = 1
    IF (l_north_nogather) THEN
      ALLOCATE(ipj_s(ipf))
      !$ACC KERNELS
      ipj = 2
      ipj_s(:) = 1
      !$ACC END KERNELS
      ALLOCATE(jj_s(ipf, 2))
      !$ACC KERNELS
      l_full_nf_update = .TRUE.
      IF (l_full_nf_update .OR. (ncom_stp == nit000 .AND. .NOT. ln_rstart)) ipj_s(:) = 2
      DO jf = 1, ipf
        SELECT CASE (npolj)
        CASE (3, 4)
          SELECT CASE (cd_nat)
          CASE ('T', 'W', 'U')
            jj_s(jf, 1) = nlcj - 2
            jj_s(jf, 2) = nlcj - 1
          CASE ('V', 'F')
            jj_s(jf, 1) = nlcj - 3
            jj_s(jf, 2) = nlcj - 2
          END SELECT
        CASE (5, 6)
          SELECT CASE (cd_nat)
          CASE ('T', 'W', 'U')
            jj_s(jf, 1) = nlcj - 1
            ipj_s(jf) = 1
          CASE ('V', 'F')
            jj_s(jf, 1) = nlcj - 2
            jj_s(jf, 2) = nlcj - 1
          END SELECT
        END SELECT
      END DO
      !$ACC END KERNELS
      ipf_j = SUM(ipj_s(:))
      ALLOCATE(znorthloc(jpimax, ipf_j, ipk, ipl, 1))
      js = 0
      DO jf = 1, ipf
        DO jj = 1, ipj_s(jf)
          js = js + 1
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              znorthloc(1 : jpi, js, jk, jl, 1) = ptab(1 : jpi, jj_s(jf, jj), jk, jl)
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      ibuffsize = jpimax * ipf_j * ipk * ipl
      ALLOCATE(zfoldwk(jpimax, ipf_j, ipk, ipl, 1))
      ALLOCATE(ztabr(jpimax * jpmaxngh, ipj, ipk, ipl, ipf))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztabr(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      IF (ln_timing) CALL tic_tac(.TRUE.)
      DO jr = 1, nsndto
        IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
          CALL mppsend(5, znorthloc, ibuffsize, nfipproc(isendto(jr), jpnj), ml_req_nf(jr))
        END IF
      END DO
      DO jr = 1, nsndto
        iproc = nfipproc(isendto(jr), jpnj)
        IF (iproc /= - 1) THEN
          iilb = nimppt(iproc + 1)
          ilci = nlcit(iproc + 1)
          ildi = nldit(iproc + 1)
          ilei = nleit(iproc + 1)
          IF (iilb == 1) ildi = 1
          IF (iilb + ilci - 1 == jpiglo) ilei = ilci
          iilb = nfiimpp(isendto(jr), jpnj) - nfiimpp(isendto(1), jpnj)
        END IF
        IF (iproc /= narea - 1 .AND. iproc /= - 1) THEN
          CALL mpprecv(5, zfoldwk, ibuffsize, iproc)
          js = 0
          DO jf = 1, ipf
            DO jj = 1, ipj_s(jf)
              js = js + 1
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = zfoldwk(ji, js, jk, jl, 1)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END DO
          END DO
        ELSE IF (iproc == narea - 1) THEN
          DO jf = 1, ipf
            DO jj = 1, ipj_s(jf)
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = ptab(ji, jj_s(jf, jj), jk, jl)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END DO
          END DO
        END IF
      END DO
      IF (l_isend) THEN
        DO jr = 1, nsndto
          IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
            CALL mpi_wait(ml_req_nf(jr), ml_stat, ml_err)
          END IF
        END DO
      END IF
      IF (ln_timing) CALL tic_tac(.FALSE.)
      DO jf = 1, ipf
        CALL lbc_nfd_nogather(ptab(:, :, :, :), ztabr(:, 1 : ipj_s(jf), :, :, jf), cd_nat, psgn)
      END DO
      DEALLOCATE(zfoldwk)
      DEALLOCATE(ztabr)
      DEALLOCATE(jj_s)
      DEALLOCATE(ipj_s)
    ELSE
      ipj = 4
      ALLOCATE(znorthloc(jpimax, ipj, ipk, ipl, ipf))
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              znorthloc(1 : jpi, ij, jk, jl, jf) = ptab(1 : jpi, jj, jk, jl)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      ibuffsize = jpimax * ipj * ipk * ipl * ipf
      ALLOCATE(ztab(jpiglo, ipj, ipk, ipl, ipf))
      ALLOCATE(znorthgloio(jpimax, ipj, ipk, ipl, ipf, jpni))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztab(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      IF (ln_timing) CALL tic_tac(.TRUE.)
      CALL MPI_ALLGATHER(znorthloc, ibuffsize, MPI_DOUBLE_PRECISION, znorthgloio, ibuffsize, MPI_DOUBLE_PRECISION, ncomm_north, &
&ierr)
      IF (ln_timing) CALL tic_tac(.FALSE.)
      DO jr = 1, ndim_rank_north
        iproc = nrank_north(jr) + 1
        iilb = nimppt(iproc)
        ilci = nlcit(iproc)
        ildi = nldit(iproc)
        ilei = nleit(iproc)
        IF (iilb == 1) ildi = 1
        IF (iilb + ilci - 1 == jpiglo) ilei = ilci
        DO jf = 1, ipf
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              !$ACC LOOP INDEPENDENT COLLAPSE(2)
              DO jj = 1, ipj
                DO ji = ildi, ilei
                  ztab(ji + iilb - 1, jj, jk, jl, jf) = znorthgloio(ji, jj, jk, jl, jf, jr)
                END DO
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      DO jf = 1, ipf
        CALL lbc_nfd(ztab(:, :, :, :, jf), cd_nat, psgn)
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(ji, jj, jk, jl) = ztab(ji + nimpp - 1, ij, jk, jl, jf)
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      DEALLOCATE(ztab)
      DEALLOCATE(znorthgloio)
    END IF
    DEALLOCATE(znorthloc)
  END SUBROUTINE mpp_nfd_4d
  SUBROUTINE mpp_nfd_4d_ptr(ptab, cd_nat, psgn, kfld)
    TYPE(PTR_4D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf, jr
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ierr, ibuffsize, ilci, ildi, ilei, iilb
    INTEGER :: ij, iproc
    INTEGER, DIMENSION(jpmaxngh) :: ml_req_nf
    INTEGER :: ml_err
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    INTEGER :: ipf_j
    INTEGER :: js
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: jj_s
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipj_s
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ztabl
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: ztab, ztabr
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: znorthloc, zfoldwk
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: znorthgloio
    ipk = SIZE(ptab(1) % pt4d, 3)
    ipl = SIZE(ptab(1) % pt4d, 4)
    ipf = kfld
    IF (l_north_nogather) THEN
      ALLOCATE(ipj_s(ipf))
      !$ACC KERNELS
      ipj = 2
      ipj_s(:) = 1
      !$ACC END KERNELS
      ALLOCATE(jj_s(ipf, 2))
      !$ACC KERNELS
      l_full_nf_update = .TRUE.
      IF (l_full_nf_update .OR. (ncom_stp == nit000 .AND. .NOT. ln_rstart)) ipj_s(:) = 2
      !$ACC END KERNELS
      DO jf = 1, ipf
        SELECT CASE (npolj)
        CASE (3, 4)
          SELECT CASE (cd_nat(jf))
          CASE ('T', 'W', 'U')
            jj_s(jf, 1) = nlcj - 2
            jj_s(jf, 2) = nlcj - 1
          CASE ('V', 'F')
            jj_s(jf, 1) = nlcj - 3
            jj_s(jf, 2) = nlcj - 2
          END SELECT
        CASE (5, 6)
          SELECT CASE (cd_nat(jf))
          CASE ('T', 'W', 'U')
            jj_s(jf, 1) = nlcj - 1
            ipj_s(jf) = 1
          CASE ('V', 'F')
            jj_s(jf, 1) = nlcj - 2
            jj_s(jf, 2) = nlcj - 1
          END SELECT
        END SELECT
      END DO
      ipf_j = SUM(ipj_s(:))
      ALLOCATE(znorthloc(jpimax, ipf_j, ipk, ipl, 1))
      js = 0
      DO jf = 1, ipf
        DO jj = 1, ipj_s(jf)
          js = js + 1
          DO jl = 1, ipl
            DO jk = 1, ipk
              znorthloc(1 : jpi, js, jk, jl, 1) = ptab(jf) % pt4d(1 : jpi, jj_s(jf, jj), jk, jl)
            END DO
          END DO
        END DO
      END DO
      ibuffsize = jpimax * ipf_j * ipk * ipl
      ALLOCATE(zfoldwk(jpimax, ipf_j, ipk, ipl, 1))
      ALLOCATE(ztabr(jpimax * jpmaxngh, ipj, ipk, ipl, ipf))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztabr(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      IF (ln_timing) CALL tic_tac(.TRUE.)
      DO jr = 1, nsndto
        IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
          CALL mppsend(5, znorthloc, ibuffsize, nfipproc(isendto(jr), jpnj), ml_req_nf(jr))
        END IF
      END DO
      DO jr = 1, nsndto
        iproc = nfipproc(isendto(jr), jpnj)
        IF (iproc /= - 1) THEN
          iilb = nimppt(iproc + 1)
          ilci = nlcit(iproc + 1)
          ildi = nldit(iproc + 1)
          ilei = nleit(iproc + 1)
          IF (iilb == 1) ildi = 1
          IF (iilb + ilci - 1 == jpiglo) ilei = ilci
          iilb = nfiimpp(isendto(jr), jpnj) - nfiimpp(isendto(1), jpnj)
        END IF
        IF (iproc /= narea - 1 .AND. iproc /= - 1) THEN
          CALL mpprecv(5, zfoldwk, ibuffsize, iproc)
          js = 0
          DO jf = 1, ipf
            DO jj = 1, ipj_s(jf)
              js = js + 1
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = zfoldwk(ji, js, jk, jl, 1)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END DO
          END DO
        ELSE IF (iproc == narea - 1) THEN
          DO jf = 1, ipf
            DO jj = 1, ipj_s(jf)
              DO jl = 1, ipl
                DO jk = 1, ipk
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = ptab(jf) % pt4d(ji, jj_s(jf, jj), jk, jl)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END IF
      END DO
      IF (l_isend) THEN
        DO jr = 1, nsndto
          IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
            CALL mpi_wait(ml_req_nf(jr), ml_stat, ml_err)
          END IF
        END DO
      END IF
      IF (ln_timing) CALL tic_tac(.FALSE.)
      DO jf = 1, ipf
        CALL lbc_nfd_nogather(ptab(jf) % pt4d(:, :, :, :), ztabr(:, 1 : ipj_s(jf), :, :, jf), cd_nat(jf), psgn(jf))
      END DO
      DEALLOCATE(zfoldwk)
      DEALLOCATE(ztabr)
      DEALLOCATE(jj_s)
      DEALLOCATE(ipj_s)
    ELSE
      ipj = 4
      ALLOCATE(znorthloc(jpimax, ipj, ipk, ipl, ipf))
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              znorthloc(1 : jpi, ij, jk, jl, jf) = ptab(jf) % pt4d(1 : jpi, jj, jk, jl)
            END DO
          END DO
        END DO
      END DO
      ibuffsize = jpimax * ipj * ipk * ipl * ipf
      ALLOCATE(ztab(jpiglo, ipj, ipk, ipl, ipf))
      ALLOCATE(znorthgloio(jpimax, ipj, ipk, ipl, ipf, jpni))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztab(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      IF (ln_timing) CALL tic_tac(.TRUE.)
      CALL MPI_ALLGATHER(znorthloc, ibuffsize, MPI_DOUBLE_PRECISION, znorthgloio, ibuffsize, MPI_DOUBLE_PRECISION, ncomm_north, &
&ierr)
      IF (ln_timing) CALL tic_tac(.FALSE.)
      DO jr = 1, ndim_rank_north
        iproc = nrank_north(jr) + 1
        iilb = nimppt(iproc)
        ilci = nlcit(iproc)
        ildi = nldit(iproc)
        ilei = nleit(iproc)
        IF (iilb == 1) ildi = 1
        IF (iilb + ilci - 1 == jpiglo) ilei = ilci
        DO jf = 1, ipf
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              !$ACC LOOP INDEPENDENT COLLAPSE(2)
              DO jj = 1, ipj
                DO ji = ildi, ilei
                  ztab(ji + iilb - 1, jj, jk, jl, jf) = znorthgloio(ji, jj, jk, jl, jf, jr)
                END DO
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      DO jf = 1, ipf
        CALL lbc_nfd(ztab(:, :, :, :, jf), cd_nat(jf), psgn(jf))
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(jf) % pt4d(ji, jj, jk, jl) = ztab(ji + nimpp - 1, ij, jk, jl, jf)
              END DO
            END DO
          END DO
        END DO
      END DO
      DEALLOCATE(ztab)
      DEALLOCATE(znorthgloio)
    END IF
    DEALLOCATE(znorthloc)
  END SUBROUTINE mpp_nfd_4d_ptr
  SUBROUTINE mpp_lnk_bdy_2d(cdname, ptab, cd_nat, psgn, kb_bdy)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, INTENT(IN) :: kb_bdy
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = 1
    ipl = 1
    ipf = 1
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ipk, ipl, ipf, ld_lbc = .TRUE.)
    ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2), zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), &
&zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    zland = 0._wp
    DO jf = 1, ipf
      IF (nbondi == 2) THEN
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(1, :) = ptab(jpim1, :)
          ptab(jpi, :) = ptab(2, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :) = zland
          ptab(nlci - nn_hls + 1 : jpi, :) = zland
          !$ACC END KERNELS
        END IF
      ELSE IF (nbondi == - 1) THEN
        !$ACC KERNELS
        IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :) = zland
        !$ACC END KERNELS
      ELSE IF (nbondi == 1) THEN
        ptab(nlci - nn_hls + 1 : jpi, :) = zland
      END IF
      IF (nbondj == 2) THEN
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(:, 1) = ptab(:, jpjm1)
          ptab(:, jpj) = ptab(:, 2)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls) = zland
          ptab(:, nlcj - nn_hls + 1 : jpj) = zland
          !$ACC END KERNELS
        END IF
      ELSE IF (nbondj == - 1) THEN
        !$ACC KERNELS
        IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls) = zland
        !$ACC END KERNELS
      ELSE IF (nbondj == 1) THEN
        ptab(:, nlcj - nn_hls + 1 : jpj) = zland
      END IF
    END DO
    DO jf = 1, ipf
      SELECT CASE (nbondi_bdy(kb_bdy))
      CASE (- 1, 0, 1)
        iihom = nlci - nreci
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(nn_hls + jh, :)
              zt3we(:, jh, jk, jl, jf, 1) = ptab(iihom + jh, :)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
      imigr = nn_hls * jpj * ipk * ipl
      IF (ln_timing) CALL tic_tac(.TRUE.)
      SELECT CASE (nbondi_bdy(kb_bdy))
      CASE (- 1)
        CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
      CASE (0)
        CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
        CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
      CASE (1)
        CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      END SELECT
      SELECT CASE (nbondi_bdy_b(kb_bdy))
      CASE (- 1)
        CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      CASE (0)
        CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
        CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      CASE (1)
        CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      END SELECT
      SELECT CASE (nbondi_bdy(kb_bdy))
      CASE (- 1)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE (0)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
        IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE (1)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      IF (ln_timing) CALL tic_tac(.FALSE.)
      iihom = nlci - nn_hls
      SELECT CASE (nbondi_bdy_b(kb_bdy))
      CASE (- 1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(iihom + jh, :) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (0)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :) = zt3we(:, jh, jk, jl, jf, 2)
              ptab(iihom + jh, :) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :) = zt3we(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
    END DO
    IF (npolj /= 0) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat, psgn)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat, psgn)
      END SELECT
    END IF
    DO jf = 1, ipf
      IF (nbondj_bdy(kb_bdy) /= 2) THEN
        ijhom = nlcj - nrecj
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(:, ijhom + jh)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(:, nn_hls + jh)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END IF
      imigr = nn_hls * jpi * ipk * ipl
      IF (ln_timing) CALL tic_tac(.TRUE.)
      SELECT CASE (nbondj_bdy(kb_bdy))
      CASE (- 1)
        CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
      CASE (0)
        CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
        CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
      CASE (1)
        CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      END SELECT
      SELECT CASE (nbondj_bdy_b(kb_bdy))
      CASE (- 1)
        CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      CASE (0)
        CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
        CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      CASE (1)
        CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      END SELECT
      SELECT CASE (nbondj_bdy(kb_bdy))
      CASE (- 1)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE (0)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
        IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE (1)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      IF (ln_timing) CALL tic_tac(.FALSE.)
      ijhom = nlcj - nn_hls
      SELECT CASE (nbondj_bdy_b(kb_bdy))
      CASE (- 1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, ijhom + jh) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (0)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(:, ijhom + jh) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh) = zt3sn(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
    END DO
    DEALLOCATE(zt3ns, zt3sn, zt3ew, zt3we)
  END SUBROUTINE mpp_lnk_bdy_2d
  SUBROUTINE mpp_lnk_bdy_3d(cdname, ptab, cd_nat, psgn, kb_bdy)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, INTENT(IN) :: kb_bdy
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = SIZE(ptab, 3)
    ipl = 1
    ipf = 1
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ipk, ipl, ipf, ld_lbc = .TRUE.)
    ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2), zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), &
&zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    zland = 0._wp
    DO jf = 1, ipf
      IF (nbondi == 2) THEN
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(1, :, :) = ptab(jpim1, :, :)
          ptab(jpi, :, :) = ptab(2, :, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :, :) = zland
          ptab(nlci - nn_hls + 1 : jpi, :, :) = zland
          !$ACC END KERNELS
        END IF
      ELSE IF (nbondi == - 1) THEN
        !$ACC KERNELS
        IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :, :) = zland
        !$ACC END KERNELS
      ELSE IF (nbondi == 1) THEN
        ptab(nlci - nn_hls + 1 : jpi, :, :) = zland
      END IF
      IF (nbondj == 2) THEN
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(:, 1, :) = ptab(:, jpjm1, :)
          ptab(:, jpj, :) = ptab(:, 2, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls, :) = zland
          ptab(:, nlcj - nn_hls + 1 : jpj, :) = zland
          !$ACC END KERNELS
        END IF
      ELSE IF (nbondj == - 1) THEN
        !$ACC KERNELS
        IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls, :) = zland
        !$ACC END KERNELS
      ELSE IF (nbondj == 1) THEN
        ptab(:, nlcj - nn_hls + 1 : jpj, :) = zland
      END IF
    END DO
    DO jf = 1, ipf
      SELECT CASE (nbondi_bdy(kb_bdy))
      CASE (- 1, 0, 1)
        iihom = nlci - nreci
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(nn_hls + jh, :, jk)
              zt3we(:, jh, jk, jl, jf, 1) = ptab(iihom + jh, :, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
      imigr = nn_hls * jpj * ipk * ipl
      IF (ln_timing) CALL tic_tac(.TRUE.)
      SELECT CASE (nbondi_bdy(kb_bdy))
      CASE (- 1)
        CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
      CASE (0)
        CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
        CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
      CASE (1)
        CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      END SELECT
      SELECT CASE (nbondi_bdy_b(kb_bdy))
      CASE (- 1)
        CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      CASE (0)
        CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
        CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      CASE (1)
        CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      END SELECT
      SELECT CASE (nbondi_bdy(kb_bdy))
      CASE (- 1)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE (0)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
        IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE (1)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      IF (ln_timing) CALL tic_tac(.FALSE.)
      iihom = nlci - nn_hls
      SELECT CASE (nbondi_bdy_b(kb_bdy))
      CASE (- 1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(iihom + jh, :, jk) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (0)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :, jk) = zt3we(:, jh, jk, jl, jf, 2)
              ptab(iihom + jh, :, jk) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :, jk) = zt3we(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
    END DO
    IF (npolj /= 0) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat, psgn)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat, psgn)
      END SELECT
    END IF
    DO jf = 1, ipf
      IF (nbondj_bdy(kb_bdy) /= 2) THEN
        ijhom = nlcj - nrecj
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(:, ijhom + jh, jk)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(:, nn_hls + jh, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END IF
      imigr = nn_hls * jpi * ipk * ipl
      IF (ln_timing) CALL tic_tac(.TRUE.)
      SELECT CASE (nbondj_bdy(kb_bdy))
      CASE (- 1)
        CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
      CASE (0)
        CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
        CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
      CASE (1)
        CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      END SELECT
      SELECT CASE (nbondj_bdy_b(kb_bdy))
      CASE (- 1)
        CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      CASE (0)
        CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
        CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      CASE (1)
        CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      END SELECT
      SELECT CASE (nbondj_bdy(kb_bdy))
      CASE (- 1)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE (0)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
        IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE (1)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      IF (ln_timing) CALL tic_tac(.FALSE.)
      ijhom = nlcj - nn_hls
      SELECT CASE (nbondj_bdy_b(kb_bdy))
      CASE (- 1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, ijhom + jh, jk) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (0)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh, jk) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(:, ijhom + jh, jk) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh, jk) = zt3sn(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
    END DO
    DEALLOCATE(zt3ns, zt3sn, zt3ew, zt3we)
  END SUBROUTINE mpp_lnk_bdy_3d
  SUBROUTINE mpp_lnk_bdy_4d(cdname, ptab, cd_nat, psgn, kb_bdy)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, INTENT(IN) :: kb_bdy
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = SIZE(ptab, 3)
    ipl = SIZE(ptab, 4)
    ipf = 1
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ipk, ipl, ipf, ld_lbc = .TRUE.)
    ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2), zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), &
&zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    zland = 0._wp
    DO jf = 1, ipf
      IF (nbondi == 2) THEN
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(1, :, :, :) = ptab(jpim1, :, :, :)
          ptab(jpi, :, :, :) = ptab(2, :, :, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :, :, :) = zland
          ptab(nlci - nn_hls + 1 : jpi, :, :, :) = zland
          !$ACC END KERNELS
        END IF
      ELSE IF (nbondi == - 1) THEN
        !$ACC KERNELS
        IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :, :, :) = zland
        !$ACC END KERNELS
      ELSE IF (nbondi == 1) THEN
        ptab(nlci - nn_hls + 1 : jpi, :, :, :) = zland
      END IF
      IF (nbondj == 2) THEN
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(:, 1, :, :) = ptab(:, jpjm1, :, :)
          ptab(:, jpj, :, :) = ptab(:, 2, :, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls, :, :) = zland
          ptab(:, nlcj - nn_hls + 1 : jpj, :, :) = zland
          !$ACC END KERNELS
        END IF
      ELSE IF (nbondj == - 1) THEN
        !$ACC KERNELS
        IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls, :, :) = zland
        !$ACC END KERNELS
      ELSE IF (nbondj == 1) THEN
        ptab(:, nlcj - nn_hls + 1 : jpj, :, :) = zland
      END IF
    END DO
    DO jf = 1, ipf
      SELECT CASE (nbondi_bdy(kb_bdy))
      CASE (- 1, 0, 1)
        iihom = nlci - nreci
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(nn_hls + jh, :, jk, jl)
              zt3we(:, jh, jk, jl, jf, 1) = ptab(iihom + jh, :, jk, jl)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
      imigr = nn_hls * jpj * ipk * ipl
      IF (ln_timing) CALL tic_tac(.TRUE.)
      SELECT CASE (nbondi_bdy(kb_bdy))
      CASE (- 1)
        CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
      CASE (0)
        CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
        CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
      CASE (1)
        CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      END SELECT
      SELECT CASE (nbondi_bdy_b(kb_bdy))
      CASE (- 1)
        CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      CASE (0)
        CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
        CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      CASE (1)
        CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      END SELECT
      SELECT CASE (nbondi_bdy(kb_bdy))
      CASE (- 1)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE (0)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
        IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE (1)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      IF (ln_timing) CALL tic_tac(.FALSE.)
      iihom = nlci - nn_hls
      SELECT CASE (nbondi_bdy_b(kb_bdy))
      CASE (- 1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(iihom + jh, :, jk, jl) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (0)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :, jk, jl) = zt3we(:, jh, jk, jl, jf, 2)
              ptab(iihom + jh, :, jk, jl) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :, jk, jl) = zt3we(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
    END DO
    IF (npolj /= 0) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat, psgn)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat, psgn)
      END SELECT
    END IF
    DO jf = 1, ipf
      IF (nbondj_bdy(kb_bdy) /= 2) THEN
        ijhom = nlcj - nrecj
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(:, ijhom + jh, jk, jl)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(:, nn_hls + jh, jk, jl)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END IF
      imigr = nn_hls * jpi * ipk * ipl
      IF (ln_timing) CALL tic_tac(.TRUE.)
      SELECT CASE (nbondj_bdy(kb_bdy))
      CASE (- 1)
        CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
      CASE (0)
        CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
        CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
      CASE (1)
        CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      END SELECT
      SELECT CASE (nbondj_bdy_b(kb_bdy))
      CASE (- 1)
        CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      CASE (0)
        CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
        CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      CASE (1)
        CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      END SELECT
      SELECT CASE (nbondj_bdy(kb_bdy))
      CASE (- 1)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE (0)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
        IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE (1)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      IF (ln_timing) CALL tic_tac(.FALSE.)
      ijhom = nlcj - nn_hls
      SELECT CASE (nbondj_bdy_b(kb_bdy))
      CASE (- 1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, ijhom + jh, jk, jl) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (0)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh, jk, jl) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(:, ijhom + jh, jk, jl) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh, jk, jl) = zt3sn(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
    END DO
    DEALLOCATE(zt3ns, zt3sn, zt3ew, zt3we)
  END SUBROUTINE mpp_lnk_bdy_4d
  SUBROUTINE mppsend(ktyp, pmess, kbytes, kdest, md_req)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(INOUT) :: pmess(*)
    INTEGER, INTENT(IN) :: kbytes
    INTEGER, INTENT(IN) :: kdest
    INTEGER, INTENT(IN) :: ktyp
    INTEGER, INTENT(IN) :: md_req
    INTEGER :: iflag
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppsend', 'r0', 0, 0)
    SELECT CASE (cn_mpi_send)
    CASE ('S')
      CALL mpi_send(pmess, kbytes, mpi_double_precision, kdest, ktyp, mpi_comm_oce, iflag)
    CASE ('B')
      CALL mpi_bsend(pmess, kbytes, mpi_double_precision, kdest, ktyp, mpi_comm_oce, iflag)
    CASE ('I')
      CALL mpi_isend(pmess, kbytes, mpi_double_precision, kdest, ktyp, mpi_comm_oce, md_req, iflag)
    END SELECT
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE mppsend
  SUBROUTINE mpprecv(ktyp, pmess, kbytes, ksource)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(INOUT) :: pmess(*)
    INTEGER, INTENT(IN) :: kbytes
    INTEGER, INTENT(IN) :: ktyp
    INTEGER, OPTIONAL, INTENT(IN) :: ksource
    INTEGER :: istatus(mpi_status_size)
    INTEGER :: iflag
    INTEGER :: use_source
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mpprecv', 'r0', 0, 0)
    use_source = mpi_any_source
    IF (PRESENT(ksource)) use_source = ksource
    CALL mpi_recv(pmess, kbytes, mpi_double_precision, use_source, ktyp, mpi_comm_oce, istatus, iflag)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE mpprecv
  SUBROUTINE mppgather(ptab, kp, pio)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: ptab
    INTEGER, INTENT(IN) :: kp
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpnij), INTENT(OUT) :: pio
    INTEGER :: itaille, ierror
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppgather', 'r0', 0, 0)
    itaille = jpi * jpj
    CALL mpi_gather(ptab, itaille, mpi_double_precision, pio, itaille, mpi_double_precision, kp, mpi_comm_oce, ierror)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE mppgather
  SUBROUTINE mppscatter(pio, kp, ptab)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpnij) :: pio
    INTEGER :: kp
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ptab
    INTEGER :: itaille, ierror
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppscatter', 'r0', 0, 0)
    itaille = jpi * jpj
    CALL mpi_scatter(pio, itaille, mpi_double_precision, ptab, itaille, mpi_double_precision, kp, mpi_comm_oce, ierror)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE mppscatter
  SUBROUTINE mpp_delay_sum(cdname, cdelay, y_in, pout, ldlast, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    CHARACTER(LEN = *), INTENT(IN) :: cdelay
    COMPLEX(KIND = wp), INTENT(IN), DIMENSION(:) :: y_in
    REAL(KIND = wp), INTENT(OUT), DIMENSION(:) :: pout
    LOGICAL, INTENT(IN) :: ldlast
    INTEGER, INTENT(IN), OPTIONAL :: kcom
    INTEGER :: ji, isz
    INTEGER :: idvar
    INTEGER :: ierr, ilocalcomm
    COMPLEX(KIND = wp), ALLOCATABLE, DIMENSION(:) :: ytmp
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    isz = SIZE(y_in)
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_dlg = .TRUE.)
    idvar = - 1
    DO ji = 1, nbdelay
      IF (TRIM(cdelay) == TRIM(c_delaylist(ji))) idvar = ji
    END DO
    IF (idvar == - 1) CALL ctl_stop('STOP', ' mpp_delay_sum : please add a new delayed exchange for ' // TRIM(cdname))
    IF (ndelayid(idvar) == 0) THEN
      IF (SIZE(todelay(idvar) % z1d) /= isz) THEN
        IF (lwp) WRITE(numout, FMT = *) ' WARNING: the nb of delayed variables in restart file is not the model one'
        DEALLOCATE(todelay(idvar) % z1d)
        ndelayid(idvar) = - 1
      ELSE
        ALLOCATE(todelay(idvar) % y1d(isz))
        todelay(idvar) % y1d(:) = CMPLX(todelay(idvar) % z1d(:), 0., wp)
      END IF
    END IF
    IF (ndelayid(idvar) == - 1) THEN
      ALLOCATE(todelay(idvar) % z1d(isz), todelay(idvar) % y1d(isz))
      CALL mpi_allreduce(y_in(:), todelay(idvar) % y1d(:), isz, MPI_DOUBLE_COMPLEX, mpi_sumdd, ilocalcomm, ierr)
      todelay(idvar) % z1d(:) = REAL(todelay(idvar) % y1d(:), wp)
    END IF
    IF (ndelayid(idvar) > 0) CALL mpp_delay_rcv(idvar)
    pout(:) = todelay(idvar) % z1d(:)
    CALL mpi_iallreduce(y_in(:), todelay(idvar) % y1d(:), isz, MPI_DOUBLE_COMPLEX, mpi_sumdd, ilocalcomm, ndelayid(idvar), ierr)
  END SUBROUTINE mpp_delay_sum
  SUBROUTINE mpp_delay_max(cdname, cdelay, p_in, pout, ldlast, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    CHARACTER(LEN = *), INTENT(IN) :: cdelay
    REAL(KIND = wp), INTENT(IN), DIMENSION(:) :: p_in
    REAL(KIND = wp), INTENT(OUT), DIMENSION(:) :: pout
    LOGICAL, INTENT(IN) :: ldlast
    INTEGER, INTENT(IN), OPTIONAL :: kcom
    INTEGER :: ji, isz
    INTEGER :: idvar
    INTEGER :: ierr, ilocalcomm
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    isz = SIZE(p_in)
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_dlg = .TRUE.)
    idvar = - 1
    DO ji = 1, nbdelay
      IF (TRIM(cdelay) == TRIM(c_delaylist(ji))) idvar = ji
    END DO
    IF (idvar == - 1) CALL ctl_stop('STOP', ' mpp_delay_max : please add a new delayed exchange for ' // TRIM(cdname))
    IF (ndelayid(idvar) == 0) THEN
      IF (SIZE(todelay(idvar) % z1d) /= isz) THEN
        IF (lwp) WRITE(numout, FMT = *) ' WARNING: the nb of delayed variables in restart file is not the model one'
        DEALLOCATE(todelay(idvar) % z1d)
        ndelayid(idvar) = - 1
      END IF
    END IF
    IF (ndelayid(idvar) == - 1) THEN
      ALLOCATE(todelay(idvar) % z1d(isz))
      CALL mpi_allreduce(p_in(:), todelay(idvar) % z1d(:), isz, MPI_DOUBLE_PRECISION, mpi_max, ilocalcomm, ierr)
    END IF
    IF (ndelayid(idvar) > 0) CALL mpp_delay_rcv(idvar)
    pout(:) = todelay(idvar) % z1d(:)
    CALL mpi_iallreduce(p_in(:), todelay(idvar) % z1d(:), isz, MPI_DOUBLE_PRECISION, mpi_max, ilocalcomm, ndelayid(idvar), ierr)
  END SUBROUTINE mpp_delay_max
  SUBROUTINE mpp_delay_rcv(kid)
    INTEGER, INTENT(IN) :: kid
    INTEGER :: ierr
    IF (ndelayid(kid) /= - 2) THEN
      IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
      CALL mpi_wait(ndelayid(kid), MPI_STATUS_IGNORE, ierr)
      IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
      IF (ASSOCIATED(todelay(kid) % y1d)) todelay(kid) % z1d(:) = REAL(todelay(kid) % y1d(:), wp)
      ndelayid(kid) = - 2
    END IF
  END SUBROUTINE mpp_delay_rcv
  SUBROUTINE mppmax_int(cdname, ptab, kdim, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    INTEGER, INTENT(INOUT) :: ptab
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ipi, ii, ierr
    INTEGER :: ierror, ilocalcomm
    INTEGER, ALLOCATABLE :: work(:)
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmax_int', 'r0', 0, 0)
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    IF (PRESENT(kdim)) THEN
      ipi = kdim
    ELSE
      ipi = 1
    END IF
    ALLOCATE(work(ipi))
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL mpi_allreduce(ptab, work, ipi, mpi_integer, mpi_max, ilocalcomm, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    DO ii = 1, ipi
      ptab = work(ii)
    END DO
    !$ACC END KERNELS
    DEALLOCATE(work)
  END SUBROUTINE mppmax_int
  SUBROUTINE mppmax_a_int(cdname, ptab, kdim, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    INTEGER, INTENT(INOUT) :: ptab(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ipi, ii, ierr
    INTEGER :: ierror, ilocalcomm
    INTEGER, ALLOCATABLE :: work(:)
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmax_a_int', 'r0', 0, 0)
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    IF (PRESENT(kdim)) THEN
      ipi = kdim
    ELSE
      ipi = SIZE(ptab, 1)
    END IF
    ALLOCATE(work(ipi))
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL mpi_allreduce(ptab(:), work, ipi, mpi_integer, mpi_max, ilocalcomm, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    DO ii = 1, ipi
      ptab(ii) = work(ii)
    END DO
    !$ACC END KERNELS
    DEALLOCATE(work)
  END SUBROUTINE mppmax_a_int
  SUBROUTINE mppmax_real(cdname, ptab, kdim, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ipi, ii, ierr
    INTEGER :: ierror, ilocalcomm
    REAL(KIND = wp), ALLOCATABLE :: work(:)
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmax_real', 'r0', 0, 0)
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    IF (PRESENT(kdim)) THEN
      ipi = kdim
    ELSE
      ipi = 1
    END IF
    ALLOCATE(work(ipi))
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL mpi_allreduce(ptab, work, ipi, mpi_double_precision, mpi_max, ilocalcomm, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    DO ii = 1, ipi
      ptab = work(ii)
    END DO
    !$ACC END KERNELS
    DEALLOCATE(work)
  END SUBROUTINE mppmax_real
  SUBROUTINE mppmax_a_real(cdname, ptab, kdim, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ipi, ii, ierr
    INTEGER :: ierror, ilocalcomm
    REAL(KIND = wp), ALLOCATABLE :: work(:)
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmax_a_real', 'r0', 0, 0)
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    IF (PRESENT(kdim)) THEN
      ipi = kdim
    ELSE
      ipi = SIZE(ptab, 1)
    END IF
    ALLOCATE(work(ipi))
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL mpi_allreduce(ptab(:), work, ipi, mpi_double_precision, mpi_max, ilocalcomm, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    DO ii = 1, ipi
      ptab(ii) = work(ii)
    END DO
    !$ACC END KERNELS
    DEALLOCATE(work)
  END SUBROUTINE mppmax_a_real
  SUBROUTINE mppmin_int(cdname, ptab, kdim, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    INTEGER, INTENT(INOUT) :: ptab
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ipi, ii, ierr
    INTEGER :: ierror, ilocalcomm
    INTEGER, ALLOCATABLE :: work(:)
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmin_int', 'r0', 0, 0)
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    IF (PRESENT(kdim)) THEN
      ipi = kdim
    ELSE
      ipi = 1
    END IF
    ALLOCATE(work(ipi))
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL mpi_allreduce(ptab, work, ipi, mpi_integer, mpi_min, ilocalcomm, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    DO ii = 1, ipi
      ptab = work(ii)
    END DO
    !$ACC END KERNELS
    DEALLOCATE(work)
  END SUBROUTINE mppmin_int
  SUBROUTINE mppmin_a_int(cdname, ptab, kdim, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    INTEGER, INTENT(INOUT) :: ptab(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ipi, ii, ierr
    INTEGER :: ierror, ilocalcomm
    INTEGER, ALLOCATABLE :: work(:)
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmin_a_int', 'r0', 0, 0)
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    IF (PRESENT(kdim)) THEN
      ipi = kdim
    ELSE
      ipi = SIZE(ptab, 1)
    END IF
    ALLOCATE(work(ipi))
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL mpi_allreduce(ptab(:), work, ipi, mpi_integer, mpi_min, ilocalcomm, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    DO ii = 1, ipi
      ptab(ii) = work(ii)
    END DO
    !$ACC END KERNELS
    DEALLOCATE(work)
  END SUBROUTINE mppmin_a_int
  SUBROUTINE mppmin_real(cdname, ptab, kdim, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ipi, ii, ierr
    INTEGER :: ierror, ilocalcomm
    REAL(KIND = wp), ALLOCATABLE :: work(:)
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmin_real', 'r0', 0, 0)
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    IF (PRESENT(kdim)) THEN
      ipi = kdim
    ELSE
      ipi = 1
    END IF
    ALLOCATE(work(ipi))
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL mpi_allreduce(ptab, work, ipi, mpi_double_precision, mpi_min, ilocalcomm, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    DO ii = 1, ipi
      ptab = work(ii)
    END DO
    !$ACC END KERNELS
    DEALLOCATE(work)
  END SUBROUTINE mppmin_real
  SUBROUTINE mppmin_a_real(cdname, ptab, kdim, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ipi, ii, ierr
    INTEGER :: ierror, ilocalcomm
    REAL(KIND = wp), ALLOCATABLE :: work(:)
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmin_a_real', 'r0', 0, 0)
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    IF (PRESENT(kdim)) THEN
      ipi = kdim
    ELSE
      ipi = SIZE(ptab, 1)
    END IF
    ALLOCATE(work(ipi))
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL mpi_allreduce(ptab(:), work, ipi, mpi_double_precision, mpi_min, ilocalcomm, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    DO ii = 1, ipi
      ptab(ii) = work(ii)
    END DO
    !$ACC END KERNELS
    DEALLOCATE(work)
  END SUBROUTINE mppmin_a_real
  SUBROUTINE mppsum_int(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    INTEGER, INTENT(INOUT) :: ptab
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ipi, ii, ierr
    INTEGER :: ierror, ilocalcomm
    INTEGER, ALLOCATABLE :: work(:)
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    IF (PRESENT(kdim)) THEN
      ipi = kdim
    ELSE
      ipi = 1
    END IF
    ALLOCATE(work(ipi))
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL mpi_allreduce(ptab, work, ipi, mpi_integer, mpi_sum, ilocalcomm, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    !$ACC KERNELS
    DO ii = 1, ipi
      ptab = work(ii)
    END DO
    !$ACC END KERNELS
    DEALLOCATE(work)
  END SUBROUTINE mppsum_int
  SUBROUTINE mppsum_a_int(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    INTEGER, INTENT(INOUT) :: ptab(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ipi, ii, ierr
    INTEGER :: ierror, ilocalcomm
    INTEGER, ALLOCATABLE :: work(:)
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    IF (PRESENT(kdim)) THEN
      ipi = kdim
    ELSE
      ipi = SIZE(ptab, 1)
    END IF
    ALLOCATE(work(ipi))
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL mpi_allreduce(ptab(:), work, ipi, mpi_integer, mpi_sum, ilocalcomm, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    !$ACC KERNELS
    DO ii = 1, ipi
      ptab(ii) = work(ii)
    END DO
    !$ACC END KERNELS
    DEALLOCATE(work)
  END SUBROUTINE mppsum_a_int
  SUBROUTINE mppsum_real(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ipi, ii, ierr
    INTEGER :: ierror, ilocalcomm
    REAL(KIND = wp), ALLOCATABLE :: work(:)
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    IF (PRESENT(kdim)) THEN
      ipi = kdim
    ELSE
      ipi = 1
    END IF
    ALLOCATE(work(ipi))
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL mpi_allreduce(ptab, work, ipi, mpi_double_precision, mpi_sum, ilocalcomm, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    !$ACC KERNELS
    DO ii = 1, ipi
      ptab = work(ii)
    END DO
    !$ACC END KERNELS
    DEALLOCATE(work)
  END SUBROUTINE mppsum_real
  SUBROUTINE mppsum_a_real(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ipi, ii, ierr
    INTEGER :: ierror, ilocalcomm
    REAL(KIND = wp), ALLOCATABLE :: work(:)
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    IF (PRESENT(kdim)) THEN
      ipi = kdim
    ELSE
      ipi = SIZE(ptab, 1)
    END IF
    ALLOCATE(work(ipi))
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL mpi_allreduce(ptab(:), work, ipi, mpi_double_precision, mpi_sum, ilocalcomm, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    !$ACC KERNELS
    DO ii = 1, ipi
      ptab(ii) = work(ii)
    END DO
    !$ACC END KERNELS
    DEALLOCATE(work)
  END SUBROUTINE mppsum_a_real
  SUBROUTINE mppsum_realdd(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    COMPLEX, INTENT(INOUT) :: ptab
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ipi, ii, ierr
    INTEGER :: ierror, ilocalcomm
    COMPLEX, ALLOCATABLE :: work(:)
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    IF (PRESENT(kdim)) THEN
      ipi = kdim
    ELSE
      ipi = 1
    END IF
    ALLOCATE(work(ipi))
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL mpi_allreduce(ptab, work, ipi, mpi_double_complex, mpi_sumdd, ilocalcomm, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    !$ACC KERNELS
    DO ii = 1, ipi
      ptab = work(ii)
    END DO
    !$ACC END KERNELS
    DEALLOCATE(work)
  END SUBROUTINE mppsum_realdd
  SUBROUTINE mppsum_a_realdd(cdname, ptab, kdim, kcom)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    COMPLEX, INTENT(INOUT) :: ptab(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ipi, ii, ierr
    INTEGER :: ierror, ilocalcomm
    COMPLEX, ALLOCATABLE :: work(:)
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    IF (PRESENT(kdim)) THEN
      ipi = kdim
    ELSE
      ipi = SIZE(ptab, 1)
    END IF
    ALLOCATE(work(ipi))
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL mpi_allreduce(ptab(:), work, ipi, mpi_double_complex, mpi_sumdd, ilocalcomm, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    !$ACC KERNELS
    DO ii = 1, ipi
      ptab(ii) = work(ii)
    END DO
    !$ACC END KERNELS
    DEALLOCATE(work)
  END SUBROUTINE mppsum_a_realdd
  SUBROUTINE mpp_minloc2d(cdname, ptab, pmask, pmin, kindex)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :)
    REAL(KIND = wp), INTENT(IN) :: pmask(:, :)
    REAL(KIND = wp), INTENT(OUT) :: pmin
    INTEGER, INTENT(OUT) :: kindex(2)
    INTEGER :: ierror, ii, idim
    INTEGER :: index0
    REAL(KIND = wp) :: zmin
    INTEGER, DIMENSION(:), ALLOCATABLE :: ilocs
    REAL(KIND = wp), DIMENSION(2, 1) :: zain, zaout
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    idim = SIZE(kindex)
    IF (ALL(pmask(:, :) /= 1._wp)) THEN
      zmin = HUGE(zmin)
      index0 = 0
    ELSE
      ALLOCATE(ilocs(idim))
      ilocs = MINLOC(ptab(:, :), mask = pmask(:, :) == 1._wp)
      zmin = ptab(ilocs(1), ilocs(2))
      kindex(1) = mig(ilocs(1))
      kindex(2) = mjg(ilocs(2))
      DEALLOCATE(ilocs)
      index0 = kindex(1) - 1
      index0 = index0 + jpiglo * (kindex(2) - 1)
    END IF
    !$ACC KERNELS
    zain(1, :) = zmin
    zain(2, :) = REAL(index0, wp)
    !$ACC END KERNELS
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL MPI_ALLREDUCE(zain, zaout, 1, MPI_2DOUBLE_PRECISION, mpi_minloc, MPI_COMM_OCE, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    !$ACC KERNELS
    pmin = zaout(1, 1)
    index0 = NINT(zaout(2, 1))
    kindex(2) = index0 / jpiglo
    index0 = index0 - kindex(2) * jpiglo
    kindex(1) = index0
    kindex(:) = kindex(:) + 1
    !$ACC END KERNELS
  END SUBROUTINE mpp_minloc2d
  SUBROUTINE mpp_minloc3d(cdname, ptab, pmask, pmin, kindex)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :, :)
    REAL(KIND = wp), INTENT(IN) :: pmask(:, :, :)
    REAL(KIND = wp), INTENT(OUT) :: pmin
    INTEGER, INTENT(OUT) :: kindex(3)
    INTEGER :: ierror, ii, idim
    INTEGER :: index0
    REAL(KIND = wp) :: zmin
    INTEGER, DIMENSION(:), ALLOCATABLE :: ilocs
    REAL(KIND = wp), DIMENSION(2, 1) :: zain, zaout
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    idim = SIZE(kindex)
    IF (ALL(pmask(:, :, :) /= 1._wp)) THEN
      zmin = HUGE(zmin)
      index0 = 0
    ELSE
      ALLOCATE(ilocs(idim))
      ilocs = MINLOC(ptab(:, :, :), mask = pmask(:, :, :) == 1._wp)
      zmin = ptab(ilocs(1), ilocs(2), ilocs(3))
      kindex(1) = mig(ilocs(1))
      kindex(2) = mjg(ilocs(2))
      kindex(3) = ilocs(3)
      DEALLOCATE(ilocs)
      index0 = kindex(1) - 1
      index0 = index0 + jpiglo * (kindex(2) - 1)
      index0 = index0 + jpiglo * jpjglo * (kindex(3) - 1)
    END IF
    !$ACC KERNELS
    zain(1, :) = zmin
    zain(2, :) = REAL(index0, wp)
    !$ACC END KERNELS
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL MPI_ALLREDUCE(zain, zaout, 1, MPI_2DOUBLE_PRECISION, mpi_minloc, MPI_COMM_OCE, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    !$ACC KERNELS
    pmin = zaout(1, 1)
    index0 = NINT(zaout(2, 1))
    kindex(3) = index0 / (jpiglo * jpjglo)
    index0 = index0 - kindex(3) * (jpiglo * jpjglo)
    kindex(2) = index0 / jpiglo
    index0 = index0 - kindex(2) * jpiglo
    kindex(1) = index0
    kindex(:) = kindex(:) + 1
    !$ACC END KERNELS
  END SUBROUTINE mpp_minloc3d
  SUBROUTINE mpp_maxloc2d(cdname, ptab, pmask, pmin, kindex)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :)
    REAL(KIND = wp), INTENT(IN) :: pmask(:, :)
    REAL(KIND = wp), INTENT(OUT) :: pmin
    INTEGER, INTENT(OUT) :: kindex(2)
    INTEGER :: ierror, ii, idim
    INTEGER :: index0
    REAL(KIND = wp) :: zmin
    INTEGER, DIMENSION(:), ALLOCATABLE :: ilocs
    REAL(KIND = wp), DIMENSION(2, 1) :: zain, zaout
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    idim = SIZE(kindex)
    IF (ALL(pmask(:, :) /= 1._wp)) THEN
      zmin = - HUGE(zmin)
      index0 = 0
    ELSE
      ALLOCATE(ilocs(idim))
      ilocs = MAXLOC(ptab(:, :), mask = pmask(:, :) == 1._wp)
      zmin = ptab(ilocs(1), ilocs(2))
      kindex(1) = mig(ilocs(1))
      kindex(2) = mjg(ilocs(2))
      DEALLOCATE(ilocs)
      index0 = kindex(1) - 1
      index0 = index0 + jpiglo * (kindex(2) - 1)
    END IF
    !$ACC KERNELS
    zain(1, :) = zmin
    zain(2, :) = REAL(index0, wp)
    !$ACC END KERNELS
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL MPI_ALLREDUCE(zain, zaout, 1, MPI_2DOUBLE_PRECISION, mpi_maxloc, MPI_COMM_OCE, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    !$ACC KERNELS
    pmin = zaout(1, 1)
    index0 = NINT(zaout(2, 1))
    kindex(2) = index0 / jpiglo
    index0 = index0 - kindex(2) * jpiglo
    kindex(1) = index0
    kindex(:) = kindex(:) + 1
    !$ACC END KERNELS
  END SUBROUTINE mpp_maxloc2d
  SUBROUTINE mpp_maxloc3d(cdname, ptab, pmask, pmin, kindex)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :, :)
    REAL(KIND = wp), INTENT(IN) :: pmask(:, :, :)
    REAL(KIND = wp), INTENT(OUT) :: pmin
    INTEGER, INTENT(OUT) :: kindex(3)
    INTEGER :: ierror, ii, idim
    INTEGER :: index0
    REAL(KIND = wp) :: zmin
    INTEGER, DIMENSION(:), ALLOCATABLE :: ilocs
    REAL(KIND = wp), DIMENSION(2, 1) :: zain, zaout
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, ld_glb = .TRUE.)
    idim = SIZE(kindex)
    IF (ALL(pmask(:, :, :) /= 1._wp)) THEN
      zmin = - HUGE(zmin)
      index0 = 0
    ELSE
      ALLOCATE(ilocs(idim))
      ilocs = MAXLOC(ptab(:, :, :), mask = pmask(:, :, :) == 1._wp)
      zmin = ptab(ilocs(1), ilocs(2), ilocs(3))
      kindex(1) = mig(ilocs(1))
      kindex(2) = mjg(ilocs(2))
      kindex(3) = ilocs(3)
      DEALLOCATE(ilocs)
      index0 = kindex(1) - 1
      index0 = index0 + jpiglo * (kindex(2) - 1)
      index0 = index0 + jpiglo * jpjglo * (kindex(3) - 1)
    END IF
    !$ACC KERNELS
    zain(1, :) = zmin
    zain(2, :) = REAL(index0, wp)
    !$ACC END KERNELS
    IF (ln_timing) CALL tic_tac(.TRUE., ld_global = .TRUE.)
    CALL MPI_ALLREDUCE(zain, zaout, 1, MPI_2DOUBLE_PRECISION, mpi_maxloc, MPI_COMM_OCE, ierror)
    IF (ln_timing) CALL tic_tac(.FALSE., ld_global = .TRUE.)
    !$ACC KERNELS
    pmin = zaout(1, 1)
    index0 = NINT(zaout(2, 1))
    kindex(3) = index0 / (jpiglo * jpjglo)
    index0 = index0 - kindex(3) * (jpiglo * jpjglo)
    kindex(2) = index0 / jpiglo
    index0 = index0 - kindex(2) * jpiglo
    kindex(1) = index0
    kindex(:) = kindex(:) + 1
    !$ACC END KERNELS
  END SUBROUTINE mpp_maxloc3d
  SUBROUTINE mppsync
    INTEGER :: ierror
    CALL mpi_barrier(mpi_comm_oce, ierror)
  END SUBROUTINE mppsync
  SUBROUTINE mppstop(ldfinal, ld_force_abort)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    LOGICAL, OPTIONAL, INTENT(IN) :: ldfinal
    LOGICAL, OPTIONAL, INTENT(IN) :: ld_force_abort
    LOGICAL :: llfinal, ll_force_abort
    INTEGER :: info
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppstop', 'r0', 0, 0)
    llfinal = .FALSE.
    IF (PRESENT(ldfinal)) llfinal = ldfinal
    ll_force_abort = .FALSE.
    IF (PRESENT(ld_force_abort)) ll_force_abort = ld_force_abort
    IF (ll_force_abort) THEN
      CALL mpi_abort(MPI_COMM_WORLD)
    ELSE
      CALL mppsync
      CALL mpi_finalize(info)
    END IF
    IF (.NOT. llfinal) STOP 123
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE mppstop
  SUBROUTINE mpp_comm_free(kcom)
    INTEGER, INTENT(IN) :: kcom
    INTEGER :: ierr
    CALL MPI_COMM_FREE(kcom, ierr)
  END SUBROUTINE mpp_comm_free
  SUBROUTINE mpp_ini_znl(kumout)
    INTEGER, INTENT(IN) :: kumout
    INTEGER :: jproc
    INTEGER :: ierr, ii
    INTEGER, ALLOCATABLE, DIMENSION(:) :: kwork
    ALLOCATE(kwork(jpnij), STAT = ierr)
    IF (ierr /= 0) THEN
      WRITE(kumout, cform_err)
      WRITE(kumout, FMT = *) 'mpp_ini_znl : failed to allocate 1D array of length jpnij'
      CALL mppstop
    END IF
    IF (jpnj == 1) THEN
      ngrp_znl = ngrp_world
      ncomm_znl = mpi_comm_oce
    ELSE
      CALL MPI_ALLGATHER(njmpp, 1, mpi_integer, kwork, 1, mpi_integer, mpi_comm_oce, ierr)
      ndim_rank_znl = 0
      DO jproc = 1, jpnij
        IF (kwork(jproc) == njmpp) THEN
          ndim_rank_znl = ndim_rank_znl + 1
        END IF
      END DO
      IF (ALLOCATED(nrank_znl)) DEALLOCATE(nrank_znl)
      ALLOCATE(nrank_znl(ndim_rank_znl))
      !$ACC KERNELS
      ii = 0
      nrank_znl(:) = 0
      !$ACC END KERNELS
      DO jproc = 1, jpnij
        IF (kwork(jproc) == njmpp) THEN
          ii = ii + 1
          nrank_znl(ii) = jproc - 1
        END IF
      END DO
      CALL MPI_COMM_GROUP(mpi_comm_oce, ngrp_opa, ierr)
      CALL MPI_GROUP_INCL(ngrp_opa, ndim_rank_znl, nrank_znl, ngrp_znl, ierr)
      CALL MPI_COMM_CREATE(mpi_comm_oce, ngrp_znl, ncomm_znl, ierr)
    END IF
    IF (jpni == 1) THEN
      l_znl_root = .TRUE.
    ELSE
      l_znl_root = .FALSE.
      kwork(1) = nimpp
      CALL mpp_min('lib_mpp', kwork(1), kcom = ncomm_znl)
      IF (nimpp == kwork(1)) l_znl_root = .TRUE.
    END IF
    DEALLOCATE(kwork)
  END SUBROUTINE mpp_ini_znl
  SUBROUTINE mpp_ini_north
    INTEGER :: ierr
    INTEGER :: jjproc
    INTEGER :: ii, ji
    njmppmax = MAXVAL(njmppt)
    ndim_rank_north = 0
    DO jjproc = 1, jpnij
      IF (njmppt(jjproc) == njmppmax) ndim_rank_north = ndim_rank_north + 1
    END DO
    IF (ALLOCATED(nrank_north)) DEALLOCATE(nrank_north)
    ALLOCATE(nrank_north(ndim_rank_north))
    ii = 0
    DO ji = 1, jpnij
      IF (njmppt(ji) == njmppmax) THEN
        ii = ii + 1
        nrank_north(ii) = ji - 1
      END IF
    END DO
    CALL MPI_COMM_GROUP(mpi_comm_oce, ngrp_world, ierr)
    CALL MPI_GROUP_INCL(ngrp_world, ndim_rank_north, nrank_north, ngrp_north, ierr)
    CALL MPI_COMM_CREATE(mpi_comm_oce, ngrp_north, ncomm_north, ierr)
  END SUBROUTINE mpp_ini_north
  SUBROUTINE mpi_init_oce(ldtxt, ksft, code)
    CHARACTER(LEN = *), DIMENSION(:), INTENT(OUT) :: ldtxt
    INTEGER, INTENT(INOUT) :: ksft
    INTEGER, INTENT(OUT) :: code
    INTEGER :: ierr, ji
    LOGICAL :: mpi_was_called
    CALL mpi_initialized(mpi_was_called, code)
    IF (code /= MPI_SUCCESS) THEN
      DO ji = 1, SIZE(ldtxt)
        IF (TRIM(ldtxt(ji)) /= '') WRITE(*, FMT = *) ldtxt(ji)
      END DO
      WRITE(*, cform_err)
      WRITE(*, FMT = *) ' lib_mpp: Error in routine mpi_initialized'
      CALL mpi_abort(mpi_comm_world, code, ierr)
    END IF
    IF (.NOT. mpi_was_called) THEN
      CALL mpi_init(code)
      CALL mpi_comm_dup(mpi_comm_world, mpi_comm_oce, code)
      IF (code /= MPI_SUCCESS) THEN
        DO ji = 1, SIZE(ldtxt)
          IF (TRIM(ldtxt(ji)) /= '') WRITE(*, FMT = *) ldtxt(ji)
        END DO
        WRITE(*, cform_err)
        WRITE(*, FMT = *) ' lib_mpp: Error in routine mpi_comm_dup'
        CALL mpi_abort(mpi_comm_world, code, ierr)
      END IF
    END IF
    IF (nn_buffer > 0) THEN
      WRITE(ldtxt(ksft), FMT = *) 'mpi_bsend, buffer allocation of  : ', nn_buffer
      ksft = ksft + 1
      ALLOCATE(tampon(nn_buffer), STAT = ierr)
      IF (ierr /= 0) THEN
        DO ji = 1, SIZE(ldtxt)
          IF (TRIM(ldtxt(ji)) /= '') WRITE(*, FMT = *) ldtxt(ji)
        END DO
        WRITE(*, cform_err)
        WRITE(*, FMT = *) ' lib_mpp: Error in ALLOCATE', ierr
        CALL mpi_abort(mpi_comm_world, code, ierr)
      END IF
      CALL mpi_buffer_attach(tampon, nn_buffer, code)
    END IF
  END SUBROUTINE mpi_init_oce
  SUBROUTINE DDPDD_MPI(ydda, yddb, ilen, itype)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: ilen, itype
    COMPLEX(KIND = wp), DIMENSION(ilen), INTENT(IN) :: ydda
    COMPLEX(KIND = wp), DIMENSION(ilen), INTENT(INOUT) :: yddb
    REAL(KIND = wp) :: zerr, zt1, zt2
    INTEGER :: ji, ztmp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ddpdd_mpi', 'r0', 0, 0)
    ztmp = itype
    DO ji = 1, ilen
      zt1 = REAL(ydda(ji)) + REAL(yddb(ji))
      zerr = zt1 - REAL(ydda(ji))
      zt2 = ((REAL(yddb(ji)) - zerr) + (REAL(ydda(ji)) - (zt1 - zerr))) + AIMAG(ydda(ji)) + AIMAG(yddb(ji))
      yddb(ji) = CMPLX(zt1 + zt2, zt2 - ((zt1 + zt2) - zt1), wp)
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE DDPDD_MPI
  SUBROUTINE mpp_lbc_north_icb(pt2d, cd_type, psgn, kextj)
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: pt2d
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, INTENT(IN) :: kextj
    INTEGER :: ji, jj, jr
    INTEGER :: ierr, itaille, ildi, ilei, iilb
    INTEGER :: ipj, ij, iproc
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: ztab_e, znorthloc_e
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: znorthgloio_e
    ipj = 4
    ALLOCATE(ztab_e(jpiglo, 1 - kextj : ipj + kextj), znorthloc_e(jpimax, 1 - kextj : ipj + kextj), znorthgloio_e(jpimax, 1 - &
&kextj : ipj + kextj, jpni))
    !$ACC KERNELS
    ztab_e(:, :) = 0._wp
    znorthloc_e(:, :) = 0._wp
    ij = 1 - kextj
    DO jj = jpj - ipj + 1 - kextj, jpj + kextj
      znorthloc_e(1 : jpi, ij) = pt2d(1 : jpi, jj)
      ij = ij + 1
    END DO
    itaille = jpimax * (ipj + 2 * kextj)
    !$ACC END KERNELS
    IF (ln_timing) CALL tic_tac(.TRUE.)
    CALL MPI_ALLGATHER(znorthloc_e(1, 1 - kextj), itaille, MPI_DOUBLE_PRECISION, znorthgloio_e(1, 1 - kextj, 1), itaille, &
&MPI_DOUBLE_PRECISION, ncomm_north, ierr)
    IF (ln_timing) CALL tic_tac(.FALSE.)
    !$ACC KERNELS
    DO jr = 1, ndim_rank_north
      iproc = nrank_north(jr) + 1
      ildi = nldit(iproc)
      ilei = nleit(iproc)
      iilb = nimppt(iproc)
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1 - kextj, ipj + kextj
        DO ji = ildi, ilei
          ztab_e(ji + iilb - 1, jj) = znorthgloio_e(ji, jj, jr)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_nfd(ztab_e(:, 1 - kextj : ipj + kextj), cd_type, psgn, kextj)
    !$ACC KERNELS
    ij = 1 - kextj
    DO jj = jpj - ipj + 1 - kextj, jpj + kextj
      DO ji = 1, jpi
        pt2d(ji, jj) = ztab_e(ji + nimpp - 1, ij)
      END DO
      ij = ij + 1
    END DO
    !$ACC END KERNELS
    DEALLOCATE(ztab_e, znorthloc_e, znorthgloio_e)
  END SUBROUTINE mpp_lbc_north_icb
  SUBROUTINE mpp_lnk_2d_icb(cdname, pt2d, cd_type, psgn, kexti, kextj)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), DIMENSION(1 - kexti : jpi + kexti, 1 - kextj : jpj + kextj), INTENT(INOUT) :: pt2d
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, INTENT(IN) :: kexti
    INTEGER, INTENT(IN) :: kextj
    INTEGER :: jl
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ipreci, iprecj
    INTEGER :: ml_req1, ml_req2, ml_err
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(1 - kexti : jpi + kexti, nn_hls + kextj, 2) :: r2dns, r2dsn
    REAL(KIND = wp), DIMENSION(1 - kextj : jpj + kextj, nn_hls + kexti, 2) :: r2dwe, r2dew
    ipreci = nn_hls + kexti
    iprecj = nn_hls + kextj
    IF (narea == 1 .AND. numcom == - 1) CALL mpp_report(cdname, 1, 1, 1, ld_lbc = .TRUE.)
    IF (l_Iperio) THEN
      !$ACC KERNELS
      pt2d(1 - kexti : 1, :) = pt2d(jpim1 - kexti : jpim1, :)
      pt2d(jpi : jpi + kexti, :) = pt2d(2 : 2 + kexti, :)
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      IF (.NOT. cd_type == 'F') pt2d(1 - kexti : nn_hls, :) = 0._wp
      pt2d(jpi - nn_hls + 1 : jpi + kexti, :) = 0._wp
      !$ACC END KERNELS
    END IF
    IF (l_Jperio) THEN
      !$ACC KERNELS
      pt2d(:, 1 - kextj : 1) = pt2d(:, jpjm1 - kextj : jpjm1)
      pt2d(:, jpj : jpj + kextj) = pt2d(:, 2 : 2 + kextj)
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      IF (.NOT. cd_type == 'F') pt2d(:, 1 - kextj : nn_hls) = 0._wp
      pt2d(:, jpj - nn_hls + 1 : jpj + kextj) = 0._wp
      !$ACC END KERNELS
    END IF
    IF (npolj /= 0) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(pt2d(1 : jpi, 1 : jpj + kextj), cd_type, psgn, kextj)
      CASE DEFAULT
        CALL mpp_lbc_north_icb(pt2d(1 : jpi, 1 : jpj + kextj), cd_type, psgn, kextj)
      END SELECT
    END IF
    !$ACC KERNELS
    SELECT CASE (nbondi)
    CASE (- 1, 0, 1)
      iihom = jpi - nreci - kexti
      DO jl = 1, ipreci
        r2dew(:, jl, 1) = pt2d(nn_hls + jl, :)
        r2dwe(:, jl, 1) = pt2d(iihom + jl, :)
      END DO
    END SELECT
    imigr = ipreci * (jpj + 2 * kextj)
    !$ACC END KERNELS
    IF (ln_timing) CALL tic_tac(.TRUE.)
    SELECT CASE (nbondi)
    CASE (- 1)
      CALL mppsend(2, r2dwe(1 - kextj, 1, 1), imigr, noea, ml_req1)
      CALL mpprecv(1, r2dew(1 - kextj, 1, 2), imigr, noea)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(1, r2dew(1 - kextj, 1, 1), imigr, nowe, ml_req1)
      CALL mppsend(2, r2dwe(1 - kextj, 1, 1), imigr, noea, ml_req2)
      CALL mpprecv(1, r2dew(1 - kextj, 1, 2), imigr, noea)
      CALL mpprecv(2, r2dwe(1 - kextj, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(1, r2dew(1 - kextj, 1, 1), imigr, nowe, ml_req1)
      CALL mpprecv(2, r2dwe(1 - kextj, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    IF (ln_timing) CALL tic_tac(.FALSE.)
    !$ACC KERNELS
    iihom = jpi - nn_hls
    SELECT CASE (nbondi)
    CASE (- 1)
      DO jl = 1, ipreci
        pt2d(iihom + jl, :) = r2dew(:, jl, 2)
      END DO
    CASE (0)
      DO jl = 1, ipreci
        pt2d(jl - kexti, :) = r2dwe(:, jl, 2)
        pt2d(iihom + jl, :) = r2dew(:, jl, 2)
      END DO
    CASE (1)
      DO jl = 1, ipreci
        pt2d(jl - kexti, :) = r2dwe(:, jl, 2)
      END DO
    END SELECT
    IF (nbondj /= 2) THEN
      ijhom = jpj - nrecj - kextj
      DO jl = 1, iprecj
        r2dsn(:, jl, 1) = pt2d(:, ijhom + jl)
        r2dns(:, jl, 1) = pt2d(:, nn_hls + jl)
      END DO
    END IF
    imigr = iprecj * (jpi + 2 * kexti)
    !$ACC END KERNELS
    IF (ln_timing) CALL tic_tac(.TRUE.)
    SELECT CASE (nbondj)
    CASE (- 1)
      CALL mppsend(4, r2dsn(1 - kexti, 1, 1), imigr, nono, ml_req1)
      CALL mpprecv(3, r2dns(1 - kexti, 1, 2), imigr, nono)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(3, r2dns(1 - kexti, 1, 1), imigr, noso, ml_req1)
      CALL mppsend(4, r2dsn(1 - kexti, 1, 1), imigr, nono, ml_req2)
      CALL mpprecv(3, r2dns(1 - kexti, 1, 2), imigr, nono)
      CALL mpprecv(4, r2dsn(1 - kexti, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(3, r2dns(1 - kexti, 1, 1), imigr, noso, ml_req1)
      CALL mpprecv(4, r2dsn(1 - kexti, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    IF (ln_timing) CALL tic_tac(.FALSE.)
    !$ACC KERNELS
    ijhom = jpj - nn_hls
    SELECT CASE (nbondj)
    CASE (- 1)
      DO jl = 1, iprecj
        pt2d(:, ijhom + jl) = r2dns(:, jl, 2)
      END DO
    CASE (0)
      DO jl = 1, iprecj
        pt2d(:, jl - kextj) = r2dsn(:, jl, 2)
        pt2d(:, ijhom + jl) = r2dns(:, jl, 2)
      END DO
    CASE (1)
      DO jl = 1, iprecj
        pt2d(:, jl - kextj) = r2dsn(:, jl, 2)
      END DO
    END SELECT
    !$ACC END KERNELS
  END SUBROUTINE mpp_lnk_2d_icb
  SUBROUTINE mpp_report(cdname, kpk, kpl, kpf, ld_lbc, ld_glb, ld_dlg)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    INTEGER, OPTIONAL, INTENT(IN) :: kpk, kpl, kpf
    LOGICAL, OPTIONAL, INTENT(IN) :: ld_lbc, ld_glb, ld_dlg
    CHARACTER(LEN = 128) :: ccountname
    LOGICAL :: ll_lbc, ll_glb, ll_dlg
    INTEGER :: ji, jj, jk, jh, jf, jcount
    ll_lbc = .FALSE.
    IF (PRESENT(ld_lbc)) ll_lbc = ld_lbc
    ll_glb = .FALSE.
    IF (PRESENT(ld_glb)) ll_glb = ld_glb
    ll_dlg = .FALSE.
    IF (PRESENT(ld_dlg)) ll_dlg = ld_dlg
    IF (ncom_dttrc /= 1) CALL ctl_stop('STOP', 'mpp_report, ncom_dttrc /= 1 not coded...')
    ncom_freq = ncom_fsbc
    IF (ncom_stp == nit000 + ncom_freq) THEN
      IF (ll_lbc) THEN
        IF (.NOT. ALLOCATED(ncomm_sequence)) ALLOCATE(ncomm_sequence(ncom_rec_max, 2))
        IF (.NOT. ALLOCATED(crname_lbc)) ALLOCATE(crname_lbc(ncom_rec_max))
        n_sequence_lbc = n_sequence_lbc + 1
        IF (n_sequence_lbc > ncom_rec_max) CALL ctl_stop('STOP', 'lib_mpp, increase ncom_rec_max')
        crname_lbc(n_sequence_lbc) = cdname
        ncomm_sequence(n_sequence_lbc, 1) = kpk * kpl
        ncomm_sequence(n_sequence_lbc, 2) = kpf
      END IF
      IF (ll_glb) THEN
        IF (.NOT. ALLOCATED(crname_glb)) ALLOCATE(crname_glb(ncom_rec_max))
        n_sequence_glb = n_sequence_glb + 1
        IF (n_sequence_glb > ncom_rec_max) CALL ctl_stop('STOP', 'lib_mpp, increase ncom_rec_max')
        crname_glb(n_sequence_glb) = cdname
      END IF
      IF (ll_dlg) THEN
        IF (.NOT. ALLOCATED(crname_dlg)) ALLOCATE(crname_dlg(ncom_rec_max))
        n_sequence_dlg = n_sequence_dlg + 1
        IF (n_sequence_dlg > ncom_rec_max) CALL ctl_stop('STOP', 'lib_mpp, increase ncom_rec_max')
        crname_dlg(n_sequence_dlg) = cdname
      END IF
    ELSE IF (ncom_stp == nit000 + 2 * ncom_freq) THEN
      CALL ctl_opn(numcom, 'communication_report.txt', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, numout, .FALSE., narea)
      WRITE(numcom, FMT = *) ' '
      WRITE(numcom, FMT = *) ' ------------------------------------------------------------'
      WRITE(numcom, FMT = *) ' Communication pattern report (second oce+sbc+top time step):'
      WRITE(numcom, FMT = *) ' ------------------------------------------------------------'
      WRITE(numcom, FMT = *) ' '
      WRITE(numcom, FMT = '(A,I4)') ' Exchanged halos : ', n_sequence_lbc
      !$ACC KERNELS
      jj = 0
      jk = 0
      jf = 0
      jh = 0
      DO ji = 1, n_sequence_lbc
        IF (ncomm_sequence(ji, 1) .GT. 1) jk = jk + 1
        IF (ncomm_sequence(ji, 2) .GT. 1) jf = jf + 1
        IF (ncomm_sequence(ji, 1) .GT. 1 .AND. ncomm_sequence(ji, 2) .GT. 1) jj = jj + 1
        jh = MAX(jh, ncomm_sequence(ji, 1) * ncomm_sequence(ji, 2))
      END DO
      !$ACC END KERNELS
      WRITE(numcom, FMT = '(A,I3)') ' 3D Exchanged halos : ', jk
      WRITE(numcom, FMT = '(A,I3)') ' Multi arrays exchanged halos : ', jf
      WRITE(numcom, FMT = '(A,I3)') '   from which 3D : ', jj
      WRITE(numcom, FMT = '(A,I10)') ' Array max size : ', jh * jpi * jpj
      WRITE(numcom, FMT = *) ' '
      WRITE(numcom, FMT = *) ' lbc_lnk called'
      DO ji = 1, n_sequence_lbc - 1
        IF (crname_lbc(ji) /= 'already counted') THEN
          ccountname = crname_lbc(ji)
          crname_lbc(ji) = 'already counted'
          jcount = 1
          DO jj = ji + 1, n_sequence_lbc
            IF (ccountname == crname_lbc(jj)) THEN
              jcount = jcount + 1
              crname_lbc(jj) = 'already counted'
            END IF
          END DO
          WRITE(numcom, FMT = '(A, I4, A, A)') ' - ', jcount, ' times by subroutine ', TRIM(ccountname)
        END IF
      END DO
      IF (crname_lbc(n_sequence_lbc) /= 'already counted') THEN
        WRITE(numcom, FMT = '(A, I4, A, A)') ' - ', 1, ' times by subroutine ', TRIM(crname_lbc(ncom_rec_max))
      END IF
      WRITE(numcom, FMT = *) ' '
      IF (n_sequence_glb > 0) THEN
        WRITE(numcom, FMT = '(A,I4)') ' Global communications : ', n_sequence_glb
        jj = 1
        DO ji = 2, n_sequence_glb
          IF (crname_glb(ji - 1) /= crname_glb(ji)) THEN
            WRITE(numcom, FMT = '(A, I4, A, A)') ' - ', jj, ' times by subroutine ', TRIM(crname_glb(ji - 1))
            jj = 0
          END IF
          jj = jj + 1
        END DO
        WRITE(numcom, FMT = '(A, I4, A, A)') ' - ', jj, ' times by subroutine ', TRIM(crname_glb(n_sequence_glb))
        DEALLOCATE(crname_glb)
      ELSE
        WRITE(numcom, FMT = *) ' No MPI global communication '
      END IF
      WRITE(numcom, FMT = *) ' '
      IF (n_sequence_dlg > 0) THEN
        WRITE(numcom, FMT = '(A,I4)') ' Delayed global communications : ', n_sequence_dlg
        jj = 1
        DO ji = 2, n_sequence_dlg
          IF (crname_dlg(ji - 1) /= crname_dlg(ji)) THEN
            WRITE(numcom, FMT = '(A, I4, A, A)') ' - ', jj, ' times by subroutine ', TRIM(crname_dlg(ji - 1))
            jj = 0
          END IF
          jj = jj + 1
        END DO
        WRITE(numcom, FMT = '(A, I4, A, A)') ' - ', jj, ' times by subroutine ', TRIM(crname_dlg(n_sequence_dlg))
        DEALLOCATE(crname_dlg)
      ELSE
        WRITE(numcom, FMT = *) ' No MPI delayed global communication '
      END IF
      WRITE(numcom, FMT = *) ' '
      WRITE(numcom, FMT = *) ' -----------------------------------------------'
      WRITE(numcom, FMT = *) ' '
      DEALLOCATE(ncomm_sequence)
      DEALLOCATE(crname_lbc)
    END IF
  END SUBROUTINE mpp_report
  SUBROUTINE tic_tac(ld_tic, ld_global)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    LOGICAL, INTENT(IN) :: ld_tic
    LOGICAL, OPTIONAL, INTENT(IN) :: ld_global
    REAL(KIND = wp), DIMENSION(2), SAVE :: tic_wt
    REAL(KIND = wp), SAVE :: tic_ct = 0._wp
    INTEGER :: ii
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    IF (ncom_stp <= nit000) RETURN
    IF (ncom_stp == nitend) RETURN
    CALL profile_psy_data0 % PreStart('tic_tac', 'r0', 0, 0)
    ii = 1
    IF (PRESENT(ld_global)) THEN
      IF (ld_global) ii = 2
    END IF
    IF (ld_tic) THEN
      tic_wt(ii) = MPI_Wtime()
      IF (tic_ct > 0.0_wp) compute_time = compute_time + MPI_Wtime() - tic_ct
    ELSE
      waiting_time(ii) = waiting_time(ii) + MPI_Wtime() - tic_wt(ii)
      tic_ct = MPI_Wtime()
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE tic_tac
  SUBROUTINE ctl_stop(cd1, cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cd1, cd2, cd3, cd4, cd5
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cd6, cd7, cd8, cd9, cd10
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ctl_stop', 'r0', 0, 0)
    nstop = nstop + 1
    IF (numout == 6) CALL ctl_opn(numout, 'ocean.output', 'APPEND', 'FORMATTED', 'SEQUENTIAL', - 1, 6, .FALSE.)
    WRITE(numout, cform_err)
    IF (PRESENT(cd1)) WRITE(numout, FMT = *) TRIM(cd1)
    IF (PRESENT(cd2)) WRITE(numout, FMT = *) TRIM(cd2)
    IF (PRESENT(cd3)) WRITE(numout, FMT = *) TRIM(cd3)
    IF (PRESENT(cd4)) WRITE(numout, FMT = *) TRIM(cd4)
    IF (PRESENT(cd5)) WRITE(numout, FMT = *) TRIM(cd5)
    IF (PRESENT(cd6)) WRITE(numout, FMT = *) TRIM(cd6)
    IF (PRESENT(cd7)) WRITE(numout, FMT = *) TRIM(cd7)
    IF (PRESENT(cd8)) WRITE(numout, FMT = *) TRIM(cd8)
    IF (PRESENT(cd9)) WRITE(numout, FMT = *) TRIM(cd9)
    IF (PRESENT(cd10)) WRITE(numout, FMT = *) TRIM(cd10)
    CALL FLUSH(numout)
    IF (numstp /= - 1) CALL FLUSH(numstp)
    IF (numrun /= - 1) CALL FLUSH(numrun)
    IF (numevo_ice /= - 1) CALL FLUSH(numevo_ice)
    IF (cd1 == 'STOP') THEN
      WRITE(numout, FMT = *) 'huge E-R-R-O-R : immediate stop'
      CALL mppstop(ld_force_abort = .TRUE.)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ctl_stop
  SUBROUTINE ctl_warn(cd1, cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cd1, cd2, cd3, cd4, cd5
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cd6, cd7, cd8, cd9, cd10
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ctl_warn', 'r0', 0, 0)
    nwarn = nwarn + 1
    IF (lwp) THEN
      WRITE(numout, cform_war)
      IF (PRESENT(cd1)) WRITE(numout, FMT = *) TRIM(cd1)
      IF (PRESENT(cd2)) WRITE(numout, FMT = *) TRIM(cd2)
      IF (PRESENT(cd3)) WRITE(numout, FMT = *) TRIM(cd3)
      IF (PRESENT(cd4)) WRITE(numout, FMT = *) TRIM(cd4)
      IF (PRESENT(cd5)) WRITE(numout, FMT = *) TRIM(cd5)
      IF (PRESENT(cd6)) WRITE(numout, FMT = *) TRIM(cd6)
      IF (PRESENT(cd7)) WRITE(numout, FMT = *) TRIM(cd7)
      IF (PRESENT(cd8)) WRITE(numout, FMT = *) TRIM(cd8)
      IF (PRESENT(cd9)) WRITE(numout, FMT = *) TRIM(cd9)
      IF (PRESENT(cd10)) WRITE(numout, FMT = *) TRIM(cd10)
    END IF
    CALL FLUSH(numout)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ctl_warn
  SUBROUTINE ctl_opn(knum, cdfile, cdstat, cdform, cdacce, klengh, kout, ldwp, karea)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(OUT) :: knum
    CHARACTER(LEN = *), INTENT(IN) :: cdfile
    CHARACTER(LEN = *), INTENT(IN) :: cdstat
    CHARACTER(LEN = *), INTENT(IN) :: cdform
    CHARACTER(LEN = *), INTENT(IN) :: cdacce
    INTEGER, INTENT(IN) :: klengh
    INTEGER, INTENT(IN) :: kout
    LOGICAL, INTENT(IN) :: ldwp
    INTEGER, OPTIONAL, INTENT(IN) :: karea
    CHARACTER(LEN = 80) :: clfile
    INTEGER :: iost
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ctl_opn', 'r0', 0, 0)
    clfile = TRIM(cdfile)
    IF (PRESENT(karea)) THEN
      IF (karea > 1) WRITE(clfile, FMT = "(a,'_',i6.6)") TRIM(clfile), karea - 1
    END IF
    knum = get_unit()
    IF (TRIM(cdfile) == '/dev/null') clfile = TRIM(cdfile)
    iost = 0
    IF (cdacce(1 : 6) == 'DIRECT') THEN
      OPEN(UNIT = knum, FILE = clfile, FORM = cdform, ACCESS = cdacce, STATUS = cdstat, RECL = klengh, ERR = 100, IOSTAT = iost)
    ELSE IF (TRIM(cdstat) == 'APPEND') THEN
      OPEN(UNIT = knum, FILE = clfile, FORM = cdform, ACCESS = cdacce, STATUS = 'UNKNOWN', POSITION = 'APPEND', ERR = 100, IOSTAT &
&= iost)
    ELSE
      OPEN(UNIT = knum, FILE = clfile, FORM = cdform, ACCESS = cdacce, STATUS = cdstat, ERR = 100, IOSTAT = iost)
    END IF
    IF (iost /= 0 .AND. TRIM(clfile) == '/dev/null') OPEN(UNIT = knum, FILE = 'NUL', FORM = cdform, ACCESS = cdacce, STATUS = &
&cdstat, ERR = 100, IOSTAT = iost)
    IF (iost == 0) THEN
      IF (ldwp) THEN
        WRITE(kout, FMT = *) '     file   : ', TRIM(clfile), ' open ok'
        WRITE(kout, FMT = *) '     unit   = ', knum
        WRITE(kout, FMT = *) '     status = ', cdstat
        WRITE(kout, FMT = *) '     form   = ', cdform
        WRITE(kout, FMT = *) '     access = ', cdacce
        WRITE(kout, FMT = *)
      END IF
    END IF
100 CONTINUE
    IF (iost /= 0) THEN
      IF (ldwp) THEN
        WRITE(kout, FMT = *)
        WRITE(kout, FMT = *) ' ===>>>> : bad opening file: ', TRIM(clfile)
        WRITE(kout, FMT = *) ' =======   ===  '
        WRITE(kout, FMT = *) '           unit   = ', knum
        WRITE(kout, FMT = *) '           status = ', cdstat
        WRITE(kout, FMT = *) '           form   = ', cdform
        WRITE(kout, FMT = *) '           access = ', cdacce
        WRITE(kout, FMT = *) '           iostat = ', iost
        WRITE(kout, FMT = *) '           we stop. verify the file '
        WRITE(kout, FMT = *)
      ELSE
        WRITE(*, FMT = *)
        WRITE(*, FMT = *) ' ===>>>> : bad opening file: ', TRIM(clfile)
        WRITE(*, FMT = *) ' =======   ===  '
        WRITE(*, FMT = *) '           unit   = ', knum
        WRITE(*, FMT = *) '           status = ', cdstat
        WRITE(*, FMT = *) '           form   = ', cdform
        WRITE(*, FMT = *) '           access = ', cdacce
        WRITE(*, FMT = *) '           iostat = ', iost
        WRITE(*, FMT = *) '           we stop. verify the file '
        WRITE(*, FMT = *)
      END IF
      CALL FLUSH(kout)
      STOP 'ctl_opn bad opening'
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ctl_opn
  SUBROUTINE ctl_nam(kios, cdnam, ldwp)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(INOUT) :: kios
    CHARACTER(LEN = *), INTENT(IN) :: cdnam
    CHARACTER(LEN = 5) :: clios
    LOGICAL, INTENT(IN) :: ldwp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ctl_nam', 'r0', 0, 0)
    WRITE(clios, FMT = '(I5.0)') kios
    IF (kios < 0) THEN
      CALL ctl_warn('end of record or file while reading namelist ' // TRIM(cdnam) // ' iostat = ' // TRIM(clios))
    END IF
    IF (kios > 0) THEN
      CALL ctl_stop('misspelled variable in namelist ' // TRIM(cdnam) // ' iostat = ' // TRIM(clios))
    END IF
    kios = 0
    CALL profile_psy_data0 % PostEnd
    RETURN
  END SUBROUTINE ctl_nam
  INTEGER FUNCTION get_unit()
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    LOGICAL :: llopn
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('get_unit', 'r0', 0, 0)
    get_unit = 15
    llopn = .TRUE.
    DO WHILE ((get_unit < 998) .AND. llopn)
      get_unit = get_unit + 1
      INQUIRE(UNIT = get_unit, OPENED = llopn)
    END DO
    IF ((get_unit == 999) .AND. llopn) THEN
      CALL ctl_stop('get_unit: All logical units until 999 are used...')
      get_unit = - 1
    END IF
    CALL profile_psy_data0 % PostEnd
  END FUNCTION get_unit
END MODULE lib_mpp
