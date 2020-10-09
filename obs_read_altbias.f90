MODULE obs_read_altbias
  USE par_kind, ONLY: wp, dp, sp
  USE par_oce, ONLY: jpi, jpj, jpim1
  USE in_out_manager, ONLY: lwp, numout
  USE obs_surf_def
  USE dom_oce, ONLY: tmask, tmask_i, e1t, e2t, gphit
  USE oce, ONLY: sshn
  USE obs_inter_h2d
  USE obs_utils
  USE obs_inter_sup
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: obs_rea_altbias
  CONTAINS
  SUBROUTINE obs_rea_altbias(sladata, k2dint, bias_file)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE iom
    TYPE(obs_surf), INTENT(INOUT) :: sladata
    INTEGER, INTENT(IN) :: k2dint
    CHARACTER(LEN = 128) :: bias_file
    CHARACTER(LEN = 12), PARAMETER :: cpname = 'obs_rea_altbias'
    INTEGER :: jobs
    INTEGER :: jpialtbias
    INTEGER :: jpjaltbias
    INTEGER :: iico
    INTEGER :: ijco
    INTEGER :: i_nx_id
    INTEGER :: i_ny_id
    INTEGER :: i_file_id
    INTEGER :: i_var_id
    REAL(KIND = wp), DIMENSION(1) :: zext, zobsmask
    REAL(KIND = wp), DIMENSION(2, 2, 1) :: zweig
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: zmask, zbias, zglam, zgphi
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z_altbias
    REAL(KIND = wp) :: zlam
    REAL(KIND = wp) :: zphi
    INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: igrdi, igrdj
    INTEGER :: numaltbias
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('obs_rea_altbias', 'r0', 0, 0)
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) ' obs_rea_altbias : '
    IF (lwp) WRITE(numout, FMT = *) ' ------------- '
    IF (lwp) WRITE(numout, FMT = *) '   Read altimeter bias'
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    z_altbias(:, :) = 0.0_wp
    numaltbias = 0
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('obs_rea_altbias', 'r1', 0, 0)
    IF (lwp) WRITE(numout, FMT = *) 'Opening ', bias_file
    CALL iom_open(bias_file, numaltbias, ldstop = .FALSE.)
    IF (numaltbias .GT. 0) THEN
      CALL iom_get(numaltbias, jpdom_data, 'altbias', z_altbias(:, :), 1)
      CALL iom_close(numaltbias)
    ELSE
      IF (lwp) WRITE(numout, FMT = *) 'no file found'
    END IF
    ALLOCATE(igrdi(2, 2, sladata % nsurf), igrdj(2, 2, sladata % nsurf), zglam(2, 2, sladata % nsurf), zgphi(2, 2, sladata % &
&nsurf), zmask(2, 2, sladata % nsurf), zbias(2, 2, sladata % nsurf))
    DO jobs = 1, sladata % nsurf
      igrdi(1, 1, jobs) = sladata % mi(jobs) - 1
      igrdj(1, 1, jobs) = sladata % mj(jobs) - 1
      igrdi(1, 2, jobs) = sladata % mi(jobs) - 1
      igrdj(1, 2, jobs) = sladata % mj(jobs)
      igrdi(2, 1, jobs) = sladata % mi(jobs)
      igrdj(2, 1, jobs) = sladata % mj(jobs) - 1
      igrdi(2, 2, jobs) = sladata % mi(jobs)
      igrdj(2, 2, jobs) = sladata % mj(jobs)
    END DO
    CALL obs_int_comm_2d(2, 2, sladata % nsurf, jpi, jpj, igrdi, igrdj, glamt, zglam)
    CALL obs_int_comm_2d(2, 2, sladata % nsurf, jpi, jpj, igrdi, igrdj, gphit, zgphi)
    CALL obs_int_comm_2d(2, 2, sladata % nsurf, jpi, jpj, igrdi, igrdj, tmask(:, :, 1), zmask)
    CALL obs_int_comm_2d(2, 2, sladata % nsurf, jpi, jpj, igrdi, igrdj, z_altbias, zbias)
    DO jobs = 1, sladata % nsurf
      zlam = sladata % rlam(jobs)
      zphi = sladata % rphi(jobs)
      iico = sladata % mi(jobs)
      ijco = sladata % mj(jobs)
      CALL obs_int_h2d_init(1, 1, k2dint, zlam, zphi, zglam(:, :, jobs), zgphi(:, :, jobs), zmask(:, :, jobs), zweig, zobsmask)
      CALL obs_int_h2d(1, 1, zweig, zbias(:, :, jobs), zext)
      sladata % rext(jobs, 2) = sladata % rext(jobs, 2) - zext(1)
    END DO
    DEALLOCATE(igrdi, igrdj, zglam, zgphi, zmask, zbias)
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE obs_rea_altbias
END MODULE obs_read_altbias