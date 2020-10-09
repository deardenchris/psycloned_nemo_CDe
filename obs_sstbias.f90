MODULE obs_sstbias
  USE par_kind, ONLY: wp, dp, sp
  USE par_oce, ONLY: jpi, jpj, jpim1
  USE in_out_manager, ONLY: lwp, numout
  USE obs_surf_def
  USE dom_oce, ONLY: tmask, tmask_i, e1t, e2t, gphit, glamt
  USE oce, ONLY: sshn
  USE obs_inter_h2d
  USE obs_utils
  USE obs_inter_sup
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: obs_app_sstbias
  CONTAINS
  SUBROUTINE obs_app_sstbias(sstdata, k2dint, knumtypes, cl_bias_files)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE iom
    USE netcdf
    TYPE(obs_surf), INTENT(INOUT) :: sstdata
    INTEGER, INTENT(IN) :: k2dint
    INTEGER, INTENT(IN) :: knumtypes
    CHARACTER(LEN = 128), DIMENSION(knumtypes), INTENT(IN) :: cl_bias_files
    INTEGER :: jobs
    INTEGER :: jpisstbias
    INTEGER :: jpjsstbias
    INTEGER :: iico
    INTEGER :: ijco
    INTEGER :: jt
    INTEGER :: i_nx_id
    INTEGER :: i_ny_id
    INTEGER :: i_file_id
    INTEGER :: i_var_id
    INTEGER, DIMENSION(knumtypes) :: ibiastypes
    REAL(KIND = wp), DIMENSION(jpi, jpj, knumtypes) :: z_sstbias
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z_sstbias_2d
    REAL(KIND = wp), DIMENSION(1) :: zext, zobsmask
    REAL(KIND = wp), DIMENSION(2, 2, 1) :: zweig
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: zmask, zglam, zgphi
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: zmask_tmp, zglam_tmp, zgphi_tmp
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zbias
    REAL(KIND = wp) :: zlam
    REAL(KIND = wp) :: zphi
    INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: igrdi, igrdj
    INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: igrdi_tmp, igrdj_tmp
    INTEGER :: numsstbias
    INTEGER(KIND = NF90_INT) :: ifile_source
    INTEGER :: incfile
    INTEGER :: jtype
    INTEGER :: iret
    INTEGER :: inumtype
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    CALL profile_psy_data0 % PreStart('obs_app_sstbias', 'r0', 0, 0)
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) 'obs_rea_sstbias : '
    IF (lwp) WRITE(numout, FMT = *) '----------------- '
    IF (lwp) WRITE(numout, FMT = *) 'Read SST bias '
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    z_sstbias(:, :, :) = 0.0_wp
    !$ACC END KERNELS
    DO jtype = 1, knumtypes
      CALL profile_psy_data1 % PreStart('obs_app_sstbias', 'r1', 0, 0)
      numsstbias = 0
      IF (lwp) WRITE(numout, FMT = *) 'Opening ', cl_bias_files(jtype)
      CALL iom_open(cl_bias_files(jtype), numsstbias, ldstop = .FALSE.)
      CALL profile_psy_data1 % PostEnd
      IF (numsstbias > 0) THEN
        CALL profile_psy_data2 % PreStart('obs_app_sstbias', 'r2', 0, 0)
        iret = NF90_OPEN(TRIM(cl_bias_files(jtype)), NF90_NOWRITE, incfile)
        iret = NF90_GET_ATT(incfile, NF90_GLOBAL, "SST_source", ifile_source)
        ibiastypes(jtype) = ifile_source
        iret = NF90_CLOSE(incfile)
        IF (iret /= 0) CALL ctl_stop('obs_rea_sstbias : Cannot read bias type from file ' // cl_bias_files(jtype))
        CALL iom_get(numsstbias, jpdom_data, 'tn', z_sstbias_2d(:, :), 1)
        CALL profile_psy_data2 % PostEnd
        !$ACC KERNELS
        z_sstbias(:, :, jtype) = z_sstbias_2d(:, :)
        !$ACC END KERNELS
        CALL iom_close(numsstbias)
      ELSE
        CALL ctl_stop('obs_read_sstbias: File ' // TRIM(cl_bias_files(jtype)) // ' Not found')
      END IF
    END DO
    CALL profile_psy_data3 % PreStart('obs_app_sstbias', 'r3', 0, 0)
    ALLOCATE(igrdi(2, 2, sstdata % nsurf), igrdj(2, 2, sstdata % nsurf), zglam(2, 2, sstdata % nsurf), zgphi(2, 2, sstdata % &
&nsurf), zmask(2, 2, sstdata % nsurf))
    DO jobs = 1, sstdata % nsurf
      igrdi(1, 1, jobs) = sstdata % mi(jobs) - 1
      igrdj(1, 1, jobs) = sstdata % mj(jobs) - 1
      igrdi(1, 2, jobs) = sstdata % mi(jobs) - 1
      igrdj(1, 2, jobs) = sstdata % mj(jobs)
      igrdi(2, 1, jobs) = sstdata % mi(jobs)
      igrdj(2, 1, jobs) = sstdata % mj(jobs) - 1
      igrdi(2, 2, jobs) = sstdata % mi(jobs)
      igrdj(2, 2, jobs) = sstdata % mj(jobs)
    END DO
    CALL obs_int_comm_2d(2, 2, sstdata % nsurf, jpi, jpj, igrdi, igrdj, glamt, zglam)
    CALL obs_int_comm_2d(2, 2, sstdata % nsurf, jpi, jpj, igrdi, igrdj, gphit, zgphi)
    CALL obs_int_comm_2d(2, 2, sstdata % nsurf, jpi, jpj, igrdi, igrdj, tmask(:, :, 1), zmask)
    CALL profile_psy_data3 % PostEnd
    DO jtype = 1, knumtypes
      CALL profile_psy_data4 % PreStart('obs_app_sstbias', 'r4', 0, 0)
      inumtype = COUNT(sstdata % ntyp(:) == ibiastypes(jtype))
      ALLOCATE(igrdi_tmp(2, 2, inumtype), igrdj_tmp(2, 2, inumtype), zglam_tmp(2, 2, inumtype), zgphi_tmp(2, 2, inumtype), &
&zmask_tmp(2, 2, inumtype), zbias(2, 2, inumtype))
      jt = 1
      CALL profile_psy_data4 % PostEnd
      DO jobs = 1, sstdata % nsurf
        IF (sstdata % ntyp(jobs) == ibiastypes(jtype)) THEN
          !$ACC KERNELS
          igrdi_tmp(:, :, jt) = igrdi(:, :, jobs)
          igrdj_tmp(:, :, jt) = igrdj(:, :, jobs)
          zglam_tmp(:, :, jt) = zglam(:, :, jobs)
          zgphi_tmp(:, :, jt) = zgphi(:, :, jobs)
          zmask_tmp(:, :, jt) = zmask(:, :, jobs)
          jt = jt + 1
          !$ACC END KERNELS
        END IF
      END DO
      CALL profile_psy_data5 % PreStart('obs_app_sstbias', 'r5', 0, 0)
      CALL obs_int_comm_2d(2, 2, inumtype, jpi, jpj, igrdi_tmp(:, :, :), igrdj_tmp(:, :, :), z_sstbias(:, :, jtype), zbias(:, :, :))
      jt = 1
      DO jobs = 1, sstdata % nsurf
        IF (sstdata % ntyp(jobs) == ibiastypes(jtype)) THEN
          zlam = sstdata % rlam(jobs)
          zphi = sstdata % rphi(jobs)
          iico = sstdata % mi(jobs)
          ijco = sstdata % mj(jobs)
          CALL obs_int_h2d_init(1, 1, k2dint, zlam, zphi, zglam_tmp(:, :, jt), zgphi_tmp(:, :, jt), zmask_tmp(:, :, jt), zweig, &
&zobsmask)
          CALL obs_int_h2d(1, 1, zweig, zbias(:, :, jt), zext)
          sstdata % robs(jobs, 1) = sstdata % robs(jobs, 1) - zext(1)
          jt = jt + 1
        END IF
      END DO
      DEALLOCATE(igrdi_tmp, igrdj_tmp, zglam_tmp, zgphi_tmp, zmask_tmp, zbias)
      CALL profile_psy_data5 % PostEnd
    END DO
    CALL profile_psy_data6 % PreStart('obs_app_sstbias', 'r6', 0, 0)
    DEALLOCATE(igrdi, igrdj, zglam, zgphi, zmask)
    IF (lwp) THEN
      WRITE(numout, FMT = *) " "
      WRITE(numout, FMT = *) "SST bias correction applied successfully"
      WRITE(numout, FMT = *) "Obs types: ", ibiastypes(:), " Have all been bias corrected\n"
    END IF
    CALL profile_psy_data6 % PostEnd
  END SUBROUTINE obs_app_sstbias
END MODULE obs_sstbias