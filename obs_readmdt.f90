MODULE obs_readmdt
  USE par_kind
  USE par_oce
  USE in_out_manager
  USE obs_surf_def
  USE obs_inter_sup
  USE obs_inter_h2d
  USE obs_utils
  USE iom_nf90
  USE netcdf
  USE lib_mpp
  USE dom_oce, ONLY: tmask, tmask_i, e1e2t, gphit, glamt
  USE obs_const, ONLY: obfillflt
  USE oce, ONLY: sshn
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: obs_rea_mdt
  PUBLIC :: obs_offset_mdt
  INTEGER, PUBLIC :: nn_msshc = 1
  REAL(KIND = wp), PUBLIC :: rn_mdtcorr = 1.61_wp
  REAL(KIND = wp), PUBLIC :: rn_mdtcutoff = 65.0_wp
  CONTAINS
  SUBROUTINE obs_rea_mdt(sladata, k2dint)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE iom
    TYPE(obs_surf), INTENT(INOUT) :: sladata
    INTEGER, INTENT(IN) :: k2dint
    CHARACTER(LEN = 12), PARAMETER :: cpname = 'obs_rea_mdt'
    CHARACTER(LEN = 20), PARAMETER :: mdtname = 'slaReferenceLevel.nc'
    INTEGER :: jobs
    INTEGER :: jpimdt, jpjmdt
    INTEGER :: iico, ijco
    INTEGER :: i_nx_id, i_ny_id, i_file_id, i_var_id, i_stat
    INTEGER :: nummdt
    REAL(KIND = wp), DIMENSION(1) :: zext, zobsmask
    REAL(KIND = wp), DIMENSION(2, 2, 1) :: zweig
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: zmask, zmdtl, zglam, zgphi
    INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: igrdi, igrdj
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z_mdt, mdtmask
    REAL(KIND = wp) :: zlam, zphi, zfill, zinfill
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('obs_rea_mdt', 'r0', 0, 0)
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) ' obs_rea_mdt : Read MDT for referencing altimeter anomalies'
    IF (lwp) WRITE(numout, FMT = *) ' ------------- '
    CALL FLUSH(numout)
    CALL iom_open(mdtname, nummdt)
    CALL iom_get(nummdt, jpdom_data, 'sossheig', z_mdt(:, :), 1)
    CALL iom_close(nummdt)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zinfill = 0.0
    i_stat = nf90_open(mdtname, nf90_nowrite, nummdt)
    i_stat = nf90_inq_varid(nummdt, 'sossheig', i_var_id)
    i_stat = nf90_get_att(nummdt, i_var_id, "_FillValue", zinfill)
    zfill = zinfill
    i_stat = nf90_close(nummdt)
    WHERE (z_mdt(:, :) /= zfill)
      mdtmask(:, :) = tmask(:, :, 1)
    ELSEWHERE
      mdtmask(:, :) = 0
    END WHERE
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('obs_rea_mdt', 'r1', 0, 0)
    IF (nn_msshc == 1 .OR. nn_msshc == 2) CALL obs_offset_mdt(jpi, jpj, z_mdt, zfill)
    ALLOCATE(igrdi(2, 2, sladata % nsurf), igrdj(2, 2, sladata % nsurf), zglam(2, 2, sladata % nsurf), zgphi(2, 2, sladata % &
&nsurf), zmask(2, 2, sladata % nsurf), zmdtl(2, 2, sladata % nsurf))
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
    CALL obs_int_comm_2d(2, 2, sladata % nsurf, jpi, jpj, igrdi, igrdj, mdtmask, zmask)
    CALL obs_int_comm_2d(2, 2, sladata % nsurf, jpi, jpj, igrdi, igrdj, z_mdt, zmdtl)
    DO jobs = 1, sladata % nsurf
      zlam = sladata % rlam(jobs)
      zphi = sladata % rphi(jobs)
      CALL obs_int_h2d_init(1, 1, k2dint, zlam, zphi, zglam(:, :, jobs), zgphi(:, :, jobs), zmask(:, :, jobs), zweig, zobsmask)
      CALL obs_int_h2d(1, 1, zweig, zmdtl(:, :, jobs), zext)
      sladata % rext(jobs, 2) = zext(1)
      IF (zobsmask(1) == 0) sladata % nqc(jobs) = IBSET(sladata % nqc(jobs), 15)
    END DO
    DEALLOCATE(igrdi, igrdj, zglam, zgphi, zmask, zmdtl)
    IF (lwp) WRITE(numout, FMT = *) ' ------------- '
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE obs_rea_mdt
  SUBROUTINE obs_offset_mdt(kpi, kpj, mdt, zfill)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kpi, kpj
    REAL(KIND = wp), DIMENSION(kpi, kpj), INTENT(INOUT) :: mdt
    REAL(KIND = wp), INTENT(IN) :: zfill
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zdxdy, zarea, zeta1, zeta2, zcorr_mdt, zcorr_bcketa, zcorr
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zpromsk
    CHARACTER(LEN = 14), PARAMETER :: cpname = 'obs_offset_mdt'
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    !$ACC KERNELS
    DO ji = 1, jpi
      DO jj = 1, jpj
        zpromsk(ji, jj) = tmask_i(ji, jj)
        IF ((gphit(ji, jj) .GT. rn_mdtcutoff) .OR. (gphit(ji, jj) .LT. - rn_mdtcutoff) .OR. (mdt(ji, jj) .EQ. zfill)) zpromsk(ji, &
&jj) = 0.0
      END DO
    END DO
    zarea = 0.0
    zeta1 = 0.0
    zeta2 = 0.0
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zdxdy = e1e2t(ji, jj) * zpromsk(ji, jj)
        zarea = zarea + zdxdy
        zeta1 = zeta1 + mdt(ji, jj) * zdxdy
        zeta2 = zeta2 + sshn(ji, jj) * zdxdy
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('obs_offset_mdt', 'r0', 0, 0)
    CALL mpp_sum('obs_readmdt', zeta1)
    CALL mpp_sum('obs_readmdt', zeta2)
    CALL mpp_sum('obs_readmdt', zarea)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zcorr_mdt = zeta1 / zarea
    zcorr_bcketa = zeta2 / zarea
    zcorr = zcorr_mdt - zcorr_bcketa
    IF (nn_msshc == 1) mdt(:, :) = mdt(:, :) - zcorr
    IF (nn_msshc == 2) mdt(:, :) = mdt(:, :) - rn_mdtcorr
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('obs_offset_mdt', 'r1', 0, 0)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) ' obs_readmdt : rn_mdtcutoff     = ', rn_mdtcutoff
      WRITE(numout, FMT = *) ' -----------   zcorr_mdt     = ', zcorr_mdt
      WRITE(numout, FMT = *) '               zcorr_bcketa  = ', zcorr_bcketa
      WRITE(numout, FMT = *) '               zcorr         = ', zcorr
      WRITE(numout, FMT = *) '               nn_msshc        = ', nn_msshc
    END IF
    IF (nn_msshc == 0) WRITE(numout, FMT = *) '           MSSH correction is not applied'
    IF (nn_msshc == 1) WRITE(numout, FMT = *) '           MSSH correction is applied'
    IF (nn_msshc == 2) WRITE(numout, FMT = *) '           User defined MSSH correction'
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE obs_offset_mdt
END MODULE obs_readmdt