MODULE obs_oper
  USE obs_inter_sup
  USE obs_inter_h2d, ONLY: obs_int_h2d, obs_int_h2d_init
  USE obs_averg_h2d, ONLY: obs_avg_h2d, obs_avg_h2d_init, obs_max_fpsize
  USE obs_inter_z1d, ONLY: obs_int_z1d, obs_int_z1d_spl
  USE obs_const, ONLY: obfillflt
  USE dom_oce, ONLY: glamt, glamf, gphit, gphif
  USE lib_mpp, ONLY: ctl_warn, ctl_stop
  USE sbcdcy, ONLY: sbc_dcy, nday_qsr
  USE obs_grid, ONLY: obs_level_search
  USE par_kind, ONLY: wp
  USE in_out_manager
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: obs_prof_opt
  PUBLIC :: obs_surf_opt
  INTEGER, PARAMETER, PUBLIC :: imaxavtypes = 20
  CONTAINS
  SUBROUTINE obs_prof_opt(prodatqc, kt, kpi, kpj, kpk, kit000, kdaystp, pvar1, pvar2, pgdept, pgdepw, pmask1, pmask2, plam1, &
&plam2, pphi1, pphi2, k1dint, k2dint, kdailyavtypes)
    USE obs_profiles_def
    IMPLICIT NONE
    TYPE(obs_prof), INTENT(INOUT) :: prodatqc
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kpi, kpj, kpk
    INTEGER, INTENT(IN) :: kit000
    INTEGER, INTENT(IN) :: k1dint
    INTEGER, INTENT(IN) :: k2dint
    INTEGER, INTENT(IN) :: kdaystp
    REAL(KIND = wp), INTENT(IN), DIMENSION(kpi, kpj, kpk) :: pvar1, pvar2
    REAL(KIND = wp), INTENT(IN), DIMENSION(kpi, kpj, kpk) :: pmask1, pmask2
    REAL(KIND = wp), INTENT(IN), DIMENSION(kpi, kpj) :: plam1, plam2
    REAL(KIND = wp), INTENT(IN), DIMENSION(kpi, kpj) :: pphi1, pphi2
    REAL(KIND = wp), INTENT(IN), DIMENSION(kpi, kpj, kpk) :: pgdept, pgdepw
    INTEGER, DIMENSION(imaxavtypes), OPTIONAL :: kdailyavtypes
    INTEGER :: ji
    INTEGER :: jj
    INTEGER :: jk
    INTEGER :: jobs
    INTEGER :: inrc
    INTEGER :: ipro
    INTEGER :: idayend
    INTEGER :: ista
    INTEGER :: iend
    INTEGER :: iobs
    INTEGER :: iin, ijn, ikn, ik
    INTEGER :: inum_obs
    INTEGER, DIMENSION(imaxavtypes) :: idailyavtypes
    INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: igrdi1, igrdi2, igrdj1, igrdj2
    INTEGER, ALLOCATABLE, DIMENSION(:) :: iv_indic
    REAL(KIND = wp) :: zlam
    REAL(KIND = wp) :: zphi
    REAL(KIND = wp) :: zdaystp
    REAL(KIND = wp), DIMENSION(kpk) :: zobsmask1, zobsmask2, zobsk, zobs2k
    REAL(KIND = wp), DIMENSION(2, 2, 1) :: zweig1, zweig2, zweig
    REAL(KIND = wp), DIMENSION(:, :, :, :), ALLOCATABLE :: zmask1, zmask2, zint1, zint2, zinm1, zinm2, zgdept, zgdepw
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: zglam1, zglam2, zgphi1, zgphi2
    REAL(KIND = wp), DIMENSION(1) :: zmsk_1, zmsk_2
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: interp_corner
    LOGICAL :: ld_dailyav
    inrc = kt - kit000 + 2
    ipro = prodatqc % npstp(inrc)
    ld_dailyav = .FALSE.
    IF (PRESENT(kdailyavtypes)) THEN
      idailyavtypes(:) = kdailyavtypes(:)
      IF (ANY(idailyavtypes(:) /= - 1)) ld_dailyav = .TRUE.
    ELSE
      idailyavtypes(:) = - 1
    END IF
    idayend = MOD(kt - kit000 + 1, kdaystp)
    IF (ld_dailyav) THEN
      IF (idayend == 1 .OR. kt == 0) THEN
        DO jk = 1, jpk
          DO jj = 1, jpj
            DO ji = 1, jpi
              prodatqc % vdmean(ji, jj, jk, 1) = 0.0
              prodatqc % vdmean(ji, jj, jk, 2) = 0.0
            END DO
          END DO
        END DO
      END IF
      DO jk = 1, jpk
        DO jj = 1, jpj
          DO ji = 1, jpi
            prodatqc % vdmean(ji, jj, jk, 1) = prodatqc % vdmean(ji, jj, jk, 1) + pvar1(ji, jj, jk)
            prodatqc % vdmean(ji, jj, jk, 2) = prodatqc % vdmean(ji, jj, jk, 2) + pvar2(ji, jj, jk)
          END DO
        END DO
      END DO
      zdaystp = 1.0 / REAL(kdaystp)
      IF (idayend == 0) THEN
        IF (lwp) WRITE(numout, FMT = *) 'Calculating prodatqc%vdmean on time-step: ', kt
        CALL FLUSH(numout)
        DO jk = 1, jpk
          DO jj = 1, jpj
            DO ji = 1, jpi
              prodatqc % vdmean(ji, jj, jk, 1) = prodatqc % vdmean(ji, jj, jk, 1) * zdaystp
              prodatqc % vdmean(ji, jj, jk, 2) = prodatqc % vdmean(ji, jj, jk, 2) * zdaystp
            END DO
          END DO
        END DO
      END IF
    END IF
    ALLOCATE(igrdi1(2, 2, ipro), igrdi2(2, 2, ipro), igrdj1(2, 2, ipro), igrdj2(2, 2, ipro), zglam1(2, 2, ipro), zglam2(2, 2, &
&ipro), zgphi1(2, 2, ipro), zgphi2(2, 2, ipro), zmask1(2, 2, kpk, ipro), zmask2(2, 2, kpk, ipro), zint1(2, 2, kpk, ipro), zint2(2, &
&2, kpk, ipro), zgdept(2, 2, kpk, ipro), zgdepw(2, 2, kpk, ipro))
    DO jobs = prodatqc % nprofup + 1, prodatqc % nprofup + ipro
      iobs = jobs - prodatqc % nprofup
      igrdi1(1, 1, iobs) = prodatqc % mi(jobs, 1) - 1
      igrdj1(1, 1, iobs) = prodatqc % mj(jobs, 1) - 1
      igrdi1(1, 2, iobs) = prodatqc % mi(jobs, 1) - 1
      igrdj1(1, 2, iobs) = prodatqc % mj(jobs, 1)
      igrdi1(2, 1, iobs) = prodatqc % mi(jobs, 1)
      igrdj1(2, 1, iobs) = prodatqc % mj(jobs, 1) - 1
      igrdi1(2, 2, iobs) = prodatqc % mi(jobs, 1)
      igrdj1(2, 2, iobs) = prodatqc % mj(jobs, 1)
      igrdi2(1, 1, iobs) = prodatqc % mi(jobs, 2) - 1
      igrdj2(1, 1, iobs) = prodatqc % mj(jobs, 2) - 1
      igrdi2(1, 2, iobs) = prodatqc % mi(jobs, 2) - 1
      igrdj2(1, 2, iobs) = prodatqc % mj(jobs, 2)
      igrdi2(2, 1, iobs) = prodatqc % mi(jobs, 2)
      igrdj2(2, 1, iobs) = prodatqc % mj(jobs, 2) - 1
      igrdi2(2, 2, iobs) = prodatqc % mi(jobs, 2)
      igrdj2(2, 2, iobs) = prodatqc % mj(jobs, 2)
    END DO
    zgdept(:, :, :, :) = 0.0
    zgdepw(:, :, :, :) = 0.0
    CALL obs_int_comm_2d(2, 2, ipro, kpi, kpj, igrdi1, igrdj1, plam1, zglam1)
    CALL obs_int_comm_2d(2, 2, ipro, kpi, kpj, igrdi1, igrdj1, pphi1, zgphi1)
    CALL obs_int_comm_3d(2, 2, ipro, kpi, kpj, kpk, igrdi1, igrdj1, pmask1, zmask1)
    CALL obs_int_comm_3d(2, 2, ipro, kpi, kpj, kpk, igrdi1, igrdj1, pvar1, zint1)
    CALL obs_int_comm_2d(2, 2, ipro, kpi, kpj, igrdi2, igrdj2, plam2, zglam2)
    CALL obs_int_comm_2d(2, 2, ipro, kpi, kpj, igrdi2, igrdj2, pphi2, zgphi2)
    CALL obs_int_comm_3d(2, 2, ipro, kpi, kpj, kpk, igrdi2, igrdj2, pmask2, zmask2)
    CALL obs_int_comm_3d(2, 2, ipro, kpi, kpj, kpk, igrdi2, igrdj2, pvar2, zint2)
    CALL obs_int_comm_3d(2, 2, ipro, kpi, kpj, kpk, igrdi1, igrdj1, pgdept, zgdept)
    CALL obs_int_comm_3d(2, 2, ipro, kpi, kpj, kpk, igrdi1, igrdj1, pgdepw, zgdepw)
    IF (ld_dailyav .AND. idayend == 0) THEN
      ALLOCATE(zinm1(2, 2, kpk, ipro), zinm2(2, 2, kpk, ipro))
      CALL obs_int_comm_3d(2, 2, ipro, kpi, kpj, kpk, igrdi1, igrdj1, prodatqc % vdmean(:, :, :, 1), zinm1)
      CALL obs_int_comm_3d(2, 2, ipro, kpi, kpj, kpk, igrdi2, igrdj2, prodatqc % vdmean(:, :, :, 2), zinm2)
    END IF
    IF (ipro == 0) RETURN
    DO jobs = prodatqc % nprofup + 1, prodatqc % nprofup + ipro
      iobs = jobs - prodatqc % nprofup
      IF (kt /= prodatqc % mstp(jobs)) THEN
        IF (lwp) THEN
          WRITE(numout, FMT = *)
          WRITE(numout, FMT = *) ' E R R O R : Observation', ' time step is not consistent with the', ' model time step'
          WRITE(numout, FMT = *) ' ========='
          WRITE(numout, FMT = *)
          WRITE(numout, FMT = *) ' Record  = ', jobs, ' kt      = ', kt, ' mstp    = ', prodatqc % mstp(jobs), ' ntyp    = ', &
&prodatqc % ntyp(jobs)
        END IF
        CALL ctl_stop('obs_pro_opt', 'Inconsistent time')
      END IF
      zlam = prodatqc % rlam(jobs)
      zphi = prodatqc % rphi(jobs)
      IF (prodatqc % npvend(jobs, 1) > 0) THEN
        CALL obs_int_h2d_init(1, 1, k2dint, zlam, zphi, zglam1(:, :, iobs), zgphi1(:, :, iobs), zmask1(:, :, 1, iobs), zweig1, &
&zmsk_1)
      END IF
      IF (prodatqc % npvend(jobs, 2) > 0) THEN
        CALL obs_int_h2d_init(1, 1, k2dint, zlam, zphi, zglam2(:, :, iobs), zgphi2(:, :, iobs), zmask2(:, :, 1, iobs), zweig2, &
&zmsk_2)
      END IF
      IF (prodatqc % npvend(jobs, 1) > 0) THEN
        zobsk(:) = obfillflt
        IF (ANY(idailyavtypes(:) == prodatqc % ntyp(jobs))) THEN
          IF (idayend == 0) THEN
            ista = prodatqc % npvsta(jobs, 1)
            iend = prodatqc % npvend(jobs, 1)
            inum_obs = iend - ista + 1
            ALLOCATE(interp_corner(2, 2, inum_obs), iv_indic(inum_obs))
            DO iin = 1, 2
              DO ijn = 1, 2
                IF (k1dint == 1) THEN
                  CALL obs_int_z1d_spl(kpk, zinm1(iin, ijn, :, iobs), zobs2k, zgdept(iin, ijn, :, iobs), zmask1(iin, ijn, :, iobs))
                END IF
                CALL obs_level_search(kpk, zgdept(iin, ijn, :, iobs), inum_obs, prodatqc % var(1) % vdep(ista : iend), iv_indic)
                CALL obs_int_z1d(kpk, iv_indic, k1dint, inum_obs, prodatqc % var(1) % vdep(ista : iend), zinm1(iin, ijn, :, iobs), &
&zobs2k, interp_corner(iin, ijn, :), zgdept(iin, ijn, :, iobs), zmask1(iin, ijn, :, iobs))
              END DO
            END DO
          END IF
        ELSE
          ista = prodatqc % npvsta(jobs, 1)
          iend = prodatqc % npvend(jobs, 1)
          inum_obs = iend - ista + 1
          ALLOCATE(interp_corner(2, 2, inum_obs), iv_indic(inum_obs))
          DO iin = 1, 2
            DO ijn = 1, 2
              IF (k1dint == 1) THEN
                CALL obs_int_z1d_spl(kpk, zint1(iin, ijn, :, iobs), zobs2k, zgdept(iin, ijn, :, iobs), zmask1(iin, ijn, :, iobs))
              END IF
              CALL obs_level_search(kpk, zgdept(iin, ijn, :, iobs), inum_obs, prodatqc % var(1) % vdep(ista : iend), iv_indic)
              CALL obs_int_z1d(kpk, iv_indic, k1dint, inum_obs, prodatqc % var(1) % vdep(ista : iend), zint1(iin, ijn, :, iobs), &
&zobs2k, interp_corner(iin, ijn, :), zgdept(iin, ijn, :, iobs), zmask1(iin, ijn, :, iobs))
            END DO
          END DO
        END IF
        DO ikn = 1, inum_obs
          iend = ista + ikn - 1
          zweig(:, :, 1) = 0._wp
          DO iin = 1, 2
            DO ijn = 1, 2
              depth_loop1:DO ik = kpk, 2, - 1
                IF (zmask1(iin, ijn, ik - 1, iobs) > 0.9) THEN
                  zweig(iin, ijn, 1) = zweig1(iin, ijn, 1) * MAX(SIGN(1._wp, (zgdepw(iin, ijn, ik, iobs)) - prodatqc % var(1) % &
&vdep(iend)), 0._wp)
                  EXIT depth_loop1
                END IF
              END DO depth_loop1
            END DO
          END DO
          CALL obs_int_h2d(1, 1, zweig, interp_corner(:, :, ikn), prodatqc % var(1) % vmod(iend : iend))
          IF (SUM(zweig) == 0.0_wp) prodatqc % var(1) % nvqc(iend : iend) = 4
        END DO
        DEALLOCATE(interp_corner, iv_indic)
      END IF
      IF (prodatqc % npvend(jobs, 2) > 0) THEN
        zobsk(:) = obfillflt
        IF (ANY(idailyavtypes(:) == prodatqc % ntyp(jobs))) THEN
          IF (idayend == 0) THEN
            ista = prodatqc % npvsta(jobs, 2)
            iend = prodatqc % npvend(jobs, 2)
            inum_obs = iend - ista + 1
            ALLOCATE(interp_corner(2, 2, inum_obs), iv_indic(inum_obs))
            DO iin = 1, 2
              DO ijn = 1, 2
                IF (k1dint == 1) THEN
                  CALL obs_int_z1d_spl(kpk, zinm2(iin, ijn, :, iobs), zobs2k, zgdept(iin, ijn, :, iobs), zmask2(iin, ijn, :, iobs))
                END IF
                CALL obs_level_search(kpk, zgdept(iin, ijn, :, iobs), inum_obs, prodatqc % var(2) % vdep(ista : iend), iv_indic)
                CALL obs_int_z1d(kpk, iv_indic, k1dint, inum_obs, prodatqc % var(2) % vdep(ista : iend), zinm2(iin, ijn, :, iobs), &
&zobs2k, interp_corner(iin, ijn, :), zgdept(iin, ijn, :, iobs), zmask2(iin, ijn, :, iobs))
              END DO
            END DO
          END IF
        ELSE
          ista = prodatqc % npvsta(jobs, 2)
          iend = prodatqc % npvend(jobs, 2)
          inum_obs = iend - ista + 1
          ALLOCATE(interp_corner(2, 2, inum_obs), iv_indic(inum_obs))
          DO iin = 1, 2
            DO ijn = 1, 2
              IF (k1dint == 1) THEN
                CALL obs_int_z1d_spl(kpk, zint2(iin, ijn, :, iobs), zobs2k, zgdept(iin, ijn, :, iobs), zmask2(iin, ijn, :, iobs))
              END IF
              CALL obs_level_search(kpk, zgdept(iin, ijn, :, iobs), inum_obs, prodatqc % var(2) % vdep(ista : iend), iv_indic)
              CALL obs_int_z1d(kpk, iv_indic, k1dint, inum_obs, prodatqc % var(2) % vdep(ista : iend), zint2(iin, ijn, :, iobs), &
&zobs2k, interp_corner(iin, ijn, :), zgdept(iin, ijn, :, iobs), zmask2(iin, ijn, :, iobs))
            END DO
          END DO
        END IF
        DO ikn = 1, inum_obs
          iend = ista + ikn - 1
          zweig(:, :, 1) = 0._wp
          DO iin = 1, 2
            DO ijn = 1, 2
              depth_loop2:DO ik = kpk, 2, - 1
                IF (zmask2(iin, ijn, ik - 1, iobs) > 0.9) THEN
                  zweig(iin, ijn, 1) = zweig2(iin, ijn, 1) * MAX(SIGN(1._wp, (zgdepw(iin, ijn, ik, iobs)) - prodatqc % var(2) % &
&vdep(iend)), 0._wp)
                  EXIT depth_loop2
                END IF
              END DO depth_loop2
            END DO
          END DO
          CALL obs_int_h2d(1, 1, zweig, interp_corner(:, :, ikn), prodatqc % var(2) % vmod(iend : iend))
          IF (SUM(zweig) == 0.0_wp) prodatqc % var(2) % nvqc(iend : iend) = 4
        END DO
        DEALLOCATE(interp_corner, iv_indic)
      END IF
    END DO
    DEALLOCATE(igrdi1, igrdi2, igrdj1, igrdj2, zglam1, zglam2, zgphi1, zgphi2, zmask1, zmask2, zint1, zint2, zgdept, zgdepw)
    IF (ld_dailyav .AND. idayend == 0) THEN
      DEALLOCATE(zinm1, zinm2)
    END IF
    prodatqc % nprofup = prodatqc % nprofup + ipro
  END SUBROUTINE obs_prof_opt
  SUBROUTINE obs_surf_opt(surfdataqc, kt, kpi, kpj, kit000, kdaystp, psurf, psurfmask, k2dint, ldnightav, plamscl, pphiscl, &
&lindegrees)
    USE obs_surf_def
    IMPLICIT NONE
    TYPE(obs_surf), INTENT(INOUT) :: surfdataqc
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kpi
    INTEGER, INTENT(IN) :: kpj
    INTEGER, INTENT(IN) :: kit000
    INTEGER, INTENT(IN) :: kdaystp
    INTEGER, INTENT(IN) :: k2dint
    REAL(KIND = wp), INTENT(IN), DIMENSION(kpi, kpj) :: psurf, psurfmask
    LOGICAL, INTENT(IN) :: ldnightav
    REAL(KIND = wp), INTENT(IN) :: plamscl, pphiscl
    LOGICAL, INTENT(IN) :: lindegrees
    INTEGER :: ji
    INTEGER :: jj
    INTEGER :: jobs
    INTEGER :: inrc
    INTEGER :: isurf
    INTEGER :: iobs
    INTEGER :: imaxifp, imaxjfp
    INTEGER :: imodi, imodj
    INTEGER :: idayend
    INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: igrdi, igrdj, igrdip1, igrdjp1
    INTEGER, DIMENSION(:, :), SAVE, ALLOCATABLE :: icount_night, imask_night
    REAL(KIND = wp) :: zlam
    REAL(KIND = wp) :: zphi
    REAL(KIND = wp), DIMENSION(1) :: zext, zobsmask
    REAL(KIND = wp) :: zdaystp
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: zweig, zmask, zsurf, zsurfm, zsurftmp, zglam, zgphi, zglamf, zgphif
    REAL(KIND = wp), DIMENSION(:, :), SAVE, ALLOCATABLE :: zintmp, zouttmp, zmeanday
    inrc = kt - kit000 + 2
    isurf = surfdataqc % nsstp(inrc)
    CALL obs_max_fpsize(k2dint, plamscl, pphiscl, lindegrees, psurfmask, imaxifp, imaxjfp)
    IF (ldnightav) THEN
      IF (kt == 0) THEN
        ALLOCATE(icount_night(kpi, kpj))
        ALLOCATE(imask_night(kpi, kpj))
        ALLOCATE(zintmp(kpi, kpj))
        ALLOCATE(zouttmp(kpi, kpj))
        ALLOCATE(zmeanday(kpi, kpj))
        nday_qsr = - 1
      END IF
      idayend = MOD(kt - kit000 + 1, kdaystp)
      IF (idayend == 1 .OR. kt == 0) THEN
        DO jj = 1, jpj
          DO ji = 1, jpi
            surfdataqc % vdmean(ji, jj) = 0.0
            zmeanday(ji, jj) = 0.0
            icount_night(ji, jj) = 0
          END DO
        END DO
      END IF
      zintmp(:, :) = 0.0
      zouttmp(:, :) = sbc_dcy(zintmp(:, :), .TRUE.)
      imask_night(:, :) = INT(zouttmp(:, :))
      DO jj = 1, jpj
        DO ji = 1, jpi
          surfdataqc % vdmean(ji, jj) = surfdataqc % vdmean(ji, jj) + psurf(ji, jj) * REAL(imask_night(ji, jj))
          zmeanday(ji, jj) = zmeanday(ji, jj) + psurf(ji, jj)
          icount_night(ji, jj) = icount_night(ji, jj) + imask_night(ji, jj)
        END DO
      END DO
      zdaystp = 1.0 / REAL(kdaystp)
      IF (idayend == 0) THEN
        IF (lwp) WRITE(numout, FMT = *) 'Calculating surfdataqc%vdmean on time-step: ', kt
        DO jj = 1, jpj
          DO ji = 1, jpi
            IF (icount_night(ji, jj) > 0) THEN
              surfdataqc % vdmean(ji, jj) = surfdataqc % vdmean(ji, jj) / REAL(icount_night(ji, jj))
            ELSE
              surfdataqc % vdmean(ji, jj) = zmeanday(ji, jj) * zdaystp
            END IF
          END DO
        END DO
      END IF
    END IF
    ALLOCATE(zweig(imaxifp, imaxjfp, 1), igrdi(imaxifp, imaxjfp, isurf), igrdj(imaxifp, imaxjfp, isurf), zglam(imaxifp, imaxjfp, &
&isurf), zgphi(imaxifp, imaxjfp, isurf), zmask(imaxifp, imaxjfp, isurf), zsurf(imaxifp, imaxjfp, isurf), zsurftmp(imaxifp, &
&imaxjfp, isurf), zglamf(imaxifp + 1, imaxjfp + 1, isurf), zgphif(imaxifp + 1, imaxjfp + 1, isurf), igrdip1(imaxifp + 1, imaxjfp + &
&1, isurf), igrdjp1(imaxifp + 1, imaxjfp + 1, isurf))
    DO jobs = surfdataqc % nsurfup + 1, surfdataqc % nsurfup + isurf
      iobs = jobs - surfdataqc % nsurfup
      DO ji = 0, imaxifp
        imodi = surfdataqc % mi(jobs) - INT(imaxifp / 2) + ji - 1
        IF (imodi < 1) imodi = imodi + jpiglo
        IF (imodi > jpiglo) imodi = imodi - jpiglo
        DO jj = 0, imaxjfp
          imodj = surfdataqc % mj(jobs) - INT(imaxjfp / 2) + jj - 1
          IF (imodj < 1) imodj = 1
          IF (imodj > jpjglo) imodj = jpjglo
          igrdip1(ji + 1, jj + 1, iobs) = imodi
          igrdjp1(ji + 1, jj + 1, iobs) = imodj
          IF (ji >= 1 .AND. jj >= 1) THEN
            igrdi(ji, jj, iobs) = imodi
            igrdj(ji, jj, iobs) = imodj
          END IF
        END DO
      END DO
    END DO
    CALL obs_int_comm_2d(imaxifp, imaxjfp, isurf, kpi, kpj, igrdi, igrdj, glamt, zglam)
    CALL obs_int_comm_2d(imaxifp, imaxjfp, isurf, kpi, kpj, igrdi, igrdj, gphit, zgphi)
    CALL obs_int_comm_2d(imaxifp, imaxjfp, isurf, kpi, kpj, igrdi, igrdj, psurfmask, zmask)
    CALL obs_int_comm_2d(imaxifp, imaxjfp, isurf, kpi, kpj, igrdi, igrdj, psurf, zsurf)
    CALL obs_int_comm_2d(imaxifp + 1, imaxjfp + 1, isurf, kpi, kpj, igrdip1, igrdjp1, glamf, zglamf)
    CALL obs_int_comm_2d(imaxifp + 1, imaxjfp + 1, isurf, kpi, kpj, igrdip1, igrdjp1, gphif, zgphif)
    IF (idayend == 0 .AND. ldnightav) THEN
      ALLOCATE(zsurfm(imaxifp, imaxjfp, isurf))
      CALL obs_int_comm_2d(imaxifp, imaxjfp, isurf, kpi, kpj, igrdi, igrdj, surfdataqc % vdmean(:, :), zsurfm)
    END IF
    DO jobs = surfdataqc % nsurfup + 1, surfdataqc % nsurfup + isurf
      iobs = jobs - surfdataqc % nsurfup
      IF (kt /= surfdataqc % mstp(jobs)) THEN
        IF (lwp) THEN
          WRITE(numout, FMT = *)
          WRITE(numout, FMT = *) ' E R R O R : Observation', ' time step is not consistent with the', ' model time step'
          WRITE(numout, FMT = *) ' ========='
          WRITE(numout, FMT = *)
          WRITE(numout, FMT = *) ' Record  = ', jobs, ' kt      = ', kt, ' mstp    = ', surfdataqc % mstp(jobs), ' ntyp    = ', &
&surfdataqc % ntyp(jobs)
        END IF
        CALL ctl_stop('obs_surf_opt', 'Inconsistent time')
      END IF
      zlam = surfdataqc % rlam(jobs)
      zphi = surfdataqc % rphi(jobs)
      IF (ldnightav .AND. idayend == 0) THEN
        zsurftmp(:, :, iobs) = zsurfm(:, :, iobs)
      ELSE
        zsurftmp(:, :, iobs) = zsurf(:, :, iobs)
      END IF
      IF (k2dint <= 4) THEN
        CALL obs_int_h2d_init(1, 1, k2dint, zlam, zphi, zglam(:, :, iobs), zgphi(:, :, iobs), zmask(:, :, iobs), zweig, zobsmask)
        CALL obs_int_h2d(1, 1, zweig, zsurftmp(:, :, iobs), zext)
      ELSE
        CALL obs_avg_h2d_init(1, 1, imaxifp, imaxjfp, k2dint, zlam, zphi, zglam(:, :, iobs), zgphi(:, :, iobs), zglamf(:, :, &
&iobs), zgphif(:, :, iobs), zmask(:, :, iobs), plamscl, pphiscl, lindegrees, zweig, zobsmask)
        CALL obs_avg_h2d(1, 1, imaxifp, imaxjfp, zweig, zsurftmp(:, :, iobs), zext)
      END IF
      IF (TRIM(surfdataqc % cvars(1)) == 'SLA' .AND. surfdataqc % nextra == 2) THEN
        surfdataqc % rext(jobs, 1) = zext(1)
        surfdataqc % rmod(jobs, 1) = surfdataqc % rext(jobs, 1) - surfdataqc % rext(jobs, 2)
      ELSE
        surfdataqc % rmod(jobs, 1) = zext(1)
      END IF
      IF (zext(1) == obfillflt) THEN
        surfdataqc % nqc(jobs) = 4
      END IF
    END DO
    DEALLOCATE(zweig, igrdi, igrdj, zglam, zgphi, zmask, zsurf, zsurftmp, zglamf, zgphif, igrdip1, igrdjp1)
    IF (idayend == 0 .AND. ldnightav) THEN
      DEALLOCATE(zsurfm)
    END IF
    surfdataqc % nsurfup = surfdataqc % nsurfup + isurf
  END SUBROUTINE obs_surf_opt
END MODULE obs_oper