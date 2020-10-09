MODULE obs_surf_def
  USE par_kind, ONLY: wp
  USE obs_mpp, ONLY: obs_mpp_sum_integer
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: obs_surf, obs_surf_alloc, obs_surf_dealloc, obs_surf_compress, obs_surf_decompress
  TYPE :: obs_surf
    INTEGER :: nsurf
    INTEGER :: nsurfmpp
    INTEGER :: nvar
    INTEGER :: nextra
    INTEGER :: nstp
    INTEGER :: npi
    INTEGER :: npj
    INTEGER :: nsurfup
    INTEGER :: nrec
    INTEGER, POINTER, DIMENSION(:) :: mi, mj, mt, nsidx, nsfil, nyea, nmon, nday, nhou, nmin, mstp, nqc, ntyp
    CHARACTER(LEN = 8), POINTER, DIMENSION(:) :: cvars
    CHARACTER(LEN = 8), POINTER, DIMENSION(:) :: cwmo
    REAL(KIND = wp), POINTER, DIMENSION(:) :: rlam, rphi
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: robs, rmod
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: rext
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: vdmean
    INTEGER, POINTER, DIMENSION(:) :: nsstp, nsstpmpp
    INTEGER, POINTER, DIMENSION(:) :: mrecstp
    INTEGER, POINTER, DIMENSION(:) :: nsind
    LOGICAL :: lgrid
  END TYPE obs_surf
  CONTAINS
  SUBROUTINE obs_surf_alloc(surf, ksurf, kvar, kextra, kstp, kpi, kpj)
    TYPE(obs_surf), INTENT(INOUT) :: surf
    INTEGER, INTENT(IN) :: ksurf
    INTEGER, INTENT(IN) :: kvar
    INTEGER, INTENT(IN) :: kextra
    INTEGER, INTENT(IN) :: kstp
    INTEGER, INTENT(IN) :: kpi
    INTEGER, INTENT(IN) :: kpj
    INTEGER :: ji
    INTEGER :: jvar
    surf % nsurf = ksurf
    surf % nsurfmpp = 0
    surf % nextra = kextra
    surf % nvar = kvar
    surf % nstp = kstp
    surf % npi = kpi
    surf % npj = kpj
    ALLOCATE(surf % cvars(kvar))
    DO jvar = 1, kvar
      surf % cvars(jvar) = "NotSet"
    END DO
    ALLOCATE(surf % mi(ksurf), surf % mj(ksurf), surf % mt(ksurf), surf % nsidx(ksurf), surf % nsfil(ksurf), surf % nyea(ksurf), &
&surf % nmon(ksurf), surf % nday(ksurf), surf % nhou(ksurf), surf % nmin(ksurf), surf % mstp(ksurf), surf % nqc(ksurf), surf % &
&ntyp(ksurf), surf % cwmo(ksurf), surf % rlam(ksurf), surf % rphi(ksurf), surf % nsind(ksurf))
    surf % mt(:) = - 1
    ALLOCATE(surf % robs(ksurf, kvar), surf % rmod(ksurf, kvar))
    ALLOCATE(surf % rext(ksurf, kextra))
    surf % rext(:, :) = 0.0_wp
    ALLOCATE(surf % nsstp(kstp), surf % nsstpmpp(kstp))
    ALLOCATE(surf % vdmean(kpi, kpj))
    DO ji = 1, ksurf
      surf % nsind(ji) = ji
    END DO
    surf % nsstp(:) = 0
    surf % nsstpmpp(:) = 0
    surf % nsurfup = 0
    surf % lgrid = .FALSE.
  END SUBROUTINE obs_surf_alloc
  SUBROUTINE obs_surf_dealloc(surf)
    TYPE(obs_surf), INTENT(INOUT) :: surf
    DEALLOCATE(surf % mi, surf % mj, surf % mt, surf % nsidx, surf % nsfil, surf % nyea, surf % nmon, surf % nday, surf % nhou, &
&surf % nmin, surf % mstp, surf % nqc, surf % ntyp, surf % cwmo, surf % rlam, surf % rphi, surf % nsind)
    DEALLOCATE(surf % robs, surf % rmod)
    DEALLOCATE(surf % rext)
    DEALLOCATE(surf % vdmean)
    DEALLOCATE(surf % nsstp, surf % nsstpmpp)
    DEALLOCATE(surf % cvars)
  END SUBROUTINE obs_surf_dealloc
  SUBROUTINE obs_surf_compress(surf, newsurf, lallocate, kumout, lvalid)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(obs_surf), INTENT(IN) :: surf
    TYPE(obs_surf), INTENT(INOUT) :: newsurf
    LOGICAL :: lallocate
    INTEGER, INTENT(IN) :: kumout
    LOGICAL, OPTIONAL, INTENT(IN), DIMENSION(:) :: lvalid
    INTEGER :: insurf
    INTEGER :: ji
    INTEGER :: jk
    LOGICAL, DIMENSION(:), ALLOCATABLE :: llvalid
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('obs_surf_compress', 'r0', 0, 0)
    IF (PRESENT(lvalid)) THEN
      insurf = 0
      DO ji = 1, surf % nsurf
        IF (lvalid(ji)) THEN
          insurf = insurf + 1
        END IF
      END DO
    ELSE
      insurf = surf % nsurf
    END IF
    IF (lallocate) THEN
      CALL obs_surf_alloc(newsurf, insurf, surf % nvar, surf % nextra, surf % nstp, surf % npi, surf % npj)
    END IF
    ALLOCATE(llvalid(surf % nsurf))
    CALL profile_psy_data0 % PostEnd
    IF (PRESENT(lvalid)) THEN
      !$ACC KERNELS
      llvalid(:) = lvalid(:)
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      llvalid(:) = .TRUE.
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data1 % PreStart('obs_surf_compress', 'r1', 0, 0)
    insurf = 0
    DO ji = 1, surf % nsurf
      IF (llvalid(ji)) THEN
        insurf = insurf + 1
        newsurf % mi(insurf) = surf % mi(ji)
        newsurf % mj(insurf) = surf % mj(ji)
        newsurf % mt(insurf) = surf % mt(ji)
        newsurf % nsidx(insurf) = surf % nsidx(ji)
        newsurf % nsfil(insurf) = surf % nsfil(ji)
        newsurf % nyea(insurf) = surf % nyea(ji)
        newsurf % nmon(insurf) = surf % nmon(ji)
        newsurf % nday(insurf) = surf % nday(ji)
        newsurf % nhou(insurf) = surf % nhou(ji)
        newsurf % nmin(insurf) = surf % nmin(ji)
        newsurf % mstp(insurf) = surf % mstp(ji)
        newsurf % nqc(insurf) = surf % nqc(ji)
        newsurf % ntyp(insurf) = surf % ntyp(ji)
        newsurf % cwmo(insurf) = surf % cwmo(ji)
        newsurf % rlam(insurf) = surf % rlam(ji)
        newsurf % rphi(insurf) = surf % rphi(ji)
        DO jk = 1, surf % nvar
          newsurf % robs(insurf, jk) = surf % robs(ji, jk)
          newsurf % rmod(insurf, jk) = surf % rmod(ji, jk)
        END DO
        DO jk = 1, surf % nextra
          newsurf % rext(insurf, jk) = surf % rext(ji, jk)
        END DO
        newsurf % nsind(insurf) = ji
      END IF
    END DO
    newsurf % nsurf = insurf
    CALL obs_mpp_sum_integer(newsurf % nsurf, newsurf % nsurfmpp)
    newsurf % nstp = surf % nstp
    newsurf % cvars(:) = surf % cvars(:)
    newsurf % mt(insurf) = surf % mt(ji)
    DEALLOCATE(llvalid)
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE obs_surf_compress
  SUBROUTINE obs_surf_decompress(surf, oldsurf, ldeallocate, kumout)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(obs_surf), INTENT(INOUT) :: surf
    TYPE(obs_surf), INTENT(INOUT) :: oldsurf
    LOGICAL :: ldeallocate
    INTEGER, INTENT(IN) :: kumout
    INTEGER :: ji
    INTEGER :: jj
    INTEGER :: jk
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_surf_decompress', 'r0', 0, 0)
    DO ji = 1, surf % nsurf
      jj = surf % nsind(ji)
      oldsurf % mi(jj) = surf % mi(ji)
      oldsurf % mj(jj) = surf % mj(ji)
      oldsurf % mt(jj) = surf % mt(ji)
      oldsurf % nsidx(jj) = surf % nsidx(ji)
      oldsurf % nsfil(jj) = surf % nsfil(ji)
      oldsurf % nyea(jj) = surf % nyea(ji)
      oldsurf % nmon(jj) = surf % nmon(ji)
      oldsurf % nday(jj) = surf % nday(ji)
      oldsurf % nhou(jj) = surf % nhou(ji)
      oldsurf % nmin(jj) = surf % nmin(ji)
      oldsurf % mstp(jj) = surf % mstp(ji)
      oldsurf % nqc(jj) = surf % nqc(ji)
      oldsurf % ntyp(jj) = surf % ntyp(ji)
      oldsurf % cwmo(jj) = surf % cwmo(ji)
      oldsurf % rlam(jj) = surf % rlam(ji)
      oldsurf % rphi(jj) = surf % rphi(ji)
    END DO
    DO jk = 1, surf % nvar
      DO ji = 1, surf % nsurf
        jj = surf % nsind(ji)
        oldsurf % robs(jj, jk) = surf % robs(ji, jk)
        oldsurf % rmod(jj, jk) = surf % rmod(ji, jk)
      END DO
    END DO
    DO jk = 1, surf % nextra
      DO ji = 1, surf % nsurf
        jj = surf % nsind(ji)
        oldsurf % rext(jj, jk) = surf % rext(ji, jk)
      END DO
    END DO
    IF (ldeallocate) CALL obs_surf_dealloc(surf)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_surf_decompress
END MODULE obs_surf_def