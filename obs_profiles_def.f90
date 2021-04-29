MODULE obs_profiles_def
  USE par_kind, ONLY: wp
  USE in_out_manager
  USE obs_mpp, ONLY: obs_mpp_sum_integers
  USE obs_fbm
  USE lib_mpp, ONLY: ctl_warn, ctl_stop
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: obs_prof, obs_prof_var, obs_prof_valid, obs_prof_alloc, obs_prof_alloc_var, obs_prof_dealloc, obs_prof_compress, obs_prof_decompress, obs_prof_staend
  TYPE :: obs_prof_valid
    LOGICAL, POINTER, DIMENSION(:) :: luse
  END TYPE obs_prof_valid
  TYPE :: obs_prof_var
    INTEGER, POINTER, DIMENSION(:) :: mvk, nvpidx, nvlidx, nvqc, idqc
    REAL(KIND = wp), POINTER, DIMENSION(:) :: vdep, vobs, vmod
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: vext
    INTEGER, POINTER, DIMENSION(:) :: nvind
    INTEGER, POINTER, DIMENSION(:, :) :: idqcf, nvqcf
  END TYPE obs_prof_var
  TYPE :: obs_prof
    INTEGER :: nvar
    INTEGER :: next
    INTEGER :: nprof
    INTEGER :: nstp
    INTEGER :: npi
    INTEGER :: npj
    INTEGER :: npk
    INTEGER :: nprofup
    CHARACTER(LEN = 8), POINTER, DIMENSION(:) :: cvars
    INTEGER, POINTER, DIMENSION(:) :: nvprot, nvprotmpp
    INTEGER, POINTER, DIMENSION(:) :: npidx, npfil, nyea, nmon, nday, nhou, nmin, mstp, nqc, ntyp, ipqc, itqc
    REAL(KIND = wp), POINTER, DIMENSION(:) :: rlam, rphi
    CHARACTER(LEN = 8), POINTER, DIMENSION(:) :: cwmo
    INTEGER, POINTER, DIMENSION(:, :) :: npvsta, npvend, mi, mj, ivqc
    INTEGER, POINTER, DIMENSION(:, :) :: nqcf, ipqcf, itqcf
    INTEGER, POINTER, DIMENSION(:, :, :) :: ivqcf
    TYPE(obs_prof_var), POINTER, DIMENSION(:) :: var
    INTEGER, POINTER, DIMENSION(:) :: npstp, npstpmpp
    INTEGER, POINTER, DIMENSION(:, :) :: nvstp, nvstpmpp
    REAL(KIND = wp), POINTER, DIMENSION(:, :, :, :) :: vdmean
    INTEGER, POINTER, DIMENSION(:) :: npind
  END TYPE obs_prof
  CONTAINS
  SUBROUTINE obs_prof_alloc(prof, kvar, kext, kprof, ko3dt, kstp, kpi, kpj, kpk)
    TYPE(obs_prof), INTENT(INOUT) :: prof
    INTEGER, INTENT(IN) :: kprof
    INTEGER, INTENT(IN) :: kvar
    INTEGER, INTENT(IN) :: kext
    INTEGER, INTENT(IN), DIMENSION(kvar) :: ko3dt
    INTEGER, INTENT(IN) :: kstp
    INTEGER, INTENT(IN) :: kpi
    INTEGER, INTENT(IN) :: kpj
    INTEGER, INTENT(IN) :: kpk
    INTEGER :: jvar
    INTEGER :: ji
    prof % nvar = kvar
    prof % next = kext
    prof % nprof = kprof
    prof % nstp = kstp
    prof % npi = kpi
    prof % npj = kpj
    prof % npk = kpk
    ALLOCATE(prof % cvars(kvar), prof % nvprot(kvar), prof % nvprotmpp(kvar))
    DO jvar = 1, kvar
      prof % cvars(jvar) = "NotSet"
      prof % nvprot(jvar) = ko3dt(jvar)
      prof % nvprotmpp(jvar) = 0
    END DO
    ALLOCATE(prof % npvsta(kprof, kvar), prof % npvend(kprof, kvar), prof % mi(kprof, kvar), prof % mj(kprof, kvar), prof % ivqc(kprof, kvar))
    ALLOCATE(prof % ivqcf(idefnqcf, kprof, kvar))
    ALLOCATE(prof % npidx(kprof), prof % npfil(kprof), prof % nyea(kprof), prof % nmon(kprof), prof % nday(kprof), prof % nhou(kprof), prof % nmin(kprof), prof % mstp(kprof), prof % nqc(kprof), prof % ipqc(kprof), prof % itqc(kprof), prof % ntyp(kprof), prof % rlam(kprof), prof % rphi(kprof), prof % cwmo(kprof), prof % npind(kprof))
    ALLOCATE(prof % nqcf(idefnqcf, kprof), prof % ipqcf(idefnqcf, kprof), prof % itqcf(idefnqcf, kprof))
    ALLOCATE(prof % var(kvar))
    DO jvar = 1, kvar
      IF (ko3dt(jvar) >= 0) THEN
        CALL obs_prof_alloc_var(prof, jvar, kext, ko3dt(jvar))
      END IF
    END DO
    ALLOCATE(prof % npstp(kstp), prof % npstpmpp(kstp))
    ALLOCATE(prof % nvstp(kstp, kvar), prof % nvstpmpp(kstp, kvar))
    ALLOCATE(prof % vdmean(kpi, kpj, kpk, kvar))
    DO ji = 1, kprof
      prof % npind(ji) = ji
    END DO
    DO jvar = 1, kvar
      DO ji = 1, ko3dt(jvar)
        prof % var(jvar) % nvind(ji) = ji
      END DO
    END DO
    prof % npstp(:) = 0
    prof % npstpmpp(:) = 0
    prof % nvstp(:, :) = 0
    prof % nvstpmpp(:, :) = 0
    prof % nprofup = 0
  END SUBROUTINE obs_prof_alloc
  SUBROUTINE obs_prof_dealloc(prof)
    TYPE(obs_prof), INTENT(INOUT) :: prof
    INTEGER :: jvar
    DEALLOCATE(prof % npvsta, prof % npvend)
    DEALLOCATE(prof % mi, prof % mj, prof % ivqc, prof % ivqcf, prof % npidx, prof % npfil, prof % nyea, prof % nmon, prof % nday, prof % nhou, prof % nmin, prof % mstp, prof % nqc, prof % ipqc, prof % itqc, prof % nqcf, prof % ipqcf, prof % itqcf, prof % ntyp, prof % rlam, prof % rphi, prof % cwmo, prof % npind)
    DO jvar = 1, prof % nvar
      IF (prof % nvprot(jvar) >= 0) THEN
        CALL obs_prof_dealloc_var(prof, jvar)
      END IF
    END DO
    DEALLOCATE(prof % var)
    DEALLOCATE(prof % npstp, prof % npstpmpp)
    DEALLOCATE(prof % nvstp, prof % nvstpmpp)
    DEALLOCATE(prof % vdmean)
    DEALLOCATE(prof % cvars, prof % nvprot, prof % nvprotmpp)
  END SUBROUTINE obs_prof_dealloc
  SUBROUTINE obs_prof_alloc_var(prof, kvar, kext, kobs)
    TYPE(obs_prof), INTENT(INOUT) :: prof
    INTEGER, INTENT(IN) :: kvar
    INTEGER, INTENT(IN) :: kext
    INTEGER, INTENT(IN) :: kobs
    ALLOCATE(prof % var(kvar) % mvk(kobs), prof % var(kvar) % nvpidx(kobs), prof % var(kvar) % nvlidx(kobs), prof % var(kvar) % nvqc(kobs), prof % var(kvar) % idqc(kobs), prof % var(kvar) % vdep(kobs), prof % var(kvar) % vobs(kobs), prof % var(kvar) % vmod(kobs), prof % var(kvar) % nvind(kobs))
    ALLOCATE(prof % var(kvar) % idqcf(idefnqcf, kobs), prof % var(kvar) % nvqcf(idefnqcf, kobs))
    IF (kext > 0) THEN
      ALLOCATE(prof % var(kvar) % vext(kobs, kext))
    END IF
  END SUBROUTINE obs_prof_alloc_var
  SUBROUTINE obs_prof_dealloc_var(prof, kvar)
    TYPE(obs_prof), INTENT(INOUT) :: prof
    INTEGER, INTENT(IN) :: kvar
    DEALLOCATE(prof % var(kvar) % mvk, prof % var(kvar) % nvpidx, prof % var(kvar) % nvlidx, prof % var(kvar) % nvqc, prof % var(kvar) % idqc, prof % var(kvar) % vdep, prof % var(kvar) % vobs, prof % var(kvar) % vmod, prof % var(kvar) % nvind, prof % var(kvar) % idqcf, prof % var(kvar) % nvqcf)
    IF (prof % next > 0) THEN
      DEALLOCATE(prof % var(kvar) % vext)
    END IF
  END SUBROUTINE obs_prof_dealloc_var
  SUBROUTINE obs_prof_compress(prof, newprof, lallocate, kumout, lvalid, lvvalid)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(obs_prof), INTENT(IN) :: prof
    TYPE(obs_prof), INTENT(INOUT) :: newprof
    LOGICAL :: lallocate
    INTEGER, INTENT(IN) :: kumout
    TYPE(obs_prof_valid), OPTIONAL, INTENT(IN) :: lvalid
    TYPE(obs_prof_valid), OPTIONAL, INTENT(IN), DIMENSION(prof % nvar) :: lvvalid
    INTEGER :: inprof
    INTEGER, DIMENSION(prof % nvar) :: invpro
    INTEGER :: jvar
    INTEGER :: jext
    INTEGER :: ji
    INTEGER :: jj
    LOGICAL :: lfirst
    TYPE(obs_prof_valid) :: llvalid
    TYPE(obs_prof_valid), DIMENSION(prof % nvar) :: llvvalid
    LOGICAL :: lallpresent
    LOGICAL :: lnonepresent
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_prof_compress', 'r0', 0, 0)
    lallpresent = .FALSE.
    lnonepresent = .FALSE.
    IF (PRESENT(lvalid) .AND. PRESENT(lvvalid)) THEN
      lallpresent = .TRUE.
    ELSE IF ((.NOT. PRESENT(lvalid)) .AND. (.NOT. PRESENT(lvvalid))) THEN
      lnonepresent = .TRUE.
    ELSE
      CALL ctl_stop('Error in obs_prof_compress:', 'Either all selection variables should be set', 'or no selection variable should be set')
    END IF
    IF (lallpresent) THEN
      inprof = 0
      invpro(:) = 0
      DO ji = 1, prof % nprof
        IF (lvalid % luse(ji)) THEN
          inprof = inprof + 1
          DO jvar = 1, prof % nvar
            DO jj = prof % npvsta(ji, jvar), prof % npvend(ji, jvar)
              IF (lvvalid(jvar) % luse(jj)) invpro(jvar) = invpro(jvar) + 1
            END DO
          END DO
        END IF
      END DO
    ELSE
      inprof = prof % nprof
      invpro(:) = prof % nvprot(:)
    END IF
    IF (lallocate) THEN
      CALL obs_prof_alloc(newprof, prof % nvar, prof % next, inprof, invpro, prof % nstp, prof % npi, prof % npj, prof % npk)
    END IF
    ALLOCATE(llvalid % luse(prof % nprof))
    DO jvar = 1, prof % nvar
      ALLOCATE(llvvalid(jvar) % luse(prof % nvprot(jvar)))
    END DO
    IF (lallpresent) THEN
      llvalid % luse(:) = lvalid % luse(:)
      DO jvar = 1, prof % nvar
        llvvalid(jvar) % luse(:) = lvvalid(jvar) % luse(:)
      END DO
    ELSE
      llvalid % luse(:) = .TRUE.
      DO jvar = 1, prof % nvar
        llvvalid(jvar) % luse(:) = .TRUE.
      END DO
    END IF
    inprof = 0
    invpro(:) = 0
    newprof % npvsta(:, :) = 0
    newprof % npvend(:, :) = - 1
    DO ji = 1, prof % nprof
      IF (llvalid % luse(ji)) THEN
        inprof = inprof + 1
        newprof % mi(inprof, :) = prof % mi(ji, :)
        newprof % mj(inprof, :) = prof % mj(ji, :)
        newprof % npidx(inprof) = prof % npidx(ji)
        newprof % npfil(inprof) = prof % npfil(ji)
        newprof % nyea(inprof) = prof % nyea(ji)
        newprof % nmon(inprof) = prof % nmon(ji)
        newprof % nday(inprof) = prof % nday(ji)
        newprof % nhou(inprof) = prof % nhou(ji)
        newprof % nmin(inprof) = prof % nmin(ji)
        newprof % mstp(inprof) = prof % mstp(ji)
        newprof % nqc(inprof) = prof % nqc(ji)
        newprof % ipqc(inprof) = prof % ipqc(ji)
        newprof % itqc(inprof) = prof % itqc(ji)
        newprof % ivqc(inprof, :) = prof % ivqc(ji, :)
        newprof % ntyp(inprof) = prof % ntyp(ji)
        newprof % rlam(inprof) = prof % rlam(ji)
        newprof % rphi(inprof) = prof % rphi(ji)
        newprof % cwmo(inprof) = prof % cwmo(ji)
        newprof % nqcf(:, inprof) = prof % nqcf(:, ji)
        newprof % ipqcf(:, inprof) = prof % ipqcf(:, ji)
        newprof % itqcf(:, inprof) = prof % itqcf(:, ji)
        newprof % ivqcf(:, inprof, :) = prof % ivqcf(:, ji, :)
        newprof % npind(inprof) = ji
        DO jvar = 1, prof % nvar
          lfirst = .TRUE.
          DO jj = prof % npvsta(ji, jvar), prof % npvend(ji, jvar)
            IF (llvvalid(jvar) % luse(jj)) THEN
              invpro(jvar) = invpro(jvar) + 1
              IF (lfirst) THEN
                lfirst = .FALSE.
                newprof % npvsta(inprof, jvar) = invpro(jvar)
              END IF
              newprof % npvend(inprof, jvar) = invpro(jvar)
              newprof % var(jvar) % mvk(invpro(jvar)) = prof % var(jvar) % mvk(jj)
              newprof % var(jvar) % nvpidx(invpro(jvar)) = prof % var(jvar) % nvpidx(jj)
              newprof % var(jvar) % nvlidx(invpro(jvar)) = prof % var(jvar) % nvlidx(jj)
              newprof % var(jvar) % nvqc(invpro(jvar)) = prof % var(jvar) % nvqc(jj)
              newprof % var(jvar) % idqc(invpro(jvar)) = prof % var(jvar) % idqc(jj)
              newprof % var(jvar) % idqcf(:, invpro(jvar)) = prof % var(jvar) % idqcf(:, jj)
              newprof % var(jvar) % nvqcf(:, invpro(jvar)) = prof % var(jvar) % nvqcf(:, jj)
              newprof % var(jvar) % vdep(invpro(jvar)) = prof % var(jvar) % vdep(jj)
              newprof % var(jvar) % vobs(invpro(jvar)) = prof % var(jvar) % vobs(jj)
              newprof % var(jvar) % vmod(invpro(jvar)) = prof % var(jvar) % vmod(jj)
              DO jext = 1, prof % next
                newprof % var(jvar) % vext(invpro(jvar), jext) = prof % var(jvar) % vext(jj, jext)
              END DO
              newprof % var(jvar) % nvind(invpro(jvar)) = jj
            END IF
          END DO
        END DO
      END IF
    END DO
    DO jvar = 1, prof % nvar
      newprof % nvprot(jvar) = invpro(jvar)
    END DO
    CALL obs_mpp_sum_integers(newprof % nvprot, newprof % nvprotmpp, prof % nvar)
    newprof % nvar = prof % nvar
    newprof % next = prof % next
    newprof % nstp = prof % nstp
    newprof % npi = prof % npi
    newprof % npj = prof % npj
    newprof % npk = prof % npk
    newprof % cvars(:) = prof % cvars(:)
    DO jvar = 1, prof % nvar
      DEALLOCATE(llvvalid(jvar) % luse)
    END DO
    DEALLOCATE(llvalid % luse)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_prof_compress
  SUBROUTINE obs_prof_decompress(prof, oldprof, ldeallocate, kumout)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(obs_prof), INTENT(INOUT) :: prof
    TYPE(obs_prof), INTENT(INOUT) :: oldprof
    LOGICAL :: ldeallocate
    INTEGER, INTENT(IN) :: kumout
    INTEGER :: jvar
    INTEGER :: jext
    INTEGER :: ji
    INTEGER :: jj
    INTEGER :: jk
    INTEGER :: jl
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_prof_decompress', 'r0', 0, 0)
    DO ji = 1, prof % nprof
      jk = prof % npind(ji)
      oldprof % mi(jk, :) = prof % mi(ji, :)
      oldprof % mj(jk, :) = prof % mj(ji, :)
      oldprof % npidx(jk) = prof % npidx(ji)
      oldprof % npfil(jk) = prof % npfil(ji)
      oldprof % nyea(jk) = prof % nyea(ji)
      oldprof % nmon(jk) = prof % nmon(ji)
      oldprof % nday(jk) = prof % nday(ji)
      oldprof % nhou(jk) = prof % nhou(ji)
      oldprof % nmin(jk) = prof % nmin(ji)
      oldprof % mstp(jk) = prof % mstp(ji)
      oldprof % nqc(jk) = prof % nqc(ji)
      oldprof % ipqc(jk) = prof % ipqc(ji)
      oldprof % itqc(jk) = prof % itqc(ji)
      oldprof % ivqc(jk, :) = prof % ivqc(ji, :)
      oldprof % ntyp(jk) = prof % ntyp(ji)
      oldprof % rlam(jk) = prof % rlam(ji)
      oldprof % rphi(jk) = prof % rphi(ji)
      oldprof % cwmo(jk) = prof % cwmo(ji)
      oldprof % nqcf(:, jk) = prof % nqcf(:, ji)
      oldprof % ipqcf(:, jk) = prof % ipqcf(:, ji)
      oldprof % itqcf(:, jk) = prof % itqcf(:, ji)
      oldprof % ivqcf(:, jk, :) = prof % ivqcf(:, ji, :)
      DO jvar = 1, prof % nvar
        DO jj = prof % npvsta(ji, jvar), prof % npvend(ji, jvar)
          jl = prof % var(jvar) % nvind(jj)
          oldprof % var(jvar) % mvk(jl) = prof % var(jvar) % mvk(jj)
          oldprof % var(jvar) % nvpidx(jl) = prof % var(jvar) % nvpidx(jj)
          oldprof % var(jvar) % nvlidx(jl) = prof % var(jvar) % nvlidx(jj)
          oldprof % var(jvar) % nvqc(jl) = prof % var(jvar) % nvqc(jj)
          oldprof % var(jvar) % idqc(jl) = prof % var(jvar) % idqc(jj)
          oldprof % var(jvar) % vdep(jl) = prof % var(jvar) % vdep(jj)
          oldprof % var(jvar) % vobs(jl) = prof % var(jvar) % vobs(jj)
          oldprof % var(jvar) % vmod(jl) = prof % var(jvar) % vmod(jj)
          oldprof % var(jvar) % idqcf(:, jl) = prof % var(jvar) % idqcf(:, jj)
          oldprof % var(jvar) % nvqcf(:, jl) = prof % var(jvar) % nvqcf(:, jj)
          DO jext = 1, prof % next
            oldprof % var(jvar) % vext(jl, jext) = prof % var(jvar) % vext(jj, jext)
          END DO
        END DO
      END DO
    END DO
    IF (ldeallocate) CALL obs_prof_dealloc(prof)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_prof_decompress
  SUBROUTINE obs_prof_staend(prof, kvarno)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(obs_prof), INTENT(INOUT) :: prof
    INTEGER, INTENT(IN) :: kvarno
    INTEGER :: ji
    INTEGER :: iprofno
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_prof_staend', 'r0', 0, 0)
    prof % npvsta(:, kvarno) = prof % nvprot(kvarno) + 1
    prof % npvend(:, kvarno) = - 1
    DO ji = 1, prof % nvprot(kvarno)
      iprofno = prof % var(kvarno) % nvpidx(ji)
      prof % npvsta(iprofno, kvarno) = MIN(ji, prof % npvsta(iprofno, kvarno))
      prof % npvend(iprofno, kvarno) = MAX(ji, prof % npvend(iprofno, kvarno))
    END DO
    DO ji = 1, prof % nprof
      IF (prof % npvsta(ji, kvarno) == (prof % nvprot(kvarno) + 1)) prof % npvsta(ji, kvarno) = 0
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_prof_staend
END MODULE obs_profiles_def