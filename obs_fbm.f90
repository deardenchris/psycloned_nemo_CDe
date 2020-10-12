MODULE obs_fbm
  USE netcdf
  USE obs_utils
  IMPLICIT NONE
  PUBLIC
  INTEGER, PARAMETER :: fbsp = SELECTED_REAL_KIND(6, 37)
  INTEGER, PARAMETER :: fbdp = SELECTED_REAL_KIND(12, 307)
  INTEGER, PARAMETER :: ilenwmo = 8
  INTEGER, PARAMETER :: ilentyp = 4
  INTEGER, PARAMETER :: ilenname = 8
  INTEGER, PARAMETER :: ilengrid = 1
  INTEGER, PARAMETER :: ilenjuld = 14
  INTEGER, PARAMETER :: idefnqcf = 2
  INTEGER, PARAMETER :: ilenlong = 128
  INTEGER, PARAMETER :: ilenunit = 32
  INTEGER, PARAMETER :: fbimdi = - 99999
  REAL(KIND = fbsp), PARAMETER :: fbrmdi = 99999
  TYPE :: obfbdata
    LOGICAL :: lalloc
    LOGICAL :: lgrid
    INTEGER :: nvar
    INTEGER :: nobs
    INTEGER :: nlev
    INTEGER :: nadd
    INTEGER :: next
    INTEGER :: nqcf
    CHARACTER(LEN = ilenwmo), DIMENSION(:), POINTER :: cdwmo
    CHARACTER(LEN = ilentyp), DIMENSION(:), POINTER :: cdtyp
    CHARACTER(LEN = ilenjuld) :: cdjuldref
    INTEGER, DIMENSION(:), POINTER :: kindex
    INTEGER, DIMENSION(:), POINTER :: ioqc, ipqc, itqc
    INTEGER, DIMENSION(:, :), POINTER :: ioqcf, ipqcf, itqcf
    INTEGER, DIMENSION(:, :), POINTER :: idqc
    INTEGER, DIMENSION(:, :, :), POINTER :: idqcf
    REAL(KIND = fbdp), DIMENSION(:), POINTER :: plam, pphi, ptim
    REAL(KIND = fbsp), DIMENSION(:, :), POINTER :: pdep
    CHARACTER(LEN = ilenname), DIMENSION(:), POINTER :: cname
    REAL(KIND = fbsp), DIMENSION(:, :, :), POINTER :: pob
    CHARACTER(LEN = ilenlong), DIMENSION(:), POINTER :: coblong
    CHARACTER(LEN = ilenunit), DIMENSION(:), POINTER :: cobunit
    INTEGER, DIMENSION(:, :), POINTER :: ivqc
    INTEGER, DIMENSION(:, :, :), POINTER :: ivqcf
    INTEGER, DIMENSION(:, :, :), POINTER :: ivlqc
    INTEGER, DIMENSION(:, :, :, :), POINTER :: ivlqcf
    INTEGER, DIMENSION(:, :), POINTER :: iproc, iobsi, iobsj
    INTEGER, DIMENSION(:, :, :), POINTER :: iobsk
    CHARACTER(LEN = ilengrid), DIMENSION(:), POINTER :: cgrid
    CHARACTER(LEN = ilenname), DIMENSION(:), POINTER :: caddname
    CHARACTER(LEN = ilenlong), DIMENSION(:, :), POINTER :: caddlong
    CHARACTER(LEN = ilenunit), DIMENSION(:, :), POINTER :: caddunit
    REAL(KIND = fbsp), DIMENSION(:, :, :, :), POINTER :: padd
    CHARACTER(LEN = ilenname), DIMENSION(:), POINTER :: cextname
    CHARACTER(LEN = ilenlong), DIMENSION(:), POINTER :: cextlong
    CHARACTER(LEN = ilenunit), DIMENSION(:), POINTER :: cextunit
    REAL(KIND = fbsp), DIMENSION(:, :, :), POINTER :: pext
  END TYPE obfbdata
  PRIVATE :: putvaratt_obfbdata
  CONTAINS
  SUBROUTINE init_obfbdata(fbdata)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(obfbdata) :: fbdata
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('init_obfbdata', 'r0', 0, 0)
    fbdata % nvar = 0
    fbdata % nobs = 0
    fbdata % nlev = 0
    fbdata % nadd = 0
    fbdata % next = 0
    fbdata % nqcf = idefnqcf
    fbdata % lalloc = .FALSE.
    fbdata % lgrid = .FALSE.
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE init_obfbdata
  SUBROUTINE alloc_obfbdata(fbdata, kvar, kobs, klev, kadd, kext, lgrid, kqcf)
    TYPE(obfbdata) :: fbdata
    INTEGER, INTENT(IN) :: kvar
    INTEGER, INTENT(IN) :: kobs
    INTEGER, INTENT(IN) :: klev
    INTEGER, INTENT(IN) :: kadd
    INTEGER, INTENT(IN) :: kext
    LOGICAL, INTENT(IN) :: lgrid
    INTEGER, OPTIONAL :: kqcf
    INTEGER :: ji
    INTEGER :: jv
    IF (fbdata % lalloc) THEN
      CALL dealloc_obfbdata(fbdata)
    END IF
    fbdata % lalloc = .TRUE.
    fbdata % nvar = kvar
    fbdata % nobs = kobs
    fbdata % nlev = MAX(klev, 1)
    fbdata % nadd = kadd
    fbdata % next = kext
    IF (PRESENT(kqcf)) THEN
      fbdata % nqcf = kqcf
    ELSE
      fbdata % nqcf = idefnqcf
    END IF
    fbdata % cdjuldref = REPEAT('X', ilenjuld)
    ALLOCATE(fbdata % cname(fbdata % nvar), fbdata % coblong(fbdata % nvar), fbdata % cobunit(fbdata % nvar))
    DO ji = 1, fbdata % nvar
      WRITE(fbdata % cname(ji), FMT = '(A,I2.2)') 'V_', ji
      fbdata % coblong(ji) = REPEAT(' ', ilenlong)
      fbdata % cobunit(ji) = REPEAT(' ', ilenunit)
    END DO
    IF (lgrid) THEN
      ALLOCATE(fbdata % cgrid(fbdata % nvar))
      fbdata % cgrid(:) = REPEAT('X', ilengrid)
      fbdata % lgrid = .TRUE.
    END IF
    IF (fbdata % nadd > 0) THEN
      ALLOCATE(fbdata % caddname(fbdata % nadd), fbdata % caddlong(fbdata % nadd, fbdata % nvar), fbdata % caddunit(fbdata % nadd, &
&fbdata % nvar))
      DO ji = 1, fbdata % nadd
        WRITE(fbdata % caddname(ji), FMT = '(A,I2.2)') 'A', ji
      END DO
      DO jv = 1, fbdata % nvar
        DO ji = 1, fbdata % nadd
          fbdata % caddlong(ji, jv) = REPEAT(' ', ilenlong)
          fbdata % caddunit(ji, jv) = REPEAT(' ', ilenunit)
        END DO
      END DO
    END IF
    IF (fbdata % next > 0) THEN
      ALLOCATE(fbdata % cextname(fbdata % next), fbdata % cextlong(fbdata % next), fbdata % cextunit(fbdata % next))
      DO ji = 1, fbdata % next
        WRITE(fbdata % cextname(ji), FMT = '(A,I2.2)') 'E_', ji
        fbdata % cextlong(ji) = REPEAT(' ', ilenlong)
        fbdata % cextunit(ji) = REPEAT(' ', ilenunit)
      END DO
    END IF
    IF (fbdata % nobs > 0) THEN
      ALLOCATE(fbdata % cdwmo(fbdata % nobs), fbdata % cdtyp(fbdata % nobs), fbdata % ioqc(fbdata % nobs), fbdata % ioqcf(fbdata % &
&nqcf, fbdata % nobs), fbdata % ipqc(fbdata % nobs), fbdata % ipqcf(fbdata % nqcf, fbdata % nobs), fbdata % itqc(fbdata % nobs), &
&fbdata % itqcf(fbdata % nqcf, fbdata % nobs), fbdata % idqc(fbdata % nlev, fbdata % nobs), fbdata % idqcf(fbdata % nqcf, fbdata % &
&nlev, fbdata % nobs), fbdata % plam(fbdata % nobs), fbdata % pphi(fbdata % nobs), fbdata % pdep(fbdata % nlev, fbdata % nobs), &
&fbdata % ptim(fbdata % nobs), fbdata % kindex(fbdata % nobs), fbdata % ivqc(fbdata % nobs, fbdata % nvar), fbdata % ivqcf(fbdata &
&% nqcf, fbdata % nobs, fbdata % nvar), fbdata % ivlqc(fbdata % nlev, fbdata % nobs, fbdata % nvar), fbdata % ivlqcf(fbdata % &
&nqcf, fbdata % nlev, fbdata % nobs, fbdata % nvar), fbdata % pob(fbdata % nlev, fbdata % nobs, fbdata % nvar))
      fbdata % kindex(:) = fbimdi
      fbdata % cdwmo(:) = REPEAT('X', ilenwmo)
      fbdata % cdtyp(:) = REPEAT('X', ilentyp)
      fbdata % ioqc(:) = fbimdi
      fbdata % ioqcf(:, :) = fbimdi
      fbdata % ipqc(:) = fbimdi
      fbdata % ipqcf(:, :) = fbimdi
      fbdata % itqc(:) = fbimdi
      fbdata % itqcf(:, :) = fbimdi
      fbdata % idqc(:, :) = fbimdi
      fbdata % idqcf(:, :, :) = fbimdi
      fbdata % plam(:) = fbrmdi
      fbdata % pphi(:) = fbrmdi
      fbdata % pdep(:, :) = fbrmdi
      fbdata % ptim(:) = fbrmdi
      fbdata % ivqc(:, :) = fbimdi
      fbdata % ivqcf(:, :, :) = fbimdi
      fbdata % ivlqc(:, :, :) = fbimdi
      fbdata % ivlqcf(:, :, :, :) = fbimdi
      fbdata % pob(:, :, :) = fbrmdi
      IF (lgrid) THEN
        ALLOCATE(fbdata % iproc(fbdata % nobs, fbdata % nvar), fbdata % iobsi(fbdata % nobs, fbdata % nvar), fbdata % iobsj(fbdata &
&% nobs, fbdata % nvar), fbdata % iobsk(fbdata % nlev, fbdata % nobs, fbdata % nvar))
        fbdata % iproc(:, :) = fbimdi
        fbdata % iobsi(:, :) = fbimdi
        fbdata % iobsj(:, :) = fbimdi
        fbdata % iobsk(:, :, :) = fbimdi
        fbdata % lgrid = .TRUE.
      END IF
      IF (fbdata % nadd > 0) THEN
        ALLOCATE(fbdata % padd(fbdata % nlev, fbdata % nobs, fbdata % nadd, fbdata % nvar))
        fbdata % padd(:, :, :, :) = fbrmdi
      END IF
      IF (fbdata % next > 0) THEN
        ALLOCATE(fbdata % pext(fbdata % nlev, fbdata % nobs, fbdata % next))
        fbdata % pext(:, :, :) = fbrmdi
      END IF
    END IF
  END SUBROUTINE alloc_obfbdata
  SUBROUTINE dealloc_obfbdata(fbdata)
    TYPE(obfbdata) :: fbdata
    DEALLOCATE(fbdata % cname, fbdata % coblong, fbdata % cobunit)
    IF (fbdata % lgrid) THEN
      DEALLOCATE(fbdata % cgrid)
    END IF
    IF (fbdata % nadd > 0) THEN
      DEALLOCATE(fbdata % caddname, fbdata % caddlong, fbdata % caddunit)
    END IF
    IF (fbdata % next > 0) THEN
      DEALLOCATE(fbdata % cextname, fbdata % cextlong, fbdata % cextunit)
    END IF
    IF (fbdata % nobs > 0) THEN
      DEALLOCATE(fbdata % cdwmo, fbdata % cdtyp, fbdata % ioqc, fbdata % ioqcf, fbdata % ipqc, fbdata % ipqcf, fbdata % itqc, &
&fbdata % itqcf, fbdata % idqc, fbdata % idqcf, fbdata % plam, fbdata % pphi, fbdata % pdep, fbdata % ptim, fbdata % kindex, &
&fbdata % ivqc, fbdata % ivqcf, fbdata % ivlqc, fbdata % ivlqcf, fbdata % pob)
      IF (fbdata % lgrid) THEN
        DEALLOCATE(fbdata % iproc, fbdata % iobsi, fbdata % iobsj, fbdata % iobsk)
      END IF
      IF (fbdata % nadd > 0) THEN
        DEALLOCATE(fbdata % padd)
      END IF
      IF (fbdata % next > 0) THEN
        DEALLOCATE(fbdata % pext)
      END IF
    END IF
    fbdata % lalloc = .FALSE.
    fbdata % lgrid = .FALSE.
    fbdata % nvar = 0
    fbdata % nobs = 0
    fbdata % nlev = 0
    fbdata % nadd = 0
    fbdata % next = 0
  END SUBROUTINE dealloc_obfbdata
  SUBROUTINE copy_obfbdata(fbdata1, fbdata2, kadd, kext, lgrid, kqcf)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(obfbdata) :: fbdata1
    TYPE(obfbdata) :: fbdata2
    INTEGER, INTENT(IN), OPTIONAL :: kadd
    INTEGER, INTENT(IN), OPTIONAL :: kext
    INTEGER, INTENT(IN), OPTIONAL :: kqcf
    LOGICAL, OPTIONAL :: lgrid
    INTEGER :: nadd
    INTEGER :: next
    INTEGER :: nqcf
    LOGICAL :: llgrid
    INTEGER :: jv
    INTEGER :: je
    INTEGER :: ji
    INTEGER :: jk
    INTEGER :: jq
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('copy_obfbdata', 'r0', 0, 0)
    IF (.NOT. fbdata1 % lalloc) THEN
      CALL fatal_error('copy_obfbdata: input data not allocated', 511)
    END IF
    IF (PRESENT(kadd)) THEN
      nadd = kadd
      IF (nadd < fbdata1 % nadd) THEN
        CALL warning('copy_obfbdata: ' // 'nadd smaller than input nadd', 521)
      END IF
    ELSE
      nadd = fbdata1 % nadd
    END IF
    IF (PRESENT(kext)) THEN
      next = kext
      IF (next < fbdata1 % next) THEN
        CALL fatal_error('copy_obfbdata: ' // 'next smaller than input next', 530)
      END IF
    ELSE
      next = fbdata1 % next
    END IF
    IF (PRESENT(lgrid)) THEN
      llgrid = lgrid
      IF (fbdata1 % lgrid .AND. (.NOT. llgrid)) THEN
        CALL fatal_error('copy_obfbdata: ' // 'switching off grid info not possible', 540)
      END IF
    ELSE
      llgrid = fbdata1 % lgrid
    END IF
    IF (PRESENT(kqcf)) THEN
      nqcf = kqcf
      IF (nqcf < fbdata1 % nqcf) THEN
        CALL fatal_error('copy_obfbdata: ' // 'nqcf smaller than input nqcf', 549)
      END IF
    ELSE
      nqcf = fbdata1 % nqcf
    END IF
    IF (fbdata2 % lalloc) THEN
      IF (fbdata1 % nvar > fbdata2 % nvar) THEN
        CALL fatal_error('copy_obfbdata: ' // 'output kvar smaller than input kvar', 562)
      END IF
      IF (fbdata1 % nobs > fbdata2 % nobs) THEN
        CALL fatal_error('copy_obfbdata: ' // 'output kobs smaller than input kobs', 566)
      END IF
      IF (fbdata1 % nlev > fbdata2 % nlev) THEN
        CALL fatal_error('copy_obfbdata: ' // 'output klev smaller than input klev', 570)
      END IF
      IF (fbdata1 % nadd > fbdata2 % nadd) THEN
        CALL warning('copy_obfbdata: ' // 'output nadd smaller than input nadd', 574)
      END IF
      IF (fbdata1 % next > fbdata2 % next) THEN
        CALL fatal_error('copy_obfbdata: ' // 'output next smaller than input next', 578)
      END IF
      IF (fbdata1 % lgrid .NEQV. fbdata2 % lgrid) THEN
        CALL fatal_error('copy_obfbdata: ' // 'lgrid inconsistent', 582)
      END IF
      IF (fbdata1 % next > fbdata2 % next) THEN
        CALL fatal_error('copy_obfbdata: ' // 'output next smaller than input next', 586)
      END IF
      IF (fbdata1 % nqcf > fbdata2 % nqcf) THEN
        CALL fatal_error('copy_obfbdata: ' // 'output  smaller than input kext', 590)
      END IF
    ELSE
      CALL alloc_obfbdata(fbdata2, fbdata1 % nvar, fbdata1 % nobs, fbdata1 % nlev, nadd, next, llgrid, kqcf = nqcf)
    END IF
    fbdata2 % cdjuldref = fbdata1 % cdjuldref
    DO ji = 1, fbdata1 % nobs
      fbdata2 % cdwmo(ji) = fbdata1 % cdwmo(ji)
      fbdata2 % cdtyp(ji) = fbdata1 % cdtyp(ji)
      fbdata2 % ioqc(ji) = fbdata1 % ioqc(ji)
      fbdata2 % ipqc(ji) = fbdata1 % ipqc(ji)
      fbdata2 % itqc(ji) = fbdata1 % itqc(ji)
      fbdata2 % plam(ji) = fbdata1 % plam(ji)
      fbdata2 % pphi(ji) = fbdata1 % pphi(ji)
      fbdata2 % ptim(ji) = fbdata1 % ptim(ji)
      fbdata2 % kindex(ji) = fbdata1 % kindex(ji)
      DO jq = 1, fbdata1 % nqcf
        fbdata2 % ioqcf(jq, ji) = fbdata1 % ioqcf(jq, ji)
        fbdata2 % ipqcf(jq, ji) = fbdata1 % ipqcf(jq, ji)
        fbdata2 % itqcf(jq, ji) = fbdata1 % itqcf(jq, ji)
      END DO
      DO jk = 1, fbdata1 % nlev
        fbdata2 % idqc(jk, ji) = fbdata1 % idqc(jk, ji)
        fbdata2 % pdep(jk, ji) = fbdata1 % pdep(jk, ji)
        DO jq = 1, fbdata1 % nqcf
          fbdata2 % idqcf(jq, jk, ji) = fbdata1 % idqcf(jq, jk, ji)
        END DO
      END DO
    END DO
    DO jv = 1, fbdata1 % nvar
      fbdata2 % cname(jv) = fbdata1 % cname(jv)
      fbdata2 % coblong(jv) = fbdata1 % coblong(jv)
      fbdata2 % cobunit(jv) = fbdata1 % cobunit(jv)
      DO ji = 1, fbdata1 % nobs
        fbdata2 % ivqc(ji, jv) = fbdata1 % ivqc(ji, jv)
        DO jq = 1, fbdata1 % nqcf
          fbdata2 % ivqcf(jq, ji, jv) = fbdata1 % ivqcf(jq, ji, jv)
        END DO
        DO jk = 1, fbdata1 % nlev
          fbdata2 % ivlqc(jk, ji, jv) = fbdata1 % ivlqc(jk, ji, jv)
          fbdata2 % pob(jk, ji, jv) = fbdata1 % pob(jk, ji, jv)
          DO jq = 1, fbdata1 % nqcf
            fbdata2 % ivlqcf(jq, jk, ji, jv) = fbdata1 % ivlqcf(jq, jk, ji, jv)
          END DO
        END DO
      END DO
    END DO
    IF (fbdata1 % lgrid) THEN
      DO jv = 1, fbdata1 % nvar
        fbdata2 % cgrid(jv) = fbdata1 % cgrid(jv)
        DO ji = 1, fbdata1 % nobs
          fbdata2 % iproc(ji, jv) = fbdata1 % iproc(ji, jv)
          fbdata2 % iobsi(ji, jv) = fbdata1 % iobsi(ji, jv)
          fbdata2 % iobsj(ji, jv) = fbdata1 % iobsj(ji, jv)
          DO jk = 1, fbdata1 % nlev
            fbdata2 % iobsk(jk, ji, jv) = fbdata1 % iobsk(jk, ji, jv)
          END DO
        END DO
      END DO
    END IF
    DO je = 1, MIN(fbdata1 % nadd, fbdata2 % nadd)
      fbdata2 % caddname(je) = fbdata1 % caddname(je)
    END DO
    DO jv = 1, fbdata1 % nvar
      DO je = 1, MIN(fbdata1 % nadd, fbdata2 % nadd)
        fbdata2 % caddlong(je, jv) = fbdata1 % caddlong(je, jv)
        fbdata2 % caddunit(je, jv) = fbdata1 % caddunit(je, jv)
        DO ji = 1, fbdata1 % nobs
          DO jk = 1, fbdata1 % nlev
            fbdata2 % padd(jk, ji, je, jv) = fbdata1 % padd(jk, ji, je, jv)
          END DO
        END DO
      END DO
    END DO
    DO je = 1, fbdata1 % next
      fbdata2 % cextname(je) = fbdata1 % cextname(je)
      fbdata2 % cextlong(je) = fbdata1 % cextlong(je)
      fbdata2 % cextunit(je) = fbdata1 % cextunit(je)
    END DO
    DO je = 1, fbdata1 % next
      DO ji = 1, fbdata1 % nobs
        DO jk = 1, fbdata1 % nlev
          fbdata2 % pext(jk, ji, je) = fbdata1 % pext(jk, ji, je)
        END DO
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE copy_obfbdata
  SUBROUTINE subsamp_obfbdata(fbdata1, fbdata2, llvalid)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(obfbdata) :: fbdata1
    TYPE(obfbdata) :: fbdata2
    LOGICAL, DIMENSION(fbdata1 % nobs) :: llvalid
    INTEGER :: nobs
    INTEGER :: jv
    INTEGER :: je
    INTEGER :: ji
    INTEGER :: jk
    INTEGER :: jq
    INTEGER :: ij
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('subsamp_obfbdata', 'r0', 0, 0)
    IF (.NOT. fbdata1 % lalloc) THEN
      CALL fatal_error('copy_obfbdata: input data not allocated', 726)
    END IF
    IF (fbdata2 % lalloc) THEN
      CALL fatal_error('subsample_obfbdata: ' // 'fbdata2 already allocated', 733)
    END IF
    nobs = COUNT(llvalid)
    CALL alloc_obfbdata(fbdata2, fbdata1 % nvar, nobs, fbdata1 % nlev, fbdata1 % nadd, fbdata1 % next, fbdata1 % lgrid, &
&kqcf = fbdata1 % nqcf)
    fbdata2 % cdjuldref = fbdata1 % cdjuldref
    ij = 0
    DO ji = 1, fbdata1 % nobs
      IF (llvalid(ji)) THEN
        ij = ij + 1
        fbdata2 % cdwmo(ij) = fbdata1 % cdwmo(ji)
        fbdata2 % cdtyp(ij) = fbdata1 % cdtyp(ji)
        fbdata2 % ioqc(ij) = fbdata1 % ioqc(ji)
        fbdata2 % ipqc(ij) = fbdata1 % ipqc(ji)
        fbdata2 % itqc(ij) = fbdata1 % itqc(ji)
        fbdata2 % plam(ij) = fbdata1 % plam(ji)
        fbdata2 % pphi(ij) = fbdata1 % pphi(ji)
        fbdata2 % ptim(ij) = fbdata1 % ptim(ji)
        fbdata2 % kindex(ij) = fbdata1 % kindex(ji)
        DO jq = 1, fbdata1 % nqcf
          fbdata2 % ioqcf(jq, ij) = fbdata1 % ioqcf(jq, ji)
          fbdata2 % ipqcf(jq, ij) = fbdata1 % ipqcf(jq, ji)
          fbdata2 % itqcf(jq, ij) = fbdata1 % itqcf(jq, ji)
        END DO
        DO jk = 1, fbdata1 % nlev
          fbdata2 % idqc(jk, ij) = fbdata1 % idqc(jk, ji)
          fbdata2 % pdep(jk, ij) = fbdata1 % pdep(jk, ji)
          DO jq = 1, fbdata1 % nqcf
            fbdata2 % idqcf(jq, jk, ij) = fbdata1 % idqcf(jq, jk, ji)
          END DO
        END DO
      END IF
    END DO
    DO jv = 1, fbdata1 % nvar
      fbdata2 % cname(jv) = fbdata1 % cname(jv)
      fbdata2 % coblong(jv) = fbdata1 % coblong(jv)
      fbdata2 % cobunit(jv) = fbdata1 % cobunit(jv)
      ij = 0
      DO ji = 1, fbdata1 % nobs
        IF (llvalid(ji)) THEN
          ij = ij + 1
          fbdata2 % ivqc(ij, jv) = fbdata1 % ivqc(ji, jv)
          DO jq = 1, fbdata1 % nqcf
            fbdata2 % ivqcf(jq, ij, jv) = fbdata1 % ivqcf(jq, ji, jv)
          END DO
          DO jk = 1, fbdata1 % nlev
            fbdata2 % ivlqc(jk, ij, jv) = fbdata1 % ivlqc(jk, ji, jv)
            fbdata2 % pob(jk, ij, jv) = fbdata1 % pob(jk, ji, jv)
            DO jq = 1, fbdata1 % nqcf
              fbdata2 % ivlqcf(jq, jk, ij, jv) = fbdata1 % ivlqcf(jq, jk, ji, jv)
            END DO
          END DO
        END IF
      END DO
    END DO
    IF (fbdata1 % lgrid) THEN
      DO jv = 1, fbdata1 % nvar
        fbdata2 % cgrid(jv) = fbdata1 % cgrid(jv)
        ij = 0
        DO ji = 1, fbdata1 % nobs
          IF (llvalid(ji)) THEN
            ij = ij + 1
            fbdata2 % iproc(ij, jv) = fbdata1 % iproc(ji, jv)
            fbdata2 % iobsi(ij, jv) = fbdata1 % iobsi(ji, jv)
            fbdata2 % iobsj(ij, jv) = fbdata1 % iobsj(ji, jv)
            DO jk = 1, fbdata1 % nlev
              fbdata2 % iobsk(jk, ij, jv) = fbdata1 % iobsk(jk, ji, jv)
            END DO
          END IF
        END DO
      END DO
    END IF
    DO je = 1, fbdata1 % nadd
      fbdata2 % caddname(je) = fbdata1 % caddname(je)
    END DO
    DO jv = 1, fbdata1 % nvar
      DO je = 1, fbdata1 % nadd
        fbdata2 % caddlong(je, jv) = fbdata1 % caddlong(je, jv)
        fbdata2 % caddunit(je, jv) = fbdata1 % caddunit(je, jv)
        ij = 0
        DO ji = 1, fbdata1 % nobs
          IF (llvalid(ji)) THEN
            ij = ij + 1
            DO jk = 1, fbdata1 % nlev
              fbdata2 % padd(jk, ij, je, jv) = fbdata1 % padd(jk, ji, je, jv)
            END DO
          END IF
        END DO
      END DO
    END DO
    DO je = 1, fbdata1 % next
      fbdata2 % cextname(je) = fbdata1 % cextname(je)
      fbdata2 % cextlong(je) = fbdata1 % cextlong(je)
      fbdata2 % cextunit(je) = fbdata1 % cextunit(je)
    END DO
    DO je = 1, fbdata1 % next
      ij = 0
      DO ji = 1, fbdata1 % nobs
        IF (llvalid(ji)) THEN
          ij = ij + 1
          DO jk = 1, fbdata1 % nlev
            fbdata2 % pext(jk, ij, je) = fbdata1 % pext(jk, ji, je)
          END DO
        END IF
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE subsamp_obfbdata
  SUBROUTINE merge_obfbdata(nsets, fbdatain, fbdataout, iset, inum, iind)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: nsets
    TYPE(obfbdata), DIMENSION(nsets) :: fbdatain
    TYPE(obfbdata) :: fbdataout
    INTEGER, INTENT(IN), DIMENSION(fbdataout % nobs) :: iset
    INTEGER, INTENT(IN), DIMENSION(fbdataout % nobs) :: inum
    INTEGER, INTENT(IN), DIMENSION(fbdataout % nobs) :: iind
    INTEGER :: js
    INTEGER :: jo
    INTEGER :: jv
    INTEGER :: je
    INTEGER :: ji
    INTEGER :: jk
    INTEGER :: jq
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('merge_obfbdata', 'r0', 0, 0)
    DO js = 1, nsets
      IF (.NOT. fbdatain(js) % lalloc) THEN
        CALL fatal_error('merge_obfbdata: input data not allocated', 906)
      END IF
    END DO
    IF (.NOT. fbdataout % lalloc) THEN
      CALL fatal_error('merge_obfbdata: output data not allocated', 914)
    END IF
    DO jv = 1, fbdatain(1) % nvar
      fbdataout % cname(jv) = fbdatain(1) % cname(jv)
      fbdataout % coblong(jv) = fbdatain(1) % coblong(jv)
      fbdataout % cobunit(jv) = fbdatain(1) % cobunit(jv)
      IF (fbdatain(1) % lgrid) THEN
        fbdataout % cgrid(jv) = fbdatain(1) % cgrid(jv)
      END IF
    END DO
    DO jv = 1, fbdatain(1) % nadd
      fbdataout % caddname(jv) = fbdatain(1) % caddname(jv)
    END DO
    DO jv = 1, fbdatain(1) % nvar
      DO je = 1, fbdatain(1) % nadd
        fbdataout % caddlong(je, jv) = fbdatain(1) % caddlong(je, jv)
        fbdataout % caddunit(je, jv) = fbdatain(1) % caddunit(je, jv)
      END DO
    END DO
    DO jv = 1, fbdatain(1) % next
      fbdataout % cextname(jv) = fbdatain(1) % cextname(jv)
      fbdataout % cextlong(jv) = fbdatain(1) % cextlong(jv)
      fbdataout % cextunit(jv) = fbdatain(1) % cextunit(jv)
    END DO
    fbdataout % cdjuldref = fbdatain(1) % cdjuldref
    DO jo = 1, fbdataout % nobs
      js = iset(iind(jo))
      ji = inum(iind(jo))
      fbdataout % cdwmo(jo) = fbdatain(js) % cdwmo(ji)
      fbdataout % cdtyp(jo) = fbdatain(js) % cdtyp(ji)
      fbdataout % ioqc(jo) = fbdatain(js) % ioqc(ji)
      fbdataout % ipqc(jo) = fbdatain(js) % ipqc(ji)
      fbdataout % itqc(jo) = fbdatain(js) % itqc(ji)
      fbdataout % plam(jo) = fbdatain(js) % plam(ji)
      fbdataout % pphi(jo) = fbdatain(js) % pphi(ji)
      fbdataout % ptim(jo) = fbdatain(js) % ptim(ji)
      fbdataout % kindex(jo) = fbdatain(js) % kindex(ji)
      DO jq = 1, fbdatain(js) % nqcf
        fbdataout % ioqcf(jq, jo) = fbdatain(js) % ioqcf(jq, ji)
        fbdataout % ipqcf(jq, jo) = fbdatain(js) % ipqcf(jq, ji)
        fbdataout % itqcf(jq, jo) = fbdatain(js) % itqcf(jq, ji)
      END DO
      DO jk = 1, fbdatain(js) % nlev
        fbdataout % pdep(jk, jo) = fbdatain(js) % pdep(jk, ji)
        fbdataout % idqc(jk, jo) = fbdatain(js) % idqc(jk, ji)
        DO jq = 1, fbdatain(js) % nqcf
          fbdataout % idqcf(jq, jk, jo) = fbdatain(js) % idqcf(jq, jk, ji)
        END DO
      END DO
      DO jv = 1, fbdatain(js) % nvar
        fbdataout % ivqc(jo, jv) = fbdatain(js) % ivqc(ji, jv)
        DO jq = 1, fbdatain(js) % nqcf
          fbdataout % ivqcf(jq, jo, jv) = fbdatain(js) % ivqcf(jq, ji, jv)
        END DO
        DO jk = 1, fbdatain(js) % nlev
          fbdataout % ivlqc(jk, jo, jv) = fbdatain(js) % ivlqc(jk, ji, jv)
          fbdataout % pob(jk, jo, jv) = fbdatain(js) % pob(jk, ji, jv)
          DO jq = 1, fbdatain(js) % nqcf
            fbdataout % ivlqcf(jq, jk, jo, jv) = fbdatain(js) % ivlqcf(jq, jk, ji, jv)
          END DO
        END DO
      END DO
      IF (fbdatain(js) % lgrid) THEN
        DO jv = 1, fbdatain(js) % nvar
          fbdataout % cgrid(jv) = fbdatain(js) % cgrid(jv)
          fbdataout % iproc(jo, jv) = fbdatain(js) % iproc(ji, jv)
          fbdataout % iobsi(jo, jv) = fbdatain(js) % iobsi(ji, jv)
          fbdataout % iobsj(jo, jv) = fbdatain(js) % iobsj(ji, jv)
          DO jk = 1, fbdatain(js) % nlev
            fbdataout % iobsk(jk, jo, jv) = fbdatain(js) % iobsk(jk, ji, jv)
          END DO
        END DO
      END IF
      DO jv = 1, fbdatain(js) % nvar
        DO je = 1, fbdatain(js) % nadd
          DO jk = 1, fbdatain(js) % nlev
            fbdataout % padd(jk, jo, je, jv) = fbdatain(js) % padd(jk, ji, je, jv)
          END DO
        END DO
      END DO
      DO je = 1, fbdatain(js) % next
        DO jk = 1, fbdatain(js) % nlev
          fbdataout % pext(jk, jo, je) = fbdatain(js) % pext(jk, ji, je)
        END DO
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE merge_obfbdata
  SUBROUTINE write_obfbdata(cdfilename, fbdata)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *) :: cdfilename
    TYPE(obfbdata) :: fbdata
    CHARACTER(LEN = 14), PARAMETER :: cpname = 'write_obfbdata'
    INTEGER :: idfile
    INTEGER :: idodim
    INTEGER :: idldim
    INTEGER :: idvdim
    INTEGER :: idadim
    INTEGER :: idedim
    INTEGER :: idsndim
    INTEGER :: idsgdim
    INTEGER :: idswdim
    INTEGER :: idstdim
    INTEGER :: idjddim
    INTEGER :: idqcdim
    INTEGER :: idvard
    INTEGER :: idaddd
    INTEGER :: idextd
    INTEGER :: idcdwmo
    INTEGER :: idcdtyp
    INTEGER :: idplam
    INTEGER :: idpphi
    INTEGER :: idpdep
    INTEGER :: idptim
    INTEGER :: idptimr
    INTEGER :: idioqc
    INTEGER :: idioqcf
    INTEGER :: idipqc
    INTEGER :: idipqcf
    INTEGER :: iditqc
    INTEGER :: iditqcf
    INTEGER :: ididqc
    INTEGER :: ididqcf
    INTEGER :: idkindex
    INTEGER, DIMENSION(fbdata % nvar) :: idpob, idivqc, idivqcf, idivlqc, idivlqcf, idiobsi, idiobsj, idiobsk, idcgrid
    INTEGER, DIMENSION(fbdata % nadd, fbdata % nvar) :: idpadd
    INTEGER, DIMENSION(fbdata % next) :: idpext
    INTEGER, DIMENSION(1) :: incdim1
    INTEGER, DIMENSION(2) :: incdim2
    INTEGER, DIMENSION(3) :: incdim3
    INTEGER, DIMENSION(4) :: incdim4
    INTEGER :: jv
    INTEGER :: je
    INTEGER :: ioldfill
    CHARACTER(LEN = nf90_max_name) :: cdtmp
    CHARACTER(LEN = 16), PARAMETER :: cdqcconv = 'q where q =[0,9]'
    CHARACTER(LEN = 24), PARAMETER :: cdqcfconv = 'NEMOVAR flag conventions'
    CHARACTER(LEN = ilenlong) :: cdltmp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('write_obfbdata', 'r0', 0, 0)
    CALL chkerr(nf90_create(TRIM(cdfilename), nf90_clobber, idfile), cpname, 1107)
    CALL chkerr(nf90_set_fill(idfile, nf90_nofill, ioldfill), cpname, 1109)
    CALL chkerr(nf90_put_att(idfile, nf90_global, 'title', 'NEMO observation operator output'), cpname, 1112)
    CALL chkerr(nf90_put_att(idfile, nf90_global, 'Convention', 'NEMO unified observation operator output'), cpname, 1115)
    CALL chkerr(nf90_def_dim(idfile, 'N_OBS', fbdata % nobs, idodim), cpname, 1120)
    CALL chkerr(nf90_def_dim(idfile, 'N_LEVELS', fbdata % nlev, idldim), cpname, 1122)
    CALL chkerr(nf90_def_dim(idfile, 'N_VARS', fbdata % nvar, idvdim), cpname, 1124)
    CALL chkerr(nf90_def_dim(idfile, 'N_QCF', fbdata % nqcf, idqcdim), cpname, 1126)
    IF (fbdata % nadd > 0) THEN
      CALL chkerr(nf90_def_dim(idfile, 'N_ENTRIES', fbdata % nadd, idadim), cpname, 1129)
    END IF
    IF (fbdata % next > 0) THEN
      CALL chkerr(nf90_def_dim(idfile, 'N_EXTRA', fbdata % next, idedim), cpname, 1133)
    END IF
    CALL chkerr(nf90_def_dim(idfile, 'STRINGNAM', ilenname, idsndim), cpname, 1136)
    IF (fbdata % lgrid) THEN
      CALL chkerr(nf90_def_dim(idfile, 'STRINGGRID', ilengrid, idsgdim), cpname, 1139)
    END IF
    CALL chkerr(nf90_def_dim(idfile, 'STRINGWMO', ilenwmo, idswdim), cpname, 1142)
    CALL chkerr(nf90_def_dim(idfile, 'STRINGTYP', ilentyp, idstdim), cpname, 1144)
    CALL chkerr(nf90_def_dim(idfile, 'STRINGJULD', ilenjuld, idjddim), cpname, 1146)
    incdim2(1) = idsndim
    incdim2(2) = idvdim
    CALL chkerr(nf90_def_var(idfile, 'VARIABLES', nf90_char, incdim2, idvard), cpname, 1154)
    CALL putvaratt_obfbdata(idfile, idvard, 'List of variables in feedback files')
    IF (fbdata % nadd > 0) THEN
      incdim2(1) = idsndim
      incdim2(2) = idadim
      CALL chkerr(nf90_def_var(idfile, 'ENTRIES', nf90_char, incdim2, idaddd), cpname, 1162)
      CALL putvaratt_obfbdata(idfile, idaddd, 'List of additional entries for each ' // 'variable in feedback files')
    END IF
    IF (fbdata % next > 0) THEN
      incdim2(1) = idsndim
      incdim2(2) = idedim
      CALL chkerr(nf90_def_var(idfile, 'EXTRA', nf90_char, incdim2, idextd), cpname, 1172)
      CALL putvaratt_obfbdata(idfile, idextd, 'List of extra variables')
    END IF
    incdim2(1) = idswdim
    incdim2(2) = idodim
    CALL chkerr(nf90_def_var(idfile, 'STATION_IDENTIFIER', nf90_char, incdim2, idcdwmo), cpname, 1181)
    CALL putvaratt_obfbdata(idfile, idcdwmo, 'Station identifier')
    incdim2(1) = idstdim
    incdim2(2) = idodim
    CALL chkerr(nf90_def_var(idfile, 'STATION_TYPE', nf90_char, incdim2, idcdtyp), cpname, 1188)
    CALL putvaratt_obfbdata(idfile, idcdtyp, 'Code instrument type')
    incdim1(1) = idodim
    CALL chkerr(nf90_def_var(idfile, 'LONGITUDE', nf90_double, incdim1, idplam), cpname, 1194)
    CALL putvaratt_obfbdata(idfile, idplam, 'Longitude', cdunits = 'degrees_east', rfillvalue = fbrmdi)
    CALL chkerr(nf90_def_var(idfile, 'LATITUDE', nf90_double, incdim1, idpphi), cpname, 1200)
    CALL putvaratt_obfbdata(idfile, idpphi, 'Latitude', cdunits = 'degrees_north', rfillvalue = fbrmdi)
    incdim2(1) = idldim
    incdim2(2) = idodim
    CALL chkerr(nf90_def_var(idfile, 'DEPTH', nf90_double, incdim2, idpdep), cpname, 1208)
    CALL putvaratt_obfbdata(idfile, idpdep, 'Depth', cdunits = 'metre', rfillvalue = fbrmdi)
    incdim3(1) = idqcdim
    incdim3(2) = idldim
    incdim3(3) = idodim
    CALL chkerr(nf90_def_var(idfile, 'DEPTH_QC', nf90_int, incdim2, ididqc), cpname, 1217)
    CALL putvaratt_obfbdata(idfile, ididqc, 'Quality on depth', conventions = cdqcconv, ifillvalue = 0)
    CALL chkerr(nf90_def_var(idfile, 'DEPTH_QC_FLAGS', nf90_int, incdim3, ididqcf), cpname, 1224)
    CALL putvaratt_obfbdata(idfile, ididqcf, 'Quality flags on depth', conventions = cdqcfconv)
    CALL chkerr(nf90_def_var(idfile, 'JULD', nf90_double, incdim1, idptim), cpname, 1230)
    CALL putvaratt_obfbdata(idfile, idptim, 'Julian day', cdunits = 'days since JULD_REFERENCE', &
&conventions = 'relative julian days with ' // 'decimal part (as parts of day)', rfillvalue = fbrmdi)
    incdim1(1) = idjddim
    CALL chkerr(nf90_def_var(idfile, 'JULD_REFERENCE', nf90_char, incdim1, idptimr), cpname, 1240)
    CALL putvaratt_obfbdata(idfile, idptimr, 'Date of reference for julian days ', conventions = 'YYYYMMDDHHMMSS')
    incdim1(1) = idodim
    CALL chkerr(nf90_def_var(idfile, 'OBSERVATION_QC', nf90_int, incdim1, idioqc), cpname, 1247)
    CALL putvaratt_obfbdata(idfile, idioqc, 'Quality on observation', conventions = cdqcconv, ifillvalue = 0)
    incdim2(1) = idqcdim
    incdim2(2) = idodim
    CALL chkerr(nf90_def_var(idfile, 'OBSERVATION_QC_FLAGS', nf90_int, incdim2, idioqcf), cpname, 1256)
    CALL putvaratt_obfbdata(idfile, idioqcf, 'Quality flags on observation', conventions = cdqcfconv, ifillvalue = 0)
    CALL chkerr(nf90_def_var(idfile, 'POSITION_QC', nf90_int, incdim1, idipqc), cpname, 1263)
    CALL putvaratt_obfbdata(idfile, idipqc, 'Quality on position (latitude and longitude)', conventions = cdqcconv, ifillvalue = 0)
    CALL chkerr(nf90_def_var(idfile, 'POSITION_QC_FLAGS', nf90_int, incdim2, idipqcf), cpname, 1270)
    CALL putvaratt_obfbdata(idfile, idipqcf, 'Quality flags on position', conventions = cdqcfconv, ifillvalue = 0)
    CALL chkerr(nf90_def_var(idfile, 'JULD_QC', nf90_int, incdim1, iditqc), cpname, 1277)
    CALL putvaratt_obfbdata(idfile, iditqc, 'Quality on date and time', conventions = cdqcconv, ifillvalue = 0)
    CALL chkerr(nf90_def_var(idfile, 'JULD_QC_FLAGS', nf90_int, incdim2, iditqcf), cpname, 1284)
    CALL putvaratt_obfbdata(idfile, iditqcf, 'Quality flags on date and time', conventions = cdqcfconv, ifillvalue = 0)
    CALL chkerr(nf90_def_var(idfile, 'ORIGINAL_FILE_INDEX', nf90_int, incdim1, idkindex), cpname, 1291)
    CALL putvaratt_obfbdata(idfile, idkindex, 'Index in original data file', ifillvalue = fbimdi)
    CALL profile_psy_data0 % PostEnd
    DO jv = 1, fbdata % nvar
      CALL profile_psy_data1 % PreStart('write_obfbdata', 'r1', 0, 0)
      incdim1(1) = idodim
      incdim2(1) = idldim
      incdim2(2) = idodim
      WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_OBS'
      CALL chkerr(nf90_def_var(idfile, cdtmp, nf90_float, incdim2, idpob(jv)), cpname, 1306)
      CALL putvaratt_obfbdata(idfile, idpob(jv), fbdata % coblong(jv), cdunits = fbdata % cobunit(jv), rfillvalue = fbrmdi)
      IF (fbdata % nadd > 0) THEN
        DO je = 1, fbdata % nadd
          WRITE(cdtmp, FMT = '(3A)') TRIM(fbdata % cname(jv)), '_', TRIM(fbdata % caddname(je))
          CALL chkerr(nf90_def_var(idfile, cdtmp, nf90_float, incdim2, idpadd(je, jv)), cpname, 1318)
          CALL putvaratt_obfbdata(idfile, idpadd(je, jv), fbdata % caddlong(je, jv), cdunits = fbdata % caddunit(je, jv), &
&rfillvalue = fbrmdi)
        END DO
      END IF
      cdltmp = fbdata % coblong(jv)
      CALL profile_psy_data1 % PostEnd
      IF ((cdltmp(1 : 1) >= 'A') .AND. (cdltmp(1 : 1) <= 'Z')) cdltmp(1 : 1) = ACHAR(IACHAR(cdltmp(1 : 1)) + 32)
      CALL profile_psy_data2 % PreStart('write_obfbdata', 'r2', 0, 0)
      WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_QC'
      CALL chkerr(nf90_def_var(idfile, cdtmp, nf90_int, incdim1, idivqc(jv)), cpname, 1332)
      CALL putvaratt_obfbdata(idfile, idivqc(jv), 'Quality on ' // cdltmp, conventions = cdqcconv, ifillvalue = 0)
      incdim2(1) = idqcdim
      incdim2(2) = idodim
      WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_QC_FLAGS'
      CALL chkerr(nf90_def_var(idfile, cdtmp, nf90_int, incdim2, idivqcf(jv)), cpname, 1342)
      CALL putvaratt_obfbdata(idfile, idivqcf(jv), 'Quality flags on ' // cdltmp, conventions = cdqcfconv, ifillvalue = 0)
      incdim2(1) = idldim
      incdim2(2) = idodim
      WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_LEVEL_QC'
      CALL chkerr(nf90_def_var(idfile, cdtmp, nf90_int, incdim2, idivlqc(jv)), cpname, 1352)
      CALL putvaratt_obfbdata(idfile, idivlqc(jv), 'Quality for each level on ' // cdltmp, conventions = cdqcconv, ifillvalue = 0)
      incdim3(1) = idqcdim
      incdim3(2) = idldim
      incdim3(3) = idodim
      WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_LEVEL_QC_FLAGS'
      CALL chkerr(nf90_def_var(idfile, cdtmp, nf90_int, incdim3, idivlqcf(jv)), cpname, 1363)
      CALL putvaratt_obfbdata(idfile, idivlqcf(jv), 'Quality flags for each level on ' // cdltmp, conventions = cdqcfconv, &
&ifillvalue = 0)
      IF (fbdata % lgrid) THEN
        incdim2(1) = idldim
        incdim2(2) = idodim
        WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_IOBSI'
        CALL chkerr(nf90_def_var(idfile, cdtmp, nf90_int, incdim1, idiobsi(jv)), cpname, 1376)
        CALL putvaratt_obfbdata(idfile, idiobsi(jv), 'ORCA grid search I coordinate')
        WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_IOBSJ'
        CALL chkerr(nf90_def_var(idfile, cdtmp, nf90_int, incdim1, idiobsj(jv)), cpname, 1382)
        CALL putvaratt_obfbdata(idfile, idiobsj(jv), 'ORCA grid search J coordinate')
        WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_IOBSK'
        CALL chkerr(nf90_def_var(idfile, cdtmp, nf90_int, incdim2, idiobsk(jv)), cpname, 1388)
        CALL putvaratt_obfbdata(idfile, idiobsk(jv), 'ORCA grid search K coordinate')
        incdim1(1) = idsgdim
        WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_GRID'
        CALL chkerr(nf90_def_var(idfile, cdtmp, nf90_char, incdim1, idcgrid(jv)), cpname, 1394)
        CALL putvaratt_obfbdata(idfile, idcgrid(jv), 'ORCA grid search grid (T,U,V)')
      END IF
      CALL profile_psy_data2 % PostEnd
    END DO
    CALL profile_psy_data3 % PreStart('write_obfbdata', 'r3', 0, 0)
    IF (fbdata % next > 0) THEN
      DO je = 1, fbdata % next
        incdim2(1) = idldim
        incdim2(2) = idodim
        WRITE(cdtmp, FMT = '(A)') TRIM(fbdata % cextname(je))
        CALL chkerr(nf90_def_var(idfile, cdtmp, nf90_float, incdim2, idpext(je)), cpname, 1408)
        CALL putvaratt_obfbdata(idfile, idpext(je), fbdata % cextlong(je), cdunits = fbdata % cextunit(je), rfillvalue = fbrmdi)
      END DO
    END IF
    CALL chkerr(nf90_enddef(idfile), cpname, 1418)
    CALL chkerr(nf90_put_var(idfile, idvard, fbdata % cname), cpname, 1423)
    IF (fbdata % nadd > 0) THEN
      CALL chkerr(nf90_put_var(idfile, idaddd, fbdata % caddname), cpname, 1427)
    END IF
    IF (fbdata % next > 0) THEN
      CALL chkerr(nf90_put_var(idfile, idextd, fbdata % cextname), cpname, 1432)
    END IF
    CALL chkerr(nf90_put_var(idfile, idptimr, fbdata % cdjuldref), cpname, 1436)
    IF (fbdata % nobs > 0) THEN
      CALL chkerr(nf90_put_var(idfile, idcdwmo, fbdata % cdwmo), cpname, 1443)
      CALL chkerr(nf90_put_var(idfile, idcdtyp, fbdata % cdtyp), cpname, 1445)
      CALL chkerr(nf90_put_var(idfile, idplam, fbdata % plam), cpname, 1447)
      CALL chkerr(nf90_put_var(idfile, idpphi, fbdata % pphi), cpname, 1449)
      CALL chkerr(nf90_put_var(idfile, idpdep, fbdata % pdep), cpname, 1451)
      CALL chkerr(nf90_put_var(idfile, idptim, fbdata % ptim), cpname, 1453)
      CALL chkerr(nf90_put_var(idfile, idioqc, fbdata % ioqc), cpname, 1455)
      CALL chkerr(nf90_put_var(idfile, idioqcf, fbdata % ioqcf), cpname, 1457)
      CALL chkerr(nf90_put_var(idfile, idipqc, fbdata % ipqc), cpname, 1459)
      CALL chkerr(nf90_put_var(idfile, idipqcf, fbdata % ipqcf), cpname, 1461)
      CALL chkerr(nf90_put_var(idfile, iditqc, fbdata % itqc), cpname, 1463)
      CALL chkerr(nf90_put_var(idfile, iditqcf, fbdata % itqcf), cpname, 1465)
      CALL chkerr(nf90_put_var(idfile, ididqc, fbdata % idqc), cpname, 1467)
      CALL chkerr(nf90_put_var(idfile, ididqcf, fbdata % idqcf), cpname, 1469)
      CALL chkerr(nf90_put_var(idfile, idkindex, fbdata % kindex), cpname, 1471)
      DO jv = 1, fbdata % nvar
        CALL chkerr(nf90_put_var(idfile, idpob(jv), fbdata % pob(:, :, jv)), cpname, 1475)
        IF (fbdata % nadd > 0) THEN
          DO je = 1, fbdata % nadd
            CALL chkerr(nf90_put_var(idfile, idpadd(je, jv), fbdata % padd(:, :, je, jv)), cpname, 1480)
          END DO
        END IF
        CALL chkerr(nf90_put_var(idfile, idivqc(jv), fbdata % ivqc(:, jv)), cpname, 1485)
        CALL chkerr(nf90_put_var(idfile, idivqcf(jv), fbdata % ivqcf(:, :, jv)), cpname, 1488)
        CALL chkerr(nf90_put_var(idfile, idivlqc(jv), fbdata % ivlqc(:, :, jv)), cpname, 1491)
        CALL chkerr(nf90_put_var(idfile, idivlqcf(jv), fbdata % ivlqcf(:, :, :, jv)), cpname, 1494)
        IF (fbdata % lgrid) THEN
          CALL chkerr(nf90_put_var(idfile, idiobsi(jv), fbdata % iobsi(:, jv)), cpname, 1498)
          CALL chkerr(nf90_put_var(idfile, idiobsj(jv), fbdata % iobsj(:, jv)), cpname, 1501)
          CALL chkerr(nf90_put_var(idfile, idiobsk(jv), fbdata % iobsk(:, :, jv)), cpname, 1504)
          CALL chkerr(nf90_put_var(idfile, idcgrid(jv), fbdata % cgrid(jv)), cpname, 1507)
        END IF
      END DO
      IF (fbdata % next > 0) THEN
        DO je = 1, fbdata % next
          CALL chkerr(nf90_put_var(idfile, idpext(je), fbdata % pext(:, :, je)), cpname, 1515)
        END DO
      END IF
    END IF
    CALL chkerr(nf90_close(idfile), cpname, 1523)
    CALL profile_psy_data3 % PostEnd
  END SUBROUTINE write_obfbdata
  SUBROUTINE putvaratt_obfbdata(idfile, idvar, cdlongname, cdunits, conventions, cfillvalue, ifillvalue, rfillvalue)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: idfile
    INTEGER :: idvar
    CHARACTER(LEN = *) :: cdlongname
    CHARACTER(LEN = *), OPTIONAL :: cdunits
    CHARACTER(LEN = *), OPTIONAL :: cfillvalue
    INTEGER, OPTIONAL :: ifillvalue
    REAL(KIND = fbsp), OPTIONAL :: rfillvalue
    CHARACTER(LEN = *), OPTIONAL :: conventions
    CHARACTER(LEN = 18), PARAMETER :: cpname = 'putvaratt_obfbdata'
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('putvaratt_obfbdata', 'r0', 0, 0)
    CALL chkerr(nf90_put_att(idfile, idvar, 'long_name', TRIM(cdlongname)), cpname, 1556)
    IF (PRESENT(cdunits)) THEN
      CALL chkerr(nf90_put_att(idfile, idvar, 'units', TRIM(cdunits)), cpname, 1562)
    END IF
    IF (PRESENT(conventions)) THEN
      CALL chkerr(nf90_put_att(idfile, idvar, 'Conventions', TRIM(conventions)), cpname, 1570)
    END IF
    IF (PRESENT(cfillvalue)) THEN
      CALL chkerr(nf90_put_att(idfile, idvar, '_Fillvalue', TRIM(cfillvalue)), cpname, 1578)
    END IF
    IF (PRESENT(ifillvalue)) THEN
      CALL chkerr(nf90_put_att(idfile, idvar, '_Fillvalue', ifillvalue), cpname, 1586)
    END IF
    IF (PRESENT(rfillvalue)) THEN
      CALL chkerr(nf90_put_att(idfile, idvar, '_Fillvalue', rfillvalue), cpname, 1594)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE putvaratt_obfbdata
  SUBROUTINE read_obfbdata(cdfilename, fbdata, ldgrid)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *) :: cdfilename
    TYPE(obfbdata) :: fbdata
    LOGICAL, OPTIONAL :: ldgrid
    CHARACTER(LEN = 14), PARAMETER :: cpname = 'read_obfbdata'
    INTEGER :: idfile
    INTEGER :: idodim
    INTEGER :: idldim
    INTEGER :: idvdim
    INTEGER :: idadim
    INTEGER :: idedim
    INTEGER :: idgdim
    INTEGER :: idvard
    INTEGER :: idaddd
    INTEGER :: idextd
    INTEGER :: idcdwmo
    INTEGER :: idcdtyp
    INTEGER :: idplam
    INTEGER :: idpphi
    INTEGER :: idpdep
    INTEGER :: idptim
    INTEGER :: idptimr
    INTEGER :: idioqc
    INTEGER :: idioqcf
    INTEGER :: idipqc
    INTEGER :: idipqcf
    INTEGER :: ididqc
    INTEGER :: ididqcf
    INTEGER :: iditqc
    INTEGER :: iditqcf
    INTEGER :: idkindex
    INTEGER, DIMENSION(:), ALLOCATABLE :: idpob, idivqc, idivqcf, idivlqc, idivlqcf, idiobsi, idiobsj, idiobsk, idcgrid, idpext
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: idpadd
    INTEGER :: jv
    INTEGER :: je
    INTEGER :: nvar
    INTEGER :: nobs
    INTEGER :: nlev
    INTEGER :: nadd
    INTEGER :: next
    LOGICAL :: lgrid
    CHARACTER(LEN = NF90_MAX_NAME) :: cdtmp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('read_obfbdata', 'r0', 0, 0)
    IF (fbdata % lalloc) THEN
      CALL dealloc_obfbdata(fbdata)
    END IF
    CALL chkerr(nf90_open(TRIM(cdfilename), nf90_nowrite, idfile), cpname, 1675)
    CALL chkerr(nf90_inq_dimid(idfile, 'N_OBS', idodim), cpname, 1680)
    CALL chkerr(nf90_inquire_dimension(idfile, idodim, len = nobs), cpname, 1682)
    CALL chkerr(nf90_inq_dimid(idfile, 'N_LEVELS', idldim), cpname, 1684)
    CALL chkerr(nf90_inquire_dimension(idfile, idldim, len = nlev), cpname, 1686)
    CALL chkerr(nf90_inq_dimid(idfile, 'N_VARS', idvdim), cpname, 1688)
    CALL chkerr(nf90_inquire_dimension(idfile, idvdim, len = nvar), cpname, 1690)
    IF (nf90_inq_dimid(idfile, 'N_ENTRIES', idadim) == 0) THEN
      CALL chkerr(nf90_inquire_dimension(idfile, idadim, len = nadd), cpname, 1693)
    ELSE
      nadd = 0
    END IF
    IF (nf90_inq_dimid(idfile, 'N_EXTRA', idedim) == 0) THEN
      CALL chkerr(nf90_inquire_dimension(idfile, idedim, len = next), cpname, 1699)
    ELSE
      next = 0
    END IF
    lgrid = (nf90_inq_dimid(idfile, 'STRINGGRID', idgdim) == 0)
    IF (PRESENT(ldgrid)) THEN
      CALL alloc_obfbdata(fbdata, nvar, nobs, nlev, nadd, next, lgrid .OR. ldgrid)
    ELSE
      CALL alloc_obfbdata(fbdata, nvar, nobs, nlev, nadd, next, lgrid)
    END IF
    ALLOCATE(idpob(fbdata % nvar), idivqc(fbdata % nvar), idivqcf(fbdata % nvar), idivlqc(fbdata % nvar), idivlqcf(fbdata % nvar), &
&idiobsi(fbdata % nvar), idiobsj(fbdata % nvar), idiobsk(fbdata % nvar), idcgrid(fbdata % nvar))
    IF (fbdata % nadd > 0) THEN
      ALLOCATE(idpadd(fbdata % nadd, fbdata % nvar))
    END IF
    IF (fbdata % next > 0) THEN
      ALLOCATE(idpext(fbdata % next))
    END IF
    CALL chkerr(nf90_inq_varid(idfile, 'VARIABLES', idvard), cpname, 1745)
    CALL chkerr(nf90_get_var(idfile, idvard, fbdata % cname), cpname, 1747)
    IF (fbdata % nadd > 0) THEN
      CALL chkerr(nf90_inq_varid(idfile, 'ENTRIES', idaddd), cpname, 1750)
      CALL chkerr(nf90_get_var(idfile, idaddd, fbdata % caddname), cpname, 1752)
    END IF
    IF (fbdata % next > 0) THEN
      CALL chkerr(nf90_inq_varid(idfile, 'EXTRA', idextd), cpname, 1756)
      CALL chkerr(nf90_get_var(idfile, idextd, fbdata % cextname), cpname, 1758)
    END IF
    CALL chkerr(nf90_inq_varid(idfile, 'JULD_REFERENCE', idptimr), cpname, 1762)
    CALL chkerr(nf90_get_var(idfile, idptimr, fbdata % cdjuldref), cpname, 1764)
    IF (fbdata % nobs > 0) THEN
      CALL chkerr(nf90_inq_varid(idfile, 'STATION_IDENTIFIER', idcdwmo), cpname, 1769)
      CALL chkerr(nf90_get_var(idfile, idcdwmo, fbdata % cdwmo), cpname, 1771)
      CALL chkerr(nf90_inq_varid(idfile, 'STATION_TYPE', idcdtyp), cpname, 1773)
      CALL chkerr(nf90_get_var(idfile, idcdtyp, fbdata % cdtyp), cpname, 1775)
      CALL chkerr(nf90_inq_varid(idfile, 'LONGITUDE', idplam), cpname, 1777)
      CALL chkerr(nf90_get_var(idfile, idplam, fbdata % plam), cpname, 1779)
      CALL chkerr(nf90_inq_varid(idfile, 'LATITUDE', idpphi), cpname, 1781)
      CALL chkerr(nf90_get_var(idfile, idpphi, fbdata % pphi), cpname, 1783)
      CALL chkerr(nf90_inq_varid(idfile, 'DEPTH', idpdep), cpname, 1785)
      CALL chkerr(nf90_get_var(idfile, idpdep, fbdata % pdep), cpname, 1787)
      CALL chkerr(nf90_inq_varid(idfile, 'JULD', idptim), cpname, 1789)
      CALL chkerr(nf90_get_var(idfile, idptim, fbdata % ptim), cpname, 1791)
      CALL chkerr(nf90_inq_varid(idfile, 'OBSERVATION_QC', idioqc), cpname, 1793)
      CALL chkerr(nf90_get_var(idfile, idioqc, fbdata % ioqc), cpname, 1795)
      CALL chkerr(nf90_inq_varid(idfile, 'OBSERVATION_QC_FLAGS', idioqcf), cpname, 1797)
      CALL chkerr(nf90_get_var(idfile, idioqcf, fbdata % ioqcf), cpname, 1799)
      CALL chkerr(nf90_inq_varid(idfile, 'POSITION_QC', idipqc), cpname, 1801)
      CALL chkerr(nf90_get_var(idfile, idipqc, fbdata % ipqc), cpname, 1803)
      CALL chkerr(nf90_inq_varid(idfile, 'POSITION_QC_FLAGS', idipqcf), cpname, 1805)
      CALL chkerr(nf90_get_var(idfile, idipqcf, fbdata % ipqcf), cpname, 1807)
      CALL chkerr(nf90_inq_varid(idfile, 'DEPTH_QC', ididqc), cpname, 1809)
      CALL chkerr(nf90_get_var(idfile, ididqc, fbdata % idqc), cpname, 1811)
      CALL chkerr(nf90_inq_varid(idfile, 'DEPTH_QC_FLAGS', ididqcf), cpname, 1813)
      CALL chkerr(nf90_get_var(idfile, ididqcf, fbdata % idqcf), cpname, 1815)
      CALL chkerr(nf90_inq_varid(idfile, 'JULD_QC', iditqc), cpname, 1817)
      CALL chkerr(nf90_get_var(idfile, iditqc, fbdata % itqc), cpname, 1819)
      CALL chkerr(nf90_inq_varid(idfile, 'JULD_QC_FLAGS', iditqcf), cpname, 1821)
      CALL chkerr(nf90_get_var(idfile, iditqcf, fbdata % itqcf), cpname, 1823)
      CALL chkerr(nf90_inq_varid(idfile, 'ORIGINAL_FILE_INDEX', idkindex), cpname, 1825)
      CALL chkerr(nf90_get_var(idfile, idkindex, fbdata % kindex), cpname, 1827)
      DO jv = 1, fbdata % nvar
        WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_OBS'
        CALL chkerr(nf90_inq_varid(idfile, cdtmp, idpob(jv)), cpname, 1835)
        CALL chkerr(nf90_get_var(idfile, idpob(jv), fbdata % pob(:, :, jv)), cpname, 1838)
        CALL getvaratt_obfbdata(idfile, idpob(jv), fbdata % coblong(jv), fbdata % cobunit(jv))
        IF (fbdata % nadd > 0) THEN
          DO je = 1, fbdata % nadd
            WRITE(cdtmp, FMT = '(3A)') TRIM(fbdata % cname(jv)), '_', TRIM(fbdata % caddname(je))
            CALL chkerr(nf90_inq_varid(idfile, cdtmp, idpadd(je, jv)), cpname, 1848)
            CALL chkerr(nf90_get_var(idfile, idpadd(je, jv), fbdata % padd(:, :, je, jv)), cpname, 1851)
            CALL getvaratt_obfbdata(idfile, idpadd(je, jv), fbdata % caddlong(je, jv), fbdata % caddunit(je, jv))
          END DO
        END IF
        WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_QC'
        CALL chkerr(nf90_inq_varid(idfile, cdtmp, idivqc(jv)), cpname, 1860)
        CALL chkerr(nf90_get_var(idfile, idivqc(jv), fbdata % ivqc(:, jv)), cpname, 1863)
        WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_QC_FLAGS'
        CALL chkerr(nf90_inq_varid(idfile, cdtmp, idivqcf(jv)), cpname, 1866)
        CALL chkerr(nf90_get_var(idfile, idivqcf(jv), fbdata % ivqcf(:, :, jv)), cpname, 1869)
        WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_LEVEL_QC'
        CALL chkerr(nf90_inq_varid(idfile, cdtmp, idivlqc(jv)), cpname, 1872)
        CALL chkerr(nf90_get_var(idfile, idivlqc(jv), fbdata % ivlqc(:, :, jv)), cpname, 1875)
        WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_LEVEL_QC_FLAGS'
        CALL chkerr(nf90_inq_varid(idfile, cdtmp, idivlqcf(jv)), cpname, 1878)
        CALL chkerr(nf90_get_var(idfile, idivlqcf(jv), fbdata % ivlqcf(:, :, :, jv)), cpname, 1881)
        IF (lgrid) THEN
          WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_IOBSI'
          CALL chkerr(nf90_inq_varid(idfile, cdtmp, idiobsi(jv)), cpname, 1885)
          CALL chkerr(nf90_get_var(idfile, idiobsi(jv), fbdata % iobsi(:, jv)), cpname, 1888)
          WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_IOBSJ'
          CALL chkerr(nf90_inq_varid(idfile, cdtmp, idiobsj(jv)), cpname, 1891)
          CALL chkerr(nf90_get_var(idfile, idiobsj(jv), fbdata % iobsj(:, jv)), cpname, 1894)
          WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_IOBSK'
          CALL chkerr(nf90_inq_varid(idfile, cdtmp, idiobsk(jv)), cpname, 1897)
          CALL chkerr(nf90_get_var(idfile, idiobsk(jv), fbdata % iobsk(:, :, jv)), cpname, 1900)
          WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_GRID'
          CALL chkerr(nf90_inq_varid(idfile, cdtmp, idcgrid(jv)), cpname, 1903)
          CALL chkerr(nf90_get_var(idfile, idcgrid(jv), fbdata % cgrid(jv)), cpname, 1906)
        END IF
      END DO
      IF (fbdata % next > 0) THEN
        DO je = 1, fbdata % next
          WRITE(cdtmp, FMT = '(A)') TRIM(fbdata % cextname(je))
          CALL chkerr(nf90_inq_varid(idfile, cdtmp, idpext(je)), cpname, 1915)
          CALL chkerr(nf90_get_var(idfile, idpext(je), fbdata % pext(:, :, je)), cpname, 1918)
          CALL getvaratt_obfbdata(idfile, idpext(je), fbdata % cextlong(je), fbdata % cextunit(je))
        END DO
      END IF
    ELSE
      DO jv = 1, fbdata % nvar
        WRITE(cdtmp, FMT = '(2A)') TRIM(fbdata % cname(jv)), '_OBS'
        CALL chkerr(nf90_inq_varid(idfile, cdtmp, idpob(jv)), cpname, 1931)
        CALL getvaratt_obfbdata(idfile, idpob(jv), fbdata % coblong(jv), fbdata % cobunit(jv))
        IF (fbdata % nadd > 0) THEN
          DO je = 1, fbdata % nadd
            WRITE(cdtmp, FMT = '(3A)') TRIM(fbdata % cname(jv)), '_', TRIM(fbdata % caddname(je))
            CALL chkerr(nf90_inq_varid(idfile, cdtmp, idpadd(je, jv)), cpname, 1941)
            CALL getvaratt_obfbdata(idfile, idpadd(je, jv), fbdata % caddlong(je, jv), fbdata % caddunit(je, jv))
          END DO
        END IF
      END DO
      IF (fbdata % next > 0) THEN
        DO je = 1, fbdata % next
          WRITE(cdtmp, FMT = '(A)') TRIM(fbdata % cextname(je))
          CALL chkerr(nf90_inq_varid(idfile, cdtmp, idpext(je)), cpname, 1954)
          CALL getvaratt_obfbdata(idfile, idpext(je), fbdata % cextlong(je), fbdata % cextunit(je))
        END DO
      END IF
    END IF
    CALL chkerr(nf90_close(idfile), cpname, 1965)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE read_obfbdata
  SUBROUTINE getvaratt_obfbdata(idfile, idvar, cdlongname, cdunits)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: idfile
    INTEGER :: idvar
    CHARACTER(LEN = *) :: cdlongname
    CHARACTER(LEN = *) :: cdunits
    CHARACTER(LEN = 18), PARAMETER :: cpname = 'getvaratt_obfbdata'
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('getvaratt_obfbdata', 'r0', 0, 0)
    CALL chkerr(nf90_get_att(idfile, idvar, 'long_name', cdlongname), cpname, 1990)
    CALL chkerr(nf90_get_att(idfile, idvar, 'units', cdunits), cpname, 1994)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE getvaratt_obfbdata
END MODULE obs_fbm
