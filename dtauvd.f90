MODULE dtauvd
  USE oce
  USE phycst
  USE dom_oce
  USE in_out_manager
  USE fldread
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dta_uvd_init
  PUBLIC :: dta_uvd
  LOGICAL, PUBLIC :: ln_uvd_init
  LOGICAL, PUBLIC :: ln_uvd_dyndmp
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_uvd
  CONTAINS
  SUBROUTINE dta_uvd_init(ld_dyndmp)
    LOGICAL, INTENT(IN), OPTIONAL :: ld_dyndmp
    INTEGER :: ios, ierr0, ierr1, ierr2, ierr3
    CHARACTER(LEN = 100) :: cn_dir
    TYPE(FLD_N), DIMENSION(2) :: suv_i
    TYPE(FLD_N) :: sn_ucur, sn_vcur
    NAMELIST /namc1d_uvd/ ln_uvd_init, ln_uvd_dyndmp, cn_dir, sn_ucur, sn_vcur
    ierr0 = 0
    ierr1 = 0
    ierr2 = 0
    ierr3 = 0
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namc1d_uvd, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namc1d_uvd in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namc1d_uvd, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namc1d_uvd in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namc1d_uvd)
    IF (PRESENT(ld_dyndmp)) ln_uvd_dyndmp = .TRUE.
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dta_uvd_init : U & V current data '
      WRITE(numout, FMT = *) '~~~~~~~~~~~~ '
      WRITE(numout, FMT = *) '   Namelist namc1d_uvd : Set flags'
      WRITE(numout, FMT = *) '      Initialization of ocean U & V current with input data   ln_uvd_init   = ', ln_uvd_init
      WRITE(numout, FMT = *) '      Damping of ocean U & V current toward input data        ln_uvd_dyndmp = ', ln_uvd_dyndmp
      WRITE(numout, FMT = *)
      IF (.NOT. ln_uvd_init .AND. .NOT. ln_uvd_dyndmp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) '   U & V current data not used'
      END IF
    END IF
    IF (ln_rstart .AND. ln_uvd_init) THEN
      CALL ctl_warn('dta_uvd_init: ocean restart and U & V current data initialization, ', &
&'we keep the restart U & V current values and set ln_uvd_init to FALSE')
      ln_uvd_init = .FALSE.
    END IF
    IF (ln_uvd_init .OR. ln_uvd_dyndmp) THEN
      ALLOCATE(sf_uvd(2), STAT = ierr0)
      IF (ierr0 > 0) THEN
        CALL ctl_stop('dta_uvd_init: unable to allocate sf_uvd structure')
        RETURN
      END IF
      ALLOCATE(sf_uvd(1) % fnow(jpi, jpj, jpk), STAT = ierr0)
      IF (sn_ucur % ln_tint) ALLOCATE(sf_uvd(1) % fdta(jpi, jpj, jpk, 2), STAT = ierr1)
      ALLOCATE(sf_uvd(2) % fnow(jpi, jpj, jpk), STAT = ierr2)
      IF (sn_vcur % ln_tint) ALLOCATE(sf_uvd(2) % fdta(jpi, jpj, jpk, 2), STAT = ierr3)
      IF (ierr0 + ierr1 + ierr2 + ierr3 > 0) THEN
        CALL ctl_stop('dta_uvd_init : unable to allocate U & V current data arrays')
        RETURN
      END IF
      suv_i(1) = sn_ucur
      suv_i(2) = sn_vcur
      CALL fld_fill(sf_uvd, suv_i, cn_dir, 'dta_uvd', 'U & V current data', 'namc1d_uvd')
    END IF
  END SUBROUTINE dta_uvd_init
  SUBROUTINE dta_uvd(kt, puvd)
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, 2), INTENT(OUT) :: puvd
    INTEGER :: ji, jj, jk, jl, jkk
    INTEGER :: ik, il0, il1, ii0, ii1, ij0, ij1
    REAL(KIND = wp) :: zl, zi
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:) :: zup, zvp
    IF (ln_timing) CALL timing_start('dta_uvd')
    CALL fld_read(kt, 1, sf_uvd)
    puvd(:, :, :, 1) = sf_uvd(1) % fnow(:, :, :)
    puvd(:, :, :, 2) = sf_uvd(2) % fnow(:, :, :)
    IF (ln_sco) THEN
      ALLOCATE(zup(jpk), zvp(jpk))
      IF (kt == nit000 .AND. lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) 'dta_uvd: interpolate U & V current data onto the s- or mixed s-z-coordinate mesh'
      END IF
      DO jj = 1, jpj
        DO ji = 1, jpi
          DO jk = 1, jpk
            zl = gdept_n(ji, jj, jk)
            IF (zl < gdept_1d(1)) THEN
              zup(jk) = puvd(ji, jj, 1, 1)
              zvp(jk) = puvd(ji, jj, 1, 2)
            ELSE IF (zl > gdept_1d(jpk)) THEN
              zup(jk) = puvd(ji, jj, jpkm1, 1)
              zvp(jk) = puvd(ji, jj, jpkm1, 2)
            ELSE
              DO jkk = 1, jpkm1
                IF ((zl - gdept_1d(jkk)) * (zl - gdept_1d(jkk + 1)) <= 0._wp) THEN
                  zi = (zl - gdept_1d(jkk)) / (gdept_1d(jkk + 1) - gdept_1d(jkk))
                  zup(jk) = puvd(ji, jj, jkk, 1) + (puvd(ji, jj, jkk + 1, 1) - puvd(ji, jj, jkk, 1)) * zi
                  zvp(jk) = puvd(ji, jj, jkk, 2) + (puvd(ji, jj, jkk + 1, 2) - puvd(ji, jj, jkk, 2)) * zi
                END IF
              END DO
            END IF
          END DO
          DO jk = 1, jpkm1
            puvd(ji, jj, jk, 1) = zup(jk) * umask(ji, jj, jk)
            puvd(ji, jj, jk, 2) = zvp(jk) * vmask(ji, jj, jk)
          END DO
          puvd(ji, jj, jpk, 1) = 0._wp
          puvd(ji, jj, jpk, 2) = 0._wp
        END DO
      END DO
      DEALLOCATE(zup, zvp)
    ELSE
      puvd(:, :, :, 1) = puvd(:, :, :, 1) * umask(:, :, :)
      puvd(:, :, :, 2) = puvd(:, :, :, 2) * vmask(:, :, :)
      IF (ln_zps) THEN
        DO jj = 1, jpj
          DO ji = 1, jpi
            ik = mbkt(ji, jj)
            IF (ik > 1) THEN
              zl = (gdept_1d(ik) - gdept_0(ji, jj, ik)) / (gdept_1d(ik) - gdept_1d(ik - 1))
              puvd(ji, jj, ik, 1) = (1. - zl) * puvd(ji, jj, ik, 1) + zl * puvd(ji, jj, ik - 1, 1)
              puvd(ji, jj, ik, 2) = (1. - zl) * puvd(ji, jj, ik, 2) + zl * puvd(ji, jj, ik - 1, 2)
            END IF
          END DO
        END DO
      END IF
    END IF
    IF (.NOT. ln_uvd_dyndmp) THEN
      IF (lwp) WRITE(numout, FMT = *) 'dta_uvd: deallocate U & V current arrays as they are only used to initialize the run'
      DEALLOCATE(sf_uvd(1) % fnow)
      IF (sf_uvd(1) % ln_tint) DEALLOCATE(sf_uvd(1) % fdta)
      DEALLOCATE(sf_uvd(2) % fnow)
      IF (sf_uvd(2) % ln_tint) DEALLOCATE(sf_uvd(2) % fdta)
      DEALLOCATE(sf_uvd)
    END IF
    IF (ln_timing) CALL timing_stop('dta_uvd')
  END SUBROUTINE dta_uvd
END MODULE dtauvd