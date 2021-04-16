MODULE diaptr
  USE oce
  USE dom_oce
  USE phycst
  USE iom
  USE in_out_manager
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  INTERFACE ptr_sj
    MODULE PROCEDURE ptr_sj_3d, ptr_sj_2d
  END INTERFACE
  PUBLIC :: ptr_sj
  PUBLIC :: ptr_sjk
  PUBLIC :: dia_ptr_init
  PUBLIC :: dia_ptr
  PUBLIC :: dia_ptr_hst
  REAL(KIND = wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:, :) :: htr_adv, htr_ldf, htr_eiv
  REAL(KIND = wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:, :) :: str_adv, str_ldf, str_eiv
  REAL(KIND = wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:, :) :: htr_ove, str_ove
  REAL(KIND = wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:, :) :: htr_btr, str_btr
  LOGICAL, PUBLIC :: ln_diaptr
  LOGICAL, PUBLIC :: ln_subbas
  INTEGER, PUBLIC :: nptr
  REAL(KIND = wp) :: rc_sv = 1.E-6_wp
  REAL(KIND = wp) :: rc_pwatt = 1.E-15_wp
  REAL(KIND = wp) :: rc_ggram = 1.E-6_wp
  CHARACTER(LEN = 3), ALLOCATABLE, SAVE, DIMENSION(:) :: clsubb
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: btmsk
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: btm30
  REAL(KIND = wp), TARGET, ALLOCATABLE, SAVE, DIMENSION(:) :: p_fval1d
  REAL(KIND = wp), TARGET, ALLOCATABLE, SAVE, DIMENSION(:, :) :: p_fval2d
  CONTAINS
  SUBROUTINE dia_ptr(pvtr)
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN), OPTIONAL :: pvtr
    INTEGER :: ji, jj, jk, jn
    REAL(KIND = wp) :: zsfc, zvfc
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z2d
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: z3d
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zmask
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, jpts) :: zts
    REAL(KIND = wp), DIMENSION(jpj) :: vsum
    REAL(KIND = wp), DIMENSION(jpj, jpts) :: tssum
    REAL(KIND = wp), DIMENSION(jpj, jpk, nptr) :: sjk, r1_sjk
    REAL(KIND = wp), DIMENSION(jpj, jpk, nptr) :: v_msf, sn_jk, tn_jk
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zvn
    CHARACTER(LEN = 12) :: cl1
    IF (ln_timing) CALL timing_start('dia_ptr')
    IF (PRESENT(pvtr)) THEN
      IF (iom_use("zomsfglo")) THEN
        z3d(1, :, :) = ptr_sjk(pvtr(:, :, :))
        ! !$OMP parallel default(shared), private(jk)
        ! !$OMP do schedule(static) ! CDe race condition
        DO jk = 2, jpkm1
          z3d(1, :, jk) = z3d(1, :, jk - 1) + z3d(1, :, jk)
        END DO
        ! !$OMP end do
        ! !$OMP end parallel
        DO ji = 1, jpi
          z3d(ji, :, :) = z3d(1, :, :)
        END DO
        cl1 = TRIM('zomsf' // clsubb(1))
        CALL iom_put(cl1, z3d * rc_sv)
        DO jn = 2, nptr
          z3d(1, :, :) = ptr_sjk(pvtr(:, :, :), btmsk(:, :, jn) * btm30(:, :))
          ! !$OMP parallel default(shared), private(jk)
          ! !$OMP do schedule(static) ! CDe race condition
          DO jk = 2, jpkm1
            z3d(1, :, jk) = z3d(1, :, jk - 1) + z3d(1, :, jk)
          END DO
          ! !$OMP end do
          ! !$OMP end parallel
          DO ji = 1, jpi
            z3d(ji, :, :) = z3d(1, :, :)
          END DO
          cl1 = TRIM('zomsf' // clsubb(jn))
          CALL iom_put(cl1, z3d * rc_sv)
        END DO
      END IF
      IF (iom_use("sopstove") .OR. iom_use("sophtove") .OR. iom_use("sopstbtr") .OR. iom_use("sophtbtr")) THEN
        zmask(:, :, :) = 0._wp
        zts(:, :, :, :) = 0._wp
        zvn(:, :, :) = 0._wp
        !$OMP parallel default(shared), private(ji,jj,jk,zvfc)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 1, jpjm1
            DO ji = 1, jpi
              zvfc = e1v(ji, jj) * e3v_n(ji, jj, jk)
              zmask(ji, jj, jk) = vmask(ji, jj, jk) * zvfc
              zts(ji, jj, jk, jp_tem) = (tsn(ji, jj, jk, jp_tem) + tsn(ji, jj + 1, jk, jp_tem)) * 0.5 * zvfc
              zts(ji, jj, jk, jp_sal) = (tsn(ji, jj, jk, jp_sal) + tsn(ji, jj + 1, jk, jp_sal)) * 0.5 * zvfc
              zvn(ji, jj, jk) = vn(ji, jj, jk) * zvfc
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
      END IF
      IF (iom_use("sopstove") .OR. iom_use("sophtove")) THEN
        sjk(:, :, 1) = ptr_sjk(zmask(:, :, :), btmsk(:, :, 1))
        r1_sjk(:, :, 1) = 0._wp
        WHERE (sjk(:, :, 1) /= 0._wp) r1_sjk(:, :, 1) = 1._wp / sjk(:, :, 1)
        tn_jk(:, :, 1) = ptr_sjk(zts(:, :, :, jp_tem)) * r1_sjk(:, :, 1)
        sn_jk(:, :, 1) = ptr_sjk(zts(:, :, :, jp_sal)) * r1_sjk(:, :, 1)
        v_msf(:, :, 1) = ptr_sjk(zvn(:, :, :))
        htr_ove(:, 1) = SUM(v_msf(:, :, 1) * tn_jk(:, :, 1), 2)
        str_ove(:, 1) = SUM(v_msf(:, :, 1) * sn_jk(:, :, 1), 2)
        z2d(1, :) = htr_ove(:, 1) * rc_pwatt
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sophtove'
        CALL iom_put(TRIM(cl1), z2d)
        z2d(1, :) = str_ove(:, 1) * rc_ggram
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sopstove'
        CALL iom_put(TRIM(cl1), z2d)
        IF (ln_subbas) THEN
          DO jn = 2, nptr
            sjk(:, :, jn) = ptr_sjk(zmask(:, :, :), btmsk(:, :, jn))
            r1_sjk(:, :, jn) = 0._wp
            WHERE (sjk(:, :, jn) /= 0._wp) r1_sjk(:, :, jn) = 1._wp / sjk(:, :, jn)
            tn_jk(:, :, jn) = ptr_sjk(zts(:, :, :, jp_tem), btmsk(:, :, jn)) * r1_sjk(:, :, jn)
            sn_jk(:, :, jn) = ptr_sjk(zts(:, :, :, jp_sal), btmsk(:, :, jn)) * r1_sjk(:, :, jn)
            v_msf(:, :, jn) = ptr_sjk(zvn(:, :, :), btmsk(:, :, jn))
            htr_ove(:, jn) = SUM(v_msf(:, :, jn) * tn_jk(:, :, jn), 2)
            str_ove(:, jn) = SUM(v_msf(:, :, jn) * sn_jk(:, :, jn), 2)
            z2d(1, :) = htr_ove(:, jn) * rc_pwatt
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            cl1 = TRIM('sophtove_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            z2d(1, :) = str_ove(:, jn) * rc_ggram
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            cl1 = TRIM('sopstove_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
          END DO
        END IF
      END IF
      IF (iom_use("sopstbtr") .OR. iom_use("sophtbtr")) THEN
        sjk(:, 1, 1) = ptr_sj(zmask(:, :, :), btmsk(:, :, 1))
        r1_sjk(:, 1, 1) = 0._wp
        WHERE (sjk(:, 1, 1) /= 0._wp) r1_sjk(:, 1, 1) = 1._wp / sjk(:, 1, 1)
        vsum = ptr_sj(zvn(:, :, :), btmsk(:, :, 1))
        tssum(:, jp_tem) = ptr_sj(zts(:, :, :, jp_tem), btmsk(:, :, 1))
        tssum(:, jp_sal) = ptr_sj(zts(:, :, :, jp_sal), btmsk(:, :, 1))
        htr_btr(:, 1) = vsum * tssum(:, jp_tem) * r1_sjk(:, 1, 1)
        str_btr(:, 1) = vsum * tssum(:, jp_sal) * r1_sjk(:, 1, 1)
        z2d(1, :) = htr_btr(:, 1) * rc_pwatt
        DO ji = 2, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sophtbtr'
        CALL iom_put(TRIM(cl1), z2d)
        z2d(1, :) = str_btr(:, 1) * rc_ggram
        DO ji = 2, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sopstbtr'
        CALL iom_put(TRIM(cl1), z2d)
        IF (ln_subbas) THEN
          DO jn = 2, nptr
            sjk(:, 1, jn) = ptr_sj(zmask(:, :, :), btmsk(:, :, jn))
            r1_sjk(:, 1, jn) = 0._wp
            WHERE (sjk(:, 1, jn) /= 0._wp) r1_sjk(:, 1, jn) = 1._wp / sjk(:, 1, jn)
            vsum = ptr_sj(zvn(:, :, :), btmsk(:, :, jn))
            tssum(:, jp_tem) = ptr_sj(zts(:, :, :, jp_tem), btmsk(:, :, jn))
            tssum(:, jp_sal) = ptr_sj(zts(:, :, :, jp_sal), btmsk(:, :, jn))
            htr_btr(:, jn) = vsum * tssum(:, jp_tem) * r1_sjk(:, 1, jn)
            str_btr(:, jn) = vsum * tssum(:, jp_sal) * r1_sjk(:, 1, jn)
            z2d(1, :) = htr_btr(:, jn) * rc_pwatt
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            cl1 = TRIM('sophtbtr_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            z2d(1, :) = str_btr(:, jn) * rc_ggram
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            cl1 = TRIM('sopstbtr_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
          END DO
        END IF
      END IF
    ELSE
      IF (iom_use("zotemglo")) THEN
        !$OMP parallel default(shared), private(ji,jj,jk,zsfc)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 1, jpj
            DO ji = 1, jpi
              zsfc = e1t(ji, jj) * e3t_n(ji, jj, jk)
              zmask(ji, jj, jk) = tmask(ji, jj, jk) * zsfc
              zts(ji, jj, jk, jp_tem) = tsn(ji, jj, jk, jp_tem) * zsfc
              zts(ji, jj, jk, jp_sal) = tsn(ji, jj, jk, jp_sal) * zsfc
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
        DO jn = 1, nptr
          zmask(1, :, :) = ptr_sjk(zmask(:, :, :), btmsk(:, :, jn))
          cl1 = TRIM('zosrf' // clsubb(jn))
          CALL iom_put(cl1, zmask)
          z3d(1, :, :) = ptr_sjk(zts(:, :, :, jp_tem), btmsk(:, :, jn)) / MAX(zmask(1, :, :), 10.E-15)
          DO ji = 1, jpi
            z3d(ji, :, :) = z3d(1, :, :)
          END DO
          cl1 = TRIM('zotem' // clsubb(jn))
          CALL iom_put(cl1, z3d)
          z3d(1, :, :) = ptr_sjk(zts(:, :, :, jp_sal), btmsk(:, :, jn)) / MAX(zmask(1, :, :), 10.E-15)
          DO ji = 1, jpi
            z3d(ji, :, :) = z3d(1, :, :)
          END DO
          cl1 = TRIM('zosal' // clsubb(jn))
          CALL iom_put(cl1, z3d)
        END DO
      END IF
      IF (iom_use("sophtadv") .OR. iom_use("sopstadv")) THEN
        z2d(1, :) = htr_adv(:, 1) * rc_pwatt
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sophtadv'
        CALL iom_put(TRIM(cl1), z2d)
        z2d(1, :) = str_adv(:, 1) * rc_ggram
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sopstadv'
        CALL iom_put(TRIM(cl1), z2d)
        IF (ln_subbas) THEN
          DO jn = 2, nptr
            z2d(1, :) = htr_adv(:, jn) * rc_pwatt
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            cl1 = TRIM('sophtadv_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            z2d(1, :) = str_adv(:, jn) * rc_ggram
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            cl1 = TRIM('sopstadv_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
          END DO
        END IF
      END IF
      IF (iom_use("sophtldf") .OR. iom_use("sopstldf")) THEN
        z2d(1, :) = htr_ldf(:, 1) * rc_pwatt
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sophtldf'
        CALL iom_put(TRIM(cl1), z2d)
        z2d(1, :) = str_ldf(:, 1) * rc_ggram
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sopstldf'
        CALL iom_put(TRIM(cl1), z2d)
        IF (ln_subbas) THEN
          DO jn = 2, nptr
            z2d(1, :) = htr_ldf(:, jn) * rc_pwatt
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            cl1 = TRIM('sophtldf_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            z2d(1, :) = str_ldf(:, jn) * rc_ggram
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            cl1 = TRIM('sopstldf_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
          END DO
        END IF
      END IF
      IF (iom_use("sophteiv") .OR. iom_use("sopsteiv")) THEN
        z2d(1, :) = htr_eiv(:, 1) * rc_pwatt
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sophteiv'
        CALL iom_put(TRIM(cl1), z2d)
        z2d(1, :) = str_eiv(:, 1) * rc_ggram
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sopsteiv'
        CALL iom_put(TRIM(cl1), z2d)
        IF (ln_subbas) THEN
          DO jn = 2, nptr
            z2d(1, :) = htr_eiv(:, jn) * rc_pwatt
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            cl1 = TRIM('sophteiv_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            z2d(1, :) = str_eiv(:, jn) * rc_ggram
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            cl1 = TRIM('sopsteiv_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
          END DO
        END IF
      END IF
    END IF
    IF (ln_timing) CALL timing_stop('dia_ptr')
  END SUBROUTINE dia_ptr
  SUBROUTINE dia_ptr_init
    INTEGER :: jn
    INTEGER :: inum, ierr
    INTEGER :: ios
    NAMELIST /namptr/ ln_diaptr, ln_subbas
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namptr, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namptr in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namptr, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namptr in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namptr)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dia_ptr_init : poleward transport and msf initialization'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namptr : set ptr parameters'
      WRITE(numout, FMT = *) '      Poleward heat & salt transport (T) or not (F)      ln_diaptr  = ', ln_diaptr
      WRITE(numout, FMT = *) '      Global (F) or glo/Atl/Pac/Ind/Indo-Pac basins      ln_subbas  = ', ln_subbas
    END IF
    IF (ln_diaptr) THEN
      IF (ln_subbas) THEN
        nptr = 5
        ALLOCATE(clsubb(nptr))
        clsubb(1) = 'glo'
        clsubb(2) = 'atl'
        clsubb(3) = 'pac'
        clsubb(4) = 'ind'
        clsubb(5) = 'ipc'
      ELSE
        nptr = 1
        ALLOCATE(clsubb(nptr))
        clsubb(1) = 'glo'
      END IF
      IF (dia_ptr_alloc() /= 0) CALL ctl_stop('STOP', 'dia_ptr_init : unable to allocate arrays')
      rc_pwatt = rc_pwatt * rau0_rcp
      IF (lk_mpp) CALL mpp_ini_znl(numout)
      IF (ln_subbas) THEN
        CALL iom_open('subbasins', inum, ldstop = .FALSE.)
        CALL iom_get(inum, jpdom_data, 'atlmsk', btmsk(:, :, 2))
        CALL iom_get(inum, jpdom_data, 'pacmsk', btmsk(:, :, 3))
        CALL iom_get(inum, jpdom_data, 'indmsk', btmsk(:, :, 4))
        CALL iom_close(inum)
        btmsk(:, :, 5) = MAX(btmsk(:, :, 3), btmsk(:, :, 4))
        WHERE (gphit(:, :) < - 30._wp)
          btm30(:, :) = 0._wp
        ELSEWHERE
          btm30(:, :) = ssmask(:, :)
        END WHERE
      END IF
      btmsk(:, :, 1) = tmask_i(:, :)
      DO jn = 1, nptr
        btmsk(:, :, jn) = btmsk(:, :, jn) * tmask_i(:, :)
      END DO
      htr_adv(:, :) = 0._wp
      str_adv(:, :) = 0._wp
      htr_ldf(:, :) = 0._wp
      str_ldf(:, :) = 0._wp
      htr_eiv(:, :) = 0._wp
      str_eiv(:, :) = 0._wp
      htr_ove(:, :) = 0._wp
      str_ove(:, :) = 0._wp
      htr_btr(:, :) = 0._wp
      str_btr(:, :) = 0._wp
    END IF
  END SUBROUTINE dia_ptr_init
  SUBROUTINE dia_ptr_hst(ktra, cptr, pva)
    INTEGER, INTENT(IN) :: ktra
    CHARACTER(LEN = 3), INTENT(IN) :: cptr
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pva
    INTEGER :: jn
    IF (cptr == 'adv') THEN
      IF (ktra == jp_tem) htr_adv(:, 1) = ptr_sj(pva(:, :, :))
      IF (ktra == jp_sal) str_adv(:, 1) = ptr_sj(pva(:, :, :))
    END IF
    IF (cptr == 'ldf') THEN
      IF (ktra == jp_tem) htr_ldf(:, 1) = ptr_sj(pva(:, :, :))
      IF (ktra == jp_sal) str_ldf(:, 1) = ptr_sj(pva(:, :, :))
    END IF
    IF (cptr == 'eiv') THEN
      IF (ktra == jp_tem) htr_eiv(:, 1) = ptr_sj(pva(:, :, :))
      IF (ktra == jp_sal) str_eiv(:, 1) = ptr_sj(pva(:, :, :))
    END IF
    IF (ln_subbas) THEN
      IF (cptr == 'adv') THEN
        IF (ktra == jp_tem) THEN
          DO jn = 2, nptr
            htr_adv(:, jn) = ptr_sj(pva(:, :, :), btmsk(:, :, jn))
          END DO
        END IF
        IF (ktra == jp_sal) THEN
          DO jn = 2, nptr
            str_adv(:, jn) = ptr_sj(pva(:, :, :), btmsk(:, :, jn))
          END DO
        END IF
      END IF
      IF (cptr == 'ldf') THEN
        IF (ktra == jp_tem) THEN
          DO jn = 2, nptr
            htr_ldf(:, jn) = ptr_sj(pva(:, :, :), btmsk(:, :, jn))
          END DO
        END IF
        IF (ktra == jp_sal) THEN
          DO jn = 2, nptr
            str_ldf(:, jn) = ptr_sj(pva(:, :, :), btmsk(:, :, jn))
          END DO
        END IF
      END IF
      IF (cptr == 'eiv') THEN
        IF (ktra == jp_tem) THEN
          DO jn = 2, nptr
            htr_eiv(:, jn) = ptr_sj(pva(:, :, :), btmsk(:, :, jn))
          END DO
        END IF
        IF (ktra == jp_sal) THEN
          DO jn = 2, nptr
            str_eiv(:, jn) = ptr_sj(pva(:, :, :), btmsk(:, :, jn))
          END DO
        END IF
      END IF
    END IF
  END SUBROUTINE dia_ptr_hst
  FUNCTION dia_ptr_alloc()
    INTEGER :: dia_ptr_alloc
    INTEGER, DIMENSION(3) :: ierr
    ierr(:) = 0
    ALLOCATE(btmsk(jpi, jpj, nptr), htr_adv(jpj, nptr), str_adv(jpj, nptr), htr_eiv(jpj, nptr), str_eiv(jpj, nptr), htr_ove(jpj, nptr), str_ove(jpj, nptr), htr_btr(jpj, nptr), str_btr(jpj, nptr), htr_ldf(jpj, nptr), str_ldf(jpj, nptr), STAT = ierr(1))
    ALLOCATE(p_fval1d(jpj), p_fval2d(jpj, jpk), STAT = ierr(2))
    ALLOCATE(btm30(jpi, jpj), STAT = ierr(3))
    dia_ptr_alloc = MAXVAL(ierr)
    CALL mpp_sum('diaptr', dia_ptr_alloc)
  END FUNCTION dia_ptr_alloc
  FUNCTION ptr_sj_3d(pva, pmsk) RESULT(p_fval)
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj, jpk) :: pva
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj), OPTIONAL :: pmsk
    INTEGER :: ji, jj, jk
    INTEGER :: ijpj
    REAL(KIND = wp), POINTER, DIMENSION(:) :: p_fval
    p_fval => p_fval1d
    ijpj = jpj
    p_fval(:) = 0._wp
    IF (PRESENT(pmsk)) THEN
      ! !$OMP parallel default(shared), private(ji,jj,jk)
      ! !$OMP do schedule(static) ! CDe race condition
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            p_fval(jj) = p_fval(jj) + pva(ji, jj, jk) * tmask_i(ji, jj) * pmsk(ji, jj)
          END DO
        END DO
      END DO
      ! !$OMP end do
      ! !$OMP end parallel
    ELSE
      ! !$OMP parallel default(shared), private(ji,jj,jk)
      ! !$OMP do schedule(static) ! CDe race condition
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            p_fval(jj) = p_fval(jj) + pva(ji, jj, jk) * tmask_i(ji, jj)
          END DO
        END DO
      END DO
      ! !$OMP end do
      ! !$OMP end parallel
    END IF
    CALL mpp_sum('diaptr', p_fval, ijpj, ncomm_znl)
  END FUNCTION ptr_sj_3d
  FUNCTION ptr_sj_2d(pva, pmsk) RESULT(p_fval)
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj) :: pva
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj), OPTIONAL :: pmsk
    INTEGER :: ji, jj
    INTEGER :: ijpj
    REAL(KIND = wp), POINTER, DIMENSION(:) :: p_fval
    p_fval => p_fval1d
    ijpj = jpj
    p_fval(:) = 0._wp
    IF (PRESENT(pmsk)) THEN
      DO jj = 2, jpjm1
        DO ji = nldi, nlei
          p_fval(jj) = p_fval(jj) + pva(ji, jj) * tmask_i(ji, jj) * pmsk(ji, jj)
        END DO
      END DO
    ELSE
      DO jj = 2, jpjm1
        DO ji = nldi, nlei
          p_fval(jj) = p_fval(jj) + pva(ji, jj) * tmask_i(ji, jj)
        END DO
      END DO
    END IF
    CALL mpp_sum('diaptr', p_fval, ijpj, ncomm_znl)
  END FUNCTION ptr_sj_2d
  FUNCTION ptr_sjk(pta, pmsk) RESULT(p_fval)
    IMPLICIT NONE
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj, jpk) :: pta
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj), OPTIONAL :: pmsk
    INTEGER :: ji, jj, jk
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: p_fval
    INTEGER, DIMENSION(1) :: ish
    INTEGER, DIMENSION(2) :: ish2
    INTEGER :: ijpjjpk
    REAL(KIND = wp), DIMENSION(jpj * jpk) :: zwork
    p_fval => p_fval2d
    p_fval(:, :) = 0._wp
    IF (PRESENT(pmsk)) THEN
      !$OMP parallel default(shared), private(ji,jj,jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = nldi, nlei
            p_fval(jj, jk) = p_fval(jj, jk) + pta(ji, jj, jk) * pmsk(ji, jj)
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
    ELSE
      !$OMP parallel default(shared), private(ji,jj,jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = nldi, nlei
            p_fval(jj, jk) = p_fval(jj, jk) + pta(ji, jj, jk) * tmask_i(ji, jj)
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
    END IF
    ijpjjpk = jpj * jpk
    ish(1) = ijpjjpk
    ish2(1) = jpj
    ish2(2) = jpk
    zwork(1 : ijpjjpk) = RESHAPE(p_fval, ish)
    CALL mpp_sum('diaptr', zwork, ijpjjpk, ncomm_znl)
    p_fval(:, :) = RESHAPE(zwork, ish2)
  END FUNCTION ptr_sjk
END MODULE diaptr
