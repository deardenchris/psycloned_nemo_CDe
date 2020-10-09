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
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
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
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data7
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data8
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data9
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data10
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data11
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data12
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data13
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data14
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data15
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data16
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data17
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data18
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data19
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data20
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data21
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data22
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data23
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data24
    IF (ln_timing) CALL timing_start('dia_ptr')
    IF (PRESENT(pvtr)) THEN
      IF (iom_use("zomsfglo")) THEN
        CALL profile_psy_data0 % PreStart('dia_ptr', 'r0', 0, 0)
        z3d(1, :, :) = ptr_sjk(pvtr(:, :, :))
        CALL profile_psy_data0 % PostEnd
        !$ACC KERNELS
        DO jk = 2, jpkm1
          z3d(1, :, jk) = z3d(1, :, jk - 1) + z3d(1, :, jk)
        END DO
        DO ji = 1, jpi
          z3d(ji, :, :) = z3d(1, :, :)
        END DO
        !$ACC END KERNELS
        CALL profile_psy_data1 % PreStart('dia_ptr', 'r1', 0, 0)
        cl1 = TRIM('zomsf' // clsubb(1))
        CALL iom_put(cl1, z3d * rc_sv)
        CALL profile_psy_data1 % PostEnd
        DO jn = 2, nptr
          CALL profile_psy_data2 % PreStart('dia_ptr', 'r2', 0, 0)
          z3d(1, :, :) = ptr_sjk(pvtr(:, :, :), btmsk(:, :, jn) * btm30(:, :))
          CALL profile_psy_data2 % PostEnd
          !$ACC KERNELS
          DO jk = 2, jpkm1
            z3d(1, :, jk) = z3d(1, :, jk - 1) + z3d(1, :, jk)
          END DO
          DO ji = 1, jpi
            z3d(ji, :, :) = z3d(1, :, :)
          END DO
          !$ACC END KERNELS
          CALL profile_psy_data3 % PreStart('dia_ptr', 'r3', 0, 0)
          cl1 = TRIM('zomsf' // clsubb(jn))
          CALL iom_put(cl1, z3d * rc_sv)
          CALL profile_psy_data3 % PostEnd
        END DO
      END IF
      IF (iom_use("sopstove") .OR. iom_use("sophtove") .OR. iom_use("sopstbtr") .OR. iom_use("sophtbtr")) THEN
        !$ACC KERNELS
        zmask(:, :, :) = 0._wp
        zts(:, :, :, :) = 0._wp
        zvn(:, :, :) = 0._wp
        DO jk = 1, jpkm1
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
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
        !$ACC END KERNELS
      END IF
      IF (iom_use("sopstove") .OR. iom_use("sophtove")) THEN
        CALL profile_psy_data4 % PreStart('dia_ptr', 'r4', 0, 0)
        sjk(:, :, 1) = ptr_sjk(zmask(:, :, :), btmsk(:, :, 1))
        CALL profile_psy_data4 % PostEnd
        !$ACC KERNELS
        r1_sjk(:, :, 1) = 0._wp
        WHERE (sjk(:, :, 1) /= 0._wp) r1_sjk(:, :, 1) = 1._wp / sjk(:, :, 1)
        !$ACC END KERNELS
        CALL profile_psy_data5 % PreStart('dia_ptr', 'r5', 0, 0)
        tn_jk(:, :, 1) = ptr_sjk(zts(:, :, :, jp_tem)) * r1_sjk(:, :, 1)
        sn_jk(:, :, 1) = ptr_sjk(zts(:, :, :, jp_sal)) * r1_sjk(:, :, 1)
        v_msf(:, :, 1) = ptr_sjk(zvn(:, :, :))
        CALL profile_psy_data5 % PostEnd
        !$ACC KERNELS
        htr_ove(:, 1) = SUM(v_msf(:, :, 1) * tn_jk(:, :, 1), 2)
        str_ove(:, 1) = SUM(v_msf(:, :, 1) * sn_jk(:, :, 1), 2)
        z2d(1, :) = htr_ove(:, 1) * rc_pwatt
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sophtove'
        !$ACC END KERNELS
        CALL iom_put(TRIM(cl1), z2d)
        !$ACC KERNELS
        z2d(1, :) = str_ove(:, 1) * rc_ggram
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sopstove'
        !$ACC END KERNELS
        CALL iom_put(TRIM(cl1), z2d)
        IF (ln_subbas) THEN
          DO jn = 2, nptr
            CALL profile_psy_data6 % PreStart('dia_ptr', 'r6', 0, 0)
            sjk(:, :, jn) = ptr_sjk(zmask(:, :, :), btmsk(:, :, jn))
            CALL profile_psy_data6 % PostEnd
            !$ACC KERNELS
            r1_sjk(:, :, jn) = 0._wp
            WHERE (sjk(:, :, jn) /= 0._wp) r1_sjk(:, :, jn) = 1._wp / sjk(:, :, jn)
            !$ACC END KERNELS
            CALL profile_psy_data7 % PreStart('dia_ptr', 'r7', 0, 0)
            tn_jk(:, :, jn) = ptr_sjk(zts(:, :, :, jp_tem), btmsk(:, :, jn)) * r1_sjk(:, :, jn)
            sn_jk(:, :, jn) = ptr_sjk(zts(:, :, :, jp_sal), btmsk(:, :, jn)) * r1_sjk(:, :, jn)
            v_msf(:, :, jn) = ptr_sjk(zvn(:, :, :), btmsk(:, :, jn))
            CALL profile_psy_data7 % PostEnd
            !$ACC KERNELS
            htr_ove(:, jn) = SUM(v_msf(:, :, jn) * tn_jk(:, :, jn), 2)
            str_ove(:, jn) = SUM(v_msf(:, :, jn) * sn_jk(:, :, jn), 2)
            z2d(1, :) = htr_ove(:, jn) * rc_pwatt
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL profile_psy_data8 % PreStart('dia_ptr', 'r8', 0, 0)
            cl1 = TRIM('sophtove_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL profile_psy_data8 % PostEnd
            !$ACC KERNELS
            z2d(1, :) = str_ove(:, jn) * rc_ggram
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL profile_psy_data9 % PreStart('dia_ptr', 'r9', 0, 0)
            cl1 = TRIM('sopstove_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL profile_psy_data9 % PostEnd
          END DO
        END IF
      END IF
      IF (iom_use("sopstbtr") .OR. iom_use("sophtbtr")) THEN
        CALL profile_psy_data10 % PreStart('dia_ptr', 'r10', 0, 0)
        sjk(:, 1, 1) = ptr_sj(zmask(:, :, :), btmsk(:, :, 1))
        CALL profile_psy_data10 % PostEnd
        !$ACC KERNELS
        r1_sjk(:, 1, 1) = 0._wp
        WHERE (sjk(:, 1, 1) /= 0._wp) r1_sjk(:, 1, 1) = 1._wp / sjk(:, 1, 1)
        !$ACC END KERNELS
        CALL profile_psy_data11 % PreStart('dia_ptr', 'r11', 0, 0)
        vsum = ptr_sj(zvn(:, :, :), btmsk(:, :, 1))
        tssum(:, jp_tem) = ptr_sj(zts(:, :, :, jp_tem), btmsk(:, :, 1))
        tssum(:, jp_sal) = ptr_sj(zts(:, :, :, jp_sal), btmsk(:, :, 1))
        CALL profile_psy_data11 % PostEnd
        !$ACC KERNELS
        htr_btr(:, 1) = vsum * tssum(:, jp_tem) * r1_sjk(:, 1, 1)
        str_btr(:, 1) = vsum * tssum(:, jp_sal) * r1_sjk(:, 1, 1)
        z2d(1, :) = htr_btr(:, 1) * rc_pwatt
        DO ji = 2, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sophtbtr'
        !$ACC END KERNELS
        CALL iom_put(TRIM(cl1), z2d)
        !$ACC KERNELS
        z2d(1, :) = str_btr(:, 1) * rc_ggram
        DO ji = 2, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sopstbtr'
        !$ACC END KERNELS
        CALL iom_put(TRIM(cl1), z2d)
        IF (ln_subbas) THEN
          DO jn = 2, nptr
            CALL profile_psy_data12 % PreStart('dia_ptr', 'r12', 0, 0)
            sjk(:, 1, jn) = ptr_sj(zmask(:, :, :), btmsk(:, :, jn))
            CALL profile_psy_data12 % PostEnd
            !$ACC KERNELS
            r1_sjk(:, 1, jn) = 0._wp
            WHERE (sjk(:, 1, jn) /= 0._wp) r1_sjk(:, 1, jn) = 1._wp / sjk(:, 1, jn)
            !$ACC END KERNELS
            CALL profile_psy_data13 % PreStart('dia_ptr', 'r13', 0, 0)
            vsum = ptr_sj(zvn(:, :, :), btmsk(:, :, jn))
            tssum(:, jp_tem) = ptr_sj(zts(:, :, :, jp_tem), btmsk(:, :, jn))
            tssum(:, jp_sal) = ptr_sj(zts(:, :, :, jp_sal), btmsk(:, :, jn))
            CALL profile_psy_data13 % PostEnd
            !$ACC KERNELS
            htr_btr(:, jn) = vsum * tssum(:, jp_tem) * r1_sjk(:, 1, jn)
            str_btr(:, jn) = vsum * tssum(:, jp_sal) * r1_sjk(:, 1, jn)
            z2d(1, :) = htr_btr(:, jn) * rc_pwatt
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL profile_psy_data14 % PreStart('dia_ptr', 'r14', 0, 0)
            cl1 = TRIM('sophtbtr_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL profile_psy_data14 % PostEnd
            !$ACC KERNELS
            z2d(1, :) = str_btr(:, jn) * rc_ggram
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL profile_psy_data15 % PreStart('dia_ptr', 'r15', 0, 0)
            cl1 = TRIM('sopstbtr_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL profile_psy_data15 % PostEnd
          END DO
        END IF
      END IF
    ELSE
      IF (iom_use("zotemglo")) THEN
        !$ACC KERNELS
        DO jk = 1, jpkm1
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 1, jpj
            DO ji = 1, jpi
              zsfc = e1t(ji, jj) * e3t_n(ji, jj, jk)
              zmask(ji, jj, jk) = tmask(ji, jj, jk) * zsfc
              zts(ji, jj, jk, jp_tem) = tsn(ji, jj, jk, jp_tem) * zsfc
              zts(ji, jj, jk, jp_sal) = tsn(ji, jj, jk, jp_sal) * zsfc
            END DO
          END DO
        END DO
        !$ACC END KERNELS
        DO jn = 1, nptr
          CALL profile_psy_data16 % PreStart('dia_ptr', 'r16', 0, 0)
          zmask(1, :, :) = ptr_sjk(zmask(:, :, :), btmsk(:, :, jn))
          cl1 = TRIM('zosrf' // clsubb(jn))
          CALL iom_put(cl1, zmask)
          z3d(1, :, :) = ptr_sjk(zts(:, :, :, jp_tem), btmsk(:, :, jn)) / MAX(zmask(1, :, :), 10.E-15)
          CALL profile_psy_data16 % PostEnd
          !$ACC KERNELS
          DO ji = 1, jpi
            z3d(ji, :, :) = z3d(1, :, :)
          END DO
          !$ACC END KERNELS
          CALL profile_psy_data17 % PreStart('dia_ptr', 'r17', 0, 0)
          cl1 = TRIM('zotem' // clsubb(jn))
          CALL iom_put(cl1, z3d)
          z3d(1, :, :) = ptr_sjk(zts(:, :, :, jp_sal), btmsk(:, :, jn)) / MAX(zmask(1, :, :), 10.E-15)
          CALL profile_psy_data17 % PostEnd
          !$ACC KERNELS
          DO ji = 1, jpi
            z3d(ji, :, :) = z3d(1, :, :)
          END DO
          !$ACC END KERNELS
          CALL profile_psy_data18 % PreStart('dia_ptr', 'r18', 0, 0)
          cl1 = TRIM('zosal' // clsubb(jn))
          CALL iom_put(cl1, z3d)
          CALL profile_psy_data18 % PostEnd
        END DO
      END IF
      IF (iom_use("sophtadv") .OR. iom_use("sopstadv")) THEN
        !$ACC KERNELS
        z2d(1, :) = htr_adv(:, 1) * rc_pwatt
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sophtadv'
        !$ACC END KERNELS
        CALL iom_put(TRIM(cl1), z2d)
        !$ACC KERNELS
        z2d(1, :) = str_adv(:, 1) * rc_ggram
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sopstadv'
        !$ACC END KERNELS
        CALL iom_put(TRIM(cl1), z2d)
        IF (ln_subbas) THEN
          DO jn = 2, nptr
            !$ACC KERNELS
            z2d(1, :) = htr_adv(:, jn) * rc_pwatt
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL profile_psy_data19 % PreStart('dia_ptr', 'r19', 0, 0)
            cl1 = TRIM('sophtadv_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL profile_psy_data19 % PostEnd
            !$ACC KERNELS
            z2d(1, :) = str_adv(:, jn) * rc_ggram
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL profile_psy_data20 % PreStart('dia_ptr', 'r20', 0, 0)
            cl1 = TRIM('sopstadv_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL profile_psy_data20 % PostEnd
          END DO
        END IF
      END IF
      IF (iom_use("sophtldf") .OR. iom_use("sopstldf")) THEN
        !$ACC KERNELS
        z2d(1, :) = htr_ldf(:, 1) * rc_pwatt
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sophtldf'
        !$ACC END KERNELS
        CALL iom_put(TRIM(cl1), z2d)
        !$ACC KERNELS
        z2d(1, :) = str_ldf(:, 1) * rc_ggram
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sopstldf'
        !$ACC END KERNELS
        CALL iom_put(TRIM(cl1), z2d)
        IF (ln_subbas) THEN
          DO jn = 2, nptr
            !$ACC KERNELS
            z2d(1, :) = htr_ldf(:, jn) * rc_pwatt
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL profile_psy_data21 % PreStart('dia_ptr', 'r21', 0, 0)
            cl1 = TRIM('sophtldf_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL profile_psy_data21 % PostEnd
            !$ACC KERNELS
            z2d(1, :) = str_ldf(:, jn) * rc_ggram
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL profile_psy_data22 % PreStart('dia_ptr', 'r22', 0, 0)
            cl1 = TRIM('sopstldf_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL profile_psy_data22 % PostEnd
          END DO
        END IF
      END IF
      IF (iom_use("sophteiv") .OR. iom_use("sopsteiv")) THEN
        !$ACC KERNELS
        z2d(1, :) = htr_eiv(:, 1) * rc_pwatt
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sophteiv'
        !$ACC END KERNELS
        CALL iom_put(TRIM(cl1), z2d)
        !$ACC KERNELS
        z2d(1, :) = str_eiv(:, 1) * rc_ggram
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        cl1 = 'sopsteiv'
        !$ACC END KERNELS
        CALL iom_put(TRIM(cl1), z2d)
        IF (ln_subbas) THEN
          DO jn = 2, nptr
            !$ACC KERNELS
            z2d(1, :) = htr_eiv(:, jn) * rc_pwatt
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL profile_psy_data23 % PreStart('dia_ptr', 'r23', 0, 0)
            cl1 = TRIM('sophteiv_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL profile_psy_data23 % PostEnd
            !$ACC KERNELS
            z2d(1, :) = str_eiv(:, jn) * rc_ggram
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL profile_psy_data24 % PreStart('dia_ptr', 'r24', 0, 0)
            cl1 = TRIM('sopsteiv_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL profile_psy_data24 % PostEnd
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
        !$ACC KERNELS
        btmsk(:, :, 5) = MAX(btmsk(:, :, 3), btmsk(:, :, 4))
        WHERE (gphit(:, :) < - 30._wp)
          btm30(:, :) = 0._wp
        ELSEWHERE
          btm30(:, :) = ssmask(:, :)
        END WHERE
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
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
      !$ACC END KERNELS
    END IF
  END SUBROUTINE dia_ptr_init
  SUBROUTINE dia_ptr_hst(ktra, cptr, pva)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: ktra
    CHARACTER(LEN = 3), INTENT(IN) :: cptr
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pva
    INTEGER :: jn
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('dia_ptr_hst', 'r0', 0, 0)
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
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE dia_ptr_hst
  FUNCTION dia_ptr_alloc()
    INTEGER :: dia_ptr_alloc
    INTEGER, DIMENSION(3) :: ierr
    !$ACC KERNELS
    ierr(:) = 0
    !$ACC END KERNELS
    ALLOCATE(btmsk(jpi, jpj, nptr), htr_adv(jpj, nptr), str_adv(jpj, nptr), htr_eiv(jpj, nptr), str_eiv(jpj, nptr), htr_ove(jpj, &
&nptr), str_ove(jpj, nptr), htr_btr(jpj, nptr), str_btr(jpj, nptr), htr_ldf(jpj, nptr), str_ldf(jpj, nptr), STAT = ierr(1))
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
    !$ACC KERNELS
    ijpj = jpj
    p_fval(:) = 0._wp
    !$ACC END KERNELS
    IF (PRESENT(pmsk)) THEN
      !$ACC KERNELS
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            p_fval(jj) = p_fval(jj) + pva(ji, jj, jk) * tmask_i(ji, jj) * pmsk(ji, jj)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            p_fval(jj) = p_fval(jj) + pva(ji, jj, jk) * tmask_i(ji, jj)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END IF
  END FUNCTION ptr_sj_3d
  FUNCTION ptr_sj_2d(pva, pmsk) RESULT(p_fval)
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj) :: pva
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj), OPTIONAL :: pmsk
    INTEGER :: ji, jj
    INTEGER :: ijpj
    REAL(KIND = wp), POINTER, DIMENSION(:) :: p_fval
    p_fval => p_fval1d
    !$ACC KERNELS
    ijpj = jpj
    p_fval(:) = 0._wp
    !$ACC END KERNELS
    IF (PRESENT(pmsk)) THEN
      !$ACC KERNELS
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = nldi, nlei
          p_fval(jj) = p_fval(jj) + pva(ji, jj) * tmask_i(ji, jj) * pmsk(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = nldi, nlei
          p_fval(jj) = p_fval(jj) + pva(ji, jj) * tmask_i(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
    END IF
  END FUNCTION ptr_sj_2d
  FUNCTION ptr_sjk(pta, pmsk) RESULT(p_fval)
    IMPLICIT NONE
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj, jpk) :: pta
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj), OPTIONAL :: pmsk
    INTEGER :: ji, jj, jk
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: p_fval
    p_fval => p_fval2d
    !$ACC KERNELS
    p_fval(:, :) = 0._wp
    !$ACC END KERNELS
    IF (PRESENT(pmsk)) THEN
      !$ACC KERNELS
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = nldi, nlei
            p_fval(jj, jk) = p_fval(jj, jk) + pta(ji, jj, jk) * pmsk(ji, jj)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = nldi, nlei
            p_fval(jj, jk) = p_fval(jj, jk) + pta(ji, jj, jk) * tmask_i(ji, jj)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END IF
  END FUNCTION ptr_sjk
END MODULE diaptr