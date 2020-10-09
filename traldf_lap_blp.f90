MODULE traldf_lap_blp
  USE oce
  USE dom_oce
  USE ldftra
  USE traldf_iso
  USE traldf_triad
  USE diaptr
  USE diaar5
  USE trc_oce
  USE zpshde
  USE in_out_manager
  USE iom
  USE lbclnk
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_ldf_lap
  PUBLIC :: tra_ldf_blp
  LOGICAL :: l_ptr
  LOGICAL :: l_hst
  CONTAINS
  SUBROUTINE tra_ldf_lap(kt, kit000, cdtype, pahu, pahv, pgu, pgv, pgui, pgvi, ptb, pta, kjpt, kpass)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kit000
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    INTEGER, INTENT(IN) :: kjpt
    INTEGER, INTENT(IN) :: kpass
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pahu, pahv
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt), INTENT(IN) :: pgu, pgv
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt), INTENT(IN) :: pgui, pgvi
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: ptb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    INTEGER :: ji, jj, jk, jn
    REAL(KIND = wp) :: zsign
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: ztu, ztv, zaheeu, zaheev
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('tra_ldf_lap', 'r0', 0, 0)
    IF (kt == nit000 .AND. lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'tra_ldf_lap : iso-level laplacian diffusion on ', cdtype, ', pass=', kpass
      WRITE(numout, FMT = *) '~~~~~~~~~~~ '
    END IF
    l_hst = .FALSE.
    l_ptr = .FALSE.
    IF (cdtype == 'TRA' .AND. ln_diaptr) l_ptr = .TRUE.
    IF (cdtype == 'TRA' .AND. (iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR. iom_use("uadv_salttr") .OR. &
&iom_use("vadv_salttr"))) l_hst = .TRUE.
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    IF (kpass == 1) THEN
      zsign = 1._wp
    ELSE
      zsign = - 1._wp
    END IF
    DO jk = 1, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zaheeu(ji, jj, jk) = zsign * pahu(ji, jj, jk) * e2_e1u(ji, jj) * e3u_n(ji, jj, jk)
          zaheev(ji, jj, jk) = zsign * pahv(ji, jj, jk) * e1_e2v(ji, jj) * e3v_n(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    DO jn = 1, kjpt
      !$ACC KERNELS
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            ztu(ji, jj, jk) = zaheeu(ji, jj, jk) * (ptb(ji + 1, jj, jk, jn) - ptb(ji, jj, jk, jn))
            ztv(ji, jj, jk) = zaheev(ji, jj, jk) * (ptb(ji, jj + 1, jk, jn) - ptb(ji, jj, jk, jn))
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      IF (ln_zps) THEN
        !$ACC KERNELS
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            ztu(ji, jj, mbku(ji, jj)) = zaheeu(ji, jj, mbku(ji, jj)) * pgu(ji, jj, jn)
            ztv(ji, jj, mbkv(ji, jj)) = zaheev(ji, jj, mbkv(ji, jj)) * pgv(ji, jj, jn)
          END DO
        END DO
        !$ACC END KERNELS
        IF (ln_isfcav) THEN
          !$ACC KERNELS
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              IF (miku(ji, jj) > 1) ztu(ji, jj, miku(ji, jj)) = zaheeu(ji, jj, miku(ji, jj)) * pgui(ji, jj, jn)
              IF (mikv(ji, jj) > 1) ztv(ji, jj, mikv(ji, jj)) = zaheev(ji, jj, mikv(ji, jj)) * pgvi(ji, jj, jn)
            END DO
          END DO
          !$ACC END KERNELS
        END IF
      END IF
      !$ACC KERNELS
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            pta(ji, jj, jk, jn) = pta(ji, jj, jk, jn) + (ztu(ji, jj, jk) - ztu(ji - 1, jj, jk) + ztv(ji, jj, jk) - ztv(ji, jj - 1, &
&jk)) / (e1e2t(ji, jj) * e3t_n(ji, jj, jk))
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('tra_ldf_lap', 'r1', 0, 0)
      IF ((kpass == 1 .AND. .NOT. ln_traldf_blp) .OR. (kpass == 2 .AND. ln_traldf_blp)) THEN
        IF (l_ptr) CALL dia_ptr_hst(jn, 'ldf', - ztv(:, :, :))
        IF (l_hst) CALL dia_ar5_hst(jn, 'ldf', - ztu(:, :, :), - ztv(:, :, :))
      END IF
      CALL profile_psy_data1 % PostEnd
    END DO
  END SUBROUTINE tra_ldf_lap
  SUBROUTINE tra_ldf_blp(kt, kit000, cdtype, pahu, pahv, pgu, pgv, pgui, pgvi, ptb, pta, kjpt, kldf)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kit000
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    INTEGER, INTENT(IN) :: kjpt
    INTEGER, INTENT(IN) :: kldf
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pahu, pahv
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt), INTENT(IN) :: pgu, pgv
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt), INTENT(IN) :: pgui, pgvi
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: ptb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    INTEGER :: ji, jj, jk, jn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt) :: zlap
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt) :: zglu, zglv
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt) :: zgui, zgvi
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('tra_ldf_blp', 'r0', 0, 0)
    IF (kt == kit000 .AND. lwp) THEN
      WRITE(numout, FMT = *)
      SELECT CASE (kldf)
      CASE (np_blp)
        WRITE(numout, FMT = *) 'tra_ldf_blp : iso-level   bilaplacian operator on ', cdtype
      CASE (np_blp_i)
        WRITE(numout, FMT = *) 'tra_ldf_blp : iso-neutral bilaplacian operator on ', cdtype, ' (Standard)'
      CASE (np_blp_it)
        WRITE(numout, FMT = *) 'tra_ldf_blp : iso-neutral bilaplacian operator on ', cdtype, ' (triad)'
      END SELECT
      WRITE(numout, FMT = *) '~~~~~~~~~~~'
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zlap(:, :, :, :) = 0._wp
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('tra_ldf_blp', 'r1', 0, 0)
    SELECT CASE (kldf)
    CASE (np_blp)
      CALL tra_ldf_lap(kt, kit000, cdtype, pahu, pahv, pgu, pgv, pgui, pgvi, ptb, zlap, kjpt, 1)
    CASE (np_blp_i)
      CALL tra_ldf_iso(kt, kit000, cdtype, pahu, pahv, pgu, pgv, pgui, pgvi, ptb, ptb, zlap, kjpt, 1)
    CASE (np_blp_it)
      CALL tra_ldf_triad(kt, kit000, cdtype, pahu, pahv, pgu, pgv, pgui, pgvi, ptb, ptb, zlap, kjpt, 1)
    END SELECT
    CALL lbc_lnk('traldf_lap_blp', zlap(:, :, :, :), 'T', 1.)
    IF (ln_isfcav .AND. ln_zps) THEN
      CALL zps_hde_isf(kt, kjpt, zlap, zglu, zglv, zgui, zgvi)
    ELSE IF (ln_zps) THEN
      CALL zps_hde(kt, kjpt, zlap, zglu, zglv)
    END IF
    SELECT CASE (kldf)
    CASE (np_blp)
      CALL tra_ldf_lap(kt, kit000, cdtype, pahu, pahv, zglu, zglv, zgui, zgvi, zlap, pta, kjpt, 2)
    CASE (np_blp_i)
      CALL tra_ldf_iso(kt, kit000, cdtype, pahu, pahv, zglu, zglv, zgui, zgvi, zlap, ptb, pta, kjpt, 2)
    CASE (np_blp_it)
      CALL tra_ldf_triad(kt, kit000, cdtype, pahu, pahv, zglu, zglv, zgui, zgvi, zlap, ptb, pta, kjpt, 2)
    END SELECT
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE tra_ldf_blp
END MODULE traldf_lap_blp