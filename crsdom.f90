MODULE crsdom
  USE dom_oce
  USE crs
  USE in_out_manager
  USE par_kind
  USE crslbclnk
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: crs_dom_ope
  PUBLIC :: crs_dom_e3, crs_dom_sfc, crs_dom_msk, crs_dom_hgr, crs_dom_coordinates
  PUBLIC :: crs_dom_facvol, crs_dom_def, crs_dom_bat
  INTERFACE crs_dom_ope
    MODULE PROCEDURE crs_dom_ope_3d, crs_dom_ope_2d
  END INTERFACE
  REAL(KIND = wp) :: r_inf = 1E+36
  CONTAINS
  SUBROUTINE crs_dom_msk
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jj, jk
    INTEGER :: ijie, ijis, ijje, ijjs, ij, je_2
    REAL(KIND = wp) :: zmask
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    !$ACC KERNELS
    tmask_crs(:, :, :) = 0.0
    vmask_crs(:, :, :) = 0.0
    umask_crs(:, :, :) = 0.0
    fmask_crs(:, :, :) = 0.0
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('crs_dom_msk', 'r0', 0, 0)
    IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
      IF (mje_crs(2) - mjs_crs(2) == 1) THEN
        je_2 = mje_crs(2)
        ij = je_2
      END IF
    ELSE
      je_2 = mje_crs(2)
      ij = mjs_crs(2)
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    DO jk = 1, jpkm1
      DO ji = 2, nlei_crs
        ijis = mis_crs(ji)
        ijie = mie_crs(ji)
        zmask = 0.0
        zmask = SUM(tmask(ijis : ijie, ij : je_2, jk))
        IF (zmask > 0.0) tmask_crs(ji, 2, jk) = 1.0
        zmask = 0.0
        zmask = SUM(vmask(ijis : ijie, je_2, jk))
        IF (zmask > 0.0) vmask_crs(ji, 2, jk) = 1.0
        zmask = 0.0
        zmask = SUM(umask(ijie, ij : je_2, jk))
        IF (zmask > 0.0) umask_crs(ji, 2, jk) = 1.0
        fmask_crs(ji, je_2, jk) = fmask(ijie, 2, jk)
      END DO
    END DO
    DO jk = 1, jpkm1
      DO ji = 2, nlei_crs
        ijis = mis_crs(ji)
        ijie = mie_crs(ji)
        DO jj = 3, nlej_crs
          ijjs = mjs_crs(jj)
          ijje = mje_crs(jj)
          zmask = 0.0
          zmask = SUM(tmask(ijis : ijie, ijjs : ijje, jk))
          IF (zmask > 0.0) tmask_crs(ji, jj, jk) = 1.0
          zmask = 0.0
          zmask = SUM(vmask(ijis : ijie, ijje, jk))
          IF (zmask > 0.0) vmask_crs(ji, jj, jk) = 1.0
          zmask = 0.0
          zmask = SUM(umask(ijie, ijjs : ijje, jk))
          IF (zmask > 0.0) umask_crs(ji, jj, jk) = 1.0
          fmask_crs(ji, jj, jk) = fmask(ijie, ijje, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('crs_dom_msk', 'r1', 0, 0)
    CALL crs_lbc_lnk(tmask_crs, 'T', 1.0)
    CALL crs_lbc_lnk(vmask_crs, 'V', 1.0)
    CALL crs_lbc_lnk(umask_crs, 'U', 1.0)
    CALL crs_lbc_lnk(fmask_crs, 'F', 1.0)
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE crs_dom_msk
  SUBROUTINE crs_dom_coordinates(p_gphi, p_glam, cd_type, p_gphi_crs, p_glam_crs)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: p_gphi
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: p_glam
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs), INTENT(OUT) :: p_gphi_crs
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs), INTENT(OUT) :: p_glam_crs
    INTEGER :: ji, jj, jk
    INTEGER :: ijis, ijjs
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    SELECT CASE (cd_type)
    CASE ('T')
      DO jj = nldj_crs, nlej_crs
        ijjs = mjs_crs(jj) + mybinctr
        DO ji = 2, nlei_crs
          ijis = mis_crs(ji) + mxbinctr
          p_gphi_crs(ji, jj) = p_gphi(ijis, ijjs)
          p_glam_crs(ji, jj) = p_glam(ijis, ijjs)
        END DO
      END DO
    CASE ('U')
      DO jj = nldj_crs, nlej_crs
        ijjs = mjs_crs(jj) + mybinctr
        DO ji = 2, nlei_crs
          ijis = mis_crs(ji)
          p_gphi_crs(ji, jj) = p_gphi(ijis, ijjs)
          p_glam_crs(ji, jj) = p_glam(ijis, ijjs)
        END DO
      END DO
    CASE ('V')
      DO jj = nldj_crs, nlej_crs
        ijjs = mjs_crs(jj)
        DO ji = 2, nlei_crs
          ijis = mis_crs(ji) + mxbinctr
          p_gphi_crs(ji, jj) = p_gphi(ijis, ijjs)
          p_glam_crs(ji, jj) = p_glam(ijis, ijjs)
        END DO
      END DO
    CASE ('F')
      DO jj = nldj_crs, nlej_crs
        ijjs = mjs_crs(jj)
        DO ji = 2, nlei_crs
          ijis = mis_crs(ji)
          p_gphi_crs(ji, jj) = p_gphi(ijis, ijjs)
          p_glam_crs(ji, jj) = p_glam(ijis, ijjs)
        END DO
      END DO
    END SELECT
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('crs_dom_coordinates', 'r0', 0, 0)
    CALL crs_lbc_lnk(p_gphi_crs, cd_type, 1.0)
    CALL crs_lbc_lnk(p_glam_crs, cd_type, 1.0)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    SELECT CASE (cd_type)
    CASE ('T', 'V')
      DO ji = 2, nlei_crs
        ijis = mis_crs(ji) + mxbinctr
        p_gphi_crs(ji, 1) = p_gphi(ijis, 1)
        p_glam_crs(ji, 1) = p_glam(ijis, 1)
      END DO
    CASE ('U', 'F')
      DO ji = 2, nlei_crs
        ijis = mis_crs(ji)
        p_gphi_crs(ji, 1) = p_gphi(ijis, 1)
        p_glam_crs(ji, 1) = p_glam(ijis, 1)
      END DO
    END SELECT
    !$ACC END KERNELS
  END SUBROUTINE crs_dom_coordinates
  SUBROUTINE crs_dom_hgr(p_e1, p_e2, cd_type, p_e1_crs, p_e2_crs)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: p_e1
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: p_e2
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs), INTENT(OUT) :: p_e1_crs
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs), INTENT(OUT) :: p_e2_crs
    INTEGER :: ji, jj, jk
    INTEGER :: ijie, ijje, ijrs
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    DO jk = 1, jpk
      DO ji = 2, nlei_crs
        ijie = mie_crs(ji)
        DO jj = nldj_crs, nlej_crs
          ijje = mje_crs(jj)
          ijrs = mje_crs(jj) - mjs_crs(jj)
          SELECT CASE (cd_type)
          CASE ('T')
            IF (ijrs == 0 .OR. ijrs == 1) THEN
              p_e1_crs(ji, jj) = p_e1(ijie - 1, ijje) * nn_factx
              p_e2_crs(ji, jj) = p_e2(ijie - 1, ijje) * nn_facty
            ELSE
              p_e1_crs(ji, jj) = p_e1(ijie - 1, ijje - 1) * nn_factx
              p_e2_crs(ji, jj) = p_e2(ijie - 1, ijje - 1) * nn_facty
            END IF
          CASE ('U')
            IF (ijrs == 0 .OR. ijrs == 1) THEN
              p_e1_crs(ji, jj) = p_e1(ijie, ijje) * nn_factx
              p_e2_crs(ji, jj) = p_e2(ijie, ijje) * nn_facty
            ELSE
              p_e1_crs(ji, jj) = p_e1(ijie, ijje - 1) * nn_factx
              p_e2_crs(ji, jj) = p_e2(ijie, ijje - 1) * nn_facty
            END IF
          CASE ('V')
            p_e1_crs(ji, jj) = p_e1(ijie - 1, ijje) * nn_factx
            p_e2_crs(ji, jj) = p_e2(ijie - 1, ijje) * nn_facty
          CASE ('F')
            p_e1_crs(ji, jj) = p_e1(ijie, ijje) * nn_factx
            p_e2_crs(ji, jj) = p_e2(ijie, ijje) * nn_facty
          END SELECT
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('crs_dom_hgr', 'r0', 0, 0)
    CALL crs_lbc_lnk(p_e1_crs, cd_type, 1.0, pval = 1.0)
    CALL crs_lbc_lnk(p_e2_crs, cd_type, 1.0, pval = 1.0)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE crs_dom_hgr
  SUBROUTINE crs_dom_facvol(p_mask, cd_type, p_e1, p_e2, p_e3, p_fld1_crs, p_fld2_crs)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: p_mask
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: p_e1
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: p_e2
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: p_e3
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs, jpk), INTENT(OUT) :: p_fld1_crs
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs, jpk), INTENT(OUT) :: p_fld2_crs
    INTEGER :: ji, jj, jk, ii, ij, je_2
    REAL(KIND = wp) :: zdAm
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zvol, zmask
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    p_fld1_crs(:, :, :) = 0._wp
    p_fld2_crs(:, :, :) = 0._wp
    DO jk = 1, jpk
      zvol(:, :, jk) = p_e1(:, :) * p_e2(:, :) * p_e3(:, :, jk)
    END DO
    zmask(:, :, :) = 0._wp
    IF (cd_type == 'W') THEN
      zmask(:, :, 1) = p_mask(:, :, 1)
      DO jk = 2, jpk
        zmask(:, :, jk) = p_mask(:, :, jk - 1)
      END DO
    ELSE
      DO jk = 1, jpk
        zmask(:, :, jk) = p_mask(:, :, jk)
      END DO
    END IF
    !$ACC END KERNELS
    IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
      IF (mje_crs(2) - mjs_crs(2) == 1) THEN
        !$ACC KERNELS
        je_2 = mje_crs(2)
        DO jk = 1, jpk
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            p_fld1_crs(ii, 2, jk) = zvol(ji, je_2, jk) + zvol(ji + 1, je_2, jk) + zvol(ji + 2, je_2, jk) + zvol(ji, je_2 - 1, jk) + zvol(ji + 1, je_2 - 1, jk) + zvol(ji + 2, je_2 - 1, jk)
            zdAm = zvol(ji, je_2, jk) * zmask(ji, je_2, jk) + zvol(ji + 1, je_2, jk) * zmask(ji + 1, je_2, jk) + zvol(ji + 2, je_2, jk) * zmask(ji + 2, je_2, jk)
            p_fld2_crs(ii, 2, jk) = zdAm / p_fld1_crs(ii, 2, jk)
          END DO
        END DO
        !$ACC END KERNELS
      END IF
    ELSE
      !$ACC KERNELS
      je_2 = mjs_crs(2)
      DO jk = 1, jpk
        DO ji = nistr, niend, nn_factx
          ii = (ji - mis_crs(2)) * rfactx_r + 2
          p_fld1_crs(ii, 2, jk) = zvol(ji, je_2, jk) + zvol(ji + 1, je_2, jk) + zvol(ji + 2, je_2, jk) + zvol(ji, je_2 + 1, jk) + zvol(ji + 1, je_2 + 1, jk) + zvol(ji + 2, je_2 + 1, jk) + zvol(ji, je_2 + 2, jk) + zvol(ji + 1, je_2 + 2, jk) + zvol(ji + 2, je_2 + 2, jk)
          zdAm = zvol(ji, je_2, jk) * zmask(ji, je_2, jk) + zvol(ji + 1, je_2, jk) * zmask(ji + 1, je_2, jk) + zvol(ji + 2, je_2, jk) * zmask(ji + 2, je_2, jk) + zvol(ji, je_2 + 1, jk) * zmask(ji, je_2 + 1, jk) + zvol(ji + 1, je_2 + 1, jk) * zmask(ji + 1, je_2 + 1, jk) + zvol(ji + 2, je_2 + 1, jk) * zmask(ji + 2, je_2 + 1, jk) + zvol(ji, je_2 + 2, jk) * zmask(ji, je_2 + 2, jk) + zvol(ji + 1, je_2 + 2, jk) * zmask(ji + 1, je_2 + 2, jk) + zvol(ji + 2, je_2 + 2, jk) * zmask(ji + 2, je_2 + 2, jk)
          p_fld2_crs(ii, 2, jk) = zdAm / p_fld1_crs(ii, 2, jk)
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    DO jk = 1, jpk
      DO jj = njstr, njend, nn_facty
        DO ji = nistr, niend, nn_factx
          ii = (ji - mis_crs(2)) * rfactx_r + 2
          ij = (jj - njstr) * rfacty_r + 3
          p_fld1_crs(ii, ij, jk) = zvol(ji, jj, jk) + zvol(ji + 1, jj, jk) + zvol(ji + 2, jj, jk) + zvol(ji, jj + 1, jk) + zvol(ji + 1, jj + 1, jk) + zvol(ji + 2, jj + 1, jk) + zvol(ji, jj + 2, jk) + zvol(ji + 1, jj + 2, jk) + zvol(ji + 2, jj + 2, jk)
          zdAm = zvol(ji, jj, jk) * zmask(ji, jj, jk) + zvol(ji + 1, jj, jk) * zmask(ji + 1, jj, jk) + zvol(ji + 2, jj, jk) * zmask(ji + 2, jj, jk) + zvol(ji, jj + 1, jk) * zmask(ji, jj + 1, jk) + zvol(ji + 1, jj + 1, jk) * zmask(ji + 1, jj + 1, jk) + zvol(ji + 2, jj + 1, jk) * zmask(ji + 2, jj + 1, jk) + zvol(ji, jj + 2, jk) * zmask(ji, jj + 2, jk) + zvol(ji + 1, jj + 2, jk) * zmask(ji + 1, jj + 2, jk) + zvol(ji + 2, jj + 2, jk) * zmask(ji + 2, jj + 2, jk)
          p_fld2_crs(ii, ij, jk) = zdAm / p_fld1_crs(ii, ij, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('crs_dom_facvol', 'r0', 0, 0)
    CALL crs_lbc_lnk(p_fld1_crs, cd_type, 1.0)
    CALL crs_lbc_lnk(p_fld2_crs, cd_type, 1.0)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE crs_dom_facvol
  SUBROUTINE crs_dom_ope_3d(p_fld, cd_op, cd_type, p_mask, p_fld_crs, p_e12, p_e3, p_surf_crs, p_mask_crs, psgn)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: p_fld
    CHARACTER(LEN = 3), INTENT(IN) :: cd_op
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: p_mask
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN), OPTIONAL :: p_e12
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN), OPTIONAL :: p_e3
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs, jpk), INTENT(IN), OPTIONAL :: p_surf_crs
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs, jpk), INTENT(IN), OPTIONAL :: p_mask_crs
    REAL(KIND = wp), INTENT(IN) :: psgn
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs, jpk), INTENT(OUT) :: p_fld_crs
    INTEGER :: ji, jj, jk
    INTEGER :: ii, ij, ijie, ijje, je_2
    REAL(KIND = wp) :: zflcrs, zsfcrs
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: zsurf, zsurfmsk, zmask
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    !$ACC KERNELS
    p_fld_crs(:, :, :) = 0._wp
    !$ACC END KERNELS
    SELECT CASE (cd_op)
    CASE ('VOL')
      ALLOCATE(zsurf(jpi, jpj, jpk), zsurfmsk(jpi, jpj, jpk))
      SELECT CASE (cd_type)
      CASE ('T', 'W')
        !$ACC KERNELS
        IF (cd_type == 'T') THEN
          DO jk = 1, jpk
            zsurf(:, :, jk) = p_e12(:, :) * p_e3(:, :, jk) * p_mask(:, :, jk)
            zsurfmsk(:, :, jk) = zsurf(:, :, jk)
          END DO
        ELSE
          zsurf(:, :, 1) = p_e12(:, :) * p_e3(:, :, 1)
          zsurfmsk(:, :, 1) = zsurf(:, :, 1) * p_mask(:, :, 1)
          DO jk = 2, jpk
            zsurf(:, :, jk) = p_e12(:, :) * p_e3(:, :, jk)
            zsurfmsk(:, :, jk) = zsurf(:, :, jk) * p_mask(:, :, jk - 1)
          END DO
        END IF
        !$ACC END KERNELS
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            !$ACC KERNELS
            je_2 = mje_crs(2)
            DO jk = 1, jpk
              DO ji = nistr, niend, nn_factx
                ii = (ji - mis_crs(2)) * rfactx_r + 2
                zflcrs = p_fld(ji, je_2, jk) * zsurfmsk(ji, je_2, jk) + p_fld(ji + 1, je_2, jk) * zsurfmsk(ji + 1, je_2, jk) + p_fld(ji + 2, je_2, jk) * zsurfmsk(ji + 2, je_2, jk)
                zsfcrs = zsurf(ji, je_2, jk) + zsurf(ji + 1, je_2, jk) + zsurf(ji + 2, je_2, jk)
                p_fld_crs(ii, 2, jk) = zflcrs
                IF (zsfcrs /= 0.0) p_fld_crs(ii, 2, jk) = zflcrs / zsfcrs
              END DO
            END DO
            !$ACC END KERNELS
          END IF
        ELSE
          !$ACC KERNELS
          je_2 = mjs_crs(2)
          DO jk = 1, jpk
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              zflcrs = p_fld(ji, je_2, jk) * zsurfmsk(ji, je_2, jk) + p_fld(ji + 1, je_2, jk) * zsurfmsk(ji + 1, je_2, jk) + p_fld(ji + 2, je_2, jk) * zsurfmsk(ji + 2, je_2, jk) + p_fld(ji, je_2 + 1, jk) * zsurfmsk(ji, je_2 + 1, jk) + p_fld(ji + 1, je_2 + 1, jk) * zsurfmsk(ji + 1, je_2 + 1, jk) + p_fld(ji + 2, je_2 + 1, jk) * zsurfmsk(ji + 2, je_2 + 1, jk) + p_fld(ji, je_2 + 2, jk) * zsurfmsk(ji, je_2 + 2, jk) + p_fld(ji + 1, je_2 + 2, jk) * zsurfmsk(ji + 1, je_2 + 2, jk) + p_fld(ji + 2, je_2 + 2, jk) * zsurfmsk(ji + 2, je_2 + 2, jk)
              zsfcrs = zsurf(ji, je_2, jk) + zsurf(ji + 1, je_2, jk) + zsurf(ji + 2, je_2, jk) + zsurf(ji, je_2 + 1, jk) + zsurf(ji + 1, je_2 + 1, jk) + zsurf(ji + 2, je_2 + 1, jk) + zsurf(ji, je_2 + 2, jk) + zsurf(ji + 1, je_2 + 2, jk) + zsurf(ji + 2, je_2 + 2, jk)
              p_fld_crs(ii, 2, jk) = zflcrs
              IF (zsfcrs /= 0.0) p_fld_crs(ii, 2, jk) = zflcrs / zsfcrs
            END DO
          END DO
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        DO jk = 1, jpk
          DO jj = njstr, njend, nn_facty
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              ij = (jj - njstr) * rfacty_r + 3
              zflcrs = p_fld(ji, jj, jk) * zsurfmsk(ji, jj, jk) + p_fld(ji + 1, jj, jk) * zsurfmsk(ji + 1, jj, jk) + p_fld(ji + 2, jj, jk) * zsurfmsk(ji + 2, jj, jk) + p_fld(ji, jj + 1, jk) * zsurfmsk(ji, jj + 1, jk) + p_fld(ji + 1, jj + 1, jk) * zsurfmsk(ji + 1, jj + 1, jk) + p_fld(ji + 2, jj + 1, jk) * zsurfmsk(ji + 2, jj + 1, jk) + p_fld(ji, jj + 2, jk) * zsurfmsk(ji, jj + 2, jk) + p_fld(ji + 1, jj + 2, jk) * zsurfmsk(ji + 1, jj + 2, jk) + p_fld(ji + 2, jj + 2, jk) * zsurfmsk(ji + 2, jj + 2, jk)
              zsfcrs = zsurf(ji, jj, jk) + zsurf(ji + 1, jj, jk) + zsurf(ji + 2, jj, jk) + zsurf(ji, jj + 1, jk) + zsurf(ji + 1, jj + 1, jk) + zsurf(ji + 2, jj + 1, jk) + zsurf(ji, jj + 2, jk) + zsurf(ji + 1, jj + 2, jk) + zsurf(ji + 2, jj + 2, jk)
              p_fld_crs(ii, ij, jk) = zflcrs
              IF (zsfcrs /= 0.0) p_fld_crs(ii, ij, jk) = zflcrs / zsfcrs
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      CASE DEFAULT
        STOP
      END SELECT
      DEALLOCATE(zsurf, zsurfmsk)
    CASE ('SUM')
      ALLOCATE(zsurfmsk(jpi, jpj, jpk))
      SELECT CASE (cd_type)
      CASE ('W')
        IF (PRESENT(p_e3)) THEN
          !$ACC KERNELS
          zsurfmsk(:, :, 1) = p_e12(:, :) * p_e3(:, :, 1) * p_mask(:, :, 1)
          DO jk = 2, jpk
            zsurfmsk(:, :, jk) = p_e12(:, :) * p_e3(:, :, jk) * p_mask(:, :, jk - 1)
          END DO
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          zsurfmsk(:, :, 1) = p_e12(:, :) * p_mask(:, :, 1)
          DO jk = 2, jpk
            zsurfmsk(:, :, jk) = p_e12(:, :) * p_mask(:, :, jk - 1)
          END DO
          !$ACC END KERNELS
        END IF
      CASE DEFAULT
        IF (PRESENT(p_e3)) THEN
          !$ACC KERNELS
          DO jk = 1, jpk
            zsurfmsk(:, :, jk) = p_e12(:, :) * p_e3(:, :, jk) * p_mask(:, :, jk)
          END DO
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          DO jk = 1, jpk
            zsurfmsk(:, :, jk) = p_e12(:, :) * p_mask(:, :, jk)
          END DO
          !$ACC END KERNELS
        END IF
      END SELECT
      SELECT CASE (cd_type)
      CASE ('T', 'W')
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            !$ACC KERNELS
            je_2 = mje_crs(2)
            DO jk = 1, jpk
              DO ji = nistr, niend, nn_factx
                ii = (ji - mis_crs(2)) * rfactx_r + 2
                zflcrs = p_fld(ji, je_2, jk) * zsurfmsk(ji, je_2, jk) + p_fld(ji + 1, je_2, jk) * zsurfmsk(ji + 1, je_2, jk) + p_fld(ji + 2, je_2, jk) * zsurfmsk(ji + 2, je_2, jk)
                p_fld_crs(ii, 2, jk) = zflcrs
              END DO
            END DO
            !$ACC END KERNELS
          END IF
        ELSE
          !$ACC KERNELS
          je_2 = mjs_crs(2)
          DO jk = 1, jpk
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              zflcrs = p_fld(ji, je_2, jk) * zsurfmsk(ji, je_2, jk) + p_fld(ji + 1, je_2, jk) * zsurfmsk(ji + 1, je_2, jk) + p_fld(ji + 2, je_2, jk) * zsurfmsk(ji + 2, je_2, jk) + p_fld(ji, je_2 + 1, jk) * zsurfmsk(ji, je_2 + 1, jk) + p_fld(ji + 1, je_2 + 1, jk) * zsurfmsk(ji + 1, je_2 + 1, jk) + p_fld(ji + 2, je_2 + 1, jk) * zsurfmsk(ji + 2, je_2 + 1, jk) + p_fld(ji, je_2 + 2, jk) * zsurfmsk(ji, je_2 + 2, jk) + p_fld(ji + 1, je_2 + 2, jk) * zsurfmsk(ji + 1, je_2 + 2, jk) + p_fld(ji + 2, je_2 + 2, jk) * zsurfmsk(ji + 2, je_2 + 2, jk)
              p_fld_crs(ii, 2, jk) = zflcrs
            END DO
          END DO
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        DO jk = 1, jpk
          DO jj = njstr, njend, nn_facty
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              ij = (jj - njstr) * rfacty_r + 3
              zflcrs = p_fld(ji, jj, jk) * zsurfmsk(ji, jj, jk) + p_fld(ji + 1, jj, jk) * zsurfmsk(ji + 1, jj, jk) + p_fld(ji + 2, jj, jk) * zsurfmsk(ji + 2, jj, jk) + p_fld(ji, jj + 1, jk) * zsurfmsk(ji, jj + 1, jk) + p_fld(ji + 1, jj + 1, jk) * zsurfmsk(ji + 1, jj + 1, jk) + p_fld(ji + 2, jj + 1, jk) * zsurfmsk(ji + 2, jj + 1, jk) + p_fld(ji, jj + 2, jk) * zsurfmsk(ji, jj + 2, jk) + p_fld(ji + 1, jj + 2, jk) * zsurfmsk(ji + 1, jj + 2, jk) + p_fld(ji + 2, jj + 2, jk) * zsurfmsk(ji + 2, jj + 2, jk)
              p_fld_crs(ii, ij, jk) = zflcrs
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      CASE ('V')
        CALL profile_psy_data0 % PreStart('crs_dom_ope_3d', 'r0', 0, 0)
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            ijje = mje_crs(2)
          END IF
        ELSE
          ijje = mjs_crs(2)
        END IF
        CALL profile_psy_data0 % PostEnd
        !$ACC KERNELS
        DO jk = 1, jpk
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            zflcrs = p_fld(ji, ijje, jk) * zsurfmsk(ji, ijje, jk) + p_fld(ji + 1, ijje, jk) * zsurfmsk(ji + 1, ijje, jk) + p_fld(ji + 2, ijje, jk) * zsurfmsk(ji + 2, ijje, jk)
            p_fld_crs(ii, 2, jk) = zflcrs
          END DO
        END DO
        DO jk = 1, jpk
          DO jj = njstr, njend, nn_facty
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              ij = (jj - njstr) * rfacty_r + 3
              ijje = mje_crs(ij)
              zflcrs = p_fld(ji, ijje, jk) * zsurfmsk(ji, ijje, jk) + p_fld(ji + 1, ijje, jk) * zsurfmsk(ji + 1, ijje, jk) + p_fld(ji + 2, ijje, jk) * zsurfmsk(ji + 2, ijje, jk)
              p_fld_crs(ii, ij, jk) = zflcrs
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      CASE ('U')
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            !$ACC KERNELS
            je_2 = mje_crs(2)
            DO jk = 1, jpk
              DO ji = nistr, niend, nn_factx
                ii = (ji - mis_crs(2)) * rfactx_r + 2
                ijie = mie_crs(ii)
                zflcrs = p_fld(ijie, je_2, jk) * zsurfmsk(ijie, je_2, jk)
                p_fld_crs(ii, 2, jk) = zflcrs
              END DO
            END DO
            !$ACC END KERNELS
          END IF
        ELSE
          !$ACC KERNELS
          je_2 = mjs_crs(2)
          DO jk = 1, jpk
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              ijie = mie_crs(ii)
              zflcrs = p_fld(ijie, je_2, jk) * zsurfmsk(ijie, je_2, jk) + p_fld(ijie, je_2 + 1, jk) * zsurfmsk(ijie, je_2 + 1, jk) + p_fld(ijie, je_2 + 2, jk) * zsurfmsk(ijie, je_2 + 2, jk)
              p_fld_crs(ii, 2, jk) = zflcrs
            END DO
          END DO
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        DO jk = 1, jpk
          DO jj = njstr, njend, nn_facty
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              ij = (jj - njstr) * rfacty_r + 3
              ijie = mie_crs(ii)
              zflcrs = p_fld(ijie, jj, jk) * zsurfmsk(ijie, jj, jk) + p_fld(ijie, jj + 1, jk) * zsurfmsk(ijie, jj + 1, jk) + p_fld(ijie, jj + 2, jk) * zsurfmsk(ijie, jj + 2, jk)
              p_fld_crs(ii, ij, jk) = zflcrs
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      END SELECT
      CALL profile_psy_data1 % PreStart('crs_dom_ope_3d', 'r1', 0, 0)
      IF (PRESENT(p_surf_crs)) THEN
        WHERE (p_surf_crs /= 0.0) p_fld_crs(:, :, :) = p_fld_crs(:, :, :) / p_surf_crs(:, :, :)
      END IF
      DEALLOCATE(zsurfmsk)
      CALL profile_psy_data1 % PostEnd
    CASE ('MAX')
      ALLOCATE(zmask(jpi, jpj, jpk))
      SELECT CASE (cd_type)
      CASE ('W')
        !$ACC KERNELS
        zmask(:, :, 1) = p_mask(:, :, 1)
        DO jk = 2, jpk
          zmask(:, :, jk) = p_mask(:, :, jk - 1)
        END DO
        !$ACC END KERNELS
      CASE ('T')
        !$ACC KERNELS
        DO jk = 1, jpk
          zmask(:, :, jk) = p_mask(:, :, jk)
        END DO
        !$ACC END KERNELS
      END SELECT
      SELECT CASE (cd_type)
      CASE ('T', 'W')
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            !$ACC KERNELS
            je_2 = mje_crs(2)
            DO jk = 1, jpk
              DO ji = nistr, niend, nn_factx
                ii = (ji - mis_crs(2)) * rfactx_r + 2
                zflcrs = MAX(p_fld(ji, je_2, jk) * zmask(ji, je_2, jk) - (1. - zmask(ji, je_2, jk)) * r_inf, p_fld(ji + 1, je_2, jk) * zmask(ji + 1, je_2, jk) - (1. - zmask(ji + 1, je_2, jk)) * r_inf, p_fld(ji + 2, je_2, jk) * zmask(ji + 2, je_2, jk) - (1. - zmask(ji + 2, je_2, jk)) * r_inf)
                p_fld_crs(ii, 2, jk) = zflcrs
              END DO
            END DO
            !$ACC END KERNELS
          END IF
        ELSE
          !$ACC KERNELS
          je_2 = mjs_crs(2)
          DO jk = 1, jpk
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              zflcrs = MAX(p_fld(ji, je_2, jk) * zmask(ji, je_2, jk) - (1. - zmask(ji, je_2, jk)) * r_inf, p_fld(ji + 1, je_2, jk) * zmask(ji + 1, je_2, jk) - (1. - zmask(ji + 1, je_2, jk)) * r_inf, p_fld(ji + 2, je_2, jk) * zmask(ji + 2, je_2, jk) - (1. - zmask(ji + 2, je_2, jk)) * r_inf, p_fld(ji, je_2 + 1, jk) * zmask(ji, je_2 + 1, jk) - (1. - zmask(ji, je_2 + 1, jk)) * r_inf, p_fld(ji + 1, je_2 + 1, jk) * zmask(ji + 1, je_2 + 1, jk) - (1. - zmask(ji + 1, je_2 + 1, jk)) * r_inf, p_fld(ji + 2, je_2 + 1, jk) * zmask(ji + 2, je_2 + 1, jk) - (1. - zmask(ji + 2, je_2 + 1, jk)) * r_inf, p_fld(ji, je_2 + 2, jk) * zmask(ji, je_2 + 2, jk) - (1. - zmask(ji, je_2 + 2, jk)) * r_inf, p_fld(ji + 1, je_2 + 2, jk) * zmask(ji + 1, je_2 + 2, jk) - (1. - zmask(ji + 1, je_2 + 2, jk)) * r_inf, p_fld(ji + 2, je_2 + 2, jk) * zmask(ji + 2, je_2 + 2, jk) - (1. - zmask(ji + 2, je_2 + 2, jk)) * r_inf)
              p_fld_crs(ii, 2, jk) = zflcrs
            END DO
          END DO
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        DO jk = 1, jpk
          DO jj = njstr, njend, nn_facty
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              ij = (jj - njstr) * rfacty_r + 3
              zflcrs = MAX(p_fld(ji, jj, jk) * zmask(ji, jj, jk) - (1. - zmask(ji, jj, jk)) * r_inf, p_fld(ji + 1, jj, jk) * zmask(ji + 1, jj, jk) - (1. - zmask(ji + 1, jj, jk)) * r_inf, p_fld(ji + 2, jj, jk) * zmask(ji + 2, jj, jk) - (1. - zmask(ji + 2, jj, jk)) * r_inf, p_fld(ji, jj + 1, jk) * zmask(ji, jj + 1, jk) - (1. - zmask(ji, jj + 1, jk)) * r_inf, p_fld(ji + 1, jj + 1, jk) * zmask(ji + 1, jj + 1, jk) - (1. - zmask(ji + 1, jj + 1, jk)) * r_inf, p_fld(ji + 2, jj + 1, jk) * zmask(ji + 2, jj + 1, jk) - (1. - zmask(ji + 2, jj + 1, jk)) * r_inf, p_fld(ji, jj + 2, jk) * zmask(ji, jj + 2, jk) - (1. - zmask(ji, jj + 2, jk)) * r_inf, p_fld(ji + 1, jj + 2, jk) * zmask(ji + 1, jj + 2, jk) - (1. - zmask(ji + 1, jj + 2, jk)) * r_inf, p_fld(ji + 2, jj + 2, jk) * zmask(ji + 2, jj + 2, jk) - (1. - zmask(ji + 2, jj + 2, jk)) * r_inf)
              p_fld_crs(ii, ij, jk) = zflcrs
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      CASE ('V')
        CALL profile_psy_data2 % PreStart('crs_dom_ope_3d', 'r2', 0, 0)
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            ijje = mje_crs(2)
          END IF
        ELSE
          ijje = mjs_crs(2)
        END IF
        CALL profile_psy_data2 % PostEnd
        !$ACC KERNELS
        DO jk = 1, jpk
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            zflcrs = MAX(p_fld(ji, ijje, jk) * p_mask(ji, ijje, jk) - (1. - p_mask(ji, ijje, jk)) * r_inf, p_fld(ji + 1, ijje, jk) * p_mask(ji + 1, ijje, jk) - (1. - p_mask(ji, ijje, jk)) * r_inf, p_fld(ji + 2, ijje, jk) * p_mask(ji + 2, ijje, jk) - (1. - p_mask(ji, ijje, jk)) * r_inf)
            p_fld_crs(ii, 2, jk) = zflcrs
          END DO
        END DO
        DO jk = 1, jpk
          DO jj = njstr, njend, nn_facty
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              ij = (jj - njstr) * rfacty_r + 3
              ijje = mje_crs(ij)
              zflcrs = MAX(p_fld(ji, ijje, jk) * p_mask(ji, ijje, jk) - (1. - p_mask(ji, ijje, jk)) * r_inf, p_fld(ji + 1, ijje, jk) * p_mask(ji + 1, ijje, jk) - (1. - p_mask(ji, ijje, jk)) * r_inf, p_fld(ji + 2, ijje, jk) * p_mask(ji + 2, ijje, jk) - (1. - p_mask(ji, ijje, jk)) * r_inf)
              p_fld_crs(ii, ij, jk) = zflcrs
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      CASE ('U')
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            !$ACC KERNELS
            je_2 = mje_crs(2)
            DO jk = 1, jpk
              DO ji = nistr, niend, nn_factx
                ii = (ji - mis_crs(2)) * rfactx_r + 2
                ijie = mie_crs(ii)
                zflcrs = p_fld(ijie, je_2, jk) * p_mask(ijie, je_2, jk) - (1. - p_mask(ijie, je_2, jk)) * r_inf
                p_fld_crs(ii, 2, jk) = zflcrs
              END DO
            END DO
            !$ACC END KERNELS
          END IF
        ELSE
          !$ACC KERNELS
          je_2 = mjs_crs(2)
          DO jk = 1, jpk
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              ijie = mie_crs(ii)
              zflcrs = MAX(p_fld(ijie, je_2, jk) * p_mask(ijie, je_2, jk) - (1. - p_mask(ijie, je_2, jk)) * r_inf, p_fld(ijie, je_2 + 1, jk) * p_mask(ijie, je_2 + 1, jk) - (1. - p_mask(ijie, je_2, jk)) * r_inf, p_fld(ijie, je_2 + 2, jk) * p_mask(ijie, je_2 + 2, jk) - (1. - p_mask(ijie, je_2, jk)) * r_inf)
              p_fld_crs(ii, 2, jk) = zflcrs
            END DO
          END DO
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        DO jk = 1, jpk
          DO jj = njstr, njend, nn_facty
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              ij = (jj - njstr) * rfacty_r + 3
              ijie = mie_crs(ii)
              zflcrs = MAX(p_fld(ijie, jj, jk) * p_mask(ijie, jj, jk) - (1. - p_mask(ijie, jj, jk)) * r_inf, p_fld(ijie, jj + 1, jk) * p_mask(ijie, jj + 1, jk) - (1. - p_mask(ijie, jj, jk)) * r_inf, p_fld(ijie, jj + 2, jk) * p_mask(ijie, jj + 2, jk) - (1. - p_mask(ijie, jj, jk)) * r_inf)
              p_fld_crs(ii, ij, jk) = zflcrs
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      END SELECT
      DEALLOCATE(zmask)
    CASE ('MIN')
      ALLOCATE(zmask(jpi, jpj, jpk))
      SELECT CASE (cd_type)
      CASE ('W')
        !$ACC KERNELS
        zmask(:, :, 1) = p_mask(:, :, 1)
        DO jk = 2, jpk
          zmask(:, :, jk) = p_mask(:, :, jk - 1)
        END DO
        !$ACC END KERNELS
      CASE ('T')
        !$ACC KERNELS
        DO jk = 1, jpk
          zmask(:, :, jk) = p_mask(:, :, jk)
        END DO
        !$ACC END KERNELS
      END SELECT
      SELECT CASE (cd_type)
      CASE ('T', 'W')
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            !$ACC KERNELS
            je_2 = mje_crs(2)
            DO jk = 1, jpk
              DO ji = nistr, niend, nn_factx
                ii = (ji - mis_crs(2)) * rfactx_r + 2
                zflcrs = MIN(p_fld(ji, je_2, jk) * zmask(ji, je_2, jk) + (1. - zmask(ji, je_2, jk)) * r_inf, p_fld(ji + 1, je_2, jk) * zmask(ji + 1, je_2, jk) + (1. - zmask(ji + 1, je_2, jk)) * r_inf, p_fld(ji + 2, je_2, jk) * zmask(ji + 2, je_2, jk) + (1. - zmask(ji + 2, je_2, jk)) * r_inf)
                p_fld_crs(ii, 2, jk) = zflcrs
              END DO
            END DO
            !$ACC END KERNELS
          END IF
        ELSE
          !$ACC KERNELS
          je_2 = mjs_crs(2)
          DO jk = 1, jpk
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              zflcrs = MIN(p_fld(ji, je_2, jk) * zmask(ji, je_2, jk) + (1. - zmask(ji, je_2, jk)) * r_inf, p_fld(ji + 1, je_2, jk) * zmask(ji + 1, je_2, jk) + (1. - zmask(ji + 1, je_2, jk)) * r_inf, p_fld(ji + 2, je_2, jk) * zmask(ji + 2, je_2, jk) + (1. - zmask(ji + 2, je_2, jk)) * r_inf, p_fld(ji, je_2 + 1, jk) * zmask(ji, je_2 + 1, jk) + (1. - zmask(ji, je_2 + 1, jk)) * r_inf, p_fld(ji + 1, je_2 + 1, jk) * zmask(ji + 1, je_2 + 1, jk) + (1. - zmask(ji + 1, je_2 + 1, jk)) * r_inf, p_fld(ji + 2, je_2 + 1, jk) * zmask(ji + 2, je_2 + 1, jk) + (1. - zmask(ji + 2, je_2 + 1, jk)) * r_inf, p_fld(ji, je_2 + 2, jk) * zmask(ji, je_2 + 2, jk) + (1. - zmask(ji, je_2 + 2, jk)) * r_inf, p_fld(ji + 1, je_2 + 2, jk) * zmask(ji + 1, je_2 + 2, jk) + (1. - zmask(ji + 1, je_2 + 2, jk)) * r_inf, p_fld(ji + 2, je_2 + 2, jk) * zmask(ji + 2, je_2 + 2, jk) + (1. - zmask(ji + 2, je_2 + 2, jk)) * r_inf)
              p_fld_crs(ii, 2, jk) = zflcrs
            END DO
          END DO
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        DO jk = 1, jpk
          DO jj = njstr, njend, nn_facty
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              ij = (jj - njstr) * rfacty_r + 3
              zflcrs = MIN(p_fld(ji, jj, jk) * zmask(ji, jj, jk) + (1. - zmask(ji, jj, jk)) * r_inf, p_fld(ji + 1, jj, jk) * zmask(ji + 1, jj, jk) + (1. - zmask(ji + 1, jj, jk)) * r_inf, p_fld(ji + 2, jj, jk) * zmask(ji + 2, jj, jk) + (1. - zmask(ji + 2, jj, jk)) * r_inf, p_fld(ji, jj + 1, jk) * zmask(ji, jj + 1, jk) + (1. - zmask(ji, jj + 1, jk)) * r_inf, p_fld(ji + 1, jj + 1, jk) * zmask(ji + 1, jj + 1, jk) + (1. - zmask(ji + 1, jj + 1, jk)) * r_inf, p_fld(ji + 2, jj + 1, jk) * zmask(ji + 2, jj + 1, jk) + (1. - zmask(ji + 2, jj + 1, jk)) * r_inf, p_fld(ji, jj + 2, jk) * zmask(ji, jj + 2, jk) + (1. - zmask(ji, jj + 2, jk)) * r_inf, p_fld(ji + 1, jj + 2, jk) * zmask(ji + 1, jj + 2, jk) + (1. - zmask(ji + 1, jj + 2, jk)) * r_inf, p_fld(ji + 2, jj + 2, jk) * zmask(ji + 2, jj + 2, jk) + (1. - zmask(ji + 2, jj + 2, jk)) * r_inf)
              p_fld_crs(ii, ij, jk) = zflcrs
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      CASE ('V')
        CALL profile_psy_data3 % PreStart('crs_dom_ope_3d', 'r3', 0, 0)
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            ijje = mje_crs(2)
          END IF
        ELSE
          ijje = mjs_crs(2)
        END IF
        CALL profile_psy_data3 % PostEnd
        !$ACC KERNELS
        DO jk = 1, jpk
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            zflcrs = MIN(p_fld(ji, ijje, jk) * p_mask(ji, ijje, jk) + (1. - p_mask(ji, ijje, jk)) * r_inf, p_fld(ji + 1, ijje, jk) * p_mask(ji + 1, ijje, jk) + (1. - p_mask(ji, ijje, jk)) * r_inf, p_fld(ji + 2, ijje, jk) * p_mask(ji + 2, ijje, jk) + (1. - p_mask(ji, ijje, jk)) * r_inf)
            p_fld_crs(ii, 2, jk) = zflcrs
          END DO
        END DO
        DO jk = 1, jpk
          DO jj = njstr, njend, nn_facty
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              ij = (jj - njstr) * rfacty_r + 3
              ijje = mje_crs(ij)
              zflcrs = MIN(p_fld(ji, ijje, jk) * p_mask(ji, ijje, jk) + (1. - p_mask(ji, ijje, jk)) * r_inf, p_fld(ji + 1, ijje, jk) * p_mask(ji + 1, ijje, jk) + (1. - p_mask(ji, ijje, jk)) * r_inf, p_fld(ji + 2, ijje, jk) * p_mask(ji + 2, ijje, jk) + (1. - p_mask(ji, ijje, jk)) * r_inf)
              p_fld_crs(ii, ij, jk) = zflcrs
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      CASE ('U')
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            !$ACC KERNELS
            je_2 = mje_crs(2)
            DO jk = 1, jpk
              DO ji = nistr, niend, nn_factx
                ii = (ji - mis_crs(2)) * rfactx_r + 2
                ijie = mie_crs(ii)
                zflcrs = p_fld(ijie, je_2, jk) * p_mask(ijie, je_2, jk) + (1. - p_mask(ijie, je_2, jk)) * r_inf
                p_fld_crs(ii, 2, jk) = zflcrs
              END DO
            END DO
            !$ACC END KERNELS
          END IF
        ELSE
          !$ACC KERNELS
          je_2 = mjs_crs(2)
          DO jk = 1, jpk
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              ijie = mie_crs(ii)
              zflcrs = MIN(p_fld(ijie, je_2, jk) * p_mask(ijie, je_2, jk) + (1. - p_mask(ijie, je_2, jk)) * r_inf, p_fld(ijie, je_2 + 1, jk) * p_mask(ijie, je_2 + 1, jk) + (1. - p_mask(ijie, je_2, jk)) * r_inf, p_fld(ijie, je_2 + 2, jk) * p_mask(ijie, je_2 + 2, jk) + (1. - p_mask(ijie, je_2, jk)) * r_inf)
              p_fld_crs(ii, 2, jk) = zflcrs
            END DO
          END DO
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        DO jk = 1, jpk
          DO jj = njstr, njend, nn_facty
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              ij = (jj - njstr) * rfacty_r + 3
              ijie = mie_crs(ii)
              zflcrs = MIN(p_fld(ijie, jj, jk) * p_mask(ijie, jj, jk) + (1. - p_mask(ijie, jj, jk)) * r_inf, p_fld(ijie, jj + 1, jk) * p_mask(ijie, jj + 1, jk) + (1. - p_mask(ijie, jj, jk)) * r_inf, p_fld(ijie, jj + 2, jk) * p_mask(ijie, jj + 2, jk) + (1. - p_mask(ijie, jj, jk)) * r_inf)
              p_fld_crs(ii, ij, jk) = zflcrs
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      END SELECT
      DEALLOCATE(zmask)
    END SELECT
    CALL profile_psy_data4 % PreStart('crs_dom_ope_3d', 'r4', 0, 0)
    CALL crs_lbc_lnk(p_fld_crs, cd_type, psgn)
    CALL profile_psy_data4 % PostEnd
  END SUBROUTINE crs_dom_ope_3d
  SUBROUTINE crs_dom_ope_2d(p_fld, cd_op, cd_type, p_mask, p_fld_crs, p_e12, p_e3, p_surf_crs, p_mask_crs, psgn)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: p_fld
    CHARACTER(LEN = 3), INTENT(IN) :: cd_op
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: p_mask
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN), OPTIONAL :: p_e12
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN), OPTIONAL :: p_e3
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs), INTENT(IN), OPTIONAL :: p_surf_crs
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs, jpk), INTENT(IN), OPTIONAL :: p_mask_crs
    REAL(KIND = wp), INTENT(IN) :: psgn
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs), INTENT(OUT) :: p_fld_crs
    INTEGER :: ji, jj, jk
    INTEGER :: ijie, ijje, ii, ij, je_2
    REAL(KIND = wp) :: zflcrs, zsfcrs
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: zsurfmsk
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    !$ACC KERNELS
    p_fld_crs(:, :) = 0._wp
    !$ACC END KERNELS
    SELECT CASE (cd_op)
    CASE ('VOL')
      ALLOCATE(zsurfmsk(jpi, jpj))
      !$ACC KERNELS
      zsurfmsk(:, :) = p_e12(:, :) * p_e3(:, :, 1) * p_mask(:, :, 1)
      !$ACC END KERNELS
      IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
        IF (mje_crs(2) - mjs_crs(2) == 1) THEN
          !$ACC KERNELS
          je_2 = mje_crs(2)
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            zflcrs = p_fld(ji, je_2) * zsurfmsk(ji, je_2) + p_fld(ji + 1, je_2) * zsurfmsk(ji + 1, je_2) + p_fld(ji + 2, je_2) * zsurfmsk(ji + 2, je_2)
            zsfcrs = zsurfmsk(ji, je_2) + zsurfmsk(ji + 1, je_2) + zsurfmsk(ji + 2, je_2)
            p_fld_crs(ii, 2) = zflcrs
            IF (zsfcrs /= 0.0) p_fld_crs(ii, 2) = zflcrs / zsfcrs
          END DO
          !$ACC END KERNELS
        END IF
      ELSE
        !$ACC KERNELS
        je_2 = mjs_crs(2)
        DO ji = nistr, niend, nn_factx
          ii = (ji - mis_crs(2)) * rfactx_r + 2
          zflcrs = p_fld(ji, je_2) * zsurfmsk(ji, je_2) + p_fld(ji + 1, je_2) * zsurfmsk(ji + 1, je_2) + p_fld(ji + 2, je_2) * zsurfmsk(ji + 2, je_2) + p_fld(ji, je_2 + 1) * zsurfmsk(ji, je_2 + 1) + p_fld(ji + 1, je_2 + 1) * zsurfmsk(ji + 1, je_2 + 1) + p_fld(ji + 2, je_2 + 1) * zsurfmsk(ji + 2, je_2 + 1) + p_fld(ji, je_2 + 2) * zsurfmsk(ji, je_2 + 2) + p_fld(ji + 1, je_2 + 2) * zsurfmsk(ji + 1, je_2 + 2) + p_fld(ji + 2, je_2 + 2) * zsurfmsk(ji + 2, je_2 + 2)
          zsfcrs = zsurfmsk(ji, je_2) + zsurfmsk(ji + 1, je_2) + zsurfmsk(ji + 2, je_2) + zsurfmsk(ji, je_2 + 1) + zsurfmsk(ji + 1, je_2 + 1) + zsurfmsk(ji + 2, je_2 + 1) + zsurfmsk(ji, je_2 + 2) + zsurfmsk(ji + 1, je_2 + 2) + zsurfmsk(ji + 2, je_2 + 2)
          p_fld_crs(ii, 2) = zflcrs
          IF (zsfcrs /= 0.0) p_fld_crs(ii, 2) = zflcrs / zsfcrs
        END DO
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      !$ACC LOOP INDEPENDENT ! CDe compiler bug prevents use of COLLAPSE(2) here
      DO jj = njstr, njend, nn_facty
        !$ACC LOOP INDEPENDENT ! CDe 
        DO ji = nistr, niend, nn_factx
          ii = (ji - mis_crs(2)) * rfactx_r + 2
          ij = (jj - njstr) * rfacty_r + 3
          zflcrs = p_fld(ji, jj) * zsurfmsk(ji, jj) + p_fld(ji + 1, jj) * zsurfmsk(ji + 1, jj) + p_fld(ji + 2, jj) * zsurfmsk(ji + 2, jj) + p_fld(ji, jj + 1) * zsurfmsk(ji, jj + 1) + p_fld(ji + 1, jj + 1) * zsurfmsk(ji + 1, jj + 1) + p_fld(ji + 2, jj + 1) * zsurfmsk(ji + 2, jj + 1) + p_fld(ji, jj + 2) * zsurfmsk(ji, jj + 2) + p_fld(ji + 1, jj + 2) * zsurfmsk(ji + 1, jj + 2) + p_fld(ji + 2, jj + 2) * zsurfmsk(ji + 2, jj + 2)
          zsfcrs = zsurfmsk(ji, jj) + zsurfmsk(ji + 1, jj) + zsurfmsk(ji + 2, jj) + zsurfmsk(ji, jj + 1) + zsurfmsk(ji + 1, jj + 1) + zsurfmsk(ji + 2, jj + 1) + zsurfmsk(ji, jj + 2) + zsurfmsk(ji + 1, jj + 2) + zsurfmsk(ji + 2, jj + 2)
          p_fld_crs(ii, ij) = zflcrs
          IF (zsfcrs /= 0.0) p_fld_crs(ii, ij) = zflcrs / zsfcrs
        END DO
      END DO
      !$ACC END KERNELS
      DEALLOCATE(zsurfmsk)
    CASE ('SUM')
      ALLOCATE(zsurfmsk(jpi, jpj))
      IF (PRESENT(p_e3)) THEN
        !$ACC KERNELS
        zsurfmsk(:, :) = p_e12(:, :) * p_e3(:, :, 1) * p_mask(:, :, 1)
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        zsurfmsk(:, :) = p_e12(:, :) * p_mask(:, :, 1)
        !$ACC END KERNELS
      END IF
      SELECT CASE (cd_type)
      CASE ('T', 'W')
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            !$ACC KERNELS
            je_2 = mje_crs(2)
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              zflcrs = p_fld(ji, je_2) * zsurfmsk(ji, je_2) + p_fld(ji + 1, je_2) * zsurfmsk(ji + 1, je_2) + p_fld(ji + 2, je_2) * zsurfmsk(ji + 2, je_2)
              p_fld_crs(ii, 2) = zflcrs
            END DO
            !$ACC END KERNELS
          END IF
        ELSE
          !$ACC KERNELS
          je_2 = mjs_crs(2)
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            zflcrs = p_fld(ji, je_2) * zsurfmsk(ji, je_2) + p_fld(ji + 1, je_2) * zsurfmsk(ji + 1, je_2) + p_fld(ji + 2, je_2) * zsurfmsk(ji + 2, je_2) + p_fld(ji, je_2 + 1) * zsurfmsk(ji, je_2 + 1) + p_fld(ji + 1, je_2 + 1) * zsurfmsk(ji + 1, je_2 + 1) + p_fld(ji + 2, je_2 + 1) * zsurfmsk(ji + 2, je_2 + 1) + p_fld(ji, je_2 + 2) * zsurfmsk(ji, je_2 + 2) + p_fld(ji + 1, je_2 + 2) * zsurfmsk(ji + 1, je_2 + 2) + p_fld(ji + 2, je_2 + 2) * zsurfmsk(ji + 2, je_2 + 2)
            p_fld_crs(ii, 2) = zflcrs
          END DO
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        !$ACC LOOP INDEPENDENT ! CDe
        DO jj = njstr, njend, nn_facty
          !$ACC LOOP INDEPENDENT ! CDe
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            ij = (jj - njstr) * rfacty_r + 3
            zflcrs = p_fld(ji, jj) * zsurfmsk(ji, jj) + p_fld(ji + 1, jj) * zsurfmsk(ji + 1, jj) + p_fld(ji + 2, jj) * zsurfmsk(ji + 2, jj) + p_fld(ji, jj + 1) * zsurfmsk(ji, jj + 1) + p_fld(ji + 1, jj + 1) * zsurfmsk(ji + 1, jj + 1) + p_fld(ji + 2, jj + 1) * zsurfmsk(ji + 2, jj + 1) + p_fld(ji, jj + 2) * zsurfmsk(ji, jj + 2) + p_fld(ji + 1, jj + 2) * zsurfmsk(ji + 1, jj + 2) + p_fld(ji + 2, jj + 2) * zsurfmsk(ji + 2, jj + 2)
            p_fld_crs(ii, ij) = zflcrs
          END DO
        END DO
        !$ACC END KERNELS
      CASE ('V')
        CALL profile_psy_data0 % PreStart('crs_dom_ope_2d', 'r0', 0, 0)
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            ijje = mje_crs(2)
          END IF
        ELSE
          ijje = mjs_crs(2)
        END IF
        CALL profile_psy_data0 % PostEnd
        !$ACC KERNELS
        DO ji = nistr, niend, nn_factx
          ii = (ji - mis_crs(2)) * rfactx_r + 2
          zflcrs = p_fld(ji, ijje) * zsurfmsk(ji, ijje) + p_fld(ji + 1, ijje) * zsurfmsk(ji + 1, ijje) + p_fld(ji + 2, ijje) * zsurfmsk(ji + 2, ijje)
          p_fld_crs(ii, 2) = zflcrs
        END DO
        !$ACC LOOP INDEPENDENT ! CDe
        DO jj = njstr, njend, nn_facty
          !$ACC LOOP INDEPENDENT ! CDe
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            ij = (jj - njstr) * rfacty_r + 3
            ijje = mje_crs(ij)
            zflcrs = p_fld(ji, ijje) * zsurfmsk(ji, ijje) + p_fld(ji + 1, ijje) * zsurfmsk(ji + 1, ijje) + p_fld(ji + 2, ijje) * zsurfmsk(ji + 2, ijje)
            p_fld_crs(ii, ij) = zflcrs
          END DO
        END DO
        !$ACC END KERNELS
      CASE ('U')
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            !$ACC KERNELS
            je_2 = mje_crs(2)
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              ijie = mie_crs(ii)
              zflcrs = p_fld(ijie, je_2) * zsurfmsk(ijie, je_2)
              p_fld_crs(ii, 2) = zflcrs
            END DO
            !$ACC END KERNELS
          END IF
        ELSE
          !$ACC KERNELS
          je_2 = mjs_crs(2)
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            ijie = mie_crs(ii)
            zflcrs = p_fld(ijie, je_2) * zsurfmsk(ijie, je_2) + p_fld(ijie, je_2 + 1) * zsurfmsk(ijie, je_2 + 1) + p_fld(ijie, je_2 + 2) * zsurfmsk(ijie, je_2 + 2)
            p_fld_crs(ii, 2) = zflcrs
          END DO
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        !$ACC LOOP INDEPENDENT ! CDe
        DO jj = njstr, njend, nn_facty
          !$ACC LOOP INDEPENDENT ! CDe
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            ij = (jj - njstr) * rfacty_r + 3
            ijie = mie_crs(ii)
            zflcrs = p_fld(ijie, jj) * zsurfmsk(ijie, jj) + p_fld(ijie, jj + 1) * zsurfmsk(ijie, jj + 1) + p_fld(ijie, jj + 2) * zsurfmsk(ijie, jj + 2)
            p_fld_crs(ii, ij) = zflcrs
          END DO
        END DO
        !$ACC END KERNELS
      END SELECT
      CALL profile_psy_data1 % PreStart('crs_dom_ope_2d', 'r1', 0, 0)
      IF (PRESENT(p_surf_crs)) THEN
        WHERE (p_surf_crs /= 0.0) p_fld_crs(:, :) = p_fld_crs(:, :) / p_surf_crs(:, :)
      END IF
      DEALLOCATE(zsurfmsk)
      CALL profile_psy_data1 % PostEnd
    CASE ('MAX')
      SELECT CASE (cd_type)
      CASE ('T', 'W')
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            !$ACC KERNELS
            je_2 = mje_crs(2)
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              zflcrs = MAX(p_fld(ji, je_2) * p_mask(ji, je_2, 1) - (1. - p_mask(ji, je_2, 1)) * r_inf, p_fld(ji + 1, je_2) * p_mask(ji + 1, je_2, 1) - (1. - p_mask(ji + 1, je_2, 1)) * r_inf, p_fld(ji + 2, je_2) * p_mask(ji + 2, je_2, 1) - (1. - p_mask(ji + 2, je_2, 1)) * r_inf)
              p_fld_crs(ii, 2) = zflcrs
            END DO
            !$ACC END KERNELS
          END IF
        ELSE
          CALL profile_psy_data2 % PreStart('crs_dom_ope_2d', 'r2', 0, 0)
          je_2 = mjs_crs(2)
          zflcrs = MAX(p_fld(ji, je_2) * p_mask(ji, je_2, 1) - (1. - p_mask(ji, je_2, 1)) * r_inf, p_fld(ji + 1, je_2) * p_mask(ji + 1, je_2, 1) - (1. - p_mask(ji + 1, je_2, 1)) * r_inf, p_fld(ji + 2, je_2) * p_mask(ji + 2, je_2, 1) - (1. - p_mask(ji + 2, je_2, 1)) * r_inf, p_fld(ji, je_2 + 1) * p_mask(ji, je_2 + 1, 1) - (1. - p_mask(ji, je_2 + 1, 1)) * r_inf, p_fld(ji + 1, je_2 + 1) * p_mask(ji + 1, je_2 + 1, 1) - (1. - p_mask(ji + 1, je_2 + 1, 1)) * r_inf, p_fld(ji + 2, je_2 + 1) * p_mask(ji + 2, je_2 + 1, 1) - (1. - p_mask(ji + 2, je_2 + 1, 1)) * r_inf, p_fld(ji, je_2 + 2) * p_mask(ji, je_2 + 2, 1) - (1. - p_mask(ji, je_2 + 2, 1)) * r_inf, p_fld(ji + 1, je_2 + 2) * p_mask(ji + 1, je_2 + 2, 1) - (1. - p_mask(ji + 1, je_2 + 2, 1)) * r_inf, p_fld(ji + 2, je_2 + 2) * p_mask(ji + 2, je_2 + 2, 1) - (1. - p_mask(ji + 2, je_2 + 2, 1)) * r_inf)
          p_fld_crs(ii, 2) = zflcrs
          CALL profile_psy_data2 % PostEnd
        END IF
        !$ACC KERNELS
        !$ACC LOOP INDEPENDENT ! CDe
        DO jj = njstr, njend, nn_facty
          !$ACC LOOP INDEPENDENT ! CDe
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            ij = (jj - njstr) * rfacty_r + 3
            zflcrs = MAX(p_fld(ji, jj) * p_mask(ji, jj, 1) - (1. - p_mask(ji, jj, 1)) * r_inf, p_fld(ji + 1, jj) * p_mask(ji + 1, jj, 1) - (1. - p_mask(ji + 1, jj, 1)) * r_inf, p_fld(ji + 2, jj) * p_mask(ji + 2, jj, 1) - (1. - p_mask(ji + 2, jj, 1)) * r_inf, p_fld(ji, jj + 1) * p_mask(ji, jj + 1, 1) - (1. - p_mask(ji, jj + 1, 1)) * r_inf, p_fld(ji + 1, jj + 1) * p_mask(ji + 1, jj + 1, 1) - (1. - p_mask(ji + 1, jj + 1, 1)) * r_inf, p_fld(ji + 2, jj + 1) * p_mask(ji + 2, jj + 1, 1) - (1. - p_mask(ji + 2, jj + 1, 1)) * r_inf, p_fld(ji, jj + 2) * p_mask(ji, jj + 2, 1) - (1. - p_mask(ji, jj + 2, 1)) * r_inf, p_fld(ji + 1, jj + 2) * p_mask(ji + 1, jj + 2, 1) - (1. - p_mask(ji + 1, jj + 2, 1)) * r_inf, p_fld(ji + 2, jj + 2) * p_mask(ji + 2, jj + 2, 1) - (1. - p_mask(ji + 2, jj + 2, 1)) * r_inf)
            p_fld_crs(ii, ij) = zflcrs
          END DO
        END DO
        !$ACC END KERNELS
      CASE ('V')
        CALL profile_psy_data3 % PreStart('crs_dom_ope_2d', 'r3', 0, 0)
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            ijje = mje_crs(2)
          END IF
        ELSE
          ijje = mjs_crs(2)
        END IF
        CALL profile_psy_data3 % PostEnd
        !$ACC KERNELS
        DO ji = nistr, niend, nn_factx
          ii = (ji - mis_crs(2)) * rfactx_r + 2
          zflcrs = MAX(p_fld(ji, ijje) * p_mask(ji, ijje, 1) - (1. - p_mask(ji, ijje, 1)) * r_inf, p_fld(ji + 1, ijje) * p_mask(ji + 1, ijje, 1) - (1. - p_mask(ji, ijje, 1)) * r_inf, p_fld(ji + 2, ijje) * p_mask(ji + 2, ijje, 1) - (1. - p_mask(ji, ijje, 1)) * r_inf)
          p_fld_crs(ii, 2) = zflcrs
        END DO
        !$ACC LOOP INDEPENDENT ! CDe
        DO jj = njstr, njend, nn_facty
          !$ACC LOOP INDEPENDENT ! CDe
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            ij = (jj - njstr) * rfacty_r + 3
            ijje = mje_crs(ij)
            zflcrs = MAX(p_fld(ji, ijje) * p_mask(ji, ijje, 1) - (1. - p_mask(ji, ijje, 1)) * r_inf, p_fld(ji + 1, ijje) * p_mask(ji + 1, ijje, 1) - (1. - p_mask(ji, ijje, 1)) * r_inf, p_fld(ji + 2, ijje) * p_mask(ji + 2, ijje, 1) - (1. - p_mask(ji, ijje, 1)) * r_inf)
            p_fld_crs(ii, ij) = zflcrs
          END DO
        END DO
        !$ACC END KERNELS
      CASE ('U')
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            !$ACC KERNELS
            je_2 = mje_crs(2)
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              ijie = mie_crs(ii)
              zflcrs = p_fld(ijie, je_2) * p_mask(ijie, je_2, 1) - (1. - p_mask(ijie, je_2, 1)) * r_inf
              p_fld_crs(ii, 2) = zflcrs
            END DO
            !$ACC END KERNELS
          END IF
        ELSE
          !$ACC KERNELS
          je_2 = mjs_crs(2)
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            ijie = mie_crs(ii)
            zflcrs = MAX(p_fld(ijie, je_2) * p_mask(ijie, je_2, 1) - (1. - p_mask(ijie, je_2, 1)) * r_inf, p_fld(ijie, je_2 + 1) * p_mask(ijie, je_2 + 1, 1) - (1. - p_mask(ijie, je_2, 1)) * r_inf, p_fld(ijie, je_2 + 2) * p_mask(ijie, je_2 + 2, 1) - (1. - p_mask(ijie, je_2, 1)) * r_inf)
            p_fld_crs(ii, 2) = zflcrs
          END DO
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        !$ACC LOOP INDEPENDENT ! CDe
        DO jj = njstr, njend, nn_facty
          !$ACC LOOP INDEPENDENT ! CDe
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            ij = (jj - njstr) * rfacty_r + 3
            ijie = mie_crs(ii)
            zflcrs = MAX(p_fld(ijie, jj) * p_mask(ijie, jj, 1) - (1. - p_mask(ijie, jj, 1)) * r_inf, p_fld(ijie, jj + 1) * p_mask(ijie, jj + 1, 1) - (1. - p_mask(ijie, jj, 1)) * r_inf, p_fld(ijie, jj + 2) * p_mask(ijie, jj + 2, 1) - (1. - p_mask(ijie, jj, 1)) * r_inf)
            p_fld_crs(ii, ij) = zflcrs
          END DO
        END DO
        !$ACC END KERNELS
      END SELECT
    CASE ('MIN')
      SELECT CASE (cd_type)
      CASE ('T', 'W')
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            !$ACC KERNELS
            je_2 = mje_crs(2)
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              zflcrs = MIN(p_fld(ji, je_2) * p_mask(ji, je_2, 1) + (1. - p_mask(ji, je_2, 1)) * r_inf, p_fld(ji + 1, je_2) * p_mask(ji + 1, je_2, 1) + (1. - p_mask(ji + 1, je_2, 1)) * r_inf, p_fld(ji + 2, je_2) * p_mask(ji + 2, je_2, 1) + (1. - p_mask(ji + 2, je_2, 1)) * r_inf)
              p_fld_crs(ii, 2) = zflcrs
            END DO
            !$ACC END KERNELS
          END IF
        ELSE
          CALL profile_psy_data4 % PreStart('crs_dom_ope_2d', 'r4', 0, 0)
          je_2 = mjs_crs(2)
          zflcrs = MIN(p_fld(ji, je_2) * p_mask(ji, je_2, 1) + (1. - p_mask(ji, je_2, 1)) * r_inf, p_fld(ji + 1, je_2) * p_mask(ji + 1, je_2, 1) + (1. - p_mask(ji + 1, je_2, 1)) * r_inf, p_fld(ji + 2, je_2) * p_mask(ji + 2, je_2, 1) + (1. - p_mask(ji + 2, je_2, 1)) * r_inf, p_fld(ji, je_2 + 1) * p_mask(ji, je_2 + 1, 1) + (1. - p_mask(ji, je_2 + 1, 1)) * r_inf, p_fld(ji + 1, je_2 + 1) * p_mask(ji + 1, je_2 + 1, 1) + (1. - p_mask(ji + 1, je_2 + 1, 1)) * r_inf, p_fld(ji + 2, je_2 + 1) * p_mask(ji + 2, je_2 + 1, 1) + (1. - p_mask(ji + 2, je_2 + 1, 1)) * r_inf, p_fld(ji, je_2 + 2) * p_mask(ji, je_2 + 2, 1) + (1. - p_mask(ji, je_2 + 2, 1)) * r_inf, p_fld(ji + 1, je_2 + 2) * p_mask(ji + 1, je_2 + 2, 1) + (1. - p_mask(ji + 1, je_2 + 2, 1)) * r_inf, p_fld(ji + 2, je_2 + 2) * p_mask(ji + 2, je_2 + 2, 1) + (1. - p_mask(ji + 2, je_2 + 2, 1)) * r_inf)
          p_fld_crs(ii, 2) = zflcrs
          CALL profile_psy_data4 % PostEnd
        END IF
        !$ACC KERNELS
        !$ACC LOOP INDEPENDENT ! CDe
        DO jj = njstr, njend, nn_facty
          !$ACC LOOP INDEPENDENT ! CDe
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            ij = (jj - njstr) * rfacty_r + 3
            zflcrs = MIN(p_fld(ji, jj) * p_mask(ji, jj, 1) + (1. - p_mask(ji, jj, 1)) * r_inf, p_fld(ji + 1, jj) * p_mask(ji + 1, jj, 1) + (1. - p_mask(ji + 1, jj, 1)) * r_inf, p_fld(ji + 2, jj) * p_mask(ji + 2, jj, 1) + (1. - p_mask(ji + 2, jj, 1)) * r_inf, p_fld(ji, jj + 1) * p_mask(ji, jj + 1, 1) + (1. - p_mask(ji, jj + 1, 1)) * r_inf, p_fld(ji + 1, jj + 1) * p_mask(ji + 1, jj + 1, 1) + (1. - p_mask(ji + 1, jj + 1, 1)) * r_inf, p_fld(ji + 2, jj + 1) * p_mask(ji + 2, jj + 1, 1) + (1. - p_mask(ji + 2, jj + 1, 1)) * r_inf, p_fld(ji, jj + 2) * p_mask(ji, jj + 2, 1) + (1. - p_mask(ji, jj + 2, 1)) * r_inf, p_fld(ji + 1, jj + 2) * p_mask(ji + 1, jj + 2, 1) + (1. - p_mask(ji + 1, jj + 2, 1)) * r_inf, p_fld(ji + 2, jj + 2) * p_mask(ji + 2, jj + 2, 1) + (1. - p_mask(ji + 2, jj + 2, 1)) * r_inf)
            p_fld_crs(ii, ij) = zflcrs
          END DO
        END DO
        !$ACC END KERNELS
      CASE ('V')
        CALL profile_psy_data5 % PreStart('crs_dom_ope_2d', 'r5', 0, 0)
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            ijje = mje_crs(2)
          END IF
        ELSE
          ijje = mjs_crs(2)
        END IF
        CALL profile_psy_data5 % PostEnd
        !$ACC KERNELS
        DO ji = nistr, niend, nn_factx
          ii = (ji - mis_crs(2)) * rfactx_r + 2
          zflcrs = MIN(p_fld(ji, ijje) * p_mask(ji, ijje, 1) + (1. - p_mask(ji, ijje, 1)) * r_inf, p_fld(ji + 1, ijje) * p_mask(ji + 1, ijje, 1) + (1. - p_mask(ji, ijje, 1)) * r_inf, p_fld(ji + 2, ijje) * p_mask(ji + 2, ijje, 1) + (1. - p_mask(ji, ijje, 1)) * r_inf)
          p_fld_crs(ii, 2) = zflcrs
        END DO
        !$ACC LOOP INDEPENDENT ! CDe
        DO jj = njstr, njend, nn_facty
          !$ACC LOOP INDEPENDENT ! CDe
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            ij = (jj - njstr) * rfacty_r + 3
            ijje = mje_crs(ij)
            zflcrs = MIN(p_fld(ji, ijje) * p_mask(ji, ijje, 1) + (1. - p_mask(ji, ijje, 1)) * r_inf, p_fld(ji + 1, ijje) * p_mask(ji + 1, ijje, 1) + (1. - p_mask(ji, ijje, 1)) * r_inf, p_fld(ji + 2, ijje) * p_mask(ji + 2, ijje, 1) + (1. - p_mask(ji, ijje, 1)) * r_inf)
            p_fld_crs(ii, ij) = zflcrs
          END DO
        END DO
        !$ACC END KERNELS
      CASE ('U')
        IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
          IF (mje_crs(2) - mjs_crs(2) == 1) THEN
            !$ACC KERNELS
            je_2 = mje_crs(2)
            DO ji = nistr, niend, nn_factx
              ii = (ji - mis_crs(2)) * rfactx_r + 2
              ijie = mie_crs(ii)
              zflcrs = p_fld(ijie, je_2) * p_mask(ijie, je_2, 1) + (1. - p_mask(ijie, je_2, 1)) * r_inf
              p_fld_crs(ii, 2) = zflcrs
            END DO
            !$ACC END KERNELS
          END IF
        ELSE
          !$ACC KERNELS
          je_2 = mjs_crs(2)
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            ijie = mie_crs(ii)
            zflcrs = MIN(p_fld(ijie, je_2) * p_mask(ijie, je_2, 1) + (1. - p_mask(ijie, je_2, 1)) * r_inf, p_fld(ijie, je_2 + 1) * p_mask(ijie, je_2 + 1, 1) + (1. - p_mask(ijie, je_2, 1)) * r_inf, p_fld(ijie, je_2 + 2) * p_mask(ijie, je_2 + 2, 1) + (1. - p_mask(ijie, je_2, 1)) * r_inf)
            p_fld_crs(ii, 2) = zflcrs
          END DO
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        !$ACC LOOP INDEPENDENT ! CDe
        DO jj = njstr, njend, nn_facty
          !$ACC LOOP INDEPENDENT ! CDe
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            ij = (jj - njstr) * rfacty_r + 3
            ijie = mie_crs(ii)
            zflcrs = MIN(p_fld(ijie, jj) * p_mask(ijie, jj, 1) + (1. - p_mask(ijie, jj, 1)) * r_inf, p_fld(ijie, jj + 1) * p_mask(ijie, jj + 1, 1) + (1. - p_mask(ijie, jj, 1)) * r_inf, p_fld(ijie, jj + 2) * p_mask(ijie, jj + 2, 1) + (1. - p_mask(ijie, jj, 1)) * r_inf)
            p_fld_crs(ii, ij) = zflcrs
          END DO
        END DO
        !$ACC END KERNELS
      END SELECT
    END SELECT
    CALL profile_psy_data6 % PreStart('crs_dom_ope_2d', 'r6', 0, 0)
    CALL crs_lbc_lnk(p_fld_crs, cd_type, psgn)
    CALL profile_psy_data6 % PostEnd
  END SUBROUTINE crs_dom_ope_2d
  SUBROUTINE crs_dom_e3(p_e1, p_e2, p_e3, p_sfc_crs, cd_type, p_mask, p_e3_crs, p_e3_max_crs)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: p_mask
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: p_e1, p_e2
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: p_e3
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs, jpk), INTENT(IN) :: p_sfc_crs
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs, jpk), INTENT(INOUT) :: p_e3_crs
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs, jpk), INTENT(INOUT) :: p_e3_max_crs
    INTEGER :: ji, jj, jk
    INTEGER :: ijie, ijje, ii, ij, je_2
    REAL(KIND = wp) :: ze3crs
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zmask, zsurf
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    p_e3_crs(:, :, :) = 0.
    p_e3_max_crs(:, :, :) = 1.
    SELECT CASE (cd_type)
    CASE ('W')
      zmask(:, :, 1) = p_mask(:, :, 1)
      DO jk = 2, jpk
        zmask(:, :, jk) = p_mask(:, :, jk - 1)
      END DO
    CASE DEFAULT
      DO jk = 1, jpk
        zmask(:, :, jk) = p_mask(:, :, jk)
      END DO
    END SELECT
    DO jk = 1, jpk
      zsurf(:, :, jk) = p_e1(:, :) * p_e2(:, :) * p_e3(:, :, jk)
    END DO
    !$ACC END KERNELS
    IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
      IF (mje_crs(2) - mjs_crs(2) == 1) THEN
        !$ACC KERNELS
        je_2 = mje_crs(2)
        DO jk = 1, jpk
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            ze3crs = zsurf(ji, je_2, jk) * zmask(ji, je_2, jk) + zsurf(ji + 1, je_2, jk) * zmask(ji + 1, je_2, jk) + zsurf(ji + 2, je_2, jk) * zmask(ji + 2, je_2, jk)
            p_e3_crs(ii, 2, jk) = ze3crs / p_sfc_crs(ii, ij, jk)
            ze3crs = MAX(p_e3(ji, je_2, jk) * zmask(ji, je_2, jk), p_e3(ji + 1, je_2, jk) * zmask(ji + 1, je_2, jk), p_e3(ji + 2, je_2, jk) * zmask(ji + 2, je_2, jk))
            p_e3_max_crs(ii, 2, jk) = ze3crs
          END DO
        END DO
        !$ACC END KERNELS
      END IF
    ELSE
      !$ACC KERNELS
      je_2 = mjs_crs(2)
      DO jk = 1, jpk
        DO ji = nistr, niend, nn_factx
          ii = (ji - mis_crs(2)) * rfactx_r + 2
          ze3crs = zsurf(ji, je_2, jk) * zmask(ji, je_2, jk) + zsurf(ji + 1, je_2, jk) * zmask(ji + 1, je_2, jk) + zsurf(ji + 2, je_2, jk) * zmask(ji + 2, je_2, jk) + zsurf(ji, je_2 + 1, jk) * zmask(ji, je_2 + 1, jk) + zsurf(ji + 1, je_2 + 1, jk) * zmask(ji + 1, je_2 + 1, jk) + zsurf(ji + 2, je_2 + 1, jk) * zmask(ji + 2, je_2 + 1, jk) + zsurf(ji, je_2 + 2, jk) * zmask(ji, je_2 + 2, jk) + zsurf(ji + 1, je_2 + 2, jk) * zmask(ji + 1, je_2 + 2, jk) + zsurf(ji + 2, je_2 + 2, jk) * zmask(ji + 2, je_2 + 2, jk)
          p_e3_crs(ii, 2, jk) = ze3crs / p_sfc_crs(ii, 2, jk)
          ze3crs = MAX(p_e3(ji, je_2, jk) * zmask(ji, je_2, jk), p_e3(ji + 1, je_2, jk) * zmask(ji + 1, je_2, jk), p_e3(ji + 2, je_2, jk) * zmask(ji + 2, je_2, jk), p_e3(ji, je_2 + 1, jk) * zmask(ji, je_2 + 1, jk), p_e3(ji + 1, je_2 + 1, jk) * zmask(ji + 1, je_2 + 1, jk), p_e3(ji + 2, je_2 + 1, jk) * zmask(ji + 2, je_2 + 1, jk), p_e3(ji, je_2 + 2, jk) * zmask(ji, je_2 + 2, jk), p_e3(ji + 1, je_2 + 2, jk) * zmask(ji + 1, je_2 + 2, jk), p_e3(ji + 2, je_2 + 2, jk) * zmask(ji + 2, je_2 + 2, jk))
          p_e3_max_crs(ii, 2, jk) = ze3crs
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    DO jk = 1, jpk
      DO jj = njstr, njend, nn_facty
        DO ji = nistr, niend, nn_factx
          ii = (ji - mis_crs(2)) * rfactx_r + 2
          ij = (jj - njstr) * rfacty_r + 3
          ze3crs = zsurf(ji, jj, jk) * zmask(ji, jj, jk) + zsurf(ji + 1, jj, jk) * zmask(ji + 1, jj, jk) + zsurf(ji + 2, jj, jk) * zmask(ji + 2, jj, jk) + zsurf(ji, jj + 1, jk) * zmask(ji, jj + 1, jk) + zsurf(ji + 1, jj + 1, jk) * zmask(ji + 1, jj + 1, jk) + zsurf(ji + 2, jj + 1, jk) * zmask(ji + 2, jj + 1, jk) + zsurf(ji, jj + 2, jk) * zmask(ji, jj + 2, jk) + zsurf(ji + 1, jj + 2, jk) * zmask(ji + 1, jj + 2, jk) + zsurf(ji + 2, jj + 2, jk) * zmask(ji + 2, jj + 2, jk)
          p_e3_crs(ii, ij, jk) = ze3crs / p_sfc_crs(ii, ij, jk)
          ze3crs = MAX(p_e3(ji, jj, jk) * zmask(ji, jj, jk), p_e3(ji + 1, jj, jk) * zmask(ji + 1, jj, jk), p_e3(ji + 2, jj, jk) * zmask(ji + 2, jj, jk), p_e3(ji, jj + 1, jk) * zmask(ji, jj + 1, jk), p_e3(ji + 1, jj + 1, jk) * zmask(ji + 1, jj + 1, jk), p_e3(ji + 2, jj + 1, jk) * zmask(ji + 2, jj + 1, jk), p_e3(ji, jj + 2, jk) * zmask(ji, jj + 2, jk), p_e3(ji + 1, jj + 2, jk) * zmask(ji + 1, jj + 2, jk), p_e3(ji + 2, jj + 2, jk) * zmask(ji + 2, jj + 2, jk))
          p_e3_max_crs(ii, ij, jk) = ze3crs
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('crs_dom_e3', 'r0', 0, 0)
    CALL crs_lbc_lnk(p_e3_crs, cd_type, 1.0, pval = 1.0)
    CALL crs_lbc_lnk(p_e3_max_crs, cd_type, 1.0, pval = 1.0)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE crs_dom_e3
  SUBROUTINE crs_dom_sfc(p_mask, cd_type, p_surf_crs, p_surf_crs_msk, p_e1, p_e2, p_e3)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: p_mask
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN), OPTIONAL :: p_e1, p_e2
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN), OPTIONAL :: p_e3
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs, jpk), INTENT(OUT) :: p_surf_crs
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs, jpk), INTENT(OUT) :: p_surf_crs_msk
    INTEGER :: ji, jj, jk
    INTEGER :: ii, ij, je_2
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zsurf, zsurfmsk
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    SELECT CASE (cd_type)
    CASE ('W')
      DO jk = 1, jpk
        zsurf(:, :, jk) = p_e1(:, :) * p_e2(:, :)
      END DO
      zsurfmsk(:, :, 1) = zsurf(:, :, 1) * p_mask(:, :, 1)
      DO jk = 2, jpk
        zsurfmsk(:, :, jk) = zsurf(:, :, jk) * p_mask(:, :, jk - 1)
      END DO
    CASE ('V')
      DO jk = 1, jpk
        zsurf(:, :, jk) = p_e1(:, :) * p_e3(:, :, jk)
      END DO
      DO jk = 1, jpk
        zsurfmsk(:, :, jk) = zsurf(:, :, jk) * p_mask(:, :, jk)
      END DO
    CASE ('U')
      DO jk = 1, jpk
        zsurf(:, :, jk) = p_e2(:, :) * p_e3(:, :, jk)
      END DO
      DO jk = 1, jpk
        zsurfmsk(:, :, jk) = zsurf(:, :, jk) * p_mask(:, :, jk)
      END DO
    CASE DEFAULT
      DO jk = 1, jpk
        zsurf(:, :, jk) = p_e1(:, :) * p_e2(:, :)
      END DO
      DO jk = 1, jpk
        zsurfmsk(:, :, jk) = zsurf(:, :, jk) * p_mask(:, :, jk)
      END DO
    END SELECT
    !$ACC END KERNELS
    IF (nldj_crs == 1 .AND. ((mje_crs(2) - mjs_crs(2)) < 2)) THEN
      IF (mje_crs(2) - mjs_crs(2) == 1) THEN
        !$ACC KERNELS
        je_2 = mje_crs(2)
        DO jk = 1, jpk
          DO ji = nistr, niend, nn_factx
            ii = (ji - mis_crs(2)) * rfactx_r + 2
            p_surf_crs(ii, 2, jk) = zsurf(ji, je_2, jk) + zsurf(ji + 1, je_2, jk) + zsurf(ji + 2, je_2, jk) + zsurf(ji, je_2 - 1, jk) + zsurf(ji + 1, je_2 - 1, jk) + zsurf(ji + 2, je_2 - 1, jk)
            p_surf_crs_msk(ii, 2, jk) = zsurfmsk(ji, je_2, jk) + zsurfmsk(ji + 1, je_2, jk) + zsurfmsk(ji + 2, je_2, jk)
          END DO
        END DO
        !$ACC END KERNELS
      END IF
    ELSE
      !$ACC KERNELS
      je_2 = mjs_crs(2)
      DO jk = 1, jpk
        DO ji = nistr, niend, nn_factx
          ii = (ji - mis_crs(2)) * rfactx_r + 2
          p_surf_crs(ii, 2, jk) = zsurf(ji, je_2, jk) + zsurf(ji + 1, je_2, jk) + zsurf(ji + 2, je_2, jk) + zsurf(ji, je_2 + 1, jk) + zsurf(ji + 1, je_2 + 1, jk) + zsurf(ji + 2, je_2 + 1, jk) + zsurf(ji, je_2 + 2, jk) + zsurf(ji + 1, je_2 + 2, jk) + zsurf(ji + 2, je_2 + 2, jk)
          p_surf_crs_msk(ii, 2, jk) = zsurfmsk(ji, je_2, jk) + zsurfmsk(ji + 1, je_2, jk) + zsurfmsk(ji + 2, je_2, jk) + zsurfmsk(ji, je_2 + 1, jk) + zsurfmsk(ji + 1, je_2 + 1, jk) + zsurfmsk(ji + 2, je_2 + 1, jk) + zsurfmsk(ji, je_2 + 2, jk) + zsurfmsk(ji + 1, je_2 + 2, jk) + zsurfmsk(ji + 2, je_2 + 2, jk)
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    DO jk = 1, jpk
      DO jj = njstr, njend, nn_facty
        DO ji = nistr, niend, nn_factx
          ii = (ji - mis_crs(2)) * rfactx_r + 2
          ij = (jj - njstr) * rfacty_r + 3
          p_surf_crs(ii, ij, jk) = zsurf(ji, jj, jk) + zsurf(ji + 1, jj, jk) + zsurf(ji + 2, jj, jk) + zsurf(ji, jj + 1, jk) + zsurf(ji + 1, jj + 1, jk) + zsurf(ji + 2, jj + 1, jk) + zsurf(ji, jj + 2, jk) + zsurf(ji + 1, jj + 2, jk) + zsurf(ji + 2, jj + 2, jk)
          p_surf_crs_msk(ii, ij, jk) = zsurfmsk(ji, jj, jk) + zsurfmsk(ji + 1, jj, jk) + zsurfmsk(ji + 2, jj, jk) + zsurfmsk(ji, jj + 1, jk) + zsurfmsk(ji + 1, jj + 1, jk) + zsurfmsk(ji + 2, jj + 1, jk) + zsurfmsk(ji, jj + 2, jk) + zsurfmsk(ji + 1, jj + 2, jk) + zsurfmsk(ji + 2, jj + 2, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('crs_dom_sfc', 'r0', 0, 0)
    CALL crs_lbc_lnk(p_surf_crs, cd_type, 1.0, pval = 1.0)
    CALL crs_lbc_lnk(p_surf_crs_msk, cd_type, 1.0, pval = 1.0)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE crs_dom_sfc
  SUBROUTINE crs_dom_def
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jj, jk, ijjgloT, ijis, ijie, ijjs, ijje, jn
    INTEGER :: ierr
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
    CALL profile_psy_data0 % PreStart('crs_dom_def', 'r0', 0, 0)
    jpiglo_crs = INT((jpiglo - 2) / nn_factx) + 2
    jpjglo_crs = INT((jpjglo - MOD(jpjglo, nn_facty)) / nn_facty) + 3
    jpiglo_crsm1 = jpiglo_crs - 1
    jpjglo_crsm1 = jpjglo_crs - 1
    jpi_crs = (jpiglo_crs - 2 * nn_hls + (jpni - 1)) / jpni + 2 * nn_hls
    jpj_crs = (jpjglo_crsm1 - 2 * nn_hls + (jpnj - 1)) / jpnj + 2 * nn_hls
    IF (noso < 0) jpj_crs = jpj_crs + 1
    jpi_crsm1 = jpi_crs - 1
    jpj_crsm1 = jpj_crs - 1
    nperio_crs = jperio
    npolj_crs = npolj
    ierr = crs_dom_alloc()
    CALL profile_psy_data0 % PostEnd
    IF (.NOT. lk_mpp) THEN
      CALL profile_psy_data1 % PreStart('crs_dom_def', 'r1', 0, 0)
      nimpp_crs = 1
      njmpp_crs = 1
      nlci_crs = jpi_crs
      nlcj_crs = jpj_crs
      nldi_crs = 1
      nldj_crs = 1
      nlei_crs = jpi_crs
      nlej_crs = jpj_crs
      CALL profile_psy_data1 % PostEnd
    ELSE
      CALL profile_psy_data2 % PreStart('crs_dom_def', 'r2', 0, 0)
      nimpp_crs = 1
      njmpp_crs = 1
      nlci_crs = jpi_crs
      nlcj_crs = jpj_crs
      nldi_crs = 1
      nldj_crs = 1
      nlei_crs = jpi_crs
      nlej_crs = jpj_crs
      CALL profile_psy_data2 % PostEnd
      DO jn = 1, jpnij, jpni
        CALL profile_psy_data3 % PreStart('crs_dom_def', 'r3', 0, 0)
        IF (jn < (jpnij - jpni + 1)) THEN
          nlejt_crs(jn) = AINT(REAL((jpjglo - (njmppt(jn) - 1)) / nn_facty, wp)) - AINT(REAL((jpjglo - (njmppt(jn + jpni) - 1)) / nn_facty, wp))
        ELSE
          nlejt_crs(jn) = AINT(REAL(nlejt(jn) / nn_facty, wp)) + 1
        END IF
        IF (noso < 0) nlejt_crs(jn) = nlejt_crs(jn) + 1
        SELECT CASE (ibonjt(jn))
        CASE (- 1)
          IF (MOD(jpjglo - njmppt(jn), nn_facty) > 0) nlejt_crs(jn) = nlejt_crs(jn) + 1
          nlcjt_crs(jn) = nlejt_crs(jn) + nn_hls
          nldjt_crs(jn) = nldjt(jn)
        CASE (0)
          nldjt_crs(jn) = nldjt(jn)
          IF (nldjt(jn) == 1) nlejt_crs(jn) = nlejt_crs(jn) + 1
          nlejt_crs(jn) = nlejt_crs(jn) + nn_hls
          nlcjt_crs(jn) = nlejt_crs(jn) + nn_hls
        CASE (1, 2)
          nlejt_crs(jn) = nlejt_crs(jn) + nn_hls
          nlcjt_crs(jn) = nlejt_crs(jn)
          nldjt_crs(jn) = nldjt(jn)
        CASE DEFAULT
          STOP
        END SELECT
        IF (nlcjt_crs(jn) > jpj_crs) jpj_crs = jpj_crs + 1
        IF (nldjt_crs(jn) == 1) THEN
          njmppt_crs(jn) = 1
        ELSE
          njmppt_crs(jn) = 2 + ANINT(REAL((njmppt(jn) + 1 - MOD(jpjglo, nn_facty)) / nn_facty, wp))
        END IF
        CALL profile_psy_data3 % PostEnd
        !$ACC KERNELS
        DO jj = jn + 1, jn + jpni - 1
          nlejt_crs(jj) = nlejt_crs(jn)
          nlcjt_crs(jj) = nlcjt_crs(jn)
          nldjt_crs(jj) = nldjt_crs(jn)
          njmppt_crs(jj) = njmppt_crs(jn)
        END DO
        !$ACC END KERNELS
      END DO
      CALL profile_psy_data4 % PreStart('crs_dom_def', 'r4', 0, 0)
      nlej_crs = nlejt_crs(nproc + 1)
      nlcj_crs = nlcjt_crs(nproc + 1)
      nldj_crs = nldjt_crs(nproc + 1)
      njmpp_crs = njmppt_crs(nproc + 1)
      CALL profile_psy_data4 % PostEnd
      DO jn = 1, jpni
        CALL profile_psy_data5 % PreStart('crs_dom_def', 'r5', 0, 0)
        IF (jn == 1) THEN
          nleit_crs(jn) = AINT(REAL((nimppt(jn) - 1 + nlcit(jn)) / nn_factx, wp))
        ELSE
          nleit_crs(jn) = AINT(REAL((nimppt(jn) - 1 + nlcit(jn)) / nn_factx, wp)) - AINT(REAL((nimppt(jn - 1) - 1 + nlcit(jn - 1)) / nn_factx, wp))
        END IF
        SELECT CASE (ibonit(jn))
        CASE (- 1)
          nleit_crs(jn) = nleit_crs(jn) + nn_hls
          nlcit_crs(jn) = nleit_crs(jn) + nn_hls
          nldit_crs(jn) = nldit(jn)
        CASE (0)
          nleit_crs(jn) = nleit_crs(jn) + nn_hls
          nlcit_crs(jn) = nleit_crs(jn) + nn_hls
          nldit_crs(jn) = nldit(jn)
        CASE (1, 2)
          IF (MOD(jpiglo - nimppt(jn), nn_factx) > 0) nleit_crs(jn) = nleit_crs(jn) + 1
          nleit_crs(jn) = nleit_crs(jn) + nn_hls
          nlcit_crs(jn) = nleit_crs(jn)
          nldit_crs(jn) = nldit(jn)
        CASE DEFAULT
          STOP
        END SELECT
        nimppt_crs(jn) = ANINT(REAL((nimppt(jn) + 1) / nn_factx, wp)) + 1
        CALL profile_psy_data5 % PostEnd
        !$ACC KERNELS
        DO jj = jn + jpni, jpnij, jpni
          nleit_crs(jj) = nleit_crs(jn)
          nlcit_crs(jj) = nlcit_crs(jn)
          nldit_crs(jj) = nldit_crs(jn)
          nimppt_crs(jj) = nimppt_crs(jn)
        END DO
        !$ACC END KERNELS
      END DO
      !$ACC KERNELS
      nlei_crs = nleit_crs(nproc + 1)
      nlci_crs = nlcit_crs(nproc + 1)
      nldi_crs = nldit_crs(nproc + 1)
      nimpp_crs = nimppt_crs(nproc + 1)
      DO ji = 1, jpi_crs
        mig_crs(ji) = ji + nimpp_crs - 1
      END DO
      DO jj = 1, jpj_crs
        mjg_crs(jj) = jj + njmpp_crs - 1
      END DO
      DO ji = 1, jpiglo_crs
        mi0_crs(ji) = MAX(1, MIN(ji - nimpp_crs + 1, jpi_crs + 1))
        mi1_crs(ji) = MAX(0, MIN(ji - nimpp_crs + 1, jpi_crs))
      END DO
      DO jj = 1, jpjglo_crs
        mj0_crs(jj) = MAX(1, MIN(jj - njmpp_crs + 1, jpj_crs + 1))
        mj1_crs(jj) = MAX(0, MIN(jj - njmpp_crs + 1, jpj_crs))
      END DO
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    jpi_full = jpi
    jpj_full = jpj
    jpim1_full = jpim1
    jpjm1_full = jpjm1
    nperio_full = jperio
    npolj_full = npolj
    jpiglo_full = jpiglo
    jpjglo_full = jpjglo
    nlcj_full = nlcj
    nlci_full = nlci
    nldi_full = nldi
    nldj_full = nldj
    nlei_full = nlei
    nlej_full = nlej
    nimpp_full = nimpp
    njmpp_full = njmpp
    nlcit_full(:) = nlcit(:)
    nldit_full(:) = nldit(:)
    nleit_full(:) = nleit(:)
    nimppt_full(:) = nimppt(:)
    nlcjt_full(:) = nlcjt(:)
    nldjt_full(:) = nldjt(:)
    nlejt_full(:) = nlejt(:)
    njmppt_full(:) = njmppt(:)
    !$ACC END KERNELS
    CALL profile_psy_data6 % PreStart('crs_dom_def', 'r6', 0, 0)
    CALL dom_grid_crs
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'crs_init : coarse grid dimensions'
      WRITE(numout, *) '~~~~~~~   coarse domain global j-dimension           jpjglo = ', jpjglo
      WRITE(numout, *) '~~~~~~~   coarse domain global i-dimension           jpiglo = ', jpiglo
      WRITE(numout, *) '~~~~~~~   coarse domain local  i-dimension              jpi = ', jpi
      WRITE(numout, *) '~~~~~~~   coarse domain local  j-dimension              jpj = ', jpj
      WRITE(numout, *)
      WRITE(numout, *) ' nproc  = ', nproc
      WRITE(numout, *) ' nlci   = ', nlci
      WRITE(numout, *) ' nlcj   = ', nlcj
      WRITE(numout, *) ' nldi   = ', nldi
      WRITE(numout, *) ' nldj   = ', nldj
      WRITE(numout, *) ' nlei   = ', nlei
      WRITE(numout, *) ' nlej   = ', nlej
      WRITE(numout, *) ' nlei_full=', nlei_full
      WRITE(numout, *) ' nldi_full=', nldi_full
      WRITE(numout, *) ' nimpp  = ', nimpp
      WRITE(numout, *) ' njmpp  = ', njmpp
      WRITE(numout, *) ' njmpp_full  = ', njmpp_full
      WRITE(numout, *)
    END IF
    CALL dom_grid_glo
    mxbinctr = INT(nn_factx * 0.5)
    mybinctr = INT(nn_facty * 0.5)
    nrestx = MOD(nn_factx, 2)
    nresty = MOD(nn_facty, 2)
    IF (nrestx == 0) THEN
      mxbinctr = mxbinctr - 1
    END IF
    IF (nresty == 0) THEN
      mybinctr = mybinctr - 1
      IF (jperio == 3 .OR. jperio == 4) nperio_crs = jperio + 2
      IF (jperio == 5 .OR. jperio == 6) nperio_crs = jperio - 2
      IF (npolj == 3) npolj_crs = 5
      IF (npolj == 5) npolj_crs = 3
    END IF
    rfactxy = nn_factx * nn_facty
    ierr = crs_dom_alloc2()
    CALL profile_psy_data6 % PostEnd
    !$ACC KERNELS
    mis2_crs(:) = 0
    mie2_crs(:) = 0
    mjs2_crs(:) = 0
    mje2_crs(:) = 0
    !$ACC END KERNELS
    SELECT CASE (nn_binref)
    CASE (0)
      SELECT CASE (jperio)
      CASE (0, 1, 3, 4)
        !$ACC KERNELS
        DO ji = 2, jpiglo_crsm1
          ijie = (ji * nn_factx) - nn_factx
          ijis = ijie - nn_factx + 1
          mis2_crs(ji) = ijis
          mie2_crs(ji) = ijie
        END DO
        !$ACC END KERNELS
        CALL profile_psy_data7 % PreStart('crs_dom_def', 'r7', 0, 0)
        IF (jpiglo - 1 - mie2_crs(jpiglo_crsm1) <= nn_factx) mie2_crs(jpiglo_crsm1) = jpiglo - 2
        CALL profile_psy_data7 % PostEnd
        !$ACC KERNELS
        IF (nn_facty == 2) THEN
          ijjglot = jpjglo - 1
        ELSE
          ijjglot = jpjglo
        END IF
        DO jj = 2, jpjglo_crs
          ijje = ijjgloT - nn_facty * (jj - 3)
          ijjs = ijje - nn_facty + 1
          mjs2_crs(jpjglo_crs - jj + 2) = ijjs
          mje2_crs(jpjglo_crs - jj + 2) = ijje
        END DO
        !$ACC END KERNELS
      CASE (2)
        WRITE(numout, *) 'crs_init, jperio=2 not supported'
      CASE (5, 6)
        !$ACC KERNELS
        DO ji = 2, jpiglo_crsm1
          ijie = (ji * nn_factx) - nn_factx
          ijis = ijie - nn_factx + 1
          mis2_crs(ji) = ijis
          mie2_crs(ji) = ijie
        END DO
        !$ACC END KERNELS
        CALL profile_psy_data8 % PreStart('crs_dom_def', 'r8', 0, 0)
        IF (jpiglo - 1 - mie2_crs(jpiglo_crsm1) <= nn_factx) mie_crs(jpiglo_crsm1) = jpiglo - 2
        CALL profile_psy_data8 % PostEnd
        !$ACC KERNELS
        jj = 2
        ijje = jpj - nn_facty * (jj - 2)
        IF (nn_facty == 3) THEN
          ijjs = ijje - 1
        ELSE
          ijjs = ijje - nn_facty + 1
        END IF
        mjs2_crs(jpj_crs - jj + 1) = ijjs
        mje2_crs(jpj_crs - jj + 1) = ijje
        DO jj = 3, jpjglo_crsm1
          ijje = jpjglo - nn_facty * (jj - 2)
          ijjs = ijje - nn_facty + 1
          IF (ijjs <= nn_facty) ijjs = 2
          mjs2_crs(jpj_crs - jj + 1) = ijjs
          mje2_crs(jpj_crs - jj + 1) = ijje
        END DO
        !$ACC END KERNELS
      CASE DEFAULT
        WRITE(numout, *) 'crs_init. Only jperio = 0, 1, 3, 4, 5, 6 supported'
      END SELECT
    CASE (1)
      CALL profile_psy_data9 % PreStart('crs_dom_def', 'r9', 0, 0)
      WRITE(numout, *) 'crs_init.  Equator-centered bins option not yet available'
      CALL profile_psy_data9 % PostEnd
    END SELECT
    CALL profile_psy_data10 % PreStart('crs_dom_def', 'r10', 0, 0)
    mis2_crs(2) = 1
    mis2_crs(jpiglo_crs) = mie2_crs(jpiglo_crs - 1) + 1
    mie2_crs(2) = nn_factx
    mie2_crs(jpiglo_crs) = jpiglo
    mjs2_crs(1) = 1
    mje2_crs(1) = 1
    mje2_crs(2) = mjs2_crs(3) - 1
    mje2_crs(jpjglo_crs) = jpjglo
    mjs2_crs(2) = 1
    mjs2_crs(jpjglo_crs) = mje2_crs(jpjglo_crs) - nn_facty + 1
    CALL profile_psy_data10 % PostEnd
    IF (.NOT. lk_mpp) THEN
      !$ACC KERNELS
      mis_crs(:) = mis2_crs(:)
      mie_crs(:) = mie2_crs(:)
      mjs_crs(:) = mjs2_crs(:)
      mje_crs(:) = mje2_crs(:)
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      DO jj = 1, nlej_crs
        mjs_crs(jj) = mjs2_crs(mjg_crs(jj)) - njmpp + 1
        mje_crs(jj) = mje2_crs(mjg_crs(jj)) - njmpp + 1
      END DO
      DO ji = 1, nlei_crs
        mis_crs(ji) = mis2_crs(mig_crs(ji)) - nimpp + 1
        mie_crs(ji) = mie2_crs(mig_crs(ji)) - nimpp + 1
      END DO
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data11 % PreStart('crs_dom_def', 'r11', 0, 0)
    nistr = mis_crs(2)
    niend = mis_crs(nlci_crs - 1)
    njstr = mjs_crs(3)
    njend = mjs_crs(nlcj_crs - 1)
    CALL profile_psy_data11 % PostEnd
  END SUBROUTINE crs_dom_def
  SUBROUTINE crs_dom_bat
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jj, jk
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs) :: zmbk
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    !$ACC KERNELS
    mbathy_crs(:, :) = jpkm1
    mbkt_crs(:, :) = 1
    mbku_crs(:, :) = 1
    mbkv_crs(:, :) = 1
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('crs_dom_bat', 'r0', 0, 0)
    DO jj = 1, jpj_crs
      DO ji = 1, jpi_crs
        jk = 0
        DO WHILE (tmask_crs(ji, jj, jk + 1) > 0.)
          jk = jk + 1
        END DO
        mbathy_crs(ji, jj) = FLOAT(jk)
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zmbk(:, :) = 0.0
    zmbk(:, :) = REAL(mbathy_crs(:, :), wp)
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('crs_dom_bat', 'r1', 0, 0)
    CALL crs_lbc_lnk(zmbk, 'T', 1.0)
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    mbathy_crs(:, :) = NINT(zmbk(:, :))
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('crs_dom_bat', 'r2', 0, 0)
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) '    crsini : mbkt is ocean bottom k-index of T-, U-, V- and W-levels '
    IF (lwp) WRITE(numout, *) '    ~~~~~~~~~~~~~'
    CALL profile_psy_data2 % PostEnd
    !$ACC KERNELS
    mbkt_crs(:, :) = MAX(mbathy_crs(:, :), 1)
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj_crsm1
      DO ji = 1, jpi_crsm1
        mbku_crs(ji, jj) = MIN(mbkt_crs(ji + 1, jj), mbkt_crs(ji, jj))
        mbkv_crs(ji, jj) = MIN(mbkt_crs(ji, jj + 1), mbkt_crs(ji, jj))
      END DO
    END DO
    zmbk(:, :) = 1.E0
    zmbk(:, :) = REAL(mbku_crs(:, :), wp)
    !$ACC END KERNELS
    CALL profile_psy_data3 % PreStart('crs_dom_bat', 'r3', 0, 0)
    CALL crs_lbc_lnk(zmbk, 'U', 1.0)
    CALL profile_psy_data3 % PostEnd
    !$ACC KERNELS
    mbku_crs(:, :) = MAX(NINT(zmbk(:, :)), 1)
    zmbk(:, :) = REAL(mbkv_crs(:, :), wp)
    !$ACC END KERNELS
    CALL profile_psy_data4 % PreStart('crs_dom_bat', 'r4', 0, 0)
    CALL crs_lbc_lnk(zmbk, 'V', 1.0)
    CALL profile_psy_data4 % PostEnd
    !$ACC KERNELS
    mbkv_crs(:, :) = MAX(NINT(zmbk(:, :)), 1)
    !$ACC END KERNELS
  END SUBROUTINE crs_dom_bat
END MODULE crsdom
