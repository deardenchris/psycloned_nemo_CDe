MODULE crslbclnk
  USE par_kind, ONLY: wp
  USE dom_oce
  USE crs
  USE lbclnk
  USE in_out_manager
  INTERFACE crs_lbc_lnk
    MODULE PROCEDURE crs_lbc_lnk_3d, crs_lbc_lnk_2d
  END INTERFACE
  PUBLIC :: crs_lbc_lnk
  CONTAINS
  SUBROUTINE crs_lbc_lnk_3d(pt3d1, cd_type1, psgn, cd_mpp, pval)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type1
    REAL(KIND = wp), INTENT(IN) :: psgn
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs, jpk), INTENT(INOUT) :: pt3d1
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    LOGICAL :: ll_grid_crs
    REAL(KIND = wp) :: zval
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('crs_lbc_lnk_3d', 'r0', 0, 0)
    ll_grid_crs = (jpi == jpi_crs)
    IF (PRESENT(pval)) THEN
      zval = pval
    ELSE
      zval = 0._wp
    END IF
    IF (.NOT. ll_grid_crs) CALL dom_grid_crs
    IF (PRESENT(cd_mpp)) THEN
      CALL lbc_lnk('crslbclnk', pt3d1, cd_type1, psgn, cd_mpp, pval = zval)
    ELSE
      CALL lbc_lnk('crslbclnk', pt3d1, cd_type1, psgn, pval = zval)
    END IF
    IF (.NOT. ll_grid_crs) CALL dom_grid_glo
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE crs_lbc_lnk_3d
  SUBROUTINE crs_lbc_lnk_2d(pt2d, cd_type, psgn, cd_mpp, pval)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    REAL(KIND = wp), INTENT(IN) :: psgn
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs), INTENT(INOUT) :: pt2d
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    LOGICAL :: ll_grid_crs
    REAL(KIND = wp) :: zval
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('crs_lbc_lnk_2d', 'r0', 0, 0)
    ll_grid_crs = (jpi == jpi_crs)
    IF (PRESENT(pval)) THEN
      zval = pval
    ELSE
      zval = 0._wp
    END IF
    IF (.NOT. ll_grid_crs) CALL dom_grid_crs
    IF (PRESENT(cd_mpp)) THEN
      CALL lbc_lnk('crslbclnk', pt2d, cd_type, psgn, cd_mpp, pval = zval)
    ELSE
      CALL lbc_lnk('crslbclnk', pt2d, cd_type, psgn, pval = zval)
    END IF
    IF (.NOT. ll_grid_crs) CALL dom_grid_glo
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE crs_lbc_lnk_2d
END MODULE crslbclnk