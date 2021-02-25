MODULE lbclnk
  USE par_oce
  USE lib_mpp
  USE lbcnfd
  INTERFACE lbc_lnk
    MODULE PROCEDURE mpp_lnk_2d, mpp_lnk_3d, mpp_lnk_4d
  END INTERFACE
  INTERFACE lbc_lnk_ptr
    MODULE PROCEDURE mpp_lnk_2d_ptr, mpp_lnk_3d_ptr, mpp_lnk_4d_ptr
  END INTERFACE
  INTERFACE lbc_lnk_multi
    MODULE PROCEDURE lbc_lnk_2d_multi, lbc_lnk_3d_multi, lbc_lnk_4d_multi
  END INTERFACE
  INTERFACE lbc_bdy_lnk
    MODULE PROCEDURE mpp_lnk_bdy_2d, mpp_lnk_bdy_3d, mpp_lnk_bdy_4d
  END INTERFACE
  INTERFACE lbc_lnk_icb
    MODULE PROCEDURE mpp_lnk_2d_icb
  END INTERFACE
  PUBLIC :: lbc_lnk
  PUBLIC :: lbc_lnk_multi
  PUBLIC :: lbc_bdy_lnk
  PUBLIC :: lbc_lnk_icb
  CONTAINS
  SUBROUTINE lbc_lnk_2d_multi(cdname, pt1, cdna1, psgn1, pt2, cdna2, psgn2, pt3, cdna3, psgn3, pt4, cdna4, psgn4, pt5, cdna5, &
&psgn5, pt6, cdna6, psgn6, pt7, cdna7, psgn7, pt8, cdna8, psgn8, pt9, cdna9, psgn9, cd_mpp, pval)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), DIMENSION(:, :), TARGET, INTENT(INOUT) :: pt1
    REAL(KIND = wp), DIMENSION(:, :), OPTIONAL, TARGET, INTENT(INOUT) :: pt2, pt3, pt4, pt5, pt6, pt7, pt8, pt9
    CHARACTER(LEN = 1), INTENT(IN) :: cdna1
    CHARACTER(LEN = 1), OPTIONAL, INTENT(IN) :: cdna2, cdna3, cdna4, cdna5, cdna6, cdna7, cdna8, cdna9
    REAL(KIND = wp), INTENT(IN) :: psgn1
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: psgn2, psgn3, psgn4, psgn5, psgn6, psgn7, psgn8, psgn9
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: kfld
    TYPE(PTR_2D), DIMENSION(9) :: ptab_ptr
    CHARACTER(LEN = 1), DIMENSION(9) :: cdna_ptr
    REAL(KIND = wp), DIMENSION(9) :: psgn_ptr
    kfld = 0
    CALL load_ptr_2d(pt1, cdna1, psgn1, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn2)) CALL load_ptr_2d(pt2, cdna2, psgn2, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn3)) CALL load_ptr_2d(pt3, cdna3, psgn3, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn4)) CALL load_ptr_2d(pt4, cdna4, psgn4, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn5)) CALL load_ptr_2d(pt5, cdna5, psgn5, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn6)) CALL load_ptr_2d(pt6, cdna6, psgn6, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn7)) CALL load_ptr_2d(pt7, cdna7, psgn7, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn8)) CALL load_ptr_2d(pt8, cdna8, psgn8, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn9)) CALL load_ptr_2d(pt9, cdna9, psgn9, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    CALL lbc_lnk_ptr(cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, cd_mpp, pval)
  END SUBROUTINE lbc_lnk_2d_multi
  SUBROUTINE load_ptr_2d(ptab, cdna, psgn, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    REAL(KIND = wp), DIMENSION(:, :), TARGET, INTENT(INOUT) :: ptab
    CHARACTER(LEN = 1), INTENT(IN) :: cdna
    REAL(KIND = wp), INTENT(IN) :: psgn
    TYPE(PTR_2D), DIMENSION(:), INTENT(INOUT) :: ptab_ptr
    CHARACTER(LEN = 1), DIMENSION(:), INTENT(INOUT) :: cdna_ptr
    REAL(KIND = wp), DIMENSION(:), INTENT(INOUT) :: psgn_ptr
    INTEGER, INTENT(INOUT) :: kfld
    kfld = kfld + 1
    ptab_ptr(kfld) % pt2d => ptab
    cdna_ptr(kfld) = cdna
    psgn_ptr(kfld) = psgn
  END SUBROUTINE load_ptr_2d
  SUBROUTINE lbc_lnk_3d_multi(cdname, pt1, cdna1, psgn1, pt2, cdna2, psgn2, pt3, cdna3, psgn3, pt4, cdna4, psgn4, pt5, cdna5, &
&psgn5, pt6, cdna6, psgn6, pt7, cdna7, psgn7, pt8, cdna8, psgn8, pt9, cdna9, psgn9, cd_mpp, pval)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), DIMENSION(:, :, :), TARGET, INTENT(INOUT) :: pt1
    REAL(KIND = wp), DIMENSION(:, :, :), OPTIONAL, TARGET, INTENT(INOUT) :: pt2, pt3, pt4, pt5, pt6, pt7, pt8, pt9
    CHARACTER(LEN = 1), INTENT(IN) :: cdna1
    CHARACTER(LEN = 1), OPTIONAL, INTENT(IN) :: cdna2, cdna3, cdna4, cdna5, cdna6, cdna7, cdna8, cdna9
    REAL(KIND = wp), INTENT(IN) :: psgn1
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: psgn2, psgn3, psgn4, psgn5, psgn6, psgn7, psgn8, psgn9
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: kfld
    TYPE(PTR_3D), DIMENSION(9) :: ptab_ptr
    CHARACTER(LEN = 1), DIMENSION(9) :: cdna_ptr
    REAL(KIND = wp), DIMENSION(9) :: psgn_ptr
    kfld = 0
    CALL load_ptr_3d(pt1, cdna1, psgn1, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn2)) CALL load_ptr_3d(pt2, cdna2, psgn2, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn3)) CALL load_ptr_3d(pt3, cdna3, psgn3, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn4)) CALL load_ptr_3d(pt4, cdna4, psgn4, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn5)) CALL load_ptr_3d(pt5, cdna5, psgn5, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn6)) CALL load_ptr_3d(pt6, cdna6, psgn6, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn7)) CALL load_ptr_3d(pt7, cdna7, psgn7, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn8)) CALL load_ptr_3d(pt8, cdna8, psgn8, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn9)) CALL load_ptr_3d(pt9, cdna9, psgn9, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    CALL lbc_lnk_ptr(cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, cd_mpp, pval)
  END SUBROUTINE lbc_lnk_3d_multi
  SUBROUTINE load_ptr_3d(ptab, cdna, psgn, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    REAL(KIND = wp), DIMENSION(:, :, :), TARGET, INTENT(INOUT) :: ptab
    CHARACTER(LEN = 1), INTENT(IN) :: cdna
    REAL(KIND = wp), INTENT(IN) :: psgn
    TYPE(PTR_3D), DIMENSION(:), INTENT(INOUT) :: ptab_ptr
    CHARACTER(LEN = 1), DIMENSION(:), INTENT(INOUT) :: cdna_ptr
    REAL(KIND = wp), DIMENSION(:), INTENT(INOUT) :: psgn_ptr
    INTEGER, INTENT(INOUT) :: kfld
    kfld = kfld + 1
    ptab_ptr(kfld) % pt3d => ptab
    cdna_ptr(kfld) = cdna
    psgn_ptr(kfld) = psgn
  END SUBROUTINE load_ptr_3d
  SUBROUTINE lbc_lnk_4d_multi(cdname, pt1, cdna1, psgn1, pt2, cdna2, psgn2, pt3, cdna3, psgn3, pt4, cdna4, psgn4, pt5, cdna5, &
&psgn5, pt6, cdna6, psgn6, pt7, cdna7, psgn7, pt8, cdna8, psgn8, pt9, cdna9, psgn9, cd_mpp, pval)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), DIMENSION(:, :, :, :), TARGET, INTENT(INOUT) :: pt1
    REAL(KIND = wp), DIMENSION(:, :, :, :), OPTIONAL, TARGET, INTENT(INOUT) :: pt2, pt3, pt4, pt5, pt6, pt7, pt8, pt9
    CHARACTER(LEN = 1), INTENT(IN) :: cdna1
    CHARACTER(LEN = 1), OPTIONAL, INTENT(IN) :: cdna2, cdna3, cdna4, cdna5, cdna6, cdna7, cdna8, cdna9
    REAL(KIND = wp), INTENT(IN) :: psgn1
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: psgn2, psgn3, psgn4, psgn5, psgn6, psgn7, psgn8, psgn9
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: kfld
    TYPE(PTR_4D), DIMENSION(9) :: ptab_ptr
    CHARACTER(LEN = 1), DIMENSION(9) :: cdna_ptr
    REAL(KIND = wp), DIMENSION(9) :: psgn_ptr
    kfld = 0
    CALL load_ptr_4d(pt1, cdna1, psgn1, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn2)) CALL load_ptr_4d(pt2, cdna2, psgn2, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn3)) CALL load_ptr_4d(pt3, cdna3, psgn3, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn4)) CALL load_ptr_4d(pt4, cdna4, psgn4, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn5)) CALL load_ptr_4d(pt5, cdna5, psgn5, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn6)) CALL load_ptr_4d(pt6, cdna6, psgn6, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn7)) CALL load_ptr_4d(pt7, cdna7, psgn7, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn8)) CALL load_ptr_4d(pt8, cdna8, psgn8, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    IF (PRESENT(psgn9)) CALL load_ptr_4d(pt9, cdna9, psgn9, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    CALL lbc_lnk_ptr(cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, cd_mpp, pval)
  END SUBROUTINE lbc_lnk_4d_multi
  SUBROUTINE load_ptr_4d(ptab, cdna, psgn, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    REAL(KIND = wp), DIMENSION(:, :, :, :), TARGET, INTENT(INOUT) :: ptab
    CHARACTER(LEN = 1), INTENT(IN) :: cdna
    REAL(KIND = wp), INTENT(IN) :: psgn
    TYPE(PTR_4D), DIMENSION(:), INTENT(INOUT) :: ptab_ptr
    CHARACTER(LEN = 1), DIMENSION(:), INTENT(INOUT) :: cdna_ptr
    REAL(KIND = wp), DIMENSION(:), INTENT(INOUT) :: psgn_ptr
    INTEGER, INTENT(INOUT) :: kfld
    kfld = kfld + 1
    ptab_ptr(kfld) % pt4d => ptab
    cdna_ptr(kfld) = cdna
    psgn_ptr(kfld) = psgn
  END SUBROUTINE load_ptr_4d
END MODULE lbclnk