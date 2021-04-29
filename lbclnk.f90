MODULE lbclnk
  USE oce
  USE dom_oce
  USE in_out_manager
  USE lbcnfd
  IMPLICIT NONE
  PRIVATE
  INTERFACE lbc_lnk
    MODULE PROCEDURE lbc_lnk_2d, lbc_lnk_3d, lbc_lnk_4d
  END INTERFACE
  INTERFACE lbc_lnk_ptr
    MODULE PROCEDURE lbc_lnk_2d_ptr, lbc_lnk_3d_ptr, lbc_lnk_4d_ptr
  END INTERFACE
  INTERFACE lbc_lnk_multi
    MODULE PROCEDURE lbc_lnk_2d_multi, lbc_lnk_3d_multi, lbc_lnk_4d_multi
  END INTERFACE
  INTERFACE lbc_bdy_lnk
    MODULE PROCEDURE lbc_bdy_lnk_2d, lbc_bdy_lnk_3d, lbc_bdy_lnk_4d
  END INTERFACE
  INTERFACE lbc_lnk_icb
    MODULE PROCEDURE lbc_lnk_2d_icb
  END INTERFACE
  PUBLIC :: lbc_lnk
  PUBLIC :: lbc_lnk_multi
  PUBLIC :: lbc_bdy_lnk
  PUBLIC :: lbc_lnk_icb
  CONTAINS
  SUBROUTINE lbc_lnk_2d(cdname, ptab, cd_nat, psgn, cd_mpp, pval)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    REAL(KIND = wp) :: zland
    LOGICAL :: ll_nfd
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('lbc_lnk_2d', 'r0', 0, 0)
    ipk = 1
    ipl = 1
    ipf = 1
    ll_nfd = jperio == 3 .OR. jperio == 4 .OR. jperio == 5 .OR. jperio == 6
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (.NOT. PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(1, :) = ptab(jpim1, :)
          ptab(jpi, :) = ptab(2, :)
          !$ACC END KERNELS
        ELSE
          CALL profile_psy_data1 % PreStart('lbc_lnk_2d', 'r1', 0, 0)
          IF (.NOT. cd_nat == 'F') ptab(1, :) = zland
          CALL profile_psy_data1 % PostEnd
          !$ACC KERNELS
          ptab(jpi, :) = zland
          !$ACC END KERNELS
        END IF
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(:, 1) = ptab(:, jpjm1)
          ptab(:, jpj) = ptab(:, 2)
          !$ACC END KERNELS
        ELSE IF (ll_nfd) THEN
          CALL profile_psy_data2 % PreStart('lbc_lnk_2d', 'r2', 0, 0)
          IF (.NOT. cd_nat == 'F') ptab(:, 1) = zland
          CALL lbc_nfd(ptab, cd_nat, psgn)
          CALL profile_psy_data2 % PostEnd
        ELSE
          CALL profile_psy_data3 % PreStart('lbc_lnk_2d', 'r3', 0, 0)
          IF (.NOT. cd_nat == 'F') ptab(:, 1) = zland
          CALL profile_psy_data3 % PostEnd
          !$ACC KERNELS
          ptab(:, jpj) = zland
          !$ACC END KERNELS
        END IF
      END DO
    END IF
  END SUBROUTINE lbc_lnk_2d
  SUBROUTINE lbc_lnk_2d_ptr(cdname, ptab, cd_nat, psgn, kfld, cd_mpp, pval)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kfld
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    TYPE(PTR_2D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    REAL(KIND = wp) :: zland
    LOGICAL :: ll_nfd
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('lbc_lnk_2d_ptr', 'r0', 0, 0)
    ipk = 1
    ipl = 1
    ipf = kfld
    ll_nfd = jperio == 3 .OR. jperio == 4 .OR. jperio == 5 .OR. jperio == 6
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (.NOT. PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(jf) % pt2d(1, :) = ptab(jf) % pt2d(jpim1, :)
          ptab(jf) % pt2d(jpi, :) = ptab(jf) % pt2d(2, :)
          !$ACC END KERNELS
        ELSE
          CALL profile_psy_data1 % PreStart('lbc_lnk_2d_ptr', 'r1', 0, 0)
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt2d(1, :) = zland
          CALL profile_psy_data1 % PostEnd
          !$ACC KERNELS
          ptab(jf) % pt2d(jpi, :) = zland
          !$ACC END KERNELS
        END IF
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(jf) % pt2d(:, 1) = ptab(jf) % pt2d(:, jpjm1)
          ptab(jf) % pt2d(:, jpj) = ptab(jf) % pt2d(:, 2)
          !$ACC END KERNELS
        ELSE IF (ll_nfd) THEN
          CALL profile_psy_data2 % PreStart('lbc_lnk_2d_ptr', 'r2', 0, 0)
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt2d(:, 1) = zland
          CALL lbc_nfd(ptab, cd_nat(:), psgn(:), ipf)
          CALL profile_psy_data2 % PostEnd
        ELSE
          CALL profile_psy_data3 % PreStart('lbc_lnk_2d_ptr', 'r3', 0, 0)
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt2d(:, 1) = zland
          CALL profile_psy_data3 % PostEnd
          !$ACC KERNELS
          ptab(jf) % pt2d(:, jpj) = zland
          !$ACC END KERNELS
        END IF
      END DO
    END IF
  END SUBROUTINE lbc_lnk_2d_ptr
  SUBROUTINE lbc_lnk_3d(cdname, ptab, cd_nat, psgn, cd_mpp, pval)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    REAL(KIND = wp) :: zland
    LOGICAL :: ll_nfd
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('lbc_lnk_3d', 'r0', 0, 0)
    ipk = SIZE(ptab, 3)
    ipl = 1
    ipf = 1
    ll_nfd = jperio == 3 .OR. jperio == 4 .OR. jperio == 5 .OR. jperio == 6
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (.NOT. PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(1, :, :) = ptab(jpim1, :, :)
          ptab(jpi, :, :) = ptab(2, :, :)
          !$ACC END KERNELS
        ELSE
          CALL profile_psy_data1 % PreStart('lbc_lnk_3d', 'r1', 0, 0)
          IF (.NOT. cd_nat == 'F') ptab(1, :, :) = zland
          CALL profile_psy_data1 % PostEnd
          !$ACC KERNELS
          ptab(jpi, :, :) = zland
          !$ACC END KERNELS
        END IF
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(:, 1, :) = ptab(:, jpjm1, :)
          ptab(:, jpj, :) = ptab(:, 2, :)
          !$ACC END KERNELS
        ELSE IF (ll_nfd) THEN
          CALL profile_psy_data2 % PreStart('lbc_lnk_3d', 'r2', 0, 0)
          IF (.NOT. cd_nat == 'F') ptab(:, 1, :) = zland
          CALL lbc_nfd(ptab, cd_nat, psgn)
          CALL profile_psy_data2 % PostEnd
        ELSE
          CALL profile_psy_data3 % PreStart('lbc_lnk_3d', 'r3', 0, 0)
          IF (.NOT. cd_nat == 'F') ptab(:, 1, :) = zland
          CALL profile_psy_data3 % PostEnd
          !$ACC KERNELS
          ptab(:, jpj, :) = zland
          !$ACC END KERNELS
        END IF
      END DO
    END IF
  END SUBROUTINE lbc_lnk_3d
  SUBROUTINE lbc_lnk_3d_ptr(cdname, ptab, cd_nat, psgn, kfld, cd_mpp, pval)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kfld
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    TYPE(PTR_3D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    REAL(KIND = wp) :: zland
    LOGICAL :: ll_nfd
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('lbc_lnk_3d_ptr', 'r0', 0, 0)
    ipk = SIZE(ptab(1) % pt3d, 3)
    ipl = 1
    ipf = kfld
    ll_nfd = jperio == 3 .OR. jperio == 4 .OR. jperio == 5 .OR. jperio == 6
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (.NOT. PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(jf) % pt3d(1, :, :) = ptab(jf) % pt3d(jpim1, :, :)
          ptab(jf) % pt3d(jpi, :, :) = ptab(jf) % pt3d(2, :, :)
          !$ACC END KERNELS
        ELSE
          CALL profile_psy_data1 % PreStart('lbc_lnk_3d_ptr', 'r1', 0, 0)
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt3d(1, :, :) = zland
          CALL profile_psy_data1 % PostEnd
          !$ACC KERNELS
          ptab(jf) % pt3d(jpi, :, :) = zland
          !$ACC END KERNELS
        END IF
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(jf) % pt3d(:, 1, :) = ptab(jf) % pt3d(:, jpjm1, :)
          ptab(jf) % pt3d(:, jpj, :) = ptab(jf) % pt3d(:, 2, :)
          !$ACC END KERNELS
        ELSE IF (ll_nfd) THEN
          CALL profile_psy_data2 % PreStart('lbc_lnk_3d_ptr', 'r2', 0, 0)
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt3d(:, 1, :) = zland
          CALL lbc_nfd(ptab, cd_nat(:), psgn(:), ipf)
          CALL profile_psy_data2 % PostEnd
        ELSE
          CALL profile_psy_data3 % PreStart('lbc_lnk_3d_ptr', 'r3', 0, 0)
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt3d(:, 1, :) = zland
          CALL profile_psy_data3 % PostEnd
          !$ACC KERNELS
          ptab(jf) % pt3d(:, jpj, :) = zland
          !$ACC END KERNELS
        END IF
      END DO
    END IF
  END SUBROUTINE lbc_lnk_3d_ptr
  SUBROUTINE lbc_lnk_4d(cdname, ptab, cd_nat, psgn, cd_mpp, pval)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    REAL(KIND = wp) :: zland
    LOGICAL :: ll_nfd
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('lbc_lnk_4d', 'r0', 0, 0)
    ipk = SIZE(ptab, 3)
    ipl = SIZE(ptab, 4)
    ipf = 1
    ll_nfd = jperio == 3 .OR. jperio == 4 .OR. jperio == 5 .OR. jperio == 6
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (.NOT. PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(1, :, :, :) = ptab(jpim1, :, :, :)
          ptab(jpi, :, :, :) = ptab(2, :, :, :)
          !$ACC END KERNELS
        ELSE
          CALL profile_psy_data1 % PreStart('lbc_lnk_4d', 'r1', 0, 0)
          IF (.NOT. cd_nat == 'F') ptab(1, :, :, :) = zland
          CALL profile_psy_data1 % PostEnd
          !$ACC KERNELS
          ptab(jpi, :, :, :) = zland
          !$ACC END KERNELS
        END IF
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(:, 1, :, :) = ptab(:, jpjm1, :, :)
          ptab(:, jpj, :, :) = ptab(:, 2, :, :)
          !$ACC END KERNELS
        ELSE IF (ll_nfd) THEN
          CALL profile_psy_data2 % PreStart('lbc_lnk_4d', 'r2', 0, 0)
          IF (.NOT. cd_nat == 'F') ptab(:, 1, :, :) = zland
          CALL lbc_nfd(ptab, cd_nat, psgn)
          CALL profile_psy_data2 % PostEnd
        ELSE
          CALL profile_psy_data3 % PreStart('lbc_lnk_4d', 'r3', 0, 0)
          IF (.NOT. cd_nat == 'F') ptab(:, 1, :, :) = zland
          CALL profile_psy_data3 % PostEnd
          !$ACC KERNELS
          ptab(:, jpj, :, :) = zland
          !$ACC END KERNELS
        END IF
      END DO
    END IF
  END SUBROUTINE lbc_lnk_4d
  SUBROUTINE lbc_lnk_4d_ptr(cdname, ptab, cd_nat, psgn, kfld, cd_mpp, pval)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kfld
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    TYPE(PTR_4D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    REAL(KIND = wp) :: zland
    LOGICAL :: ll_nfd
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('lbc_lnk_4d_ptr', 'r0', 0, 0)
    ipk = SIZE(ptab(1) % pt4d, 3)
    ipl = SIZE(ptab(1) % pt4d, 4)
    ipf = kfld
    ll_nfd = jperio == 3 .OR. jperio == 4 .OR. jperio == 5 .OR. jperio == 6
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (.NOT. PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(jf) % pt4d(1, :, :, :) = ptab(jf) % pt4d(jpim1, :, :, :)
          ptab(jf) % pt4d(jpi, :, :, :) = ptab(jf) % pt4d(2, :, :, :)
          !$ACC END KERNELS
        ELSE
          CALL profile_psy_data1 % PreStart('lbc_lnk_4d_ptr', 'r1', 0, 0)
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt4d(1, :, :, :) = zland
          CALL profile_psy_data1 % PostEnd
          !$ACC KERNELS
          ptab(jf) % pt4d(jpi, :, :, :) = zland
          !$ACC END KERNELS
        END IF
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(jf) % pt4d(:, 1, :, :) = ptab(jf) % pt4d(:, jpjm1, :, :)
          ptab(jf) % pt4d(:, jpj, :, :) = ptab(jf) % pt4d(:, 2, :, :)
          !$ACC END KERNELS
        ELSE IF (ll_nfd) THEN
          CALL profile_psy_data2 % PreStart('lbc_lnk_4d_ptr', 'r2', 0, 0)
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt4d(:, 1, :, :) = zland
          CALL lbc_nfd(ptab, cd_nat(:), psgn(:), ipf)
          CALL profile_psy_data2 % PostEnd
        ELSE
          CALL profile_psy_data3 % PreStart('lbc_lnk_4d_ptr', 'r3', 0, 0)
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt4d(:, 1, :, :) = zland
          CALL profile_psy_data3 % PostEnd
          !$ACC KERNELS
          ptab(jf) % pt4d(:, jpj, :, :) = zland
          !$ACC END KERNELS
        END IF
      END DO
    END IF
  END SUBROUTINE lbc_lnk_4d_ptr
  SUBROUTINE lbc_bdy_lnk_4d(cdname, pt4d, cd_type, psgn, ib_bdy)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), DIMENSION(:, :, :, :), INTENT(INOUT) :: pt4d
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, INTENT(IN) :: ib_bdy
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('lbc_bdy_lnk_4d', 'r0', 0, 0)
    CALL lbc_lnk_4d(cdname, pt4d, cd_type, psgn)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE lbc_bdy_lnk_4d
  SUBROUTINE lbc_bdy_lnk_3d(cdname, pt3d, cd_type, psgn, ib_bdy)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pt3d
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, INTENT(IN) :: ib_bdy
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('lbc_bdy_lnk_3d', 'r0', 0, 0)
    CALL lbc_lnk_3d(cdname, pt3d, cd_type, psgn)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE lbc_bdy_lnk_3d
  SUBROUTINE lbc_bdy_lnk_2d(cdname, pt2d, cd_type, psgn, ib_bdy)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: pt2d
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, INTENT(IN) :: ib_bdy
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('lbc_bdy_lnk_2d', 'r0', 0, 0)
    CALL lbc_lnk_2d(cdname, pt2d, cd_type, psgn)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE lbc_bdy_lnk_2d
  SUBROUTINE lbc_lnk_2d_icb(cdname, pt2d, cd_type, psgn, ki, kj)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: pt2d
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, INTENT(IN) :: ki, kj
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('lbc_lnk_2d_icb', 'r0', 0, 0)
    CALL lbc_lnk_2d(cdname, pt2d, cd_type, psgn)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE lbc_lnk_2d_icb
  SUBROUTINE lbc_lnk_2d_multi(cdname, pt1, cdna1, psgn1, pt2, cdna2, psgn2, pt3, cdna3, psgn3, pt4, cdna4, psgn4, pt5, cdna5, psgn5, pt6, cdna6, psgn6, pt7, cdna7, psgn7, pt8, cdna8, psgn8, pt9, cdna9, psgn9, cd_mpp, pval)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
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
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('lbc_lnk_2d_multi', 'r0', 0, 0)
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
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE lbc_lnk_2d_multi
  SUBROUTINE load_ptr_2d(ptab, cdna, psgn, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :), TARGET, INTENT(INOUT) :: ptab
    CHARACTER(LEN = 1), INTENT(IN) :: cdna
    REAL(KIND = wp), INTENT(IN) :: psgn
    TYPE(PTR_2D), DIMENSION(:), INTENT(INOUT) :: ptab_ptr
    CHARACTER(LEN = 1), DIMENSION(:), INTENT(INOUT) :: cdna_ptr
    REAL(KIND = wp), DIMENSION(:), INTENT(INOUT) :: psgn_ptr
    INTEGER, INTENT(INOUT) :: kfld
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('load_ptr_2d', 'r0', 0, 0)
    kfld = kfld + 1
    ptab_ptr(kfld) % pt2d => ptab
    cdna_ptr(kfld) = cdna
    psgn_ptr(kfld) = psgn
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE load_ptr_2d
  SUBROUTINE lbc_lnk_3d_multi(cdname, pt1, cdna1, psgn1, pt2, cdna2, psgn2, pt3, cdna3, psgn3, pt4, cdna4, psgn4, pt5, cdna5, psgn5, pt6, cdna6, psgn6, pt7, cdna7, psgn7, pt8, cdna8, psgn8, pt9, cdna9, psgn9, cd_mpp, pval)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
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
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('lbc_lnk_3d_multi', 'r0', 0, 0)
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
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE lbc_lnk_3d_multi
  SUBROUTINE load_ptr_3d(ptab, cdna, psgn, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :), TARGET, INTENT(INOUT) :: ptab
    CHARACTER(LEN = 1), INTENT(IN) :: cdna
    REAL(KIND = wp), INTENT(IN) :: psgn
    TYPE(PTR_3D), DIMENSION(:), INTENT(INOUT) :: ptab_ptr
    CHARACTER(LEN = 1), DIMENSION(:), INTENT(INOUT) :: cdna_ptr
    REAL(KIND = wp), DIMENSION(:), INTENT(INOUT) :: psgn_ptr
    INTEGER, INTENT(INOUT) :: kfld
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('load_ptr_3d', 'r0', 0, 0)
    kfld = kfld + 1
    ptab_ptr(kfld) % pt3d => ptab
    cdna_ptr(kfld) = cdna
    psgn_ptr(kfld) = psgn
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE load_ptr_3d
  SUBROUTINE lbc_lnk_4d_multi(cdname, pt1, cdna1, psgn1, pt2, cdna2, psgn2, pt3, cdna3, psgn3, pt4, cdna4, psgn4, pt5, cdna5, psgn5, pt6, cdna6, psgn6, pt7, cdna7, psgn7, pt8, cdna8, psgn8, pt9, cdna9, psgn9, cd_mpp, pval)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
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
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('lbc_lnk_4d_multi', 'r0', 0, 0)
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
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE lbc_lnk_4d_multi
  SUBROUTINE load_ptr_4d(ptab, cdna, psgn, ptab_ptr, cdna_ptr, psgn_ptr, kfld)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :, :), TARGET, INTENT(INOUT) :: ptab
    CHARACTER(LEN = 1), INTENT(IN) :: cdna
    REAL(KIND = wp), INTENT(IN) :: psgn
    TYPE(PTR_4D), DIMENSION(:), INTENT(INOUT) :: ptab_ptr
    CHARACTER(LEN = 1), DIMENSION(:), INTENT(INOUT) :: cdna_ptr
    REAL(KIND = wp), DIMENSION(:), INTENT(INOUT) :: psgn_ptr
    INTEGER, INTENT(INOUT) :: kfld
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('load_ptr_4d', 'r0', 0, 0)
    kfld = kfld + 1
    ptab_ptr(kfld) % pt4d => ptab
    cdna_ptr(kfld) = cdna
    psgn_ptr(kfld) = psgn
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE load_ptr_4d
END MODULE lbclnk