MODULE lbcnfd
  USE dom_oce
  USE in_out_manager
  IMPLICIT NONE
  PRIVATE
  INTERFACE lbc_nfd
    MODULE PROCEDURE lbc_nfd_2d, lbc_nfd_3d, lbc_nfd_4d
    MODULE PROCEDURE lbc_nfd_2d_ptr, lbc_nfd_3d_ptr, lbc_nfd_4d_ptr
    MODULE PROCEDURE lbc_nfd_2d_ext
  END INTERFACE
  INTERFACE lbc_nfd_nogather
    MODULE PROCEDURE lbc_nfd_nogather_2d, lbc_nfd_nogather_3d
    MODULE PROCEDURE lbc_nfd_nogather_4d
    MODULE PROCEDURE lbc_nfd_nogather_2d_ptr, lbc_nfd_nogather_3d_ptr
  END INTERFACE
  TYPE, PUBLIC :: PTR_2D
    REAL(KIND = wp), DIMENSION(:, :), POINTER :: pt2d
  END TYPE PTR_2D
  TYPE, PUBLIC :: PTR_3D
    REAL(KIND = wp), DIMENSION(:, :, :), POINTER :: pt3d
  END TYPE PTR_3D
  TYPE, PUBLIC :: PTR_4D
    REAL(KIND = wp), DIMENSION(:, :, :, :), POINTER :: pt4d
  END TYPE PTR_4D
  PUBLIC :: lbc_nfd
  PUBLIC :: lbc_nfd_nogather
  INTEGER, PUBLIC, PARAMETER :: jpmaxngh = 3
  INTEGER, PUBLIC :: nsndto, nfsloop, nfeloop
  INTEGER, PUBLIC, DIMENSION(jpmaxngh) :: isendto
  CONTAINS
  SUBROUTINE lbc_nfd_2d(ptab, cd_nat, psgn)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ipjm1
    !$ACC KERNELS
    ipk = 1
    ipl = 1
    ipf = 1
    SELECT CASE (jpni)
    CASE (1)
      ipj = nlcj
    CASE DEFAULT
      ipj = 4
    END SELECT
    ipjm1 = ipj - 1
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipj) = psgn * ptab(ijt, ipj - 2)
          END DO
          ptab(1, ipj) = psgn * ptab(3, ipj - 2)
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipjm1) = psgn * ptab(ijt, ipjm1)
          END DO
        CASE ('U')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipj) = psgn * ptab(iju, ipj - 2)
          END DO
          ptab(1, ipj) = psgn * ptab(2, ipj - 2)
          ptab(jpiglo, ipj) = psgn * ptab(jpiglo - 1, ipj - 2)
          DO ji = jpiglo / 2, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipjm1) = psgn * ptab(iju, ipjm1)
          END DO
        CASE ('V')
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipj - 1) = psgn * ptab(ijt, ipj - 2)
            ptab(ji, ipj) = psgn * ptab(ijt, ipj - 3)
          END DO
          ptab(1, ipj) = psgn * ptab(3, ipj - 3)
        CASE ('F')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipj - 1) = psgn * ptab(iju, ipj - 2)
            ptab(ji, ipj) = psgn * ptab(iju, ipj - 3)
          END DO
          ptab(1, ipj) = psgn * ptab(2, ipj - 3)
          ptab(jpiglo, ipj) = psgn * ptab(jpiglo - 1, ipj - 3)
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipj) = psgn * ptab(ijt, ipj - 1)
          END DO
        CASE ('U')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipj) = psgn * ptab(iju, ipj - 1)
          END DO
          ptab(jpiglo, ipj) = psgn * ptab(jpiglo - 2, ipj - 1)
        CASE ('V')
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipj) = psgn * ptab(ijt, ipj - 2)
          END DO
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipjm1) = psgn * ptab(ijt, ipjm1)
          END DO
        CASE ('F')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipj) = psgn * ptab(iju, ipj - 2)
          END DO
          ptab(jpiglo, ipj) = psgn * ptab(jpiglo - 2, ipj - 2)
          DO ji = jpiglo / 2 + 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipjm1) = psgn * ptab(iju, ipjm1)
          END DO
        END SELECT
      CASE DEFAULT
        SELECT CASE (cd_nat)
        CASE ('T', 'U', 'V', 'W')
          ptab(:, 1) = 0._wp
          ptab(:, ipj) = 0._wp
        CASE ('F')
          ptab(:, ipj) = 0._wp
        END SELECT
      END SELECT
    END DO
    !$ACC END KERNELS
  END SUBROUTINE lbc_nfd_2d
  SUBROUTINE lbc_nfd_2d_ptr(ptab, cd_nat, psgn, kfld)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kfld
    TYPE(PTR_2D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ipjm1
!    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
!    CALL profile_psy_data0 % PreStart('lbc_nfd_2d_ptr', 'r0', 0, 0)
    ipk = 1
    ipl = 1
    ipf = kfld
    SELECT CASE (jpni)
    CASE (1)
      ipj = nlcj
    CASE DEFAULT
      ipj = 4
    END SELECT
    ipjm1 = ipj - 1
    CALL profile_psy_data0 % PostEnd
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'W')
          !$ACC KERNELS
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt2d(ji, ipj) = psgn(jf) * ptab(jf) % pt2d(ijt, ipj - 2)
          END DO
          ptab(jf) % pt2d(1, ipj) = psgn(jf) * ptab(jf) % pt2d(3, ipj - 2)
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt2d(ji, ipjm1) = psgn(jf) * ptab(jf) % pt2d(ijt, ipjm1)
          END DO
          !$ACC END KERNELS
        CASE ('U')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt2d(ji, ipj) = psgn(jf) * ptab(jf) % pt2d(iju, ipj - 2)
          END DO
          ptab(jf) % pt2d(1, ipj) = psgn(jf) * ptab(jf) % pt2d(2, ipj - 2)
          ptab(jf) % pt2d(jpiglo, ipj) = psgn(jf) * ptab(jf) % pt2d(jpiglo - 1, ipj - 2)
          DO ji = jpiglo / 2, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt2d(ji, ipjm1) = psgn(jf) * ptab(jf) % pt2d(iju, ipjm1)
          END DO
          !$ACC END KERNELS
        CASE ('V')
          !$ACC KERNELS ! CDe      
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt2d(ji, ipj - 1) = psgn(jf) * ptab(jf) % pt2d(ijt, ipj - 2)
            ptab(jf) % pt2d(ji, ipj) = psgn(jf) * ptab(jf) % pt2d(ijt, ipj - 3)
          END DO
          ptab(jf) % pt2d(1, ipj) = psgn(jf) * ptab(jf) % pt2d(3, ipj - 3)
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS ! CDe      
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt2d(ji, ipj - 1) = psgn(jf) * ptab(jf) % pt2d(iju, ipj - 2)
            ptab(jf) % pt2d(ji, ipj) = psgn(jf) * ptab(jf) % pt2d(iju, ipj - 3)
          END DO
          ptab(jf) % pt2d(1, ipj) = psgn(jf) * ptab(jf) % pt2d(2, ipj - 3)
          ptab(jf) % pt2d(jpiglo, ipj) = psgn(jf) * ptab(jf) % pt2d(jpiglo - 1, ipj - 3)
          !$ACC END KERNELS
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'W')
          !$ACC KERNELS
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(jf) % pt2d(ji, ipj) = psgn(jf) * ptab(jf) % pt2d(ijt, ipj - 1)
          END DO
          !$ACC END KERNELS
        CASE ('U')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(jf) % pt2d(ji, ipj) = psgn(jf) * ptab(jf) % pt2d(iju, ipj - 1)
          END DO
          ptab(jf) % pt2d(jpiglo, ipj) = psgn(jf) * ptab(jf) % pt2d(jpiglo - 2, ipj - 1)
          !$ACC END KERNELS
        CASE ('V')
          !$ACC KERNELS
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(jf) % pt2d(ji, ipj) = psgn(jf) * ptab(jf) % pt2d(ijt, ipj - 2)
          END DO
          !$ACC LOOP INDEPENDENT ! CDe forces parallelism
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(jf) % pt2d(ji, ipjm1) = psgn(jf) * ptab(jf) % pt2d(ijt, ipjm1)
          END DO
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(jf) % pt2d(ji, ipj) = psgn(jf) * ptab(jf) % pt2d(iju, ipj - 2)
          END DO
          ptab(jf) % pt2d(jpiglo, ipj) = psgn(jf) * ptab(jf) % pt2d(jpiglo - 2, ipj - 2)
          DO ji = jpiglo / 2 + 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(jf) % pt2d(ji, ipjm1) = psgn(jf) * ptab(jf) % pt2d(iju, ipjm1)
          END DO
          !$ACC END KERNELS
        END SELECT
      CASE DEFAULT
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'U', 'V', 'W')
          !$ACC KERNELS
          ptab(jf) % pt2d(:, 1) = 0._wp
          ptab(jf) % pt2d(:, ipj) = 0._wp
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          ptab(jf) % pt2d(:, ipj) = 0._wp
          !$ACC END KERNELS
        END SELECT
      END SELECT
    END DO
!    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE lbc_nfd_2d_ptr
  SUBROUTINE lbc_nfd_2d_ext(ptab, cd_nat, psgn, kextj)
    INTEGER, INTENT(IN) :: kextj
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, 1 - kextj :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ipjm1
    !$ACC KERNELS
    ipk = 1
    ipl = 1
    ipf = 1
    SELECT CASE (jpni)
    CASE (1)
      ipj = nlcj
    CASE DEFAULT
      ipj = 4
    END SELECT
    ipjm1 = ipj - 1
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          DO jh = 0, kextj
            DO ji = 2, jpiglo
              ijt = jpiglo - ji + 2
              ptab(ji, ipj + jh) = psgn * ptab(ijt, ipj - 2 - jh)
            END DO
            ptab(1, ipj + jh) = psgn * ptab(3, ipj - 2 - jh)
          END DO
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipjm1) = psgn * ptab(ijt, ipjm1)
          END DO
        CASE ('U')
          DO jh = 0, kextj
            DO ji = 2, jpiglo - 1
              iju = jpiglo - ji + 1
              ptab(ji, ipj + jh) = psgn * ptab(iju, ipj - 2 - jh)
            END DO
            ptab(1, ipj + jh) = psgn * ptab(2, ipj - 2 - jh)
            ptab(jpiglo, ipj + jh) = psgn * ptab(jpiglo - 1, ipj - 2 - jh)
          END DO
          DO ji = jpiglo / 2, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipjm1) = psgn * ptab(iju, ipjm1)
          END DO
        CASE ('V')
          DO jh = 0, kextj
            DO ji = 2, jpiglo
              ijt = jpiglo - ji + 2
              ptab(ji, ipj - 1 + jh) = psgn * ptab(ijt, ipj - 2 - jh)
              ptab(ji, ipj + jh) = psgn * ptab(ijt, ipj - 3 - jh)
            END DO
            ptab(1, ipj + jh) = psgn * ptab(3, ipj - 3 - jh)
          END DO
        CASE ('F')
          DO jh = 0, kextj
            DO ji = 1, jpiglo - 1
              iju = jpiglo - ji + 1
              ptab(ji, ipj - 1 + jh) = psgn * ptab(iju, ipj - 2 - jh)
              ptab(ji, ipj + jh) = psgn * ptab(iju, ipj - 3 - jh)
            END DO
          END DO
          DO jh = 0, kextj
            ptab(1, ipj + jh) = psgn * ptab(2, ipj - 3 - jh)
            ptab(jpiglo, ipj + jh) = psgn * ptab(jpiglo - 1, ipj - 3 - jh)
          END DO
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          DO jh = 0, kextj
            DO ji = 1, jpiglo
              ijt = jpiglo - ji + 1
              ptab(ji, ipj + jh) = psgn * ptab(ijt, ipj - 1 - jh)
            END DO
          END DO
        CASE ('U')
          DO jh = 0, kextj
            DO ji = 1, jpiglo - 1
              iju = jpiglo - ji
              ptab(ji, ipj + jh) = psgn * ptab(iju, ipj - 1 - jh)
            END DO
            ptab(jpiglo, ipj + jh) = psgn * ptab(jpiglo - 2, ipj - 1 - jh)
          END DO
        CASE ('V')
          DO jh = 0, kextj
            DO ji = 1, jpiglo
              ijt = jpiglo - ji + 1
              ptab(ji, ipj + jh) = psgn * ptab(ijt, ipj - 2 - jh)
            END DO
          END DO
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipjm1) = psgn * ptab(ijt, ipjm1)
          END DO
        CASE ('F')
          DO jh = 0, kextj
            DO ji = 1, jpiglo - 1
              iju = jpiglo - ji
              ptab(ji, ipj + jh) = psgn * ptab(iju, ipj - 2 - jh)
            END DO
            ptab(jpiglo, ipj + jh) = psgn * ptab(jpiglo - 2, ipj - 2 - jh)
          END DO
          DO ji = jpiglo / 2 + 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipjm1) = psgn * ptab(iju, ipjm1)
          END DO
        END SELECT
      CASE DEFAULT
        SELECT CASE (cd_nat)
        CASE ('T', 'U', 'V', 'W')
          ptab(:, 1 : 1 - kextj) = 0._wp
          ptab(:, ipj : ipj + kextj) = 0._wp
        CASE ('F')
          ptab(:, ipj : ipj + kextj) = 0._wp
        END SELECT
      END SELECT
    END DO
    !$ACC END KERNELS
  END SUBROUTINE lbc_nfd_2d_ext
  SUBROUTINE lbc_nfd_3d(ptab, cd_nat, psgn)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ipjm1
    !$ACC KERNELS
    ipk = SIZE(ptab, 3)
    ipl = 1
    ipf = 1
    SELECT CASE (jpni)
    CASE (1)
      ipj = nlcj
    CASE DEFAULT
      ipj = 4
    END SELECT
    ipjm1 = ipj - 1
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipj, :) = psgn * ptab(ijt, ipj - 2, :)
          END DO
          ptab(1, ipj, :) = psgn * ptab(3, ipj - 2, :)
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipjm1, :) = psgn * ptab(ijt, ipjm1, :)
          END DO
        CASE ('U')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipj, :) = psgn * ptab(iju, ipj - 2, :)
          END DO
          ptab(1, ipj, :) = psgn * ptab(2, ipj - 2, :)
          ptab(jpiglo, ipj, :) = psgn * ptab(jpiglo - 1, ipj - 2, :)
          DO ji = jpiglo / 2, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipjm1, :) = psgn * ptab(iju, ipjm1, :)
          END DO
        CASE ('V')
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipj - 1, :) = psgn * ptab(ijt, ipj - 2, :)
            ptab(ji, ipj, :) = psgn * ptab(ijt, ipj - 3, :)
          END DO
          ptab(1, ipj, :) = psgn * ptab(3, ipj - 3, :)
        CASE ('F')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipj - 1, :) = psgn * ptab(iju, ipj - 2, :)
            ptab(ji, ipj, :) = psgn * ptab(iju, ipj - 3, :)
          END DO
          ptab(1, ipj, :) = psgn * ptab(2, ipj - 3, :)
          ptab(jpiglo, ipj, :) = psgn * ptab(jpiglo - 1, ipj - 3, :)
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipj, :) = psgn * ptab(ijt, ipj - 1, :)
          END DO
        CASE ('U')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipj, :) = psgn * ptab(iju, ipj - 1, :)
          END DO
          ptab(jpiglo, ipj, :) = psgn * ptab(jpiglo - 2, ipj - 1, :)
        CASE ('V')
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipj, :) = psgn * ptab(ijt, ipj - 2, :)
          END DO
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipjm1, :) = psgn * ptab(ijt, ipjm1, :)
          END DO
        CASE ('F')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipj, :) = psgn * ptab(iju, ipj - 2, :)
          END DO
          ptab(jpiglo, ipj, :) = psgn * ptab(jpiglo - 2, ipj - 2, :)
          DO ji = jpiglo / 2 + 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipjm1, :) = psgn * ptab(iju, ipjm1, :)
          END DO
        END SELECT
      CASE DEFAULT
        SELECT CASE (cd_nat)
        CASE ('T', 'U', 'V', 'W')
          ptab(:, 1, :) = 0._wp
          ptab(:, ipj, :) = 0._wp
        CASE ('F')
          ptab(:, ipj, :) = 0._wp
        END SELECT
      END SELECT
    END DO
    !$ACC END KERNELS
  END SUBROUTINE lbc_nfd_3d
  SUBROUTINE lbc_nfd_3d_ptr(ptab, cd_nat, psgn, kfld)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kfld
    TYPE(PTR_3D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ipjm1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('lbc_nfd_3d_ptr', 'r0', 0, 0)
    ipk = SIZE(ptab(1) % pt3d, 3)
    ipl = 1
    ipf = kfld
    SELECT CASE (jpni)
    CASE (1)
      ipj = nlcj
    CASE DEFAULT
      ipj = 4
    END SELECT
    ipjm1 = ipj - 1
    CALL profile_psy_data0 % PostEnd
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'W')
          !$ACC KERNELS
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt3d(ji, ipj, :) = psgn(jf) * ptab(jf) % pt3d(ijt, ipj - 2, :)
          END DO
          ptab(jf) % pt3d(1, ipj, :) = psgn(jf) * ptab(jf) % pt3d(3, ipj - 2, :)
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt3d(ji, ipjm1, :) = psgn(jf) * ptab(jf) % pt3d(ijt, ipjm1, :)
          END DO
          !$ACC END KERNELS
        CASE ('U')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt3d(ji, ipj, :) = psgn(jf) * ptab(jf) % pt3d(iju, ipj - 2, :)
          END DO
          ptab(jf) % pt3d(1, ipj, :) = psgn(jf) * ptab(jf) % pt3d(2, ipj - 2, :)
          ptab(jf) % pt3d(jpiglo, ipj, :) = psgn(jf) * ptab(jf) % pt3d(jpiglo - 1, ipj - 2, :)
          DO ji = jpiglo / 2, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt3d(ji, ipjm1, :) = psgn(jf) * ptab(jf) % pt3d(iju, ipjm1, :)
          END DO
          !$ACC END KERNELS
        CASE ('V')
          !$ACC KERNELS
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt3d(ji, ipj - 1, :) = psgn(jf) * ptab(jf) % pt3d(ijt, ipj - 2, :)
            ptab(jf) % pt3d(ji, ipj, :) = psgn(jf) * ptab(jf) % pt3d(ijt, ipj - 3, :)
          END DO
          ptab(jf) % pt3d(1, ipj, :) = psgn(jf) * ptab(jf) % pt3d(3, ipj - 3, :)
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt3d(ji, ipj - 1, :) = psgn(jf) * ptab(jf) % pt3d(iju, ipj - 2, :)
            ptab(jf) % pt3d(ji, ipj, :) = psgn(jf) * ptab(jf) % pt3d(iju, ipj - 3, :)
          END DO
          ptab(jf) % pt3d(1, ipj, :) = psgn(jf) * ptab(jf) % pt3d(2, ipj - 3, :)
          ptab(jf) % pt3d(jpiglo, ipj, :) = psgn(jf) * ptab(jf) % pt3d(jpiglo - 1, ipj - 3, :)
          !$ACC END KERNELS
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'W')
          !$ACC KERNELS ! CDe      
          DO jk = 1, jpk
            !$ACC LOOP INDEPENDENT ! CDe This runs OK and verifies
            DO ji = 1, jpiglo
              ijt = jpiglo - ji + 1
              ptab(jf) % pt3d(ji, ipj, jk) = psgn(jf) * ptab(jf) % pt3d(ijt, ipj - 1, jk)
            END DO
          END DO
          !$ACC END KERNELS
        CASE ('U')
          !$ACC KERNELS ! CDe
          DO jk = 1, jpk
            !$ACC LOOP INDEPENDENT ! CDe verifies
            DO ji = 1, jpiglo - 1
              iju = jpiglo - ji
              ptab(jf) % pt3d(ji, ipj, jk) = psgn(jf) * ptab(jf) % pt3d(iju, ipj - 1, jk)
            END DO
          END DO
          DO jk = 1, jpk
            ptab(jf) % pt3d(jpiglo, ipj, jk) = psgn(jf) * ptab(jf) % pt3d(jpiglo - 2, ipj - 1, jk)
          END DO
          !$ACC END KERNELS
        CASE ('V')
          !$ACC KERNELS ! CDe
          DO jk = 1, jpk
            !$ACC LOOP INDEPENDENT ! CDe verifies
            DO ji = 1, jpiglo
              ijt = jpiglo - ji + 1
              ptab(jf) % pt3d(ji, ipj, jk) = psgn(jf) * ptab(jf) % pt3d(ijt, ipj - 2, jk)
            END DO
            !$ACC LOOP INDEPENDENT ! CDe verifies
            DO ji = jpiglo / 2 + 1, jpiglo
              ijt = jpiglo - ji + 1
              ptab(jf) % pt3d(ji, ipjm1, jk) = psgn(jf) * ptab(jf) % pt3d(ijt, ipjm1, jk)
            END DO
          END DO
          !$ACC END KERNELS
        CASE ('F')      
          !$ACC KERNELS
          DO jk =1, jpk      
            DO ji = 1, jpiglo - 1
              iju = jpiglo - ji
              ptab(jf) % pt3d(ji, ipj, jk) = psgn(jf) * ptab(jf) % pt3d(iju, ipj - 2, jk)
            END DO
            ptab(jf) % pt3d(jpiglo, ipj, jk) = psgn(jf) * ptab(jf) % pt3d(jpiglo - 2, ipj - 2, jk)
            DO ji = jpiglo / 2 + 1, jpiglo - 1
              iju = jpiglo - ji
              ptab(jf) % pt3d(ji, ipjm1, jk) = psgn(jf) * ptab(jf) % pt3d(iju, ipjm1, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END SELECT
      CASE DEFAULT
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'U', 'V', 'W')
          !$ACC KERNELS
          ptab(jf) % pt3d(:, 1, :) = 0._wp
          ptab(jf) % pt3d(:, ipj, :) = 0._wp
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          ptab(jf) % pt3d(:, ipj, :) = 0._wp
          !$ACC END KERNELS
        END SELECT
      END SELECT
    END DO
  END SUBROUTINE lbc_nfd_3d_ptr
  SUBROUTINE lbc_nfd_4d(ptab, cd_nat, psgn)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ipjm1
    !$ACC KERNELS
    ipk = SIZE(ptab, 3)
    ipl = SIZE(ptab, 4)
    ipf = 1
    SELECT CASE (jpni)
    CASE (1)
      ipj = nlcj
    CASE DEFAULT
      ipj = 4
    END SELECT
    ipjm1 = ipj - 1
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipj, :, :) = psgn * ptab(ijt, ipj - 2, :, :)
          END DO
          ptab(1, ipj, :, :) = psgn * ptab(3, ipj - 2, :, :)
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipjm1, :, :) = psgn * ptab(ijt, ipjm1, :, :)
          END DO
        CASE ('U')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipj, :, :) = psgn * ptab(iju, ipj - 2, :, :)
          END DO
          ptab(1, ipj, :, :) = psgn * ptab(2, ipj - 2, :, :)
          ptab(jpiglo, ipj, :, :) = psgn * ptab(jpiglo - 1, ipj - 2, :, :)
          DO ji = jpiglo / 2, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipjm1, :, :) = psgn * ptab(iju, ipjm1, :, :)
          END DO
        CASE ('V')
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipj - 1, :, :) = psgn * ptab(ijt, ipj - 2, :, :)
            ptab(ji, ipj, :, :) = psgn * ptab(ijt, ipj - 3, :, :)
          END DO
          ptab(1, ipj, :, :) = psgn * ptab(3, ipj - 3, :, :)
        CASE ('F')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipj - 1, :, :) = psgn * ptab(iju, ipj - 2, :, :)
            ptab(ji, ipj, :, :) = psgn * ptab(iju, ipj - 3, :, :)
          END DO
          ptab(1, ipj, :, :) = psgn * ptab(2, ipj - 3, :, :)
          ptab(jpiglo, ipj, :, :) = psgn * ptab(jpiglo - 1, ipj - 3, :, :)
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipj, :, :) = psgn * ptab(ijt, ipj - 1, :, :)
          END DO
        CASE ('U')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipj, :, :) = psgn * ptab(iju, ipj - 1, :, :)
          END DO
          ptab(jpiglo, ipj, :, :) = psgn * ptab(jpiglo - 2, ipj - 1, :, :)
        CASE ('V')
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipj, :, :) = psgn * ptab(ijt, ipj - 2, :, :)
          END DO
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipjm1, :, :) = psgn * ptab(ijt, ipjm1, :, :)
          END DO
        CASE ('F')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipj, :, :) = psgn * ptab(iju, ipj - 2, :, :)
          END DO
          ptab(jpiglo, ipj, :, :) = psgn * ptab(jpiglo - 2, ipj - 2, :, :)
          DO ji = jpiglo / 2 + 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipjm1, :, :) = psgn * ptab(iju, ipjm1, :, :)
          END DO
        END SELECT
      CASE DEFAULT
        SELECT CASE (cd_nat)
        CASE ('T', 'U', 'V', 'W')
          ptab(:, 1, :, :) = 0._wp
          ptab(:, ipj, :, :) = 0._wp
        CASE ('F')
          ptab(:, ipj, :, :) = 0._wp
        END SELECT
      END SELECT
    END DO
    !$ACC END KERNELS
  END SUBROUTINE lbc_nfd_4d
  SUBROUTINE lbc_nfd_4d_ptr(ptab, cd_nat, psgn, kfld)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kfld
    TYPE(PTR_4D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ipjm1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('lbc_nfd_4d_ptr', 'r0', 0, 0)
    ipk = SIZE(ptab(1) % pt4d, 3)
    ipl = SIZE(ptab(1) % pt4d, 4)
    ipf = kfld
    SELECT CASE (jpni)
    CASE (1)
      ipj = nlcj
    CASE DEFAULT
      ipj = 4
    END SELECT
    ipjm1 = ipj - 1
    CALL profile_psy_data0 % PostEnd
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'W')
          !$ACC KERNELS
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt4d(ji, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(ijt, ipj - 2, :, :)
          END DO
          ptab(jf) % pt4d(1, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(3, ipj - 2, :, :)
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt4d(ji, ipjm1, :, :) = psgn(jf) * ptab(jf) % pt4d(ijt, ipjm1, :, :)
          END DO
          !$ACC END KERNELS
        CASE ('U')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt4d(ji, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(iju, ipj - 2, :, :)
          END DO
          ptab(jf) % pt4d(1, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(2, ipj - 2, :, :)
          ptab(jf) % pt4d(jpiglo, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(jpiglo - 1, ipj - 2, :, :)
          DO ji = jpiglo / 2, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt4d(ji, ipjm1, :, :) = psgn(jf) * ptab(jf) % pt4d(iju, ipjm1, :, :)
          END DO
          !$ACC END KERNELS
        CASE ('V')
          !$ACC KERNELS
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt4d(ji, ipj - 1, :, :) = psgn(jf) * ptab(jf) % pt4d(ijt, ipj - 2, :, :)
            ptab(jf) % pt4d(ji, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(ijt, ipj - 3, :, :)
          END DO
          ptab(jf) % pt4d(1, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(3, ipj - 3, :, :)
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt4d(ji, ipj - 1, :, :) = psgn(jf) * ptab(jf) % pt4d(iju, ipj - 2, :, :)
            ptab(jf) % pt4d(ji, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(iju, ipj - 3, :, :)
          END DO
          ptab(jf) % pt4d(1, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(2, ipj - 3, :, :)
          ptab(jf) % pt4d(jpiglo, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(jpiglo - 1, ipj - 3, :, :)
          !$ACC END KERNELS
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'W')
          !$ACC KERNELS
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(jf) % pt4d(ji, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(ijt, ipj - 1, :, :)
          END DO
          !$ACC END KERNELS
        CASE ('U')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(jf) % pt4d(ji, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(iju, ipj - 1, :, :)
          END DO
          ptab(jf) % pt4d(jpiglo, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(jpiglo - 2, ipj - 1, :, :)
          !$ACC END KERNELS
        CASE ('V')
          !$ACC KERNELS
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(jf) % pt4d(ji, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(ijt, ipj - 2, :, :)
          END DO
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(jf) % pt4d(ji, ipjm1, :, :) = psgn(jf) * ptab(jf) % pt4d(ijt, ipjm1, :, :)
          END DO
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(jf) % pt4d(ji, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(iju, ipj - 2, :, :)
          END DO
          ptab(jf) % pt4d(jpiglo, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(jpiglo - 2, ipj - 2, :, :)
          DO ji = jpiglo / 2 + 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(jf) % pt4d(ji, ipjm1, :, :) = psgn(jf) * ptab(jf) % pt4d(iju, ipjm1, :, :)
          END DO
          !$ACC END KERNELS
        END SELECT
      CASE DEFAULT
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'U', 'V', 'W')
          !$ACC KERNELS
          ptab(jf) % pt4d(:, 1, :, :) = 0._wp
          ptab(jf) % pt4d(:, ipj, :, :) = 0._wp
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          ptab(jf) % pt4d(:, ipj, :, :) = 0._wp
          !$ACC END KERNELS
        END SELECT
      END SELECT
    END DO
  END SUBROUTINE lbc_nfd_4d_ptr
  SUBROUTINE lbc_nfd_nogather_2d(ptab, ptab2, cd_nat, psgn, kfld)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :)
    REAL(KIND = wp), INTENT(INOUT) :: ptab2(:, :, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ijpj, ijpjp1, ijta, ijua, jia, startloop, endloop
    LOGICAL :: l_fast_exchanges
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
    CALL profile_psy_data0 % PreStart('lbc_nfd_nogather_2d', 'r0', 0, 0)
    ipj = SIZE(ptab2, 2)
    ipk = 1
    ipl = 1
    ipf = 1
    IF (ipf > 1) THEN
      WRITE(6, *) 'lbc_nfd_nogather: multiple fields not allowed. Revise implementation'
      WRITE(6, *) 'You should not be there...'
      STOP
    END IF
    ijpj = 1
    ijpjp1 = 2
    IF (ipj == 1) THEN
      l_fast_exchanges = .TRUE.
    ELSE
      l_fast_exchanges = .FALSE.
    END IF
    CALL profile_psy_data0 % PostEnd
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          CALL profile_psy_data1 % PreStart('lbc_nfd_nogather_2d', 'r1', 0, 0)
          IF (nimpp /= 1) THEN
            startloop = 1
          ELSE
            startloop = 2
          END IF
          CALL profile_psy_data1 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = startloop, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                ptab(ji, nlcj) = psgn * ptab2(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF (nimpp == 1) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                ptab(1, nlcj) = psgn * ptab(3, nlcj - 2)
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data2 % PreStart('lbc_nfd_nogather_2d', 'r2', 0, 0)
            IF (nimpp >= jpiglo / 2 + 1) THEN
              startloop = 1
            ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
              startloop = jpiglo / 2 + 1 - nimpp + 1
            ELSE
              startloop = nlci + 1
            END IF
            CALL profile_psy_data2 % PostEnd
            IF (startloop <= nlci) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, nlci
                    ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                    jia = ji + nimpp - 1
                    ijta = jpiglo - jia + 2
                    IF (ijta >= startloop + nimpp - 1 .AND. ijta < jia) THEN
                      ptab(ji, nlcj - 1) = psgn * ptab(ijta - nimpp + 1, nlcj - 1)
                    ELSE
                      ptab(ji, nlcj - 1) = psgn * ptab2(ijt, ijpjp1, jk, jl)
                    END IF
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        CASE ('U')
          CALL profile_psy_data3 % PreStart('lbc_nfd_nogather_2d', 'r3', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data3 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(ji, nlcj) = psgn * ptab2(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          CALL profile_psy_data4 % PreStart('lbc_nfd_nogather_2d', 'r4', 0, 0)
          IF (nimpp .EQ. 1) THEN
            ptab(1, nlcj) = psgn * ptab(2, nlcj - 2)
          END IF
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            ptab(nlci, nlcj) = psgn * ptab(nlci - 1, nlcj - 2)
          END IF
          CALL profile_psy_data4 % PostEnd
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data5 % PreStart('lbc_nfd_nogather_2d', 'r5', 0, 0)
            IF (nimpp + nlci - 1 /= jpiglo) THEN
              endloop = nlci
            ELSE
              endloop = nlci - 1
            END IF
            IF (nimpp >= jpiglo / 2) THEN
              startloop = 1
            ELSE IF ((nimpp + nlci - 1 >= jpiglo / 2) .AND. (nimpp < jpiglo / 2)) THEN
              startloop = jpiglo / 2 - nimpp + 1
            ELSE
              startloop = endloop + 1
            END IF
            CALL profile_psy_data5 % PostEnd
            IF (startloop <= endloop) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, endloop
                    iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                    jia = ji + nimpp - 1
                    ijua = jpiglo - jia + 1
                    IF (ijua >= startloop + nimpp - 1 .AND. ijua < jia) THEN
                      ptab(ji, nlcj - 1) = psgn * ptab(ijua - nimpp + 1, nlcj - 1)
                    ELSE
                      ptab(ji, nlcj - 1) = psgn * ptab2(iju, ijpjp1, jk, jl)
                    END IF
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        CASE ('V')
          CALL profile_psy_data6 % PreStart('lbc_nfd_nogather_2d', 'r6', 0, 0)
          IF (nimpp /= 1) THEN
            startloop = 1
          ELSE
            startloop = 2
          END IF
          CALL profile_psy_data6 % PostEnd
          IF (.NOT. l_fast_exchanges) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                DO ji = startloop, nlci
                  ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                  ptab(ji, nlcj - 1) = psgn * ptab2(ijt, ijpjp1, jk, jl)
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = startloop, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                ptab(ji, nlcj) = psgn * ptab2(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          CALL profile_psy_data7 % PreStart('lbc_nfd_nogather_2d', 'r7', 0, 0)
          IF (nimpp .EQ. 1) THEN
            ptab(1, nlcj) = psgn * ptab(3, nlcj - 3)
          END IF
          CALL profile_psy_data7 % PostEnd
        CASE ('F')
          CALL profile_psy_data8 % PreStart('lbc_nfd_nogather_2d', 'r8', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data8 % PostEnd
          IF (.NOT. l_fast_exchanges) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                DO ji = 1, endloop
                  iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                  ptab(ji, nlcj - 1) = psgn * ptab2(iju, ijpjp1, jk, jl)
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(ji, nlcj) = psgn * ptab2(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          CALL profile_psy_data9 % PreStart('lbc_nfd_nogather_2d', 'r9', 0, 0)
          IF (nimpp .EQ. 1) THEN
            ptab(1, nlcj) = psgn * ptab(2, nlcj - 3)
            IF (.NOT. l_fast_exchanges) ptab(1, nlcj - 1) = psgn * ptab(2, nlcj - 2)
          END IF
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            ptab(nlci, nlcj) = psgn * ptab(nlci - 1, nlcj - 3)
            IF (.NOT. l_fast_exchanges) ptab(nlci, nlcj - 1) = psgn * ptab(nlci - 1, nlcj - 2)
          END IF
          CALL profile_psy_data9 % PostEnd
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(ji, nlcj) = psgn * ptab2(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        CASE ('U')
          CALL profile_psy_data10 % PreStart('lbc_nfd_nogather_2d', 'r10', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data10 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
                ptab(ji, nlcj) = psgn * ptab2(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                ptab(nlci, nlcj) = psgn * ptab(nlci - 2, nlcj - 1)
              END DO
              !$ACC END KERNELS
            END DO
          END IF
        CASE ('V')
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(ji, nlcj) = psgn * ptab2(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data11 % PreStart('lbc_nfd_nogather_2d', 'r11', 0, 0)
            IF (nimpp >= jpiglo / 2 + 1) THEN
              startloop = 1
            ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
              startloop = jpiglo / 2 + 1 - nimpp + 1
            ELSE
              startloop = nlci + 1
            END IF
            CALL profile_psy_data11 % PostEnd
            IF (startloop <= nlci) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, nlci
                    ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                    ptab(ji, nlcj - 1) = psgn * ptab2(ijt, ijpjp1, jk, jl)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        CASE ('F')
          CALL profile_psy_data12 % PreStart('lbc_nfd_nogather_2d', 'r12', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data12 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
                ptab(ji, nlcj) = psgn * ptab2(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                ptab(nlci, nlcj) = psgn * ptab(nlci - 2, nlcj - 2)
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data13 % PreStart('lbc_nfd_nogather_2d', 'r13', 0, 0)
            IF (nimpp + nlci - 1 /= jpiglo) THEN
              endloop = nlci
            ELSE
              endloop = nlci - 1
            END IF
            IF (nimpp >= jpiglo / 2 + 1) THEN
              startloop = 1
            ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
              startloop = jpiglo / 2 + 1 - nimpp + 1
            ELSE
              startloop = endloop + 1
            END IF
            CALL profile_psy_data13 % PostEnd
            IF (startloop <= endloop) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, endloop
                    iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
                    ptab(ji, nlcj - 1) = psgn * ptab2(iju, ijpjp1, jk, jl)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        END SELECT
      CASE DEFAULT
        WRITE(*, *) 'lbc_nfd_nogather_generic: You should not have seen this print! error?', npolj
      END SELECT
    END DO
  END SUBROUTINE lbc_nfd_nogather_2d
  SUBROUTINE lbc_nfd_nogather_2d_ptr(ptab, ptab2, cd_nat, psgn, kfld)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(PTR_2D), INTENT(INOUT) :: ptab(:)
    TYPE(PTR_4D), INTENT(INOUT) :: ptab2(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ijpj, ijpjp1, ijta, ijua, jia, startloop, endloop
    LOGICAL :: l_fast_exchanges
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
    CALL profile_psy_data0 % PreStart('lbc_nfd_nogather_2d_ptr', 'r0', 0, 0)
    ipj = SIZE(ptab2(1) % pt4d, 2)
    ipk = 1
    ipl = 1
    ipf = kfld
    IF (ipf > 1) THEN
      WRITE(6, *) 'lbc_nfd_nogather: multiple fields not allowed. Revise implementation'
      WRITE(6, *) 'You should not be there...'
      STOP
    END IF
    ijpj = 1
    ijpjp1 = 2
    IF (ipj == 1) THEN
      l_fast_exchanges = .TRUE.
    ELSE
      l_fast_exchanges = .FALSE.
    END IF
    CALL profile_psy_data0 % PostEnd
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'W')
          CALL profile_psy_data1 % PreStart('lbc_nfd_nogather_2d_ptr', 'r1', 0, 0)
          IF (nimpp /= 1) THEN
            startloop = 1
          ELSE
            startloop = 2
          END IF
          CALL profile_psy_data1 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = startloop, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                ptab(jf) % pt2d(ji, nlcj) = psgn(jf) * ptab2(jf) % pt4d(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF (nimpp == 1) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                ptab(jf) % pt2d(1, nlcj) = psgn(jf) * ptab(jf) % pt2d(3, nlcj - 2)
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data2 % PreStart('lbc_nfd_nogather_2d_ptr', 'r2', 0, 0)
            IF (nimpp >= jpiglo / 2 + 1) THEN
              startloop = 1
            ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
              startloop = jpiglo / 2 + 1 - nimpp + 1
            ELSE
              startloop = nlci + 1
            END IF
            CALL profile_psy_data2 % PostEnd
            IF (startloop <= nlci) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, nlci
                    ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                    jia = ji + nimpp - 1
                    ijta = jpiglo - jia + 2
                    IF (ijta >= startloop + nimpp - 1 .AND. ijta < jia) THEN
                      ptab(jf) % pt2d(ji, nlcj - 1) = psgn(jf) * ptab(jf) % pt2d(ijta - nimpp + 1, nlcj - 1)
                    ELSE
                      ptab(jf) % pt2d(ji, nlcj - 1) = psgn(jf) * ptab2(jf) % pt4d(ijt, ijpjp1, jk, jl)
                    END IF
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        CASE ('U')
          CALL profile_psy_data3 % PreStart('lbc_nfd_nogather_2d_ptr', 'r3', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data3 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(jf) % pt2d(ji, nlcj) = psgn(jf) * ptab2(jf) % pt4d(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          CALL profile_psy_data4 % PreStart('lbc_nfd_nogather_2d_ptr', 'r4', 0, 0)
          IF (nimpp .EQ. 1) THEN
            ptab(jf) % pt2d(1, nlcj) = psgn(jf) * ptab(jf) % pt2d(2, nlcj - 2)
          END IF
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            ptab(jf) % pt2d(nlci, nlcj) = psgn(jf) * ptab(jf) % pt2d(nlci - 1, nlcj - 2)
          END IF
          CALL profile_psy_data4 % PostEnd
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data5 % PreStart('lbc_nfd_nogather_2d_ptr', 'r5', 0, 0)
            IF (nimpp + nlci - 1 /= jpiglo) THEN
              endloop = nlci
            ELSE
              endloop = nlci - 1
            END IF
            IF (nimpp >= jpiglo / 2) THEN
              startloop = 1
            ELSE IF ((nimpp + nlci - 1 >= jpiglo / 2) .AND. (nimpp < jpiglo / 2)) THEN
              startloop = jpiglo / 2 - nimpp + 1
            ELSE
              startloop = endloop + 1
            END IF
            CALL profile_psy_data5 % PostEnd
            IF (startloop <= endloop) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, endloop
                    iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                    jia = ji + nimpp - 1
                    ijua = jpiglo - jia + 1
                    IF (ijua >= startloop + nimpp - 1 .AND. ijua < jia) THEN
                      ptab(jf) % pt2d(ji, nlcj - 1) = psgn(jf) * ptab(jf) % pt2d(ijua - nimpp + 1, nlcj - 1)
                    ELSE
                      ptab(jf) % pt2d(ji, nlcj - 1) = psgn(jf) * ptab2(jf) % pt4d(iju, ijpjp1, jk, jl)
                    END IF
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        CASE ('V')
          CALL profile_psy_data6 % PreStart('lbc_nfd_nogather_2d_ptr', 'r6', 0, 0)
          IF (nimpp /= 1) THEN
            startloop = 1
          ELSE
            startloop = 2
          END IF
          CALL profile_psy_data6 % PostEnd
          IF (.NOT. l_fast_exchanges) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                DO ji = startloop, nlci
                  ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                  ptab(jf) % pt2d(ji, nlcj - 1) = psgn(jf) * ptab2(jf) % pt4d(ijt, ijpjp1, jk, jl)
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = startloop, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                ptab(jf) % pt2d(ji, nlcj) = psgn(jf) * ptab2(jf) % pt4d(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          CALL profile_psy_data7 % PreStart('lbc_nfd_nogather_2d_ptr', 'r7', 0, 0)
          IF (nimpp .EQ. 1) THEN
            ptab(jf) % pt2d(1, nlcj) = psgn(jf) * ptab(jf) % pt2d(3, nlcj - 3)
          END IF
          CALL profile_psy_data7 % PostEnd
        CASE ('F')
          CALL profile_psy_data8 % PreStart('lbc_nfd_nogather_2d_ptr', 'r8', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data8 % PostEnd
          IF (.NOT. l_fast_exchanges) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                DO ji = 1, endloop
                  iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                  ptab(jf) % pt2d(ji, nlcj - 1) = psgn(jf) * ptab2(jf) % pt4d(iju, ijpjp1, jk, jl)
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(jf) % pt2d(ji, nlcj) = psgn(jf) * ptab2(jf) % pt4d(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          CALL profile_psy_data9 % PreStart('lbc_nfd_nogather_2d_ptr', 'r9', 0, 0)
          IF (nimpp .EQ. 1) THEN
            ptab(jf) % pt2d(1, nlcj) = psgn(jf) * ptab(jf) % pt2d(2, nlcj - 3)
            IF (.NOT. l_fast_exchanges) ptab(jf) % pt2d(1, nlcj - 1) = psgn(jf) * ptab(jf) % pt2d(2, nlcj - 2)
          END IF
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            ptab(jf) % pt2d(nlci, nlcj) = psgn(jf) * ptab(jf) % pt2d(nlci - 1, nlcj - 3)
            IF (.NOT. l_fast_exchanges) ptab(jf) % pt2d(nlci, nlcj - 1) = psgn(jf) * ptab(jf) % pt2d(nlci - 1, nlcj - 2)
          END IF
          CALL profile_psy_data9 % PostEnd
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'W')
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(jf) % pt2d(ji, nlcj) = psgn(jf) * ptab2(jf) % pt4d(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        CASE ('U')
          CALL profile_psy_data10 % PreStart('lbc_nfd_nogather_2d_ptr', 'r10', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data10 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
                ptab(jf) % pt2d(ji, nlcj) = psgn(jf) * ptab2(jf) % pt4d(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                ptab(jf) % pt2d(nlci, nlcj) = psgn(jf) * ptab(jf) % pt2d(nlci - 2, nlcj - 1)
              END DO
              !$ACC END KERNELS
            END DO
          END IF
        CASE ('V')
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(jf) % pt2d(ji, nlcj) = psgn(jf) * ptab2(jf) % pt4d(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data11 % PreStart('lbc_nfd_nogather_2d_ptr', 'r11', 0, 0)
            IF (nimpp >= jpiglo / 2 + 1) THEN
              startloop = 1
            ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
              startloop = jpiglo / 2 + 1 - nimpp + 1
            ELSE
              startloop = nlci + 1
            END IF
            CALL profile_psy_data11 % PostEnd
            IF (startloop <= nlci) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, nlci
                    ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                    ptab(jf) % pt2d(ji, nlcj - 1) = psgn(jf) * ptab2(jf) % pt4d(ijt, ijpjp1, jk, jl)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        CASE ('F')
          CALL profile_psy_data12 % PreStart('lbc_nfd_nogather_2d_ptr', 'r12', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data12 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
                ptab(jf) % pt2d(ji, nlcj) = psgn(jf) * ptab2(jf) % pt4d(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                ptab(jf) % pt2d(nlci, nlcj) = psgn(jf) * ptab(jf) % pt2d(nlci - 2, nlcj - 2)
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data13 % PreStart('lbc_nfd_nogather_2d_ptr', 'r13', 0, 0)
            IF (nimpp + nlci - 1 /= jpiglo) THEN
              endloop = nlci
            ELSE
              endloop = nlci - 1
            END IF
            IF (nimpp >= jpiglo / 2 + 1) THEN
              startloop = 1
            ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
              startloop = jpiglo / 2 + 1 - nimpp + 1
            ELSE
              startloop = endloop + 1
            END IF
            CALL profile_psy_data13 % PostEnd
            IF (startloop <= endloop) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, endloop
                    iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
                    ptab(jf) % pt2d(ji, nlcj - 1) = psgn(jf) * ptab2(jf) % pt4d(iju, ijpjp1, jk, jl)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        END SELECT
      CASE DEFAULT
        WRITE(*, *) 'lbc_nfd_nogather_generic: You should not have seen this print! error?', npolj
      END SELECT
    END DO
  END SUBROUTINE lbc_nfd_nogather_2d_ptr
  SUBROUTINE lbc_nfd_nogather_3d(ptab, ptab2, cd_nat, psgn, kfld)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :)
    REAL(KIND = wp), INTENT(INOUT) :: ptab2(:, :, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ijpj, ijpjp1, ijta, ijua, jia, startloop, endloop
    LOGICAL :: l_fast_exchanges
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
    CALL profile_psy_data0 % PreStart('lbc_nfd_nogather_3d', 'r0', 0, 0)
    ipj = SIZE(ptab2, 2)
    ipk = SIZE(ptab, 3)
    ipl = 1
    ipf = 1
    IF (ipf > 1) THEN
      WRITE(6, *) 'lbc_nfd_nogather: multiple fields not allowed. Revise implementation'
      WRITE(6, *) 'You should not be there...'
      STOP
    END IF
    ijpj = 1
    ijpjp1 = 2
    IF (ipj == 1) THEN
      l_fast_exchanges = .TRUE.
    ELSE
      l_fast_exchanges = .FALSE.
    END IF
    CALL profile_psy_data0 % PostEnd
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          CALL profile_psy_data1 % PreStart('lbc_nfd_nogather_3d', 'r1', 0, 0)
          IF (nimpp /= 1) THEN
            startloop = 1
          ELSE
            startloop = 2
          END IF
          CALL profile_psy_data1 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = startloop, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                ptab(ji, nlcj, jk) = psgn * ptab2(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF (nimpp == 1) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                ptab(1, nlcj, jk) = psgn * ptab(3, nlcj - 2, jk)
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data2 % PreStart('lbc_nfd_nogather_3d', 'r2', 0, 0)
            IF (nimpp >= jpiglo / 2 + 1) THEN
              startloop = 1
            ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
              startloop = jpiglo / 2 + 1 - nimpp + 1
            ELSE
              startloop = nlci + 1
            END IF
            CALL profile_psy_data2 % PostEnd
            IF (startloop <= nlci) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, nlci
                    ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                    jia = ji + nimpp - 1
                    ijta = jpiglo - jia + 2
                    IF (ijta >= startloop + nimpp - 1 .AND. ijta < jia) THEN
                      ptab(ji, nlcj - 1, jk) = psgn * ptab(ijta - nimpp + 1, nlcj - 1, jk)
                    ELSE
                      ptab(ji, nlcj - 1, jk) = psgn * ptab2(ijt, ijpjp1, jk, jl)
                    END IF
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        CASE ('U')
          CALL profile_psy_data3 % PreStart('lbc_nfd_nogather_3d', 'r3', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data3 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(ji, nlcj, jk) = psgn * ptab2(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          !$ACC KERNELS
          IF (nimpp .EQ. 1) THEN
            ptab(1, nlcj, :) = psgn * ptab(2, nlcj - 2, :)
          END IF
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            ptab(nlci, nlcj, :) = psgn * ptab(nlci - 1, nlcj - 2, :)
          END IF
          !$ACC END KERNELS
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data4 % PreStart('lbc_nfd_nogather_3d', 'r4', 0, 0)
            IF (nimpp + nlci - 1 /= jpiglo) THEN
              endloop = nlci
            ELSE
              endloop = nlci - 1
            END IF
            IF (nimpp >= jpiglo / 2) THEN
              startloop = 1
            ELSE IF ((nimpp + nlci - 1 >= jpiglo / 2) .AND. (nimpp < jpiglo / 2)) THEN
              startloop = jpiglo / 2 - nimpp + 1
            ELSE
              startloop = endloop + 1
            END IF
            CALL profile_psy_data4 % PostEnd
            IF (startloop <= endloop) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, endloop
                    iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                    jia = ji + nimpp - 1
                    ijua = jpiglo - jia + 1
                    IF (ijua >= startloop + nimpp - 1 .AND. ijua < jia) THEN
                      ptab(ji, nlcj - 1, jk) = psgn * ptab(ijua - nimpp + 1, nlcj - 1, jk)
                    ELSE
                      ptab(ji, nlcj - 1, jk) = psgn * ptab2(iju, ijpjp1, jk, jl)
                    END IF
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        CASE ('V')
          CALL profile_psy_data5 % PreStart('lbc_nfd_nogather_3d', 'r5', 0, 0)
          IF (nimpp /= 1) THEN
            startloop = 1
          ELSE
            startloop = 2
          END IF
          CALL profile_psy_data5 % PostEnd
          IF (.NOT. l_fast_exchanges) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                DO ji = startloop, nlci
                  ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                  ptab(ji, nlcj - 1, jk) = psgn * ptab2(ijt, ijpjp1, jk, jl)
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = startloop, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                ptab(ji, nlcj, jk) = psgn * ptab2(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          !$ACC KERNELS
          IF (nimpp .EQ. 1) THEN
            ptab(1, nlcj, :) = psgn * ptab(3, nlcj - 3, :)
          END IF
          !$ACC END KERNELS
        CASE ('F')
          CALL profile_psy_data6 % PreStart('lbc_nfd_nogather_3d', 'r6', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data6 % PostEnd
          IF (.NOT. l_fast_exchanges) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                DO ji = 1, endloop
                  iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                  ptab(ji, nlcj - 1, jk) = psgn * ptab2(iju, ijpjp1, jk, jl)
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(ji, nlcj, jk) = psgn * ptab2(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF (nimpp .EQ. 1) THEN
            !$ACC KERNELS
            ptab(1, nlcj, :) = psgn * ptab(2, nlcj - 3, :)
            !$ACC END KERNELS
            CALL profile_psy_data7 % PreStart('lbc_nfd_nogather_3d', 'r7', 0, 0)
            IF (.NOT. l_fast_exchanges) ptab(1, nlcj - 1, :) = psgn * ptab(2, nlcj - 2, :)
            CALL profile_psy_data7 % PostEnd
          END IF
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            !$ACC KERNELS
            ptab(nlci, nlcj, :) = psgn * ptab(nlci - 1, nlcj - 3, :)
            !$ACC END KERNELS
            CALL profile_psy_data8 % PreStart('lbc_nfd_nogather_3d', 'r8', 0, 0)
            IF (.NOT. l_fast_exchanges) ptab(nlci, nlcj - 1, :) = psgn * ptab(nlci - 1, nlcj - 2, :)
            CALL profile_psy_data8 % PostEnd
          END IF
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(ji, nlcj, jk) = psgn * ptab2(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        CASE ('U')
          CALL profile_psy_data9 % PreStart('lbc_nfd_nogather_3d', 'r9', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data9 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
                ptab(ji, nlcj, jk) = psgn * ptab2(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                ptab(nlci, nlcj, jk) = psgn * ptab(nlci - 2, nlcj - 1, jk)
              END DO
              !$ACC END KERNELS
            END DO
          END IF
        CASE ('V')
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(ji, nlcj, jk) = psgn * ptab2(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data10 % PreStart('lbc_nfd_nogather_3d', 'r10', 0, 0)
            IF (nimpp >= jpiglo / 2 + 1) THEN
              startloop = 1
            ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
              startloop = jpiglo / 2 + 1 - nimpp + 1
            ELSE
              startloop = nlci + 1
            END IF
            CALL profile_psy_data10 % PostEnd
            IF (startloop <= nlci) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, nlci
                    ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                    ptab(ji, nlcj - 1, jk) = psgn * ptab2(ijt, ijpjp1, jk, jl)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        CASE ('F')
          CALL profile_psy_data11 % PreStart('lbc_nfd_nogather_3d', 'r11', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data11 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
                ptab(ji, nlcj, jk) = psgn * ptab2(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                ptab(nlci, nlcj, jk) = psgn * ptab(nlci - 2, nlcj - 2, jk)
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data12 % PreStart('lbc_nfd_nogather_3d', 'r12', 0, 0)
            IF (nimpp + nlci - 1 /= jpiglo) THEN
              endloop = nlci
            ELSE
              endloop = nlci - 1
            END IF
            IF (nimpp >= jpiglo / 2 + 1) THEN
              startloop = 1
            ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
              startloop = jpiglo / 2 + 1 - nimpp + 1
            ELSE
              startloop = endloop + 1
            END IF
            CALL profile_psy_data12 % PostEnd
            IF (startloop <= endloop) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, endloop
                    iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
                    ptab(ji, nlcj - 1, jk) = psgn * ptab2(iju, ijpjp1, jk, jl)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        END SELECT
      CASE DEFAULT
        WRITE(*, *) 'lbc_nfd_nogather_generic: You should not have seen this print! error?', npolj
      END SELECT
    END DO
  END SUBROUTINE lbc_nfd_nogather_3d
  SUBROUTINE lbc_nfd_nogather_3d_ptr(ptab, ptab2, cd_nat, psgn, kfld)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(PTR_3D), INTENT(INOUT) :: ptab(:)
    TYPE(PTR_4D), INTENT(INOUT) :: ptab2(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ijpj, ijpjp1, ijta, ijua, jia, startloop, endloop
    LOGICAL :: l_fast_exchanges
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
    CALL profile_psy_data0 % PreStart('lbc_nfd_nogather_3d_ptr', 'r0', 0, 0)
    ipj = SIZE(ptab2(1) % pt4d, 2)
    ipk = SIZE(ptab(1) % pt3d, 3)
    ipl = 1
    ipf = kfld
    IF (ipf > 1) THEN
      WRITE(6, *) 'lbc_nfd_nogather: multiple fields not allowed. Revise implementation'
      WRITE(6, *) 'You should not be there...'
      STOP
    END IF
    ijpj = 1
    ijpjp1 = 2
    IF (ipj == 1) THEN
      l_fast_exchanges = .TRUE.
    ELSE
      l_fast_exchanges = .FALSE.
    END IF
    CALL profile_psy_data0 % PostEnd
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'W')
          CALL profile_psy_data1 % PreStart('lbc_nfd_nogather_3d_ptr', 'r1', 0, 0)
          IF (nimpp /= 1) THEN
            startloop = 1
          ELSE
            startloop = 2
          END IF
          CALL profile_psy_data1 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = startloop, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                ptab(jf) % pt3d(ji, nlcj, jk) = psgn(jf) * ptab2(jf) % pt4d(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF (nimpp == 1) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                ptab(jf) % pt3d(1, nlcj, jk) = psgn(jf) * ptab(jf) % pt3d(3, nlcj - 2, jk)
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data2 % PreStart('lbc_nfd_nogather_3d_ptr', 'r2', 0, 0)
            IF (nimpp >= jpiglo / 2 + 1) THEN
              startloop = 1
            ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
              startloop = jpiglo / 2 + 1 - nimpp + 1
            ELSE
              startloop = nlci + 1
            END IF
            CALL profile_psy_data2 % PostEnd
            IF (startloop <= nlci) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, nlci
                    ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                    jia = ji + nimpp - 1
                    ijta = jpiglo - jia + 2
                    IF (ijta >= startloop + nimpp - 1 .AND. ijta < jia) THEN
                      ptab(jf) % pt3d(ji, nlcj - 1, jk) = psgn(jf) * ptab(jf) % pt3d(ijta - nimpp + 1, nlcj - 1, jk)
                    ELSE
                      ptab(jf) % pt3d(ji, nlcj - 1, jk) = psgn(jf) * ptab2(jf) % pt4d(ijt, ijpjp1, jk, jl)
                    END IF
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        CASE ('U')
          CALL profile_psy_data3 % PreStart('lbc_nfd_nogather_3d_ptr', 'r3', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data3 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(jf) % pt3d(ji, nlcj, jk) = psgn(jf) * ptab2(jf) % pt4d(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          !$ACC KERNELS
          IF (nimpp .EQ. 1) THEN
            ptab(jf) % pt3d(1, nlcj, :) = psgn(jf) * ptab(jf) % pt3d(2, nlcj - 2, :)
          END IF
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            ptab(jf) % pt3d(nlci, nlcj, :) = psgn(jf) * ptab(jf) % pt3d(nlci - 1, nlcj - 2, :)
          END IF
          !$ACC END KERNELS
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data4 % PreStart('lbc_nfd_nogather_3d_ptr', 'r4', 0, 0)
            IF (nimpp + nlci - 1 /= jpiglo) THEN
              endloop = nlci
            ELSE
              endloop = nlci - 1
            END IF
            IF (nimpp >= jpiglo / 2) THEN
              startloop = 1
            ELSE IF ((nimpp + nlci - 1 >= jpiglo / 2) .AND. (nimpp < jpiglo / 2)) THEN
              startloop = jpiglo / 2 - nimpp + 1
            ELSE
              startloop = endloop + 1
            END IF
            CALL profile_psy_data4 % PostEnd
            IF (startloop <= endloop) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, endloop
                    iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                    jia = ji + nimpp - 1
                    ijua = jpiglo - jia + 1
                    IF (ijua >= startloop + nimpp - 1 .AND. ijua < jia) THEN
                      ptab(jf) % pt3d(ji, nlcj - 1, jk) = psgn(jf) * ptab(jf) % pt3d(ijua - nimpp + 1, nlcj - 1, jk)
                    ELSE
                      ptab(jf) % pt3d(ji, nlcj - 1, jk) = psgn(jf) * ptab2(jf) % pt4d(iju, ijpjp1, jk, jl)
                    END IF
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        CASE ('V')
          CALL profile_psy_data5 % PreStart('lbc_nfd_nogather_3d_ptr', 'r5', 0, 0)
          IF (nimpp /= 1) THEN
            startloop = 1
          ELSE
            startloop = 2
          END IF
          CALL profile_psy_data5 % PostEnd
          IF (.NOT. l_fast_exchanges) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                DO ji = startloop, nlci
                  ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                  ptab(jf) % pt3d(ji, nlcj - 1, jk) = psgn(jf) * ptab2(jf) % pt4d(ijt, ijpjp1, jk, jl)
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = startloop, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                ptab(jf) % pt3d(ji, nlcj, jk) = psgn(jf) * ptab2(jf) % pt4d(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          !$ACC KERNELS
          IF (nimpp .EQ. 1) THEN
            ptab(jf) % pt3d(1, nlcj, :) = psgn(jf) * ptab(jf) % pt3d(3, nlcj - 3, :)
          END IF
          !$ACC END KERNELS
        CASE ('F')
          CALL profile_psy_data6 % PreStart('lbc_nfd_nogather_3d_ptr', 'r6', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data6 % PostEnd
          IF (.NOT. l_fast_exchanges) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                DO ji = 1, endloop
                  iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                  ptab(jf) % pt3d(ji, nlcj - 1, jk) = psgn(jf) * ptab2(jf) % pt4d(iju, ijpjp1, jk, jl)
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(jf) % pt3d(ji, nlcj, jk) = psgn(jf) * ptab2(jf) % pt4d(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF (nimpp .EQ. 1) THEN
            !$ACC KERNELS
            ptab(jf) % pt3d(1, nlcj, :) = psgn(jf) * ptab(jf) % pt3d(2, nlcj - 3, :)
            !$ACC END KERNELS
            CALL profile_psy_data7 % PreStart('lbc_nfd_nogather_3d_ptr', 'r7', 0, 0)
            IF (.NOT. l_fast_exchanges) ptab(jf) % pt3d(1, nlcj - 1, :) = psgn(jf) * ptab(jf) % pt3d(2, nlcj - 2, :)
            CALL profile_psy_data7 % PostEnd
          END IF
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            !$ACC KERNELS
            ptab(jf) % pt3d(nlci, nlcj, :) = psgn(jf) * ptab(jf) % pt3d(nlci - 1, nlcj - 3, :)
            !$ACC END KERNELS
            CALL profile_psy_data8 % PreStart('lbc_nfd_nogather_3d_ptr', 'r8', 0, 0)
            IF (.NOT. l_fast_exchanges) ptab(jf) % pt3d(nlci, nlcj - 1, :) = psgn(jf) * ptab(jf) % pt3d(nlci - 1, nlcj - 2, :)
            CALL profile_psy_data8 % PostEnd
          END IF
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'W')
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(jf) % pt3d(ji, nlcj, jk) = psgn(jf) * ptab2(jf) % pt4d(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        CASE ('U')
          CALL profile_psy_data9 % PreStart('lbc_nfd_nogather_3d_ptr', 'r9', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data9 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
                ptab(jf) % pt3d(ji, nlcj, jk) = psgn(jf) * ptab2(jf) % pt4d(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                ptab(jf) % pt3d(nlci, nlcj, jk) = psgn(jf) * ptab(jf) % pt3d(nlci - 2, nlcj - 1, jk)
              END DO
              !$ACC END KERNELS
            END DO
          END IF
        CASE ('V')
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(jf) % pt3d(ji, nlcj, jk) = psgn(jf) * ptab2(jf) % pt4d(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data10 % PreStart('lbc_nfd_nogather_3d_ptr', 'r10', 0, 0)
            IF (nimpp >= jpiglo / 2 + 1) THEN
              startloop = 1
            ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
              startloop = jpiglo / 2 + 1 - nimpp + 1
            ELSE
              startloop = nlci + 1
            END IF
            CALL profile_psy_data10 % PostEnd
            IF (startloop <= nlci) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, nlci
                    ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                    ptab(jf) % pt3d(ji, nlcj - 1, jk) = psgn(jf) * ptab2(jf) % pt4d(ijt, ijpjp1, jk, jl)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        CASE ('F')
          CALL profile_psy_data11 % PreStart('lbc_nfd_nogather_3d_ptr', 'r11', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data11 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
                ptab(jf) % pt3d(ji, nlcj, jk) = psgn(jf) * ptab2(jf) % pt4d(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                ptab(jf) % pt3d(nlci, nlcj, jk) = psgn(jf) * ptab(jf) % pt3d(nlci - 2, nlcj - 2, jk)
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data12 % PreStart('lbc_nfd_nogather_3d_ptr', 'r12', 0, 0)
            IF (nimpp + nlci - 1 /= jpiglo) THEN
              endloop = nlci
            ELSE
              endloop = nlci - 1
            END IF
            IF (nimpp >= jpiglo / 2 + 1) THEN
              startloop = 1
            ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
              startloop = jpiglo / 2 + 1 - nimpp + 1
            ELSE
              startloop = endloop + 1
            END IF
            CALL profile_psy_data12 % PostEnd
            IF (startloop <= endloop) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, endloop
                    iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
                    ptab(jf) % pt3d(ji, nlcj - 1, jk) = psgn(jf) * ptab2(jf) % pt4d(iju, ijpjp1, jk, jl)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        END SELECT
      CASE DEFAULT
        WRITE(*, *) 'lbc_nfd_nogather_generic: You should not have seen this print! error?', npolj
      END SELECT
    END DO
  END SUBROUTINE lbc_nfd_nogather_3d_ptr
  SUBROUTINE lbc_nfd_nogather_4d(ptab, ptab2, cd_nat, psgn, kfld)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :, :)
    REAL(KIND = wp), INTENT(INOUT) :: ptab2(:, :, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ijpj, ijpjp1, ijta, ijua, jia, startloop, endloop
    LOGICAL :: l_fast_exchanges
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
    CALL profile_psy_data0 % PreStart('lbc_nfd_nogather_4d', 'r0', 0, 0)
    ipj = SIZE(ptab2, 2)
    ipk = SIZE(ptab, 3)
    ipl = SIZE(ptab, 4)
    ipf = 1
    IF (ipf > 1) THEN
      WRITE(6, *) 'lbc_nfd_nogather: multiple fields not allowed. Revise implementation'
      WRITE(6, *) 'You should not be there...'
      STOP
    END IF
    ijpj = 1
    ijpjp1 = 2
    IF (ipj == 1) THEN
      l_fast_exchanges = .TRUE.
    ELSE
      l_fast_exchanges = .FALSE.
    END IF
    CALL profile_psy_data0 % PostEnd
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          CALL profile_psy_data1 % PreStart('lbc_nfd_nogather_4d', 'r1', 0, 0)
          IF (nimpp /= 1) THEN
            startloop = 1
          ELSE
            startloop = 2
          END IF
          CALL profile_psy_data1 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = startloop, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                ptab(ji, nlcj, jk, jl) = psgn * ptab2(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF (nimpp == 1) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                ptab(1, nlcj, jk, jl) = psgn * ptab(3, nlcj - 2, jk, jl)
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data2 % PreStart('lbc_nfd_nogather_4d', 'r2', 0, 0)
            IF (nimpp >= jpiglo / 2 + 1) THEN
              startloop = 1
            ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
              startloop = jpiglo / 2 + 1 - nimpp + 1
            ELSE
              startloop = nlci + 1
            END IF
            CALL profile_psy_data2 % PostEnd
            IF (startloop <= nlci) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, nlci
                    ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                    jia = ji + nimpp - 1
                    ijta = jpiglo - jia + 2
                    IF (ijta >= startloop + nimpp - 1 .AND. ijta < jia) THEN
                      ptab(ji, nlcj - 1, jk, jl) = psgn * ptab(ijta - nimpp + 1, nlcj - 1, jk, jl)
                    ELSE
                      ptab(ji, nlcj - 1, jk, jl) = psgn * ptab2(ijt, ijpjp1, jk, jl)
                    END IF
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        CASE ('U')
          CALL profile_psy_data3 % PreStart('lbc_nfd_nogather_4d', 'r3', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data3 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(ji, nlcj, jk, jl) = psgn * ptab2(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          !$ACC KERNELS
          IF (nimpp .EQ. 1) THEN
            ptab(1, nlcj, :, :) = psgn * ptab(2, nlcj - 2, :, :)
          END IF
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            ptab(nlci, nlcj, :, :) = psgn * ptab(nlci - 1, nlcj - 2, :, :)
          END IF
          !$ACC END KERNELS
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data4 % PreStart('lbc_nfd_nogather_4d', 'r4', 0, 0)
            IF (nimpp + nlci - 1 /= jpiglo) THEN
              endloop = nlci
            ELSE
              endloop = nlci - 1
            END IF
            IF (nimpp >= jpiglo / 2) THEN
              startloop = 1
            ELSE IF ((nimpp + nlci - 1 >= jpiglo / 2) .AND. (nimpp < jpiglo / 2)) THEN
              startloop = jpiglo / 2 - nimpp + 1
            ELSE
              startloop = endloop + 1
            END IF
            CALL profile_psy_data4 % PostEnd
            IF (startloop <= endloop) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, endloop
                    iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                    jia = ji + nimpp - 1
                    ijua = jpiglo - jia + 1
                    IF (ijua >= startloop + nimpp - 1 .AND. ijua < jia) THEN
                      ptab(ji, nlcj - 1, jk, jl) = psgn * ptab(ijua - nimpp + 1, nlcj - 1, jk, jl)
                    ELSE
                      ptab(ji, nlcj - 1, jk, jl) = psgn * ptab2(iju, ijpjp1, jk, jl)
                    END IF
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        CASE ('V')
          CALL profile_psy_data5 % PreStart('lbc_nfd_nogather_4d', 'r5', 0, 0)
          IF (nimpp /= 1) THEN
            startloop = 1
          ELSE
            startloop = 2
          END IF
          CALL profile_psy_data5 % PostEnd
          IF (.NOT. l_fast_exchanges) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                DO ji = startloop, nlci
                  ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                  ptab(ji, nlcj - 1, jk, jl) = psgn * ptab2(ijt, ijpjp1, jk, jl)
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = startloop, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
                ptab(ji, nlcj, jk, jl) = psgn * ptab2(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          !$ACC KERNELS
          IF (nimpp .EQ. 1) THEN
            ptab(1, nlcj, :, :) = psgn * ptab(3, nlcj - 3, :, :)
          END IF
          !$ACC END KERNELS
        CASE ('F')
          CALL profile_psy_data6 % PreStart('lbc_nfd_nogather_4d', 'r6', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data6 % PostEnd
          IF (.NOT. l_fast_exchanges) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                DO ji = 1, endloop
                  iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                  ptab(ji, nlcj - 1, jk, jl) = psgn * ptab2(iju, ijpjp1, jk, jl)
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(ji, nlcj, jk, jl) = psgn * ptab2(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF (nimpp .EQ. 1) THEN
            !$ACC KERNELS
            ptab(1, nlcj, :, :) = psgn * ptab(2, nlcj - 3, :, :)
            !$ACC END KERNELS
            CALL profile_psy_data7 % PreStart('lbc_nfd_nogather_4d', 'r7', 0, 0)
            IF (.NOT. l_fast_exchanges) ptab(1, nlcj - 1, :, :) = psgn * ptab(2, nlcj - 2, :, :)
            CALL profile_psy_data7 % PostEnd
          END IF
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            !$ACC KERNELS
            ptab(nlci, nlcj, :, :) = psgn * ptab(nlci - 1, nlcj - 3, :, :)
            !$ACC END KERNELS
            CALL profile_psy_data8 % PreStart('lbc_nfd_nogather_4d', 'r8', 0, 0)
            IF (.NOT. l_fast_exchanges) ptab(nlci, nlcj - 1, :, :) = psgn * ptab(nlci - 1, nlcj - 2, :, :)
            CALL profile_psy_data8 % PostEnd
          END IF
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(ji, nlcj, jk, jl) = psgn * ptab2(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        CASE ('U')
          CALL profile_psy_data9 % PreStart('lbc_nfd_nogather_4d', 'r9', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data9 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
                ptab(ji, nlcj, jk, jl) = psgn * ptab2(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                ptab(nlci, nlcj, jk, jl) = psgn * ptab(nlci - 2, nlcj - 1, jk, jl)
              END DO
              !$ACC END KERNELS
            END DO
          END IF
        CASE ('V')
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, nlci
                ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                ptab(ji, nlcj, jk, jl) = psgn * ptab2(ijt, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data10 % PreStart('lbc_nfd_nogather_4d', 'r10', 0, 0)
            IF (nimpp >= jpiglo / 2 + 1) THEN
              startloop = 1
            ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
              startloop = jpiglo / 2 + 1 - nimpp + 1
            ELSE
              startloop = nlci + 1
            END IF
            CALL profile_psy_data10 % PostEnd
            IF (startloop <= nlci) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, nlci
                    ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
                    ptab(ji, nlcj - 1, jk, jl) = psgn * ptab2(ijt, ijpjp1, jk, jl)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        CASE ('F')
          CALL profile_psy_data11 % PreStart('lbc_nfd_nogather_4d', 'r11', 0, 0)
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          CALL profile_psy_data11 % PostEnd
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              DO ji = 1, endloop
                iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
                ptab(ji, nlcj, jk, jl) = psgn * ptab2(iju, ijpj, jk, jl)
              END DO
            END DO
            !$ACC END KERNELS
          END DO
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                ptab(nlci, nlcj, jk, jl) = psgn * ptab(nlci - 2, nlcj - 2, jk, jl)
              END DO
              !$ACC END KERNELS
            END DO
          END IF
          IF (.NOT. l_fast_exchanges) THEN
            CALL profile_psy_data12 % PreStart('lbc_nfd_nogather_4d', 'r12', 0, 0)
            IF (nimpp + nlci - 1 /= jpiglo) THEN
              endloop = nlci
            ELSE
              endloop = nlci - 1
            END IF
            IF (nimpp >= jpiglo / 2 + 1) THEN
              startloop = 1
            ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
              startloop = jpiglo / 2 + 1 - nimpp + 1
            ELSE
              startloop = endloop + 1
            END IF
            CALL profile_psy_data12 % PostEnd
            IF (startloop <= endloop) THEN
              DO jl = 1, ipl
                !$ACC KERNELS
                DO jk = 1, ipk
                  DO ji = startloop, endloop
                    iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
                    ptab(ji, nlcj - 1, jk, jl) = psgn * ptab2(iju, ijpjp1, jk, jl)
                  END DO
                END DO
                !$ACC END KERNELS
              END DO
            END IF
          END IF
        END SELECT
      CASE DEFAULT
        WRITE(*, *) 'lbc_nfd_nogather_generic: You should not have seen this print! error?', npolj
      END SELECT
    END DO
  END SUBROUTINE lbc_nfd_nogather_4d
END MODULE lbcnfd
