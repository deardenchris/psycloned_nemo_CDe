MODULE storng
  USE par_kind
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: kiss, kiss_seed, kiss_state, kiss_reset
  PUBLIC :: kiss_uniform, kiss_gaussian, kiss_gamma, kiss_sample
  INTEGER(KIND = i8) :: x = 1234567890987654321_8
  INTEGER(KIND = i8) :: y = 362436362436362436_8
  INTEGER(KIND = i8) :: z = 1066149217761810_8
  INTEGER(KIND = i8) :: w = 123456123456123456_8
  REAL(KIND = wp), PARAMETER :: huge64 = 9223372036854775808.0
  REAL(KIND = wp), PARAMETER :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0
  INTEGER(KIND = i8), SAVE :: ig = 1
  REAL(KIND = wp), SAVE :: gran1, gran2
  CONTAINS
  FUNCTION kiss()
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    IMPLICIT NONE
    INTEGER(KIND = i8) :: kiss, t
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('kiss', 'r0', 0, 0)
    t = ISHFT(x, 58) + w
    IF (s(x) .EQ. s(t)) THEN
      w = ISHFT(x, - 6) + s(x)
    ELSE
      w = ISHFT(x, - 6) + 1 - s(x + t)
    END IF
    x = t + x
    y = m(m(m(y, 13_8), - 17_8), 43_8)
    z = 6906969069_8 * z + 1234567_8
    kiss = x + y + z
    CALL profile_psy_data0 % PostEnd
    CONTAINS
    FUNCTION s(k)
      USE profile_psy_data_mod, ONLY: profile_PSyDataType
      INTEGER(KIND = i8) :: s, k
      TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
      CALL profile_psy_data0 % PreStart('s', 'r0', 0, 0)
      s = ISHFT(k, - 63)
      CALL profile_psy_data0 % PostEnd
    END FUNCTION s
    FUNCTION m(k, n)
      USE profile_psy_data_mod, ONLY: profile_PSyDataType
      INTEGER(KIND = i8) :: m, k, n
      TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
      CALL profile_psy_data0 % PreStart('m', 'r0', 0, 0)
      m = IEOR(k, ISHFT(k, n))
      CALL profile_psy_data0 % PostEnd
    END FUNCTION m
  END FUNCTION kiss
  SUBROUTINE kiss_seed(ix, iy, iz, iw)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    IMPLICIT NONE
    INTEGER(KIND = i8) :: ix, iy, iz, iw
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('kiss_seed', 'r0', 0, 0)
    x = ix
    y = iy
    z = iz
    w = iw
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE kiss_seed
  SUBROUTINE kiss_state(ix, iy, iz, iw)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    IMPLICIT NONE
    INTEGER(KIND = i8) :: ix, iy, iz, iw
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('kiss_state', 'r0', 0, 0)
    ix = x
    iy = y
    iz = z
    iw = w
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE kiss_state
  SUBROUTINE kiss_reset
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    IMPLICIT NONE
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('kiss_reset', 'r0', 0, 0)
    x = 1234567890987654321_8
    y = 362436362436362436_8
    z = 1066149217761810_8
    w = 123456123456123456_8
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE kiss_reset
  SUBROUTINE kiss_uniform(uran)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    IMPLICIT NONE
    REAL(KIND = wp) :: uran
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('kiss_uniform', 'r0', 0, 0)
    uran = half * (one + REAL(kiss(), wp) / huge64)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE kiss_uniform
  SUBROUTINE kiss_gaussian(gran)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    IMPLICIT NONE
    REAL(KIND = wp) :: gran, u1, u2, rsq, fac
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('kiss_gaussian', 'r0', 0, 0)
    IF (ig .EQ. 1) THEN
      rsq = two
      DO WHILE ((rsq .GE. one) .OR. (rsq .EQ. zero))
        u1 = REAL(kiss(), wp) / huge64
        u2 = REAL(kiss(), wp) / huge64
        rsq = u1 * u1 + u2 * u2
      END DO
      fac = SQRT(- two * LOG(rsq) / rsq)
      gran1 = u1 * fac
      gran2 = u2 * fac
    END IF
    IF (ig .EQ. 1) THEN
      gran = gran1
      ig = 2
    ELSE
      gran = gran2
      ig = 1
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE kiss_gaussian
  SUBROUTINE kiss_gamma(gamr, k)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    IMPLICIT NONE
    REAL(KIND = wp), PARAMETER :: p1 = 4.5_8
    REAL(KIND = wp), PARAMETER :: p2 = 2.50407739677627_8
    REAL(KIND = wp), PARAMETER :: p3 = 1.38629436111989_8
    REAL(KIND = wp) :: gamr, k, u1, u2, b, c, d, xx, yy, zz, rr, ee
    LOGICAL :: accepted
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('kiss_gamma', 'r0', 0, 0)
    IF (k .GT. one) THEN
      b = k - p3
      d = SQRT(two * k - one)
      c = k + d
      accepted = .FALSE.
      DO WHILE (.NOT. accepted)
        CALL kiss_uniform(u1)
        yy = LOG(u1 / (one - u1)) / d
        xx = k * EXP(yy)
        rr = b + c * yy - xx
        CALL kiss_uniform(u2)
        zz = u1 * u1 * u2
        accepted = rr .GE. (zz * p1 - p2)
        IF (.NOT. accepted) accepted = rr .GE. LOG(zz)
      END DO
      gamr = xx
    ELSE IF (k .LT. one) THEN
      c = one / k
      d = (one - k) * EXP((k / (one - k)) * LOG(k))
      accepted = .FALSE.
      DO WHILE (.NOT. accepted)
        CALL kiss_uniform(u1)
        zz = - LOG(u1)
        xx = EXP(c * LOG(zz))
        CALL kiss_uniform(u2)
        ee = - LOG(u2)
        accepted = (zz + ee) .GE. (d + xx)
      END DO
      gamr = xx
    ELSE
      CALL kiss_uniform(u1)
      gamr = - LOG(u1)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE kiss_gamma
  SUBROUTINE kiss_sample(a, n, k)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    IMPLICIT NONE
    INTEGER(KIND = i8), DIMENSION(:) :: a
    INTEGER(KIND = i8) :: n, k, i, j, atmp
    REAL(KIND = wp) :: uran
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('kiss_sample', 'r0', 0, 0)
    DO i = 1, k
      CALL kiss_uniform(uran)
      j = i - 1 + CEILING(REAL(n - i + 1, 8) * uran)
      atmp = a(i)
      a(i) = a(j)
      a(j) = atmp
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE kiss_sample
END MODULE storng