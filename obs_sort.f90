MODULE obs_sort
  USE par_kind, ONLY: dp
  IMPLICIT NONE
  PRIVATE :: index_sort
  PUBLIC :: sort_dp_indx
  CONTAINS
  SUBROUTINE sort_dp_indx(kvals, pvals, kindx)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kvals
    REAL(KIND = dp), DIMENSION(kvals), INTENT(IN) :: pvals
    INTEGER, DIMENSION(kvals), INTENT(OUT) :: kindx
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('sort_dp_indx', 'r0', 0, 0)
    IF (kvals >= 1) THEN
      CALL index_sort(pvals, kindx, kvals)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE sort_dp_indx
  SUBROUTINE index_sort(pval, kindx, kvals)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kvals
    REAL(KIND = dp), DIMENSION(kvals), INTENT(IN) :: pval
    INTEGER, DIMENSION(kvals), INTENT(INOUT) :: kindx
    INTEGER :: ji
    INTEGER :: jj
    INTEGER :: jt
    INTEGER :: jn
    INTEGER :: jparent
    INTEGER :: jchild
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('index_sort', 'r0', 0, 0)
    DO ji = 1, kvals
      kindx(ji) = ji
    END DO
    ji = kvals / 2 + 1
    jn = kvals
    main_loop:DO
      IF (ji > 1) THEN
        ji = ji - 1
        jt = kindx(ji)
      ELSE
        jt = kindx(jn)
        kindx(jn) = kindx(1)
        jn = jn - 1
        IF (jn <= 1) THEN
          kindx(1) = jt
          EXIT main_loop
        END IF
      END IF
      jparent = ji
      jchild = 2 * ji
      inner_loop:DO
        IF (jchild > jn) EXIT inner_loop
        IF (jchild < jn) THEN
          IF (pval(kindx(jchild)) < pval(kindx(jchild + 1))) THEN
            jchild = jchild + 1
          END IF
        END IF
        IF (pval(jt) < pval(kindx(jchild))) THEN
          kindx(jparent) = kindx(jchild)
          jparent = jchild
          jchild = jchild * 2
        ELSE
          jchild = jn + 1
        END IF
      END DO inner_loop
      kindx(jparent) = jt
    END DO main_loop
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE index_sort
END MODULE obs_sort