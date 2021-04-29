MODULE usrdef_istate
  USE par_oce
  USE phycst
  USE in_out_manager
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: usr_def_istate
  CONTAINS
  SUBROUTINE usr_def_istate(pdept, ptmask, pts, pu, pv, pssh)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pdept
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: ptmask
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, jpts), INTENT(OUT) :: pts
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(OUT) :: pu
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(OUT) :: pv
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(OUT) :: pssh
    INTEGER :: ji, jj, jk
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('usr_def_istate', 'r0', 0, 0)
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) 'usr_def_istate : analytical definition of initial state '
    IF (lwp) WRITE(numout, *) '~~~~~~~~~~~~~~   Ocean at rest, with an horizontally uniform T and S profiles'
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    pu(:, :, :) = 0._wp
    pv(:, :, :) = 0._wp
    pssh(:, :) = 0._wp
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('usr_def_istate', 'r1', 0, 0)
    DO jk = 1, jpk
      DO jj = 1, jpj
        DO ji = 1, jpi
          pts(ji, jj, jk, jp_tem) = ((16. - 12. * TANH((pdept(ji, jj, jk) - 400) / 700)) * (- TANH((500. - pdept(ji, jj, jk)) / 150.) + 1.) / 2. + (15. * (1. - TANH((pdept(ji, jj, jk) - 50.) / 1500.)) - 1.4 * TANH((pdept(ji, jj, jk) - 100.) / 100.) + 7. * (1500. - pdept(ji, jj, jk)) / 1500.) * (- TANH((pdept(ji, jj, jk) - 500.) / 150.) + 1.) / 2.) * ptmask(ji, jj, jk)
          pts(ji, jj, jk, jp_sal) = ((36.25 - 1.13 * TANH((pdept(ji, jj, jk) - 305) / 460)) * (- TANH((500. - pdept(ji, jj, jk)) / 150.) + 1.) / 2 + (35.55 + 1.25 * (5000. - pdept(ji, jj, jk)) / 5000. - 1.62 * TANH((pdept(ji, jj, jk) - 60.) / 650.) + 0.2 * TANH((pdept(ji, jj, jk) - 35.) / 100.) + 0.2 * TANH((pdept(ji, jj, jk) - 1000.) / 5000.)) * (- TANH((pdept(ji, jj, jk) - 500.) / 150.) + 1.) / 2) * ptmask(ji, jj, jk)
        END DO
      END DO
    END DO
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE usr_def_istate
END MODULE usrdef_istate