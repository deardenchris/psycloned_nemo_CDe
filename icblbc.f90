MODULE icblbc
  USE par_oce
  USE dom_oce
  USE in_out_manager
  USE lib_mpp
  USE icb_oce
  USE icbutl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: icb_lbc
  PUBLIC :: icb_lbc_mpp
  CONTAINS
  SUBROUTINE icb_lbc
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(iceberg), POINTER :: this
    TYPE(point), POINTER :: pt
    INTEGER :: iine
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_lbc', 'r0', 0, 0)
    IF (l_Iperio) THEN
      this => first_berg
      DO WHILE (ASSOCIATED(this))
        pt => this % current_point
        iine = INT(pt % xi + 0.5)
        IF (iine > mig(nicbei)) THEN
          pt % xi = ricb_right + MOD(pt % xi, 1._wp) - 1._wp
        ELSE IF (iine < mig(nicbdi)) THEN
          pt % xi = ricb_left + MOD(pt % xi, 1._wp)
        END IF
        this => this % next
      END DO
    END IF
    IF (l_Jperio) CALL ctl_stop(' north-south periodicity not implemented for icebergs')
    IF (npolj /= 0) CALL icb_lbc_nfld
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_lbc
  SUBROUTINE icb_lbc_nfld
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(iceberg), POINTER :: this
    TYPE(point), POINTER :: pt
    INTEGER :: iine, ijne, ipts
    INTEGER :: iiglo, ijglo
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_lbc_nfld', 'r0', 0, 0)
    this => first_berg
    DO WHILE (ASSOCIATED(this))
      pt => this % current_point
      ijne = INT(pt % yj + 0.5)
      IF (ijne .GT. mjg(nicbej)) THEN
        iine = INT(pt % xi + 0.5)
        ipts = nicbfldpts(mi1(iine))
        ijglo = INT(ipts / nicbpack)
        iiglo = ipts - nicbpack * ijglo
        pt % xi = iiglo - (pt % xi - REAL(iine, wp))
        pt % yj = ijglo - (pt % yj - REAL(ijne, wp))
        pt % uvel = - 1._wp * pt % uvel
        pt % vvel = - 1._wp * pt % vvel
      END IF
      this => this % next
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_lbc_nfld
  SUBROUTINE icb_lbc_mpp
    WRITE(numout, *) 'icb_lbc_mpp: You should not have seen this message!!'
  END SUBROUTINE icb_lbc_mpp
END MODULE icblbc