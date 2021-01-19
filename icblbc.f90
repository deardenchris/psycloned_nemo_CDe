MODULE icblbc
  USE par_oce
  USE dom_oce
  USE in_out_manager
  USE lib_mpp
  USE icb_oce
  USE icbutl
  IMPLICIT NONE
  PRIVATE
  INCLUDE 'mpif.h'
  TYPE, PUBLIC :: buffer
    INTEGER :: size = 0
    REAL(KIND = wp), DIMENSION(:, :), POINTER :: data
  END TYPE buffer
  TYPE(buffer), POINTER :: obuffer_n => NULL(), ibuffer_n => NULL()
  TYPE(buffer), POINTER :: obuffer_s => NULL(), ibuffer_s => NULL()
  TYPE(buffer), POINTER :: obuffer_e => NULL(), ibuffer_e => NULL()
  TYPE(buffer), POINTER :: obuffer_w => NULL(), ibuffer_w => NULL()
  TYPE(buffer), POINTER :: obuffer_f => NULL(), ibuffer_f => NULL()
  INTEGER, PARAMETER, PRIVATE :: jp_delta_buf = 25
  INTEGER, PARAMETER, PRIVATE :: jp_buffer_width = 15 + nkounts
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
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(iceberg), POINTER :: tmpberg, this
    TYPE(point), POINTER :: pt
    INTEGER :: ibergs_to_send_e, ibergs_to_send_w
    INTEGER :: ibergs_to_send_n, ibergs_to_send_s
    INTEGER :: ibergs_rcvd_from_e, ibergs_rcvd_from_w
    INTEGER :: ibergs_rcvd_from_n, ibergs_rcvd_from_s
    INTEGER :: i, ibergs_start, ibergs_end
    INTEGER :: iine, ijne
    INTEGER :: ipe_N, ipe_S, ipe_W, ipe_E
    REAL(KIND = wp), DIMENSION(2) :: zewbergs, zwebergs, znsbergs, zsnbergs
    INTEGER :: iml_req1, iml_req2, iml_req3, iml_req4
    INTEGER :: iml_req5, iml_req6, iml_req7, iml_req8, iml_err
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: iml_stat
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_lbc_mpp', 'r0', 0, 0)
    ipe_N = - 1
    ipe_S = - 1
    ipe_W = - 1
    ipe_E = - 1
    IF (nbondi .EQ. 0 .OR. nbondi .EQ. 1) ipe_W = nowe
    IF (nbondi .EQ. - 1 .OR. nbondi .EQ. 0) ipe_E = noea
    IF (nbondj .EQ. 0 .OR. nbondj .EQ. 1) ipe_S = noso
    IF (nbondj .EQ. - 1 .OR. nbondj .EQ. 0) ipe_N = nono
    IF (npolj > 0) ipe_N = - 1
    IF (jpni == 1) THEN
      ipe_E = - 1
      ipe_W = - 1
    END IF
    IF (nn_verbose_level >= 2) THEN
      WRITE(numicb, FMT = *) 'processor west  : ', ipe_W
      WRITE(numicb, FMT = *) 'processor east  : ', ipe_E
      WRITE(numicb, FMT = *) 'processor north : ', ipe_N
      WRITE(numicb, FMT = *) 'processor south : ', ipe_S
      WRITE(numicb, FMT = *) 'processor nimpp : ', nimpp
      WRITE(numicb, FMT = *) 'processor njmpp : ', njmpp
      WRITE(numicb, FMT = *) 'processor nbondi: ', nbondi
      WRITE(numicb, FMT = *) 'processor nbondj: ', nbondj
      CALL flush(numicb)
    END IF
    IF (jpni == 1) CALL icb_lbc
    IF (nn_verbose_level > 0) THEN
      ibergs_start = icb_utl_count()
    END IF
    ibergs_to_send_e = 0
    ibergs_to_send_w = 0
    ibergs_to_send_n = 0
    ibergs_to_send_s = 0
    ibergs_rcvd_from_e = 0
    ibergs_rcvd_from_w = 0
    ibergs_rcvd_from_n = 0
    ibergs_rcvd_from_s = 0
    IF (ASSOCIATED(first_berg)) THEN
      this => first_berg
      DO WHILE (ASSOCIATED(this))
        pt => this % current_point
        iine = INT(pt % xi + 0.5)
        IF (ipe_E >= 0 .AND. iine > mig(nicbei)) THEN
          tmpberg => this
          this => this % next
          ibergs_to_send_e = ibergs_to_send_e + 1
          IF (nn_verbose_level >= 4) THEN
            WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' packing berg ', tmpberg % number(:), ' for transfer to east'
            CALL flush(numicb)
          END IF
          tmpberg % current_point % xi = ricb_right + MOD(tmpberg % current_point % xi, 1._wp) - 1._wp
          CALL icb_pack_into_buffer(tmpberg, obuffer_e, ibergs_to_send_e)
          CALL icb_utl_delete(first_berg, tmpberg)
        ELSE IF (ipe_W >= 0 .AND. iine < mig(nicbdi)) THEN
          tmpberg => this
          this => this % next
          ibergs_to_send_w = ibergs_to_send_w + 1
          IF (nn_verbose_level >= 4) THEN
            WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' packing berg ', tmpberg % number(:), ' for transfer to west'
            CALL flush(numicb)
          END IF
          tmpberg % current_point % xi = ricb_left + MOD(tmpberg % current_point % xi, 1._wp)
          CALL icb_pack_into_buffer(tmpberg, obuffer_w, ibergs_to_send_w)
          CALL icb_utl_delete(first_berg, tmpberg)
        ELSE
          this => this % next
        END IF
      END DO
    END IF
    IF (nn_verbose_level >= 3) THEN
      WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' send ew: ', ibergs_to_send_e, ibergs_to_send_w
      CALL flush(numicb)
    END IF
    SELECT CASE (nbondi)
    CASE (- 1)
      zwebergs(1) = ibergs_to_send_e
      CALL mppsend(12, zwebergs(1), 1, ipe_E, iml_req1)
      CALL mpprecv(11, zewbergs(2), 1, ipe_E)
      IF (l_isend) CALL mpi_wait(iml_req1, iml_stat, iml_err)
      ibergs_rcvd_from_e = INT(zewbergs(2))
    CASE (0)
      zewbergs(1) = ibergs_to_send_w
      zwebergs(1) = ibergs_to_send_e
      CALL mppsend(11, zewbergs(1), 1, ipe_W, iml_req2)
      CALL mppsend(12, zwebergs(1), 1, ipe_E, iml_req3)
      CALL mpprecv(11, zewbergs(2), 1, ipe_E)
      CALL mpprecv(12, zwebergs(2), 1, ipe_W)
      IF (l_isend) CALL mpi_wait(iml_req2, iml_stat, iml_err)
      IF (l_isend) CALL mpi_wait(iml_req3, iml_stat, iml_err)
      ibergs_rcvd_from_e = INT(zewbergs(2))
      ibergs_rcvd_from_w = INT(zwebergs(2))
    CASE (1)
      zewbergs(1) = ibergs_to_send_w
      CALL mppsend(11, zewbergs(1), 1, ipe_W, iml_req4)
      CALL mpprecv(12, zwebergs(2), 1, ipe_W)
      IF (l_isend) CALL mpi_wait(iml_req4, iml_stat, iml_err)
      ibergs_rcvd_from_w = INT(zwebergs(2))
    END SELECT
    IF (nn_verbose_level >= 3) THEN
      WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' recv ew: ', ibergs_rcvd_from_w, ibergs_rcvd_from_e
      CALL flush(numicb)
    END IF
    SELECT CASE (nbondi)
    CASE (- 1)
      IF (ibergs_to_send_e > 0) CALL mppsend(14, obuffer_e % data, ibergs_to_send_e * jp_buffer_width, ipe_E, iml_req1)
      IF (ibergs_rcvd_from_e > 0) THEN
        CALL icb_increase_ibuffer(ibuffer_e, ibergs_rcvd_from_e)
        CALL mpprecv(13, ibuffer_e % data, ibergs_rcvd_from_e * jp_buffer_width)
      END IF
      IF (ibergs_to_send_e > 0 .AND. l_isend) CALL mpi_wait(iml_req1, iml_stat, iml_err)
      DO i = 1, ibergs_rcvd_from_e
        IF (nn_verbose_level >= 4) THEN
          WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' unpacking berg ', INT(ibuffer_e % data(16, i)), ' from east'
          CALL flush(numicb)
        END IF
        CALL icb_unpack_from_buffer(first_berg, ibuffer_e, i)
      END DO
    CASE (0)
      IF (ibergs_to_send_w > 0) CALL mppsend(13, obuffer_w % data, ibergs_to_send_w * jp_buffer_width, ipe_W, iml_req2)
      IF (ibergs_to_send_e > 0) CALL mppsend(14, obuffer_e % data, ibergs_to_send_e * jp_buffer_width, ipe_E, iml_req3)
      IF (ibergs_rcvd_from_e > 0) THEN
        CALL icb_increase_ibuffer(ibuffer_e, ibergs_rcvd_from_e)
        CALL mpprecv(13, ibuffer_e % data, ibergs_rcvd_from_e * jp_buffer_width)
      END IF
      IF (ibergs_rcvd_from_w > 0) THEN
        CALL icb_increase_ibuffer(ibuffer_w, ibergs_rcvd_from_w)
        CALL mpprecv(14, ibuffer_w % data, ibergs_rcvd_from_w * jp_buffer_width)
      END IF
      IF (ibergs_to_send_w > 0 .AND. l_isend) CALL mpi_wait(iml_req2, iml_stat, iml_err)
      IF (ibergs_to_send_e > 0 .AND. l_isend) CALL mpi_wait(iml_req3, iml_stat, iml_err)
      DO i = 1, ibergs_rcvd_from_e
        IF (nn_verbose_level >= 4) THEN
          WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' unpacking berg ', INT(ibuffer_e % data(16, i)), ' from east'
          CALL flush(numicb)
        END IF
        CALL icb_unpack_from_buffer(first_berg, ibuffer_e, i)
      END DO
      DO i = 1, ibergs_rcvd_from_w
        IF (nn_verbose_level >= 4) THEN
          WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' unpacking berg ', INT(ibuffer_w % data(16, i)), ' from west'
          CALL flush(numicb)
        END IF
        CALL icb_unpack_from_buffer(first_berg, ibuffer_w, i)
      END DO
    CASE (1)
      IF (ibergs_to_send_w > 0) CALL mppsend(13, obuffer_w % data, ibergs_to_send_w * jp_buffer_width, ipe_W, iml_req4)
      IF (ibergs_rcvd_from_w > 0) THEN
        CALL icb_increase_ibuffer(ibuffer_w, ibergs_rcvd_from_w)
        CALL mpprecv(14, ibuffer_w % data, ibergs_rcvd_from_w * jp_buffer_width)
      END IF
      IF (ibergs_to_send_w > 0 .AND. l_isend) CALL mpi_wait(iml_req4, iml_stat, iml_err)
      DO i = 1, ibergs_rcvd_from_w
        IF (nn_verbose_level >= 4) THEN
          WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' unpacking berg ', INT(ibuffer_w % data(16, i)), ' from west'
          CALL flush(numicb)
        END IF
        CALL icb_unpack_from_buffer(first_berg, ibuffer_w, i)
      END DO
    END SELECT
    IF (ASSOCIATED(first_berg)) THEN
      this => first_berg
      DO WHILE (ASSOCIATED(this))
        pt => this % current_point
        ijne = INT(pt % yj + 0.5)
        IF (ipe_N >= 0 .AND. ijne .GT. mjg(nicbej)) THEN
          tmpberg => this
          this => this % next
          ibergs_to_send_n = ibergs_to_send_n + 1
          IF (nn_verbose_level >= 4) THEN
            WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' packing berg ', tmpberg % number(:), ' for transfer to north'
            CALL flush(numicb)
          END IF
          CALL icb_pack_into_buffer(tmpberg, obuffer_n, ibergs_to_send_n)
          CALL icb_utl_delete(first_berg, tmpberg)
        ELSE IF (ipe_S >= 0 .AND. ijne .LT. mjg(nicbdj)) THEN
          tmpberg => this
          this => this % next
          ibergs_to_send_s = ibergs_to_send_s + 1
          IF (nn_verbose_level >= 4) THEN
            WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' packing berg ', tmpberg % number(:), ' for transfer to south'
            CALL flush(numicb)
          END IF
          CALL icb_pack_into_buffer(tmpberg, obuffer_s, ibergs_to_send_s)
          CALL icb_utl_delete(first_berg, tmpberg)
        ELSE
          this => this % next
        END IF
      END DO
    END IF
    IF (nn_verbose_level >= 3) THEN
      WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' send ns: ', ibergs_to_send_n, ibergs_to_send_s
      CALL flush(numicb)
    END IF
    SELECT CASE (nbondj)
    CASE (- 1)
      zsnbergs(1) = ibergs_to_send_n
      CALL mppsend(16, zsnbergs(1), 1, ipe_N, iml_req1)
      CALL mpprecv(15, znsbergs(2), 1, ipe_N)
      IF (l_isend) CALL mpi_wait(iml_req1, iml_stat, iml_err)
      ibergs_rcvd_from_n = INT(znsbergs(2))
    CASE (0)
      znsbergs(1) = ibergs_to_send_s
      zsnbergs(1) = ibergs_to_send_n
      CALL mppsend(15, znsbergs(1), 1, ipe_S, iml_req2)
      CALL mppsend(16, zsnbergs(1), 1, ipe_N, iml_req3)
      CALL mpprecv(15, znsbergs(2), 1, ipe_N)
      CALL mpprecv(16, zsnbergs(2), 1, ipe_S)
      IF (l_isend) CALL mpi_wait(iml_req2, iml_stat, iml_err)
      IF (l_isend) CALL mpi_wait(iml_req3, iml_stat, iml_err)
      ibergs_rcvd_from_n = INT(znsbergs(2))
      ibergs_rcvd_from_s = INT(zsnbergs(2))
    CASE (1)
      znsbergs(1) = ibergs_to_send_s
      CALL mppsend(15, znsbergs(1), 1, ipe_S, iml_req4)
      CALL mpprecv(16, zsnbergs(2), 1, ipe_S)
      IF (l_isend) CALL mpi_wait(iml_req4, iml_stat, iml_err)
      ibergs_rcvd_from_s = INT(zsnbergs(2))
    END SELECT
    IF (nn_verbose_level >= 3) THEN
      WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' recv ns: ', ibergs_rcvd_from_s, ibergs_rcvd_from_n
      CALL flush(numicb)
    END IF
    SELECT CASE (nbondj)
    CASE (- 1)
      IF (ibergs_to_send_n > 0) CALL mppsend(18, obuffer_n % data, ibergs_to_send_n * jp_buffer_width, ipe_N, iml_req1)
      IF (ibergs_rcvd_from_n > 0) THEN
        CALL icb_increase_ibuffer(ibuffer_n, ibergs_rcvd_from_n)
        CALL mpprecv(17, ibuffer_n % data, ibergs_rcvd_from_n * jp_buffer_width)
      END IF
      IF (ibergs_to_send_n > 0 .AND. l_isend) CALL mpi_wait(iml_req1, iml_stat, iml_err)
      DO i = 1, ibergs_rcvd_from_n
        IF (nn_verbose_level >= 4) THEN
          WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' unpacking berg ', INT(ibuffer_n % data(16, i)), ' from north'
          CALL flush(numicb)
        END IF
        CALL icb_unpack_from_buffer(first_berg, ibuffer_n, i)
      END DO
    CASE (0)
      IF (ibergs_to_send_s > 0) CALL mppsend(17, obuffer_s % data, ibergs_to_send_s * jp_buffer_width, ipe_S, iml_req2)
      IF (ibergs_to_send_n > 0) CALL mppsend(18, obuffer_n % data, ibergs_to_send_n * jp_buffer_width, ipe_N, iml_req3)
      IF (ibergs_rcvd_from_n > 0) THEN
        CALL icb_increase_ibuffer(ibuffer_n, ibergs_rcvd_from_n)
        CALL mpprecv(17, ibuffer_n % data, ibergs_rcvd_from_n * jp_buffer_width)
      END IF
      IF (ibergs_rcvd_from_s > 0) THEN
        CALL icb_increase_ibuffer(ibuffer_s, ibergs_rcvd_from_s)
        CALL mpprecv(18, ibuffer_s % data, ibergs_rcvd_from_s * jp_buffer_width)
      END IF
      IF (ibergs_to_send_s > 0 .AND. l_isend) CALL mpi_wait(iml_req2, iml_stat, iml_err)
      IF (ibergs_to_send_n > 0 .AND. l_isend) CALL mpi_wait(iml_req3, iml_stat, iml_err)
      DO i = 1, ibergs_rcvd_from_n
        IF (nn_verbose_level >= 4) THEN
          WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' unpacking berg ', INT(ibuffer_n % data(16, i)), ' from north'
          CALL flush(numicb)
        END IF
        CALL icb_unpack_from_buffer(first_berg, ibuffer_n, i)
      END DO
      DO i = 1, ibergs_rcvd_from_s
        IF (nn_verbose_level >= 4) THEN
          WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' unpacking berg ', INT(ibuffer_s % data(16, i)), ' from south'
          CALL flush(numicb)
        END IF
        CALL icb_unpack_from_buffer(first_berg, ibuffer_s, i)
      END DO
    CASE (1)
      IF (ibergs_to_send_s > 0) CALL mppsend(17, obuffer_s % data, ibergs_to_send_s * jp_buffer_width, ipe_S, iml_req4)
      IF (ibergs_rcvd_from_s > 0) THEN
        CALL icb_increase_ibuffer(ibuffer_s, ibergs_rcvd_from_s)
        CALL mpprecv(18, ibuffer_s % data, ibergs_rcvd_from_s * jp_buffer_width)
      END IF
      IF (ibergs_to_send_s > 0 .AND. l_isend) CALL mpi_wait(iml_req4, iml_stat, iml_err)
      DO i = 1, ibergs_rcvd_from_s
        IF (nn_verbose_level >= 4) THEN
          WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' unpacking berg ', INT(ibuffer_s % data(16, i)), ' from south'
          CALL flush(numicb)
        END IF
        CALL icb_unpack_from_buffer(first_berg, ibuffer_s, i)
      END DO
    END SELECT
    IF (nn_verbose_level > 0) THEN
      ibergs_end = icb_utl_count()
      i = (ibergs_rcvd_from_n + ibergs_rcvd_from_s + ibergs_rcvd_from_e + ibergs_rcvd_from_w) - (ibergs_to_send_n + &
&ibergs_to_send_s + ibergs_to_send_e + ibergs_to_send_w)
      IF (ibergs_end - (ibergs_start + i) .NE. 0) THEN
        WRITE(numicb, FMT = *) 'send_bergs_to_other_pes: net change in number of icebergs'
        WRITE(numicb, 1000) 'send_bergs_to_other_pes: ibergs_end=', ibergs_end, ' on PE', narea
        WRITE(numicb, 1000) 'send_bergs_to_other_pes: ibergs_start=', ibergs_start, ' on PE', narea
        WRITE(numicb, 1000) 'send_bergs_to_other_pes: delta=', i, ' on PE', narea
        WRITE(numicb, 1000) 'send_bergs_to_other_pes: error=', ibergs_end - (ibergs_start + i), ' on PE', narea
        WRITE(numicb, 1000) 'send_bergs_to_other_pes: ibergs_to_send_n=', ibergs_to_send_n, ' on PE', narea
        WRITE(numicb, 1000) 'send_bergs_to_other_pes: ibergs_to_send_s=', ibergs_to_send_s, ' on PE', narea
        WRITE(numicb, 1000) 'send_bergs_to_other_pes: ibergs_to_send_e=', ibergs_to_send_e, ' on PE', narea
        WRITE(numicb, 1000) 'send_bergs_to_other_pes: ibergs_to_send_w=', ibergs_to_send_w, ' on PE', narea
        WRITE(numicb, 1000) 'send_bergs_to_other_pes: ibergs_rcvd_from_n=', ibergs_rcvd_from_n, ' on PE', narea
        WRITE(numicb, 1000) 'send_bergs_to_other_pes: ibergs_rcvd_from_s=', ibergs_rcvd_from_s, ' on PE', narea
        WRITE(numicb, 1000) 'send_bergs_to_other_pes: ibergs_rcvd_from_e=', ibergs_rcvd_from_e, ' on PE', narea
        WRITE(numicb, 1000) 'send_bergs_to_other_pes: ibergs_rcvd_from_w=', ibergs_rcvd_from_w, ' on PE', narea
1000    FORMAT(A, I5, A, I4)
        CALL ctl_stop('send_bergs_to_other_pes: lost or gained an iceberg or two')
      END IF
    END IF
    IF (npolj /= 0 .AND. jpni > 1) CALL icb_lbc_mpp_nfld
    IF (nn_verbose_level > 0) THEN
      i = 0
      this => first_berg
      DO WHILE (ASSOCIATED(this))
        pt => this % current_point
        iine = INT(pt % xi + 0.5)
        ijne = INT(pt % yj + 0.5)
        IF (iine .LT. mig(nicbdi) .OR. iine .GT. mig(nicbei) .OR. ijne .LT. mjg(nicbdj) .OR. ijne .GT. mjg(nicbej)) THEN
          i = i + 1
          WRITE(numicb, FMT = *) 'berg lost in halo: ', this % number(:), iine, ijne
          WRITE(numicb, FMT = *) '                   ', nimpp, njmpp
          WRITE(numicb, FMT = *) '                   ', nicbdi, nicbei, nicbdj, nicbej
          CALL flush(numicb)
        END IF
        this => this % next
      END DO
      CALL mpp_sum('icblbc', i)
      IF (i .GT. 0) THEN
        WRITE(numicb, FMT = '(a,i4)') 'send_bergs_to_other_pes: # of bergs outside computational domain = ', i
        CALL ctl_stop('send_bergs_to_other_pes:  there are bergs still in halos!')
      END IF
    END IF
    CALL mppsync
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_lbc_mpp
  SUBROUTINE icb_lbc_mpp_nfld
    TYPE(iceberg), POINTER :: tmpberg, this
    TYPE(point), POINTER :: pt
    INTEGER :: ibergs_to_send
    INTEGER :: ibergs_to_rcv
    INTEGER :: iiglo, ijglo, jk, jn
    INTEGER :: ifldproc, iproc, ipts
    INTEGER :: iine, ijne
    INTEGER :: jjn
    REAL(KIND = wp), DIMENSION(0 : 3) :: zsbergs, znbergs
    INTEGER :: iml_req1, iml_req2, iml_err
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: iml_stat
    nicbfldnsend(:) = 0
    nicbfldexpect(:) = 0
    nicbfldreq(:) = 0
    DO jn = 1, jpni
      IF (nicbfldproc(jn) /= - 1) THEN
        ifldproc = nicbfldproc(jn)
        nicbfldnsend(jn) = 0
        IF (ASSOCIATED(first_berg)) THEN
          this => first_berg
          DO WHILE (ASSOCIATED(this))
            pt => this % current_point
            iine = INT(pt % xi + 0.5)
            ijne = INT(pt % yj + 0.5)
            iproc = nicbflddest(mi1(iine))
            IF (ijne .GT. mjg(nicbej)) THEN
              IF (iproc == ifldproc) THEN
                IF (iproc /= narea) THEN
                  tmpberg => this
                  nicbfldnsend(jn) = nicbfldnsend(jn) + 1
                END IF
              END IF
            END IF
            this => this % next
          END DO
        END IF
      END IF
    END DO
    DO jn = 1, jpni
      IF (nicbfldproc(jn) /= - 1) THEN
        ifldproc = nicbfldproc(jn)
        IF (ifldproc == narea) CYCLE
        zsbergs(0) = narea
        zsbergs(1) = nicbfldnsend(jn)
        CALL mppsend(21, zsbergs(0 : 1), 2, ifldproc - 1, nicbfldreq(jn))
      END IF
    END DO
    DO jn = 1, jpni
      IF (nicbfldproc(jn) /= - 1) THEN
        ifldproc = nicbfldproc(jn)
        IF (ifldproc == narea) CYCLE
        CALL mpprecv(21, znbergs(1 : 2), 2)
        DO jjn = 1, jpni
          IF (nicbfldproc(jjn) .EQ. INT(znbergs(1))) EXIT
        END DO
        IF (jjn .GT. jpni .AND. nn_verbose_level > 0) WRITE(numicb, FMT = *) 'ICB ERROR'
        nicbfldexpect(jjn) = INT(znbergs(2))
      END IF
    END DO
    DO jn = 1, jpni
      IF (nicbfldproc(jn) /= - 1) THEN
        ifldproc = nicbfldproc(jn)
        IF (ifldproc == narea) CYCLE
        IF (l_isend) CALL mpi_wait(nicbfldreq(jn), iml_stat, iml_err)
      END IF
    END DO
    DO jn = 1, jpni
      IF (nicbfldproc(jn) /= - 1) THEN
        ifldproc = nicbfldproc(jn)
        ibergs_to_send = 0
        IF (ASSOCIATED(first_berg)) THEN
          this => first_berg
          DO WHILE (ASSOCIATED(this))
            pt => this % current_point
            iine = INT(pt % xi + 0.5)
            ijne = INT(pt % yj + 0.5)
            ipts = nicbfldpts(mi1(iine))
            iproc = nicbflddest(mi1(iine))
            IF (ijne .GT. mjg(nicbej)) THEN
              IF (iproc == ifldproc) THEN
                ijglo = INT(ipts / nicbpack)
                iiglo = ipts - nicbpack * ijglo
                pt % xi = iiglo - (pt % xi - REAL(iine, wp))
                pt % yj = ijglo - (pt % yj - REAL(ijne, wp))
                pt % uvel = - 1._wp * pt % uvel
                pt % vvel = - 1._wp * pt % vvel
                IF (iproc /= narea) THEN
                  tmpberg => this
                  ibergs_to_send = ibergs_to_send + 1
                  IF (nn_verbose_level >= 4) THEN
                    WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' packing berg ', tmpberg % number(:), ' for north fold'
                    CALL flush(numicb)
                  END IF
                  CALL icb_pack_into_buffer(tmpberg, obuffer_f, ibergs_to_send)
                  CALL icb_utl_delete(first_berg, tmpberg)
                END IF
              END IF
            END IF
            this => this % next
          END DO
        END IF
        IF (nn_verbose_level >= 3) THEN
          WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' send nfld: ', ibergs_to_send
          CALL flush(numicb)
        END IF
        IF (ifldproc == narea) CYCLE
        IF (ibergs_to_send > 0) CALL mppsend(12, obuffer_f % data, ibergs_to_send * jp_buffer_width, ifldproc - 1, nicbfldreq(jn))
      END IF
    END DO
    DO jn = 1, jpni
      IF (nicbfldproc(jn) /= - 1) THEN
        ifldproc = nicbfldproc(jn)
        IF (ifldproc == narea) CYCLE
        ibergs_to_rcv = nicbfldexpect(jn)
        IF (ibergs_to_rcv > 0) THEN
          CALL icb_increase_ibuffer(ibuffer_f, ibergs_to_rcv)
          CALL mpprecv(12, ibuffer_f % data, ibergs_to_rcv * jp_buffer_width, ifldproc - 1)
        END IF
        DO jk = 1, ibergs_to_rcv
          IF (nn_verbose_level >= 4) THEN
            WRITE(numicb, FMT = *) 'bergstep ', nktberg, ' unpacking berg ', INT(ibuffer_f % data(16, jk)), ' from north fold'
            CALL flush(numicb)
          END IF
          CALL icb_unpack_from_buffer(first_berg, ibuffer_f, jk)
        END DO
      END IF
    END DO
    DO jn = 1, jpni
      IF (nicbfldproc(jn) /= - 1) THEN
        ifldproc = nicbfldproc(jn)
        IF (ifldproc == narea) CYCLE
        IF (l_isend) CALL mpi_wait(nicbfldreq(jn), iml_stat, iml_err)
      END IF
    END DO
  END SUBROUTINE icb_lbc_mpp_nfld
  SUBROUTINE icb_pack_into_buffer(berg, pbuff, kb)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(iceberg), POINTER :: berg
    TYPE(buffer), POINTER :: pbuff
    INTEGER, INTENT(IN) :: kb
    INTEGER :: k
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_pack_into_buffer', 'r0', 0, 0)
    IF (.NOT. ASSOCIATED(pbuff)) CALL icb_increase_buffer(pbuff, jp_delta_buf)
    IF (kb .GT. pbuff % size) CALL icb_increase_buffer(pbuff, jp_delta_buf)
    pbuff % data(1, kb) = berg % current_point % lon
    pbuff % data(2, kb) = berg % current_point % lat
    pbuff % data(3, kb) = berg % current_point % uvel
    pbuff % data(4, kb) = berg % current_point % vvel
    pbuff % data(5, kb) = berg % current_point % xi
    pbuff % data(6, kb) = berg % current_point % yj
    pbuff % data(7, kb) = FLOAT(berg % current_point % year)
    pbuff % data(8, kb) = berg % current_point % day
    pbuff % data(9, kb) = berg % current_point % mass
    pbuff % data(10, kb) = berg % current_point % thickness
    pbuff % data(11, kb) = berg % current_point % width
    pbuff % data(12, kb) = berg % current_point % length
    pbuff % data(13, kb) = berg % current_point % mass_of_bits
    pbuff % data(14, kb) = berg % current_point % heat_density
    pbuff % data(15, kb) = berg % mass_scaling
    DO k = 1, nkounts
      pbuff % data(15 + k, kb) = REAL(berg % number(k), wp)
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_pack_into_buffer
  SUBROUTINE icb_unpack_from_buffer(first, pbuff, kb)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(iceberg), POINTER :: first
    TYPE(buffer), POINTER :: pbuff
    INTEGER, INTENT(IN) :: kb
    TYPE(iceberg) :: currentberg
    TYPE(point) :: pt
    INTEGER :: ik
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_unpack_from_buffer', 'r0', 0, 0)
    pt % lon = pbuff % data(1, kb)
    pt % lat = pbuff % data(2, kb)
    pt % uvel = pbuff % data(3, kb)
    pt % vvel = pbuff % data(4, kb)
    pt % xi = pbuff % data(5, kb)
    pt % yj = pbuff % data(6, kb)
    pt % year = INT(pbuff % data(7, kb))
    pt % day = pbuff % data(8, kb)
    pt % mass = pbuff % data(9, kb)
    pt % thickness = pbuff % data(10, kb)
    pt % width = pbuff % data(11, kb)
    pt % length = pbuff % data(12, kb)
    pt % mass_of_bits = pbuff % data(13, kb)
    pt % heat_density = pbuff % data(14, kb)
    currentberg % mass_scaling = pbuff % data(15, kb)
    DO ik = 1, nkounts
      currentberg % number(ik) = INT(pbuff % data(15 + ik, kb))
    END DO
    CALL icb_utl_add(currentberg, pt)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_unpack_from_buffer
  SUBROUTINE icb_increase_buffer(old, kdelta)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(buffer), POINTER :: old
    INTEGER, INTENT(IN) :: kdelta
    TYPE(buffer), POINTER :: new
    INTEGER :: inew_size
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_increase_buffer', 'r0', 0, 0)
    IF (.NOT. ASSOCIATED(old)) THEN
      inew_size = kdelta
    ELSE
      inew_size = old % size + kdelta
    END IF
    ALLOCATE(new)
    ALLOCATE(new % data(jp_buffer_width, inew_size))
    new % size = inew_size
    IF (ASSOCIATED(old)) THEN
      new % data(:, 1 : old % size) = old % data(:, 1 : old % size)
      DEALLOCATE(old % data)
      DEALLOCATE(old)
    END IF
    old => new
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_increase_buffer
  SUBROUTINE icb_increase_ibuffer(old, kdelta)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(buffer), POINTER :: old
    INTEGER, INTENT(IN) :: kdelta
    TYPE(buffer), POINTER :: new
    INTEGER :: inew_size, iold_size
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_increase_ibuffer', 'r0', 0, 0)
    IF (.NOT. ASSOCIATED(old)) THEN
      inew_size = kdelta + jp_delta_buf
      iold_size = 0
    ELSE
      iold_size = old % size
      IF (kdelta .LT. old % size) THEN
        inew_size = old % size + kdelta
      ELSE
        inew_size = kdelta + jp_delta_buf
      END IF
    END IF
    IF (iold_size .NE. inew_size) THEN
      ALLOCATE(new)
      ALLOCATE(new % data(jp_buffer_width, inew_size))
      new % size = inew_size
      IF (ASSOCIATED(old)) THEN
        new % data(:, 1 : old % size) = old % data(:, 1 : old % size)
        DEALLOCATE(old % data)
        DEALLOCATE(old)
      END IF
      old => new
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_increase_ibuffer
END MODULE icblbc