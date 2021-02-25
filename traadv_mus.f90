MODULE traadv_mus
  USE oce
  USE trc_oce
  USE dom_oce
  USE trd_oce
  USE trdtra
  USE sbcrnf
  USE diaptr
  USE diaar5
  USE iom
  USE in_out_manager
  USE lib_mpp
  USE lbclnk
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_adv_mus
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: upsmsk
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: xind
  LOGICAL :: l_trd
  LOGICAL :: l_ptr
  LOGICAL :: l_hst
  CONTAINS
  SUBROUTINE tra_adv_mus(kt, kit000, cdtype, p2dt, pun, pvn, pwn, ptb, pta, kjpt, ld_msc_ups)
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kit000
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    INTEGER, INTENT(IN) :: kjpt
    LOGICAL, INTENT(IN) :: ld_msc_ups
    REAL(KIND = wp), INTENT(IN) :: p2dt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pun, pvn, pwn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: ptb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    INTEGER :: ji, jj, jk, jn
    INTEGER :: ierr
    REAL(KIND = wp) :: zu, z0u, zzwx, zw, zalpha
    REAL(KIND = wp) :: zv, z0v, zzwy, z0w
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwx, zslpx
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwy, zslpy
    IF (kt == kit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'tra_adv : MUSCL advection scheme on ', cdtype
      IF (lwp) WRITE(numout, FMT = *) '        : mixed up-stream           ', ld_msc_ups
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~'
      IF (lwp) WRITE(numout, FMT = *)
      ALLOCATE(xind(jpi, jpj, jpk), STAT = ierr)
      xind(:, :, :) = 1._wp
      IF (ld_msc_ups) THEN
        ALLOCATE(upsmsk(jpi, jpj), STAT = ierr)
        upsmsk(:, :) = 0._wp
        !$OMP parallel default(shared), private(jk)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          xind(:, :, jk) = 1._wp - MAX(rnfmsk(:, :) * rnfmsk_z(jk), upsmsk(:, :)) * tmask(:, :, jk)
        END DO
        !$OMP end do
        !$OMP end parallel
      END IF
    END IF
    l_trd = .FALSE.
    l_hst = .FALSE.
    l_ptr = .FALSE.
    IF ((cdtype == 'TRA' .AND. l_trdtra) .OR. (cdtype == 'TRC' .AND. l_trdtrc)) l_trd = .TRUE.
    IF (cdtype == 'TRA' .AND. ln_diaptr) l_ptr = .TRUE.
    IF (cdtype == 'TRA' .AND. (iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR. iom_use("uadv_salttr") .OR. &
&iom_use("vadv_salttr"))) l_hst = .TRUE.
    DO jn = 1, kjpt
      zwx(:, :, jpk) = 0._wp
      zwy(:, :, jpk) = 0._wp
      !$OMP parallel default(shared), private(ji,jj,jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zwx(ji, jj, jk) = umask(ji, jj, jk) * (ptb(ji + 1, jj, jk, jn) - ptb(ji, jj, jk, jn))
            zwy(ji, jj, jk) = vmask(ji, jj, jk) * (ptb(ji, jj + 1, jk, jn) - ptb(ji, jj, jk, jn))
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
      CALL lbc_lnk_multi('traadv_mus', zwx, 'U', - 1., zwy, 'V', - 1.)
      zslpx(:, :, jpk) = 0._wp
      zslpy(:, :, jpk) = 0._wp
      !$OMP parallel default(shared), private(ji,jj,jk,z0u,z0v,zalpha,zu,zv,zzwx,zzwy)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 2, jpj
          DO ji = 2, jpi
            zslpx(ji, jj, jk) = (zwx(ji, jj, jk) + zwx(ji - 1, jj, jk)) * (0.25 + SIGN(0.25, zwx(ji, jj, jk) * zwx(ji - 1, jj, jk)))
            zslpy(ji, jj, jk) = (zwy(ji, jj, jk) + zwy(ji, jj - 1, jk)) * (0.25 + SIGN(0.25, zwy(ji, jj, jk) * zwy(ji, jj - 1, jk)))
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 2, jpj
          DO ji = 2, jpi
            zslpx(ji, jj, jk) = SIGN(1., zslpx(ji, jj, jk)) * MIN(ABS(zslpx(ji, jj, jk)), 2. * ABS(zwx(ji - 1, jj, jk)), 2. * &
&ABS(zwx(ji, jj, jk)))
            zslpy(ji, jj, jk) = SIGN(1., zslpy(ji, jj, jk)) * MIN(ABS(zslpy(ji, jj, jk)), 2. * ABS(zwy(ji, jj - 1, jk)), 2. * &
&ABS(zwy(ji, jj, jk)))
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            z0u = SIGN(0.5, pun(ji, jj, jk))
            zalpha = 0.5 - z0u
            zu = z0u - 0.5 * pun(ji, jj, jk) * p2dt * r1_e1e2u(ji, jj) / e3u_n(ji, jj, jk)
            zzwx = ptb(ji + 1, jj, jk, jn) + xind(ji, jj, jk) * zu * zslpx(ji + 1, jj, jk)
            zzwy = ptb(ji, jj, jk, jn) + xind(ji, jj, jk) * zu * zslpx(ji, jj, jk)
            zwx(ji, jj, jk) = pun(ji, jj, jk) * (zalpha * zzwx + (1. - zalpha) * zzwy)
            z0v = SIGN(0.5, pvn(ji, jj, jk))
            zalpha = 0.5 - z0v
            zv = z0v - 0.5 * pvn(ji, jj, jk) * p2dt * r1_e1e2v(ji, jj) / e3v_n(ji, jj, jk)
            zzwx = ptb(ji, jj + 1, jk, jn) + xind(ji, jj, jk) * zv * zslpy(ji, jj + 1, jk)
            zzwy = ptb(ji, jj, jk, jn) + xind(ji, jj, jk) * zv * zslpy(ji, jj, jk)
            zwy(ji, jj, jk) = pvn(ji, jj, jk) * (zalpha * zzwx + (1. - zalpha) * zzwy)
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
      CALL lbc_lnk_multi('traadv_mus', zwx, 'U', - 1., zwy, 'V', - 1.)
      !$OMP parallel default(shared), private(ji,jj,jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            pta(ji, jj, jk, jn) = pta(ji, jj, jk, jn) - (zwx(ji, jj, jk) - zwx(ji - 1, jj, jk) + zwy(ji, jj, jk) - zwy(ji, jj - 1, &
&jk)) * r1_e1e2t(ji, jj) / e3t_n(ji, jj, jk)
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
      IF (l_trd) THEN
        CALL trd_tra(kt, cdtype, jn, jptra_xad, zwx, pun, ptb(:, :, :, jn))
        CALL trd_tra(kt, cdtype, jn, jptra_yad, zwy, pvn, ptb(:, :, :, jn))
      END IF
      IF (l_ptr) CALL dia_ptr_hst(jn, 'adv', zwy(:, :, :))
      IF (l_hst) CALL dia_ar5_hst(jn, 'adv', zwx(:, :, :), zwy(:, :, :))
      zwx(:, :, 1) = 0._wp
      zwx(:, :, jpk) = 0._wp
      !$OMP parallel default(shared), private(jk)
      !$OMP do schedule(static)
      DO jk = 2, jpkm1
        zwx(:, :, jk) = tmask(:, :, jk) * (ptb(:, :, jk - 1, jn) - ptb(:, :, jk, jn))
      END DO
      !$OMP end do
      !$OMP end parallel
      zslpx(:, :, 1) = 0._wp
      !$OMP parallel default(shared), private(ji,jj,jk,z0w,zalpha,zw,zzwx,zzwy)
      !$OMP do schedule(static)
      DO jk = 2, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            zslpx(ji, jj, jk) = (zwx(ji, jj, jk) + zwx(ji, jj, jk + 1)) * (0.25 + SIGN(0.25, zwx(ji, jj, jk) * zwx(ji, jj, jk + 1)))
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP do schedule(static)
      DO jk = 2, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            zslpx(ji, jj, jk) = SIGN(1., zslpx(ji, jj, jk)) * MIN(ABS(zslpx(ji, jj, jk)), 2. * ABS(zwx(ji, jj, jk + 1)), 2. * &
&ABS(zwx(ji, jj, jk)))
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP do schedule(static)
      DO jk = 1, jpk - 2
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            z0w = SIGN(0.5, pwn(ji, jj, jk + 1))
            zalpha = 0.5 + z0w
            zw = z0w - 0.5 * pwn(ji, jj, jk + 1) * p2dt * r1_e1e2t(ji, jj) / e3w_n(ji, jj, jk + 1)
            zzwx = ptb(ji, jj, jk + 1, jn) + xind(ji, jj, jk) * zw * zslpx(ji, jj, jk + 1)
            zzwy = ptb(ji, jj, jk, jn) + xind(ji, jj, jk) * zw * zslpx(ji, jj, jk)
            zwx(ji, jj, jk + 1) = pwn(ji, jj, jk + 1) * (zalpha * zzwx + (1. - zalpha) * zzwy) * wmask(ji, jj, jk)
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
      IF (ln_linssh) THEN
        IF (ln_isfcav) THEN
          DO jj = 1, jpj
            DO ji = 1, jpi
              zwx(ji, jj, mikt(ji, jj)) = pwn(ji, jj, mikt(ji, jj)) * ptb(ji, jj, mikt(ji, jj), jn)
            END DO
          END DO
        ELSE
          zwx(:, :, 1) = pwn(:, :, 1) * ptb(:, :, 1, jn)
        END IF
      END IF
      !$OMP parallel default(shared), private(ji,jj,jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            pta(ji, jj, jk, jn) = pta(ji, jj, jk, jn) - (zwx(ji, jj, jk) - zwx(ji, jj, jk + 1)) * r1_e1e2t(ji, jj) / e3t_n(ji, jj, &
&jk)
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
      IF (l_trd) CALL trd_tra(kt, cdtype, jn, jptra_zad, zwx, pwn, ptb(:, :, :, jn))
    END DO
  END SUBROUTINE tra_adv_mus
END MODULE traadv_mus