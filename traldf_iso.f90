MODULE traldf_iso
  USE oce
  USE dom_oce
  USE trc_oce
  USE zdf_oce
  USE ldftra
  USE ldfslp
  USE diaptr
  USE diaar5
  USE in_out_manager
  USE iom
  USE phycst
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_ldf_iso
  LOGICAL :: l_ptr
  LOGICAL :: l_hst
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: zdkt, zdk1t, z2d
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: zdit, zdjt, zftu, zftv, ztfw
  CONTAINS
  SUBROUTINE tra_ldf_iso(kt, kit000, cdtype, pahu, pahv, pgu, pgv, pgui, pgvi, ptb, ptbb, pta, kjpt, kpass)
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kit000
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    INTEGER, INTENT(IN) :: kjpt
    INTEGER, INTENT(IN) :: kpass
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pahu, pahv
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt), INTENT(IN) :: pgu, pgv
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt), INTENT(IN) :: pgui, pgvi
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: ptb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: ptbb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    INTEGER :: ji, jj, jk, jn
    INTEGER :: ikt
    INTEGER :: ierr
    REAL(KIND = wp) :: zmsku, zahu_w, zabe1, zcof1, zcoef3
    REAL(KIND = wp) :: zmskv, zahv_w, zabe2, zcof2, zcoef4
    REAL(KIND = wp) :: zcoef0, ze3w_2, zsign, z2dt, z1_2dt
    IF (kpass == 1 .AND. kt == kit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'tra_ldf_iso : rotated laplacian diffusion operator on ', cdtype
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~'
      akz(:, :, :) = 0._wp
      ah_wslp2(:, :, :) = 0._wp
      ALLOCATE(zdit(jpi, jpj, jpk), zdjt(jpi, jpj, jpk), zftu(jpi, jpj, jpk), zftv(jpi, jpj, jpk), ztfw(jpi, jpj, jpk))
      ALLOCATE(zdkt(jpi, jpj), zdk1t(jpi, jpj), z2d(jpi, jpj))
    END IF
    l_hst = .FALSE.
    l_ptr = .FALSE.
    IF (cdtype == 'TRA' .AND. ln_diaptr) l_ptr = .TRUE.
    IF (cdtype == 'TRA' .AND. (iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR. iom_use("uadv_salttr") .OR. &
&iom_use("vadv_salttr"))) l_hst = .TRUE.
    IF (neuler == 0 .AND. kt == nit000) THEN
      z2dt = rdt
    ELSE
      z2dt = 2. * rdt
    END IF
    z1_2dt = 1._wp / z2dt
    IF (kpass == 1) THEN
      zsign = 1._wp
    ELSE
      zsign = - 1._wp
    END IF
    IF (kpass == 1) THEN
      !$OMP parallel default(shared), private(ji,jj,jk,zahu_w,zahv_w,zmsku,zmskv)
      !$OMP do schedule(static)
      DO jk = 2, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zmsku = wmask(ji, jj, jk) / MAX(umask(ji, jj, jk - 1) + umask(ji - 1, jj, jk) + umask(ji - 1, jj, jk - 1) + umask(ji, &
&jj, jk), 1._wp)
            zmskv = wmask(ji, jj, jk) / MAX(vmask(ji, jj, jk - 1) + vmask(ji, jj - 1, jk) + vmask(ji, jj - 1, jk - 1) + vmask(ji, &
&jj, jk), 1._wp)
            zahu_w = (pahu(ji, jj, jk - 1) + pahu(ji - 1, jj, jk) + pahu(ji - 1, jj, jk - 1) + pahu(ji, jj, jk)) * zmsku
            zahv_w = (pahv(ji, jj, jk - 1) + pahv(ji, jj - 1, jk) + pahv(ji, jj - 1, jk - 1) + pahv(ji, jj, jk)) * zmskv
            ah_wslp2(ji, jj, jk) = zahu_w * wslpi(ji, jj, jk) * wslpi(ji, jj, jk) + zahv_w * wslpj(ji, jj, jk) * wslpj(ji, jj, jk)
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
      IF (ln_traldf_msc) THEN
        !$OMP parallel default(shared), private(ji,jj,jk)
        !$OMP do schedule(static)
        DO jk = 2, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              akz(ji, jj, jk) = 0.25_wp * ((pahu(ji, jj, jk) + pahu(ji, jj, jk - 1)) / (e1u(ji, jj) * e1u(ji, jj)) + (pahu(ji - 1, &
&jj, jk) + pahu(ji - 1, jj, jk - 1)) / (e1u(ji - 1, jj) * e1u(ji - 1, jj)) + (pahv(ji, jj, jk) + pahv(ji, jj, jk - 1)) / (e2v(ji, &
&jj) * e2v(ji, jj)) + (pahv(ji, jj - 1, jk) + pahv(ji, jj - 1, jk - 1)) / (e2v(ji, jj - 1) * e2v(ji, jj - 1)))
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
        IF (ln_traldf_blp) THEN
          !$OMP parallel default(shared), private(ji,jj,jk)
          !$OMP do schedule(static)
          DO jk = 2, jpkm1
            DO jj = 1, jpjm1
              DO ji = 1, jpim1
                akz(ji, jj, jk) = 16._wp * ah_wslp2(ji, jj, jk) * (akz(ji, jj, jk) + ah_wslp2(ji, jj, jk) / (e3w_n(ji, jj, jk) * &
&e3w_n(ji, jj, jk)))
              END DO
            END DO
          END DO
          !$OMP end do
          !$OMP end parallel
        ELSE IF (ln_traldf_lap) THEN
          !$OMP parallel default(shared), private(ji,jj,jk,zcoef0,ze3w_2)
          !$OMP do schedule(static)
          DO jk = 2, jpkm1
            DO jj = 1, jpjm1
              DO ji = 1, jpim1
                ze3w_2 = e3w_n(ji, jj, jk) * e3w_n(ji, jj, jk)
                zcoef0 = z2dt * (akz(ji, jj, jk) + ah_wslp2(ji, jj, jk) / ze3w_2)
                akz(ji, jj, jk) = MAX(zcoef0 - 0.5_wp, 0._wp) * ze3w_2 * z1_2dt
              END DO
            END DO
          END DO
          !$OMP end do
          !$OMP end parallel
        END IF
      ELSE
        akz(:, :, :) = ah_wslp2(:, :, :)
      END IF
    END IF
    DO jn = 1, kjpt
      zdit(1, :, :) = 0._wp
      zdit(jpi, :, :) = 0._wp
      zdjt(1, :, :) = 0._wp
      zdjt(jpi, :, :) = 0._wp
      !$OMP parallel default(shared), private(ji,jj,jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zdit(ji, jj, jk) = (ptb(ji + 1, jj, jk, jn) - ptb(ji, jj, jk, jn)) * umask(ji, jj, jk)
            zdjt(ji, jj, jk) = (ptb(ji, jj + 1, jk, jn) - ptb(ji, jj, jk, jn)) * vmask(ji, jj, jk)
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
      IF (ln_zps) THEN
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zdit(ji, jj, mbku(ji, jj)) = pgu(ji, jj, jn)
            zdjt(ji, jj, mbkv(ji, jj)) = pgv(ji, jj, jn)
          END DO
        END DO
        IF (ln_isfcav) THEN
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              IF (miku(ji, jj) > 1) zdit(ji, jj, miku(ji, jj)) = pgui(ji, jj, jn)
              IF (mikv(ji, jj) > 1) zdjt(ji, jj, mikv(ji, jj)) = pgvi(ji, jj, jn)
            END DO
          END DO
        END IF
      END IF
      ! !$OMP parallel default(shared), private(ji,jj,jk,zabe1,zabe2,zcof1,zcof2,zmsku,zmskv)
      !  CDe possible race condition wth zdk1t and zdkt being shared??
      ! !$OMP do schedule(static)
      DO jk = 1, jpkm1
        zdk1t(:, :) = (ptb(:, :, jk, jn) - ptb(:, :, jk + 1, jn)) * wmask(:, :, jk + 1)
        IF (jk == 1) THEN
          zdkt(:, :) = zdk1t(:, :)
        ELSE
          zdkt(:, :) = (ptb(:, :, jk - 1, jn) - ptb(:, :, jk, jn)) * wmask(:, :, jk)
        END IF
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zabe1 = pahu(ji, jj, jk) * e2_e1u(ji, jj) * e3u_n(ji, jj, jk)
            zabe2 = pahv(ji, jj, jk) * e1_e2v(ji, jj) * e3v_n(ji, jj, jk)
            zmsku = 1. / MAX(wmask(ji + 1, jj, jk) + wmask(ji, jj, jk + 1) + wmask(ji + 1, jj, jk + 1) + wmask(ji, jj, jk), 1.)
            zmskv = 1. / MAX(wmask(ji, jj + 1, jk) + wmask(ji, jj, jk + 1) + wmask(ji, jj + 1, jk + 1) + wmask(ji, jj, jk), 1.)
            zcof1 = - pahu(ji, jj, jk) * e2u(ji, jj) * uslp(ji, jj, jk) * zmsku
            zcof2 = - pahv(ji, jj, jk) * e1v(ji, jj) * vslp(ji, jj, jk) * zmskv
            zftu(ji, jj, jk) = (zabe1 * zdit(ji, jj, jk) + zcof1 * (zdkt(ji + 1, jj) + zdk1t(ji, jj) + zdk1t(ji + 1, jj) + &
&zdkt(ji, jj))) * umask(ji, jj, jk)
            zftv(ji, jj, jk) = (zabe2 * zdjt(ji, jj, jk) + zcof2 * (zdkt(ji, jj + 1) + zdk1t(ji, jj) + zdk1t(ji, jj + 1) + &
&zdkt(ji, jj))) * vmask(ji, jj, jk)
          END DO
        END DO
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            pta(ji, jj, jk, jn) = pta(ji, jj, jk, jn) + zsign * (zftu(ji, jj, jk) - zftu(ji - 1, jj, jk) + zftv(ji, jj, jk) - &
&zftv(ji, jj - 1, jk)) * r1_e1e2t(ji, jj) / e3t_n(ji, jj, jk)
          END DO
        END DO
      END DO
      ! !$OMP end do
      ! !$OMP end parallel
      ztfw(1, :, :) = 0._wp
      ztfw(jpi, :, :) = 0._wp
      ztfw(:, :, 1) = 0._wp
      ztfw(:, :, jpk) = 0._wp
      !$OMP parallel default(shared), private(ji,jj,jk,zahu_w,zahv_w,zcoef3,zcoef4,zmsku,zmskv)
      !$OMP do schedule(static)
      DO jk = 2, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zmsku = wmask(ji, jj, jk) / MAX(umask(ji, jj, jk - 1) + umask(ji - 1, jj, jk) + umask(ji - 1, jj, jk - 1) + umask(ji, &
&jj, jk), 1._wp)
            zmskv = wmask(ji, jj, jk) / MAX(vmask(ji, jj, jk - 1) + vmask(ji, jj - 1, jk) + vmask(ji, jj - 1, jk - 1) + vmask(ji, &
&jj, jk), 1._wp)
            zahu_w = (pahu(ji, jj, jk - 1) + pahu(ji - 1, jj, jk) + pahu(ji - 1, jj, jk - 1) + pahu(ji, jj, jk)) * zmsku
            zahv_w = (pahv(ji, jj, jk - 1) + pahv(ji, jj - 1, jk) + pahv(ji, jj - 1, jk - 1) + pahv(ji, jj, jk)) * zmskv
            zcoef3 = - zahu_w * e2t(ji, jj) * zmsku * wslpi(ji, jj, jk)
            zcoef4 = - zahv_w * e1t(ji, jj) * zmskv * wslpj(ji, jj, jk)
            ztfw(ji, jj, jk) = zcoef3 * (zdit(ji, jj, jk - 1) + zdit(ji - 1, jj, jk) + zdit(ji - 1, jj, jk - 1) + zdit(ji, jj, &
&jk)) + zcoef4 * (zdjt(ji, jj, jk - 1) + zdjt(ji, jj - 1, jk) + zdjt(ji, jj - 1, jk - 1) + zdjt(ji, jj, jk))
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
      IF (ln_traldf_lap) THEN
        !$OMP parallel default(shared), private(ji,jj,jk)
        !$OMP do schedule(static)
        DO jk = 2, jpkm1
          DO jj = 1, jpjm1
            DO ji = 2, jpim1
              ztfw(ji, jj, jk) = ztfw(ji, jj, jk) + e1e2t(ji, jj) / e3w_n(ji, jj, jk) * wmask(ji, jj, jk) * (ah_wslp2(ji, jj, jk) &
&- akz(ji, jj, jk)) * (ptb(ji, jj, jk - 1, jn) - ptb(ji, jj, jk, jn))
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
      ELSE
        SELECT CASE (kpass)
        CASE (1)
          !$OMP parallel default(shared), private(ji,jj,jk)
          !$OMP do schedule(static)
          DO jk = 2, jpkm1
            DO jj = 1, jpjm1
              DO ji = 2, jpim1
                ztfw(ji, jj, jk) = ztfw(ji, jj, jk) + ah_wslp2(ji, jj, jk) * e1e2t(ji, jj) * (ptb(ji, jj, jk - 1, jn) - ptb(ji, &
&jj, jk, jn)) / e3w_n(ji, jj, jk) * wmask(ji, jj, jk)
              END DO
            END DO
          END DO
          !$OMP end do
          !$OMP end parallel
        CASE (2)
          !$OMP parallel default(shared), private(ji,jj,jk)
          !$OMP do schedule(static)
          DO jk = 2, jpkm1
            DO jj = 1, jpjm1
              DO ji = 2, jpim1
                ztfw(ji, jj, jk) = ztfw(ji, jj, jk) + e1e2t(ji, jj) / e3w_n(ji, jj, jk) * wmask(ji, jj, jk) * (ah_wslp2(ji, jj, &
&jk) * (ptb(ji, jj, jk - 1, jn) - ptb(ji, jj, jk, jn)) + akz(ji, jj, jk) * (ptbb(ji, jj, jk - 1, jn) - ptbb(ji, jj, jk, jn)))
              END DO
            END DO
          END DO
          !$OMP end do
          !$OMP end parallel
        END SELECT
      END IF
      !$OMP parallel default(shared), private(ji,jj,jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            pta(ji, jj, jk, jn) = pta(ji, jj, jk, jn) + zsign * (ztfw(ji, jj, jk) - ztfw(ji, jj, jk + 1)) * r1_e1e2t(ji, jj) / &
&e3t_n(ji, jj, jk)
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
      IF ((kpass == 1 .AND. ln_traldf_lap) .OR. (kpass == 2 .AND. ln_traldf_blp)) THEN
        IF (l_ptr) CALL dia_ptr_hst(jn, 'ldf', - zftv(:, :, :))
        IF (l_hst) CALL dia_ar5_hst(jn, 'ldf', - zftu(:, :, :), - zftv(:, :, :))
      END IF
    END DO
  END SUBROUTINE tra_ldf_iso
END MODULE traldf_iso
