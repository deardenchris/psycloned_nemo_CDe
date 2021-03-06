MODULE geo2ocean
  USE dom_oce
  USE phycst
  USE in_out_manager
  USE lbclnk
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rot_rep
  PUBLIC :: geo2oce
  PUBLIC :: oce2geo
  PUBLIC :: obs_rot
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: gsint, gcost
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: gsinu, gcosu
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: gsinv, gcosv
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: gsinf, gcosf
  LOGICAL, SAVE, DIMENSION(4) :: linit = .FALSE.
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: gsinlon, gcoslon, gsinlat, gcoslat
  LOGICAL :: lmust_init = .TRUE.
  CONTAINS
  SUBROUTINE rot_rep(pxin, pyin, cd_type, cdtodo, prot)
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN   ) :: pxin, pyin
    CHARACTER(LEN = 1), INTENT(IN   ) :: cd_type
    CHARACTER(LEN = 5), INTENT(IN   ) :: cdtodo
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(  OUT) :: prot
    IF (lmust_init) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) ' rot_rep: coordinate transformation : geographic <==> model (i,j)-components'
      IF (lwp) WRITE(numout, FMT = *) ' ~~~~~~~~    '
      CALL angle
      lmust_init = .FALSE.
    END IF
    SELECT CASE (cdtodo)
    CASE ('en->i')
      SELECT CASE (cd_type)
      CASE ('T')
        !$ACC KERNELS
        prot(:, :) = pxin(:, :) * gcost(:, :) + pyin(:, :) * gsint(:, :)
        !$ACC END KERNELS
      CASE ('U')
        !$ACC KERNELS
        prot(:, :) = pxin(:, :) * gcosu(:, :) + pyin(:, :) * gsinu(:, :)
        !$ACC END KERNELS
      CASE ('V')
        !$ACC KERNELS
        prot(:, :) = pxin(:, :) * gcosv(:, :) + pyin(:, :) * gsinv(:, :)
        !$ACC END KERNELS
      CASE ('F')
        !$ACC KERNELS
        prot(:, :) = pxin(:, :) * gcosf(:, :) + pyin(:, :) * gsinf(:, :)
        !$ACC END KERNELS
      CASE DEFAULT
        CALL ctl_stop('Only T, U, V and F grid points are coded')
      END SELECT
    CASE ('en->j')
      SELECT CASE (cd_type)
      CASE ('T')
        !$ACC KERNELS
        prot(:, :) = pyin(:, :) * gcost(:, :) - pxin(:, :) * gsint(:, :)
        !$ACC END KERNELS
      CASE ('U')
        !$ACC KERNELS
        prot(:, :) = pyin(:, :) * gcosu(:, :) - pxin(:, :) * gsinu(:, :)
        !$ACC END KERNELS
      CASE ('V')
        !$ACC KERNELS
        prot(:, :) = pyin(:, :) * gcosv(:, :) - pxin(:, :) * gsinv(:, :)
        !$ACC END KERNELS
      CASE ('F')
        !$ACC KERNELS
        prot(:, :) = pyin(:, :) * gcosf(:, :) - pxin(:, :) * gsinf(:, :)
        !$ACC END KERNELS
      CASE DEFAULT
        CALL ctl_stop('Only T, U, V and F grid points are coded')
      END SELECT
    CASE ('ij->e')
      SELECT CASE (cd_type)
      CASE ('T')
        !$ACC KERNELS
        prot(:, :) = pxin(:, :) * gcost(:, :) - pyin(:, :) * gsint(:, :)
        !$ACC END KERNELS
      CASE ('U')
        !$ACC KERNELS
        prot(:, :) = pxin(:, :) * gcosu(:, :) - pyin(:, :) * gsinu(:, :)
        !$ACC END KERNELS
      CASE ('V')
        !$ACC KERNELS
        prot(:, :) = pxin(:, :) * gcosv(:, :) - pyin(:, :) * gsinv(:, :)
        !$ACC END KERNELS
      CASE ('F')
        !$ACC KERNELS
        prot(:, :) = pxin(:, :) * gcosf(:, :) - pyin(:, :) * gsinf(:, :)
        !$ACC END KERNELS
      CASE DEFAULT
        CALL ctl_stop('Only T, U, V and F grid points are coded')
      END SELECT
    CASE ('ij->n')
      SELECT CASE (cd_type)
      CASE ('T')
        !$ACC KERNELS
        prot(:, :) = pyin(:, :) * gcost(:, :) + pxin(:, :) * gsint(:, :)
        !$ACC END KERNELS
      CASE ('U')
        !$ACC KERNELS
        prot(:, :) = pyin(:, :) * gcosu(:, :) + pxin(:, :) * gsinu(:, :)
        !$ACC END KERNELS
      CASE ('V')
        !$ACC KERNELS
        prot(:, :) = pyin(:, :) * gcosv(:, :) + pxin(:, :) * gsinv(:, :)
        !$ACC END KERNELS
      CASE ('F')
        !$ACC KERNELS
        prot(:, :) = pyin(:, :) * gcosf(:, :) + pxin(:, :) * gsinf(:, :)
        !$ACC END KERNELS
      CASE DEFAULT
        CALL ctl_stop('Only T, U, V and F grid points are coded')
      END SELECT
    CASE DEFAULT
      CALL ctl_stop('rot_rep: Syntax Error in the definition of cdtodo')
    END SELECT
  END SUBROUTINE rot_rep
  SUBROUTINE angle
    INTEGER :: ji, jj
    INTEGER :: ierr
    REAL(KIND = wp) :: zlam, zphi
    REAL(KIND = wp) :: zlan, zphh
    REAL(KIND = wp) :: zxnpt, zynpt, znnpt
    REAL(KIND = wp) :: zxnpu, zynpu, znnpu
    REAL(KIND = wp) :: zxnpv, zynpv, znnpv
    REAL(KIND = wp) :: zxnpf, zynpf, znnpf
    REAL(KIND = wp) :: zxvvt, zyvvt, znvvt
    REAL(KIND = wp) :: zxffu, zyffu, znffu
    REAL(KIND = wp) :: zxffv, zyffv, znffv
    REAL(KIND = wp) :: zxuuf, zyuuf, znuuf
    ALLOCATE(gsint(jpi, jpj), gcost(jpi, jpj), gsinu(jpi, jpj), gcosu(jpi, jpj), gsinv(jpi, jpj), gcosv(jpi, jpj), gsinf(jpi, jpj), gcosf(jpi, jpj), STAT = ierr)
    IF (lk_mpp) CALL mpp_sum(ierr)
    IF (ierr /= 0) CALL ctl_stop('angle: unable to allocate arrays')
    !$ACC KERNELS
    DO jj = 2, jpjm1
      DO ji = 2, jpi
        zlam = glamt(ji, jj)
        zphi = gphit(ji, jj)
        zxnpt = 0. - 2. * COS(rad * zlam) * TAN(rpi / 4. - rad * zphi / 2.)
        zynpt = 0. - 2. * SIN(rad * zlam) * TAN(rpi / 4. - rad * zphi / 2.)
        znnpt = zxnpt * zxnpt + zynpt * zynpt
        zlam = glamu(ji, jj)
        zphi = gphiu(ji, jj)
        zxnpu = 0. - 2. * COS(rad * zlam) * TAN(rpi / 4. - rad * zphi / 2.)
        zynpu = 0. - 2. * SIN(rad * zlam) * TAN(rpi / 4. - rad * zphi / 2.)
        znnpu = zxnpu * zxnpu + zynpu * zynpu
        zlam = glamv(ji, jj)
        zphi = gphiv(ji, jj)
        zxnpv = 0. - 2. * COS(rad * zlam) * TAN(rpi / 4. - rad * zphi / 2.)
        zynpv = 0. - 2. * SIN(rad * zlam) * TAN(rpi / 4. - rad * zphi / 2.)
        znnpv = zxnpv * zxnpv + zynpv * zynpv
        zlam = glamf(ji, jj)
        zphi = gphif(ji, jj)
        zxnpf = 0. - 2. * COS(rad * zlam) * TAN(rpi / 4. - rad * zphi / 2.)
        zynpf = 0. - 2. * SIN(rad * zlam) * TAN(rpi / 4. - rad * zphi / 2.)
        znnpf = zxnpf * zxnpf + zynpf * zynpf
        zlam = glamv(ji, jj)
        zphi = gphiv(ji, jj)
        zlan = glamv(ji, jj - 1)
        zphh = gphiv(ji, jj - 1)
        zxvvt = 2. * COS(rad * zlam) * TAN(rpi / 4. - rad * zphi / 2.) - 2. * COS(rad * zlan) * TAN(rpi / 4. - rad * zphh / 2.)
        zyvvt = 2. * SIN(rad * zlam) * TAN(rpi / 4. - rad * zphi / 2.) - 2. * SIN(rad * zlan) * TAN(rpi / 4. - rad * zphh / 2.)
        znvvt = SQRT(znnpt * (zxvvt * zxvvt + zyvvt * zyvvt))
        znvvt = MAX(znvvt, 1.E-14)
        zlam = glamf(ji, jj)
        zphi = gphif(ji, jj)
        zlan = glamf(ji, jj - 1)
        zphh = gphif(ji, jj - 1)
        zxffu = 2. * COS(rad * zlam) * TAN(rpi / 4. - rad * zphi / 2.) - 2. * COS(rad * zlan) * TAN(rpi / 4. - rad * zphh / 2.)
        zyffu = 2. * SIN(rad * zlam) * TAN(rpi / 4. - rad * zphi / 2.) - 2. * SIN(rad * zlan) * TAN(rpi / 4. - rad * zphh / 2.)
        znffu = SQRT(znnpu * (zxffu * zxffu + zyffu * zyffu))
        znffu = MAX(znffu, 1.E-14)
        zlam = glamf(ji, jj)
        zphi = gphif(ji, jj)
        zlan = glamf(ji - 1, jj)
        zphh = gphif(ji - 1, jj)
        zxffv = 2. * COS(rad * zlam) * TAN(rpi / 4. - rad * zphi / 2.) - 2. * COS(rad * zlan) * TAN(rpi / 4. - rad * zphh / 2.)
        zyffv = 2. * SIN(rad * zlam) * TAN(rpi / 4. - rad * zphi / 2.) - 2. * SIN(rad * zlan) * TAN(rpi / 4. - rad * zphh / 2.)
        znffv = SQRT(znnpv * (zxffv * zxffv + zyffv * zyffv))
        znffv = MAX(znffv, 1.E-14)
        zlam = glamu(ji, jj + 1)
        zphi = gphiu(ji, jj + 1)
        zlan = glamu(ji, jj)
        zphh = gphiu(ji, jj)
        zxuuf = 2. * COS(rad * zlam) * TAN(rpi / 4. - rad * zphi / 2.) - 2. * COS(rad * zlan) * TAN(rpi / 4. - rad * zphh / 2.)
        zyuuf = 2. * SIN(rad * zlam) * TAN(rpi / 4. - rad * zphi / 2.) - 2. * SIN(rad * zlan) * TAN(rpi / 4. - rad * zphh / 2.)
        znuuf = SQRT(znnpf * (zxuuf * zxuuf + zyuuf * zyuuf))
        znuuf = MAX(znuuf, 1.E-14)
        gsint(ji, jj) = (zxnpt * zyvvt - zynpt * zxvvt) / znvvt
        gcost(ji, jj) = (zxnpt * zxvvt + zynpt * zyvvt) / znvvt
        gsinu(ji, jj) = (zxnpu * zyffu - zynpu * zxffu) / znffu
        gcosu(ji, jj) = (zxnpu * zxffu + zynpu * zyffu) / znffu
        gsinf(ji, jj) = (zxnpf * zyuuf - zynpf * zxuuf) / znuuf
        gcosf(ji, jj) = (zxnpf * zxuuf + zynpf * zyuuf) / znuuf
        gsinv(ji, jj) = (zxnpv * zxffv + zynpv * zyffv) / znffv
        gcosv(ji, jj) = - (zxnpv * zyffv - zynpv * zxffv) / znffv
      END DO
    END DO
    DO jj = 2, jpjm1
      DO ji = 2, jpi
        IF (MOD(ABS(glamv(ji, jj) - glamv(ji, jj - 1)), 360.) < 1.E-8) THEN
          gsint(ji, jj) = 0.
          gcost(ji, jj) = 1.
        END IF
        IF (MOD(ABS(glamf(ji, jj) - glamf(ji, jj - 1)), 360.) < 1.E-8) THEN
          gsinu(ji, jj) = 0.
          gcosu(ji, jj) = 1.
        END IF
        IF (ABS(gphif(ji, jj) - gphif(ji - 1, jj)) < 1.E-8) THEN
          gsinv(ji, jj) = 0.
          gcosv(ji, jj) = 1.
        END IF
        IF (MOD(ABS(glamu(ji, jj) - glamu(ji, jj + 1)), 360.) < 1.E-8) THEN
          gsinf(ji, jj) = 0.
          gcosf(ji, jj) = 1.
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk_multi(gcost, 'T', - 1., gsint, 'T', - 1., gcosu, 'U', - 1., gsinu, 'U', - 1., gcosv, 'V', - 1., gsinv, 'V', - 1., gcosf, 'F', - 1., gsinf, 'F', - 1.)
  END SUBROUTINE angle
  SUBROUTINE geo2oce(pxx, pyy, pzz, cgrid, pte, ptn)
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN   ) :: pxx, pyy, pzz
    CHARACTER(LEN = 1), INTENT(IN   ) :: cgrid
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(  OUT) :: pte, ptn
    REAL(KIND = wp), PARAMETER :: rpi = 3.141592653E0
    REAL(KIND = wp), PARAMETER :: rad = rpi / 180.E0
    INTEGER :: ig
    INTEGER :: ierr
    IF (.NOT. ALLOCATED(gsinlon)) THEN
      ALLOCATE(gsinlon(jpi, jpj, 4), gcoslon(jpi, jpj, 4), gsinlat(jpi, jpj, 4), gcoslat(jpi, jpj, 4), STAT = ierr)
      IF (lk_mpp) CALL mpp_sum(ierr)
      IF (ierr /= 0) CALL ctl_stop('geo2oce: unable to allocate arrays')
    END IF
    SELECT CASE (cgrid)
    CASE ('T')
      ig = 1
      IF (.NOT. linit(ig)) THEN
        !$ACC KERNELS
        gsinlon(:, :, ig) = SIN(rad * glamt(:, :))
        gcoslon(:, :, ig) = COS(rad * glamt(:, :))
        gsinlat(:, :, ig) = SIN(rad * gphit(:, :))
        gcoslat(:, :, ig) = COS(rad * gphit(:, :))
        linit(ig) = .TRUE.
        !$ACC END KERNELS
      END IF
    CASE ('U')
      ig = 2
      IF (.NOT. linit(ig)) THEN
        !$ACC KERNELS
        gsinlon(:, :, ig) = SIN(rad * glamu(:, :))
        gcoslon(:, :, ig) = COS(rad * glamu(:, :))
        gsinlat(:, :, ig) = SIN(rad * gphiu(:, :))
        gcoslat(:, :, ig) = COS(rad * gphiu(:, :))
        linit(ig) = .TRUE.
        !$ACC END KERNELS
      END IF
    CASE ('V')
      ig = 3
      IF (.NOT. linit(ig)) THEN
        !$ACC KERNELS
        gsinlon(:, :, ig) = SIN(rad * glamv(:, :))
        gcoslon(:, :, ig) = COS(rad * glamv(:, :))
        gsinlat(:, :, ig) = SIN(rad * gphiv(:, :))
        gcoslat(:, :, ig) = COS(rad * gphiv(:, :))
        linit(ig) = .TRUE.
        !$ACC END KERNELS
      END IF
    CASE ('F')
      ig = 4
      IF (.NOT. linit(ig)) THEN
        !$ACC KERNELS
        gsinlon(:, :, ig) = SIN(rad * glamf(:, :))
        gcoslon(:, :, ig) = COS(rad * glamf(:, :))
        gsinlat(:, :, ig) = SIN(rad * gphif(:, :))
        gcoslat(:, :, ig) = COS(rad * gphif(:, :))
        linit(ig) = .TRUE.
        !$ACC END KERNELS
      END IF
    CASE DEFAULT
      WRITE(ctmp1, FMT = *) 'geo2oce : bad grid argument : ', cgrid
      CALL ctl_stop(ctmp1)
    END SELECT
    pte = - gsinlon(:, :, ig) * pxx + gcoslon(:, :, ig) * pyy
    ptn = - gcoslon(:, :, ig) * gsinlat(:, :, ig) * pxx - gsinlon(:, :, ig) * gsinlat(:, :, ig) * pyy + gcoslat(:, :, ig) * pzz
  END SUBROUTINE geo2oce
  SUBROUTINE oce2geo(pte, ptn, cgrid, pxx, pyy, pzz)
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT( IN    ) :: pte, ptn
    CHARACTER(LEN = 1), INTENT( IN    ) :: cgrid
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(   OUT ) :: pxx, pyy, pzz
    REAL(KIND = wp), PARAMETER :: rpi = 3.141592653E0
    REAL(KIND = wp), PARAMETER :: rad = rpi / 180.E0
    INTEGER :: ig
    INTEGER :: ierr
    IF (.NOT. ALLOCATED(gsinlon)) THEN
      ALLOCATE(gsinlon(jpi, jpj, 4), gcoslon(jpi, jpj, 4), gsinlat(jpi, jpj, 4), gcoslat(jpi, jpj, 4), STAT = ierr)
      IF (lk_mpp) CALL mpp_sum(ierr)
      IF (ierr /= 0) CALL ctl_stop('oce2geo: unable to allocate arrays')
    END IF
    SELECT CASE (cgrid)
    CASE ('T')
      ig = 1
      IF (.NOT. linit(ig)) THEN
        !$ACC KERNELS
        gsinlon(:, :, ig) = SIN(rad * glamt(:, :))
        gcoslon(:, :, ig) = COS(rad * glamt(:, :))
        gsinlat(:, :, ig) = SIN(rad * gphit(:, :))
        gcoslat(:, :, ig) = COS(rad * gphit(:, :))
        linit(ig) = .TRUE.
        !$ACC END KERNELS
      END IF
    CASE ('U')
      ig = 2
      IF (.NOT. linit(ig)) THEN
        !$ACC KERNELS
        gsinlon(:, :, ig) = SIN(rad * glamu(:, :))
        gcoslon(:, :, ig) = COS(rad * glamu(:, :))
        gsinlat(:, :, ig) = SIN(rad * gphiu(:, :))
        gcoslat(:, :, ig) = COS(rad * gphiu(:, :))
        linit(ig) = .TRUE.
        !$ACC END KERNELS
      END IF
    CASE ('V')
      ig = 3
      IF (.NOT. linit(ig)) THEN
        !$ACC KERNELS
        gsinlon(:, :, ig) = SIN(rad * glamv(:, :))
        gcoslon(:, :, ig) = COS(rad * glamv(:, :))
        gsinlat(:, :, ig) = SIN(rad * gphiv(:, :))
        gcoslat(:, :, ig) = COS(rad * gphiv(:, :))
        linit(ig) = .TRUE.
        !$ACC END KERNELS
      END IF
    CASE ('F')
      ig = 4
      IF (.NOT. linit(ig)) THEN
        !$ACC KERNELS
        gsinlon(:, :, ig) = SIN(rad * glamf(:, :))
        gcoslon(:, :, ig) = COS(rad * glamf(:, :))
        gsinlat(:, :, ig) = SIN(rad * gphif(:, :))
        gcoslat(:, :, ig) = COS(rad * gphif(:, :))
        linit(ig) = .TRUE.
        !$ACC END KERNELS
      END IF
    CASE DEFAULT
      WRITE(ctmp1, FMT = *) 'geo2oce : bad grid argument : ', cgrid
      CALL ctl_stop(ctmp1)
    END SELECT
    pxx = - gsinlon(:, :, ig) * pte - gcoslon(:, :, ig) * gsinlat(:, :, ig) * ptn
    pyy = gcoslon(:, :, ig) * pte - gsinlon(:, :, ig) * gsinlat(:, :, ig) * ptn
    pzz = gcoslat(:, :, ig) * ptn
  END SUBROUTINE oce2geo
  SUBROUTINE obs_rot(psinu, pcosu, psinv, pcosv)
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT( OUT ) :: psinu, pcosu, psinv, pcosv
    IF (lmust_init) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) ' obs_rot : geographic <--> stretched'
      IF (lwp) WRITE(numout, FMT = *) ' ~~~~~~~   coordinate transformation'
      CALL angle
      lmust_init = .FALSE.
    END IF
    !$ACC KERNELS
    psinu(:, :) = gsinu(:, :)
    pcosu(:, :) = gcosu(:, :)
    psinv(:, :) = gsinv(:, :)
    pcosv(:, :) = gcosv(:, :)
    !$ACC END KERNELS
  END SUBROUTINE obs_rot
END MODULE geo2ocean