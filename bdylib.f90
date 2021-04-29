MODULE bdylib
  USE oce
  USE dom_oce
  USE bdy_oce
  USE phycst
  USE in_out_manager
  USE lbclnk
  USE lib_mpp, ONLY: ctl_stop
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: bdy_frs, bdy_spe, bdy_nmn, bdy_orl
  PUBLIC :: bdy_orlanski_2d
  PUBLIC :: bdy_orlanski_3d
  CONTAINS
  SUBROUTINE bdy_frs(idx, pta, dta)
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: dta
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pta
    REAL(KIND = wp) :: zwgt
    INTEGER :: ib, ik, igrd
    INTEGER :: ii, ij
    !$ACC KERNELS
    igrd = 1
    DO ib = 1, idx % nblen(igrd)
      DO ik = 1, jpkm1
        ii = idx % nbi(ib, igrd)
        ij = idx % nbj(ib, igrd)
        zwgt = idx % nbw(ib, igrd)
        pta(ii, ij, ik) = (pta(ii, ij, ik) + zwgt * (dta(ib, ik) - pta(ii, ij, ik))) * tmask(ii, ij, ik)
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE bdy_frs
  SUBROUTINE bdy_spe(idx, pta, dta)
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: dta
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pta
    REAL(KIND = wp) :: zwgt
    INTEGER :: ib, ik, igrd
    INTEGER :: ii, ij
    !$ACC KERNELS
    igrd = 1
    DO ib = 1, idx % nblenrim(igrd)
      ii = idx % nbi(ib, igrd)
      ij = idx % nbj(ib, igrd)
      DO ik = 1, jpkm1
        pta(ii, ij, ik) = dta(ib, ik) * tmask(ii, ij, ik)
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE bdy_spe
  SUBROUTINE bdy_orl(idx, ptb, pta, dta, ll_npo)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: dta
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: ptb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pta
    LOGICAL, INTENT(IN) :: ll_npo
    INTEGER :: igrd
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bdy_orl', 'r0', 0, 0)
    igrd = 1
    CALL bdy_orlanski_3d(idx, igrd, ptb(:, :, :), pta(:, :, :), dta, ll_npo)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE bdy_orl
  SUBROUTINE bdy_orlanski_2d(idx, igrd, phib, phia, phi_ext, ll_npo)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    INTEGER, INTENT(IN) :: igrd
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: phib
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: phia
    REAL(KIND = wp), DIMENSION(:), INTENT(IN) :: phi_ext
    LOGICAL, INTENT(IN) :: ll_npo
    INTEGER :: jb
    INTEGER :: ii, ij, iibm1, iibm2, ijbm1, ijbm2
    INTEGER :: iijm1, iijp1, ijjm1, ijjp1
    INTEGER :: iibm1jp1, iibm1jm1, ijbm1jp1, ijbm1jm1
    INTEGER :: ii_offset, ij_offset
    INTEGER :: flagu, flagv
    REAL(KIND = wp) :: zmask_x, zmask_y1, zmask_y2
    REAL(KIND = wp) :: zex1, zex2, zey, zey1, zey2
    REAL(KIND = wp) :: zdt, zdx, zdy, znor2, zrx, zry
    REAL(KIND = wp) :: zout, zwgt, zdy_centred
    REAL(KIND = wp) :: zdy_1, zdy_2, zsign_ups
    REAL(KIND = wp), PARAMETER :: zepsilon = 1.E-30
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: pmask
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: pmask_xdif
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: pmask_ydif
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: pe_xdif
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: pe_ydif
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bdy_orlanski_2d', 'r0', 0, 0)
    SELECT CASE (igrd)
    CASE (1)
      pmask => tmask(:, :, 1)
      pmask_xdif => umask(:, :, 1)
      pmask_ydif => vmask(:, :, 1)
      pe_xdif => e1u(:, :)
      pe_ydif => e2v(:, :)
      ii_offset = 0
      ij_offset = 0
    CASE (2)
      pmask => umask(:, :, 1)
      pmask_xdif => tmask(:, :, 1)
      pmask_ydif => fmask(:, :, 1)
      pe_xdif => e1t(:, :)
      pe_ydif => e2f(:, :)
      ii_offset = 1
      ij_offset = 0
    CASE (3)
      pmask => vmask(:, :, 1)
      pmask_xdif => fmask(:, :, 1)
      pmask_ydif => tmask(:, :, 1)
      pe_xdif => e1f(:, :)
      pe_ydif => e2t(:, :)
      ii_offset = 0
      ij_offset = 1
    CASE DEFAULT
      CALL ctl_stop('unrecognised value for igrd in bdy_orlanksi_2d')
    END SELECT
    DO jb = 1, idx % nblenrim(igrd)
      ii = idx % nbi(jb, igrd)
      ij = idx % nbj(jb, igrd)
      flagu = INT(idx % flagu(jb, igrd))
      flagv = INT(idx % flagv(jb, igrd))
      iibm1 = ii + flagu
      iibm2 = ii + 2 * flagu
      ijbm1 = ij + flagv
      ijbm2 = ij + 2 * flagv
      iijm1 = ii - ABS(flagv)
      iijp1 = ii + ABS(flagv)
      ijjm1 = ij - ABS(flagu)
      ijjp1 = ij + ABS(flagu)
      iibm1jm1 = ii + flagu - ABS(flagv)
      iibm1jp1 = ii + flagu + ABS(flagv)
      ijbm1jm1 = ij + flagv - ABS(flagu)
      ijbm1jp1 = ij + flagv + ABS(flagu)
      zex1 = (ABS(iibm1 - iibm2) * pe_xdif(iibm1 + ii_offset, ijbm1) + ABS(ijbm1 - ijbm2) * pe_ydif(iibm1, ijbm1 + ij_offset))
      zex2 = (ABS(iibm1 - iibm2) * pe_xdif(iibm2 + ii_offset, ijbm2) + ABS(ijbm1 - ijbm2) * pe_ydif(iibm2, ijbm2 + ij_offset))
      zey1 = ((iibm1 - iibm1jm1) * pe_xdif(iibm1jm1 + ii_offset, ijbm1jm1) + (ijbm1 - ijbm1jm1) * pe_ydif(iibm1jm1, ijbm1jm1 + ij_offset))
      zey2 = ((iibm1jp1 - iibm1) * pe_xdif(iibm1 + ii_offset, ijbm1) + (ijbm1jp1 - ijbm1) * pe_ydif(iibm1, ijbm1 + ij_offset))
      IF (zey1 .LT. rsmall) zey1 = zey2
      IF (zey2 .LT. rsmall) zey2 = zey1
      zex1 = MAX(zex1, rsmall)
      zex2 = MAX(zex2, rsmall)
      zey1 = MAX(zey1, rsmall)
      zey2 = MAX(zey2, rsmall)
      zmask_x = (ABS(iibm1 - iibm2) * pmask_xdif(iibm2 + ii_offset, ijbm2) + ABS(ijbm1 - ijbm2) * pmask_ydif(iibm2, ijbm2 + ij_offset))
      zmask_y1 = ((iibm1 - iibm1jm1) * pmask_xdif(iibm1jm1 + ii_offset, ijbm1jm1) + (ijbm1 - ijbm1jm1) * pmask_ydif(iibm1jm1, ijbm1jm1 + ij_offset))
      zmask_y2 = ((iibm1jp1 - iibm1) * pmask_xdif(iibm1 + ii_offset, ijbm1) + (ijbm1jp1 - ijbm1) * pmask_ydif(iibm1, ijbm1 + ij_offset))
      zdt = phia(iibm1, ijbm1) - phib(iibm1, ijbm1)
      zdx = ((phia(iibm1, ijbm1) - phia(iibm2, ijbm2)) / zex2) * zmask_x
      zdy_1 = ((phib(iibm1, ijbm1) - phib(iibm1jm1, ijbm1jm1)) / zey1) * zmask_y1
      zdy_2 = ((phib(iibm1jp1, ijbm1jp1) - phib(iibm1, ijbm1)) / zey2) * zmask_y2
      zdy_centred = 0.5 * (zdy_1 + zdy_2)
      zsign_ups = SIGN(1., zdt * zdy_centred)
      zsign_ups = 0.5 * (zsign_ups + ABS(zsign_ups))
      zdy = zsign_ups * zdy_1 + (1. - zsign_ups) * zdy_2
      znor2 = zdx * zdx + zdy * zdy
      znor2 = MAX(znor2, zepsilon)
      zrx = zdt * zdx / (zex1 * znor2)
      zout = SIGN(1., zrx)
      zout = 0.5 * (zout + ABS(zout))
      zwgt = 2. * rdt * ((1. - zout) * idx % nbd(jb, igrd) + zout * idx % nbdout(jb, igrd))
      IF (ll_npo) THEN
        phia(ii, ij) = (1. - zout) * (phib(ii, ij) + zwgt * (phi_ext(jb) - phib(ii, ij))) + zout * (phib(ii, ij) + zrx * phia(iibm1, ijbm1) + zwgt * (phi_ext(jb) - phib(ii, ij))) / (1. + zrx)
      ELSE
        zsign_ups = SIGN(1., zdt * zdy)
        zsign_ups = 0.5 * (zsign_ups + ABS(zsign_ups))
        zey = zsign_ups * zey1 + (1. - zsign_ups) * zey2
        zry = zdt * zdy / (zey * znor2)
        phia(ii, ij) = (1. - zout) * (phib(ii, ij) + zwgt * (phi_ext(jb) - phib(ii, ij))) + zout * (phib(ii, ij) + zrx * phia(iibm1, ijbm1) - zsign_ups * zry * (phib(ii, ij) - phib(iijm1, ijjm1)) - (1. - zsign_ups) * zry * (phib(iijp1, ijjp1) - phib(ii, ij)) + zwgt * (phi_ext(jb) - phib(ii, ij))) / (1. + zrx)
      END IF
      phia(ii, ij) = phia(ii, ij) * pmask(ii, ij)
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE bdy_orlanski_2d
  SUBROUTINE bdy_orlanski_3d(idx, igrd, phib, phia, phi_ext, ll_npo)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    INTEGER, INTENT(IN) :: igrd
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: phib
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: phia
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: phi_ext
    LOGICAL, INTENT(IN) :: ll_npo
    INTEGER :: jb, jk
    INTEGER :: ii, ij, iibm1, iibm2, ijbm1, ijbm2
    INTEGER :: iijm1, iijp1, ijjm1, ijjp1
    INTEGER :: iibm1jp1, iibm1jm1, ijbm1jp1, ijbm1jm1
    INTEGER :: ii_offset, ij_offset
    INTEGER :: flagu, flagv
    REAL(KIND = wp) :: zmask_x, zmask_y1, zmask_y2
    REAL(KIND = wp) :: zex1, zex2, zey, zey1, zey2
    REAL(KIND = wp) :: zdt, zdx, zdy, znor2, zrx, zry
    REAL(KIND = wp) :: zout, zwgt, zdy_centred
    REAL(KIND = wp) :: zdy_1, zdy_2, zsign_ups
    REAL(KIND = wp), PARAMETER :: zepsilon = 1.E-30
    REAL(KIND = wp), POINTER, DIMENSION(:, :, :) :: pmask
    REAL(KIND = wp), POINTER, DIMENSION(:, :, :) :: pmask_xdif
    REAL(KIND = wp), POINTER, DIMENSION(:, :, :) :: pmask_ydif
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: pe_xdif
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: pe_ydif
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bdy_orlanski_3d', 'r0', 0, 0)
    SELECT CASE (igrd)
    CASE (1)
      pmask => tmask(:, :, :)
      pmask_xdif => umask(:, :, :)
      pmask_ydif => vmask(:, :, :)
      pe_xdif => e1u(:, :)
      pe_ydif => e2v(:, :)
      ii_offset = 0
      ij_offset = 0
    CASE (2)
      pmask => umask(:, :, :)
      pmask_xdif => tmask(:, :, :)
      pmask_ydif => fmask(:, :, :)
      pe_xdif => e1t(:, :)
      pe_ydif => e2f(:, :)
      ii_offset = 1
      ij_offset = 0
    CASE (3)
      pmask => vmask(:, :, :)
      pmask_xdif => fmask(:, :, :)
      pmask_ydif => tmask(:, :, :)
      pe_xdif => e1f(:, :)
      pe_ydif => e2t(:, :)
      ii_offset = 0
      ij_offset = 1
    CASE DEFAULT
      CALL ctl_stop('unrecognised value for igrd in bdy_orlanksi_2d')
    END SELECT
    DO jk = 1, jpk
      DO jb = 1, idx % nblenrim(igrd)
        ii = idx % nbi(jb, igrd)
        ij = idx % nbj(jb, igrd)
        flagu = INT(idx % flagu(jb, igrd))
        flagv = INT(idx % flagv(jb, igrd))
        iibm1 = ii + flagu
        iibm2 = ii + 2 * flagu
        ijbm1 = ij + flagv
        ijbm2 = ij + 2 * flagv
        iijm1 = ii - ABS(flagv)
        iijp1 = ii + ABS(flagv)
        ijjm1 = ij - ABS(flagu)
        ijjp1 = ij + ABS(flagu)
        iibm1jm1 = ii + flagu - ABS(flagv)
        iibm1jp1 = ii + flagu + ABS(flagv)
        ijbm1jm1 = ij + flagv - ABS(flagu)
        ijbm1jp1 = ij + flagv + ABS(flagu)
        zex1 = (ABS(iibm1 - iibm2) * pe_xdif(iibm1 + ii_offset, ijbm1) + ABS(ijbm1 - ijbm2) * pe_ydif(iibm1, ijbm1 + ij_offset))
        zex2 = (ABS(iibm1 - iibm2) * pe_xdif(iibm2 + ii_offset, ijbm2) + ABS(ijbm1 - ijbm2) * pe_ydif(iibm2, ijbm2 + ij_offset))
        zey1 = ((iibm1 - iibm1jm1) * pe_xdif(iibm1jm1 + ii_offset, ijbm1jm1) + (ijbm1 - ijbm1jm1) * pe_ydif(iibm1jm1, ijbm1jm1 + ij_offset))
        zey2 = ((iibm1jp1 - iibm1) * pe_xdif(iibm1 + ii_offset, ijbm1) + (ijbm1jp1 - ijbm1) * pe_ydif(iibm1, ijbm1 + ij_offset))
        IF (zey1 .LT. rsmall) zey1 = zey2
        IF (zey2 .LT. rsmall) zey2 = zey1
        zex1 = MAX(zex1, rsmall)
        zex2 = MAX(zex2, rsmall)
        zey1 = MAX(zey1, rsmall)
        zey2 = MAX(zey2, rsmall)
        zmask_x = (ABS(iibm1 - iibm2) * pmask_xdif(iibm2 + ii_offset, ijbm2, jk) + ABS(ijbm1 - ijbm2) * pmask_ydif(iibm2, ijbm2 + ij_offset, jk))
        zmask_y1 = ((iibm1 - iibm1jm1) * pmask_xdif(iibm1jm1 + ii_offset, ijbm1jm1, jk) + (ijbm1 - ijbm1jm1) * pmask_ydif(iibm1jm1, ijbm1jm1 + ij_offset, jk))
        zmask_y2 = ((iibm1jp1 - iibm1) * pmask_xdif(iibm1 + ii_offset, ijbm1, jk) + (ijbm1jp1 - ijbm1) * pmask_ydif(iibm1, ijbm1 + ij_offset, jk))
        zdt = phia(iibm1, ijbm1, jk) - phib(iibm1, ijbm1, jk)
        zdx = ((phia(iibm1, ijbm1, jk) - phia(iibm2, ijbm2, jk)) / zex2) * zmask_x
        zdy_1 = ((phib(iibm1, ijbm1, jk) - phib(iibm1jm1, ijbm1jm1, jk)) / zey1) * zmask_y1
        zdy_2 = ((phib(iibm1jp1, ijbm1jp1, jk) - phib(iibm1, ijbm1, jk)) / zey2) * zmask_y2
        zdy_centred = 0.5 * (zdy_1 + zdy_2)
        zsign_ups = SIGN(1., zdt * zdy_centred)
        zsign_ups = 0.5 * (zsign_ups + ABS(zsign_ups))
        zdy = zsign_ups * zdy_1 + (1. - zsign_ups) * zdy_2
        znor2 = zdx * zdx + zdy * zdy
        znor2 = MAX(znor2, zepsilon)
        zrx = zdt * zdx / (zex1 * znor2)
        zout = SIGN(1., zrx)
        zout = 0.5 * (zout + ABS(zout))
        zwgt = 2. * rdt * ((1. - zout) * idx % nbd(jb, igrd) + zout * idx % nbdout(jb, igrd))
        IF (ll_npo) THEN
          phia(ii, ij, jk) = (1. - zout) * (phib(ii, ij, jk) + zwgt * (phi_ext(jb, jk) - phib(ii, ij, jk))) + zout * (phib(ii, ij, jk) + zrx * phia(iibm1, ijbm1, jk) + zwgt * (phi_ext(jb, jk) - phib(ii, ij, jk))) / (1. + zrx)
        ELSE
          zsign_ups = SIGN(1., zdt * zdy)
          zsign_ups = 0.5 * (zsign_ups + ABS(zsign_ups))
          zey = zsign_ups * zey1 + (1. - zsign_ups) * zey2
          zry = zdt * zdy / (zey * znor2)
          phia(ii, ij, jk) = (1. - zout) * (phib(ii, ij, jk) + zwgt * (phi_ext(jb, jk) - phib(ii, ij, jk))) + zout * (phib(ii, ij, jk) + zrx * phia(iibm1, ijbm1, jk) - zsign_ups * zry * (phib(ii, ij, jk) - phib(iijm1, ijjm1, jk)) - (1. - zsign_ups) * zry * (phib(iijp1, ijjp1, jk) - phib(ii, ij, jk)) + zwgt * (phi_ext(jb, jk) - phib(ii, ij, jk))) / (1. + zrx)
        END IF
        phia(ii, ij, jk) = phia(ii, ij, jk) * pmask(ii, ij, jk)
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE bdy_orlanski_3d
  SUBROUTINE bdy_nmn(idx, igrd, phia)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: igrd
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: phia
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    REAL(KIND = wp) :: zcoef, zcoef1, zcoef2
    REAL(KIND = wp), POINTER, DIMENSION(:, :, :) :: pmask
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: bdypmask
    INTEGER :: ib, ik
    INTEGER :: ii, ij, ip, jp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bdy_nmn', 'r0', 0, 0)
    SELECT CASE (igrd)
    CASE (1)
      pmask => tmask(:, :, :)
      bdypmask => bdytmask(:, :)
    CASE (2)
      pmask => umask(:, :, :)
      bdypmask => bdyumask(:, :)
    CASE (3)
      pmask => vmask(:, :, :)
      bdypmask => bdyvmask(:, :)
    CASE DEFAULT
      CALL ctl_stop('unrecognised value for igrd in bdy_nmn')
    END SELECT
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    DO ib = 1, idx % nblenrim(igrd)
      ii = idx % nbi(ib, igrd)
      ij = idx % nbj(ib, igrd)
      DO ik = 1, jpkm1
        zcoef1 = bdypmask(ii - 1, ij) * pmask(ii - 1, ij, ik) + bdypmask(ii + 1, ij) * pmask(ii + 1, ij, ik)
        zcoef2 = bdypmask(ii, ij - 1) * pmask(ii, ij - 1, ik) + bdypmask(ii, ij + 1) * pmask(ii, ij + 1, ik)
        IF (NINT(zcoef1 + zcoef2) == 0) THEN
          zcoef = pmask(ii - 1, ij, ik) + pmask(ii + 1, ij, ik) + pmask(ii, ij - 1, ik) + pmask(ii, ij + 1, ik)
          IF (zcoef > .5_wp) THEN
            phia(ii, ij, ik) = phia(ii - 1, ij, ik) * pmask(ii - 1, ij, ik) + phia(ii + 1, ij, ik) * pmask(ii + 1, ij, ik) + phia(ii, ij - 1, ik) * pmask(ii, ij - 1, ik) + phia(ii, ij + 1, ik) * pmask(ii, ij + 1, ik)
            phia(ii, ij, ik) = (phia(ii, ij, ik) / zcoef) * pmask(ii, ij, ik)
          ELSE
            phia(ii, ij, ik) = phia(ii, ij, ik) * pmask(ii, ij, ik)
          END IF
        ELSE IF (NINT(zcoef1 + zcoef2) == 2) THEN
          zcoef = pmask(ii - 1, ij, ik) * bdypmask(ii - 1, ij) + pmask(ii + 1, ij, ik) * bdypmask(ii + 1, ij) + pmask(ii, ij - 1, ik) * bdypmask(ii, ij - 1) + pmask(ii, ij + 1, ik) * bdypmask(ii, ij + 1)
          phia(ii, ij, ik) = phia(ii - 1, ij, ik) * pmask(ii - 1, ij, ik) * bdypmask(ii - 1, ij) + phia(ii + 1, ij, ik) * pmask(ii + 1, ij, ik) * bdypmask(ii + 1, ij) + phia(ii, ij - 1, ik) * pmask(ii, ij - 1, ik) * bdypmask(ii, ij - 1) + phia(ii, ij + 1, ik) * pmask(ii, ij + 1, ik) * bdypmask(ii, ij + 1)
          phia(ii, ij, ik) = (phia(ii, ij, ik) / MAX(1._wp, zcoef)) * pmask(ii, ij, ik)
        ELSE
          ip = NINT(bdypmask(ii + 1, ij) * pmask(ii + 1, ij, ik) - bdypmask(ii - 1, ij) * pmask(ii - 1, ij, ik))
          jp = NINT(bdypmask(ii, ij + 1) * pmask(ii, ij + 1, ik) - bdypmask(ii, ij - 1) * pmask(ii, ij - 1, ik))
          phia(ii, ij, ik) = phia(ii + ip, ij + jp, ik) * pmask(ii + ip, ij + jp, ik)
        END IF
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE bdy_nmn
END MODULE bdylib