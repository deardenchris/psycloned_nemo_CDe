MODULE bdydyn2d
  USE oce
  USE dom_oce
  USE bdy_oce
  USE bdylib
  USE phycst
  USE lbclnk
  USE wet_dry
  USE in_out_manager
  USE lib_mpp, ONLY: ctl_stop
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: bdy_dyn2d
  PUBLIC :: bdy_ssh
  CONTAINS
  SUBROUTINE bdy_dyn2d(kt, pua2d, pva2d, pub2d, pvb2d, phur, phvr, pssh)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: pua2d, pva2d
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pub2d, pvb2d
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: phur, phvr
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pssh
    INTEGER :: ib_bdy
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bdy_dyn2d', 'r0', 0, 0)
    DO ib_bdy = 1, nb_bdy
      SELECT CASE (cn_dyn2d(ib_bdy))
      CASE ('none')
        CYCLE
      CASE ('frs')
        CALL bdy_dyn2d_frs(idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, pua2d, pva2d)
      CASE ('flather')
        CALL bdy_dyn2d_fla(idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, pua2d, pva2d, pssh, phur, phvr)
      CASE ('orlanski')
        CALL bdy_dyn2d_orlanski(idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, pua2d, pva2d, pub2d, pvb2d, ll_npo = .FALSE.)
      CASE ('orlanski_npo')
        CALL bdy_dyn2d_orlanski(idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, pua2d, pva2d, pub2d, pvb2d, ll_npo = .TRUE.)
      CASE DEFAULT
        CALL ctl_stop('bdy_dyn2d : unrecognised option for open boundaries for barotropic variables')
      END SELECT
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE bdy_dyn2d
  SUBROUTINE bdy_dyn2d_frs(idx, dta, ib_bdy, pua2d, pva2d)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    TYPE(OBC_DATA), INTENT(IN) :: dta
    INTEGER, INTENT(IN) :: ib_bdy
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: pua2d, pva2d
    INTEGER :: jb, jk
    INTEGER :: ii, ij, igrd
    REAL(KIND = wp) :: zwgt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bdy_dyn2d_frs', 'r0', 0, 0)
    igrd = 2
    DO jb = 1, idx % nblen(igrd)
      ii = idx % nbi(jb, igrd)
      ij = idx % nbj(jb, igrd)
      zwgt = idx % nbw(jb, igrd)
      pua2d(ii, ij) = (pua2d(ii, ij) + zwgt * (dta % u2d(jb) - pua2d(ii, ij))) * umask(ii, ij, 1)
    END DO
    igrd = 3
    DO jb = 1, idx % nblen(igrd)
      ii = idx % nbi(jb, igrd)
      ij = idx % nbj(jb, igrd)
      zwgt = idx % nbw(jb, igrd)
      pva2d(ii, ij) = (pva2d(ii, ij) + zwgt * (dta % v2d(jb) - pva2d(ii, ij))) * vmask(ii, ij, 1)
    END DO
    CALL lbc_bdy_lnk('bdydyn2d', pua2d, 'U', - 1., ib_bdy)
    CALL lbc_bdy_lnk('bdydyn2d', pva2d, 'V', - 1., ib_bdy)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE bdy_dyn2d_frs
  SUBROUTINE bdy_dyn2d_fla(idx, dta, ib_bdy, pua2d, pva2d, pssh, phur, phvr)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    TYPE(OBC_DATA), INTENT(IN) :: dta
    INTEGER, INTENT(IN) :: ib_bdy
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: pua2d, pva2d
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pssh, phur, phvr
    INTEGER :: jb, igrd
    INTEGER :: ii, ij, iim1, iip1, ijm1, ijp1
    REAL(KIND = wp), POINTER :: flagu, flagv
    REAL(KIND = wp) :: zcorr
    REAL(KIND = wp) :: zforc
    REAL(KIND = wp) :: zflag, z1_2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    z1_2 = 0.5_wp
    igrd = 1
    spgu(:, :) = 0.0
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('bdy_dyn2d_fla', 'r0', 0, 0)
    DO jb = 1, idx % nblenrim(igrd)
      ii = idx % nbi(jb, igrd)
      ij = idx % nbj(jb, igrd)
      IF (ll_wd) THEN
        spgu(ii, ij) = dta % ssh(jb) - ssh_ref
      ELSE
        spgu(ii, ij) = dta % ssh(jb)
      END IF
    END DO
    CALL lbc_bdy_lnk('bdydyn2d', spgu(:, :), 'T', 1., ib_bdy)
    igrd = 2
    DO jb = 1, idx % nblenrim(igrd)
      ii = idx % nbi(jb, igrd)
      ij = idx % nbj(jb, igrd)
      flagu => idx % flagu(jb, igrd)
      iim1 = ii + MAX(0, INT(flagu))
      iip1 = ii - MIN(0, INT(flagu))
      zcorr = - flagu * SQRT(grav * phur(ii, ij)) * (pssh(iim1, ij) - spgu(iip1, ij))
      zflag = ABS(flagu)
      zforc = dta % u2d(jb) * (1._wp - z1_2 * zflag) + z1_2 * zflag * pua2d(ii + NINT(flagu), ij)
      pua2d(ii, ij) = zforc + (1._wp - z1_2 * zflag) * zcorr * umask(ii, ij, 1)
    END DO
    igrd = 3
    DO jb = 1, idx % nblenrim(igrd)
      ii = idx % nbi(jb, igrd)
      ij = idx % nbj(jb, igrd)
      flagv => idx % flagv(jb, igrd)
      ijm1 = ij + MAX(0, INT(flagv))
      ijp1 = ij - MIN(0, INT(flagv))
      zcorr = - flagv * SQRT(grav * phvr(ii, ij)) * (pssh(ii, ijm1) - spgu(ii, ijp1))
      zflag = ABS(flagv)
      zforc = dta % v2d(jb) * (1._wp - z1_2 * zflag) + z1_2 * zflag * pva2d(ii, ij + NINT(flagv))
      pva2d(ii, ij) = zforc + (1._wp - z1_2 * zflag) * zcorr * vmask(ii, ij, 1)
    END DO
    CALL lbc_bdy_lnk('bdydyn2d', pua2d, 'U', - 1., ib_bdy)
    CALL lbc_bdy_lnk('bdydyn2d', pva2d, 'V', - 1., ib_bdy)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE bdy_dyn2d_fla
  SUBROUTINE bdy_dyn2d_orlanski(idx, dta, ib_bdy, pua2d, pva2d, pub2d, pvb2d, ll_npo)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    TYPE(OBC_DATA), INTENT(IN) :: dta
    INTEGER, INTENT(IN) :: ib_bdy
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: pua2d, pva2d
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pub2d, pvb2d
    LOGICAL, INTENT(IN) :: ll_npo
    INTEGER :: ib, igrd
    INTEGER :: ii, ij, iibm1, ijbm1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bdy_dyn2d_orlanski', 'r0', 0, 0)
    igrd = 2
    CALL bdy_orlanski_2d(idx, igrd, pub2d, pua2d, dta % u2d, ll_npo)
    igrd = 3
    CALL bdy_orlanski_2d(idx, igrd, pvb2d, pva2d, dta % v2d, ll_npo)
    CALL lbc_bdy_lnk('bdydyn2d', pua2d, 'U', - 1., ib_bdy)
    CALL lbc_bdy_lnk('bdydyn2d', pva2d, 'V', - 1., ib_bdy)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE bdy_dyn2d_orlanski
  SUBROUTINE bdy_ssh(zssh)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: zssh
    INTEGER :: ib_bdy, ib, igrd
    INTEGER :: ii, ij, zcoef, zcoef1, zcoef2, ip, jp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bdy_ssh', 'r0', 0, 0)
    igrd = 1
    DO ib_bdy = 1, nb_bdy
      DO ib = 1, idx_bdy(ib_bdy) % nblenrim(igrd)
        ii = idx_bdy(ib_bdy) % nbi(ib, igrd)
        ij = idx_bdy(ib_bdy) % nbj(ib, igrd)
        zcoef1 = bdytmask(ii - 1, ij) + bdytmask(ii + 1, ij)
        zcoef2 = bdytmask(ii, ij - 1) + bdytmask(ii, ij + 1)
        IF (zcoef1 + zcoef2 == 0) THEN
          zcoef = bdytmask(ii - 1, ij - 1) + bdytmask(ii + 1, ij + 1) + bdytmask(ii + 1, ij - 1) + bdytmask(ii - 1, ij + 1)
          zssh(ii, ij) = zssh(ii - 1, ij - 1) * bdytmask(ii - 1, ij - 1) + zssh(ii + 1, ij + 1) * bdytmask(ii + 1, ij + 1) + &
&zssh(ii + 1, ij - 1) * bdytmask(ii + 1, ij - 1) + zssh(ii - 1, ij + 1) * bdytmask(ii - 1, ij + 1)
          zssh(ii, ij) = (zssh(ii, ij) / MAX(1, zcoef)) * tmask(ii, ij, 1)
        ELSE
          ip = bdytmask(ii + 1, ij) - bdytmask(ii - 1, ij)
          jp = bdytmask(ii, ij + 1) - bdytmask(ii, ij - 1)
          zssh(ii, ij) = zssh(ii + ip, ij + jp) * tmask(ii + ip, ij + jp, 1)
        END IF
      END DO
      CALL lbc_bdy_lnk('bdydyn2d', zssh(:, :), 'T', 1., ib_bdy)
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE bdy_ssh
END MODULE bdydyn2d