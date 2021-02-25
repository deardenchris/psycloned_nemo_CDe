MODULE bdydyn3d
  USE timing
  USE oce
  USE dom_oce
  USE bdy_oce
  USE bdylib
  USE lbclnk
  USE in_out_manager
  USE lib_mpp, ONLY: ctl_stop
  USE phycst
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: bdy_dyn3d
  PUBLIC :: bdy_dyn3d_dmp
  CONTAINS
  SUBROUTINE bdy_dyn3d(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ib_bdy
    DO ib_bdy = 1, nb_bdy
      SELECT CASE (cn_dyn3d(ib_bdy))
      CASE ('none')
        CYCLE
      CASE ('frs')
        CALL bdy_dyn3d_frs(idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt, ib_bdy)
      CASE ('specified')
        CALL bdy_dyn3d_spe(idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt, ib_bdy)
      CASE ('zero')
        CALL bdy_dyn3d_zro(idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt, ib_bdy)
      CASE ('orlanski')
        CALL bdy_dyn3d_orlanski(idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, ll_npo = .FALSE.)
      CASE ('orlanski_npo')
        CALL bdy_dyn3d_orlanski(idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, ll_npo = .TRUE.)
      CASE ('zerograd')
        CALL bdy_dyn3d_zgrad(idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt, ib_bdy)
      CASE ('neumann')
        CALL bdy_dyn3d_nmn(idx_bdy(ib_bdy), ib_bdy)
      CASE DEFAULT
        CALL ctl_stop('bdy_dyn3d : unrecognised option for open boundaries for baroclinic velocities')
      END SELECT
    END DO
  END SUBROUTINE bdy_dyn3d
  SUBROUTINE bdy_dyn3d_spe(idx, dta, kt, ib_bdy)
    INTEGER, INTENT(IN) :: kt
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    TYPE(OBC_DATA), INTENT(IN) :: dta
    INTEGER, INTENT(IN) :: ib_bdy
    INTEGER :: jb, jk
    INTEGER :: ii, ij, igrd
    REAL(KIND = wp) :: zwgt
    igrd = 2
    DO jb = 1, idx % nblenrim(igrd)
      DO jk = 1, jpkm1
        ii = idx % nbi(jb, igrd)
        ij = idx % nbj(jb, igrd)
        ua(ii, ij, jk) = dta % u3d(jb, jk) * umask(ii, ij, jk)
      END DO
    END DO
    igrd = 3
    DO jb = 1, idx % nblenrim(igrd)
      DO jk = 1, jpkm1
        ii = idx % nbi(jb, igrd)
        ij = idx % nbj(jb, igrd)
        va(ii, ij, jk) = dta % v3d(jb, jk) * vmask(ii, ij, jk)
      END DO
    END DO
    CALL lbc_bdy_lnk('bdydyn3d', ua, 'U', - 1., ib_bdy)
    CALL lbc_bdy_lnk('bdydyn3d', va, 'V', - 1., ib_bdy)
    IF (kt == nit000) CLOSE(UNIT = 102)
  END SUBROUTINE bdy_dyn3d_spe
  SUBROUTINE bdy_dyn3d_zgrad(idx, dta, kt, ib_bdy)
    INTEGER :: kt
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    TYPE(OBC_DATA), INTENT(IN) :: dta
    INTEGER, INTENT(IN) :: ib_bdy
    INTEGER :: jb, jk
    INTEGER :: ii, ij, igrd
    REAL(KIND = wp) :: zwgt
    INTEGER :: fu, fv
    igrd = 2
    DO jb = 1, idx % nblenrim(igrd)
      DO jk = 1, jpkm1
        ii = idx % nbi(jb, igrd)
        ij = idx % nbj(jb, igrd)
        fu = ABS(ABS(NINT(idx % flagu(jb, igrd))) - 1)
        ua(ii, ij, jk) = ua(ii, ij, jk) * REAL(1 - fu) + (ua(ii, ij + fu, jk) * umask(ii, ij + fu, jk) + ua(ii, ij - fu, jk) * &
&umask(ii, ij - fu, jk)) * umask(ii, ij, jk) * REAL(fu)
      END DO
    END DO
    igrd = 3
    DO jb = 1, idx % nblenrim(igrd)
      DO jk = 1, jpkm1
        ii = idx % nbi(jb, igrd)
        ij = idx % nbj(jb, igrd)
        fv = ABS(ABS(NINT(idx % flagv(jb, igrd))) - 1)
        va(ii, ij, jk) = va(ii, ij, jk) * REAL(1 - fv) + (va(ii + fv, ij, jk) * vmask(ii + fv, ij, jk) + va(ii - fv, ij, jk) * &
&vmask(ii - fv, ij, jk)) * vmask(ii, ij, jk) * REAL(fv)
      END DO
    END DO
    CALL lbc_bdy_lnk('bdydyn3d', ua, 'U', - 1., ib_bdy)
    CALL lbc_bdy_lnk('bdydyn3d', va, 'V', - 1., ib_bdy)
    IF (kt == nit000) CLOSE(UNIT = 102)
  END SUBROUTINE bdy_dyn3d_zgrad
  SUBROUTINE bdy_dyn3d_zro(idx, dta, kt, ib_bdy)
    INTEGER, INTENT(IN) :: kt
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    TYPE(OBC_DATA), INTENT(IN) :: dta
    INTEGER, INTENT(IN) :: ib_bdy
    INTEGER :: ib, ik
    INTEGER :: ii, ij, igrd
    REAL(KIND = wp) :: zwgt
    igrd = 2
    DO ib = 1, idx % nblenrim(igrd)
      ii = idx % nbi(ib, igrd)
      ij = idx % nbj(ib, igrd)
      DO ik = 1, jpkm1
        ua(ii, ij, ik) = 0._wp
      END DO
    END DO
    igrd = 3
    DO ib = 1, idx % nblenrim(igrd)
      ii = idx % nbi(ib, igrd)
      ij = idx % nbj(ib, igrd)
      DO ik = 1, jpkm1
        va(ii, ij, ik) = 0._wp
      END DO
    END DO
    CALL lbc_bdy_lnk('bdydyn3d', ua, 'U', - 1., ib_bdy)
    CALL lbc_bdy_lnk('bdydyn3d', va, 'V', - 1., ib_bdy)
    IF (kt == nit000) CLOSE(UNIT = 102)
  END SUBROUTINE bdy_dyn3d_zro
  SUBROUTINE bdy_dyn3d_frs(idx, dta, kt, ib_bdy)
    INTEGER, INTENT(IN) :: kt
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    TYPE(OBC_DATA), INTENT(IN) :: dta
    INTEGER, INTENT(IN) :: ib_bdy
    INTEGER :: jb, jk
    INTEGER :: ii, ij, igrd
    REAL(KIND = wp) :: zwgt
    igrd = 2
    DO jb = 1, idx % nblen(igrd)
      DO jk = 1, jpkm1
        ii = idx % nbi(jb, igrd)
        ij = idx % nbj(jb, igrd)
        zwgt = idx % nbw(jb, igrd)
        ua(ii, ij, jk) = (ua(ii, ij, jk) + zwgt * (dta % u3d(jb, jk) - ua(ii, ij, jk))) * umask(ii, ij, jk)
      END DO
    END DO
    igrd = 3
    DO jb = 1, idx % nblen(igrd)
      DO jk = 1, jpkm1
        ii = idx % nbi(jb, igrd)
        ij = idx % nbj(jb, igrd)
        zwgt = idx % nbw(jb, igrd)
        va(ii, ij, jk) = (va(ii, ij, jk) + zwgt * (dta % v3d(jb, jk) - va(ii, ij, jk))) * vmask(ii, ij, jk)
      END DO
    END DO
    CALL lbc_bdy_lnk('bdydyn3d', ua, 'U', - 1., ib_bdy)
    CALL lbc_bdy_lnk('bdydyn3d', va, 'V', - 1., ib_bdy)
    IF (kt == nit000) CLOSE(UNIT = 102)
  END SUBROUTINE bdy_dyn3d_frs
  SUBROUTINE bdy_dyn3d_orlanski(idx, dta, ib_bdy, ll_npo)
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    TYPE(OBC_DATA), INTENT(IN) :: dta
    INTEGER, INTENT(IN) :: ib_bdy
    LOGICAL, INTENT(IN) :: ll_npo
    INTEGER :: jb, igrd
    igrd = 2
    CALL bdy_orlanski_3d(idx, igrd, ub, ua, dta % u3d, ll_npo)
    igrd = 3
    CALL bdy_orlanski_3d(idx, igrd, vb, va, dta % v3d, ll_npo)
    CALL lbc_bdy_lnk('bdydyn3d', ua, 'U', - 1., ib_bdy)
    CALL lbc_bdy_lnk('bdydyn3d', va, 'V', - 1., ib_bdy)
  END SUBROUTINE bdy_dyn3d_orlanski
  SUBROUTINE bdy_dyn3d_dmp(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: jb, jk
    INTEGER :: ib_bdy
    INTEGER :: ii, ij, igrd
    REAL(KIND = wp) :: zwgt
    IF (ln_timing) CALL timing_start('bdy_dyn3d_dmp')
    DO ib_bdy = 1, nb_bdy
      IF (ln_dyn3d_dmp(ib_bdy) .AND. cn_dyn3d(ib_bdy) /= 'none') THEN
        igrd = 2
        DO jb = 1, idx_bdy(ib_bdy) % nblen(igrd)
          ii = idx_bdy(ib_bdy) % nbi(jb, igrd)
          ij = idx_bdy(ib_bdy) % nbj(jb, igrd)
          zwgt = idx_bdy(ib_bdy) % nbd(jb, igrd)
          DO jk = 1, jpkm1
            ua(ii, ij, jk) = (ua(ii, ij, jk) + zwgt * (dta_bdy(ib_bdy) % u3d(jb, jk) - ub(ii, ij, jk) + ub_b(ii, ij))) * umask(ii, &
&ij, jk)
          END DO
        END DO
        igrd = 3
        DO jb = 1, idx_bdy(ib_bdy) % nblen(igrd)
          ii = idx_bdy(ib_bdy) % nbi(jb, igrd)
          ij = idx_bdy(ib_bdy) % nbj(jb, igrd)
          zwgt = idx_bdy(ib_bdy) % nbd(jb, igrd)
          DO jk = 1, jpkm1
            va(ii, ij, jk) = (va(ii, ij, jk) + zwgt * (dta_bdy(ib_bdy) % v3d(jb, jk) - vb(ii, ij, jk) + vb_b(ii, ij))) * vmask(ii, &
&ij, jk)
          END DO
        END DO
      END IF
    END DO
    CALL lbc_lnk_multi('bdydyn3d', ua, 'U', - 1., va, 'V', - 1.)
    IF (ln_timing) CALL timing_stop('bdy_dyn3d_dmp')
  END SUBROUTINE bdy_dyn3d_dmp
  SUBROUTINE bdy_dyn3d_nmn(idx, ib_bdy)
    TYPE(OBC_INDEX), INTENT(IN) :: idx
    INTEGER, INTENT(IN) :: ib_bdy
    INTEGER :: jb, igrd
    igrd = 2
    CALL bdy_nmn(idx, igrd, ua)
    igrd = 3
    CALL bdy_nmn(idx, igrd, va)
    CALL lbc_bdy_lnk('bdydyn3d', ua, 'U', - 1., ib_bdy)
    CALL lbc_bdy_lnk('bdydyn3d', va, 'V', - 1., ib_bdy)
  END SUBROUTINE bdy_dyn3d_nmn
END MODULE bdydyn3d