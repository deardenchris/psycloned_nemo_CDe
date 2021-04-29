MODULE bdyvol
  USE oce
  USE bdy_oce
  USE sbc_oce
  USE dom_oce
  USE phycst
  USE sbcisf
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: bdy_vol2d
  CONTAINS
  SUBROUTINE bdy_vol2d(kt, kc, pua2d, pva2d, phu, phv)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt, kc
    INTEGER :: ji, jj, jk, jb, jgrd
    INTEGER :: ib_bdy, ii, ij
    REAL(KIND = wp) :: zubtpecor, ztranst
    REAL(KIND = wp), SAVE :: z_cflxemp
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: pua2d, pva2d
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: phu, phv
    TYPE(OBC_INDEX), POINTER :: idx
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bdy_vol2d', 'r0', 0, 0)
    IF (kc == 1) z_cflxemp = glob_sum('bdyvol', (emp(:, :) - rnf(:, :) + fwfisf(:, :)) * bdytmask(:, :) * e1e2t(:, :)) / rau0
    IF (.NOT. ln_linssh) THEN
      bdysurftot = bdy_segs_surf(phu, phv)
    ELSE
      IF ((kt == nit000) .AND. (kc == 1)) bdysurftot = bdy_segs_surf(phu, phv)
    END IF
    zubtpecor = 0._wp
    DO ib_bdy = 1, nb_bdy
      idx => idx_bdy(ib_bdy)
      jgrd = 2
      DO jb = 1, idx % nblenrim(jgrd)
        ii = idx % nbi(jb, jgrd)
        ij = idx % nbj(jb, jgrd)
        zubtpecor = zubtpecor + idx % flagu(jb, jgrd) * pua2d(ii, ij) * e2u(ii, ij) * phu(ii, ij) * tmask_i(ii, ij) * tmask_i(ii + 1, ij)
      END DO
      jgrd = 3
      DO jb = 1, idx % nblenrim(jgrd)
        ii = idx % nbi(jb, jgrd)
        ij = idx % nbj(jb, jgrd)
        zubtpecor = zubtpecor + idx % flagv(jb, jgrd) * pva2d(ii, ij) * e1v(ii, ij) * phv(ii, ij) * tmask_i(ii, ij) * tmask_i(ii, ij + 1)
      END DO
    END DO
    IF (lk_mpp) CALL mpp_sum('bdyvol', zubtpecor)
    IF (nn_volctl == 1) THEN
      zubtpecor = (zubtpecor - z_cflxemp) / bdysurftot
    ELSE
      zubtpecor = zubtpecor / bdysurftot
    END IF
    DO ib_bdy = 1, nb_bdy
      idx => idx_bdy(ib_bdy)
      jgrd = 2
      DO jb = 1, idx % nblenrim(jgrd)
        ii = idx % nbi(jb, jgrd)
        ij = idx % nbj(jb, jgrd)
        pua2d(ii, ij) = pua2d(ii, ij) - idx % flagu(jb, jgrd) * zubtpecor * tmask_i(ii, ij) * tmask_i(ii + 1, ij)
      END DO
      jgrd = 3
      DO jb = 1, idx % nblenrim(jgrd)
        ii = idx % nbi(jb, jgrd)
        ij = idx % nbj(jb, jgrd)
        pva2d(ii, ij) = pva2d(ii, ij) - idx % flagv(jb, jgrd) * zubtpecor * tmask_i(ii, ij) * tmask_i(ii, ij + 1)
      END DO
    END DO
    IF (MOD(kt, nwrite) == 0 .AND. (kc == 1)) THEN
      ztranst = 0._wp
      DO ib_bdy = 1, nb_bdy
        idx => idx_bdy(ib_bdy)
        jgrd = 2
        DO jb = 1, idx % nblenrim(jgrd)
          ii = idx % nbi(jb, jgrd)
          ij = idx % nbj(jb, jgrd)
          ztranst = ztranst + idx % flagu(jb, jgrd) * pua2d(ii, ij) * e2u(ii, ij) * phu(ii, ij) * tmask_i(ii, ij) * tmask_i(ii + 1, ij)
        END DO
        jgrd = 3
        DO jb = 1, idx % nblenrim(jgrd)
          ii = idx % nbi(jb, jgrd)
          ij = idx % nbj(jb, jgrd)
          ztranst = ztranst + idx % flagv(jb, jgrd) * pva2d(ii, ij) * e1v(ii, ij) * phv(ii, ij) * tmask_i(ii, ij) * tmask_i(ii, ij + 1)
        END DO
      END DO
      IF (lk_mpp) CALL mpp_sum('bdyvol', ztranst)
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'bdy_vol : time step :', kt
      IF (lwp) WRITE(numout, *) '~~~~~~~ '
      IF (lwp) WRITE(numout, *) '          cumulate flux EMP             =', z_cflxemp, ' (m3/s)'
      IF (lwp) WRITE(numout, *) '          total lateral surface of OBC  =', bdysurftot, '(m2)'
      IF (lwp) WRITE(numout, *) '          correction velocity zubtpecor =', zubtpecor, '(m/s)'
      IF (lwp) WRITE(numout, *) '          cumulated transport ztranst   =', ztranst, '(m3/s)'
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE bdy_vol2d
  REAL(KIND = wp) FUNCTION bdy_segs_surf(phu, phv)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: phu, phv
    INTEGER :: igrd, ib_bdy, ib
    INTEGER, POINTER :: nbi, nbj
    REAL(KIND = wp), POINTER :: zflagu, zflagv
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wp', 'r0', 0, 0)
    bdy_segs_surf = 0._wp
    igrd = 2
    DO ib_bdy = 1, nb_bdy
      DO ib = 1, idx_bdy(ib_bdy) % nblenrim(igrd)
        nbi => idx_bdy(ib_bdy) % nbi(ib, igrd)
        nbj => idx_bdy(ib_bdy) % nbj(ib, igrd)
        zflagu => idx_bdy(ib_bdy) % flagu(ib, igrd)
        bdy_segs_surf = bdy_segs_surf + phu(nbi, nbj) * e2u(nbi, nbj) * ABS(zflagu) * tmask_i(nbi, nbj) * tmask_i(nbi + 1, nbj)
      END DO
    END DO
    igrd = 3
    DO ib_bdy = 1, nb_bdy
      DO ib = 1, idx_bdy(ib_bdy) % nblenrim(igrd)
        nbi => idx_bdy(ib_bdy) % nbi(ib, igrd)
        nbj => idx_bdy(ib_bdy) % nbj(ib, igrd)
        zflagv => idx_bdy(ib_bdy) % flagv(ib, igrd)
        bdy_segs_surf = bdy_segs_surf + phv(nbi, nbj) * e1v(nbi, nbj) * ABS(zflagv) * tmask_i(nbi, nbj) * tmask_i(nbi, nbj + 1)
      END DO
    END DO
    IF (lk_mpp) CALL mpp_sum('bdyvol', bdy_segs_surf)
    CALL profile_psy_data0 % PostEnd
  END FUNCTION bdy_segs_surf
END MODULE bdyvol