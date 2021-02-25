MODULE diacfl
  USE oce
  USE dom_oce
  USE domvvl
  USE lib_mpp
  USE lbclnk
  USE in_out_manager
  USE iom
  USE timing
  IMPLICIT NONE
  PRIVATE
  CHARACTER(LEN = 50) :: clname = "cfl_diagnostics.ascii"
  INTEGER :: numcfl
  INTEGER, DIMENSION(3) :: nCu_loc, nCv_loc, nCw_loc
  REAL(KIND = wp) :: rCu_max, rCv_max, rCw_max
  PUBLIC :: dia_cfl
  PUBLIC :: dia_cfl_init
  CONTAINS
  SUBROUTINE dia_cfl(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: z2dt, zCu_max, zCv_max, zCw_max
    INTEGER, DIMENSION(3) :: iloc_u, iloc_v, iloc_w, iloc
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zCu_cfl, zCv_cfl, zCw_cfl
    IF (ln_timing) CALL timing_start('dia_cfl')
    IF (neuler == 0 .AND. kt == nit000) THEN
      z2dt = rdt
    ELSE
      z2dt = rdt * 2._wp
    END IF
    !$OMP parallel default(shared), private(ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = 1, jpk
      DO jj = 1, jpj
        DO ji = 1, jpi
          zCu_cfl(ji, jj, jk) = ABS(un(ji, jj, jk)) * z2dt / e1u(ji, jj)
          zCv_cfl(ji, jj, jk) = ABS(vn(ji, jj, jk)) * z2dt / e2v(ji, jj)
          zCw_cfl(ji, jj, jk) = ABS(wn(ji, jj, jk)) * z2dt / e3w_n(ji, jj, jk)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    IF (iom_use('cfl_cu')) CALL iom_put('cfl_cu', MAXVAL(zCu_cfl, dim = 3))
    IF (iom_use('cfl_cv')) CALL iom_put('cfl_cv', MAXVAL(zCv_cfl, dim = 3))
    IF (iom_use('cfl_cw')) CALL iom_put('cfl_cw', MAXVAL(zCw_cfl, dim = 3))
    IF (lk_mpp) THEN
      CALL mpp_maxloc('diacfl', zCu_cfl, umask, zCu_max, iloc_u)
      CALL mpp_maxloc('diacfl', zCv_cfl, vmask, zCv_max, iloc_v)
      CALL mpp_maxloc('diacfl', zCw_cfl, wmask, zCw_max, iloc_w)
    ELSE
      iloc = MAXLOC(ABS(zcu_cfl(:, :, :)))
      iloc_u(1) = iloc(1) + nimpp - 1
      iloc_u(2) = iloc(2) + njmpp - 1
      iloc_u(3) = iloc(3)
      zCu_max = zCu_cfl(iloc(1), iloc(2), iloc(3))
      iloc = MAXLOC(ABS(zcv_cfl(:, :, :)))
      iloc_v(1) = iloc(1) + nimpp - 1
      iloc_v(2) = iloc(2) + njmpp - 1
      iloc_v(3) = iloc(3)
      zCv_max = zCv_cfl(iloc(1), iloc(2), iloc(3))
      iloc = MAXLOC(ABS(zcw_cfl(:, :, :)))
      iloc_w(1) = iloc(1) + nimpp - 1
      iloc_w(2) = iloc(2) + njmpp - 1
      iloc_w(3) = iloc(3)
      zCw_max = zCw_cfl(iloc(1), iloc(2), iloc(3))
    END IF
    IF (lwp) THEN
      WRITE(numcfl, FMT = '(2x,i6,3x,a6,4x,f7.4,1x,i4,1x,i4,1x,i4)') kt, 'Max Cu', zCu_max, iloc_u(1), iloc_u(2), iloc_u(3)
      WRITE(numcfl, FMT = '(11x,     a6,4x,f7.4,1x,i4,1x,i4,1x,i4)') 'Max Cv', zCv_max, iloc_v(1), iloc_v(2), iloc_v(3)
      WRITE(numcfl, FMT = '(11x,     a6,4x,f7.4,1x,i4,1x,i4,1x,i4)') 'Max Cw', zCw_max, iloc_w(1), iloc_w(2), iloc_w(3)
    END IF
    IF (zcu_max > rcu_max) THEN
      rcu_max = zcu_max
      ncu_loc(:) = iloc_u(:)
    END IF
    IF (zcv_max > rcv_max) THEN
      rcv_max = zcv_max
      ncv_loc(:) = iloc_v(:)
    END IF
    IF (zcw_max > rcw_max) THEN
      rcw_max = zcw_max
      ncw_loc(:) = iloc_w(:)
    END IF
    IF (kt == nitend .AND. lwp) THEN
      WRITE(numcfl, FMT = *) '******************************************'
      WRITE(numcfl, FMT = '(3x,a12,6x,f7.4,1x,i4,1x,i4,1x,i4)') 'Run Max Cu', rCu_max, nCu_loc(1), nCu_loc(2), nCu_loc(3)
      WRITE(numcfl, FMT = '(3x,a8,11x,f15.1)') ' => dt/C', z2dt / rCu_max
      WRITE(numcfl, FMT = *) '******************************************'
      WRITE(numcfl, FMT = '(3x,a12,6x,f7.4,1x,i4,1x,i4,1x,i4)') 'Run Max Cv', rCv_max, nCv_loc(1), nCv_loc(2), nCv_loc(3)
      WRITE(numcfl, FMT = '(3x,a8,11x,f15.1)') ' => dt/C', z2dt / rCv_max
      WRITE(numcfl, FMT = *) '******************************************'
      WRITE(numcfl, FMT = '(3x,a12,6x,f7.4,1x,i4,1x,i4,1x,i4)') 'Run Max Cw', rCw_max, nCw_loc(1), nCw_loc(2), nCw_loc(3)
      WRITE(numcfl, FMT = '(3x,a8,11x,f15.1)') ' => dt/C', z2dt / rCw_max
      CLOSE(UNIT = numcfl)
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dia_cfl : Maximum Courant number information for the run '
      WRITE(numout, FMT = *) '~~~~~~~'
      WRITE(numout, FMT = *) '   Max Cu = ', rCu_max, ' at (i,j,k) = (', nCu_loc(1), nCu_loc(2), nCu_loc(3), ') => dt/C = ', z2dt &
&/ rCu_max
      WRITE(numout, FMT = *) '   Max Cv = ', rCv_max, ' at (i,j,k) = (', nCv_loc(1), nCv_loc(2), nCv_loc(3), ') => dt/C = ', z2dt &
&/ rCv_max
      WRITE(numout, FMT = *) '   Max Cw = ', rCw_max, ' at (i,j,k) = (', nCw_loc(1), nCw_loc(2), nCw_loc(3), ') => dt/C = ', z2dt &
&/ rCw_max
    END IF
    IF (ln_timing) CALL timing_stop('dia_cfl')
  END SUBROUTINE dia_cfl
  SUBROUTINE dia_cfl_init
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dia_cfl : Outputting CFL diagnostics to ', TRIM(clname), ' file'
      WRITE(numout, FMT = *) '~~~~~~~'
      WRITE(numout, FMT = *)
      CALL ctl_opn(numcfl, clname, 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', 1, numout, lwp, 1)
      WRITE(numcfl, FMT = *) 'Timestep  Direction  Max C     i    j    k'
      WRITE(numcfl, FMT = *) '******************************************'
    END IF
    rCu_max = 0._wp
    rCv_max = 0._wp
    rCw_max = 0._wp
  END SUBROUTINE dia_cfl_init
END MODULE diacfl