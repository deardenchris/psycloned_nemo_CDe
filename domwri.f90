MODULE domwri
  USE dom_oce
  USE phycst, ONLY: rsmall
  USE wet_dry, ONLY: ll_wd
  USE in_out_manager
  USE iom
  USE lbclnk
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dom_wri
  PUBLIC :: dom_stiff
  CONTAINS
  SUBROUTINE dom_wri
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: inum
    CHARACTER(LEN = 21) :: clnam
    INTEGER :: ji, jj, jk
    INTEGER :: izco, izps, isco, icav
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zprt, zprw
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zdepu, zdepv
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    CALL profile_psy_data0 % PreStart('dom_wri', 'r0', 0, 0)
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) 'dom_wri : create NetCDF mesh and mask information file(s)'
    IF (lwp) WRITE(numout, *) '~~~~~~~'
    clnam = 'mesh_mask'
    CALL iom_open(TRIM(clnam), inum, ldwrt = .TRUE.)
    CALL iom_rstput(0, 0, inum, 'jpiglo', REAL(jpiglo, wp), ktype = jp_i4)
    CALL iom_rstput(0, 0, inum, 'jpjglo', REAL(jpjglo, wp), ktype = jp_i4)
    CALL iom_rstput(0, 0, inum, 'jpkglo', REAL(jpkglo, wp), ktype = jp_i4)
    CALL iom_rstput(0, 0, inum, 'jperio', REAL(jperio, wp), ktype = jp_i4)
    IF (ln_zco) THEN
      izco = 1
    ELSE
      izco = 0
    END IF
    IF (ln_zps) THEN
      izps = 1
    ELSE
      izps = 0
    END IF
    IF (ln_sco) THEN
      isco = 1
    ELSE
      isco = 0
    END IF
    CALL iom_rstput(0, 0, inum, 'ln_zco', REAL(izco, wp), ktype = jp_i4)
    CALL iom_rstput(0, 0, inum, 'ln_zps', REAL(izps, wp), ktype = jp_i4)
    CALL iom_rstput(0, 0, inum, 'ln_sco', REAL(isco, wp), ktype = jp_i4)
    IF (ln_isfcav) THEN
      icav = 1
    ELSE
      icav = 0
    END IF
    CALL iom_rstput(0, 0, inum, 'ln_isfcav', REAL(icav, wp), ktype = jp_i4)
    CALL iom_rstput(0, 0, inum, 'tmask', tmask, ktype = jp_i1)
    CALL iom_rstput(0, 0, inum, 'umask', umask, ktype = jp_i1)
    CALL iom_rstput(0, 0, inum, 'vmask', vmask, ktype = jp_i1)
    CALL iom_rstput(0, 0, inum, 'fmask', fmask, ktype = jp_i1)
    CALL dom_uniq(zprw, 'T')
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zprt(ji, jj) = ssmask(ji, jj) * zprw(ji, jj)
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('dom_wri', 'r1', 0, 0)
    CALL iom_rstput(0, 0, inum, 'tmaskutil', zprt, ktype = jp_i1)
    CALL dom_uniq(zprw, 'U')
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zprt(ji, jj) = ssumask(ji, jj) * zprw(ji, jj)
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('dom_wri', 'r2', 0, 0)
    CALL iom_rstput(0, 0, inum, 'umaskutil', zprt, ktype = jp_i1)
    CALL dom_uniq(zprw, 'V')
    CALL profile_psy_data2 % PostEnd
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zprt(ji, jj) = ssvmask(ji, jj) * zprw(ji, jj)
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data3 % PreStart('dom_wri', 'r3', 0, 0)
    CALL iom_rstput(0, 0, inum, 'vmaskutil', zprt, ktype = jp_i1)
    CALL iom_rstput(0, 0, inum, 'glamt', glamt, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'glamu', glamu, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'glamv', glamv, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'glamf', glamf, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'gphit', gphit, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'gphiu', gphiu, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'gphiv', gphiv, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'gphif', gphif, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e1t', e1t, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e1u', e1u, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e1v', e1v, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e1f', e1f, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e2t', e2t, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e2u', e2u, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e2v', e2v, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e2f', e2f, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'ff_f', ff_f, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'ff_t', ff_t, ktype = jp_r8)
    CALL profile_psy_data3 % PostEnd
    !$ACC KERNELS
    zprt(:, :) = ssmask(:, :) * REAL(mbkt(:, :), wp)
    !$ACC END KERNELS
    CALL profile_psy_data4 % PreStart('dom_wri', 'r4', 0, 0)
    CALL iom_rstput(0, 0, inum, 'mbathy', zprt, ktype = jp_i4)
    CALL profile_psy_data4 % PostEnd
    !$ACC KERNELS
    zprt(:, :) = ssmask(:, :) * REAL(mikt(:, :), wp)
    !$ACC END KERNELS
    CALL profile_psy_data5 % PreStart('dom_wri', 'r5', 0, 0)
    CALL iom_rstput(0, 0, inum, 'misf', zprt, ktype = jp_i4)
    CALL profile_psy_data5 % PostEnd
    !$ACC KERNELS
    zprt(:, :) = ssmask(:, :) * REAL(risfdep(:, :), wp)
    !$ACC END KERNELS
    CALL profile_psy_data6 % PreStart('dom_wri', 'r6', 0, 0)
    CALL iom_rstput(0, 0, inum, 'isfdraft', zprt, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e3t_0', e3t_0, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e3u_0', e3u_0, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e3v_0', e3v_0, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e3w_0', e3w_0, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'gdept_1d', gdept_1d, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'gdepw_1d', gdepw_1d, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'gdept_0', gdept_0, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'gdepw_0', gdepw_0, ktype = jp_r8)
    IF (ln_sco) THEN
      CALL dom_stiff(zprt)
      CALL iom_rstput(0, 0, inum, 'stiffness', zprt)
    END IF
    IF (ll_wd) CALL iom_rstput(0, 0, inum, 'ht_0', ht_0, ktype = jp_r8)
    CALL iom_close(inum)
    CALL profile_psy_data6 % PostEnd
  END SUBROUTINE dom_wri
  SUBROUTINE dom_uniq(puniq, cdgrd)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = 1), INTENT(IN) :: cdgrd
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: puniq
    REAL(KIND = wp) :: zshift
    INTEGER :: ji
    LOGICAL, DIMENSION(SIZE(puniq, 1), SIZE(puniq, 2), 1) :: lldbl
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztstref
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('dom_uniq', 'r0', 0, 0)
    zshift = jpi * jpj * (narea - 1)
    ztstref(:, :) = RESHAPE((/(zshift + REAL(ji, wp), ji = 1, jpi * jpj)/), (/jpi, jpj/))
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    puniq(:, :) = ztstref(:, :)
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('dom_uniq', 'r1', 0, 0)
    CALL lbc_lnk('domwri', puniq, cdgrd, 1.)
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    lldbl(:, :, 1) = puniq(:, :) == ztstref(:, :)
    puniq(:, :) = 1.
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('dom_uniq', 'r2', 0, 0)
    puniq(nldi : nlei, nldj : nlej) = REAL(COUNT(lldbl(nldi : nlei, nldj : nlej, :), dim = 3), wp)
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE dom_uniq
  SUBROUTINE dom_stiff(px1)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT), OPTIONAL :: px1
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zrxmax
    REAL(KIND = wp), DIMENSION(4) :: zr1
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zx1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    zx1(:, :) = 0._wp
    zrxmax = 0._wp
    zr1(:) = 0._wp
    !$ACC END KERNELS
    DO ji = 2, jpim1
      DO jj = 2, jpjm1
        !$ACC KERNELS
        DO jk = 1, jpkm1
          zr1(1) = ABS((gdepw_0(ji, jj, jk) - gdepw_0(ji - 1, jj, jk) + gdepw_0(ji, jj, jk + 1) - gdepw_0(ji - 1, jj, jk + 1)) / (gdepw_0(ji, jj, jk) + gdepw_0(ji - 1, jj, jk) - gdepw_0(ji, jj, jk + 1) - gdepw_0(ji - 1, jj, jk + 1) + rsmall)) * umask(ji - 1, jj, jk)
          zr1(2) = ABS((gdepw_0(ji + 1, jj, jk) - gdepw_0(ji, jj, jk) + gdepw_0(ji + 1, jj, jk + 1) - gdepw_0(ji, jj, jk + 1)) / (gdepw_0(ji + 1, jj, jk) + gdepw_0(ji, jj, jk) - gdepw_0(ji + 1, jj, jk + 1) - gdepw_0(ji, jj, jk + 1) + rsmall)) * umask(ji, jj, jk)
          zr1(3) = ABS((gdepw_0(ji, jj + 1, jk) - gdepw_0(ji, jj, jk) + gdepw_0(ji, jj + 1, jk + 1) - gdepw_0(ji, jj, jk + 1)) / (gdepw_0(ji, jj + 1, jk) + gdepw_0(ji, jj, jk) - gdepw_0(ji, jj + 1, jk + 1) - gdepw_0(ji, jj, jk + 1) + rsmall)) * vmask(ji, jj, jk)
          zr1(4) = ABS((gdepw_0(ji, jj, jk) - gdepw_0(ji, jj - 1, jk) + gdepw_0(ji, jj, jk + 1) - gdepw_0(ji, jj - 1, jk + 1)) / (gdepw_0(ji, jj, jk) + gdepw_0(ji, jj - 1, jk) - gdepw_0(ji, jj, jk + 1) - gdepw_0(ji, jj - 1, jk + 1) + rsmall)) * vmask(ji, jj - 1, jk)
          zrxmax = MAXVAL(zr1(1 : 4))
          zx1(ji, jj) = MAX(zx1(ji, jj), zrxmax)
        END DO
        !$ACC END KERNELS
      END DO
    END DO
    CALL profile_psy_data0 % PreStart('dom_stiff', 'r0', 0, 0)
    CALL lbc_lnk('domwri', zx1, 'T', 1.)
    IF (PRESENT(px1)) px1 = zx1
    zrxmax = MAXVAL(zx1)
    CALL mpp_max('domwri', zrxmax)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'dom_stiff : maximum grid stiffness ratio: ', zrxmax
      WRITE(numout, *) '~~~~~~~~~'
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE dom_stiff
END MODULE domwri