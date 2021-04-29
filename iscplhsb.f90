MODULE iscplhsb
  USE oce
  USE dom_oce
  USE domwri
  USE domngb
  USE phycst
  USE sbc_oce
  USE iscplini
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: iscpl_div
  PUBLIC :: iscpl_cons
  CONTAINS
  SUBROUTINE iscpl_cons(ptmask_b, psmask_b, pe3t_b, pts_flx, pvol_flx, prdt_iscpl)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: ptmask_b
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pe3t_b
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: psmask_b
    REAL(KIND = wp), DIMENSION(:, :, :, :), INTENT(OUT) :: pts_flx
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: pvol_flx
    REAL(KIND = wp), INTENT(IN) :: prdt_iscpl
    INTEGER :: ji, jj, jk
    INTEGER :: jip1, jim1, jjp1, jjm1
    REAL(KIND = wp) :: summsk, zsum, zsumn, zjip1_ratio, zjim1_ratio, zdtem, zde3t, z1_rdtiscpl
    REAL(KIND = wp) :: zarea, zsum1, zsumb, zjjp1_ratio, zjjm1_ratio, zdsal
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdssh
    REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: zlon, zlat
    REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: zcorr_vol, zcorr_tem, zcorr_sal
    INTEGER, DIMENSION(:), ALLOCATABLE :: ixpts, iypts, izpts, inpts
    INTEGER :: jpts, npts
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    !$ACC KERNELS
    zde3t = 0._wp
    zdsal = 0._wp
    zdtem = 0._wp
    pvol_flx(:, :, :) = 0._wp
    pts_flx(:, :, :, :) = 0._wp
    z1_rdtiscpl = 1._wp / prdt_iscpl
    tsb(:, :, :, jp_tem) = tsb(:, :, :, jp_tem) * ptmask_b(:, :, :)
    tsn(:, :, :, jp_tem) = tsn(:, :, :, jp_tem) * tmask(:, :, :)
    tsb(:, :, :, jp_sal) = tsb(:, :, :, jp_sal) * ptmask_b(:, :, :)
    tsn(:, :, :, jp_sal) = tsn(:, :, :, jp_sal) * tmask(:, :, :)
    zdssh(:, :) = sshn(:, :) * ssmask(:, :) - sshb(:, :) * psmask_b(:, :)
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('iscpl_cons', 'r0', 0, 0)
    IF (.NOT. ln_linssh) zdssh = 0.0_wp
    DO jk = 1, jpk - 1
      DO jj = 2, jpj - 1
        DO ji = 2, jpim1
          IF (tmask_h(ji, jj) == 1._wp) THEN
            zde3t = e3t_n(ji, jj, jk) * tmask(ji, jj, jk) - pe3t_b(ji, jj, jk) * ptmask_b(ji, jj, jk)
            zdtem = tsn(ji, jj, jk, jp_tem) * e3t_n(ji, jj, jk) * tmask(ji, jj, jk) - tsb(ji, jj, jk, jp_tem) * pe3t_b(ji, jj, jk) * ptmask_b(ji, jj, jk)
            zdsal = tsn(ji, jj, jk, jp_sal) * e3t_n(ji, jj, jk) * tmask(ji, jj, jk) - tsb(ji, jj, jk, jp_sal) * pe3t_b(ji, jj, jk) * ptmask_b(ji, jj, jk)
            IF (ptmask_b(ji, jj, jk) == 1._wp .OR. tmask(ji, jj, jk) == 1._wp) THEN
              zde3t = zde3t + zdssh(ji, jj)
              zdssh(ji, jj) = 0._wp
            END IF
            pvol_flx(ji, jj, jk) = pvol_flx(ji, jj, jk) + zde3t * z1_rdtiscpl
            pts_flx(ji, jj, jk, jp_sal) = pts_flx(ji, jj, jk, jp_sal) + zdsal * z1_rdtiscpl
            pts_flx(ji, jj, jk, jp_tem) = pts_flx(ji, jj, jk, jp_tem) + zdtem * z1_rdtiscpl
            IF (tmask(ji, jj, jk) == 0._wp .AND. ptmask_b(ji, jj, jk) == 1._wp) THEN
              jip1 = ji + 1
              jim1 = ji - 1
              jjp1 = jj + 1
              jjm1 = jj - 1
              zsum = e1e2t(ji, jjp1) * tmask(ji, jjp1, jk) + e1e2t(ji, jjm1) * tmask(ji, jjm1, jk) + e1e2t(jim1, jj) * tmask(jim1, jj, jk) + e1e2t(jip1, jj) * tmask(jip1, jj, jk)
              IF (zsum /= 0._wp) THEN
                zjip1_ratio = e1e2t(jip1, jj) * tmask(jip1, jj, jk) / zsum
                zjim1_ratio = e1e2t(jim1, jj) * tmask(jim1, jj, jk) / zsum
                zjjp1_ratio = e1e2t(ji, jjp1) * tmask(ji, jjp1, jk) / zsum
                zjjm1_ratio = e1e2t(ji, jjm1) * tmask(ji, jjm1, jk) / zsum
                pvol_flx(ji, jjp1, jk) = pvol_flx(ji, jjp1, jk) + pvol_flx(ji, jj, jk) * zjjp1_ratio
                pvol_flx(ji, jjm1, jk) = pvol_flx(ji, jjm1, jk) + pvol_flx(ji, jj, jk) * zjjm1_ratio
                pvol_flx(jip1, jj, jk) = pvol_flx(jip1, jj, jk) + pvol_flx(ji, jj, jk) * zjip1_ratio
                pvol_flx(jim1, jj, jk) = pvol_flx(jim1, jj, jk) + pvol_flx(ji, jj, jk) * zjim1_ratio
                pts_flx(ji, jjp1, jk, jp_sal) = pts_flx(ji, jjp1, jk, jp_sal) + pts_flx(ji, jj, jk, jp_sal) * zjjp1_ratio
                pts_flx(ji, jjm1, jk, jp_sal) = pts_flx(ji, jjm1, jk, jp_sal) + pts_flx(ji, jj, jk, jp_sal) * zjjm1_ratio
                pts_flx(jip1, jj, jk, jp_sal) = pts_flx(jip1, jj, jk, jp_sal) + pts_flx(ji, jj, jk, jp_sal) * zjip1_ratio
                pts_flx(jim1, jj, jk, jp_sal) = pts_flx(jim1, jj, jk, jp_sal) + pts_flx(ji, jj, jk, jp_sal) * zjim1_ratio
                pts_flx(ji, jjp1, jk, jp_tem) = pts_flx(ji, jjp1, jk, jp_tem) + pts_flx(ji, jj, jk, jp_tem) * zjjp1_ratio
                pts_flx(ji, jjm1, jk, jp_tem) = pts_flx(ji, jjm1, jk, jp_tem) + pts_flx(ji, jj, jk, jp_tem) * zjjm1_ratio
                pts_flx(jip1, jj, jk, jp_tem) = pts_flx(jip1, jj, jk, jp_tem) + pts_flx(ji, jj, jk, jp_tem) * zjip1_ratio
                pts_flx(jim1, jj, jk, jp_tem) = pts_flx(jim1, jj, jk, jp_tem) + pts_flx(ji, jj, jk, jp_tem) * zjim1_ratio
                pvol_flx(ji, jj, jk) = 0._wp
                pts_flx(ji, jj, jk, jp_sal) = 0._wp
                pts_flx(ji, jj, jk, jp_tem) = 0._wp
              ELSE IF (zsum == 0._wp) THEN
                IF (tmask(ji, jj, jk + 1) == 1._wp) THEN
                  pvol_flx(ji, jj, jk + 1) = pvol_flx(ji, jj, jk + 1) + pvol_flx(ji, jj, jk)
                  pts_flx(ji, jj, jk + 1, jp_sal) = pts_flx(ji, jj, jk + 1, jp_sal) + pts_flx(ji, jj, jk, jp_sal)
                  pts_flx(ji, jj, jk + 1, jp_tem) = pts_flx(ji, jj, jk + 1, jp_tem) + pts_flx(ji, jj, jk, jp_tem)
                  pvol_flx(ji, jj, jk) = 0._wp
                  pts_flx(ji, jj, jk, jp_sal) = 0._wp
                  pts_flx(ji, jj, jk, jp_tem) = 0._wp
                ELSE
                  IF (lwp) THEN
                    WRITE(numout, *) 'W A R N I N G iscpl: no adjacent cell on the vertical and horizontal'
                    WRITE(numout, *) '                     ', mig(ji), ' ', mjg(jj), ' ', jk
                    WRITE(numout, *) '                     ', ji, ' ', jj, ' ', jk, ' ', narea
                    WRITE(numout, *) ' we are now looking for the closest wet cell on the horizontal '
                  END IF
                END IF
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    STOP ' iscpl_cons:   please modify this module !'
    ALLOCATE(inpts(jpnij))
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    inpts(:) = 0
    DO jk = 1, jpk - 1
      !$ACC loop independent collapse(2)
      DO jj = 2, jpj - 1
        DO ji = 2, jpim1
          IF (ptmask_b(ji, jj, jk) == 1._wp .AND. tmask(ji, jj, jk + 1) == 0._wp .AND. tmask_h(ji, jj) == 1._wp .AND. SUM(tmask(ji - 1 : ji + 1, jj, jk)) + SUM(tmask(ji, jj - 1 : jj + 1, jk)) == 0._wp) THEN
            inpts(narea) = inpts(narea) + 1
          END IF
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('iscpl_cons', 'r1', 0, 0)
    CALL mpp_max('iscplhsb', inpts, jpnij)
    npts = SUM(inpts)
    ALLOCATE(ixpts(npts), iypts(npts), izpts(npts), zcorr_vol(npts), zcorr_sal(npts), zcorr_tem(npts), zlon(npts), zlat(npts))
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    ixpts(:) = - 9999
    iypts(:) = - 9999
    izpts(:) = - 9999
    zlon(:) = - 1.0E20_wp
    zlat(:) = - 1.0E20_wp
    zcorr_vol(:) = - 1.0E20_wp
    zcorr_sal(:) = - 1.0E20_wp
    zcorr_tem(:) = - 1.0E20_wp
    jpts = SUM(inpts(1 : narea - 1))
    DO jk = 1, jpk - 1
      !$ACC loop independent collapse(2)
      DO jj = 2, jpj - 1
        DO ji = 2, jpim1
          IF (ptmask_b(ji, jj, jk) == 1._wp .AND. tmask(ji, jj, jk + 1) == 0._wp .AND. tmask_h(ji, jj) == 1._wp .AND. SUM(tmask(ji - 1 : ji + 1, jj, jk)) + SUM(tmask(ji, jj - 1 : jj + 1, jk)) == 0._wp) THEN
            jpts = jpts + 1
            ixpts(jpts) = ji
            iypts(jpts) = jj
            izpts(jpts) = jk
            zlon(jpts) = glamt(ji, jj)
            zlat(jpts) = gphit(ji, jj)
            zcorr_vol(jpts) = pvol_flx(ji, jj, jk)
            zcorr_sal(jpts) = pts_flx(ji, jj, jk, jp_sal)
            zcorr_tem(jpts) = pts_flx(ji, jj, jk, jp_tem)
            pvol_flx(ji, jj, jk) = 0.0_wp
            pts_flx(ji, jj, jk, jp_sal) = 0.0_wp
            pts_flx(ji, jj, jk, jp_tem) = 0.0_wp
          END IF
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('iscpl_cons', 'r2', 0, 0)
    CALL mpp_max('iscplhsb', zlat, npts)
    CALL mpp_max('iscplhsb', zlon, npts)
    CALL mpp_max('iscplhsb', izpts, npts)
    CALL mpp_max('iscplhsb', zcorr_vol, npts)
    CALL mpp_max('iscplhsb', zcorr_sal, npts)
    CALL mpp_max('iscplhsb', zcorr_tem, npts)
    CALL profile_psy_data2 % PostEnd
    DO jpts = 1, npts
      CALL profile_psy_data3 % PreStart('iscpl_cons', 'r3', 0, 0)
      CALL dom_ngb(zlon(jpts), zlat(jpts), ixpts(jpts), iypts(jpts), 'T', izpts(jpts))
      CALL profile_psy_data3 % PostEnd
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = mj0(iypts(jpts)), mj1(iypts(jpts))
        DO ji = mi0(ixpts(jpts)), mi1(ixpts(jpts))
          jk = izpts(jpts)
          IF (tmask_h(ji, jj) == 1._wp) THEN
            pvol_flx(ji, jj, jk) = pvol_flx(ji, jj, jk) + zcorr_vol(jpts)
            pts_flx(ji, jj, jk, jp_sal) = pts_flx(ji, jj, jk, jp_sal) + zcorr_sal(jpts)
            pts_flx(ji, jj, jk, jp_tem) = pts_flx(ji, jj, jk, jp_tem) + zcorr_tem(jpts)
            zcorr_vol(jpts) = 0.0_wp
            zcorr_sal(jpts) = 0.0_wp
            zcorr_tem(jpts) = 0.0_wp
          END IF
        END DO
      END DO
      !$ACC END KERNELS
    END DO
    CALL profile_psy_data4 % PreStart('iscpl_cons', 'r4', 0, 0)
    DEALLOCATE(inpts)
    DEALLOCATE(ixpts, iypts, izpts, zcorr_vol, zcorr_sal, zcorr_tem, zlon, zlat)
    CALL profile_psy_data4 % PostEnd
    !$ACC KERNELS
    pvol_flx(:, :, :) = pvol_flx(:, :, :) * tmask(:, :, :)
    pts_flx(:, :, :, jp_sal) = pts_flx(:, :, :, jp_sal) * tmask(:, :, :)
    pts_flx(:, :, :, jp_tem) = pts_flx(:, :, :, jp_tem) * tmask(:, :, :)
    !$ACC END KERNELS
  END SUBROUTINE iscpl_cons
  SUBROUTINE iscpl_div(phdivn)
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: phdivn
    INTEGER :: ji, jj, jk
    !$ACC KERNELS
    DO jk = 1, jpk
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          phdivn(ji, jj, jk) = phdivn(ji, jj, jk) + hdiv_iscpl(ji, jj, jk) / e3t_n(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE iscpl_div
END MODULE iscplhsb