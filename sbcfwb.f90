MODULE sbcfwb
  USE oce
  USE dom_oce
  USE sbc_oce
  USE sbc_ice, ONLY: snwice_mass, snwice_mass_b, snwice_fmass
  USE phycst
  USE sbcrnf
  USE sbcisf
  USE sbcssr
  USE in_out_manager
  USE lib_mpp
  USE timing
  USE lbclnk
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sbc_fwb
  REAL(KIND = wp) :: a_fwb_b
  REAL(KIND = wp) :: a_fwb
  REAL(KIND = wp) :: fwfold
  REAL(KIND = wp) :: area
  CONTAINS
  SUBROUTINE sbc_fwb(kt, kn_fwb, kn_fsbc)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kn_fsbc
    INTEGER, INTENT(IN) :: kn_fwb
    INTEGER :: inum, ikty, iyear
    REAL(KIND = wp) :: z_fwf, z_fwf_nsrf, zsum_fwf, zsum_erp
    REAL(KIND = wp) :: zsurf_neg, zsurf_pos, zsurf_tospread, zcoef
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: ztmsk_neg, ztmsk_pos, z_wgt
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: ztmsk_tospread, zerp_cor
    REAL(KIND = wp), DIMENSION(1) :: z_fwfprv
    COMPLEX(KIND = wp), DIMENSION(1) :: y_fwfnow
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data7
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data8
    IF (kt == nit000) THEN
      CALL profile_psy_data0 % PreStart('sbc_fwb', 'r0', 0, 0)
      IF (lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) 'sbc_fwb : FreshWater Budget correction'
        WRITE(numout, FMT = *) '~~~~~~~'
        IF (kn_fwb == 1) WRITE(numout, FMT = *) '          instantaneously set to zero'
        IF (kn_fwb == 2) WRITE(numout, FMT = *) '          adjusted from previous year budget'
        IF (kn_fwb == 3) WRITE(numout, FMT = *) '          fwf set to zero and spread out over erp area'
      END IF
      IF (kn_fwb == 3 .AND. nn_sssr /= 2) CALL ctl_stop('sbc_fwb: nn_fwb = 3 requires nn_sssr = 2, we stop ')
      IF (kn_fwb == 3 .AND. ln_isfcav) CALL ctl_stop('sbc_fwb: nn_fwb = 3 with ln_isfcav = .TRUE. not working, we stop ')
      area = glob_sum('sbcfwb', e1e2t(:, :) * tmask(:, :, 1))
      CALL profile_psy_data0 % PostEnd
      !$ACC KERNELS
      snwice_mass_b(:, :) = 0.E0
      snwice_mass(:, :) = 0.E0
      !$ACC END KERNELS
    END IF
    SELECT CASE (kn_fwb)
    CASE (1)
      IF (MOD(kt - 1, kn_fsbc) == 0) THEN
        CALL profile_psy_data1 % PreStart('sbc_fwb', 'r1', 0, 0)
        y_fwfnow(1) = local_sum(e1e2t(:, :) * (emp(:, :) - rnf(:, :) + fwfisf(:, :) - snwice_fmass(:, :)))
        CALL mpp_delay_sum('sbcfwb', 'fwb', y_fwfnow(:), z_fwfprv(:), kt == nitend - nn_fsbc + 1)
        CALL profile_psy_data1 % PostEnd
        !$ACC KERNELS
        z_fwfprv(1) = z_fwfprv(1) / area
        zcoef = z_fwfprv(1) * rcp
        emp(:, :) = emp(:, :) - z_fwfprv(1) * tmask(:, :, 1)
        qns(:, :) = qns(:, :) + zcoef * sst_m(:, :) * tmask(:, :, 1)
        !$ACC END KERNELS
      END IF
    CASE (2)
      CALL profile_psy_data2 % PreStart('sbc_fwb', 'r2', 0, 0)
      IF (kt == nit000) THEN
        CALL ctl_opn(inum, 'EMPave_old.dat', 'OLD', 'FORMATTED', 'SEQUENTIAL', - 1, numout, .FALSE.)
        READ(inum, FMT = "(24X,I8,2ES24.16)") iyear, a_fwb_b, a_fwb
        CLOSE(UNIT = inum)
        fwfold = a_fwb
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) 'sbc_fwb : year = ', iyear, ' freshwater budget correction = ', fwfold
        IF (lwp) WRITE(numout, FMT = *) '          year = ', iyear - 1, ' freshwater budget read       = ', a_fwb
        IF (lwp) WRITE(numout, FMT = *) '          year = ', iyear - 2, ' freshwater budget read       = ', a_fwb_b
      END IF
      ikty = 365 * 86400 / rdt
      IF (MOD(kt, ikty) == 0) THEN
        a_fwb_b = a_fwb
        a_fwb = glob_sum('sbcfwb', e1e2t(:, :) * (sshn(:, :) + snwice_mass(:, :) * r1_rau0))
        a_fwb = a_fwb * 1.E+3 / (area * rday * 365.)
        fwfold = a_fwb
      END IF
      CALL profile_psy_data2 % PostEnd
      !$ACC KERNELS
      IF (MOD(kt - 1, kn_fsbc) == 0) THEN
        zcoef = fwfold * rcp
        emp(:, :) = emp(:, :) + fwfold * tmask(:, :, 1)
        qns(:, :) = qns(:, :) - zcoef * sst_m(:, :) * tmask(:, :, 1)
      END IF
      !$ACC END KERNELS
      CALL profile_psy_data3 % PreStart('sbc_fwb', 'r3', 0, 0)
      IF (kt == nitend .AND. lwm) THEN
        CALL ctl_opn(inum, 'EMPave.dat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, numout, .FALSE., narea)
        WRITE(inum, FMT = "(24X,I8,2ES24.16)") nyear, a_fwb_b, a_fwb
        CLOSE(UNIT = inum)
      END IF
      CALL profile_psy_data3 % PostEnd
    CASE (3)
      ALLOCATE(ztmsk_neg(jpi, jpj), ztmsk_pos(jpi, jpj), ztmsk_tospread(jpi, jpj), z_wgt(jpi, jpj), zerp_cor(jpi, jpj))
      IF (MOD(kt - 1, kn_fsbc) == 0) THEN
        !$ACC KERNELS
        ztmsk_pos(:, :) = tmask_i(:, :)
        !$ACC END KERNELS
        WHERE (erp < 0._wp) ztmsk_pos = 0._wp
        !$ACC KERNELS
        ztmsk_neg(:, :) = tmask_i(:, :) - ztmsk_pos(:, :)
        !$ACC END KERNELS
        CALL profile_psy_data4 % PreStart('sbc_fwb', 'r4', 0, 0)
        z_fwf = glob_sum('sbcfwb', e1e2t(:, :) * (emp(:, :) - rnf(:, :) + fwfisf(:, :) - snwice_fmass(:, :))) / area
        CALL profile_psy_data4 % PostEnd
        IF (z_fwf < 0._wp) THEN
          CALL profile_psy_data5 % PreStart('sbc_fwb', 'r5', 0, 0)
          zsurf_pos = glob_sum('sbcfwb', e1e2t(:, :) * ztmsk_pos(:, :))
          CALL profile_psy_data5 % PostEnd
          !$ACC KERNELS
          zsurf_tospread = zsurf_pos
          ztmsk_tospread(:, :) = ztmsk_pos(:, :)
          !$ACC END KERNELS
        ELSE
          CALL profile_psy_data6 % PreStart('sbc_fwb', 'r6', 0, 0)
          zsurf_neg = glob_sum('sbcfwb', e1e2t(:, :) * ztmsk_neg(:, :))
          CALL profile_psy_data6 % PostEnd
          !$ACC KERNELS
          zsurf_tospread = zsurf_neg
          ztmsk_tospread(:, :) = ztmsk_neg(:, :)
          !$ACC END KERNELS
        END IF
        CALL profile_psy_data7 % PreStart('sbc_fwb', 'r7', 0, 0)
        zsum_fwf = glob_sum('sbcfwb', e1e2t(:, :) * z_fwf)
        z_fwf_nsrf = zsum_fwf / (zsurf_tospread + rsmall)
        zsum_erp = glob_sum('sbcfwb', ztmsk_tospread(:, :) * erp(:, :) * e1e2t(:, :))
        CALL profile_psy_data7 % PostEnd
        !$ACC KERNELS
        z_wgt(:, :) = ztmsk_tospread(:, :) * erp(:, :) / (zsum_erp + rsmall)
        zerp_cor(:, :) = - 1. * z_fwf_nsrf * zsurf_tospread * z_wgt(:, :)
        !$ACC END KERNELS
        CALL lbc_lnk('sbcfwb', zerp_cor, 'T', 1.)
        !$ACC KERNELS
        emp(:, :) = emp(:, :) + zerp_cor(:, :)
        qns(:, :) = qns(:, :) - zerp_cor(:, :) * rcp * sst_m(:, :)
        erp(:, :) = erp(:, :) + zerp_cor(:, :)
        !$ACC END KERNELS
        CALL profile_psy_data8 % PreStart('sbc_fwb', 'r8', 0, 0)
        IF (nprint == 1 .AND. lwp) THEN
          IF (z_fwf < 0._wp) THEN
            WRITE(numout, FMT = *) '   z_fwf < 0'
            WRITE(numout, FMT = *) '   SUM(erp+)     = ', SUM(ztmsk_tospread(:, :) * erp(:, :) * e1e2t(:, :)) * 1.E-9, ' Sv'
          ELSE
            WRITE(numout, FMT = *) '   z_fwf >= 0'
            WRITE(numout, FMT = *) '   SUM(erp-)     = ', SUM(ztmsk_tospread(:, :) * erp(:, :) * e1e2t(:, :)) * 1.E-9, ' Sv'
          END IF
          WRITE(numout, FMT = *) '   SUM(empG)     = ', SUM(z_fwf * e1e2t(:, :)) * 1.E-9, ' Sv'
          WRITE(numout, FMT = *) '   z_fwf         = ', z_fwf, ' Kg/m2/s'
          WRITE(numout, FMT = *) '   z_fwf_nsrf    = ', z_fwf_nsrf, ' Kg/m2/s'
          WRITE(numout, FMT = *) '   MIN(zerp_cor) = ', MINVAL(zerp_cor)
          WRITE(numout, FMT = *) '   MAX(zerp_cor) = ', MAXVAL(zerp_cor)
        END IF
        CALL profile_psy_data8 % PostEnd
      END IF
      DEALLOCATE(ztmsk_neg, ztmsk_pos, ztmsk_tospread, z_wgt, zerp_cor)
    CASE DEFAULT
      CALL ctl_stop('sbc_fwb : wrong nn_fwb value for the FreshWater Budget correction, choose either 1, 2 or 3')
    END SELECT
  END SUBROUTINE sbc_fwb
END MODULE sbcfwb