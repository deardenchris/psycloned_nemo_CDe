MODULE trdtra
  USE oce
  USE dom_oce
  USE sbc_oce
  USE zdf_oce
  USE trd_oce
  USE trdtrc
  USE trdglo
  USE trdpen
  USE trdmxl
  USE ldftra
  USE ldfslp
  USE zdfddm
  USE phycst
  USE in_out_manager
  USE iom
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: trd_tra
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: trdtx, trdty, trdt
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: avt_evd
  CONTAINS
  INTEGER FUNCTION trd_tra_alloc()
    ALLOCATE(trdtx(jpi, jpj, jpk), trdty(jpi, jpj, jpk), trdt(jpi, jpj, jpk), avt_evd(jpi, jpj, jpk), STAT = trd_tra_alloc)
    CALL mpp_sum('trdtra', trd_tra_alloc)
    IF (trd_tra_alloc /= 0) CALL ctl_stop('STOP', 'trd_tra_alloc: failed to allocate arrays')
  END FUNCTION trd_tra_alloc
  SUBROUTINE trd_tra(kt, ctype, ktra, ktrd, ptrd, pun, ptra)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    CHARACTER(LEN = 3), INTENT(IN) :: ctype
    INTEGER, INTENT(IN) :: ktra
    INTEGER, INTENT(IN) :: ktrd
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: ptrd
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN), OPTIONAL :: pun
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN), OPTIONAL :: ptra
    INTEGER :: jk
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: ztrds
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zwt, zws, ztrdt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data7
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data8
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data9
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data10
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data11
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data12
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data13
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data14
    CALL profile_psy_data0 % PreStart('trd_tra', 'r0', 0, 0)
    IF (.NOT. ALLOCATED(trdtx)) THEN
      IF (trd_tra_alloc() /= 0) CALL ctl_stop('STOP', 'trd_tra : unable to allocate arrays')
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (ctype == 'TRA' .AND. ktra == jp_tem) THEN
      SELECT CASE (ktrd)
      CASE (jptra_xad)
        CALL profile_psy_data1 % PreStart('trd_tra', 'r1', 0, 0)
        CALL trd_tra_adv(ptrd, pun, ptra, 'X', trdtx)
        CALL profile_psy_data1 % PostEnd
      CASE (jptra_yad)
        CALL profile_psy_data2 % PreStart('trd_tra', 'r2', 0, 0)
        CALL trd_tra_adv(ptrd, pun, ptra, 'Y', trdty)
        CALL profile_psy_data2 % PostEnd
      CASE (jptra_zad)
        CALL profile_psy_data3 % PreStart('trd_tra', 'r3', 0, 0)
        CALL trd_tra_adv(ptrd, pun, ptra, 'Z', trdt)
        CALL profile_psy_data3 % PostEnd
      CASE (jptra_bbc, jptra_qsr)
        !$ACC KERNELS
        trdt(:, :, :) = ptrd(:, :, :) * tmask(:, :, :)
        ztrds(:, :, :) = 0._wp
        !$ACC END KERNELS
        CALL profile_psy_data4 % PreStart('trd_tra', 'r4', 0, 0)
        CALL trd_tra_mng(trdt, ztrds, ktrd, kt)
        CALL profile_psy_data4 % PostEnd
      CASE (jptra_evd)
        !$ACC KERNELS
        avt_evd(:, :, :) = ptrd(:, :, :) * tmask(:, :, :)
        !$ACC END KERNELS
      CASE DEFAULT
        !$ACC KERNELS
        trdt(:, :, :) = ptrd(:, :, :) * tmask(:, :, :)
        !$ACC END KERNELS
      END SELECT
    END IF
    IF (ctype == 'TRA' .AND. ktra == jp_sal) THEN
      SELECT CASE (ktrd)
      CASE (jptra_xad)
        CALL profile_psy_data5 % PreStart('trd_tra', 'r5', 0, 0)
        CALL trd_tra_adv(ptrd, pun, ptra, 'X', ztrds)
        CALL trd_tra_mng(trdtx, ztrds, ktrd, kt)
        CALL profile_psy_data5 % PostEnd
      CASE (jptra_yad)
        CALL profile_psy_data6 % PreStart('trd_tra', 'r6', 0, 0)
        CALL trd_tra_adv(ptrd, pun, ptra, 'Y', ztrds)
        CALL trd_tra_mng(trdty, ztrds, ktrd, kt)
        CALL profile_psy_data6 % PostEnd
      CASE (jptra_zad)
        CALL profile_psy_data7 % PreStart('trd_tra', 'r7', 0, 0)
        CALL trd_tra_adv(ptrd, pun, ptra, 'Z', ztrds)
        CALL trd_tra_mng(trdt, ztrds, ktrd, kt)
        CALL profile_psy_data7 % PostEnd
      CASE (jptra_zdfp)
        ALLOCATE(zwt(jpi, jpj, jpk), zws(jpi, jpj, jpk), ztrdt(jpi, jpj, jpk))
        !$ACC KERNELS
        zwt(:, :, 1) = 0._wp
        zws(:, :, 1) = 0._wp
        zwt(:, :, jpk) = 0._wp
        zws(:, :, jpk) = 0._wp
        DO jk = 2, jpk
          zwt(:, :, jk) = avt(:, :, jk) * (tsa(:, :, jk - 1, jp_tem) - tsa(:, :, jk, jp_tem)) / e3w_n(:, :, jk) * tmask(:, :, jk)
          zws(:, :, jk) = avs(:, :, jk) * (tsa(:, :, jk - 1, jp_sal) - tsa(:, :, jk, jp_sal)) / e3w_n(:, :, jk) * tmask(:, :, jk)
        END DO
        ztrdt(:, :, jpk) = 0._wp
        ztrds(:, :, jpk) = 0._wp
        DO jk = 1, jpkm1
          ztrdt(:, :, jk) = (zwt(:, :, jk) - zwt(:, :, jk + 1)) / e3t_n(:, :, jk)
          ztrds(:, :, jk) = (zws(:, :, jk) - zws(:, :, jk + 1)) / e3t_n(:, :, jk)
        END DO
        !$ACC END KERNELS
        CALL profile_psy_data8 % PreStart('trd_tra', 'r8', 0, 0)
        CALL trd_tra_mng(ztrdt, ztrds, jptra_zdfp, kt)
        CALL profile_psy_data8 % PostEnd
        !$ACC KERNELS
        zwt(:, :, :) = 0._wp
        zws(:, :, :) = 0._wp
        DO jk = 2, jpk
          zwt(:, :, jk) = avt_evd(:, :, jk) * (tsa(:, :, jk - 1, jp_tem) - tsa(:, :, jk, jp_tem)) / e3w_n(:, :, jk) * tmask(:, :, jk)
          zws(:, :, jk) = avt_evd(:, :, jk) * (tsa(:, :, jk - 1, jp_sal) - tsa(:, :, jk, jp_sal)) / e3w_n(:, :, jk) * tmask(:, :, jk)
        END DO
        ztrdt(:, :, jpk) = 0._wp
        ztrds(:, :, jpk) = 0._wp
        DO jk = 1, jpkm1
          ztrdt(:, :, jk) = (zwt(:, :, jk) - zwt(:, :, jk + 1)) / e3t_n(:, :, jk)
          ztrds(:, :, jk) = (zws(:, :, jk) - zws(:, :, jk + 1)) / e3t_n(:, :, jk)
        END DO
        !$ACC END KERNELS
        CALL profile_psy_data9 % PreStart('trd_tra', 'r9', 0, 0)
        CALL trd_tra_mng(ztrdt, ztrds, jptra_evd, kt)
        DEALLOCATE(zwt, zws, ztrdt)
        CALL profile_psy_data9 % PostEnd
      CASE DEFAULT
        !$ACC KERNELS
        ztrds(:, :, :) = ptrd(:, :, :) * tmask(:, :, :)
        !$ACC END KERNELS
        CALL profile_psy_data10 % PreStart('trd_tra', 'r10', 0, 0)
        CALL trd_tra_mng(trdt, ztrds, ktrd, kt)
        CALL profile_psy_data10 % PostEnd
      END SELECT
    END IF
    IF (ctype == 'TRC') THEN
      SELECT CASE (ktrd)
      CASE (jptra_xad)
        CALL profile_psy_data11 % PreStart('trd_tra', 'r11', 0, 0)
        CALL trd_tra_adv(ptrd, pun, ptra, 'X', ztrds)
        CALL profile_psy_data11 % PostEnd
      CASE (jptra_yad)
        CALL profile_psy_data12 % PreStart('trd_tra', 'r12', 0, 0)
        CALL trd_tra_adv(ptrd, pun, ptra, 'Y', ztrds)
        CALL profile_psy_data12 % PostEnd
      CASE (jptra_zad)
        CALL profile_psy_data13 % PreStart('trd_tra', 'r13', 0, 0)
        CALL trd_tra_adv(ptrd, pun, ptra, 'Z', ztrds)
        CALL profile_psy_data13 % PostEnd
      CASE DEFAULT
        !$ACC KERNELS
        ztrds(:, :, :) = ptrd(:, :, :) * tmask(:, :, :)
        !$ACC END KERNELS
      END SELECT
      CALL profile_psy_data14 % PreStart('trd_tra', 'r14', 0, 0)
      CALL trd_trc(ztrds, ktra, ktrd, kt)
      CALL profile_psy_data14 % PostEnd
    END IF
  END SUBROUTINE trd_tra
  SUBROUTINE trd_tra_adv(pf, pun, ptn, cdir, ptrd)
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pf
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pun
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: ptn
    CHARACTER(LEN = 1), INTENT(IN) :: cdir
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(OUT) :: ptrd
    INTEGER :: ji, jj, jk
    INTEGER :: ii, ij, ik
    !$ACC KERNELS
    SELECT CASE (cdir)
    CASE ('X')
      ii = 1
      ij = 0
      ik = 0
    CASE ('Y')
      ii = 0
      ij = 1
      ik = 0
    CASE ('Z')
      ii = 0
      ij = 0
      ik = - 1
    END SELECT
    ptrd(jpi, :, :) = 0._wp
    ptrd(1, :, :) = 0._wp
    ptrd(:, jpj, :) = 0._wp
    ptrd(:, 1, :) = 0._wp
    ptrd(:, :, jpk) = 0._wp
    DO jk = 1, jpkm1
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ptrd(ji, jj, jk) = - (pf(ji, jj, jk) - pf(ji - ii, jj - ij, jk - ik) - (pun(ji, jj, jk) - pun(ji - ii, jj - ij, jk - ik)) * ptn(ji, jj, jk)) * r1_e1e2t(ji, jj) / e3t_n(ji, jj, jk) * tmask(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE trd_tra_adv
  SUBROUTINE trd_tra_mng(ptrdx, ptrdy, ktrd, kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: ptrdx
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: ptrdy
    INTEGER, INTENT(IN) :: ktrd
    INTEGER, INTENT(IN) :: kt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data7
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data8
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data9
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data10
    CALL profile_psy_data0 % PreStart('trd_tra_mng', 'r0', 0, 0)
    IF (neuler == 0 .AND. kt == nit000) THEN
      r2dt = rdt
    ELSE IF (kt <= nit000 + 1) THEN
      r2dt = 2. * rdt
    END IF
    IF (ln_tra_trd) CALL trd_tra_iom(ptrdx, ptrdy, ktrd, kt)
    IF (ln_glo_trd) CALL trd_glo(ptrdx, ptrdy, ktrd, 'TRA', kt)
    IF (ln_PE_trd) CALL trd_pen(ptrdx, ptrdy, ktrd, kt, r2dt)
    CALL profile_psy_data0 % PostEnd
    IF (ln_tra_mxl) THEN
      SELECT CASE (ktrd)
      CASE (jptra_xad)
        CALL profile_psy_data1 % PreStart('trd_tra_mng', 'r1', 0, 0)
        CALL trd_mxl_zint(ptrdx, ptrdy, jpmxl_xad, '3D')
        CALL profile_psy_data1 % PostEnd
      CASE (jptra_yad)
        CALL profile_psy_data2 % PreStart('trd_tra_mng', 'r2', 0, 0)
        CALL trd_mxl_zint(ptrdx, ptrdy, jpmxl_yad, '3D')
        CALL profile_psy_data2 % PostEnd
      CASE (jptra_zad)
        CALL profile_psy_data3 % PreStart('trd_tra_mng', 'r3', 0, 0)
        CALL trd_mxl_zint(ptrdx, ptrdy, jpmxl_zad, '3D')
        CALL profile_psy_data3 % PostEnd
      CASE (jptra_ldf)
        CALL profile_psy_data4 % PreStart('trd_tra_mng', 'r4', 0, 0)
        CALL trd_mxl_zint(ptrdx, ptrdy, jpmxl_ldf, '3D')
        CALL profile_psy_data4 % PostEnd
      CASE (jptra_bbl)
        CALL profile_psy_data5 % PreStart('trd_tra_mng', 'r5', 0, 0)
        CALL trd_mxl_zint(ptrdx, ptrdy, jpmxl_bbl, '3D')
        CALL profile_psy_data5 % PostEnd
      CASE (jptra_zdf)
        CALL profile_psy_data6 % PreStart('trd_tra_mng', 'r6', 0, 0)
        IF (ln_traldf_iso) THEN
          CALL trd_mxl_zint(ptrdx, ptrdy, jpmxl_ldf, '3D')
        ELSE
          CALL trd_mxl_zint(ptrdx, ptrdy, jpmxl_zdf, '3D')
        END IF
        CALL profile_psy_data6 % PostEnd
      CASE (jptra_dmp)
        CALL profile_psy_data7 % PreStart('trd_tra_mng', 'r7', 0, 0)
        CALL trd_mxl_zint(ptrdx, ptrdy, jpmxl_dmp, '3D')
        CALL profile_psy_data7 % PostEnd
      CASE (jptra_qsr)
        CALL profile_psy_data8 % PreStart('trd_tra_mng', 'r8', 0, 0)
        CALL trd_mxl_zint(ptrdx, ptrdy, jpmxl_for, '3D')
        CALL profile_psy_data8 % PostEnd
      CASE (jptra_nsr)
        !$ACC KERNELS
        ptrdx(:, :, 2 : jpk) = 0._wp
        ptrdy(:, :, 2 : jpk) = 0._wp
        !$ACC END KERNELS
        CALL profile_psy_data9 % PreStart('trd_tra_mng', 'r9', 0, 0)
        CALL trd_mxl_zint(ptrdx, ptrdy, jpmxl_for, '2D')
        CALL profile_psy_data9 % PostEnd
      CASE (jptra_bbc)
        CALL profile_psy_data10 % PreStart('trd_tra_mng', 'r10', 0, 0)
        CALL trd_mxl_zint(ptrdx, ptrdy, jpmxl_bbc, '3D')
        CALL profile_psy_data10 % PostEnd
      CASE (jptra_npc)
        CALL trd_mxl_zint(ptrdx, ptrdy, jpmxl_npc, '3D')
      CASE (jptra_atf)
        CALL trd_mxl_zint(ptrdx, ptrdy, jpmxl_atf, '3D')
        CALL trd_mxl(kt, r2dt)
      END SELECT
    END IF
  END SUBROUTINE trd_tra_mng
  SUBROUTINE trd_tra_iom(ptrdx, ptrdy, ktrd, kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: ptrdx
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: ptrdy
    INTEGER, INTENT(IN) :: ktrd
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: ikbu, ikbv
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: z2dx, z2dy
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    CALL profile_psy_data0 % PreStart('trd_tra_iom', 'r0', 0, 0)
    SELECT CASE (ktrd)
    CASE (jptra_tot)
      CALL iom_put("ttrd_tot", ptrdx)
      CALL iom_put("strd_tot", ptrdy)
    END SELECT
    CALL profile_psy_data0 % PostEnd
    IF (MOD(kt, 2) == 0) THEN
      SELECT CASE (ktrd)
      CASE (jptra_xad)
        CALL profile_psy_data1 % PreStart('trd_tra_iom', 'r1', 0, 0)
        CALL iom_put("ttrd_xad", ptrdx)
        CALL iom_put("strd_xad", ptrdy)
        CALL profile_psy_data1 % PostEnd
      CASE (jptra_yad)
        CALL profile_psy_data2 % PreStart('trd_tra_iom', 'r2', 0, 0)
        CALL iom_put("ttrd_yad", ptrdx)
        CALL iom_put("strd_yad", ptrdy)
        CALL profile_psy_data2 % PostEnd
      CASE (jptra_zad)
        CALL profile_psy_data3 % PreStart('trd_tra_iom', 'r3', 0, 0)
        CALL iom_put("ttrd_zad", ptrdx)
        CALL iom_put("strd_zad", ptrdy)
        CALL profile_psy_data3 % PostEnd
        IF (ln_linssh) THEN
          ALLOCATE(z2dx(jpi, jpj), z2dy(jpi, jpj))
          !$ACC KERNELS
          z2dx(:, :) = wn(:, :, 1) * tsn(:, :, 1, jp_tem) / e3t_n(:, :, 1)
          z2dy(:, :) = wn(:, :, 1) * tsn(:, :, 1, jp_sal) / e3t_n(:, :, 1)
          !$ACC END KERNELS
          CALL profile_psy_data4 % PreStart('trd_tra_iom', 'r4', 0, 0)
          CALL iom_put("ttrd_sad", z2dx)
          CALL iom_put("strd_sad", z2dy)
          DEALLOCATE(z2dx, z2dy)
          CALL profile_psy_data4 % PostEnd
        END IF
      CASE (jptra_totad)
        CALL profile_psy_data5 % PreStart('trd_tra_iom', 'r5', 0, 0)
        CALL iom_put("ttrd_totad", ptrdx)
        CALL iom_put("strd_totad", ptrdy)
        CALL profile_psy_data5 % PostEnd
      CASE (jptra_ldf)
        CALL iom_put("ttrd_ldf", ptrdx)
        CALL iom_put("strd_ldf", ptrdy)
      CASE (jptra_zdf)
        CALL iom_put("ttrd_zdf", ptrdx)
        CALL iom_put("strd_zdf", ptrdy)
      CASE (jptra_zdfp)
        CALL iom_put("ttrd_zdfp", ptrdx)
        CALL iom_put("strd_zdfp", ptrdy)
      CASE (jptra_evd)
        CALL iom_put("ttrd_evd", ptrdx)
        CALL iom_put("strd_evd", ptrdy)
      CASE (jptra_dmp)
        CALL iom_put("ttrd_dmp", ptrdx)
        CALL iom_put("strd_dmp", ptrdy)
      CASE (jptra_bbl)
        CALL iom_put("ttrd_bbl", ptrdx)
        CALL iom_put("strd_bbl", ptrdy)
      CASE (jptra_npc)
        CALL iom_put("ttrd_npc", ptrdx)
        CALL iom_put("strd_npc", ptrdy)
      CASE (jptra_bbc)
        CALL iom_put("ttrd_bbc", ptrdx)
      CASE (jptra_nsr)
        CALL iom_put("ttrd_qns", ptrdx(:, :, 1))
        CALL iom_put("strd_cdt", ptrdy(:, :, 1))
      CASE (jptra_qsr)
        CALL iom_put("ttrd_qsr", ptrdx)
      END SELECT
    ELSE IF (MOD(kt, 2) == 1) THEN
      CALL profile_psy_data6 % PreStart('trd_tra_iom', 'r6', 0, 0)
      SELECT CASE (ktrd)
      CASE (jptra_atf)
        CALL iom_put("ttrd_atf", ptrdx)
        CALL iom_put("strd_atf", ptrdy)
      END SELECT
      CALL profile_psy_data6 % PostEnd
    END IF
  END SUBROUTINE trd_tra_iom
END MODULE trdtra