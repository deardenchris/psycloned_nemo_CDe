MODULE icedyn_adv_umx
  USE phycst
  USE dom_oce
  USE sbc_oce, ONLY: nn_fsbc
  USE ice
  USE icevar
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_dyn_adv_umx
  INTEGER, PARAMETER :: np_advS = 1
  INTEGER, PARAMETER :: np_limiter = 1
  LOGICAL :: ll_upsxy = .TRUE.
  LOGICAL :: ll_hoxy = .TRUE.
  LOGICAL :: ll_neg = .TRUE.
  LOGICAL :: ll_prelim = .FALSE.
  REAL(KIND = wp) :: z1_6 = 1._wp / 6._wp
  REAL(KIND = wp) :: z1_120 = 1._wp / 120._wp
  INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: imsk_small, jmsk_small
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zfu_ho, zfv_ho, zpt
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zfu_ups, zfv_ups, zt_ups
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ztu1, ztu2, ztu3, ztu4
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ztv1, ztv2, ztv3, ztv4
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zswitch
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zslpy
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zslpx
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zbup, zbdo
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zbetup, zbetdo, zti_ups, ztj_ups
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zt_u, zt_v
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zudy, zvdx, zcu_box, zcv_box
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zati1, zati2
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zu_cat, zv_cat
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zua_ho, zva_ho, zua_ups, zva_ups
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: z1_ai, z1_aip, zhvar
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zhi_max, zhs_max, zhip_max
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zuv_ho, zvv_ho, zuv_ups, zvv_ups, z1_vi, z1_vs
  CONTAINS
  SUBROUTINE ice_dyn_adv_umx(kn_umx, kt, pu_ice, pv_ice, ph_i, ph_s, ph_ip, pato_i, pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, &
&pe_s, pe_i)
    INTEGER, INTENT(IN) :: kn_umx
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pu_ice
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pv_ice
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: ph_i
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: ph_s
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: ph_ip
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: pato_i
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pv_i
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pv_s
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: psv_i
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: poa_i
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pa_i
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pa_ip
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pv_ip
    REAL(KIND = wp), DIMENSION(:, :, :, :), INTENT(INOUT) :: pe_s
    REAL(KIND = wp), DIMENSION(:, :, :, :), INTENT(INOUT) :: pe_i
    INTEGER :: ji, jj, jk, jl, jt
    INTEGER :: icycle
    REAL(KIND = wp) :: zamsk
    REAL(KIND = wp) :: zdt, zvi_cen
    REAL(KIND = wp), DIMENSION(1) :: zcflprv, zcflnow
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *) '-- ice_dyn_adv_umx: Ultimate-Macho advection scheme'
      ALLOCATE(zudy(jpi, jpj), zvdx(jpi, jpj), zcu_box(jpi, jpj), zcv_box(jpi, jpj))
      ALLOCATE(zati1(jpi, jpj), zati2(jpi, jpj))
      ALLOCATE(zu_cat(jpi, jpj, jpl), zv_cat(jpi, jpj, jpl))
      ALLOCATE(zua_ho(jpi, jpj, jpl), zva_ho(jpi, jpj, jpl), zua_ups(jpi, jpj, jpl), zva_ups(jpi, jpj, jpl))
      ALLOCATE(z1_ai(jpi, jpj, jpl), z1_aip(jpi, jpj, jpl), zhvar(jpi, jpj, jpl))
      ALLOCATE(zhi_max(jpi, jpj, jpl), zhs_max(jpi, jpj, jpl), zhip_max(jpi, jpj, jpl))
      ALLOCATE(zfu_ho(jpi, jpj, jpl), zfv_ho(jpi, jpj, jpl), zpt(jpi, jpj, jpl))
      ALLOCATE(zfu_ups(jpi, jpj, jpl), zfv_ups(jpi, jpj, jpl), zt_ups(jpi, jpj, jpl))
      ALLOCATE(ztu1(jpi, jpj, jpl), ztu2(jpi, jpj, jpl), ztu3(jpi, jpj, jpl), ztu4(jpi, jpj, jpl))
      ALLOCATE(ztv1(jpi, jpj, jpl), ztv2(jpi, jpj, jpl), ztv3(jpi, jpj, jpl), ztv4(jpi, jpj, jpl))
      ALLOCATE(zswitch(jpi, jpj))
      ALLOCATE(zslpy(jpi, jpj, jpl))
      ALLOCATE(zslpx(jpi, jpj, jpl))
      ALLOCATE(zbup(jpi, jpj), zbdo(jpi, jpj))
      ALLOCATE(zbetup(jpi, jpj, jpl), zbetdo(jpi, jpj, jpl), zti_ups(jpi, jpj, jpl), ztj_ups(jpi, jpj, jpl))
      ALLOCATE(zt_u(jpi, jpj, jpl), zt_v(jpi, jpj, jpl))
    END IF
    DO jl = 1, jpl
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zhip_max(ji, jj, jl) = MAX(epsi20, ph_ip(ji, jj, jl), ph_ip(ji + 1, jj, jl), ph_ip(ji, jj + 1, jl), ph_ip(ji - 1, jj, &
&jl), ph_ip(ji, jj - 1, jl), ph_ip(ji + 1, jj + 1, jl), ph_ip(ji - 1, jj - 1, jl), ph_ip(ji + 1, jj - 1, jl), ph_ip(ji - 1, jj + &
&1, jl))
          zhi_max(ji, jj, jl) = MAX(epsi20, ph_i(ji, jj, jl), ph_i(ji + 1, jj, jl), ph_i(ji, jj + 1, jl), ph_i(ji - 1, jj, jl), &
&ph_i(ji, jj - 1, jl), ph_i(ji + 1, jj + 1, jl), ph_i(ji - 1, jj - 1, jl), ph_i(ji + 1, jj - 1, jl), ph_i(ji - 1, jj + 1, jl))
          zhs_max(ji, jj, jl) = MAX(epsi20, ph_s(ji, jj, jl), ph_s(ji + 1, jj, jl), ph_s(ji, jj + 1, jl), ph_s(ji - 1, jj, jl), &
&ph_s(ji, jj - 1, jl), ph_s(ji + 1, jj + 1, jl), ph_s(ji - 1, jj - 1, jl), ph_s(ji + 1, jj - 1, jl), ph_s(ji - 1, jj + 1, jl))
        END DO
      END DO
    END DO
    CALL lbc_lnk_multi('icedyn_adv_umx', zhi_max, 'T', 1., zhs_max, 'T', 1., zhip_max, 'T', 1.)
    zcflnow(1) = MAXVAL(ABS(pu_ice(:, :)) * rdt_ice * r1_e1u(:, :))
    zcflnow(1) = MAX(zcflnow(1), MAXVAL(ABS(pv_ice(:, :)) * rdt_ice * r1_e2v(:, :)))
    CALL mpp_delay_max('icedyn_adv_umx', 'cflice', zcflnow(:), zcflprv(:), kt == nitend - nn_fsbc + 1)
    IF (zcflprv(1) > .5) THEN
      icycle = 2
    ELSE
      icycle = 1
    END IF
    zdt = rdt_ice / REAL(icycle)
    zudy(:, :) = pu_ice(:, :) * e2u(:, :)
    zvdx(:, :) = pv_ice(:, :) * e1v(:, :)
    DO jl = 1, jpl
      zu_cat(:, :, jl) = zudy(:, :)
      zv_cat(:, :, jl) = zvdx(:, :)
    END DO
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        IF (pu_ice(ji, jj) * pu_ice(ji - 1, jj) <= 0._wp) THEN
          zcu_box(ji, jj) = 0._wp
        ELSE IF (pu_ice(ji, jj) > 0._wp) THEN
          zcu_box(ji, jj) = pu_ice(ji - 1, jj)
        ELSE
          zcu_box(ji, jj) = pu_ice(ji, jj)
        END IF
        IF (pv_ice(ji, jj) * pv_ice(ji, jj - 1) <= 0._wp) THEN
          zcv_box(ji, jj) = 0._wp
        ELSE IF (pv_ice(ji, jj) > 0._wp) THEN
          zcv_box(ji, jj) = pv_ice(ji, jj - 1)
        ELSE
          zcv_box(ji, jj) = pv_ice(ji, jj)
        END IF
      END DO
    END DO
    DO jt = 1, icycle
      zati1(:, :) = SUM(pa_i(:, :, :), dim = 3)
      WHERE (pa_i(:, :, :) >= epsi20)
        z1_ai(:, :, :) = 1._wp / pa_i(:, :, :)
      ELSEWHERE
        z1_ai(:, :, :) = 0.
      END WHERE
      WHERE (pa_ip(:, :, :) >= epsi20)
        z1_aip(:, :, :) = 1._wp / pa_ip(:, :, :)
      ELSEWHERE
        z1_aip(:, :, :) = 0.
      END WHERE
      IF (ll_neg) THEN
        IF (.NOT. ALLOCATED(imsk_small)) ALLOCATE(imsk_small(jpi, jpj, jpl))
        IF (.NOT. ALLOCATED(jmsk_small)) ALLOCATE(jmsk_small(jpi, jpj, jpl))
        DO jl = 1, jpl
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              zvi_cen = 0.5_wp * (pv_i(ji + 1, jj, jl) + pv_i(ji, jj, jl))
              IF (zvi_cen < epsi06) THEN
                imsk_small(ji, jj, jl) = 0
              ELSE
                imsk_small(ji, jj, jl) = 1
              END IF
              zvi_cen = 0.5_wp * (pv_i(ji, jj + 1, jl) + pv_i(ji, jj, jl))
              IF (zvi_cen < epsi06) THEN
                jmsk_small(ji, jj, jl) = 0
              ELSE
                jmsk_small(ji, jj, jl) = 1
              END IF
            END DO
          END DO
        END DO
      END IF
      zamsk = 1._wp
      CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zu_cat, zv_cat, zcu_box, zcv_box, pa_i, pa_i, zua_ups, zva_ups, zua_ho, &
&zva_ho)
      IF (np_advS == 1) THEN
        zamsk = 0._wp
        zhvar(:, :, :) = pv_i(:, :, :) * z1_ai(:, :, :)
        CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, pv_i, zua_ups, zva_ups)
        zhvar(:, :, :) = pv_s(:, :, :) * z1_ai(:, :, :)
        CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, pv_s, zua_ups, zva_ups)
        zamsk = 1._wp
        CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zu_cat, zv_cat, zcu_box, zcv_box, psv_i, psv_i)
        DO jk = 1, nlay_i
          CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zu_cat, zv_cat, zcu_box, zcv_box, pe_i(:, :, jk, :), pe_i(:, :, jk, &
&:))
        END DO
        DO jk = 1, nlay_s
          CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zu_cat, zv_cat, zcu_box, zcv_box, pe_s(:, :, jk, :), pe_s(:, :, jk, &
&:))
        END DO
      ELSE IF (np_advS == 2) THEN
        zamsk = 0._wp
        zhvar(:, :, :) = pv_i(:, :, :) * z1_ai(:, :, :)
        CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, pv_i, zua_ups, zva_ups)
        zhvar(:, :, :) = pv_s(:, :, :) * z1_ai(:, :, :)
        CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, pv_s, zua_ups, zva_ups)
        zhvar(:, :, :) = psv_i(:, :, :) * z1_ai(:, :, :)
        CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, psv_i, zua_ups, zva_ups)
        DO jk = 1, nlay_i
          zhvar(:, :, :) = pe_i(:, :, jk, :) * z1_ai(:, :, :)
          CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, pe_i(:, :, jk, :), &
&zua_ups, zva_ups)
        END DO
        DO jk = 1, nlay_s
          zhvar(:, :, :) = pe_s(:, :, jk, :) * z1_ai(:, :, :)
          CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, pe_s(:, :, jk, :), &
&zua_ups, zva_ups)
        END DO
      ELSE IF (np_advS == 3) THEN
        zamsk = 0._wp
        ALLOCATE(zuv_ho(jpi, jpj, jpl), zvv_ho(jpi, jpj, jpl), zuv_ups(jpi, jpj, jpl), zvv_ups(jpi, jpj, jpl), z1_vi(jpi, jpj, &
&jpl), z1_vs(jpi, jpj, jpl))
        WHERE (pv_i(:, :, :) >= epsi20)
          z1_vi(:, :, :) = 1._wp / pv_i(:, :, :)
        ELSEWHERE
          z1_vi(:, :, :) = 0.
        END WHERE
        WHERE (pv_s(:, :, :) >= epsi20)
          z1_vs(:, :, :) = 1._wp / pv_s(:, :, :)
        ELSEWHERE
          z1_vs(:, :, :) = 0.
        END WHERE
        zuv_ups = zua_ups
        zvv_ups = zva_ups
        zhvar(:, :, :) = pv_i(:, :, :) * z1_ai(:, :, :)
        CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, pv_i, zuv_ups, zvv_ups, &
&zuv_ho, zvv_ho)
        zhvar(:, :, :) = psv_i(:, :, :) * z1_vi(:, :, :)
        CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zuv_ho, zvv_ho, zcu_box, zcv_box, zhvar, psv_i, zuv_ups, zvv_ups)
        DO jk = 1, nlay_i
          zhvar(:, :, :) = pe_i(:, :, jk, :) * z1_vi(:, :, :)
          CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zuv_ho, zvv_ho, zcu_box, zcv_box, zhvar, pe_i(:, :, jk, :), &
&zuv_ups, zvv_ups)
        END DO
        zuv_ups = zua_ups
        zvv_ups = zva_ups
        zhvar(:, :, :) = pv_s(:, :, :) * z1_ai(:, :, :)
        CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, pv_s, zuv_ups, zvv_ups, &
&zuv_ho, zvv_ho)
        DO jk = 1, nlay_s
          zhvar(:, :, :) = pe_s(:, :, jk, :) * z1_vs(:, :, :)
          CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zuv_ho, zvv_ho, zcu_box, zcv_box, zhvar, pe_s(:, :, jk, :), &
&zuv_ups, zvv_ups)
        END DO
        DEALLOCATE(zuv_ho, zvv_ho, zuv_ups, zvv_ups, z1_vi, z1_vs)
      END IF
      IF (iom_use('iceage') .OR. iom_use('iceage_cat')) THEN
        zamsk = 1._wp
        CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zu_cat, zv_cat, zcu_box, zcv_box, poa_i, poa_i)
      END IF
      IF (ln_pnd_H12) THEN
        zamsk = 1._wp
        CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zu_cat, zv_cat, zcu_box, zcv_box, pa_ip, pa_ip, zua_ups, zva_ups, &
&zua_ho, zva_ho)
        zamsk = 0._wp
        zhvar(:, :, :) = pv_ip(:, :, :) * z1_aip(:, :, :)
        CALL adv_umx(zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, pv_ip, zua_ups, zva_ups)
      END IF
      zati2(:, :) = SUM(pa_i(:, :, :), dim = 3)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          pato_i(ji, jj) = pato_i(ji, jj) - (zati2(ji, jj) - zati1(ji, jj)) - (zudy(ji, jj) - zudy(ji - 1, jj) + zvdx(ji, jj) - &
&zvdx(ji, jj - 1)) * r1_e1e2t(ji, jj) * zdt
        END DO
      END DO
      CALL lbc_lnk('icedyn_adv_umx', pato_i, 'T', 1.)
      CALL ice_var_zapneg(zdt, pato_i, pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, pe_s, pe_i)
      CALL Hbig(zdt, zhi_max, zhs_max, zhip_max, pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, pe_s, pe_i)
    END DO
  END SUBROUTINE ice_dyn_adv_umx
  SUBROUTINE adv_umx(pamsk, kn_umx, jt, kt, pdt, pu, pv, puc, pvc, pubox, pvbox, pt, ptc, pua_ups, pva_ups, pua_ho, pva_ho)
    REAL(KIND = wp), INTENT(IN) :: pamsk
    INTEGER, INTENT(IN) :: kn_umx
    INTEGER, INTENT(IN) :: jt
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), INTENT(IN) :: pdt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pu, pv
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: puc, pvc
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pubox, pvbox
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: ptc
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl), INTENT(INOUT), OPTIONAL :: pua_ups, pva_ups
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl), INTENT(OUT), OPTIONAL :: pua_ho, pva_ho
    INTEGER :: ji, jj, jl
    REAL(KIND = wp) :: ztra
    CALL upstream(pamsk, jt, kt, pdt, pt, pu, pv, zt_ups, zfu_ups, zfv_ups)
    SELECT CASE (kn_umx)
    CASE (20)
      CALL cen2(pamsk, jt, kt, pdt, pt, pu, pv, zt_ups, zfu_ups, zfv_ups, zfu_ho, zfv_ho)
    CASE (1 : 5)
      CALL macho(pamsk, kn_umx, jt, kt, pdt, pt, pu, pv, pubox, pvbox, zt_ups, zfu_ups, zfv_ups, zfu_ho, zfv_ho)
    END SELECT
    IF (pamsk == 0._wp) THEN
      DO jl = 1, jpl
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            IF (ABS(pu(ji, jj)) > epsi10) THEN
              zfu_ho(ji, jj, jl) = zfu_ho(ji, jj, jl) * puc(ji, jj, jl) / pu(ji, jj)
              zfu_ups(ji, jj, jl) = zfu_ups(ji, jj, jl) * pua_ups(ji, jj, jl) / pu(ji, jj)
            ELSE
              zfu_ho(ji, jj, jl) = 0._wp
              zfu_ups(ji, jj, jl) = 0._wp
            END IF
            IF (ABS(pv(ji, jj)) > epsi10) THEN
              zfv_ho(ji, jj, jl) = zfv_ho(ji, jj, jl) * pvc(ji, jj, jl) / pv(ji, jj)
              zfv_ups(ji, jj, jl) = zfv_ups(ji, jj, jl) * pva_ups(ji, jj, jl) / pv(ji, jj)
            ELSE
              zfv_ho(ji, jj, jl) = 0._wp
              zfv_ups(ji, jj, jl) = 0._wp
            END IF
          END DO
        END DO
      END DO
      DO jl = 1, jpl
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            ztra = - (zfu_ups(ji, jj, jl) - zfu_ups(ji - 1, jj, jl) + zfv_ups(ji, jj, jl) - zfv_ups(ji, jj - 1, jl))
            zt_ups(ji, jj, jl) = (ptc(ji, jj, jl) + ztra * r1_e1e2t(ji, jj) * pdt) * tmask(ji, jj, 1)
          END DO
        END DO
      END DO
      CALL lbc_lnk('icedyn_adv_umx', zt_ups, 'T', 1.)
      IF (np_limiter == 1) THEN
        CALL nonosc_ice(1._wp, pdt, pu, pv, ptc, zt_ups, zfu_ups, zfv_ups, zfu_ho, zfv_ho)
      ELSE IF (np_limiter == 2 .OR. np_limiter == 3) THEN
        CALL limiter_x(pdt, pu, ptc, zfu_ups, zfu_ho)
        CALL limiter_y(pdt, pv, ptc, zfv_ups, zfv_ho)
      END IF
    END IF
    IF (PRESENT(pua_ho)) THEN
      DO jl = 1, jpl
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            pua_ho(ji, jj, jl) = zfu_ho(ji, jj, jl)
            pva_ho(ji, jj, jl) = zfv_ho(ji, jj, jl)
            pua_ups(ji, jj, jl) = zfu_ups(ji, jj, jl)
            pva_ups(ji, jj, jl) = zfv_ups(ji, jj, jl)
          END DO
        END DO
      END DO
    END IF
    DO jl = 1, jpl
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ztra = - (zfu_ho(ji, jj, jl) - zfu_ho(ji - 1, jj, jl) + zfv_ho(ji, jj, jl) - zfv_ho(ji, jj - 1, jl))
          ptc(ji, jj, jl) = (ptc(ji, jj, jl) + ztra * r1_e1e2t(ji, jj) * pdt) * tmask(ji, jj, 1)
        END DO
      END DO
    END DO
    CALL lbc_lnk('icedyn_adv_umx', ptc, 'T', 1.)
  END SUBROUTINE adv_umx
  SUBROUTINE upstream(pamsk, jt, kt, pdt, pt, pu, pv, pt_ups, pfu_ups, pfv_ups)
    REAL(KIND = wp), INTENT(IN) :: pamsk
    INTEGER, INTENT(IN) :: jt
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), INTENT(IN) :: pdt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pu, pv
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl), INTENT(OUT) :: pt_ups
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl), INTENT(OUT) :: pfu_ups, pfv_ups
    INTEGER :: ji, jj, jl
    REAL(KIND = wp) :: ztra
    IF (.NOT. ll_upsxy) THEN
      DO jl = 1, jpl
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            pfu_ups(ji, jj, jl) = MAX(pu(ji, jj), 0._wp) * pt(ji, jj, jl) + MIN(pu(ji, jj), 0._wp) * pt(ji + 1, jj, jl)
            pfv_ups(ji, jj, jl) = MAX(pv(ji, jj), 0._wp) * pt(ji, jj, jl) + MIN(pv(ji, jj), 0._wp) * pt(ji, jj + 1, jl)
          END DO
        END DO
      END DO
    ELSE
      IF (MOD((kt - 1) / nn_fsbc, 2) == MOD((jt - 1), 2)) THEN
        DO jl = 1, jpl
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              pfu_ups(ji, jj, jl) = MAX(pu(ji, jj), 0._wp) * pt(ji, jj, jl) + MIN(pu(ji, jj), 0._wp) * pt(ji + 1, jj, jl)
            END DO
          END DO
        END DO
        DO jl = 1, jpl
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              ztra = - (pfu_ups(ji, jj, jl) - pfu_ups(ji - 1, jj, jl)) + (pu(ji, jj) - pu(ji - 1, jj)) * pt(ji, jj, jl) * (1. - &
&pamsk)
              zpt(ji, jj, jl) = (pt(ji, jj, jl) + ztra * pdt * r1_e1e2t(ji, jj)) * tmask(ji, jj, 1)
            END DO
          END DO
        END DO
        CALL lbc_lnk('icedyn_adv_umx', zpt, 'T', 1.)
        DO jl = 1, jpl
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              pfv_ups(ji, jj, jl) = MAX(pv(ji, jj), 0._wp) * zpt(ji, jj, jl) + MIN(pv(ji, jj), 0._wp) * zpt(ji, jj + 1, jl)
            END DO
          END DO
        END DO
      ELSE
        DO jl = 1, jpl
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              pfv_ups(ji, jj, jl) = MAX(pv(ji, jj), 0._wp) * pt(ji, jj, jl) + MIN(pv(ji, jj), 0._wp) * pt(ji, jj + 1, jl)
            END DO
          END DO
        END DO
        DO jl = 1, jpl
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              ztra = - (pfv_ups(ji, jj, jl) - pfv_ups(ji, jj - 1, jl)) + (pv(ji, jj) - pv(ji, jj - 1)) * pt(ji, jj, jl) * (1. - &
&pamsk)
              zpt(ji, jj, jl) = (pt(ji, jj, jl) + ztra * pdt * r1_e1e2t(ji, jj)) * tmask(ji, jj, 1)
            END DO
          END DO
        END DO
        CALL lbc_lnk('icedyn_adv_umx', zpt, 'T', 1.)
        DO jl = 1, jpl
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              pfu_ups(ji, jj, jl) = MAX(pu(ji, jj), 0._wp) * zpt(ji, jj, jl) + MIN(pu(ji, jj), 0._wp) * zpt(ji + 1, jj, jl)
            END DO
          END DO
        END DO
      END IF
    END IF
    DO jl = 1, jpl
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ztra = - (pfu_ups(ji, jj, jl) - pfu_ups(ji - 1, jj, jl) + pfv_ups(ji, jj, jl) - pfv_ups(ji, jj - 1, jl)) + (pu(ji, jj) - &
&pu(ji - 1, jj) + pv(ji, jj) - pv(ji, jj - 1)) * pt(ji, jj, jl) * (1. - pamsk)
          pt_ups(ji, jj, jl) = (pt(ji, jj, jl) + ztra * pdt * r1_e1e2t(ji, jj)) * tmask(ji, jj, 1)
        END DO
      END DO
    END DO
    CALL lbc_lnk('icedyn_adv_umx', pt_ups, 'T', 1.)
  END SUBROUTINE upstream
  SUBROUTINE cen2(pamsk, jt, kt, pdt, pt, pu, pv, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho)
    REAL(KIND = wp), INTENT(IN) :: pamsk
    INTEGER, INTENT(IN) :: jt
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), INTENT(IN) :: pdt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pu, pv
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pt_ups
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pfu_ups, pfv_ups
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl), INTENT(OUT) :: pfu_ho, pfv_ho
    INTEGER :: ji, jj, jl
    REAL(KIND = wp) :: ztra
    IF (.NOT. ll_hoxy) THEN
      DO jl = 1, jpl
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            pfu_ho(ji, jj, jl) = 0.5_wp * pu(ji, jj) * (pt(ji, jj, jl) + pt(ji + 1, jj, jl))
            pfv_ho(ji, jj, jl) = 0.5_wp * pv(ji, jj) * (pt(ji, jj, jl) + pt(ji, jj + 1, jl))
          END DO
        END DO
      END DO
      IF (np_limiter == 1) THEN
        CALL nonosc_ice(pamsk, pdt, pu, pv, pt, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho)
      ELSE IF (np_limiter == 2 .OR. np_limiter == 3) THEN
        CALL limiter_x(pdt, pu, pt, pfu_ups, pfu_ho)
        CALL limiter_y(pdt, pv, pt, pfv_ups, pfv_ho)
      END IF
    ELSE
      IF (MOD((kt - 1) / nn_fsbc, 2) == MOD((jt - 1), 2)) THEN
        DO jl = 1, jpl
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              pfu_ho(ji, jj, jl) = 0.5_wp * pu(ji, jj) * (pt(ji, jj, jl) + pt(ji + 1, jj, jl))
            END DO
          END DO
        END DO
        IF (np_limiter == 2 .OR. np_limiter == 3) CALL limiter_x(pdt, pu, pt, pfu_ups, pfu_ho)
        DO jl = 1, jpl
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              ztra = - (pfu_ho(ji, jj, jl) - pfu_ho(ji - 1, jj, jl)) + (pu(ji, jj) - pu(ji - 1, jj)) * pt(ji, jj, jl) * (1. - pamsk)
              zpt(ji, jj, jl) = (pt(ji, jj, jl) + ztra * pdt * r1_e1e2t(ji, jj)) * tmask(ji, jj, 1)
            END DO
          END DO
        END DO
        CALL lbc_lnk('icedyn_adv_umx', zpt, 'T', 1.)
        DO jl = 1, jpl
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              pfv_ho(ji, jj, jl) = 0.5_wp * pv(ji, jj) * (zpt(ji, jj, jl) + zpt(ji, jj + 1, jl))
            END DO
          END DO
        END DO
        IF (np_limiter == 2 .OR. np_limiter == 3) CALL limiter_y(pdt, pv, pt, pfv_ups, pfv_ho)
      ELSE
        DO jl = 1, jpl
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              pfv_ho(ji, jj, jl) = 0.5_wp * pv(ji, jj) * (pt(ji, jj, jl) + pt(ji, jj + 1, jl))
            END DO
          END DO
        END DO
        IF (np_limiter == 2 .OR. np_limiter == 3) CALL limiter_y(pdt, pv, pt, pfv_ups, pfv_ho)
        DO jl = 1, jpl
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              ztra = - (pfv_ho(ji, jj, jl) - pfv_ho(ji, jj - 1, jl)) + (pv(ji, jj) - pv(ji, jj - 1)) * pt(ji, jj, jl) * (1. - pamsk)
              zpt(ji, jj, jl) = (pt(ji, jj, jl) + ztra * pdt * r1_e1e2t(ji, jj)) * tmask(ji, jj, 1)
            END DO
          END DO
        END DO
        CALL lbc_lnk('icedyn_adv_umx', zpt, 'T', 1.)
        DO jl = 1, jpl
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              pfu_ho(ji, jj, jl) = 0.5_wp * pu(ji, jj) * (zpt(ji, jj, jl) + zpt(ji + 1, jj, jl))
            END DO
          END DO
        END DO
        IF (np_limiter == 2 .OR. np_limiter == 3) CALL limiter_x(pdt, pu, pt, pfu_ups, pfu_ho)
      END IF
      IF (np_limiter == 1) CALL nonosc_ice(pamsk, pdt, pu, pv, pt, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho)
    END IF
  END SUBROUTINE cen2
  SUBROUTINE macho(pamsk, kn_umx, jt, kt, pdt, pt, pu, pv, pubox, pvbox, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho)
    REAL(KIND = wp), INTENT(IN) :: pamsk
    INTEGER, INTENT(IN) :: kn_umx
    INTEGER, INTENT(IN) :: jt
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), INTENT(IN) :: pdt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pu, pv
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pubox, pvbox
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pt_ups
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pfu_ups, pfv_ups
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl), INTENT(OUT) :: pfu_ho, pfv_ho
    INTEGER :: ji, jj, jl
    IF (MOD((kt - 1) / nn_fsbc, 2) == MOD((jt - 1), 2)) THEN
      CALL ultimate_x(pamsk, kn_umx, pdt, pt, pu, zt_u, pfu_ho)
      IF (np_limiter == 2 .OR. np_limiter == 3) CALL limiter_x(pdt, pu, pt, pfu_ups, pfu_ho)
      DO jl = 1, jpl
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zpt(ji, jj, jl) = (pt(ji, jj, jl) - (pubox(ji, jj) * (zt_u(ji, jj, jl) - zt_u(ji - 1, jj, jl)) * r1_e1t(ji, jj) + &
&pt(ji, jj, jl) * (pu(ji, jj) - pu(ji - 1, jj)) * r1_e1e2t(ji, jj) * pamsk) * pdt) * tmask(ji, jj, 1)
          END DO
        END DO
      END DO
      CALL lbc_lnk('icedyn_adv_umx', zpt, 'T', 1.)
      IF (ll_hoxy) THEN
        CALL ultimate_y(pamsk, kn_umx, pdt, zpt, pv, zt_v, pfv_ho)
      ELSE
        CALL ultimate_y(pamsk, kn_umx, pdt, pt, pv, zt_v, pfv_ho)
      END IF
      IF (np_limiter == 2 .OR. np_limiter == 3) CALL limiter_y(pdt, pv, pt, pfv_ups, pfv_ho)
    ELSE
      CALL ultimate_y(pamsk, kn_umx, pdt, pt, pv, zt_v, pfv_ho)
      IF (np_limiter == 2 .OR. np_limiter == 3) CALL limiter_y(pdt, pv, pt, pfv_ups, pfv_ho)
      DO jl = 1, jpl
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zpt(ji, jj, jl) = (pt(ji, jj, jl) - (pvbox(ji, jj) * (zt_v(ji, jj, jl) - zt_v(ji, jj - 1, jl)) * r1_e2t(ji, jj) + &
&pt(ji, jj, jl) * (pv(ji, jj) - pv(ji, jj - 1)) * r1_e1e2t(ji, jj) * pamsk) * pdt) * tmask(ji, jj, 1)
          END DO
        END DO
      END DO
      CALL lbc_lnk('icedyn_adv_umx', zpt, 'T', 1.)
      IF (ll_hoxy) THEN
        CALL ultimate_x(pamsk, kn_umx, pdt, zpt, pu, zt_u, pfu_ho)
      ELSE
        CALL ultimate_x(pamsk, kn_umx, pdt, pt, pu, zt_u, pfu_ho)
      END IF
      IF (np_limiter == 2 .OR. np_limiter == 3) CALL limiter_x(pdt, pu, pt, pfu_ups, pfu_ho)
    END IF
    IF (np_limiter == 1) CALL nonosc_ice(pamsk, pdt, pu, pv, pt, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho)
  END SUBROUTINE macho
  SUBROUTINE ultimate_x(pamsk, kn_umx, pdt, pt, pu, pt_u, pfu_ho)
    REAL(KIND = wp), INTENT(IN) :: pamsk
    INTEGER, INTENT(IN) :: kn_umx
    REAL(KIND = wp), INTENT(IN) :: pdt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pu
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl), INTENT(OUT) :: pt_u
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl), INTENT(OUT) :: pfu_ho
    INTEGER :: ji, jj, jl
    REAL(KIND = wp) :: zcu, zdx2, zdx4
    DO jl = 1, jpl
      DO jj = 2, jpjm1
        DO ji = 1, jpim1
          ztu1(ji, jj, jl) = (pt(ji + 1, jj, jl) - pt(ji, jj, jl)) * r1_e1u(ji, jj) * umask(ji, jj, 1)
        END DO
        DO ji = 2, jpim1
          ztu2(ji, jj, jl) = (ztu1(ji, jj, jl) - ztu1(ji - 1, jj, jl)) * r1_e1t(ji, jj)
        END DO
      END DO
    END DO
    CALL lbc_lnk('icedyn_adv_umx', ztu2, 'T', 1.)
    DO jl = 1, jpl
      DO jj = 2, jpjm1
        DO ji = 1, jpim1
          ztu3(ji, jj, jl) = (ztu2(ji + 1, jj, jl) - ztu2(ji, jj, jl)) * r1_e1u(ji, jj) * umask(ji, jj, 1)
        END DO
        DO ji = 2, jpim1
          ztu4(ji, jj, jl) = (ztu3(ji, jj, jl) - ztu3(ji - 1, jj, jl)) * r1_e1t(ji, jj)
        END DO
      END DO
    END DO
    CALL lbc_lnk('icedyn_adv_umx', ztu4, 'T', 1.)
    SELECT CASE (kn_umx)
    CASE (1)
      DO jl = 1, jpl
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            pt_u(ji, jj, jl) = 0.5_wp * umask(ji, jj, 1) * (pt(ji + 1, jj, jl) + pt(ji, jj, jl) - SIGN(1._wp, pu(ji, jj)) * (pt(ji &
&+ 1, jj, jl) - pt(ji, jj, jl)))
          END DO
        END DO
      END DO
    CASE (2)
      DO jl = 1, jpl
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zcu = pu(ji, jj) * r1_e2u(ji, jj) * pdt * r1_e1u(ji, jj)
            pt_u(ji, jj, jl) = 0.5_wp * umask(ji, jj, 1) * (pt(ji + 1, jj, jl) + pt(ji, jj, jl) - zcu * (pt(ji + 1, jj, jl) - &
&pt(ji, jj, jl)))
          END DO
        END DO
      END DO
    CASE (3)
      DO jl = 1, jpl
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zcu = pu(ji, jj) * r1_e2u(ji, jj) * pdt * r1_e1u(ji, jj)
            zdx2 = e1u(ji, jj) * e1u(ji, jj)
            pt_u(ji, jj, jl) = 0.5_wp * umask(ji, jj, 1) * ((pt(ji + 1, jj, jl) + pt(ji, jj, jl) - zcu * (pt(ji + 1, jj, jl) - &
&pt(ji, jj, jl))) + z1_6 * zdx2 * (zcu * zcu - 1._wp) * (ztu2(ji + 1, jj, jl) + ztu2(ji, jj, jl) - SIGN(1._wp, zcu) * (ztu2(ji + &
&1, jj, jl) - ztu2(ji, jj, jl))))
          END DO
        END DO
      END DO
    CASE (4)
      DO jl = 1, jpl
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zcu = pu(ji, jj) * r1_e2u(ji, jj) * pdt * r1_e1u(ji, jj)
            zdx2 = e1u(ji, jj) * e1u(ji, jj)
            pt_u(ji, jj, jl) = 0.5_wp * umask(ji, jj, 1) * ((pt(ji + 1, jj, jl) + pt(ji, jj, jl) - zcu * (pt(ji + 1, jj, jl) - &
&pt(ji, jj, jl))) + z1_6 * zdx2 * (zcu * zcu - 1._wp) * (ztu2(ji + 1, jj, jl) + ztu2(ji, jj, jl) - 0.5_wp * zcu * (ztu2(ji + 1, &
&jj, jl) - ztu2(ji, jj, jl))))
          END DO
        END DO
      END DO
    CASE (5)
      DO jl = 1, jpl
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zcu = pu(ji, jj) * r1_e2u(ji, jj) * pdt * r1_e1u(ji, jj)
            zdx2 = e1u(ji, jj) * e1u(ji, jj)
            zdx4 = zdx2 * zdx2
            pt_u(ji, jj, jl) = 0.5_wp * umask(ji, jj, 1) * ((pt(ji + 1, jj, jl) + pt(ji, jj, jl) - zcu * (pt(ji + 1, jj, jl) - &
&pt(ji, jj, jl))) + z1_6 * zdx2 * (zcu * zcu - 1._wp) * (ztu2(ji + 1, jj, jl) + ztu2(ji, jj, jl) - 0.5_wp * zcu * (ztu2(ji + 1, &
&jj, jl) - ztu2(ji, jj, jl))) + z1_120 * zdx4 * (zcu * zcu - 1._wp) * (zcu * zcu - 4._wp) * (ztu4(ji + 1, jj, jl) + ztu4(ji, jj, &
&jl) - SIGN(1._wp, zcu) * (ztu4(ji + 1, jj, jl) - ztu4(ji, jj, jl))))
          END DO
        END DO
      END DO
    END SELECT
    IF (ll_neg) THEN
      DO jl = 1, jpl
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            IF (pt_u(ji, jj, jl) < 0._wp .OR. (imsk_small(ji, jj, jl) == 0 .AND. pamsk == 0.)) THEN
              pt_u(ji, jj, jl) = 0.5_wp * umask(ji, jj, 1) * (pt(ji + 1, jj, jl) + pt(ji, jj, jl) - SIGN(1._wp, pu(ji, jj)) * &
&(pt(ji + 1, jj, jl) - pt(ji, jj, jl)))
            END IF
          END DO
        END DO
      END DO
    END IF
    DO jl = 1, jpl
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          pfu_ho(ji, jj, jl) = pu(ji, jj) * pt_u(ji, jj, jl)
        END DO
      END DO
    END DO
  END SUBROUTINE ultimate_x
  SUBROUTINE ultimate_y(pamsk, kn_umx, pdt, pt, pv, pt_v, pfv_ho)
    REAL(KIND = wp), INTENT(IN) :: pamsk
    INTEGER, INTENT(IN) :: kn_umx
    REAL(KIND = wp), INTENT(IN) :: pdt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pv
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl), INTENT(OUT) :: pt_v
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl), INTENT(OUT) :: pfv_ho
    INTEGER :: ji, jj, jl
    REAL(KIND = wp) :: zcv, zdy2, zdy4
    DO jl = 1, jpl
      DO jj = 1, jpjm1
        DO ji = 2, jpim1
          ztv1(ji, jj, jl) = (pt(ji, jj + 1, jl) - pt(ji, jj, jl)) * r1_e2v(ji, jj) * vmask(ji, jj, 1)
        END DO
      END DO
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ztv2(ji, jj, jl) = (ztv1(ji, jj, jl) - ztv1(ji, jj - 1, jl)) * r1_e2t(ji, jj)
        END DO
      END DO
    END DO
    CALL lbc_lnk('icedyn_adv_umx', ztv2, 'T', 1.)
    DO jl = 1, jpl
      DO jj = 1, jpjm1
        DO ji = 2, jpim1
          ztv3(ji, jj, jl) = (ztv2(ji, jj + 1, jl) - ztv2(ji, jj, jl)) * r1_e2v(ji, jj) * vmask(ji, jj, 1)
        END DO
      END DO
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ztv4(ji, jj, jl) = (ztv3(ji, jj, jl) - ztv3(ji, jj - 1, jl)) * r1_e2t(ji, jj)
        END DO
      END DO
    END DO
    CALL lbc_lnk('icedyn_adv_umx', ztv4, 'T', 1.)
    SELECT CASE (kn_umx)
    CASE (1)
      DO jl = 1, jpl
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            pt_v(ji, jj, jl) = 0.5_wp * vmask(ji, jj, 1) * (pt(ji, jj + 1, jl) + pt(ji, jj, jl) - SIGN(1._wp, pv(ji, jj)) * &
&(pt(ji, jj + 1, jl) - pt(ji, jj, jl)))
          END DO
        END DO
      END DO
    CASE (2)
      DO jl = 1, jpl
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zcv = pv(ji, jj) * r1_e1v(ji, jj) * pdt * r1_e2v(ji, jj)
            pt_v(ji, jj, jl) = 0.5_wp * vmask(ji, jj, 1) * (pt(ji, jj + 1, jl) + pt(ji, jj, jl) - zcv * (pt(ji, jj + 1, jl) - &
&pt(ji, jj, jl)))
          END DO
        END DO
      END DO
    CASE (3)
      DO jl = 1, jpl
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zcv = pv(ji, jj) * r1_e1v(ji, jj) * pdt * r1_e2v(ji, jj)
            zdy2 = e2v(ji, jj) * e2v(ji, jj)
            pt_v(ji, jj, jl) = 0.5_wp * vmask(ji, jj, 1) * ((pt(ji, jj + 1, jl) + pt(ji, jj, jl) - zcv * (pt(ji, jj + 1, jl) - &
&pt(ji, jj, jl))) + z1_6 * zdy2 * (zcv * zcv - 1._wp) * (ztv2(ji, jj + 1, jl) + ztv2(ji, jj, jl) - SIGN(1._wp, zcv) * (ztv2(ji, jj &
&+ 1, jl) - ztv2(ji, jj, jl))))
          END DO
        END DO
      END DO
    CASE (4)
      DO jl = 1, jpl
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zcv = pv(ji, jj) * r1_e1v(ji, jj) * pdt * r1_e2v(ji, jj)
            zdy2 = e2v(ji, jj) * e2v(ji, jj)
            pt_v(ji, jj, jl) = 0.5_wp * vmask(ji, jj, 1) * ((pt(ji, jj + 1, jl) + pt(ji, jj, jl) - zcv * (pt(ji, jj + 1, jl) - &
&pt(ji, jj, jl))) + z1_6 * zdy2 * (zcv * zcv - 1._wp) * (ztv2(ji, jj + 1, jl) + ztv2(ji, jj, jl) - 0.5_wp * zcv * (ztv2(ji, jj + &
&1, jl) - ztv2(ji, jj, jl))))
          END DO
        END DO
      END DO
    CASE (5)
      DO jl = 1, jpl
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zcv = pv(ji, jj) * r1_e1v(ji, jj) * pdt * r1_e2v(ji, jj)
            zdy2 = e2v(ji, jj) * e2v(ji, jj)
            zdy4 = zdy2 * zdy2
            pt_v(ji, jj, jl) = 0.5_wp * vmask(ji, jj, 1) * ((pt(ji, jj + 1, jl) + pt(ji, jj, jl) - zcv * (pt(ji, jj + 1, jl) - &
&pt(ji, jj, jl))) + z1_6 * zdy2 * (zcv * zcv - 1._wp) * (ztv2(ji, jj + 1, jl) + ztv2(ji, jj, jl) - 0.5_wp * zcv * (ztv2(ji, jj + &
&1, jl) - ztv2(ji, jj, jl))) + z1_120 * zdy4 * (zcv * zcv - 1._wp) * (zcv * zcv - 4._wp) * (ztv4(ji, jj + 1, jl) + ztv4(ji, jj, &
&jl) - SIGN(1._wp, zcv) * (ztv4(ji, jj + 1, jl) - ztv4(ji, jj, jl))))
          END DO
        END DO
      END DO
    END SELECT
    IF (ll_neg) THEN
      DO jl = 1, jpl
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            IF (pt_v(ji, jj, jl) < 0._wp .OR. (jmsk_small(ji, jj, jl) == 0 .AND. pamsk == 0.)) THEN
              pt_v(ji, jj, jl) = 0.5_wp * vmask(ji, jj, 1) * ((pt(ji, jj + 1, jl) + pt(ji, jj, jl)) - SIGN(1._wp, pv(ji, jj)) * &
&(pt(ji, jj + 1, jl) - pt(ji, jj, jl)))
            END IF
          END DO
        END DO
      END DO
    END IF
    DO jl = 1, jpl
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          pfv_ho(ji, jj, jl) = pv(ji, jj) * pt_v(ji, jj, jl)
        END DO
      END DO
    END DO
  END SUBROUTINE ultimate_y
  SUBROUTINE nonosc_ice(pamsk, pdt, pu, pv, pt, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho)
    REAL(KIND = wp), INTENT(IN) :: pamsk
    REAL(KIND = wp), INTENT(IN) :: pdt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pu
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pv
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pt, pt_ups
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pfv_ups, pfu_ups
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pfv_ho, pfu_ho
    INTEGER :: ji, jj, jl
    REAL(KIND = wp) :: zpos, zneg, zbig, zup, zdo, z1_dt
    REAL(KIND = wp) :: zau, zbu, zcu, zav, zbv, zcv, zcoef, zzt
    zbig = 1.E+40_wp
    DO jl = 1, jpl
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          pfu_ho(ji, jj, jl) = pfu_ho(ji, jj, jl) - pfu_ups(ji, jj, jl)
          pfv_ho(ji, jj, jl) = pfv_ho(ji, jj, jl) - pfv_ups(ji, jj, jl)
        END DO
      END DO
    END DO
    IF (ll_prelim) THEN
      DO jl = 1, jpl
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zti_ups(ji, jj, jl) = pt_ups(ji + 1, jj, jl)
            ztj_ups(ji, jj, jl) = pt_ups(ji, jj + 1, jl)
          END DO
        END DO
      END DO
      CALL lbc_lnk_multi('icedyn_adv_umx', zti_ups, 'T', 1., ztj_ups, 'T', 1.)
      DO jl = 1, jpl
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            IF (pfu_ho(ji, jj, jl) * (pt_ups(ji + 1, jj, jl) - pt_ups(ji, jj, jl)) <= 0._wp .AND. pfv_ho(ji, jj, jl) * (pt_ups(ji, &
&jj + 1, jl) - pt_ups(ji, jj, jl)) <= 0._wp) THEN
              IF (pfu_ho(ji, jj, jl) * (zti_ups(ji + 1, jj, jl) - zti_ups(ji, jj, jl)) <= 0._wp .AND. pfv_ho(ji, jj, jl) * &
&(ztj_ups(ji, jj + 1, jl) - ztj_ups(ji, jj, jl)) <= 0._wp) THEN
                pfu_ho(ji, jj, jl) = 0._wp
                pfv_ho(ji, jj, jl) = 0._wp
              END IF
              IF (pfu_ho(ji, jj, jl) * (pt_ups(ji, jj, jl) - pt_ups(ji - 1, jj, jl)) <= 0._wp .AND. pfv_ho(ji, jj, jl) * &
&(pt_ups(ji, jj, jl) - pt_ups(ji, jj - 1, jl)) <= 0._wp) THEN
                pfu_ho(ji, jj, jl) = 0._wp
                pfv_ho(ji, jj, jl) = 0._wp
              END IF
            END IF
          END DO
        END DO
      END DO
      CALL lbc_lnk_multi('icedyn_adv_umx', pfu_ho, 'U', - 1., pfv_ho, 'V', - 1.)
    END IF
    z1_dt = 1._wp / pdt
    DO jl = 1, jpl
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (pt(ji, jj, jl) <= 0._wp .AND. pt_ups(ji, jj, jl) <= 0._wp) THEN
            zbup(ji, jj) = - zbig
            zbdo(ji, jj) = zbig
          ELSE IF (pt(ji, jj, jl) <= 0._wp .AND. pt_ups(ji, jj, jl) > 0._wp) THEN
            zbup(ji, jj) = pt_ups(ji, jj, jl)
            zbdo(ji, jj) = pt_ups(ji, jj, jl)
          ELSE IF (pt(ji, jj, jl) > 0._wp .AND. pt_ups(ji, jj, jl) <= 0._wp) THEN
            zbup(ji, jj) = pt(ji, jj, jl)
            zbdo(ji, jj) = pt(ji, jj, jl)
          ELSE
            zbup(ji, jj) = MAX(pt(ji, jj, jl), pt_ups(ji, jj, jl))
            zbdo(ji, jj) = MIN(pt(ji, jj, jl), pt_ups(ji, jj, jl))
          END IF
        END DO
      END DO
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zup = MAX(zbup(ji, jj), zbup(ji - 1, jj), zbup(ji + 1, jj), zbup(ji, jj - 1), zbup(ji, jj + 1))
          zdo = MIN(zbdo(ji, jj), zbdo(ji - 1, jj), zbdo(ji + 1, jj), zbdo(ji, jj - 1), zbdo(ji, jj + 1))
          zpos = MAX(0._wp, pfu_ho(ji - 1, jj, jl)) - MIN(0._wp, pfu_ho(ji, jj, jl)) + MAX(0._wp, pfv_ho(ji, jj - 1, jl)) - &
&MIN(0._wp, pfv_ho(ji, jj, jl))
          zneg = MAX(0._wp, pfu_ho(ji, jj, jl)) - MIN(0._wp, pfu_ho(ji - 1, jj, jl)) + MAX(0._wp, pfv_ho(ji, jj, jl)) - MIN(0._wp, &
&pfv_ho(ji, jj - 1, jl))
          zpos = zpos - (pt(ji, jj, jl) * MIN(0., pu(ji, jj) - pu(ji - 1, jj)) + pt(ji, jj, jl) * MIN(0., pv(ji, jj) - pv(ji, jj - &
&1))) * (1. - pamsk)
          zneg = zneg + (pt(ji, jj, jl) * MAX(0., pu(ji, jj) - pu(ji - 1, jj)) + pt(ji, jj, jl) * MAX(0., pv(ji, jj) - pv(ji, jj - &
&1))) * (1. - pamsk)
          IF (zpos > epsi10) THEN
            zbetup(ji, jj, jl) = MAX(0._wp, zup - pt_ups(ji, jj, jl)) / zpos * e1e2t(ji, jj) * z1_dt
          ELSE
            zbetup(ji, jj, jl) = 0._wp
          END IF
          IF (zneg > epsi10) THEN
            zbetdo(ji, jj, jl) = MAX(0._wp, pt_ups(ji, jj, jl) - zdo) / zneg * e1e2t(ji, jj) * z1_dt
          ELSE
            zbetdo(ji, jj, jl) = 0._wp
          END IF
          IF (zup == - zbig) zbetup(ji, jj, jl) = 0._wp
          IF (zdo == zbig) zbetdo(ji, jj, jl) = 0._wp
        END DO
      END DO
    END DO
    CALL lbc_lnk_multi('icedyn_adv_umx', zbetup, 'T', 1., zbetdo, 'T', 1.)
    DO jl = 1, jpl
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zau = MIN(1._wp, zbetdo(ji, jj, jl), zbetup(ji + 1, jj, jl))
          zbu = MIN(1._wp, zbetup(ji, jj, jl), zbetdo(ji + 1, jj, jl))
          zcu = 0.5_wp + SIGN(0.5_wp, pfu_ho(ji, jj, jl))
          zcoef = (zcu * zau + (1._wp - zcu) * zbu)
          pfu_ho(ji, jj, jl) = pfu_ho(ji, jj, jl) * zcoef + pfu_ups(ji, jj, jl)
        END DO
      END DO
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zav = MIN(1._wp, zbetdo(ji, jj, jl), zbetup(ji, jj + 1, jl))
          zbv = MIN(1._wp, zbetup(ji, jj, jl), zbetdo(ji, jj + 1, jl))
          zcv = 0.5_wp + SIGN(0.5_wp, pfv_ho(ji, jj, jl))
          zcoef = (zcv * zav + (1._wp - zcv) * zbv)
          pfv_ho(ji, jj, jl) = pfv_ho(ji, jj, jl) * zcoef + pfv_ups(ji, jj, jl)
        END DO
      END DO
    END DO
  END SUBROUTINE nonosc_ice
  SUBROUTINE limiter_x(pdt, pu, pt, pfu_ups, pfu_ho)
    REAL(KIND = wp), INTENT(IN) :: pdt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pu
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pfu_ups
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pfu_ho
    REAL(KIND = wp) :: Cr, Rjm, Rj, Rjp, uCFL, zpsi, zh3, zlimiter, Rr
    INTEGER :: ji, jj, jl
    DO jl = 1, jpl
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zslpx(ji, jj, jl) = (pt(ji + 1, jj, jl) - pt(ji, jj, jl)) * umask(ji, jj, 1)
        END DO
      END DO
    END DO
    CALL lbc_lnk('icedyn_adv_umx', zslpx, 'U', - 1.)
    DO jl = 1, jpl
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          uCFL = pdt * ABS(pu(ji, jj)) * r1_e1e2t(ji, jj)
          Rjm = zslpx(ji - 1, jj, jl)
          Rj = zslpx(ji, jj, jl)
          Rjp = zslpx(ji + 1, jj, jl)
          IF (np_limiter == 3) THEN
            IF (pu(ji, jj) > 0.) THEN
              rr = rjm
            ELSE
              rr = rjp
            END IF
            zh3 = pfu_ho(ji, jj, jl) - pfu_ups(ji, jj, jl)
            IF (Rj > 0.) THEN
              zlimiter = MAX(0., MIN(zh3, MAX(- Rr * 0.5 * ABS(pu(ji, jj)), MIN(2. * Rr * 0.5 * ABS(pu(ji, jj)), zh3, 1.5 * Rj * &
&0.5 * ABS(pu(ji, jj))))))
            ELSE
              zlimiter = - MAX(0., MIN(- zh3, MAX(Rr * 0.5 * ABS(pu(ji, jj)), MIN(- 2. * Rr * 0.5 * ABS(pu(ji, jj)), - zh3, - 1.5 &
&* Rj * 0.5 * ABS(pu(ji, jj))))))
            END IF
            pfu_ho(ji, jj, jl) = pfu_ups(ji, jj, jl) + zlimiter
          ELSE IF (np_limiter == 2) THEN
            IF (Rj /= 0.) THEN
              IF (pu(ji, jj) > 0.) THEN
                cr = rjm / rj
              ELSE
                cr = rjp / rj
              END IF
            ELSE
              Cr = 0.
            END IF
            zpsi = MAX(0., MAX(MIN(1., 2. * Cr), MIN(2., Cr)))
            pfu_ho(ji, jj, jl) = pfu_ho(ji, jj, jl) - ABS(pu(ji, jj)) * ((1. - zpsi) + uCFL * zpsi) * Rj * 0.5
          END IF
        END DO
      END DO
    END DO
    CALL lbc_lnk('icedyn_adv_umx', pfu_ho, 'U', - 1.)
  END SUBROUTINE limiter_x
  SUBROUTINE limiter_y(pdt, pv, pt, pfv_ups, pfv_ho)
    REAL(KIND = wp), INTENT(IN) :: pdt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pv
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pfv_ups
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pfv_ho
    REAL(KIND = wp) :: Cr, Rjm, Rj, Rjp, vCFL, zpsi, zh3, zlimiter, Rr
    INTEGER :: ji, jj, jl
    DO jl = 1, jpl
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zslpy(ji, jj, jl) = (pt(ji, jj + 1, jl) - pt(ji, jj, jl)) * vmask(ji, jj, 1)
        END DO
      END DO
    END DO
    CALL lbc_lnk('icedyn_adv_umx', zslpy, 'V', - 1.)
    DO jl = 1, jpl
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          vCFL = pdt * ABS(pv(ji, jj)) * r1_e1e2t(ji, jj)
          Rjm = zslpy(ji, jj - 1, jl)
          Rj = zslpy(ji, jj, jl)
          Rjp = zslpy(ji, jj + 1, jl)
          IF (np_limiter == 3) THEN
            IF (pv(ji, jj) > 0.) THEN
              rr = rjm
            ELSE
              rr = rjp
            END IF
            zh3 = pfv_ho(ji, jj, jl) - pfv_ups(ji, jj, jl)
            IF (Rj > 0.) THEN
              zlimiter = MAX(0., MIN(zh3, MAX(- Rr * 0.5 * ABS(pv(ji, jj)), MIN(2. * Rr * 0.5 * ABS(pv(ji, jj)), zh3, 1.5 * Rj * &
&0.5 * ABS(pv(ji, jj))))))
            ELSE
              zlimiter = - MAX(0., MIN(- zh3, MAX(Rr * 0.5 * ABS(pv(ji, jj)), MIN(- 2. * Rr * 0.5 * ABS(pv(ji, jj)), - zh3, - 1.5 &
&* Rj * 0.5 * ABS(pv(ji, jj))))))
            END IF
            pfv_ho(ji, jj, jl) = pfv_ups(ji, jj, jl) + zlimiter
          ELSE IF (np_limiter == 2) THEN
            IF (Rj /= 0.) THEN
              IF (pv(ji, jj) > 0.) THEN
                cr = rjm / rj
              ELSE
                cr = rjp / rj
              END IF
            ELSE
              Cr = 0.
            END IF
            zpsi = MAX(0., MAX(MIN(1., 2. * Cr), MIN(2., Cr)))
            pfv_ho(ji, jj, jl) = pfv_ho(ji, jj, jl) - ABS(pv(ji, jj)) * ((1. - zpsi) + vCFL * zpsi) * Rj * 0.5
          END IF
        END DO
      END DO
    END DO
    CALL lbc_lnk('icedyn_adv_umx', pfv_ho, 'V', - 1.)
  END SUBROUTINE limiter_y
  SUBROUTINE Hbig(pdt, phi_max, phs_max, phip_max, pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, pe_s, pe_i)
    REAL(KIND = wp), INTENT(IN) :: pdt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: phi_max, phs_max, phip_max
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip
    REAL(KIND = wp), DIMENSION(:, :, :, :), INTENT(INOUT) :: pe_s
    REAL(KIND = wp), DIMENSION(:, :, :, :), INTENT(INOUT) :: pe_i
    INTEGER :: ji, jj, jk, jl
    REAL(KIND = wp) :: z1_dt, zhip, zhi, zhs, zvs_excess, zfra
    z1_dt = 1._wp / pdt
    DO jl = 1, jpl
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (pv_i(ji, jj, jl) > 0._wp) THEN
            IF (ln_pnd_H12 .AND. pv_ip(ji, jj, jl) > 0._wp) THEN
              zhip = pv_ip(ji, jj, jl) / MAX(epsi20, pa_ip(ji, jj, jl))
              IF (zhip > phip_max(ji, jj, jl) .AND. pa_ip(ji, jj, jl) < 0.15) THEN
                pa_ip(ji, jj, jl) = pv_ip(ji, jj, jl) / phip_max(ji, jj, jl)
              END IF
            END IF
            zhi = pv_i(ji, jj, jl) / pa_i(ji, jj, jl)
            IF (zhi > phi_max(ji, jj, jl) .AND. pa_i(ji, jj, jl) < 0.15) THEN
              pa_i(ji, jj, jl) = pv_i(ji, jj, jl) / MIN(phi_max(ji, jj, jl), hi_max(jpl))
            END IF
            zhs = pv_s(ji, jj, jl) / pa_i(ji, jj, jl)
            IF (pv_s(ji, jj, jl) > 0._wp .AND. zhs > phs_max(ji, jj, jl) .AND. pa_i(ji, jj, jl) < 0.15) THEN
              zfra = phs_max(ji, jj, jl) / MAX(zhs, epsi20)
              wfx_res(ji, jj) = wfx_res(ji, jj) + (pv_s(ji, jj, jl) - pa_i(ji, jj, jl) * phs_max(ji, jj, jl)) * rhos * z1_dt
              hfx_res(ji, jj) = hfx_res(ji, jj) - SUM(pe_s(ji, jj, 1 : nlay_s, jl)) * (1._wp - zfra) * z1_dt
              pe_s(ji, jj, 1 : nlay_s, jl) = pe_s(ji, jj, 1 : nlay_s, jl) * zfra
              pv_s(ji, jj, jl) = pa_i(ji, jj, jl) * phs_max(ji, jj, jl)
            END IF
            zvs_excess = MAX(0._wp, pv_s(ji, jj, jl) - pv_i(ji, jj, jl) * (rau0 - rhoi) * r1_rhos)
            IF (zvs_excess > 0._wp) THEN
              zfra = (pv_s(ji, jj, jl) - zvs_excess) / MAX(pv_s(ji, jj, jl), epsi20)
              wfx_res(ji, jj) = wfx_res(ji, jj) + zvs_excess * rhos * z1_dt
              hfx_res(ji, jj) = hfx_res(ji, jj) - SUM(pe_s(ji, jj, 1 : nlay_s, jl)) * (1._wp - zfra) * z1_dt
              pe_s(ji, jj, 1 : nlay_s, jl) = pe_s(ji, jj, 1 : nlay_s, jl) * zfra
              pv_s(ji, jj, jl) = pv_s(ji, jj, jl) - zvs_excess
            END IF
          END IF
        END DO
      END DO
    END DO
    WHERE (pa_ip(:, :, :) > pa_i(:, :, :)) pa_ip(:, :, :) = pa_i(:, :, :)
  END SUBROUTINE Hbig
END MODULE icedyn_adv_umx