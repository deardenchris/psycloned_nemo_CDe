MODULE icedyn_adv_pra
  USE dom_oce
  USE ice
  USE sbc_oce, ONLY: nn_fsbc
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  USE lbclnk
  USE prtctl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_dyn_adv_pra
  PUBLIC :: adv_pra_init
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: sxice, syice, sxxice, syyice, sxyice
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: sxsn, sysn, sxxsn, syysn, sxysn
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: sxa, sya, sxxa, syya, sxya
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: sxsal, sysal, sxxsal, syysal, sxysal
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: sxage, syage, sxxage, syyage, sxyage
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: sxopw, syopw, sxxopw, syyopw, sxyopw
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :, :) :: sxc0, syc0, sxxc0, syyc0, sxyc0
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :, :) :: sxe, sye, sxxe, syye, sxye
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: sxap, syap, sxxap, syyap, sxyap
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: sxvp, syvp, sxxvp, syyvp, sxyvp
  CONTAINS
  SUBROUTINE ice_dyn_adv_pra(kt, pu_ice, pv_ice, pato_i, pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, pe_s, pe_i)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pu_ice
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pv_ice
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
    INTEGER :: jk, jl, jt
    INTEGER :: initad
    REAL(KIND = wp) :: zcfl, zusnit
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zarea
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: z0opw
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: z0ice, z0snw, z0ai, z0smi, z0oi
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: z0ap, z0vp
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :, :) :: z0es
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :, :) :: z0ei
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ice_dyn_adv_pra', 'r0', 0, 0)
    IF (kt == nit000 .AND. lwp) WRITE(numout, *) '-- ice_dyn_adv_pra: Prather advection scheme'
    ALLOCATE(zarea(jpi, jpj), z0opw(jpi, jpj, 1), z0ice(jpi, jpj, jpl), z0snw(jpi, jpj, jpl), z0ai(jpi, jpj, jpl), z0smi(jpi, jpj, jpl), z0oi(jpi, jpj, jpl), z0ap(jpi, jpj, jpl), z0vp(jpi, jpj, jpl), z0es(jpi, jpj, nlay_s, jpl), z0ei(jpi, jpj, nlay_i, jpl))
    zcfl = MAXVAL(ABS(pu_ice(:, :)) * rdt_ice * r1_e1u(:, :))
    zcfl = MAX(zcfl, MAXVAL(ABS(pv_ice(:, :)) * rdt_ice * r1_e2v(:, :)))
    CALL mpp_max('icedyn_adv_pra', zcfl)
    IF (zcfl > 0.5) THEN
      initad = 2
      zusnit = 0.5_wp
    ELSE
      initad = 1
      zusnit = 1.0_wp
    END IF
    zarea(:, :) = e1e2t(:, :)
    z0opw(:, :, 1) = pato_i(:, :) * e1e2t(:, :)
    DO jl = 1, jpl
      z0snw(:, :, jl) = pv_s(:, :, jl) * e1e2t(:, :)
      z0ice(:, :, jl) = pv_i(:, :, jl) * e1e2t(:, :)
      z0ai(:, :, jl) = pa_i(:, :, jl) * e1e2t(:, :)
      z0smi(:, :, jl) = psv_i(:, :, jl) * e1e2t(:, :)
      z0oi(:, :, jl) = poa_i(:, :, jl) * e1e2t(:, :)
      DO jk = 1, nlay_s
        z0es(:, :, jk, jl) = pe_s(:, :, jk, jl) * e1e2t(:, :)
      END DO
      DO jk = 1, nlay_i
        z0ei(:, :, jk, jl) = pe_i(:, :, jk, jl) * e1e2t(:, :)
      END DO
      IF (ln_pnd_H12) THEN
        z0ap(:, :, jl) = pa_ip(:, :, jl) * e1e2t(:, :)
        z0vp(:, :, jl) = pv_ip(:, :, jl) * e1e2t(:, :)
      END IF
    END DO
    IF (MOD((kt - 1) / nn_fsbc, 2) == 0) THEN
      DO jt = 1, initad
        CALL adv_x(zusnit, pu_ice, 1._wp, zarea, z0opw(:, :, 1), sxopw(:, :), sxxopw(:, :), syopw(:, :), syyopw(:, :), sxyopw(:, :))
        CALL adv_y(zusnit, pv_ice, 0._wp, zarea, z0opw(:, :, 1), sxopw(:, :), sxxopw(:, :), syopw(:, :), syyopw(:, :), sxyopw(:, :))
        DO jl = 1, jpl
          CALL adv_x(zusnit, pu_ice, 1._wp, zarea, z0ice(:, :, jl), sxice(:, :, jl), sxxice(:, :, jl), syice(:, :, jl), syyice(:, :, jl), sxyice(:, :, jl))
          CALL adv_y(zusnit, pv_ice, 0._wp, zarea, z0ice(:, :, jl), sxice(:, :, jl), sxxice(:, :, jl), syice(:, :, jl), syyice(:, :, jl), sxyice(:, :, jl))
          CALL adv_x(zusnit, pu_ice, 1._wp, zarea, z0snw(:, :, jl), sxsn(:, :, jl), sxxsn(:, :, jl), sysn(:, :, jl), syysn(:, :, jl), sxysn(:, :, jl))
          CALL adv_y(zusnit, pv_ice, 0._wp, zarea, z0snw(:, :, jl), sxsn(:, :, jl), sxxsn(:, :, jl), sysn(:, :, jl), syysn(:, :, jl), sxysn(:, :, jl))
          CALL adv_x(zusnit, pu_ice, 1._wp, zarea, z0smi(:, :, jl), sxsal(:, :, jl), sxxsal(:, :, jl), sysal(:, :, jl), syysal(:, :, jl), sxysal(:, :, jl))
          CALL adv_y(zusnit, pv_ice, 0._wp, zarea, z0smi(:, :, jl), sxsal(:, :, jl), sxxsal(:, :, jl), sysal(:, :, jl), syysal(:, :, jl), sxysal(:, :, jl))
          CALL adv_x(zusnit, pu_ice, 1._wp, zarea, z0oi(:, :, jl), sxage(:, :, jl), sxxage(:, :, jl), syage(:, :, jl), syyage(:, :, jl), sxyage(:, :, jl))
          CALL adv_y(zusnit, pv_ice, 0._wp, zarea, z0oi(:, :, jl), sxage(:, :, jl), sxxage(:, :, jl), syage(:, :, jl), syyage(:, :, jl), sxyage(:, :, jl))
          CALL adv_x(zusnit, pu_ice, 1._wp, zarea, z0ai(:, :, jl), sxa(:, :, jl), sxxa(:, :, jl), sya(:, :, jl), syya(:, :, jl), sxya(:, :, jl))
          CALL adv_y(zusnit, pv_ice, 0._wp, zarea, z0ai(:, :, jl), sxa(:, :, jl), sxxa(:, :, jl), sya(:, :, jl), syya(:, :, jl), sxya(:, :, jl))
          DO jk = 1, nlay_s
            CALL adv_x(zusnit, pu_ice, 1._wp, zarea, z0es(:, :, jk, jl), sxc0(:, :, jk, jl), sxxc0(:, :, jk, jl), syc0(:, :, jk, jl), syyc0(:, :, jk, jl), sxyc0(:, :, jk, jl))
            CALL adv_y(zusnit, pv_ice, 0._wp, zarea, z0es(:, :, jk, jl), sxc0(:, :, jk, jl), sxxc0(:, :, jk, jl), syc0(:, :, jk, jl), syyc0(:, :, jk, jl), sxyc0(:, :, jk, jl))
          END DO
          DO jk = 1, nlay_i
            CALL adv_x(zusnit, pu_ice, 1._wp, zarea, z0ei(:, :, jk, jl), sxe(:, :, jk, jl), sxxe(:, :, jk, jl), sye(:, :, jk, jl), syye(:, :, jk, jl), sxye(:, :, jk, jl))
            CALL adv_y(zusnit, pv_ice, 0._wp, zarea, z0ei(:, :, jk, jl), sxe(:, :, jk, jl), sxxe(:, :, jk, jl), sye(:, :, jk, jl), syye(:, :, jk, jl), sxye(:, :, jk, jl))
          END DO
          IF (ln_pnd_H12) THEN
            CALL adv_x(zusnit, pu_ice, 1._wp, zarea, z0ap(:, :, jl), sxap(:, :, jl), sxxap(:, :, jl), syap(:, :, jl), syyap(:, :, jl), sxyap(:, :, jl))
            CALL adv_y(zusnit, pv_ice, 0._wp, zarea, z0ap(:, :, jl), sxap(:, :, jl), sxxap(:, :, jl), syap(:, :, jl), syyap(:, :, jl), sxyap(:, :, jl))
            CALL adv_x(zusnit, pu_ice, 1._wp, zarea, z0vp(:, :, jl), sxvp(:, :, jl), sxxvp(:, :, jl), syvp(:, :, jl), syyvp(:, :, jl), sxyvp(:, :, jl))
            CALL adv_y(zusnit, pv_ice, 0._wp, zarea, z0vp(:, :, jl), sxvp(:, :, jl), sxxvp(:, :, jl), syvp(:, :, jl), syyvp(:, :, jl), sxyvp(:, :, jl))
          END IF
        END DO
      END DO
    ELSE
      DO jt = 1, initad
        CALL adv_y(zusnit, pv_ice, 1._wp, zarea, z0opw(:, :, 1), sxopw(:, :), sxxopw(:, :), syopw(:, :), syyopw(:, :), sxyopw(:, :))
        CALL adv_x(zusnit, pu_ice, 0._wp, zarea, z0opw(:, :, 1), sxopw(:, :), sxxopw(:, :), syopw(:, :), syyopw(:, :), sxyopw(:, :))
        DO jl = 1, jpl
          CALL adv_y(zusnit, pv_ice, 1._wp, zarea, z0ice(:, :, jl), sxice(:, :, jl), sxxice(:, :, jl), syice(:, :, jl), syyice(:, :, jl), sxyice(:, :, jl))
          CALL adv_x(zusnit, pu_ice, 0._wp, zarea, z0ice(:, :, jl), sxice(:, :, jl), sxxice(:, :, jl), syice(:, :, jl), syyice(:, :, jl), sxyice(:, :, jl))
          CALL adv_y(zusnit, pv_ice, 1._wp, zarea, z0snw(:, :, jl), sxsn(:, :, jl), sxxsn(:, :, jl), sysn(:, :, jl), syysn(:, :, jl), sxysn(:, :, jl))
          CALL adv_x(zusnit, pu_ice, 0._wp, zarea, z0snw(:, :, jl), sxsn(:, :, jl), sxxsn(:, :, jl), sysn(:, :, jl), syysn(:, :, jl), sxysn(:, :, jl))
          CALL adv_y(zusnit, pv_ice, 1._wp, zarea, z0smi(:, :, jl), sxsal(:, :, jl), sxxsal(:, :, jl), sysal(:, :, jl), syysal(:, :, jl), sxysal(:, :, jl))
          CALL adv_x(zusnit, pu_ice, 0._wp, zarea, z0smi(:, :, jl), sxsal(:, :, jl), sxxsal(:, :, jl), sysal(:, :, jl), syysal(:, :, jl), sxysal(:, :, jl))
          CALL adv_y(zusnit, pv_ice, 1._wp, zarea, z0oi(:, :, jl), sxage(:, :, jl), sxxage(:, :, jl), syage(:, :, jl), syyage(:, :, jl), sxyage(:, :, jl))
          CALL adv_x(zusnit, pu_ice, 0._wp, zarea, z0oi(:, :, jl), sxage(:, :, jl), sxxage(:, :, jl), syage(:, :, jl), syyage(:, :, jl), sxyage(:, :, jl))
          CALL adv_y(zusnit, pv_ice, 1._wp, zarea, z0ai(:, :, jl), sxa(:, :, jl), sxxa(:, :, jl), sya(:, :, jl), syya(:, :, jl), sxya(:, :, jl))
          CALL adv_x(zusnit, pu_ice, 0._wp, zarea, z0ai(:, :, jl), sxa(:, :, jl), sxxa(:, :, jl), sya(:, :, jl), syya(:, :, jl), sxya(:, :, jl))
          DO jk = 1, nlay_s
            CALL adv_y(zusnit, pv_ice, 1._wp, zarea, z0es(:, :, jk, jl), sxc0(:, :, jk, jl), sxxc0(:, :, jk, jl), syc0(:, :, jk, jl), syyc0(:, :, jk, jl), sxyc0(:, :, jk, jl))
            CALL adv_x(zusnit, pu_ice, 0._wp, zarea, z0es(:, :, jk, jl), sxc0(:, :, jk, jl), sxxc0(:, :, jk, jl), syc0(:, :, jk, jl), syyc0(:, :, jk, jl), sxyc0(:, :, jk, jl))
          END DO
          DO jk = 1, nlay_i
            CALL adv_y(zusnit, pv_ice, 1._wp, zarea, z0ei(:, :, jk, jl), sxe(:, :, jk, jl), sxxe(:, :, jk, jl), sye(:, :, jk, jl), syye(:, :, jk, jl), sxye(:, :, jk, jl))
            CALL adv_x(zusnit, pu_ice, 0._wp, zarea, z0ei(:, :, jk, jl), sxe(:, :, jk, jl), sxxe(:, :, jk, jl), sye(:, :, jk, jl), syye(:, :, jk, jl), sxye(:, :, jk, jl))
          END DO
          IF (ln_pnd_H12) THEN
            CALL adv_y(zusnit, pv_ice, 1._wp, zarea, z0ap(:, :, jl), sxap(:, :, jl), sxxap(:, :, jl), syap(:, :, jl), syyap(:, :, jl), sxyap(:, :, jl))
            CALL adv_x(zusnit, pu_ice, 0._wp, zarea, z0ap(:, :, jl), sxap(:, :, jl), sxxap(:, :, jl), syap(:, :, jl), syyap(:, :, jl), sxyap(:, :, jl))
            CALL adv_y(zusnit, pv_ice, 1._wp, zarea, z0vp(:, :, jl), sxvp(:, :, jl), sxxvp(:, :, jl), syvp(:, :, jl), syyvp(:, :, jl), sxyvp(:, :, jl))
            CALL adv_x(zusnit, pu_ice, 0._wp, zarea, z0vp(:, :, jl), sxvp(:, :, jl), sxxvp(:, :, jl), syvp(:, :, jl), syyvp(:, :, jl), sxyvp(:, :, jl))
          END IF
        END DO
      END DO
    END IF
    pato_i(:, :) = z0opw(:, :, 1) * r1_e1e2t(:, :) * tmask(:, :, 1)
    DO jl = 1, jpl
      pv_i(:, :, jl) = z0ice(:, :, jl) * r1_e1e2t(:, :) * tmask(:, :, 1)
      pv_s(:, :, jl) = z0snw(:, :, jl) * r1_e1e2t(:, :) * tmask(:, :, 1)
      psv_i(:, :, jl) = z0smi(:, :, jl) * r1_e1e2t(:, :) * tmask(:, :, 1)
      poa_i(:, :, jl) = z0oi(:, :, jl) * r1_e1e2t(:, :) * tmask(:, :, 1)
      pa_i(:, :, jl) = z0ai(:, :, jl) * r1_e1e2t(:, :) * tmask(:, :, 1)
      DO jk = 1, nlay_s
        pe_s(:, :, jk, jl) = z0es(:, :, jk, jl) * r1_e1e2t(:, :) * tmask(:, :, 1)
      END DO
      DO jk = 1, nlay_i
        pe_i(:, :, jk, jl) = z0ei(:, :, jk, jl) * r1_e1e2t(:, :) * tmask(:, :, 1)
      END DO
      IF (ln_pnd_H12) THEN
        pa_ip(:, :, jl) = z0ap(:, :, jl) * r1_e1e2t(:, :) * tmask(:, :, 1)
        pv_ip(:, :, jl) = z0vp(:, :, jl) * r1_e1e2t(:, :) * tmask(:, :, 1)
      END IF
    END DO
    DEALLOCATE(zarea, z0opw, z0ice, z0snw, z0ai, z0smi, z0oi, z0ap, z0vp, z0es, z0ei)
    IF (lrst_ice) CALL adv_pra_rst('WRITE', kt)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ice_dyn_adv_pra
  SUBROUTINE adv_x(pdf, put, pcrh, psm, ps0, psx, psxx, psy, psyy, psxy)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN) :: pdf
    REAL(KIND = wp), INTENT(IN) :: pcrh
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: put
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: psm
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: ps0
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: psx, psy
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: psxx, psyy, psxy
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zs1max, zrdt, zslpmax, ztemp
    REAL(KIND = wp) :: zs1new, zalf, zalfq, zbt
    REAL(KIND = wp) :: zs2new, zalf1, zalf1q, zbt1
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zf0, zfx, zfy, zbet
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zfm, zfxx, zfyy, zfxy
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zalg, zalg1, zalg1q
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('adv_x', 'r0', 0, 0)
    zrdt = rdt_ice * pdf
    DO jj = 1, jpj
      DO ji = 1, jpi
        zslpmax = MAX(0._wp, ps0(ji, jj))
        zs1max = 1.5 * zslpmax
        zs1new = MIN(zs1max, MAX(- zs1max, psx(ji, jj)))
        zs2new = MIN(2.0 * zslpmax - 0.3334 * ABS(zs1new), MAX(ABS(zs1new) - zslpmax, psxx(ji, jj)))
        rswitch = (1.0 - MAX(0._wp, SIGN(1._wp, - zslpmax))) * tmask(ji, jj, 1)
        ps0(ji, jj) = zslpmax
        psx(ji, jj) = zs1new * rswitch
        psxx(ji, jj) = zs2new * rswitch
        psy(ji, jj) = psy(ji, jj) * rswitch
        psyy(ji, jj) = psyy(ji, jj) * rswitch
        psxy(ji, jj) = MIN(zslpmax, MAX(- zslpmax, psxy(ji, jj))) * rswitch
      END DO
    END DO
    psm(:, :) = MAX(pcrh * e1e2t(:, :) + (1.0 - pcrh) * psm(:, :), epsi20)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zbet(ji, jj) = MAX(0._wp, SIGN(1._wp, put(ji, jj)))
        zalf = MAX(0._wp, put(ji, jj)) * zrdt * e2u(ji, jj) / psm(ji, jj)
        zalfq = zalf * zalf
        zalf1 = 1.0 - zalf
        zalf1q = zalf1 * zalf1
        zfm(ji, jj) = zalf * psm(ji, jj)
        zf0(ji, jj) = zalf * (ps0(ji, jj) + zalf1 * (psx(ji, jj) + (zalf1 - zalf) * psxx(ji, jj)))
        zfx(ji, jj) = zalfq * (psx(ji, jj) + 3.0 * zalf1 * psxx(ji, jj))
        zfxx(ji, jj) = zalf * psxx(ji, jj) * zalfq
        zfy(ji, jj) = zalf * (psy(ji, jj) + zalf1 * psxy(ji, jj))
        zfxy(ji, jj) = zalfq * psxy(ji, jj)
        zfyy(ji, jj) = zalf * psyy(ji, jj)
        psm(ji, jj) = psm(ji, jj) - zfm(ji, jj)
        ps0(ji, jj) = ps0(ji, jj) - zf0(ji, jj)
        psx(ji, jj) = zalf1q * (psx(ji, jj) - 3.0 * zalf * psxx(ji, jj))
        psxx(ji, jj) = zalf1 * zalf1q * psxx(ji, jj)
        psy(ji, jj) = psy(ji, jj) - zfy(ji, jj)
        psyy(ji, jj) = psyy(ji, jj) - zfyy(ji, jj)
        psxy(ji, jj) = zalf1q * psxy(ji, jj)
      END DO
    END DO
    DO jj = 1, jpjm1
      DO ji = 1, jpim1
        zalf = MAX(0._wp, - put(ji, jj)) * zrdt * e2u(ji, jj) / psm(ji + 1, jj)
        zalg(ji, jj) = zalf
        zalfq = zalf * zalf
        zalf1 = 1.0 - zalf
        zalg1(ji, jj) = zalf1
        zalf1q = zalf1 * zalf1
        zalg1q(ji, jj) = zalf1q
        zfm(ji, jj) = zfm(ji, jj) + zalf * psm(ji + 1, jj)
        zf0(ji, jj) = zf0(ji, jj) + zalf * (ps0(ji + 1, jj) - zalf1 * (psx(ji + 1, jj) - (zalf1 - zalf) * psxx(ji + 1, jj)))
        zfx(ji, jj) = zfx(ji, jj) + zalfq * (psx(ji + 1, jj) - 3.0 * zalf1 * psxx(ji + 1, jj))
        zfxx(ji, jj) = zfxx(ji, jj) + zalf * psxx(ji + 1, jj) * zalfq
        zfy(ji, jj) = zfy(ji, jj) + zalf * (psy(ji + 1, jj) - zalf1 * psxy(ji + 1, jj))
        zfxy(ji, jj) = zfxy(ji, jj) + zalfq * psxy(ji + 1, jj)
        zfyy(ji, jj) = zfyy(ji, jj) + zalf * psyy(ji + 1, jj)
      END DO
    END DO
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zbt = zbet(ji - 1, jj)
        zbt1 = 1.0 - zbet(ji - 1, jj)
        psm(ji, jj) = zbt * psm(ji, jj) + zbt1 * (psm(ji, jj) - zfm(ji - 1, jj))
        ps0(ji, jj) = zbt * ps0(ji, jj) + zbt1 * (ps0(ji, jj) - zf0(ji - 1, jj))
        psx(ji, jj) = zalg1q(ji - 1, jj) * (psx(ji, jj) + 3.0 * zalg(ji - 1, jj) * psxx(ji, jj))
        psxx(ji, jj) = zalg1(ji - 1, jj) * zalg1q(ji - 1, jj) * psxx(ji, jj)
        psy(ji, jj) = zbt * psy(ji, jj) + zbt1 * (psy(ji, jj) - zfy(ji - 1, jj))
        psyy(ji, jj) = zbt * psyy(ji, jj) + zbt1 * (psyy(ji, jj) - zfyy(ji - 1, jj))
        psxy(ji, jj) = zalg1q(ji - 1, jj) * psxy(ji, jj)
      END DO
    END DO
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zbt = zbet(ji - 1, jj)
        zbt1 = 1.0 - zbet(ji - 1, jj)
        psm(ji, jj) = zbt * (psm(ji, jj) + zfm(ji - 1, jj)) + zbt1 * psm(ji, jj)
        zalf = zbt * zfm(ji - 1, jj) / psm(ji, jj)
        zalf1 = 1.0 - zalf
        ztemp = zalf * ps0(ji, jj) - zalf1 * zf0(ji - 1, jj)
        ps0(ji, jj) = zbt * (ps0(ji, jj) + zf0(ji - 1, jj)) + zbt1 * ps0(ji, jj)
        psx(ji, jj) = zbt * (zalf * zfx(ji - 1, jj) + zalf1 * psx(ji, jj) + 3.0 * ztemp) + zbt1 * psx(ji, jj)
        psxx(ji, jj) = zbt * (zalf * zalf * zfxx(ji - 1, jj) + zalf1 * zalf1 * psxx(ji, jj) + 5.0 * (zalf * zalf1 * (psx(ji, jj) - zfx(ji - 1, jj)) - (zalf1 - zalf) * ztemp)) + zbt1 * psxx(ji, jj)
        psxy(ji, jj) = zbt * (zalf * zfxy(ji - 1, jj) + zalf1 * psxy(ji, jj) + 3.0 * (- zalf1 * zfy(ji - 1, jj) + zalf * psy(ji, jj))) + zbt1 * psxy(ji, jj)
        psy(ji, jj) = zbt * (psy(ji, jj) + zfy(ji - 1, jj)) + zbt1 * psy(ji, jj)
        psyy(ji, jj) = zbt * (psyy(ji, jj) + zfyy(ji - 1, jj)) + zbt1 * psyy(ji, jj)
      END DO
    END DO
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zbt = zbet(ji, jj)
        zbt1 = 1.0 - zbet(ji, jj)
        psm(ji, jj) = zbt * psm(ji, jj) + zbt1 * (psm(ji, jj) + zfm(ji, jj))
        zalf = zbt1 * zfm(ji, jj) / psm(ji, jj)
        zalf1 = 1.0 - zalf
        ztemp = - zalf * ps0(ji, jj) + zalf1 * zf0(ji, jj)
        ps0(ji, jj) = zbt * ps0(ji, jj) + zbt1 * (ps0(ji, jj) + zf0(ji, jj))
        psx(ji, jj) = zbt * psx(ji, jj) + zbt1 * (zalf * zfx(ji, jj) + zalf1 * psx(ji, jj) + 3.0 * ztemp)
        psxx(ji, jj) = zbt * psxx(ji, jj) + zbt1 * (zalf * zalf * zfxx(ji, jj) + zalf1 * zalf1 * psxx(ji, jj) + 5.0 * (zalf * zalf1 * (- psx(ji, jj) + zfx(ji, jj)) + (zalf1 - zalf) * ztemp))
        psxy(ji, jj) = zbt * psxy(ji, jj) + zbt1 * (zalf * zfxy(ji, jj) + zalf1 * psxy(ji, jj) + 3.0 * (zalf1 * zfy(ji, jj) - zalf * psy(ji, jj)))
        psy(ji, jj) = zbt * psy(ji, jj) + zbt1 * (psy(ji, jj) + zfy(ji, jj))
        psyy(ji, jj) = zbt * psyy(ji, jj) + zbt1 * (psyy(ji, jj) + zfyy(ji, jj))
      END DO
    END DO
    CALL lbc_lnk_multi('icedyn_adv_pra', psm, 'T', 1., ps0, 'T', 1., psx, 'T', - 1., psy, 'T', - 1., psxx, 'T', 1., psyy, 'T', 1., psxy, 'T', 1.)
    IF (ln_ctl) THEN
      CALL prt_ctl(tab2d_1 = psm, clinfo1 = ' adv_x: psm  :', tab2d_2 = ps0, clinfo2 = ' ps0  : ')
      CALL prt_ctl(tab2d_1 = psx, clinfo1 = ' adv_x: psx  :', tab2d_2 = psxx, clinfo2 = ' psxx : ')
      CALL prt_ctl(tab2d_1 = psy, clinfo1 = ' adv_x: psy  :', tab2d_2 = psyy, clinfo2 = ' psyy : ')
      CALL prt_ctl(tab2d_1 = psxy, clinfo1 = ' adv_x: psxy :')
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE adv_x
  SUBROUTINE adv_y(pdf, pvt, pcrh, psm, ps0, psx, psxx, psy, psyy, psxy)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN) :: pdf
    REAL(KIND = wp), INTENT(IN) :: pcrh
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pvt
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: psm
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: ps0
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: psx, psy
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: psxx, psyy, psxy
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zs1max, zrdt, zslpmax, ztemp
    REAL(KIND = wp) :: zs1new, zalf, zalfq, zbt
    REAL(KIND = wp) :: zs2new, zalf1, zalf1q, zbt1
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zf0, zfx, zfy, zbet
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zfm, zfxx, zfyy, zfxy
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zalg, zalg1, zalg1q
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('adv_y', 'r0', 0, 0)
    zrdt = rdt_ice * pdf
    DO jj = 1, jpj
      DO ji = 1, jpi
        zslpmax = MAX(0._wp, ps0(ji, jj))
        zs1max = 1.5 * zslpmax
        zs1new = MIN(zs1max, MAX(- zs1max, psy(ji, jj)))
        zs2new = MIN((2.0 * zslpmax - 0.3334 * ABS(zs1new)), MAX(ABS(zs1new) - zslpmax, psyy(ji, jj)))
        rswitch = (1.0 - MAX(0._wp, SIGN(1._wp, - zslpmax))) * tmask(ji, jj, 1)
        ps0(ji, jj) = zslpmax
        psx(ji, jj) = psx(ji, jj) * rswitch
        psxx(ji, jj) = psxx(ji, jj) * rswitch
        psy(ji, jj) = zs1new * rswitch
        psyy(ji, jj) = zs2new * rswitch
        psxy(ji, jj) = MIN(zslpmax, MAX(- zslpmax, psxy(ji, jj))) * rswitch
      END DO
    END DO
    psm(:, :) = MAX(pcrh * e1e2t(:, :) + (1.0 - pcrh) * psm(:, :), epsi20)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zbet(ji, jj) = MAX(0._wp, SIGN(1._wp, pvt(ji, jj)))
        zalf = MAX(0._wp, pvt(ji, jj)) * zrdt * e1v(ji, jj) / psm(ji, jj)
        zalfq = zalf * zalf
        zalf1 = 1.0 - zalf
        zalf1q = zalf1 * zalf1
        zfm(ji, jj) = zalf * psm(ji, jj)
        zf0(ji, jj) = zalf * (ps0(ji, jj) + zalf1 * (psy(ji, jj) + (zalf1 - zalf) * psyy(ji, jj)))
        zfy(ji, jj) = zalfq * (psy(ji, jj) + 3.0 * zalf1 * psyy(ji, jj))
        zfyy(ji, jj) = zalf * zalfq * psyy(ji, jj)
        zfx(ji, jj) = zalf * (psx(ji, jj) + zalf1 * psxy(ji, jj))
        zfxy(ji, jj) = zalfq * psxy(ji, jj)
        zfxx(ji, jj) = zalf * psxx(ji, jj)
        psm(ji, jj) = psm(ji, jj) - zfm(ji, jj)
        ps0(ji, jj) = ps0(ji, jj) - zf0(ji, jj)
        psy(ji, jj) = zalf1q * (psy(ji, jj) - 3.0 * zalf * psyy(ji, jj))
        psyy(ji, jj) = zalf1 * zalf1q * psyy(ji, jj)
        psx(ji, jj) = psx(ji, jj) - zfx(ji, jj)
        psxx(ji, jj) = psxx(ji, jj) - zfxx(ji, jj)
        psxy(ji, jj) = zalf1q * psxy(ji, jj)
      END DO
    END DO
    DO jj = 1, jpjm1
      DO ji = 1, jpi
        zalf = (MAX(0._wp, - pvt(ji, jj)) * zrdt * e1v(ji, jj)) / psm(ji, jj + 1)
        zalg(ji, jj) = zalf
        zalfq = zalf * zalf
        zalf1 = 1.0 - zalf
        zalg1(ji, jj) = zalf1
        zalf1q = zalf1 * zalf1
        zalg1q(ji, jj) = zalf1q
        zfm(ji, jj) = zfm(ji, jj) + zalf * psm(ji, jj + 1)
        zf0(ji, jj) = zf0(ji, jj) + zalf * (ps0(ji, jj + 1) - zalf1 * (psy(ji, jj + 1) - (zalf1 - zalf) * psyy(ji, jj + 1)))
        zfy(ji, jj) = zfy(ji, jj) + zalfq * (psy(ji, jj + 1) - 3.0 * zalf1 * psyy(ji, jj + 1))
        zfyy(ji, jj) = zfyy(ji, jj) + zalf * psyy(ji, jj + 1) * zalfq
        zfx(ji, jj) = zfx(ji, jj) + zalf * (psx(ji, jj + 1) - zalf1 * psxy(ji, jj + 1))
        zfxy(ji, jj) = zfxy(ji, jj) + zalfq * psxy(ji, jj + 1)
        zfxx(ji, jj) = zfxx(ji, jj) + zalf * psxx(ji, jj + 1)
      END DO
    END DO
    DO jj = 2, jpj
      DO ji = 1, jpi
        zbt = zbet(ji, jj - 1)
        zbt1 = (1.0 - zbet(ji, jj - 1))
        psm(ji, jj) = zbt * psm(ji, jj) + zbt1 * (psm(ji, jj) - zfm(ji, jj - 1))
        ps0(ji, jj) = zbt * ps0(ji, jj) + zbt1 * (ps0(ji, jj) - zf0(ji, jj - 1))
        psy(ji, jj) = zalg1q(ji, jj - 1) * (psy(ji, jj) + 3.0 * zalg(ji, jj - 1) * psyy(ji, jj))
        psyy(ji, jj) = zalg1(ji, jj - 1) * zalg1q(ji, jj - 1) * psyy(ji, jj)
        psx(ji, jj) = zbt * psx(ji, jj) + zbt1 * (psx(ji, jj) - zfx(ji, jj - 1))
        psxx(ji, jj) = zbt * psxx(ji, jj) + zbt1 * (psxx(ji, jj) - zfxx(ji, jj - 1))
        psxy(ji, jj) = zalg1q(ji, jj - 1) * psxy(ji, jj)
      END DO
    END DO
    DO jj = 2, jpjm1
      DO ji = 1, jpi
        zbt = zbet(ji, jj - 1)
        zbt1 = (1.0 - zbet(ji, jj - 1))
        psm(ji, jj) = zbt * (psm(ji, jj) + zfm(ji, jj - 1)) + zbt1 * psm(ji, jj)
        zalf = zbt * zfm(ji, jj - 1) / psm(ji, jj)
        zalf1 = 1.0 - zalf
        ztemp = zalf * ps0(ji, jj) - zalf1 * zf0(ji, jj - 1)
        ps0(ji, jj) = zbt * (ps0(ji, jj) + zf0(ji, jj - 1)) + zbt1 * ps0(ji, jj)
        psy(ji, jj) = zbt * (zalf * zfy(ji, jj - 1) + zalf1 * psy(ji, jj) + 3.0 * ztemp) + zbt1 * psy(ji, jj)
        psyy(ji, jj) = zbt * (zalf * zalf * zfyy(ji, jj - 1) + zalf1 * zalf1 * psyy(ji, jj) + 5.0 * (zalf * zalf1 * (psy(ji, jj) - zfy(ji, jj - 1)) - (zalf1 - zalf) * ztemp)) + zbt1 * psyy(ji, jj)
        psxy(ji, jj) = zbt * (zalf * zfxy(ji, jj - 1) + zalf1 * psxy(ji, jj) + 3.0 * (- zalf1 * zfx(ji, jj - 1) + zalf * psx(ji, jj))) + zbt1 * psxy(ji, jj)
        psx(ji, jj) = zbt * (psx(ji, jj) + zfx(ji, jj - 1)) + zbt1 * psx(ji, jj)
        psxx(ji, jj) = zbt * (psxx(ji, jj) + zfxx(ji, jj - 1)) + zbt1 * psxx(ji, jj)
      END DO
    END DO
    DO jj = 2, jpjm1
      DO ji = 1, jpi
        zbt = zbet(ji, jj)
        zbt1 = (1.0 - zbet(ji, jj))
        psm(ji, jj) = zbt * psm(ji, jj) + zbt1 * (psm(ji, jj) + zfm(ji, jj))
        zalf = zbt1 * zfm(ji, jj) / psm(ji, jj)
        zalf1 = 1.0 - zalf
        ztemp = - zalf * ps0(ji, jj) + zalf1 * zf0(ji, jj)
        ps0(ji, jj) = zbt * ps0(ji, jj) + zbt1 * (ps0(ji, jj) + zf0(ji, jj))
        psy(ji, jj) = zbt * psy(ji, jj) + zbt1 * (zalf * zfy(ji, jj) + zalf1 * psy(ji, jj) + 3.0 * ztemp)
        psyy(ji, jj) = zbt * psyy(ji, jj) + zbt1 * (zalf * zalf * zfyy(ji, jj) + zalf1 * zalf1 * psyy(ji, jj) + 5.0 * (zalf * zalf1 * (- psy(ji, jj) + zfy(ji, jj)) + (zalf1 - zalf) * ztemp))
        psxy(ji, jj) = zbt * psxy(ji, jj) + zbt1 * (zalf * zfxy(ji, jj) + zalf1 * psxy(ji, jj) + 3.0 * (zalf1 * zfx(ji, jj) - zalf * psx(ji, jj)))
        psx(ji, jj) = zbt * psx(ji, jj) + zbt1 * (psx(ji, jj) + zfx(ji, jj))
        psxx(ji, jj) = zbt * psxx(ji, jj) + zbt1 * (psxx(ji, jj) + zfxx(ji, jj))
      END DO
    END DO
    CALL lbc_lnk_multi('icedyn_adv_pra', psm, 'T', 1., ps0, 'T', 1., psx, 'T', - 1., psy, 'T', - 1., psxx, 'T', 1., psyy, 'T', 1., psxy, 'T', 1.)
    IF (ln_ctl) THEN
      CALL prt_ctl(tab2d_1 = psm, clinfo1 = ' adv_y: psm  :', tab2d_2 = ps0, clinfo2 = ' ps0  : ')
      CALL prt_ctl(tab2d_1 = psx, clinfo1 = ' adv_y: psx  :', tab2d_2 = psxx, clinfo2 = ' psxx : ')
      CALL prt_ctl(tab2d_1 = psy, clinfo1 = ' adv_y: psy  :', tab2d_2 = psyy, clinfo2 = ' psyy : ')
      CALL prt_ctl(tab2d_1 = psxy, clinfo1 = ' adv_y: psxy :')
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE adv_y
  SUBROUTINE adv_pra_init
    INTEGER :: ierr
    ALLOCATE(sxopw(jpi, jpj), syopw(jpi, jpj), sxxopw(jpi, jpj), syyopw(jpi, jpj), sxyopw(jpi, jpj), sxice(jpi, jpj, jpl), syice(jpi, jpj, jpl), sxxice(jpi, jpj, jpl), &
            syyice(jpi, jpj, jpl), sxyice(jpi, jpj, jpl), sxsn(jpi, jpj, jpl), sysn(jpi, jpj, jpl), sxxsn(jpi, jpj, jpl), syysn(jpi,jpj, jpl), sxysn(jpi, jpj, jpl), &
            sxa(jpi, jpj, jpl), sya(jpi, jpj, jpl), sxxa(jpi, jpj, jpl), syya(jpi, jpj, jpl), sxya(jpi, jpj, jpl), sxsal(jpi, jpj,jpl), sysal(jpi, jpj, jpl), &
    sxxsal(jpi, jpj, jpl), syysal(jpi, jpj, jpl), sxysal(jpi, jpj, jpl), sxage(jpi, jpj, jpl), syage(jpi, jpj, jpl), sxxage(jpi,jpj, jpl), syyage(jpi, jpj, jpl), &
        sxyage(jpi, jpj, jpl), sxap(jpi, jpj, jpl), syap(jpi, jpj, jpl), sxxap(jpi, jpj, jpl), syyap(jpi, jpj, jpl), sxyap(jpi, jpj,jpl), sxvp(jpi, jpj, jpl), &
        syvp(jpi, jpj, jpl), sxxvp(jpi, jpj, jpl), syyvp(jpi, jpj, jpl), sxyvp(jpi, jpj, jpl), sxc0(jpi, jpj, nlay_s, jpl),syc0(jpi, jpj, nlay_s, jpl), &
        sxxc0(jpi, jpj, nlay_s, jpl), syyc0(jpi, jpj, nlay_s, jpl), sxyc0(jpi, jpj, nlay_s, jpl), sxe(jpi, jpj, nlay_i, jpl),sye(jpi, jpj, nlay_i, jpl), &
        sxxe(jpi, jpj, nlay_i, jpl), syye(jpi, jpj, nlay_i, jpl), sxye(jpi, jpj, nlay_i, jpl), STAT = ierr)
    CALL mpp_sum('icedyn_adv_pra', ierr)
    IF (ierr /= 0) CALL ctl_stop('STOP', 'adv_pra_init : unable to allocate ice arrays for Prather advection scheme')
    CALL adv_pra_rst('READ')
  END SUBROUTINE adv_pra_init
  SUBROUTINE adv_pra_rst(cdrw, kt)
    CHARACTER(LEN = *), INTENT(IN) :: cdrw
    INTEGER, OPTIONAL, INTENT(IN) :: kt
    INTEGER :: jk, jl
    INTEGER :: iter
    INTEGER :: id1
    CHARACTER(LEN = 25) :: znam
    CHARACTER(LEN = 2) :: zchar, zchar1
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: z3d
    IF (TRIM(cdrw) == 'READ') THEN
      IF (ln_rstart) THEN
        id1 = iom_varid(numrir, 'sxopw', ldstop = .FALSE.)
      ELSE
        id1 = 0
      END IF
      IF (id1 > 0) THEN
        CALL iom_get(numrir, jpdom_autoglo, 'sxice', sxice)
        CALL iom_get(numrir, jpdom_autoglo, 'syice', syice)
        CALL iom_get(numrir, jpdom_autoglo, 'sxxice', sxxice)
        CALL iom_get(numrir, jpdom_autoglo, 'syyice', syyice)
        CALL iom_get(numrir, jpdom_autoglo, 'sxyice', sxyice)
        CALL iom_get(numrir, jpdom_autoglo, 'sxsn', sxsn)
        CALL iom_get(numrir, jpdom_autoglo, 'sysn', sysn)
        CALL iom_get(numrir, jpdom_autoglo, 'sxxsn', sxxsn)
        CALL iom_get(numrir, jpdom_autoglo, 'syysn', syysn)
        CALL iom_get(numrir, jpdom_autoglo, 'sxysn', sxysn)
        CALL iom_get(numrir, jpdom_autoglo, 'sxa', sxa)
        CALL iom_get(numrir, jpdom_autoglo, 'sya', sya)
        CALL iom_get(numrir, jpdom_autoglo, 'sxxa', sxxa)
        CALL iom_get(numrir, jpdom_autoglo, 'syya', syya)
        CALL iom_get(numrir, jpdom_autoglo, 'sxya', sxya)
        CALL iom_get(numrir, jpdom_autoglo, 'sxsal', sxsal)
        CALL iom_get(numrir, jpdom_autoglo, 'sysal', sysal)
        CALL iom_get(numrir, jpdom_autoglo, 'sxxsal', sxxsal)
        CALL iom_get(numrir, jpdom_autoglo, 'syysal', syysal)
        CALL iom_get(numrir, jpdom_autoglo, 'sxysal', sxysal)
        CALL iom_get(numrir, jpdom_autoglo, 'sxage', sxage)
        CALL iom_get(numrir, jpdom_autoglo, 'syage', syage)
        CALL iom_get(numrir, jpdom_autoglo, 'sxxage', sxxage)
        CALL iom_get(numrir, jpdom_autoglo, 'syyage', syyage)
        CALL iom_get(numrir, jpdom_autoglo, 'sxyage', sxyage)
        CALL iom_get(numrir, jpdom_autoglo, 'sxopw', sxopw)
        CALL iom_get(numrir, jpdom_autoglo, 'syopw', syopw)
        CALL iom_get(numrir, jpdom_autoglo, 'sxxopw', sxxopw)
        CALL iom_get(numrir, jpdom_autoglo, 'syyopw', syyopw)
        CALL iom_get(numrir, jpdom_autoglo, 'sxyopw', sxyopw)
        DO jk = 1, nlay_s
          WRITE(zchar1, '(I2.2)') jk
          znam = 'sxc0' // '_l' // zchar1
          CALL iom_get(numrir, jpdom_autoglo, znam, z3d)
          sxc0(:, :, jk, :) = z3d(:, :, :)
          znam = 'syc0' // '_l' // zchar1
          CALL iom_get(numrir, jpdom_autoglo, znam, z3d)
          syc0(:, :, jk, :) = z3d(:, :, :)
          znam = 'sxxc0' // '_l' // zchar1
          CALL iom_get(numrir, jpdom_autoglo, znam, z3d)
          sxxc0(:, :, jk, :) = z3d(:, :, :)
          znam = 'syyc0' // '_l' // zchar1
          CALL iom_get(numrir, jpdom_autoglo, znam, z3d)
          syyc0(:, :, jk, :) = z3d(:, :, :)
          znam = 'sxyc0' // '_l' // zchar1
          CALL iom_get(numrir, jpdom_autoglo, znam, z3d)
          sxyc0(:, :, jk, :) = z3d(:, :, :)
        END DO
        DO jk = 1, nlay_i
          WRITE(zchar1, '(I2.2)') jk
          znam = 'sxe' // '_l' // zchar1
          CALL iom_get(numrir, jpdom_autoglo, znam, z3d)
          sxe(:, :, jk, :) = z3d(:, :, :)
          znam = 'sye' // '_l' // zchar1
          CALL iom_get(numrir, jpdom_autoglo, znam, z3d)
          sye(:, :, jk, :) = z3d(:, :, :)
          znam = 'sxxe' // '_l' // zchar1
          CALL iom_get(numrir, jpdom_autoglo, znam, z3d)
          sxxe(:, :, jk, :) = z3d(:, :, :)
          znam = 'syye' // '_l' // zchar1
          CALL iom_get(numrir, jpdom_autoglo, znam, z3d)
          syye(:, :, jk, :) = z3d(:, :, :)
          znam = 'sxye' // '_l' // zchar1
          CALL iom_get(numrir, jpdom_autoglo, znam, z3d)
          sxye(:, :, jk, :) = z3d(:, :, :)
        END DO
        IF (ln_pnd_H12) THEN
          CALL iom_get(numrir, jpdom_autoglo, 'sxap', sxap)
          CALL iom_get(numrir, jpdom_autoglo, 'syap', syap)
          CALL iom_get(numrir, jpdom_autoglo, 'sxxap', sxxap)
          CALL iom_get(numrir, jpdom_autoglo, 'syyap', syyap)
          CALL iom_get(numrir, jpdom_autoglo, 'sxyap', sxyap)
          CALL iom_get(numrir, jpdom_autoglo, 'sxvp', sxvp)
          CALL iom_get(numrir, jpdom_autoglo, 'syvp', syvp)
          CALL iom_get(numrir, jpdom_autoglo, 'sxxvp', sxxvp)
          CALL iom_get(numrir, jpdom_autoglo, 'syyvp', syyvp)
          CALL iom_get(numrir, jpdom_autoglo, 'sxyvp', sxyvp)
        END IF
      ELSE
        IF (lwp) WRITE(numout, *) '   ==>>   start from rest OR previous run without Prather, set moments to 0'
        sxice = 0._wp
        syice = 0._wp
        sxxice = 0._wp
        syyice = 0._wp
        sxyice = 0._wp
        sxsn = 0._wp
        sysn = 0._wp
        sxxsn = 0._wp
        syysn = 0._wp
        sxysn = 0._wp
        sxa = 0._wp
        sya = 0._wp
        sxxa = 0._wp
        syya = 0._wp
        sxya = 0._wp
        sxsal = 0._wp
        sysal = 0._wp
        sxxsal = 0._wp
        syysal = 0._wp
        sxysal = 0._wp
        sxage = 0._wp
        syage = 0._wp
        sxxage = 0._wp
        syyage = 0._wp
        sxyage = 0._wp
        sxopw = 0._wp
        syopw = 0._wp
        sxxopw = 0._wp
        syyopw = 0._wp
        sxyopw = 0._wp
        sxc0 = 0._wp
        syc0 = 0._wp
        sxxc0 = 0._wp
        syyc0 = 0._wp
        sxyc0 = 0._wp
        sxe = 0._wp
        sye = 0._wp
        sxxe = 0._wp
        syye = 0._wp
        sxye = 0._wp
        IF (ln_pnd_H12) THEN
          sxap = 0._wp
          syap = 0._wp
          sxxap = 0._wp
          syyap = 0._wp
          sxyap = 0._wp
          sxvp = 0._wp
          syvp = 0._wp
          sxxvp = 0._wp
          syyvp = 0._wp
          sxyvp = 0._wp
        END IF
      END IF
    ELSE IF (TRIM(cdrw) == 'WRITE') THEN
      IF (lwp) WRITE(numout, *) '----  ice-adv-rst  ----'
      iter = kt + nn_fsbc - 1
      CALL iom_rstput(iter, nitrst, numriw, 'sxice', sxice)
      CALL iom_rstput(iter, nitrst, numriw, 'syice', syice)
      CALL iom_rstput(iter, nitrst, numriw, 'sxxice', sxxice)
      CALL iom_rstput(iter, nitrst, numriw, 'syyice', syyice)
      CALL iom_rstput(iter, nitrst, numriw, 'sxyice', sxyice)
      CALL iom_rstput(iter, nitrst, numriw, 'sxsn', sxsn)
      CALL iom_rstput(iter, nitrst, numriw, 'sysn', sysn)
      CALL iom_rstput(iter, nitrst, numriw, 'sxxsn', sxxsn)
      CALL iom_rstput(iter, nitrst, numriw, 'syysn', syysn)
      CALL iom_rstput(iter, nitrst, numriw, 'sxysn', sxysn)
      CALL iom_rstput(iter, nitrst, numriw, 'sxa', sxa)
      CALL iom_rstput(iter, nitrst, numriw, 'sya', sya)
      CALL iom_rstput(iter, nitrst, numriw, 'sxxa', sxxa)
      CALL iom_rstput(iter, nitrst, numriw, 'syya', syya)
      CALL iom_rstput(iter, nitrst, numriw, 'sxya', sxya)
      CALL iom_rstput(iter, nitrst, numriw, 'sxsal', sxsal)
      CALL iom_rstput(iter, nitrst, numriw, 'sysal', sysal)
      CALL iom_rstput(iter, nitrst, numriw, 'sxxsal', sxxsal)
      CALL iom_rstput(iter, nitrst, numriw, 'syysal', syysal)
      CALL iom_rstput(iter, nitrst, numriw, 'sxysal', sxysal)
      CALL iom_rstput(iter, nitrst, numriw, 'sxage', sxage)
      CALL iom_rstput(iter, nitrst, numriw, 'syage', syage)
      CALL iom_rstput(iter, nitrst, numriw, 'sxxage', sxxage)
      CALL iom_rstput(iter, nitrst, numriw, 'syyage', syyage)
      CALL iom_rstput(iter, nitrst, numriw, 'sxyage', sxyage)
      CALL iom_rstput(iter, nitrst, numriw, 'sxopw', sxopw)
      CALL iom_rstput(iter, nitrst, numriw, 'syopw', syopw)
      CALL iom_rstput(iter, nitrst, numriw, 'sxxopw', sxxopw)
      CALL iom_rstput(iter, nitrst, numriw, 'syyopw', syyopw)
      CALL iom_rstput(iter, nitrst, numriw, 'sxyopw', sxyopw)
      DO jk = 1, nlay_s
        WRITE(zchar1, '(I2.2)') jk
        znam = 'sxc0' // '_l' // zchar1
        z3d(:, :, :) = sxc0(:, :, jk, :)
        CALL iom_rstput(iter, nitrst, numriw, znam, z3d)
        znam = 'syc0' // '_l' // zchar1
        z3d(:, :, :) = syc0(:, :, jk, :)
        CALL iom_rstput(iter, nitrst, numriw, znam, z3d)
        znam = 'sxxc0' // '_l' // zchar1
        z3d(:, :, :) = sxxc0(:, :, jk, :)
        CALL iom_rstput(iter, nitrst, numriw, znam, z3d)
        znam = 'syyc0' // '_l' // zchar1
        z3d(:, :, :) = syyc0(:, :, jk, :)
        CALL iom_rstput(iter, nitrst, numriw, znam, z3d)
        znam = 'sxyc0' // '_l' // zchar1
        z3d(:, :, :) = sxyc0(:, :, jk, :)
        CALL iom_rstput(iter, nitrst, numriw, znam, z3d)
      END DO
      DO jk = 1, nlay_i
        WRITE(zchar1, '(I2.2)') jk
        znam = 'sxe' // '_l' // zchar1
        z3d(:, :, :) = sxe(:, :, jk, :)
        CALL iom_rstput(iter, nitrst, numriw, znam, z3d)
        znam = 'sye' // '_l' // zchar1
        z3d(:, :, :) = sye(:, :, jk, :)
        CALL iom_rstput(iter, nitrst, numriw, znam, z3d)
        znam = 'sxxe' // '_l' // zchar1
        z3d(:, :, :) = sxxe(:, :, jk, :)
        CALL iom_rstput(iter, nitrst, numriw, znam, z3d)
        znam = 'syye' // '_l' // zchar1
        z3d(:, :, :) = syye(:, :, jk, :)
        CALL iom_rstput(iter, nitrst, numriw, znam, z3d)
        znam = 'sxye' // '_l' // zchar1
        z3d(:, :, :) = sxye(:, :, jk, :)
        CALL iom_rstput(iter, nitrst, numriw, znam, z3d)
      END DO
      IF (ln_pnd_H12) THEN
        CALL iom_rstput(iter, nitrst, numriw, 'sxap', sxap)
        CALL iom_rstput(iter, nitrst, numriw, 'syap', syap)
        CALL iom_rstput(iter, nitrst, numriw, 'sxxap', sxxap)
        CALL iom_rstput(iter, nitrst, numriw, 'syyap', syyap)
        CALL iom_rstput(iter, nitrst, numriw, 'sxyap', sxyap)
        CALL iom_rstput(iter, nitrst, numriw, 'sxvp', sxvp)
        CALL iom_rstput(iter, nitrst, numriw, 'syvp', syvp)
        CALL iom_rstput(iter, nitrst, numriw, 'sxxvp', sxxvp)
        CALL iom_rstput(iter, nitrst, numriw, 'syyvp', syyvp)
        CALL iom_rstput(iter, nitrst, numriw, 'sxyvp', sxyvp)
      END IF
    END IF
  END SUBROUTINE adv_pra_rst
END MODULE icedyn_adv_pra
