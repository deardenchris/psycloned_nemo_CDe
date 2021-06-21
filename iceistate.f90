MODULE iceistate
  USE phycst
  USE oce
  USE dom_oce
  USE sbc_oce, ONLY: sst_m, sss_m, ln_ice_embd
  USE sbc_ice, ONLY: tn_ice, snwice_mass, snwice_mass_b
  USE eosbn2
  USE domvvl
  USE ice
  USE icevar
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  USE fldread
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_istate
  PUBLIC :: ice_istate_init
  INTEGER, PARAMETER :: jpfldi = 6
  INTEGER, PARAMETER :: jp_hti = 1
  INTEGER, PARAMETER :: jp_hts = 2
  INTEGER, PARAMETER :: jp_ati = 3
  INTEGER, PARAMETER :: jp_tsu = 4
  INTEGER, PARAMETER :: jp_tmi = 5
  INTEGER, PARAMETER :: jp_smi = 6
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: si
  LOGICAL :: ln_iceini
  LOGICAL :: ln_iceini_file
  REAL(KIND = wp) :: rn_thres_sst
  REAL(KIND = wp) :: rn_hts_ini_n
  REAL(KIND = wp) :: rn_hts_ini_s
  REAL(KIND = wp) :: rn_hti_ini_n
  REAL(KIND = wp) :: rn_hti_ini_s
  REAL(KIND = wp) :: rn_ati_ini_n
  REAL(KIND = wp) :: rn_ati_ini_s
  REAL(KIND = wp) :: rn_smi_ini_n
  REAL(KIND = wp) :: rn_smi_ini_s
  REAL(KIND = wp) :: rn_tmi_ini_n
  REAL(KIND = wp) :: rn_tmi_ini_s
  CONTAINS
  SUBROUTINE ice_istate
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jj, jk, jl
    INTEGER :: i_hemis, i_fill, jl0
    REAL(KIND = wp) :: ztmelts, zdh
    REAL(KIND = wp) :: zarg, zV, zconv, zdv, zfac
    INTEGER, DIMENSION(4) :: itest
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z2d
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zswitch
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zht_i_ini, zat_i_ini, zvt_i_ini
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zts_u_ini, zht_s_ini, zsm_i_ini, ztm_i_ini
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: zh_i_ini, za_i_ini
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ice_istate', 'r0', 0, 0)
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) 'ice_istate: sea-ice initialization '
    IF (lwp) WRITE(numout, *) '~~~~~~~~~~'
    DO jl = 1, jpl
      t_su(:, :, jl) = rt0 * tmask(:, :, 1)
      cnd_ice(:, :, jl) = 0._wp
    END DO
    CALL eos_fzp(sss_m(:, :), t_bo(:, :))
    t_bo(:, :) = (t_bo(:, :) + rt0) * tmask(:, :, 1)
    IF (ln_iceini) THEN
      IF (ln_iceini_file) THEN
        zht_i_ini(:, :) = si(jp_hti) % fnow(:, :, 1)
        zht_s_ini(:, :) = si(jp_hts) % fnow(:, :, 1)
        zat_i_ini(:, :) = si(jp_ati) % fnow(:, :, 1)
        zts_u_ini(:, :) = si(jp_tsu) % fnow(:, :, 1)
        ztm_i_ini(:, :) = si(jp_tmi) % fnow(:, :, 1)
        zsm_i_ini(:, :) = si(jp_smi) % fnow(:, :, 1)
        WHERE (zat_i_ini(:, :) > 0._wp)
          zswitch(:, :) = tmask(:, :, 1)
        ELSEWHERE
          zswitch(:, :) = 0._wp
        END WHERE
        zvt_i_ini(:, :) = zht_i_ini(:, :) * zat_i_ini(:, :)
      ELSE
        WHERE ((sst_m(:, :) - (t_bo(:, :) - rt0)) * tmask(:, :, 1) >= rn_thres_sst)
          zswitch(:, :) = 0._wp
        ELSEWHERE
          zswitch(:, :) = tmask(:, :, 1)
        END WHERE
        WHERE (ff_t(:, :) >= 0._wp)
          zht_i_ini(:, :) = rn_hti_ini_n * zswitch(:, :)
          zht_s_ini(:, :) = rn_hts_ini_n * zswitch(:, :)
          zat_i_ini(:, :) = rn_ati_ini_n * zswitch(:, :)
          zts_u_ini(:, :) = rn_tmi_ini_n * zswitch(:, :)
          zsm_i_ini(:, :) = rn_smi_ini_n * zswitch(:, :)
          ztm_i_ini(:, :) = rn_tmi_ini_n * zswitch(:, :)
        ELSEWHERE
          zht_i_ini(:, :) = rn_hti_ini_s * zswitch(:, :)
          zht_s_ini(:, :) = rn_hts_ini_s * zswitch(:, :)
          zat_i_ini(:, :) = rn_ati_ini_s * zswitch(:, :)
          zts_u_ini(:, :) = rn_tmi_ini_s * zswitch(:, :)
          zsm_i_ini(:, :) = rn_smi_ini_s * zswitch(:, :)
          ztm_i_ini(:, :) = rn_tmi_ini_s * zswitch(:, :)
        END WHERE
        zvt_i_ini(:, :) = zht_i_ini(:, :) * zat_i_ini(:, :)
      END IF
      IF (jpl == 1) THEN
        zh_i_ini(:, :, 1) = zht_i_ini(:, :)
        za_i_ini(:, :, 1) = zat_i_ini(:, :)
      ELSE
        zh_i_ini(:, :, :) = 0._wp
        za_i_ini(:, :, :) = 0._wp
        DO jj = 1, jpj
          DO ji = 1, jpi
            IF (zat_i_ini(ji, jj) > 0._wp .AND. zht_i_ini(ji, jj) > 0._wp) THEN
              jl0 = jpl
              DO jl = 1, jpl
                IF ((zht_i_ini(ji, jj) > hi_max(jl - 1)) .AND. (zht_i_ini(ji, jj) <= hi_max(jl))) THEN
                  jl0 = jl
                  CYCLE
                END IF
              END DO
              itest(:) = 0
              i_fill = jpl + 1
              DO WHILE ((SUM(itest(:)) /= 4) .AND. (i_fill >= 2))
                i_fill = i_fill - 1
                zh_i_ini(ji, jj, :) = 0._wp
                za_i_ini(ji, jj, :) = 0._wp
                itest(:) = 0
                IF (i_fill == 1) THEN
                  zh_i_ini(ji, jj, 1) = zht_i_ini(ji, jj)
                  za_i_ini(ji, jj, 1) = zat_i_ini(ji, jj)
                ELSE
                  DO jl = 1, i_fill - 1
                    zh_i_ini(ji, jj, jl) = hi_mean(jl)
                  END DO
                  za_i_ini(ji, jj, jl0) = zat_i_ini(ji, jj) / SQRT(REAL(jpl))
                  DO jl = 1, i_fill - 1
                    IF (jl /= jl0) THEN
                      zarg = (zh_i_ini(ji, jj, jl) - zht_i_ini(ji, jj)) / (0.5_wp * zht_i_ini(ji, jj))
                      za_i_ini(ji, jj, jl) = za_i_ini(ji, jj, jl0) * EXP(- zarg ** 2)
                    END IF
                  END DO
                  za_i_ini(ji, jj, i_fill) = zat_i_ini(ji, jj) - SUM(za_i_ini(ji, jj, 1 : i_fill - 1))
                  zV = SUM(za_i_ini(ji, jj, 1 : i_fill - 1) * zh_i_ini(ji, jj, 1 : i_fill - 1))
                  zh_i_ini(ji, jj, i_fill) = (zvt_i_ini(ji, jj) - zV) / MAX(za_i_ini(ji, jj, i_fill), epsi10)
                  IF (jl0 /= jpl) THEN
                    DO jl = jpl, jl0 + 1, - 1
                      IF (za_i_ini(ji, jj, jl) > za_i_ini(ji, jj, jl - 1)) THEN
                        zdv = zh_i_ini(ji, jj, jl) * za_i_ini(ji, jj, jl)
                        zh_i_ini(ji, jj, jl) = 0._wp
                        za_i_ini(ji, jj, jl) = 0._wp
                        za_i_ini(ji, jj, 1 : jl - 1) = za_i_ini(ji, jj, 1 : jl - 1) + zdv / MAX(REAL(jl - 1) * zht_i_ini(ji, jj), epsi10)
                      END IF
                    END DO
                  END IF
                END IF
                zconv = ABS(zat_i_ini(ji, jj) - SUM(za_i_ini(ji, jj, 1 : jpl)))
                IF (zconv < epsi06) itest(1) = 1
                zconv = ABS(zat_i_ini(ji, jj) * zht_i_ini(ji, jj) - SUM(za_i_ini(ji, jj, 1 : jpl) * zh_i_ini(ji, jj, 1 : jpl)))
                IF (zconv < epsi06) itest(2) = 1
                IF (zh_i_ini(ji, jj, i_fill) >= hi_max(i_fill - 1)) itest(3) = 1
                itest(4) = 1
                DO jl = 1, i_fill
                  IF (za_i_ini(ji, jj, jl) < 0._wp) itest(4) = 0
                END DO
              END DO
              IF (lwp .AND. SUM(itest) /= 4) THEN
                WRITE(numout, *)
                WRITE(numout, *) ' !!!! ALERT itest is not equal to 4      !!! '
                WRITE(numout, *) ' !!!! Something is wrong in the SI3 initialization procedure '
                WRITE(numout, *)
                WRITE(numout, *) ' *** itest_i (i=1,4) = ', itest(:)
                WRITE(numout, *) ' zat_i_ini : ', zat_i_ini(ji, jj)
                WRITE(numout, *) ' zht_i_ini : ', zht_i_ini(ji, jj)
              END IF
            END IF
          END DO
        END DO
      END IF
      DO jl = 1, jpl
        DO jj = 1, jpj
          DO ji = 1, jpi
            a_i(ji, jj, jl) = zswitch(ji, jj) * za_i_ini(ji, jj, jl)
            h_i(ji, jj, jl) = zswitch(ji, jj) * zh_i_ini(ji, jj, jl)
            s_i(ji, jj, jl) = zswitch(ji, jj) * zsm_i_ini(ji, jj)
            o_i(ji, jj, jl) = 0._wp
            t_su(ji, jj, jl) = zswitch(ji, jj) * zts_u_ini(ji, jj) + (1._wp - zswitch(ji, jj)) * rt0
            IF (zht_i_ini(ji, jj) > 0._wp) THEN
              h_s(ji, jj, jl) = h_i(ji, jj, jl) * (zht_s_ini(ji, jj) / zht_i_ini(ji, jj))
            ELSE
              h_s(ji, jj, jl) = 0._wp
            END IF
            zdh = MAX(0._wp, (rhos * h_s(ji, jj, jl) + (rhoi - rau0) * h_i(ji, jj, jl)) * r1_rau0)
            h_i(ji, jj, jl) = MIN(hi_max(jl), h_i(ji, jj, jl) + zdh)
            h_s(ji, jj, jl) = MAX(0._wp, h_s(ji, jj, jl) - zdh * rhoi * r1_rhos)
            v_i(ji, jj, jl) = h_i(ji, jj, jl) * a_i(ji, jj, jl)
            v_s(ji, jj, jl) = h_s(ji, jj, jl) * a_i(ji, jj, jl)
            sv_i(ji, jj, jl) = MIN(s_i(ji, jj, jl), sss_m(ji, jj)) * v_i(ji, jj, jl)
            oa_i(ji, jj, jl) = o_i(ji, jj, jl) * a_i(ji, jj, jl)
          END DO
        END DO
      END DO
      IF (nn_icesal /= 2) THEN
        CALL ice_var_salprof
        sv_i = s_i * v_i
      END IF
      DO jk = 1, nlay_s
        DO jl = 1, jpl
          DO jj = 1, jpj
            DO ji = 1, jpi
              t_s(ji, jj, jk, jl) = zswitch(ji, jj) * ztm_i_ini(ji, jj) + (1._wp - zswitch(ji, jj)) * rt0
              e_s(ji, jj, jk, jl) = zswitch(ji, jj) * rhos * (rcpi * (rt0 - t_s(ji, jj, jk, jl)) + rLfus)
              e_s(ji, jj, jk, jl) = e_s(ji, jj, jk, jl) * v_s(ji, jj, jl) * r1_nlay_s
            END DO
          END DO
        END DO
      END DO
      DO jk = 1, nlay_i
        DO jl = 1, jpl
          DO jj = 1, jpj
            DO ji = 1, jpi
              t_i(ji, jj, jk, jl) = zswitch(ji, jj) * ztm_i_ini(ji, jj) + (1._wp - zswitch(ji, jj)) * rt0
              sz_i(ji, jj, jk, jl) = zswitch(ji, jj) * zsm_i_ini(ji, jj) + (1._wp - zswitch(ji, jj)) * rn_simin
              ztmelts = - rTmlt * sz_i(ji, jj, jk, jl) + rt0
              e_i(ji, jj, jk, jl) = zswitch(ji, jj) * rhoi * (rcpi * (ztmelts - t_i(ji, jj, jk, jl)) + rLfus * (1._wp - (ztmelts - rt0) / MIN((t_i(ji, jj, jk, jl) - rt0), - epsi20)) - rcp * (ztmelts - rt0))
              e_i(ji, jj, jk, jl) = e_i(ji, jj, jk, jl) * v_i(ji, jj, jl) * r1_nlay_i
            END DO
          END DO
        END DO
      END DO
      tn_ice(:, :, :) = t_su(:, :, :)
      t1_ice(:, :, :) = t_i(:, :, 1, :)
      IF (ln_pnd_cst .OR. ln_pnd_h12) THEN
        zfac = 1._wp
      ELSE
        zfac = 0._wp
      END IF
      DO jl = 1, jpl
        a_ip_frac(:, :, jl) = rn_apnd * zswitch(:, :) * zfac
        h_ip(:, :, jl) = rn_hpnd * zswitch(:, :) * zfac
      END DO
      a_ip(:, :, :) = a_ip_frac(:, :, :) * a_i(:, :, :)
      v_ip(:, :, :) = h_ip(:, :, :) * a_ip(:, :, :)
    ELSE
      a_i(:, :, :) = 0._wp
      v_i(:, :, :) = 0._wp
      v_s(:, :, :) = 0._wp
      sv_i(:, :, :) = 0._wp
      oa_i(:, :, :) = 0._wp
      h_i(:, :, :) = 0._wp
      h_s(:, :, :) = 0._wp
      s_i(:, :, :) = 0._wp
      o_i(:, :, :) = 0._wp
      e_i(:, :, :, :) = 0._wp
      e_s(:, :, :, :) = 0._wp
      DO jl = 1, jpl
        DO jk = 1, nlay_i
          t_i(:, :, jk, jl) = rt0 * tmask(:, :, 1)
        END DO
        DO jk = 1, nlay_s
          t_s(:, :, jk, jl) = rt0 * tmask(:, :, 1)
        END DO
      END DO
      tn_ice(:, :, :) = t_i(:, :, 1, :)
      t1_ice(:, :, :) = t_i(:, :, 1, :)
      a_ip(:, :, :) = 0._wp
      v_ip(:, :, :) = 0._wp
      a_ip_frac(:, :, :) = 0._wp
      h_ip(:, :, :) = 0._wp
    END IF
    at_i(:, :) = 0.0_wp
    DO jl = 1, jpl
      at_i(:, :) = at_i(:, :) + a_i(:, :, jl)
    END DO
    u_ice(:, :) = 0._wp
    v_ice(:, :) = 0._wp
    l_split_advumx(1) = .FALSE.
    snwice_mass(:, :) = tmask(:, :, 1) * SUM(rhos * v_s(:, :, :) + rhoi * v_i(:, :, :), dim = 3)
    snwice_mass_b(:, :) = snwice_mass(:, :)
    IF (ln_ice_embd) THEN
      sshn(:, :) = sshn(:, :) - snwice_mass(:, :) * r1_rau0
      sshb(:, :) = sshb(:, :) - snwice_mass(:, :) * r1_rau0
      IF (.NOT. ln_linssh) THEN
        WHERE (ht_0(:, :) > 0)
          z2d(:, :) = 1._wp + sshn(:, :) * tmask(:, :, 1) / ht_0(:, :)
        ELSEWHERE
          z2d(:, :) = 1._wp
        END WHERE
        DO jk = 1, jpkm1
          e3t_n(:, :, jk) = e3t_0(:, :, jk) * z2d(:, :)
          e3t_b(:, :, jk) = e3t_n(:, :, jk)
          e3t_a(:, :, jk) = e3t_n(:, :, jk)
        END DO
        CALL dom_vvl_interpol(e3t_b(:, :, :), e3u_b(:, :, :), 'U')
        CALL dom_vvl_interpol(e3t_b(:, :, :), e3v_b(:, :, :), 'V')
        CALL dom_vvl_interpol(e3t_n(:, :, :), e3u_n(:, :, :), 'U')
        CALL dom_vvl_interpol(e3t_n(:, :, :), e3v_n(:, :, :), 'V')
        CALL dom_vvl_interpol(e3u_n(:, :, :), e3f_n(:, :, :), 'F')
        CALL dom_vvl_interpol(e3t_n(:, :, :), e3w_n(:, :, :), 'W')
        CALL dom_vvl_interpol(e3u_n(:, :, :), e3uw_n(:, :, :), 'UW')
        CALL dom_vvl_interpol(e3v_n(:, :, :), e3vw_n(:, :, :), 'VW')
        CALL dom_vvl_interpol(e3u_b(:, :, :), e3uw_b(:, :, :), 'UW')
        CALL dom_vvl_interpol(e3v_b(:, :, :), e3vw_b(:, :, :), 'VW')
        gdept_n(:, :, 1) = 0.5_wp * e3w_n(:, :, 1)
        gdepw_n(:, :, 1) = 0.0_wp
        gde3w_n(:, :, 1) = gdept_n(:, :, 1) - sshn(:, :)
        DO jk = 2, jpk
          gdept_n(:, :, jk) = gdept_n(:, :, jk - 1) + e3w_n(:, :, jk)
          gdepw_n(:, :, jk) = gdepw_n(:, :, jk - 1) + e3t_n(:, :, jk - 1)
          gde3w_n(:, :, jk) = gdept_n(:, :, jk) - sshn(:, :)
        END DO
      END IF
    END IF
    a_i_b(:, :, :) = a_i(:, :, :)
    e_i_b(:, :, :, :) = e_i(:, :, :, :)
    v_i_b(:, :, :) = v_i(:, :, :)
    v_s_b(:, :, :) = v_s(:, :, :)
    e_s_b(:, :, :, :) = e_s(:, :, :, :)
    sv_i_b(:, :, :) = sv_i(:, :, :)
    oa_i_b(:, :, :) = oa_i(:, :, :)
    u_ice_b(:, :) = u_ice(:, :)
    v_ice_b(:, :) = v_ice(:, :)
    at_i_b(:, :) = at_i(:, :)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ice_istate
  SUBROUTINE ice_istate_init
    INTEGER :: ji, jj
    INTEGER :: ios, ierr, inum_ice
    INTEGER :: ifpr, ierror
    CHARACTER(LEN = 256) :: cn_dir
    TYPE(FLD_N) :: sn_hti, sn_hts, sn_ati, sn_tsu, sn_tmi, sn_smi
    TYPE(FLD_N), DIMENSION(jpfldi) :: slf_i
    NAMELIST /namini/ ln_iceini, ln_iceini_file, rn_thres_sst, rn_hts_ini_n, rn_hts_ini_s, rn_hti_ini_n, rn_hti_ini_s, rn_ati_ini_n, rn_ati_ini_s, rn_smi_ini_n, rn_smi_ini_s, rn_tmi_ini_n, rn_tmi_ini_s, sn_hti, sn_hts, sn_ati, sn_tsu, sn_tmi, sn_smi, cn_dir
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, namini, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namini in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, namini, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namini in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, namini)
    slf_i(jp_hti) = sn_hti
    slf_i(jp_hts) = sn_hts
    slf_i(jp_ati) = sn_ati
    slf_i(jp_tsu) = sn_tsu
    slf_i(jp_tmi) = sn_tmi
    slf_i(jp_smi) = sn_smi
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ice_istate_init: ice parameters inititialisation '
      WRITE(numout, *) '~~~~~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist namini:'
      WRITE(numout, *) '      initialization with ice (T) or not (F)                 ln_iceini       = ', ln_iceini
      WRITE(numout, *) '      ice initialization from a netcdf file                  ln_iceini_file  = ', ln_iceini_file
      WRITE(numout, *) '      max delta ocean temp. above Tfreeze with initial ice   rn_thres_sst    = ', rn_thres_sst
      WRITE(numout, *) '      initial snow thickness in the north                    rn_hts_ini_n    = ', rn_hts_ini_n
      WRITE(numout, *) '      initial snow thickness in the south                    rn_hts_ini_s    = ', rn_hts_ini_s
      WRITE(numout, *) '      initial ice thickness  in the north                    rn_hti_ini_n    = ', rn_hti_ini_n
      WRITE(numout, *) '      initial ice thickness  in the south                    rn_hti_ini_s    = ', rn_hti_ini_s
      WRITE(numout, *) '      initial ice concentr.  in the north                    rn_ati_ini_n    = ', rn_ati_ini_n
      WRITE(numout, *) '      initial ice concentr.  in the north                    rn_ati_ini_s    = ', rn_ati_ini_s
      WRITE(numout, *) '      initial  ice salinity  in the north                    rn_smi_ini_n    = ', rn_smi_ini_n
      WRITE(numout, *) '      initial  ice salinity  in the south                    rn_smi_ini_s    = ', rn_smi_ini_s
      WRITE(numout, *) '      initial  ice/snw temp  in the north                    rn_tmi_ini_n    = ', rn_tmi_ini_n
      WRITE(numout, *) '      initial  ice/snw temp  in the south                    rn_tmi_ini_s    = ', rn_tmi_ini_s
    END IF
    IF (ln_iceini_file) THEN
      ALLOCATE(si(jpfldi), STAT = ierror)
      IF (ierror > 0) THEN
        CALL ctl_stop('Ice_ini in iceistate: unable to allocate si structure')
        RETURN
      END IF
      DO ifpr = 1, jpfldi
        ALLOCATE(si(ifpr) % fnow(jpi, jpj, 1))
        ALLOCATE(si(ifpr) % fdta(jpi, jpj, 1, 2))
      END DO
      CALL fld_fill(si, slf_i, cn_dir, 'ice_istate', 'ice istate ini', 'numnam_ice')
      CALL fld_read(nit000, 1, si)
    END IF
  END SUBROUTINE ice_istate_init
END MODULE iceistate