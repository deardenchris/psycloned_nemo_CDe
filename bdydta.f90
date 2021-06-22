MODULE bdydta
  USE oce
  USE dom_oce
  USE phycst
  USE sbcapr
  USE sbctide
  USE bdy_oce
  USE bdytides
  USE ice
  USE icevar
  USE lib_mpp, ONLY: ctl_stop, ctl_nam
  USE fldread
  USE iom
  USE in_out_manager
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: bdy_dta
  PUBLIC :: bdy_dta_init
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nb_bdy_fld
  INTEGER :: nb_bdy_fld_sum
  LOGICAL, DIMENSION(jp_bdy) :: ln_full_vel_array
  TYPE(FLD), PUBLIC, ALLOCATABLE, DIMENSION(:), TARGET :: bf
  TYPE(MAP_POINTER), ALLOCATABLE, DIMENSION(:) :: nbmap_ptr
  INTEGER :: nice_cat
  INTEGER :: jfld_hti, jfld_hts, jfld_ai
  INTEGER, DIMENSION(jp_bdy) :: jfld_htit, jfld_htst, jfld_ait
  CONTAINS
  SUBROUTINE bdy_dta(kt, jit, time_offset)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN), OPTIONAL :: jit
    INTEGER, INTENT(IN), OPTIONAL :: time_offset
    INTEGER :: jbdy, jfld, jstart, jend, ib, jl
    INTEGER :: ii, ij, ik, igrd
    INTEGER, DIMENSION(jpbgrd) :: ilen1
    INTEGER, POINTER, DIMENSION(:) :: nblen, nblenrim
    TYPE(OBC_DATA), POINTER :: dta
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('bdy_dta', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('bdy_dta')
    IF (kt == nit000 .AND. .NOT. PRESENT(jit)) THEN
      DO jbdy = 1, nb_bdy
        nblen => idx_bdy(jbdy) % nblen
        nblenrim => idx_bdy(jbdy) % nblenrim
        dta => dta_bdy(jbdy)
        IF (nn_dyn2d_dta(jbdy) == 0) THEN
          ilen1(:) = nblen(:)
          IF (dta % ll_ssh) THEN
            igrd = 1
            DO ib = 1, ilen1(igrd)
              ii = idx_bdy(jbdy) % nbi(ib, igrd)
              ij = idx_bdy(jbdy) % nbj(ib, igrd)
              dta_bdy(jbdy) % ssh(ib) = sshn(ii, ij) * tmask(ii, ij, 1)
            END DO
          END IF
          IF (dta % ll_u2d) THEN
            igrd = 2
            DO ib = 1, ilen1(igrd)
              ii = idx_bdy(jbdy) % nbi(ib, igrd)
              ij = idx_bdy(jbdy) % nbj(ib, igrd)
              dta_bdy(jbdy) % u2d(ib) = un_b(ii, ij) * umask(ii, ij, 1)
            END DO
          END IF
          IF (dta % ll_v2d) THEN
            igrd = 3
            DO ib = 1, ilen1(igrd)
              ii = idx_bdy(jbdy) % nbi(ib, igrd)
              ij = idx_bdy(jbdy) % nbj(ib, igrd)
              dta_bdy(jbdy) % v2d(ib) = vn_b(ii, ij) * vmask(ii, ij, 1)
            END DO
          END IF
        END IF
        IF (nn_dyn3d_dta(jbdy) == 0) THEN
          ilen1(:) = nblen(:)
          IF (dta % ll_u3d) THEN
            igrd = 2
            DO ib = 1, ilen1(igrd)
              DO ik = 1, jpkm1
                ii = idx_bdy(jbdy) % nbi(ib, igrd)
                ij = idx_bdy(jbdy) % nbj(ib, igrd)
                dta_bdy(jbdy) % u3d(ib, ik) = (un(ii, ij, ik) - un_b(ii, ij)) * umask(ii, ij, ik)
              END DO
            END DO
          END IF
          IF (dta % ll_v3d) THEN
            igrd = 3
            DO ib = 1, ilen1(igrd)
              DO ik = 1, jpkm1
                ii = idx_bdy(jbdy) % nbi(ib, igrd)
                ij = idx_bdy(jbdy) % nbj(ib, igrd)
                dta_bdy(jbdy) % v3d(ib, ik) = (vn(ii, ij, ik) - vn_b(ii, ij)) * vmask(ii, ij, ik)
              END DO
            END DO
          END IF
        END IF
        IF (nn_tra_dta(jbdy) == 0) THEN
          ilen1(:) = nblen(:)
          IF (dta % ll_tem) THEN
            igrd = 1
            DO ib = 1, ilen1(igrd)
              DO ik = 1, jpkm1
                ii = idx_bdy(jbdy) % nbi(ib, igrd)
                ij = idx_bdy(jbdy) % nbj(ib, igrd)
                dta_bdy(jbdy) % tem(ib, ik) = tsn(ii, ij, ik, jp_tem) * tmask(ii, ij, ik)
              END DO
            END DO
          END IF
          IF (dta % ll_sal) THEN
            igrd = 1
            DO ib = 1, ilen1(igrd)
              DO ik = 1, jpkm1
                ii = idx_bdy(jbdy) % nbi(ib, igrd)
                ij = idx_bdy(jbdy) % nbj(ib, igrd)
                dta_bdy(jbdy) % sal(ib, ik) = tsn(ii, ij, ik, jp_sal) * tmask(ii, ij, ik)
              END DO
            END DO
          END IF
        END IF
        IF (nn_ice_dta(jbdy) == 0) THEN
          ilen1(:) = nblen(:)
          IF (dta % ll_a_i) THEN
            igrd = 1
            DO jl = 1, jpl
              DO ib = 1, ilen1(igrd)
                ii = idx_bdy(jbdy) % nbi(ib, igrd)
                ij = idx_bdy(jbdy) % nbj(ib, igrd)
                dta_bdy(jbdy) % a_i(ib, jl) = a_i(ii, ij, jl) * tmask(ii, ij, 1)
              END DO
            END DO
          END IF
          IF (dta % ll_h_i) THEN
            igrd = 1
            DO jl = 1, jpl
              DO ib = 1, ilen1(igrd)
                ii = idx_bdy(jbdy) % nbi(ib, igrd)
                ij = idx_bdy(jbdy) % nbj(ib, igrd)
                dta_bdy(jbdy) % h_i(ib, jl) = h_i(ii, ij, jl) * tmask(ii, ij, 1)
              END DO
            END DO
          END IF
          IF (dta % ll_h_s) THEN
            igrd = 1
            DO jl = 1, jpl
              DO ib = 1, ilen1(igrd)
                ii = idx_bdy(jbdy) % nbi(ib, igrd)
                ij = idx_bdy(jbdy) % nbj(ib, igrd)
                dta_bdy(jbdy) % h_s(ib, jl) = h_s(ii, ij, jl) * tmask(ii, ij, 1)
              END DO
            END DO
          END IF
        END IF
      END DO
    END IF
    jstart = 1
    DO jbdy = 1, nb_bdy
      dta => dta_bdy(jbdy)
      IF (nn_dta(jbdy) == 1) THEN
        IF (PRESENT(jit)) THEN
          IF (cn_dyn2d(jbdy) /= 'none') THEN
            IF (nn_dyn2d_dta(jbdy) == 2) THEN
              IF (dta % ll_ssh) dta % ssh(:) = 0._wp
              IF (dta % ll_u2d) dta % u2d(:) = 0._wp
              IF (dta % ll_u3d) dta % v2d(:) = 0._wp
            END IF
            IF (cn_tra(jbdy) /= 'runoff') THEN
              IF (nn_dyn2d_dta(jbdy) == 1 .OR. nn_dyn2d_dta(jbdy) == 3) THEN
                jend = jstart + dta % nread(2) - 1
                IF (ln_full_vel_array(jbdy)) THEN
                  CALL fld_read(kt = kt, kn_fsbc = 1, sd = bf(jstart : jend), map = nbmap_ptr(jstart : jend), kit = jit, kt_offset = time_offset, jpk_bdy = nb_jpk_bdy(jbdy), fvl = ln_full_vel_array(jbdy))
                ELSE
                  CALL fld_read(kt = kt, kn_fsbc = 1, sd = bf(jstart : jend), map = nbmap_ptr(jstart : jend), kit = jit, kt_offset = time_offset)
                END IF
                IF (ln_full_vel_array(jbdy) .AND. (nn_dyn2d_dta(jbdy) == 1 .OR. nn_dyn2d_dta(jbdy) == 3 .OR. nn_dyn3d_dta(jbdy) == 1)) THEN
                  igrd = 2
                  dta % u2d(:) = 0._wp
                  DO ib = 1, idx_bdy(jbdy) % nblen(igrd)
                    ii = idx_bdy(jbdy) % nbi(ib, igrd)
                    ij = idx_bdy(jbdy) % nbj(ib, igrd)
                    DO ik = 1, jpkm1
                      dta % u2d(ib) = dta % u2d(ib) + e3u_n(ii, ij, ik) * umask(ii, ij, ik) * dta % u3d(ib, ik)
                    END DO
                    dta % u2d(ib) = dta % u2d(ib) * r1_hu_n(ii, ij)
                  END DO
                  igrd = 3
                  dta % v2d(:) = 0._wp
                  DO ib = 1, idx_bdy(jbdy) % nblen(igrd)
                    ii = idx_bdy(jbdy) % nbi(ib, igrd)
                    ij = idx_bdy(jbdy) % nbj(ib, igrd)
                    DO ik = 1, jpkm1
                      dta % v2d(ib) = dta % v2d(ib) + e3v_n(ii, ij, ik) * vmask(ii, ij, ik) * dta % v3d(ib, ik)
                    END DO
                    dta % v2d(ib) = dta % v2d(ib) * r1_hv_n(ii, ij)
                  END DO
                END IF
              END IF
              IF (nn_dyn2d_dta(jbdy) .GE. 2) THEN
                CALL bdytide_update(kt = kt, idx = idx_bdy(jbdy), dta = dta, td = tides(jbdy), jit = jit, time_offset = time_offset)
              END IF
            END IF
          END IF
        ELSE
          IF (cn_tra(jbdy) == 'runoff') THEN
            jend = nb_bdy_fld(jbdy)
            CALL fld_read(kt = kt, kn_fsbc = 1, sd = bf(jstart : jend), map = nbmap_ptr(jstart : jend), kt_offset = time_offset)
            igrd = 2
            DO ib = 1, idx_bdy(jbdy) % nblen(igrd)
              ii = idx_bdy(jbdy) % nbi(ib, igrd)
              ij = idx_bdy(jbdy) % nbj(ib, igrd)
              dta % u2d(ib) = dta % u2d(ib) / (e2u(ii, ij) * hu_0(ii, ij))
            END DO
            igrd = 3
            DO ib = 1, idx_bdy(jbdy) % nblen(igrd)
              ii = idx_bdy(jbdy) % nbi(ib, igrd)
              ij = idx_bdy(jbdy) % nbj(ib, igrd)
              dta % v2d(ib) = dta % v2d(ib) / (e1v(ii, ij) * hv_0(ii, ij))
            END DO
          ELSE
            IF (nn_dyn2d_dta(jbdy) == 2) THEN
              IF (dta % ll_ssh) dta % ssh(:) = 0._wp
              IF (dta % ll_u2d) dta % u2d(:) = 0._wp
              IF (dta % ll_v2d) dta % v2d(:) = 0._wp
            END IF
            IF (dta % nread(1) .GT. 0) THEN
              jend = jstart + dta % nread(1) - 1
              CALL fld_read(kt = kt, kn_fsbc = 1, sd = bf(jstart : jend), map = nbmap_ptr(jstart : jend), kt_offset = time_offset, jpk_bdy = nb_jpk_bdy(jbdy), fvl = ln_full_vel_array(jbdy))
            END IF
            IF (ln_full_vel_array(jbdy) .AND. (nn_dyn2d_dta(jbdy) == 1 .OR. nn_dyn2d_dta(jbdy) == 3 .OR. nn_dyn3d_dta(jbdy) == 1)) THEN
              igrd = 2
              dta % u2d(:) = 0._wp
              DO ib = 1, idx_bdy(jbdy) % nblen(igrd)
                ii = idx_bdy(jbdy) % nbi(ib, igrd)
                ij = idx_bdy(jbdy) % nbj(ib, igrd)
                DO ik = 1, jpkm1
                  dta % u2d(ib) = dta % u2d(ib) + e3u_n(ii, ij, ik) * umask(ii, ij, ik) * dta % u3d(ib, ik)
                END DO
                dta % u2d(ib) = dta % u2d(ib) * r1_hu_n(ii, ij)
                DO ik = 1, jpkm1
                  dta % u3d(ib, ik) = dta % u3d(ib, ik) - dta % u2d(ib)
                END DO
              END DO
              igrd = 3
              dta % v2d(:) = 0._wp
              DO ib = 1, idx_bdy(jbdy) % nblen(igrd)
                ii = idx_bdy(jbdy) % nbi(ib, igrd)
                ij = idx_bdy(jbdy) % nbj(ib, igrd)
                DO ik = 1, jpkm1
                  dta % v2d(ib) = dta % v2d(ib) + e3v_n(ii, ij, ik) * vmask(ii, ij, ik) * dta % v3d(ib, ik)
                END DO
                dta % v2d(ib) = dta % v2d(ib) * r1_hv_n(ii, ij)
                DO ik = 1, jpkm1
                  dta % v3d(ib, ik) = dta % v3d(ib, ik) - dta % v2d(ib)
                END DO
              END DO
            END IF
          END IF
          IF (cn_ice(jbdy) /= 'none' .AND. nn_ice_dta(jbdy) == 1) THEN
            jfld_hti = jfld_htit(jbdy)
            jfld_hts = jfld_htst(jbdy)
            jfld_ai = jfld_ait(jbdy)
            IF (jpl /= 1 .AND. nice_cat == 1) THEN
              CALL ice_var_itd(bf(jfld_hti) % fnow(:, 1, 1), bf(jfld_hts) % fnow(:, 1, 1), bf(jfld_ai) % fnow(:, 1, 1), dta_bdy(jbdy) % h_i, dta_bdy(jbdy) % h_s, dta_bdy(jbdy) % a_i)
            ELSE IF (jpl /= 1 .AND. nice_cat /= 1 .AND. nice_cat /= jpl) THEN
              CALL ice_var_itd2(bf(jfld_hti) % fnow(:, 1, :), bf(jfld_hts) % fnow(:, 1, :), bf(jfld_ai) % fnow(:, 1, :), dta_bdy(jbdy) % h_i, dta_bdy(jbdy) % h_s, dta_bdy(jbdy) % a_i)
            END IF
          END IF
        END IF
        jstart = jstart + dta % nread(1)
      END IF
    END DO
    IF (ln_apr_obc) THEN
      DO jbdy = 1, nb_bdy
        IF (cn_tra(jbdy) /= 'runoff') THEN
          igrd = 1
          DO ib = 1, idx_bdy(jbdy) % nblenrim(igrd)
            ii = idx_bdy(jbdy) % nbi(ib, igrd)
            ij = idx_bdy(jbdy) % nbj(ib, igrd)
            dta_bdy(jbdy) % ssh(ib) = dta_bdy(jbdy) % ssh(ib) + ssh_ib(ii, ij)
          END DO
        END IF
      END DO
    END IF
    IF (ln_tide) THEN
      IF (ln_dynspg_ts) THEN
        DO jbdy = 1, nb_bdy
          IF (nn_dyn2d_dta(jbdy) .GE. 2) THEN
            nblen => idx_bdy(jbdy) % nblen
            nblenrim => idx_bdy(jbdy) % nblenrim
            IF (cn_dyn2d(jbdy) == 'frs') THEN
              ilen1(:) = nblen(:)
            ELSE
              ilen1(:) = nblenrim(:)
            END IF
            IF (dta_bdy(jbdy) % ll_ssh) dta_bdy_s(jbdy) % ssh(1 : ilen1(1)) = dta_bdy(jbdy) % ssh(1 : ilen1(1))
            IF (dta_bdy(jbdy) % ll_u2d) dta_bdy_s(jbdy) % u2d(1 : ilen1(2)) = dta_bdy(jbdy) % u2d(1 : ilen1(2))
            IF (dta_bdy(jbdy) % ll_v2d) dta_bdy_s(jbdy) % v2d(1 : ilen1(3)) = dta_bdy(jbdy) % v2d(1 : ilen1(3))
          END IF
        END DO
      ELSE
        CALL bdy_dta_tides(kt = kt, time_offset = time_offset)
      END IF
    END IF
    IF (ln_timing) CALL timing_stop('bdy_dta')
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE bdy_dta
  SUBROUTINE bdy_dta_init
    INTEGER :: jbdy, jfld, jstart, jend, ierror, ios
    CHARACTER(LEN = 100) :: cn_dir
    CHARACTER(LEN = 100), DIMENSION(nb_bdy) :: cn_dir_array
    CHARACTER(LEN = 256) :: clname
    LOGICAL :: ln_full_vel
    INTEGER :: ilen_global
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ilen1, ilen3
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ibdy
    INTEGER, ALLOCATABLE, DIMENSION(:) :: igrid
    INTEGER, POINTER, DIMENSION(:) :: nblen, nblenrim
    TYPE(OBC_DATA), POINTER :: dta
    INTEGER :: kndims
    INTEGER, DIMENSION(4) :: kdimsz
    INTEGER :: inum, id1
    TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) :: blf_i
    TYPE(FLD_N) :: bn_tem, bn_sal, bn_u3d, bn_v3d
    TYPE(FLD_N) :: bn_ssh, bn_u2d, bn_v2d
    TYPE(FLD_N) :: bn_a_i, bn_h_i, bn_h_s
    NAMELIST /nambdy_dta/ cn_dir, bn_tem, bn_sal, bn_u3d, bn_v3d, bn_ssh, bn_u2d, bn_v2d
    NAMELIST /nambdy_dta/ bn_a_i, bn_h_i, bn_h_s
    NAMELIST /nambdy_dta/ ln_full_vel
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) 'bdy_dta_ini : initialization of data at the open boundaries'
    IF (lwp) WRITE(numout, *) '~~~~~~~~~~'
    IF (lwp) WRITE(numout, *) ''
    DO jbdy = 1, nb_bdy
      nn_dta(jbdy) = MAX(nn_dyn2d_dta(jbdy), nn_dyn3d_dta(jbdy), nn_tra_dta(jbdy), nn_ice_dta(jbdy))
      IF (nn_dta(jbdy) > 1) nn_dta(jbdy) = 1
    END DO
    ALLOCATE(nb_bdy_fld(nb_bdy))
    nb_bdy_fld(:) = 0
    DO jbdy = 1, nb_bdy
      IF (cn_dyn2d(jbdy) /= 'none' .AND. (nn_dyn2d_dta(jbdy) == 1 .OR. nn_dyn2d_dta(jbdy) == 3)) THEN
        nb_bdy_fld(jbdy) = nb_bdy_fld(jbdy) + 3
      END IF
      IF (cn_dyn3d(jbdy) /= 'none' .AND. nn_dyn3d_dta(jbdy) == 1) THEN
        nb_bdy_fld(jbdy) = nb_bdy_fld(jbdy) + 2
      END IF
      IF (cn_tra(jbdy) /= 'none' .AND. nn_tra_dta(jbdy) == 1) THEN
        nb_bdy_fld(jbdy) = nb_bdy_fld(jbdy) + 2
      END IF
      IF (cn_ice(jbdy) /= 'none' .AND. nn_ice_dta(jbdy) == 1) THEN
        nb_bdy_fld(jbdy) = nb_bdy_fld(jbdy) + 3
      END IF
      IF (lwp) WRITE(numout, *) 'Maximum number of files to open =', nb_bdy_fld(jbdy)
    END DO
    nb_bdy_fld_sum = SUM(nb_bdy_fld)
    ALLOCATE(bf(nb_bdy_fld_sum), STAT = ierror)
    IF (ierror > 0) THEN
      CALL ctl_stop('bdy_dta: unable to allocate bf structure')
      RETURN
    END IF
    ALLOCATE(blf_i(nb_bdy_fld_sum), STAT = ierror)
    IF (ierror > 0) THEN
      CALL ctl_stop('bdy_dta: unable to allocate blf_i structure')
      RETURN
    END IF
    ALLOCATE(nbmap_ptr(nb_bdy_fld_sum), STAT = ierror)
    IF (ierror > 0) THEN
      CALL ctl_stop('bdy_dta: unable to allocate nbmap_ptr structure')
      RETURN
    END IF
    ALLOCATE(ilen1(nb_bdy_fld_sum), ilen3(nb_bdy_fld_sum))
    ALLOCATE(ibdy(nb_bdy_fld_sum))
    ALLOCATE(igrid(nb_bdy_fld_sum))
    REWIND(UNIT = numnam_cfg)
    jfld = 0
    DO jbdy = 1, nb_bdy
      IF (nn_dta(jbdy) == 1) THEN
        REWIND(UNIT = numnam_ref)
        READ(numnam_ref, nambdy_dta, IOSTAT = ios, ERR = 901)
901     IF (ios /= 0) CALL ctl_nam(ios, 'nambdy_dta in reference namelist', lwp)
        READ(numnam_cfg, nambdy_dta, IOSTAT = ios, ERR = 902)
902     IF (ios > 0) CALL ctl_nam(ios, 'nambdy_dta in configuration namelist', lwp)
        IF (lwm) WRITE(numond, nambdy_dta)
        cn_dir_array(jbdy) = cn_dir
        ln_full_vel_array(jbdy) = ln_full_vel
        nblen => idx_bdy(jbdy) % nblen
        nblenrim => idx_bdy(jbdy) % nblenrim
        dta => dta_bdy(jbdy)
        dta % nread(2) = 0
        IF (nn_dyn2d_dta(jbdy) == 1 .OR. nn_dyn2d_dta(jbdy) == 3) THEN
          IF (dta % ll_ssh) THEN
            IF (lwp) WRITE(numout, *) '++++++ reading in ssh field'
            jfld = jfld + 1
            blf_i(jfld) = bn_ssh
            ibdy(jfld) = jbdy
            igrid(jfld) = 1
            ilen1(jfld) = nblen(igrid(jfld))
            ilen3(jfld) = 1
            dta % nread(2) = dta % nread(2) + 1
          END IF
          IF (dta % ll_u2d .AND. .NOT. ln_full_vel_array(jbdy)) THEN
            IF (lwp) WRITE(numout, *) '++++++ reading in u2d field'
            jfld = jfld + 1
            blf_i(jfld) = bn_u2d
            ibdy(jfld) = jbdy
            igrid(jfld) = 2
            ilen1(jfld) = nblen(igrid(jfld))
            ilen3(jfld) = 1
            dta % nread(2) = dta % nread(2) + 1
          END IF
          IF (dta % ll_v2d .AND. .NOT. ln_full_vel_array(jbdy)) THEN
            IF (lwp) WRITE(numout, *) '++++++ reading in v2d field'
            jfld = jfld + 1
            blf_i(jfld) = bn_v2d
            ibdy(jfld) = jbdy
            igrid(jfld) = 3
            ilen1(jfld) = nblen(igrid(jfld))
            ilen3(jfld) = 1
            dta % nread(2) = dta % nread(2) + 1
          END IF
        END IF
        IF (nn_dyn3d_dta(jbdy) == 1 .OR. (ln_full_vel_array(jbdy) .AND. (nn_dyn2d_dta(jbdy) == 1 .OR. nn_dyn2d_dta(jbdy) == 3))) THEN
          IF (dta % ll_u3d .OR. (ln_full_vel_array(jbdy) .AND. dta % ll_u2d)) THEN
            IF (lwp) WRITE(numout, *) '++++++ reading in u3d field'
            jfld = jfld + 1
            blf_i(jfld) = bn_u3d
            ibdy(jfld) = jbdy
            igrid(jfld) = 2
            ilen1(jfld) = nblen(igrid(jfld))
            ilen3(jfld) = jpk
            IF (ln_full_vel_array(jbdy) .AND. dta % ll_u2d) dta % nread(2) = dta % nread(2) + 1
          END IF
          IF (dta % ll_v3d .OR. (ln_full_vel_array(jbdy) .AND. dta % ll_v2d)) THEN
            IF (lwp) WRITE(numout, *) '++++++ reading in v3d field'
            jfld = jfld + 1
            blf_i(jfld) = bn_v3d
            ibdy(jfld) = jbdy
            igrid(jfld) = 3
            ilen1(jfld) = nblen(igrid(jfld))
            ilen3(jfld) = jpk
            IF (ln_full_vel_array(jbdy) .AND. dta % ll_v2d) dta % nread(2) = dta % nread(2) + 1
          END IF
        END IF
        IF (nn_tra_dta(jbdy) == 1) THEN
          IF (dta % ll_tem) THEN
            IF (lwp) WRITE(numout, *) '++++++ reading in tem field'
            jfld = jfld + 1
            blf_i(jfld) = bn_tem
            ibdy(jfld) = jbdy
            igrid(jfld) = 1
            ilen1(jfld) = nblen(igrid(jfld))
            ilen3(jfld) = jpk
          END IF
          IF (dta % ll_sal) THEN
            IF (lwp) WRITE(numout, *) '++++++ reading in sal field'
            jfld = jfld + 1
            blf_i(jfld) = bn_sal
            ibdy(jfld) = jbdy
            igrid(jfld) = 1
            ilen1(jfld) = nblen(igrid(jfld))
            ilen3(jfld) = jpk
          END IF
        END IF
        IF (nn_ice_dta(jbdy) == 1) THEN
          clname = TRIM(cn_dir) // TRIM(bn_a_i % clname)
          IF (.NOT. bn_a_i % ln_clim) THEN
            WRITE(clname, '(a,"_y",i4.4)') TRIM(clname), nyear
            IF (bn_a_i % cltype /= 'yearly') WRITE(clname, '(a,"m" ,i2.2)') TRIM(clname), nmonth
          ELSE
            IF (bn_a_i % cltype /= 'yearly') WRITE(clname, '(a,"_m",i2.2)') TRIM(clname), nmonth
          END IF
          IF (bn_a_i % cltype == 'daily' .OR. bn_a_i % cltype(1 : 4) == 'week') WRITE(clname, '(a,"d" ,i2.2)') TRIM(clname), nday
          CALL iom_open(clname, inum)
          id1 = iom_varid(inum, bn_a_i % clvar, kdimsz = kdimsz, kndims = kndims, ldstop = .FALSE.)
          CALL iom_close(inum)
          IF (kndims == 4) THEN
            nice_cat = kdimsz(4)
          ELSE
            nice_cat = 1
          END IF
          IF (dta % ll_a_i) THEN
            jfld = jfld + 1
            blf_i(jfld) = bn_a_i
            ibdy(jfld) = jbdy
            igrid(jfld) = 1
            ilen1(jfld) = nblen(igrid(jfld))
            ilen3(jfld) = nice_cat
          END IF
          IF (dta % ll_h_i) THEN
            jfld = jfld + 1
            blf_i(jfld) = bn_h_i
            ibdy(jfld) = jbdy
            igrid(jfld) = 1
            ilen1(jfld) = nblen(igrid(jfld))
            ilen3(jfld) = nice_cat
          END IF
          IF (dta % ll_h_s) THEN
            jfld = jfld + 1
            blf_i(jfld) = bn_h_s
            ibdy(jfld) = jbdy
            igrid(jfld) = 1
            ilen1(jfld) = nblen(igrid(jfld))
            ilen3(jfld) = nice_cat
          END IF
        END IF
        IF (jbdy == 1) THEN
          nb_bdy_fld_sum = 0
          nb_bdy_fld(jbdy) = jfld
          nb_bdy_fld_sum = jfld
        ELSE
          nb_bdy_fld(jbdy) = jfld - nb_bdy_fld_sum
          nb_bdy_fld_sum = nb_bdy_fld_sum + nb_bdy_fld(jbdy)
        END IF
        dta % nread(1) = nb_bdy_fld(jbdy)
      END IF
    END DO
    DO jfld = 1, nb_bdy_fld_sum
      ALLOCATE(bf(jfld) % fnow(ilen1(jfld), 1, ilen3(jfld)))
      IF (blf_i(jfld) % ln_tint) ALLOCATE(bf(jfld) % fdta(ilen1(jfld), 1, ilen3(jfld), 2))
      nbmap_ptr(jfld) % ptr => idx_bdy(ibdy(jfld)) % nbmap(:, igrid(jfld))
      nbmap_ptr(jfld) % ll_unstruc = ln_coords_file(ibdy(jfld))
    END DO
    jstart = 1
    DO jbdy = 1, nb_bdy
      jend = jstart - 1 + nb_bdy_fld(jbdy)
      CALL fld_fill(bf(jstart : jend), blf_i(jstart : jend), cn_dir_array(jbdy), 'bdy_dta', 'open boundary conditions', 'nambdy_dta')
      jstart = jend + 1
    END DO
    DO jfld = 1, nb_bdy_fld_sum
      bf(jfld) % igrd = igrid(jfld)
      bf(jfld) % ibdy = ibdy(jfld)
    END DO
    jfld = 0
    DO jbdy = 1, nb_bdy
      nblen => idx_bdy(jbdy) % nblen
      dta => dta_bdy(jbdy)
      IF (lwp) THEN
        WRITE(numout, *) '++++++ dta%ll_ssh = ', dta % ll_ssh
        WRITE(numout, *) '++++++ dta%ll_u2d = ', dta % ll_u2d
        WRITE(numout, *) '++++++ dta%ll_v2d = ', dta % ll_v2d
        WRITE(numout, *) '++++++ dta%ll_u3d = ', dta % ll_u3d
        WRITE(numout, *) '++++++ dta%ll_v3d = ', dta % ll_v3d
        WRITE(numout, *) '++++++ dta%ll_tem = ', dta % ll_tem
        WRITE(numout, *) '++++++ dta%ll_sal = ', dta % ll_sal
      END IF
      IF (nn_dyn2d_dta(jbdy) == 0 .OR. nn_dyn2d_dta(jbdy) == 2) THEN
        IF (lwp) WRITE(numout, *) '++++++ dta%ssh/u2d/u3d allocated space'
        IF (dta % ll_ssh) ALLOCATE(dta % ssh(nblen(1)))
        IF (dta % ll_u2d) ALLOCATE(dta % u2d(nblen(2)))
        IF (dta % ll_v2d) ALLOCATE(dta % v2d(nblen(3)))
      END IF
      IF (nn_dyn2d_dta(jbdy) == 1 .OR. nn_dyn2d_dta(jbdy) == 3) THEN
        IF (dta % ll_ssh) THEN
          IF (lwp) WRITE(numout, *) '++++++ dta%ssh pointing to fnow'
          jfld = jfld + 1
          dta % ssh => bf(jfld) % fnow(:, 1, 1)
        END IF
        IF (dta % ll_u2d) THEN
          IF (ln_full_vel_array(jbdy)) THEN
            IF (lwp) WRITE(numout, *) '++++++ dta%u2d allocated space'
            ALLOCATE(dta % u2d(nblen(2)))
          ELSE
            IF (lwp) WRITE(numout, *) '++++++ dta%u2d pointing to fnow'
            jfld = jfld + 1
            dta % u2d => bf(jfld) % fnow(:, 1, 1)
          END IF
        END IF
        IF (dta % ll_v2d) THEN
          IF (ln_full_vel_array(jbdy)) THEN
            IF (lwp) WRITE(numout, *) '++++++ dta%v2d allocated space'
            ALLOCATE(dta % v2d(nblen(3)))
          ELSE
            IF (lwp) WRITE(numout, *) '++++++ dta%v2d pointing to fnow'
            jfld = jfld + 1
            dta % v2d => bf(jfld) % fnow(:, 1, 1)
          END IF
        END IF
      END IF
      IF (nn_dyn3d_dta(jbdy) == 0) THEN
        IF (lwp) WRITE(numout, *) '++++++ dta%u3d/v3d allocated space'
        IF (dta % ll_u3d) ALLOCATE(dta_bdy(jbdy) % u3d(nblen(2), jpk))
        IF (dta % ll_v3d) ALLOCATE(dta_bdy(jbdy) % v3d(nblen(3), jpk))
      END IF
      IF (nn_dyn3d_dta(jbdy) == 1 .OR. (ln_full_vel_array(jbdy) .AND. (nn_dyn2d_dta(jbdy) == 1 .OR. nn_dyn2d_dta(jbdy) == 3))) THEN
        IF (dta % ll_u3d .OR. (ln_full_vel_array(jbdy) .AND. dta % ll_u2d)) THEN
          IF (lwp) WRITE(numout, *) '++++++ dta%u3d pointing to fnow'
          jfld = jfld + 1
          dta_bdy(jbdy) % u3d => bf(jfld) % fnow(:, 1, :)
        END IF
        IF (dta % ll_v3d .OR. (ln_full_vel_array(jbdy) .AND. dta % ll_v2d)) THEN
          IF (lwp) WRITE(numout, *) '++++++ dta%v3d pointing to fnow'
          jfld = jfld + 1
          dta_bdy(jbdy) % v3d => bf(jfld) % fnow(:, 1, :)
        END IF
      END IF
      IF (nn_tra_dta(jbdy) == 0) THEN
        IF (lwp) WRITE(numout, *) '++++++ dta%tem/sal allocated space'
        IF (dta % ll_tem) ALLOCATE(dta_bdy(jbdy) % tem(nblen(1), jpk))
        IF (dta % ll_sal) ALLOCATE(dta_bdy(jbdy) % sal(nblen(1), jpk))
      ELSE
        IF (dta % ll_tem) THEN
          IF (lwp) WRITE(numout, *) '++++++ dta%tem pointing to fnow'
          jfld = jfld + 1
          dta_bdy(jbdy) % tem => bf(jfld) % fnow(:, 1, :)
        END IF
        IF (dta % ll_sal) THEN
          IF (lwp) WRITE(numout, *) '++++++ dta%sal pointing to fnow'
          jfld = jfld + 1
          dta_bdy(jbdy) % sal => bf(jfld) % fnow(:, 1, :)
        END IF
      END IF
      IF (cn_ice(jbdy) /= 'none') THEN
        IF (nn_ice_dta(jbdy) == 0) THEN
          ALLOCATE(dta_bdy(jbdy) % a_i(nblen(1), jpl))
          ALLOCATE(dta_bdy(jbdy) % h_i(nblen(1), jpl))
          ALLOCATE(dta_bdy(jbdy) % h_s(nblen(1), jpl))
        ELSE
          IF (nice_cat == jpl) THEN
            jfld = jfld + 1
            dta_bdy(jbdy) % a_i => bf(jfld) % fnow(:, 1, :)
            jfld = jfld + 1
            dta_bdy(jbdy) % h_i => bf(jfld) % fnow(:, 1, :)
            jfld = jfld + 1
            dta_bdy(jbdy) % h_s => bf(jfld) % fnow(:, 1, :)
          ELSE
            jfld_ait(jbdy) = jfld + 1
            jfld_htit(jbdy) = jfld + 2
            jfld_htst(jbdy) = jfld + 3
            jfld = jfld + 3
            ALLOCATE(dta_bdy(jbdy) % a_i(nblen(1), jpl))
            ALLOCATE(dta_bdy(jbdy) % h_i(nblen(1), jpl))
            ALLOCATE(dta_bdy(jbdy) % h_s(nblen(1), jpl))
            dta_bdy(jbdy) % a_i(:, :) = 0._wp
            dta_bdy(jbdy) % h_i(:, :) = 0._wp
            dta_bdy(jbdy) % h_s(:, :) = 0._wp
          END IF
        END IF
      END IF
    END DO
  END SUBROUTINE bdy_dta_init
END MODULE bdydta