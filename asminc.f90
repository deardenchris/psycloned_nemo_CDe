MODULE asminc
  USE oce
  USE par_oce
  USE dom_oce
  USE domvvl
  USE ldfdyn
  USE eosbn2
  USE zpshde
  USE asmpar
  USE asmbkg
  USE c1d
  USE sbc_oce
  USE diaobs, ONLY: calc_date
  USE ice, ONLY: hm_i, at_i, at_i_b
  USE in_out_manager
  USE iom
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: asm_inc_init
  PUBLIC :: tra_asm_inc
  PUBLIC :: dyn_asm_inc
  PUBLIC :: ssh_asm_inc
  PUBLIC :: ssh_asm_div
  PUBLIC :: seaice_asm_inc
  LOGICAL, PUBLIC, PARAMETER :: lk_asminc = .FALSE.
  LOGICAL, PUBLIC :: ln_bkgwri
  LOGICAL, PUBLIC :: ln_asmiau
  LOGICAL, PUBLIC :: ln_asmdin
  LOGICAL, PUBLIC :: ln_trainc
  LOGICAL, PUBLIC :: ln_dyninc
  LOGICAL, PUBLIC :: ln_sshinc
  LOGICAL, PUBLIC :: ln_seaiceinc
  LOGICAL, PUBLIC :: ln_salfix
  LOGICAL, PUBLIC :: ln_temnofreeze = .FALSE.
  INTEGER, PUBLIC :: nn_divdmp
  REAL(KIND = wp), PUBLIC, DIMENSION(:, :, :), ALLOCATABLE :: t_bkg, s_bkg
  REAL(KIND = wp), PUBLIC, DIMENSION(:, :, :), ALLOCATABLE :: u_bkg, v_bkg
  REAL(KIND = wp), PUBLIC, DIMENSION(:, :, :), ALLOCATABLE :: t_bkginc, s_bkginc
  REAL(KIND = wp), PUBLIC, DIMENSION(:, :, :), ALLOCATABLE :: u_bkginc, v_bkginc
  REAL(KIND = wp), PUBLIC, DIMENSION(:), ALLOCATABLE :: wgtiau
  INTEGER, PUBLIC :: nitbkg
  INTEGER, PUBLIC :: nitdin
  INTEGER, PUBLIC :: nitiaustr
  INTEGER, PUBLIC :: nitiaufin
  INTEGER, PUBLIC :: niaufn
  REAL(KIND = wp), PUBLIC :: salfixmin
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: ssh_bkg, ssh_bkginc
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: seaice_bkginc
  CONTAINS
  SUBROUTINE asm_inc_init
    INTEGER :: ji, jj, jk, jt
    INTEGER :: imid, inum
    INTEGER :: ios
    INTEGER :: iiauper
    INTEGER :: icycper
    REAL(KIND = dp) :: ditend_date
    REAL(KIND = dp) :: ditbkg_date
    REAL(KIND = dp) :: ditdin_date
    REAL(KIND = dp) :: ditiaustr_date
    REAL(KIND = dp) :: ditiaufin_date
    REAL(KIND = wp) :: znorm
    REAL(KIND = wp) :: ztotwgt
    REAL(KIND = wp) :: z_inc_dateb
    REAL(KIND = wp) :: z_inc_datef
    REAL(KIND = wp) :: zdate_bkg
    REAL(KIND = wp) :: zdate_inc
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zhdiv
    NAMELIST /nam_asminc/ ln_bkgwri, ln_trainc, ln_dyninc, ln_sshinc, ln_asmdin, ln_asmiau, nitbkg, nitdin, nitiaustr, nitiaufin, niaufn, ln_salfix, salfixmin, nn_divdmp
    ln_seaiceinc = .FALSE.
    ln_temnofreeze = .FALSE.
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, nam_asminc, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'nam_asminc in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, nam_asminc, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'nam_asminc in configuration namelist', lwp)
    IF (lwm) WRITE(numond, nam_asminc)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'asm_inc_init : Assimilation increment initialization :'
      WRITE(numout, *) '~~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist namasm : set assimilation increment parameters'
      WRITE(numout, *) '      Logical switch for writing out background state          ln_bkgwri = ', ln_bkgwri
      WRITE(numout, *) '      Logical switch for applying tracer increments            ln_trainc = ', ln_trainc
      WRITE(numout, *) '      Logical switch for applying velocity increments          ln_dyninc = ', ln_dyninc
      WRITE(numout, *) '      Logical switch for applying SSH increments               ln_sshinc = ', ln_sshinc
      WRITE(numout, *) '      Logical switch for Direct Initialization (DI)            ln_asmdin = ', ln_asmdin
      WRITE(numout, *) '      Logical switch for applying sea ice increments        ln_seaiceinc = ', ln_seaiceinc
      WRITE(numout, *) '      Logical switch for Incremental Analysis Updating (IAU)   ln_asmiau = ', ln_asmiau
      WRITE(numout, *) '      Timestep of background in [0,nitend-nit000-1]            nitbkg    = ', nitbkg
      WRITE(numout, *) '      Timestep of background for DI in [0,nitend-nit000-1]     nitdin    = ', nitdin
      WRITE(numout, *) '      Timestep of start of IAU interval in [0,nitend-nit000-1] nitiaustr = ', nitiaustr
      WRITE(numout, *) '      Timestep of end of IAU interval in [0,nitend-nit000-1]   nitiaufin = ', nitiaufin
      WRITE(numout, *) '      Type of IAU weighting function                           niaufn    = ', niaufn
      WRITE(numout, *) '      Logical switch for ensuring that the sa > salfixmin      ln_salfix = ', ln_salfix
      WRITE(numout, *) '      Minimum salinity after applying the increments           salfixmin = ', salfixmin
    END IF
    nitbkg_r = nitbkg + nit000 - 1
    nitdin_r = nitdin + nit000 - 1
    nitiaustr_r = nitiaustr + nit000 - 1
    nitiaufin_r = nitiaufin + nit000 - 1
    iiauper = nitiaufin_r - nitiaustr_r + 1
    icycper = nitend - nit000 + 1
    CALL calc_date(nitend, ditend_date)
    CALL calc_date(nitbkg_r, ditbkg_date)
    CALL calc_date(nitdin_r, ditdin_date)
    CALL calc_date(nitiaustr_r, ditiaustr_date)
    CALL calc_date(nitiaufin_r, ditiaufin_date)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) '   Time steps referenced to current cycle:'
      WRITE(numout, *) '       iitrst      = ', nit000 - 1
      WRITE(numout, *) '       nit000      = ', nit000
      WRITE(numout, *) '       nitend      = ', nitend
      WRITE(numout, *) '       nitbkg_r    = ', nitbkg_r
      WRITE(numout, *) '       nitdin_r    = ', nitdin_r
      WRITE(numout, *) '       nitiaustr_r = ', nitiaustr_r
      WRITE(numout, *) '       nitiaufin_r = ', nitiaufin_r
      WRITE(numout, *)
      WRITE(numout, *) '   Dates referenced to current cycle:'
      WRITE(numout, *) '       ndastp         = ', ndastp
      WRITE(numout, *) '       ndate0         = ', ndate0
      WRITE(numout, *) '       nn_time0       = ', nn_time0
      WRITE(numout, *) '       ditend_date    = ', ditend_date
      WRITE(numout, *) '       ditbkg_date    = ', ditbkg_date
      WRITE(numout, *) '       ditdin_date    = ', ditdin_date
      WRITE(numout, *) '       ditiaustr_date = ', ditiaustr_date
      WRITE(numout, *) '       ditiaufin_date = ', ditiaufin_date
    END IF
    IF ((ln_asmdin) .AND. (ln_asmiau)) CALL ctl_stop(' ln_asmdin and ln_asmiau :', ' Choose Direct Initialization OR Incremental Analysis Updating')
    IF (((.NOT. ln_asmdin) .AND. (.NOT. ln_asmiau)) .AND. ((ln_trainc) .OR. (ln_dyninc) .OR. (ln_sshinc) .OR. (ln_seaiceinc))) CALL ctl_stop(' One or more of ln_trainc, ln_dyninc, ln_sshinc and ln_seaiceinc is set to .true.', ' but ln_asmdin and ln_asmiau are both set to .false. :', ' Inconsistent options')
    IF ((niaufn /= 0) .AND. (niaufn /= 1)) CALL ctl_stop(' niaufn /= 0 or niaufn /=1 :', ' Type IAU weighting function is invalid')
    IF ((.NOT. ln_trainc) .AND. (.NOT. ln_dyninc) .AND. (.NOT. ln_sshinc) .AND. (.NOT. ln_seaiceinc)) CALL ctl_warn(' ln_trainc, ln_dyninc, ln_sshinc and ln_seaiceinc are set to .false. :', ' The assimilation increments are not applied')
    IF ((ln_asmiau) .AND. (nitiaustr == nitiaufin)) CALL ctl_stop(' nitiaustr = nitiaufin :', ' IAU interval is of zero length')
    IF ((ln_asmiau) .AND. ((nitiaustr_r < nit000) .OR. (nitiaufin_r > nitend))) CALL ctl_stop(' nitiaustr or nitiaufin :', ' IAU starting or final time step is outside the cycle interval', ' Valid range nit000 to nitend')
    IF ((nitbkg_r < nit000 - 1) .OR. (nitbkg_r > nitend)) CALL ctl_stop(' nitbkg :', ' Background time step is outside the cycle interval')
    IF ((nitdin_r < nit000 - 1) .OR. (nitdin_r > nitend)) CALL ctl_stop(' nitdin :', ' Background time step for Direct Initialization is outside', ' the cycle interval')
    IF (nstop > 0) RETURN
    IF (ln_asmiau) THEN
      ALLOCATE(wgtiau(icycper))
      wgtiau(:) = 0._wp
      IF (niaufn == 0) THEN
        DO jt = 1, iiauper
          wgtiau(jt + nitiaustr - 1) = 1.0 / REAL(iiauper)
        END DO
      ELSE IF (niaufn == 1) THEN
        znorm = 0._wp
        IF (MOD(iiauper, 2) == 0) THEN
          imid = iiauper / 2
          DO jt = 1, imid
            znorm = znorm + REAL(jt)
          END DO
          znorm = 2.0 * znorm
        ELSE
          imid = (iiauper + 1) / 2
          DO jt = 1, imid - 1
            znorm = znorm + REAL(jt)
          END DO
          znorm = 2.0 * znorm + REAL(imid)
        END IF
        znorm = 1.0 / znorm
        DO jt = 1, imid - 1
          wgtiau(jt + nitiaustr - 1) = REAL(jt) * znorm
        END DO
        DO jt = imid, iiauper
          wgtiau(jt + nitiaustr - 1) = REAL(iiauper - jt + 1) * znorm
        END DO
      END IF
      IF (lwp) THEN
        WRITE(numout, *)
        WRITE(numout, *) 'asm_inc_init : IAU weights'
        WRITE(numout, *) '~~~~~~~~~~~~'
        WRITE(numout, *) '             time step         IAU  weight'
        WRITE(numout, *) '             =========     ====================='
        ztotwgt = 0.0
        DO jt = 1, icycper
          ztotwgt = ztotwgt + wgtiau(jt)
          WRITE(numout, *) '         ', jt, '       ', wgtiau(jt)
        END DO
        WRITE(numout, *) '         ==================================='
        WRITE(numout, *) '         Time-integrated weight = ', ztotwgt
        WRITE(numout, *) '         ==================================='
      END IF
    END IF
    ALLOCATE(t_bkginc(jpi, jpj, jpk))
    t_bkginc(:, :, :) = 0._wp
    ALLOCATE(s_bkginc(jpi, jpj, jpk))
    s_bkginc(:, :, :) = 0._wp
    ALLOCATE(u_bkginc(jpi, jpj, jpk))
    u_bkginc(:, :, :) = 0._wp
    ALLOCATE(v_bkginc(jpi, jpj, jpk))
    v_bkginc(:, :, :) = 0._wp
    ALLOCATE(ssh_bkginc(jpi, jpj))
    ssh_bkginc(:, :) = 0._wp
    ALLOCATE(seaice_bkginc(jpi, jpj))
    seaice_bkginc(:, :) = 0._wp
    IF (ln_trainc .OR. ln_dyninc .OR. ln_sshinc .OR. ln_seaiceinc) THEN
      CALL iom_open(c_asminc, inum)
      CALL iom_get(inum, 'time', zdate_inc)
      CALL iom_get(inum, 'z_inc_dateb', z_inc_dateb)
      CALL iom_get(inum, 'z_inc_datef', z_inc_datef)
      z_inc_dateb = zdate_inc
      z_inc_datef = zdate_inc
      IF (lwp) THEN
        WRITE(numout, *)
        WRITE(numout, *) 'asm_inc_init : Assimilation increments valid between dates ', z_inc_dateb, ' and ', z_inc_datef
        WRITE(numout, *) '~~~~~~~~~~~~'
      END IF
      IF ((z_inc_dateb < ndastp + nn_time0 * 0.0001_wp) .OR. (z_inc_datef > ditend_date)) CALL ctl_warn(' Validity time of assimilation increments is ', ' outside the assimilation interval')
      IF ((ln_asmdin) .AND. (zdate_inc /= ditdin_date)) CALL ctl_warn(' Validity time of assimilation increments does ', ' not agree with Direct Initialization time')
      IF (ln_trainc) THEN
        CALL iom_get(inum, jpdom_autoglo, 'bckint', t_bkginc, 1)
        CALL iom_get(inum, jpdom_autoglo, 'bckins', s_bkginc, 1)
        t_bkginc(:, :, :) = t_bkginc(:, :, :) * tmask(:, :, :)
        s_bkginc(:, :, :) = s_bkginc(:, :, :) * tmask(:, :, :)
        WHERE (ABS(t_bkginc(:, :, :)) > 1.0E+10) t_bkginc(:, :, :) = 0.0
        WHERE (ABS(s_bkginc(:, :, :)) > 1.0E+10) s_bkginc(:, :, :) = 0.0
      END IF
      IF (ln_dyninc) THEN
        CALL iom_get(inum, jpdom_autoglo, 'bckinu', u_bkginc, 1)
        CALL iom_get(inum, jpdom_autoglo, 'bckinv', v_bkginc, 1)
        u_bkginc(:, :, :) = u_bkginc(:, :, :) * umask(:, :, :)
        v_bkginc(:, :, :) = v_bkginc(:, :, :) * vmask(:, :, :)
        WHERE (ABS(u_bkginc(:, :, :)) > 1.0E+10) u_bkginc(:, :, :) = 0.0
        WHERE (ABS(v_bkginc(:, :, :)) > 1.0E+10) v_bkginc(:, :, :) = 0.0
      END IF
      IF (ln_sshinc) THEN
        CALL iom_get(inum, jpdom_autoglo, 'bckineta', ssh_bkginc, 1)
        ssh_bkginc(:, :) = ssh_bkginc(:, :) * tmask(:, :, 1)
        WHERE (ABS(ssh_bkginc(:, :)) > 1.0E+10) ssh_bkginc(:, :) = 0.0
      END IF
      IF (ln_seaiceinc) THEN
        CALL iom_get(inum, jpdom_autoglo, 'bckinseaice', seaice_bkginc, 1)
        seaice_bkginc(:, :) = seaice_bkginc(:, :) * tmask(:, :, 1)
        WHERE (ABS(seaice_bkginc(:, :)) > 1.0E+10) seaice_bkginc(:, :) = 0.0
      END IF
      CALL iom_close(inum)
    END IF
    IF (ln_dyninc .AND. nn_divdmp > 0) THEN
      ALLOCATE(zhdiv(jpi, jpj))
      DO jt = 1, nn_divdmp
        DO jk = 1, jpkm1
          zhdiv(:, :) = 0._wp
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zhdiv(ji, jj) = (e2u(ji, jj) * e3u_n(ji, jj, jk) * u_bkginc(ji, jj, jk) - e2u(ji - 1, jj) * e3u_n(ji - 1, jj, jk) * u_bkginc(ji - 1, jj, jk) + e1v(ji, jj) * e3v_n(ji, jj, jk) * v_bkginc(ji, jj, jk) - e1v(ji, jj - 1) * e3v_n(ji, jj - 1, jk) * v_bkginc(ji, jj - 1, jk)) / e3t_n(ji, jj, jk)
            END DO
          END DO
          CALL lbc_lnk('asminc', zhdiv, 'T', 1.)
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              u_bkginc(ji, jj, jk) = u_bkginc(ji, jj, jk) + 0.2_wp * (zhdiv(ji + 1, jj) - zhdiv(ji, jj)) * r1_e1u(ji, jj) * umask(ji, jj, jk)
              v_bkginc(ji, jj, jk) = v_bkginc(ji, jj, jk) + 0.2_wp * (zhdiv(ji, jj + 1) - zhdiv(ji, jj)) * r1_e2v(ji, jj) * vmask(ji, jj, jk)
            END DO
          END DO
        END DO
      END DO
      DEALLOCATE(zhdiv)
    END IF
    IF (ln_asmdin) THEN
      ALLOCATE(t_bkg(jpi, jpj, jpk))
      t_bkg(:, :, :) = 0._wp
      ALLOCATE(s_bkg(jpi, jpj, jpk))
      s_bkg(:, :, :) = 0._wp
      ALLOCATE(u_bkg(jpi, jpj, jpk))
      u_bkg(:, :, :) = 0._wp
      ALLOCATE(v_bkg(jpi, jpj, jpk))
      v_bkg(:, :, :) = 0._wp
      ALLOCATE(ssh_bkg(jpi, jpj))
      ssh_bkg(:, :) = 0._wp
      CALL iom_open(c_asmdin, inum)
      CALL iom_get(inum, 'rdastp', zdate_bkg)
      IF (lwp) THEN
        WRITE(numout, *)
        WRITE(numout, *) '   ==>>>  Assimilation background state valid at : ', zdate_bkg
        WRITE(numout, *)
      END IF
      IF (zdate_bkg /= ditdin_date) CALL ctl_warn(' Validity time of assimilation background state does', ' not agree with Direct Initialization time')
      IF (ln_trainc) THEN
        CALL iom_get(inum, jpdom_autoglo, 'tn', t_bkg)
        CALL iom_get(inum, jpdom_autoglo, 'sn', s_bkg)
        t_bkg(:, :, :) = t_bkg(:, :, :) * tmask(:, :, :)
        s_bkg(:, :, :) = s_bkg(:, :, :) * tmask(:, :, :)
      END IF
      IF (ln_dyninc) THEN
        CALL iom_get(inum, jpdom_autoglo, 'un', u_bkg)
        CALL iom_get(inum, jpdom_autoglo, 'vn', v_bkg)
        u_bkg(:, :, :) = u_bkg(:, :, :) * umask(:, :, :)
        v_bkg(:, :, :) = v_bkg(:, :, :) * vmask(:, :, :)
      END IF
      IF (ln_sshinc) THEN
        CALL iom_get(inum, jpdom_autoglo, 'sshn', ssh_bkg)
        ssh_bkg(:, :) = ssh_bkg(:, :) * tmask(:, :, 1)
      END IF
      CALL iom_close(inum)
    END IF
    IF (lwp) WRITE(numout, *) '   ==>>>   Euler time step switch is ', neuler
    IF (lk_asminc) THEN
      IF (ln_bkgwri) CALL asm_bkg_wri(nit000 - 1)
      IF (ln_asmdin) THEN
        IF (ln_trainc) CALL tra_asm_inc(nit000 - 1)
        IF (ln_dyninc) CALL dyn_asm_inc(nit000 - 1)
        IF (ln_sshinc) CALL ssh_asm_inc(nit000 - 1)
      END IF
    END IF
  END SUBROUTINE asm_inc_init
  SUBROUTINE tra_asm_inc(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: it
    REAL(KIND = wp) :: zincwgt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: fzptnz
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    CALL profile_psy_data0 % PreStart('tra_asm_inc', 'r0', 0, 0)
    DO jk = 1, jpkm1
      CALL eos_fzp(tsn(:, :, jk, jp_sal), fzptnz(:, :, jk), gdept_n(:, :, jk))
    END DO
    CALL profile_psy_data0 % PostEnd
    IF (ln_asmiau) THEN
      IF ((kt >= nitiaustr_r) .AND. (kt <= nitiaufin_r)) THEN
        CALL profile_psy_data1 % PreStart('tra_asm_inc', 'r1', 0, 0)
        it = kt - nit000 + 1
        zincwgt = wgtiau(it) / rdt
        IF (lwp) THEN
          WRITE(numout, *)
          WRITE(numout, *) 'tra_asm_inc : Tracer IAU at time step = ', kt, ' with IAU weight = ', wgtiau(it)
          WRITE(numout, *) '~~~~~~~~~~~~'
        END IF
        CALL profile_psy_data1 % PostEnd
        DO jk = 1, jpkm1
          IF (ln_temnofreeze) THEN
            WHERE (t_bkginc(:, :, jk) > 0.0_wp .OR. tsn(:, :, jk, jp_tem) + tsa(:, :, jk, jp_tem) + t_bkginc(:, :, jk) * wgtiau(it) > fzptnz(:, :, jk))
              tsa(:, :, jk, jp_tem) = tsa(:, :, jk, jp_tem) + t_bkginc(:, :, jk) * zincwgt
            END WHERE
          ELSE
            !$ACC KERNELS
            tsa(:, :, jk, jp_tem) = tsa(:, :, jk, jp_tem) + t_bkginc(:, :, jk) * zincwgt
            !$ACC END KERNELS
          END IF
          IF (ln_salfix) THEN
            WHERE (s_bkginc(:, :, jk) > 0.0_wp .OR. tsn(:, :, jk, jp_sal) + tsa(:, :, jk, jp_sal) + s_bkginc(:, :, jk) * wgtiau(it) > salfixmin)
              tsa(:, :, jk, jp_sal) = tsa(:, :, jk, jp_sal) + s_bkginc(:, :, jk) * zincwgt
            END WHERE
          ELSE
            !$ACC KERNELS
            tsa(:, :, jk, jp_sal) = tsa(:, :, jk, jp_sal) + s_bkginc(:, :, jk) * zincwgt
            !$ACC END KERNELS
          END IF
        END DO
      END IF
      CALL profile_psy_data2 % PreStart('tra_asm_inc', 'r2', 0, 0)
      IF (kt == nitiaufin_r + 1) THEN
        DEALLOCATE(t_bkginc)
        DEALLOCATE(s_bkginc)
      END IF
      CALL profile_psy_data2 % PostEnd
    ELSE IF (ln_asmdin) THEN
      IF (kt == nitdin_r) THEN
        CALL profile_psy_data3 % PreStart('tra_asm_inc', 'r3', 0, 0)
        neuler = 0
        CALL profile_psy_data3 % PostEnd
        IF (ln_temnofreeze) THEN
          !$ACC KERNELS
          WHERE (t_bkginc(:, :, :) > 0.0_wp .OR. tsn(:, :, :, jp_tem) + t_bkginc(:, :, :) > fzptnz(:, :, :))
            tsn(:, :, :, jp_tem) = t_bkg(:, :, :) + t_bkginc(:, :, :)
          END WHERE
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          tsn(:, :, :, jp_tem) = t_bkg(:, :, :) + t_bkginc(:, :, :)
          !$ACC END KERNELS
        END IF
        IF (ln_salfix) THEN
          !$ACC KERNELS
          WHERE (s_bkginc(:, :, :) > 0.0_wp .OR. tsn(:, :, :, jp_sal) + s_bkginc(:, :, :) > salfixmin)
            tsn(:, :, :, jp_sal) = s_bkg(:, :, :) + s_bkginc(:, :, :)
          END WHERE
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          tsn(:, :, :, jp_sal) = s_bkg(:, :, :) + s_bkginc(:, :, :)
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        tsb(:, :, :, :) = tsn(:, :, :, :)
        !$ACC END KERNELS
        CALL profile_psy_data4 % PreStart('tra_asm_inc', 'r4', 0, 0)
        CALL eos(tsb, rhd, rhop, gdept_0(:, :, :))
        IF (ln_zps .AND. .NOT. lk_c1d .AND. .NOT. ln_isfcav) CALL zps_hde(kt, jpts, tsb, gtsu, gtsv, rhd, gru, grv)
        IF (ln_zps .AND. .NOT. lk_c1d .AND. ln_isfcav) CALL zps_hde_isf(nit000, jpts, tsb, gtsu, gtsv, gtui, gtvi, rhd, gru, grv, grui, grvi)
        DEALLOCATE(t_bkginc)
        DEALLOCATE(s_bkginc)
        DEALLOCATE(t_bkg)
        DEALLOCATE(s_bkg)
        CALL profile_psy_data4 % PostEnd
      END IF
    END IF
    CALL profile_psy_data5 % PreStart('tra_asm_inc', 'r5', 0, 0)
    IF (ln_seaiceinc) CALL seaice_asm_inc(kt)
    CALL profile_psy_data5 % PostEnd
  END SUBROUTINE tra_asm_inc
  SUBROUTINE dyn_asm_inc(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: jk
    INTEGER :: it
    REAL(KIND = wp) :: zincwgt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    IF (ln_asmiau) THEN
      IF ((kt >= nitiaustr_r) .AND. (kt <= nitiaufin_r)) THEN
        CALL profile_psy_data0 % PreStart('dyn_asm_inc', 'r0', 0, 0)
        it = kt - nit000 + 1
        zincwgt = wgtiau(it) / rdt
        IF (lwp) THEN
          WRITE(numout, *)
          WRITE(numout, *) 'dyn_asm_inc : Dynamics IAU at time step = ', kt, ' with IAU weight = ', wgtiau(it)
          WRITE(numout, *) '~~~~~~~~~~~~'
        END IF
        CALL profile_psy_data0 % PostEnd
        !$ACC KERNELS
        DO jk = 1, jpkm1
          ua(:, :, jk) = ua(:, :, jk) + u_bkginc(:, :, jk) * zincwgt
          va(:, :, jk) = va(:, :, jk) + v_bkginc(:, :, jk) * zincwgt
        END DO
        !$ACC END KERNELS
        CALL profile_psy_data1 % PreStart('dyn_asm_inc', 'r1', 0, 0)
        IF (kt == nitiaufin_r) THEN
          DEALLOCATE(u_bkginc)
          DEALLOCATE(v_bkginc)
        END IF
        CALL profile_psy_data1 % PostEnd
      END IF
    ELSE IF (ln_asmdin) THEN
      IF (kt == nitdin_r) THEN
        !$ACC KERNELS
        neuler = 0
        un(:, :, :) = u_bkg(:, :, :) + u_bkginc(:, :, :)
        vn(:, :, :) = v_bkg(:, :, :) + v_bkginc(:, :, :)
        ub(:, :, :) = un(:, :, :)
        vb(:, :, :) = vn(:, :, :)
        !$ACC END KERNELS
        CALL profile_psy_data2 % PreStart('dyn_asm_inc', 'r2', 0, 0)
        DEALLOCATE(u_bkg)
        DEALLOCATE(v_bkg)
        DEALLOCATE(u_bkginc)
        DEALLOCATE(v_bkginc)
        CALL profile_psy_data2 % PostEnd
      END IF
    END IF
  END SUBROUTINE dyn_asm_inc
  SUBROUTINE ssh_asm_inc(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: it
    INTEGER :: jk
    REAL(KIND = wp) :: zincwgt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    IF (ln_asmiau) THEN
      CALL profile_psy_data0 % PreStart('ssh_asm_inc', 'r0', 0, 0)
      IF ((kt >= nitiaustr_r) .AND. (kt <= nitiaufin_r)) THEN
        it = kt - nit000 + 1
        zincwgt = wgtiau(it) / rdt
        IF (lwp) THEN
          WRITE(numout, *)
          WRITE(numout, *) 'ssh_asm_inc : SSH IAU at time step = ', kt, ' with IAU weight = ', wgtiau(it)
          WRITE(numout, *) '~~~~~~~~~~~~'
        END IF
      ELSE IF (kt == nitiaufin_r + 1) THEN
        IF (ALLOCATED(ssh_bkginc)) DEALLOCATE(ssh_bkginc)
      END IF
      CALL profile_psy_data0 % PostEnd
    ELSE IF (ln_asmdin) THEN
      IF (kt == nitdin_r) THEN
        !$ACC KERNELS
        neuler = 0
        sshn(:, :) = ssh_bkg(:, :) + ssh_bkginc(:, :)
        sshb(:, :) = sshn(:, :)
        e3t_b(:, :, :) = e3t_n(:, :, :)
        !$ACC END KERNELS
        CALL profile_psy_data1 % PreStart('ssh_asm_inc', 'r1', 0, 0)
        DEALLOCATE(ssh_bkg)
        DEALLOCATE(ssh_bkginc)
        CALL profile_psy_data1 % PostEnd
      END IF
    END IF
  END SUBROUTINE ssh_asm_inc
  SUBROUTINE ssh_asm_div(kt, phdivn)
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: phdivn
    INTEGER :: jk
    REAL(KIND = wp), DIMENSION(:, :), POINTER :: ztim
  END SUBROUTINE ssh_asm_div
  SUBROUTINE seaice_asm_inc(kt, kindic)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN), OPTIONAL :: kindic
    INTEGER :: it
    REAL(KIND = wp) :: zincwgt
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zofrld, zohicif, zseaicendg, zhicifinc
    REAL(KIND = wp) :: zhicifmin = 0.5_wp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    IF (ln_asmiau) THEN
      IF ((kt >= nitiaustr_r) .AND. (kt <= nitiaufin_r)) THEN
        CALL profile_psy_data0 % PreStart('seaice_asm_inc', 'r0', 0, 0)
        it = kt - nit000 + 1
        zincwgt = wgtiau(it)
        IF (lwp) THEN
          WRITE(numout, *)
          WRITE(numout, *) 'seaice_asm_inc : sea ice conc IAU at time step = ', kt, ' with IAU weight = ', wgtiau(it)
          WRITE(numout, *) '~~~~~~~~~~~~'
        END IF
        CALL profile_psy_data0 % PostEnd
        !$ACC KERNELS
        zofrld(:, :) = 1._wp - at_i(:, :)
        zohicif(:, :) = hm_i(:, :)
        at_i(:, :) = 1. - MIN(MAX(1. - at_i(:, :) - seaice_bkginc(:, :) * zincwgt, 0.0_wp), 1.0_wp)
        at_i_b(:, :) = 1. - MIN(MAX(1. - at_i_b(:, :) - seaice_bkginc(:, :) * zincwgt, 0.0_wp), 1.0_wp)
        fr_i(:, :) = at_i(:, :)
        zseaicendg(:, :) = zofrld(:, :) - (1. - at_i(:, :))
        WHERE (zseaicendg(:, :) > 0.0_wp .AND. hm_i(:, :) < zhicifmin)
          zhicifinc(:, :) = (zhicifmin - hm_i(:, :)) * zincwgt
        ELSEWHERE
          zhicifinc(:, :) = 0.0_wp
        END WHERE
        hm_i(:, :) = hm_i(:, :) + zhicifinc(:, :)
        !$ACC END KERNELS
        CALL profile_psy_data1 % PreStart('seaice_asm_inc', 'r1', 0, 0)
        IF (kt == nitiaufin_r) THEN
          DEALLOCATE(seaice_bkginc)
        END IF
        CALL profile_psy_data1 % PostEnd
      ELSE
      END IF
    ELSE IF (ln_asmdin) THEN
      IF (kt == nitdin_r) THEN
        !$ACC KERNELS
        neuler = 0
        zofrld(:, :) = 1._wp - at_i(:, :)
        zohicif(:, :) = hm_i(:, :)
        at_i(:, :) = 1. - MIN(MAX(1. - at_i(:, :) - seaice_bkginc(:, :), 0.0_wp), 1.0_wp)
        at_i_b(:, :) = at_i(:, :)
        fr_i(:, :) = at_i(:, :)
        zseaicendg(:, :) = zofrld(:, :) - (1. - at_i(:, :))
        WHERE (zseaicendg(:, :) > 0.0_wp .AND. hm_i(:, :) < zhicifmin)
          zhicifinc(:, :) = zhicifmin - hm_i(:, :)
        ELSEWHERE
          zhicifinc(:, :) = 0.0_wp
        END WHERE
        hm_i(:, :) = hm_i(:, :) + zhicifinc(:, :)
        !$ACC END KERNELS
        CALL profile_psy_data2 % PreStart('seaice_asm_inc', 'r2', 0, 0)
        IF (.NOT. PRESENT(kindic)) THEN
          DEALLOCATE(seaice_bkginc)
        END IF
        CALL profile_psy_data2 % PostEnd
      ELSE
      END IF
    END IF
  END SUBROUTINE seaice_asm_inc
END MODULE asminc