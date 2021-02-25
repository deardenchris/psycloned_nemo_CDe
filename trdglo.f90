MODULE trdglo
  USE oce
  USE dom_oce
  USE sbc_oce
  USE trd_oce
  USE phycst
  USE ldftra
  USE ldfdyn
  USE zdf_oce
  USE zdfddm
  USE eosbn2
  USE phycst
  USE lib_mpp
  USE in_out_manager
  USE iom
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: trd_glo
  PUBLIC :: trd_glo_init
  REAL(KIND = wp) :: tvolt
  REAL(KIND = wp) :: tvolu
  REAL(KIND = wp) :: tvolv
  REAL(KIND = wp) :: rpktrd
  REAL(KIND = wp) :: peke
  REAL(KIND = wp), DIMENSION(jptot_tra) :: tmo, smo
  REAL(KIND = wp), DIMENSION(jptot_tra) :: t2, s2
  REAL(KIND = wp), DIMENSION(jptot_dyn) :: umo, vmo
  REAL(KIND = wp), DIMENSION(jptot_dyn) :: hke
  CONTAINS
  SUBROUTINE trd_glo(ptrdx, ptrdy, ktrd, ctype, kt)
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: ptrdx
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: ptrdy
    INTEGER, INTENT(IN) :: ktrd
    CHARACTER(LEN = 3), INTENT(IN) :: ctype
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: ikbu, ikbv
    REAL(KIND = wp) :: zvm, zvt, zvs, z1_2rau0
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztswu, ztswv, z2dx, z2dy
    IF (MOD(kt, nn_trd) == 0 .OR. kt == nit000 .OR. kt == nitend) THEN
      SELECT CASE (ctype)
      CASE ('TRA')
        ! !$OMP parallel default(shared), private(ji,jj,jk,zvm,zvs,zvt) ! CDe race condition??
        ! !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 1, jpj
            DO ji = 1, jpi
              zvm = e1e2t(ji, jj) * e3t_n(ji, jj, jk) * tmask(ji, jj, jk) * tmask_i(ji, jj)
              zvt = ptrdx(ji, jj, jk) * zvm
              zvs = ptrdy(ji, jj, jk) * zvm
              tmo(ktrd) = tmo(ktrd) + zvt
              smo(ktrd) = smo(ktrd) + zvs
              t2(ktrd) = t2(ktrd) + zvt * tsn(ji, jj, jk, jp_tem)
              s2(ktrd) = s2(ktrd) + zvs * tsn(ji, jj, jk, jp_sal)
            END DO
          END DO
        END DO
        ! !$OMP end do
        ! !$OMP end parallel
        IF (ln_linssh .AND. ktrd == jptra_zad) THEN
          tmo(jptra_sad) = SUM(wn(:, :, 1) * tsn(:, :, 1, jp_tem) * e1e2t(:, :) * tmask_i(:, :))
          smo(jptra_sad) = SUM(wn(:, :, 1) * tsn(:, :, 1, jp_sal) * e1e2t(:, :) * tmask_i(:, :))
          t2(jptra_sad) = SUM(wn(:, :, 1) * tsn(:, :, 1, jp_tem) * tsn(:, :, 1, jp_tem) * e1e2t(:, :) * tmask_i(:, :))
          s2(jptra_sad) = SUM(wn(:, :, 1) * tsn(:, :, 1, jp_sal) * tsn(:, :, 1, jp_sal) * e1e2t(:, :) * tmask_i(:, :))
        END IF
        IF (ktrd == jptra_atf) THEN
          CALL glo_tra_wri(kt)
          tmo(:) = 0._wp
          smo(:) = 0._wp
          t2(:) = 0._wp
          s2(:) = 0._wp
        END IF
      CASE ('DYN')
        ! !$OMP parallel default(shared), private(ji,jj,jk,zvs,zvt) ! CDe race condition??
        ! !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              zvt = ptrdx(ji, jj, jk) * tmask_i(ji + 1, jj) * tmask_i(ji, jj) * umask(ji, jj, jk) * e1e2u(ji, jj) * e3u_n(ji, jj, &
&jk)
              zvs = ptrdy(ji, jj, jk) * tmask_i(ji, jj + 1) * tmask_i(ji, jj) * vmask(ji, jj, jk) * e1e2v(ji, jj) * e3u_n(ji, jj, &
&jk)
              umo(ktrd) = umo(ktrd) + zvt
              vmo(ktrd) = vmo(ktrd) + zvs
              hke(ktrd) = hke(ktrd) + un(ji, jj, jk) * zvt + vn(ji, jj, jk) * zvs
            END DO
          END DO
        END DO
        ! !$OMP end do
        ! !$OMP end parallel
        IF (ktrd == jpdyn_zdf) THEN
          z1_2rau0 = 0.5_wp / rau0
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              zvt = (utau_b(ji, jj) + utau(ji, jj)) * tmask_i(ji + 1, jj) * tmask_i(ji, jj) * umask(ji, jj, jk) * z1_2rau0 * &
&e1e2u(ji, jj)
              zvs = (vtau_b(ji, jj) + vtau(ji, jj)) * tmask_i(ji, jj + 1) * tmask_i(ji, jj) * vmask(ji, jj, jk) * z1_2rau0 * &
&e1e2v(ji, jj)
              umo(jpdyn_tau) = umo(jpdyn_tau) + zvt
              vmo(jpdyn_tau) = vmo(jpdyn_tau) + zvs
              hke(jpdyn_tau) = hke(jpdyn_tau) + un(ji, jj, 1) * zvt + vn(ji, jj, 1) * zvs
            END DO
          END DO
        END IF
      END SELECT
    END IF
  END SUBROUTINE trd_glo
  SUBROUTINE glo_dyn_wri(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zcof
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zkx, zky, zkz, zkepe
    IF (MOD(kt, nn_trd) == 0 .OR. kt == nit000 .OR. kt == nitend) THEN
      zkx(:, :, :) = 0._wp
      zky(:, :, :) = 0._wp
      zkz(:, :, :) = 0._wp
      zkepe(:, :, :) = 0._wp
      CALL eos(tsn, rhd, rhop)
      zcof = 0.5_wp / rau0
      zkz(:, :, 1) = 0._wp
      DO jk = 2, jpk
        zkz(:, :, jk) = zcof * e1e2t(:, :) * wn(:, :, jk) * (rhop(:, :, jk) + rhop(:, :, jk - 1)) * tmask_i(:, :)
      END DO
      zcof = 0.5_wp / rau0
      DO jk = 1, jpkm1
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zkx(ji, jj, jk) = zcof * e2u(ji, jj) * e3u_n(ji, jj, jk) * un(ji, jj, jk) * (rhop(ji, jj, jk) + rhop(ji + 1, jj, jk))
            zky(ji, jj, jk) = zcof * e1v(ji, jj) * e3v_n(ji, jj, jk) * vn(ji, jj, jk) * (rhop(ji, jj, jk) + rhop(ji, jj + 1, jk))
          END DO
        END DO
      END DO
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zkepe(ji, jj, jk) = - (zkz(ji, jj, jk) - zkz(ji, jj, jk + 1) + zkx(ji, jj, jk) - zkx(ji - 1, jj, jk) + zky(ji, jj, jk) &
&- zky(ji, jj - 1, jk)) / (e1e2t(ji, jj) * e3t_n(ji, jj, jk)) * tmask(ji, jj, jk) * tmask_i(ji, jj)
          END DO
        END DO
      END DO
      peke = 0._wp
      DO jk = 1, jpkm1
        peke = peke + SUM(zkepe(:, :, jk) * gdept_n(:, :, jk) * e1e2t(:, :) * e3t_n(:, :, jk))
      END DO
      peke = grav * peke
      IF (lk_mpp) THEN
        CALL mpp_sum('trdglo', peke)
        CALL mpp_sum('trdglo', umo, jptot_dyn)
        CALL mpp_sum('trdglo', vmo, jptot_dyn)
        CALL mpp_sum('trdglo', hke, jptot_dyn)
      END IF
      IF (lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *)
        WRITE(numout, 9500) kt
        WRITE(numout, 9501) umo(jpdyn_hpg) / tvolu, vmo(jpdyn_hpg) / tvolv
        WRITE(numout, 9502) umo(jpdyn_keg) / tvolu, vmo(jpdyn_keg) / tvolv
        WRITE(numout, 9503) umo(jpdyn_rvo) / tvolu, vmo(jpdyn_rvo) / tvolv
        WRITE(numout, 9504) umo(jpdyn_pvo) / tvolu, vmo(jpdyn_pvo) / tvolv
        WRITE(numout, 9505) umo(jpdyn_zad) / tvolu, vmo(jpdyn_zad) / tvolv
        WRITE(numout, 9506) umo(jpdyn_ldf) / tvolu, vmo(jpdyn_ldf) / tvolv
        WRITE(numout, 9507) umo(jpdyn_zdf) / tvolu, vmo(jpdyn_zdf) / tvolv
        WRITE(numout, 9508) umo(jpdyn_spg) / tvolu, vmo(jpdyn_spg) / tvolv
        WRITE(numout, 9509) umo(jpdyn_bfr) / tvolu, vmo(jpdyn_bfr) / tvolv
        WRITE(numout, 9510) umo(jpdyn_atf) / tvolu, vmo(jpdyn_atf) / tvolv
        WRITE(numout, 9511)
        WRITE(numout, 9512) (umo(jpdyn_hpg) + umo(jpdyn_keg) + umo(jpdyn_rvo) + umo(jpdyn_pvo) + umo(jpdyn_zad) + umo(jpdyn_ldf) + &
&umo(jpdyn_zdf) + umo(jpdyn_spg) + umo(jpdyn_bfr) + umo(jpdyn_atf)) / tvolu, (vmo(jpdyn_hpg) + vmo(jpdyn_keg) + vmo(jpdyn_rvo) + &
&vmo(jpdyn_pvo) + vmo(jpdyn_zad) + vmo(jpdyn_ldf) + vmo(jpdyn_zdf) + vmo(jpdyn_spg) + vmo(jpdyn_bfr) + vmo(jpdyn_atf)) / tvolv
        WRITE(numout, 9513) umo(jpdyn_tau) / tvolu, vmo(jpdyn_tau) / tvolv
      END IF
9500  FORMAT(' momentum trend at it= ', I6, ' :', /, ' ==============================')
9501  FORMAT(' hydro pressure gradient    u= ', E20.13, '    v= ', E20.13)
9502  FORMAT(' ke gradient                u= ', E20.13, '    v= ', E20.13)
9503  FORMAT(' relative vorticity term    u= ', E20.13, '    v= ', E20.13)
9504  FORMAT(' planetary vorticity term   u= ', E20.13, '    v= ', E20.13)
9505  FORMAT(' vertical advection         u= ', E20.13, '    v= ', E20.13)
9506  FORMAT(' horizontal diffusion       u= ', E20.13, '    v= ', E20.13)
9507  FORMAT(' vertical diffusion         u= ', E20.13, '    v= ', E20.13)
9508  FORMAT(' surface pressure gradient  u= ', E20.13, '    v= ', E20.13)
9509  FORMAT(' explicit bottom friction   u= ', E20.13, '    v= ', E20.13)
9510  FORMAT(' Asselin time filter        u= ', E20.13, '    v= ', E20.13)
9511  FORMAT(' -----------------------------------------------------------------------------')
9512  FORMAT(' total trend                u= ', E20.13, '    v= ', E20.13)
9513  FORMAT(' incl. surface wind stress  u= ', E20.13, '    v= ', E20.13)
9514  FORMAT('       bottom stress        u= ', E20.13, '    v= ', E20.13)
      IF (lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *)
        WRITE(numout, 9520) kt
        WRITE(numout, 9521) hke(jpdyn_hpg) / tvolt
        WRITE(numout, 9522) hke(jpdyn_keg) / tvolt
        WRITE(numout, 9523) hke(jpdyn_rvo) / tvolt
        WRITE(numout, 9524) hke(jpdyn_pvo) / tvolt
        WRITE(numout, 9525) hke(jpdyn_zad) / tvolt
        WRITE(numout, 9526) hke(jpdyn_ldf) / tvolt
        WRITE(numout, 9527) hke(jpdyn_zdf) / tvolt
        WRITE(numout, 9528) hke(jpdyn_spg) / tvolt
        WRITE(numout, 9529) hke(jpdyn_bfr) / tvolt
        WRITE(numout, 9530) hke(jpdyn_atf) / tvolt
        WRITE(numout, 9531)
        WRITE(numout, 9532) (hke(jpdyn_hpg) + hke(jpdyn_keg) + hke(jpdyn_rvo) + hke(jpdyn_pvo) + hke(jpdyn_zad) + hke(jpdyn_ldf) + &
&hke(jpdyn_zdf) + hke(jpdyn_spg) + hke(jpdyn_bfr) + hke(jpdyn_atf)) / tvolt
        WRITE(numout, 9533) hke(jpdyn_tau) / tvolt
      END IF
9520  FORMAT(' kinetic energy trend at it= ', I6, ' :', /, ' ====================================')
9521  FORMAT(' hydro pressure gradient   u2= ', E20.13)
9522  FORMAT(' ke gradient               u2= ', E20.13)
9523  FORMAT(' relative vorticity term   u2= ', E20.13)
9524  FORMAT(' planetary vorticity term  u2= ', E20.13)
9525  FORMAT(' vertical advection        u2= ', E20.13)
9526  FORMAT(' horizontal diffusion      u2= ', E20.13)
9527  FORMAT(' vertical diffusion        u2= ', E20.13)
9528  FORMAT(' surface pressure gradient u2= ', E20.13)
9529  FORMAT(' explicit bottom friction  u2= ', E20.13)
9530  FORMAT(' Asselin time filter       u2= ', E20.13)
9531  FORMAT(' --------------------------------------------------')
9532  FORMAT(' total trend               u2= ', E20.13)
9533  FORMAT(' incl. surface wind stress u2= ', E20.13)
9534  FORMAT('       bottom stress       u2= ', E20.13)
      IF (lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *)
        WRITE(numout, 9540) kt
        WRITE(numout, 9541) (hke(jpdyn_keg) + hke(jpdyn_rvo) + hke(jpdyn_zad)) / tvolt
        WRITE(numout, 9542) (hke(jpdyn_keg) + hke(jpdyn_zad)) / tvolt
        WRITE(numout, 9543) (hke(jpdyn_pvo)) / tvolt
        WRITE(numout, 9544) (hke(jpdyn_rvo)) / tvolt
        WRITE(numout, 9545) (hke(jpdyn_spg)) / tvolt
        WRITE(numout, 9546) (hke(jpdyn_ldf)) / tvolt
        WRITE(numout, 9547) (hke(jpdyn_zdf)) / tvolt
        WRITE(numout, 9548) (hke(jpdyn_hpg)) / tvolt, rpktrd / tvolt
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *)
      END IF
9540  FORMAT(' energetic consistency at it= ', I6, ' :', /, ' =========================================')
9541  FORMAT(' 0 = non linear term (true if KE conserved)                : ', E20.13)
9542  FORMAT(' 0 = ke gradient + vertical advection                      : ', E20.13)
9543  FORMAT(' 0 = coriolis term  (true if KE conserving scheme)         : ', E20.13)
9544  FORMAT(' 0 = vorticity term (true if KE conserving scheme)         : ', E20.13)
9545  FORMAT(' 0 = surface pressure gradient  ???                        : ', E20.13)
9546  FORMAT(' 0 < horizontal diffusion                                  : ', E20.13)
9547  FORMAT(' 0 < vertical diffusion                                    : ', E20.13)
9548  FORMAT(' pressure gradient u2 = - 1/rau0 u.dz(rhop)                : ', E20.13, '  u.dz(rhop) =', E20.13)
      rpktrd = peke
    END IF
  END SUBROUTINE glo_dyn_wri
  SUBROUTINE glo_tra_wri(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: jk
    IF (MOD(kt, nn_trd) == 0 .OR. kt == nit000 .OR. kt == nitend) THEN
      IF (lk_mpp) THEN
        CALL mpp_sum('trdglo', tmo, jptot_tra)
        CALL mpp_sum('trdglo', smo, jptot_tra)
        CALL mpp_sum('trdglo', t2, jptot_tra)
        CALL mpp_sum('trdglo', s2, jptot_tra)
      END IF
      IF (lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *)
        WRITE(numout, 9400) kt
        WRITE(numout, 9401) tmo(jptra_xad) / tvolt, smo(jptra_xad) / tvolt
        WRITE(numout, 9411) tmo(jptra_yad) / tvolt, smo(jptra_yad) / tvolt
        WRITE(numout, 9402) tmo(jptra_zad) / tvolt, smo(jptra_zad) / tvolt
        WRITE(numout, 9403) tmo(jptra_ldf) / tvolt, smo(jptra_ldf) / tvolt
        WRITE(numout, 9404) tmo(jptra_zdf) / tvolt, smo(jptra_zdf) / tvolt
        WRITE(numout, 9405) tmo(jptra_npc) / tvolt, smo(jptra_npc) / tvolt
        WRITE(numout, 9406) tmo(jptra_dmp) / tvolt, smo(jptra_dmp) / tvolt
        WRITE(numout, 9407) tmo(jptra_qsr) / tvolt
        WRITE(numout, 9408) tmo(jptra_nsr) / tvolt, smo(jptra_nsr) / tvolt
        WRITE(numout, 9409)
        WRITE(numout, 9410) (tmo(jptra_xad) + tmo(jptra_yad) + tmo(jptra_zad) + tmo(jptra_ldf) + tmo(jptra_zdf) + tmo(jptra_npc) + &
&tmo(jptra_dmp) + tmo(jptra_qsr) + tmo(jptra_nsr)) / tvolt, (smo(jptra_xad) + smo(jptra_yad) + smo(jptra_zad) + smo(jptra_ldf) + &
&smo(jptra_zdf) + smo(jptra_npc) + smo(jptra_dmp) + smo(jptra_nsr)) / tvolt
      END IF
9400  FORMAT(' tracer trend at it= ', I6, ' :     temperature', '              salinity', /, ' ============================')
9401  FORMAT(' zonal      advection        ', E20.13, '     ', E20.13)
9411  FORMAT(' meridional advection        ', E20.13, '     ', E20.13)
9402  FORMAT(' vertical advection          ', E20.13, '     ', E20.13)
9403  FORMAT(' horizontal diffusion        ', E20.13, '     ', E20.13)
9404  FORMAT(' vertical diffusion          ', E20.13, '     ', E20.13)
9405  FORMAT(' static instability mixing   ', E20.13, '     ', E20.13)
9406  FORMAT(' damping term                ', E20.13, '     ', E20.13)
9407  FORMAT(' penetrative qsr             ', E20.13)
9408  FORMAT(' non solar radiation         ', E20.13, '     ', E20.13)
9409  FORMAT(' -------------------------------------------------------------------------')
9410  FORMAT(' total trend                 ', E20.13, '     ', E20.13)
      IF (lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *)
        WRITE(numout, 9420) kt
        WRITE(numout, 9421) t2(jptra_xad) / tvolt, s2(jptra_xad) / tvolt
        WRITE(numout, 9431) t2(jptra_yad) / tvolt, s2(jptra_yad) / tvolt
        WRITE(numout, 9422) t2(jptra_zad) / tvolt, s2(jptra_zad) / tvolt
        WRITE(numout, 9423) t2(jptra_ldf) / tvolt, s2(jptra_ldf) / tvolt
        WRITE(numout, 9424) t2(jptra_zdf) / tvolt, s2(jptra_zdf) / tvolt
        WRITE(numout, 9425) t2(jptra_npc) / tvolt, s2(jptra_npc) / tvolt
        WRITE(numout, 9426) t2(jptra_dmp) / tvolt, s2(jptra_dmp) / tvolt
        WRITE(numout, 9427) t2(jptra_qsr) / tvolt
        WRITE(numout, 9428) t2(jptra_nsr) / tvolt, s2(jptra_nsr) / tvolt
        WRITE(numout, 9429)
        WRITE(numout, 9430) (t2(jptra_xad) + t2(jptra_yad) + t2(jptra_zad) + t2(jptra_ldf) + t2(jptra_zdf) + t2(jptra_npc) + &
&t2(jptra_dmp) + t2(jptra_qsr) + t2(jptra_nsr)) / tvolt, (s2(jptra_xad) + s2(jptra_yad) + s2(jptra_zad) + s2(jptra_ldf) + &
&s2(jptra_zdf) + s2(jptra_npc) + s2(jptra_dmp) + s2(jptra_nsr)) / tvolt
      END IF
9420  FORMAT(' tracer**2 trend at it= ', I6, ' :      temperature', '               salinity', /, ' &
&===============================')
9421  FORMAT(' zonal      advection      * t   ', E20.13, '     ', E20.13)
9431  FORMAT(' meridional advection      * t   ', E20.13, '     ', E20.13)
9422  FORMAT(' vertical advection        * t   ', E20.13, '     ', E20.13)
9423  FORMAT(' horizontal diffusion      * t   ', E20.13, '     ', E20.13)
9424  FORMAT(' vertical diffusion        * t   ', E20.13, '     ', E20.13)
9425  FORMAT(' static instability mixing * t   ', E20.13, '     ', E20.13)
9426  FORMAT(' damping term              * t   ', E20.13, '     ', E20.13)
9427  FORMAT(' penetrative qsr           * t   ', E20.13)
9428  FORMAT(' non solar radiation       * t   ', E20.13, '     ', E20.13)
9429  FORMAT(' -----------------------------------------------------------------------------')
9430  FORMAT(' total trend                *t = ', E20.13, '  *s = ', E20.13)
      IF (lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *)
        WRITE(numout, 9440) kt
        WRITE(numout, 9441) (tmo(jptra_xad) + tmo(jptra_yad) + tmo(jptra_zad)) / tvolt, (smo(jptra_xad) + smo(jptra_yad) + &
&smo(jptra_zad)) / tvolt
        WRITE(numout, 9442) tmo(jptra_sad) / tvolt, smo(jptra_sad) / tvolt
        WRITE(numout, 9443) tmo(jptra_ldf) / tvolt, smo(jptra_ldf) / tvolt
        WRITE(numout, 9444) tmo(jptra_zdf) / tvolt, smo(jptra_zdf) / tvolt
        WRITE(numout, 9445) tmo(jptra_npc) / tvolt, smo(jptra_npc) / tvolt
        WRITE(numout, 9446) (t2(jptra_xad) + t2(jptra_yad) + t2(jptra_zad)) / tvolt, (s2(jptra_xad) + s2(jptra_yad) + &
&s2(jptra_zad)) / tvolt
        WRITE(numout, 9447) t2(jptra_ldf) / tvolt, s2(jptra_ldf) / tvolt
        WRITE(numout, 9448) t2(jptra_zdf) / tvolt, s2(jptra_zdf) / tvolt
        WRITE(numout, 9449) t2(jptra_npc) / tvolt, s2(jptra_npc) / tvolt
      END IF
9440  FORMAT(' tracer consistency at it= ', I6, ' :         temperature', '                salinity', /, ' &
&==================================')
9441  FORMAT(' 0 = horizontal+vertical advection +    ', E20.13, '       ', E20.13)
9442  FORMAT('     1st lev vertical advection         ', E20.13, '       ', E20.13)
9443  FORMAT(' 0 = horizontal diffusion               ', E20.13, '       ', E20.13)
9444  FORMAT(' 0 = vertical diffusion                 ', E20.13, '       ', E20.13)
9445  FORMAT(' 0 = static instability mixing          ', E20.13, '       ', E20.13)
9446  FORMAT(' 0 = horizontal+vertical advection * t  ', E20.13, '       ', E20.13)
9447  FORMAT(' 0 > horizontal diffusion          * t  ', E20.13, '       ', E20.13)
9448  FORMAT(' 0 > vertical diffusion            * t  ', E20.13, '       ', E20.13)
9449  FORMAT(' 0 > static instability mixing     * t  ', E20.13, '       ', E20.13)
    END IF
  END SUBROUTINE glo_tra_wri
  SUBROUTINE trd_glo_init
    INTEGER :: ji, jj, jk
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'trd_glo_init : integral constraints properties trends'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~~'
    END IF
    tvolt = 0._wp
    DO jk = 1, jpkm1
      tvolt = tvolt + SUM(e1e2t(:, :) * e3t_n(:, :, jk) * tmask(:, :, jk) * tmask_i(:, :))
    END DO
    CALL mpp_sum('trdglo', tvolt)
    IF (lwp) WRITE(numout, FMT = *) '                total ocean volume at T-point   tvolt = ', tvolt
    rpktrd = 0._wp
    tvolu = 0._wp
    tvolv = 0._wp
    DO jk = 1, jpk
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          tvolu = tvolu + e1u(ji, jj) * e2u(ji, jj) * e3u_n(ji, jj, jk) * tmask_i(ji + 1, jj) * tmask_i(ji, jj) * umask(ji, jj, jk)
          tvolv = tvolv + e1v(ji, jj) * e2v(ji, jj) * e3v_n(ji, jj, jk) * tmask_i(ji, jj + 1) * tmask_i(ji, jj) * vmask(ji, jj, jk)
        END DO
      END DO
    END DO
    CALL mpp_sum('trdglo', tvolu)
    CALL mpp_sum('trdglo', tvolv)
    IF (lwp) THEN
      WRITE(numout, FMT = *) '                total ocean volume at U-point   tvolu = ', tvolu
      WRITE(numout, FMT = *) '                total ocean volume at V-point   tvolv = ', tvolv
    END IF
  END SUBROUTINE trd_glo_init
END MODULE trdglo
