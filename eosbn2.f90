MODULE eosbn2
  USE dom_oce
  USE phycst
  USE stopar
  USE stopts
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  USE prtctl
  USE lbclnk
  USE timing
  IMPLICIT NONE
  PRIVATE
  INTERFACE eos
    MODULE PROCEDURE eos_insitu, eos_insitu_pot, eos_insitu_2d
  END INTERFACE
  INTERFACE eos_rab
    MODULE PROCEDURE rab_3d, rab_2d, rab_0d
  END INTERFACE
  INTERFACE eos_fzp
    MODULE PROCEDURE eos_fzp_2d, eos_fzp_0d
  END INTERFACE
  PUBLIC :: eos
  PUBLIC :: bn2
  PUBLIC :: eos_rab
  PUBLIC :: eos_pt_from_ct
  PUBLIC :: eos_fzp
  PUBLIC :: eos_pen
  PUBLIC :: eos_init
  LOGICAL, PUBLIC :: ln_TEOS10
  LOGICAL, PUBLIC :: ln_EOS80
  LOGICAL, PUBLIC :: ln_SEOS
  LOGICAL, PUBLIC :: l_useCT
  INTEGER, PUBLIC :: neos
  INTEGER, PARAMETER :: np_teos10 = - 1
  INTEGER, PARAMETER :: np_eos80 = 0
  INTEGER, PARAMETER :: np_seos = 1
  REAL(KIND = wp) :: rn_a0 = 1.6550E-1_wp
  REAL(KIND = wp) :: rn_b0 = 7.6554E-1_wp
  REAL(KIND = wp) :: rn_lambda1 = 5.9520E-2_wp
  REAL(KIND = wp) :: rn_lambda2 = 5.4914E-4_wp
  REAL(KIND = wp) :: rn_mu1 = 1.4970E-4_wp
  REAL(KIND = wp) :: rn_mu2 = 1.1090E-5_wp
  REAL(KIND = wp) :: rn_nu = 2.4341E-3_wp
  REAL(KIND = wp) :: r1_S0, r1_T0, r1_Z0, rdeltaS
  REAL(KIND = wp) :: EOS000, EOS100, EOS200, EOS300, EOS400, EOS500, EOS600
  REAL(KIND = wp) :: EOS010, EOS110, EOS210, EOS310, EOS410, EOS510
  REAL(KIND = wp) :: EOS020, EOS120, EOS220, EOS320, EOS420
  REAL(KIND = wp) :: EOS030, EOS130, EOS230, EOS330
  REAL(KIND = wp) :: EOS040, EOS140, EOS240
  REAL(KIND = wp) :: EOS050, EOS150
  REAL(KIND = wp) :: EOS060
  REAL(KIND = wp) :: EOS001, EOS101, EOS201, EOS301, EOS401
  REAL(KIND = wp) :: EOS011, EOS111, EOS211, EOS311
  REAL(KIND = wp) :: EOS021, EOS121, EOS221
  REAL(KIND = wp) :: EOS031, EOS131
  REAL(KIND = wp) :: EOS041
  REAL(KIND = wp) :: EOS002, EOS102, EOS202
  REAL(KIND = wp) :: EOS012, EOS112
  REAL(KIND = wp) :: EOS022
  REAL(KIND = wp) :: EOS003, EOS103
  REAL(KIND = wp) :: EOS013
  REAL(KIND = wp) :: ALP000, ALP100, ALP200, ALP300, ALP400, ALP500
  REAL(KIND = wp) :: ALP010, ALP110, ALP210, ALP310, ALP410
  REAL(KIND = wp) :: ALP020, ALP120, ALP220, ALP320
  REAL(KIND = wp) :: ALP030, ALP130, ALP230
  REAL(KIND = wp) :: ALP040, ALP140
  REAL(KIND = wp) :: ALP050
  REAL(KIND = wp) :: ALP001, ALP101, ALP201, ALP301
  REAL(KIND = wp) :: ALP011, ALP111, ALP211
  REAL(KIND = wp) :: ALP021, ALP121
  REAL(KIND = wp) :: ALP031
  REAL(KIND = wp) :: ALP002, ALP102
  REAL(KIND = wp) :: ALP012
  REAL(KIND = wp) :: ALP003
  REAL(KIND = wp) :: BET000, BET100, BET200, BET300, BET400, BET500
  REAL(KIND = wp) :: BET010, BET110, BET210, BET310, BET410
  REAL(KIND = wp) :: BET020, BET120, BET220, BET320
  REAL(KIND = wp) :: BET030, BET130, BET230
  REAL(KIND = wp) :: BET040, BET140
  REAL(KIND = wp) :: BET050
  REAL(KIND = wp) :: BET001, BET101, BET201, BET301
  REAL(KIND = wp) :: BET011, BET111, BET211
  REAL(KIND = wp) :: BET021, BET121
  REAL(KIND = wp) :: BET031
  REAL(KIND = wp) :: BET002, BET102
  REAL(KIND = wp) :: BET012
  REAL(KIND = wp) :: BET003
  REAL(KIND = wp) :: PEN000, PEN100, PEN200, PEN300, PEN400
  REAL(KIND = wp) :: PEN010, PEN110, PEN210, PEN310
  REAL(KIND = wp) :: PEN020, PEN120, PEN220
  REAL(KIND = wp) :: PEN030, PEN130
  REAL(KIND = wp) :: PEN040
  REAL(KIND = wp) :: PEN001, PEN101, PEN201
  REAL(KIND = wp) :: PEN011, PEN111
  REAL(KIND = wp) :: PEN021
  REAL(KIND = wp) :: PEN002, PEN102
  REAL(KIND = wp) :: PEN012
  REAL(KIND = wp) :: APE000, APE100, APE200, APE300
  REAL(KIND = wp) :: APE010, APE110, APE210
  REAL(KIND = wp) :: APE020, APE120
  REAL(KIND = wp) :: APE030
  REAL(KIND = wp) :: APE001, APE101
  REAL(KIND = wp) :: APE011
  REAL(KIND = wp) :: APE002
  REAL(KIND = wp) :: BPE000, BPE100, BPE200, BPE300
  REAL(KIND = wp) :: BPE010, BPE110, BPE210
  REAL(KIND = wp) :: BPE020, BPE120
  REAL(KIND = wp) :: BPE030
  REAL(KIND = wp) :: BPE001, BPE101
  REAL(KIND = wp) :: BPE011
  REAL(KIND = wp) :: BPE002
  CONTAINS
  SUBROUTINE eos_insitu(pts, prd, pdep)
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, jpts), INTENT(IN) :: pts
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(OUT) :: prd
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pdep
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zt, zh, zs, ztm
    REAL(KIND = wp) :: zn, zn0, zn1, zn2, zn3
    IF (ln_timing) CALL timing_start('eos-insitu')
    SELECT CASE (neos)
    CASE (np_teos10, np_eos80)
      !$OMP parallel default(shared), private(ji,jj,jk,zh,zn,zn0,zn1,zn2,zn3,zs,zt,ztm)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            zh = pdep(ji, jj, jk) * r1_Z0
            zt = pts(ji, jj, jk, jp_tem) * r1_T0
            zs = SQRT(ABS(pts(ji, jj, jk, jp_sal) + rdeltaS) * r1_S0)
            ztm = tmask(ji, jj, jk)
            zn3 = EOS013 * zt + EOS103 * zs + EOS003
            zn2 = (EOS022 * zt + EOS112 * zs + EOS012) * zt + (EOS202 * zs + EOS102) * zs + EOS002
            zn1 = (((EOS041 * zt + EOS131 * zs + EOS031) * zt + (EOS221 * zs + EOS121) * zs + EOS021) * zt + ((EOS311 * zs + &
&EOS211) * zs + EOS111) * zs + EOS011) * zt + (((EOS401 * zs + EOS301) * zs + EOS201) * zs + EOS101) * zs + EOS001
            zn0 = (((((EOS060 * zt + EOS150 * zs + EOS050) * zt + (EOS240 * zs + EOS140) * zs + EOS040) * zt + ((EOS330 * zs + &
&EOS230) * zs + EOS130) * zs + EOS030) * zt + (((EOS420 * zs + EOS320) * zs + EOS220) * zs + EOS120) * zs + EOS020) * zt + &
&((((EOS510 * zs + EOS410) * zs + EOS310) * zs + EOS210) * zs + EOS110) * zs + EOS010) * zt + (((((EOS600 * zs + EOS500) * zs + &
&EOS400) * zs + EOS300) * zs + EOS200) * zs + EOS100) * zs + EOS000
            zn = ((zn3 * zh + zn2) * zh + zn1) * zh + zn0
            prd(ji, jj, jk) = (zn * r1_rau0 - 1._wp) * ztm
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
    CASE (np_seos)
      !$OMP parallel default(shared), private(ji,jj,jk,zh,zn,zs,zt,ztm)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            zt = pts(ji, jj, jk, jp_tem) - 10._wp
            zs = pts(ji, jj, jk, jp_sal) - 35._wp
            zh = pdep(ji, jj, jk)
            ztm = tmask(ji, jj, jk)
            zn = - rn_a0 * (1._wp + 0.5_wp * rn_lambda1 * zt + rn_mu1 * zh) * zt + rn_b0 * (1._wp - 0.5_wp * rn_lambda2 * zs - &
&rn_mu2 * zh) * zs - rn_nu * zt * zs
            prd(ji, jj, jk) = zn * r1_rau0 * ztm
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
    END SELECT
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = prd, clinfo1 = ' eos-insitu  : ', kdim = jpk)
    IF (ln_timing) CALL timing_stop('eos-insitu')
  END SUBROUTINE eos_insitu
  SUBROUTINE eos_insitu_pot(pts, prd, prhop, pdep)
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, jpts), INTENT(IN) :: pts
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(OUT) :: prd
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(OUT) :: prhop
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pdep
    INTEGER :: ji, jj, jk, jsmp
    INTEGER :: jdof
    REAL(KIND = wp) :: zt, zh, zstemp, zs, ztm
    REAL(KIND = wp) :: zn, zn0, zn1, zn2, zn3
    REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: zn0_sto, zn_sto, zsign
    IF (ln_timing) CALL timing_start('eos-pot')
    SELECT CASE (neos)
    CASE (np_teos10, np_eos80)
      IF (ln_sto_eos) THEN
        ALLOCATE(zn0_sto(1 : 2 * nn_sto_eos))
        ALLOCATE(zn_sto(1 : 2 * nn_sto_eos))
        ALLOCATE(zsign(1 : 2 * nn_sto_eos))
        DO jsmp = 1, 2 * nn_sto_eos, 2
          zsign(jsmp) = 1._wp
          zsign(jsmp + 1) = - 1._wp
        END DO
        DO jk = 1, jpkm1
          DO jj = 1, jpj
            DO ji = 1, jpi
              DO jsmp = 1, nn_sto_eos * 2
                jdof = (jsmp + 1) / 2
                zh = pdep(ji, jj, jk) * r1_Z0
                zt = (pts(ji, jj, jk, jp_tem) + pts_ran(ji, jj, jk, jp_tem, jdof) * zsign(jsmp)) * r1_T0
                zstemp = pts(ji, jj, jk, jp_sal) + pts_ran(ji, jj, jk, jp_sal, jdof) * zsign(jsmp)
                zs = SQRT(ABS(zstemp + rdeltaS) * r1_S0)
                ztm = tmask(ji, jj, jk)
                zn3 = EOS013 * zt + EOS103 * zs + EOS003
                zn2 = (EOS022 * zt + EOS112 * zs + EOS012) * zt + (EOS202 * zs + EOS102) * zs + EOS002
                zn1 = (((EOS041 * zt + EOS131 * zs + EOS031) * zt + (EOS221 * zs + EOS121) * zs + EOS021) * zt + ((EOS311 * zs + &
&EOS211) * zs + EOS111) * zs + EOS011) * zt + (((EOS401 * zs + EOS301) * zs + EOS201) * zs + EOS101) * zs + EOS001
                zn0_sto(jsmp) = (((((EOS060 * zt + EOS150 * zs + EOS050) * zt + (EOS240 * zs + EOS140) * zs + EOS040) * zt + &
&((EOS330 * zs + EOS230) * zs + EOS130) * zs + EOS030) * zt + (((EOS420 * zs + EOS320) * zs + EOS220) * zs + EOS120) * zs + &
&EOS020) * zt + ((((EOS510 * zs + EOS410) * zs + EOS310) * zs + EOS210) * zs + EOS110) * zs + EOS010) * zt + (((((EOS600 * zs + &
&EOS500) * zs + EOS400) * zs + EOS300) * zs + EOS200) * zs + EOS100) * zs + EOS000
                zn_sto(jsmp) = ((zn3 * zh + zn2) * zh + zn1) * zh + zn0_sto(jsmp)
              END DO
              prhop(ji, jj, jk) = 0._wp
              prd(ji, jj, jk) = 0._wp
              DO jsmp = 1, nn_sto_eos * 2
                prhop(ji, jj, jk) = prhop(ji, jj, jk) + zn0_sto(jsmp)
                prd(ji, jj, jk) = prd(ji, jj, jk) + (zn_sto(jsmp) * r1_rau0 - 1._wp)
              END DO
              prhop(ji, jj, jk) = 0.5_wp * prhop(ji, jj, jk) * ztm / nn_sto_eos
              prd(ji, jj, jk) = 0.5_wp * prd(ji, jj, jk) * ztm / nn_sto_eos
            END DO
          END DO
        END DO
        DEALLOCATE(zn0_sto, zn_sto, zsign)
      ELSE
        DO jk = 1, jpkm1
          DO jj = 1, jpj
            DO ji = 1, jpi
              zh = pdep(ji, jj, jk) * r1_Z0
              zt = pts(ji, jj, jk, jp_tem) * r1_T0
              zs = SQRT(ABS(pts(ji, jj, jk, jp_sal) + rdeltaS) * r1_S0)
              ztm = tmask(ji, jj, jk)
              zn3 = EOS013 * zt + EOS103 * zs + EOS003
              zn2 = (EOS022 * zt + EOS112 * zs + EOS012) * zt + (EOS202 * zs + EOS102) * zs + EOS002
              zn1 = (((EOS041 * zt + EOS131 * zs + EOS031) * zt + (EOS221 * zs + EOS121) * zs + EOS021) * zt + ((EOS311 * zs + &
&EOS211) * zs + EOS111) * zs + EOS011) * zt + (((EOS401 * zs + EOS301) * zs + EOS201) * zs + EOS101) * zs + EOS001
              zn0 = (((((EOS060 * zt + EOS150 * zs + EOS050) * zt + (EOS240 * zs + EOS140) * zs + EOS040) * zt + ((EOS330 * zs + &
&EOS230) * zs + EOS130) * zs + EOS030) * zt + (((EOS420 * zs + EOS320) * zs + EOS220) * zs + EOS120) * zs + EOS020) * zt + &
&((((EOS510 * zs + EOS410) * zs + EOS310) * zs + EOS210) * zs + EOS110) * zs + EOS010) * zt + (((((EOS600 * zs + EOS500) * zs + &
&EOS400) * zs + EOS300) * zs + EOS200) * zs + EOS100) * zs + EOS000
              zn = ((zn3 * zh + zn2) * zh + zn1) * zh + zn0
              prhop(ji, jj, jk) = zn0 * ztm
              prd(ji, jj, jk) = (zn * r1_rau0 - 1._wp) * ztm
            END DO
          END DO
        END DO
      END IF
    CASE (np_seos)
      DO jk = 1, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            zt = pts(ji, jj, jk, jp_tem) - 10._wp
            zs = pts(ji, jj, jk, jp_sal) - 35._wp
            zh = pdep(ji, jj, jk)
            ztm = tmask(ji, jj, jk)
            zn = - rn_a0 * (1._wp + 0.5_wp * rn_lambda1 * zt) * zt + rn_b0 * (1._wp - 0.5_wp * rn_lambda2 * zs) * zs - rn_nu * zt &
&* zs
            prhop(ji, jj, jk) = (rau0 + zn) * ztm
            zn = zn - (rn_a0 * rn_mu1 * zt + rn_b0 * rn_mu2 * zs) * zh
            prd(ji, jj, jk) = zn * r1_rau0 * ztm
          END DO
        END DO
      END DO
    END SELECT
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = prd, clinfo1 = ' eos-pot: ', tab3d_2 = prhop, clinfo2 = ' pot : ', kdim = jpk)
    IF (ln_timing) CALL timing_stop('eos-pot')
  END SUBROUTINE eos_insitu_pot
  SUBROUTINE eos_insitu_2d(pts, pdep, prd)
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpts), INTENT(IN) :: pts
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pdep
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(OUT) :: prd
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zt, zh, zs
    REAL(KIND = wp) :: zn, zn0, zn1, zn2, zn3
    IF (ln_timing) CALL timing_start('eos2d')
    prd(:, :) = 0._wp
    SELECT CASE (neos)
    CASE (np_teos10, np_eos80)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zh = pdep(ji, jj) * r1_Z0
          zt = pts(ji, jj, jp_tem) * r1_T0
          zs = SQRT(ABS(pts(ji, jj, jp_sal) + rdeltaS) * r1_S0)
          zn3 = EOS013 * zt + EOS103 * zs + EOS003
          zn2 = (EOS022 * zt + EOS112 * zs + EOS012) * zt + (EOS202 * zs + EOS102) * zs + EOS002
          zn1 = (((EOS041 * zt + EOS131 * zs + EOS031) * zt + (EOS221 * zs + EOS121) * zs + EOS021) * zt + ((EOS311 * zs + EOS211) &
&* zs + EOS111) * zs + EOS011) * zt + (((EOS401 * zs + EOS301) * zs + EOS201) * zs + EOS101) * zs + EOS001
          zn0 = (((((EOS060 * zt + EOS150 * zs + EOS050) * zt + (EOS240 * zs + EOS140) * zs + EOS040) * zt + ((EOS330 * zs + &
&EOS230) * zs + EOS130) * zs + EOS030) * zt + (((EOS420 * zs + EOS320) * zs + EOS220) * zs + EOS120) * zs + EOS020) * zt + &
&((((EOS510 * zs + EOS410) * zs + EOS310) * zs + EOS210) * zs + EOS110) * zs + EOS010) * zt + (((((EOS600 * zs + EOS500) * zs + &
&EOS400) * zs + EOS300) * zs + EOS200) * zs + EOS100) * zs + EOS000
          zn = ((zn3 * zh + zn2) * zh + zn1) * zh + zn0
          prd(ji, jj) = zn * r1_rau0 - 1._wp
        END DO
      END DO
      CALL lbc_lnk('eosbn2', prd, 'T', 1.)
    CASE (np_seos)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zt = pts(ji, jj, jp_tem) - 10._wp
          zs = pts(ji, jj, jp_sal) - 35._wp
          zh = pdep(ji, jj)
          zn = - rn_a0 * (1._wp + 0.5_wp * rn_lambda1 * zt + rn_mu1 * zh) * zt + rn_b0 * (1._wp - 0.5_wp * rn_lambda2 * zs - &
&rn_mu2 * zh) * zs - rn_nu * zt * zs
          prd(ji, jj) = zn * r1_rau0
        END DO
      END DO
      CALL lbc_lnk('eosbn2', prd, 'T', 1.)
    END SELECT
    IF (ln_ctl) CALL prt_ctl(tab2d_1 = prd, clinfo1 = ' eos2d: ')
    IF (ln_timing) CALL timing_stop('eos2d')
  END SUBROUTINE eos_insitu_2d
  SUBROUTINE rab_3d(pts, pab)
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, jpts), INTENT(IN) :: pts
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, jpts), INTENT(OUT) :: pab
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zt, zh, zs, ztm
    REAL(KIND = wp) :: zn, zn0, zn1, zn2, zn3
    IF (ln_timing) CALL timing_start('rab_3d')
    SELECT CASE (neos)
    CASE (np_teos10, np_eos80)
      DO jk = 1, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            zh = gdept_n(ji, jj, jk) * r1_Z0
            zt = pts(ji, jj, jk, jp_tem) * r1_T0
            zs = SQRT(ABS(pts(ji, jj, jk, jp_sal) + rdeltaS) * r1_S0)
            ztm = tmask(ji, jj, jk)
            zn3 = ALP003
            zn2 = ALP012 * zt + ALP102 * zs + ALP002
            zn1 = ((ALP031 * zt + ALP121 * zs + ALP021) * zt + (ALP211 * zs + ALP111) * zs + ALP011) * zt + ((ALP301 * zs + &
&ALP201) * zs + ALP101) * zs + ALP001
            zn0 = ((((ALP050 * zt + ALP140 * zs + ALP040) * zt + (ALP230 * zs + ALP130) * zs + ALP030) * zt + ((ALP320 * zs + &
&ALP220) * zs + ALP120) * zs + ALP020) * zt + (((ALP410 * zs + ALP310) * zs + ALP210) * zs + ALP110) * zs + ALP010) * zt + &
&((((ALP500 * zs + ALP400) * zs + ALP300) * zs + ALP200) * zs + ALP100) * zs + ALP000
            zn = ((zn3 * zh + zn2) * zh + zn1) * zh + zn0
            pab(ji, jj, jk, jp_tem) = zn * r1_rau0 * ztm
            zn3 = BET003
            zn2 = BET012 * zt + BET102 * zs + BET002
            zn1 = ((BET031 * zt + BET121 * zs + BET021) * zt + (BET211 * zs + BET111) * zs + BET011) * zt + ((BET301 * zs + &
&BET201) * zs + BET101) * zs + BET001
            zn0 = ((((BET050 * zt + BET140 * zs + BET040) * zt + (BET230 * zs + BET130) * zs + BET030) * zt + ((BET320 * zs + &
&BET220) * zs + BET120) * zs + BET020) * zt + (((BET410 * zs + BET310) * zs + BET210) * zs + BET110) * zs + BET010) * zt + &
&((((BET500 * zs + BET400) * zs + BET300) * zs + BET200) * zs + BET100) * zs + BET000
            zn = ((zn3 * zh + zn2) * zh + zn1) * zh + zn0
            pab(ji, jj, jk, jp_sal) = zn / zs * r1_rau0 * ztm
          END DO
        END DO
      END DO
    CASE (np_seos)
      DO jk = 1, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            zt = pts(ji, jj, jk, jp_tem) - 10._wp
            zs = pts(ji, jj, jk, jp_sal) - 35._wp
            zh = gdept_n(ji, jj, jk)
            ztm = tmask(ji, jj, jk)
            zn = rn_a0 * (1._wp + rn_lambda1 * zt + rn_mu1 * zh) + rn_nu * zs
            pab(ji, jj, jk, jp_tem) = zn * r1_rau0 * ztm
            zn = rn_b0 * (1._wp - rn_lambda2 * zs - rn_mu2 * zh) - rn_nu * zt
            pab(ji, jj, jk, jp_sal) = zn * r1_rau0 * ztm
          END DO
        END DO
      END DO
    CASE DEFAULT
      WRITE(ctmp1, FMT = *) '          bad flag value for neos = ', neos
      CALL ctl_stop('rab_3d:', ctmp1)
    END SELECT
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = pab(:, :, :, jp_tem), clinfo1 = ' rab_3d_t: ', tab3d_2 = pab(:, :, :, jp_sal), clinfo2 = ' &
&rab_3d_s : ', kdim = jpk)
    IF (ln_timing) CALL timing_stop('rab_3d')
  END SUBROUTINE rab_3d
  SUBROUTINE rab_2d(pts, pdep, pab)
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpts), INTENT(IN) :: pts
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pdep
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpts), INTENT(OUT) :: pab
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zt, zh, zs
    REAL(KIND = wp) :: zn, zn0, zn1, zn2, zn3
    IF (ln_timing) CALL timing_start('rab_2d')
    pab(:, :, :) = 0._wp
    SELECT CASE (neos)
    CASE (np_teos10, np_eos80)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zh = pdep(ji, jj) * r1_Z0
          zt = pts(ji, jj, jp_tem) * r1_T0
          zs = SQRT(ABS(pts(ji, jj, jp_sal) + rdeltaS) * r1_S0)
          zn3 = ALP003
          zn2 = ALP012 * zt + ALP102 * zs + ALP002
          zn1 = ((ALP031 * zt + ALP121 * zs + ALP021) * zt + (ALP211 * zs + ALP111) * zs + ALP011) * zt + ((ALP301 * zs + ALP201) &
&* zs + ALP101) * zs + ALP001
          zn0 = ((((ALP050 * zt + ALP140 * zs + ALP040) * zt + (ALP230 * zs + ALP130) * zs + ALP030) * zt + ((ALP320 * zs + &
&ALP220) * zs + ALP120) * zs + ALP020) * zt + (((ALP410 * zs + ALP310) * zs + ALP210) * zs + ALP110) * zs + ALP010) * zt + &
&((((ALP500 * zs + ALP400) * zs + ALP300) * zs + ALP200) * zs + ALP100) * zs + ALP000
          zn = ((zn3 * zh + zn2) * zh + zn1) * zh + zn0
          pab(ji, jj, jp_tem) = zn * r1_rau0
          zn3 = BET003
          zn2 = BET012 * zt + BET102 * zs + BET002
          zn1 = ((BET031 * zt + BET121 * zs + BET021) * zt + (BET211 * zs + BET111) * zs + BET011) * zt + ((BET301 * zs + BET201) &
&* zs + BET101) * zs + BET001
          zn0 = ((((BET050 * zt + BET140 * zs + BET040) * zt + (BET230 * zs + BET130) * zs + BET030) * zt + ((BET320 * zs + &
&BET220) * zs + BET120) * zs + BET020) * zt + (((BET410 * zs + BET310) * zs + BET210) * zs + BET110) * zs + BET010) * zt + &
&((((BET500 * zs + BET400) * zs + BET300) * zs + BET200) * zs + BET100) * zs + BET000
          zn = ((zn3 * zh + zn2) * zh + zn1) * zh + zn0
          pab(ji, jj, jp_sal) = zn / zs * r1_rau0
        END DO
      END DO
      CALL lbc_lnk_multi('eosbn2', pab(:, :, jp_tem), 'T', 1., pab(:, :, jp_sal), 'T', 1.)
    CASE (np_seos)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zt = pts(ji, jj, jp_tem) - 10._wp
          zs = pts(ji, jj, jp_sal) - 35._wp
          zh = pdep(ji, jj)
          zn = rn_a0 * (1._wp + rn_lambda1 * zt + rn_mu1 * zh) + rn_nu * zs
          pab(ji, jj, jp_tem) = zn * r1_rau0
          zn = rn_b0 * (1._wp - rn_lambda2 * zs - rn_mu2 * zh) - rn_nu * zt
          pab(ji, jj, jp_sal) = zn * r1_rau0
        END DO
      END DO
      CALL lbc_lnk_multi('eosbn2', pab(:, :, jp_tem), 'T', 1., pab(:, :, jp_sal), 'T', 1.)
    CASE DEFAULT
      WRITE(ctmp1, FMT = *) '          bad flag value for neos = ', neos
      CALL ctl_stop('rab_2d:', ctmp1)
    END SELECT
    IF (ln_ctl) CALL prt_ctl(tab2d_1 = pab(:, :, jp_tem), clinfo1 = ' rab_2d_t: ', tab2d_2 = pab(:, :, jp_sal), clinfo2 = ' &
&rab_2d_s : ')
    IF (ln_timing) CALL timing_stop('rab_2d')
  END SUBROUTINE rab_2d
  SUBROUTINE rab_0d(pts, pdep, pab)
    REAL(KIND = wp), DIMENSION(jpts), INTENT(IN) :: pts
    REAL(KIND = wp), INTENT(IN) :: pdep
    REAL(KIND = wp), DIMENSION(jpts), INTENT(OUT) :: pab
    REAL(KIND = wp) :: zt, zh, zs
    REAL(KIND = wp) :: zn, zn0, zn1, zn2, zn3
    IF (ln_timing) CALL timing_start('rab_0d')
    pab(:) = 0._wp
    SELECT CASE (neos)
    CASE (np_teos10, np_eos80)
      zh = pdep * r1_Z0
      zt = pts(jp_tem) * r1_T0
      zs = SQRT(ABS(pts(jp_sal) + rdeltaS) * r1_S0)
      zn3 = ALP003
      zn2 = ALP012 * zt + ALP102 * zs + ALP002
      zn1 = ((ALP031 * zt + ALP121 * zs + ALP021) * zt + (ALP211 * zs + ALP111) * zs + ALP011) * zt + ((ALP301 * zs + ALP201) * zs &
&+ ALP101) * zs + ALP001
      zn0 = ((((ALP050 * zt + ALP140 * zs + ALP040) * zt + (ALP230 * zs + ALP130) * zs + ALP030) * zt + ((ALP320 * zs + ALP220) * &
&zs + ALP120) * zs + ALP020) * zt + (((ALP410 * zs + ALP310) * zs + ALP210) * zs + ALP110) * zs + ALP010) * zt + ((((ALP500 * zs + &
&ALP400) * zs + ALP300) * zs + ALP200) * zs + ALP100) * zs + ALP000
      zn = ((zn3 * zh + zn2) * zh + zn1) * zh + zn0
      pab(jp_tem) = zn * r1_rau0
      zn3 = BET003
      zn2 = BET012 * zt + BET102 * zs + BET002
      zn1 = ((BET031 * zt + BET121 * zs + BET021) * zt + (BET211 * zs + BET111) * zs + BET011) * zt + ((BET301 * zs + BET201) * zs &
&+ BET101) * zs + BET001
      zn0 = ((((BET050 * zt + BET140 * zs + BET040) * zt + (BET230 * zs + BET130) * zs + BET030) * zt + ((BET320 * zs + BET220) * &
&zs + BET120) * zs + BET020) * zt + (((BET410 * zs + BET310) * zs + BET210) * zs + BET110) * zs + BET010) * zt + ((((BET500 * zs + &
&BET400) * zs + BET300) * zs + BET200) * zs + BET100) * zs + BET000
      zn = ((zn3 * zh + zn2) * zh + zn1) * zh + zn0
      pab(jp_sal) = zn / zs * r1_rau0
    CASE (np_seos)
      zt = pts(jp_tem) - 10._wp
      zs = pts(jp_sal) - 35._wp
      zh = pdep
      zn = rn_a0 * (1._wp + rn_lambda1 * zt + rn_mu1 * zh) + rn_nu * zs
      pab(jp_tem) = zn * r1_rau0
      zn = rn_b0 * (1._wp - rn_lambda2 * zs - rn_mu2 * zh) - rn_nu * zt
      pab(jp_sal) = zn * r1_rau0
    CASE DEFAULT
      WRITE(ctmp1, FMT = *) '          bad flag value for neos = ', neos
      CALL ctl_stop('rab_0d:', ctmp1)
    END SELECT
    IF (ln_timing) CALL timing_stop('rab_0d')
  END SUBROUTINE rab_0d
  SUBROUTINE bn2(pts, pab, pn2)
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, jpts), INTENT(IN) :: pts
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, jpts), INTENT(IN) :: pab
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(OUT) :: pn2
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zaw, zbw, zrw
    IF (ln_timing) CALL timing_start('bn2')
    DO jk = 2, jpkm1
      DO jj = 1, jpj
        DO ji = 1, jpi
          zrw = (gdepw_n(ji, jj, jk) - gdept_n(ji, jj, jk)) / (gdept_n(ji, jj, jk - 1) - gdept_n(ji, jj, jk))
          zaw = pab(ji, jj, jk, jp_tem) * (1. - zrw) + pab(ji, jj, jk - 1, jp_tem) * zrw
          zbw = pab(ji, jj, jk, jp_sal) * (1. - zrw) + pab(ji, jj, jk - 1, jp_sal) * zrw
          pn2(ji, jj, jk) = grav * (zaw * (pts(ji, jj, jk - 1, jp_tem) - pts(ji, jj, jk, jp_tem)) - zbw * (pts(ji, jj, jk - 1, &
&jp_sal) - pts(ji, jj, jk, jp_sal))) / e3w_n(ji, jj, jk) * wmask(ji, jj, jk)
        END DO
      END DO
    END DO
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = pn2, clinfo1 = ' bn2  : ', kdim = jpk)
    IF (ln_timing) CALL timing_stop('bn2')
  END SUBROUTINE bn2
  FUNCTION eos_pt_from_ct(ctmp, psal) RESULT(ptmp)
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: ctmp
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: psal
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ptmp
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zt, zs, ztm
    REAL(KIND = wp) :: zn, zd
    REAL(KIND = wp) :: zdeltaS, z1_S0, z1_T0
    IF (ln_timing) CALL timing_start('eos_pt_from_ct')
    zdeltaS = 5._wp
    z1_S0 = 0.875_wp / 35.16504_wp
    z1_T0 = 1._wp / 40._wp
    DO jj = 1, jpj
      DO ji = 1, jpi
        zt = ctmp(ji, jj) * z1_T0
        zs = SQRT(ABS(psal(ji, jj) + zdeltaS) * r1_S0)
        ztm = tmask(ji, jj, 1)
        zn = ((((- 2.1385727895E-01_wp * zt - 2.7674419971E-01_wp * zs + 1.0728094330_wp) * zt + (2.6366564313_wp * zs + &
&3.3546960647_wp) * zs - 7.8012209473_wp) * zt + ((1.8835586562_wp * zs + 7.3949191679_wp) * zs - 3.3937395875_wp) * zs - &
&5.6414948432_wp) * zt + (((3.5737370589_wp * zs - 1.5512427389E+01_wp) * zs + 2.4625741105E+01_wp) * zs + 1.9912291000E+01_wp) * &
&zs - 3.2191146312E+01_wp) * zt + ((((5.7153204649E-01_wp * zs - 3.0943149543_wp) * zs + 9.3052495181_wp) * zs - 9.4528934807_wp) &
&* zs + 3.1066408996_wp) * zs - 4.3504021262E-01_wp
        zd = (2.0035003456_wp * zt - 3.4570358592E-01_wp * zs + 5.6471810638_wp) * zt + (1.5393993508_wp * zs - 6.9394762624_wp) * &
&zs + 1.2750522650E+01_wp
        ptmp(ji, jj) = (zt / z1_T0 + zn / zd) * ztm
      END DO
    END DO
    IF (ln_timing) CALL timing_stop('eos_pt_from_ct')
  END FUNCTION eos_pt_from_ct
  SUBROUTINE eos_fzp_2d(psal, ptf, pdep)
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: psal
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN), OPTIONAL :: pdep
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(OUT) :: ptf
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zt, zs, z1_S0
    SELECT CASE (neos)
    CASE (np_teos10, np_seos)
      z1_S0 = 1._wp / 35.16504_wp
      DO jj = 1, jpj
        DO ji = 1, jpi
          zs = SQRT(ABS(psal(ji, jj)) * z1_S0)
          ptf(ji, jj) = ((((1.46873E-03_wp * zs - 9.64972E-03_wp) * zs + 2.28348E-02_wp) * zs - 3.12775E-02_wp) * zs + &
&2.07679E-02_wp) * zs - 5.87701E-02_wp
        END DO
      END DO
      ptf(:, :) = ptf(:, :) * psal(:, :)
      IF (PRESENT(pdep)) ptf(:, :) = ptf(:, :) - 7.53E-4 * pdep(:, :)
    CASE (np_eos80)
      ptf(:, :) = (- 0.0575_wp + 1.710523E-3_wp * SQRT(psal(:, :)) - 2.154996E-4_wp * psal(:, :)) * psal(:, :)
      IF (PRESENT(pdep)) ptf(:, :) = ptf(:, :) - 7.53E-4 * pdep(:, :)
    CASE DEFAULT
      WRITE(ctmp1, FMT = *) '          bad flag value for neos = ', neos
      CALL ctl_stop('eos_fzp_2d:', ctmp1)
    END SELECT
  END SUBROUTINE eos_fzp_2d
  SUBROUTINE eos_fzp_0d(psal, ptf, pdep)
    REAL(KIND = wp), INTENT(IN) :: psal
    REAL(KIND = wp), INTENT(IN), OPTIONAL :: pdep
    REAL(KIND = wp), INTENT(OUT) :: ptf
    REAL(KIND = wp) :: zs
    SELECT CASE (neos)
    CASE (np_teos10, np_seos)
      zs = SQRT(ABS(psal) / 35.16504_wp)
      ptf = ((((1.46873E-03_wp * zs - 9.64972E-03_wp) * zs + 2.28348E-02_wp) * zs - 3.12775E-02_wp) * zs + 2.07679E-02_wp) * zs - &
&5.87701E-02_wp
      ptf = ptf * psal
      IF (PRESENT(pdep)) ptf = ptf - 7.53E-4 * pdep
    CASE (np_eos80)
      ptf = (- 0.0575_wp + 1.710523E-3_wp * SQRT(psal) - 2.154996E-4_wp * psal) * psal
      IF (PRESENT(pdep)) ptf = ptf - 7.53E-4 * pdep
    CASE DEFAULT
      WRITE(ctmp1, FMT = *) '          bad flag value for neos = ', neos
      CALL ctl_stop('eos_fzp_0d:', ctmp1)
    END SELECT
  END SUBROUTINE eos_fzp_0d
  SUBROUTINE eos_pen(pts, pab_pe, ppen)
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, jpts), INTENT(IN) :: pts
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, jpts), INTENT(OUT) :: pab_pe
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(OUT) :: ppen
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zt, zh, zs, ztm
    REAL(KIND = wp) :: zn, zn0, zn1, zn2
    IF (ln_timing) CALL timing_start('eos_pen')
    SELECT CASE (neos)
    CASE (np_teos10, np_eos80)
      DO jk = 1, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            zh = gdept_n(ji, jj, jk) * r1_Z0
            zt = pts(ji, jj, jk, jp_tem) * r1_T0
            zs = SQRT(ABS(pts(ji, jj, jk, jp_sal) + rdeltaS) * r1_S0)
            ztm = tmask(ji, jj, jk)
            zn2 = (PEN012) * zt + PEN102 * zs + PEN002
            zn1 = ((PEN021) * zt + PEN111 * zs + PEN011) * zt + (PEN201 * zs + PEN101) * zs + PEN001
            zn0 = ((((PEN040) * zt + PEN130 * zs + PEN030) * zt + (PEN220 * zs + PEN120) * zs + PEN020) * zt + ((PEN310 * zs + &
&PEN210) * zs + PEN110) * zs + PEN010) * zt + (((PEN400 * zs + PEN300) * zs + PEN200) * zs + PEN100) * zs + PEN000
            zn = (zn2 * zh + zn1) * zh + zn0
            ppen(ji, jj, jk) = zn * zh * r1_rau0 * ztm
            zn2 = APE002
            zn1 = (APE011) * zt + APE101 * zs + APE001
            zn0 = (((APE030) * zt + APE120 * zs + APE020) * zt + (APE210 * zs + APE110) * zs + APE010) * zt + ((APE300 * zs + &
&APE200) * zs + APE100) * zs + APE000
            zn = (zn2 * zh + zn1) * zh + zn0
            pab_pe(ji, jj, jk, jp_tem) = zn * zh * r1_rau0 * ztm
            zn2 = BPE002
            zn1 = (BPE011) * zt + BPE101 * zs + BPE001
            zn0 = (((BPE030) * zt + BPE120 * zs + BPE020) * zt + (BPE210 * zs + BPE110) * zs + BPE010) * zt + ((BPE300 * zs + &
&BPE200) * zs + BPE100) * zs + BPE000
            zn = (zn2 * zh + zn1) * zh + zn0
            pab_pe(ji, jj, jk, jp_sal) = zn / zs * zh * r1_rau0 * ztm
          END DO
        END DO
      END DO
    CASE (np_seos)
      DO jk = 1, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            zt = pts(ji, jj, jk, jp_tem) - 10._wp
            zs = pts(ji, jj, jk, jp_sal) - 35._wp
            zh = gdept_n(ji, jj, jk)
            ztm = tmask(ji, jj, jk)
            zn = 0.5_wp * zh * r1_rau0 * ztm
            ppen(ji, jj, jk) = (rn_a0 * rn_mu1 * zt + rn_b0 * rn_mu2 * zs) * zn
            pab_pe(ji, jj, jk, jp_tem) = - rn_a0 * rn_mu1 * zn
            pab_pe(ji, jj, jk, jp_sal) = rn_b0 * rn_mu2 * zn
          END DO
        END DO
      END DO
    CASE DEFAULT
      WRITE(ctmp1, FMT = *) '          bad flag value for neos = ', neos
      CALL ctl_stop('eos_pen:', ctmp1)
    END SELECT
    IF (ln_timing) CALL timing_stop('eos_pen')
  END SUBROUTINE eos_pen
  SUBROUTINE eos_init
    INTEGER :: ios
    INTEGER :: ioptio
    NAMELIST /nameos/ ln_TEOS10, ln_EOS80, ln_SEOS, rn_a0, rn_b0, rn_lambda1, rn_mu1, rn_lambda2, rn_mu2, rn_nu
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, nameos, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'nameos in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, nameos, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'nameos in configuration namelist', lwp)
    IF (lwm) WRITE(numond, nameos)
    rau0 = 1026._wp
    rcp = 3991.86795711963_wp
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'eos_init : equation of state'
      WRITE(numout, FMT = *) '~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist nameos : Chosen the Equation Of Seawater (EOS)'
      WRITE(numout, FMT = *) '      TEOS-10 : rho=F(Conservative Temperature, Absolute  Salinity, depth)   ln_TEOS10 = ', ln_TEOS10
      WRITE(numout, FMT = *) '      EOS-80  : rho=F(Potential    Temperature, Practical Salinity, depth)   ln_EOS80  = ', ln_EOS80
      WRITE(numout, FMT = *) '      S-EOS   : rho=F(Conservative Temperature, Absolute  Salinity, depth)   ln_SEOS   = ', ln_SEOS
    END IF
    ioptio = 0
    IF (ln_teos10) THEN
      ioptio = ioptio + 1
      neos = np_teos10
    END IF
    IF (ln_eos80) THEN
      ioptio = ioptio + 1
      neos = np_eos80
    END IF
    IF (ln_seos) THEN
      ioptio = ioptio + 1
      neos = np_seos
    END IF
    IF (ioptio /= 1) CALL ctl_stop("Exactly one equation of state option must be selected")
    SELECT CASE (neos)
    CASE (np_teos10)
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   use of TEOS-10 equation of state (cons. temp. and abs. salinity)'
      l_useCT = .TRUE.
      rdeltaS = 32._wp
      r1_S0 = 0.875_wp / 35.16504_wp
      r1_T0 = 1._wp / 40._wp
      r1_Z0 = 1.E-4_wp
      EOS000 = 8.0189615746E+02_wp
      EOS100 = 8.6672408165E+02_wp
      EOS200 = - 1.7864682637E+03_wp
      EOS300 = 2.0375295546E+03_wp
      EOS400 = - 1.2849161071E+03_wp
      EOS500 = 4.3227585684E+02_wp
      EOS600 = - 6.0579916612E+01_wp
      EOS010 = 2.6010145068E+01_wp
      EOS110 = - 6.5281885265E+01_wp
      EOS210 = 8.1770425108E+01_wp
      EOS310 = - 5.6888046321E+01_wp
      EOS410 = 1.7681814114E+01_wp
      EOS510 = - 1.9193502195_wp
      EOS020 = - 3.7074170417E+01_wp
      EOS120 = 6.1548258127E+01_wp
      EOS220 = - 6.0362551501E+01_wp
      EOS320 = 2.9130021253E+01_wp
      EOS420 = - 5.4723692739_wp
      EOS030 = 2.1661789529E+01_wp
      EOS130 = - 3.3449108469E+01_wp
      EOS230 = 1.9717078466E+01_wp
      EOS330 = - 3.1742946532_wp
      EOS040 = - 8.3627885467_wp
      EOS140 = 1.1311538584E+01_wp
      EOS240 = - 5.3563304045_wp
      EOS050 = 5.4048723791E-01_wp
      EOS150 = 4.8169980163E-01_wp
      EOS060 = - 1.9083568888E-01_wp
      EOS001 = 1.9681925209E+01_wp
      EOS101 = - 4.2549998214E+01_wp
      EOS201 = 5.0774768218E+01_wp
      EOS301 = - 3.0938076334E+01_wp
      EOS401 = 6.6051753097_wp
      EOS011 = - 1.3336301113E+01_wp
      EOS111 = - 4.4870114575_wp
      EOS211 = 5.0042598061_wp
      EOS311 = - 6.5399043664E-01_wp
      EOS021 = 6.7080479603_wp
      EOS121 = 3.5063081279_wp
      EOS221 = - 1.8795372996_wp
      EOS031 = - 2.4649669534_wp
      EOS131 = - 5.5077101279E-01_wp
      EOS041 = 5.5927935970E-01_wp
      EOS002 = 2.0660924175_wp
      EOS102 = - 4.9527603989_wp
      EOS202 = 2.5019633244_wp
      EOS012 = 2.0564311499_wp
      EOS112 = - 2.1311365518E-01_wp
      EOS022 = - 1.2419983026_wp
      EOS003 = - 2.3342758797E-02_wp
      EOS103 = - 1.8507636718E-02_wp
      EOS013 = 3.7969820455E-01_wp
      ALP000 = - 6.5025362670E-01_wp
      ALP100 = 1.6320471316_wp
      ALP200 = - 2.0442606277_wp
      ALP300 = 1.4222011580_wp
      ALP400 = - 4.4204535284E-01_wp
      ALP500 = 4.7983755487E-02_wp
      ALP010 = 1.8537085209_wp
      ALP110 = - 3.0774129064_wp
      ALP210 = 3.0181275751_wp
      ALP310 = - 1.4565010626_wp
      ALP410 = 2.7361846370E-01_wp
      ALP020 = - 1.6246342147_wp
      ALP120 = 2.5086831352_wp
      ALP220 = - 1.4787808849_wp
      ALP320 = 2.3807209899E-01_wp
      ALP030 = 8.3627885467E-01_wp
      ALP130 = - 1.1311538584_wp
      ALP230 = 5.3563304045E-01_wp
      ALP040 = - 6.7560904739E-02_wp
      ALP140 = - 6.0212475204E-02_wp
      ALP050 = 2.8625353333E-02_wp
      ALP001 = 3.3340752782E-01_wp
      ALP101 = 1.1217528644E-01_wp
      ALP201 = - 1.2510649515E-01_wp
      ALP301 = 1.6349760916E-02_wp
      ALP011 = - 3.3540239802E-01_wp
      ALP111 = - 1.7531540640E-01_wp
      ALP211 = 9.3976864981E-02_wp
      ALP021 = 1.8487252150E-01_wp
      ALP121 = 4.1307825959E-02_wp
      ALP031 = - 5.5927935970E-02_wp
      ALP002 = - 5.1410778748E-02_wp
      ALP102 = 5.3278413794E-03_wp
      ALP012 = 6.2099915132E-02_wp
      ALP003 = - 9.4924551138E-03_wp
      BET000 = 1.0783203594E+01_wp
      BET100 = - 4.4452095908E+01_wp
      BET200 = 7.6048755820E+01_wp
      BET300 = - 6.3944280668E+01_wp
      BET400 = 2.6890441098E+01_wp
      BET500 = - 4.5221697773_wp
      BET010 = - 8.1219372432E-01_wp
      BET110 = 2.0346663041_wp
      BET210 = - 2.1232895170_wp
      BET310 = 8.7994140485E-01_wp
      BET410 = - 1.1939638360E-01_wp
      BET020 = 7.6574242289E-01_wp
      BET120 = - 1.5019813020_wp
      BET220 = 1.0872489522_wp
      BET320 = - 2.7233429080E-01_wp
      BET030 = - 4.1615152308E-01_wp
      BET130 = 4.9061350869E-01_wp
      BET230 = - 1.1847737788E-01_wp
      BET040 = 1.4073062708E-01_wp
      BET140 = - 1.3327978879E-01_wp
      BET050 = 5.9929880134E-03_wp
      BET001 = - 5.2937873009E-01_wp
      BET101 = 1.2634116779_wp
      BET201 = - 1.1547328025_wp
      BET301 = 3.2870876279E-01_wp
      BET011 = - 5.5824407214E-02_wp
      BET111 = 1.2451933313E-01_wp
      BET211 = - 2.4409539932E-02_wp
      BET021 = 4.3623149752E-02_wp
      BET121 = - 4.6767901790E-02_wp
      BET031 = - 6.8523260060E-03_wp
      BET002 = - 6.1618945251E-02_wp
      BET102 = 6.2255521644E-02_wp
      BET012 = - 2.6514181169E-03_wp
      BET003 = - 2.3025968587E-04_wp
      PEN000 = - 9.8409626043_wp
      PEN100 = 2.1274999107E+01_wp
      PEN200 = - 2.5387384109E+01_wp
      PEN300 = 1.5469038167E+01_wp
      PEN400 = - 3.3025876549_wp
      PEN010 = 6.6681505563_wp
      PEN110 = 2.2435057288_wp
      PEN210 = - 2.5021299030_wp
      PEN310 = 3.2699521832E-01_wp
      PEN020 = - 3.3540239802_wp
      PEN120 = - 1.7531540640_wp
      PEN220 = 9.3976864981E-01_wp
      PEN030 = 1.2324834767_wp
      PEN130 = 2.7538550639E-01_wp
      PEN040 = - 2.7963967985E-01_wp
      PEN001 = - 1.3773949450_wp
      PEN101 = 3.3018402659_wp
      PEN201 = - 1.6679755496_wp
      PEN011 = - 1.3709540999_wp
      PEN111 = 1.4207577012E-01_wp
      PEN021 = 8.2799886843E-01_wp
      PEN002 = 1.7507069098E-02_wp
      PEN102 = 1.3880727538E-02_wp
      PEN012 = - 2.8477365341E-01_wp
      APE000 = - 1.6670376391E-01_wp
      APE100 = - 5.6087643219E-02_wp
      APE200 = 6.2553247576E-02_wp
      APE300 = - 8.1748804580E-03_wp
      APE010 = 1.6770119901E-01_wp
      APE110 = 8.7657703198E-02_wp
      APE210 = - 4.6988432490E-02_wp
      APE020 = - 9.2436260751E-02_wp
      APE120 = - 2.0653912979E-02_wp
      APE030 = 2.7963967985E-02_wp
      APE001 = 3.4273852498E-02_wp
      APE101 = - 3.5518942529E-03_wp
      APE011 = - 4.1399943421E-02_wp
      APE002 = 7.1193413354E-03_wp
      BPE000 = 2.6468936504E-01_wp
      BPE100 = - 6.3170583896E-01_wp
      BPE200 = 5.7736640125E-01_wp
      BPE300 = - 1.6435438140E-01_wp
      BPE010 = 2.7912203607E-02_wp
      BPE110 = - 6.2259666565E-02_wp
      BPE210 = 1.2204769966E-02_wp
      BPE020 = - 2.1811574876E-02_wp
      BPE120 = 2.3383950895E-02_wp
      BPE030 = 3.4261630030E-03_wp
      BPE001 = 4.1079296834E-02_wp
      BPE101 = - 4.1503681096E-02_wp
      BPE011 = 1.7676120780E-03_wp
      BPE002 = 1.7269476440E-04_wp
    CASE (np_eos80)
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   use of EOS-80 equation of state (pot. temp. and pract. salinity)'
      l_useCT = .FALSE.
      rdeltaS = 20._wp
      r1_S0 = 1._wp / 40._wp
      r1_T0 = 1._wp / 40._wp
      r1_Z0 = 1.E-4_wp
      EOS000 = 9.5356891948E+02_wp
      EOS100 = 1.7136499189E+02_wp
      EOS200 = - 3.7501039454E+02_wp
      EOS300 = 5.1856810420E+02_wp
      EOS400 = - 3.7264470465E+02_wp
      EOS500 = 1.4302533998E+02_wp
      EOS600 = - 2.2856621162E+01_wp
      EOS010 = 1.0087518651E+01_wp
      EOS110 = - 1.3647741861E+01_wp
      EOS210 = 8.8478359933_wp
      EOS310 = - 7.2329388377_wp
      EOS410 = 1.4774410611_wp
      EOS510 = 2.0036720553E-01_wp
      EOS020 = - 2.5579830599E+01_wp
      EOS120 = 2.4043512327E+01_wp
      EOS220 = - 1.6807503990E+01_wp
      EOS320 = 8.3811577084_wp
      EOS420 = - 1.9771060192_wp
      EOS030 = 1.6846451198E+01_wp
      EOS130 = - 2.1482926901E+01_wp
      EOS230 = 1.0108954054E+01_wp
      EOS330 = - 6.2675951440E-01_wp
      EOS040 = - 8.0812310102_wp
      EOS140 = 1.0102374985E+01_wp
      EOS240 = - 4.8340368631_wp
      EOS050 = 1.2079167803_wp
      EOS150 = 1.1515380987E-01_wp
      EOS060 = - 2.4520288837E-01_wp
      EOS001 = 1.0748601068E+01_wp
      EOS101 = - 1.7817043500E+01_wp
      EOS201 = 2.2181366768E+01_wp
      EOS301 = - 1.6750916338E+01_wp
      EOS401 = 4.1202230403_wp
      EOS011 = - 1.5852644587E+01_wp
      EOS111 = - 7.6639383522E-01_wp
      EOS211 = 4.1144627302_wp
      EOS311 = - 6.6955877448E-01_wp
      EOS021 = 9.9994861860_wp
      EOS121 = - 1.9467067787E-01_wp
      EOS221 = - 1.2177554330_wp
      EOS031 = - 3.4866102017_wp
      EOS131 = 2.2229155620E-01_wp
      EOS041 = 5.9503008642E-01_wp
      EOS002 = 1.0375676547_wp
      EOS102 = - 3.4249470629_wp
      EOS202 = 2.0542026429_wp
      EOS012 = 2.1836324814_wp
      EOS112 = - 3.4453674320E-01_wp
      EOS022 = - 1.2548163097_wp
      EOS003 = 1.8729078427E-02_wp
      EOS103 = - 5.7238495240E-02_wp
      EOS013 = 3.8306136687E-01_wp
      ALP000 = - 2.5218796628E-01_wp
      ALP100 = 3.4119354654E-01_wp
      ALP200 = - 2.2119589983E-01_wp
      ALP300 = 1.8082347094E-01_wp
      ALP400 = - 3.6936026529E-02_wp
      ALP500 = - 5.0091801383E-03_wp
      ALP010 = 1.2789915300_wp
      ALP110 = - 1.2021756164_wp
      ALP210 = 8.4037519952E-01_wp
      ALP310 = - 4.1905788542E-01_wp
      ALP410 = 9.8855300959E-02_wp
      ALP020 = - 1.2634838399_wp
      ALP120 = 1.6112195176_wp
      ALP220 = - 7.5817155402E-01_wp
      ALP320 = 4.7006963580E-02_wp
      ALP030 = 8.0812310102E-01_wp
      ALP130 = - 1.0102374985_wp
      ALP230 = 4.8340368631E-01_wp
      ALP040 = - 1.5098959754E-01_wp
      ALP140 = - 1.4394226233E-02_wp
      ALP050 = 3.6780433255E-02_wp
      ALP001 = 3.9631611467E-01_wp
      ALP101 = 1.9159845880E-02_wp
      ALP201 = - 1.0286156825E-01_wp
      ALP301 = 1.6738969362E-02_wp
      ALP011 = - 4.9997430930E-01_wp
      ALP111 = 9.7335338937E-03_wp
      ALP211 = 6.0887771651E-02_wp
      ALP021 = 2.6149576513E-01_wp
      ALP121 = - 1.6671866715E-02_wp
      ALP031 = - 5.9503008642E-02_wp
      ALP002 = - 5.4590812035E-02_wp
      ALP102 = 8.6134185799E-03_wp
      ALP012 = 6.2740815484E-02_wp
      ALP003 = - 9.5765341718E-03_wp
      BET000 = 2.1420623987_wp
      BET100 = - 9.3752598635_wp
      BET200 = 1.9446303907E+01_wp
      BET300 = - 1.8632235232E+01_wp
      BET400 = 8.9390837485_wp
      BET500 = - 1.7142465871_wp
      BET010 = - 1.7059677327E-01_wp
      BET110 = 2.2119589983E-01_wp
      BET210 = - 2.7123520642E-01_wp
      BET310 = 7.3872053057E-02_wp
      BET410 = 1.2522950346E-02_wp
      BET020 = 3.0054390409E-01_wp
      BET120 = - 4.2018759976E-01_wp
      BET220 = 3.1429341406E-01_wp
      BET320 = - 9.8855300959E-02_wp
      BET030 = - 2.6853658626E-01_wp
      BET130 = 2.5272385134E-01_wp
      BET230 = - 2.3503481790E-02_wp
      BET040 = 1.2627968731E-01_wp
      BET140 = - 1.2085092158E-01_wp
      BET050 = 1.4394226233E-03_wp
      BET001 = - 2.2271304375E-01_wp
      BET101 = 5.5453416919E-01_wp
      BET201 = - 6.2815936268E-01_wp
      BET301 = 2.0601115202E-01_wp
      BET011 = - 9.5799229402E-03_wp
      BET111 = 1.0286156825E-01_wp
      BET211 = - 2.5108454043E-02_wp
      BET021 = - 2.4333834734E-03_wp
      BET121 = - 3.0443885826E-02_wp
      BET031 = 2.7786444526E-03_wp
      BET002 = - 4.2811838287E-02_wp
      BET102 = 5.1355066072E-02_wp
      BET012 = - 4.3067092900E-03_wp
      BET003 = - 7.1548119050E-04_wp
      PEN000 = - 5.3743005340_wp
      PEN100 = 8.9085217499_wp
      PEN200 = - 1.1090683384E+01_wp
      PEN300 = 8.3754581690_wp
      PEN400 = - 2.0601115202_wp
      PEN010 = 7.9263222935_wp
      PEN110 = 3.8319691761E-01_wp
      PEN210 = - 2.0572313651_wp
      PEN310 = 3.3477938724E-01_wp
      PEN020 = - 4.9997430930_wp
      PEN120 = 9.7335338937E-02_wp
      PEN220 = 6.0887771651E-01_wp
      PEN030 = 1.7433051009_wp
      PEN130 = - 1.1114577810E-01_wp
      PEN040 = - 2.9751504321E-01_wp
      PEN001 = - 6.9171176978E-01_wp
      PEN101 = 2.2832980419_wp
      PEN201 = - 1.3694684286_wp
      PEN011 = - 1.4557549876_wp
      PEN111 = 2.2969116213E-01_wp
      PEN021 = 8.3654420645E-01_wp
      PEN002 = - 1.4046808820E-02_wp
      PEN102 = 4.2928871430E-02_wp
      PEN012 = - 2.8729602515E-01_wp
      APE000 = - 1.9815805734E-01_wp
      APE100 = - 9.5799229402E-03_wp
      APE200 = 5.1430784127E-02_wp
      APE300 = - 8.3694846809E-03_wp
      APE010 = 2.4998715465E-01_wp
      APE110 = - 4.8667669469E-03_wp
      APE210 = - 3.0443885826E-02_wp
      APE020 = - 1.3074788257E-01_wp
      APE120 = 8.3359333577E-03_wp
      APE030 = 2.9751504321E-02_wp
      APE001 = 3.6393874690E-02_wp
      APE101 = - 5.7422790533E-03_wp
      APE011 = - 4.1827210323E-02_wp
      APE002 = 7.1824006288E-03_wp
      BPE000 = 1.1135652187E-01_wp
      BPE100 = - 2.7726708459E-01_wp
      BPE200 = 3.1407968134E-01_wp
      BPE300 = - 1.0300557601E-01_wp
      BPE010 = 4.7899614701E-03_wp
      BPE110 = - 5.1430784127E-02_wp
      BPE210 = 1.2554227021E-02_wp
      BPE020 = 1.2166917367E-03_wp
      BPE120 = 1.5221942913E-02_wp
      BPE030 = - 1.3893222263E-03_wp
      BPE001 = 2.8541225524E-02_wp
      BPE101 = - 3.4236710714E-02_wp
      BPE011 = 2.8711395266E-03_wp
      BPE002 = 5.3661089288E-04_wp
    CASE (np_seos)
      IF (lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) '   ==>>>   use of simplified eos:    '
        WRITE(numout, FMT = *) '              rhd(dT=T-10,dS=S-35,Z) = [-a0*(1+lambda1/2*dT+mu1*Z)*dT '
        WRITE(numout, FMT = *) '                                       + b0*(1+lambda2/2*dT+mu2*Z)*dS - nu*dT*dS] / rau0'
        WRITE(numout, FMT = *) '              with the following coefficients :'
        WRITE(numout, FMT = *) '                 thermal exp. coef.    rn_a0      = ', rn_a0
        WRITE(numout, FMT = *) '                 saline  cont. coef.   rn_b0      = ', rn_b0
        WRITE(numout, FMT = *) '                 cabbeling coef.       rn_lambda1 = ', rn_lambda1
        WRITE(numout, FMT = *) '                 cabbeling coef.       rn_lambda2 = ', rn_lambda2
        WRITE(numout, FMT = *) '                 thermobar. coef.      rn_mu1     = ', rn_mu1
        WRITE(numout, FMT = *) '                 thermobar. coef.      rn_mu2     = ', rn_mu2
        WRITE(numout, FMT = *) '                 2nd cabbel. coef.     rn_nu      = ', rn_nu
        WRITE(numout, FMT = *) '              Caution: rn_beta0=0 incompatible with ddm parameterization '
      END IF
      l_useCT = .TRUE.
    CASE DEFAULT
      WRITE(ctmp1, FMT = *) '          bad flag value for neos = ', neos, '. You should never see this error'
      CALL ctl_stop(ctmp1)
    END SELECT
    rau0_rcp = rau0 * rcp
    r1_rau0 = 1._wp / rau0
    r1_rcp = 1._wp / rcp
    r1_rau0_rcp = 1._wp / rau0_rcp
    IF (lwp) THEN
      IF (l_useCT) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) '   ==>>>   model uses Conservative Temperature'
        WRITE(numout, FMT = *) '           Important: model must be initialized with CT and SA fields'
      ELSE
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) '   ==>>>   model does not use Conservative Temperature'
      END IF
    END IF
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) '   Associated physical constant'
    IF (lwp) WRITE(numout, FMT = *) '      volumic mass of reference           rau0  = ', rau0, ' kg/m^3'
    IF (lwp) WRITE(numout, FMT = *) '      1. / rau0                        r1_rau0  = ', r1_rau0, ' m^3/kg'
    IF (lwp) WRITE(numout, FMT = *) '      ocean specific heat                 rcp   = ', rcp, ' J/Kelvin'
    IF (lwp) WRITE(numout, FMT = *) '      rau0 * rcp                       rau0_rcp = ', rau0_rcp
    IF (lwp) WRITE(numout, FMT = *) '      1. / ( rau0 * rcp )           r1_rau0_rcp = ', r1_rau0_rcp
  END SUBROUTINE eos_init
END MODULE eosbn2