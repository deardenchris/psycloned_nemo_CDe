MODULE zdf_oce
  USE par_oce
  USE in_out_manager
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: zdf_oce_alloc
  LOGICAL, PUBLIC :: ln_zad_Aimp
  LOGICAL, PUBLIC :: ln_zdfcst
  LOGICAL, PUBLIC :: ln_zdfric
  LOGICAL, PUBLIC :: ln_zdftke
  LOGICAL, PUBLIC :: ln_zdfgls
  LOGICAL, PUBLIC :: ln_zdfosm
  LOGICAL, PUBLIC :: ln_zdfevd
  INTEGER, PUBLIC :: nn_evdm
  REAL(KIND = wp), PUBLIC :: rn_evd
  LOGICAL, PUBLIC :: ln_zdfnpc
  INTEGER, PUBLIC :: nn_npc
  INTEGER, PUBLIC :: nn_npcp
  LOGICAL, PUBLIC :: ln_zdfddm
  REAL(KIND = wp), PUBLIC :: rn_avts
  REAL(KIND = wp), PUBLIC :: rn_hsbfr
  LOGICAL, PUBLIC :: ln_zdfswm
  LOGICAL, PUBLIC :: ln_zdfiwm
  LOGICAL, PUBLIC :: ln_zdftmx
  REAL(KIND = wp), PUBLIC :: rn_avm0
  REAL(KIND = wp), PUBLIC :: rn_avt0
  INTEGER, PUBLIC :: nn_avb
  INTEGER, PUBLIC :: nn_havtb
  REAL(KIND = wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:, :, :) :: avm, avt, avs
  REAL(KIND = wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:, :, :) :: avm_k, avt_k
  REAL(KIND = wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:, :, :) :: en
  REAL(KIND = wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:) :: avmb, avtb
  REAL(KIND = wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:, :) :: avtb_2d
  CONTAINS
  INTEGER FUNCTION zdf_oce_alloc()
    ALLOCATE(avm(jpi, jpj, jpk), avm_k(jpi, jpj, jpk), avs(jpi, jpj, jpk), avt(jpi, jpj, jpk), avt_k(jpi, jpj, jpk), en(jpi, jpj, jpk), avmb(jpk), avtb(jpk), avtb_2d(jpi, jpj), STAT = zdf_oce_alloc)
    IF (zdf_oce_alloc /= 0) CALL ctl_stop('STOP', 'zdf_oce_alloc: failed to allocate arrays')
  END FUNCTION zdf_oce_alloc
END MODULE zdf_oce