MODULE par_oce
  USE par_kind
  IMPLICIT NONE
  PUBLIC
  LOGICAL :: ln_read_cfg
  CHARACTER(LEN = lc) :: cn_domcfg
  LOGICAL :: ln_write_cfg
  CHARACTER(LEN = lc) :: cn_domcfg_out
  LOGICAL :: ln_use_jattr
  CHARACTER(LEN = lc) :: cn_cfg
  INTEGER :: nn_cfg
  INTEGER :: jpiglo
  INTEGER :: jpjglo
  INTEGER :: jpkglo
  INTEGER, PUBLIC :: nbug_in_agrif_conv_do_not_remove_or_modify = 1 - 1
  INTEGER, PUBLIC, PARAMETER :: nbghostcells = 3
  INTEGER, PUBLIC :: nbcellsx
  INTEGER, PUBLIC :: nbcellsy
  INTEGER, PUBLIC :: jpi
  INTEGER, PUBLIC :: jpj
  INTEGER, PUBLIC :: jpk
  INTEGER, PUBLIC :: jpim1
  INTEGER, PUBLIC :: jpjm1
  INTEGER, PUBLIC :: jpkm1
  INTEGER, PUBLIC :: jpij
  INTEGER, PUBLIC :: jpimax
  INTEGER, PUBLIC :: jpjmax
  INTEGER, PUBLIC, PARAMETER :: jpts = 2
  INTEGER, PUBLIC, PARAMETER :: jp_tem = 1
  INTEGER, PUBLIC, PARAMETER :: jp_sal = 2
  INTEGER, PUBLIC :: jpni
  INTEGER, PUBLIC :: jpnj
  INTEGER, PUBLIC :: jpnij
  INTEGER, PUBLIC, PARAMETER :: jpr2di = 0
  INTEGER, PUBLIC, PARAMETER :: jpr2dj = 0
  INTEGER, PUBLIC, PARAMETER :: nn_hls = 1
END MODULE par_oce