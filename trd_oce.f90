MODULE trd_oce
  USE par_oce
  USE trdmxl_oce
  USE trdvor_oce
  IMPLICIT NONE
  PUBLIC
  LOGICAL, PUBLIC :: ln_dyn_trd = .FALSE.
  LOGICAL, PUBLIC :: ln_tra_trd = .FALSE.
  LOGICAL, PUBLIC :: ln_KE_trd = .FALSE.
  LOGICAL, PUBLIC :: ln_PE_trd = .FALSE.
  LOGICAL, PUBLIC :: ln_vor_trd = .FALSE.
  LOGICAL, PUBLIC :: ln_glo_trd = .FALSE.
  LOGICAL, PUBLIC :: ln_dyn_mxl = .FALSE.
  LOGICAL, PUBLIC :: ln_tra_mxl = .FALSE.
  INTEGER, PUBLIC :: nn_trd = 10
  LOGICAL, PUBLIC :: l_trdtra
  LOGICAL, PUBLIC :: l_trddyn
  LOGICAL, PUBLIC :: l_trdtrc = .FALSE.
  INTEGER, PUBLIC, PARAMETER :: jptot_tra = 20
  INTEGER, PUBLIC, PARAMETER :: jptra_xad = 1
  INTEGER, PUBLIC, PARAMETER :: jptra_yad = 2
  INTEGER, PUBLIC, PARAMETER :: jptra_zad = 3
  INTEGER, PUBLIC, PARAMETER :: jptra_sad = 4
  INTEGER, PUBLIC, PARAMETER :: jptra_totad = 5
  INTEGER, PUBLIC, PARAMETER :: jptra_ldf = 6
  INTEGER, PUBLIC, PARAMETER :: jptra_zdf = 7
  INTEGER, PUBLIC, PARAMETER :: jptra_zdfp = 8
  INTEGER, PUBLIC, PARAMETER :: jptra_evd = 9
  INTEGER, PUBLIC, PARAMETER :: jptra_bbc = 10
  INTEGER, PUBLIC, PARAMETER :: jptra_bbl = 11
  INTEGER, PUBLIC, PARAMETER :: jptra_npc = 12
  INTEGER, PUBLIC, PARAMETER :: jptra_dmp = 13
  INTEGER, PUBLIC, PARAMETER :: jptra_qsr = 14
  INTEGER, PUBLIC, PARAMETER :: jptra_nsr = 15
  INTEGER, PUBLIC, PARAMETER :: jptra_atf = 16
  INTEGER, PUBLIC, PARAMETER :: jptra_tot = 17
  INTEGER, PUBLIC, PARAMETER :: jptra_sms = 18
  INTEGER, PUBLIC, PARAMETER :: jptra_radn = 19
  INTEGER, PUBLIC, PARAMETER :: jptra_radb = 20
  INTEGER, PUBLIC, PARAMETER :: jptot_dyn = 13
  INTEGER, PUBLIC, PARAMETER :: jpdyn_hpg = 1
  INTEGER, PUBLIC, PARAMETER :: jpdyn_spg = 2
  INTEGER, PUBLIC, PARAMETER :: jpdyn_keg = 3
  INTEGER, PUBLIC, PARAMETER :: jpdyn_rvo = 4
  INTEGER, PUBLIC, PARAMETER :: jpdyn_pvo = 5
  INTEGER, PUBLIC, PARAMETER :: jpdyn_zad = 6
  INTEGER, PUBLIC, PARAMETER :: jpdyn_ldf = 7
  INTEGER, PUBLIC, PARAMETER :: jpdyn_zdf = 8
  INTEGER, PUBLIC, PARAMETER :: jpdyn_bfr = 9
  INTEGER, PUBLIC, PARAMETER :: jpdyn_atf = 10
  INTEGER, PUBLIC, PARAMETER :: jpdyn_tau = 11
  INTEGER, PUBLIC, PARAMETER :: jpdyn_bfri = 12
  INTEGER, PUBLIC, PARAMETER :: jpdyn_ken = 13
END MODULE trd_oce