MODULE trdvor_oce
  USE par_oce
  IMPLICIT NONE
  PRIVATE
  INTEGER, PUBLIC, PARAMETER :: jpltot_vor = 11
  INTEGER, PUBLIC, PARAMETER :: jpvor_prg = 1
  INTEGER, PUBLIC, PARAMETER :: jpvor_keg = 2
  INTEGER, PUBLIC, PARAMETER :: jpvor_rvo = 3
  INTEGER, PUBLIC, PARAMETER :: jpvor_pvo = 4
  INTEGER, PUBLIC, PARAMETER :: jpvor_ldf = 5
  INTEGER, PUBLIC, PARAMETER :: jpvor_zad = 6
  INTEGER, PUBLIC, PARAMETER :: jpvor_zdf = 7
  INTEGER, PUBLIC, PARAMETER :: jpvor_spg = 8
  INTEGER, PUBLIC, PARAMETER :: jpvor_bev = 9
  INTEGER, PUBLIC, PARAMETER :: jpvor_swf = 10
  INTEGER, PUBLIC, PARAMETER :: jpvor_bfr = 11
END MODULE trdvor_oce