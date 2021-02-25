MODULE in_out_manager
  USE par_oce
  USE nc4interface
  IMPLICIT NONE
  PUBLIC
  CHARACTER(LEN = lc) :: cn_exp
  CHARACTER(LEN = lc) :: cn_ocerst_in
  CHARACTER(LEN = lc) :: cn_ocerst_indir
  CHARACTER(LEN = lc) :: cn_ocerst_out
  CHARACTER(LEN = lc) :: cn_ocerst_outdir
  LOGICAL :: ln_rstart
  LOGICAL :: ln_rst_list
  INTEGER :: nn_rstctl
  INTEGER :: nn_rstssh = 0
  INTEGER :: nn_it000
  INTEGER :: nn_itend
  INTEGER :: nn_date0
  INTEGER :: nn_time0
  INTEGER :: nn_leapy
  INTEGER :: nn_istate
  INTEGER :: nn_write
  INTEGER :: nn_stock
  INTEGER, DIMENSION(10) :: nn_stocklist
  LOGICAL :: ln_mskland
  LOGICAL :: ln_rstdate
  LOGICAL :: ln_cfmeta
  LOGICAL :: ln_clobber
  INTEGER :: nn_chunksz
  LOGICAL :: ln_xios_read
  INTEGER :: nn_wxios
  INTEGER :: nn_no
  TYPE(snc4_ctl) :: snc4set
  CHARACTER(LEN = lc) :: cexper
  INTEGER :: nrstdt
  INTEGER :: nit000
  INTEGER :: nitend
  INTEGER :: ndate0
  INTEGER :: nleapy
  INTEGER :: ninist
  INTEGER :: nwrite
  INTEGER :: nstock
  INTEGER, DIMENSION(10) :: nstocklist
  INTEGER :: nitrst
  LOGICAL :: lrst_oce
  LOGICAL :: lrst_ice
  INTEGER :: numror = 0
  INTEGER :: numrir
  INTEGER :: numrow
  INTEGER :: numriw
  INTEGER :: nrst_lst
  LOGICAL :: ln_ctl
  TYPE :: sn_ctl
    LOGICAL :: l_config = .FALSE.
    LOGICAL :: l_runstat = .FALSE.
    LOGICAL :: l_trcstat = .FALSE.
    LOGICAL :: l_oceout = .FALSE.
    LOGICAL :: l_layout = .FALSE.
    LOGICAL :: l_mppout = .FALSE.
    LOGICAL :: l_mpptop = .FALSE.
    INTEGER :: procmin = 0
    INTEGER :: procmax = 1000000
    INTEGER :: procincr = 1
    INTEGER :: ptimincr = 1
  END TYPE
  TYPE(sn_ctl), SAVE :: sn_cfctl
  LOGICAL :: ln_timing
  LOGICAL :: ln_diacfl
  INTEGER :: nn_print
  INTEGER :: nn_ictls
  INTEGER :: nn_ictle
  INTEGER :: nn_jctls
  INTEGER :: nn_jctle
  INTEGER :: nn_isplt
  INTEGER :: nn_jsplt
  INTEGER :: nprint, nictls, nictle, njctls, njctle, isplt, jsplt
  INTEGER :: ijsplt = 1
  INTEGER :: numstp = - 1
  INTEGER :: numtime = - 1
  INTEGER :: numout = 6
  INTEGER :: numnul = - 1
  INTEGER :: numnam_ref = - 1
  INTEGER :: numnam_cfg = - 1
  INTEGER :: numond = - 1
  INTEGER :: numnam_ice_ref = - 1
  INTEGER :: numnam_ice_cfg = - 1
  INTEGER :: numoni = - 1
  INTEGER :: numevo_ice = - 1
  INTEGER :: numrun = - 1
  INTEGER :: numdct_in = - 1
  INTEGER :: numdct_vol = - 1
  INTEGER :: numdct_heat = - 1
  INTEGER :: numdct_salt = - 1
  INTEGER :: numfl = - 1
  INTEGER :: numflo = - 1
  INTEGER :: no_print = 0
  INTEGER :: nstop = 0
  INTEGER :: nwarn = 0
  CHARACTER(LEN = lc) :: ctmp1, ctmp2, ctmp3
  CHARACTER(LEN = lc) :: ctmp4, ctmp5, ctmp6
  CHARACTER(LEN = lc) :: ctmp7, ctmp8, ctmp9
  CHARACTER(LEN = lc) :: ctmp10
  CHARACTER(LEN = lc) :: cform_err = "(/,' ===>>> : E R R O R',     /,'         ===========',/)"
  CHARACTER(LEN = lc) :: cform_war = "(/,' ===>>> : W A R N I N G', /,'         ===============',/)"
  LOGICAL :: lwm = .FALSE.
  LOGICAL :: lwp = .FALSE.
  LOGICAL :: lsp_area = .TRUE.
  CHARACTER(LEN = lc) :: cxios_context
  CHARACTER(LEN = lc) :: crxios_context
  CHARACTER(LEN = lc) :: cwxios_context
END MODULE in_out_manager