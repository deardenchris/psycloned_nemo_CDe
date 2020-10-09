MODULE dom_oce
  USE par_oce
  IMPLICIT NONE
  PUBLIC
  PUBLIC :: dom_oce_alloc
  LOGICAL, PUBLIC :: ln_linssh
  LOGICAL, PUBLIC :: ln_meshmask
  REAL(KIND = wp), PUBLIC :: rn_isfhmin
  REAL(KIND = wp), PUBLIC :: rn_rdt
  REAL(KIND = wp), PUBLIC :: rn_atfp
  INTEGER, PUBLIC :: nn_euler
  LOGICAL, PUBLIC :: ln_iscpl
  LOGICAL, PUBLIC :: ln_crs
  LOGICAL, PUBLIC :: ln_dynspg_exp
  LOGICAL, PUBLIC :: ln_dynspg_ts
  LOGICAL, PUBLIC :: ln_bt_fw
  LOGICAL, PUBLIC :: ln_bt_av
  LOGICAL, PUBLIC :: ln_bt_auto
  INTEGER, PUBLIC :: nn_bt_flt
  INTEGER, PUBLIC :: nn_baro
  REAL(KIND = wp), PUBLIC :: rn_bt_cmax
  REAL(KIND = wp), PUBLIC :: rn_bt_alpha
  REAL(KIND = wp), PUBLIC :: atfp
  REAL(KIND = wp), PUBLIC :: rdt
  INTEGER, PUBLIC :: neuler
  REAL(KIND = wp), PUBLIC :: r2dt
  INTEGER, PUBLIC :: jperio
  LOGICAL, PUBLIC :: l_Iperio, l_Jperio
  INTEGER, PUBLIC :: nimpp, njmpp
  INTEGER, PUBLIC :: nreci, nrecj
  INTEGER, PUBLIC :: nproc
  INTEGER, PUBLIC :: narea
  INTEGER, PUBLIC :: nbondi, nbondj
  INTEGER, ALLOCATABLE, PUBLIC :: nbondi_bdy(:)
  INTEGER, ALLOCATABLE, PUBLIC :: nbondj_bdy(:)
  INTEGER, ALLOCATABLE, PUBLIC :: nbondi_bdy_b(:)
  INTEGER, ALLOCATABLE, PUBLIC :: nbondj_bdy_b(:)
  INTEGER, PUBLIC :: npolj
  INTEGER, PUBLIC :: nlci, nldi, nlei
  INTEGER, PUBLIC :: nlcj, nldj, nlej
  INTEGER, PUBLIC :: noea, nowe
  INTEGER, PUBLIC :: noso, nono
  INTEGER, PUBLIC :: nidom
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: mig
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: mjg
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: mi0, mi1
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: mj0, mj1
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: nimppt, njmppt
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: ibonit, ibonjt
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: nlcit, nlcjt
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: nldit, nldjt
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: nleit, nlejt
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: nfiimpp, nfipproc, nfilcit
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: glamt, glamu, glamv, glamf
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: gphit, gphiu, gphiv, gphif
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:, :) :: e1t, e2t, r1_e1t, r1_e2t
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:, :) :: e1u, e2u, r1_e1u, r1_e2u
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:, :) :: e1v, e2v, r1_e1v, r1_e2v
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:, :) :: e1f, e2f, r1_e1f, r1_e2f
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: e1e2t, r1_e1e2t
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: e1e2u, r1_e1e2u, e2_e1u
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: e1e2v, r1_e1e2v, e1_e2v
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: e1e2f, r1_e1e2f
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ff_f, ff_t
  LOGICAL, PUBLIC :: ln_zco
  LOGICAL, PUBLIC :: ln_zps
  LOGICAL, PUBLIC :: ln_sco
  LOGICAL, PUBLIC :: ln_isfcav
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: e3t_0, e3t_b, e3t_n, e3t_a
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: e3u_0, e3u_b, e3u_n, e3u_a
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: e3v_0, e3v_b, e3v_n, e3v_a
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: e3f_0, e3f_n
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: e3w_0, e3w_b, e3w_n
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: e3uw_0, e3uw_b, e3uw_n
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: e3vw_0, e3vw_b, e3vw_n
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: gdept_0, gdept_b, gdept_n
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: gdepw_0, gdepw_b, gdepw_n
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: gde3w_0, gde3w_n
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ht_0, ht_n
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hu_0, hu_b, hu_n, hu_a
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hv_0, hv_b, hv_n, hv_a
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: r1_hu_b, r1_hu_n, r1_hu_a
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: r1_hv_b, r1_hv_n, r1_hv_a
  INTEGER, PUBLIC :: nla10
  INTEGER, PUBLIC :: nlb10
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: gdept_1d, gdepw_1d
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: e3t_1d, e3w_1d
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: mbkt, mbku, mbkv
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: tmask_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: tmask_h
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: misfdep
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: mikt, miku, mikv, mikf
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: risfdep
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ssmask, ssumask, ssvmask
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :), TARGET :: tmask, umask, vmask, fmask
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :), TARGET :: wmask, wumask, wvmask
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: tpol, fpol
  INTEGER, PUBLIC :: nyear
  INTEGER, PUBLIC :: nmonth
  INTEGER, PUBLIC :: nday
  INTEGER, PUBLIC :: nhour
  INTEGER, PUBLIC :: nminute
  INTEGER, PUBLIC :: ndastp
  INTEGER, PUBLIC :: nday_year
  INTEGER, PUBLIC :: nsec_year
  INTEGER, PUBLIC :: nsec_month
  INTEGER, PUBLIC :: nsec_week
  INTEGER, PUBLIC :: nsec_day
  REAL(KIND = wp), PUBLIC :: fjulday
  REAL(KIND = wp), PUBLIC :: fjulstartyear
  REAL(KIND = wp), PUBLIC :: adatrj
  INTEGER, PUBLIC, DIMENSION(0 : 2) :: nyear_len
  INTEGER, PUBLIC, DIMENSION(0 : 13) :: nmonth_len
  INTEGER, PUBLIC, DIMENSION(0 : 13) :: nmonth_half
  INTEGER, PUBLIC, DIMENSION(0 : 13) :: nmonth_end
  INTEGER, PUBLIC :: nsec1jan000
  LOGICAL, PUBLIC, PARAMETER :: lk_agrif = .FALSE.
  CONTAINS
  LOGICAL FUNCTION Agrif_Root()
    Agrif_Root = .TRUE.
  END FUNCTION Agrif_Root
  CHARACTER(LEN = 3) FUNCTION Agrif_CFixed()
    Agrif_CFixed = '0'
  END FUNCTION Agrif_CFixed
  INTEGER FUNCTION dom_oce_alloc()
    INTEGER, DIMENSION(12) :: ierr
    !$ACC KERNELS
    ierr(:) = 0
    !$ACC END KERNELS
    ALLOCATE(mig(jpi), mjg(jpj), STAT = ierr(1))
    ALLOCATE(mi0(jpiglo), mi1(jpiglo), mj0(jpjglo), mj1(jpjglo), tpol(jpiglo), fpol(jpiglo), STAT = ierr(2))
    ALLOCATE(glamt(jpi, jpj), glamu(jpi, jpj), glamv(jpi, jpj), glamf(jpi, jpj), gphit(jpi, jpj), gphiu(jpi, jpj), gphiv(jpi, &
&jpj), gphif(jpi, jpj), e1t(jpi, jpj), e2t(jpi, jpj), r1_e1t(jpi, jpj), r1_e2t(jpi, jpj), e1u(jpi, jpj), e2u(jpi, jpj), &
&r1_e1u(jpi, jpj), r1_e2u(jpi, jpj), e1v(jpi, jpj), e2v(jpi, jpj), r1_e1v(jpi, jpj), r1_e2v(jpi, jpj), e1f(jpi, jpj), e2f(jpi, &
&jpj), r1_e1f(jpi, jpj), r1_e2f(jpi, jpj), e1e2t(jpi, jpj), r1_e1e2t(jpi, jpj), e1e2u(jpi, jpj), r1_e1e2u(jpi, jpj), e2_e1u(jpi, &
&jpj), e1e2v(jpi, jpj), r1_e1e2v(jpi, jpj), e1_e2v(jpi, jpj), e1e2f(jpi, jpj), r1_e1e2f(jpi, jpj), ff_f(jpi, jpj), ff_t(jpi, jpj), &
&STAT = ierr(3))
    ALLOCATE(gdept_0(jpi, jpj, jpk), gdepw_0(jpi, jpj, jpk), gde3w_0(jpi, jpj, jpk), gdept_b(jpi, jpj, jpk), gdepw_b(jpi, jpj, &
&jpk), gdept_n(jpi, jpj, jpk), gdepw_n(jpi, jpj, jpk), gde3w_n(jpi, jpj, jpk), STAT = ierr(4))
    ALLOCATE(e3t_0(jpi, jpj, jpk), e3u_0(jpi, jpj, jpk), e3v_0(jpi, jpj, jpk), e3f_0(jpi, jpj, jpk), e3w_0(jpi, jpj, jpk), &
&e3t_b(jpi, jpj, jpk), e3u_b(jpi, jpj, jpk), e3v_b(jpi, jpj, jpk), e3w_b(jpi, jpj, jpk), e3t_n(jpi, jpj, jpk), e3u_n(jpi, jpj, &
&jpk), e3v_n(jpi, jpj, jpk), e3f_n(jpi, jpj, jpk), e3w_n(jpi, jpj, jpk), e3t_a(jpi, jpj, jpk), e3u_a(jpi, jpj, jpk), e3v_a(jpi, &
&jpj, jpk), e3uw_0(jpi, jpj, jpk), e3vw_0(jpi, jpj, jpk), e3uw_b(jpi, jpj, jpk), e3vw_b(jpi, jpj, jpk), e3uw_n(jpi, jpj, jpk), &
&e3vw_n(jpi, jpj, jpk), STAT = ierr(5))
    ALLOCATE(ht_0(jpi, jpj), hu_0(jpi, jpj), hv_0(jpi, jpj), hu_b(jpi, jpj), hv_b(jpi, jpj), r1_hu_b(jpi, jpj), r1_hv_b(jpi, jpj), &
&ht_n(jpi, jpj), hu_n(jpi, jpj), hv_n(jpi, jpj), r1_hu_n(jpi, jpj), r1_hv_n(jpi, jpj), hu_a(jpi, jpj), hv_a(jpi, jpj), &
&r1_hu_a(jpi, jpj), r1_hv_a(jpi, jpj), STAT = ierr(6))
    ALLOCATE(gdept_1d(jpk), gdepw_1d(jpk), e3t_1d(jpk), e3w_1d(jpk), STAT = ierr(7))
    ALLOCATE(tmask_i(jpi, jpj), tmask_h(jpi, jpj), ssmask(jpi, jpj), ssumask(jpi, jpj), ssvmask(jpi, jpj), mbkt(jpi, jpj), &
&mbku(jpi, jpj), mbkv(jpi, jpj), STAT = ierr(9))
    ALLOCATE(misfdep(jpi, jpj), mikt(jpi, jpj), miku(jpi, jpj), risfdep(jpi, jpj), mikv(jpi, jpj), mikf(jpi, jpj), STAT = ierr(10))
    ALLOCATE(tmask(jpi, jpj, jpk), umask(jpi, jpj, jpk), vmask(jpi, jpj, jpk), fmask(jpi, jpj, jpk), STAT = ierr(11))
    ALLOCATE(wmask(jpi, jpj, jpk), wumask(jpi, jpj, jpk), wvmask(jpi, jpj, jpk), STAT = ierr(12))
    dom_oce_alloc = MAXVAL(ierr)
  END FUNCTION dom_oce_alloc
END MODULE dom_oce