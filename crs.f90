MODULE crs
  USE par_oce
  USE dom_oce
  USE in_out_manager
  IMPLICIT NONE
  PUBLIC
  PUBLIC :: crs_dom_alloc
  PUBLIC :: crs_dom_alloc2
  PUBLIC :: dom_grid_glo
  PUBLIC :: dom_grid_crs
  INTEGER :: jpiglo_crs, jpjglo_crs
  INTEGER :: jpi_crs, jpj_crs
  INTEGER :: jpi_full, jpj_full
  INTEGER :: nistr, njstr
  INTEGER :: niend, njend
  INTEGER :: jpi_crsm1, jpj_crsm1
  INTEGER :: jpiglo_crsm1, jpjglo_crsm1
  INTEGER :: nperio_full, nperio_crs
  INTEGER :: npolj_full, npolj_crs
  INTEGER :: jpiglo_full, jpjglo_full
  INTEGER :: npiglo, npjglo
  INTEGER :: nlci_full, nlcj_full
  INTEGER :: nldi_full, nldj_full
  INTEGER :: nlei_full, nlej_full
  INTEGER :: nlci_crs, nlcj_crs
  INTEGER :: nldi_crs, nldj_crs
  INTEGER :: nlei_crs, nlej_crs
  INTEGER :: narea_full, narea_crs
  INTEGER :: jpnij_full, jpnij_crs
  INTEGER :: jpim1_full, jpjm1_full
  INTEGER :: nimpp_full, njmpp_full
  INTEGER :: nimpp_crs, njmpp_crs
  INTEGER :: nreci_full, nrecj_full
  INTEGER :: nreci_crs, nrecj_crs
  INTEGER :: noea_full, nowe_full
  INTEGER :: noso_full, nono_full
  INTEGER :: npne_full, npnw_full
  INTEGER :: npse_full, npsw_full
  INTEGER :: nbne_full, nbnw_full
  INTEGER :: nbse_full, nbsw_full
  INTEGER :: nidom_full
  INTEGER :: nproc_full
  INTEGER :: nbondi_full, nbondj_full
  INTEGER :: noea_crs, nowe_crs
  INTEGER :: noso_crs, nono_crs
  INTEGER :: npne_crs, npnw_crs
  INTEGER :: npse_crs, npsw_crs
  INTEGER :: nbne_crs, nbnw_crs
  INTEGER :: nbse_crs, nbsw_crs
  INTEGER :: nidom_crs
  INTEGER :: nproc_crs
  INTEGER :: nbondi_crs, nbondj_crs
  INTEGER, DIMENSION(:), ALLOCATABLE :: mis_crs, mie_crs, mis2_crs, mie2_crs
  INTEGER, DIMENSION(:), ALLOCATABLE :: mjs_crs, mje_crs, mjs2_crs, mje2_crs
  INTEGER, DIMENSION(:), ALLOCATABLE :: mjg_crs, mig_crs
  INTEGER, DIMENSION(:), ALLOCATABLE :: mi0_crs, mi1_crs, mj0_crs, mj1_crs
  INTEGER :: mxbinctr, mybinctr
  INTEGER, DIMENSION(:), ALLOCATABLE :: nlcit_crs, nlcit_full
  INTEGER, DIMENSION(:), ALLOCATABLE :: nldit_crs, nldit_full
  INTEGER, DIMENSION(:), ALLOCATABLE :: nleit_crs, nleit_full
  INTEGER, DIMENSION(:), ALLOCATABLE :: nimppt_crs, nimppt_full
  INTEGER, DIMENSION(:), ALLOCATABLE :: nlcjt_crs, nlcjt_full
  INTEGER, DIMENSION(:), ALLOCATABLE :: nldjt_crs, nldjt_full
  INTEGER, DIMENSION(:), ALLOCATABLE :: nlejt_crs, nlejt_full
  INTEGER, DIMENSION(:), ALLOCATABLE :: njmppt_crs, njmppt_full
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: tmask_crs, umask_crs, vmask_crs, fmask_crs
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: tmask_i_crs, rnfmsk_crs, tpol_crs, fpol_crs
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: e1t_crs, e2t_crs, e1e2t_crs
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: e1u_crs, e2u_crs
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: e1v_crs, e2v_crs
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: e1f_crs, e2f_crs
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: e3t_crs, e3u_crs, e3v_crs, e3f_crs, e3w_crs
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: e3t_max_crs, e3u_max_crs, e3v_max_crs, e3f_max_crs, e3w_max_crs
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: e1e2w_crs, e2e3u_crs, e1e3v_crs
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: e1e2w_msk, e2e3u_msk, e1e3v_msk
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: gphit_crs, glamt_crs, gphif_crs, glamf_crs
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: gphiu_crs, glamu_crs, gphiv_crs, glamv_crs
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: ff_crs
  INTEGER, DIMENSION(:, :), ALLOCATABLE :: mbathy_crs, mbkt_crs, mbku_crs, mbkv_crs
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: gdept_crs, gdepu_crs, gdepv_crs, gdepw_crs
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: facsurfv, facsurfu, facvol_t, facvol_w
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ocean_volume_crs_t, ocean_volume_crs_w, bt_crs, r1_bt_crs
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: crs_surfu_wgt, crs_surfv_wgt, crs_surfw_wgt, crs_volt_wgt
  INTEGER :: nn_factx = 3
  INTEGER :: nn_facty = 3
  INTEGER :: nn_binref = 0
  LOGICAL :: ln_msh_crs = .TRUE.
  INTEGER :: nn_crs_kz = 0
  LOGICAL :: ln_crs_wn = .FALSE.
  INTEGER :: nrestx, nresty
  REAL(KIND = wp) :: rfactx_r
  REAL(KIND = wp) :: rfacty_r
  REAL(KIND = wp) :: rfactxy
  REAL(KIND = wp), DIMENSION(:, :, :, :), ALLOCATABLE :: tsn_crs
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: un_crs, vn_crs, wn_crs
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: hdivn_crs
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: sshn_crs
  REAL(KIND = wp), PUBLIC, DIMENSION(:, :), ALLOCATABLE :: qsr_crs, fr_i_crs, wndm_crs
  REAL(KIND = wp), PUBLIC, DIMENSION(:, :), ALLOCATABLE :: emp_crs, emp_b_crs, sfx_crs
  REAL(KIND = wp), PUBLIC, DIMENSION(:, :), ALLOCATABLE :: utau_crs, vtau_crs
  REAL(KIND = wp), PUBLIC, DIMENSION(:, :), ALLOCATABLE :: rnf_crs
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :, :) :: avt_crs
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :, :) :: avs_crs
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: nmln_crs, hmld_crs, hmlp_crs, hmlpt_crs
  CONTAINS
  INTEGER FUNCTION crs_dom_alloc()
    INTEGER, DIMENSION(17) :: ierr
    ierr(:) = 0
    ALLOCATE(mis2_crs(jpiglo_crs), mie2_crs(jpiglo_crs), mjs2_crs(jpjglo_crs), mje2_crs(jpjglo_crs), mi0_crs(jpiglo_crs), &
&mi1_crs(jpiglo_crs), mj0_crs(jpjglo_crs), mj1_crs(jpjglo_crs), mig_crs(jpi_crs), mjg_crs(jpj_crs), STAT = ierr(1))
    ALLOCATE(tmask_crs(jpi_crs, jpj_crs, jpk), fmask_crs(jpi_crs, jpj_crs, jpk), umask_crs(jpi_crs, jpj_crs, jpk), &
&vmask_crs(jpi_crs, jpj_crs, jpk), STAT = ierr(2))
    ALLOCATE(tmask_i_crs(jpi_crs, jpj_crs), rnfmsk_crs(jpi_crs, jpj_crs), tpol_crs(jpiglo_crs, jpjglo_crs), fpol_crs(jpiglo_crs, &
&jpjglo_crs), STAT = ierr(3))
    ALLOCATE(gphit_crs(jpi_crs, jpj_crs), glamt_crs(jpi_crs, jpj_crs), gphiu_crs(jpi_crs, jpj_crs), glamu_crs(jpi_crs, jpj_crs), &
&gphiv_crs(jpi_crs, jpj_crs), glamv_crs(jpi_crs, jpj_crs), gphif_crs(jpi_crs, jpj_crs), glamf_crs(jpi_crs, jpj_crs), &
&ff_crs(jpi_crs, jpj_crs), STAT = ierr(4))
    ALLOCATE(e1t_crs(jpi_crs, jpj_crs), e2t_crs(jpi_crs, jpj_crs), e1u_crs(jpi_crs, jpj_crs), e2u_crs(jpi_crs, jpj_crs), &
&e1v_crs(jpi_crs, jpj_crs), e2v_crs(jpi_crs, jpj_crs), e1f_crs(jpi_crs, jpj_crs), e2f_crs(jpi_crs, jpj_crs), e1e2t_crs(jpi_crs, &
&jpj_crs), STAT = ierr(5))
    ALLOCATE(e3t_crs(jpi_crs, jpj_crs, jpk), e3w_crs(jpi_crs, jpj_crs, jpk), e3u_crs(jpi_crs, jpj_crs, jpk), e3v_crs(jpi_crs, &
&jpj_crs, jpk), e3f_crs(jpi_crs, jpj_crs, jpk), e1e2w_msk(jpi_crs, jpj_crs, jpk), e2e3u_msk(jpi_crs, jpj_crs, jpk), &
&e1e3v_msk(jpi_crs, jpj_crs, jpk), e1e2w_crs(jpi_crs, jpj_crs, jpk), e2e3u_crs(jpi_crs, jpj_crs, jpk), e1e3v_crs(jpi_crs, jpj_crs, &
&jpk), e3t_max_crs(jpi_crs, jpj_crs, jpk), e3w_max_crs(jpi_crs, jpj_crs, jpk), e3u_max_crs(jpi_crs, jpj_crs, jpk), &
&e3v_max_crs(jpi_crs, jpj_crs, jpk), STAT = ierr(6))
    ALLOCATE(facsurfv(jpi_crs, jpj_crs, jpk), facsurfu(jpi_crs, jpj_crs, jpk), facvol_t(jpi_crs, jpj_crs, jpk), facvol_w(jpi_crs, &
&jpj_crs, jpk), ocean_volume_crs_t(jpi_crs, jpj_crs, jpk), ocean_volume_crs_w(jpi_crs, jpj_crs, jpk), bt_crs(jpi_crs, jpj_crs, &
&jpk), r1_bt_crs(jpi_crs, jpj_crs, jpk), STAT = ierr(7))
    ALLOCATE(crs_surfu_wgt(jpi_crs, jpj_crs, jpk), crs_surfv_wgt(jpi_crs, jpj_crs, jpk), crs_surfw_wgt(jpi_crs, jpj_crs, jpk), &
&crs_volt_wgt(jpi_crs, jpj_crs, jpk), STAT = ierr(8))
    ALLOCATE(mbathy_crs(jpi_crs, jpj_crs), mbkt_crs(jpi_crs, jpj_crs), mbku_crs(jpi_crs, jpj_crs), mbkv_crs(jpi_crs, jpj_crs), &
&STAT = ierr(9))
    ALLOCATE(gdept_crs(jpi_crs, jpj_crs, jpk), gdepu_crs(jpi_crs, jpj_crs, jpk), gdepv_crs(jpi_crs, jpj_crs, jpk), &
&gdepw_crs(jpi_crs, jpj_crs, jpk), STAT = ierr(10))
    ALLOCATE(un_crs(jpi_crs, jpj_crs, jpk), vn_crs(jpi_crs, jpj_crs, jpk), wn_crs(jpi_crs, jpj_crs, jpk), hdivn_crs(jpi_crs, &
&jpj_crs, jpk), STAT = ierr(11))
    ALLOCATE(sshn_crs(jpi_crs, jpj_crs), emp_crs(jpi_crs, jpj_crs), emp_b_crs(jpi_crs, jpj_crs), qsr_crs(jpi_crs, jpj_crs), &
&wndm_crs(jpi_crs, jpj_crs), utau_crs(jpi_crs, jpj_crs), vtau_crs(jpi_crs, jpj_crs), rnf_crs(jpi_crs, jpj_crs), fr_i_crs(jpi_crs, &
&jpj_crs), sfx_crs(jpi_crs, jpj_crs), STAT = ierr(12))
    ALLOCATE(tsn_crs(jpi_crs, jpj_crs, jpk, jpts), avt_crs(jpi_crs, jpj_crs, jpk), avs_crs(jpi_crs, jpj_crs, jpk), STAT = ierr(13))
    ALLOCATE(nmln_crs(jpi_crs, jpj_crs), hmld_crs(jpi_crs, jpj_crs), hmlp_crs(jpi_crs, jpj_crs), hmlpt_crs(jpi_crs, jpj_crs), STAT &
&= ierr(14))
    ALLOCATE(nimppt_crs(jpnij), nlcit_crs(jpnij), nldit_crs(jpnij), nleit_crs(jpnij), nimppt_full(jpnij), nlcit_full(jpnij), &
&nldit_full(jpnij), nleit_full(jpnij), njmppt_crs(jpnij), nlcjt_crs(jpnij), nldjt_crs(jpnij), nlejt_crs(jpnij), &
&njmppt_full(jpnij), nlcjt_full(jpnij), nldjt_full(jpnij), nlejt_full(jpnij), STAT = ierr(15))
    crs_dom_alloc = MAXVAL(ierr)
  END FUNCTION crs_dom_alloc
  INTEGER FUNCTION crs_dom_alloc2()
    INTEGER, DIMENSION(1) :: ierr
    ierr(:) = 0
    ALLOCATE(mjs_crs(nlej_crs), mje_crs(nlej_crs), mis_crs(nlei_crs), mie_crs(nlei_crs), STAT = ierr(1))
    crs_dom_alloc2 = MAXVAL(ierr)
  END FUNCTION crs_dom_alloc2
  SUBROUTINE dom_grid_glo
    jpi = jpi_full
    jpj = jpj_full
    jpim1 = jpim1_full
    jpjm1 = jpjm1_full
    jperio = nperio_full
    npolj = npolj_full
    jpiglo = jpiglo_full
    jpjglo = jpjglo_full
    nlci = nlci_full
    nlcj = nlcj_full
    nldi = nldi_full
    nldj = nldj_full
    nlei = nlei_full
    nlej = nlej_full
    nimpp = nimpp_full
    njmpp = njmpp_full
    nlcit(:) = nlcit_full(:)
    nldit(:) = nldit_full(:)
    nleit(:) = nleit_full(:)
    nimppt(:) = nimppt_full(:)
    nlcjt(:) = nlcjt_full(:)
    nldjt(:) = nldjt_full(:)
    nlejt(:) = nlejt_full(:)
    njmppt(:) = njmppt_full(:)
  END SUBROUTINE dom_grid_glo
  SUBROUTINE dom_grid_crs
    jpi = jpi_crs
    jpj = jpj_crs
    jpim1 = jpi_crsm1
    jpjm1 = jpj_crsm1
    jperio = nperio_crs
    npolj = npolj_crs
    jpiglo = jpiglo_crs
    jpjglo = jpjglo_crs
    nlci = nlci_crs
    nlcj = nlcj_crs
    nldi = nldi_crs
    nlei = nlei_crs
    nlej = nlej_crs
    nldj = nldj_crs
    nimpp = nimpp_crs
    njmpp = njmpp_crs
    nlcit(:) = nlcit_crs(:)
    nldit(:) = nldit_crs(:)
    nleit(:) = nleit_crs(:)
    nimppt(:) = nimppt_crs(:)
    nlcjt(:) = nlcjt_crs(:)
    nldjt(:) = nldjt_crs(:)
    nlejt(:) = nlejt_crs(:)
    njmppt(:) = njmppt_crs(:)
  END SUBROUTINE dom_grid_crs
END MODULE crs