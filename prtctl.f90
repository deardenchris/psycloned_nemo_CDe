MODULE prtctl
  USE dom_oce
  USE in_out_manager
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: numid
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: nlditl, nldjtl
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: nleitl, nlejtl
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: nimpptl, njmpptl
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: nlcitl, nlcjtl
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ibonitl, ibonjtl
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE, SAVE :: t_ctll, s_ctll
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE, SAVE :: u_ctll, v_ctll
  INTEGER :: ktime
  PUBLIC :: prt_ctl
  PUBLIC :: prt_ctl_info
  PUBLIC :: prt_ctl_init
  PUBLIC :: sub_dom
  CONTAINS
  SUBROUTINE prt_ctl(tab2d_1, tab3d_1, mask1, clinfo1, tab2d_2, tab3d_2, mask2, clinfo2, kdim, clinfo3)
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN), OPTIONAL :: tab2d_1
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN), OPTIONAL :: tab3d_1
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN), OPTIONAL :: mask1
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: clinfo1
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN), OPTIONAL :: tab2d_2
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN), OPTIONAL :: tab3d_2
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN), OPTIONAL :: mask2
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: clinfo2
    INTEGER, INTENT(IN), OPTIONAL :: kdim
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: clinfo3
    CHARACTER(LEN = 15) :: cl2
    INTEGER :: jn, sind, eind, kdir, j_id
    REAL(KIND = wp) :: zsum1, zsum2, zvctl1, zvctl2
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztab2d_1, ztab2d_2
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zmask1, zmask2, ztab3d_1, ztab3d_2
    kdir = jpkm1
    cl2 = ''
    zsum1 = 0.E0
    zsum2 = 0.E0
    zvctl1 = 0.E0
    zvctl2 = 0.E0
    ztab2d_1(:, :) = 0.E0
    ztab2d_2(:, :) = 0.E0
    ztab3d_1(:, :, :) = 0.E0
    ztab3d_2(:, :, :) = 0.E0
    zmask1(:, :, :) = 1.E0
    zmask2(:, :, :) = 1.E0
    IF (PRESENT(clinfo2)) cl2 = clinfo2
    IF (PRESENT(kdim)) kdir = kdim
    IF (PRESENT(tab2d_1)) ztab2d_1(:, :) = tab2d_1(:, :)
    IF (PRESENT(tab2d_2)) ztab2d_2(:, :) = tab2d_2(:, :)
    IF (PRESENT(tab3d_1)) ztab3d_1(:, :, 1 : kdir) = tab3d_1(:, :, 1 : kdir)
    IF (PRESENT(tab3d_2)) ztab3d_2(:, :, 1 : kdir) = tab3d_2(:, :, 1 : kdir)
    IF (PRESENT(mask1)) zmask1(:, :, :) = mask1(:, :, :)
    IF (PRESENT(mask2)) zmask2(:, :, :) = mask2(:, :, :)
    IF (lk_mpp .AND. jpnij > 1) THEN
      sind = narea
      eind = narea
    ELSE
      sind = 1
      eind = ijsplt
    END IF
    DO jn = sind, eind
      j_id = numid(jn - narea + 1)
      IF (.NOT. lsp_area) THEN
        IF (lk_mpp .AND. jpnij > 1) THEN
          nictls = MAX(1, nlditl(jn))
          nictle = MIN(jpi, nleitl(jn))
          njctls = MAX(1, nldjtl(jn))
          njctle = MIN(jpj, nlejtl(jn))
          IF (ibonitl(jn) == - 1 .OR. ibonitl(jn) == 2) nictls = MAX(2, nictls)
          IF (ibonjtl(jn) == - 1 .OR. ibonjtl(jn) == 2) njctls = MAX(2, njctls)
          IF (ibonitl(jn) == 1 .OR. ibonitl(jn) == 2) nictle = MIN(nictle, nleitl(jn) - 1)
          IF (ibonjtl(jn) == 1 .OR. ibonjtl(jn) == 2) njctle = MIN(njctle, nlejtl(jn) - 1)
        ELSE
          nictls = MAX(1, nimpptl(jn) - 1 + nlditl(jn))
          nictle = MIN(jpi, nimpptl(jn) - 1 + nleitl(jn))
          njctls = MAX(1, njmpptl(jn) - 1 + nldjtl(jn))
          njctle = MIN(jpj, njmpptl(jn) - 1 + nlejtl(jn))
          IF (ibonitl(jn) == - 1 .OR. ibonitl(jn) == 2) nictls = MAX(2, nictls)
          IF (ibonjtl(jn) == - 1 .OR. ibonjtl(jn) == 2) njctls = MAX(2, njctls)
          IF (ibonitl(jn) == 1 .OR. ibonitl(jn) == 2) nictle = MIN(nictle, nimpptl(jn) + nleitl(jn) - 2)
          IF (ibonjtl(jn) == 1 .OR. ibonjtl(jn) == 2) njctle = MIN(njctle, njmpptl(jn) + nlejtl(jn) - 2)
        END IF
      END IF
      IF (PRESENT(clinfo3)) THEN
        IF (clinfo3 == 'tra') THEN
          zvctl1 = t_ctll(jn)
          zvctl2 = s_ctll(jn)
        ELSE IF (clinfo3 == 'dyn') THEN
          zvctl1 = u_ctll(jn)
          zvctl2 = v_ctll(jn)
        END IF
      END IF
      IF (PRESENT(tab2d_1)) THEN
        zsum1 = SUM(ztab2d_1(nictls : nictle, njctls : njctle) * zmask1(nictls : nictle, njctls : njctle, 1))
        zsum2 = SUM(ztab2d_2(nictls : nictle, njctls : njctle) * zmask2(nictls : nictle, njctls : njctle, 1))
      END IF
      IF (PRESENT(tab3d_1)) THEN
        zsum1 = SUM(ztab3d_1(nictls : nictle, njctls : njctle, 1 : kdir) * zmask1(nictls : nictle, njctls : njctle, 1 : kdir))
        zsum2 = SUM(ztab3d_2(nictls : nictle, njctls : njctle, 1 : kdir) * zmask2(nictls : nictle, njctls : njctle, 1 : kdir))
      END IF
      IF (PRESENT(clinfo3)) THEN
        WRITE(j_id, FMT = '(a,D23.16,3x,a,D23.16)') clinfo1, zsum1 - zvctl1, cl2, zsum2 - zvctl2
        SELECT CASE (clinfo3)
        CASE ('tra-ta')
          t_ctll(jn) = zsum1
        CASE ('tra')
          t_ctll(jn) = zsum1
          s_ctll(jn) = zsum2
        CASE ('dyn')
          u_ctll(jn) = zsum1
          v_ctll(jn) = zsum2
        END SELECT
      ELSE IF (PRESENT(clinfo2) .OR. PRESENT(tab2d_2) .OR. PRESENT(tab3d_2)) THEN
        WRITE(j_id, FMT = '(a,D23.16,3x,a,D23.16)') clinfo1, zsum1, cl2, zsum2
      ELSE
        WRITE(j_id, FMT = '(a,D23.16)') clinfo1, zsum1
      END IF
    END DO
  END SUBROUTINE prt_ctl
  SUBROUTINE prt_ctl_info(clinfo1, ivar1, clinfo2, ivar2, itime)
    CHARACTER(LEN = *), INTENT(IN) :: clinfo1
    INTEGER, INTENT(IN), OPTIONAL :: ivar1
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: clinfo2
    INTEGER, INTENT(IN), OPTIONAL :: ivar2
    INTEGER, INTENT(IN), OPTIONAL :: itime
    INTEGER :: jn, sind, eind, iltime, j_id
    IF (lk_mpp .AND. jpnij > 1) THEN
      sind = narea
      eind = narea
    ELSE
      sind = 1
      eind = ijsplt
    END IF
    IF (PRESENT(itime)) THEN
      iltime = itime
      IF (iltime > ktime) THEN
        t_ctll(:) = 0.E0
        s_ctll(:) = 0.E0
        u_ctll(:) = 0.E0
        v_ctll(:) = 0.E0
        ktime = iltime
      END IF
    END IF
    DO jn = sind, eind
      j_id = numid(jn - narea + 1)
      IF (PRESENT(ivar1) .AND. PRESENT(clinfo2) .AND. PRESENT(ivar2)) THEN
        WRITE(j_id, FMT = *) clinfo1, ivar1, clinfo2, ivar2
      ELSE IF (PRESENT(ivar1) .AND. PRESENT(clinfo2) .AND. .NOT. PRESENT(ivar2)) THEN
        WRITE(j_id, FMT = *) clinfo1, ivar1, clinfo2
      ELSE IF (PRESENT(ivar1) .AND. .NOT. PRESENT(clinfo2) .AND. PRESENT(ivar2)) THEN
        WRITE(j_id, FMT = *) clinfo1, ivar1, ivar2
      ELSE IF (PRESENT(ivar1) .AND. .NOT. PRESENT(clinfo2) .AND. .NOT. PRESENT(ivar2)) THEN
        WRITE(j_id, FMT = *) clinfo1, ivar1
      ELSE
        WRITE(j_id, FMT = *) clinfo1
      END IF
    END DO
  END SUBROUTINE prt_ctl_info
  SUBROUTINE prt_ctl_init
    INTEGER :: jn, sind, eind, j_id
    CHARACTER(LEN = 28) :: clfile_out
    CHARACTER(LEN = 23) :: clb_name
    CHARACTER(LEN = 19) :: cl_run
    ALLOCATE(nlditl(ijsplt), nleitl(ijsplt), nimpptl(ijsplt), ibonitl(ijsplt), nldjtl(ijsplt), nlejtl(ijsplt), njmpptl(ijsplt), &
&ibonjtl(ijsplt), nlcitl(ijsplt), t_ctll(ijsplt), u_ctll(ijsplt), nlcjtl(ijsplt), s_ctll(ijsplt), v_ctll(ijsplt))
    t_ctll(:) = 0.E0
    s_ctll(:) = 0.E0
    u_ctll(:) = 0.E0
    v_ctll(:) = 0.E0
    ktime = 1
    IF (lk_mpp .AND. jpnij > 1) THEN
      sind = narea
      eind = narea
      clb_name = "('mpp.output_',I4.4)"
      cl_run = 'MULTI processor run'
      nlditl(1 : jpnij) = nldit(:)
      nleitl(1 : jpnij) = nleit(:)
      nldjtl(1 : jpnij) = nldjt(:)
      nlejtl(1 : jpnij) = nlejt(:)
      nimpptl(1 : jpnij) = nimppt(:)
      njmpptl(1 : jpnij) = njmppt(:)
      nlcitl(1 : jpnij) = nlcit(:)
      nlcjtl(1 : jpnij) = nlcjt(:)
      ibonitl(1 : jpnij) = ibonit(:)
      ibonjtl(1 : jpnij) = ibonjt(:)
    ELSE
      sind = 1
      eind = ijsplt
      clb_name = "('mono.output_',I4.4)"
      cl_run = 'MONO processor run '
      CALL sub_dom
    END IF
    ALLOCATE(numid(eind - sind + 1))
    DO jn = sind, eind
      WRITE(clfile_out, FMT = clb_name) jn - 1
      CALL ctl_opn(numid(jn - narea + 1), clfile_out, 'REPLACE', 'FORMATTED', 'SEQUENTIAL', 1, numout, .FALSE.)
      j_id = numid(jn - narea + 1)
      WRITE(j_id, FMT = *)
      WRITE(j_id, FMT = *) '                 L O D Y C - I P S L'
      WRITE(j_id, FMT = *) '                     O P A model'
      WRITE(j_id, FMT = *) '            Ocean General Circulation Model'
      WRITE(j_id, FMT = *) '               version OPA 9.0  (2005) '
      WRITE(j_id, FMT = *)
      WRITE(j_id, FMT = *) '                   PROC number: ', jn
      WRITE(j_id, FMT = *)
      WRITE(j_id, FMT = "(19x,a20)") cl_run
      IF (.NOT. lsp_area) THEN
        nictls = nimpptl(jn) + nlditl(jn) - 1
        nictle = nimpptl(jn) + nleitl(jn) - 1
        njctls = njmpptl(jn) + nldjtl(jn) - 1
        njctle = njmpptl(jn) + nlejtl(jn) - 1
      END IF
      WRITE(j_id, FMT = *)
      WRITE(j_id, FMT = *) 'prt_ctl :  Sum control indices'
      WRITE(j_id, FMT = *) '~~~~~~~'
      WRITE(j_id, FMT = *)
      WRITE(j_id, 9000) '                                nlej   = ', nlejtl(jn), '              '
      WRITE(j_id, 9000) '                  ------------- njctle = ', njctle, ' -------------'
      WRITE(j_id, 9001) '                  |                                       |'
      WRITE(j_id, 9001) '                  |                                       |'
      WRITE(j_id, 9001) '                  |                                       |'
      WRITE(j_id, 9002) '           nictls = ', nictls, '                           nictle = ', nictle
      WRITE(j_id, 9002) '           nldi   = ', nlditl(jn), '                           nlei   = ', nleitl(jn)
      WRITE(j_id, 9001) '                  |                                       |'
      WRITE(j_id, 9001) '                  |                                       |'
      WRITE(j_id, 9001) '                  |                                       |'
      WRITE(j_id, 9004) '  njmpp  = ', njmpptl(jn), '   ------------- njctls = ', njctls, ' -------------'
      WRITE(j_id, 9003) '           nimpp  = ', nimpptl(jn), '        nldj   = ', nldjtl(jn), '              '
      WRITE(j_id, FMT = *)
      WRITE(j_id, FMT = *)
9000  FORMAT(A41, I4.4, A14)
9001  FORMAT(A59)
9002  FORMAT(A20, I4.4, A36, I3.3)
9003  FORMAT(A20, I4.4, A17, I4.4)
9004  FORMAT(A11, I4.4, A26, I4.4, A14)
    END DO
  END SUBROUTINE prt_ctl_init
  SUBROUTINE sub_dom
    INTEGER :: ji, jj, jn
    INTEGER :: ii, ij, irestil, irestjl, ijpi, ijpj, nlcil, nlcjl, nbondil, nbondjl, nrecil, nrecjl, nldil, nleil, nldjl, nlejl
    INTEGER, DIMENSION(jpi, jpj) :: iimpptl, ijmpptl, ilcitl, ilcjtl
    REAL(KIND = wp) :: zidom, zjdom
    INTEGER :: inum
    ijpi = (jpiglo - 2 * nn_hls + (isplt - 1)) / isplt + 2 * nn_hls
    ijpj = (jpjglo - 2 * nn_hls + (jsplt - 1)) / jsplt + 2 * nn_hls
    nrecil = 2 * nn_hls
    nrecjl = 2 * nn_hls
    irestil = MOD(jpiglo - nrecil, isplt)
    irestjl = MOD(jpjglo - nrecjl, jsplt)
    IF (irestil == 0) irestil = isplt
    DO jj = 1, jsplt
      DO ji = 1, irestil
        ilcitl(ji, jj) = ijpi
      END DO
      DO ji = irestil + 1, isplt
        ilcitl(ji, jj) = ijpi - 1
      END DO
    END DO
    IF (irestjl == 0) irestjl = jsplt
    DO ji = 1, isplt
      DO jj = 1, irestjl
        ilcjtl(ji, jj) = ijpj
      END DO
      DO jj = irestjl + 1, jsplt
        ilcjtl(ji, jj) = ijpj - 1
      END DO
    END DO
    zidom = nrecil
    DO ji = 1, isplt
      zidom = zidom + ilcitl(ji, 1) - nrecil
    END DO
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) ' sum ilcitl(i,1) = ', zidom, ' jpiglo = ', jpiglo
    zjdom = nrecjl
    DO jj = 1, jsplt
      zjdom = zjdom + ilcjtl(1, jj) - nrecjl
    END DO
    IF (lwp) WRITE(numout, FMT = *) ' sum ilcitl(1,j) = ', zjdom, ' jpjglo = ', jpjglo
    IF (lwp) WRITE(numout, FMT = *)
    iimpptl(:, :) = 1
    ijmpptl(:, :) = 1
    IF (isplt > 1) THEN
      DO jj = 1, jsplt
        DO ji = 2, isplt
          iimpptl(ji, jj) = iimpptl(ji - 1, jj) + ilcitl(ji - 1, jj) - nrecil
        END DO
      END DO
    END IF
    IF (jsplt > 1) THEN
      DO jj = 2, jsplt
        DO ji = 1, isplt
          ijmpptl(ji, jj) = ijmpptl(ji, jj - 1) + ilcjtl(ji, jj - 1) - nrecjl
        END DO
      END DO
    END IF
    DO jn = 1, ijsplt
      ii = 1 + MOD(jn - 1, isplt)
      ij = 1 + (jn - 1) / isplt
      nimpptl(jn) = iimpptl(ii, ij)
      njmpptl(jn) = ijmpptl(ii, ij)
      nlcitl(jn) = ilcitl(ii, ij)
      nlcil = nlcitl(jn)
      nlcjtl(jn) = ilcjtl(ii, ij)
      nlcjl = nlcjtl(jn)
      nbondjl = - 1
      IF (jn > isplt) nbondjl = 0
      IF (jn > (jsplt - 1) * isplt) nbondjl = 1
      IF (jsplt == 1) nbondjl = 2
      ibonjtl(jn) = nbondjl
      nbondil = 0
      IF (MOD(jn, isplt) == 1) nbondil = - 1
      IF (MOD(jn, isplt) == 0) nbondil = 1
      IF (isplt == 1) nbondil = 2
      ibonitl(jn) = nbondil
      nldil = 1 + nn_hls
      nleil = nlcil - nn_hls
      IF (nbondil == - 1 .OR. nbondil == 2) nldil = 1
      IF (nbondil == 1 .OR. nbondil == 2) nleil = nlcil
      nldjl = 1 + nn_hls
      nlejl = nlcjl - nn_hls
      IF (nbondjl == - 1 .OR. nbondjl == 2) nldjl = 1
      IF (nbondjl == 1 .OR. nbondjl == 2) nlejl = nlcjl
      nlditl(jn) = nldil
      nleitl(jn) = nleil
      nldjtl(jn) = nldjl
      nlejtl(jn) = nlejl
    END DO
    IF (lwp) THEN
      CALL ctl_opn(inum, 'layout_prtctl.dat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, numout, .FALSE., narea)
      WRITE(inum, FMT = '(a)') 'nproc nlcil nlcjl nldil nldjl nleil nlejl nimpptl njmpptl ibonitl ibonjtl'
      DO jn = 1, ijsplt
        WRITE(inum, FMT = '(i5,6i6,4i8)') jn - 1, nlcitl(jn), nlcjtl(jn), nlditl(jn), nldjtl(jn), nleitl(jn), nlejtl(jn), &
&nimpptl(jn), njmpptl(jn), ibonitl(jn), ibonjtl(jn)
      END DO
      CLOSE(UNIT = inum)
    END IF
  END SUBROUTINE sub_dom
END MODULE prtctl