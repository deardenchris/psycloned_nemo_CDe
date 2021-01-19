MODULE mppini
  USE dom_oce
  USE bdy_oce
  USE lbcnfd, ONLY: isendto, nsndto, nfsloop, nfeloop
  USE lib_mpp
  USE iom
  USE ioipsl
  USE in_out_manager
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: mpp_init
  INTEGER :: numbot = - 1
  INTEGER :: numbdy = - 1
  CONTAINS
  SUBROUTINE mpp_init
    INTEGER :: ji, jj, jn, jproc, jarea
    INTEGER :: inijmin
    INTEGER :: i2add
    INTEGER :: inum
    INTEGER :: idir, ifreq, icont
    INTEGER :: ii, il1, ili, imil
    INTEGER :: ij, il2, ilj, ijm1
    INTEGER :: iino, ijno, iiso, ijso
    INTEGER :: iiea, ijea, iiwe, ijwe
    INTEGER :: iarea0
    INTEGER :: ierr, ios
    INTEGER :: inbi, inbj, iimax, ijmax, icnt1, icnt2
    LOGICAL :: llbest, llauto
    LOGICAL :: llwrtlay
    INTEGER, ALLOCATABLE, DIMENSION(:) :: iin, ii_nono, ii_noea
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ijn, ii_noso, ii_nowe
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: iimppt, ilci, ibondi, ipproc
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: ijmppt, ilcj, ibondj, ipolj
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: ilei, ildi, iono, ioea
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: ilej, ildj, ioso, iowe
    LOGICAL, ALLOCATABLE, DIMENSION(:, :) :: llisoce
    NAMELIST /nambdy/ ln_bdy, nb_bdy, ln_coords_file, cn_coords_file, ln_mask_file, cn_mask_file, cn_dyn2d, nn_dyn2d_dta, &
&cn_dyn3d, nn_dyn3d_dta, cn_tra, nn_tra_dta, ln_tra_dmp, ln_dyn3d_dmp, rn_time_dmp, rn_time_dmp_out, cn_ice, nn_ice_dta, &
&rn_ice_tem, rn_ice_sal, rn_ice_age, ln_vol, nn_volctl, nn_rimwidth, nb_jpk_bdy
    llwrtlay = lwp .OR. ln_ctl .OR. sn_cfctl % l_layout
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, nambdy, IOSTAT = ios, ERR = 903)
903 IF (ios /= 0) CALL ctl_nam(ios, 'nambdy in reference namelist (mppini)', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, nambdy, IOSTAT = ios, ERR = 904)
904 IF (ios > 0) CALL ctl_nam(ios, 'nambdy in configuration namelist (mppini)', lwp)
    IF (ln_read_cfg) CALL iom_open(cn_domcfg, numbot)
    IF (ln_bdy .AND. ln_mask_file) CALL iom_open(cn_mask_file, numbdy)
    IF (lwp) THEN
      WRITE(numout, FMT = *) 'mpp_init:'
      WRITE(numout, FMT = *) '~~~~~~~~ '
      WRITE(numout, FMT = *)
    END IF
    IF (jpni < 1 .OR. jpnj < 1) THEN
      CALL mpp_init_bestpartition(mppsize, jpni, jpnj)
      llauto = .TRUE.
      llbest = .TRUE.
    ELSE
      llauto = .FALSE.
      CALL mpp_init_bestpartition(mppsize, inbi, inbj, icnt2)
      CALL mpp_basic_decomposition(jpni, jpnj, jpimax, jpjmax)
      CALL mpp_basic_decomposition(inbi, inbj, iimax, ijmax)
      icnt1 = jpni * jpnj - mppsize
      IF (lwp) THEN
        WRITE(numout, 9000) '   The chosen domain decomposition ', jpni, ' x ', jpnj, ' with ', icnt1, ' land subdomains'
        WRITE(numout, 9002) '      - uses a total of ', mppsize, ' mpi process'
        WRITE(numout, 9000) '      - has mpi subdomains with a maximum size of (jpi = ', jpimax, ', jpj = ', jpjmax, ', jpi*jpj = &
&', jpimax * jpjmax, ')'
        WRITE(numout, 9000) '   The best domain decompostion ', inbi, ' x ', inbj, ' with ', icnt2, ' land subdomains'
        WRITE(numout, 9002) '      - uses a total of ', inbi * inbj - icnt2, ' mpi process'
        WRITE(numout, 9000) '      - has mpi subdomains with a maximum size of (jpi = ', iimax, ', jpj = ', ijmax, ', jpi*jpj = ', &
&iimax * ijmax, ')'
      END IF
      IF (iimax * ijmax < jpimax * jpjmax) THEN
        llbest = .FALSE.
        IF (inbi * inbj - icnt2 < mppsize) THEN
          WRITE(ctmp1, FMT = *) '   ==> You could therefore have smaller mpi subdomains with less mpi processes'
        ELSE
          WRITE(ctmp1, FMT = *) '   ==> You could therefore have smaller mpi subdomains with the same number of mpi processes'
        END IF
        CALL ctl_warn(' ', ctmp1, ' ', '    ---   YOU ARE WASTING CPU...   ---', ' ')
      ELSE IF (iimax * ijmax == jpimax * jpjmax .AND. (inbi * inbj - icnt2) < mppsize) THEN
        llbest = .FALSE.
        WRITE(ctmp1, FMT = *) '   ==> You could therefore have the same mpi subdomains size with less mpi processes'
        CALL ctl_warn(' ', ctmp1, ' ', '    ---   YOU ARE WASTING CPU...   ---', ' ')
      ELSE
        llbest = .TRUE.
      END IF
    END IF
    ALLOCATE(llisoce(jpni, jpnj))
    CALL mpp_init_isoce(jpni, jpnj, llisoce)
    inijmin = COUNT(llisoce)
    IF (mppsize < inijmin) THEN
      WRITE(ctmp1, 9001) '   With this specified domain decomposition: jpni = ', jpni, ' jpnj = ', jpnj
      WRITE(ctmp2, 9002) '   we can eliminate only ', jpni * jpnj - inijmin, ' land mpi subdomains therefore '
      WRITE(ctmp3, 9001) '   the number of ocean mpi subdomains (', inijmin, ') exceed the number of MPI processes:', mppsize
      WRITE(ctmp4, FMT = *) '   ==>>> There is the list of best domain decompositions you should use: '
      CALL ctl_stop(ctmp1, ctmp2, ctmp3, ' ', ctmp4, ' ')
      CALL mpp_init_bestpartition(mppsize, ldlist = .TRUE.)
      CALL ctl_stop('STOP')
    END IF
    IF (mppsize > jpni * jpnj) THEN
      IF (lwp) THEN
        WRITE(numout, 9003) '   The number of mpi processes: ', mppsize
        WRITE(numout, 9003) '   exceeds the maximum number of subdomains (ocean+land) = ', jpni * jpnj
        WRITE(numout, 9001) '   defined by the following domain decomposition: jpni = ', jpni, ' jpnj = ', jpnj
        WRITE(numout, FMT = *) '   You should: '
        IF (llauto) THEN
          WRITE(numout, FMT = *) '     - either prescribe your domain decomposition with the namelist variables'
          WRITE(numout, FMT = *) '       jpni and jpnj to match the number of mpi process you want to use, '
          WRITE(numout, FMT = *) '       even IF it not the best choice...'
          WRITE(numout, FMT = *) '     - or keep the automatic and optimal domain decomposition by picking up one'
          WRITE(numout, FMT = *) '       of the number of mpi process proposed in the list bellow'
        ELSE
          WRITE(numout, FMT = *) '     - either properly prescribe your domain decomposition with jpni and jpnj'
          WRITE(numout, FMT = *) '       in order to be consistent with the number of mpi process you want to use'
          WRITE(numout, FMT = *) '       even IF it not the best choice...'
          WRITE(numout, FMT = *) '     - or use the automatic and optimal domain decomposition and pick up one of'
          WRITE(numout, FMT = *) '       the domain decomposition proposed in the list bellow'
        END IF
        WRITE(numout, FMT = *)
      END IF
      CALL mpp_init_bestpartition(mppsize, ldlist = .TRUE.)
      CALL ctl_stop('STOP')
    END IF
    jpnij = mppsize
    IF (mppsize > inijmin) THEN
      WRITE(ctmp1, 9003) '   The number of mpi processes: ', mppsize
      WRITE(ctmp2, 9003) '   exceeds the maximum number of ocean subdomains = ', inijmin
      WRITE(ctmp3, 9002) '   we suppressed ', jpni * jpnj - mppsize, ' land subdomains '
      WRITE(ctmp4, 9002) '   BUT we had to keep ', mppsize - inijmin, ' land subdomains that are useless...'
      CALL ctl_warn(ctmp1, ctmp2, ctmp3, ctmp4, ' ', '    --- YOU ARE WASTING CPU... ---', ' ')
    ELSE
      IF (lwp) THEN
        IF (llbest) WRITE(numout, FMT = *) '   ==> you use the best mpi decomposition'
        WRITE(numout, FMT = *)
        WRITE(numout, 9003) '   Number of mpi processes: ', mppsize
        WRITE(numout, 9003) '   Number of ocean subdomains = ', inijmin
        WRITE(numout, 9003) '   Number of suppressed land subdomains = ', jpni * jpnj - inijmin
        WRITE(numout, FMT = *)
      END IF
    END IF
9000 FORMAT(A, I4, A, I4, A, I7, A)
9001 FORMAT(A, I4, A, I4)
9002 FORMAT(A, I4, A)
9003 FORMAT(A, I5)
    IF (numbot /= - 1) CALL iom_close(numbot)
    IF (numbdy /= - 1) CALL iom_close(numbdy)
    ALLOCATE(nfiimpp(jpni, jpnj), nfipproc(jpni, jpnj), nfilcit(jpni, jpnj), nimppt(jpnij), ibonit(jpnij), nlcit(jpnij), &
&nlcjt(jpnij), njmppt(jpnij), ibonjt(jpnij), nldit(jpnij), nldjt(jpnij), nleit(jpnij), nlejt(jpnij), iin(jpnij), ii_nono(jpnij), &
&ii_noea(jpnij), ijn(jpnij), ii_noso(jpnij), ii_nowe(jpnij), iimppt(jpni, jpnj), ilci(jpni, jpnj), ibondi(jpni, jpnj), &
&ipproc(jpni, jpnj), ijmppt(jpni, jpnj), ilcj(jpni, jpnj), ibondj(jpni, jpnj), ipolj(jpni, jpnj), ilei(jpni, jpnj), ildi(jpni, &
&jpnj), iono(jpni, jpnj), ioea(jpni, jpnj), ilej(jpni, jpnj), ildj(jpni, jpnj), ioso(jpni, jpnj), iowe(jpni, jpnj), STAT = ierr)
    CALL mpp_sum('mppini', ierr)
    IF (ierr /= 0) CALL ctl_stop('STOP', 'mpp_init: unable to allocate standard ocean arrays')
    nreci = 2 * nn_hls
    nrecj = 2 * nn_hls
    CALL mpp_basic_decomposition(jpni, jpnj, jpimax, jpjmax, iimppt, ijmppt, ilci, ilcj)
    nfiimpp(:, :) = iimppt(:, :)
    nfilcit(:, :) = ilci(:, :)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'MPI Message Passing MPI - domain lay out over processors'
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) '   defines mpp subdomains'
      WRITE(numout, FMT = *) '      jpni = ', jpni
      WRITE(numout, FMT = *) '      jpnj = ', jpnj
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) '      sum ilci(i,1) = ', SUM(ilci(:, 1)), ' jpiglo = ', jpiglo
      WRITE(numout, FMT = *) '      sum ilcj(1,j) = ', SUM(ilcj(1, :)), ' jpjglo = ', jpjglo
    END IF
    l_Iperio = jpni == 1 .AND. (jperio == 1 .OR. jperio == 4 .OR. jperio == 6 .OR. jperio == 7)
    l_Jperio = jpnj == 1 .AND. (jperio == 2 .OR. jperio == 7)
    DO jarea = 1, jpni * jpnj
      iarea0 = jarea - 1
      ii = 1 + MOD(iarea0, jpni)
      ij = 1 + iarea0 / jpni
      ili = ilci(ii, ij)
      ilj = ilcj(ii, ij)
      ibondi(ii, ij) = 0
      IF (ii == 1) ibondi(ii, ij) = - 1
      IF (ii == jpni) ibondi(ii, ij) = 1
      IF (jpni == 1) ibondi(ii, ij) = 2
      ibondj(ii, ij) = 0
      IF (ij == 1) ibondj(ii, ij) = - 1
      IF (ij == jpnj) ibondj(ii, ij) = 1
      IF (jpnj == 1) ibondj(ii, ij) = 2
      ioso(ii, ij) = iarea0 - jpni
      iowe(ii, ij) = iarea0 - 1
      ioea(ii, ij) = iarea0 + 1
      iono(ii, ij) = iarea0 + jpni
      ildi(ii, ij) = 1 + nn_hls
      ilei(ii, ij) = ili - nn_hls
      ildj(ii, ij) = 1 + nn_hls
      ilej(ii, ij) = ilj - nn_hls
      IF (jperio == 1 .OR. jperio == 4 .OR. jperio == 6 .OR. jperio == 7) THEN
        IF (jpni /= 1) ibondi(ii, ij) = 0
        IF (ii == 1) iowe(ii, ij) = iarea0 + (jpni - 1)
        IF (ii == jpni) ioea(ii, ij) = iarea0 - (jpni - 1)
      END IF
      IF (jperio == 2 .OR. jperio == 7) THEN
        IF (jpnj /= 1) ibondj(ii, ij) = 0
        IF (ij == 1) ioso(ii, ij) = iarea0 + jpni * (jpnj - 1)
        IF (ij == jpnj) iono(ii, ij) = iarea0 - jpni * (jpnj - 1)
      END IF
      ipolj(ii, ij) = 0
      IF (jperio == 3 .OR. jperio == 4) THEN
        ijm1 = jpni * (jpnj - 1)
        imil = ijm1 + (jpni + 1) / 2
        IF (jarea > ijm1) ipolj(ii, ij) = 3
        IF (MOD(jpni, 2) == 1 .AND. jarea == imil) ipolj(ii, ij) = 4
        IF (ipolj(ii, ij) == 3) iono(ii, ij) = jpni * jpnj - jarea + ijm1
      END IF
      IF (jperio == 5 .OR. jperio == 6) THEN
        ijm1 = jpni * (jpnj - 1)
        imil = ijm1 + (jpni + 1) / 2
        IF (jarea > ijm1) ipolj(ii, ij) = 5
        IF (MOD(jpni, 2) == 1 .AND. jarea == imil) ipolj(ii, ij) = 6
        IF (ipolj(ii, ij) == 5) iono(ii, ij) = jpni * jpnj - jarea + ijm1
      END IF
    END DO
    ipproc(:, :) = - 1
    icont = - 1
    DO jarea = 1, jpni * jpnj
      iarea0 = jarea - 1
      ii = 1 + MOD(iarea0, jpni)
      ij = 1 + iarea0 / jpni
      IF (llisoce(ii, ij)) THEN
        icont = icont + 1
        ipproc(ii, ij) = icont
        iin(icont + 1) = ii
        ijn(icont + 1) = ij
      END IF
    END DO
    i2add = jpnij - inijmin
    DO jarea = 1, jpni * jpnj
      iarea0 = jarea - 1
      ii = 1 + MOD(iarea0, jpni)
      ij = 1 + iarea0 / jpni
      IF (.NOT. llisoce(ii, ij) .AND. i2add > 0) THEN
        icont = icont + 1
        ipproc(ii, ij) = icont
        iin(icont + 1) = ii
        ijn(icont + 1) = ij
        i2add = i2add - 1
      END IF
    END DO
    nfipproc(:, :) = ipproc(:, :)
    DO jarea = 1, jpni * jpnj
      ii = 1 + MOD(jarea - 1, jpni)
      ij = 1 + (jarea - 1) / jpni
      IF (ipproc(ii, ij) == - 1 .AND. 0 <= iono(ii, ij) .AND. iono(ii, ij) <= jpni * jpnj - 1) THEN
        iino = 1 + MOD(iono(ii, ij), jpni)
        ijno = 1 + iono(ii, ij) / jpni
        idir = 1
        IF (ij == jpnj .AND. ijno == jpnj) idir = - 1
        IF (ibondj(iino, ijno) == idir) ibondj(iino, ijno) = 2
        IF (ibondj(iino, ijno) == 0) ibondj(iino, ijno) = - idir
      END IF
      IF (ipproc(ii, ij) == - 1 .AND. 0 <= ioso(ii, ij) .AND. ioso(ii, ij) <= jpni * jpnj - 1) THEN
        iiso = 1 + MOD(ioso(ii, ij), jpni)
        ijso = 1 + ioso(ii, ij) / jpni
        IF (ibondj(iiso, ijso) == - 1) ibondj(iiso, ijso) = 2
        IF (ibondj(iiso, ijso) == 0) ibondj(iiso, ijso) = 1
      END IF
      IF (ipproc(ii, ij) == - 1 .AND. 0 <= ioea(ii, ij) .AND. ioea(ii, ij) <= jpni * jpnj - 1) THEN
        iiea = 1 + MOD(ioea(ii, ij), jpni)
        ijea = 1 + ioea(ii, ij) / jpni
        IF (ibondi(iiea, ijea) == 1) ibondi(iiea, ijea) = 2
        IF (ibondi(iiea, ijea) == 0) ibondi(iiea, ijea) = - 1
      END IF
      IF (ipproc(ii, ij) == - 1 .AND. 0 <= iowe(ii, ij) .AND. iowe(ii, ij) <= jpni * jpnj - 1) THEN
        iiwe = 1 + MOD(iowe(ii, ij), jpni)
        ijwe = 1 + iowe(ii, ij) / jpni
        IF (ibondi(iiwe, ijwe) == - 1) ibondi(iiwe, ijwe) = 2
        IF (ibondi(iiwe, ijwe) == 0) ibondi(iiwe, ijwe) = 1
      END IF
    END DO
    DO jproc = 1, jpnij
      ii = iin(jproc)
      ij = ijn(jproc)
      IF (ibondi(ii, ij) == - 1 .OR. ibondi(ii, ij) == 2) ildi(ii, ij) = 1
      IF (ibondi(ii, ij) == 1 .OR. ibondi(ii, ij) == 2) ilei(ii, ij) = ilci(ii, ij)
      IF (ibondj(ii, ij) == - 1 .OR. ibondj(ii, ij) == 2) ildj(ii, ij) = 1
      IF (ibondj(ii, ij) == 1 .OR. ibondj(ii, ij) == 2) ilej(ii, ij) = ilcj(ii, ij)
    END DO
    IF (lwp) THEN
      ifreq = 4
      il1 = 1
      DO jn = 1, (jpni - 1) / ifreq + 1
        il2 = MIN(jpni, il1 + ifreq - 1)
        WRITE(numout, FMT = *)
        WRITE(numout, 9400) ('***', ji = il1, il2 - 1)
        DO jj = jpnj, 1, - 1
          WRITE(numout, 9403) ('   ', ji = il1, il2 - 1)
          WRITE(numout, 9402) jj, (ilci(ji, jj), ilcj(ji, jj), ji = il1, il2)
          WRITE(numout, 9404) (ipproc(ji, jj), ji = il1, il2)
          WRITE(numout, 9403) ('   ', ji = il1, il2 - 1)
          WRITE(numout, 9400) ('***', ji = il1, il2 - 1)
        END DO
        WRITE(numout, 9401) (ji, ji = il1, il2)
        il1 = il1 + ifreq
      END DO
9400  FORMAT('           ***', 20('*************', A3))
9403  FORMAT('           *     ', 20('         *   ', A3))
9401  FORMAT('              ', 20('   ', I3, '          '))
9402  FORMAT('       ', I3, ' *  ', 20(I3, '  x', I3, '   *   '))
9404  FORMAT('           *  ', 20('      ', I3, '   *   '))
    END IF
    ii_noso(:) = - 1
    ii_nono(:) = - 1
    ii_noea(:) = - 1
    ii_nowe(:) = - 1
    DO jproc = 1, jpnij
      ii = iin(jproc)
      ij = ijn(jproc)
      IF (0 <= ioso(ii, ij) .AND. ioso(ii, ij) <= (jpni * jpnj - 1)) THEN
        iiso = 1 + MOD(ioso(ii, ij), jpni)
        ijso = 1 + ioso(ii, ij) / jpni
        ii_noso(jproc) = ipproc(iiso, ijso)
      END IF
      IF (0 <= iowe(ii, ij) .AND. iowe(ii, ij) <= (jpni * jpnj - 1)) THEN
        iiwe = 1 + MOD(iowe(ii, ij), jpni)
        ijwe = 1 + iowe(ii, ij) / jpni
        ii_nowe(jproc) = ipproc(iiwe, ijwe)
      END IF
      IF (0 <= ioea(ii, ij) .AND. ioea(ii, ij) <= (jpni * jpnj - 1)) THEN
        iiea = 1 + MOD(ioea(ii, ij), jpni)
        ijea = 1 + ioea(ii, ij) / jpni
        ii_noea(jproc) = ipproc(iiea, ijea)
      END IF
      IF (0 <= iono(ii, ij) .AND. iono(ii, ij) <= (jpni * jpnj - 1)) THEN
        iino = 1 + MOD(iono(ii, ij), jpni)
        ijno = 1 + iono(ii, ij) / jpni
        ii_nono(jproc) = ipproc(iino, ijno)
      END IF
    END DO
    ii = iin(narea)
    ij = ijn(narea)
    noso = ii_noso(narea)
    nowe = ii_nowe(narea)
    noea = ii_noea(narea)
    nono = ii_nono(narea)
    nlci = ilci(ii, ij)
    nldi = ildi(ii, ij)
    nlei = ilei(ii, ij)
    nlcj = ilcj(ii, ij)
    nldj = ildj(ii, ij)
    nlej = ilej(ii, ij)
    nbondi = ibondi(ii, ij)
    nbondj = ibondj(ii, ij)
    nimpp = iimppt(ii, ij)
    njmpp = ijmppt(ii, ij)
    jpi = nlci
    jpj = nlcj
    jpk = jpkglo
    jpim1 = jpi - 1
    jpjm1 = jpj - 1
    jpkm1 = MAX(1, jpk - 1)
    jpij = jpi * jpj
    DO jproc = 1, jpnij
      ii = iin(jproc)
      ij = ijn(jproc)
      nlcit(jproc) = ilci(ii, ij)
      nldit(jproc) = ildi(ii, ij)
      nleit(jproc) = ilei(ii, ij)
      nlcjt(jproc) = ilcj(ii, ij)
      nldjt(jproc) = ildj(ii, ij)
      nlejt(jproc) = ilej(ii, ij)
      ibonit(jproc) = ibondi(ii, ij)
      ibonjt(jproc) = ibondj(ii, ij)
      nimppt(jproc) = iimppt(ii, ij)
      njmppt(jproc) = ijmppt(ii, ij)
    END DO
    IF (llwrtlay) THEN
      CALL ctl_opn(inum, 'layout.dat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, numout, .FALSE., narea)
      WRITE(inum, FMT = '(a)') '   jpnij   jpimax  jpjmax    jpk  jpiglo  jpjglo' // ' ( local:    narea     jpi     jpj )'
      WRITE(inum, FMT = '(6i8,a,3i8,a)') jpnij, jpimax, jpjmax, jpk, jpiglo, jpjglo, ' ( local: ', narea, jpi, jpj, ' )'
      WRITE(inum, FMT = '(a)') 'nproc nlci nlcj nldi nldj nlei nlej nimp njmp nono noso nowe noea nbondi nbondj '
      DO jproc = 1, jpnij
        WRITE(inum, FMT = '(13i5,2i7)') jproc - 1, nlcit(jproc), nlcjt(jproc), nldit(jproc), nldjt(jproc), nleit(jproc), &
&nlejt(jproc), nimppt(jproc), njmppt(jproc), ii_nono(jproc), ii_noso(jproc), ii_nowe(jproc), ii_noea(jproc), ibonit(jproc), &
&ibonjt(jproc)
      END DO
    END IF
    npolj = 0
    ij = ijn(narea)
    IF (jperio == 3 .OR. jperio == 4) THEN
      IF (ij == jpnj) npolj = 3
    END IF
    IF (jperio == 5 .OR. jperio == 6) THEN
      IF (ij == jpnj) npolj = 5
    END IF
    nproc = narea - 1
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) '   resulting internal parameters : '
      WRITE(numout, FMT = *) '      nproc  = ', nproc
      WRITE(numout, FMT = *) '      nowe   = ', nowe, '   noea  =  ', noea
      WRITE(numout, FMT = *) '      nono   = ', nono, '   noso  =  ', noso
      WRITE(numout, FMT = *) '      nbondi = ', nbondi
      WRITE(numout, FMT = *) '      nbondj = ', nbondj
      WRITE(numout, FMT = *) '      npolj  = ', npolj
      WRITE(numout, FMT = *) '    l_Iperio = ', l_Iperio
      WRITE(numout, FMT = *) '    l_Jperio = ', l_Jperio
      WRITE(numout, FMT = *) '      nlci   = ', nlci
      WRITE(numout, FMT = *) '      nlcj   = ', nlcj
      WRITE(numout, FMT = *) '      nimpp  = ', nimpp
      WRITE(numout, FMT = *) '      njmpp  = ', njmpp
      WRITE(numout, FMT = *) '      nreci  = ', nreci
      WRITE(numout, FMT = *) '      nrecj  = ', nrecj
      WRITE(numout, FMT = *) '      nn_hls = ', nn_hls
    END IF
    IF (jperio >= 3 .AND. jperio <= 6 .AND. jpni > 1) THEN
      CALL mpp_ini_north
      IF (lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) '   ==>>>   North fold boundary prepared for jpni >1'
      END IF
      IF (llwrtlay) THEN
        WRITE(inum, FMT = *)
        WRITE(inum, FMT = *)
        WRITE(inum, FMT = *) 'number of subdomains located along the north fold : ', ndim_rank_north
        WRITE(inum, FMT = *) 'Rank of the subdomains located along the north fold : ', ndim_rank_north
        DO jproc = 1, ndim_rank_north, 5
          WRITE(inum, FMT = *) nrank_north(jproc : MINVAL((/jproc + 4, ndim_rank_north/)))
        END DO
      END IF
    END IF
    CALL mpp_init_ioipsl
    IF (ln_nnogather) THEN
      CALL mpp_init_nfdcom
      IF (llwrtlay) THEN
        WRITE(inum, FMT = *)
        WRITE(inum, FMT = *)
        WRITE(inum, FMT = *) 'north fold exchanges with explicit point-to-point messaging :'
        WRITE(inum, FMT = *) 'nfsloop : ', nfsloop
        WRITE(inum, FMT = *) 'nfeloop : ', nfeloop
        WRITE(inum, FMT = *) 'nsndto : ', nsndto
        WRITE(inum, FMT = *) 'isendto : ', isendto
      END IF
    END IF
    IF (llwrtlay) CLOSE(UNIT = inum)
    DEALLOCATE(iin, ijn, ii_nono, ii_noea, ii_noso, ii_nowe, iimppt, ijmppt, ibondi, ibondj, ipproc, ipolj, ilci, ilcj, ilei, &
&ilej, ildi, ildj, iono, ioea, ioso, iowe, llisoce)
  END SUBROUTINE mpp_init
  SUBROUTINE mpp_basic_decomposition(knbi, knbj, kimax, kjmax, kimppt, kjmppt, klci, klcj)
    INTEGER, INTENT(IN) :: knbi, knbj
    INTEGER, INTENT(OUT) :: kimax, kjmax
    INTEGER, DIMENSION(knbi, knbj), OPTIONAL, INTENT(OUT) :: kimppt, kjmppt
    INTEGER, DIMENSION(knbi, knbj), OPTIONAL, INTENT(OUT) :: klci, klcj
    INTEGER :: ji, jj
    INTEGER :: iresti, irestj, irm, ijpjmin
    INTEGER :: ireci, irecj
    kimax = (jpiglo - 2 * nn_hls + (knbi - 1)) / knbi + 2 * nn_hls
    kjmax = (jpjglo - 2 * nn_hls + (knbj - 1)) / knbj + 2 * nn_hls
    IF (.NOT. PRESENT(kimppt)) RETURN
    ireci = 2 * nn_hls
    irecj = 2 * nn_hls
    iresti = 1 + MOD(jpiglo - ireci - 1, knbi)
    irestj = 1 + MOD(jpjglo - irecj - 1, knbj)
    klci(1 : iresti, :) = kimax
    klci(iresti + 1 : knbi, :) = kimax - 1
    IF (MINVAL(klci) < 3) THEN
      WRITE(ctmp1, FMT = *) '   mpp_basic_decomposition: minimum value of jpi must be >= 3'
      WRITE(ctmp2, FMT = *) '   We have ', MINVAL(klci)
      CALL ctl_stop('STOP', ctmp1, ctmp2)
    END IF
    IF (jperio == 3 .OR. jperio == 4 .OR. jperio == 5 .OR. jperio == 6) THEN
      IF (jperio == 3 .OR. jperio == 4) ijpjmin = 5
      IF (jperio == 5 .OR. jperio == 6) ijpjmin = 4
      irm = knbj - irestj
      klcj(:, knbj) = MAX(ijpjmin, kjmax - irm)
      irm = irm - (kjmax - klcj(1, knbj))
      irestj = knbj - 1 - irm
      klcj(:, 1 : irestj) = kjmax
      klcj(:, irestj + 1 : knbj - 1) = kjmax - 1
    ELSE
      ijpjmin = 3
      klcj(:, 1 : irestj) = kjmax
      klcj(:, irestj + 1 : knbj) = kjmax - 1
    END IF
    IF (MINVAL(klcj) < ijpjmin) THEN
      WRITE(ctmp1, FMT = *) '   mpp_basic_decomposition: minimum value of jpj must be >= ', ijpjmin
      WRITE(ctmp2, FMT = *) '   We have ', MINVAL(klcj)
      CALL ctl_stop('STOP', ctmp1, ctmp2)
    END IF
    kimppt(:, :) = 1
    kjmppt(:, :) = 1
    IF (knbi > 1) THEN
      DO jj = 1, knbj
        DO ji = 2, knbi
          kimppt(ji, jj) = kimppt(ji - 1, jj) + klci(ji - 1, jj) - ireci
        END DO
      END DO
    END IF
    IF (knbj > 1) THEN
      DO jj = 2, knbj
        DO ji = 1, knbi
          kjmppt(ji, jj) = kjmppt(ji, jj - 1) + klcj(ji, jj - 1) - irecj
        END DO
      END DO
    END IF
  END SUBROUTINE mpp_basic_decomposition
  SUBROUTINE mpp_init_bestpartition(knbij, knbi, knbj, knbcnt, ldlist)
    INTEGER, INTENT(IN) :: knbij
    INTEGER, OPTIONAL, INTENT(OUT) :: knbi, knbj
    INTEGER, OPTIONAL, INTENT(OUT) :: knbcnt
    LOGICAL, OPTIONAL, INTENT(IN) :: ldlist
    INTEGER :: ji, jj, ii, iitarget
    INTEGER :: iszitst, iszjtst
    INTEGER :: isziref, iszjref
    INTEGER :: inbij, iszij
    INTEGER :: inbimax, inbjmax, inbijmax
    INTEGER :: isz0, isz1
    INTEGER, DIMENSION(:), ALLOCATABLE :: indexok
    INTEGER, DIMENSION(:), ALLOCATABLE :: inbi0, inbj0, inbij0
    INTEGER, DIMENSION(:), ALLOCATABLE :: iszi0, iszj0, iszij0
    INTEGER, DIMENSION(:), ALLOCATABLE :: inbi1, inbj1, inbij1
    INTEGER, DIMENSION(:), ALLOCATABLE :: iszi1, iszj1, iszij1
    LOGICAL :: llist
    LOGICAL, DIMENSION(:, :), ALLOCATABLE :: llmsk2d
    LOGICAL, DIMENSION(:, :), ALLOCATABLE :: llisoce
    REAL(KIND = wp) :: zpropland
    llist = .FALSE.
    IF (PRESENT(ldlist)) llist = ldlist
    CALL mpp_init_landprop(zpropland)
    inbij = NINT(REAL(knbij, wp) / (1.0 - zpropland))
    IF (llist) THEN
      inbijmax = inbij * 2
    ELSE
      inbijmax = inbij
    END IF
    ALLOCATE(inbi0(inbijmax), inbj0(inbijmax), iszi0(inbijmax), iszj0(inbijmax))
    inbimax = 0
    inbjmax = 0
    isziref = jpiglo * jpjglo + 1
    iszjref = jpiglo * jpjglo + 1
    DO ji = 1, inbijmax
      iszitst = (jpiglo - 2 * nn_hls + (ji - 1)) / ji + 2 * nn_hls
      IF (iszitst < isziref) THEN
        isziref = iszitst
        inbimax = inbimax + 1
        inbi0(inbimax) = ji
        iszi0(inbimax) = isziref
      END IF
      iszjtst = (jpjglo - 2 * nn_hls + (ji - 1)) / ji + 2 * nn_hls
      IF (iszjtst < iszjref) THEN
        iszjref = iszjtst
        inbjmax = inbjmax + 1
        inbj0(inbjmax) = ji
        iszj0(inbjmax) = iszjref
      END IF
    END DO
    ALLOCATE(llmsk2d(inbimax, inbjmax))
    DO jj = 1, inbjmax
      DO ji = 1, inbimax
        IF (inbi0(ji) * inbj0(jj) <= inbijmax) THEN
          llmsk2d(ji, jj) = .TRUE.
        ELSE
          llmsk2d(ji, jj) = .FALSE.
        END IF
      END DO
    END DO
    isz1 = COUNT(llmsk2d)
    ALLOCATE(inbi1(isz1), inbj1(isz1), iszi1(isz1), iszj1(isz1))
    ii = 0
    DO jj = 1, inbjmax
      DO ji = 1, inbimax
        IF (llmsk2d(ji, jj) .EQV. .TRUE.) THEN
          ii = ii + 1
          inbi1(ii) = inbi0(ji)
          inbj1(ii) = inbj0(jj)
          iszi1(ii) = iszi0(ji)
          iszj1(ii) = iszj0(jj)
        END IF
      END DO
    END DO
    DEALLOCATE(inbi0, inbj0, iszi0, iszj0)
    DEALLOCATE(llmsk2d)
    ALLOCATE(inbij1(isz1), iszij1(isz1))
    inbij1(:) = inbi1(:) * inbj1(:)
    iszij1(:) = iszi1(:) * iszj1(:)
    IF (.NOT. llist .AND. numbot == - 1 .AND. numbdy == - 1) THEN
      ii = MINLOC(inbij1, mask = iszij1 == MINVAL(iszij1), dim = 1)
      knbi = inbi1(ii)
      knbj = inbj1(ii)
      IF (PRESENT(knbcnt)) knbcnt = 0
      DEALLOCATE(inbi1, inbj1, inbij1, iszi1, iszj1, iszij1)
      RETURN
    END IF
    ALLOCATE(indexok(isz1))
    isz0 = 0
    inbij = 1
    iszij = jpiglo * jpjglo + 1
    DO WHILE (inbij <= inbijmax)
      ii = MINLOC(iszij1, mask = inbij1 == inbij, dim = 1)
      IF (iszij1(ii) < iszij) THEN
        isz0 = isz0 + 1
        indexok(isz0) = ii
        iszij = iszij1(ii)
      END IF
      inbij = MINVAL(inbij1, mask = inbij1 > inbij)
    END DO
    DEALLOCATE(inbij1, iszij1)
    ALLOCATE(inbi0(isz0), inbj0(isz0), iszi0(isz0), iszj0(isz0))
    DO ji = 1, isz0
      ii = indexok(ji)
      inbi0(ji) = inbi1(ii)
      inbj0(ji) = inbj1(ii)
      iszi0(ji) = iszi1(ii)
      iszj0(ji) = iszj1(ii)
    END DO
    DEALLOCATE(indexok, inbi1, inbj1, iszi1, iszj1)
    IF (llist) THEN
      IF (lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) '                  For your information:'
        WRITE(numout, FMT = '(a,i5,a)') '  list of the best partitions around ', knbij, ' mpi processes'
        WRITE(numout, FMT = *) '  --------------------------------------', '-----', '--------------'
        WRITE(numout, FMT = *)
      END IF
      iitarget = MINLOC(inbi0(:) * inbj0(:), mask = inbi0(:) * inbj0(:) >= knbij, dim = 1)
      DO ji = MAX(1, iitarget - 10), MIN(isz0, iitarget + 10)
        ALLOCATE(llisoce(inbi0(ji), inbj0(ji)))
        CALL mpp_init_isoce(inbi0(ji), inbj0(ji), llisoce)
        inbij = COUNT(llisoce)
        DEALLOCATE(llisoce)
        IF (lwp) WRITE(numout, FMT = '(a, i5, a, i5, a, i4, a, i4, a, i9, a, i5, a, i5, a)') 'nb_cores ', inbij, ' oce + ', &
&inbi0(ji) * inbj0(ji) - inbij, ' land ( ', inbi0(ji), ' x ', inbj0(ji), ' ), nb_points ', iszi0(ji) * iszj0(ji), ' ( ', &
&iszi0(ji), ' x ', iszj0(ji), ' )'
      END DO
      DEALLOCATE(inbi0, inbj0, iszi0, iszj0)
      RETURN
    END IF
    DEALLOCATE(iszi0, iszj0)
    inbij = inbijmax + 1
    ii = isz0 + 1
    DO WHILE (inbij > knbij)
      ii = ii - 1
      ALLOCATE(llisoce(inbi0(ii), inbj0(ii)))
      CALL mpp_init_isoce(inbi0(ii), inbj0(ii), llisoce)
      inbij = COUNT(llisoce)
      DEALLOCATE(llisoce)
    END DO
    knbi = inbi0(ii)
    knbj = inbj0(ii)
    IF (PRESENT(knbcnt)) knbcnt = knbi * knbj - inbij
    DEALLOCATE(inbi0, inbj0)
  END SUBROUTINE mpp_init_bestpartition
  SUBROUTINE mpp_init_landprop(propland)
    REAL(KIND = wp), INTENT(OUT) :: propland
    INTEGER, DIMENSION(jpni * jpnj) :: kusedom_1d
    INTEGER :: inboce, iarea
    INTEGER :: iproc, idiv, ijsz
    INTEGER :: ijstr
    LOGICAL, ALLOCATABLE, DIMENSION(:, :) :: lloce
    IF (numbot == - 1 .AND. numbdy == - 1) THEN
      propland = 0.
      RETURN
    END IF
    iproc = MINVAL((/mppsize, jpjglo / 2, 100/))
    IF (iproc == 1) THEN
      idiv = mppsize
    ELSE
      idiv = (mppsize - 1) / (iproc - 1)
    END IF
    iarea = (narea - 1) / idiv
    IF (MOD(narea - 1, idiv) == 0 .AND. iarea < iproc) THEN
      ijsz = jpjglo / iproc
      IF (iarea < MOD(jpjglo, iproc)) ijsz = ijsz + 1
      ijstr = iarea * (jpjglo / iproc) + MIN(iarea, MOD(jpjglo, iproc)) + 1
      ALLOCATE(lloce(jpiglo, ijsz))
      CALL mpp_init_readbot_strip(ijstr, ijsz, lloce)
      inboce = COUNT(lloce)
      DEALLOCATE(lloce)
    ELSE
      inboce = 0
    END IF
    CALL mpp_sum('mppini', inboce)
    propland = REAL(jpiglo * jpjglo - inboce, wp) / REAL(jpiglo * jpjglo, wp)
  END SUBROUTINE mpp_init_landprop
  SUBROUTINE mpp_init_isoce(knbi, knbj, ldisoce)
    INTEGER, INTENT(IN) :: knbi, knbj
    LOGICAL, DIMENSION(knbi, knbj), INTENT(OUT) :: ldisoce
    INTEGER, DIMENSION(knbi, knbj) :: inboce
    INTEGER, DIMENSION(knbi * knbj) :: inboce_1d
    INTEGER :: idiv, iimax, ijmax, iarea
    INTEGER :: ji, jn
    LOGICAL, ALLOCATABLE, DIMENSION(:, :) :: lloce
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: iimppt, ilci
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: ijmppt, ilcj
    IF (numbot == - 1 .AND. numbdy == - 1) THEN
      ldisoce(:, :) = .TRUE.
      RETURN
    END IF
    IF (knbj == 1) THEN
      idiv = mppsize
    ELSE IF (mppsize < knbj) THEN
      idiv = 1
    ELSE
      idiv = (mppsize - 1) / (knbj - 1)
    END IF
    inboce(:, :) = 0
    DO jn = 0, (knbj - 1) / mppsize
      iarea = (narea - 1) / idiv + jn * mppsize
      IF (MOD(narea - 1, idiv) == 0 .AND. iarea < knbj) THEN
        ALLOCATE(iimppt(knbi, knbj), ijmppt(knbi, knbj), ilci(knbi, knbj), ilcj(knbi, knbj))
        CALL mpp_basic_decomposition(knbi, knbj, iimax, ijmax, iimppt, ijmppt, ilci, ilcj)
        ALLOCATE(lloce(jpiglo, ilcj(1, iarea + 1)))
        CALL mpp_init_readbot_strip(ijmppt(1, iarea + 1), ilcj(1, iarea + 1), lloce)
        DO ji = 1, knbi
          inboce(ji, iarea + 1) = COUNT(lloce(iimppt(ji, 1) : iimppt(ji, 1) + ilci(ji, 1) - 1, :))
        END DO
        DEALLOCATE(lloce)
        DEALLOCATE(iimppt, ijmppt, ilci, ilcj)
      END IF
    END DO
    inboce_1d = RESHAPE(inboce, (/knbi * knbj/))
    CALL mpp_sum('mppini', inboce_1d)
    inboce = RESHAPE(inboce_1d, (/knbi, knbj/))
    ldisoce(:, :) = inboce(:, :) /= 0
  END SUBROUTINE mpp_init_isoce
  SUBROUTINE mpp_init_readbot_strip(kjstr, kjcnt, ldoce)
    INTEGER, INTENT(IN) :: kjstr
    INTEGER, INTENT(IN) :: kjcnt
    LOGICAL, DIMENSION(jpiglo, kjcnt), INTENT(OUT) :: ldoce
    INTEGER :: inumsave
    REAL(KIND = wp), DIMENSION(jpiglo, kjcnt) :: zbot, zbdy
    inumsave = numout
    numout = numnul
    IF (numbot /= - 1) THEN
      CALL iom_get(numbot, jpdom_unknown, 'bottom_level', zbot, kstart = (/1, kjstr/), kcount = (/jpiglo, kjcnt/))
    ELSE
      zbot(:, :) = 1.
    END IF
    IF (numbdy /= - 1) THEN
      CALL iom_get(numbdy, jpdom_unknown, 'bdy_msk', zbdy, kstart = (/1, kjstr/), kcount = (/jpiglo, kjcnt/))
      zbot(:, :) = zbot(:, :) * zbdy(:, :)
    END IF
    ldoce(:, :) = zbot(:, :) > 0.
    numout = inumsave
  END SUBROUTINE mpp_init_readbot_strip
  SUBROUTINE mpp_init_ioipsl
    INTEGER, DIMENSION(2) :: iglo, iloc, iabsf, iabsl, ihals, ihale, idid
    iglo(1) = jpiglo
    iglo(2) = jpjglo
    iloc(1) = nlci
    iloc(2) = nlcj
    iabsf(1) = nimppt(narea)
    iabsf(2) = njmppt(narea)
    iabsl(:) = iabsf(:) + iloc(:) - 1
    ihals(1) = nldi - 1
    ihals(2) = nldj - 1
    ihale(1) = nlci - nlei
    ihale(2) = nlcj - nlej
    idid(1) = 1
    idid(2) = 2
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'mpp_init_ioipsl :   iloc  = ', iloc(1), iloc(2)
      WRITE(numout, FMT = *) '~~~~~~~~~~~~~~~     iabsf = ', iabsf(1), iabsf(2)
      WRITE(numout, FMT = *) '                    ihals = ', ihals(1), ihals(2)
      WRITE(numout, FMT = *) '                    ihale = ', ihale(1), ihale(2)
    END IF
    CALL flio_dom_set(jpnij, nproc, idid, iglo, iloc, iabsf, iabsl, ihals, ihale, 'BOX', nidom)
  END SUBROUTINE mpp_init_ioipsl
  SUBROUTINE mpp_init_nfdcom
    INTEGER :: sxM, dxM, sxT, dxT, jn
    INTEGER :: njmppmax
    njmppmax = MAXVAL(njmppt)
    isendto(:) = 0
    nsndto = 0
    IF (njmpp == njmppmax) THEN
      sxM = jpiglo - nimppt(narea) - nlcit(narea) + 1
      dxM = jpiglo - nimppt(narea) + 2
      DO jn = 1, jpni
        sxT = nfiimpp(jn, jpnj)
        dxT = nfiimpp(jn, jpnj) + nfilcit(jn, jpnj) - 1
        IF (sxT < sxM .AND. sxM < dxT) THEN
          nsndto = nsndto + 1
          isendto(nsndto) = jn
        ELSE IF (sxM <= sxT .AND. dxM >= dxT) THEN
          nsndto = nsndto + 1
          isendto(nsndto) = jn
        ELSE IF (dxM < dxT .AND. sxT < dxM) THEN
          nsndto = nsndto + 1
          isendto(nsndto) = jn
        END IF
      END DO
      nfsloop = 1
      nfeloop = nlci
      DO jn = 2, jpni - 1
        IF (nfipproc(jn, jpnj) == (narea - 1)) THEN
          IF (nfipproc(jn - 1, jpnj) == - 1) nfsloop = nldi
          IF (nfipproc(jn + 1, jpnj) == - 1) nfeloop = nlei
        END IF
      END DO
    END IF
    l_north_nogather = .TRUE.
  END SUBROUTINE mpp_init_nfdcom
END MODULE mppini