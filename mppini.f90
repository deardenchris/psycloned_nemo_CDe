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
    jpimax = jpiglo
    jpjmax = jpjglo
    jpi = jpiglo
    jpj = jpjglo
    jpk = jpkglo
    jpim1 = jpi - 1
    jpjm1 = jpj - 1
    jpkm1 = MAX(1, jpk - 1)
    jpij = jpi * jpj
    jpni = 1
    jpnj = 1
    jpnij = jpni * jpnj
    nimpp = 1
    njmpp = 1
    nlci = jpi
    nlcj = jpj
    nldi = 1
    nldj = 1
    nlei = jpi
    nlej = jpj
    nbondi = 2
    nbondj = 2
    nidom = FLIO_DOM_NONE
    npolj = jperio
    l_Iperio = jpni == 1 .AND. (jperio == 1 .OR. jperio == 4 .OR. jperio == 6 .OR. jperio == 7)
    l_Jperio = jpnj == 1 .AND. (jperio == 2 .OR. jperio == 7)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'mpp_init : NO massively parallel processing'
      WRITE(numout, *) '~~~~~~~~ '
      WRITE(numout, *) '   l_Iperio = ', l_Iperio, '    l_Jperio = ', l_Jperio
      WRITE(numout, *) '     npolj  = ', npolj, '      njmpp  = ', njmpp
    END IF
    IF (jpni /= 1 .OR. jpnj /= 1 .OR. jpnij /= 1) CALL ctl_stop('mpp_init: equality  jpni = jpnj = jpnij = 1 is not satisfied', 'the domain is lay out for distributed memory computing!')
  END SUBROUTINE mpp_init
END MODULE mppini