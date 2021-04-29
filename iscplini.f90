MODULE iscplini
  USE oce
  USE dom_oce
  USE lib_mpp
  USE lib_fortran
  USE in_out_manager
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: iscpl_init
  PUBLIC :: iscpl_alloc
  LOGICAL, PUBLIC :: ln_hsb
  INTEGER, PUBLIC :: nn_fiscpl
  INTEGER, PUBLIC :: nn_drown
  INTEGER, PUBLIC :: nstp_iscpl
  REAL(KIND = wp), PUBLIC :: rdt_iscpl
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: hdiv_iscpl
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :, :) :: htsc_iscpl
  CONTAINS
  INTEGER FUNCTION iscpl_alloc()
    ALLOCATE(htsc_iscpl(jpi, jpj, jpk, jpts), hdiv_iscpl(jpi, jpj, jpk), STAT = iscpl_alloc)
    CALL mpp_sum('iscplini', iscpl_alloc)
    IF (iscpl_alloc > 0) CALL ctl_warn('iscpl_alloc: allocation of arrays failed')
  END FUNCTION iscpl_alloc
  SUBROUTINE iscpl_init
    INTEGER :: ios
    NAMELIST /namsbc_iscpl/ nn_fiscpl, ln_hsb, nn_drown
    nn_fiscpl = 0
    ln_hsb = .FALSE.
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namsbc_iscpl, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namsbc_iscpl in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namsbc_iscpl, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namsbc_iscpl in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namsbc_iscpl)
    nstp_iscpl = MIN(nn_fiscpl, nitend - nit000 + 1)
    rdt_iscpl = nstp_iscpl * rn_rdt
    IF (lwp) THEN
      WRITE(numout, *) 'iscpl_rst:'
      WRITE(numout, *) '~~~~~~~~~'
      WRITE(numout, *) ' coupling     flag (ln_iscpl )            = ', ln_iscpl
      WRITE(numout, *) ' conservation flag (ln_hsb   )            = ', ln_hsb
      WRITE(numout, *) ' nb of stp for cons (rn_fiscpl)           = ', nstp_iscpl
      IF (nstp_iscpl .NE. nn_fiscpl) WRITE(numout, *) 'W A R N I N G: nb of stp for cons has been modified                                            (larger than run length)'
      WRITE(numout, *) ' coupling time step                       = ', rdt_iscpl
      WRITE(numout, *) ' number of call of the extrapolation loop = ', nn_drown
    END IF
  END SUBROUTINE iscpl_init
END MODULE iscplini