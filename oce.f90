MODULE oce
  USE par_oce
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: oce_alloc
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: ub, un, ua
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: vb, vn, va
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: wn
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: wi
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: hdivn
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :, :) :: tsb, tsn, tsa
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :, :) :: rab_b, rab_n
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: rn2b, rn2
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: rhd
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: rhop
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: Cu_adv
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ub_b, un_b, ua_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: vb_b, vn_b, va_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sshb, sshn, ssha
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ubb_e, ub_e, un_e, ua_e
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: vbb_e, vb_e, vn_e, va_e
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sshbb_e, sshb_e, sshn_e, ssha_e
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hu_e
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hv_e
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hur_e
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hvr_e
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ub2_b, vb2_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: un_bf, vn_bf
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: spgu, spgv
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: gtsu, gtsv
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: gru, grv
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: gtui, gtvi
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: grui, grvi
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: riceload
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: fraqsr_1lev
  CONTAINS
  INTEGER FUNCTION oce_alloc()
    INTEGER :: ierr(6)
    ierr(:) = 0
    ALLOCATE(ub(jpi, jpj, jpk), un(jpi, jpj, jpk), ua(jpi, jpj, jpk), vb(jpi, jpj, jpk), vn(jpi, jpj, jpk), va(jpi, jpj, jpk), wn(jpi, jpj, jpk), hdivn(jpi, jpj, jpk), tsb(jpi, jpj, jpk, jpts), tsn(jpi, jpj, jpk, jpts), tsa(jpi, jpj, jpk, jpts), rab_b(jpi, jpj, jpk, jpts), rab_n(jpi, jpj, jpk, jpts), rn2b(jpi, jpj, jpk), rn2(jpi, jpj, jpk), rhd(jpi, jpj, jpk), rhop(jpi, jpj, jpk), STAT = ierr(1))
    ALLOCATE(sshb(jpi, jpj), sshn(jpi, jpj), ssha(jpi, jpj), ub_b(jpi, jpj), un_b(jpi, jpj), ua_b(jpi, jpj), vb_b(jpi, jpj), vn_b(jpi, jpj), va_b(jpi, jpj), spgu(jpi, jpj), spgv(jpi, jpj), gtsu(jpi, jpj, jpts), gtsv(jpi, jpj, jpts), gru(jpi, jpj), grv(jpi, jpj), gtui(jpi, jpj, jpts), gtvi(jpi, jpj, jpts), grui(jpi, jpj), grvi(jpi, jpj), riceload(jpi, jpj), STAT = ierr(2))
    ALLOCATE(fraqsr_1lev(jpi, jpj), STAT = ierr(3))
    ALLOCATE(ssha_e(jpi, jpj), sshn_e(jpi, jpj), sshb_e(jpi, jpj), sshbb_e(jpi, jpj), ua_e(jpi, jpj), un_e(jpi, jpj), ub_e(jpi, jpj), ubb_e(jpi, jpj), va_e(jpi, jpj), vn_e(jpi, jpj), vb_e(jpi, jpj), vbb_e(jpi, jpj), hu_e(jpi, jpj), hur_e(jpi, jpj), hv_e(jpi, jpj), hvr_e(jpi, jpj), STAT = ierr(4))
    ALLOCATE(ub2_b(jpi, jpj), vb2_b(jpi, jpj), un_bf(jpi, jpj), vn_bf(jpi, jpj), STAT = ierr(6))
    oce_alloc = MAXVAL(ierr)
    IF (oce_alloc /= 0) CALL ctl_stop('STOP', 'oce_alloc: failed to allocate arrays')
  END FUNCTION oce_alloc
END MODULE oce