MODULE bdy_oce
  USE par_oce
  IMPLICIT NONE
  PUBLIC
  INTEGER, PUBLIC, PARAMETER :: jp_bdy = 10
  INTEGER, PUBLIC, PARAMETER :: jpbgrd = 3
  TYPE, PUBLIC :: OBC_INDEX
    INTEGER, DIMENSION(jpbgrd) :: nblen
    INTEGER, DIMENSION(jpbgrd) :: nblenrim
    INTEGER, POINTER, DIMENSION(:, :) :: nbi
    INTEGER, POINTER, DIMENSION(:, :) :: nbj
    INTEGER, POINTER, DIMENSION(:, :) :: nbr
    INTEGER, POINTER, DIMENSION(:, :) :: nbmap
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: nbw
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: nbd
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: nbdout
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: flagu
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: flagv
  END TYPE OBC_INDEX
  TYPE, PUBLIC :: OBC_DATA
    INTEGER, DIMENSION(2) :: nread
    LOGICAL :: ll_ssh
    LOGICAL :: ll_u2d
    LOGICAL :: ll_v2d
    LOGICAL :: ll_u3d
    LOGICAL :: ll_v3d
    LOGICAL :: ll_tem
    LOGICAL :: ll_sal
    LOGICAL :: ll_fvl
    REAL(KIND = wp), POINTER, DIMENSION(:) :: ssh
    REAL(KIND = wp), POINTER, DIMENSION(:) :: u2d
    REAL(KIND = wp), POINTER, DIMENSION(:) :: v2d
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: u3d
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: v3d
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: tem
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: sal
  END TYPE OBC_DATA
  LOGICAL, PUBLIC :: ln_bdy
  CHARACTER(LEN = 80), DIMENSION(jp_bdy) :: cn_coords_file
  CHARACTER(LEN = 80) :: cn_mask_file
  LOGICAL, DIMENSION(jp_bdy) :: ln_coords_file
  LOGICAL :: ln_mask_file
  LOGICAL :: ln_vol
  INTEGER :: nb_bdy
  INTEGER, DIMENSION(jp_bdy) :: nb_jpk_bdy
  INTEGER, DIMENSION(jp_bdy) :: nn_rimwidth
  INTEGER :: nn_volctl
  CHARACTER(LEN = 20), DIMENSION(jp_bdy) :: cn_dyn2d
  INTEGER, DIMENSION(jp_bdy) :: nn_dyn2d_dta
  CHARACTER(LEN = 20), DIMENSION(jp_bdy) :: cn_dyn3d
  INTEGER, DIMENSION(jp_bdy) :: nn_dyn3d_dta
  CHARACTER(LEN = 20), DIMENSION(jp_bdy) :: cn_tra
  INTEGER, DIMENSION(jp_bdy) :: nn_tra_dta
  LOGICAL, DIMENSION(jp_bdy) :: ln_tra_dmp
  LOGICAL, DIMENSION(jp_bdy) :: ln_dyn3d_dmp
  REAL(KIND = wp), DIMENSION(jp_bdy) :: rn_time_dmp
  REAL(KIND = wp), DIMENSION(jp_bdy) :: rn_time_dmp_out
  CHARACTER(LEN = 20), DIMENSION(jp_bdy) :: cn_ice
  INTEGER, DIMENSION(jp_bdy) :: nn_ice_dta
  REAL(KIND = wp), DIMENSION(jp_bdy) :: rn_ice_tem
  REAL(KIND = wp), DIMENSION(jp_bdy) :: rn_ice_sal
  REAL(KIND = wp), DIMENSION(jp_bdy) :: rn_ice_age
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :), TARGET :: bdytmask
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :), TARGET :: bdyumask
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :), TARGET :: bdyvmask
  REAL(KIND = wp) :: bdysurftot
  INTEGER, DIMENSION(jp_bdy) :: nn_dta
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :), TARGET :: dta_global
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :), TARGET :: dta_global_z
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :), TARGET :: dta_global_dz
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :), TARGET :: dta_global2
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :), TARGET :: dta_global2_z
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :), TARGET :: dta_global2_dz
  TYPE(OBC_INDEX), DIMENSION(jp_bdy), TARGET :: idx_bdy
  TYPE(OBC_DATA), DIMENSION(jp_bdy), TARGET :: dta_bdy
  CONTAINS
  FUNCTION bdy_oce_alloc()
    USE lib_mpp, ONLY: ctl_stop, mpp_sum
    INTEGER :: bdy_oce_alloc
    ALLOCATE(bdytmask(jpi, jpj), bdyumask(jpi, jpj), bdyvmask(jpi, jpj), STAT = bdy_oce_alloc)
    !$ACC KERNELS
    bdytmask(:, :) = 1._wp
    bdyumask(:, :) = 1._wp
    bdyvmask(:, :) = 1._wp
    !$ACC END KERNELS
    CALL mpp_sum('bdy_oce', bdy_oce_alloc)
    IF (bdy_oce_alloc /= 0) CALL ctl_stop('STOP', 'bdy_oce_alloc: failed to allocate arrays.')
  END FUNCTION bdy_oce_alloc
END MODULE bdy_oce