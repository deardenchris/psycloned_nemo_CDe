MODULE iom_def
  USE par_kind
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER, PUBLIC :: jpdom_data = 1
  INTEGER, PARAMETER, PUBLIC :: jpdom_global = 2
  INTEGER, PARAMETER, PUBLIC :: jpdom_local = 3
  INTEGER, PARAMETER, PUBLIC :: jpdom_local_full = 4
  INTEGER, PARAMETER, PUBLIC :: jpdom_local_noextra = 5
  INTEGER, PARAMETER, PUBLIC :: jpdom_local_noovlap = 6
  INTEGER, PARAMETER, PUBLIC :: jpdom_unknown = 7
  INTEGER, PARAMETER, PUBLIC :: jpdom_autoglo = 8
  INTEGER, PARAMETER, PUBLIC :: jpdom_autoglo_xy = 9
  INTEGER, PARAMETER, PUBLIC :: jpdom_autodta = 10
  INTEGER, PARAMETER, PUBLIC :: jp_r8 = 200
  INTEGER, PARAMETER, PUBLIC :: jp_r4 = 201
  INTEGER, PARAMETER, PUBLIC :: jp_i4 = 202
  INTEGER, PARAMETER, PUBLIC :: jp_i2 = 203
  INTEGER, PARAMETER, PUBLIC :: jp_i1 = 204
  INTEGER, PARAMETER, PUBLIC :: jpmax_files = 100
  INTEGER, PARAMETER, PUBLIC :: jpmax_vars = 1200
  INTEGER, PARAMETER, PUBLIC :: jpmax_dims = 4
  INTEGER, PARAMETER, PUBLIC :: jpmax_digits = 5
  INTEGER, PUBLIC :: iom_open_init = 0
  LOGICAL, PUBLIC :: lwxios
  INTEGER, PUBLIC :: nxioso
  LOGICAL, PUBLIC :: lrxios
  LOGICAL, PUBLIC :: lxios_sini = .FALSE.
  LOGICAL, PUBLIC :: lxios_set = .FALSE.
  TYPE, PUBLIC :: file_descriptor
    CHARACTER(LEN = 240) :: name
    INTEGER :: nfid
    INTEGER :: nvars
    INTEGER :: iduld
    INTEGER :: lenuld
    INTEGER :: irec
    CHARACTER(LEN = 32) :: uldname
    CHARACTER(LEN = 32), DIMENSION(jpmax_vars) :: cn_var
    INTEGER, DIMENSION(jpmax_vars) :: nvid
    INTEGER, DIMENSION(jpmax_vars) :: ndims
    LOGICAL, DIMENSION(jpmax_vars) :: luld
    INTEGER, DIMENSION(jpmax_dims, jpmax_vars) :: dimsz
    REAL(KIND = wp), DIMENSION(jpmax_vars) :: scf
    REAL(KIND = wp), DIMENSION(jpmax_vars) :: ofs
    INTEGER :: nlev
  END TYPE file_descriptor
  TYPE(file_descriptor), DIMENSION(jpmax_files), PUBLIC :: iom_file
  INTEGER, PARAMETER, PUBLIC :: max_rst_fields = 95
  TYPE, PUBLIC :: RST_FIELD
    CHARACTER(LEN = 30) :: vname = "NO_NAME"
    CHARACTER(LEN = 30) :: grid = "NO_GRID"
    LOGICAL :: active = .FALSE.
  END TYPE RST_FIELD
  TYPE(RST_FIELD), PUBLIC, SAVE :: rst_wfields(max_rst_fields), rst_rfields(max_rst_fields)
END MODULE iom_def