MODULE asmpar
  IMPLICIT NONE
  PRIVATE
  CHARACTER(LEN = 40), PUBLIC, PARAMETER :: c_asmbkg = 'assim_background_state_Jb'
  CHARACTER(LEN = 40), PUBLIC, PARAMETER :: c_asmdin = 'assim_background_state_DI'
  CHARACTER(LEN = 40), PUBLIC, PARAMETER :: c_asmtrj = 'assim_trj'
  CHARACTER(LEN = 40), PUBLIC, PARAMETER :: c_asminc = 'assim_background_increments'
  INTEGER, PUBLIC :: nitbkg_r
  INTEGER, PUBLIC :: nitdin_r
  INTEGER, PUBLIC :: nitiaustr_r
  INTEGER, PUBLIC :: nitiaufin_r
  INTEGER, PUBLIC :: nittrjfrq
END MODULE asmpar