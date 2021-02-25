MODULE obs_profiles
  USE obs_profiles_def
  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC :: nprofsets, nprofvars, nprofextr, profdata, prodatqc
  PUBLIC :: nvelosets, nvelovars, nveloextr, velodata, veldatqc
  INTEGER :: nprofsets
  INTEGER :: nprofvars
  INTEGER :: nprofextr
  TYPE(obs_prof), POINTER :: profdata(:)
  TYPE(obs_prof), POINTER :: prodatqc(:)
  INTEGER :: nvelosets
  INTEGER :: nvelovars
  INTEGER :: nveloextr
  TYPE(obs_prof), POINTER :: velodata(:)
  TYPE(obs_prof), POINTER :: veldatqc(:)
END MODULE obs_profiles