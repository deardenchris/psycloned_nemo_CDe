MODULE diawri
  USE oce
  USE dom_oce
  USE phycst
  USE dianam
  USE diahth
  USE dynadv, ONLY: ln_dynadv_vec
  USE icb_oce
  USE icbdia
  USE ldftra
  USE ldfdyn
  USE sbc_oce
  USE sbc_ice
  USE sbcssr
  USE sbcwave
  USE wet_dry
  USE zdf_oce
  USE zdfdrg
  USE zdfmxl
  USE lbclnk
  USE in_out_manager
  USE diatmb
  USE dia25h
  USE iom
  USE ioipsl
  USE lib_mpp
  USE timing
  USE diurnal_bulk
  USE cool_skin
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dia_wri
  PUBLIC :: dia_wri_state
  PUBLIC :: dia_wri_alloc
  INTEGER :: nid_T, nz_T, nh_T, ndim_T, ndim_hT
  INTEGER :: nb_T, ndim_bT
  INTEGER :: nid_U, nz_U, nh_U, ndim_U, ndim_hU
  INTEGER :: nid_V, nz_V, nh_V, ndim_V, ndim_hV
  INTEGER :: nid_W, nz_W, nh_W
  INTEGER :: ndex(1)
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_hT, ndex_hU, ndex_hV
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_T, ndex_U, ndex_V
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_bT
  CONTAINS
  INTEGER FUNCTION dia_wri_alloc()
    INTEGER, DIMENSION(2) :: ierr
    ierr = 0
    ALLOCATE(ndex_hT(jpi * jpj), ndex_T(jpi * jpj * jpk), ndex_hU(jpi * jpj), ndex_U(jpi * jpj * jpk), ndex_hV(jpi * jpj), ndex_V(jpi * jpj * jpk), STAT = ierr(1))
    dia_wri_alloc = MAXVAL(ierr)
    CALL mpp_sum('diawri', dia_wri_alloc)
  END FUNCTION dia_wri_alloc
  SUBROUTINE dia_wri(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    LOGICAL :: ll_print = .FALSE.
    CHARACTER(LEN = 40) :: clhstnam, clop, clmx
    INTEGER :: inum = 11
    INTEGER :: ji, jj, jk
    INTEGER :: ierr
    INTEGER :: iimi, iima, ipk, it, itmod, ijmi, ijma
    INTEGER :: jn, ierror
    REAL(KIND = wp) :: zsto, zout, zmax, zjulian
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zw2d
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zw3d
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    RETURN
    CALL profile_psy_data0 % PreStart('dia_wri', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('dia_wri')
    IF (ninist == 1) THEN
      CALL dia_wri_state('output.init')
      ninist = 0
    END IF
    ll_print = .FALSE.
    ll_print = ll_print .AND. lwp
    clop = "x"
    zsto = rdt
    clop = "ave(" // TRIM(clop) // ")"
    zout = nwrite * rdt
    zmax = (nitend - nit000 + 1) * rdt
    iimi = 1
    iima = jpi
    ijmi = 1
    ijma = jpj
    ipk = jpk
    it = kt
    itmod = kt - nit000 + 1
    CALL profile_psy_data0 % PostEnd
    IF (kt == nit000) THEN
      CALL profile_psy_data1 % PreStart('dia_wri', 'r1', 0, 0)
      CALL ymds2ju(nyear, nmonth, nday, rdt, zjulian)
      zjulian = zjulian - adatrj
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'Date 0 used :', nit000, ' YEAR ', nyear, ' MONTH ', nmonth, ' DAY ', nday, 'Julian day : ', zjulian
      IF (lwp) WRITE(numout, *) ' indexes of zoom = ', iimi, iima, ijmi, ijma, ' limit storage in depth = ', ipk
      IF (lwp) THEN
        CALL dia_nam(clhstnam, nwrite, ' ')
        CALL ctl_opn(inum, 'date.file', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, numout, lwp, narea)
        WRITE(inum, *) clhstnam
        CLOSE(UNIT = inum)
      END IF
      CALL dia_nam(clhstnam, nwrite, 'grid_T')
      IF (lwp) WRITE(numout, *) " Name of NETCDF file ", clhstnam
      CALL histbeg(clhstnam, jpi, glamt, jpj, gphit, iimi, iima - iimi + 1, ijmi, ijma - ijmi + 1, nit000 - 1, zjulian, rdt, nh_T, nid_T, domain_id = nidom, snc4chunks = snc4set)
      CALL histvert(nid_T, "deptht", "Vertical T levels", "m", ipk, gdept_1d, nz_T, "down")
      CALL wheneq(jpi * jpj * ipk, tmask, 1, 1., ndex_T, ndim_T)
      CALL wheneq(jpi * jpj, tmask, 1, 1., ndex_hT, ndim_hT)
      CALL profile_psy_data1 % PostEnd
      IF (ln_icebergs) THEN
        CALL profile_psy_data2 % PreStart('dia_wri', 'r2', 0, 0)
        ALLOCATE(ndex_bT(jpi * jpj * nclasses), STAT = ierror)
        CALL mpp_sum('diawri', ierror)
        CALL profile_psy_data2 % PostEnd
        IF (ierror /= 0) THEN
          CALL profile_psy_data3 % PreStart('dia_wri', 'r3', 0, 0)
          CALL ctl_stop('dia_wri: failed to allocate iceberg diagnostic array')
          CALL profile_psy_data3 % PostEnd
          RETURN
        END IF
        CALL profile_psy_data4 % PreStart('dia_wri', 'r4', 0, 0)
        CALL histvert(nid_T, "class", "Iceberg class", "number", nclasses, class_num, nb_T)
        ndim_bT = 3
        DO jn = 1, nclasses
          ndex_bT((jn - 1) * jpi * jpj + 1 : jn * jpi * jpj) = ndex_hT(1 : jpi * jpj)
        END DO
        CALL profile_psy_data4 % PostEnd
      END IF
      CALL profile_psy_data5 % PreStart('dia_wri', 'r5', 0, 0)
      CALL dia_nam(clhstnam, nwrite, 'grid_U')
      IF (lwp) WRITE(numout, *) " Name of NETCDF file ", clhstnam
      CALL histbeg(clhstnam, jpi, glamu, jpj, gphiu, iimi, iima - iimi + 1, ijmi, ijma - ijmi + 1, nit000 - 1, zjulian, rdt, nh_U, nid_U, domain_id = nidom, snc4chunks = snc4set)
      CALL histvert(nid_U, "depthu", "Vertical U levels", "m", ipk, gdept_1d, nz_U, "down")
      CALL wheneq(jpi * jpj * ipk, umask, 1, 1., ndex_U, ndim_U)
      CALL wheneq(jpi * jpj, umask, 1, 1., ndex_hU, ndim_hU)
      CALL dia_nam(clhstnam, nwrite, 'grid_V')
      IF (lwp) WRITE(numout, *) " Name of NETCDF file ", clhstnam
      CALL histbeg(clhstnam, jpi, glamv, jpj, gphiv, iimi, iima - iimi + 1, ijmi, ijma - ijmi + 1, nit000 - 1, zjulian, rdt, nh_V, nid_V, domain_id = nidom, snc4chunks = snc4set)
      CALL histvert(nid_V, "depthv", "Vertical V levels", "m", ipk, gdept_1d, nz_V, "down")
      CALL wheneq(jpi * jpj * ipk, vmask, 1, 1., ndex_V, ndim_V)
      CALL wheneq(jpi * jpj, vmask, 1, 1., ndex_hV, ndim_hV)
      CALL dia_nam(clhstnam, nwrite, 'grid_W')
      IF (lwp) WRITE(numout, *) " Name of NETCDF file ", clhstnam
      CALL histbeg(clhstnam, jpi, glamt, jpj, gphit, iimi, iima - iimi + 1, ijmi, ijma - ijmi + 1, nit000 - 1, zjulian, rdt, nh_W, nid_W, domain_id = nidom, snc4chunks = snc4set)
      CALL histvert(nid_W, "depthw", "Vertical W levels", "m", ipk, gdepw_1d, nz_W, "down")
      CALL histdef(nid_T, "votemper", "Temperature", "C", jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout)
      CALL histdef(nid_T, "vosaline", "Salinity", "PSU", jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout)
      IF (.NOT. ln_linssh) THEN
        CALL histdef(nid_T, "vovvle3t", "Level thickness", "m", jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout)
        CALL histdef(nid_T, "vovvldep", "T point depth", "m", jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout)
        CALL histdef(nid_T, "vovvldef", "Squared level deformation", "%^2", jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout)
      END IF
      CALL histdef(nid_T, "sosstsst", "Sea Surface temperature", "C", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "sosaline", "Sea Surface Salinity", "PSU", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "sossheig", "Sea Surface Height", "m", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "sowaflup", "Net Upward Water Flux", "Kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "sorunoff", "River runoffs", "Kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "sosfldow", "downward salt flux", "PSU/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      IF (ln_linssh) THEN
        CALL histdef(nid_T, "sosst_cd", "Concentration/Dilution term on temperature", "KgC/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
        CALL histdef(nid_T, "sosss_cd", "Concentration/Dilution term on salinity", "KgPSU/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      END IF
      CALL histdef(nid_T, "sohefldo", "Net Downward Heat Flux", "W/m2", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "soshfldo", "Shortwave Radiation", "W/m2", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "somixhgt", "Turbocline Depth", "m", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "somxl010", "Mixed Layer Depth 0.01", "m", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "soicecov", "Ice fraction", "[0,1]", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "sowindsp", "wind speed at 10m", "m/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      IF (ln_icebergs) THEN
        CALL histdef(nid_T, "calving", "calving mass input", "kg/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
        CALL histdef(nid_T, "calving_heat", "calving heat flux", "XXXX", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
        CALL histdef(nid_T, "berg_floating_melt", "Melt rate of icebergs + bits", "kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
        CALL histdef(nid_T, "berg_stored_ice", "Accumulated ice mass by class", "kg", jpi, jpj, nh_T, nclasses, 1, nclasses, nb_T, 32, clop, zsto, zout)
        IF (ln_bergdia) THEN
          CALL histdef(nid_T, "berg_melt", "Melt rate of icebergs", "kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "berg_buoy_melt", "Buoyancy component of iceberg melt rate", "kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "berg_eros_melt", "Erosion component of iceberg melt rate", "kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "berg_conv_melt", "Convective component of iceberg melt rate", "kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "berg_virtual_area", "Virtual coverage by icebergs", "m2", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "bits_src", "Mass source of bergy bits", "kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "bits_melt", "Melt rate of bergy bits", "kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "bits_mass", "Bergy bit density field", "kg/m2", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "berg_mass", "Iceberg density field", "kg/m2", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "berg_real_calving", "Calving into iceberg class", "kg/s", jpi, jpj, nh_T, nclasses, 1, nclasses, nb_T, 32, clop, zsto, zout)
        END IF
      END IF
      IF (.NOT. ln_cpl) THEN
        CALL histdef(nid_T, "sohefldp", "Surface Heat Flux: Damping", "W/m2", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
        CALL histdef(nid_T, "sowafldp", "Surface Water Flux: Damping", "Kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
        CALL histdef(nid_T, "sosafldp", "Surface salt flux: damping", "Kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      END IF
      IF (ln_cpl .AND. nn_ice <= 1) THEN
        CALL histdef(nid_T, "sohefldp", "Surface Heat Flux: Damping", "W/m2", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
        CALL histdef(nid_T, "sowafldp", "Surface Water Flux: Damping", "Kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
        CALL histdef(nid_T, "sosafldp", "Surface salt flux: Damping", "Kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      END IF
      clmx = "l_max(only(x))"
      CALL histend(nid_T, snc4chunks = snc4set)
      CALL histdef(nid_U, "vozocrtx", "Zonal Current", "m/s", jpi, jpj, nh_U, ipk, 1, ipk, nz_U, 32, clop, zsto, zout)
      IF (ln_wave .AND. ln_sdw) THEN
        CALL histdef(nid_U, "sdzocrtx", "Stokes Drift Zonal Current", "m/s", jpi, jpj, nh_U, ipk, 1, ipk, nz_U, 32, clop, zsto, zout)
      END IF
      CALL histdef(nid_U, "sozotaux", "Wind Stress along i-axis", "N/m2", jpi, jpj, nh_U, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histend(nid_U, snc4chunks = snc4set)
      CALL histdef(nid_V, "vomecrty", "Meridional Current", "m/s", jpi, jpj, nh_V, ipk, 1, ipk, nz_V, 32, clop, zsto, zout)
      IF (ln_wave .AND. ln_sdw) THEN
        CALL histdef(nid_V, "sdmecrty", "Stokes Drift Meridional Current", "m/s", jpi, jpj, nh_V, ipk, 1, ipk, nz_V, 32, clop, zsto, zout)
      END IF
      CALL histdef(nid_V, "sometauy", "Wind Stress along j-axis", "N/m2", jpi, jpj, nh_V, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histend(nid_V, snc4chunks = snc4set)
      CALL histdef(nid_W, "vovecrtz", "Vertical Velocity", "m/s", jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout)
      CALL histdef(nid_W, "votkeavt", "Vertical Eddy Diffusivity", "m2/s", jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout)
      CALL histdef(nid_W, "votkeavm", "Vertical Eddy Viscosity", "m2/s", jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout)
      IF (ln_zdfddm) THEN
        CALL histdef(nid_W, "voddmavs", "Salt Vertical Eddy Diffusivity", "m2/s", jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout)
      END IF
      IF (ln_wave .AND. ln_sdw) THEN
        CALL histdef(nid_W, "sdvecrtz", "Stokes Drift Vertical Current", "m/s", jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout)
      END IF
      CALL histend(nid_W, snc4chunks = snc4set)
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'End of NetCDF Initialization'
      IF (ll_print) CALL FLUSH(numout)
      CALL profile_psy_data5 % PostEnd
    END IF
    CALL profile_psy_data6 % PreStart('dia_wri', 'r6', 0, 0)
    IF (lwp .AND. MOD(itmod, nwrite) == 0) THEN
      WRITE(numout, *) 'dia_wri : write model outputs in NetCDF files at ', kt, 'time-step'
      WRITE(numout, *) '~~~~~~ '
    END IF
    IF (.NOT. ln_linssh) THEN
      CALL histwrite(nid_T, "votemper", it, tsn(:, :, :, jp_tem) * e3t_n(:, :, :), ndim_T, ndex_T)
      CALL histwrite(nid_T, "vosaline", it, tsn(:, :, :, jp_sal) * e3t_n(:, :, :), ndim_T, ndex_T)
      CALL histwrite(nid_T, "sosstsst", it, tsn(:, :, 1, jp_tem) * e3t_n(:, :, 1), ndim_hT, ndex_hT)
      CALL histwrite(nid_T, "sosaline", it, tsn(:, :, 1, jp_sal) * e3t_n(:, :, 1), ndim_hT, ndex_hT)
    ELSE
      CALL histwrite(nid_T, "votemper", it, tsn(:, :, :, jp_tem), ndim_T, ndex_T)
      CALL histwrite(nid_T, "vosaline", it, tsn(:, :, :, jp_sal), ndim_T, ndex_T)
      CALL histwrite(nid_T, "sosstsst", it, tsn(:, :, 1, jp_tem), ndim_hT, ndex_hT)
      CALL histwrite(nid_T, "sosaline", it, tsn(:, :, 1, jp_sal), ndim_hT, ndex_hT)
    END IF
    IF (.NOT. ln_linssh) THEN
      zw3d(:, :, :) = ((e3t_n(:, :, :) - e3t_0(:, :, :)) / e3t_0(:, :, :) * 100 * tmask(:, :, :)) ** 2
      CALL histwrite(nid_T, "vovvle3t", it, e3t_n(:, :, :), ndim_T, ndex_T)
      CALL histwrite(nid_T, "vovvldep", it, gdept_n(:, :, :), ndim_T, ndex_T)
      CALL histwrite(nid_T, "vovvldef", it, zw3d, ndim_T, ndex_T)
    END IF
    CALL histwrite(nid_T, "sossheig", it, sshn, ndim_hT, ndex_hT)
    CALL histwrite(nid_T, "sowaflup", it, (emp - rnf), ndim_hT, ndex_hT)
    CALL histwrite(nid_T, "sorunoff", it, rnf, ndim_hT, ndex_hT)
    CALL histwrite(nid_T, "sosfldow", it, sfx, ndim_hT, ndex_hT)
    IF (ln_linssh) THEN
      zw2d(:, :) = emp(:, :) * tsn(:, :, 1, jp_tem)
      CALL histwrite(nid_T, "sosst_cd", it, zw2d, ndim_hT, ndex_hT)
      zw2d(:, :) = emp(:, :) * tsn(:, :, 1, jp_sal)
      CALL histwrite(nid_T, "sosss_cd", it, zw2d, ndim_hT, ndex_hT)
    END IF
    CALL histwrite(nid_T, "sohefldo", it, qns + qsr, ndim_hT, ndex_hT)
    CALL histwrite(nid_T, "soshfldo", it, qsr, ndim_hT, ndex_hT)
    CALL histwrite(nid_T, "somixhgt", it, hmld, ndim_hT, ndex_hT)
    CALL histwrite(nid_T, "somxl010", it, hmlp, ndim_hT, ndex_hT)
    CALL histwrite(nid_T, "soicecov", it, fr_i, ndim_hT, ndex_hT)
    CALL histwrite(nid_T, "sowindsp", it, wndm, ndim_hT, ndex_hT)
    IF (ln_icebergs) THEN
      CALL histwrite(nid_T, "calving", it, berg_grid(1) % calving, ndim_hT, ndex_hT)
      CALL histwrite(nid_T, "calving_heat", it, berg_grid(1) % calving_hflx, ndim_hT, ndex_hT)
      CALL histwrite(nid_T, "berg_floating_melt", it, berg_grid(1) % floating_melt, ndim_hT, ndex_hT)
      CALL histwrite(nid_T, "berg_stored_ice", it, berg_grid(1) % stored_ice, ndim_bT, ndex_bT)
      IF (ln_bergdia) THEN
        CALL histwrite(nid_T, "berg_melt", it, berg_melt, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "berg_buoy_melt", it, buoy_melt, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "berg_eros_melt", it, eros_melt, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "berg_conv_melt", it, conv_melt, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "berg_virtual_area", it, virtual_area, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "bits_src", it, bits_src, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "bits_melt", it, bits_melt, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "bits_mass", it, bits_mass, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "berg_mass", it, berg_mass, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "berg_real_calving", it, real_calving, ndim_bT, ndex_bT)
      END IF
    END IF
    IF (.NOT. ln_cpl) THEN
      CALL histwrite(nid_T, "sohefldp", it, qrp, ndim_hT, ndex_hT)
      CALL histwrite(nid_T, "sowafldp", it, erp, ndim_hT, ndex_hT)
      IF (ln_ssr) zw2d(:, :) = erp(:, :) * tsn(:, :, 1, jp_sal) * tmask(:, :, 1)
      CALL histwrite(nid_T, "sosafldp", it, zw2d, ndim_hT, ndex_hT)
    END IF
    IF (ln_cpl .AND. nn_ice <= 1) THEN
      CALL histwrite(nid_T, "sohefldp", it, qrp, ndim_hT, ndex_hT)
      CALL histwrite(nid_T, "sowafldp", it, erp, ndim_hT, ndex_hT)
      IF (ln_ssr) zw2d(:, :) = erp(:, :) * tsn(:, :, 1, jp_sal) * tmask(:, :, 1)
      CALL histwrite(nid_T, "sosafldp", it, zw2d, ndim_hT, ndex_hT)
    END IF
    CALL histwrite(nid_U, "vozocrtx", it, un, ndim_U, ndex_U)
    CALL histwrite(nid_U, "sozotaux", it, utau, ndim_hU, ndex_hU)
    CALL histwrite(nid_V, "vomecrty", it, vn, ndim_V, ndex_V)
    CALL histwrite(nid_V, "sometauy", it, vtau, ndim_hV, ndex_hV)
    CALL histwrite(nid_W, "vovecrtz", it, wn, ndim_T, ndex_T)
    CALL histwrite(nid_W, "votkeavt", it, avt, ndim_T, ndex_T)
    CALL histwrite(nid_W, "votkeavm", it, avm, ndim_T, ndex_T)
    IF (ln_zdfddm) THEN
      CALL histwrite(nid_W, "voddmavs", it, avs, ndim_T, ndex_T)
    END IF
    IF (ln_wave .AND. ln_sdw) THEN
      CALL histwrite(nid_U, "sdzocrtx", it, usd, ndim_U, ndex_U)
      CALL histwrite(nid_V, "sdmecrty", it, vsd, ndim_V, ndex_V)
      CALL histwrite(nid_W, "sdvecrtz", it, wsd, ndim_T, ndex_T)
    END IF
    IF (kt == nitend) THEN
      CALL histclo(nid_T)
      CALL histclo(nid_U)
      CALL histclo(nid_V)
      CALL histclo(nid_W)
    END IF
    IF (ln_timing) CALL timing_stop('dia_wri')
    CALL profile_psy_data6 % PostEnd
  END SUBROUTINE dia_wri
  SUBROUTINE dia_wri_state(cdfile_name)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdfile_name
    INTEGER :: inum
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('dia_wri_state', 'r0', 0, 0)
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) 'dia_wri_state : single instantaneous ocean state'
    IF (lwp) WRITE(numout, *) '~~~~~~~~~~~~~   and forcing fields file created '
    IF (lwp) WRITE(numout, *) '                and named :', cdfile_name, '...nc'
    CALL iom_open(TRIM(cdfile_name), inum, ldwrt = .TRUE.)
    CALL iom_rstput(0, 0, inum, 'votemper', tsn(:, :, :, jp_tem))
    CALL iom_rstput(0, 0, inum, 'vosaline', tsn(:, :, :, jp_sal))
    CALL iom_rstput(0, 0, inum, 'sossheig', sshn)
    CALL iom_rstput(0, 0, inum, 'vozocrtx', un)
    CALL iom_rstput(0, 0, inum, 'vomecrty', vn)
    CALL iom_rstput(0, 0, inum, 'vovecrtz', wn)
    IF (ALLOCATED(ahtu)) THEN
      CALL iom_rstput(0, 0, inum, 'ahtu', ahtu)
      CALL iom_rstput(0, 0, inum, 'ahtv', ahtv)
    END IF
    IF (ALLOCATED(ahmt)) THEN
      CALL iom_rstput(0, 0, inum, 'ahmt', ahmt)
      CALL iom_rstput(0, 0, inum, 'ahmf', ahmf)
    END IF
    CALL iom_rstput(0, 0, inum, 'sowaflup', emp - rnf)
    CALL iom_rstput(0, 0, inum, 'sohefldo', qsr + qns)
    CALL iom_rstput(0, 0, inum, 'soshfldo', qsr)
    CALL iom_rstput(0, 0, inum, 'soicecov', fr_i)
    CALL iom_rstput(0, 0, inum, 'sozotaux', utau)
    CALL iom_rstput(0, 0, inum, 'sometauy', vtau)
    IF (.NOT. ln_linssh) THEN
      CALL iom_rstput(0, 0, inum, 'vovvldep', gdept_n)
      CALL iom_rstput(0, 0, inum, 'vovvle3t', e3t_n)
    END IF
    IF (ln_wave .AND. ln_sdw) THEN
      CALL iom_rstput(0, 0, inum, 'sdzocrtx', usd)
      CALL iom_rstput(0, 0, inum, 'sdmecrty', vsd)
      CALL iom_rstput(0, 0, inum, 'sdvecrtz', wsd)
    END IF
    CALL iom_close(inum)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE dia_wri_state
END MODULE diawri
