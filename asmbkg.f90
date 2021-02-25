MODULE asmbkg
  USE oce
  USE sbc_oce
  USE zdf_oce
  USE zdfddm
  USE ldftra
  USE ldfslp
  USE tradmp
  USE zdftke
  USE eosbn2
  USE zdfmxl
  USE dom_oce, ONLY: ndastp
  USE in_out_manager
  USE iom
  USE asmpar
  USE zdfmxl
  USE ice
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: asm_bkg_wri
  CONTAINS
  SUBROUTINE asm_bkg_wri(kt)
    INTEGER, INTENT(IN) :: kt
    CHARACTER(LEN = 50) :: cl_asmbkg
    CHARACTER(LEN = 50) :: cl_asmdin
    LOGICAL :: llok
    INTEGER :: inum
    REAL(KIND = wp) :: zdate
    IF (kt == nitbkg_r) THEN
      WRITE(cl_asmbkg, FMT = '(A,".nc")') TRIM(c_asmbkg)
      cl_asmbkg = TRIM(cl_asmbkg)
      INQUIRE(FILE = cl_asmbkg, EXIST = llok)
      IF (.NOT. llok) THEN
        IF (lwp) WRITE(numout, FMT = *) ' Setting up assimilation background file ' // TRIM(c_asmbkg)
        CALL iom_open(c_asmbkg, inum, ldwrt = .TRUE.)
        IF (nitbkg_r == nit000 - 1) THEN
          zdate = REAL(ndastp)
          IF (ln_zdftke) THEN
            IF (lwp) WRITE(numout, FMT = *) ' Reading TKE (en) from restart...'
            CALL tke_rst(nit000, 'READ')
          END IF
        ELSE
          zdate = REAL(ndastp)
        END IF
        CALL iom_rstput(kt, nitbkg_r, inum, 'rdastp', zdate)
        CALL iom_rstput(kt, nitbkg_r, inum, 'un', un)
        CALL iom_rstput(kt, nitbkg_r, inum, 'vn', vn)
        CALL iom_rstput(kt, nitbkg_r, inum, 'tn', tsn(:, :, :, jp_tem))
        CALL iom_rstput(kt, nitbkg_r, inum, 'sn', tsn(:, :, :, jp_sal))
        CALL iom_rstput(kt, nitbkg_r, inum, 'sshn', sshn)
        IF (ln_zdftke) CALL iom_rstput(kt, nitbkg_r, inum, 'en', en)
        CALL iom_close(inum)
      END IF
    END IF
    IF (kt == nitdin_r) THEN
      WRITE(cl_asmdin, FMT = '(A,".nc")') TRIM(c_asmdin)
      cl_asmdin = TRIM(cl_asmdin)
      INQUIRE(FILE = cl_asmdin, EXIST = llok)
      IF (.NOT. llok) THEN
        IF (lwp) WRITE(numout, FMT = *) ' Setting up assimilation background file ' // TRIM(c_asmdin)
        CALL iom_open(c_asmdin, inum, ldwrt = .TRUE.)
        IF (nitdin_r == nit000 - 1) THEN
          zdate = REAL(ndastp)
        ELSE
          zdate = REAL(ndastp)
        END IF
        CALL iom_rstput(kt, nitdin_r, inum, 'rdastp', zdate)
        CALL iom_rstput(kt, nitdin_r, inum, 'un', un)
        CALL iom_rstput(kt, nitdin_r, inum, 'vn', vn)
        CALL iom_rstput(kt, nitdin_r, inum, 'tn', tsn(:, :, :, jp_tem))
        CALL iom_rstput(kt, nitdin_r, inum, 'sn', tsn(:, :, :, jp_sal))
        CALL iom_rstput(kt, nitdin_r, inum, 'sshn', sshn)
        IF (nn_ice == 2) THEN
          IF (ALLOCATED(at_i)) THEN
            CALL iom_rstput(kt, nitdin_r, inum, 'iceconc', at_i(:, :))
          ELSE
            CALL ctl_warn('asm_bkg_wri: Ice concentration not written to background ', &
&'as ice variable at_i not allocated on this timestep')
          END IF
        END IF
        CALL iom_close(inum)
      END IF
    END IF
  END SUBROUTINE asm_bkg_wri
END MODULE asmbkg