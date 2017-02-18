
! ------------------------------------------------------------- SMC_Npts
      subroutine smc_npts(f_smc, npts, sps)

* Open smc file and obtain length of time series:

* Dates: 12/07/00 - Written by D. Boore
*        07/17/03 - Read all of int_head and real_head, and add sps to output
!        04/04/10 - Add trim_c, and use nc_f_smc.

      integer int_head(48)
      real real_head(50)
      character f_smc*(*)

      call get_lun(nu_smc)
      call trim_c(f_smc,nc_f_smc)
      open(unit=nu_smc,file=f_smc(1:nc_f_smc),status='unknown')
      call skip(nu_smc,11)

      read(nu_smc, '(8I10)') (int_head(i), i=1, 48)
      npts = int_head(17)

      read(nu_smc, '(5e15.7)') (real_head(i), i = 1, 50)
      sps = real_head(2)

      close(nu_smc)        

      return
      end
! ------------------------------------------------------------- SMC_Npts  
