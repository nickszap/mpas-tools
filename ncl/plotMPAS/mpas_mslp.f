C NCLFORTSTART
      subroutine mpas_mslp( nCells, mslp, thetaSurf, hSurf, pSurf)

      implicit none
      ! get pressure at h=0m

      integer nCells
      real*8  mslp(nCells), thetaSurf(nCells)
      real*8  hSurf(nCells), pSurf(nCells)
      
C NCLEND

      !  local
      
      integer iCell
      real*8 t0, p0, g, dt_dz, exp_scale, Rd, cp, fac, Rd_cp
      
      dt_dz = -6.5/1000. !K/m
      g = 9.806
      Rd = 287.04
      cp = 1004.5
      p0 = 1.e5
      Rd_cp = Rd/cp
      
      do iCell = 1, nCells
        !convert potential temperature to temperature
        fac = pSurf(iCell)/p0
        t0 = thetaSurf(iCell)*(fac**(Rd_cp))
        
        !"reference temperature" halfway between ground and 0m
        t0 = t0-.5*hSurf(iCell)*dt_dz
        
        exp_scale = -g*hSurf(iCell)/(Rd*t0)
        mslp(iCell) = pSurf(iCell)*exp(exp_scale)
      end do
      
      end
