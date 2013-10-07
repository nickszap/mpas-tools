C NCLFORTSTART
      subroutine mpas_vinterp( nCells, nLevels, interpVal,
     &                         latCell, field0, field1,field,
     &                         missingVal)

      implicit none
      ! interpolate columns of field1 to where field0 is interpVal
      ! input

      integer nCells, nLevels
      real*8  field0(nLevels,nCells), field1(nLevels,nCells)
      real*8  interpVal, missingVal
      real*8  field(nCells), latCell(nCells)
C NCLEND

      !  local
      
      integer iCell, iLev, levInd
      real*8 valh, vall, sgnh, sgnl, dv_dl, levFrac, valInterpCell

      do iCell = 1, nCells
        !starting from top, trap val if values on opposite side
        levInd = -1 !what should happen with missing values?
        levFrac = 0.0
        valInterpCell = sign(interpVal, latCell(iCell))
        do iLev = nLevels,2,-1
        !do iLev = 2,nLevels-1,1
          valh = field0(iLev,iCell)
          vall = field0(iLev-1,iCell)
          sgnh = valh-valInterpCell
          sgnl = vall-valInterpCell
          if (sgnh*sgnl<= 0.0) then
            !sandwiched value. equal in case val0 is a vals[l].
            !get linear interpolation: val0 = vals[l]+dvals/dl * dl
            !Avoid divide by 0 by just assuming value is 
            !halfway between...
     
            dv_dl = valh-vall;
            if (abs(dv_dl)<1.e-6) then
              levFrac = 0.5;
            else
              levFrac = (interpVal-vall)/dv_dl
            end if
            
            levInd = iLev-1
            
            exit
          end if
        end do ! done searching column

        !find value of field using index we just found
        if (levInd<0) then
          field(iCell) = missingVal
        else
          valh = field1(levInd+1,iCell)
          vall = field1(levInd,iCell)
        
          dv_dl = valh-vall
          field(iCell) = vall+dv_dl*levFrac
        end if
      end do
      
      return
      
      end
