# -*- coding: utf-8 -*-

# imports
import numpy as np
import output_data

def driver_tracers():
  '''
  given .nc (MPAS NETCDF) file and initial location for tracer, we want to 
  advect that tracer with the wind velocity through all time steps.
  '''
  
  #file properties
  ncfname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc' #input file
  data = output_data.open_netcdf_data(ncfname)
  
  tSteps = len(data.dimensions['Time'])
  vertLevels = len(data.dimensions['nVertLevels']);
  
  r = np.empty(3);
  r[0] = 50.; r[1]=50.; r[2] = 5000.;
  dt = 6.*60.*60.; #can parse xtime in file as actual timestep
  
  #loop through time steps
  vel = np.empty(3);
  for i in range(tSteps):
    (hCell, vCell) = findOwner_coord(r, data, vertLevels);
    vel[0]=data.variables['uReconstructX'][i,hCell,vCell]
    vel[1]=data.variables['uReconstructY'][i,hCell,vCell]
    vel[2]=data.variables['uReconstructZ'][i,hCell,vCell]
    r = integrate(r,u,dt);
    print r
    
  data.close()
  
def integrate(r, u, dt):
  '''
  given a position and velocity as numpy arrays, get new position at next time
  '''
  return(r+u*dt) #if these aren't numpy,netcdf4,..., + will append them!

def findSurfaceCell_bruteForce(r, data):
  '''
  Find the cell on a given spherical level that contains the x,y,z coordinate.
  '''
  #as a quick first try: 
  #1) we're going to assume the mesh is voronoi so the closest cell
  #  is also the cell that contains the point.
  #2) find closest by calculating 'distance' to all cell centers. a global argmin seems silly since we have a mesh!
  #   at the least, we could do the simple walk (especially since we're w/o boundaries)
  
  pass

def findOwner_coord(r, data, vLevels):
  #find the horizontal owner and then the owner in that column.
  #return (horizIndex, vertLevel)
  
  hCell = findOwner_horizNbrs(r, data, 0); #initial guess of cell 0 could be better?
  zInterface = data.variables['zgrid'][:]
  (l,dl) = output_data.calcIndexOfValue(r[2],zInterface[hCell,:], vLevels+1) #this might be buggy outside top and bottom boundaries
  return(hCell,l)

def findOwner_horizNbrs(r, data, cellId):
  '''
  find the cell that contains the xyz coordinate by walking in the direction
  of the closest neighbor (in spherical distance) to the initial guess. 
  Owner is the cell closer to point than any nbr (local min).
  '''
  
  #if we have the memory, we can avoid recomputing distances for nbrs of previous
  #nbrs by storing the whole mesh. say d=Inf until i calc it.
  
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:];
  latCell = data.variables['latCell'][:];
  lonCell = data.variables['lonCell'][:];
  
  #convert r to spherical coords
  rll = calc_xyz2Spherical(r);
  llNbr = np.empty(2);
  llNbr[0]=latCell[cellId]; llNbr[1]=lonCell[cellId];
  dCell = calc_distSphere(1., rll[1:], llNbr);
  
  flag = 1;
  while (flag==1):
    #keep going towards cell closer to point until no nbr is closer
    flag =0;
    nbrs = getNbrs_cell(cellId, nEdgesOnCell[cellId], cellsOnCell[cellId,:]);
    for i in range(nEdgesOnCell[cellId]):
      nbrId = nbrs[i];
      llNbr[0]=latCell[nbrId]; llNbr[1]=lonCell[nbrId];
      dNbr = calc_distSphere(1., rll[1:], llNbr);
      if (dNbr<dCell):
        dCell = dNbr;
        cellId = nbrId;
        flag=1;
      #running min
    #end of getting min of self and neighbors
  return(cellId);


def getNbrs_cell(cell, nEdgeOnCell, cellsOnCell):
  '''
  return the indices of all neighbors to a given cell
  '''
  #since know that have nEdgesOnCells neighbors since we don't have lateral boundaries,
  #can loop through cellsOnCell to get nbrs on level
  nNbrs = nEdgeOnCell[cell];
  #return (cellsOnCell[0:nNbrs-1]);
  return (cellsOnCell[0:nNbrs]-1); #0:6 give [0,1,...,5] I believe. We also change 1-index to 0
  pass

def calc_xyz2Spherical(r):
  '''
  return the (radius,lat,lon) of an xyz point as the mesher does
  '''
  s = np.empty(3)
  s[0] = np.sqrt(np.dot(r,r))
  s[1] = np.arcsin(r[2]/s[0]) #this is the mesher way
  #s[1] = np.arccos(r[2]/s[0]) #this is the wiki way
  
  #maybe beware divide by 0, but numpy handled some basic cases like atan2(-1,0), atan2(0,0) fine
  s[2] = np.arctan2(r[1],r[0])
  return(s)

def calc_distSphere(r, ll1, ll2):
  '''
  calculate the distance on the sphere between two latLon points on the same radius.
  could use haversine or vincenty formulas for better conditioning/more accuracy apparently.
  '''
  
  '''
  in some function of some part of mpas init...:
  arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
              cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
 sphere_distance = 2.*radius*asin(arg1)
  '''
  
  dll = ll2-ll1;
  latTerm = np.sin(.5*dll[0]); latTerm = latTerm*latTerm;
  lonTerm = np.sin(.5*dll[1]); lonTerm = lonTerm*lonTerm*cos(ll1[0])*cos(ll2[0]);
  dAngle = np.sqrt(latTerm+lonTerm)
  
  dist = 2.*r*np.arcsin(dAngle)
  
  return(dist)

def test_spherical(data):
  '''
  test whether our lat lon calcs match the lat lon of the mesh
  surprise, they don't!!! apparently lon is sometimes of by -2pi. correct that!!!!
  '''
  
  #cell locations on unit sphere
  xc = data.variables['xCell'][:]
  yc = data.variables['yCell'][:]
  zc = data.variables['zCell'][:]
  
  latc = data.variables['latCell'][:]
  lonc = data.variables['lonCell'][:]
  
  #calc rms difference for surface layer
  dims = xc.shape
  nCells = dims[0]
  rms = 0.
  xyz = np.empty(3)
  dlat = np.empty(nCells); dlon = np.empty(nCells)
  for i in range(nCells):
    xyz[0]=xc[i]; xyz[1]=yc[i]; xyz[2]=zc[i];
    rll = calc_xyz2Spherical(xyz);
    
    dlat[i] = rll[1]-latc[i]; dlon[i] = rll[2]-lonc[i];
    #it doesn't matter if we're off by multiples of 2pi
    eps = float(1e-6); #if we mod -.0001 by 2pi, it becomes 2pi-.0001. we want it to be 0
    dlat[i] = np.mod(dlat[i]+eps, 2.*np.pi); dlon[i]=np.mod(dlon[i]+eps,2.*np.pi);
    
    rms =rms + dlat[i]*dlat[i] + dlon[i]*dlon[i]
    
  rms = rms/nCells; rms = np.sqrt(rms)
  print "RMS of latLon between mesh and calc is: ",rms

  import matplotlib.pyplot as plt
  plt.plot(dlat,'b', dlon, 'g');
  plt.show()
  return (dlat,dlon)
  
'''
Routine from the module_sphere_utilities.F in mesher to check whether we're working with the same lat/lon:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE CONVERT_XL
!
! Convert (x, y, z) to a (lat, lon) location on a sphere with 
!    radius sqrt(x^2 + y^2 + z^2).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine convert_xl(x, y, z, latlon)

   use data_types

   implicit none

   real, intent(in) :: x, y, z
   type (geo_point), intent(out) :: latlon

   real :: dl, clat, pii, rtod
   real :: eps
   parameter (eps=1.e-10)

   pii = 2.*asin(1.0)
   rtod=180./pii
   dl = sqrt(x*x + y*y + z*z)

   latlon%lat = asin(z/dl)

!  check for being close to either pole

   if (abs(x) > eps) then

      if (abs(y) > eps) then

         latlon%lon = atan(abs(y/x))

         if ((x <= 0.) .and. (y >= 0.)) then
            latlon%lon = pii-latlon%lon
         else if ((x <= 0.) .and. (y < 0.)) then
            latlon%lon = latlon%lon+pii
         else if ((x >= 0.) .and. (y <= 0.)) then
            latlon%lon = 2*pii-latlon%lon
         end if

      else ! we're either on longitude 0 or 180

         if (x > 0) then
            latlon%lon = 0.
         else
            latlon%lon = pii
         end if

      end if

   else if (abs(y) > eps) then  

      if (y > 0) then 
         latlon%lon = pii/2.
      else
         latlon%lon = 3.*pii/2.
      end if

   else  ! we are at a pole

      latlon%lon = 0.

   end if

end subroutine convert_xl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION SPHERE_DISTANCE
!
! Given two (latitude, longitude) coordinates on the surface of a sphere,
!    plus the radius of the sphere, compute the great circle distance between
!    the points.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function sphere_distance(p1, p2, radius)

   use data_types

   implicit none

   type (geo_point), intent(in) :: p1, p2
   real, intent(in) :: radius

   real :: arg1

   arg1 = sqrt( sin(0.5*(p2%lat-p1%lat))**2 +  &
                cos(p1%lat)*cos(p2%lat)*sin(0.5*(p2%lon-p1%lon))**2 )
   sphere_distance = 2.*radius*asin(arg1)

end function sphere_distance
'''

if __name__=='__main__':
  driver_tracers()

