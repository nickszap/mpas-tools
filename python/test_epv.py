import netCDF4
import numpy as np

import vars_column

def calc_tendency_thetam2theta(theta, qv, dthetam_dt, dqv_dt):
  Rd = 287.04;
  Rv = 461.6;

  Rv_Rd = Rv/Rd
  dtheta_dt = (dthetam_dt-theta*Rv_Rd*dqv_dt)/(1.0+Rv_Rd*qv)
  return dtheta_dt

def form_theta_tendencies(data, t0, hCell, level, nLevels):
  #dtheta/dt from lw,sw,pbl,cu,microphysics, advection, net, and numerics

  theta = data.variables['theta'][t0,hCell,level]
  qv = data.variables['theta'][t0,hCell,level]

  #radiation ------
  #lw
  dtheta_dt_lw = data.variables['rthratenlw'][t0,hCell,level]
  #sw
  dtheta_dt_sw = data.variables['rthratensw'][t0,hCell,level]
  print "radiation: ", dtheta_dt_lw, dtheta_dt_sw

  #pbl -------------
  dtheta_dt_pbl = data.variables['rthblten'][t0,hCell,level]

  print "pbl: ", dtheta_dt_pbl

  #convection ------------
  dtheta_dt_cu = data.variables['rthcuten'][t0,hCell,level]

  print "cu: ", dtheta_dt_cu

  #microphysics -------------
  #no rd_diabatic_tend in outputs from registry. no rqvmicrotend in registry at all
  #dthetam_dt_micro = data.variables['rt_diabatic_tend'][t0,hCell,level] #thetam
  dthetam_dt_micro = 0.
  #dqv_dt_micro = data.variables['rqvmicrotend'][t0,hCell,level]
  dqv_dt_micro = 0.
  dtheta_dt_micro = calc_tendency_thetam2theta(theta, qv, dthetam_dt_micro, dqv_dt_micro)

  print "micro: ", dtheta_dt_micro

  #net dtheta/dt -------------
  dthetam_dt_net = data.variables['tend_theta'][t0,hCell,level] #thetam
  dqv_dt_net = data.variables['tend_qv'][t0,hCell,level]
  dtheta_dt_net = calc_tendency_thetam2theta(theta, qv, dthetam_dt_net, dqv_dt_net)

  print "net: ", dtheta_dt_net

  #advection for Dtheta/Dt
  dtheta_dx, dtheta_dy, dtheta_dz = vars_column.calc_varGradient_leastSquare_data(data, 'theta', hCell, level, nLevels, t0)

  ux = data.variables['uReconstructX'][t0, hCell, level]
  uy = data.variables['uReconstructY'][t0, hCell, level]
  w = .5*(data.variables['w'][t0,hCell,level]+data.variables['w'][t0,hCell,level+1])
  normalCell = data.variables['localVerticalUnitVectors'][hCell,:]
  ux += w*np.dot(normalCell, [1,0,0])
  uy += w*np.dot(normalCell, [0,1,0])
  uz = w*np.dot(normalCell,[0,0,1])
  
  dtheta_dt_adv = ux*dtheta_dx+uy*dtheta_dy+uz*dtheta_dz
  
  print "advection: ", dtheta_dt_adv

  #numerics, net w/ advection - sum(physics)
  dtheta_dt_num = dtheta_dt_net+dtheta_dt_adv- \
                  (dtheta_dt_lw+dtheta_dt_sw+dtheta_dt_pbl+dtheta_dt_cu+dtheta_dt_micro)

  print "numerical: ", dtheta_dt_num
  
def test_dtheta_dt():

  f = '/arctic1/nick/cases/vduda/x4.t.output.2006-09-19_00.00.00.nc'
  data = netCDF4.Dataset(f,'r')
  nLevels = len(data.dimensions['nVertLevels'])

  iCell = 20
  iLevel = 18
  timeInd = 7

  #cell location
  latCell = data.variables['latCell'][iCell]
  lonCell = data.variables['lonCell'][iCell]
  print "lat,lon: ", latCell*180./np.pi, lonCell*180./np.pi

  form_theta_tendencies(data, timeInd, iCell, iLevel, nLevels)

if __name__ == '__main__':
  test_dtheta_dt()