'''
Idea is that differences will align themselves along most unstable modes 
under iteration, like the power method for eigenvalues.
'''

def main():
  #compute base field ----------------------------------------
  
  #compute initial perturbation (w/ rescaling) -----------------------
  
  #evolve perturbation ------------------------------------------
  while(True):
    #compute new IC ----
    #add perturbation delta
    
    #overwrite field in file
    
    #integrate ------------
    
    #compute perturbation ------------------------------------
    pass
  
  return 0

#The following Lorenz stuff is adapted from: 
# http://matplotlib.org/examples/mplot3d/lorenz_attractor.html
# Plot of the Lorenz Attractor based on Edward Lorenz's 1963 "Deterministic
# Nonperiodic Flow" publication.
# http://journals.ametsoc.org/doi/abs/10.1175/1520-0469%281963%29020%3C0130%3ADNF%3E2.0.CO%3B2
#
# Note: Because this is a simple non-linear ODE, it would be more easily
#       done using SciPy's ode solver, but this approach depends only
#       upon NumPy.

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def lorenz(x, y, z, s=10, r=28, b=2.667) :
    x_dot = s*(y - x)
    y_dot = r*x - y - x*z
    z_dot = x*y - b*z
    return x_dot, y_dot, z_dot

def integrate_lorenz(stepCnt=10000, dt=.01, x0=0., y0=1., z0=1.05):
  
  # Need one more for the initial values
  xs = np.empty((stepCnt + 1,))
  ys = np.empty((stepCnt + 1,))
  zs = np.empty((stepCnt + 1,))

  # Setting initial values
  xs[0], ys[0], zs[0] = (x0, y0, z0)

  # Stepping through "time".
  for i in range(stepCnt) :
      # Derivatives of the X, Y, Z state
      x_dot, y_dot, z_dot = lorenz(xs[i], ys[i], zs[i])
      xs[i + 1] = xs[i] + (x_dot * dt)
      ys[i + 1] = ys[i] + (y_dot * dt)
      zs[i + 1] = zs[i] + (z_dot * dt)
  
  return (xs, ys, zs)

def run_lorenz():
  
  xs,ys,zs = integrate_lorenz()

  fig = plt.figure()
  ax = fig.gca(projection='3d')

  ax.plot(xs, ys, zs, 'b+')
  ax.set_xlabel("X Axis")
  ax.set_ylabel("Y Axis")
  ax.set_zlabel("Z Axis")
  ax.set_title("Lorenz Attractor")

  #plt.show()
  return plt
  
def breed_lorenz(nBreedIter=10, dt = .01):
  
  #nBreedIter = 10
  nIntegrateIter = 5; #dt=.0001
  
  #compute base field ----------------------------------------
  x0,y0,z0 = integrate_lorenz(0, dt)
  
  #compute initial perturbation (w/ rescaling) -----------------------
  #we'll evolve with different timesteps to get a rough direction of growing error
  perturbMag = dt/100.
  perturbVec = np.empty((nBreedIter,3))
  xb,yb,zb = integrate_lorenz(nIntegrateIter, dt)
  xs,ys,zs = integrate_lorenz(nIntegrateIter*2, dt/2)
  
  perturbVec[0,0] = xb[-1]-xs[-1] #assume the smaller step is more realistic
  perturbVec[0,1] = yb[-1]-ys[-1]
  perturbVec[0,2] = zb[-1]-zs[-1]
  
  vMag = np.sqrt(np.dot(perturbVec[0,:],perturbVec[0,:])) #2-norm
  perturbVec[0,:] *= perturbMag/vMag
  
  #evolve perturbation ------------------------------------------
  for i in xrange(1,nBreedIter):
    xb,yb,zb = integrate_lorenz(nIntegrateIter, dt, x0+perturbVec[i-1,0], y0+perturbVec[i-1,1], z0+perturbVec[i-1,2])
    
    perturbVec[i,0] = xb[-1]-x0[0] #searching instability about x0
    perturbVec[i,1] = yb[-1]-y0[0]
    perturbVec[i,2] = zb[-1]-z0[0]

    vMag = np.sqrt(np.dot(perturbVec[i,:],perturbVec[i,:])) #2-norm
    perturbVec[i,:] *= perturbMag/vMag
    
  return perturbVec/perturbMag

def run_breed_lorenz():
  
  perturbVec = breed_lorenz()
  print perturbVec #limit of timestep tends to local tangent. can eliminate that direction to get other modes?

  fig = plt.figure()
  ax = fig.gca(projection='3d')
  
  ax.plot(perturbVec[:,0], perturbVec[:,1], perturbVec[:,2], 'b+')
  ax.set_xlabel("X Axis")
  ax.set_ylabel("Y Axis")
  ax.set_zlabel("Z Axis")
  ax.set_title("Lorenz Attractor")

  #plt.show()
  return plt
  
if __name__=='__main__':
  run_lorenz()
  run_breed_lorenz()
  
  plt.show()