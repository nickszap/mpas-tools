#we're looking segment the dynamic tropopause. away from the equator,
#something like 2PVU differentiates between the tropopause and stratosphere with higher
#static stability.

import numpy as np

def findHighPV_regionGrow(seedh, seedLevel, pv, pvThreshold, hCellsOnCell, nCells, nLevels):
  #given the horizontal and vertical indices of a seed cell, return flag[cell,level]
  #with 1 if high, 0 if low
  #Reduce returning to search already visited neighbor of neighbor by only adding a
  #cell to the queue when it's first found
  
  flag = np.zeros((nCells,nLevels),dtype=int)
  wasSeed = np.zeros((nCells,nLevels),dtype=int)
  
  #make sure starting out with seed in right region
  if (pv[seedh,seedLevel]<pvThreshold):
    print "Can't region grow since seed not in proper region\n"
    flag[seedh,seedLevel] = 0
    return flag
  flag[seed[0],seed[1]]=1
  
  active = [(seedh,seedLevel)] #list of all cells above threshold
  while(len(active)>0):
    seed = active.pop() #If no index is specified, a.pop() removes and returns the last item in the list
    if (wasSeed[seed[0],seed[1]] ==1):
      continue
    wasSeed[seed[0],seed[1]] = 1
    
    #look horizontally
    hNbrs = hCellsOnCell[seed[0],:]
    for hCell in hNbrs:
      hInd = hCell; vInd = seed[1];
      if (pv[hInd,vInd]>=pvThreshold): #in right region
        if(flag[hInd,vInd]<1): #not visited
          #cell is above threshold and wasn't searched before
          flag[hInd,vInd] = 1; active.append((hInd,vInd))
      
    #look vertically
    if (seed[1]>0): #have below
      hInd = seed[0]; vInd = seed[1]-1;
      if (pv[hInd,vInd]>=pvThreshold): #in right region
        if(flag[hInd,vInd]<1): #not visited
          #cell is above threshold and wasn't searched before
          flag[hInd,vInd] = 1; active.append((hInd,vInd))
    
    if (seed[1]<nLevels-1): #have above
      hInd = seed[0]; vInd = seed[1]+1;
      if (pv[hInd,vInd]>=pvThreshold): #in right region
        if(flag[hInd,vInd]<1): #not visited
          #cell is above threshold and wasn't searched before
          flag[hInd,vInd] = 1; active.append((hInd,vInd))
          
  return flag