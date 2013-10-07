#If we're working with big data, we may not be able to load all we want into memory.
#So, we need to be able to find how much space we have available

import os
import numpy as np

szFloat = 8 #size of a float in bytes (assume 64 bit)

def findRAM():
  #return the RAM in bytes
  txt = os.popen("free -m -b").readlines()
  freeM = txt[1].split()[3] #pretty hacky but is there another option?
  
  return int(freeM)

def calcNDecomp(nVals):
  #return # of partitions needed to store all of these.
  #we need to account for python's data structures needing space too.
  
  fac = .8 #how much of free ram we can use. <1 but who knows how much...
  totMem = int(fac*findRAM())
  
  dataMem = nVals*szFloat #this is a lower bound!
  return 1+dataMem/totMem #integer division
  