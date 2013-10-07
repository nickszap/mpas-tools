#Module for running jobs on our cluster, Arctic.
#This is setup for use with passwordless ssh and mpdboot. 

#python libs
import os
import sys
from time import sleep

#my libs

#could also read these from ~/machines if they change!!!
#say, split machines line by ':'
nNodes = 10
nCorePerNode=12

def calc_numNodes(nProcs):
  #given # of processors to use, return the # of nodes needed
  nn,m = divmod(nProcs,nCorePerNode)
  if (m>0):
    nn = nn+1; #eg 20/12 gives 2. 5/12 gives 1
  #else, 24/12 gives 2
  return nn
  
def boot(np):
  #Given the total number of processors to use for the job, request the resources to run with.
  #need to boot # of physical machines+1 if booting off of Arctic. What do we need if booting off of node???
  nn = calc_numNodes(np)
  if (nn>nNodes):
    print "\n-----------------------\nUhoh. We need %s nodes to run the requested job, but only have %s\n" %(nn,nNodes)
    return()
  else:
    #cmd = "mpdboot -n %s -f ~/machines" %(nn)
    cmd = "time mpdboot -n %s -f ~/machines" %(nn+1) #machines+1 for mpdboot since head node
    print cmd
    os.system(cmd)

def runExe(exe, np):
  #Given resources from an mpdboot, run the executable
  #cmd = "mpiexec -machinefile ~/machines -n %s %s" %(np, exe)
  cmd = "time mpiexec -machinefile ~/machines -n %s %s" %(np, exe)
  print cmd
  os.system(cmd)

def main(exe, np):
  #resources and run
  boot(np);
  runExe(exe,np)
  
  #we won't return unless execution returns, right?
  print("\nJob finished!\n")
  #can notify with text message! :)
  #sleep(5)
  os.system('mpdallexit')

if __name__=='__main__':
  #call from command line as: python arctic.py exeWithPath nprocs.
  #this is quick and dirty w/o error checking...maybe use argparse
  main(sys.argv[1], int(sys.argv[2]))

  
