#*nhyd_atmos.exe creates logs for each proc run
#we can see, say, the drift in max u over time...

import matplotlib.pyplot as plt

def example():
  #fpath = '/arctic1/nick/cases/163842/r2687/testDay/log.0000.err'
  fpath = '/arctic1/nick/cases/163842/testStart/log.0000.err'
  f = open(fpath,'r')
  
  minu = []; maxu = []; minw = []; maxw = [];
  line = 'notEmpty'
  while (line!=''): #readline returns empty string at EOF
    line = f.readline();
    if (line==''):
      continue
      
    if 'max u' in line:
      txt = line.strip().split()
      min = float(txt[-2]); max = float(txt[-1])
      minu.append(min); maxu.append(max);
    elif 'max w' in line:
      txt = line.strip().split()
      min = float(txt[-2]); max = float(txt[-1])
      minw.append(min); maxw.append(max);
  #
  f.close()

  x = range(len(minu)) #Number of iterations
  
  plt.figure(1)
  plt.suptitle('Min and max u')
  plt.subplot(2,1,1) #nrows, ncolumns, plot number within this subplot
  plt.plot(x,minu,'b.')
  plt.subplot(2,1,2)
  plt.plot(x,maxu,'r.')
  
  plt.figure(2)
  plt.suptitle('Min and max w')
  plt.subplot(2,1,1)
  plt.plot(x,minw,'b.')
  plt.subplot(2,1,2)
  plt.plot(x,maxw,'r.')
  
  plt.show()

def example_track():
  fname = '/arctic1/nick/cases/v1.0/x4/forSteven/results/2006-08-01_00/tracks/locs.txt'
  f = open(fname,'r')
  lat = []; lon = []; val = []
  lines = f.readlines()
  for s in lines:
    #print s
    if (s[0] == '['):
      print s
      lat.append( s.split()[1] )
      lon.append( s.split()[2] )
      val.append( s.split()[4] )
  f.close()
  print lat
  print lon
  print val



if __name__=='__main__':
  #example()
  example_track()
