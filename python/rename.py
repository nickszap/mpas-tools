import glob
import os

def get_filenames(fpath, matchString):
  #return a list of files that match the specified format.
  #use absolute fpath so can link to file in yet another directory
  return sorted(glob.glob(fpath+matchString))
  
def nameFiles(fileList, fpathOut):
  #a file like plot1.png gets renamed as plot1-1.png
  #file is treated as name.ext where name can have '.'s as well
  
  for iFile, fpath in enumerate(fileList):
    #get the filename out of whatever path the file is in
    f = fpath.split('/')[-1]
    
    #get the extension
    splitName = f.split('.')
    ext = splitName[-1]
    lenExt = len(ext)+1 #characters +1 for .
    #fBase = f[0:-lenExt] #for full name minus .extension
    fBase = splitName[-3] #for say f=./figs/theta_2pvu_20060801_t_x7.000003.png
    
    fnameOut = fpathOut+fBase+'-'+str(iFile)+'.'+ext
    #print fnameOut
    cmd = 'ln -s '+fpath +' '+fnameOut
    
    os.system(cmd)
    
def runRename():
  
  #fileList = get_filenames('/data01/20060801/figs/', 'theta_2pvu_20060801_t_x7*')
  #fileList = get_filenames('/data01/20060801/figs/', 'theta_2pvu_20060801_kf_x7*')
  #fdir = '/arctic1/nick/cases/v1.0/x4/forSteven/results/2006-08-01_00/t200/'
  #fileList = get_filenames(fdir,'x4_t*')
  fdir = '/arctic1/nick/cases/cfsr/press/nomads/figs/'
  fileList = get_filenames(fdir,'cfsr_2pvu_*')
  print fileList
  
  nameFiles(fileList, './forLatex/')
  
if __name__=='__main__':
  runRename()
