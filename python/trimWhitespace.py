import os
#use ImageMagick to trim the whitespace, namely around ncl's png outputs
cmd = 'convert -trim +repage ./plotCFSR/cfsr_2pvu_2006.000001.png testTrim.png'
#os.system(cmd)

#For multiple files
cmd = 'mogrify -path ./ -trim +repage ./plotCFSR/cfsr_2pvu_2006.00000*.png'


