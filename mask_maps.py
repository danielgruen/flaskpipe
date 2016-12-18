import sys
import healpy
import numpy as np
import pyfits

if(len(sys.argv)!=5):
  print("syntax:",sys.argv[0],"[input map] [mask] [nextensions (1 or 3)] [output file]")
  sys.exit(1)



# read mask
mask=(healpy.read_map(sys.argv[2])==0)


# read maps

if(sys.argv[3]=="3"): # shear / T,p,q format
  map1=healpy.read_map(sys.argv[1],field=0)
  map2=healpy.read_map(sys.argv[1],field=1)
  map3=healpy.read_map(sys.argv[1],field=2)

  map1[mask]=0.
  map2[mask]=0.
  map3[mask]=0.

  healpy.write_map(sys.argv[4],[map1,map2,map3])

else:
  map1=healpy.read_map(sys.argv[1])
  map1[mask]=0.
  healpy.write_map(sys.argv[4],map1)
