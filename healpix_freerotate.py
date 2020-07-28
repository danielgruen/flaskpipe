#!/usr/bin/env python

import numpy as np
import healpy as H
from astropy.io import fits
from math import pi
import sys

if(len(sys.argv)!=6):
  print("syntax:",sys.argv[0],"[in_map] [theta_rot] [phi_rot] [out_map] [nside]")
  sys.exit(1)

print("")
print(sys.argv[1])
print("")

nside = int(sys.argv[5])
npix = H.nside2npix(nside)
hdu = fits.open(sys.argv[1])
map = np.zeros(npix)
for index in list(range(0,len(hdu[1].data["HPIX"]))):
  unmasked_pixel = hdu[1].data["HPIX"][index]
  map[unmasked_pixel] = hdu[1].data["FRACGOOD"][index]








#map = H.read_map(sys.argv[1])
#nside = H.npix2nside(len(map))
#npix = H.nside2npix(nside)
pix = np.arange(npix)
t,p = H.pix2ang(nside,pix) #theta, phi

r = H.Rotator(deg=False, rot=[float(sys.argv[2])*pi/180., float(sys.argv[3])*pi/180.])

map_rot = np.zeros(npix)

trot, prot = r(t,p)
j=H.ang2pix(nside,trot,prot)
map_rot[pix] = map[j]
H.write_map(sys.argv[4],map_rot,nest=False)

H.write_map("./m1_new.fits",map)
