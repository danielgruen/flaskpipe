#!/usr/bin/env python
import healpy as hp
import numpy as np
from math import pi
import sys

if(len(sys.argv)!=4):
  print("syntax:",sys.argv[0]," [input map] [new nside] [output map]")
  sys.exit(0)

in_map = hp.read_map(sys.argv[1])
in_nside = hp.npix2nside(len(in_map)) # nside of input map

out_nside = int(sys.argv[2])

out_map = np.ones(hp.nside2npix(out_nside))*1.e30

out_pixlist = np.arange(len(out_map),dtype='int')
out_coord   = hp.pix2ang(out_nside, out_pixlist, nest=False)

in_pixlist  = hp.ang2pix(in_nside,out_coord[0], out_coord[1])

out_map = in_map[in_pixlist]

hp.write_map(sys.argv[3],out_map)

