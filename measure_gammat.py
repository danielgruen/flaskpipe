import treecorr
import sys
import healpy
import numpy as np
from math import pi

if(len(sys.argv)!=5):
  print("syntax:",sys.argv[0],"[lens map] [flask shear map] [lensing source mask] [output file]")
  sys.exit(1)


# read lenses
lens=healpy.read_map(sys.argv[1])
lensnside=healpy.npix2nside(len(lens))
lensid=np.arange(len(lens))[lens>0]
weight=lens[lens>0]
lens=0
ltheta,lphi=healpy.pix2ang(lensnside,lensid)
lensid=0
lra=  np.degrees(lphi)
ldec= np.degrees(pi/2.-ltheta)
ltheta=0
lphi=0
print("read",len(lra),"lenses in ra=",np.amin(lra),np.amax(lra),"dec=",np.amin(ldec),np.amax(ldec))


# read sources, apply mask
shearmask=healpy.read_map(sys.argv[3])
shear1= -healpy.read_map(sys.argv[2],1)[shearmask>0]
shear2=  healpy.read_map(sys.argv[2],2)[shearmask>0]
shearnside=healpy.npix2nside(len(shearmask))
shearid=np.arange(len(shearmask))[shearmask>0]
shearmask=0
stheta,sphi=healpy.pix2ang(shearnside,shearid)
shearid=0
sra=np.degrees(sphi)
sdec=np.degrees(pi/2.-stheta)
stheta=0
sphi=0
print("read",len(sra),"sources in ra=",np.amin(sra),np.amax(sra),"dec=",np.amin(sdec),np.amax(sdec))


lens_cat = treecorr.Catalog(ra=lra, dec=ldec, ra_units='degrees', dec_units='degrees', w=weight)
src_cat = treecorr.Catalog(ra=sra, dec=sdec,
                               g1=shear1, g2=shear2,
                               ra_units='degrees', dec_units='degrees')

#####################
# Measurement parameters
thmin = 5.
thmax = 600.
nth   = 24
bslop = 0.2
#####################

ng = treecorr.NGCorrelation(nbins = nth, min_sep = thmin, max_sep = thmax,
                                bin_slop = bslop, sep_units = 'arcmin', verbose=2)
ng.process_cross(lens_cat, src_cat)

np.savetxt(sys.argv[4], np.transpose([np.exp(ng.logr), ng.xi, ng.xi_im, ng.weight,]), '%.4e')
