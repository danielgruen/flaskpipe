import treecorr
import sys
import healpy
import numpy as np
from math import pi,sqrt

if(len(sys.argv)!=6):
  print("syntax:",sys.argv[0],"[lens map] [sigma^2_gamma [arcmin^2]] [lensing source mask] [output file] [random seed]")
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


# read mask
shearmask=healpy.read_map(sys.argv[3])
shearnside=healpy.npix2nside(len(shearmask))
print("shear_nside=",shearnside)
shearid=np.arange(len(shearmask))[shearmask>0]
shearweight=shearmask[shearmask>0]
soffset=math.sqrt(healpy.nside2pixarea(healpy.npix2nside(len(shearmask)),degrees=True)) # pixel width in deg
shearmask=0
stheta,sphi=healpy.pix2ang(shearnside,shearid)
shearid=0
sra=np.degrees(sphi)
sdec=np.degrees(pi/2.-stheta)
stheta=0
sphi=0

# randomize source positions a tiny bit
sra = sra + np.random.uniform(low=-soffset, high=+soffset)/np.cos(sdec*pi/180.)
sdec = sdec + np.random.uniform(low=-soffset, high=+soffset)


sigmagamma=sqrt(float(sys.argv[2])/(3600.*healpy.pixelfunc.nside2pixarea(shearnside,degrees=True)))

print("simulating",len(sra),"sources in ra=",np.amin(sra),np.amax(sra),"dec=",np.amin(sdec),np.amax(sdec)," with sigma_gamma=",sigmagamma)
print("by the way, I think your input sigma^2_gamma should be",0.272**2/2)
print("and, if you care to know, one of your pixels has an area of",(3600.*healpy.pixelfunc.nside2pixarea(shearnside,degrees=True)),"arcmin^2")

np.random.seed(int(sys.argv[5]))
shear1=np.random.normal(scale=sigmagamma,size=len(sra))
shear2=np.random.normal(scale=sigmagamma,size=len(sra))


lens_cat = treecorr.Catalog(ra=lra, dec=ldec, ra_units='degrees', dec_units='degrees', w=weight)
src_cat = treecorr.Catalog(ra=sra, dec=sdec, w=shearweight,
                               g1=shear1, g2=shear2,
                               ra_units='degrees', dec_units='degrees')

#####################
# Measurement parameters
thmin = 5.
thmax = 600.
nth   = 24
bslop = 0.1
#####################

ng = treecorr.NGCorrelation(nbins = nth, min_sep = thmin, max_sep = thmax,
                                bin_slop = bslop, sep_units = 'arcmin', verbose=2)
ng.process_cross(lens_cat, src_cat)

np.savetxt(sys.argv[4], np.transpose([np.exp(ng.logr), ng.xi, ng.xi_im, ng.weight,]), '%.4e')
