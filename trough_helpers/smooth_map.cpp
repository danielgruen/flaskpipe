// smooth_map: top-hat smooth input map 

#include <healpix_base.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <pointing.h>
#include <arr.h>
#include <fitshandle.h>

#include <iostream>
#include <omp.h>

using namespace std;

int main(int argc, char **argv)
{
  if(argc!=5) {
    cerr << "syntax: " << argv[0] << " [input map] [trough radius in arcmin] [output Nside] [output map]" << endl;
    return 1;
  }

  double thetat = atof(argv[2])/60./180.*M_PI;

  cerr << "// read input map" << endl;
  Healpix_Map<double> in(4,RING);
  read_Healpix_map_from_fits(argv[1],in,1,2);

  cerr << "// create healpix maps for output" << endl;
  Healpix_Map<double> smoothed(log2(atoi(argv[3])),RING);

#pragma omp parallel for
  for(int i=0; i<smoothed.Npix(); i++)
  {
        const pointing ang = smoothed.pix2ang(i); // (M_PI/2.0 - dec*M_PI/180., ra*M_PI/180.);
        vector<int> listpix;
        in.query_disc(ang, thetat, listpix);
        smoothed[i]=0;
        for(int j=0; j<listpix.size(); j++)
        {
           smoothed[i]+=in[listpix[j]]/double(listpix.size());
           if(in[listpix[j]]<0.99*Healpix_undef)
           {
             smoothed[i]=Healpix_undef;
             break;
           }
           else if (smoothed[i]<-1.e27) {
#pragma omp critical
             cout << "I just added " << in[listpix[j]] << " and now my map has gone bad." << Healpix_undef << endl;
           }
        }
  }

  fitshandle fhcount;
  fhcount.create(argv[4]);
  write_Healpix_map_to_fits(fhcount, smoothed, PLANCK_FLOAT64); 

  return 0;
}

