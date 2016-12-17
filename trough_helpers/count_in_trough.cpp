// count_in_trough: count tracers in trough regions, save count as healpix maps

#include <healpix_base.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <pointing.h>
#include <arr.h>
#include <fitshandle.h>
#include <filter.h>

#include <iostream>
#include <omp.h>

using namespace std;

int main(int argc, char **argv)
{
  if(argc<4) {
    cerr << "syntax: " << argv[0] << " [tracer input catalog with 1=RA, 2=dec, redshift cut applied] [trough radius in arcmin] [output] (Nside)" << endl;
    return 1;
  }

  double thetat = atof(argv[2])/60./180.*M_PI;

  cerr << "// read tracer catalog" << endl;
  const int keys[] = {1, 2};
  const string names[] = {"ra","dec"}; 
  const Type types[] = {doubleType, doubleType};
  ifstream in(argv[1]);
  ObjectCollection tracers(in, keys, names, types, 2);

  cerr << "// create healpix maps for trough counts and mask" << endl;
  int rank=10;
  if(argc>4) rank=atoi(argv[4]);
  Healpix_Map<int> count(rank,RING);
  count.fill(0);

  if(thetat>0) {
  for(int i=0; i<tracers.size(); i++)
  {
        double ra=tracers[i]->doublePropertyValue(0);
        double dec=tracers[i]->doublePropertyValue(1);

        const pointing ang(M_PI/2.0 - dec*M_PI/180., ra*M_PI/180.);
        vector<int> listpix;
        count.query_disc(ang, thetat, listpix);
#pragma omp parallel for
        for(int j=0; j<listpix.size(); j++)
        {
           count[listpix[j]]++;
        }
        if(rand()%10000==0) cerr << "// added tracer at " << ra << " " << dec << endl;
  }
  } else {
  for(int i=0; i<tracers.size(); i++)
  {
        double ra=tracers[i]->doublePropertyValue(0);
        double dec=tracers[i]->doublePropertyValue(1);

        const pointing ang(M_PI/2.0 - dec*M_PI/180., ra*M_PI/180.);
        int pix=count.ang2pix(ang);
        count[pix]++;
        if(rand()%10000==0) cerr << "// added tracer at " << ra << " " << dec << endl;
  }
  }

  fitshandle fhcount;
  fhcount.create(argv[3]);
  write_Healpix_map_to_fits(fhcount, count, PLANCK_INT32); 

  return 0;
}

