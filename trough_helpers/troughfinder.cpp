// troughfinder: find troughs from healpix count / mask maps

#include <healpix_base.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <pointing.h>
#include <arr.h>
#include <fitshandle.h>
#include <filter.h>

#include <iostream>
#include <omp.h>

#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

double mthresh;
vector<int> countcut;

gsl_rng *r=0;

long long int maskedcounts(long long int count, double mask)
// return Binomial-masked random realization of counts at good-area fraction mthresh (must be equal or below mask)
{
  //binomial_distribution<int> d(count,1.-mask+mthresh); // chance the individual galaxy gets hit by a random mask
  //return d(rng);
  return gsl_ran_binomial(r, 1.-mask+mthresh, count);
}

int stat_index(int count)
{
  int idx=0;
  int i=0;
  while(count>countcut[i] && i<countcut.size()){
   idx+=2;
   i++;
  }
  if(i<countcut.size() && count==countcut[i]) idx++;
  return idx;
}

int main(int argc, char **argv)
{
  if(argc<6) {
    cerr << "syntax: " << argv[0] << " [healpix count map] [healpix theta_T-smoothed z-cut mask] [masking fraction threshold (above which regions are cut)] [output prefix] [limiting percentiles 1 2 ...]" << endl;
    return 1;
  }

  r = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(r, time (NULL) * getpid());

  // (1) parse command line
  
  mthresh=1.-atof(argv[3]);
  vector<double> pcut;
  for(int i=5; i<argc; i++) {
   pcut.push_back(atof(argv[i]));
   if(pcut.size()>1) assert(pcut[pcut.size()-1]>pcut[pcut.size()-2]); // make sure it's increasing at least a bit
  }

  // (2) read counts into ObjectCollection

  cerr << "// read count map" << endl;
  Healpix_Map<int> counts(4,RING);
  read_Healpix_map_from_fits(argv[1],counts,1,2);

  cerr << "// read mask" << endl;
  Healpix_Map<double> mask(4,RING);
  read_Healpix_map_from_fits(argv[2],mask,1,2);

  // just to make sure formats are the same
  assert(mask.Npix()==counts.Npix());
  assert(mask.Scheme()==counts.Scheme());

  int keys[] = {0,0,0};
  string names[] = {"pix","counts","mask"};
  Type types[] = {integerType, integerType, doubleType};

  ObjectCollection *cylinders;
  {
  ObjectPrototype cylinderPrototype(keys,names,types,3);
  cylinders = new ObjectCollection(&cylinderPrototype);
  }

  int s=0; // pixels in survey

#pragma omp parallel for reduction(+:s)
  for(int i=0; i<mask.Npix(); i++)
  {
    double m=mask[i];
    if(m>0) s=s+1;
    if(m<mthresh) continue;

    int c=counts[i];
    Object *o = new Object(cylinders->prototype);
    o->setIntProperty(0,i);
    o->setIntProperty(1,c);
    o->setDoubleProperty(0,m);

#pragma omp critical
    cylinders->appendObject(o);
  }

  cerr << "// accepted " << cylinders->size() << "/" << s << " pixels in survey" << endl;

  // (3) determine thresholds from one homogeneously masked realization
 
   
  cylinders->transformColumnNew("maskedcounts",maskedcounts,"counts","mask");
  cylinders->sort("maskedcounts");

  for(int i=0; i<pcut.size(); i++)
  {
    countcut.push_back((*cylinders)[pcut[i]*cylinders->size()]->intPropertyValue(2));
    cerr << "// percentile " << pcut[i] << " corresponds to count " << countcut[i] << endl;
    if(i>1) assert(countcut[i]>countcut[i-2]);
  }
  pcut.push_back(1.0);

  // (4) count how many times each cylinder ends up in which bin
  
  // we need to be quick and memory efficient here, so that's what we'll do

  cylinders->sort("pix"); // sorting by count causes unequal computation times below

  int **stats = new int*[cylinders->size()];
  for(int i=0; i<cylinders->size(); i++) {
    stats[i] = new int[2*countcut.size()+1];
#pragma omp parallel
    for(int j=0; j<2*countcut.size()+1; j++)
      stats[i][j]=0;
  }
  int ntrials=1000;

  cerr << "// running trials " << flush;

  for(int i=0; i<ntrials; i++)
  {
    // first do stats
#pragma omp parallel for schedule(static)
    for(int j=0; j<cylinders->size(); j++) {
      stats[j][stat_index((*cylinders)[j]->intPropertyValue(2))]++;
    }
    
    // new dice if necessary
    if(i<ntrials-1) cylinders->transformColumnNew("maskedcounts",maskedcounts,"counts","mask");

    cerr << "." << flush;
  }
  cerr << " done." << endl;

  int cbin[2*countcut.size()+1]; // # of times of count being in each bin
  double pbin[2*countcut.size()+1]; // probability of count being in each bin
#pragma omp parallel for
  for(int i=0; i<2*countcut.size()+1; i++)
  {
    cbin[i] = 0;
    for(int j=0; j<cylinders->size(); j++) {
      cbin[i] += stats[j][i];
    } 
    pbin[i] = double(cbin[i])/double(ntrials)/double(cylinders->size());
  }
 
  double psum=0.;
  unsigned long long nsum=0; 
  for(int i=0; i<2*countcut.size()+1; i++)
  {
    cerr << "// pbin[" << i << "]=" << pbin[i] << endl;
    psum += pbin[i];
    nsum += cbin[i];
  }
  cerr << "// sum(pbin)=" << psum << " sum(cbin)=" << nsum << " ntrials=" << double(ntrials)*double(cylinders->size()) << endl;

  // (5) express that in a probability the cylinder ends up in the requested percentile range and write to healpix

  double wlimit_used=0.; // weight at lower limit of bin used up previously
  double pbelow=0.;      // percentile cut below

  for(int i=0; i<countcut.size()+1; i++)
  {
    int idx1=i*2;
    int idx2=idx1+1;
    int idx0=idx1-1;

    double wbin[3];         // weight applied to each of the count bins for this percentile range
                            // the requirement is that sum(wbin[i]*pbin[i])=pmax-pbelow
    wbin[1]=1.;             // inside bin
    wbin[0]=1.-wlimit_used; // at lower boundary

    if(idx2<2*countcut.size()+1) { 
      wbin[2]=((pcut[i]-pbelow)-wbin[1]*pbin[idx1])/pbin[idx2];
      if(idx0>=0) wbin[2] -= wbin[0]*pbin[idx0]/pbin[idx2];
    }
    else wbin[2]=0.;

    cerr << "// percentile range " << i << " has wbin[0,1,2]=" << wbin[0] << "," << wbin[1] << "," << wbin[2] << endl;

    Healpix_Map<double> map;
    map.SetNside(mask.Nside(),mask.Scheme());
    double wsum=0.;

    if(i==0) {
#pragma omp parallel for reduction(+:wsum)
      for(int j=0; j<cylinders->size(); j++)
      {
        map[(*cylinders)[j]->intPropertyValue(0)] = (wbin[1]*double(stats[j][idx1]) + wbin[2]*double(stats[j][idx2]))/double(ntrials);
        wsum += map[(*cylinders)[j]->intPropertyValue(0)];
      }
    }
    else if(i==countcut.size()) {
#pragma omp parallel for reduction(+:wsum)
      for(int j=0; j<cylinders->size(); j++)
      {
        map[(*cylinders)[j]->intPropertyValue(0)] = (wbin[1]*double(stats[j][idx1]) + wbin[0]*double(stats[j][idx0]))/double(ntrials);
        wsum += map[(*cylinders)[j]->intPropertyValue(0)];
      }
    }
    else {
#pragma omp parallel for reduction(+:wsum)
      for(int j=0; j<cylinders->size(); j++)
      {
        map[(*cylinders)[j]->intPropertyValue(0)] = (wbin[1]*double(stats[j][idx1]) + wbin[0]*double(stats[j][idx0]) + wbin[2]*double(stats[j][idx2]))/double(ntrials);
        wsum += map[(*cylinders)[j]->intPropertyValue(0)];
      }
    }

    fitshandle fh;
    fh.create(string(argv[4])+"_"+FilterFunctions::NumberToString(i)+".fits");
    write_Healpix_map_to_fits(fh, map, PLANCK_FLOAT32);

    wlimit_used=wbin[2];
    cerr << "// percentile range " << i << " has total p=" << wsum/cylinders->size() << " where it should be " << pcut[i]-pbelow << endl;
    pbelow=pcut[i];
  }
  
  // (6) finally write trough mask 

  Healpix_Map<double> troughmask;
  troughmask.SetNside(mask.Nside(),mask.Scheme());

#pragma omp parallel
  for(int j=0; j<cylinders->size(); j++)
  {
    troughmask[(*cylinders)[j]->intPropertyValue(0)] = 1.;   
  }
  fitshandle fh;
  fh.create(string(argv[4])+"_troughmask.fits");
  write_Healpix_map_to_fits(fh, troughmask, PLANCK_FLOAT32);

  return 0;
}

