// pofn_bernoulli: get p(N) with Bernoulli re-masking

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
//#include <random>

#include <ctime>
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
  if(argc<4) {
    cerr << "syntax: " << argv[0] << " [healpix count map] [healpix theta_T-smoothed z-cut mask] [masking fraction threshold (above which regions are cut)]" << endl;
    return 1;
  }

  r = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(r, time (NULL) * getpid());

  // (1) parse command line
  
  mthresh=1.-atof(argv[3]);

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

  cerr << "// saving pixel counts" << endl;
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
    //cout << "adding " << c << " " << jackknife[i] << endl;

#pragma omp critical
    cylinders->appendObject(o);
  }

  // (3) get a random Bernoulli realization
  cerr << "// drawing Bernoulli" << endl;
  const int kmaskedcounts = cylinders->transformColumnNew("maskedcounts",maskedcounts,"counts","mask");

  cerr << *(cylinders->prototype) << endl;
  cerr << (*(*cylinders)[0]) << endl;

  // (4) prepare memory for p(N)
  const int Nmax  = cylinders->maximum("maskedcounts")+2;  // includes # of good pixels in last element
  cerr << "// allocating count tables out to Nmax=" << Nmax << endl;
  cerr << "// <maskedcounts> = " << cylinders->average("maskedcounts") << endl;
  cout << "# <maskedcounts> = " << cylinders->average("maskedcounts") << endl;
  int *c = new int[Nmax]; 
  for(int j=0; j<Nmax; j++) c[j] = 0;

  // (5) get counts statistics
  cerr << "// getting counts" << endl;

  c[Nmax-1] = cylinders->size();

  for(int i=0; i<c[Nmax-1]; i++)
      c[(*cylinders)[i]->intProperty[kmaskedcounts]]++;

  for(int i=0; i<Nmax-1; i++) {
    cout << i << " " << c[i] << endl;
  }

  return 0;
}

