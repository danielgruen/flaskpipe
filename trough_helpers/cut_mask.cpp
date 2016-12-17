// cut_mask: cut input catalog by healpix mask

#include <healpix_base.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <pointing.h>
#include <arr.h>
#include <fitshandle.h>
#include <filter.h>

#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
  if(argc!=5) {
    cerr << "syntax: " << argv[0] << " [input catalog] [phi col] [theta col] [mask]" << endl;
    return 1;
  }

  vector<string> columns;
  columns.push_back(argv[2]);
  columns.push_back(argv[3]);

  cerr << "// read input catalog" << endl;
  ObjectCollection input(argv[1], 1, columns);
  cerr << "// input size:" << input.size() << endl;
  cerr << "// input header:" << *(input.prototype) << endl;

  cerr << "// read FITS mask" << endl;
  Healpix_Map<double> map(4,RING);
  read_Healpix_map_from_fits(argv[4],map,1,2);


  for(int i=0; i<input.size(); i++)
  {
        double phi=input[i]->doublePropertyValue(argv[2]);
        double theta=input[i]->doublePropertyValue(argv[3]);
        double ra=phi;
        double dec=90.-theta;

        const pointing ang(theta*M_PI/180., phi*M_PI/180.);
        int j = map.ang2pix(ang);
        double p = map[j];

        if(p>1.e-10)
        {
          cout << ra << " " << dec << endl;
        }
        else if(rand()%100000==0) {
          cerr << "warning: mask at " << ra << " " << dec << " is zero; object skipped" << endl;
        }
  }

  return 0;
}

