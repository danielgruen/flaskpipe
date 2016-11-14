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
    cerr << "syntax: " << argv[0] << " [input catalog] [RA col] [DEC col] [mask]" << endl;
    return 1;
  }

  vector<string> columns;
  columns.push_back(argv[2]);
  columns.push_back(argv[3]);

  cerr << "// read input catalog" << endl;
  ObjectCollection input(argv[1], 1, columns);
  cerr << "// input size:" << input.size() << endl;

  cerr << "// read FITS mask" << endl;
  Healpix_Map<double> map(4,RING);
  read_Healpix_map_from_fits(argv[4],map,1,2);


  for(int i=0; i<input.size(); i++)
  {
        double ra=input[i]->doublePropertyValue(argv[2]);
        double dec=input[i]->doublePropertyValue(argv[3]);

        const pointing ang(M_PI/2.0 - dec*M_PI/180., ra*M_PI/180.);
        int j = map.ang2pix(ang);
        double p = map[j];

        if(p>1.e-10)
        {
          cout << ra << " " << dec << endl;
        }
        else if(rand()%1000==0) {
          cerr << "warning: mask at " << ra << " " << dec << " is zero; object skipped" << endl;
        }
  }

  return 0;
}

