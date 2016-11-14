# config

CFLAGS_MACBOOK=-L/sw/lib -I~/include -I/Users/dgruen/daniel/work/filter -L/Users/dgruen/Downloads/Healpix_3.30/src/cxx/basic_gcc/lib/ -I/Users/dgruen/Downloads/Healpix_3.30/src/cxx/basic_gcc/include/ -L/home/dgruen/lib -fopenmp -std=c++11
CFLAGS_XOC=-L/sw/lib -I/home/dgruen/include -I/home/dgruen/work/filter -L/home/dgruen/Downloads/Healpix_3.30/src/cxx/basic_gcc/lib/ -I/home/dgruen/Downloads/Healpix_3.30/src/cxx/basic_gcc/include/ -L/home/dgruen/lib -fopenmp -std=c++11
CFLAGS_KILS=-L/home/dgruen/Downloads/Healpix_3.30/src/cxx/basic_gcc/lib/ -I/home/dgruen/Downloads/Healpix_3.30/src/cxx/basic_gcc/include/ -fopenmp -std=c++11

FILTER_O_MACBOOK=/Users/dgruen/daniel/work/filter/filter.o -lcfitsio -lCCfits
FILTER_O_XOC=/home/dgruen/work/filter/filter.o -lcfitsio -lCCfits
FILTER_O_KILS=/u/ki/users/dgruen/work/filter/filter.o -lcfitsio -lCCfits

LIBHEALPIX=-lhealpix_cxx -lcxxsupport -lc_utils -lfftpack -lsharp -lcfitsio
LIBGSL=-lgsl

CFLAGS=$(CFLAGS_KILS)
FILTER_O=$(FILTER_O_KILS)

all: cut_mask

other: mkhealpixmask mkhealpixmask_cutbug cut_mask_z add_healpix count_in_trough counts_to_maskedcounts mask_trough troughfinder troughfinder_nomask troughfinder_strictmask pofn_bernoulli_jackknife pofn_bernoulli_jackknife_tworadii matchshapez smooth_map smooth_map_NEST 
		# pofn_bernoulli_jackknife pofn_bernoulli_jackknife_tworadii matchshapez

cut_mask: cut_mask.cpp
	g++ $(CFLAGS) -o cut_mask cut_mask.cpp $(FILTER_O) $(LIBHEALPIX)


mkhealpixmask: mkhealpixmask.cpp
	g++ $(CFLAGS) -o mkhealpixmask mkhealpixmask.cpp $(FILTER_O) $(LIBHEALPIX)
mkhealpixmask_cutbug: mkhealpixmask_cutbug.cpp
	g++ $(CFLAGS) -o mkhealpixmask_cutbug mkhealpixmask_cutbug.cpp $(FILTER_O) $(LIBHEALPIX)
cut_mask_z: cut_mask_z.cpp
	g++ $(CFLAGS) -o cut_mask_z cut_mask_z.cpp $(FILTER_O) $(LIBHEALPIX)
count_in_trough: count_in_trough.cpp
	g++ $(CFLAGS) -o count_in_trough count_in_trough.cpp $(FILTER_O) $(LIBHEALPIX)
counts_to_maskedcounts: counts_to_maskedcounts.cpp
	g++ $(CFLAGS) -o counts_to_maskedcounts counts_to_maskedcounts.cpp $(FILTER_O) $(LIBHEALPIX)
add_healpix: add_healpix.cpp
	g++ $(CFLAGS) -o add_healpix add_healpix.cpp $(FILTER_O) $(LIBHEALPIX)
mask_trough: mask_trough.cpp
	g++ $(CFLAGS) -o mask_trough mask_trough.cpp $(FILTER_O) $(LIBHEALPIX)
troughfinder: troughfinder.cpp
	g++ $(CFLAGS) -o troughfinder troughfinder.cpp $(FILTER_O) $(LIBHEALPIX)
smooth_map: smooth_map.cpp
	g++ $(CFLAGS) -o smooth_map smooth_map.cpp $(LIBHEALPIX)
smooth_map_NEST: smooth_map_NEST.cpp
	g++ $(CFLAGS) -o smooth_map_NEST smooth_map_NEST.cpp $(LIBHEALPIX)
troughfinder_nomask: troughfinder_nomask.cpp
	g++ $(CFLAGS) -o troughfinder_nomask troughfinder_nomask.cpp $(FILTER_O) $(LIBHEALPIX)
troughfinder_strictmask: troughfinder_strictmask.cpp
	g++ $(CFLAGS) -o troughfinder_strictmask troughfinder_strictmask.cpp $(FILTER_O) $(LIBHEALPIX)
pofn_bernoulli_jackknife: pofn_bernoulli_jackknife.cpp
	g++ $(CFLAGS) -o pofn_bernoulli_jackknife pofn_bernoulli_jackknife.cpp $(FILTER_O) $(LIBHEALPIX)
pofn_bernoulli_jackknife_tworadii: pofn_bernoulli_jackknife_tworadii.cpp
	g++ $(CFLAGS) -o pofn_bernoulli_jackknife_tworadii pofn_bernoulli_jackknife_tworadii.cpp $(FILTER_O) $(LIBHEALPIX) $(LIBGSL)
matchshapez: matchshapez.cpp
	g++ -g $(CFLAGS) -o matchshapez matchshapez.cpp $(FILTER_O)
