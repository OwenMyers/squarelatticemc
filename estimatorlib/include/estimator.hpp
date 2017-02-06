// loop_mc_tools.hpp
#define REAL 0
#define IMAG 1

//#include <iostream>
//#include <vector>
//#include <fstream>
//#include <math.h>
//#include <time.h>
//#include <boost/random.hpp>
//#include <boost/filesystem.hpp>
//
//#include <sstream>
//
#include <complex>
#include <fftw3.h>
#include <fstream>
#include "lattice.hpp"

#ifndef estimator_hpp_
#define estimator_hpp_

using namespace std;
typedef complex<double> dcmplx;



/******************************************************************************************************/
class Estimator
{
    public:

        //Estimator(Lattice &lat, string dat_f_dir);
        //Estimator(int length,int height,int seed,string dat_f_dir);
        Estimator(string dat_f_dir,const Lattice& lat,const Lattice& br_lat);

        Lattice lat;
        // background lat
        Lattice br_lat;

        // for writing to file that is just the current running average of the bins
        int number_bins_compleated_so_far;

        int length;
        int height;
        int seed;

        ofstream out_file;
        virtual void do_measurement()=0;

        // before write to file
        virtual void divide_hist_bins_by(int dev_by)=0;

        // after we have done the above (or not if N_measure = 1) we write hist_bins to file
        virtual void write()=0;

        // after we write something to file we need to clear the hist bins
        virtual void clear_hist_bins()=0;
          
        //tells us the directory in which the files will be written
        string dir_to_write;

        int file_precision;

        void close_wf(void);

};





#endif
