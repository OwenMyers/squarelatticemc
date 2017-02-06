#include "full_lat_dimer_single_line_out_estimator.hpp"

#ifndef dimer_density_hpp_
#define dimer_density_hpp_

class DimerDensity: public EstimatorFullLatDmrSingleLineOut<double> {
    public:
        DimerDensity(string dat_f_dir,const Lattice& lat,const Lattice& br_lat);

        void do_measurement();

};
#endif
