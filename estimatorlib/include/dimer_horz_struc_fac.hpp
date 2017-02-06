#include "2d_single_line_avg_out_estimator.hpp"
#include <fftw3.h>

#ifndef dimer_horz_struc_fac_hpp_
#define dimer_horz_struc_fac_hpp_

class DimerHorizontalStructureFactor: public Estimator2DSingleLineOut<double> {
    public:
        DimerHorizontalStructureFactor(string dat_f_dir,const Lattice& lat,const Lattice& br_lat);

        void do_measurement();

};
#endif
