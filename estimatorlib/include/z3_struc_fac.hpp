#include "2d_single_line_avg_out_estimator.hpp"
#include <fftw3.h>


#ifndef z3_struc_fac_hpp_
#define z3_struc_fac_hpp_

class Z3StructureFactor: public Estimator2DSingleLineOut<dcmplx> {
    public:
        Z3StructureFactor(string dat_f_dir,const Lattice& lat,const Lattice& br_lat);

        void do_measurement();

};
#endif
