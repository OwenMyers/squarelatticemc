#include "1d_estimator.hpp"

#ifndef dimer_dimer_direction_vert_hpp_
#define dimer_dimer_direction_vert_hpp_

class DimerDimerDirectionVert: public Estimator1D <double> {
    public:
        DimerDimerDirectionVert(string dat_f_dir,const Lattice& lat,const Lattice& br_lat);

        void do_measurement(void);

};
#endif
