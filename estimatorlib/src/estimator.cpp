#include "estimator.hpp"


Estimator::Estimator(string dat_f_dir,const Lattice& lat,const Lattice& br_lat):
    dir_to_write(dat_f_dir),lat(lat),br_lat(br_lat)
{
    length = lat.length;
    height = lat.height;

    file_precision=15;
}

void Estimator::close_wf()
{
    out_file.close();
}

