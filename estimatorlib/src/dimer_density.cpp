#include "dimer_density.hpp"


DimerDensity::DimerDensity(string dat_f_dir,const Lattice& lat,const Lattice& br_lat):
    EstimatorFullLatDmrSingleLineOut(dat_f_dir,lat,br_lat)
{
    outFileString = "dimer_density_data.txt";

    // Here we do not open the file. It is rewritten every time because we are just outputting the
    // current running average.
    // file precision set in Estimator class
    out_file.setf( std::ios::fixed, std:: ios::floatfield );

}

void DimerDensity::do_measurement()
{

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < length; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                int cur_link = lat.get_vertex_link(j,i,k);
                // for this it is also [y][x]
                if (cur_link==1)
                {
                    hist_bins[i][j][k]+=(double)1;
                }
            }
        }
    }
    
}
