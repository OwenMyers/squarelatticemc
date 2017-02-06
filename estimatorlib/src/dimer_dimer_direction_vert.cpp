#include "dimer_dimer_direction_vert.hpp"

DimerDimerDirectionVert::DimerDimerDirectionVert(string dat_f_dir,const Lattice& lat,const Lattice& br_lat):
    Estimator1D(dat_f_dir,lat,br_lat)
{

    // open the actual file you want to write to.
    // file precision set in Estimator class
    out_file.setf( std::ios::fixed, std:: ios::floatfield );
    out_file.open((dir_to_write + "vert_parallel_data.txt").c_str());

}

// Add to hist bins every measure -> must take care of division (dividing hist_bin entries by
// N_measure) by number of runs added somewhere else. clearing his bins after a write must be done as
// well
void DimerDimerDirectionVert::do_measurement()
{

    // sweep link for 
    int cur_link;
    int link00;
      
    double before_add;
    for (int alpha = 0; alpha < length; alpha++)
    {
        for (int beta = 0; beta < height; beta++)
        {
            // first link of first col. always using the link to the right (.links[0])
            link00 = lat.get_vertex_link(alpha,beta,0);

            // VERY IMPORTANT using two adjacent columns does not work. Leaving this line commented out
            // below as a reminder
            // first link of first col. always using the link to the right (.links[0])
            //int link01 = * lattice[0][1].links[0];

            // Just start at the same vertex every time. Choose it to be 0,0 vertex.
            // Look for parallels links and add them to the right bin
            if (link00 == 1)
            {
                for (int i = 1; i < (height); i++)
                {
                    cur_link = lat.get_vertex_link(alpha,beta+i,0);

                    if (cur_link == 1) 
                    { 
                        //before_add = ((double)1 )/((double)(loop_over_length*loop_over_height));
                        //cout << "before add" << endl;
                        //std::cout.precision(10);
                        //std::cout.setf( std::ios::fixed, std:: ios::floatfield );
                        //cout << before_add << endl;
                        hist_bins[i-1] += ((double)1 )/((double)(length*height));
                    }
                }
            }
        }
    }
}

