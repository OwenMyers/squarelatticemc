#include "dimer_horz_struc_fac.hpp"


DimerHorizontalStructureFactor::DimerHorizontalStructureFactor(string dat_f_dir,const Lattice& lat,const Lattice& br_lat):
    Estimator2DSingleLineOut(dat_f_dir,lat,br_lat)
{
    outFileString = "dmr_horz_struc_fac.txt";

    // Here we do not open the file. It is rewritten every time because we are just outputting the
    // current running average.
    // file precision set in Estimator class
    out_file.setf( std::ios::fixed, std:: ios::floatfield );

}

void DimerHorizontalStructureFactor::do_measurement()
{


    fftw_complex row_major_cur_cor[length*height];
    fftw_complex result[length*height];
    fftw_plan fftwPlan = fftw_plan_dft_2d(length,height,
                                        row_major_cur_cor,
                                        result,
                                        FFTW_FORWARD,
                                        FFTW_ESTIMATE);



    for (int alpha = 0; alpha < length; alpha++)
    {
        for (int beta = 0; beta < height; beta ++)
        {
    
            vector < vector <double> > cur_cor;
            for (int i = 0; i < length; i++)
            {
                vector <double> to_add;
                for (int j = 0; j < height; j++)
                {
                    to_add.push_back((double)0);
                }
                cur_cor.push_back(to_add);
            }


            int link00 = lat.get_vertex_link(alpha,beta,0);
            int cur_link;
            // If this is zero than we are just adding zero to everything so skip
            if (link00!=0)
            {
                for (int i = 0; i < height; i++)
                {
                    for (int j = 0; j < length; j++)
                    {
                        cur_link = lat.get_vertex_link(alpha+i,beta+j,0);
                        // for this it is also [y][x]
                        if (cur_link==1)
                        {
                            cur_cor[i][j]+=(double)1;
                        }
                    }
                }
            }


            // reshape the array into row major order
            int cur_index = 0;
            for (int i = 0; i < length; i++)
            {
                for (int j = 0; j < height; j++)
                {
                    cur_index = i*length+j;
                    row_major_cur_cor[cur_index][REAL]
                        = cur_cor[j][i];

                    row_major_cur_cor[cur_index][IMAG]
                        = 0.0;
                    
                }
            }
            fftw_execute(fftwPlan);

            cur_index = 0;
            for (int i = 0; i < length; i++)
            {
                for (int j = 0; j < height; j++)
                {
                    cur_index = i*length+j;

                    //dcmplx curFT(result[cur_index ][REAL],result[cur_index ][IMAG]); 
                    double curFT = result[cur_index ][REAL]; 
                    hist_bins[j][i] += curFT;
                }
            }
        }
    }
    fftw_destroy_plan(fftwPlan);
}
