#include "z3_struc_fac.hpp"

Z3StructureFactor::Z3StructureFactor(string dat_f_dir,const Lattice& lat,const Lattice& br_lat):
    Estimator2DSingleLineOut(dat_f_dir,lat,br_lat)
{
    outFileString = "z3_struc_fac.txt";

    // Here we do not open the file. It is rewritten every time because we are just outputting the
    // current running average.
    // file precision set in Estimator class
    out_file.setf( std::ios::fixed, std:: ios::floatfield );

}

void Z3StructureFactor::do_measurement()
{
    int phase;
    int check_phase;


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
    
            vector < vector <dcmplx> > cur_cor;
            for (int i = 0; i < length; i++)
            {
                vector <dcmplx> to_add;
                for (int j = 0; j < height; j++)
                {
                    to_add.push_back((dcmplx)0);
                }
                cur_cor.push_back(to_add);
            }

            for (int i = 0; i < length; i++)
            {
                for (int j = 0; j < height; j++)
                {

                    phase = lat.get_z3_vison_phase(alpha,beta,i+alpha,j+beta,br_lat);

                    dcmplx complex_plane_full(cos(((double) phase) * 2.0*M_PI/3.0),sin(((double) phase) * 2.0*M_PI/3.0));

                    cur_cor[j][i] += complex_plane_full; 


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
                        = cur_cor[j][i].real();

                    row_major_cur_cor[cur_index][IMAG]
                        = cur_cor[j][i].imag();
                    
                }
            }
            fftw_execute(fftwPlan);

            cur_index = 0;
            for (int i = 0; i < length; i++)
            {
                for (int j = 0; j < height; j++)
                {
                    cur_index = i*length+j;

                    dcmplx curFT(result[cur_index ][REAL],result[cur_index ][IMAG]); 
                    hist_bins[j][i] += curFT;
                }
            }
        }
    }
    fftw_destroy_plan(fftwPlan);
}
