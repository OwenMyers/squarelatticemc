#include "estimator.hpp"


#ifndef full_lat_dimer_single_line_out_estimator_hpp_
#define full_lat_dimer_single_line_out_estimator_hpp_

template <typename measureDataType>
class EstimatorFullLatDmrSingleLineOut: public Estimator
{
    public:

        EstimatorFullLatDmrSingleLineOut ( string dat_f_dir, const Lattice& lat,const Lattice& br_lat): 
            Estimator(dat_f_dir,lat,br_lat) 
        {
            number_bins_compleated_so_far = 0;
        
            hist_bins.resize(height);
            for (int i = 0; i<height; i++)
            {
                hist_bins[i].resize(length);
                for (int j = 0; j<length; j++)
                {
                    hist_bins[i][j].resize(2);
                }
            }
            // Now initialize all values to zero
            for (int i = 0; i<height; i++)
            {
                for (int j = 0; j<length; j++)
                {
                    for (int k = 0; k<2; k++)
                    {
                        hist_bins[i][j][k] = (double)0;
                    }
                }
            }

            bin_avg.resize(height);
            for (int i = 0; i<height; i++)
            {
                bin_avg[i].resize(length);
                for (int j = 0; j<length; j++)
                {
                    bin_avg[i][j].resize(2);
                }
            }
            // Now initialize all values to zero
            for (int i = 0; i<height; i++)
            {
                for (int j = 0; j<length; j++)
                {
                    for (int k = 0; k<2; k++)
                    {
                        bin_avg[i][j][k] = (double)0;
                    }
                }
            }        
        }
        
        int number_bins_compleated_so_far;
        string outFileString;
        vector < vector < vector <measureDataType> > > hist_bins;
        vector < vector < vector <measureDataType> > > bin_avg;
        
        void divide_hist_bins_by(int dev_by)
        {
            number_bins_compleated_so_far++;

            for (int i = 0; i < hist_bins.size(); i++)
            {
                for (int j = 0; j < hist_bins[0].size(); j++)
                {

                    for (int k = 0; k < 2; k++)
                    {
                        // vector of templated type. This will work even with complex hist_bin because double
                        // complex type has properly overloaded operators. 
                        hist_bins[i][j][k] = hist_bins[i][j][k]/((double) dev_by);
                        bin_avg[i][j][k] += hist_bins[i][j][k];
                    }
                }
            }
        }
        void clear_hist_bins()
        {
            for (int i = 0; i < hist_bins.size(); i++)
            {
                for (int j = 0; j < hist_bins[0].size(); j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        // vector of templated type 
                        hist_bins[i][j][k] = (double)0;
                    }
                }
            }
        }
        void write()
        {
            out_file.open((dir_to_write + outFileString).c_str(), std::ios::out);
            for (int i = 0; i < hist_bins.size(); i++)
            {
                for (int j = 0; j < hist_bins[0].size(); j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        out_file << bin_avg[i][j][k]/((double)number_bins_compleated_so_far) << " ";
                    }
                }
            }
            out_file << endl; 
            out_file.close();
        }


};


#endif

