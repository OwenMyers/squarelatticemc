#include "estimator.hpp"
#include <complex>

typedef complex<double> dcmplx;


#ifndef twod_single_line_ave_out_estimator_hpp_
#define twod_single_line_ave_out_estimator_hpp_

template <typename measureDataType>
class Estimator2DSingleLineOut: public Estimator
{
    public:

        Estimator2DSingleLineOut ( string dat_f_dir, const Lattice& lat,const Lattice& br_lat): 
            Estimator(dat_f_dir,lat,br_lat) 
        {
            number_bins_compleated_so_far = 0;
        
            for (int i = 0; i<height; i++)
            {
                vector <measureDataType> to_add;
                for (int j = 0; j<length; j++)
                {
                    to_add.push_back((measureDataType)0);
                }
                hist_bins.push_back(to_add);
            }
            for (int i = 0; i<height; i++)
            {
                vector <measureDataType> to_add;
                for (int j = 0; j<length; j++)
                {
                    to_add.push_back((measureDataType)0);
                }
                bin_avg.push_back(to_add);
            }
        
        }
        
        int number_bins_compleated_so_far;
        string outFileString;
        vector < vector <measureDataType> > hist_bins;
        vector < vector <measureDataType> > bin_avg;
        
        void divide_hist_bins_by(int dev_by)
        {
            number_bins_compleated_so_far++;

            for (int i = 0; i < hist_bins.size(); i++)
            {
                for (int j = 0; j < hist_bins[0].size(); j++)
                {
                    // vector of templated type. This will work even with complex hist_bin because double
                    // complex type has properly overloaded operators. 
                    hist_bins[i][j] = hist_bins[i][j]/((double) dev_by);
                    bin_avg[i][j] += hist_bins[i][j];
                }
            }
        }
        void clear_hist_bins()
        {
            for (int i = 0; i < hist_bins.size(); i++)
            {
                for (int j = 0; j < hist_bins[0].size(); j++)
                {
                    // vector of templated type 
                    hist_bins[i][j] = (double)0;
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
                    out_file << bin_avg[i][j]/((double)number_bins_compleated_so_far) << " ";
                }
            }
            out_file << endl; 
            out_file.close();
        }


};


#endif

