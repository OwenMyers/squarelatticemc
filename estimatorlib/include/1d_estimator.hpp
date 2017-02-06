#include "estimator.hpp"


#ifndef estimator1d_hpp_
#define estimator1d_hpp_


template <typename measureDataType>
class Estimator1D: public Estimator
{
    public:

        Estimator1D ( string dat_f_dir, const Lattice& lat,const Lattice& br_lat):
            Estimator(dat_f_dir,lat,br_lat) 
        {
        
            hist_bins.resize(height-1);
        
        }
        
        vector <measureDataType> hist_bins;
        
        void divide_hist_bins_by(int dev_by)
        {
            for (int i = 0; i < hist_bins.size(); i++)
            {
                // vector of templated type. This will work even with complex hist_bin because double
                // complex type has properly overloaded operators. 
                hist_bins[i] = hist_bins[i]/((double) dev_by);
            }
        }
        void clear_hist_bins()
        {
            for (int i = 0; i < hist_bins.size(); i++)
            {
                // vector of templated type 
                hist_bins[i] = 0;
            }
        }
        void write()
        {
            for (int i = 0; i < hist_bins.size(); i++)
            {
                out_file << hist_bins[i] << " ";
            }
            out_file << endl; 
        }

};


#endif

