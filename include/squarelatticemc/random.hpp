
#include <boost/random.hpp>

#ifndef random_hpp_
#define random_hpp_

class O_Random
{
    public:
        typedef boost::mt19937 RandomGeneratorType;

        void set_seed(int seed);
        /* This function gets a random integer in the given range "lower" - "upper".
         *
         * Even though the get_rand_bool could just be a special case of this with "lower" = 0 and "upper"=2
         * it would be more complicated to assign the "prob" weight. Here we just have an unweighted
         * probability distribution.
         * THIS DOES NOT INCLUDE RETURNING THE UPPER BOUND i.e [lower,upper)
         */
        int get_rand_int( int lowerLimit, int upperLimit );

        double random01(void);
    private:
        RandomGeneratorType rg;
};

#endif
