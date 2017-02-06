
#include "random.hpp"

void O_Random::set_seed(int seed){
    // this is the right way and should work with newer versions of the compiler 
    //rg.seed( seed );
    // this one below is for the cluster which uses and older (2008) version of either the compiler
    // or the libraries that has a bug. This type cast deals with that bug.
    rg.seed(static_cast<unsigned int>(seed) );
}

/* This function gets a random integer in the given range "lower" - "upper".
 *
 * Even though the get_rand_bool could just be a special case of this with "lower" = 0 and "upper"=2
 * it would be more complicated to assign the "prob" weight. Here we just have an unweighted
 * probability distribution.
 * THIS DOES NOT INCLUDE RETURNING THE UPPER BOUND i.e [lower,upper)
 */
int O_Random::get_rand_int( int lowerLimit, int upperLimit ) {
    boost::uniform_int<> distribution( lowerLimit, upperLimit - 1);
    boost::variate_generator< RandomGeneratorType&, boost::uniform_int<> >
    LimitedInt( rg, distribution );
    return LimitedInt();
}

double O_Random::random01(void)
{
    static boost::uniform_01<boost::mt19937> dist(rg);
    return dist();
} 

