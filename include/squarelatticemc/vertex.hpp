
#ifndef vertex_hpp_
#define vertex_hpp_

#include <iostream>
#include <vector>

using namespace std;

class Vertex
{
    public:
        void set_num_link(int num_link_in);

        // set the sublattice 
        // 0 -> sublattice A
        // 1 -> sublattice B
        void set_sublat(int to_set);
        int get_sublat();
        // this will change but I'm a little paranoid about the memory allocation.
        int sublat;

        // For self avoiding random walk we well need to keep track of whether we have visited a
        // vertex or not. 
        // 1 --> we have
        // 0 --> we have not
        void set_been_here(int been_here_in);
        
        int been_here;

        // number of links per vertex
        int num_link;

        // declare vector containing spins at links
        // 1D array
        // links[0] ---> link to right of vertex
        // links[1] ---> link to top of vertex
        // links[2] ---> link to left of vertex
        // links[3] ---> link to bottom of vertex
        vector <int*> links;

        // A function that finds the number of occupied links at a vertex
        int num_oc_links();

        // Destructor
        //~Vertex();
};

#endif
