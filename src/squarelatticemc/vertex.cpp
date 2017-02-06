#include "vertex.hpp"

// set the sublattice 
// 0 -> sublattice A
// 1 -> sublattice B
void Vertex::set_sublat(int to_set)
{
    sublat = to_set;
}
int Vertex::get_sublat()
{
    return sublat;
}

int Vertex::num_oc_links()
{
    int count = 0;
    for (int i = 0; i<4; i++)
    {
        if (* links[i] == 1) { count += 1;}
    }

    return count;
}

void Vertex::set_been_here(int been_here_in)
{
    been_here = been_here_in;
}
