#ifndef lattice_hpp_
#define lattice_hpp_

#include "vertex.hpp"
#include "random.hpp"

#include <iostream>
#include <vector>

using namespace std;

/* be really careful putting together the Vertex
 * objects so that the pointers overlap correctly in the "image"  O's are vertices \ and - are links
 * and X is a spin
 * X   X   X   X   X
 * |   |   |   |   |   
 * O-X-O-X-O-X-O-X-O-X
 * |   |   |   |   |   
 * X   X   X   X   X   
 * |   |   |   |   |   
 * O-X-O-X-O-X-*-X-O-X
 * |   |   |   |   |   
 * X   X   X   X   X   
 * |   |   |   |   |   
 * O-X-O-X-O-X-O-X-O-X
 * |   |   |   |   |  
 * X   X   X   X   X  
 * |   |   |   |   |  
 * O-X-O-X-O-X-O-X-O-X
 *
 *
 * The "*" is located at arr[2][3]
 */

class Lattice
{
    public:
        Lattice();

        Lattice(int to_set_length, int to_set_height, O_Random rdm_in);

        O_Random rdm;

        void check_for_defects(void);

        int length;
        int height;
        // total number of vertices
        int num_vert;
        int links_per_vert;
        int tot_num_links;

        int get_z3_vison_phase(int start_x, int start_y, int end_x, int end_y,const Lattice& br_lat);

        // ****************************************************************************************************************************************
        // ***************** some functions to make the access of the lattice more accessible ******************************************************
        // ****************************************************************************************************************************************
        // "get_vertex_link(<x location of vertex>,<y location>,<link number>" returns the spin (1 or 0) in the 
        // <link number> direction> of the vertex at coordinates <x location> and <y location>.
        // Remember <link number> is from 0 to 3 going counter clockwise starting with 0 to the right
        int get_vertex_link(int vertex_x, int vertex_y, int link_num) const;
        void set_vertex_link(int vertex_x, int vertex_y, int link_num, int set_to);
        // we also want the same thing as above but in the plaquette perspective. In this case each
        // plaquette is numbered with the plaquette at the origin being the (0,0) plaquette and the one to
        // the right being the (1,0) ect. link_num still goes from 0-3 counter clockwise with the
        // 0th link being the right most edge of the plaquette.
        int get_plaquette_link(int plaquette_x, int plaquette_y, int link_num) const;
        void set_plaquette_link(int plaquette_x, int plaquette_y, int link_num, int set_to);

        // For the Z3 visons we need to know the orientation of the links and this is pretty easy to
        // do if you know wich sublattice a particular vertex belongs to
        // Sublattice A starts at the origin
        // returns integer:
        //      0 -> Sublattice A
        //      1 -> Sublattice B
        // Assume arrows point away from sublat A and to sublat B
        int get_vertex_sub_lat(int vertex_x, int vertex_y);
        // Lets make a function that allows the user to specify a specific plaquet and a link from
        // that plaquett. Assuming you are moving in the direction of the link specified this
        // function will return whether you intersect an arrow going to the right of the traveling
        // path or to the left. 
        // This returns an integer value
        //      0 -> arrow to the right
        //      1 -> arrow to the left
        // Assume arrows point away from sublat A and to sublat B
        int intersect_link_direction(int plaquette_x, int plaquette_y, int link_num);
        
        // Also we want to be able to set the links to new values. reference in the same way above but
        // include one more variable which is the value we want to set the link to.

        /* An array ordered from left to right containing the spins of all the links -> left to
         * right for each row -> then up -> then lef to right again.
         * initialize it at 0 (i.e) nothing.
         * There are 2 unique lattice points per vertex so length*height*2 is number of links */
        vector <int> link_vals;
        /* vector containing a vector containing a vertex object */
        vector < vector < Vertex > > v_arr;

        int num_oc_links(int vert_x, int vert_y);

        int find_a_vertex_up_link(int vert_x, int vert_y);
        int find_a_vertex_down_link(int vert_x, int vert_y);

        // array is [number horizontal, number vertical]
        vector <int> count_number_horz_vert_dimers();

        // ****************************************************************************************************************************************
        // ***************** initialize lattice functions ******************************************************************************************
        // ****************************************************************************************************************************************
        void build_init_fully_packed_dimer();
        void build_init_fully_packed_staggard_dimer();
        void build_init_fully_packed_loops();
        void build_other_init_fully_packed();
        void build_init_horz_columns();
        void build_empty_lat();
        void increase_winding_num(int by_how_much);
        void increase_windingX_by_one(int col_index);


    private:

        int find_link(int vert_x, int vert_y, int where);

        /* to build the actual lattice this will take care of making sure all the pointers point to
         * the correct links*/ 
        void weave_links_verts();
};
#endif
