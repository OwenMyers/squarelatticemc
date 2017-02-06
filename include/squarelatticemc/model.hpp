#include "random.hpp"
#include "lattice.hpp"
#include <fstream>

#ifndef model_hpp_
#define model_hpp_

using namespace std;

class Model
{

    public:
        Model(int type_of_model,  O_Random& randomNumberGenerator,string dat_f_dir, Lattice& lat);

        Lattice lat;
        // background lattice
        Lattice br_lat;

        //random number generator
        O_Random rdm;

        int file_precision;
        string dir_to_write;
        
        // Tell the estimator which model you are using
        // 1 ---> dimer only model
        // 2 ---> star dimer model
        int model_type;
        int num_stars;

        int flip(int to_flip);

        // the less efficient but we know is right way. will need this (something like it) for the star model
        int apply_local_update();
        void apply_loop_update();

        // Tell the estimator which model you are using
        // 1 ---> dimer only model
        // 2 ---> star dimer model
        void model(int type_of_model);

        // have a variable that is the number of types of local updates, as we will be adding and
        // subtracting different types in the future. Updates are labels 0-(num_types-1).
        int num_types;

        // Some variables to keep track of the percent success of making moves
        // Vector contains the total number of updates tried for each type
        vector <double> num_updates_by_type;
        // This vector contains the number of successful updates of each type
        vector <double> num_suc_updates_by_type;

        void build_appropriate_backround_lat(Lattice& br_lat);

        bool measure_acceptance_rates;
        vector <double> hist_bins_acceptance_rates;
        ofstream wf_acceptance_rates;
        void write();
        void divide_hist_bins_by(int dev_by);
        void clear_hist_bins();
        void close_wf();


        void get_z3_vison_phase(int start_x, int start_y, int end_x, int end_y);

        // a vector of the moves of a walk that covers the box that a star pair fits into.
        // We check two dimers per vertex because checking one does not end up being an efficient
        // walk (have to recover ground).
        // Start at the upper left vertex
        vector <int> add_horz_star_walk;
        //{0,0,0,0,0,3,2,2,2,2,2,1,3,3,0,0,0,0,0,1}; 

        // Now we want the corresponding vector that is what must be found in order for the star pair
        // to fit. we are checking in the direction of the move that is going to be
        // preformed (because we are walking along the vertices). This is easy because we can just
        // check the link using the integer that is the next move.
        // This will not touch all of the vertices -> this can be dealt with by checking the number
        // of occupied links at a vertex. If it is only one per vertex this will work. If it is not
        // it fails because there must be a star.
        vector <int> add_horz_star_there;
        //{0,1,0,1,0,0,1,0,1,0,1,0,0,0,0,1,0,1,0,0};

        // we also need a vector that is the number of occupied links at each vertex that must be
        // there
        vector <int> add_horz_star_links;
        //{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

        vector <int> annihilate_horz_star_walk;
        vector <int> annihilate_horz_star_there;
        vector <int> annihilate_horz_star_links;

        vector <int> add_vert_star_walk;
        vector <int> add_vert_star_there;
        vector <int> add_vert_star_links;

        vector <int> annihilate_vert_star_walk;
        vector <int> annihilate_vert_star_there;
        vector <int> annihilate_vert_star_links;

        vector <int> move_star_right_walk;
        vector <int> move_star_right_there;
        vector <int> move_star_right_links;

        vector <int> move_star_left_walk;
        vector <int> move_star_left_there;
        vector <int> move_star_left_links;

        vector <int> move_star_up_walk;
        vector <int> move_star_up_there;
        vector <int> move_star_up_links;

        vector <int> move_star_down_walk;
        vector <int> move_star_down_there;
        vector <int> move_star_down_links;

        // this function takes:
        // 1) a starting x point on the lattice
        // 2) a starting y point on the lattice
        // 3) a walk 
        // 4) a "what should be along the walk" vector 
        // 5) a vector of the number of occupied links each vertex must have.
        //
        // returns 1 if "what should be there" is there (and right # of links) and 0 if not
        int can_this_be_done(int start_x, int start_y, vector <int> walk, vector <int> there, vector <int> num_links);


        // if we are using a dimer pentamer model then we want to defined the relative scale of the
        // dimer/pentamer terms in the Hamiltonian. This can be done with a single parameter between
        // [0-1]. 
        // 0 -> only dimers
        // 1 -> only pentamers
        double dim_pent_frac;

        void set_dim_pent_frac(double dim_pent_frac_in);
        



        // Make the simplest pair of stars
        void make_horz_star_pair(int box_corner_x, int box_corner_y);
        void annihilate_horz_star_pair(int box_corner_x, int box_corner_y);
        void make_vert_star_pair(int box_corner_x, int box_corner_y);
        void annihilate_vert_star_pair(int box_corner_x, int box_corner_y);
        void move_star_right(int box_corner_x, int box_corner_y);
        void move_star_left(int box_corner_x, int box_corner_y);
        void move_star_up(int box_corner_x, int box_corner_y);
        void move_star_down(int box_corner_x, int box_corner_y);

        // diagonal move stuff
        int can_ur_be_done(int loc_x, int loc_y);
        int can_dl_be_done(int loc_x, int loc_y);
        int can_ul_be_done(int loc_x, int loc_y);
        int can_dr_be_done(int loc_x, int loc_y);

        void move_ur(int loc_x, int loc_y);
        void move_dl(int loc_x, int loc_y);
        void move_ul(int loc_x, int loc_y);
        void move_dr(int loc_x, int loc_y);



        // (See Ivanov paper)
        // for the correlation function separated in time and space (L/2 spatial) the correlation
        // function  is F(r_ij,t-t')=<V_i(t)V_j(t')>. But by inserting 1 as [V_j(t)]^2 = 1 it can be
        // written as 
        // <V_i(t)V_j(t')> = [V_i(t)V_j(t)][V_j(t)V_j(t')].
        // The first part on the RHS is just the vison operator that we can get from the
        // get_l_div_2_z2_vison() function. The second part is more subtle but amounts to counting
        // the number of times the jth plaquette ( the one at L/2) is flipped) is flipped. every time
        // it is flipped it changes the parity so we multiply by (-1)^(number of times j is fliped).
        // We need a variable to keep track of this number.
        int num_times_l_div_2_vj_flipped;
        int num_times_origin_flipped;
        //int num_times_diag_flipped;
        //int loc_l_div_2_vj_x;
        //int loc_l_div_2_vj_y;
        //int loc_diag_vj_x;
        //int loc_diag_vj_y;



        // ********************************************************************************************************
        // ************* Keeping track of objects (stars/monomers)  ***********************************************
        // ********************************************************************************************************


        // after "diag" we will have either "ur" for "up right.
        //                                  "dl" for "down right (inverse of up right).
        //                                  "ul" for "down right.
        //                                  "dr" for "down right.

        // must set the values of the vectors above in a function
        void set_move_vectors();

        // a variable that if true then we keep track of star positions 
        bool keep_star_locs;
        // a vector of vectors. the second vector is just the positions of a star (long rectangle
        // vector) Number of stars by 2(for xy positions)
        // star_locs[i][0] would be the star "i"'s x position 
        // star_locs[i][1] would be the star "i"'s y position 
        vector <vector <int> > star_locs;
        // these function will be used to simplify the readability of implementing the "keeping track
        // of stars"
        // add a star to the star_locs. It is at location x
        void add_star(int x,int y);
        void rm_star(int x,int y);
        // add_x(y) is the amount the star has moved by. Can be positive or negative
        void mv_star(int x, int y, int add_x, int add_y);

        // file for star locations "star_loc_file" for Write Star File.
        ofstream star_loc_file;
        void open_star_loc_file();
        void close_star_loc_file();



        // function that sets keep_star_locs to true (1)
        void keep_track_of_star_locations();
        // function that writes the star locations to a file
        void write_star_locations();


        // ********* Monomers *********
        // are we including monomers?
        bool include_monomers;

        // same stuff as for keeping track of the star locations but for the monomers.
        bool keep_monomer_locs;
        vector <vector <int> > monomer_locs;
        void add_monomer(int x, int y);
        void rm_monomer(int x, int y);
        // add_x(y) is the amount the monomer has moved by. Can be positive or negative
        void mv_monomer(int x, int y, int add_x, int add_y);
 
        // file for monomer locations "monomer_loc_file" for Write Star File.
        ofstream monomer_loc_file;
        void open_monomer_loc_file();
        void close_monomer_loc_file();

        // function that sets keep_monomer_locs to true (1)
        void keep_track_of_monomer_locations();
        // function that writes the monomer locations to a file
        void write_monomer_locations();
        
        // function sees if there is a monomer on the plaquette of interest (upper right plaquette of
        // loc_x, loc_y vertex, or just the loc_x, loc_y plaquette) and then see if the plaquette is
        // "monomer flippable
        // A monomer flip is different because the monomer must remain on the same sublattice so it
        // goes diagonally across the plaquette and the dimer (if not part of a pentamer) rotates
        // around to occupy the blank space
        bool is_there_a_monomer_can_we_flip(int loc_x,int loc_y);
        // this does the actual moving
        void monomer_flip_plaquette(int loc_x,int loc_y);

                // stuff for measuring loop lengths
        int system_find_up(const vector <int> & link_vals);
        vector <int> which_vertex_from_link(int cur_link,int w,int l);
        int find_loop_length(int cur_go_to_link, int start_x, int start_y);

        // which of these is used is determined by model type and chosen correctly in the public
        // function apply_____
        void dimer_dimer_local_update();
        void apply_full_pack_walk();
        int star_dimer_local_update();

        
        // These functions will allow us to initialize in different modifications of the lattice
        
        // make a horizontal crystal lattice
        void make_horz_crystal_lattice();

        // ********************************************************************************************************
        // ************* Toric Code stuff *************************************************************************
        // ********************************************************************************************************
        // This is just a random walk. NOT self avoiding ONLY for toric code
        void apply_toric_code_rand_walk_loop();
        void toric_code_local_update();


};
#endif
