#include "model.hpp"

Model::Model(int type_of_model,  O_Random& randomNumberGenerator,string dat_f_dir, Lattice& lat):
    rdm(randomNumberGenerator),model_type(type_of_model),lat(lat)
{
    // Hard coded because this should only change if we change the way we are implementing the star
    // dimer local updates
    num_types = 13;
    num_stars = 0;
    measure_acceptance_rates = true;
    dir_to_write = dat_f_dir;

    num_updates_by_type.resize(num_types+1);
    num_suc_updates_by_type.resize(num_types+1);

    keep_star_locs = true;

    if ((model_type == 2)||(model_type==3))
    {
        set_move_vectors();
    }

    file_precision = 15;

    hist_bins_acceptance_rates.resize(num_types);
    wf_acceptance_rates.precision(file_precision);
    wf_acceptance_rates.setf( std::ios::fixed, std:: ios::floatfield );
    wf_acceptance_rates.open((dir_to_write + "acceptance_rates.txt").c_str());
}

void Model::set_dim_pent_frac(double dim_pent_frac_in)
{
    dim_pent_frac = dim_pent_frac_in;
}

void Model::build_appropriate_backround_lat(Lattice& br_lat)
{

    if ((model_type == 1)||(model_type == 2)||(model_type == 3))
    {
        br_lat.build_init_fully_packed_dimer();
        //br_lat.build_init_horz_columns();
    }
    else if (model_type == 4)
    {
        //br_lat.build_init_fully_packed_loops();
        br_lat.build_empty_lat();
    }

}

int Model::apply_local_update()
{
    int to_return = 0;

    if (model_type == 1)
    {
        dimer_dimer_local_update();
    }
    else if ((model_type == 2)||(model_type == 3))
    {
        to_return = star_dimer_local_update();
    }
    else if (model_type == 4)
    {
        toric_code_local_update();
    }
    return to_return;

}

void Model::dimer_dimer_local_update()
{
    // Choose a random vertex to orient to.
    int loc_x = rdm.get_rand_int(0,lat.length);
    // random y (height) position)
    int loc_y = rdm.get_rand_int(0,lat.height);

    // The plaquette of interest will be the upper right plaquette from the vertex selected
    //
    // We need to see if it is flippable. In both the fully packed loop and dimer models it is
    // flippable iff parallel edges of the plaquette have the same spin at the link (equivalently both
    // yes or no dimers) and also one pair of parallel edges must be different than the other pair
    // of parallel edges
    //
    // starting with the leftmost edge we will number them clockwise. This could be a vector if
    // there were many edges but there are only 4 so the lazy way is not much different.
    // These are going to be pointers to the actual spins so if we can and will flip we can just
    // work with these short names

    // spin of left edge
    int  e1 = lat.get_plaquette_link(loc_x,loc_y,2);
    // spin of top edge. upper right diagonal vertex -> must be careful with periodic BC
    int  e2 = lat.get_plaquette_link(loc_x,loc_y,1);
    // spin of right edge
    int  e3 = lat.get_plaquette_link(loc_x,loc_y,0);
    // spin of botom edge
    int  e4 = lat.get_plaquette_link(loc_x,loc_y,3);

    // check if flipable
    if ((e1==e3)&(e2==e4)&(e1!=e2))
    {
        lat.set_plaquette_link(loc_x,loc_y,2,flip(e1));
        lat.set_plaquette_link(loc_x,loc_y,1,flip(e2));
        lat.set_plaquette_link(loc_x,loc_y,0,flip(e3));
        lat.set_plaquette_link(loc_x,loc_y,3,flip(e4));

        // see if we are flipping the plaquette at the second vison (Vj). Always flipped in toric code
        //if ( measure_l_div_2_z2_vison_time_cor ){
        //    if ((loc_x==loc_l_div_2_vj_x)&&(loc_y==loc_l_div_2_vj_y))
        //    {
        //        num_times_l_div_2_vj_flipped++;
        //    }
        //}
        //if ( measure_origin_z2_vison_time_cor ){
        //    if ((loc_x==0)&&(loc_y==0))
        //    {
        //        num_times_origin_flipped++;
        //    }
        //}
        //if (measure_diag_z2_vison_time_cor ){
        //    if ((loc_x==(length/2))&&(loc_y==(height/2)))
        //    {
        //        num_times_diag_flipped++;
        //    }
        //}
    }
    // done
}
int Model::star_dimer_local_update()
{

    // if we are using a dimer pentamer model then we want to defend the relative scale of the
    // dimer/pentamer terms in the Hamiltonian. This can be done with a single parameter between
    // [0-1]. 
    // 0 -> only dimers
    // 1 -> only pentamers
    // This is the random number between 0 and 1
    double rand_num = rdm.random01();
    //std::cout << "random number: " << rand_num << std::endl;
    //std::cout << "num_types (should be 13): " << num_types << std::endl;
    //std::cout << "dim_pent_frac: " << dim_pent_frac << std::endl;

    // return if the update was a success
    int did_it_happen = 0;

    // integer that is true false of if a move will be implemented (1 =true) (0=false)
    // move yes or no -> mv_yn
    int mv_yn;
    
    // what type of update are we going to do?
    // randomly choose a number from 0 - the following
    // 0)  create horizontal star pair
    // 1)  create vertical star pair
    // 2)  annihilate horizontal star pair
    // 3)  annihilate vertical star pair
    // 4)  move star right
    // 5)  move star left
    // 6)  move star up
    // 7)  move star down
    // 8)  move star diag up right
    // 9) move star diag down left
    // 10) move star diag up left
    // 11) move star diag down right
    // 12)  monomer stuff
    // 13) flip plaquette
    int update_type = -1;

    // Depending on whether we are going to do a dimer or a pentamer type update we will choose
    // different ranges of the move types above
    if (rand_num <= dim_pent_frac)
    {
        // if we are in here do a dimer move
        update_type = 13;
        //cout << "num_types" << num_types << "\n";
        //cout << "update type " << update_type << "\n";
    }
    else if (rand_num > dim_pent_frac)
    {
        // if we are in here do a pentamer move
        // 0-9 gives us update in interval [0,9)
        // This does not include 12 so we will not accidentally get a dimer move if we chose a
        // pentamer move
        update_type = rdm.get_rand_int(0,num_types-1);
    }
//    //std::cout << "update_type: " << update_type << std::endl;
//    if (num_types==14)
//    {
//        // STOP: NO Monomer moves right now.  Need to set up additional fraction
//        cout << "NO Monomer moves right now!"<<endl;
//        exit(EXIT_FAILURE);
//    }
//    if (update_type==-1)
//    {
//        // STOP: NO Monomer moves right now.  Need to set up additional fraction
//        cout << "No update_type was defined!"<<endl;
//        exit(EXIT_FAILURE);
//    }
    

    // Choose a random vertex to orient to. always upper left of checking box
    int loc_x = rdm.get_rand_int(0,lat.length);
    // random y (height) position)
    int loc_y = rdm.get_rand_int(0,lat.height);

    vector <int> walk;
    vector <int> there;
    vector <int> num_links;

    if (update_type == 0)
    {
        walk = add_horz_star_walk; 
        there = add_horz_star_there;
        num_links =  add_horz_star_links;
    }
    if (update_type == 1)
    {
        walk = add_vert_star_walk; 
        there = add_vert_star_there;
        num_links =  add_vert_star_links;
    }
    if (update_type == 2)
    {
        walk = annihilate_horz_star_walk; 
        there = annihilate_horz_star_there;
        num_links =  annihilate_horz_star_links;
    }
    if (update_type == 3)
    {
        walk = annihilate_vert_star_walk; 
        there = annihilate_vert_star_there;
        num_links =  annihilate_vert_star_links;
    }
    if (update_type == 4)
    {
        walk = move_star_right_walk; 
        there = move_star_right_there;
        num_links =  move_star_right_links;
    }
    if (update_type == 5)
    {
        walk = move_star_left_walk; 
        there = move_star_left_there;
        num_links =  move_star_left_links;
    }
    if (update_type == 6)
    {
        walk = move_star_up_walk; 
        there = move_star_up_there;
        num_links =  move_star_up_links;
    }
    if (update_type == 7)
    {
        walk = move_star_down_walk; 
        there = move_star_down_there;
        num_links =  move_star_down_links;
    }

    if (update_type < 8)
    {
        mv_yn = can_this_be_done(loc_x, loc_y, walk, there, num_links);
        
        //cout<<"TESTING -> See what we get from can this be done \n can it be done = " << mv_yn << "\n";

        if (mv_yn) 
        {
            if (update_type == 0) 
            {   
                make_horz_star_pair(loc_x, loc_y); 
                num_stars += 2; 
                if (keep_star_locs) 
                { 
                    add_star(loc_x,loc_y-1); add_star(loc_x+1,loc_y-1);
                } 
            }
            if (update_type == 1) { make_vert_star_pair(loc_x, loc_y); num_stars += 2; if (keep_star_locs) { add_star(loc_x+1,loc_y);  add_star(loc_x+1,loc_y-1);}}
            if (update_type == 2) { annihilate_horz_star_pair(loc_x, loc_y); num_stars -= 2; if (keep_star_locs) { rm_star(loc_x,loc_y-1);  rm_star(loc_x+1,loc_y-1);}}
            if (update_type == 3) { annihilate_vert_star_pair(loc_x, loc_y); num_stars -= 2; if (keep_star_locs) { rm_star(loc_x+1,loc_y);  rm_star(loc_x+1,loc_y-1);}}
            if (update_type == 4) { move_star_right(loc_x, loc_y); if (keep_star_locs) { mv_star(loc_x,loc_y-1,2,0); }}
            if (update_type == 5) { move_star_left(loc_x, loc_y); if (keep_star_locs) { mv_star(loc_x+2,loc_y-1,-2,0); }}
            if (update_type == 6) { move_star_up(loc_x, loc_y); if (keep_star_locs) { mv_star(loc_x+1,loc_y-2,0,2); }}
            if (update_type == 7) { move_star_down(loc_x, loc_y); if (keep_star_locs) { mv_star(loc_x+1,loc_y,0,-2); }}
            
            did_it_happen = 4;

        }
    }
    // this is for the diagonal moves which can be checked very easily
    // 9)  move star diag up right
    // 10) move star diag down left
    // 11) move star diag up left
    // 12) move star diag down right
    // INCLUDE THE MOVE ITSELF HERE
    if (update_type == 8)  
    { 
        if (can_ur_be_done(loc_x, loc_y)) 
        {
            move_ur(loc_x, loc_y);
            did_it_happen = 4;
            if (keep_star_locs) { mv_star(loc_x+1,loc_y-2,1,1); }
        }
    }
    if (update_type == 9) 
    { 
        if (can_dl_be_done(loc_x, loc_y)) 
        {
            move_dl(loc_x, loc_y);
            did_it_happen = 4;
            if (keep_star_locs) { mv_star(loc_x+2,loc_y-1,-1,-1); }
        }
    }
    if (update_type == 10) 
    { 
        if (can_ul_be_done(loc_x, loc_y)) 
        {
            move_ul(loc_x, loc_y);
            did_it_happen = 4;
            if (keep_star_locs) { mv_star(loc_x+2,loc_y-2,-1,1); }
        }
    }
    if (update_type == 11) 
    { 
        if (can_dr_be_done(loc_x, loc_y)) 
        {
            move_dr(loc_x, loc_y);
            did_it_happen = 4;
            if (keep_star_locs) { mv_star(loc_x+1,loc_y-1,1,-1); }
        }
    }
    //if (update_type == 11) { if (can_ul_be_done(loc_x, loc_y)) {move_ul(loc_x, loc_y);did_it_happen = 1;}}
    //if (update_type == 12) { if (can_dr_be_done(loc_x, loc_y)) {move_dr(loc_x, loc_y);did_it_happen = 1;}}
    //
    if (update_type == 13)
    {
        // spin of left edge
        //int  e1 = lattice[loc_y][loc_x].links[1];
        //int  e1 = lat.get_vertex_link(loc_x,loc_y,1);
        int  e1 = lat.get_plaquette_link(loc_x,loc_y,2);
        // spin of top edge. upper right diagonal vertex -> must be careful with periodic BC
        //int  e2 = lattice[(loc_y + 1)%height][(loc_x + 1)%length].links[2];
        //int  e2 = lat.get_vertex_link(loc_x + 1,loc_y + 1,2);
        int  e2 = lat.get_plaquette_link(loc_x,loc_y,1);
        // spin of right edge
        //int  e3 = lattice[(loc_y + 1)%height][(loc_x + 1)%length].links[3];
        //int  e3 = lat.get_vertex_link(loc_x + 1,loc_y + 1,3);
        int  e3 = lat.get_plaquette_link(loc_x,loc_y,0);
        // spin of bottom edge
        //int  e4 = lat.get_vertex_link(loc_x,loc_y,0);
        int  e4 = lat.get_plaquette_link(loc_x,loc_y,3);

        // check if flippable
        if ((e1==e3)&(e2==e4)&(e1!=e2))
        {
            did_it_happen = 2;
            lat.set_plaquette_link(loc_x,loc_y,2,flip(e1));
            lat.set_plaquette_link(loc_x,loc_y,1,flip(e2));
            lat.set_plaquette_link(loc_x,loc_y,0,flip(e3));
            lat.set_plaquette_link(loc_x,loc_y,3,flip(e4));

            // see if we are flipping the plaquette at the second vison (Vj). Always flipped in toric code
            //if ( measure_l_div_2_z2_vison_time_cor ){
            //    if ((loc_x==loc_l_div_2_vj_x)&&(loc_y==loc_l_div_2_vj_y))
            //    {
            //        num_times_l_div_2_vj_flipped++;
            //    }
            //}
            //if ( measure_origin_z2_vison_time_cor ){
            //    if ((loc_x==0)&&(loc_y==0))
            //    {
            //        num_times_origin_flipped++;
            //    }
            //}
            //if (measure_diag_z2_vison_time_cor ){
            //    if ((loc_x==(length/2))&&(loc_y==(height/2)))
            //    {
            //        num_times_diag_flipped++;
            //    }
            //}
        }
        // done

    }
    
    // monomer moves
    if (update_type == 12)
    {
        // this returns true if it can be done IT ALSO DOES IT
        did_it_happen = is_there_a_monomer_can_we_flip(loc_x, loc_y);
        //cout << "did_it_happen" << did_it_happen << endl;
        //just for now so we can only write configs where a monomer changed location
        if (did_it_happen) {did_it_happen=3;}
    }


    
    if (measure_acceptance_rates)
    {
        //cout << "did_in_happen " << did_it_happen << "\n";
        //cout << "update_type: " << update_type << endl;
        //cout << "num_updates_by_type[update_type] before: " << num_updates_by_type[update_type] << endl;
        num_updates_by_type[update_type] += 1.0;
        //cout << "num_updates_by_type[update_type] after: " << num_updates_by_type[update_type] << endl;
        //cout << "num_suc_updates_by_type[update_type] before: " << num_suc_updates_by_type[update_type] << endl;
        if (did_it_happen){ 
            num_suc_updates_by_type[update_type] += 1.0;
        }
        //cout << "num_suc_updates_by_type[update_type] after: " << num_suc_updates_by_type[update_type] << endl;

        // get the fractions
        for (int i = 0; i<num_types; i++)
        {
            if (num_updates_by_type[i]>(double)0)
            { 
                hist_bins_acceptance_rates[i] += (double) num_suc_updates_by_type[i]/( (double)num_updates_by_type[i]);
                //cout << "update type " << i << " fracion of suc " << hist_bins_acceptance_rates[i] << "\n";
            }
        }
    }

//    if (measure_origin_z3_vison_time_cor)
//    {
//        if (has_z3_phase_changed(0,0, origin_cur_plq))
//        {
//            z3_origin_num_inc += how_has_z3_phase_changed(0,0,origin_cur_plq,update_type);
//            //cout << "the origin plaquette has changed. originally it was: " << endl;
//            //for (int i=0; i < 4; i++)
//            //{
//            //    cout << origin_cur_plq[i] << " ";
//            //}
//            //cout << endl;
//            //cout << "reseting it... " << endl;
//            set_origin_cur_plq();
//            //cout << "Now its: " << endl;
//            //for (int i=0; i < 4; i++)
//            //{
//            //    cout << origin_cur_plq[i] << " ";
//            //}
//            //cout << endl;
//        }
//    }
//    if (measure_l_div_2_z3_vison_time_cor)
//    {
//        if (has_z3_phase_changed(0,height/2, l_div_2_cur_plq))
//        {
//            z3_l_div_2_num_inc += how_has_z3_phase_changed(0,height/2,l_div_2_cur_plq,update_type);
//            set_l_div_2_cur_plq();
//        }
//    }
    //if (measure_diag_z3_vison_time_cor)
    //{
    //    if (has_z3_phase_changed(length/2,height/2, diag_cur_plq))
    //    {
    //        z3_diag_num_inc += how_has_z3_phase_changed(length/2,height/2,diag_cur_plq,update_type);
    //        set_diag_cur_plq();
    //    }
    //}
   

    //if (did_it_happen) 
    //{
    //    cout << "current number of stars " << num_stars << "\n";
    //    std::cout << "update_type = " << update_type  << std::endl;
    //}

    return did_it_happen;

}
void Model::toric_code_local_update()
{
    int length = lat.length;
    int height = lat.height;

    // Choose a random vertex to orient to.
    int loc_x = rdm.get_rand_int(0,length);
    // random y (height) position)
    int loc_y = rdm.get_rand_int(0,height);

    //cout << "loc_x: " << loc_x << "\n";
    //cout << "loc_y: " << loc_y << "\n\n";

    // see if we are flipping the plaquette at the second vison (Vj). Always flipped in toric code
    // so always add one 
    //if (measure_l_div_2_z2_vison_time_cor ){
    //    if ((loc_x==loc_l_div_2_vj_x)&&(loc_y==loc_l_div_2_vj_y))
    //    {
    //        num_times_l_div_2_vj_flipped++;
    //    }
    //}
    //if (measure_origin_z2_vison_time_cor ){
    //    if ((loc_x==0)&&(loc_y==0))
    //    {
    //        num_times_origin_flipped++;
    //        //cout << "num_times_origin_flipped: " << num_times_origin_flipped << "\n\n";
    //    }
    //}
    //if (measure_diag_z2_vison_time_cor ){
    //    if ((loc_x==(length/2))&&(loc_y==(height/2)))
    //    {
    //        num_times_diag_flipped++;
    //    }
    //}

    // In the toric code you just flip the degrees of fredom on the plaquett no mater what
    // The plaquette of interest will be the uper right plaquette from the vertex selected
    //
    // spin of left edge
    int  e1 = lat.get_plaquette_link(loc_x,loc_y,2);
    // spin of top edge. upper right diagonal vertex -> must be carefull with periodic BC
    int  e2 = lat.get_plaquette_link(loc_x,loc_y,1);
    // spin of right edge
    int  e3 = lat.get_plaquette_link(loc_x,loc_y,0);
    // spin of botom edge
    int  e4 = lat.get_plaquette_link(loc_x,loc_y,3);

    lat.set_plaquette_link(loc_x,loc_y,2,flip(e1));
    lat.set_plaquette_link(loc_x,loc_y,1,flip(e2));
    lat.set_plaquette_link(loc_x,loc_y,0,flip(e3));
    lat.set_plaquette_link(loc_x,loc_y,3,flip(e4));


}
void Model::write()
{
    for (int i = 0; i < hist_bins_acceptance_rates.size(); i++)
    {
        wf_acceptance_rates << hist_bins_acceptance_rates[i] << " ";
    }
    wf_acceptance_rates << endl; 
}
void Model::divide_hist_bins_by(int dev_by)
{
    for (int i = 0; i<hist_bins_acceptance_rates.size(); i++)
    {
        hist_bins_acceptance_rates[i] = hist_bins_acceptance_rates[i]/ ((double) dev_by);
        //std::cout << "hist_bins_acceptance_rates[i]: "<< hist_bins_acceptance_rates[i] << std::endl;
    }
}
void Model::clear_hist_bins()
{
    for (int i = 0; i < hist_bins_acceptance_rates.size(); i++)
    {
        hist_bins_acceptance_rates[i] = (double)0;
    }
}
void Model::close_wf()
{
    wf_acceptance_rates.close();
}




//******************************************************************************************************************************************************************************************
//******************************************************************************************************************************************************************************************
//******************************************************************************************************************************************************************************************
//******************************************************************************************************************************************************************************************
int Model::can_this_be_done(int start_x, int start_y, vector <int> walk, vector <int> there, vector <int> num_links)
{
    int length = lat.length;
    int height = lat.height;
    // This is true until we find something that does not work. Then it is set to false and the loop
    // breaks.

    int can_it_be_done = 1;

    int cur_x = start_x;
    int cur_y = start_y;

    // Some self explanatory variable that will make the program more readable
    int cur_actual_num_links;
    int cur_should_be_num_links;

    int cur_direction;
    
    int cur_actual_link;
    int cur_should_be_link;

    //cout << "cnfg check walk.size() " << walk.size() << "\n";

    for (int i = 0; i<= walk.size(); i ++)
    {
        // this if statement means that we can bing i to one larger than the walk.size() in order to
        // check the number of links at the final vertex.
        if (i<walk.size())
        {
            cur_direction = walk[i];

            // check that the link along the direction we are going to walk in is correct
            //cur_actual_link = *(lattice[cur_y][cur_x].links[cur_direction]);
            cur_actual_link = lat.get_vertex_link(cur_x,cur_y,cur_direction);
            cur_should_be_link = there[i];
            if(cur_actual_link != cur_should_be_link) { can_it_be_done = 0; break; }
        }


        // make sure the number of links at this vertex is correct
        // This should be a vertex object so it calls a function to count local occupied links
        cur_actual_num_links = lat.num_oc_links(cur_x,cur_y);
        cur_should_be_num_links = num_links[i];
        if(cur_actual_num_links != cur_should_be_num_links) { can_it_be_done = 0; break; }


        
        // move the current location of the vertex to the new vertex
        if (cur_direction == 0){cur_x = (cur_x + 1)%length;}
        if (cur_direction == 1){cur_y = (cur_y + 1)%height;}
        // The "+l )% l" here takes care of the boundary conditions when the location goes negative 
        if (cur_direction == 2){cur_x = (cur_x - 1 +length)%(length);}
        if (cur_direction == 3){cur_y = (cur_y - 1 +height)%(height);}

    }
    return can_it_be_done;
}
int Model::can_ur_be_done(int loc_x, int loc_y) 
{
    int first_important_dimer;

    // first see if there is a star in the right place
    if (lat.num_oc_links(loc_x+1,loc_y-2)==4)
    {
        //cout << "found a star in the corner \n";
        // now make sure there is NOT one in the far opposite corner
        if (lat.num_oc_links(loc_x+3,loc_y)==1)
        {
            //cout << "making sure there is not one in the opposite corner \n";
            // the function "find_a_vertex_up_link" returns a random occupied link attached the the
            // vertex in question. at the vertex directly diagonal to the star in the direction of
            // potential movement there is only one occupied link so "find_a_vertex_up_link" will tell
            // us where it is.
            first_important_dimer = lat.find_a_vertex_up_link(loc_x+2,loc_y-1);

            // Only 2 possibilities
            // OR IT COULD BE A MONOMER -> in which case the function lat.find_a_vertex_up_link returns -1 
            // cant do if we have a monomer
            if (first_important_dimer == -1)
            {
                return 0;
            }
            else if (first_important_dimer == 0)
            {
                // make sure the other two are in place
                if ( (lat.get_vertex_link(loc_x+1,loc_y,0)==1)&&(lat.get_vertex_link(loc_x,loc_y,3)==1))
                {
                    return 1;
                }
            }
            else if (first_important_dimer == 1)
            {
                // make sure there is one in one of the two required locations 
                if ( (lat.get_vertex_link(loc_x+3,loc_y-1,3)==1)&&(lat.get_vertex_link(loc_x+2,loc_y-3,0)==1))
                {
                    return 1;
                }
            }

        }

    }
    return 0;
}
int Model::can_ul_be_done(int loc_x, int loc_y) 
{
    int first_important_dimer;

    // first see if there is a star in the right place
    if (lat.num_oc_links(loc_x+2,loc_y-2)==4)
    {
        // now make sure there is NOT one in the far opposite corner
        if (lat.num_oc_links(loc_x,loc_y)==1)
        {
            // the function "find_a_vertex_up_link" returns a random occupied link attached the the
            // vertex in question. at the vertex directly diagonal to the star in the direction of
            // potential movement there is only one occupied link so "find_a_vertex_up_link" will tell
            // us where it is.
            first_important_dimer = lat.find_a_vertex_up_link(loc_x+1,loc_y-1);

            // Only 2 possibilities
            // OR IT COULD BE A MONOMER -> in which case the function lat.find_a_vertex_up_link returns -1 
            // cant do if we have a monomer
            if (first_important_dimer == -1)
            {
                return 0;
            }
            else if (first_important_dimer == 2)
            {
                // make sure there is one in the two required locations 
                if ( (lat.get_vertex_link(loc_x+1,loc_y,0)==1)&&(lat.get_vertex_link(loc_x+3,loc_y,3)==1))
                {
                    return 1;
                }
            }
            else if (first_important_dimer == 1)
            {
                // make sure there is one in one of the two required locations 
                if ( (lat.get_vertex_link(loc_x,loc_y-1,3)==1)&&(lat.get_vertex_link(loc_x,loc_y-3,0)==1))
                {
                    return 1;
                }
            }

        }

    }
    return 0;
}
int Model::can_dr_be_done(int loc_x, int loc_y) 
{
    int first_important_dimer;

    // first see if there is a star in the right place
    if (lat.num_oc_links(loc_x+1,loc_y-1)==4)
    {
        // now make sure there is NOT one in the far opposite corner
        if (lat.num_oc_links(loc_x+3,loc_y-3)==1)
        {
            // the function "find_a_vertex_up_link" returns a random occupied link attached the the
            // vertex in question. at the vertex directly diagonal to the star in the direction of
            // potential movement there is only one occupied link so "find_a_vertex_up_link" will tell
            // us where it is.
            first_important_dimer = lat.find_a_vertex_up_link(loc_x+2,loc_y-2);

            // Only 2 possibilities
            // OR IT COULD BE A MONOMER -> in which case the function lat.find_a_vertex_up_link returns -1 
            // cant do if we have a monomer
            if (first_important_dimer == -1)
            {
                return 0;
            }
            else if (first_important_dimer == 0)
            {
                // make sure there is one in the two required locations 
                if ( (lat.get_vertex_link(loc_x,loc_y-2,3)==1)&&(lat.get_vertex_link(loc_x+1,loc_y-3,0)==1))
                {
                    return 1;
                }
            }
            else if (first_important_dimer == 3)
            {
                // make sure there is one in one of the two required locations 
                if ( (lat.get_vertex_link(loc_x+2,loc_y,0)==1)&&(lat.get_vertex_link(loc_x+3,loc_y-1,3)==1))
                {
                    return 1;
                }
            }

        }

    }
    return 0;
}
int Model::can_dl_be_done(int loc_x, int loc_y) 
{
    int first_important_dimer;

    // first see if there is a star in the right place
    //if (lattice[(loc_y-1+height)%height][(loc_x+2)%length].num_oc_links()==4)
    if (lat.num_oc_links(loc_x+2,loc_y-1)==4)
    {
        //cout << "found a star in the corner \n";
        // now make sure there is NOT one in the far opposite corner
        //if (lattice[(loc_y-3+height)%height][loc_x].num_oc_links()==1)
        if (lat.num_oc_links(loc_x,loc_y-3)==1)
        {
            //cout << "making sure there is not one in the opposite corner \n";
            // the function "find_a_vertex_up_link" returns a random occupied link attached the the
            // vertex in question. at the vertex directly diagonal to the star in the direction of
            // potential movement there is only one occupied link so "find_a_vertex_up_link" will tell
            // us where it is.
            first_important_dimer = lat.find_a_vertex_up_link(loc_x+1,loc_y-2);

            // Only 2 possibilities
            // OR IT COULD BE A MONOMER -> in which case the function lat.find_a_vertex_up_link returns -1 
            // cant do if we have a monomer
            if (first_important_dimer == -1)
            {
                return 0;
            }
            if (first_important_dimer == 2)
            {
                // make sure there is one in one of the two required locations 
                if (( lat.get_vertex_link(loc_x+1,loc_y-3,0)==1)&&((lat.get_vertex_link(loc_x+3,loc_y-3,1)==1)))
                {
                    return 1;
                }
            }
            else if (first_important_dimer == 3)
            {
                // make sure there is one in one of the two required locations 
                if ( (lat.get_vertex_link(loc_x,loc_y-1,3)==1)&&(lat.get_vertex_link(loc_x,loc_y,0)==1))
                {
                    return 1;
                }
            }

        }

    }
    return 0;
}
void Model::move_star_up(int box_corner_x, int box_corner_y)
{
    // This is the old box corner but we can use this (because we know it is correct) by shifting box
    // corner up along the y axis by 1
    box_corner_y += 1;

    lat.set_vertex_link(box_corner_x,box_corner_y-1,0, 1);
    lat.set_vertex_link(box_corner_x,box_corner_y-1,3, 0);
    lat.set_vertex_link(box_corner_x+1,box_corner_y-1,0, 1);
    lat.set_vertex_link(box_corner_x+1,box_corner_y-1,3, 1);
    lat.set_vertex_link(box_corner_x+2,box_corner_y-1,3, 0);
    lat.set_vertex_link(box_corner_x,box_corner_y-2,3, 1);
    lat.set_vertex_link(box_corner_x+1,box_corner_y-2,3, 0);
    lat.set_vertex_link(box_corner_x+2,box_corner_y-2,3, 1);
    lat.set_vertex_link(box_corner_x,box_corner_y-3,0, 0);
    lat.set_vertex_link(box_corner_x+1,box_corner_y-3,0, 0);
}
void Model::move_star_down(int box_corner_x, int box_corner_y)
{
    // This is the old box corner but we can use this (because we know it is correct) by shifting box
    // corner up along the y axis by 1
    box_corner_y += 1;
    lat.set_vertex_link(box_corner_x,box_corner_y-1,0, 0);
    lat.set_vertex_link(box_corner_x,box_corner_y-1,3, 1);
    lat.set_vertex_link(box_corner_x+1,box_corner_y-1,0, 0);
    lat.set_vertex_link(box_corner_x+1,box_corner_y-1,3, 0);
    lat.set_vertex_link(box_corner_x+2,box_corner_y-1,3, 1);
    lat.set_vertex_link(box_corner_x,box_corner_y-2,3, 0);
    lat.set_vertex_link(box_corner_x+1,box_corner_y-2,3, 1);
    lat.set_vertex_link(box_corner_x+2,box_corner_y-2,3, 0);
    lat.set_vertex_link(box_corner_x,box_corner_y-3,0, 1);
    lat.set_vertex_link(box_corner_x+1,box_corner_y-3,0, 1);
}
void Model::move_star_left(int box_corner_x, int box_corner_y)
{

    // This is the old box corner but we can use this (because we know it is correct) by shifting box
    // corner back along the x axis by 1
    box_corner_x -= 1;

    lat.set_vertex_link(box_corner_x+1,box_corner_y,0, 0);
    lat.set_vertex_link(box_corner_x+1,box_corner_y,3, 1);
    lat.set_vertex_link(box_corner_x+2,box_corner_y,0, 1);
    lat.set_vertex_link(box_corner_x+3,box_corner_y,3, 0);
    lat.set_vertex_link(box_corner_x+1,box_corner_y-1,0, 1);
    lat.set_vertex_link(box_corner_x+1,box_corner_y-1,3, 1);
    lat.set_vertex_link(box_corner_x+2,box_corner_y-1,0, 0);
    lat.set_vertex_link(box_corner_x+3,box_corner_y-1,3, 0);
    lat.set_vertex_link(box_corner_x+1,box_corner_y-2,0, 0);
    lat.set_vertex_link(box_corner_x+2,box_corner_y-2,0, 1);
}
void Model::move_star_right(int box_corner_x, int box_corner_y)
{
    // This is the old box corner but we can use this (because we know it is correct) by shifting box
    // corner back along the x axis by 1
    box_corner_x -= 1;

    lat.set_vertex_link(box_corner_x+1,box_corner_y,0, 1);
    lat.set_vertex_link(box_corner_x+1,box_corner_y,3, 0);
    lat.set_vertex_link(box_corner_x+2,box_corner_y,0, 0);
    lat.set_vertex_link(box_corner_x+3,box_corner_y,3, 1);
    lat.set_vertex_link(box_corner_x+1,box_corner_y-1,0, 0);
    lat.set_vertex_link(box_corner_x+1,box_corner_y-1,3, 0);
    lat.set_vertex_link(box_corner_x+2,box_corner_y-1,0, 1);
    lat.set_vertex_link(box_corner_x+3,box_corner_y-1,3, 1);
    lat.set_vertex_link(box_corner_x+1,box_corner_y-2,0, 1);
    lat.set_vertex_link(box_corner_x+2,box_corner_y-2,0, 0);
}
// Pass it the location of the upper left corner of the relevant box
void Model::make_horz_star_pair(int box_corner_x, int box_corner_y)
{
    lat.set_vertex_link(box_corner_x,box_corner_y,3, 1);
    lat.set_vertex_link(box_corner_x,box_corner_y,0, 0);
    lat.set_vertex_link(box_corner_x+1,box_corner_y,3, 1);
    lat.set_vertex_link(box_corner_x+1,box_corner_y-1,3, 1);
    lat.set_vertex_link(box_corner_x,box_corner_y-1,3, 1);
    lat.set_vertex_link(box_corner_x,box_corner_y-1,0, 1);
    lat.set_vertex_link(box_corner_x,box_corner_y-2,0, 0);
}
void Model::annihilate_horz_star_pair(int box_corner_x, int box_corner_y)
{
    lat.set_vertex_link(box_corner_x,box_corner_y,3, 0);
    lat.set_vertex_link(box_corner_x,box_corner_y,0, 1);
    lat.set_vertex_link(box_corner_x+1,box_corner_y,3, 0);
    lat.set_vertex_link(box_corner_x+1,box_corner_y-1,3, 0);
    lat.set_vertex_link(box_corner_x,box_corner_y-1,3, 0);
    lat.set_vertex_link(box_corner_x,box_corner_y-1,0, 0);
    lat.set_vertex_link(box_corner_x,box_corner_y-2,0, 1);
}

void Model::make_vert_star_pair(int box_corner_x, int box_corner_y)
{
    lat.set_vertex_link(box_corner_x,box_corner_y,3, 0);
    lat.set_vertex_link(box_corner_x,box_corner_y,0, 1);

    lat.set_vertex_link(box_corner_x,box_corner_y-1,0, 1);
   
    lat.set_vertex_link(box_corner_x+1,box_corner_y-1,0, 1);

    lat.set_vertex_link(box_corner_x+1,box_corner_y,0, 1);
    lat.set_vertex_link(box_corner_x+1,box_corner_y,3, 1);

    lat.set_vertex_link(box_corner_x+2,box_corner_y,3, 0);
}

void Model::annihilate_vert_star_pair(int box_corner_x, int box_corner_y)
{
    lat.set_vertex_link(box_corner_x,box_corner_y,3, 1);
    lat.set_vertex_link(box_corner_x,box_corner_y,0, 0);

    lat.set_vertex_link(box_corner_x,box_corner_y-1,0, 0);
   
    lat.set_vertex_link(box_corner_x+1,box_corner_y-1,0, 0);

    lat.set_vertex_link(box_corner_x+1,box_corner_y,0, 0);
    lat.set_vertex_link(box_corner_x+1,box_corner_y,3, 0);

    lat.set_vertex_link(box_corner_x+2,box_corner_y,3, 1);
}
void Model::mv_monomer(int x, int y, int add_x, int add_y)
{

    int length = lat.length;
    int height = lat.height;

    //cout << "moving monomer from: " << x << "," << y << " ... to ... " << (x+add_x) <<","<< (y+add_y) << endl;
        
    // make sure we moved something
    int not_found = 1;

    x = (x+length)%length;
    y = (y+height)%height;

    int new_x;
    int new_y;

    // first search for the monomer with the location x,y
    for (int i=0; i < monomer_locs.size(); i++)
    {
        //cout << "started at " << endl;
        //cout << "monomer_locs["<<i<<"][0] " << monomer_locs[i][0] << endl;
        //cout << "monomer_locs["<<i<<"][1] " << monomer_locs[i][1] << endl;

        // if this is it get rid of it
        if ((monomer_locs[i][0]==x)&&(monomer_locs[i][1]==y))
        {
            // we want to take care of the modulus operation for the periodic lattice in here
            new_x = (x + add_x + length) % length;
            new_y = (y + add_y + height) % height;

            monomer_locs[i][0]=new_x;
            monomer_locs[i][1]=new_y;
            not_found=0;

            //cout << "moved to " << endl;
            //cout << "monomer_locs["<<i<<"][0] " << monomer_locs[i][0] << endl;
            //cout << "monomer_locs["<<i<<"][1] " << monomer_locs[i][1] << endl;
        }
    }


    if (not_found)
    {
        cout << "EXIT_FAILURE: trying to move monomer but not found"<<endl;
        exit(EXIT_FAILURE);
    }
    //else
    //{
    //    cout << " sucsessful monomer move" << endl;
    //}
}
//******************************************************************************************************************************************************************************************
//******************************************************************************************************************************************************************************************
//******************************************************************************************************************************************************************************************
//******************************************************************************************************************************************************************************************

void Model::move_ur(int loc_x, int loc_y)
{
    //int first_important_dimer = find_a_vertex_up_link(lattice[(loc_y-1+height)%height][(loc_x+2)%length]);
    int first_important_dimer = lat.find_a_vertex_up_link(loc_x+2,loc_y-1);
    
    if (first_important_dimer == 1)
    {
        lat.set_vertex_link(loc_x+3,loc_y-2,1, 0);
        lat.set_vertex_link(loc_x+3,loc_y-2,3, 1);

        lat.set_vertex_link(loc_x+2,loc_y-3,0, 0);
        lat.set_vertex_link(loc_x+2,loc_y-3,2, 1);

        lat.set_vertex_link(loc_x+1,loc_y-2,3, 0);
        lat.set_vertex_link(loc_x+1,loc_y-2,0, 0);
        lat.set_vertex_link(loc_x+1,loc_y-2,1, 0);

        lat.set_vertex_link(loc_x+2,loc_y-1,0, 1);
        lat.set_vertex_link(loc_x+2,loc_y-1,3, 1);
        lat.set_vertex_link(loc_x+2,loc_y-1,2, 1);
    }
    if (first_important_dimer == 0)
    {
        lat.set_vertex_link(loc_x+1,loc_y-0,0, 0);
        lat.set_vertex_link(loc_x+1,loc_y-0,2, 1);

        lat.set_vertex_link(loc_x+0,loc_y-1,1, 0);
        lat.set_vertex_link(loc_x+0,loc_y-1,3, 1);

        lat.set_vertex_link(loc_x+1,loc_y-2,0, 0);
        lat.set_vertex_link(loc_x+1,loc_y-2,1, 0);
        lat.set_vertex_link(loc_x+1,loc_y-2,2, 0);

        lat.set_vertex_link(loc_x+2,loc_y-1,1, 1);
        lat.set_vertex_link(loc_x+2,loc_y-1,2, 1);
        lat.set_vertex_link(loc_x+2,loc_y-1,3, 1);
    }

}
void Model::move_dl(int loc_x, int loc_y)
{

    //int first_important_dimer = find_a_vertex_up_link(lattice[(loc_y-2+height)%height][(loc_x+1)%length]);
    int first_important_dimer = lat.find_a_vertex_up_link(loc_x+1,loc_y-2);
    
    if (first_important_dimer == 2)
    {
        lat.set_vertex_link(loc_x+2,loc_y-3,2, 0);
        lat.set_vertex_link(loc_x+2,loc_y-3,0, 1);

        lat.set_vertex_link(loc_x+3,loc_y-2,3, 0);
        lat.set_vertex_link(loc_x+3,loc_y-2,1, 1);

        lat.set_vertex_link(loc_x+2,loc_y-1,2, 0);
        lat.set_vertex_link(loc_x+2,loc_y-1,3, 0);
        lat.set_vertex_link(loc_x+2,loc_y-1,0, 0);

        lat.set_vertex_link(loc_x+1,loc_y-2,3, 1);
        lat.set_vertex_link(loc_x+1,loc_y-2,0, 1);
        lat.set_vertex_link(loc_x+1,loc_y-2,1, 1);
    }
    if (first_important_dimer == 3)
    {
        lat.set_vertex_link(loc_x+1,loc_y-0,2, 0);
        lat.set_vertex_link(loc_x+1,loc_y-0,0, 1);

        lat.set_vertex_link(loc_x+0,loc_y-1,3, 0);
        lat.set_vertex_link(loc_x+0,loc_y-1,1, 1);

        lat.set_vertex_link(loc_x+2,loc_y-1,3, 0);
        lat.set_vertex_link(loc_x+2,loc_y-1,2, 0);
        lat.set_vertex_link(loc_x+2,loc_y-1,1, 0);

        lat.set_vertex_link(loc_x+1,loc_y-2,0, 1);
        lat.set_vertex_link(loc_x+1,loc_y-2,1, 1);
        lat.set_vertex_link(loc_x+1,loc_y-2,2, 1);
    }

}
void Model::move_ul(int loc_x, int loc_y)
{

    //int first_important_dimer = find_a_vertex_up_link(lattice[(loc_y-1+height)%height][(loc_x+1)%length]);
    int first_important_dimer = lat.find_a_vertex_up_link(loc_x+1,loc_y-1);
    
    if (first_important_dimer == 2)
    {
        lat.set_vertex_link(loc_x+2,loc_y-0,2, 0);
        lat.set_vertex_link(loc_x+2,loc_y-0,0, 1);

        lat.set_vertex_link(loc_x+3,loc_y-1,1, 0);
        lat.set_vertex_link(loc_x+3,loc_y-1,3, 1);

        lat.set_vertex_link(loc_x+2,loc_y-2,0, 0);
        lat.set_vertex_link(loc_x+2,loc_y-2,1, 0);
        lat.set_vertex_link(loc_x+2,loc_y-2,2, 0);

        lat.set_vertex_link(loc_x+1,loc_y-1,3, 1);
        lat.set_vertex_link(loc_x+1,loc_y-1,0, 1);
        lat.set_vertex_link(loc_x+1,loc_y-1,1, 1);
    }
    if (first_important_dimer == 1)
    {
        lat.set_vertex_link(loc_x+0,loc_y-2,1, 0);
        lat.set_vertex_link(loc_x+0,loc_y-2,3, 1);

        lat.set_vertex_link(loc_x+1,loc_y-3,2, 0);
        lat.set_vertex_link(loc_x+1,loc_y-3,0, 1);

        lat.set_vertex_link(loc_x+2,loc_y-2,3, 0);
        lat.set_vertex_link(loc_x+2,loc_y-2,2, 0);
        lat.set_vertex_link(loc_x+2,loc_y-2,1, 0);

        lat.set_vertex_link(loc_x+1,loc_y-1,2, 1);
        lat.set_vertex_link(loc_x+1,loc_y-1,3, 1);
        lat.set_vertex_link(loc_x+1,loc_y-1,0, 1);
    }

}
void Model::move_dr(int loc_x, int loc_y)
{

    //int first_important_dimer = find_a_vertex_up_link(lattice[(loc_y-2+height)%height][(loc_x+2)%length]);
    int first_important_dimer = lat.find_a_vertex_up_link(loc_x+2,loc_y-2);
    
    if (first_important_dimer == 0)
    {
        lat.set_vertex_link(loc_x+0,loc_y-2,3, 0);
        lat.set_vertex_link(loc_x+0,loc_y-2,1, 1);

        lat.set_vertex_link(loc_x+1,loc_y-3,0, 0);
        lat.set_vertex_link(loc_x+1,loc_y-3,2, 1);

        lat.set_vertex_link(loc_x+1,loc_y-1,2, 0);
        lat.set_vertex_link(loc_x+1,loc_y-1,3, 0);
        lat.set_vertex_link(loc_x+1,loc_y-1,0, 0);

        lat.set_vertex_link(loc_x+2,loc_y-2,1, 1);
        lat.set_vertex_link(loc_x+2,loc_y-2,2, 1);
        lat.set_vertex_link(loc_x+2,loc_y-2,3, 1);
    }
    if (first_important_dimer == 3)
    {
        lat.set_vertex_link(loc_x+2,loc_y-0,0, 0);
        lat.set_vertex_link(loc_x+2,loc_y-0,2, 1);

        lat.set_vertex_link(loc_x+3,loc_y-1,3, 0);
        lat.set_vertex_link(loc_x+3,loc_y-1,1, 1);

        lat.set_vertex_link(loc_x+1,loc_y-1,3, 0);
        lat.set_vertex_link(loc_x+1,loc_y-1,0, 0);
        lat.set_vertex_link(loc_x+1,loc_y-1,1, 0);

        lat.set_vertex_link(loc_x+2,loc_y-2,0, 1);
        lat.set_vertex_link(loc_x+2,loc_y-2,1, 1);
        lat.set_vertex_link(loc_x+2,loc_y-2,2, 1);
    }
}
// add_x(y) is the amount the star has moved by. Can be positive or negative
// x and y are the initial location of the star
void Model::mv_star(int x, int y, int add_x, int add_y)
{

    //cout << "moving star from: " << x << "," << y << " ... to ... " << (x+add_x) <<","<< (y+add_y) << endl;
    
    int length = lat.length;
    int height = lat.height;
        
    // make sure we moved something
    int not_found = 1;

    x = (x+length)%length;
    y = (y+height)%height;

    int new_x;
    int new_y;

    // first search for the star with the location x,y
    for (int i=0; i < star_locs.size(); i++)
    {
        // if this is it get rid of it
        if ((star_locs[i][0]==x)&&(star_locs[i][1]==y))
        {
            // we want to take care of the modulus operation for the periodic lattice in here
            new_x = (x + add_x + length) % length;
            new_y = (y + add_y + height) % height;

            star_locs[i][0]=new_x;
            star_locs[i][1]=new_y;
            not_found=0;
        }
    }

    if (not_found)
    {
        cout << "EXIT_FAILURE: trying to move star but not found"<< endl;
        exit(EXIT_FAILURE);
    }

}
// these function will be used to simplify the readability of implementing the "keeping track
// of stars"
// add a star to the star_locs. It is at location x
void Model::add_star(int x,int y)
{
    int length = lat.length;
    int height = lat.height;
    //cout << "adding star at: " << x << "," << y <<endl;
    vector <int> to_add;
    to_add.push_back((x+length)%length);
    to_add.push_back((y+height)%height);
    star_locs.push_back(to_add);
}
void Model::rm_star(int x,int y)
{
    int length = lat.length;
    int height = lat.height;

    //cout << "removing star at: " << x << "," << y <<endl;
    int not_found = 1;

    x = (x+length)%length;
    y = (y+height)%height;
    // first search for the star with the location x,y
    for (int i=0; i < star_locs.size(); i++)
    {
        // if this is it get rid of it
        if ((star_locs[i][0]==x)&&(star_locs[i][1]==y))
        {
            star_locs.erase(star_locs.begin()+i);
            not_found = 0;
        }
    }

    if (not_found)
    {
        cout << "EXIT_FAILURE: trying to remove star but not found"<< endl;
        exit(EXIT_FAILURE);
    }

}
//******************************************************************************************************************************************************************************************
//******************************************************************************************************************************************************************************************
//******************************************************************************************************************************************************************************************
//******************************************************************************************************************************************************************************************
bool Model::is_there_a_monomer_can_we_flip(int loc_x, int loc_y)
{
    int length = lat.length;
    int height = lat.height;
    //the plaquette given the that loc_x and loc_y denote the upper left plaquette corner
    int plaq_x = loc_x;
    int plaq_y = loc_y-1;
    // these are the dimer locations on the plaquet. 0 would mean there is a dimer on the right most
    // edge. 1 is a dimer on the top and around 2-> left 3->bottom.
    vector <int> dimer_locs;
    
    // loc_x and loc_y are the locations of the plaquettes upper left had corner
    // go around and see if there is a monomer
    //going clockwise
    int num_links_v0 = lat.num_oc_links(loc_x,loc_y);
    int num_links_v1 = lat.num_oc_links(loc_x+1,loc_y);
    int num_links_v2 = lat.num_oc_links(loc_x+1,loc_y-1);
    int num_links_v3 = lat.num_oc_links(loc_x,loc_y-1);

    // each one has a x and y associated with it
    vector <int> vec_v0; vec_v0.push_back(loc_x); vec_v0.push_back(loc_y);
    vector <int> vec_v1; vec_v1.push_back((loc_x+1)%length); vec_v1.push_back(loc_y);
    vector <int> vec_v2; vec_v2.push_back((loc_x+1)%length); vec_v2.push_back((loc_y-1+height)%height);
    vector <int> vec_v3; vec_v3.push_back(loc_x); vec_v3.push_back((loc_y-1+height)%height);

    vector < vector <int> > all_vec_verts;
    // The purpose of this is so that we can reference a vertex of the plaquette as counted from zero
    // in the upper left corner going up clockwise. The first index of this matrix is the vertex. The
    // second index gives you the x location of that vertex and the third gives you the y location
    // of that vertex.
    all_vec_verts.push_back(vec_v0); all_vec_verts.push_back(vec_v1); all_vec_verts.push_back(vec_v2); all_vec_verts.push_back(vec_v3);

    // make a vector with these objects
    vector <int> v_num_links;
    v_num_links.push_back(num_links_v0);
    v_num_links.push_back(num_links_v1);
    v_num_links.push_back(num_links_v2);
    v_num_links.push_back(num_links_v3);

    // which vertex(ies) has the monomer(s)
    vector <int> monomer_locs_plqt_vrts;
    for (int i=0; i<v_num_links.size(); i++)
    {
        if (v_num_links[i]==0)
        {
            monomer_locs_plqt_vrts.push_back(i);
        }
        
    }

    if ((monomer_locs_plqt_vrts.size()>=1)&&(monomer_locs_plqt_vrts.size()<=2))
    {
        // then we have a monomer
        // can only continue if there is one link on the plaquette
        for (int i=0; i<4; i++)
        {
            if (lat.get_plaquette_link(plaq_x,plaq_y,i)==1)
            {
                dimer_locs.push_back(i);
            }
            
        }
        // the monomer (of the two max) that we will move
        int monomer_to_move = 0;
        int cur_mon_x;
        int cur_mon_y;
        int new_mon_x;
        int new_mon_y;
        if (dimer_locs.size()==1)
        {
            // if there are two monomers we must randomly select one to move
            if (monomer_locs_plqt_vrts.size()==2)
            {
                monomer_to_move= rdm.get_rand_int(0,2);
            }
            if (monomer_locs_plqt_vrts.size()>=3)
            {
                cout << "MORE THAN 2 MONOMERS ??????????"<<endl;
                exit(EXIT_FAILURE);
            }

            // then move the monomer to the same sublattice opposite corner of the plaquette and
            // flip the dimer appropriately 
            // If you draw some pictures youll find that if you add the dimer location to the
            // monomer location and its odd you subtract 1 from the dimer location -> rotate
            // clockwise
            // If you add them and get an even number then you add 1 to the dimer location ->
            // counter clockwise rotate it. 
            // All we need to do is rotate the dimer.
            if ((monomer_locs_plqt_vrts[monomer_to_move]+dimer_locs[0])%2==0) // if even
            {
                // set the current dimer location to 0    
                lat.set_plaquette_link(plaq_x,plaq_y, dimer_locs[0], 0);
                //cout << "even setting plaquette ("<<plaq_x<<","<<plaq_y<<") link "<<dimer_locs[0]<<" to 0"<<endl;
                // set the old dimer location + 1 %4 to 1
                lat.set_plaquette_link(plaq_x,plaq_y,(dimer_locs[0]+1)%4, 1);
                //cout << "even setting plaquette ("<<plaq_x<<","<<plaq_y<<") link "<<(dimer_locs[0]+1)%4<<" to 1 "<<endl;

                // this is ridiculous but the monomer vertex number is the (4 -dimer location)%4 FOR
                cur_mon_x = all_vec_verts[(4-dimer_locs[0])%4][0];
                cur_mon_y = all_vec_verts[(4-dimer_locs[0])%4][1];
                new_mon_x = all_vec_verts[(4-dimer_locs[0]+2)%4][0];
                new_mon_y = all_vec_verts[(4-dimer_locs[0]+2)%4][1];

                //mv_monomer(all_vec_verts[(4-dimer_locs[0])%4][0],all_vec_verts[(4-dimer_locs[0])%4][1],all_vec_verts[(4-dimer_locs[0]+2)%4][0],all_vec_verts[(4-dimer_locs[0]+2)%4][1]);
                //cout << "inputing even: x="<<cur_mon_x<< ", y="<<cur_mon_y<<endl;
                mv_monomer(cur_mon_x,cur_mon_y,new_mon_x-cur_mon_x,new_mon_y-cur_mon_y);
                return true;
            }
            else if ((monomer_locs_plqt_vrts[monomer_to_move]+dimer_locs[0])%2==1) // if odd
            {
                // set the current dimer location to 0    
                lat.set_plaquette_link(plaq_x,plaq_y, dimer_locs[0], 0);
                //cout << "odd setting plaquette ("<<plaq_x<<","<<plaq_y<<") link "<<dimer_locs[0]<<" to 0"<<endl;
                // set the old dimer location + 1 %4 to 1
                lat.set_plaquette_link(plaq_x,plaq_y, (dimer_locs[0]-1+4)%4, 1);
                //cout << "odd setting plaquette ("<<plaq_x<<","<<plaq_y<<") link "<<(dimer_locs[0]+1)%4<<" to 1 "<<endl;

                // this one is (3-dimer locations)%4
                cur_mon_x = all_vec_verts[(3-dimer_locs[0])%4][0];
                cur_mon_y = all_vec_verts[(3-dimer_locs[0])%4][1];
                new_mon_x = all_vec_verts[(3-dimer_locs[0]+2)%4][0];
                new_mon_y = all_vec_verts[(3-dimer_locs[0]+2)%4][1];
                //cout << "inputing odd: x="<<cur_mon_x<< ", y="<<cur_mon_y<<endl;
                mv_monomer(cur_mon_x,cur_mon_y,new_mon_x-cur_mon_x,new_mon_y-cur_mon_y);
                //mv_monomer(all_vec_verts[(3-dimer_locs[0])%4][0],all_vec_verts[(3-dimer_locs[0])%4][1],all_vec_verts[(3-dimer_locs[0]+2)%4][0],all_vec_verts[(3-dimer_locs[0]+2)%4][1]);
                return true;
            }
            else 
            {
                cout << "Something wrong with moving the monomers. function: is_there_a_monomer_can_we_flip"<<endl;
                exit(EXIT_FAILURE);
            }
        }
        else { return false; }
         
    }
    else { return false; }
}
//******************************************************************************************************************************************************************************************
//******************************************************************************************************************************************************************************************
//******************************************************************************************************************************************************************************************
//******************************************************************************************************************************************************************************************
void Model::set_move_vectors()
{
    add_horz_star_walk.push_back(0);  add_horz_star_there.push_back(1); add_horz_star_links.push_back(1);
    add_horz_star_walk.push_back(3);  add_horz_star_there.push_back(0); add_horz_star_links.push_back(1);
    add_horz_star_walk.push_back(2);  add_horz_star_there.push_back(0); add_horz_star_links.push_back(1);
    add_horz_star_walk.push_back(3);  add_horz_star_there.push_back(0); add_horz_star_links.push_back(1);
    add_horz_star_walk.push_back(0);  add_horz_star_there.push_back(1); add_horz_star_links.push_back(1); add_horz_star_links.push_back(1);

    annihilate_horz_star_walk.push_back(0);  annihilate_horz_star_there.push_back(0); annihilate_horz_star_links.push_back(1);
    annihilate_horz_star_walk.push_back(3);  annihilate_horz_star_there.push_back(1); annihilate_horz_star_links.push_back(1);
    annihilate_horz_star_walk.push_back(2);  annihilate_horz_star_there.push_back(1); annihilate_horz_star_links.push_back(4);
    annihilate_horz_star_walk.push_back(3);  annihilate_horz_star_there.push_back(1); annihilate_horz_star_links.push_back(4);
    annihilate_horz_star_walk.push_back(0);  annihilate_horz_star_there.push_back(0); annihilate_horz_star_links.push_back(1); annihilate_horz_star_links.push_back(1);

    add_vert_star_walk.push_back(3);  add_vert_star_there.push_back(1); add_vert_star_links.push_back(1);
    add_vert_star_walk.push_back(0);  add_vert_star_there.push_back(0); add_vert_star_links.push_back(1);
    add_vert_star_walk.push_back(1);  add_vert_star_there.push_back(0); add_vert_star_links.push_back(1);
    add_vert_star_walk.push_back(0);  add_vert_star_there.push_back(0); add_vert_star_links.push_back(1);
    add_vert_star_walk.push_back(3);  add_vert_star_there.push_back(1); add_vert_star_links.push_back(1); add_vert_star_links.push_back(1);

    annihilate_vert_star_walk.push_back(3);  annihilate_vert_star_there.push_back(0); annihilate_vert_star_links.push_back(1);
    annihilate_vert_star_walk.push_back(0);  annihilate_vert_star_there.push_back(1); annihilate_vert_star_links.push_back(1);
    annihilate_vert_star_walk.push_back(1);  annihilate_vert_star_there.push_back(1); annihilate_vert_star_links.push_back(4);
    annihilate_vert_star_walk.push_back(0);  annihilate_vert_star_there.push_back(1); annihilate_vert_star_links.push_back(4);
    annihilate_vert_star_walk.push_back(3);  annihilate_vert_star_there.push_back(0); annihilate_vert_star_links.push_back(1); annihilate_vert_star_links.push_back(1);

    move_star_right_walk.push_back(0);  move_star_right_there.push_back(0); move_star_right_links.push_back(1);
    move_star_right_walk.push_back(0);  move_star_right_there.push_back(1); move_star_right_links.push_back(1);
    move_star_right_walk.push_back(3);  move_star_right_there.push_back(0); move_star_right_links.push_back(1);
    move_star_right_walk.push_back(3);  move_star_right_there.push_back(0); move_star_right_links.push_back(1);
    move_star_right_walk.push_back(2);  move_star_right_there.push_back(1); move_star_right_links.push_back(1);
    move_star_right_walk.push_back(1);  move_star_right_there.push_back(0); move_star_right_links.push_back(1);
    move_star_right_walk.push_back(2);  move_star_right_there.push_back(1); move_star_right_links.push_back(1); move_star_right_links.push_back(4);

    move_star_left_walk.push_back(0);  move_star_left_there.push_back(1); move_star_left_links.push_back(1);
    move_star_left_walk.push_back(3);  move_star_left_there.push_back(0); move_star_left_links.push_back(1);
    move_star_left_walk.push_back(0);  move_star_left_there.push_back(1); move_star_left_links.push_back(1);
    move_star_left_walk.push_back(3);  move_star_left_there.push_back(1); move_star_left_links.push_back(4);
    move_star_left_walk.push_back(2);  move_star_left_there.push_back(0); move_star_left_links.push_back(1);
    move_star_left_walk.push_back(2);  move_star_left_there.push_back(1); move_star_left_links.push_back(1);
    move_star_left_walk.push_back(1);  move_star_left_there.push_back(0); move_star_left_links.push_back(1); move_star_left_links.push_back(1);

    move_star_up_walk.push_back(3);  move_star_up_there.push_back(1); move_star_up_links.push_back(1);
    move_star_up_walk.push_back(0);  move_star_up_there.push_back(0); move_star_up_links.push_back(1);
    move_star_up_walk.push_back(3);  move_star_up_there.push_back(1); move_star_up_links.push_back(1);
    move_star_up_walk.push_back(0);  move_star_up_there.push_back(1); move_star_up_links.push_back(4);
    move_star_up_walk.push_back(1);  move_star_up_there.push_back(0); move_star_up_links.push_back(1);
    move_star_up_walk.push_back(1);  move_star_up_there.push_back(1); move_star_up_links.push_back(1);
    move_star_up_walk.push_back(2);  move_star_up_there.push_back(0); move_star_up_links.push_back(1); move_star_up_links.push_back(1);

    move_star_down_walk.push_back(3);  move_star_down_there.push_back(0); move_star_down_links.push_back(1);
    move_star_down_walk.push_back(3);  move_star_down_there.push_back(1); move_star_down_links.push_back(1);
    move_star_down_walk.push_back(0);  move_star_down_there.push_back(0); move_star_down_links.push_back(1);
    move_star_down_walk.push_back(0);  move_star_down_there.push_back(0); move_star_down_links.push_back(1);
    move_star_down_walk.push_back(1);  move_star_down_there.push_back(1); move_star_down_links.push_back(1);
    move_star_down_walk.push_back(2);  move_star_down_there.push_back(0); move_star_down_links.push_back(1);
    move_star_down_walk.push_back(1);  move_star_down_there.push_back(1); move_star_down_links.push_back(1); move_star_down_links.push_back(4);

    //move_star_up_right_diag_walk.push_back();
           
}

//******************************************************************************************************************************************************************************************
//******************************************************************************************************************************************************************************************
//******************************************************************************************************************************************************************************************
//******************************************************************************************************************************************************************************************
int Model::flip(int to_flip)
{
    int new_spin;

    if (to_flip == 0){new_spin = 1;}
    else if (to_flip == 1){new_spin = 0;}
    else {cout << "BIG problem: spin at link not 0 or 1\n";}

    return new_spin;
}
