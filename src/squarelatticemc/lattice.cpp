
#include "lattice.hpp"


Lattice::Lattice()
{
    length = 0;
    height = 0;
}
Lattice::Lattice(int to_set_length, int to_set_height,O_Random rdm_in)
{
    length = to_set_length;
    height = to_set_height;
    num_vert = length*height;
    links_per_vert = 4;
    tot_num_links = 2*num_vert;

    rdm = rdm_in;
    /* An array ordered from left to right containing the spins of all the links -> left to
     * right for each row -> then up -> then left to right again.
     * initialize it at 0 (i.e) nothing.
     * There are 2 unique lattice points per vertex so length*height*2 is number of links */
    for (int i = 0; i < tot_num_links ; i++)
    {
        link_vals.push_back(0);
    }
    //cout << "link_vals.size(): " << link_vals.size() << endl;

    //with the loop below we are trying to achieve 
    //    -> vector < vector < Vertex > > v_arr (height, vector <Vertex> (length));
    for (int j = 0; j < height; j++)
    {
        vector <Vertex> to_add (length);
        v_arr.push_back(to_add);
    }

    weave_links_verts();
}

/* This function takes the vertex (vert_x and vert_y), the integer describing right left up down 
 * link from it (where) and returns the integer that is the index of this link in the link_vals arr.
 *
 * where variable 1-4: 
 * set the links in order of an angle being swept out from x axis. In this case: 
 * 0 -> right link, 1 -> up link, 2 -> left link, 3 -> down link.
 */

int Lattice::find_link(int vert_x, int vert_y, int where)
{

    // make sure the numbers inputed are on the lattice with mod
    vert_y = (vert_y + height) % height;
    vert_x = (vert_x + length) % length;

    // to return
    int the_link;

    if (where == 0) { the_link = 2*vert_y*length + vert_x; }

    if (where == 1) { the_link = (2*vert_y+1)*length+vert_x; }

    if (where ==2)
    {
        /* need to be careful on left most column of vertices with boundary conditions */
        if (vert_x == 0) { the_link = (2*vert_y +1)*length -1; }
        else { the_link = 2*vert_y*length+vert_x-1 ;}
    }

    if (where ==3 )
    {
        if (vert_y ==0 ) {the_link = tot_num_links - (length-vert_x);}
        else {the_link = (2*vert_y-1)*length+vert_x ; }
    }

    return the_link;
}

/* Initialize the square matrix (lattice) */
/* Remember that v_arr is ordered like a matrix but bottom to top*/
void Lattice::weave_links_verts()
{
    //this is the index in (link_vals) of the link we want
    int the_link;

    for (int i = 0; i<height; i++)
    {
        for (int j = 0; j<length; j++)
        {

            //closed parens at end of Vertex creates temporary that gets permanently stored in vector
            v_arr[i][j] = Vertex();
            // define the origin to be sublattice A 
            // arrows point away from sublattice A to sublattice B
            if (((i+j)%2)==0)
            {
                v_arr[i][j].set_sublat(0);
            }
            else 
            {
                v_arr[i][j].set_sublat(1);
            }

            //v_arr[i][j].set_num_link(links_per_vert);
            v_arr[i][j].num_link=links_per_vert;
            //cout << "links_per_vert: " << links_per_vert << endl;

            for (int a = 0; a < links_per_vert; a++)
            {
                the_link = find_link(j, i, a);
                //cout << "the_link: " ;
                //cout << the_link <<"\n" ;
                //v_arr[i][j].links[a] =  link_vals[the_link];
                v_arr[i][j].links.push_back( & link_vals[the_link]);
                //cout << "link_vals[the_link]: " << link_vals[the_link] << endl;

            }

            // We also might need to set the "color" of each vertex so that we can distinguish between
            // the two sublattices. color is just an integer 0 or 1.
            //if ((i+j)%2 == 0 ) { (*v_arr)[i][j].set_color(0);}
            //else {(*v_arr)[i][j].set_color(1);}


        }
    }
}

// needs to return one at random
int Lattice::find_a_vertex_up_link(int vert_x, int vert_y)
{

    // make sure the numbers inputed are on the lattice with mod
    vert_y = (vert_y + height) % height;
    vert_x = (vert_x + length) % length;

    Vertex vert = v_arr[vert_y][vert_x];

    int go_to_link = -1;
    int rand;

    // need a vector to store the indices of the up link
    vector <int> index_up;
    
    for (int i = 0; i < vert.num_link; i++ ) 
    {
        if (* vert.links[i] == 1)
        {
            index_up.push_back(i);
        }
    }
    
    if (index_up.size()!=0)
    {
        rand = rdm.get_rand_int(0,index_up.size());
        go_to_link = index_up[rand];
    }
    
    return go_to_link;    
 
}
// needs to return one at random
int Lattice::find_a_vertex_down_link(int vert_x, int vert_y)
{

    // make sure the numbers inputed are on the lattice with mod
    vert_y = (vert_y + height) % height;
    vert_x = (vert_x + length) % length;

    Vertex vert = v_arr[vert_y][vert_x];
        
    int go_to_link = -1;
    int rand;

    // need a vector to store the indices of the up link
    vector <int> index_down;
    
    for (int i = 0; i < vert.num_link; i++ ) 
    {
        if (* vert.links[i] == 0)
        {
            index_down.push_back(i);
        }
    }
    
    if (index_down.size()!=0)
    {
        rand = rdm.get_rand_int(0,index_down.size());
        go_to_link = index_down[rand];
    }

    return go_to_link;    
}

int Lattice::num_oc_links(int vert_x, int vert_y){
    // make sure the numbers inputed are on the lattice with mod
    vert_y = (vert_y + height) % height;
    vert_x = (vert_x + length) % length;

    return v_arr[vert_y][vert_x].num_oc_links();
}
int Lattice::get_z3_vison_phase(int start_x, int start_y, int end_x, int end_y,const Lattice& br_lat)
{

    // working in dual lattice so lets use the plaquettes
    int cur_plaq_x = start_x;
    int cur_plaq_y = start_y;

    // using the starting and ending points we can find the directions to go. First we move along
    // the x axis -> and then along the y axis so we only need two directions. same direction
    // numbers as those used in labeling plaquette edges

    // 4 is an invalid direction
    int direction_x = 4;
    int direction_y = 4;

    if      ((end_x-start_x)<0) { direction_x = 2; }
    else if ((end_x-start_x)>0) { direction_x = 0; }
    if      ((end_y-start_y)<0) { direction_y = 3; }
    else if ((end_y-start_y)>0) { direction_y = 1; }

    int amount_x = abs(end_x-start_x);
    int amount_y = abs(end_y-start_y);

    int cur_sublat;
    int cur_link;
    // background link
    int cur_br_link;

    // This variable will keep track of N where the phase is N*2pi/3
    int vison = 0;
    // this variable will keep track of the gauge fixing that needs to be done in accordance with the
    // background/reference dimerization(See Ivanov paper). N_gauge_fix
    int gauge_fix = 0;


    // Now lets go along the path and be careful with the gauge fixing   
    // The path is just straight up. nothing crazy
    for (int i = 0; i < amount_x; ++i)
    {
        //cout << "cur_plaq_x: " << cur_plaq_x << endl;
        //cout << "cur_plaq_y: " << cur_plaq_y << endl;
        cur_link = get_plaquette_link(cur_plaq_x,cur_plaq_y,direction_x);
        // this depends on which direction you are going.
        if (direction_x == 0)
        {
            //left most (from path perspective) vertex
            cur_sublat = get_vertex_sub_lat(cur_plaq_x+1,cur_plaq_y+1);
        }
        if (direction_x == 2)
        {
            //left most (from path prespective) vertex
            cur_sublat = get_vertex_sub_lat(cur_plaq_x,cur_plaq_y);
        }
        //cout << "cur_link: " << cur_link << endl;
        // REMEMBER THE GAUGE FIXING DIMERS GO THE OTHER WAY. (a gauge fix dimer on top of a normal
        // dimer = nothing)
        if (cur_link == 1)
        {
            // need to get the sublattice to know which way the dimers are going
            if (cur_sublat == 1)
            {
                vison ++; 
            }
            else if (cur_sublat == 0)
            {
                vison --;
            }
            else
            {
                cout << "Invalid sub lattice value in get_z3_vison_phase. cur_sublat: "<< cur_sublat <<endl;
                exit(EXIT_FAILURE);
            }
        }

        //cout << "vison: " << vison << endl;

        // be carful with the gauge fixing!! this comes from "br_lat" (background lattice) 
        // link the path is crossing over
        cur_br_link = br_lat.get_plaquette_link(cur_plaq_x,cur_plaq_y,direction_x);
        //cout << "cur_br_link: " << cur_br_link << endl;
        // the sublattice does not depend on the configuration -> can use the cur_sublat that we
        // already have BUUUUUTTTT the gauge fixing changes the phase in the other directions.
        if (cur_br_link == 1)
        {
            // need to get the sublattice to know which way the dimers are going
            if (cur_sublat == 0)
            {
                gauge_fix ++; 
            }
            else if (cur_sublat == 1)
            {
                gauge_fix --;
            }
        }
        //cout << "gauge_fix : " << gauge_fix << endl;

        if (direction_x==0){
            cur_plaq_x++;
        }
        else if (direction_x==2){
            cur_plaq_x--;
        }
        else{
            cout << "problem in getting z2 vison"<<endl;
            cout << "direction_x: "<< direction_x<<endl;
            cout << "start_x: "<< start_x <<endl;
            cout << "end_x: "<< end_x <<endl;
            cout << "start_y: "<< start_y <<endl;
            cout << "end_y: "<< end_y <<endl;
            cout << "amount_x: "<< amount_x <<endl;
            cout << "amount_y: "<< amount_y <<endl;
            exit(EXIT_FAILURE);
        }
    }
    for (int i = 0; i < amount_y; ++i)
    {
        //cout << "cur_plaq_x: " << cur_plaq_x << endl;
        //cout << "cur_plaq_y: " << cur_plaq_y << endl;
        cur_link = get_plaquette_link(cur_plaq_x,cur_plaq_y,direction_y);
        // this depends on which direction you are going.
        if (direction_y == 3)
        {
            //left most (from path perspective) vertex
            cur_sublat = get_vertex_sub_lat(cur_plaq_x+1,cur_plaq_y);
        }
        if (direction_y == 1)
        {
            //left most (from path perspective) vertex
            cur_sublat = get_vertex_sub_lat(cur_plaq_x,cur_plaq_y+1);
        }
        //cout << "cur_link: " << cur_link << endl;
        // REMEMBER THE GAUGE FIXING DIMERS GO THE OTHER WAY. (a gauge fix dimer on top of a normal
        // dimer = nothing)
        if (cur_link == 1)
        {
            // need to get the sublattice to know which way the dimers are going
            if (cur_sublat == 1)
            {
                vison ++; 
            }
            else if (cur_sublat == 0)
            {
                vison --;
            }
            else
            {
                cout << "Invalid sub lattice value in get_z3_vison_phase. cur_sublat: "<< cur_sublat <<endl;
                exit(EXIT_FAILURE);
            }
        }

        //cout << "vison: " << vison << endl;

        // be careful with the gauge fixing!! this comes from "br_lat" (background lattice) 
        // link the path is crossing over
        cur_br_link = br_lat.get_plaquette_link(cur_plaq_x,cur_plaq_y,direction_y);
        //cout << "cur_br_link: " << cur_br_link << endl;
        // the sublattice does not depend on the configuration -> can use the cur_sublat that we
        // already have BUUUUUTTTT the gauge fixing changes the phase in the other directions.
        if (cur_br_link == 1)
        {

            // need to get the sublattice to know which way the dimers are going
            if (cur_sublat == 0)
            {
                gauge_fix ++; 
            }
            else if (cur_sublat == 1)
            {
                gauge_fix --;
            }
        }
        //cout << "gauge_fix : " << gauge_fix << endl;

        if (direction_y==1){
            cur_plaq_y++;
        }
        else if (direction_y==3){
            cur_plaq_y--;
        }
        else{
            cout << "problem in getting z2 vison"<<endl;
            cout << "direction_y: "<< direction_x<<endl;
            cout << "start_x: "<< start_x <<endl;
            cout << "end_x: "<< end_x <<endl;
            cout << "start_y: "<< start_y <<endl;
            cout << "end_y: "<< end_y <<endl;
            cout << "amount_x: "<< amount_x <<endl;
            cout << "amount_y: "<< amount_y <<endl;
            exit(EXIT_FAILURE);
        }
    }

    //cout << "final vison: " << vison << endl;
    //cout << "final gauge_fix: " << gauge_fix << endl;

    return (vison+gauge_fix);
}
vector<int> Lattice::count_number_horz_vert_dimers()
{
    int num_horz_link = 0;
    int num_vert_link = 0;

    for (int i=0; i<length; i++)
    {
        for (int j=0; j<height; j++)
        {
            num_horz_link += get_vertex_link(i,j,0);
            num_vert_link += get_vertex_link(i,j,1);

        }
    }
    vector <int> num_horz_vert;
    num_horz_vert.push_back(num_horz_link);
    num_horz_vert.push_back(num_vert_link);

    return num_horz_vert;

}



// ****************************************************************************************************************************************
// ***************** some functions to make the access of the lattice more accessible ******************************************************
// ****************************************************************************************************************************************
// "get_vertex_link(<x location of vertex>,<y location>,<link number>" returns the spin (1 or 0) in the 
// <link number> direction> of the vertex at coordinates <x location> and <y location>.
// Remember <link number> is from 0 to 3 going counter clockwise starting with 0 to the right
int Lattice::get_vertex_link(int vertex_x, int vertex_y, int link_num)
const
{


    // make sure link num is a valid number
    
    if ((link_num > 3)||(link_num<0))
    {
        cout << "Invalid link_number in function get_vertex_link"<<endl;
        exit(EXIT_FAILURE);
    }

    // make sure the numbers inputed are on the lattice with mod
    vertex_y = (vertex_y + height) % height;
    vertex_x = (vertex_x + length) % length;

    if ((* v_arr[vertex_y][vertex_x].links[link_num])>1)
    {
        cout << "Invalid return value get_vertex_link"<<endl;
        cout << "value from get_vertex_link: " << * v_arr[vertex_y][vertex_x].links[link_num] << endl;
        cout << "vertex_y: " << vertex_y << endl;
        cout << "vertex_x: " << vertex_x << endl;
        cout << "link_num: " << link_num << endl;
        exit(EXIT_FAILURE);
    }

    return * v_arr[vertex_y][vertex_x].links[link_num];
    
}
void Lattice::check_for_defects(void)
{
    for (int h = 0; h < height; h++ )
    {
        for (int l=0; l < length; l++ )
        {
            for (int o=0; o<4; o++)
            {
                if (*v_arr[h][l].links[o]>1)
                {
                    cout << "Defect found" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
}


// we also want the same thing as above but in the plaquette perspective. In this case each
// plaquette is numbered with the plaquette at the origin being the (0,0) plaquette and the one to
// the right being the (1,0) ect. link_num still goes from 0-3 counter clockwise with the
// 0th link being the right most edge of the plaquette.
int Lattice::get_plaquette_link(int plaquette_x, int plaquette_y, int link_num)
const
{
    //cout << "plaquette_x: " << plaquette_x << endl;
    //cout << "plaquette_y: " << plaquette_y << endl;
    //cout << "link_num: " << link_num<< endl;

    // to be really careful we can do this by using the get_vertex_link function --> This way we are
    // sure not to grab the wrong piece of data or someting off of the lattice.
    //
    // the plaquette at plaquette_x and plaquette_y are the upper right diagonal plaquette of vertex_x
    // and vertex_y when plaquette_x=vertex_x and plaquette_y=vertex_y.

    // link numbers 2 and 3 we can do
    if (link_num==2)
    {
        return get_vertex_link(plaquette_x, plaquette_y, 1);
    }
    else if (link_num==3)
    {
        return get_vertex_link(plaquette_x, plaquette_y, 0);
    }
    // for link number 0-1 we need to hop to the upper diagonal vertex
    else if (link_num==0)
    {
        return get_vertex_link(plaquette_x+1, plaquette_y+1, 3);
    }
    else if (link_num==1)
    {
        return get_vertex_link(plaquette_x+1, plaquette_y+1, 2);
    }
    // make sure link num is a valid number
    else
    {
        cout << "Invalid link_number in function get_plaquet_link"<<endl;
        exit(EXIT_FAILURE);
    }
}

// Also we want to be able to set the links to new values. reference in the same way above but
// include one more variable which is the value we want to set the link to.
void Lattice::set_vertex_link(int vertex_x, int vertex_y, int link_num, int set_to)
{
    // make sure link num is a valid number
    if ((link_num > 3)||(link_num<0))
    {
        cout << "Invalid link_number in function get_vertex_link"<<endl;
        exit(EXIT_FAILURE);
    }
    // make sure set_to is a valid number
    if ((set_to!=1)&&(set_to!=0))
    {
        cout << "You are trying to set a vertex link to an invalid spin: "<< set_to <<endl;
        exit(EXIT_FAILURE);
    }

    // make sure the numbers inputed are on the lattice with mod
    vertex_y = (vertex_y + height) % height;
    vertex_x = (vertex_x + length) % length;

    * v_arr[vertex_y][vertex_x].links[link_num]=set_to;
}
void Lattice::set_plaquette_link(int plaquette_x, int plaquette_y, int link_num, int set_to)
{
    // make sure set_to is a valid number
    if ((set_to!=1)&&(set_to!=0))
    {
        cout << "You are trying to set a plaquette link to an invalid spin: "<< set_to <<endl;
        exit(EXIT_FAILURE);
    }

    // to be really careful we can do this by using the get_vertex_link function --> This way we are
    // sure not to grab the wrong piece of data or something off of the lattice.
    //
    // the plaquette at plaquette_x and plaquette_y are the upper right diagonal plaquette of vertex_x
    // and vertex_y when plaquette_x=vertex_x and plaquette_y=vertex_y.

    // link numbers 2 and 3 we can do
    if (link_num==2)
    {
        set_vertex_link(plaquette_x, plaquette_y, 1, set_to);
    }
    else if (link_num==3)
    {
        set_vertex_link(plaquette_x, plaquette_y, 0, set_to);
    }
    // for link number 0-1 we need to hop to the upper diagonal vertex
    else if (link_num==0)
    {
        set_vertex_link(plaquette_x+1, plaquette_y+1, 3, set_to);
    }
    else if (link_num==1)
    {
        set_vertex_link(plaquette_x+1, plaquette_y+1, 2, set_to);
    }
    // make sure link num is a valid number
    else
    {
        cout << "Invalid link_number in function get_vertex_link"<<endl;
        exit(EXIT_FAILURE);
    }
}

// For the Z3 visons we need to know the orientation of the links and this is pretty easy to
// do if you know which sublattice a particular vertex belongs to
// Sublattice A starts at the origin
// Assume arrows point away from sublat A and to sublat B
// This returns an integer value
//      0 -> arrow to the right
//      1 -> arrow to the left
int Lattice::get_vertex_sub_lat(int vertex_x, int vertex_y)
{
    vertex_y = (vertex_y + height) % height;
    vertex_x = (vertex_x + length) % length;
    
    int sum_v = vertex_x + vertex_y;

    if ((sum_v%2) == 0) {
        return 0;
    }   
    else {
        return 1; 
    }

}

// Lets make a function that allows the user to specify a specific plaquette and a link from
// that plaquette. Assuming you are moving in the direction of the link specified this
// function will return whether you intersect an arrow going to the right of the traveling
// path or to the left. 
// This returns an integer value
//      0 -> arrow to the right
//      1 -> arrow to the left
int Lattice::intersect_link_direction(int plaquette_x, int plaquette_y, int link_num)
{
    if ((link_num > 3)||(link_num<0))
    {
        cout << "Invalid link_number in function get_vertex_link"<<endl;
        exit(EXIT_FAILURE);
    }

    int intersect_direc = -1;

    // Need to find a vertex that the link is attached to
    int vert_x;
    int vert_y;

    vert_y = (plaquette_y + height) % height;
    vert_x = (plaquette_x + length) % length;

    if (get_vertex_sub_lat(vert_x,vert_y)==0)
    {
        if ((link_num == 0 )||(link_num == 2))
        {
            intersect_direc = 0;
        }        
        else
        {
            intersect_direc = 1;
        }
    }
    else if (get_vertex_sub_lat(vert_x,vert_y)==1)
    {
        if ((link_num == 0 )||(link_num == 2))
        {
            intersect_direc = 1;
        }        
        else
        {
            intersect_direc = 0;
        }
    }
    else 
    {
        cout << "in Lattice::intersect_link_direction"<<endl;
        cout << "--> Didnt get sublattice A (0) or B (1)"<<endl;
        exit(EXIT_FAILURE);
    }
    if (intersect_direc == -1)
    {
        cout << "intersect_link_direction() did not properly retrieve the direction of the intersection"<<endl;
        exit(EXIT_FAILURE);
    }

    return intersect_direc;
}


// Not staggered. Just full plaquettes with their vertical sides removed. Columnar
void Lattice::build_init_fully_packed_dimer()
{


    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < length; j ++)
        {
            if (j%2 == 0)
            {

                // Now we are at the bottom left corner of the plaquette to make
                // put the link to the right to 1
                * v_arr[i][j].links[0] = 1; //right link
            }
        }
    }
}
// Same as above but columnar pattern with vertical links
void Lattice::build_init_horz_columns()
{

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < length; j ++)
        {
            if (i%2 == 0)
            {

                // Now we are at the bottom left corner of the plaquette to make
                // put the link to the right to 1
                * v_arr[i][j].links[1] = 1; //right link
            }
        }
    }
}

// goes around the initialized lattice changing the dimers in a nontrivial loop
void Lattice::increase_winding_num(int by_how_much)
{
    int cur_link_val;

    for (int i = 0; i < by_how_much*2; i+=2)
    {
        for (int j = 0; j < length; j++)
        {
            cur_link_val = get_plaquette_link(j,i,1);
            if (cur_link_val == 0)
            {
               set_plaquette_link(j,i,1,1); 
            }
            else
            {
               set_plaquette_link(j,i,1,0); 
            }
        }
    }
        
}
// goes around the horizontal dimers (vertical column) in the vertical direction to increase Wx by 1
void Lattice::increase_windingX_by_one(int col_index)
{
    int cur_link_val;

    int orig_col_index = col_index;

    bool did_left = true;

            //cout << "In here" << endl;
            //exit(EXIT_FAILURE);
    for (int j = 0; j < height; j++)
    {
        cur_link_val = get_plaquette_link(col_index,j,3);
        if (cur_link_val == 0)
        {
            if (did_left){
                col_index -=1;
                set_plaquette_link(col_index,j,3,0); 
                set_plaquette_link(col_index,j,2,1); 
                did_left = true;
            }
            else {
                col_index +=1;
                set_plaquette_link(col_index,j,3,0); 
                set_plaquette_link(col_index,j,0,1); 
                did_left = false; 
            }
            //cout << "Wrong column being changed!" << endl;
            //exit(EXIT_FAILURE);
        }
        else
        {

            while (col_index != orig_col_index){
                if (col_index < orig_col_index ){
                    if(get_plaquette_link(col_index,j,3)==1){
                        set_plaquette_link(col_index,j,3,0); 
                        col_index+=1;
                        set_plaquette_link(col_index,j,3,1); 
                        col_index+=1;
                        did_left=true;
                    }
                }
                if (col_index > orig_col_index ){
                    if(get_plaquette_link(col_index,j,3)==1){
                        set_plaquette_link(col_index,j,3,0); 
                        col_index-=1;
                        set_plaquette_link(col_index,j,3,1); 
                        col_index-=1;
                        did_left=false;
                    }
                }
            }


            set_plaquette_link(col_index,j,3,0); 
            if (did_left){
                set_plaquette_link(col_index,j,0,1); 
                did_left = false;
            }
            else{
                set_plaquette_link(col_index,j,2,1); 
                did_left = true;
            }
        }
    }
}

// The stagged dimer configuration
void Lattice::build_init_fully_packed_staggard_dimer()
{


    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < length; j ++)
        {
            // can get the staggered configuration if we only put the link to the right of the vertex
            // we are at when i + j -> even.
            if ((i+j)%2 == 0)
            {
                // Now we are at the bottom left corner of the plaquette to make
                // put the link to the right to 1
                * v_arr[i][j].links[0] = 1; //right link
            }
        }
    }
}

/* in order to work with a fully packed loop model we need the starting configuration to be a fully
 * packed loop. Remember that the fully packed loop model is not ALL spins up. It is the not trivial
 * situation of fully packed loops. 
 * As far as I know it would be difficult to write a program that initializes a random fully packed
 * configuration which is why we want to start with the simplest case and then use random walks
 * which will preserve the fully packed property and randomize it after several applications of a
 * random walk.
 * I'm also not completely sure about this but I think you need to use a self avoiding random walk
 * on the fully packed configuration to ensure it remains fully packed.
 * 
 * So what is the simplest starting config? Its going to be plaquettes where non of them cross... so
 * no vertex has 4 links up. 
 *
 * This function takes:
 * 1st var (var_arr): the reference to the actual array of all vectors.
 */
void Lattice::build_init_fully_packed_loops()
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < length; j ++)
        {
            // we can always get the left corner of the plaquette we want to build by looking only for
            // even i and even j - NOT i+j even.
            if ((i%2 == 0)&(j%2 == 0))
            {
                // Now we are at the bottom left corner of the plaquette to make
                // put the link to the right and the link above to 1
                * v_arr[i][j].links[0] = 1; //right link
                * v_arr[i][j].links[1] = 1; //above link
                
                // now move to the upper right corer of the same plaquette and set the left and bottom
                // link of this vertex to 1
                * v_arr[i+1][j+1].links[2] = 1; //left link
                * v_arr[i+1][j+1].links[3] = 1; //below link
            }
        }
    }
}

void Lattice::build_empty_lat()
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < length; j ++)
        {
            * v_arr[i][j].links[0] = 0; 
            * v_arr[i][j].links[1] = 0; 
            * v_arr[i][j].links[2] = 0; 
            * v_arr[i][j].links[3] = 0; 
        }
    }
}

// same as above but just slightly different to see if it changes any of the data
void Lattice::build_other_init_fully_packed()
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < length; j ++)
        {
            // we can always get the left corner of the plaquette we want to build by looking only for
            // even i and even j - NOT i+j even.
            if ((i%2 == 1)&(j%2 == 0))
            {
                // Now we are at the bottom left corner of the plaquette to make
                // put the link to the right and the link above to 1
                * v_arr[i][j].links[0] = 1; //right link
                * v_arr[i][j].links[1] = 1; //above link
                
                // now move to the upper right corer of the same plaquette and set the left and bottom
                // link of this vertex to 1
                if (i != height-1)
                {
                    * v_arr[i+1][j+1].links[2] = 1; //left link
                    * v_arr[i+1][j+1].links[3] = 1; //below link
                }
                else
                {
                    * v_arr[0][j+1].links[2] = 1; //left link
                    * v_arr[0][j+1].links[3] = 1; //below link
                }
            }
        }
    }
}
