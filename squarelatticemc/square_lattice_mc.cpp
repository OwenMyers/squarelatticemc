// loops.cpp

#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <time.h>
#include <sstream>
#include <iomanip>

//#include "loop_mc_tools.hpp"
#include "model.hpp"
#include "dimer_dimer_direction_vert.hpp"
#include "dimer_horz_struc_fac.hpp"
#include "z3_struc_fac.hpp"
#include "dimer_density.hpp"

using namespace std;

// Takes the file you want to write to and the configuration you want to write
// Single row for each configuration. The links are ordered from the lower left vertex, moving left to right 
// and the odd/even columns represent the horizontal/vertical links pointing right/up from each vertex.
void write_config_to_file(ofstream& config_file,vector < vector <Vertex> >& v_arr)
{
    for (int i=0; i < v_arr.size(); i++)
    {
        for (int j = 0; j<v_arr[0].size(); j++)
        {
            //cout << "in config file write loop" << endl;
            // horizontal links
            config_file << * v_arr[i][j].links[0];
            config_file << " ";

            // vertical links
            config_file << (* v_arr[i][j].links[1]);
            config_file << " ";
        }
    }
    config_file << "\n";

}

int main( int argc, char *argv[])
{
    // where will the data file be written? default (not passed in through arg) is ./
    ostringstream dat_f_dir;

    //default seed
    int seed = 0;
    if (argc != 1)
    { 
        cout << " found argument\n";
        cout << " dat_f_dir= " << dat_f_dir << " \n";
        dat_f_dir << argv[1]; 
        dat_f_dir << "/";
        cout << " dat_f_dir= " << dat_f_dir << " \n";

        if (argc == 3)
        {
            string str_seed = argv[2];

            seed= atoi(str_seed.c_str());
            cout << seed << "\n";
        }
    }
    else
    {
        dat_f_dir << "./";
    }

    O_Random rdm;
    rdm.set_seed(seed);

    int in_length = 8;

    // Make the primary lattice object
    // input in this order (<length>,<height>)
    Lattice lat(in_length,in_length,rdm);
    Lattice br_lat(in_length,in_length,rdm);

    // Tell the estimator which model you are using
    // 1 ---> dimer only model
    // 2 ---> star dimer model
    // 3 ---> star dimer with two simple monomers
    // 4 ---> even parity toric code.
    int model = 2;

    // can be 0, 1 (not 2)
    int increase_winding_by = 2;
    // can be 0, 1, 2
    int increase_windingX_by = 0;
    

    // if we are using a dimer pentamer model then we want to defined the relative scale of the
    // dimer/pentamer terms in the Hamiltonian. This can be done with a single parameter between
    // [0-1]. 
    // 0 -> only pentamers
    // 1 -> only dimers
    double dim_pent_frac = 0.5;

    // N_bins: The number of lines in the data file
    // N_measure: How many measurements to average over per bin
    // N_update: How many updated before a measurement
    int N_bins = 10;
    int N_measure = 5;

    //// This must be 1 for Estimator_Sep
    int N_update = 2*lat.length*lat.height;

    // number of updates between gap calculation bins (in time)
    // good for toric code
    int gap_N_update;
    if (model==4)
    {
        gap_N_update = 1;
    }
    else if (model==2)
    {
        gap_N_update = lat.length*lat.height;
    }
    else 
    {
        gap_N_update = lat.length*lat.height; // this is one unit of MC time (good default)
    }


    // should we write the configuration to file?
    bool write_configs = true;

    // ***********************************
    // *** Directory and file Name *******
    // ***********************************
    // make a director with the following information contained in the name
    // M<model number>_L<lattice size(assume square)_Nb<N_bins>_Nm<N_measure>_S<seed>
    dat_f_dir << "M" << model << "L" << lat.length << "Nb" << N_bins << "Nm" << N_measure << "S" << seed << "_";
    cout << "data file: " << dat_f_dir.str() << endl;

    ofstream sys_inf_file;
    sys_inf_file.open((dat_f_dir.str() + "system_info.txt").c_str());
    sys_inf_file << "gap_N_update: " << gap_N_update << endl;
    if (model==2)
    { sys_inf_file << "dimer to pentamer move fraction: " << dim_pent_frac << endl; }
    else { sys_inf_file << "dimer to pentamer move fraction: NA"<<endl; }
    sys_inf_file << "system size: " << in_length << endl;
    sys_inf_file << "number of bins: " << N_bins << endl;
    sys_inf_file << "measurements per bin: " << N_measure << endl;
    sys_inf_file << "updates between measurements: " << N_update << endl;
    sys_inf_file << "initialize winding number: " << increase_winding_by << endl;
    sys_inf_file << "initialize windingX number: " << increase_windingX_by << endl;
    sys_inf_file << "write the configs: " << write_configs << endl;
        
    sys_inf_file.close();
    
/*****************************************************/

    // Make the model object which is the thing that performs the right updates 
    Model modelObj(model,rdm,dat_f_dir.str(),lat);
    modelObj.build_appropriate_backround_lat(br_lat);

    // Tell the estimator which model you are using
    // 1 ---> dimer only model
    // 2 ---> star dimer model
    if (model == 2)
    {
        modelObj.set_dim_pent_frac(dim_pent_frac);
    }

    // make a vector of estimators
    vector <Estimator*> estimatorVector;
    
    // make estimator
    //Estimator est(lat,dat_f_dir.str());

    DimerDimerDirectionVert verticalDimerEst(dat_f_dir.str(),lat,br_lat);
    DimerHorizontalStructureFactor horizontalDimerStructureFactorEstimator(dat_f_dir.str(),lat,br_lat);
    Z3StructureFactor z3StructureFactorEstimator(dat_f_dir.str(),lat,br_lat);
    DimerDensity dimerDensityEstimator(dat_f_dir.str(),lat,br_lat);
    Estimator *verticalDimerEstPntr = &verticalDimerEst;
    Estimator *horizontalDimerStrucFaEstPntr = &horizontalDimerStructureFactorEstimator;
    Estimator *z3StrucFaEstPntr = &z3StructureFactorEstimator;
    Estimator *dimerDensityEstPntr = &dimerDensityEstimator;

    estimatorVector.push_back(verticalDimerEstPntr);
    estimatorVector.push_back(horizontalDimerStrucFaEstPntr);
    estimatorVector.push_back(z3StrucFaEstPntr);
    estimatorVector.push_back(dimerDensityEstPntr);

    
    if (model==1) {lat.build_init_fully_packed_dimer();}
    else if (model==2) {
        

            
        //build_init_fully_packed_staggered_dimer(v_arr);
        lat.build_init_fully_packed_dimer();
        //lat.build_init_fully_packed_staggered_dimer();
        //lat.build_init_horz_columns();
        // This will change the initialized lattices winding number
        if (increase_winding_by!=0){ lat.increase_winding_num(increase_winding_by); }
        // can only do this increase in x if we have "build_init_fully_packed_dimer()" above
        if (increase_windingX_by == 1){
            lat.increase_windingX_by_one(in_length/2);
        }
        else if (increase_windingX_by == 2){
            lat.increase_windingX_by_one(in_length/2);
            lat.increase_windingX_by_one(0);
        }

    }
    else if (model==3) {
        lat.build_init_fully_packed_staggard_dimer();
    }
    else if (model==4) {
        // you don't have to do this. You can start with an empty lattice
        lat.build_init_fully_packed_loops();
    }
    //pentamer cystal
    //est.make_horz_crystal_lattice();
    
    // file object for writing the configuration file if we want it (the physical link picture).
    // "cf"-> stands for config file
    ofstream cf;
    cf.open("configurations.txt");
    write_config_to_file(cf,lat.v_arr);
    //write_config_to_file(est.cf,est.br_lat.v_arr);

    // only write config if there was a successful update
    int did_it_happen;

    clock_t t = clock();
    clock_t out_t;

    int N_equib = 10000;
    for (int j = 0; j < N_equib; j++)
    {
        for (int k = 0; k < N_update; k++)
        {   
                modelObj.apply_local_update();
        }
    }
    cout << "Done equilibrating" << endl;
    


    for (int i = 0; i < N_bins; i++)
    {
        for (int j = 0; j < N_measure; j++)
        {
            for (int k = 0; k < N_update; k++)
            {
                //est.apply_loop_update();
                modelObj.apply_local_update();
                //if (do_gap_calc){update_count++;}
                //did_it_happen = est.apply_local_update();

                // Takes the file you want to write to and the configuration you want to write
                //if (did_it_happen==4){ write_config_to_file(est.cf,lat.v_arr); est.write_star_locations();}
                //if (did_it_happen >= 2){ write_config_to_file(est.cf,v_arr); est.write_monomer_locations();
                //    cout<<"We moved a monomer"<< endl;
                //}
            }
            //if (write_configs) { write_config_to_file(est.cf,lat.v_arr); est.write_star_locations(); }
            for (int alpha = 0; alpha<estimatorVector.size(); alpha++){
                estimatorVector[alpha]->do_measurement();
            }


        }

        for (int alpha = 0; alpha<estimatorVector.size(); alpha++){
            estimatorVector[alpha]->divide_hist_bins_by(N_measure);
            estimatorVector[alpha]->write();
            estimatorVector[alpha]->clear_hist_bins();
        }
        // divide by is different because each update the acceptance rates are updated
        modelObj.divide_hist_bins_by(N_measure*N_update);
        modelObj.write();
        modelObj.clear_hist_bins();

        // For now the model contains information about the acceptance rates of the different moves
        // and we have to treat it similar to an estimator
        cout << "finished ith, i is: " << i << "\n";
        //cout << "run time: " << ((float)out_t)/CLOCKS_PER_SEC << endl;
        //cout << "time remaining (min): " << ((((float)out_t)/CLOCKS_PER_SEC)/i)*((float)(N_bins-i))/60.0 << endl;

        out_t = clock() - t;

        // break if it has been close (about 28 hrs) to 30 hrs
        //printf( "total time elapsed = %f seconds.\n",((float)out_t)/CLOCKS_PER_SEC);

        //if (((float)out_t)/CLOCKS_PER_SEC >= 90000.0) 
        //{
        //    cout << "about to go over wall time. Stopping at: " << ((float)out_t)/CLOCKS_PER_SEC << endl;
        //    break;
        //}
    }

    for (int alpha = 0; alpha<estimatorVector.size(); alpha++){
        estimatorVector[alpha]->close_wf();
    }
    modelObj.close_wf();

    cout << "Made it to very end" << endl;

    return 0;
}


