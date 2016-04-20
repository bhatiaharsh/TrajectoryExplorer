/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#include "RW_VASP.h"
#include "Utils.h"

#include <iostream>
#include <fstream>

#include <cstdlib>
#include <cmath>

void VASPmetadata::print() const {
    std::cout << "    Lattice vector = (" << lattice_vectors[0][0] << "," << lattice_vectors[1][1] << "," << lattice_vectors[2][2] << ")"
              << "    Grid = " << grid_sz << " ["<< grid_dims[0] <<"x" << grid_dims[1] <<"x" << grid_dims[2] <<"]\n";
    print_mat();
    print_time();
}
void VASPmetadata::print_time() const {
    std::cout << "    # time_steps = " << num_tsteps << ", time_step = " << time_step <<" " << time_unit << std::endl;
}
void VASPmetadata::print_mat() const {
    std::cout << "    # atoms = " << num_atoms << ", # mats = "<< materials.size() << "\n";
    for(unsigned int i = 0; i < materials.size(); i++){
        std::cout << "\t material[" << i << "] = [" << materials[i].first << " - " << materials[i].second << "]\n";
    }
}

void VASPmetadata::read_lattice(std::ifstream &infile) {

    std::string line;

    std::getline(infile, sysname);

    std::getline(infile, line);
    sscanf(line.c_str(),"%f",&scaling_factor);

    for(unsigned int i = 0; i < 3; i++){
        std::getline(infile, line);
        sscanf(line.c_str(),"%f %f %f", &lattice_vectors[i][0], &lattice_vectors[i][1], &lattice_vectors[i][2]);
    }
}
void VASPmetadata::write_lattice(std::ofstream &outfile) const {

    // write lattice
    outfile << sysname << std::endl;
    outfile << "\t" << scaling_factor << std::endl;

    for(unsigned int i = 0; i < 3; i++){
        outfile <<"\t" << lattice_vectors[i][0] << "\t" << lattice_vectors[i][1] << "\t" << lattice_vectors[i][2] << std::endl;
    }
}

void VASPmetadata::read_materials(std::ifstream &infile) {

    std::string line;

    std::getline(infile, line);
    std::vector<std::string> tokens = Utils::tokenize(line);

    num_atoms = 0;
    materials.resize(tokens.size(), std::pair<std::string, unsigned int> ("?", 0));

    // check if the files provided symbols!
    if( std::isalpha( tokens.front().at(0) ) ) {
        for(unsigned int i = 0; i < tokens.size(); i++){
            materials[i].first = tokens[i];
        }

        std::getline(infile, line);
        tokens = Utils::tokenize(line);
    }

    // read the counts
    for(unsigned int i = 0; i < tokens.size(); i++){
        materials[i].second = atoi(tokens[i].c_str());
        num_atoms += materials[i].second;
    }
}
void VASPmetadata::write_materials(std::ofstream &outfile) const {

    // write materials
    for(unsigned int i = 0; i < materials.size(); i++){
        outfile << "\t" << materials[i].first;
    }
    outfile << std::endl;

    for(unsigned int i = 0; i < materials.size(); i++){
        outfile << "\t" << materials[i].second;
    }
    outfile << std::endl;
}

void VASPmetadata::read_grid(std::ifstream &infile){

    std::string line;

    std::getline(infile, line);
    if(line.length() <= 1){
        std::getline(infile, line);
    }

    sscanf(line.c_str(),"%d %d %d", &grid_dims[0], &grid_dims[1], &grid_dims[2]);
    grid_sz = grid_dims[0]*grid_dims[1]*grid_dims[2];
}


void read_positions(std::ifstream &infile, std::vector<float*> &positions, int num) {

    std::string line;
    positions.reserve(num);
    for(unsigned int i = 0; i < num; i++){

        std::getline(infile, line);

        float *read_pos = new float[3];
        sscanf(line.c_str(), "%f %f %f\n", &read_pos[0], &read_pos[1], &read_pos[2]);
        positions.push_back(read_pos);
    }

    /*for(int i = 0; i < positions.size(); i++){
        printf(" pos[%d] = %f %f %f\n", i, positions[i][0], positions[i][1], positions[i][2]);
    }*/
}
void write_positions(std::ofstream &outfile, const std::vector<float*> &positions) {

    for(unsigned int i = 0; i < positions.size(); i++){
    for(unsigned int j = 0; j < 3; j++){
        outfile << "\t" << positions[i][j];
    }
    outfile << std::endl;
    }
}

/** POSCAR file contains the lattice geometry and the ionic positions
        http://cms.mpi.univie.ac.at/vasp/vasp/POSCAR_file.html
*/
bool RW_VASP::read_POSCAR(const std::string &filename, VASPmetadata &mdata, std::vector<float*> &positions){

    std::ifstream infile(filename.c_str());
    if(!infile.is_open()){
        std::cerr << " RW_VASP::read_POSCAR -- Could not open file " << filename << std::endl;
        return false;
    }

    // -----------------------------------------------------
    std::cout << " Reading " << filename << "...";
    fflush(stdout);

    // -----------------------------------------------------
    // POSCAR header contains lattice information and materials list

    mdata.read_lattice(infile);
    mdata.read_materials(infile);

    // -----------------------------------------------------
    // POSCAR data

    std::getline(infile, mdata.coordinate_type);     // read "Direct"
    read_positions(infile, positions, mdata.num_atoms);

    // -----------------------------------------------------
    // done
    infile.close();
    std::cout << " Done!\n";
    mdata.print_mat();
    return true;
}

bool RW_VASP::write_POSCAR(const std::string &filename, const VASPmetadata &mdata, const std::vector<float*> &positions){

    std::ofstream outfile (filename.c_str());
    std::cout << " Writing " << filename << "...";

    mdata.write_lattice(outfile);
    mdata.write_materials(outfile);

    outfile << mdata.coordinate_type << std::endl;
    write_positions(outfile, positions);
    printf(" Done!\n");
    outfile.close();
}


/** CHGCAR, AECCAR, LOCPOT files contains the lattice geometry and the ionic positions
        http://cms.mpi.univie.ac.at/vasp/vasp/CHGCAR_file.html
*/
float* RW_VASP::read_CHGCAR(const std::string &filename, VASPmetadata &mdata, std::vector<float*> &positions, bool read_pos, bool read_values){

    std::ifstream infile(filename.c_str());
    if(!infile.is_open()){
        std::cerr << " RW_VASP::read_CHGCAR -- Could not open file " << filename << std::endl;
        return 0;
    }

    // -----------------------------------------------------
    std::cout << " Reading " << filename << "...";
    fflush(stdout);

    std::string line;

    // -----------------------------------------------------
    // CAR header

    mdata.read_lattice(infile);
    mdata.read_materials(infile);

    // -----------------------------------------------------
    // read positions
    std::getline(infile, mdata.coordinate_type);     // read "Direct"

    if(read_pos){
        read_positions(infile, positions, mdata.num_atoms);
    }
    else {
        //printf(" ignoring positions!\n");
        for(unsigned int i = 0; i < mdata.num_atoms; i++){
            std::getline(infile, line);
        }
    }

    // -----------------------------------------------------
    // read grid dimensions

    mdata.read_grid(infile);

    float *values = 0;

    // read values
    if(read_values){

        size_t idx = 0;

        values = new float[mdata.grid_sz];

        while(idx < mdata.grid_sz){

            std::getline(infile, line);
            if(line.empty())
                break;

            std::vector<std::string> tokens = Utils::tokenize(line);
            for(unsigned int i = 0; i < tokens.size(); i++){
                values[idx++] = atof(tokens[i].c_str());
            }
        }
    }

    infile.close();
    std::cout <<" Done!\n";
    mdata.print();
    return values;
}


/** Read XDATCAR file which contains positions for a single time-step
    http://cms.mpi.univie.ac.at/vasp/vasp/XDATCAR_file.html
 */
bool RW_VASP::read_XDATCAR(const std::string filename, VASPmetadata &mdata, std::vector<float*> &positions){

    std::ifstream infile(filename.c_str());
    if(!infile.is_open()){
        std::cerr << " RW_VASP::read_XDATCAR -- Could not open file " << filename << std::endl;
        return false;
    }

    // -----------------------------------------------------
    std::cout << " Reading " << filename << "...";
    fflush(stdout);

    // -----------------------------------------------------
    // XDATCAR header contains lattice information and materials list

    mdata.read_lattice(infile);
    mdata.read_materials(infile);

    // -----------------------------------------------------
    // read positions
    std::getline(infile, mdata.coordinate_type);     // read "Direct"
    read_positions(infile, positions, mdata.num_atoms);

    infile.close();
    printf(" Done!\n");
    mdata.print();
    return true;
}

/** Read XDATCAR_all file which contains positions for many time-steps
    http://cms.mpi.univie.ac.at/vasp/vasp/XDATCAR_file.html
 */
bool RW_VASP::read_XDATCAR_all(const std::string filename, VASPmetadata &mdata, std::vector<std::vector<float*> > &pos, int read_num){

    std::ifstream infile(filename.c_str());
    if(!infile.is_open()){
        std::cerr << " RW_VASP::read_XDATCAR -- Could not open file " << filename << std::endl;
        return false;
    }

    // -----------------------------------------------------
    std::cout << " Reading " << filename << "...";
    fflush(stdout);

    std::string line;

    // read first line and get number of atoms and time-steps
    int temp;
    std::getline(infile, line);
    sscanf(line.c_str(), "%zd %d %zd\n", &mdata.num_atoms, &temp, &mdata.num_tsteps);

    // read other settings
    std::getline(infile, line);
    {
        std::vector<std::string> tokens = Utils::tokenize(line);
        if(tokens.size() != 5){
            std::cerr << " RW_VASP::read_XDATCAR -- Incorrect metadata!\n";
            return false;
        }

        mdata.avg_atomic_vol = atof(tokens[0].c_str());
        mdata.lattice_vectors[0][0] = atof(tokens[1].c_str());
        mdata.lattice_vectors[1][1] = atof(tokens[2].c_str());
        mdata.lattice_vectors[2][2] = atof(tokens[3].c_str());
        mdata.time_step = atof(tokens[4].c_str());

        if(mdata.time_step < pow(10.0, -12.0)){      mdata.time_factor = 12;   mdata.time_unit = std::string("psec");    }
        else if(mdata.time_step < pow(10.0, -6.0)){  mdata.time_factor = 9;    mdata.time_unit = std::string("nsec");    }
        else if(mdata.time_step < pow(10.0, -9.0)){  mdata.time_factor = 6;    mdata.time_unit = std::string("musec");   }
        else if(mdata.time_step < pow(10.0, -3.0)){  mdata.time_factor = 3;    mdata.time_unit = std::string("msec");    }

        mdata.time_step *= pow(10.0, (double)mdata.time_factor);
    }

    std::getline(infile, line);
    sscanf(line.c_str(), "%f\n", &mdata.temperature);

    std::getline(infile, line); // read comment line
    std::getline(infile, line); // read comment line

    // -----------------------------------------------------
    // allocate memory in vector

    // if i am asked to read a smaller number of values
    if(read_num > 0 && read_num < mdata.num_tsteps)
        mdata.num_tsteps = read_num;

    pos.resize(mdata.num_tsteps);
    for(unsigned int t = 0; t < mdata.num_tsteps; t++){

        // one line is comment (Konfig = #) or <empty>
        std::getline(infile, line);
        read_positions(infile, pos[t], mdata.num_atoms);
    }

    mdata.num_tsteps = pos.size();

    std::cout << " Done!\n";
    mdata.print();
    infile.close();
}

/** Read XDATCAR file containing a known number of time-steps for a single group
 */
void RW_VASP::read_XDATCAR_group(const std::string &filename, VASPmetadata &mdata, std::vector<std::vector<float*> > &pos, int num_tsteps){

    std::ifstream infile(filename.c_str());
    if(!infile.is_open()){
        std::cerr << " RW_VASP::read_XDATCAR -- Could not open file " << filename << std::endl;
        return;// false;
    }

    // -----------------------------------------------------
    std::cout << " Reading " << filename << "...";
    fflush(stdout);

    // -----------------------------------------------------
    // XDATCAR header contains lattice information and materials list

    mdata.read_lattice(infile);
    mdata.read_materials(infile);

    std::string line;

    pos.resize(num_tsteps);
    for(unsigned int tstep = 0; tstep < num_tsteps; tstep++){

        // one line is comment (Direct Configuration = #)
        std::getline(infile, line);
        read_positions(infile, pos[tstep], mdata.num_atoms);
    }
    printf(" Done!\n");
    infile.close();
    return;// true;
}
