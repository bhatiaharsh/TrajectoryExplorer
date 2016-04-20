/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/


#ifndef _RW_VASP_H_
#define _RW_VASP_H_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

struct VASPmetadata {

    // name of the system
    std::string sysname;

    // information about materials
    std::vector<std::pair<std::string, unsigned int> > materials;
    size_t num_atoms;

    // information about spatial domain
    float scaling_factor;
    float lattice_vectors[3][3];
    size_t grid_dims[3];
    size_t grid_sz;
    std::string coordinate_type;    // direct ?

    // information about time
    size_t num_tsteps;
    double time_step;
    std::string time_unit;
    unsigned int time_factor;

    // other information
    float avg_atomic_vol;
    float temperature;

    VASPmetadata() : sysname("unknown"),
                    num_atoms(-1), scaling_factor(-1), grid_sz(0), coordinate_type("unknown"),
                    num_tsteps(0), time_step(-1), time_unit("unknown"),
                    avg_atomic_vol(-1), temperature(-1) {
        for(unsigned int i = 0; i < 3; i++){
            grid_dims[i] = 0;
            for(unsigned int j = 0; j < 3; j++){
                lattice_vectors[i][j] = -1;
            }
        }
    }

    void print() const;
    void print_time() const;
    void print_mat() const;

    void read_lattice(std::ifstream &infile);
    void read_grid(std::ifstream &infile);
    void read_materials(std::ifstream &infile);

    void write_lattice(std::ofstream &outfile) const;
    void write_materials(std::ofstream &outfile) const;


    std::string find_material_type(unsigned int id) const {

        unsigned int cid = 0;

        for(unsigned int i = 0; i < materials.size(); ++i){

            if( id <= cid+materials.at(i).second ){
                return materials.at(i).first;
            }

            cid += materials.at(i).second;
        }

        return ".";
    }
};

namespace RW_VASP {

    bool read_POSCAR(const std::string &filename, VASPmetadata &mdata, std::vector<float*> &positions);
    bool write_POSCAR(const std::string &filename, const VASPmetadata &mdata, const std::vector<float*> &positions);

    float* read_CHGCAR(const std::string &filename, VASPmetadata &mdata, std::vector<float*> &positions, bool read_pos = true, bool read_values = true);
    bool write_CHGCAR(const std::string &filename, const VASPmetadata &mdata, const std::vector<float*> &positions, bool write_values = true);

    bool read_XDATCAR(const std::string filename, VASPmetadata &mdata, std::vector<float*> &positions);
    bool read_XDATCAR_all(const std::string filename, VASPmetadata &mdata, std::vector<std::vector<float*> > &pos, int read_num = -1);
    void read_XDATCAR_group(const std::string &filename, VASPmetadata &mdata, std::vector<std::vector<float*> > &pos, int num_tsteps);
}

#endif
