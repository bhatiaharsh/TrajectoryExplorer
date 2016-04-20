/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#ifndef _BOOTSTRAPPING_H_
#define _BOOTSTRAPPING_H_

#include <QVector>

namespace BootStrap {

    // -------------------------------------------------------
    // generate a list of n random numbers between 0 and n-1 (duplicates allowed)

    template <typename T>
    inline std::vector<T> random(uint n){

        std::vector<T> rlist(n);
        for(uint i = 0; i < n; i++){
            rlist[i] = (T) (rand() % n);
        }
        return rlist;
    }

    void random(std::vector<bool> &included, size_t n, size_t N);

    template <typename T>
    void random_subset(const QVector<T> &set, QVector<T> &subset, size_t n){

        std::vector<bool> included;
        random(included, n, set.size());

        subset.reserve(n);
        for(uint i = 0; i < included.size(); i++){

            if(included[i]){
                subset.push_back( set[i] );
            }
        }
    }

    // -------------------------------------------------------
    bool read(QString filename, uint num_bins, float &l2max, float &l2avg, float &l0max, float &l0avg);
    void compute_and_write(QString filename, size_t num_steps, float subset_percent, const QVector<double> &data);

    void read_forconvergence(uint tscale, uint num_bins, uint subset_percent, QVector<double> &l2, QVector<double> &l0);
    void compute_and_write_forconvergence(uint tscale, uint num_bins, size_t num_steps, uint subset_percent, const QVector<double> &data);
}

#endif
