/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#include <bootstrapping.h>
#include <histogram.h>
#include <cmath>
#include <ctime>
#include <QString>

// --------------------------------------------------------------------------
void BootStrap::random(std::vector<bool> &included, size_t n, size_t N){

    if(n >= N){
        //N = N;
        printf(" random_subset: n should be smaller than N!\n");
        return;
    }

    // if i need to include less than half the elements,
    // flag the ones to be included
    if( n < 0.5*N ) {

        included.resize(N, false);      // all elements are excluded by default
        uint cnt = 0;
        while(cnt < n){

            uint idx = (uint) (rand() % N);
            if( included[idx] )
                continue;
            included[idx] = true;
            cnt++;
        }
    }

    // if more than half the elements have to be included
    // flag the ones to be excluded instead
    else {
        n = N - n;

        included.resize(N, true);       // all elements are included by default
        uint cnt = 0;
        while(cnt < n){

            uint idx = (uint) (rand() % N);
            if( !included[idx] )
                continue;
            included[idx] = false;
            cnt++;
        }
    }
}

// --------------------------------------------------------------------------
bool BootStrap::read(QString filename, uint num_bins, float &l2max, float &l2avg, float &l0max, float &l0avg){

    bool verbose = false;

    size_t sz_data = 0;
    size_t num_steps = 0;
    size_t sz_subset = 0;
    uint bin_offset = 0;

    FILE *datafile = fopen(filename.toLatin1().data(), "rb");
    if(!datafile){
        printf(" Unable to open data file : %s\n", filename.toLatin1().data());
        return false;
    }

    if(verbose){
        printf(" Reading bootstrap data file %s...", filename.toLatin1().data());
        fflush(stdout);
    }

    // read metadata
    fread(&sz_data, sizeof(size_t), 1, datafile);
    fread(&bin_offset, sizeof(uint), 1, datafile);
    fread(&num_steps, sizeof(size_t), 1, datafile);
    fread(&sz_subset, sizeof(size_t), 1, datafile);

    int idx = ceil((float)(num_bins)/(float)bin_offset) - 1;

    fseek(datafile, idx * 4 * sizeof(float), SEEK_SET);

    fread(&l2max, sizeof(float), 1, datafile);
    fread(&l2avg, sizeof(float), 1, datafile);
    fread(&l0max, sizeof(float), 1, datafile);
    fread(&l0avg, sizeof(float), 1, datafile);

    // close and done!
    fclose(datafile);
    if(verbose){
        printf(" Done!\n");
    }
    return true;
}

void BootStrap::compute_and_write(QString filename, size_t num_steps, float subset_percent, const QVector<double> &data){

    // ---------------------------------------
    // open a file and write metadata

    FILE *datafile = fopen(filename.toLatin1().data(), "wb");
    if(!datafile){
        printf(" Unable to open data file : %s\n", filename.toLatin1().data());
        return;
    }

    size_t sz_data = data.size();
    size_t sz_subset = subset_percent*sz_data;

    static const uint bin_step = 4;

    printf(" Bootstrapping histogram of %d values: %d steps with %d subset size, and writing to %s...",
           sz_data, num_steps, sz_subset, filename.toLatin1().data());
    fflush(stdout);

    clock_t begin = std::clock();

    fwrite(&sz_data, sizeof(size_t), 1, datafile);
    fwrite(&bin_step, sizeof(uint), 1, datafile);
    fwrite(&num_steps, sizeof(size_t), 1, datafile);
    fwrite(&sz_subset, sizeof(size_t), 1, datafile);

    // ---------------------------------------
    // start bootstrapping

    std::vector<double> L2errors(num_steps, 0);
    std::vector<double> L0errors(num_steps, 0);

    for(int num_bins = bin_step; num_bins <= 180; num_bins+=bin_step){

        // create true histogram
        histogram h;
        h.compute<double>(num_bins, data.data(), sz_data, 0, 180);

        #pragma parallel for
        for(uint bstep = 0; bstep < num_steps; bstep++){

            // resample
            QVector<double> subdata;
            BootStrap::random_subset<double>(data, subdata, sz_subset);

            histogram bsh;
            bsh.compute<double>(num_bins, subdata.data(), subdata.size(), 0, 180);
            L2errors[bstep] = bsh.L2error(h);
            L0errors[bstep] = bsh.L0error(h);

            //printf(" bootstrap %d -- error = %f %f\n", bstep, L2errors[bstep], L0errors[bstep]);
        }

        float l2max = *(std::max_element(L2errors.begin(), L2errors.end()));
        float l0max = *(std::max_element(L0errors.begin(), L0errors.end()));

        long double l2sum = (long double) std::accumulate(L2errors.begin(), L2errors.end(), (long double)(0.0));
        long double l0sum = (long double) std::accumulate(L0errors.begin(), L0errors.end(), (long double)(0.0));

        float l2avg = l2sum / num_steps;
        float l0avg = l0sum / num_steps;


        // write data
        fwrite(&l2max, sizeof(float), 1, datafile);
        fwrite(&l2avg, sizeof(float), 1, datafile);
        fwrite(&l0max, sizeof(float), 1, datafile);
        fwrite(&l0avg, sizeof(float), 1, datafile);
    }

    fclose(datafile);

    clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    printf(" Done! in %f sec.\n", elapsed_secs);
}


void BootStrap::read_forconvergence(uint tscale, uint num_bins, uint subset_percent, QVector<double> &l2, QVector<double> &l0){

    QString filename("bs_");
    filename.append(QString::number(tscale)).append("_")
            .append(QString::number(num_bins)).append("_")
            .append(QString::number(subset_percent)).append(".raw");

    FILE *datafile = fopen(filename.toLatin1().data(), "rb");
    if(!datafile){
        printf(" Unable to open data file : %s\n", filename.toLatin1().data());
        return;
    }

    printf(" Reading files %s...", filename.toLatin1().data());
    fflush(stdout);

    size_t sz_data = 0, sz_subset = 0;
    size_t num_steps;

    fread(&sz_data, sizeof(size_t), 1, datafile);
    fread(&num_steps, sizeof(size_t), 1, datafile);
    fread(&sz_subset, sizeof(size_t), 1, datafile);

    double *l2f = new double[num_steps];
    double *l0f = new double[num_steps];

    fread(l2f, sizeof(double), num_steps, datafile);
    fread(l0f, sizeof(double), num_steps, datafile);

    l2.resize(num_steps);
    l0.resize(num_steps);
    for(uint i = 0; i < num_steps; i++){
        l2[i] = l2f[i];
        l0[i] = l0f[i];
    }

    fclose(datafile);
    printf(" Done!\n");
}

void BootStrap::compute_and_write_forconvergence(uint tscale, uint num_bins, size_t num_steps, uint subset_percent, const QVector<double> &data){

    size_t sz_data = data.size();
    size_t sz_subset = 0.01f*(float)subset_percent*sz_data;

    // ---------------------------------------
    // open a file and write metadata


    QString filename("bs_");
    filename.append(QString::number(tscale)).append("_")
            .append(QString::number(num_bins)).append("_")
            .append(QString::number(subset_percent)).append(".raw");

    FILE *datafile = fopen(filename.toLatin1().data(), "wb");
    if(!datafile){
        printf(" Unable to open data file : %s\n", filename.toLatin1().data());
        return;
    }

    printf(" Bootstrapping histogram of %d values: %d steps with %d subset size, and writing to %s...",
           sz_data, num_steps, sz_subset, filename.toLatin1().data());
    fflush(stdout);

    clock_t begin = std::clock();

    fwrite(&sz_data, sizeof(size_t), 1, datafile);
    fwrite(&num_steps, sizeof(size_t), 1, datafile);
    fwrite(&sz_subset, sizeof(size_t), 1, datafile);

    // ---------------------------------------
    // start bootstrapping

    std::vector<double> L2errors(num_steps, 0);
    std::vector<double> L0errors(num_steps, 0);

    // create true histogram
    histogram h;
    h.compute<double>(num_bins, data.data(), sz_data, 0, 180);

    #pragma parallel for
    for(uint bstep = 0; bstep < num_steps; bstep++){

        // resample
        QVector<double> subdata;
        BootStrap::random_subset<double>(data, subdata, sz_subset);

        histogram bsh;
        bsh.compute<double>(num_bins, subdata.data(), subdata.size(), 0, 180);
        L2errors[bstep] = bsh.L2error(h);
        L0errors[bstep] = bsh.L0error(h);

        //printf(" bootstrap %d -- error = %f %f\n", bstep, L2errors[bstep], L0errors[bstep]);
    }


    // write data
    fwrite(L2errors.data(), sizeof(double), num_steps, datafile);
    fwrite(L0errors.data(), sizeof(double), num_steps, datafile);

    fclose(datafile);

    clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    printf(" Done! in %f sec.\n", elapsed_secs);
}
