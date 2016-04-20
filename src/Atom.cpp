/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#include "Utils.h"
#include "Atom.h"
#include <fstream>
#include <iostream>

using namespace std;

void Atom::write(size_t t0, size_t t1, std::string fname) const {

    if(t0 > pos.size() || t1 > pos.size()){

        printf(" Atom::write(%ud,%ud). size mismatch! %ud\n", t0, t1, pos.size());
        exit(1);
    }

    printf(" Writing trajectory for atom %d to file %s...", this->mid, fname.c_str());
    fflush(stdout);

    std::ofstream outfile(fname.c_str());
    if(!outfile.is_open()){
        std::cerr << " Atom::write -- Could not open file " << fname << std::endl;
        return;
    }

    for(uint i = t0; i < t1; i++){
        outfile << i << "," << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << std::endl;
    }

    outfile.close();
    printf(" Done!\n");
}

// ------------------------------------------------------------------------------
// find if there is a jump across the periodic domain!
void Atom::find_Jumps(bool verbose){

    if(verbose){
        printf(" Finding splits for atom %d...", mid);
        fflush(stdout);
    }

    uint sz = pos.size();

    for(uint t = 1; t < sz; t++){

        if(Utils::dist_sq(pos[t], pos[t-1]) > 0.5){
            trace_endIdx.push_back(t-1);
        }
    }
    trace_endIdx.push_back(sz-1);
    if(verbose){
        printf(" Done! Found %d splits!\n", trace_endIdx.size());
    }

    return;
    if(verbose){
        for(uint i = 0; i < trace_endIdx.size(); i++)
            printf(" split[%d] = %d\n", i, trace_endIdx[i]);
    }
}

// ------------------------------------------------------------------------------
// ------------- compute bonds
// ------------------------------------------------------------------------------
#include "RW_VASP.h"
// call for lithium atoms!
uint Atom::create_bonds(const Atom& otherAtom, VASPmetadata &mdata, double sqthreshold, bool verbose){

    if(mid == otherAtom.mid)
        return 0;

    if(verbose){
        printf(" Checking bonds between Li atom %d and atom %d (using sq threshold = %E)...", mid, otherAtom.mid, sqthreshold);
        fflush(stdout);
    }

    float lattice_const[3] = { mdata.lattice_vectors[0][0], mdata.lattice_vectors[1][1], mdata.lattice_vectors[2][2] };

    for(uint t = 0; t < mdata.num_tsteps; t++){

        float LiPos[3] = {0,0,0};
        float aPos[3] = {0,0,0};

        // convert to cartesian!
        for(uint d = 0; d < 3; d++){
            LiPos[d] = this->pos[t][d] * lattice_const[d];
            aPos[d] = otherAtom.pos[t][d] * lattice_const[d];
        }

        float sqdist = Utils::get_periodicDist_sq(aPos, LiPos, lattice_const);
        if(sqdist < sqthreshold){
            bonded_atoms[t].insert(otherAtom.mid);
            bonded_atoms_list.insert(otherAtom.mid);
        }
    }
    if(verbose){
        printf(" Done!\n");
    }
    return bonded_atoms_list.size();
}

void Atom::compute_angle_correlation(QVector<double> &angles, QVector<double> &times, size_t tau, size_t start_idx, uint sampling){

    // since I need to look two steps ahead, find the last element i can collect the data for
    size_t end_idx = pos.size() - 2*tau;

    for(uint t0 = start_idx; t0 < end_idx; t0 += sampling){

        float ang = Utils::compute_angle(pos[t0], pos[t0+tau], pos[t0+2*tau]);

        // if this a nan
        if(ang != ang){
            printf(" ang(%d,%d,%d) = %f\n", t0, t0+tau, t0+2*tau, ang);
            exit(1);
        }

        //printf("[%d] (%d %d %d) -- %f\n", angles.size(), t0, t0+tau, t0+2*tau, ang);
        angles.push_back(ang);
        times.push_back(t0);
    }
}

/// ------------------------------------------------------------------------------
/// ------------- create tubes -- data format for gle
/// ------------------------------------------------------------------------------


#ifdef USE_GLE
// convert float* positions into gle format array!
void Atom::add_newTubePoints(const uint start_idx, const uint end_idx, bool verbose){

    if(verbose){
        printf("\t Atom::add_newTubePoints(%d, %d)...", start_idx, end_idx);
        fflush(stdout);
    }

    uint tube_len = end_idx - start_idx + 1;

    // gle tube needs an extra point on both ends!
    glePoint *currtube = new glePoint[tube_len+2];

    // fill in internal points!
    for(uint t = 0; t < tube_len; t++){
        for(uint d = 0; d < 3; d++){
            currtube[t+1][d] = pos[start_idx+t][d];
        }
    }

    // add extra points at both ends!
    {
        for(uint d = 0; d < 3; d++){

            currtube[0][d] = 1.2*currtube[1][d] - 0.2*currtube[2][d];
            currtube[tube_len+1][d] = 1.2*currtube[tube_len][d] - 0.2*currtube[tube_len-1][d];
        }
    }

    tubePoints.push_back(currtube);

    if(verbose)
        printf(" Done!\n");
}

void Atom::add_newTubeColors(const QColor col1, const QColor col2, const uint start_idx, const uint end_idx, const uint total_len, bool verbose){

    if(verbose){
        printf("\t Atom::create_tubeColors(col1, col2, %d, %d, %d)...", start_idx, end_idx, total_len);
        fflush(stdout);
    }

    uint tube_len = end_idx - start_idx + 1;

    // gle tube needs an extra point on both ends!
    gleColor* currtube = new gleColor[tube_len+2];

    // fill in internal points!
    for(uint t = 0; t < tube_len; t++){

        float l = (float) (t+start_idx) / (float) total_len;
        QColor c = Utils::lerp(col1, col2, l);

        currtube[t+1][0] = c.redF();
        currtube[t+1][1] = c.greenF();
        currtube[t+1][2] = c.blueF();
    }

    // add extra points at both ends!
    {
        currtube[0][0] = col1.redF();
        currtube[0][1] = col1.greenF();
        currtube[0][2] = col1.blueF();

        currtube[tube_len+1][0] = col2.redF();
        currtube[tube_len+1][1] = col2.greenF();
        currtube[tube_len+1][2] = col2.blueF();
    }

    tubeColors.push_back(currtube);

    if(verbose)
        printf(" Done!\n");
}
#endif

void Atom::create_Tubes(bool verbose){

#ifdef USE_GLE
    if(verbose){
        printf(" Creating tubes for atom %d...", mid);
        fflush(stdout);
    }

    // ----------------------------------------------
    // create traces for this
    uint start_idx = 0;

    for(uint j = 0; j < trace_endIdx.size(); j++){

        uint end_idx = trace_endIdx[j];

        add_newTubePoints(start_idx, end_idx);
        add_newTubeColors(Qt::lightGray, draw_color, start_idx, end_idx, pos.size());

        start_idx = trace_endIdx[j]+1;
    }

    if(verbose){
        printf(" Done! Created %d tubes!\n", tubePoints.size());
    }
#else
    printf(" Cannot create tubes. GLE is not available!\n");
#endif
}





















#if 0
/// ------------------------------------------------------------------------------
/// ------------- compute histograms
/// ------------------------------------------------------------------------------
#include <numeric>
void Atom::compute_angle_correlation(uint hist_type, uint delta, vector<float> &hist, uint num_bins, bool normalize){

    //printf(" compute_angle_correlation(delta = %d, num_bins = %d)... ", delta, num_bins);

    // ----------------------------------
    hist.resize(num_bins, 0.0f);
    float bin_width = (float) num_bins / 180.0f;

    // go over the data
    uint end_idx = pos.size() - 2*delta;    // i need to look two steps ahead!

    for(uint t0 = 0; t0 < end_idx; t0++){

        uint t1 = t0 + delta;
        uint t2 = t1 + delta;

        //if(delta == 5195)
        //    printf(" %d %d %d -- %d\n", t0, t1, t2, pos.size());

        // find the angle
        float ang = Utils::compute_angle(pos[t0], pos[t1], pos[t2]);

        // if this a nan
        if(ang != ang){
            printf(" ang(%d,%d,%d) = %f\n", t0, t1, t2, ang);
            exit(1);
            continue;
        }

        //printf(" %d %f\n", t0, ang);

        // find the corresponding bin
        int bnum = floor(ang * bin_width);

        if(bnum < 0 || bnum >= num_bins){
            printf(" Atom::compute_angle_correlation -- Error -- bnum = %d, but num_bins = %d\n", bnum, num_bins);
            exit(1);
        }

        // update the histogram
        hist[bnum]++;
    }

    if(normalize){
        float sum = std::accumulate(hist.begin(), hist.end(), 0.0f);
        for(uint bnum = 0; bnum < num_bins; bnum++){
            hist[bnum] = hist[bnum] / sum;
        }
    }

    //printf(" Done!\n");
    return;
    //for(uint bnum = 0; bnum < num_bins; bnum++)
      //  printf("bin[%d] = %d\n", bnum, hist[bnum]);

    return;
    /*uint cnt = 0;
    for(uint bnum = 0; bnum < num_bins; bnum++)
        cnt += hist[bnum];

    printf(" Done. Total = %d\n", cnt);*/
}
vector<float> Atom::compute_angle_correlation(uint hist_type, uint delta, uint num_bins, bool normalize){

    vector<float> data;
    //hist.resize(num_bins, 0.0f);
//    /float bin_width = (float) num_bins / 180.0f;

    // go over the data
    uint end_idx = pos.size() - 2*delta;    // i need to look two steps ahead!

    for(uint t0 = 0; t0 < end_idx; t0++){

        uint t1 = t0 + delta;
        uint t2 = t1 + delta;

        //if(delta == 5195)
        //    printf(" %d %d %d -- %d\n", t0, t1, t2, pos.size());

        // find the angle
        float ang = Utils::compute_angle(pos[t0], pos[t1], pos[t2]);


        // if this a nan
        if(ang != ang){
            printf(" ang(%d,%d,%d) = %f\n", t0, t1, t2, ang);
            exit(1);
            continue;
        }

        data.push_back(ang);/*
        //printf(" %d %f\n", t0, ang);

        // find the corresponding bin
        int bnum = floor(ang * bin_width);

        if(bnum < 0 || bnum >= num_bins){
            printf(" Atom::compute_angle_correlation -- Error -- bnum = %d, but num_bins = %d\n", bnum, num_bins);
            exit(1);
        }

        // update the histogram
        hist[bnum]++;*/
    }

    return data;
   //return histogram(data, 0f, M_PI, num_bins, false, true);
/*
    if(normalize){
        float sum = std::accumulate(hist.begin(), hist.end(), 0.0f);
        for(uint bnum = 0; bnum < num_bins; bnum++){
            hist[bnum] = hist[bnum] / sum;
        }
    }
*/
    //printf(" Done!\n");
    //return;
    //for(uint bnum = 0; bnum < num_bins; bnum++)
      //  printf("bin[%d] = %d\n", bnum, hist[bnum]);

    //return;
    /*uint cnt = 0;
    for(uint bnum = 0; bnum < num_bins; bnum++)
        cnt += hist[bnum];

    printf(" Done. Total = %d\n", cnt);*/
}
/*
vector<float> Atom::compute_angle_correlation(size_t tau, size_t start_idx, int sampling){

    vector<float> angles;

    // since I need to look two steps ahead,
        // figure out the last element i can collect the data for
    size_t end_idx = pos.size() - 2*tau;

    //for(uint t0 = start_idx; t0 < end_idx; t0 += sampling){
    for(uint t0 = 0; t0 < end_idx; t0++){

        float ang = Utils::compute_angle(pos[t0], pos[t0+tau], pos[t0+2*tau]);

        // if this a nan
        if(ang != ang){
            printf(" ang(%d,%d,%d) = %f\n", t0, t0+tau, t0+2*tau, ang);
            exit(1);
        }
        angles.push_back(ang);
    }
    return angles;
}*/
#endif

/*void Atom::compute_angle_correlation(statistics_old &angles, size_t tau, size_t start_idx, int sampling){

    // since I need to look two steps ahead,
        // figure out the last element i can collect the data for
    size_t end_idx = pos.size() - 2*tau;

    for(uint t0 = start_idx; t0 < end_idx; t0 += sampling){

        float ang = Utils::compute_angle(pos[t0], pos[t0+tau], pos[t0+2*tau]);

        // if this a nan
        if(ang != ang){
            printf(" ang(%d,%d,%d) = %f\n", t0, t0+tau, t0+2*tau, ang);
            exit(1);
        }
        angles.add_val(ang);
    }
}*/



#if 0
void Atom::compute_2D_hist(uint hist_type, uint min_delta_, uint max_delta_, uint del_delta_, uint num_bins_){

#if 1
    /*if(!hist_2d_angle.empty()){
        printf(" 2D histogram for atom %d exists!\n", mid);
        return;
    }*/

    min_delta = min_delta_;
    max_delta = max_delta_;
    del_delta = del_delta_;
    num_bins = num_bins_;

    uint num_delta = (max_delta-min_delta)/del_delta;

    printf(" Computing 2D histogram for atom %d -- delta (%d,%d), del_delta = %d, num_delta = %d, num_bins = %d...",
                                                            mid, min_delta, max_delta, del_delta, num_delta, num_bins);
    fflush(stdout);

    //]]hist_2d_angle.resize(num_delta);

    //hist_2d_angle_new.resize(num_delta);

    for(uint d = 0; d < num_delta; d++){

        uint delta = min_delta + d*del_delta;
        //printf(" now go for d = %d, size = %d, delta = %d\n", d, hist_2d_angle.size(), delta);

        //this->compute_angle_correlation(hist_2d_angle_new[d], delta, 0, 1);

        statistics_old h;
        this->compute_angle_correlation(h, delta, 0, 1);

        /*vector<float> data = compute_angle_correlation(delta, 0, 1);
        for(int i = 0; i < data.size(); i++)
            h.add_val(data[i]);*/

        h.need_pdf(0, 180, num_bins);

        //hist_2d_angle.push_back(h);

        /*histogram h;
        compute_angle_correlation(hist_type, delta, h, num_bins, true);
        //compute_angle_correlation(hist_type, delta, hist_2d_angle[d], num_bins, true);
        hist_2d_angle.push_back(h);*/

        /*uint cnt = 0;
        for(uint i = 0; i < num_bins; i++)
            cnt+= hist_2d_angle[d][i];
        printf(" count = %d\n", cnt);*/

        //hist_2d_angle.back().print();
    }

    printf(" Done! from here as well!\n");

    return;
    printf(" printing histogram !\n");
/*
    for(uint d = 0; d < num_delta; d++){
        for(uint b = 0; b < hist_2d_angle[d].pdf.size(); b++)
            printf(" %d ", hist_2d_angle[d][b]);
        printf("\n");
    }*/
#endif
}
#endif

