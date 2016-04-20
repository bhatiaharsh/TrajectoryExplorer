/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#include <cmath>
#include <histogram.h>

void histogram::normalize(){
    long double sum = std::accumulate(value.begin(), value.end(), (long double) 0 );
    long double one_over_sum = 1.0 / sum;
    std::transform(value.begin(), value.end(), value.begin(), std::bind1st(std::multiplies<double>(), one_over_sum));
    normalized = true;
}

void histogram::merge(const QVector<histogram> &raw, uint idx0, uint idx1){

    if(idx0 > idx1){
        printf(" ---- cannot merge. %d %d\n", idx0, idx1);
        return;
    }

    uint n = idx1-idx0+1;
    uint num_bins = raw.front().num_bins();

    count.fill(0, num_bins);
    value.fill(0, num_bins);
    L2.fill(0, n);
    L0.fill(0, n);
    indices.first = idx0;
    indices.second = idx1;

    for(uint bin = 0; bin < num_bins; bin++){

        for(uint hidx = idx0; hidx <= idx1; hidx++){
            this->count[bin] += raw[hidx].get_count(bin);
        }

        count[bin] /= (float)n;
        value[bin] = count[bin];
    }

    this->normalize();

    // now compute all error
    for(uint hidx = idx0; hidx <= idx1; hidx++){

        L2[hidx-idx0] = this->L2error(raw[hidx]);
        L0[hidx-idx0] = this->L0error(raw[hidx]);
    }
}

double histogram::L2error(const histogram &h) const {

    if(value.size() != h.num_bins()){
        printf(" histogram::L2error(). Error -- size mismatch (%d != %d)\n", this->value.size(), h.num_bins());
        return -1;
    }

    if(!normalized || !h.normalized){
        printf(" histogram::L2error(). Error -- normalized histograms required! %d %d\n", normalized, h.normalized);
        return -1;
    }

    double L2 = 0;
    for(uint bin = 0; bin < value.size(); bin++){

        double diff = value[bin] - h.get_value(bin);
        L2 += (diff*diff);
    }
    return sqrt(L2) / value.size();
}

double histogram::L0error(const histogram &h) const {

    if(value.size() != h.num_bins()){
        printf(" histogram::L0error(). Error -- size mismatch (%d != %d)\n", this->value.size(), h.num_bins());
        return -1;
    }

    if(!normalized || !h.normalized){
        printf(" histogram::L0error(). Error -- normalized histograms required! %d %d\n", normalized, h.normalized);
        return -1;
    }

    double L0 = 0;
    for(uint bin = 0; bin < value.size(); bin++){

        double diff = fabs(value[bin] - h.get_value(bin));
        if(L0 < diff)
            L0 = diff;
    }
    return L0;
}
