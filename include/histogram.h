/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#ifndef _HISTOGRAM_H_
#define _HISTOGRAM_H_

#include <numeric>
#include <cmath>
#include <QVector>

class histogram {

    bool normalized;
    std::pair<int,int> indices;

    QVector<int> count;         // element counts
    QVector<double> value;      // normalized counts
    QVector<double> L2;         // histogram error: L2 for each merged histogram
    QVector<double> L0;         // histogram error: L0 for each merged histogram

public:

    void set_indices(int idx0, int idx1){
        indices.first = idx0;
        indices.second = idx1;
    }
    int get_num_merged() const {    return indices.second-indices.first+1;   }
    std::pair<int,int> get_indices() const {    return indices; }
    //double get_param() const {  return param;   }

    uint num_bins() const {             return count.size(); }
    long int num_samples() const {      return (long int) std::accumulate(count.begin(), count.end(), (long int)(0));   }

    int get_count(uint bin) const {       return count[bin];   }
    double get_value(uint bin) const {    return value[bin];   }

    const QVector<int>& get_counts() const {        return count;    }
    const QVector<double>& get_values() const {     return value;    }
    const QVector<double>& get_L2errors() const {   return L2;    }
    const QVector<double>& get_L0errors() const {   return L0;    }

    double L2max() const {            return *std::max_element(L2.begin(), L2.end());   }
    double L0max() const {            return *std::max_element(L0.begin(), L0.end());   }
    double L2avg() const {
        long double l2sum = (long double) std::accumulate(L2.begin(), L2.end(), (long double)(0.0));
        return (l2sum / L2.size());
    }
    double L0avg() const {
        long double l0sum = (long double) std::accumulate(L0.begin(), L0.end(), (long double)(0.0));
        return (l0sum / L0.size());
    }

    double L2error(const histogram &h) const;
    double L0error(const histogram &h) const;

    void normalize();
    void merge(const QVector<histogram> &raw, uint idx0, uint idx1);

    template <typename dType>
    void compute(size_t num_bins, const dType *data, size_t datasize,
                           dType minVal, dType maxVal, int start_idx = -1, int end_idx = -1, int sampling = 1) {

        count.fill(0, num_bins);
        value.fill(0, num_bins);
        L2.fill(0, 1);       // computed singly. no error!
        L0.fill(0, 1);       // computed singly. no error!
        indices.first = -1;
        indices.second = -1;

        double factor = double (num_bins) / double (maxVal - minVal);

        if(start_idx == -1)     start_idx = 0;
        if(end_idx == -1)       end_idx = datasize;
        if(sampling < 1)        sampling = 1;

        for(int i = start_idx; i < end_idx; i+=sampling){

            int bin = floor( factor * data[i] );

            if(bin >= num_bins)     bin = num_bins-1;
            else if(bin < 0)        bin = 0;

            count[bin] ++;
            value[bin] ++;
        }

        this->normalize();
    }


};
#endif
