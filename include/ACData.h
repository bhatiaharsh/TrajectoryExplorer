/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#ifndef _ACDATA_H_
#define _ACDATA_H_

#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iterator>

#include <QVector>

#include <histogram.h>
#include <qcustomplot.h>

#include "Utils.h"

class ACData {

private:
    // raw angle data (will be used to compute histograms when needed)
    // pointers to constant data. do not allow modification!
    QVector<QVector<double> > const *m_angles;
    QVector<QVector<double> > const *m_times;
    QCPRange m_rng_theta, m_rng_tscale;

    // max number of bins (180)
    size_t m_numBins_max;

    // curr number of bins in the displayed histograms
    size_t m_numBins_curr;

    // indices of unmerged histograms that correspond to displayed histograms
    std::pair<int,int> m_tscales_curr;

    QVector<histogram> m_histograms_raw;            // histograms for every time scale in the raw data
    QVector<histogram> m_histograms_curr;

    QVector<double> m_axis_theta, m_axis_tscale;
    QVector<double> m_axis_rawtscale;

public:

    // constructor
    ACData(QVector<QVector<double> > const *angles, QVector<QVector<double> > const *times, const QVector<double> &tscales);

    // use the raw angle data to compute 2D matrix!
    bool compute(size_t numTscales, uint tidx0, uint tidx1, QCPRange rng_theta, size_t numBins, bool verbose);
    bool compute(QCPRange rng_tscale, size_t numTscales, QCPRange rng_theta, size_t numBins, bool verbose);

    bool compute_histograms(size_t numBins, int idx0_tscale, int idx1_tscale, bool verbose);
    bool merge_histograms(size_t num_tscales, int idx0_tscale, int idx1_tscale, bool verbose);

    // accessor functions
    inline size_t num_maxBins() const {                 return m_numBins_max;                           }
    inline size_t num_bins() const {                    return m_numBins_curr;                          }

    inline size_t num_maxTscales() const {              return m_angles->size();                        }
    inline size_t num_tscales() const {                 return m_histograms_curr.size();                }

    inline double get_dt_raw() const {                  return m_axis_rawtscale[1]-m_axis_rawtscale[0]; }
    inline QCPRange range_rawtscale() const {           return m_rng_tscale;                            }
    inline QCPRange range_rawtheta() const {            return m_rng_theta;                             }

    inline QCPRange range_theta() const  {              return QCPRange(m_axis_theta.front(), m_axis_theta.back());     }
    inline QCPRange range_tscale() const {              return QCPRange(m_axis_tscale.front(), m_axis_tscale.back());   }

    inline double count(size_t theta, size_t tscale) const {   return m_histograms_curr[tscale].get_count(theta); }
    inline double value(size_t theta, size_t tscale) const {   return m_histograms_curr[tscale].get_value(theta); }
    inline const histogram& hist(size_t tscale) const { return m_histograms_curr[tscale];                  }
    inline const histogram& hist_raw(double tscale) const {
        int idx = Utils::get_closestIndex(axis_tscaleraw(), tscale);
        return m_histograms_raw[idx];
    }

    inline const QVector<double>& axis_tscaleraw() const {              return m_axis_rawtscale;        }
    inline const QVector<double>& axis_tscale() const {                 return m_axis_tscale;           }
    inline const QVector<double>& axis_theta() const {                  return m_axis_theta;            }

    inline const QVector<double>& times_tscale(size_t tscale) const {   return (*m_times)[tscale];      }
    inline const QVector<double>& angles_tscale(size_t tscale) const {  return (*m_angles)[tscale];     }
    inline const QVector<double>& values_tscale(size_t tscale) const {  return m_histograms_curr[tscale].get_values(); }


    QVector<double> get_tsc_ticks(int num_ticks) const {

        double dt = (m_axis_tscale.back() - m_axis_tscale.front()) / (num_ticks-1);

        QVector<double> tv;
        for(uint i = 0; i < num_ticks; i++){

            int idx = Utils::get_closestIndex(m_axis_rawtscale, dt*i + m_axis_tscale.front());
            tv.push_back( m_axis_rawtscale[idx] );
        }
        return tv;
    }

    QVector<double> get_tht_ticks(int num_ticks) const {

        int dt = (m_axis_theta.back() - m_axis_theta.front()) / (num_ticks-1);

        printf(" dt = %d\n", dt);
        QVector<double> tv;
        for(uint i = 0; i < num_ticks; i++){
            tv.push_back( m_axis_tscale.front() + dt*i );
            printf(" tv[%d] = %f\n", i, tv[i]);
        }
        return tv;
    }

#if 0
    float angles_tscale(float tscale, QVector<double> &angles, QVector<double> &times) const {

        for(uint idx = 0; idx < m_axis_rawtscale.size(); idx++){

            //printf(" matching %f with %f\n", tscale, m_axis_rawtscale[idx]);
            if( fabs(m_axis_rawtscale[idx]-tscale) < powf(10,-6) ){
                angles = angles_tscale(idx);
                times = times_tscale(idx);
              //  printf(" found!\n");
                return this->m_axis_rawtscale[idx];
            }

            if( m_axis_rawtscale[idx] > tscale ){
                idx--;
                angles = angles_tscale(idx);
                times = times_tscale(idx);
               // printf(" found greater!\n");
                return this->m_axis_rawtscale[idx];
            }

        }

        return -1;
    }
#endif
    inline QVector<double> values_theta(size_t theta) const {

        size_t numtscales = num_tscales();
        QVector<double> theta_slice(numtscales);
        for(uint tscale = 0; tscale < numtscales; ++tscale){
            theta_slice[tscale] = m_histograms_curr[tscale].get_value(theta);
        }
        return theta_slice;
    }
};

#endif
