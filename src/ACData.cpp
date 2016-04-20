/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#include <ACViewer.h>
#include <ACData.h>
#include <numeric>
#include <bootstrapping.h>

#include <ctime>

ACData::ACData(QVector<QVector<double> > const *angles, QVector<QVector<double> > const *times,
       const QVector<double> &tscales) :
    m_angles(angles), m_times(times), m_numBins_curr(0), m_numBins_max(180)
{

    m_tscales_curr = std::make_pair(-1, -1);
    m_axis_rawtscale = tscales;
    m_histograms_raw.resize(m_axis_rawtscale.size());      // always have this space!

    m_rng_tscale = QCPRange(m_axis_rawtscale.front(), m_axis_rawtscale.back());
    m_rng_theta = QCPRange(0.0, 180.0);

    printf(" Created ACData for %d time-scale values in range (%.5f, %.5f)\n", num_maxTscales(), range_rawtscale().lower, range_rawtscale().upper);

    // create histograms for the complete data!
    this->compute(num_maxTscales(), 0, num_maxTscales()-1, this->range_rawtheta(), m_numBins_max, false);
}

bool ACData::compute_histograms(size_t num_bins, int idx0_tscale, int idx1_tscale, bool verbose){

    // cannot compute more than max number of bins
    if(num_bins > m_numBins_max)        num_bins = m_numBins_max;

    // if current data already corresponds to the required num bins
    if(num_bins == m_numBins_curr){
        if(verbose){
            printf(" - Already have %d raw histograms for %d bins.\n", m_histograms_raw.size(), m_numBins_curr);
        }
        return false;
    }

    // histograms must be created for all available time-scales
    m_numBins_curr = num_bins;
    size_t num_tscales = this->num_maxTscales();

    if(verbose){
        printf(" - Computing %d histograms for %d bins...", num_tscales , m_numBins_curr);
        fflush(stdout);
    }
    m_histograms_raw.clear();       m_histograms_raw.resize(num_tscales);

    for(uint tscale = 0; tscale < num_tscales; tscale++){

        const QVector<double> &angles = angles_tscale(tscale);
        m_histograms_raw[tscale].set_indices(tscale, tscale);
        m_histograms_raw[tscale].compute(m_numBins_curr, angles.data(), angles.size(), range_rawtheta().lower, range_rawtheta().upper);
    }

    m_axis_theta.clear();
    m_axis_theta.resize(m_numBins_curr);

    double theta_width = (this->range_rawtheta().size())/(double)(this->num_bins());

    // add 0.5 to put the value at the center of the bin
    for(uint theta = 0; theta < this->num_bins(); ++theta){
        m_axis_theta[theta] = (double)(theta+0.5)*theta_width + this->range_rawtheta().lower;
    }

    if(verbose){
        printf(" Done!\n");
    }
    return true;
}

/// including both indices!
bool ACData::merge_histograms(size_t numTscales, int idx0_tscale, int idx1_tscale, bool verbose){

    size_t num_input_tscales = idx1_tscale-idx0_tscale+1;

    // cannot compute for more than the available data
    if(numTscales > num_input_tscales)
        numTscales = num_input_tscales;

    // if curr histograms already correspond to the requested merging
    if(numTscales == m_histograms_curr.size() && idx0_tscale == m_tscales_curr.first && idx1_tscale == m_tscales_curr.second){
        if(verbose){
            printf(" - Already have %d merged histograms corresponding to tscales [%d, %d].\n", m_histograms_curr.size(),
                   idx0_tscale, idx1_tscale);
        }
        return false;
    }

    if(verbose){
        printf(" - Merging %d histograms [%d, %d] into %d...",
               num_input_tscales, idx0_tscale, idx1_tscale, numTscales);
        fflush(stdout);
    }

    m_histograms_curr.clear();
    m_histograms_curr.resize(numTscales);

    // now merge
    double kernel_support = (double) num_input_tscales / (double) numTscales;

    for(uint tscale = 0; tscale < numTscales; ++tscale){

        int tscale_srcb = roundf((float)(tscale) * kernel_support) + idx0_tscale;
        int tscale_srce = roundf((float)(tscale+1) * kernel_support) + idx0_tscale;

        m_histograms_curr[tscale].merge(m_histograms_raw, tscale_srcb, tscale_srce-1);
    }

    if(verbose){
        printf(" Done!\n");
    }

    // verify merging
    // no input histograms should have been left
    for(uint tscale = 0; tscale < numTscales-1; ++tscale){

        if(m_histograms_curr[tscale].get_indices().second+1 == m_histograms_curr[tscale+1].get_indices().first)
            continue;

        printf(" \n\n --> a histogram is missed. %d = [%d %d] and %d = [%d %d]\n",
               tscale, m_histograms_curr[tscale].get_indices().first, m_histograms_curr[tscale].get_indices().second,
               tscale+1, m_histograms_curr[tscale+1].get_indices().first, m_histograms_curr[tscale+1].get_indices().second );
        exit(1);
    }
    if(m_histograms_curr.front().get_indices().first != idx0_tscale){
        printf(" \n\n --> first histogram is missed. [%d %d], but start idx = %d\n",
               m_histograms_curr.front().get_indices().first, m_histograms_curr.front().get_indices().second, idx0_tscale);
        exit(1);
    }
    if(m_histograms_curr.back().get_indices().second != idx1_tscale){
        printf(" \n\n --> last histogram is missed. [%d %d], but end idx = %d\n", m_histograms_curr.size()-1,
               m_histograms_curr.back().get_indices().first, m_histograms_curr.back().get_indices().second, idx1_tscale);
        exit(1);
    }

    m_tscales_curr = std::make_pair(idx0_tscale, idx1_tscale);

    // tscale axis
    m_axis_tscale.clear();      m_axis_tscale.resize(this->num_tscales());

    double min_tsc = m_axis_rawtscale[m_tscales_curr.first];
    double max_tsc = m_axis_rawtscale[m_tscales_curr.second];

    double tscale_width = (max_tsc - min_tsc)/(double)(m_histograms_curr.size()-1);
    for(uint tscale = 0; tscale < num_tscales(); ++tscale){
        m_axis_tscale[tscale] = (double)tscale*tscale_width + min_tsc;
    }

    return true;
}

/// include both indices [tidx0, tidx1]
bool ACData::compute(size_t numTscales, uint tidx0, uint tidx1, QCPRange rng_theta, size_t numBins, bool verbose){

    verbose = false;
    static int count = 0;

    if(verbose){
        printf("\n ==> ACData::compute() %d\n\t === %zu tscales in range [%u %u] (%f %f)\n\t === %zu theta in range (%.1f %.1f)\n", count++,
                numTscales, tidx0, tidx1, m_axis_rawtscale[tidx0], m_axis_rawtscale[tidx1],
                numBins, rng_theta.lower, rng_theta.upper);

        if(m_tscales_curr.first == -1 || m_tscales_curr.first == -2){
            printf(" - Currently contain %d [%d %d] histograms with %d bins\n", m_histograms_curr.size(),
                    m_tscales_curr.first, m_tscales_curr.second, m_numBins_curr);
        }
        else {
            printf(" - Currently contain %d [%d %d] (%f %f) histograms with %d bins\n", m_histograms_curr.size(),
                    m_tscales_curr.first, m_tscales_curr.second, m_axis_rawtscale[m_tscales_curr.first],
                    m_axis_rawtscale[m_tscales_curr.second], m_numBins_curr);
        }
    }

    // ----------------------------------------------------------
    // cannot compute more than max number of bins
    if(numBins > m_numBins_max)        numBins = m_numBins_max;

    if(numBins == m_numBins_curr && numTscales == m_histograms_curr.size() &&
            m_tscales_curr.first == tidx0 && m_tscales_curr.second == tidx1){

        if(verbose){
            printf(" Nothing to do!\n");
        }
        return false;
    }

    // ----------------------------------------------------------
    // 1. compute histograms for all tscales (if needed)
    bool new_histograms1 = compute_histograms(numBins, tidx0, tidx1, verbose);

    // if the histograms have been updated, i need to force another merging
    // for that, overwrite the values of current tscale indices
    if(new_histograms1){
        m_tscales_curr = std::make_pair(-1,-1);
    }

    // 2. merge histograms (if needed)
    bool new_histograms2 = merge_histograms(numTscales, tidx0, tidx1, verbose);


    //printf(" returning from compute %d %d\n", new_histograms1, new_histograms2);
    return (new_histograms1 || new_histograms2);
}
