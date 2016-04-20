/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#ifndef _ATOM_H
#define _ATOM_H

#include <set>
#include <vector>
#include <QVector>
#include <QColor>

struct VASPmetadata;

#ifdef USE_GLE
#include "GL/gle.h"
typedef gleDouble glePoint[3];
#endif

// ===========================================================
// Atom
// ===========================================================

class Atom {

public:

    // raw data
    uint mid;                       // material id (from simulation)
    std::vector<float*> pos;        // positions at all time-steps

    // the entire trajectory broken into connected "traces"
    std::vector<uint> trace_endIdx;

#ifdef USE_GLE
    std::vector<glePoint*> tubePoints;
    std::vector<gleColor*> tubeColors;
#endif

    // bond detection (populated only for Li)
    std::vector<std::set<uint> > bonded_atoms;
    std::set<uint> bonded_atoms_list;

    // draw color
    QColor draw_color;

    // 2D histogram related data
    QVector<QVector<double> > hist_2d_angle;
    QVector<QVector<double> > hist_2d_times;

    uint min_delta, max_delta, del_delta;
    uint num_bins;

    // --------------------------------------------

    void write(size_t t0, size_t t1, std::string fname) const;

    Atom(uint id, QColor c){
        mid = id;
        draw_color = c;
#if 0
        static int s = 0.9f*256.f;
        static int v = 0.6f*256.f;

        // based on type
        int h = (id <= 189) ? 60 :  // carbon
                (id <= 441) ? 120 : // hydrogen
                (id <= 630) ? 180 : // oxygen
                (id <= 631) ? 240 : // phos
                (id <= 637) ? 280 : // flourine
                320;
        /*
        int h = (id == 638) ? 320 :     // lithium
                (id == 631) ? 120 :     // phosphorus
                rand() % 256;*/

        draw_color = QColor::fromHsv(h,s,v);
#endif
    }

    /*void set_color(QColor c) {
        draw_color = c;
    }*/

    bool operator ==(const Atom &rhs) const {
        return mid == rhs.mid;
    }

    void clear_Jumps_and_Tubes(){

        trace_endIdx.clear();
#ifdef USE_GLE
        for(uint i = 0; i < tubePoints.size(); i++){
            delete tubePoints[i];
            delete tubeColors[i];
        }
        tubePoints.clear();
        tubeColors.clear();
#endif
    }

    // find jumps across the domain
    void find_Jumps(bool verbose = false);

    bool is_jump(uint t) const {
        return std::binary_search(trace_endIdx.begin(), trace_endIdx.end(), t);
    }

    // create tubes
#ifdef USE_GLE
    void add_newTubePoints(const uint start_idx, const uint end_idx, bool verbose = false);
    void add_newTubeColors(const QColor col1, const QColor col2, const uint start_idx, const uint end_idx, const uint total_len, bool verbose = false);
#endif
    void create_Tubes(bool verbose = false);

    // bonds
    uint create_bonds(const Atom& LiAtom, VASPmetadata &mdata, double sqthreshold, bool verbose = false);

    // angle correlation
    void compute_angle_correlation(QVector<double> &angles, QVector<double> &times, size_t tau, size_t start_idx, uint sampling);
};
#endif
