/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/


#ifndef _UTILS_H_
#define _UTILS_H_

#pragma once

#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <QColor>
#include <numeric>

/// ------------------------------------------------------------

namespace Utils {

    // -------------------------------------------------------
    // basic mathematical utilities

    inline bool equals(double a, double b){     return (fabs(a-b) < powf(10,-5));   }

    inline float dist_sq(const float p[3], const float q[3]){
        return ((p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]));
    }

    inline float dist(const float p[3], const float q[3]){  return std::sqrt(dist_sq(p, q));            }

    inline float magn_sq(const float a[3]){                 return (a[0]*a[0] + a[1]*a[1] + a[2]*a[2]); }

    inline float magn(const float a[3]){                    return std::sqrt(magn_sq(a));               }

    inline float lerp (float a, float b, float l){          return (a + (b-a)*l);                       }

    inline float dot(const float a[3], const float b[3]){   return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]); }

    inline float* cross(const float a[3], const float b[3]){
        float *cp = new float[3];
        cp[0] = a[1]*b[2] - a[2]*b[1];
        cp[1] = a[2]*b[0] - a[0]*b[2];
        cp[2] = a[0]*b[1] - a[1]*b[0];
        return cp;
    }

    template <typename T>
    typename QVector<T>::const_iterator get_closestIterator(const QVector<T> &vals, const T& value){

        // upper_bound returns an iterator pointing to the first element that is
        // greater than the value.
        typename QVector<T>::const_iterator geq = std::upper_bound(vals.begin(), vals.end(), value);
        if(geq == vals.begin())     return geq;
        if(geq == vals.end())       return --geq;

        // find if one lower element is closer
        typename QVector<T>::const_iterator prv = geq;
        prv--;

        return (*geq - value) < (value - *prv) ? geq : prv;
    }

    template <typename T>
    typename std::vector<T>::const_iterator get_lowerIterator(const std::vector<T> &vals, const T& value){

        // upper_bound returns an iterator pointing to the first element that is
        // greater than the value.
        typename std::vector<T>::const_iterator geq = std::upper_bound(vals.begin(), vals.end(), value);
        if(geq == vals.begin())     return geq;
        return --geq;
    }

    template <typename T>
    std::size_t get_closestIndex(const QVector<T> &vals, const T& value){
        return std::distance(vals.begin(), get_closestIterator(vals, value));
    }
    template <typename T>
    std::size_t get_lowerIndex(const std::vector<T> &vals, const T& value){
        return std::distance(vals.begin(), get_lowerIterator(vals, value));
    }

    template <typename T>
    inline std::size_t snap_to_axis(T val, const std::vector<T> &values){
        return get_closestIndex(values, val);
    }




    template <typename T>
    inline int snap_to_axis(T val, const QVector<T> &values){

        if(val <= values.front())        return 0;
        if(val >= values.back())         return values.size()-1;

        for(uint i = 0; i < values.size()-1; i++){

            if( val < values[i] || val > values[i+1] )
                continue;

            return (fabs(val-values[i]) < fabs(val-values[i+1])) ? i : i+1;
        }

        printf(" cannot snap %f to vector!\n", val);
        exit(1);
        return -1;
    }

    // snap a value to a grid of uniform resolution

    template <typename T>
    inline int snap_to_uniformAxis(T val, const QVector<T> &values){

        int idx = ((val - values[0]) / (values[1]-values[0]));
        return (idx < 0) ? 0 : (idx > values.size()-1) ? idx = values.size()-1 : idx;
    }
/*
    template <typename T>
    inline int snap_value_nonuniform(T val, const QVector<T> &values){

        if(val < values.front())        return 0;
        if(val > values.back())         return values.size()-1;

        for(uint i = 0; i < values.size()-1; i++){

            if( powf() values[i] < val
            if( values[i] < val && val <= values[i+1])
                return i;
        }

        printf(" cannot snap %f to vector!\n", val);
        exit(1);
        return -1;
    }*/

    // -------------------------------------------------------
    // handling periodic boundary

    // wrap around in [0,1]
    inline float wrap_around(float val){
        return (val < 0) ? (val + 1) : (val > 1) ? (val - 1) : val;
    }

    inline float get_periodicDisplacement(float p, float q, float dim) {
        float d = p-q;
        return (d > 0.5f*dim)  ? (d - dim) :     // q is on the left end, and p is on the right end
               (d < -0.5f*dim) ? (d + dim) :     // p is on the left end, and q is on the right end
                                 d;              // else!
    }

    inline float get_periodicDist_sq(const float p[3], const float q[3], const float dims[3]) {
        float r[3];
        for(unsigned int i = 0; i < 3; i++){
            r[i] = get_periodicDisplacement(p[i], q[i], dims[i]);
        }
        return (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    }

    // -------------------------------------------------------
    // angle correlation

    float compute_angle(const float v1[], const float v2[]);
    float compute_angle(const float *p0, const float *p1, const float *p2);

    // -------------------------------------------------------
    // color related

    inline QColor lerp(QColor c1, QColor c2, float l) {

        float r = lerp(c1.redF(), c2.redF(), l);
        float g = lerp(c1.greenF(), c2.greenF(), l);
        float b = lerp(c1.blueF(), c2.blueF(), l);
        //float a = lerp(c1.alphaF(), c2.alphaF(), l);

        return QColor::fromRgbF(r,g,b);//,a);
    }

    //QColor get_color(uint idx);
    //QString get_color_tag(QColor col);


    QColor get_color(uint id, std::string type, uint colortype);

    // -------------------------------------------------------
    // RW utils

    void read_atomIDs(const char* filename, std::vector<unsigned int> &ids);
    //void read_XDATCAR(const char *filename, std::vector<std::vector<float*> > &pos, VASPmetadata1 &data, int read_num = -1);

    double* read_binary(const char *filename, size_t &sz, bool verbose);
    void write_binary(const char *filename, const double *data, size_t sz);

    std::vector<std::string> tokenize(std::string line);
}

#endif
