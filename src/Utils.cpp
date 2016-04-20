/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#include <fstream>
#include <sstream>
#include <iostream>
#include "Utils.h"

#include <iterator>

std::vector<std::string> Utils::tokenize(std::string line){

    // construct a stream from the string
    std::stringstream linestream_id(line);

    // use stream iterators to copy the stream to the vector as whitespace separated strings
    std::istream_iterator<std::string> it_line(linestream_id);
    std::istream_iterator<std::string> end_line;
    std::vector<std::string> tokens(it_line, end_line);
    return tokens;
}



float Utils::compute_angle(const float v1[], const float v2[]){

    double val = dot(v1, v2) / (magn(v1)*magn(v2));

    // handle numerical instabilities
    if(val >= 1.0)      return 0;
    if(val <= -1.0)     return M_PI;

    // return angle
    return acos(val) * 180.0 / M_PI;
}

float Utils::compute_angle(const float *p0, const float *p1, const float *p2){

    // pivot at p1, and compute displacement vectors!
    float v1[3], v2[3];
    for(unsigned int d = 0; d < 3; d++){

        v1[d] = get_periodicDisplacement(p1[d], p0[d], 1.0);
        v2[d] = get_periodicDisplacement(p2[d], p1[d], 1.0);
    }
    return compute_angle(v1, v2);
}



// --------------------------------------------------------------------------
// color picking


QColor Utils::get_color(uint id, std::string type, uint colortype){

    int h = 0;
    static int s = 0.9f*256.f;
    static int v = 0.6f*256.f;

    if(type.compare("Li") == 0)              h = 320;

    // colortype == 1 : based on specie
    else if(colortype == 1){

        if(type.compare("C") == 0)           h = 60;
        else if(type.compare("H") == 0)      h = 120;
        else if(type.compare("O") == 0)      h = 180;
        else if(type.compare("P") == 0)      h = 240;
        else if(type.compare("F") == 0)      h = 280;
    }

    // colortype == 2 : random
    else if(colortype == 2){

        h = rand() % 360;
    }

    return QColor::fromHsv(h,s,v);
}

/*
QColor Utils::get_color(uint idx){

    if(idx < 0)
        return Qt::black;

    idx = idx%15;
    switch(idx){
        case 0:  return Qt::darkGray;
        case 1:  return Qt::gray;
        case 2:  return Qt::lightGray;
        case 3:  return Qt::red;
        case 4:  return Qt::green;
        case 5:  return Qt::blue;
        case 6:  return Qt::cyan;
        case 7:  return Qt::magenta;
        case 8:  return Qt::yellow;
        case 9:  return Qt::darkRed;
        case 10: return Qt::darkGreen;
        case 11: return Qt::darkBlue;
        case 12: return Qt::darkCyan;
        case 13: return Qt::darkMagenta;
        case 14: return Qt::darkYellow;
    }
}

QString Utils::get_color_tag(QColor col){


    if(col == Qt::black) return (QString("black"));
    if(col == Qt::darkGray) return QString("darkGray");
    if(col == Qt::gray) return QString("gray");
    if(col == Qt::lightGray) return QString("lightGray");
    if(col == Qt::red) return QString("red");
    if(col == Qt::green) return QString("green");
    if(col == Qt::blue) return QString("blue");
    if(col == Qt::cyan) return QString("cyan");
    if(col == Qt::magenta) return QString("magenta");
    if(col == Qt::yellow) return QString("yellow");
    if(col == Qt::darkRed) return QString("darkRed");
    if(col == Qt::darkGreen) return QString("darkGreen");
    if(col == Qt::darkBlue) return QString("darkBlue");
    if(col == Qt::darkCyan) return QString("darkCyan");
    if(col == Qt::darkMagenta) return QString("darkMagenta");
    if(col == Qt::darkYellow) return QString("darkYellow");

    printf("QString Utils::get_color_tag(QColor col) -- invalid color!");
    exit(1);
}
*/

// --------------------------------------------------------------------------
// RW

void Utils::read_atomIDs(const char* filename, std::vector<unsigned int> &ids){

    printf(" Reading atom IDs from file %s...", filename);
    fflush(stdout);

    std::ifstream infile(filename);

    if(!infile.is_open()){
        printf(" Could not open file!\n");
        return;
    }
    std::string line;

    // read the ids
    while(!infile.eof()){
        std::getline(infile, line);
        if(line.empty()){
            break;
        }
        ids.push_back(atoi(line.c_str()));  // store VASP ids!
    }

    printf(" Done! Read %d ids!\n", ids.size());
    infile.close();
}
double* Utils::read_binary(const char *filename, size_t &sz, bool verbose){

    // open the file
    FILE *datafile = fopen(filename, "rb");
    if(!datafile){
        std::cerr << " Unable to open data file : "<<filename<<"\n";
        sz = 0;
        exit(1);
        return 0;
    }

    if(verbose){
        std::cout<<" Reading binary doubles from file "<<filename<<"...";
        fflush(stdout);
    }

    if(sz == 0){
        // find the size of the file
        fseek(datafile, 0, SEEK_END);
        long fsz = ftell(datafile);
        rewind(datafile);

        if(fsz % sizeof(double) != 0){
            std::cerr << "\n\t - Invalid number of doubles in file "<<filename<<". Size of file = "<<fsz<<", Size of double = "<<sizeof(double)<<" -- mod = "<<(fsz%sizeof(double))<<std::endl;
            fclose(datafile);
            return 0;
        }
        sz = fsz / sizeof(double);
    }


    // read the data
    double *values = new double[sz];
    unsigned int rd_sz = fread(values, sizeof(double), sz, datafile);
    if(rd_sz != sz){
        std::cerr << "\n\t - Expected "<<sz<<", but read "<<rd_sz<<" values!\n";
        sz = rd_sz;
    }

    // return
    fclose(datafile);
    if(verbose)
        printf(" Done! Read %'lu doubles!\n", sz);
    return values;
}

void Utils::write_binary(const char *filename, const double *data, size_t sz){

    // open the file
    FILE *datafile = fopen(filename, "wb");
    if(!datafile){
        printf(" Unable to open data file : %s\n", filename);
        exit(1);
    }

    printf(" Writing binary data file %s...", filename);
    fflush(stdout);

    fwrite(data, sizeof(double), sz, datafile);

    fclose(datafile);
    printf(" Done! Wrote %'lu doubles!\n", sz);
}
