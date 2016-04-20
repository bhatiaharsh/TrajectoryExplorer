/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#ifndef _TRAJECTORY_VIEWER_H
#define _TRAJECTORY_VIEWER_H

#include <QGLViewer/qglviewer.h>

// forward declaration!
class Atom;
class TrajectoryExplorerApp;

/// ===========================================================
/// AtomTrajectoryViewer
/// ===========================================================

class TrajectoryViewer : public QGLViewer{

    TrajectoryExplorerApp* parentApp;
    GLuint sphereList;

    float box_x, box_y, box_z;

public:

    float bbox_dims[3];
    uint live_animation_time;
    TrajectoryViewer(QWidget *parent);

protected :

    virtual void draw();
    virtual void animate();
    virtual void keyPressEvent(QKeyEvent* event);

    virtual void init();
    virtual QString helpString() const;

private:

    void draw_sphere(double R, double NumLatitudes, double NumLongitudes);
    void draw_bbox(const float dims[], QColor col = Qt::black, float width = 1.5);

    void draw_tube(const Atom &atom, uint start_time, uint end_time);
    void draw_linestrip(const Atom &atom, uint start_time, uint end_time, uint step = 1);

    void draw_tet(const std::set<uint> &bonded_atoms);
    void draw_bond(float p[] , float q[]);
    void draw_atomTrace(const Atom &atom, uint start_time, uint end_time);
    void draw_atomPosition(const Atom &atom, const uint time);
};
#endif
