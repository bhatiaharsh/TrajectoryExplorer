/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#ifndef _ACVIEWER_H_
#define _ACVIEWER_H_

#include "qcustomplot.h"
#include "ACData.h"
#include "ACPlots.h"

class ACWindow;

class ACViewer : public QCustomPlot {

Q_OBJECT

private:

    // how should the mouse events be handled!
    enum Mode { MODE_SELECT = 0, MODE_ZOOM_H = 1, MODE_ZOOM_V = 2, MODE_ZOOM_HV = 3 };

    ACWindow *parent;

    // orientation
    Qt::Orientation m_orientation_theta;    // whether theta should take horizontal axis or vertical
                                                // time-scale will have the opposite orientation

    // setting up a 3x3 grid
    enum AxisRectID { AX_TL = 0, AX_TC = 1, AX_TR = 2, AX_CL = 3, AX_CC = 4, AX_CR = 5, AX_BL = 6, AX_BC = 7, AX_BR = 8 };
    std::vector<QCPAxisRect*> m_axisRects;

    // central 2D plot
    QCPColorMap *m_cmap_2d;             // the 2D matrix
    QCPColorGradient *m_colorGradient;  // transfer function
    QCPColorScale *m_colorScale;        // color scale (legend)

    // five 1D graphs
    QCPGraph *m_grph_thtSlice_hist, *m_grph_tscSlice_hist, *m_grph_tscSlice_angles;
    QCPGraph *m_grph_tscError_data, *m_grph_tscError_vis;

    // legends for everyone
    QCPLegend *m_lgnd_thtSlice_hist, *m_lgnd_tscSlice_hist, *m_lgnd_tscSlice_angles;
    QCPLegend *m_lgnd_tscError;

    uint min_pixel_per_value;
    bool live_show_unnormalized;        // whether to map counts or pdf to the 2d image

    // -----------------------------------------------
public:

    PlotSelector *m_pSelector;

    ACViewer(ACWindow *parent);


    inline QSize get_plotSize() const {     return m_axisRects[AX_CC]->size();  }


    void update_transferFunction(int idx);

    void configure(const ACData *data){

        this->m_axisRects[AX_CC]->axis(QCPAxis::atBottom)->setRange(data->range_rawtscale());
        this->m_axisRects[AX_CC]->axis(QCPAxis::atLeft)->setRange(data->range_rawtheta());

        float tsc = data->axis_tscale()[data->axis_tscale().size()/2];
        float tht = data->axis_theta()[data->axis_theta().size()/2];
        m_pSelector->setCoords( tsc, tht );
        m_pSelector->setVisible(false);
        //m_pSelector->print();
    }

    void updateSelector(float tscale, float theta){
        m_pSelector->setCoords( tscale, theta );
        replot();
    }
    void updateSelector(float tscale){
        m_pSelector->setTscale(tscale);
        replot();
    }

    /*void set_autoAxes(bool e){
        printf(" set_autoAxes(%d)\n", e);
    }*/

    bool get_live_show_unnormalized() const {    return live_show_unnormalized;  }

    void saveGraphs();
    void recolor();

    void update_interactionMode(Mode mode);     // interaction mode

private:

    // --------------------------------------------------
    // setup functions that are called only once
    void setup_layout();
    void setup_graphStyling();

    // --------------------------------------------------
    // setup functions for dynamic changes

    void update_orientation();                  // depending upon the orientation, decide where to put the graphs


    // --------------------------------------------------
    bool compute_axesParams(size_t &nTheta, size_t &nTscale, QCPRange &rngTscale);

    // --------------------------------------------------
    // accessor functions
    inline QCPAxis* get_keyAxis(AxisRectID rectId){    return m_axisRects[rectId]->axis(get_keyAxisType(rectId));      }
    inline QCPAxis* get_valueAxis(AxisRectID rectId){  return m_axisRects[rectId]->axis(get_valueAxisType(rectId));    }

    QCPAxis::AxisType get_valueAxisType(AxisRectID rectId) const {

        switch(rectId){
            case AX_TC:     return QCPAxis::atLeft;
            case AX_BC:     return QCPAxis::atLeft;
            case AX_CL:     return QCPAxis::atBottom;
            case AX_CR:     return QCPAxis::atBottom;
            case AX_CC:     return (m_orientation_theta == Qt::Vertical) ? QCPAxis::atLeft : QCPAxis::atBottom;
        }

        printf(" ACViewer::get_valueAxisType -- Cannot determine value axis!\n");
        return QCPAxis::atLeft;
    }
    QCPAxis::AxisType get_keyAxisType(AxisRectID rectId) const {

        switch(rectId){
            case AX_TC:     return QCPAxis::atBottom;
            case AX_BC:     return QCPAxis::atTop;
            case AX_CL:     return QCPAxis::atRight;
            case AX_CR:     return QCPAxis::atLeft;
            case AX_CC:     return (m_orientation_theta == Qt::Vertical) ? QCPAxis::atBottom : QCPAxis::atLeft;
        }

        printf(" ACViewer::get_keyAxisType -- Cannot determine key axis!\n");
        return QCPAxis::atBottom;
    }

    // --------------------------------------------------
    // create plots using m_currdata
    void create2DPlot(bool verbose);
    void createGraphs(bool verbose);

    void tsc_hist(float tscale);
    void tsc_hist(float tscale, uint tscale_idx);
    void tht_hist(float theta, uint theta_idx);
    void tsc_angles(float tscale);
    void errors();

    // convegence graph
    void convergence();

    // --------------------------------------------------
    // save graphs

    static void saveGraphAsPdf(QCPGraph *g1, QCPGraph *g2, QString xlabel, QString ylabel, bool show_legend,
                               QString filename, int X, int Y);

    // static function that works on an axis
    static void setTicks(QCPAxis *axis, bool visible) {
        axis->setTickLabels(visible);
        axis->setTicks(visible);
    }

public:
    virtual void mousePressEvent(QMouseEvent *event);
    virtual void mouseReleaseEvent(QMouseEvent *event);
    virtual void mouseMoveEvent(QMouseEvent *event);

    virtual void keyPressEvent(QKeyEvent *event);
    virtual void keyReleaseEvent(QKeyEvent *event);

signals:
    void signal_selectorMoved(uint atomid, float tscale, float theta);

protected slots:
    void replotEvent();
};
#endif
