/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#ifndef _ACWINDOW_H_
#define _ACWINDOW_H_

#include <qcustomplot.h>
#include "ACViewer.h"
#include "ui_Statistics.h"


class ACWindow : public QMainWindow {

    Q_OBJECT

private:

    uint m_atomid;
    ACData *m_currdata;
    QString m_datauq_dir;


public:
    ACViewer *viewer;

    Ui::StatisticsUI ui;

    uint get_atomId() const {           return m_atomid;    }
    const ACData* get_data() const {    return m_currdata;  }
    QString get_datadir() const {       return m_datauq_dir;}
    bool show_datauq() const {          return ui.cb_datauq->isChecked();   }


    ACWindow(uint atomid, QVector<QVector<double> > const *angles,
                          QVector<QVector<double> > const *times,
                          const QVector<double> &tscales,
                          QString dir_datauq = "");

    void reconcile_axesParams(size_t &nTheta, size_t &nTscale, QCPRange &qRngTscale);
    void reconcile_valueRange(QCPRange &rng);

    bool recompute_data(size_t nTheta, size_t nTscale, QCPRange qRngTscale);

    void updateTitle();

    /*bool show_rawPlots() const {
        return ui.rb_raw->isChecked();
    }*/


    virtual void keyPressEvent(QKeyEvent *event);
    virtual void keyReleaseEvent(QKeyEvent *event);

protected slots:

    void on_dsb_maxcolor_editingFinished() {    viewer->replot();   }
    void on_dsb_mincolor_editingFinished() {    viewer->replot();   }
    void on_dsb_maxTscale_editingFinished() {   viewer->replot();   }
    void on_dsb_minTscale_editingFinished() {   viewer->replot();   }
    void on_sb_numBins_editingFinished() {      viewer->replot();   }
    void on_sb_numTscales_editingFinished() {   viewer->replot();   }
    void on_cb_datauq_toggled() {               viewer->replot();   }
    void on_cb_noAggr_toggled() {               viewer->replot();   }

    void on_gb_axes_toggled(bool e){            viewer->replot();   }

    void on_pb_save_clicked(){                  viewer->saveGraphs();   }

    void on_cb_colormap_currentIndexChanged(int idx){
        viewer->update_transferFunction(idx);
    }
    /*void on_gb_axes_toggled(bool e){
        viewer->set_autoAxes(e);
    }*/
};
#endif
