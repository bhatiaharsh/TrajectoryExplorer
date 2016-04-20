/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#include <ACWindow.h>


ACWindow::ACWindow(uint atomid, QVector<QVector<double> > const *angles,
                      QVector<QVector<double> > const *times,
                      const QVector<double> &tscales,
                      QString dir_datauq){

    // -------------------------------------
    m_atomid = atomid;
    m_datauq_dir = dir_datauq;
    m_currdata = new ACData(angles, times, tscales);

    // -------------------------------------
    ui.setupUi(this);
    viewer = new ACViewer(this);
    setCentralWidget(viewer);

    // -------------------------------------
    // add color map options, and default to single Hue
    ui.cb_colormap->addItem("Grayscale");
    ui.cb_colormap->addItem("Hot");
    ui.cb_colormap->addItem("Cold");
    ui.cb_colormap->addItem("Night");
    ui.cb_colormap->addItem("Candy");
    ui.cb_colormap->addItem("Geography");
    ui.cb_colormap->addItem("Ion");
    ui.cb_colormap->addItem("Thermal");
    ui.cb_colormap->addItem("Polar");
    ui.cb_colormap->addItem("Spectrum");
    ui.cb_colormap->addItem("Jet");
    ui.cb_colormap->addItem("Hues");
    ui.cb_colormap->addItem("singleHue");

    // -------------------------------------
    ui.gb_color->setChecked(false);
    ui.dsb_mincolor->setDecimals(5);
    ui.dsb_mincolor->setRange(0, 1.0);
    //ui.dsb_mincolor->setSingleStep(dt);
    ui.dsb_mincolor->setValue(0);

    ui.dsb_maxcolor->setDecimals(5);
    ui.dsb_maxcolor->setRange(0, 1.0);
    //ui.dsb_maxcolor->setSingleStep(dt);
    ui.dsb_maxcolor->setValue(1);

    // -------------------------------------
    //ui.rb_agg->setChecked(true);

    // -------------------------------------
    // set minmax in the axis params
    ui.gb_axes->setChecked(false);
    ui.sb_numBins->setRange(1, 180);
    ui.sb_numBins->setValue(180);
    ui.sb_numTscales->setRange(1, m_currdata->num_maxTscales());
    ui.sb_numTscales->setValue(m_currdata->num_maxTscales());

    QCPRange r = m_currdata->range_rawtscale();
    double dt = m_currdata->get_dt_raw();

    ui.dsb_minTscale->setDecimals(5);
    ui.dsb_minTscale->setRange(r.lower, r.upper);
    ui.dsb_minTscale->setSingleStep(dt);
    ui.dsb_minTscale->setValue(r.lower);

    ui.dsb_maxTscale->setDecimals(5);
    ui.dsb_maxTscale->setRange(r.lower, r.upper);
    ui.dsb_maxTscale->setSingleStep(dt);
    ui.dsb_maxTscale->setValue(r.upper);

    // -------------------------------------
    viewer->configure(m_currdata);
    viewer->resize(600,600);

    ui.cb_colormap->setCurrentIndex(12);
    updateTitle();
}

void ACWindow::reconcile_valueRange(QCPRange &rng){

    if(viewer->get_live_show_unnormalized()){
        printf(" ACWindow::reconcile_valueRange does not work if normalized values are not used!\n");
        exit(1);

        // UI needs to know the min max values possible for the entire data
    }

    // ---------------------------------------------------
    // if axes need to be set up manually, read them from UI
    if(this->ui.gb_color->isChecked()){

        rng.lower = ui.dsb_mincolor->text().toDouble();
        rng.upper = ui.dsb_maxcolor->text().toDouble();
    }

    // automatic setting of range
    else {

        ui.dsb_mincolor->setValue(rng.lower);
        ui.dsb_maxcolor->setValue(rng.upper);
    }
}

void ACWindow::reconcile_axesParams(size_t &nTheta, size_t &nTscale, QCPRange &qRngTscale){

    //cb_noAggr
    // ---------------------------------------------------
    // if axes need to be set up manually, read them from UI
    if(this->ui.gb_axes->isChecked()){

        nTheta = ui.sb_numBins->text().toULong();
        nTscale = ui.sb_numTscales->text().toULong();

        qRngTscale.lower = ui.dsb_minTscale->text().toDouble();
        qRngTscale.upper = ui.dsb_maxTscale->text().toDouble();
    }

    // automatic setting of range
    else {

        ui.sb_numBins->setValue(nTheta);
        ui.sb_numTscales->setValue(nTscale);
        ui.dsb_minTscale->setValue(qRngTscale.lower);
        ui.dsb_maxTscale->setValue(qRngTscale.upper);
    }


}


bool ACWindow::recompute_data(size_t nTheta, size_t nTscale, QCPRange qRngTscale){

    //printf("ACWindow::recompute_data(%d %d %f %f)\n", nTheta, nTscale, qRngTscale.lower, qRngTscale.upper);

    reconcile_axesParams(nTheta, nTscale, qRngTscale);

    // ---------------------------------------------------
    // snap to range and write final values to UI
    const QVector<double> &tsc = m_currdata->axis_tscaleraw();

    size_t idx_tscale0 = Utils::get_closestIndex(tsc, qRngTscale.lower);
    size_t idx_tscale1 = Utils::get_closestIndex(tsc, qRngTscale.upper);

    if(idx_tscale0 == idx_tscale1)      idx_tscale1 = idx_tscale0+1;
    if(idx_tscale1 > tsc.size()-1)      idx_tscale1 = tsc.size()-1;

    qRngTscale.lower = tsc[idx_tscale0];
    qRngTscale.upper = tsc[idx_tscale1];

    // write the final values to UI
    ui.dsb_minTscale->setValue(qRngTscale.lower);
    ui.dsb_maxTscale->setValue(qRngTscale.upper);

    // if no aggregation

    if(ui.gb_axes->isChecked() && ui.cb_noAggr->isChecked()){

        nTscale = idx_tscale1-idx_tscale0+1;

        //printf(" no aggregarion : (%f %f) - n = %d [%d %d]\n",
        //       qRngTscale.lower, qRngTscale.upper, nTscale, idx_tscale1, idx_tscale0);
        ui.sb_numTscales->setValue(nTscale);
    }

    // ---------------------------------------------------

    QCPRange qRngTheta (0, 180);
    return m_currdata->compute(nTscale, idx_tscale0, idx_tscale1, qRngTheta, nTheta, false);
}


void ACWindow::keyPressEvent(QKeyEvent *event){
    viewer->keyPressEvent(event);
}
void ACWindow::keyReleaseEvent(QKeyEvent *event){
    viewer->keyReleaseEvent(event);
}
void ACWindow::updateTitle(){

    //QSize plotsize = viewer->get_plotSize();
    QString t = "Atom ";
    t.append(QString::number(m_atomid))
         .append(" | # of tscales = ")
         .append(QString::number(m_currdata->num_maxTscales()))
         /*.append(" | Raw Data = [")
         .append(QString::number(m_currdata->num_maxTscales())).append(" x ")
         .append(QString::number(m_currdata->num_maxBins())).append("]")*/
         .append(" | Resolution = [")
         .append(QString::number(this->viewer->get_plotSize().width())).append(" x ")
         .append(QString::number(this->viewer->get_plotSize().height())).append("]");
         /*.append(" | Vis Data = [")
         .append(QString::number(m_currdata->num_tscales())).append(" x ")
         .append(QString::number(m_currdata->num_bins())).append("]");*/

    this->setWindowTitle(t);
}
