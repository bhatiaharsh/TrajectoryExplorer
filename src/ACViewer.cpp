/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#include <numeric>
#include <QFont>

#include "ACViewer.h"
#include "ACData.h"
#include "ACWindow.h"
#include "bootstrapping.h"

// ----------------------------------------------------------------------------
ACViewer::ACViewer(ACWindow *_parent) : QCustomPlot((QWidget*)_parent) {

    parent = _parent;
    m_orientation_theta = Qt::Vertical;         // default!
    min_pixel_per_value = 3;
    live_show_unnormalized = false;

    setup_layout();
    setup_graphStyling();

    update_orientation();
    update_interactionMode(MODE_ZOOM_HV);

    this->setMinimumSize(400, 400);
}

// ----------------------------------------------------------------------------
void ACViewer::setup_layout(){

    //  ----------------------------------------------------
    m_axisRects.resize(9, 0);

    m_colorScale = 0;
    m_cmap_2d = 0;
    m_pSelector = 0;

    m_grph_thtSlice_hist = 0;
    m_grph_tscSlice_hist = 0;
    m_grph_tscSlice_angles = 0;
    m_grph_tscError_data = 0;
    m_grph_tscError_vis = 0;

    //  ----------------------------------------------------
    // set up a 3x3 layout and create all axis rectangles

    const AxisRectID graphRects[4] = {AX_CL, AX_CR, AX_TC, AX_BC};

    // only center rectangle should be fully set up
    m_axisRects[AX_CC] = new QCPAxisRect(this, true);
    for(uint i = 0; i < 4; i++){
        m_axisRects[graphRects[i]] = new QCPAxisRect(this, false);
    }

    this->plotLayout()->clear();
    this->plotLayout()->addElement(0, 1, m_axisRects[AX_TC]);
    this->plotLayout()->addElement(1, 0, m_axisRects[AX_CL]);
    this->plotLayout()->addElement(1, 1, m_axisRects[AX_CC]);
    this->plotLayout()->addElement(1, 2, m_axisRects[AX_CR]);
    this->plotLayout()->addElement(2, 1, m_axisRects[AX_BC]);

    this->plotLayout()->setRowStretchFactor(0, 0.25);
    this->plotLayout()->setRowStretchFactor(2, 0.25);
    this->plotLayout()->setColumnStretchFactor(0, 0.25);
    this->plotLayout()->setColumnStretchFactor(2, 0.25);

    //  ----------------------------------------------------
    // synchronize the layout margins

    QCPMarginGroup *mg_bot = new QCPMarginGroup(this);
    QCPMarginGroup *mg_lft = new QCPMarginGroup(this);

    m_axisRects[AX_CC]->setMarginGroup(QCP::msLeft, mg_lft);
    m_axisRects[AX_TC]->setMarginGroup(QCP::msLeft, mg_lft);
    m_axisRects[AX_BC]->setMarginGroup(QCP::msLeft, mg_lft);

    m_axisRects[AX_CC]->setMarginGroup(QCP::msBottom, mg_bot);
    m_axisRects[AX_CL]->setMarginGroup(QCP::msBottom, mg_bot);
    m_axisRects[AX_CR]->setMarginGroup(QCP::msBottom, mg_bot);

    static bool show_scale = 1;//true;//false;
    if(show_scale){
        QCPMarginGroup *mg_top = new QCPMarginGroup(this);
        QCPMarginGroup *mg_rgt = new QCPMarginGroup(this);

        // add color scale
        m_colorScale = new QCPColorScale(this);
        this->plotLayout()->addElement(0, 2, m_colorScale); // add it to the right of the main axis rect*/

        m_colorScale->setMarginGroup(QCP::msLeft, mg_rgt);
        m_axisRects[AX_CR]->setMarginGroup(QCP::msLeft, mg_rgt);

        m_colorScale->setMarginGroup(QCP::msBottom, mg_top);
        m_axisRects[AX_TC]->setMarginGroup(QCP::msBottom, mg_top);
    }

    // ----------------------------------------------------
    // create and configure the axes

    m_axisRects[AX_CC]->setupFullAxesBox(true);

    for(uint i = 0; i < 4; i++){

        AxisRectID thisRect = graphRects[i];

        // add key and value axes
        // the value axis for all other plots should show only 2 ticks
        m_axisRects[thisRect]->addAxes(get_keyAxisType(thisRect) | get_valueAxisType(thisRect));
        m_axisRects[thisRect]->axis(get_valueAxisType(thisRect))->setAutoTickCount(2);
    }

    // the value axis for left and bottom plots need to be reversed
    m_axisRects[AX_CL]->axis(get_valueAxisType(AX_CL))->setRangeReversed(true);
    m_axisRects[AX_BC]->axis(get_valueAxisType(AX_BC))->setRangeReversed(true);

    // center plot does not show any ticks
    m_axisRects[AX_CC]->axis(QCPAxis::atTop)->setVisible(true);
    m_axisRects[AX_CC]->axis(QCPAxis::atRight)->setVisible(true);
    m_axisRects[AX_CC]->axis(QCPAxis::atBottom)->setVisible(true);
    m_axisRects[AX_CC]->axis(QCPAxis::atLeft)->setVisible(true);

    setTicks(m_axisRects[AX_CC]->axis(QCPAxis::atTop), false);
    setTicks(m_axisRects[AX_CC]->axis(QCPAxis::atRight), false);
    setTicks(m_axisRects[AX_CC]->axis(QCPAxis::atBottom), false);
    setTicks(m_axisRects[AX_CC]->axis(QCPAxis::atLeft), false);

    //  ----------------------------------------------------
    // move newly created axes on "axes" layer and grids on "grid" layer:

    foreach (QCPAxisRect *rect, this->axisRects()){
    foreach (QCPAxis *axis, rect->axes()){
        axis->setLayer("axes");
        axis->grid()->setLayer("grid");
    }
    }

    //  ----------------------------------------------------
    // create the graphs

    // -- color map at the center
    m_cmap_2d = new QCPColorMap(get_keyAxis(AX_CC), get_valueAxis(AX_CC));
    this->addPlottable(m_cmap_2d);

    m_pSelector = new PlotSelector(this, 15);
    m_pSelector->setAxes(get_keyAxis(AX_CC), get_valueAxis(AX_CC));
    m_pSelector->setAxisRect(m_axisRects[AX_CC]);

    // --
    m_grph_thtSlice_hist = this->addGraph(get_keyAxis(AX_TC), get_valueAxis(AX_TC));
    m_lgnd_thtSlice_hist = new QCPLegend;
    m_axisRects[AX_TC]->insetLayout()->addElement(m_lgnd_thtSlice_hist, Qt::AlignTop|Qt::AlignLeft);
    m_lgnd_thtSlice_hist->setLayer("legend");
    m_lgnd_thtSlice_hist->addItem(new QCPPlottableLegendItem(m_lgnd_thtSlice_hist, m_grph_thtSlice_hist));

    // --
    m_grph_tscSlice_hist = this->addGraph(get_keyAxis(AX_CL), get_valueAxis(AX_CL));
    m_lgnd_tscSlice_hist = new QCPLegend;
    m_axisRects[AX_CL]->insetLayout()->addElement(m_lgnd_tscSlice_hist, Qt::AlignTop|Qt::AlignLeft);
    m_lgnd_tscSlice_hist->setLayer("legend");
    m_lgnd_tscSlice_hist->addItem(new QCPPlottableLegendItem(m_lgnd_tscSlice_hist, m_grph_tscSlice_hist));

    // --
    m_grph_tscSlice_angles = this->addGraph(get_keyAxis(AX_CR), get_valueAxis(AX_CR));
    m_lgnd_tscSlice_angles = new QCPLegend;
    m_axisRects[AX_CR]->insetLayout()->addElement(m_lgnd_tscSlice_angles, Qt::AlignTop|Qt::AlignLeft);
    m_lgnd_tscSlice_angles->setLayer("legend");
    m_lgnd_tscSlice_angles->addItem(new QCPPlottableLegendItem(m_lgnd_tscSlice_angles, m_grph_tscSlice_angles));

    // --
    m_grph_tscError_data = this->addGraph(get_keyAxis(AX_BC), get_valueAxis(AX_BC));
    m_grph_tscError_vis = this->addGraph(get_keyAxis(AX_BC), get_valueAxis(AX_BC));
    m_lgnd_tscError = new QCPLegend;
    m_axisRects[AX_BC]->insetLayout()->addElement(m_lgnd_tscError, Qt::AlignBottom|Qt::AlignRight);
    m_lgnd_tscError->setLayer("legend");
    m_lgnd_tscError->addItem(new QCPPlottableLegendItem(m_lgnd_tscError, m_grph_tscError_data));
    m_lgnd_tscError->addItem(new QCPPlottableLegendItem(m_lgnd_tscError, m_grph_tscError_vis));

    //  ----------------------------------------------------
    // synchronize the axis zoom ranges
    connect( this->get_keyAxis(AX_CC), SIGNAL(rangeChanged(QCPRange)), this->get_keyAxis(AX_BC), SLOT(setRange(QCPRange)) );
    connect( this->get_valueAxis(AX_CC), SIGNAL(rangeChanged(QCPRange)), this->get_keyAxis(AX_CL), SLOT(setRange(QCPRange)) );

    // all updating will be done before replot
    connect(this, SIGNAL(beforeReplot()), this, SLOT(replotEvent()));
}
void ACViewer::setup_graphStyling(){

    // red histogram
    if(m_grph_tscSlice_hist != 0){

        m_grph_tscSlice_hist->setLineStyle(QCPGraph::lsStepCenter);
        m_grph_tscSlice_hist->setPen(QPen(QColor(191, 92, 1), 2));
        m_grph_tscSlice_hist->setBrush(QColor(191, 92, 1, 50));
    }

    // blue histogram
    if(m_grph_thtSlice_hist != 0){

        m_grph_thtSlice_hist->setLineStyle(QCPGraph::lsStepCenter);
        m_grph_thtSlice_hist->setPen(QPen(QColor(1, 92, 191), 2));
        m_grph_thtSlice_hist->setBrush(QColor(1, 92, 191, 50));
    }

    // black curve with gray fill
    if(m_grph_tscSlice_angles != 0){

        m_grph_tscSlice_angles->setLineStyle(QCPGraph::lsLine);
        m_grph_tscSlice_angles->setPen(QPen(QColor("#222222"), 1.25));
        m_grph_tscSlice_angles->setBrush(QColor(110, 110, 110, 30));
    }

    if(m_grph_tscError_data != 0 && m_grph_tscError_vis != 0){

        m_grph_tscError_data->setName("Data Uncertainty");
        m_grph_tscError_data->setLineStyle(QCPGraph::lsStepCenter);
        m_grph_tscError_data->setPen(QPen(QColor(110, 170, 110), 1.25));
        m_grph_tscError_data->setBrush(QColor(110, 170, 110, 50));

        m_grph_tscError_vis->setName("Visualization Uncertainty");
        m_grph_tscError_vis->setLineStyle(QCPGraph::lsStepCenter);
        m_grph_tscError_vis->setPen(QPen(QColor(150, 0, 222), 1.25));
        m_grph_tscError_vis->setBrush(QColor(150, 0, 222, 70));
    }

    if(m_cmap_2d != 0){

        m_colorGradient = new QCPColorGradient();
        m_cmap_2d->setInterpolate(false);
        m_cmap_2d->setTightBoundary(true);

        //set_transferFunction(parent->ui.cb_colormap->currentIndex());
    }
}

// ----------------------------------------------------------------------------

void ACViewer::update_interactionMode(Mode mode){

    switch(mode) {

        case MODE_SELECT:
                    this->setInteractions(QCP::iSelectPlottables);
                    this->setInteraction(QCP::iRangeDrag, false);
                    this->setInteraction(QCP::iRangeZoom, false);
                    m_pSelector->setVisible(true);
                    break;

        case MODE_ZOOM_HV:
                    this->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
                    //this->m_axisRects[AX_CC]->setRangeDrag(Qt::Horizontal | Qt::Vertical);
                    //this->m_axisRects[AX_CC]->setRangeZoom(Qt::Horizontal | Qt::Vertical);
                    m_pSelector->setVisible(false);
                    //break;

        case MODE_ZOOM_V:
                    this->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
                    //this->m_axisRects[AX_CC]->setRangeDrag(Qt::Vertical);
                    //this->m_axisRects[AX_CC]->setRangeZoom(Qt::Vertical);
                    m_pSelector->setVisible(false);
                    //break;

        case MODE_ZOOM_H:
                    this->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
                    this->m_axisRects[AX_CC]->setRangeDrag(Qt::Horizontal);
                    this->m_axisRects[AX_CC]->setRangeZoom(Qt::Horizontal);
                    m_pSelector->setVisible(false);
                    break;
    }
    //printf(" zoom - %d\n", (int)this->m_axisRects[AX_CC]->rangeZoom());
}
void ACViewer::update_orientation() {

    get_keyAxis(AX_CC)->setLabel("Time-scale");
    get_valueAxis(AX_CC)->setLabel("Angle");

    if(m_orientation_theta == Qt::Vertical) {

        // histogram slices stay on left and bottom
        m_grph_thtSlice_hist->setKeyAxis(get_keyAxis(AX_BC));
        m_grph_thtSlice_hist->setValueAxis(get_valueAxis(AX_BC));
        m_axisRects[AX_BC]->insetLayout()->addElement(m_lgnd_thtSlice_hist, Qt::AlignBottom|Qt::AlignLeft);

        m_grph_tscSlice_hist->setKeyAxis(get_keyAxis(AX_CL));
        m_grph_tscSlice_hist->setValueAxis(get_valueAxis(AX_CL));
        m_axisRects[AX_CL]->insetLayout()->addElement(m_lgnd_tscSlice_hist, Qt::AlignTop|Qt::AlignLeft);

        // errors and angles on right and top
        m_grph_tscSlice_angles->setKeyAxis(get_keyAxis(AX_CR));
        m_grph_tscSlice_angles->setValueAxis(get_valueAxis(AX_CR));
        m_axisRects[AX_CR]->insetLayout()->addElement(m_lgnd_tscSlice_angles, Qt::AlignTop|Qt::AlignLeft);

        m_grph_tscError_data->setKeyAxis(get_keyAxis(AX_TC));
        m_grph_tscError_data->setValueAxis(get_valueAxis(AX_TC));
        m_grph_tscError_vis->setKeyAxis(get_keyAxis(AX_TC));
        m_grph_tscError_vis->setValueAxis(get_valueAxis(AX_TC));
        m_axisRects[AX_TC]->insetLayout()->addElement(m_lgnd_tscError, Qt::AlignTop|Qt::AlignLeft);

        connect( this->get_keyAxis(AX_CC), SIGNAL(rangeChanged(QCPRange)), this->get_keyAxis(AX_TC), SLOT(setRange(QCPRange)) );
        disconnect( this->get_valueAxis(AX_CC), SIGNAL(rangeChanged(QCPRange)), this->get_keyAxis(AX_CR), SLOT(setRange(QCPRange)) );
    }
    else {

        // histogram slices stay on left anf bottom
        m_grph_thtSlice_hist->setKeyAxis(get_keyAxis(AX_CL));
        m_grph_thtSlice_hist->setValueAxis(get_valueAxis(AX_CL));
        m_axisRects[AX_CL]->insetLayout()->addElement(m_lgnd_thtSlice_hist, Qt::AlignTop|Qt::AlignLeft);

        m_grph_tscSlice_hist->setKeyAxis(get_keyAxis(AX_BC));
        m_grph_tscSlice_hist->setValueAxis(get_valueAxis(AX_BC));
        m_axisRects[AX_BC]->insetLayout()->addElement(m_lgnd_tscSlice_hist, Qt::AlignBottom|Qt::AlignLeft);

        // errors and angles on right and top
        m_grph_tscSlice_angles->setKeyAxis(get_keyAxis(AX_TC));
        m_grph_tscSlice_angles->setValueAxis(get_valueAxis(AX_TC));
        m_axisRects[AX_TC]->insetLayout()->addElement(m_lgnd_tscSlice_angles, Qt::AlignTop|Qt::AlignRight);

        m_grph_tscError_data->setKeyAxis(get_keyAxis(AX_CR));
        m_grph_tscError_data->setValueAxis(get_valueAxis(AX_CR));
        m_grph_tscError_vis->setKeyAxis(get_keyAxis(AX_CR));
        m_grph_tscError_vis->setValueAxis(get_valueAxis(AX_CR));
        m_axisRects[AX_CR]->insetLayout()->addElement(m_lgnd_tscError, Qt::AlignTop|Qt::AlignLeft);

        connect( this->get_keyAxis(AX_CC), SIGNAL(rangeChanged(QCPRange)), this->get_keyAxis(AX_CR), SLOT(setRange(QCPRange)) );
        disconnect( this->get_valueAxis(AX_CC), SIGNAL(rangeChanged(QCPRange)), this->get_keyAxis(AX_TC), SLOT(setRange(QCPRange)) );
    }
}
void ACViewer::update_transferFunction(int idx){

    if(idx < 0 || idx > 12){
        printf(" well, what happened? %d\n", idx);
    }
    //  QCP presets
    if(idx < 12){
        m_colorGradient->loadPreset( QCPColorGradient::GradientPreset(idx));
    }
    else {
        m_colorGradient->clearColorStops();
        m_colorGradient->setColorStopAt(0.0, QColor(240,240,240,240));
        m_colorGradient->setColorStopAt(1.0, QColor("#8070B8"));
    }

    if(m_colorScale){
        m_colorScale->setGradient(*m_colorGradient);
        m_colorScale->setType(QCPAxis::atRight); // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)
        m_cmap_2d->setColorScale(m_colorScale); // associate the color map with the color scale
    }

    else {
        m_cmap_2d->setGradient(*m_colorGradient);
    }
    replot();
}

// ----------------------------------------------------------------------------
/// save current state of the graphs
void ACViewer::saveGraphAsPdf(QCPGraph *g1, QCPGraph *g2, QString xlabel, QString ylabel, bool show_legend,
                              QString filename, int X, int Y) {

    QCustomPlot *customPlot = new QCustomPlot();

    customPlot->setLocale(QLocale(QLocale::English, QLocale::UnitedKingdom)); // period as decimal separator and comma as thousand separator
    customPlot->legend->setVisible(show_legend);

    // by default, the legend is in the inset layout of the main axis rect. So this is how we access it to change legend placement:
    customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignLeft);
    customPlot->xAxis->grid()->setPen(Qt::NoPen);
    customPlot->yAxis->grid()->setPen(Qt::NoPen);
    customPlot->xAxis->setLabel(xlabel);
    customPlot->yAxis->setLabel(ylabel);

    if(g1 != 0 && g2 != 0)
        customPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);

    if(g1 != 0){
        QCPGraph *graph = customPlot->addGraph(customPlot->xAxis, customPlot->yAxis);

        graph->setLineStyle(g1->lineStyle());
        graph->setPen(g1->pen());
        graph->setBrush(g1->brush());
        graph->setData(g1->data(), true);
        graph->setName(g1->name());
        graph->rescaleAxes();

        // use the same ticks
        if(!g1->keyAxis()->autoTicks()){
            customPlot->xAxis->setAutoTicks(false);
            customPlot->xAxis->setTickVector( g1->keyAxis()->tickVector() );
        }
    }
    if(g2 != 0){
        QCPGraph *graph = customPlot->addGraph(customPlot->xAxis, customPlot->yAxis);

        graph->setLineStyle(g2->lineStyle());
        graph->setPen(g2->pen());
        graph->setBrush(g2->brush());
        graph->setData(g2->data(), true);
        graph->setName(g2->name());
        graph->rescaleAxes(true);

        // use the same ticks
        if(!g2->keyAxis()->autoTicks()){
            customPlot->xAxis->setAutoTicks(false);
            customPlot->xAxis->setTickVector( g2->keyAxis()->tickVector() );
        }
    }



    customPlot->resize(X,Y);
    customPlot->show();
    customPlot->savePdf(filename);
    delete customPlot;
}

/// save current state of the graphs
void ACViewer::saveGraphs(){

    // file size
    int X = 500;
    int Y = 210;

    // save 1d plots
    saveGraphAsPdf(m_grph_tscSlice_angles, 0, "Time (ps)", "Angle (Degrees)", false, "graph_tsc_angles.pdf", X, Y);
    saveGraphAsPdf(m_grph_tscSlice_hist, 0, "Angle (Degrees)", "Frequency", false, "graph_tsc_hist.pdf", X, Y);
    saveGraphAsPdf(m_grph_thtSlice_hist, 0, "Time-scale (ps)", "Frequency", false, "graph_tht_hist.pdf", X, Y);
    saveGraphAsPdf(m_grph_tscError_vis, m_grph_tscError_data, "Time-scale (ps)", "Mean L2 Error", true, "graph_errors.pdf", X, Y);

    // now save the 2d plot
    QCustomPlot *customPlot = new QCustomPlot();

    customPlot->setLocale(QLocale(QLocale::English, QLocale::UnitedKingdom)); // period as decimal separator and comma as thousand separator
    //customPlot->legend->setVisible(true);

    // by default, the legend is in the inset layout of the main axis rect. So this is how we access it to change legend placement:
    customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignRight);
    customPlot->xAxis->grid()->setPen(Qt::NoPen);
    customPlot->yAxis->grid()->setPen(Qt::NoPen);

    QCPColorMap *colorMap = new QCPColorMap(customPlot->xAxis, customPlot->yAxis);
    customPlot->addPlottable(colorMap);

    colorMap->setData(m_cmap_2d->data(), true);
    colorMap->keyAxis()->setRange(m_cmap_2d->keyAxis()->range());
    colorMap->valueAxis()->setRange(m_cmap_2d->valueAxis()->range());
    colorMap->setDataRange(m_cmap_2d->dataRange());

    colorMap->setGradient(m_cmap_2d->gradient());
    colorMap->setInterpolate(m_cmap_2d->interpolate());
    colorMap->setTightBoundary(m_cmap_2d->tightBoundary());
    colorMap->setVisible(true);

    if(!m_grph_thtSlice_hist->keyAxis()->autoTicks()){
        colorMap->keyAxis()->setAutoTicks(false);
        colorMap->keyAxis()->setTickVector(m_grph_thtSlice_hist->keyAxis()->tickVector());
    }

    /*QCPColorScale *colorscale = new QCPColorScale(this);
    colorscale->setGradient(m_cmap_2d->gradient());
    colorscale->setType(QCPAxis::atRight); // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)
    colorMap->setColorScale(colorscale); // associate the color map with the color scale
    colorMap->colorScale()->setVisible(true);*/

    QFont f = customPlot->yAxis->labelFont();
    f.setPointSize(14);
    customPlot->xAxis->setLabelFont(f);
    customPlot->yAxis->setLabelFont(f);
    customPlot->yAxis->setLabel("Angle (Degrees)");
    customPlot->xAxis->setLabel("Time-Scale (ps)");

    //int w = m_axisRects[AX_CC]->width();
    //int h = 3*w/4;
    int w = 600;
    int h = 600;
    customPlot->resize(w,h);
    customPlot->show();
    customPlot->savePdf("colormap.pdf");
    customPlot->savePng("colormap.png");
    delete customPlot;

}

/// create graphs
void ACViewer::tsc_hist(float tscale, uint tscale_idx){

    /*if(parent->show_rawPlots()){
        return tsc_hist(tscale);
    }*/

    QString tag("Histogram for tscale = ");
    const QVector<double> &axisTheta = parent->get_data()->axis_theta();
    const histogram &histo = parent->get_data()->hist(tscale_idx);

    m_grph_tscSlice_hist->setData(axisTheta, histo.get_values());

    tag.append(QString::number(tscale));
    m_grph_tscSlice_hist->setName(tag);
    m_grph_tscSlice_hist->rescaleAxes();

    return;
    int num_ticks = 5;
    QVector<double> tv = parent->get_data()->get_tht_ticks(num_ticks);

    m_grph_tscSlice_hist->keyAxis()->setAutoTicks(false);
    m_grph_tscSlice_hist->keyAxis()->setTickVector(tv);
}

void ACViewer::tsc_hist(float tscale){

    const QVector<double> &axisTheta = parent->get_data()->axis_theta();
    QString tag("Histogram for tscale = ");


    const QVector<double> &rawtscales = parent->get_data()->axis_tscaleraw();
    size_t rawtscale_idx = Utils::get_closestIndex(rawtscales, (double)tscale);

    tscale = rawtscales[rawtscale_idx];

    const histogram &histo = parent->get_data()->hist_raw(tscale);
    m_grph_tscSlice_hist->setData(axisTheta, histo.get_values());

    tag.append(QString::number(tscale));
    m_grph_tscSlice_hist->setName(tag);
    m_grph_tscSlice_hist->rescaleAxes();


    return;
    int num_ticks = 5;
    QVector<double> tv = parent->get_data()->get_tht_ticks(num_ticks);

    m_grph_tscSlice_hist->keyAxis()->setAutoTicks(false);
    m_grph_tscSlice_hist->keyAxis()->setTickVector(tv);
}

void ACViewer::tht_hist(float theta, uint theta_idx){

    const QVector<double> &axisTScale = parent->get_data()->axis_tscale();
    QVector<double> histo = parent->get_data()->values_theta(theta_idx);

    QString tag("Histogram for angle = ");
    tag.append(QString::number(theta));

    m_grph_thtSlice_hist->setData(axisTScale, histo);
    m_grph_thtSlice_hist->setName(tag);
    m_grph_thtSlice_hist->rescaleAxes();



    int num_ticks = 5;
    QVector<double> tv = parent->get_data()->get_tsc_ticks(num_ticks);

    m_grph_thtSlice_hist->keyAxis()->setAutoTicks(false);
    m_grph_thtSlice_hist->keyAxis()->setTickVector(tv);
}

void ACViewer::tsc_angles(float tscale){

    // angle vs time
        // angles are available wrt raw time
        // but tscale_idx comes from reduce time

    size_t idx = Utils::get_closestIndex(this->parent->get_data()->axis_tscaleraw(), (double)tscale);
    //int idx = Utils::snap_to_uniformAxis((double)tscale, this->parent->get_data()->axis_tscaleraw());

    QVector<double> angles = parent->get_data()->angles_tscale(idx);
    QVector<double> axisTime = parent->get_data()->times_tscale(idx);

    QString tag("Angles for tscale = ");
    tag.append(QString::number(parent->get_data()->axis_tscaleraw()[idx]));
    m_grph_tscSlice_angles->setData(axisTime, angles);
    m_grph_tscSlice_angles->setName(tag);

    m_grph_tscSlice_angles->rescaleAxes();
}

void ACViewer::errors(){

    // vis error
    const QVector<double> &viserror_axis = parent->get_data()->axis_tscale();
    size_t sz_tscale = this->parent->get_data()->num_tscales();
    QVector<double> viserror(sz_tscale);
    for(uint tscale = 0; tscale < sz_tscale; ++tscale){

        viserror[tscale] = this->parent->get_data()->hist(tscale).L2avg();
        viserror[tscale] *= 100.0 /(float)this->parent->get_data()->hist(tscale).num_samples();
    }

    m_grph_tscError_vis->setData(viserror_axis, viserror);
    m_grph_tscError_vis->valueAxis()->setScaleType(QCPAxis::stLogarithmic);
    m_grph_tscError_vis->rescaleAxes();

    // data error
    if(parent->show_datauq() && !parent->get_datadir().isEmpty()){

        QString basefilename = parent->get_datadir();
        if(!basefilename.endsWith('/'))
            basefilename.append("/");
        basefilename.append(QString::number(this->parent->get_atomId())).append("_bs_");

        //printf(" basefilename = %s\n", basefilename.toLatin1().data());

        const QVector<double> &dataerror_axis = parent->get_data()->axis_tscaleraw();
        static size_t sz_rawtscale = this->parent->get_data()->num_maxTscales();
        static QVector<double> dataerror(sz_rawtscale);

        static uint dataerror_bins = 0;

        // do not recompute this if i have not changed the number of bins
        if(dataerror_bins != this->parent->get_data()->num_bins()){

            dataerror_bins = this->parent->get_data()->num_bins();
            printf(" Reading data errors for %d bins...", dataerror_bins);
            fflush(stdout);

            for(uint tscale = 0; tscale < sz_rawtscale; ++tscale){

                uint delta = 10 + tscale*10;    // starting tidx and jump

                QString filename = basefilename;
                filename.append(QString::number(delta)).append(".raw");

                float l2max = 0, l2avg = 0, l0max = 0, l0avg = 0;
                bool found = BootStrap::read(filename, dataerror_bins, l2max, l2avg, l0max, l0avg);
                dataerror[tscale] = l2avg;
                dataerror[tscale] *= 100.0/(float)parent->get_data()->angles_tscale(tscale).size();
            }
            printf(" Done!\n");
        }

        m_grph_tscError_data->setData(dataerror_axis, dataerror);
        //m_grph_tscError_data->valueAxis()->rescale();// >rescaleAxes(true);
        m_grph_tscError_data->rescaleValueAxis(true);
        //m_grph_tscError_vis->rescaleValueAxis()
    }
    else {
        m_grph_tscError_data->clearData();
    }

    int num_ticks = 5;
    QVector<double> tv = parent->get_data()->get_tsc_ticks(num_ticks);

    m_grph_tscError_vis->keyAxis()->setAutoTicks(false);
    m_grph_tscError_vis->keyAxis()->setTickVector(tv);
}

void ACViewer::createGraphs(bool verbose){

    verbose = false;
    QPointF sPoint = this->m_pSelector->getCoords();
    float tscale = sPoint.x();
    float theta = sPoint.y();

    if(verbose){
        printf("\n\n -- createGraphs for tscale %f, theta %f\n", tscale, theta);
    }

    // use the colormap to convert coords to pixels
    int tscale_idx, theta_idx;
    m_cmap_2d->data()->coordToCell(tscale, theta, &tscale_idx, &theta_idx);

    if(verbose){
        printf("\n\n -- createGraphs for tscale %f (%d), theta %f (%d)\n", tscale, tscale_idx, theta, theta_idx);
    }


    size_t nTscale = parent->get_data()->num_tscales();

    if(tscale_idx >= nTscale){

        printf(" ACViewer::set_graphs_theta -- Error. requested tscale %d, num num_tscale = %d\n", tscale_idx, nTscale);

        m_grph_tscSlice_hist->clearData();
        m_grph_tscSlice_angles->clearData();
    }
    else {
        tsc_hist(tscale, tscale_idx);
        tsc_angles(tscale);
    }

    size_t nTheta = parent->get_data()->num_bins();
    if(theta_idx >= nTheta){
        printf(" ACViewer::create_graphs_tscale -- Error. requested theta = %d, but num_theta = %d\n", theta_idx, nTheta);
        m_grph_thtSlice_hist->clearData();
    }
    else {

        tht_hist(theta, theta_idx);
    }
    errors();
    //need_1dplots = false;
}

void ACViewer::create2DPlot(bool verbose){

    //verbose = true;
    if(verbose){
        printf("\n -- create2Dplot\n");
    }

    if(parent->get_data() == 0){
        printf(" cannot create texture. currdata is null!\n");
        return;
    }
    m_cmap_2d->clearData();

    m_cmap_2d->setKeyAxis(get_keyAxis(AX_CC));
    m_cmap_2d->setValueAxis(get_valueAxis(AX_CC));
    m_pSelector->setAxes(get_keyAxis(AX_CC), get_valueAxis(AX_CC));

    size_t nTheta = parent->get_data()->num_bins();
    size_t nTscale = parent->get_data()->num_tscales();

    if(verbose){
        printf(" ntscale = %d, keyaxisrange = (%f %f)\n", nTscale, get_keyAxis(AX_CC)->range().lower, get_keyAxis(AX_CC)->range().upper);
        printf(" ntheta = %d, valaxisrange = (%f %f)\n", nTheta, get_valueAxis(AX_CC)->range().lower, get_valueAxis(AX_CC)->range().upper);
        printf(" data range = (%f %f)\n", parent->get_data()->range_tscale().lower, parent->get_data()->range_tscale().upper);
    }

    QCPColorMapData *cdata = new QCPColorMapData( nTscale, nTheta, parent->get_data()->range_tscale(), parent->get_data()->range_theta() );
    for(size_t theta = 0; theta < nTheta; ++theta){
    for(size_t tscale = 0; tscale  < nTscale; ++tscale ){

        if(live_show_unnormalized)      cdata->setCell(tscale, theta, parent->get_data()->count(theta, tscale));
        else                            cdata->setCell(tscale, theta, parent->get_data()->value(theta, tscale));
    }
    }

    m_cmap_2d->setData(cdata);


    // recompute the color range
    QCPRange vrange = cdata->dataBounds();
    parent->reconcile_valueRange(vrange);
    //m_cmap_2d->rescaleDataRange();
    m_cmap_2d->setDataRange(vrange);

    m_cmap_2d->rescaleAxes();

    //need_2dplot = false;

    if(verbose){
        printf("\n ntscale = %d, keyaxisrange = (%f %f)\n", nTscale, get_keyAxis(AX_CC)->range().lower, get_keyAxis(AX_CC)->range().upper);
        printf(" ntheta = %d, valaxisrange = (%f %f)\n", nTheta, get_valueAxis(AX_CC)->range().lower, get_valueAxis(AX_CC)->range().upper);
        printf(" data range = (%f %f)\n", parent->get_data()->range_tscale().lower, parent->get_data()->range_tscale().upper);
    }
}


/// ============================================================================

void ACViewer::convergence(){

    QCustomPlot *customPlot = new QCustomPlot();

    customPlot->setLocale(QLocale(QLocale::English, QLocale::UnitedKingdom)); // period as decimal separator and comma as thousand separator
    customPlot->legend->setVisible(true);

    //QFont legendFont = font();  // start out with MainWindow's font..
    //legendFont.setPointSize(12); // and make a bit smaller for legend
    //customPlot->legend->setFont(legendFont);
    //customPlot->legend->setBrush(QBrush(QColor(255,255,255,230)));

    // by default, the legend is in the inset layout of the main axis rect. So this is how we access it to change legend placement:
    customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignHCenter|Qt::AlignRight);

    customPlot->xAxis->grid()->setPen(Qt::NoPen);
    customPlot->yAxis->grid()->setPen(Qt::NoPen);
    customPlot->xAxis->setLabel("Number of Resampling Steps");
    customPlot->yAxis->setLabel("Error Metric");

    customPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    customPlot->yAxis->setScaleLogBase(10);
    customPlot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    customPlot->xAxis->setScaleLogBase(10);

    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);


    // ------------
    for(uint ssp = 40; ssp < 100; ssp+=10){

        QVector<double> l2, l0;
        BootStrap::read_forconvergence(30, 180, ssp, l2, l0);

        uint sz = l2.size();

        QVector<double> l2max(sz), l2avg(sz);
        QVector<double> l0max(sz), l0avg(sz);

        QVector<double> ax(sz);

        ax[0] = 0;
        l2max[0] = l2[0];       l2avg[0] = l2[0];
        l0max[0] = l0[0];       l0avg[0] = l0[0];

        long double l2sum = l2[0];
        long double l0sum = l0[0];

        for(uint i = 1; i < sz; i++){

            ax[i] = i;

            l2max[i] = (l2max[i-1] > l2[i]) ? l2max[i-1] : l2[i];
            l0max[i] = (l0max[i-1] > l0[i]) ? l0max[i-1] : l0[i];

            l2sum += l2[i];
            l2avg[i] = (l2sum / (i+1));

            l0sum += l0[i];
            l0avg[i] = (l0sum / (i+1));
        }

        QCPGraph *gl2avg = customPlot->addGraph(customPlot->xAxis, customPlot->yAxis);
        QCPGraph *gl0avg = customPlot->addGraph(customPlot->xAxis, customPlot->yAxis);

        QColor l2color(rand()%256, rand()%256, rand()%256, 255);

        //QString tl2("50 %");
        //QString tl0("L0");
        QString t = QString::number(ssp).append("% resampling");

        QPen l2pen(l2color);
        l2pen.setWidth(2);

        gl2avg->setLineStyle(QCPGraph::lsLine);
        gl2avg->setPen(l2pen);
        //QString t = tl2;
        //gl2avg->setName(t.append(" avg ").append(QString::number(ssp)));
        gl2avg->setName(t);

        l2pen.setStyle(Qt::DashLine);

        gl0avg->setLineStyle(QCPGraph::lsLine);
        gl0avg->setPen(l2pen);
        //t = tl0;
        //gl0avg->setName(t.append(" avg ").append(QString::number(ssp)));

        gl2avg->setData(ax, l2avg);
        gl0avg->setData(ax, l0avg);

        customPlot->legend->removeItem(customPlot->legend->itemCount()-1); // don't show two confidence band graphs in legend
    }


    customPlot->xAxis->rescale();
    customPlot->yAxis->rescale();
    customPlot->resize(600, 350);
    customPlot->savePdf("convergence.pdf");
    delete customPlot;
}

#include <iostream>
// compute the axis parameters for display
bool ACViewer::compute_axesParams(size_t &nTheta, size_t &nTscale, QCPRange &rngTscale){

    // findout the max size of the data for this texture
    QSize plotsize = m_axisRects[AX_CC]->size();

    if(plotsize.width() == 0 || plotsize.width() == 0){
        return false;
    }

    if(plotsize.width() < 0 || plotsize.width() < 0){
        plotsize = QSize(400,400);
    }

    if(m_orientation_theta == Qt::Vertical){
        nTscale = plotsize.width() / min_pixel_per_value;
        nTheta = plotsize.height() / min_pixel_per_value;
        rngTscale = m_axisRects[AX_CC]->axis(QCPAxis::atBottom)->range();
    }
    else {
        nTscale = plotsize.height() / min_pixel_per_value;
        nTheta = plotsize.width() / min_pixel_per_value;
        rngTscale = m_axisRects[AX_CC]->axis(QCPAxis::atLeft)->range();
    }
    return true;
}

/// ============================================================================
/// event handlers

/// before every replot, check if anything needs to be updated
void ACViewer::replotEvent() {

    static int count = 0;

    // find axes range and number of samples that I can display
    size_t nTheta, nTscale;
    QCPRange qRngTscale;

    bool is_valid = compute_axesParams(nTheta, nTscale, qRngTscale);
    if(!is_valid){
        //printf(" invalid!\n");
        return;
    }

    // request the parent ACWindow to compute the data
    bool datachanged = parent->recompute_data(nTheta, nTscale, qRngTscale);

    if(datachanged){
        create2DPlot(false);
        parent->updateTitle();
    }
    else {
        QCPRange vrange = m_cmap_2d->dataRange();
        parent->reconcile_valueRange(vrange);
        m_cmap_2d->setDataRange(vrange);
        m_cmap_2d->rescaleAxes();
    }

    //printf(" before snapping to axis");
    //this->m_pSelector->print();

    // due to zoom and drag, the selector may have moved out of range
    // i need to fix that
    this->m_pSelector->snap_to_axes(parent->get_data());

    //printf(" after snapping to axis");
    //this->m_pSelector->print();

    createGraphs(false);


    // set axes!
    int num_ticks = 5;
    QVector<double> tv = parent->get_data()->get_tsc_ticks(num_ticks);

    this->m_grph_thtSlice_hist->keyAxis()->setAutoTicks(false);
    this->m_grph_thtSlice_hist->keyAxis()->setTickVector(tv);


    count++;
}

#if 0
// before addressing Mitchell's request (02.16.16)
void ACViewer::replotEvent() {

    static int count = 0;

    // findout the max size of the data for this texture
    QSize plotsize = this->m_axisRects[AX_CC]->size();
    int X = plotsize.width() / min_pixel_per_value;
    int Y = plotsize.height() / min_pixel_per_value;

    if(X < 1 || Y < 1)
        return;

    //m_axisRects[AX_CR]->setMinimumSize(m_axisRects[AX_CR]->size());


    // find the range of theta and tscale values
    // and the number of samples i can display
    QCPRange qRngTheta, qRngTscale;
    size_t nTheta, nTscale;

    if(m_orientation_theta == Qt::Vertical){
        nTheta = Y;
        nTscale = X;

        qRngTheta = m_axisRects[AX_CC]->axis(QCPAxis::atLeft)->range();
        qRngTscale = m_axisRects[AX_CC]->axis(QCPAxis::atBottom)->range();
    }
    else {
        nTheta = X;
        nTscale = Y;

        qRngTscale = m_axisRects[AX_CC]->axis(QCPAxis::atLeft)->range();
        qRngTheta = m_axisRects[AX_CC]->axis(QCPAxis::atBottom)->range();
    }

    // snap range to data
    //printf(" replot() %d - %p %p %p\n", count, parent->get_data(), m_axisRects[AX_CC]->axis(QCPAxis::atLeft), m_axisRects[AX_CC]->axis(QCPAxis::atBottom));

    const QVector<double> &tsc = parent->get_data()->axis_tscaleraw();
    int idx_tscale0 = Utils::snap_to_uniformAxis(qRngTscale.lower, tsc);
    int idx_tscale1 = Utils::snap_to_uniformAxis(qRngTscale.upper, tsc)+1;

    if(idx_tscale0 == idx_tscale1)      idx_tscale1 = idx_tscale0+1;
    if(idx_tscale1 > tsc.size()-1)      idx_tscale1 = tsc.size()-1;

    qRngTscale.lower = tsc[idx_tscale0];
    qRngTscale.upper = tsc[idx_tscale1];

    if(qRngTheta.lower < parent->get_data()->range_rawtheta().lower){
        qRngTheta.lower = parent->get_data()->range_rawtheta().lower;
    }
    if(qRngTheta.upper > parent->get_data()->range_rawtheta().upper){
        qRngTheta.upper = parent->get_data()->range_rawtheta().upper;
    }

    //printf("\n\n -- before replot %d! texture = (%d %d). nTscale = %d (%f %f), nTheta = %d (%f %f)\n", count,
   //     plotsize.width(), plotsize.height(), nTscale, qRngTscale.lower, qRngTscale.upper, nTheta, qRngTheta.lower, qRngTheta.upper);

    bool datachanged = parent->get_data()->compute(nTscale, idx_tscale0, idx_tscale1, qRngTheta, nTheta, false);

    if(datachanged){

        create2DPlot(false);
        parent->updateTitle();
    }
    else {
        m_cmap_2d->rescaleAxes();
    }

    //printf(" before snapping to axis");
    //this->m_pSelector->print();


    // due to zoom and drag, the selector may have moved out of range
    // i need to fix that
    this->m_pSelector->snap_to_axes(parent->get_data());

    //printf(" after snapping to axis");
    //this->m_pSelector->print();

    createGraphs(false);

    count++;
}
#endif
/// start selecting if SHIFT has been pressed
void ACViewer::mousePressEvent(QMouseEvent *event) {

    if (event->button() == Qt::LeftButton && event->modifiers() == Qt::ShiftModifier) {

        if( m_pSelector->testPointInside(event->localPos()) && m_pSelector->visible() ){
            m_pSelector->setSelected(true);
            setCursor(Qt::OpenHandCursor);
            this->update_interactionMode(MODE_SELECT);
        }
    }
    QCustomPlot::mousePressEvent(event);
}


/// move the selector if it is in select mode
void ACViewer::mouseMoveEvent(QMouseEvent *event){

    if (m_pSelector->selected()) {

        QPointF pix_widgt = event->localPos();

        float pix_theta = 0, pix_tscale = 0;
        if(m_orientation_theta == Qt::Vertical) {
            pix_theta = pix_widgt.y();      pix_tscale = pix_widgt.x();
        }
        else {
            pix_theta = pix_widgt.x();      pix_tscale = pix_widgt.y();
        }

        // convert the pixel to value (on continuous axis)
        double tscale_val = m_cmap_2d->keyAxis()->pixelToCoord(pix_tscale);
        double theta_val = m_cmap_2d->valueAxis()->pixelToCoord(pix_theta);

        m_pSelector->snap_to_axes(parent->get_data(), tscale_val, theta_val);
        replot();       // needed to show the selector moving alongwith the cursor
    }
    QCustomPlot::mouseMoveEvent(event);
}


/// unselect the selector and recompute the 1D plots
void ACViewer::mouseReleaseEvent(QMouseEvent *event){

    if( m_pSelector->selected() ){

        m_pSelector->setSelected(false);
        unsetCursor();
        //need_1dplots = true;

        QPointF sPoint = this->m_pSelector->getCoords();
        signal_selectorMoved(this->parent->get_atomId(), sPoint.x(), sPoint.y());
        replot();
    }
    QCustomPlot::mouseReleaseEvent(event);
}

/// change orientation, switch interpolation, and reset view
void ACViewer::keyPressEvent(QKeyEvent *event){

    if(event->modifiers() == Qt::ShiftModifier){
        update_interactionMode(MODE_SELECT);
        replot();
        return;
    }

    switch(event->key()){

    /*case 'n':
    case 'N':
            live_show_unnormalized = !live_show_unnormalized;
            break;

    case 'h':
    case 'H':
                set_interactionMode(MODE_ZOOM_H);       break;

    case 'v':
    case 'V':   set_interactionMode(MODE_ZOOM_V);       break;

    case 'b':
    case 'B':
                set_interactionMode(MODE_ZOOM_HV);       break;*/

   /* case 'c':
    case 'C':
            convergence();
            break;
            */
    /*case 'b':
      case 'B':   live_bootstrap = !live_bootstrap;
               printf(" swithcing live_bootstrap %d\n", live_bootstrap);
                    break;*/

    /*case 't':
    case 'T':
            //m_angleTracer->setVisible(!m_angleTracer->visible());
            break;*/
    /*case 'm':
    case 'M':
    if(m_pSelector->visible()){
    m_pSelector->setVisible(false);
    m_pSelector->linesVisible(false);
    }
    else {
    m_pSelector->setVisible(true);
    m_pSelector->linesVisible(true);
    }

    break;
    */

#if BLOCKED
    case 'o':
    case 'O':
                m_orientation_theta = (m_orientation_theta == Qt::Vertical) ? Qt::Horizontal : Qt::Vertical;
                update_orientation();
                break;
#endif
        case 'i':
        case 'I':
                m_cmap_2d->setInterpolate(!m_cmap_2d->interpolate());
                break;

        case 'r':
        case 'R':
                {
                QCPRange x = (m_orientation_theta == Qt::Vertical) ? parent->get_data()->range_rawtscale() : parent->get_data()->range_rawtheta();
                QCPRange y = (m_orientation_theta == Qt::Horizontal) ? parent->get_data()->range_rawtscale() : parent->get_data()->range_rawtheta();
                m_axisRects[AX_CC]->axis(QCPAxis::atBottom)->setRange(x);
                m_axisRects[AX_CC]->axis(QCPAxis::atLeft)->setRange(y);
                }
                break;


        case 'l':
        case 'L':
                m_lgnd_tscError->setVisible( !m_lgnd_tscError->visible() );
                m_lgnd_tscSlice_angles->setVisible( !m_lgnd_tscSlice_angles->visible() );
                m_lgnd_tscSlice_hist->setVisible( !m_lgnd_tscSlice_hist->visible() );
                m_lgnd_thtSlice_hist->setVisible( !m_lgnd_thtSlice_hist->visible() );
                if(m_colorScale)
                    m_colorScale->setVisible( !m_colorScale->visible());
                break;

    }

    QCustomPlot::keyPressEvent(event);
    replot();
}
/// unselect the selector (SHIFT has ben released)
void ACViewer::keyReleaseEvent(QKeyEvent *event){

    if(m_pSelector->visible()){
        update_interactionMode(MODE_ZOOM_H);
    }

    QCustomPlot::keyReleaseEvent(event);
    replot();
}

// ============================================================================
