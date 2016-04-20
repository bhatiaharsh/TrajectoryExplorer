/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#ifndef _ACPLOTS_H_
#define _ACPLOTS_H_

#include <qcustomplot.h>
#include <ACData.h>


class PlotSelector : public QCPItemEllipse {

    Q_OBJECT

private:
    float halfSize;

    QCPItemTracer *mCenterTracer;
    QCPItemStraightLine *mHorizontal, *mVertical;
    QCPItemText* mText;

    void anchorToTracer(QCPItemPosition *pos, const QPointF &coords){

        // the coords are defined with respect to the anchor
        pos->setParentAnchor(mCenterTracer->position);
        pos->setType(QCPItemPosition::ptAbsolute);
        pos->setCoords(coords);
    }

public:

    explicit PlotSelector(QCustomPlot *parentPlot, float halfSize_ = 5)
                : QCPItemEllipse(parentPlot), halfSize(halfSize_){

        // do not show the tracer (it will be at the center of the ellipse)
        mCenterTracer = new QCPItemTracer(parentPlot);
        mCenterTracer->setStyle(QCPItemTracer::tsNone);

        // ellipse is drawn using two points (top-left and bot-right)
        anchorToTracer(topLeft, QPointF(-halfSize_, -halfSize_));
        anchorToTracer(bottomRight, QPointF(halfSize_, halfSize_));

        // line is drawn using two points (point1 and point2)
        mHorizontal = new QCPItemStraightLine(parentPlot);
        mVertical = new QCPItemStraightLine(parentPlot);

        anchorToTracer(mHorizontal->point1, QPointF(0,0));
        anchorToTracer(mHorizontal->point2, QPointF(1,0));
        anchorToTracer(mVertical->point1, QPointF(0,0));
        anchorToTracer(mVertical->point2, QPointF(0,1));

        mText = new QCPItemText(parentPlot);
        anchorToTracer(mText->position, QPointF(halfSize_, -1.5*halfSize_));
        mText->setClipToAxisRect(0);
        mText->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
        //mText->position->setType(QCPItemPosition::ptAxisRectRatio);
        //mText->setFont(QFont(font().family(), 12));

        // styline and coloring
        setPen(QPen(Qt::black));
        setSelectedPen(QPen(Qt::black, 2));
        setBrush(QBrush(QColor(qrand()%256, qrand()%256, qrand()%256, 100)));

        QPen linepen(QBrush(Qt::black), 1.75, Qt::DashDotDotLine);

        mHorizontal->setPen(linepen);
        mHorizontal->setAntialiased(false);
        mVertical->setPen(linepen);
        mVertical->setAntialiased(false);

        this->setLinesVisible(true);

        mText->setText("");
        mText->setVisible(true);

        // enable drawing and selecting
        parentPlot->addItem(this);
        parentPlot->addItem(mText);
        parentPlot->addItem(mHorizontal);
        parentPlot->addItem(mVertical);

        //setVisible(true);
        setSelectable(true);
    }

    void setVisible(bool visible){
        mText->setVisible(visible);
        QCPItemEllipse::setVisible(visible);
    }

    void setValueRange(QCPRange axisRng){

        QPointF p = mCenterTracer->position->coords();
        float py = p.y();
        if(py < axisRng.lower)           py = axisRng.lower;
        else if(py > axisRng.upper)      py = axisRng.upper;

        mCenterTracer->position->setCoords(p.x(), py);
    }
    void setKeyRange(QCPRange axisRng){

        QPointF p = mCenterTracer->position->coords();
        float px = p.x();
        if(px < axisRng.lower)           px = axisRng.lower;
        else if(px > axisRng.upper)      px = axisRng.upper;

        mCenterTracer->position->setCoords(px, p.y());
    }

    void setAxes(QCPAxis *keyAxis, QCPAxis *valueAxis){

        mCenterTracer->position->setAxes(keyAxis, valueAxis);
    }

    void setAxisRect(QCPAxisRect *rect){

        this->setClipAxisRect(rect);
        mHorizontal->setClipAxisRect(rect);
        mVertical->setClipAxisRect(rect);
    }

    void setLinesVisible(bool vis){
        this->mHorizontal->setVisible(vis);
        this->mVertical->setVisible(vis);
    }

    void setTscale(float tscale){
        mCenterTracer->position->setCoords(tscale, mCenterTracer->position->value());
        update_text();
    }
    void setTheta(float theta){
        mCenterTracer->position->setCoords(mCenterTracer->position->value(), theta);
        update_text();
    }
    void setCoords(float tscale, float theta){
        mCenterTracer->position->setCoords(tscale, theta);
        update_text();
    }

    void update_text(){
        QPointF p = mCenterTracer->position->coords();
        QString s("T = ");
        s.append(QString::number(p.x())).append(", A = ").append(QString::number(p.y()));
        mText->setText(s);
    }

    void snap_to_axes(const ACData *data, double px, double py){

        // coords can only take one of the values in the axis
        const QVector<double> &tsc = data->axis_tscale();
        const QVector<double> &tht = data->axis_theta();

        int x = Utils::get_closestIndex(tsc, px);
        int y = Utils::get_closestIndex(tht, py);

        //int x = Utils::snap_to_axis(px, data->axis_tscale());
        //int y = Utils::snap_to_axis(py, data->axis_theta());

        setCoords(tsc[x], tht[y]);
    }

    void snap_to_axes(const ACData *data){

        QPointF sPoint = this->getCoords();
        snap_to_axes(data, sPoint.x(), sPoint.y());
    }

    bool testPointInside(const QPointF &p) const{
        const QPointF &c = mCenterTracer->position->pixelPoint();
        return (sqrt( (p.x()-c.x())*(p.x()-c.x()) + (p.y()-c.y())*(p.y()-c.y()) ) < halfSize);
    }

    QPointF getCoords() const {    return mCenterTracer->position->coords();   }

    void print() const {

        const QPointF &p = mCenterTracer->position->pixelPoint();
        const QPointF &c = mCenterTracer->position->coords();
        //printf(" p = (%.1f, %.1f), c = (%.1f, %.1f)\n", p.x(), p.y(), c.x(), c.y());
        printf(" p = (%f, %f), c = (%f, %f)\n", p.x(), p.y(), c.x(), c.y());
    }
};


class StatisticalPlot  {

    QCustomPlot *m_parentPlot;
    QCPAxis *m_keyAxis, *m_valueAxis;
    QBrush *m_brush;

    std::vector<QCPStatisticalBox *> m_data;

public:
    StatisticalPlot(QCustomPlot *parentPlot, QCPAxis *keyAxis, QCPAxis *valueAxis, QBrush *brush = 0)
                    : m_parentPlot (parentPlot), m_brush(brush),
                      m_keyAxis(keyAxis), m_valueAxis(valueAxis)
    {}

    void setBrush(QBrush *brush){   m_brush = brush;    }

    void setAxes(QCPAxis *keyAxis, QCPAxis *valueAxis){

        m_keyAxis = keyAxis;
        m_valueAxis = valueAxis;

        for(uint i = 0; i < m_data.size(); i++){
            QCPStatisticalBox *sBox = m_data[i];
            sBox->setKeyAxis(m_keyAxis);
            sBox->setValueAxis(m_valueAxis);
        }
    }

    void rescaleValueAxis() {
        m_valueAxis->rescale();
    }
    void rescaleAxes() {
        m_keyAxis->rescale();
        m_valueAxis->rescale();
    }

    void clearData(){

        for(uint i = 0; i < m_data.size(); i++){
            m_parentPlot->removePlottable(m_data[i]);
        }
        m_data.clear();
    }

    void setData(const QVector<double> x, const QVector<double> y, const QVector<double> err){

        clearData();
        size_t sz = std::min( std::min(x.size(), y.size()), err.size() );

        double width = 1.25 * (x.back()-x.front()) / (double)(x.size());

        for(uint i = 0; i < sz; i++){

            QCPStatisticalBox *sBox = new QCPStatisticalBox(m_keyAxis, m_valueAxis);

            sBox->setBrush(*m_brush);
            sBox->setKey(x[i]);
            sBox->setMedian(y[i]);
            sBox->setUpperQuartile(y[i] + err[i]);
            sBox->setLowerQuartile(y[i] - err[i]);
            sBox->setWidth(width);

            m_parentPlot->addPlottable(sBox);
            m_data.push_back(sBox);
        }
    }
};

#endif
