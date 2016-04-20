#ifndef _HISTOGRAM_VIEWER_H
#define _HISTOGRAM_VIEWER_H

#include <QGLViewer/qglviewer.h>

// forward declaration!
class Atom;
class TransferFunction;
//class Histogram2D;

class QPaintEvent;
class QPainter;

#if 0
/// ===========================================================
/// new using QPainter (2D)
class QPaintEvent;
class QPainter;
class HViewer : public QGLViewer {

    uint atom_id;
    GLuint text_id;

    uint min_del, max_del;
    float max_val;

    const TransferFunction *tfunc;

public :

    HViewer(uint id, QWidget* parent = 0)
            : QGLViewer(QGLFormat(QGL::SampleBuffers), parent) {

        //printf("hv::hv()\n");
        setAttribute(Qt::WA_NoSystemBackground);

        atom_id = id;
        text_id = 0;
    }

    void create_texture(const Atom &atom, const TransferFunction *tf);
    //void draw_scale(const TransferFunction *tf, float min_v, float max_v);

protected :
    virtual void init();
    virtual void draw();

    virtual void paintGL() { update(); }
    virtual void paintEvent(QPaintEvent *event);

    //virtual QString helpString() const;
};
#endif
/// ===========================================================
/// HistogramViewer
/// ===========================================================

class HistogramViewer : public QGLViewer{

    uint atom_id;
    GLuint text_id;

    uint min_del, max_del;
    float max_val;

    const TransferFunction *tfunc;

public:

    HistogramViewer(uint id, QWidget* parent = 0)
            : QGLViewer(QGLFormat(QGL::SampleBuffers), parent) {

        //printf("hv::hv()\n");
        setAttribute(Qt::WA_NoSystemBackground);

        atom_id = id;
        text_id = 0;
    }

    void create_texture(const char *filename);

    void create_texture(const Atom &atom, const TransferFunction *tf);
    //void create_texture(const Histogram2D &hist2d, const TransferFunction *tf);
    void draw_scale(const TransferFunction *tf, float min_v, float max_v);

    //void drawOverpaint(QPainter *painter);
protected :

    virtual void draw();
    virtual void init();

    //virtual void paintGL() { update(); }
    //virtual void paintEvent(QPaintEvent *event);

    //virtual QString helpString() const;
};

#endif
