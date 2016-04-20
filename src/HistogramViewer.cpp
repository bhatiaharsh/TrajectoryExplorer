#include "Atom.h"
#include "HistogramViewer.h"
#include "transferFunction.h"
#include "Utils.h"

/// -----------------------------------------------------------------
using namespace std;

#include <QPainter>


#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE  0x809D
#endif

#if 0
void HViewer::paintEvent(QPaintEvent *event) {

    Q_UNUSED(event)

    QPainter painter;
    painter.begin(this);
    painter.setRenderHint(QPainter::Antialiasing);


    // Save current OpenGL state
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    // Reset OpenGL parameters
    glShadeModel(GL_SMOOTH);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_MULTISAMPLE);
    static GLfloat lightPosition[4] = { 1.0, 5.0, 5.0, 1.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    qglClearColor(backgroundColor());

    // Classical 3D drawing, usually performed by paintGL().
    preDraw();
    draw();
    postDraw();


    // Restore OpenGL state
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glPopAttrib();


    // --------------------------
    painter.save();

    //printf(" size = (%d, %d)\n", width(), height());
    painter.translate(width()/2, height()/2);

    QRadialGradient radialGrad(QPointF(-40, -40), 100);
    radialGrad.setColorAt(0, QColor(255, 255, 255, 100));
    radialGrad.setColorAt(1, QColor(200, 200, 0, 100));

    painter.setBrush(QBrush(radialGrad));
    painter.drawRoundRect(-100, -100, 200, 200);


    painter.restore();
    painter.end();
}
#endif
/// -----------------------------------------------------------------

void HistogramViewer::init(){


#if 01
    //printf("hv::init()\n");
    this->setWindowTitle("Angle Correlation Histogram for Atom "+QString::number(atom_id));

    this->camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);

    this->setSceneBoundingBox(qglviewer::Vec(0,0,0), qglviewer::Vec(1,1,1));
    this->setSceneCenter(qglviewer::Vec(0.5,0.5,0.5));
    this->camera()->centerScene();
    this->camera()->showEntireScene();

    this->setFPSIsDisplayed(true);
    this->setTextIsEnabled(true);

    this->setBackgroundColor(Qt::gray);
    this->setForegroundColor(Qt::black);

    this->setMouseBinding(Qt::NoModifier, Qt::LeftButton, NO_CLICK_ACTION);
    //this->setMouseBinding(Qt::NoModifier, Qt::RightButton, NO_CLICK_ACTION);

    // init the viewer
    glEnable(GL_NORMALIZE);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    glDisable(GL_CULL_FACE);

    glEnable(GL_DEPTH_TEST);
    glDepthMask(true);
    glDepthFunc(GL_LESS);

    glDisable(GL_LIGHTING);

    glEnable(GL_TEXTURE_2D);
    glDisable(GL_COLOR_MATERIAL);
#endif
}

#define GLVERTEX2POLAR(r,t)     glVertex2f(r*cos(t), r*sin(t));

void draw_horseshoe(float rad_in, float rad_out, float ang_min, float ang_max, uint ang_div){

    // for perspective correction!
    const float fact_out = sqrt(rad_out*rad_out*2);
    const float fact_in = sqrt(rad_in*rad_in*2);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //glOrtho(-1.1*rad_out,1.1*rad_out,-0.1*rad_out,1.1*rad_out,-1,1);
    glOrtho(-1.1*rad_out,1.1*rad_out,-1.1*rad_out,1.1*rad_out,-1,1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glBegin(GL_QUAD_STRIP);
    glColor3f(0,0,0);

    // first two points (angle = 0)
    glTexCoord4f(0, 0, 0, fact_in);
    GLVERTEX2POLAR(rad_in, ang_min);
    glTexCoord4f(fact_out, 0, 0, fact_out);
    GLVERTEX2POLAR(rad_out, ang_min);

    const float ang_diff = (ang_max-ang_min);
    for(uint i = 1; i < ang_div; i++){

        float tx_ang = (float)i/(float)(ang_div-1);
        float vt_ang = ang_min + tx_ang*ang_diff;

        glTexCoord4f(0, tx_ang*fact_in, 0, fact_in);
        GLVERTEX2POLAR(rad_in, vt_ang);
        glTexCoord4f(fact_out, tx_ang*fact_out, 0, fact_out);
        GLVERTEX2POLAR(rad_out, vt_ang);
    }
    glEnd();
}

void draw_horseshoe_old(float rad_in, float rad_out, float ang_min, float ang_max, uint ang_div){

    glBegin(GL_QUAD_STRIP);
    glColor3f(0,0,0);

    // first two points (angle = 0)
    glTexCoord2f(0, 0);
    GLVERTEX2POLAR(rad_in, ang_min);
    glTexCoord2f(1, 0);
    GLVERTEX2POLAR(rad_out, ang_min);

    const float ang_diff = (ang_max-ang_min);
    for(uint i = 1; i < ang_div; i++){

        float tx_ang = (float)i/(float)(ang_div-1);
        float vt_ang = ang_min + tx_ang*ang_diff;

        glTexCoord2f(0, tx_ang);
        GLVERTEX2POLAR(rad_in, vt_ang);
        glTexCoord2f(1, tx_ang);
        GLVERTEX2POLAR(rad_out, vt_ang);
    }
    glEnd();
}

#include <QPainter>

using namespace std;

#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE  0x809D
#endif

void HistogramViewer::draw(){
#if 0
    const float nbSteps = 200.0;

      glBegin(GL_QUAD_STRIP);
      for (int i=0; i<nbSteps; ++i)
        {
          const float ratio = i/nbSteps;
          const float angle = 21.0*ratio;
          const float c = cos(angle);
          const float s = sin(angle);
          const float r1 = 1.0 - 0.8f*ratio;
          const float r2 = 0.8f - 0.8f*ratio;
          const float alt = ratio - 0.5f;
          const float nor = 0.5f;
          const float up = sqrt(1.0-nor*nor);
          glColor3f(1.0-ratio, 0.2f , ratio);
          glNormal3f(nor*c, up, nor*s);
          glVertex3f(r1*c, alt, r1*s);
          glVertex3f(r2*c, alt+0.05f, r2*s);
        }
      glEnd();

#else
    //printf(" hv::draw()\n");

    if(text_id == 0)
        return;

    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glColor3f(0,1,0);

    glDisable(GL_TEXTURE_2D);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, text_id);

    //draw_horseshoe_old(0.2, 1.0, 0.0, M_PI, 36);
    //draw_horseshoe(0.2, 1.0, 0, M_PI, 36);
    //draw_horseshoe_old(0.2, 1.0, -M_PI, M_PI, 9);
    draw_horseshoe(0.2, 1.0, -M_PI/2, M_PI/2, 69);

    //return;
    /*
        // -- draw a plane!
    glBegin(GL_QUADS);

    glTexCoord2f(0,0);  glVertex2f(0,0);
    glTexCoord2f(0,1);  glVertex2f(0,1);
    glTexCoord2f(1,1);  glVertex2f(1,1);
    glTexCoord2f(1,0);  glVertex2f(1,0);

    glEnd();
*/
    {
        glColor3f(0.1,0.1,0.1);


        float mindel = 37.5*min_del/75000;
        float maxdel = 37.5*max_del/75000;

        //printf(" mindel = %f, maxdel = %f\n", mindel, maxdel);
        char tag[50], tag2[50];

        sprintf(tag, "Delta range [%.3f, %.3f]", mindel, maxdel);
        this->renderText(10,38,tag);

        sprintf(tag2, "Value range [0, %.3f]", max_val);
        this->renderText(10,52,tag2);
    }

    draw_scale(tfunc, 0, max_val);
#endif
}

#if 0
void HistogramViewer::drawOverpaint(QPainter *painter) {

    painter->save();
    painter->translate(width()/2, height()/2);

    QRadialGradient radialGrad(QPointF(-40, -40), 100);
    radialGrad.setColorAt(0, QColor(255, 255, 255, 100));
    radialGrad.setColorAt(1, QColor(200, 200, 0, 100));

    painter->setBrush(QBrush(radialGrad));

    painter->drawRoundRect(-100, -100, 200, 200);

    int n = 100;

    float x_start = -2*M_PI;
    float x_end = 2*M_PI;
    float x_step = (x_end - x_start) / (n-1);

    QPointF *line = new QPointF[n];
    for(int i = 0; i < n; i++){
        float x = x_start + (float)i*x_step;
        line[i] = QPointF(x, sin(x));
    }
    painter->drawPolyline(line, n);
    /*static const QPointF points[3] = {
        QPointF(10.0, 80.0),
        QPointF(20.0, 10.0),
        QPointF(80.0, 30.0),
    };
    painter->drawPolyline(points, 3);*/

    painter->restore();
}
#endif
/*
void HistogramViewer::paintEvent(QPaintEvent *event)
{
    Q_UNUSED(event)
    QPainter painter;
    painter.begin(this);
    painter.setRenderHint(QPainter::Antialiasing);

    // Save current OpenGL state
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    // Reset OpenGL parameters
    glShadeModel(GL_SMOOTH);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_MULTISAMPLE);
    static GLfloat lightPosition[4] = { 1.0, 5.0, 5.0, 1.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    qglClearColor(backgroundColor());

    // Classical 3D drawing, usually performed by paintGL().
    preDraw();
    draw();
    postDraw();

    // Restore OpenGL state
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glPopAttrib();

    drawOverpaint(&painter);

    painter.end();
}
*/
#include <QImage>

#if 0
void HistogramViewer::create_texture(const Atom &atom, const TransferFunction *tf){

    printf(" Creating texture...");
    fflush(stdout);

    this->tfunc = tf;
    const vector<vector<float> > &hist_2d_angle = atom.hist_2d_angle;

    uint H = 0;
    uint W = hist_2d_angle.size();

    for(uint w = 0; w < W; w++){

        float thismax = *std::max_element(hist_2d_angle[w].begin(), hist_2d_angle[w].end());
        this->max_val = std::max(max_val, thismax);
        H = std::max(H, (uint)hist_2d_angle[w].size());
    }

    this->min_del = atom.min_delta;
    this->max_del = atom.max_delta;

    // create a colored texture
    GLfloat *pat = new GLfloat[3*W*H];

    for(uint h = 0; h < H; h++){
    for(uint w = 0; w < W; w++){

        int pix = 3*(h*W + w);

        // scale the function value to bring it in [0,1]
        float grayscale_val = hist_2d_angle[w][h] / max_val;
        //float grayscale_val = (float)pix / (3.0*(float)W*H);

        tf->query(grayscale_val, pat[pix], pat[pix+1], pat[pix+2]);
    }
    }

    // now create opengl texture
    glEnable(GL_TEXTURE_2D);
    glGenTextures(1, &text_id);
    glBindTexture(GL_TEXTURE_2D, text_id);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

/*
    QImage img;
    if(!img.load("checkerboard.png")){
        std::cerr << "ERROR in loading image" << std::endl;
    }
    QImage t = QGLWidget::convertToGLFormat(img);

    //glTexImage2D(GL_TEXTURE_2D, 0, 3, t.width(), t.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, t.bits());
    //gluBuild2DMipmaps(GL_TEXTURE_2D, 3, t.width(), t.height(), GL_RGBA, GL_UNSIGNED_BYTE, t.bits());
*/
    gluBuild2DMipmaps(GL_TEXTURE_2D, 3, W, H, GL_RGB, GL_FLOAT, pat);
    //glTexImage2D(GL_TEXTURE_2D, 0, 3, W, H, 0, GL_RGB, GL_FLOAT, pat);

    printf(" Done! Texture id = %d.\n", text_id);
}
#endif

void HistogramViewer::create_texture(const char *filename){

    printf(" Creating texture...");
    fflush(stdout);

    // now create opengl texture
    glEnable(GL_TEXTURE_2D);
    glGenTextures(1, &text_id);
    glBindTexture(GL_TEXTURE_2D, text_id);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

    QImage img;
    if(!img.load(filename)){
        std::cerr << "ERROR in loading image" << std::endl;
    }
    QImage t = QGLWidget::convertToGLFormat(img);

    //glTexImage2D(GL_TEXTURE_2D, 0, 3, t.width(), t.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, t.bits());
    gluBuild2DMipmaps(GL_TEXTURE_2D, 3, t.width(), t.height(), GL_RGBA, GL_UNSIGNED_BYTE, t.bits());
}

void HistogramViewer::create_texture(const Atom &atom, const TransferFunction *tf){

    printf(" Creating texture...");
    fflush(stdout);

    this->tfunc = tf;
    const vector<statistics> &hist_2d_angle = atom.hist_2d_angle;

    max_val = 0;
    uint H = 0;
    uint W = hist_2d_angle.size();
    printf(" W = %d, H = %d, maxval = %f\n", W, H, max_val);

    for(uint w = 0; w < W; w++){

        const vector<float>& hist = hist_2d_angle[w].get_pdf();

        float thismax = *std::max_element( hist.begin(), hist.end() );
        //float thismax = *std::max_element( hist.freq.begin(), hist.freq.end() );
        this->max_val = std::max(max_val, thismax);
        H = std::max(H, (uint)hist.size());
    }

    this->min_del = atom.min_delta;
    this->max_del = atom.max_delta;

    printf("-- W = %d, H = %d, maxval = %f\n", W, H, max_val);

    // create a colored texture
    GLfloat *pat = new GLfloat[3*W*H];

    for(uint w = 0; w < W; w++){

        const vector<float> &hist = hist_2d_angle[w].get_pdf();

        for(uint h = 0; h < H; h++){

            int pix = 3*(h*W + w);

            // scale the function value to bring it in [0,1]
            float grayscale_val = (float)hist[h] / (float)max_val;
            //float grayscale_val = (float)pix / (3.0*(float)W*H);

            //printf(" %d -- grayscale = %f\n", hist.pdf[h], grayscale_val);
            tf->query(grayscale_val, pat[pix], pat[pix+1], pat[pix+2]);
        }
    }

    // now create opengl texture
    glEnable(GL_TEXTURE_2D);
    glGenTextures(1, &text_id);
    glBindTexture(GL_TEXTURE_2D, text_id);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

    gluBuild2DMipmaps(GL_TEXTURE_2D, 3, W, H, GL_RGB, GL_FLOAT, pat);
    //glTexImage2D(GL_TEXTURE_2D, 0, 3, W, H, 0, GL_RGB, GL_FLOAT, pat);

    printf(" Done! Texture id = %d.\n", text_id);
}

void HistogramViewer::draw_scale(const TransferFunction *tf, float min_v, float max_v){

    if(tfunc == 0)
        return;

    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_TEXTURE_2D);

    /// start drawing the box
    glPushAttrib( GL_ALL_ATTRIB_BITS );

    glMatrixMode( GL_PROJECTION );
    glPushMatrix();
    glLoadIdentity();

    glMatrixMode( GL_MODELVIEW );
    glPushMatrix();
    glLoadIdentity();

    glDisable(GL_DEPTH_TEST);

    static float v0 = 0.300f, v1 = 0.9f;
    static float u0 = 0.725f, u1 = 0.8f;

    static int steps = 20;

    /// -------------------------------------------------

    //double vrange = max_v - min_v;

    glBegin(GL_QUAD_STRIP);
    for(uint i = 0; i <= (steps); i++){

        float t = (float)i / (float)steps;
        //float grayscale = (float)i * vrange / (float)steps + min_v;
        float r, g, b;

        tf->query(t, r, g, b);

        glColor3f(r, g, b);
        glVertex3f( u0, Utils::lerp(v0,v1,t), 0 );
        glVertex3f( u1, Utils::lerp(v0,v1,t), 0 );
    }
    glEnd();

    /// ---------------------------------------------------
    glColor3f(0,0,0);
    glBegin(GL_LINE_LOOP);
        glVertex3f( u0, v0, 0 );
        glVertex3f( u1, v0, 0 );
        glVertex3f( u1, v1, 0 );
        glVertex3f( u0, v1, 0 );
    glEnd();

    float font_size = 12.0f;// / 750.0f * (float)this->size().height();

    //printf(" function range [%f, %f]\n", min_v, max_v);
    QFont font(QString("Courier New"),(int)font_size);

    glColor3f(0.1,0.1,0.1);

    for(int i = 0; i <= steps; i+=4){

        float t0 = (float)i/(float)(steps);
        float t1 = (float)i/(float)(steps+1);
        float fvalue = Utils::lerp( min_v, max_v, t0 );

        //cout<<" i = "<<i<<", t0 = "<<t0<<", t1 = "<<t1<<", val = "<<fvalue<<endl;
        this->renderText( u1+0.025f,  Utils::lerp(v0,v1,t1), 0, QString::number(fvalue), font );
        //this->renderText( u1+0.025f, lerp(v0,v1,t1), 0, QString("%1").arg(fvalue,0,'f',2), font );
    }

    glMatrixMode( GL_PROJECTION );
    glPopMatrix();

    glMatrixMode( GL_MODELVIEW );
    glPopMatrix();

    glPopAttrib();

    glEnable(GL_DEPTH_TEST);
    //if(live_enableLighting)
      //  glEnable(GL_LIGHTING);
}
