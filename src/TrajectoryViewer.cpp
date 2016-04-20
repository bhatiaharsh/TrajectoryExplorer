/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#include "Utils.h"
#include "TrajectoryExplorerApp.h"

/// -----------------------------------------------------------------
using namespace std;

/// -----------------------------------------------------------------

TrajectoryViewer::TrajectoryViewer(QWidget *parent) : QGLViewer(parent){

    parentApp = (TrajectoryExplorerApp*)parent;
    this->live_animation_time = 0;

    bbox_dims[0] = 1;
    bbox_dims[1] = 1;
    bbox_dims[2] = 1;

    if(parentApp != 0)
    {
        const float &lx = parentApp->mdata.lattice_vectors[0][0];
        const float &ly = parentApp->mdata.lattice_vectors[1][1];
        const float &lz = parentApp->mdata.lattice_vectors[2][2];

        box_x = 1, box_y = 1, box_z = 1;

        if(lx < ly && lx < lz){
            box_x = 1;          box_y = ly / lx;        box_z = lz / lx;
        }
        else if(ly < lx && ly < lz){
            box_x = lx / ly;    box_y = 1;              box_z = lz / ly;
        }
        else if(lz < lx && lz < ly){
            box_x = lx / lz;    box_y = ly / lz;        box_z = 1;
        }
        printf(" -- lattice = (%E, %E, %E), draw_box = (%f, %f, %f)\n", lx, ly, lz, box_x, box_y, box_z);
    }
}

void TrajectoryViewer::init(){

    this->setSceneBoundingBox(qglviewer::Vec(0,0,0), qglviewer::Vec(box_x, box_y, box_z));
    this->setSceneCenter(qglviewer::Vec(0.5*box_x,0.5*box_y,0.5*box_z));
    this->camera()->centerScene();
    this->camera()->showEntireScene();

    this->setFPSIsDisplayed(true);
    this->setTextIsEnabled(true);

    this->setBackgroundColor(Qt::white);
    this->setForegroundColor(Qt::black);

    this->setAnimationPeriod(1);        // milliseconds!

    cout << "\n -------------------------------------------------------- \n";
    cout << " | OpenGL Renderer: " << (char*) glGetString(GL_RENDERER) << "\n";
    cout << " | OpenGL Version: " << (char*) glGetString(GL_VERSION) << "\n";
    cout << " -------------------------------------------------------- \n\n";

    // create a precompiled list to draw spheres
    sphereList = glGenLists(1);
    glNewList(sphereList, GL_COMPILE);
        draw_sphere(1.0, 40, 40);
    glEndList();

    // init the viewer
    glEnable(GL_NORMALIZE);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    glDisable(GL_CULL_FACE);

    glEnable(GL_DEPTH_TEST);
    glDepthMask(true);
    glDepthFunc(GL_LESS);

    // lighting
    GLfloat ambientLight[] = {0.5f, 0.5f, 0.5f, 1.0f};
    GLfloat lightColor[] = {0.8f, 0.8f, 0.8f, 1.0f};

    GLfloat lightPos[] = {1.0, 1.0, 1.0, 1.0 };

    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientLight);

    //Diffuse (non-shiny) light component
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);

    //Specular (shiny) light component
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightColor);
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos);

    //Disable color materials, so that glMaterial calls work
    glDisable(GL_COLOR_MATERIAL);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

  // Restore previous viewer state.
  //restoreStateFromFile();

    // key bindings for TrajectoryViewer
    setKeyDescription(Qt::Key_Plus, "Trajectory Viewer: Increment Animation Time");
    setKeyDescription(Qt::Key_Minus, "Trajectory Viewer: Decrement Animation Time");

    // key bindings for ACViewer
    setKeyDescription(Qt::Key_S, "Statistics Viewer: Save Graphs");
    setKeyDescription(Qt::Key_E, "Statistics Viewer: Toggle Data Error");
    setKeyDescription(Qt::Key_O, "Statistics Viewer: Toggle Orientation of 2D plot");
    setKeyDescription(Qt::Key_I, "Statistics Viewer: Toggle Interpolation in 2D plot");
    setKeyDescription(Qt::Key_L, "Statistics Viewer: Toggle Legend");
    setKeyDescription(Qt::Key_R, "Statistics Viewer: Reset Axes");


    // Add custom mouse bindings description (see mousePressEvent())
    setMouseBindingDescription(Qt::ShiftModifier, Qt::LeftButton, "Statistics Viewer: Move Selector in 2D plot");

    // Opens help window
    help();
}

QString TrajectoryViewer::helpString() const {

    QString hs ("AtomTrajectoryViewer v1.0 contains of two components:\n");
    hs.append("(a) Trajectory Viewer displays static or animated trajectories,\n");
    hs.append("(b) Statistics Viewer displays relative-angle statistics.\n\n");
    hs.append(" Both (a) and (b) can be synchronized through the checkbox on the main UI.\n");
    hs.append(" Please note the keyboard and mouse bindings in the help window.");
    hs.append(" Main Menu (top) allows changing the transfer function in Statistics Visualization.\n");
    hs.append(" Main Menu (top) allows setting 3D visualthe transfer function in Statistics Visualization.\n");
    return hs;
}

void TrajectoryViewer::keyPressEvent(QKeyEvent* event){

    switch(event->key()){
       /* case 'N':   parentApp->ui.groupBox_anim->setChecked( !parentApp->ui.groupBox_anim->isChecked() );
                    break;*/

        case '+':   live_animation_time = (live_animation_time == parentApp->mdata.num_tsteps - 1) ? 0 : live_animation_time+1;
                    break;

        case '-':   live_animation_time = (live_animation_time == 0) ? parentApp->mdata.num_tsteps - 1 : live_animation_time-1;
                    break;

    case 'R':
        this->camera()->centerScene();
        this->camera()->showEntireScene();
    }
    updateGL();
    QGLViewer::keyPressEvent(event);
}

void TrajectoryViewer::animate(){
    live_animation_time = (live_animation_time + parentApp->ui.sb_animSpeed->value()) % parentApp->mdata.num_tsteps;
}

/// -----------------------------------------------------------------

#define PI      3.14159265359
#define PIOVER2 1.570796326795
#define DE2RA   0.01745329252
#define RA2DE   57.2957795129

void getGlobeMapping(double *s, double *t, double *q, double phi[2], double theta[2]){
   s[0] = (180.0 - phi[0]*RA2DE)/360.0;   s[1] = (theta[0] + PIOVER2)*RA2DE/180.0;
   t[0] = s[0];                           t[1] = (theta[1] + PIOVER2)*RA2DE/180.0;
   q[0] = (180.0 - phi[1]*RA2DE)/360.0;   q[1] = t[1];
}
void get3DCoords(double *a, double rad, double phi, double theta){
   a[0] = rad * cos(phi) * cos(theta);    //x
   a[1] = rad * sin(theta);               //y
   a[2] = rad * sin(phi) * cos(theta);    //z
}
void getNormalVector(double *u, double *v, double *w, double *n){
   n[0] = (u[0] + v[0] + w[0]) * 0.33;// - live_object_xz_trans[0] ;
   n[1] = (u[1] + v[1] + w[1]) * 0.33;// - live_object_y_trans ;
   n[2] = (u[2] + v[2] + w[2]) * 0.33;// + live_object_xz_trans[1];
}

void TrajectoryViewer::draw_sphere(double R, double NumLatitudes, double NumLongitudes){

    glBegin(GL_TRIANGLES);

    double start_lat = -90;
    double start_lon = 0.0;
    double lat_incr = 180.0 / NumLatitudes;
    double lon_incr = 360.0 / NumLongitudes;

    double u[3], v[3], w[3], n[3];
    double s[2], t[2], q[2];

    double phi[2], theta[2];

    for (int col = 0; col < NumLongitudes; col++){

        phi[0] = (start_lon + col * lon_incr) * DE2RA;
        phi[1] = (start_lon + (col + 1) * lon_incr) * DE2RA;

        for (int row = 0; row < NumLatitudes; row++){

            theta[0] = (start_lat + row * lat_incr) * DE2RA;
            theta[1] = (start_lat + (row + 1) * lat_incr) * DE2RA;

            get3DCoords(u, R, phi[0], theta[0]);
            get3DCoords(v, R, phi[0], theta[1]);
            get3DCoords(w, R, phi[1], theta[1]);

            getGlobeMapping(s, t, q, phi, theta);
            getNormalVector(u, v, w, n);

            glNormal3dv(n);

            glVertex3dv(u);
            glVertex3dv(v);
            glVertex3dv(w);

            get3DCoords(v, R, phi[1], theta[0]);

            getGlobeMapping(s, t, q, phi, theta);
            getNormalVector(u, v, w, n);

            // correction from the function!
            t[0] = (180.0 - phi[1]*RA2DE)/360.0;
            q[1] = (theta[0] + PIOVER2)*RA2DE/180.0;

            glNormal3dv(n);

            glVertex3dv(u);
            glVertex3dv(w);
            glVertex3dv(v);
        }
    }

    glEnd();
}
void TrajectoryViewer::draw_bbox(const float dims[], QColor col, float width){

    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);

    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    glLineWidth(width);
    glColor3f(col.redF(), col.greenF(), col.blueF());

    glBegin(GL_LINE_STRIP);

    glVertex3f(0.0f, 0.0f, 0.0f);           // xmin, ymin, zmin (0)
    glVertex3f(dims[0], 0.0f, 0.0f);        // xmax, ymin, zmin (1)
    glVertex3f(dims[0], dims[1], 0.0f);     // xmax, ymax, zmin (2)
    glVertex3f(0.0f, dims[1], 0.0f);        // xmin, ymax, zmin (3)
    glVertex3f(0.0f, 0.0f, 0.0f);           // xmin, ymin, zmin (0)

    glVertex3f(0.0f, 0.0f, dims[2]);        // xmin, ymin, zmax (4)
    glVertex3f(dims[0], 0.0f, dims[2]);     // xmax, ymin, zmax (5)
    glVertex3f(dims[0], dims[1], dims[2]);  // xmax, ymax, zmax (6)
    glVertex3f(0.0f, dims[1], dims[2]);     // xmin, ymax, zmax (7)
    glVertex3f(0.0f, 0.0f, dims[2]);        // xmin, ymin, zmax (4)

    glEnd();

    glBegin(GL_LINES);

    glVertex3f(dims[0], 0.0f, 0.0f);        // xmax, ymin, zmin (1)
    glVertex3f(dims[0], 0.0f, dims[2]);     // xmax, ymin, zmax (5)

    glVertex3f(dims[0], dims[1], 0.0f);     // xmax, ymax, zmin (2)
    glVertex3f(dims[0], dims[1], dims[2]);  // xmax, ymax, zmax (6)

    glVertex3f(0.0f, dims[1], 0.0f);        // xmin, ymax, zmin (1)
    glVertex3f(0.0f, dims[1], dims[2]);     // xmin, ymax, zmax (5)

    glEnd();

    glLineWidth(1);
    glPopAttrib();
}

void TrajectoryViewer::draw_tet(const std::set<uint> &bonded_atoms){

    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);

#if 0
    glEnable(GL_LIGHTING);

    std::vector<const float*> oatoms;
    set<uint>::const_iterator it = bonded_atoms.begin();
    for(; it != bonded_atoms.end(); it++){

        int atom1_idx = parentApp->find_atom(*it);
        if(atom1_idx == -1)
            continue;
        oatoms.push_back( parentApp->atoms[atom1_idx].pos[live_animation_time] );
    }

    if(oatoms.size() != 4)
        return;

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);

    glColor4f(0.5, 0.5, 0.5, 0.5);
    glBegin(GL_TRIANGLES);

        for(uint i = 0; i < 4; i++){

            glVertex3f(oatoms[i][0], oatoms[i][1], oatoms[i][2]);
            glVertex3f(oatoms[ (i+1)%4 ][0], oatoms[ (i+1)%4 ][1], oatoms[ (i+1)%4 ][2]);
            glVertex3f(oatoms[ (i+2)%4 ][0], oatoms[ (i+2)%4 ][1], oatoms[ (i+2)%4 ][2]);
        }

    glEnd();

    glDisable(GL_BLEND);
#else
    glLineStipple(2, 0xAAAA);
    //glEnable( GL_LINE_STIPPLE );

    glLineWidth(1);

    QColor col = Qt::darkGray;
    glColor3f(col.redF(), col.greenF(), col.blueF());

    set<uint>::const_iterator it = bonded_atoms.begin();
    for(; it != bonded_atoms.end(); it++){

        int atom1_idx = parentApp->find_atom(*it);
        if(atom1_idx == -1)
            continue;

        const float *p = parentApp->atoms[atom1_idx].pos[live_animation_time];

        set<uint>::const_iterator it2 = it;
        it2++;
        for(; it2 != bonded_atoms.end(); it2++){

            int atom2_idx = parentApp->find_atom(*it2);
            if(atom2_idx == -1)
                continue;

            const float *q = parentApp->atoms[atom2_idx].pos[live_animation_time];

            glBegin(GL_LINES);
            glVertex3f(p[0],p[1],p[2]);
            glVertex3f(q[0],q[1],q[2]);
            glEnd();
        }
    }
    glDisable( GL_LINE_STIPPLE );
#endif
    glPopAttrib();
}

void TrajectoryViewer::draw_bond(float p[], float q[]){

    glPushAttrib(GL_LIGHTING_BIT);
    QColor col = Qt::darkGray;

#ifdef USE_GLE
    // draw bonds as tubes!
    if(parentApp->ui.cb_gle_bonds->isChecked() && parentApp->ui.cb_gle_bonds->isEnabled()){

        gleDouble linearray[4][3];
        gleColor colorarray[4];

        for(uint d = 0; d < 3; d++){

            float dir_vec = q[d]-p[d];

            linearray[0][d] = p[d] - 0.2*dir_vec;        // first end point
            linearray[3][d] = q[d] + 0.2*dir_vec;        // end point

            linearray[1][d] = p[d];     // first bond point
            linearray[2][d] = q[d];     // second bond point
        }

        for(uint i = 0; i < 4; i++){
            colorarray[i][0] = col.redF();
            colorarray[i][1] = col.greenF();
            colorarray[i][2] = col.blueF();
        }

        glEnable(GL_LIGHTING);
        glEnable(GL_COLOR_MATERIAL);

        glePolyCylinder(4, linearray, colorarray, 0.002);

        glDisable(GL_LIGHTING);
        glDisable(GL_COLOR_MATERIAL);

        glPopAttrib();
        return;
    }
#endif
    // ------------------------------------------
    // draw line!

    glDisable(GL_LIGHTING);

    glLineWidth(2.5);
    glColor3f(col.redF(), col.greenF(), col.blueF());
    glBegin(GL_LINES);
    glVertex3f(p[0],p[1],p[2]);
    glVertex3f(q[0],q[1],q[2]);
    glEnd();

    glPopAttrib();
}


void TrajectoryViewer::draw(){

    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_COLOR_MATERIAL);

    // -----------------------------------------------------------
    glDisable(GL_LIGHTING);

    glScalef(box_x, box_y, box_z);

    if(parentApp->ui.cb_bbox->isChecked()){
        draw_bbox(bbox_dims);
    }

    // -----------------------------------------------------------
    set<int> atomsToDisplay;
#if 0
    for(uint i = 1; i <= 638; ++i)
        atomsToDisplay.insert(i);
#else
    QList<QListWidgetItem*> selectedItems = parentApp->ui.lw_atoms_li->selectedItems();
    for(uint i = 0; i < selectedItems.size(); i++){
        atomsToDisplay.insert( TrajectoryExplorerApp::get_atomID_from_UI(selectedItems[i]) );
    }
    selectedItems = parentApp->ui.lw_atoms_ox->selectedItems();
    for(uint i = 0; i < selectedItems.size(); i++){
        atomsToDisplay.insert( TrajectoryExplorerApp::get_atomID_from_UI(selectedItems[i]) );
    }
#endif
    // -----------------------------------------------------------
    // draw the static trajectory!
    if(parentApp->ui.groupBox_static->isChecked()){

        for(uint i = 0; i < parentApp->atoms.size(); i++){

            const Atom &curratom = parentApp->atoms[i];
            if(atomsToDisplay.find(curratom.mid) == atomsToDisplay.end())
                continue;

            int delta_idx = (parentApp->ui.sb_delta->value() / parentApp->mdata.time_step);

            if(delta_idx == 1){
                draw_atomTrace(curratom, parentApp->ui.sb_offset->value(), parentApp->mdata.num_tsteps);
            }
            else{
                draw_linestrip(curratom, parentApp->ui.sb_offset->value(), parentApp->mdata.num_tsteps, delta_idx);
            }
        }
        return;
    }

    // -----------------------------------------------------------
    // draw the atoms at current animation time

    // -------------------------------
    for(uint i = 0; i < parentApp->atoms.size(); i++){

        const Atom &curratom = parentApp->atoms[i];
        if(atomsToDisplay.find(curratom.mid) == atomsToDisplay.end())
            continue;

        draw_atomPosition(curratom, live_animation_time);

        int tailLength = parentApp->ui.sb_tailLength->value();
        int start_time = (tailLength == -1) ? 0 : std::max(0, (int)live_animation_time-(int)tailLength);

        draw_atomTrace(curratom, start_time, live_animation_time);
    }

    // draw bonds
    if(parentApp->ui.gb_bonds->isChecked()){

        for(uint i = 0; i < parentApp->track_li.size(); i++){

            const Atom &liatom = parentApp->atoms[i];

            if(atomsToDisplay.find(liatom.mid) == atomsToDisplay.end())
                continue;

            if(liatom.bonded_atoms.empty())
                continue;

            set<uint>::const_iterator it = liatom.bonded_atoms[live_animation_time].begin();
            for(; it != liatom.bonded_atoms[live_animation_time].end(); it++){

                if(atomsToDisplay.find(*it) == atomsToDisplay.end())
                    continue;

                int atom_idx = parentApp->find_atom(*it);
                if(atom_idx == -1)
                    continue;

                const Atom &otheratom = parentApp->atoms[atom_idx];
                draw_bond(liatom.pos[live_animation_time], otheratom.pos[live_animation_time]);
            }
        }
    }

    // draw tets
    if(parentApp->ui.cb_tets->isChecked()){

        for(uint i = 0; i < parentApp->track_li.size(); i++){

            const Atom &liatom = parentApp->atoms[i];

            if(atomsToDisplay.find(liatom.mid) == atomsToDisplay.end())
                continue;

            if(liatom.bonded_atoms.empty())
                continue;

            draw_tet(liatom.bonded_atoms[live_animation_time]);
        }
    }


    {   // ----- write the current time

        glColor3f(0.1,0.1,0.1);

        QString tag("Phys Time: ");
        QString tag2("Time Indx: ");

        tag2.append(QString::number(live_animation_time));

        double time_ps = (double)live_animation_time * parentApp->mdata.time_step;

        char time[25];
        sprintf(time, "%2.3f %s", time_ps, parentApp->mdata.time_unit.c_str());
        tag.append(time);

        this->renderText(10,38,tag2);
        this->renderText(10,54,tag);
    }
}



/// ------------------------------------------------------------------
/// draw atom trace
void TrajectoryViewer::draw_atomTrace(const Atom &atom, uint start_time, uint end_time){

    if(parentApp->ui.cb_gle_tubes->isChecked() && parentApp->ui.cb_gle_tubes->isEnabled()) {
        draw_tube(atom, start_time, end_time);
    }
    else {
        draw_linestrip(atom, start_time, end_time);
    }
}


/// ------------------------------------------------------------------
/// draw line strip

void TrajectoryViewer::draw_linestrip(const Atom &atom, uint start_time, uint end_time, uint step){

    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(2.0);

    QColor col1 = Qt::lightGray;
    QColor col2 = atom.draw_color;

    glBegin(GL_LINE_STRIP);
    for(uint t = start_time; t < end_time; t += step){

        float l = (float) t / (float) atom.pos.size();

        QColor col = Utils::lerp(col1, col2, l);

        glColor3f(col.redF(), col.greenF(), col.blueF());
        glVertex3f(atom.pos[t][0], atom.pos[t][1], atom.pos[t][2]);

        if(atom.is_jump(t)){
            glEnd();
            glBegin(GL_LINE_STRIP);
        }
    }

    glEnd();
    glLineWidth(1.0);
    glPopAttrib();
}

/// ------------------------------------------------------------------
/// draw tube

void TrajectoryViewer::draw_tube(const Atom &atom,  uint start_time, uint end_time){

    glPushAttrib(GL_LIGHTING_BIT);

#ifdef USE_GLE
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);

    float rad = parentApp->ui.dsb_tubeThickness->value();

    uint trc_start_time = 0;
    for(uint trc_num = 0; trc_num < atom.trace_endIdx.size(); trc_num++){

        uint trc_end_time = atom.trace_endIdx[trc_num];
        uint trc_len = trc_end_time-trc_start_time+1;

        // check if there is no overlap between this trace and the range to be drawn!
        if(end_time < trc_start_time)
            break;

        if(start_time > trc_end_time){
            trc_start_time = trc_end_time;
            continue;
        }

        // now, find the overlap!
        int draw_start_idx = std::max(start_time, trc_start_time) - trc_start_time;
        int draw_end_idx = std::min(end_time, trc_end_time) - trc_start_time;

        if(draw_end_idx < 0 || draw_start_idx < 0 || draw_start_idx >= trc_len || draw_end_idx >= trc_len){
            printf(" something wrong!\n");
            continue;
        }

        // now, draw!
        glePoint* points = atom.tubePoints[trc_num];
        gleColor* colors = atom.tubeColors[trc_num];

        glePolyCylinder( (draw_end_idx-draw_start_idx+1), points+draw_start_idx, colors+draw_start_idx, rad);

        trc_start_time = trc_end_time;
    }

    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
#else

    draw_linestrip(atom, start_time, end_time);
#endif
    glPopAttrib();
}

/// --------------------------------------------------------------
/// draw atom position

void TrajectoryViewer::draw_atomPosition(const Atom &atom, const uint time){

    glPushAttrib(GL_LIGHTING_BIT);

    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);

    // material color
    QColor col = atom.draw_color;
#if 0
    //The color of the sphere
    GLfloat materialColor[] = {col.redF(), col.greenF(), col.blueF(), 1.0f};

    //The specular (shiny) component of the material
    //GLfloat materialSpecular[] = {1.5*col.redF(), 1.5*col.greenF(), 1.5*col.blueF()};
    //GLfloat materialSpecular[] = {0,0,1,1};

    //The color emitted by the material
    //GLfloat materialEmission[] = {0.0f,1.0f,0, 1.0f};

    //GLfloat shininess = 15;

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, materialColor);
    //glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, materialSpecular);
    //glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, materialEmission);
    //glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess); //The shininess parameter
#endif

    glColor3f(col.redF(), col.greenF(), col.blueF());

    // --- draw now!
    float sphere_rad = parentApp->ui.dsb_atomRadius->value();
    const float *pos = atom.pos[time];

    glPushMatrix();
    glTranslatef(pos[0], pos[1], pos[2]);
    glScalef(sphere_rad/box_x,sphere_rad/box_y,sphere_rad/box_z);
    glCallList(sphereList);
    glPopMatrix();

    glPopAttrib();
}
