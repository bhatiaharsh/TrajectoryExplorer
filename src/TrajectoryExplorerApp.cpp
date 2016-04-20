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
#include "ACData.h"
#include <QVector>
#include "bootstrapping.h"
#include "ACWindow.h"

using namespace std;

TrajectoryExplorerApp::TrajectoryExplorerApp(QSettings &settings) {

    appsettings = &settings;

    // -------------------------------------------------------------
    // --- read the atom ids

    QString lfile = appsettings->value("data/LI_ID").toString();
    QString ofile = appsettings->value("data/OX_ID").toString();

    if(!lfile.isEmpty()) {      Utils::read_atomIDs(lfile.toLatin1().data(), track_li);   }
    if(!ofile.isEmpty()) {      Utils::read_atomIDs(ofile.toLatin1().data(), track_ox);   }

    if(track_li.empty()){
        printf(" No Lithiums specified in the system!\n");
        exit(1);
    }
    if(track_ox.empty()){
        printf(" No Oxygen specified in the system!\n");
        exit(1);
    }

    // --- read the atomic trajectories!
    QString pfile = appsettings->value("data/TRAJECTORIES").toString();
    if(pfile.isEmpty()){
        printf(" AtomTrajectoryApp needs a XDATCAR file!\n");
        exit(1);
    }

    vector<vector<float*> > rd_pos;
    RW_VASP::read_XDATCAR_all(pfile.toStdString(), mdata, rd_pos);//, 10);//, 20000);
    if(rd_pos.empty()){
        printf(" No positions read!\n");
        exit(1);
    }

    if(mdata.materials.size() == 0){
    mdata.materials.push_back(std::make_pair("C", 189));
    mdata.materials.push_back(std::make_pair("H", 252));
    mdata.materials.push_back(std::make_pair("O", 189));
    mdata.materials.push_back(std::make_pair("P", 1));
    mdata.materials.push_back(std::make_pair("F", 6));
    mdata.materials.push_back(std::make_pair("Li", 1));
    }

    //track_li.push_back(631);
    // -------------------------------------------------------------
    printf(" Processing atoms to track... ");
    fflush(stdout);

    // track all
    if(appsettings->value("data/TRACK_ALL").toBool()){
        for(uint i = 0; i < mdata.num_atoms; i++){

            uint atomid = i+1;
            atoms.push_back( Atom(atomid, Utils::get_color(atomid, mdata.find_material_type(atomid), 1)) );
        }
    }
    // track only specified atoms
    else {
        for(uint i = 0; i < track_li.size(); i++){

            uint atomid = track_li[i];
            atoms.push_back( Atom(atomid, Utils::get_color(atomid, mdata.find_material_type(atomid), 1)) );
        }
        for(uint i = 0; i < track_ox.size(); i++){

            uint atomid = track_ox[i];

            QColor c = Utils::get_color(atomid, mdata.find_material_type(atomid), 2);   // random color
#if 1
            switch(atomid){
                case 504:   c = QColor(Qt::darkRed);        break;
                case 516:   c = QColor(Qt::darkGreen);      break;
                case 525:   c = QColor(Qt::darkGray);       break;
                case 552:   c = QColor(Qt::darkYellow);     break;
                case 582:   c = QColor(Qt::darkCyan);       break;
                case 603:   c = QColor(Qt::darkBlue);       break;
            }
#endif
            atoms.push_back( Atom(atomid, c) );
        }
    }

    // -------------------------------------------------------------
    // copy positions!
    for(uint i = 0; i < atoms.size(); i++){

        Atom &atom = atoms[i];
        uint mid = atom.mid;
        atom.pos.resize(mdata.num_tsteps);

        for(uint t = 0; t < mdata.num_tsteps; t++){

            // VASP ids start from 1, but C++ ids start from 0
            atom.pos[t] = rd_pos[t][mid-1];
            for(uint d = 0; d < 3; d++){
                atom.pos[t][d] = Utils::wrap_around(atom.pos[t][d]);
            }
        }
    }
    printf(" Done! Tracking %d atoms (%d lithiums and %d others)!\n", atoms.size(), track_li.size(), track_ox.size());

    // ------------------------------------------------------------
    // set the center of the box and compute corresponding "jumps" across the domain
    current_center.first = -1;
    current_center.second = -1;
    this->find_jumps();
    this->create_tubes();


    // -------------------------------------------------------------
    // initialize ui and viewer
    ui.setupUi(this);
    viewer = new TrajectoryViewer(this);
    setCentralWidget(viewer);
    setWindowTitle("AtomTrajectoryViewer");

    // -------------------------------------------------------------

    float dist_thresh = appsettings->value("bonds/DIST_THRESH").toFloat();
    ui.dsb_bondDist->setValue(dist_thresh);
    ui.pb_computeBonds->click();

    // -------------------------------------------------------------
    // add Li ids to the UI to be selected
    {
        for(uint i = 0; i < track_li.size(); i++){
            ui.lw_atoms_li->addItem(QString::number(track_li[i]));
        }
        ui.lw_atoms_li->selectAll();
        char tag[20];
        sprintf(tag, "Select Li Atoms (%d)", ui.lw_atoms_li->count());
        this->ui.label_5->setText(tag);
    }

    // add Ox ids to the UI to be selected
    this->update_LiOxLabels();

    // -------------------------------------------------------------
    // initially, animation is off
    ui.groupBox_anim->setChecked(true);
    ui.groupBox_static->setChecked(false);

    ui.sb_time->setMaximum(mdata.num_tsteps-1);
    ui.sb_time->setValue(0);

    ui.sb_delta->setDecimals(4);
    ui.sb_delta->setSingleStep(this->mdata.time_step);
    ui.sb_delta->setMinimum(mdata.time_step);
    ui.sb_delta->setMaximum(0.5*mdata.time_step*(mdata.num_tsteps-1));

    ui.cb_center->addItem(QString::number(-1));
    for(uint i = 0; i < atoms.size(); i++){
        ui.cb_center->addItem(QString::number(atoms[i].mid));
        ui.cb_hist->addItem(QString::number(atoms[i].mid));
    }

#ifdef USE_GLE
    ui.cb_gle_bonds->setEnabled(true);
    ui.cb_gle_tubes->setEnabled(true);
    ui.dsb_tubeThickness->setEnabled(true);
#else
    ui.cb_gle_bonds->setEnabled(false);
    ui.cb_gle_tubes->setEnabled(false);
    ui.dsb_tubeThickness->setEnabled(false);
#endif

#if 0
    {
        center_at(638, 0);
        this->find_jumps();
        this->create_tubes();

        const Atom &atom = atoms[find_atom(638)];

        atom.write(0, 75000, "LiPositions.txt");
        atom.write(30000, 31226, "LiPositions-sub.txt");
    }
    exit(1);
#endif
    show();


}


void TrajectoryExplorerApp::update_LiOxLabels(){

    ui.lw_atoms_ox->clear();

    QList<QListWidgetItem*> selectedItems = ui.lw_atoms_li->selectedItems();
    uint count = 0;

    for(uint i = 0; i < selectedItems.size(); i++){

        uint liId = get_atomID_from_UI(selectedItems[i]);
        const Atom &liAtom = atoms[find_atom(liId)];

        for(set<uint>::const_iterator it = liAtom.bonded_atoms_list.begin(); it != liAtom.bonded_atoms_list.end(); it++){

            uint mid = *it;
            QColor col = atoms[find_atom(mid)].draw_color;
            //ui.lw_atoms_ox->addItem(QString::number(*it).append(" (").append(Utils::get_color_tag(col)).append(")"));
            ui.lw_atoms_ox->addItem(QString::number(*it));
            count++;
        }
    }

    ui.lw_atoms_ox->selectAll();
    char tag[20];
    sprintf(tag, "Select Ox Atoms (%d)", ui.lw_atoms_ox->count());
    this->ui.label_8->setText(tag);
}

void TrajectoryExplorerApp::compute_distance_bonds(float dist_thresh){

    printf(" Computing bonds using distance threshold = %f Angs...", dist_thresh);
    fflush(stdout);

    // get distance in absolute cartesian coordinates
    double sq_dist_thresh = dist_thresh*dist_thresh*powf(10,-20);

    set<uint> bonded_atoms;
    for(uint l = 0; l < track_li.size(); l++){

        Atom &liatom = atoms[this->find_atom(track_li[l])];

        liatom.bonded_atoms.clear();
        liatom.bonded_atoms.resize(mdata.num_tsteps);

        liatom.bonded_atoms_list.clear();

        for(uint i = 0; i < track_ox.size(); i++){

            const Atom &oxatom = atoms[this->find_atom(track_ox[i])];
            liatom.create_bonds(oxatom, mdata, sq_dist_thresh);//, 1);
        }

        bonded_atoms.insert(liatom.bonded_atoms_list.begin(), liatom.bonded_atoms_list.end());
    }

    printf(" Done! Found %d bonded atoms!\n", bonded_atoms.size());
}

bool TrajectoryExplorerApp::center_at(int mid, int time){

    if(current_center.first == mid && current_center.second == time){
        printf(" Box is already centered at atom %d time %d!\n");
        return false;
    }

    printf(" Centering the box at atom %d, time %d...", mid, time);
    fflush(stdout);

    // transvector for this shift
    float transVector[3] = {0,0,0};

    // for every shift, just keep accumulating to the transVector
    // so I can reset the box whenever I want
    static float accumulated_transVector[3] = {0,0,0};

    // i need to reset!
    if(mid == -1){

        for(uint d = 0; d < 3; d++){
            transVector[d] = -1*accumulated_transVector[d];
            accumulated_transVector[d] = 0;
        }
    }

    // transvector for centering this atom
    else {

        const Atom &atom = atoms[find_atom(mid)];

        for(uint d = 0; d < 3; d++){
            transVector[d] = 0.5*viewer->bbox_dims[d] - atom.pos[time][d];
            accumulated_transVector[d] += transVector[d];
        }
    }

    //printf(" transvector = (%f,%f,%f)\n", transVector[0], transVector[1], transVector[2]);
    //printf(" accumulated transvector = (%f,%f,%f)\n", accumulated_transVector[0], accumulated_transVector[1], accumulated_transVector[2]);

    // now shift the atoms to center!
    for(uint i = 0; i < atoms.size(); i++){

        Atom &atom = atoms[i];
        for(uint t = 0; t < mdata.num_tsteps; t++){

            for(uint d = 0; d < 3; d++){
                atom.pos[t][d] = Utils::wrap_around(atom.pos[t][d] + transVector[d]);
            }
        }

        atom.clear_Jumps_and_Tubes();
    }

    current_center.first = mid;
    current_center.second = time;

    printf(" Done!\n");
    return true;
}

void TrajectoryExplorerApp::find_jumps(){

    printf(" Finding jumps...");
    fflush(stdout);

    for(uint i = 0; i < atoms.size(); i++){
        atoms[i].find_Jumps();
    }
    printf(" Done!\n");
}

void TrajectoryExplorerApp::create_tubes(){

#ifdef USE_GLE
    // now create the tubes!
    printf(" Creating GLE tubes...");
    for(uint i = 0; i < atoms.size(); i++){
        atoms[i].create_Tubes(0);
    }
    printf(" Done!\n");
#endif
}

void TrajectoryExplorerApp::on_pb_hists_clicked(){

    uint atomid = ui.cb_hist->currentText().toInt();

    // -------------------------------------------------
    // if a histviewer corresponding to this atom exists, just show it!

    map<uint,ACViewer*>::iterator iter = histViewers.find(atomid);
    if(iter != histViewers.end()){
        iter->second->show();
        return;
    }

    // -------------------------------------------------
    // else, I need to create raw data

    Atom &atom = atoms[find_atom(atomid)];

    /*static uint offset = 15000;

    // min and max delta
    static uint min_delta = 10;
    static uint max_delta = 0.5*(this->mdata.num_tsteps-offset);
    static uint del_delta = 10;*/


    // skip equilibration
    static uint offset = appsettings->value("angles/OFFSET").toUInt();

    if(offset >= this->mdata.num_tsteps){
        printf(" Error: Offset (%d) is greater than the number of time steps (%d)\n", offset, this->mdata.num_tsteps);
        printf(" Cannot create angle-correlation plots!\n");
        return;
    }

    // min and max delta
    static uint min_delta = appsettings->value("angles/MIN_TSCALE").toUInt();
    static uint del_delta = appsettings->value("angles/DEL_TSCALE").toUInt();
    static uint max_delta = 0.5*(this->mdata.num_tsteps-offset);


    uint num_delta = (max_delta-min_delta)/del_delta;

    printf(" Computing angles for atom %d...", atomid);
    fflush(stdout);

    atom.hist_2d_angle.resize(num_delta);
    atom.hist_2d_times.resize(num_delta);
    QVector<double> tscales(num_delta);

    for(uint d = 0; d < num_delta; d++){

        uint delta = min_delta + d*del_delta;
        tscales[d] = (float)delta*mdata.time_step;

        atom.compute_angle_correlation(atom.hist_2d_angle[d], atom.hist_2d_times[d], delta, offset, 1);
        for(uint t = 0; t < atom.hist_2d_times[d].size(); t++)
            atom.hist_2d_times[d][t] *= mdata.time_step;

        //printf("[%d] = %d %f -- %d\n", d, delta, (float)delta*mdata.time_step, atom.hist_2d_times[d].size());
#if 0
        if(d > 0.14*num_delta)
            continue;

        //if(d > 0.14*num_delta)
        {
            QString filename = QString::number(atom.mid);
            filename.append("_bs_")
                    .append(QString::number(delta))
                    .append(".raw");

            BootStrap::compute_and_write(filename, 1000, 0.7, atom.hist_2d_angle[d]);
        }

        /*if(d == 0){
            for(uint ssp = 40; ssp < 100; ssp +=10){
                BootStrap::compute_and_write_forconvergence(delta, 180, 10000, ssp, atom.hist_2d_angle_new[d]);
            }
        }*/

        //printf(" delta = %d : sz = %d\n", delta, atom.hist_2d_angle_new[d].size());
#endif
    }

    printf(" Done! %d tscales in range (%.5f, %.5f) with res %.5f -- [%d %d) %d\n", num_delta,
           tscales.front(), tscales.back(), tscales[1]-tscales[0], min_delta, max_delta, del_delta);


    //return;

    //QVector<double> &h = atom.hist_2d_angle_new.front();
    //Utils::write_binary("atomtrajectory0.raw", h.data(), h.size());
    //return;

    // -------------------------------------------------------
    // now, create the plots!

    ACWindow *nac = new ACWindow(atomid, &atom.hist_2d_angle, &atom.hist_2d_times, tscales,
                                 appsettings->value("error_vis/DIRECTORY").toString());

    nac->resize(800,800);
    nac->show();

    QObject::connect(nac->viewer, SIGNAL(signal_selectorMoved(uint,float,float)), this, SLOT(slot_selectorChanged(uint,float,float)));

    // color map
    //ACViewer::TFUNC tf = get_selectedTransferFunction();

    //ACWindow *nac = new ACWindow(appsettings->value("error_vis/DIRECTORY").toString());
    //nac->viewer->set_graphStyling(tf);
    //nac->viewer->set_rawData(atomid, &atom.hist_2d_angle, &atom.hist_2d_times, tscales);

    return;
/*
    ACViewer *acviewer = new ACViewer(appsettings->value("error_vis/DIRECTORY").toString(), 0);
    acviewer->set_graphStyling(tf);
    acviewer->set_rawData(atomid, &atom.hist_2d_angle, &atom.hist_2d_times, tscales);

    acviewer->resize(800, 800);
    acviewer->show();
    histViewers.insert(std::pair<uint,ACViewer*> (atomid, acviewer));

    // -------------------------------------------------------
    // create a linked view
    QObject::connect(acviewer, SIGNAL(signal_selectorMoved(uint,float,float)), this, SLOT(slot_selectorChanged(uint,float,float)));
    */
}

// ---------------------------------------------
// #TODO! i shouldnt be using global variable
bool acviewer_just_updated = false;

void TrajectoryExplorerApp::on_sb_delta_valueChanged(double val){

    //printf("\n -- on_sb_delta_valueChanged -- value changed %f\n", val);

    if(!ui.cb_sync->isChecked())
        return;



    if(acviewer_just_updated)
        return;

    for(map<uint,ACViewer*>::iterator it = histViewers.begin(); it != histViewers.end(); ++it){
        it->second->updateSelector(val);
    }
    viewer->updateGL();
}
void TrajectoryExplorerApp::slot_selectorChanged(uint atomid, float tscale, float theta){

    //printf(" slot -- %d %f %f -- %d\n", atomid, tscale, theta, ui.cb_sync->isChecked());

    if(!ui.cb_sync->isChecked())
        return;

    //printf("\n -- updating all plots!\n");

    for(map<uint,ACViewer*>::iterator it = histViewers.begin(); it != histViewers.end(); ++it){
        if(it->first != atomid)
            it->second->updateSelector(tscale, theta);
    }

    acviewer_just_updated = true;

    this->ui.sb_delta->setValue(tscale);

    acviewer_just_updated = false;
    this->viewer->updateGL();
}
