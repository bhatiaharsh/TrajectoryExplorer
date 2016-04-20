/*
Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by Harsh Bhatia
(bhatia4@llnl.gov). CODE-682837. All rights reserved.

This file is part of TrajViewer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

For more details on the Licence, please read LICENCE file.
*/

#ifndef _TRAJECTORYEXPLORER_APP_H
#define _TRAJECTORYEXPLORER_APP_H

#include <set>
#include <vector>

#include <QMainWindow>
#include <QKeyEvent>

#include "ui_TrajectoryViewer.h"

#include "Atom.h"
#include "ACViewer.h"
#include "TrajectoryViewer.h"
#include "Utils.h"
#include "RW_VASP.h"

/// ===========================================================
/// TrajectoryExplorerApp
/// ===========================================================

class TrajectoryExplorerApp : public QMainWindow {
    Q_OBJECT

private:
    TrajectoryViewer *viewer;          // pointer to the viewer
    std::map<uint,ACViewer*> histViewers;

    void update_LiOxLabels();

    bool center_at(int mid, int time);
    void find_jumps();
    void create_tubes();

    void compute_distance_bonds(float dist_thresh);


    QSettings *appsettings;
public:

    Ui::TrajectoryViewerUI ui;
    //QActionGroup* tfunction;

    VASPmetadata mdata;

    std::vector<Atom> atoms;

    std::vector<unsigned int> track_li;
    std::vector<unsigned int> track_ox;

    // current center of the viewing box: atom idx, and time idx
    std::pair<int,int> current_center;

    // ----------------------------------------------------
    // constructor and destructor

    TrajectoryExplorerApp(QSettings &settings);
    ~TrajectoryExplorerApp(){}

    // ----------------------------------------------------

    int find_atom(uint mid) const {
        for(uint i = 0; i < atoms.size(); i++){
            if(atoms[i].mid == mid)
                return i;
        }
        return -1;
    }

    static inline int get_atomID_from_UI(const QListWidgetItem* selectedItem) {
        QString tag = selectedItem->text();
        uint spc_indx = tag.indexOf(' ');
        return tag.mid(0,spc_indx).toInt();     // convert from VASP to C++ id
    }

public slots:
    void slot_selectorChanged(uint atomid, float tscale, float theta);

protected slots:

    void on_groupBox_static_toggled(bool state){

        ui.groupBox_anim->setChecked(!state);
        if(state && viewer->animationIsStarted())
                viewer->toggleAnimation();
        viewer->updateGL();
    }

    void on_groupBox_anim_toggled(bool state){
        ui.groupBox_static->setChecked(!state);
        viewer->updateGL();
    }

    void on_pb_animPause_clicked(){
        if(viewer->animationIsStarted())
            viewer->toggleAnimation();
    }

    void on_pb_animResume_clicked(){
        if(!viewer->animationIsStarted())
            viewer->toggleAnimation();
    }

    void on_pb_animSetTime_clicked(){
        viewer->live_animation_time = ui.sb_time->text().toInt();
        viewer->updateGL();
    }

    void on_cb_gle_tubes_toggled(bool state){
        ui.dsb_tubeThickness->setEnabled(state);
        viewer->updateGL();
    }

    void on_cb_gle_bonds_toggled(bool){   viewer->updateGL(); }
    void on_cb_bbox_toggled(bool){        viewer->updateGL(); }
    void on_cb_axes_toggled(bool){        viewer->updateGL(); }
    void on_cb_tets_toggled(bool){        viewer->updateGL(); }
    void on_gb_bonds_toggled(bool){       viewer->updateGL(); }

    void on_pb_allOx_clicked(){
        this->ui.lw_atoms_ox->clearSelection();
        viewer->updateGL();
    }

    void on_pb_center_clicked(){
        if(center_at(ui.cb_center->currentText().toInt(), 0)){
            find_jumps();
            create_tubes();
        }
        viewer->updateGL();
    }

    void on_sb_delta_valueChanged(double);
    void on_sb_offset_valueChanged(int){        viewer->updateGL(); }

    void on_lw_atoms_li_itemSelectionChanged (){        update_LiOxLabels();   viewer->updateGL(); }
    void on_lw_atoms_ox_itemSelectionChanged (){        viewer->updateGL(); }
    void on_sb_animSpeed_valueChanged(int){             viewer->updateGL(); }
    void on_sb_tailLength_valueChanged(int){            viewer->updateGL(); }
    void on_dsb_tubeThickness_valueChanged(double){     viewer->updateGL(); }

    //void on_dsb_bondDist_valueChanged(double){
    void on_pb_computeBonds_clicked(){
        compute_distance_bonds(ui.dsb_bondDist->value());
        update_LiOxLabels();
        viewer->updateGL();
    }

    void on_pb_hists_clicked();
};
#endif
