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
#include <QSettings>

/// -----------------------------------------------------------------
using namespace std;

int main(int argc, char **argv) {

    if(argc != 2){
        printf(" Incorrect parameters!\n");
        printf(" Usage: %s <configuration file (*.ini)>\n", argv[0]);
        exit(1);
    }
    if(!QFile(argv[1]).exists()){
        printf(" File (%s) not found!\n", argv[1]);
        exit(1);
    }

    QSettings settings(argv[1], QSettings::IniFormat);
    printf(" Initializing the tool using file %s\n", settings.fileName().toLatin1().data());

    // start the application
    QApplication application(argc,argv);
    TrajectoryExplorerApp atApp(settings);

    // Run main loop.
    return application.exec();
}


#if 0

/// ===========================================================
/// main function
/// ===========================================================
#if 1
#if 0
#include "QDir"
#include "assert.h"
#include "QMessageBox"
#endif

#include <QSettings>

std::string posfile;
std::string track_li_file;
std::string track_ox_file;

void print_usage(const char *com){
    printf(" Incorrect parameters!\n");
    printf(" Usage: %s --pos <positions_file> --li <tracking_Li_IDs_file> --ox <tracking_ox_IDs_file>\n", com);
    exit(1);
}

void parse_arguments(int argc, char **argv){

    if(argc != 7)
        print_usage(argv[0]);

    for(uint i = 1; i < argc; i++){

        std::string tag(argv[i]);

        if(tag.compare("--pos") == 0){            posfile = std::string(argv[++i]); }
        else if(tag.compare("--li") == 0){  track_li_file = std::string(argv[++i]); }
        else if(tag.compare("--ox") == 0){  track_ox_file = std::string(argv[++i]); }
    }

    if(posfile.empty() || track_li_file.empty() || track_li_file.empty())
        print_usage(argv[0]);
}

#include <qcustomplot.h>
#include "Utils.h"







void bootstrap(const vector<int> &data, uint n_resamples) {

    uint n_data = data.size();
    uint num_bins = 100;
    int min_val = 0;
    int max_val = data.size();
#if 0
    vector<histogram_c<int> > bhists (n_resamples);
    vector<estimator<int> > bstats (n_resamples);

    vector<float> means(n_resamples), stdevs(n_resamples);

    for(uint i = 0; i < n_resamples; i++){

        vector<uint> indices = Utils::random<uint>(n_data);
        vector<int> bData = get_data<int>(data, indices);

        bstats[i] = estimator<int> (bData);
        bhists[i] = histogram_c<int> (bData, num_bins, min_val, max_val);

        means[i] = bstats[i].mean;
        stdevs[i] = bstats[i].stdev;
        //printf(" resample %d: ", i);   bstats[i].print();
    }

    // collect data on per bucket basis
    for(uint b = 0; b < num_bins; b++){

        vector<double> hvals (n_resamples);
        for(uint i = 0; i < n_resamples; i++){
            hvals[i] = bhists[i].hist[b];

            //printf(" b = %d, h = %d, val = %f\n", b, i, hvals[i]);
        }

        estimator<double> he(hvals);
        printf(" bucket_estimator[%d] ", b);
        he.print();
    }
    /*
    vector< vector<double> > histvals (num_)

    histogram<int> summary_histogram;
    summary_histogram.hist.resize(num_bins);
    for(uint i = 0; i < num_bins; i++){

    }*/


    histogram_c<int> hist_data (data, num_bins, min_val, max_val);
    hist_data.normalize();
    //hist_data.print("hmean");

    return;
    estimator<int> est_data (data);
    est_data.print("estimator of data ");

    estimator<float> est_mean (means);
    estimator<float> est_stdev (stdevs);
    est_mean.print("estimater of mean: ");
    est_stdev.print("estimater of stdev: ");
#endif
}

#if 0
void test_rw(){

    char filename[] = "testfile.raw";

    // --------------------------------------------------
    // write a file
    {

        uint num_resamples = 1000;

        FILE *datafile = fopen(filename, "wb");
        if(!datafile){
            printf(" Unable to open data file : %s\n", filename);
            exit(1);
        }

        printf(" Writing binary data file %s...", filename);
        fflush(stdout);

        fwrite(&num_resamples, sizeof(uint), 1, datafile);

        for(uint i = 0; i < num_resamples; i++){

            float d = (float)i;
            fwrite(&d, sizeof(float), 1, datafile);
        }

        fclose(datafile);
        printf(" Done!");
    }

    // --------------------------------------------------
    // read it
    {

        uint idx = 3;

        size_t sz = 0;

        // open the file
        FILE *datafile = fopen(filename, "rb");
        if(!datafile){
            std::cerr << " Unable to open data file : "<<filename<<"\n";
            sz = 0;
            return;
        }

        //if(verbose){
            std::cout<<" Reading binary doubles from file "<<filename<<"...";
            fflush(stdout);
        //}

        // read number of samples
        uint num_resamples;
        uint rd_sz = fread(&num_resamples, sizeof(uint), 1, datafile);
        printf(" num resamples = %d\n", num_resamples);

        // i want to read
        //fseek(datafile, idx, SEEK_CUR);

        fseek(datafile, idx * sizeof(float), SEEK_SET);
        //float result;
        //fread(&result, sizeof(float), 1, f);

        float n;
        unsigned int rd_sz = fread(&n, sizeof(float), 1, datafile);

        printf(" read %d %f\n", rd_sz, n);



        /*if(sz == 0){
            // find the size of the file
            fseek(datafile, 0, SEEK_END);
            long fsz = ftell(datafile);
            rewind(datafile);

            printf(" file size")
            if(fsz % sizeof(double) != 0){
                std::cerr << "\n\t - Invalid number of doubles in file "<<filename<<". Size of file = "<<fsz<<", Size of double = "<<sizeof(double)<<" -- mod = "<<(fsz%sizeof(double))<<std::endl;
                fclose(datafile);
                return 0;
            }
            sz = fsz / sizeof(double);
        }


        // read the data
        double *values = new double[sz];
        unsigned int rd_sz = fread(values, sizeof(double), sz, datafile);
        if(rd_sz != sz){
            std::cerr << "\n\t - Expected "<<sz<<", but read "<<rd_sz<<" values!\n";
            sz = rd_sz;
        }*/

        // return
        fclose(datafile);
       // if(verbose)
            printf(" Done! Read %'lu doubles!\n", sz);
        //return values;
    }
}
#endif



int main2(int argc, char** argv){

#if 0
    uint n_data = 100000;
    uint n_resamples = 100;
    vector<int> data = Utils::random<int>(n_data);

    /*int *cdata = new int[n_data];
    QVector<int> qdata(n_data);

    for(uint i = 0; i < n_data; i++){
        cdata[i] = data[i];
        qdata[i] = data[i];
    }

    estimator<int> e (data);
    e.print("estimator of data ");

    estimator<int> ce (cdata, n_data);
    e.print("estimator of cdata ");

    estimator<int> qe (qdata);
    e.print("estimator of qdata ");*/

    bootstrap(data, n_resamples);

    exit(1);
#endif

/*
        uint sz_src = 600;//rawdata->size_tscale();
        uint sz_tgt = 40;//this->size_tscale();

        float frac = (float) sz_src / (float) sz_tgt;

        for(uint tscale_tgt = 0; tscale_tgt < sz_tgt; tscale_tgt++){

            uint tscale_srcb = floorf((float)tscale_tgt * frac);
            uint tscale_srce = floorf((float)(tscale_tgt+1) * frac);
            uint n = tscale_srce-tscale_srcb;

            if(n <= 1){

                printf(" t = %d: s = %d\n", tscale_tgt, tscale_srcb);
                //this->value(theta, tscale_tgt) = rawdata->value(theta, tscale_srcb);
                //this->stddev(theta, tscale_tgt) = rawdata->stddev(theta, tscale_srcb);
            }
            else {

                printf(" t = %d: s = [%d %d] (%d)\n", tscale_tgt, tscale_srcb, tscale_srce, n);

                for(uint i = tscale_srcb; i < tscale_srce; i++){
                    printf("%d %d\n", i, i-tscale_srcb);
                }
                // average
                //this->value(theta, tscale_tgt) = average( rawdata->values_theta(theta), tscale_srcb, tscale_srce );
                //this->stddev(theta, tscale_tgt) = this->value(theta, tscale_tgt) *  (float)(abs(rand())%100)/100.0f;
            }
        }

        exit(1);
*/
/*
    float o[3] = {0,0,0};
    float p[3] = {1,0,0};
    float q[3] = {1,-1,0};

    float v1[3] = { p[0]-o[0], p[1]-o[1], p[2]-o[2] };
    float v2[3] = { q[0]-p[0], q[1]-p[1], q[2]-p[2] };

    float v3[3] = { q[0]-o[0], q[1]-o[1], q[2]-o[2] };

    float v = Utils::dot(v1, v2) / (Utils::magn(v1)*Utils::magn(v2));
    float *c = Utils::cross(v1, v2);
    printf(" magn of cross product = %f (%f,%f,%f)\n", Utils::magn(c), c[0], c[1], c[2]);

    float *d = Utils::cross(c, v3);
    printf(" magn of cross product 2 = %f (%f,%f,%f)\n", Utils::magn(d), d[0], d[1], d[2]);

    float result = acos (v) * 180.0 / M_PI;
    printf ("The arc cosine of %f is %f degrees.\n", v, result);


      exit(1);
  */
    // -----------
#if 0
    QDir dir = QDir::current();
    //printf(" dir = -%s-\n", dir.absolutePath().toLatin1().data());

    dir.cdUp();
    //printf(" dir = -%s-\n", dir.absolutePath().toLatin1().data());

    dir.cd("PlugIns");
    //printf(" dir = -%s-\n", dir.absolutePath().toLatin1().data());

    QCoreApplication::setLibraryPaths(QStringList(dir.absolutePath()));

    //printf("after change, libraryPaths=(%s)\n", QCoreApplication::libraryPaths().join(",").toUtf8().data());
#endif

    parse_arguments(argc, argv);

    // Read command lines arguments.
    QApplication application(argc,argv);

    AtomTrajectoryApp atApp(posfile, track_li_file, track_ox_file);

    // Run main loop.
    return application.exec();

}




#else

// call_function.c - A sample of calling
// python functions from C code
//
#include <Python/Python.h>

int main(int argc, char *argv[])
{
    PyObject *pName, *pModule, *pDict, *pFunc, *pValue;

    /*if (argc < 3)
    {
        printf("Usage: exe_name python_source function_name\n");
        return 1;
    }*/

    printf("1\n");

    // Initialize the Python Interpreter
    Py_Initialize();

    printf("2\n");

    // Build the name object
    //pName = PyString_FromString(argv[1]);
    pName = PyString_FromString("py_function");

    printf("3\n");

    // Load the module object
    pModule = PyImport_Import(pName);


    printf("4\n");
    // pDict is a borrowed reference
    pDict = PyModule_GetDict(pModule);

    const char* s = PyString_AsString(pDict);
    cout<<s<<endl;



    printf("5\n");
    // pFunc is also a borrowed reference
    //pFunc = PyDict_GetItemString(pDict, argv[2]);
    pFunc = PyDict_GetItemString(pDict, "multiply");

    printf("6\n");
    if (PyCallable_Check(pFunc))
    {
        PyObject_CallObject(pFunc, NULL);
    } else
    {
        PyErr_Print();
    }

    printf("7\n");
    // Clean up
    Py_DECREF(pModule);
    Py_DECREF(pName);

    // Finish the Python Interpreter
    Py_Finalize();

    return 0;
}
#endif



#if 0
#include "LithiumApp.h"
#include "LithiumViewer.h"

#include <qapplication.h>

#include "carRW.h"


/// convert Force/Force.# files into binary
void convert_forceFiles(){

    std::string mdir("/Users/bhatia4/WORK/data/BES-data");

    uint num_grp = 120;
    uint num_mats = 638;
    uint num_ts = 625;

    for(uint grp = 1; grp <= num_grp; grp++){

        // read!
        char fname[mdir.length()+20];
        sprintf(fname, "%s/Force/Force.%d", mdir.c_str(), grp);

        std::vector<std::vector<float*> > fpos;     // position read from force file (in Angstroms)
        std::vector<std::vector<float*> > fforce;    // force in Angstroms
        std::vector<std::vector<float> > lastterm;

        RW::read_forceFile(fname, fpos, fforce, lastterm);

        // now write!
        sprintf(fname, "%s.raw", fname);
        RW::write_forceFile_bin(fname, fpos, fforce, lastterm);
/*


        std::vector<std::vector<float*> > fpos2;     // position read from force file (in Angstroms)
        std::vector<std::vector<float*> > fforce2;    // force in Angstroms
        std::vector<std::vector<float> > lastterm2;

        RW::read_forceFile_bin(fname, fpos2, fforce2, lastterm2);

        for(uint ts = 0; ts < num_ts; ts++){
        for(uint mat = 0; mat < num_mats; mat++){
        for(uint k = 0; k < 3; k++){

            if(fabs(fpos[ts][mat][k] - fpos2[ts][mat][k]) > powf(10,-5)){
                printf(" mismatch in values!\n");
            }
        }
        }
        }*/

        //break;
    }
}

#if 0
/*
std::string mdir("/Users/bhatia4/WORK/data/BES-data");

uint num_mats = 638;
uint num_ts = 625;
uint num_grp = 125;

//float fact = 19.283;

for(uint grp = 1; grp <= num_grp; grp++){

    char ffname[mdir.length()+20];
    //char xfname[mdir.length()+20];

    sprintf(ffname, "%s/Force/Force.%d", mdir.c_str(), grp);
    //sprintf(xfname, "%s/xdat/XDATCAR.%d", mdir.c_str(), grp);

    //std::vector<std::vector<float*> > xpos;     // position read from xdat car file (in direct coordinates)
    std::vector<std::vector<float*> > fpos;     // position read from force file (in Angstroms)
    std::vector<std::vector<float*> > fforce;    // force in Angstroms
    std::vector<std::vector<float> > lastterm;

    //RW::read_xdatFile(xfname, xpos);
    RW::read_forceFile(ffname, fpos, fforce, lastterm);

    for(uint ts = 0; ts < num_ts; ts++){
    for(uint mat = 0; mat < num_mats; mat++){

        const float *txpos = xpos[ts][mat];
        const float *tfpos = fpos[ts][mat];

        for(uint k = 0; k < 3; k++){

            if( fabs(fact*txpos[k] - tfpos[k]) > powf(10,-5) )
                printf(" mismatch in position %.6f != %.6f! grp = %d, ts = %d, mat = %d, k = %d\n",
                       fact*txpos[k], tfpos[k], grp, ts, mat, k);
        }
    }
    }
}



/// i will write data in binary
unsigned long int sz = num_ts*num_grp;

float *pos = new float[3*sz];
float *force = new float[3*sz];
uint ind = 0;

char xffname[mdir.length()+20];
sprintf(ffname, "%s/XDATForce.%d", mdir.c_str(), grp);

// write for each material
for(uint mat = 0; mat < num_mats; mat++){

    for(uint grp = 1; grp <= num_grp; grp++){
    for(uint ts = 0; ts < num_ts; ts++){

        pos[3*ind] =



        for(uint k = 0; k < 3; k++){

            if( fabs(fact*txpos[k] - tfpos[k]) > powf(10,-5) )
                printf(" mismatch in position %.6f != %.6f! grp = %d, ts = %d, mat = %d, k = %d\n",
                       fact*txpos[k], tfpos[k], grp, ts, mat, k);
        }
        }
        }
    }
    ind++:
}
exit(1);
*/
/*
char *mainfilename = "/Users/bhatia4/WORK/data/BES-data/Force/Force.38";

std::vector<std::vector<float*> > pos;
std::vector<std::vector<float*> > force;
std::vector<std::vector<float> > lastterm;

RW::read_forceFile(mainfilename, pos, force, lastterm);

for(int i = 0; i <= 625; i+=25){

    uint ind = (i == 0) ? 0 : i-1;

    char newfilename[100];
    if(i < 10)
        sprintf(newfilename, "%s.00%d", mainfilename, (ind+1));
    else if(i < 100)
        sprintf(newfilename, "%s.0%d", mainfilename, (ind+1));
    else
        sprintf(newfilename, "%s.%d", mainfilename, (ind+1));

    //char *tempfilename = "/Users/bhatia4/WORK/data/BES-data/AECCAR-38/Force.38.001";

    //printf(" writing index %d to file %s\n", ind, newfilename );
    RW::write_forceFile(newfilename, pos[ind], force[ind], lastterm[ind]);

    std::vector<float*> tpos;
    std::vector<float*> tforce;
    std::vector<float> tlastterm;

    RW::read_forceFile(newfilename, tpos, tforce, tlastterm);

    //exit(1);

    //break;
}
exit(1);
//*/
/*for(float ts = 38; ts < 38.626; ts += 0.025){
    printf(" ts = %.3f, pts = %.4f\n", ts, RW::get_physicaltime(ts));
}
for(float ts = 39; ts < 39.626; ts += 0.025){
    printf(" ts = %.3f, pts = %.4f\n", ts, RW::get_physicaltime(ts));
}
for(float ts = 48; ts < 48.626; ts += 0.025){
    printf(" ts = %f, pts = %f\n", ts, RW::get_physicaltime(ts));
}
/*for(float ts = 49; ts < 49.626; ts += 0.025){
    printf(" ts = %.3f, pts = %.4f\n", ts, RW::get_physicaltime(ts));
}*/
//exit(1);
#endif

bool USE_LOG = false;

/// parse the command, to read the parameters
bool parse_command(int argc, char *argv[], std::vector<std::string> &args, bool display_usage = true){

    if(argc < 2 || argc > 3){
        std::cerr << " Usage: "<<argv[0]<<" [input-filename] [-use_log]\n";
        return false;
    }

    for(uint i = 1; i < argc; i++){

        std::string arg(argv[i]);

        //cout<<"\n arg["<<i<<"] = "<<arg<<endl;

        // -----
        if(arg.compare("-use_log") == 0){
            USE_LOG = true;
            continue;
        }

        // -----
        if(arg.at(0) == '-'){
            std::cerr <<" Incorrect argument: " << arg << ". The accepted command line arguments are - \n";
            std::cerr << "\t-use_log           : use log of the function\n";
            return false;
        }

        // -----
        else args.push_back(arg);

        //printf(" Ignoring unidentified argument %s\n", arg.c_str());
    }
return true;
}

#include "utils.h"

/// -----------------------
///  generating poscar for new start

using namespace std;

// ts : this is where I want to change the simulation
// percent :  how much of the velocity to add ?

void newsimulation(uint ts, float percent){

    printf(" start at %f\n", RW::getTime_fn_to_physical(39, ts));

    vector<vector<float* > > pos;
    RW::read_xdatFile("/Users/bhatia4/WORK/data/BES-data/xdat/XDATCAR.39", pos);

    // the ids of relevent atoms
    uint LiID = 637;
    //uint redOxID = 551;
    uint redOxID = 515;     // blue!

    float nudge_factor = 0.005;

    // compute the velocities
    vector<float*> m_vels(pos[ts].size());

    for(uint i = 0; i < m_vels.size(); i++){

        // these positions are in direct coordinates
        // velocity is in Ang/fs. So, velocity should be multiplied by lattice constant/timestep (Mitchell;s email)
        const float* m_currpos = pos[ts][i];
        const float* m_prevpos = pos[ts-1][i];

        // compute the velocity
        m_vels[i] = new float[3];

        for(uint k = 0; k < 3; k++){
            m_vels[i][k] = Utils::get_torusDisplacement(m_currpos[k], m_prevpos[k], 1);
            m_vels[i][k] *= (19.283/0.5);   // proper scaling
        }

        //printf(" cp = [%f,%f,%f]\n, op = [%f,%f,%f]\n", m_currpos[0],m_currpos[1],m_currpos[2],m_prevpos[0],m_prevpos[1],m_prevpos[2]);
        //printf(" vel = [%f,%f,%f]\n", m_vels[i][0],m_vels[i][1],m_vels[i][2]);

        //exit(1);
        // change the velocity for Lithium only
        if(i != LiID)
            continue;

        // compute the vectors to red and blu oxygens
        const float* red_currpos = pos[ts][redOxID];
        float li_red_vel[3];

        // get the displacement vector Li-Ox!
        for(uint k = 0; k < 3; k++){
            li_red_vel[k] = Utils::get_torusDisplacement(red_currpos[k], m_currpos[k], 1);
            li_red_vel[k] *= (19.283/0.5);  // proper scaling // why 0.5 ?
        }

        // normalize the displacement vector
        float mredvel = Utils::magn(li_red_vel);
        for(uint k = 0; k < 3; k++){
            li_red_vel[k] /= mredvel;
        }

        printf("\n li_red_vel = (%f,%f,%f) %f\n", li_red_vel[0], li_red_vel[1], li_red_vel[2], Utils::magn(li_red_vel));
        printf(" li_vel = (%f,%f,%f) %f\n", m_vels[i][0], m_vels[i][1], m_vels[i][2], Utils::magn(m_vels[i]));

        float changevel[3];
        for(uint k = 0; k < 3; k++){
            changevel[k] = nudge_factor*li_red_vel[k];
            m_vels[i][k] += changevel[k];
        }

/*
        //exit(1);
        // the li's displacement is orders of magnitude smaller
        // to change the velocity, i should scale down the li-red vector

        float orderdiff = Utils::magn(m_vels[i]) / Utils::magn(li_red_vel);
        printf("\n order diff = %f\n", orderdiff);

        float changevel[3];
        for(uint k = 0; k < 3; k++){
            changevel[k] = percent*orderdiff*li_red_vel[k];
            m_vels[i][k] += percent*orderdiff*li_red_vel[k];
        }
*/
        float changevelmagn = Utils::magn(changevel);

        printf("\n change_vel in (Ang/fs) = (%f,%f,%f) -- %f\n", changevel[0], changevel[1], changevel[2], changevelmagn);
        printf("\n new_li_vel = (%f,%f,%f)\n", m_vels[i][0], m_vels[i][1], m_vels[i][2]);
    }

    return;

    char outfilenm[100];
    sprintf(outfilenm, "/Users/bhatia4/WORK/data/BES-data/new-MD_input/toBluePOSCAR.39.%d-%.1f", ts, percent);
    RW::write_poscar(pos[ts], m_vels, outfilenm);

/*
    // --------------------------------------------
    // grid coordinates!

    printf(" \n grid coordinates!\n");
    float ligcurrpos[3];
    float ligprevpos[3];
    float redcurrpos[3];
    float blucurrpos[3];

    uint dims[3] = {280,280,280};
    for(uint k = 0; k < 3; k++){

        ligcurrpos[k] = RW::getCoord_direct_to_grid(li_currpos[k], 280);
        ligprevpos[k] = RW::getCoord_direct_to_grid(li_prevpos[k], 280);
        redcurrpos[k] = RW::getCoord_direct_to_grid(red_currpos[k], 280);
        blucurrpos[k] = RW::getCoord_direct_to_grid(blu_currpos[k], 280);
    }


    printf(" redcurrpos = (%f,%f,%f)\n", red_currpos[0], red_currpos[1], red_currpos[2]);
    printf(" bluprevpos = (%f,%f,%f)\n", blu_currpos[0], blu_currpos[1], blu_currpos[2]);

    printf(" distance from red = %f\n", sqrt(Utils::get_sq_torusdist(ligcurrpos, redcurrpos, dims)));
    printf(" distance from blu = %f\n", sqrt(Utils::get_sq_torusdist(ligcurrpos, blucurrpos, dims)));
*/
}

/*
 * ligpos[j] = RW::getCoord_direct_to_grid(pos[i][637][j], dims[j]);
                oxgpos[j] = RW::getCoord_direct_to_grid(pos[i][matid][j], dims[j]);
                */

/// -----------------------

int main(int argc, char** argv){

    /*
    vector<vector<float*> > pos;
    std::string filename("/Users/bhatia4/WORK/data/BES-data/XDATCAR_all.bin");
    RW::read_xdatAllFile_binary(filename.c_str(), pos);


    for(uint i = 0; i < pos.size(); i++){

        float *LiPos = pos[i][637];
        printf(" LiPos[%d] = (%f,%f,%f)\n", i, LiPos[0], LiPos[1], LiPos[2]);

        break;
    }

    exit(1);
*/
    //newsimulation(400, 0.5);
    //exit(1);

    std::vector<std::string> args;
    if(!parse_command(argc, argv, args))
        exit(1);

    // Read command lines arguments.
    QApplication application(argc,argv);

    LithiumApp lithApp(argv[1], true, USE_LOG);

    // Run main loop.
    return application.exec();
}
#endif
#endif
