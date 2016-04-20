# TrajectoryExplorer

Copyright (c) 2016, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. Written by *Harsh Bhatia
(bhatia4@llnl.gov)*. CODE-682837. All rights reserved.

This file is part of TrajectoryExplorer, Version: 1.0. For details, see
https://github.com/bhatiaharsh/TrajectoryExplorer.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License (as published by the Free Software
Foundation) version 2.1 dated February 1999.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
License for more details.

#### Versions

1.0: released on Apr 20, 2016. This version corresponds to the paper presented at IEEE Pacific Visualization 2016. 

* computation and visualization of relative angle distribution

#### Dependencies

The TrajectoryExplorer has the following dependencies:

1. Qt (tested with Qt 5.5.1 and Qmake 3.0)

2. The OpenGL Extension Wrangler Library (tested with GLEW 1.13)

3. QGLViewer -- A Qt-based 3D viewer: http://libqglviewer.com/ (tested with
QGLViewer 2.6.3)

4. GLE Tubing and Extrusion Library : http://www.linas.org/gle/ (tested with
GLE 3.1.0)

The code assumes that Qt and GLEW are already available, whereas the build
script downloads and installs QGLViewer and GLE automatically in a local path
(see step 2 in INSTALLATION).

To specify path for GLEW, or to correct/modify paths for any of these libraries, please see TrajectoryExplorer/TrajectoryExplorer.pro.

#### Installation & Execution

-- pwd = TrajectoryExplorer

* cd extlibs
* sh buildlibs.sh
	* This downloads and builds QGLviewer and GLE libraries
* cd ..
* mkdir build
* cd build
* qmake ../
* make
	* This builds the source code 
* ./TrajectoryExplorer <config_file>
	* The program requires a config file to specify input parameters. A sample config file can be found in TrajectoryExplorer/configs