#! /bin/bash

IN_DIR=$PWD

# --------------------------------------
# install QGLViewer

if [ ! -d "libQGLViewer-2.6.3" ]; then
    wget http://www.libqglviewer.com/src/libQGLViewer-2.6.3.tar.gz
    tar -xzf libQGLViewer-2.6.3.tar.gz
fi

cd libQGLViewer-2.6.3/QGLViewer
qmake -spec macx-g++ PREFIX=$IN_DIR QGLVIEWER_STATIC=YES
make -j8
make install

cd $IN_DIR
rm -r libQGLViewer-2.6.3

# --------------------------------------
# install GLE

if [ ! -d "gle-3.1.0" ]; then
    wget http://www.linas.org/gle/pub/gle-3.1.0.tar.gz
    tar -xvf gle-3.1.0.tar.gz
fi

cd gle-3.1.0
./configure --prefix=$IN_DIR --x-includes=/usr/X11/include --x-libraries=/usr/X11/lib

ln -s /usr/include/malloc/malloc.h ./malloc.h

make
make install

cd $IN_DIR
rm -r gle-3.1.0
