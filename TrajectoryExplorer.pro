TARGET = TrajectoryExplorer
TEMPLATE = app

CONFIG *= qt
CONFIG *= release #debug_and_release
CONFIG *= thread warn_off

QT *= xml opengl printsupport

## -------------------------------------------------
## depending upon the platform, set appropriate vals
## -------------------------------------------------
macx {
    OS = macx
    DEFINES *= MAC_OSX
    CONFIG -= app_bundle
}
else {
    OS = unix
    DEFINES *= LINUX
}

MOC_DIR = $$PWD/build/obj/
OBJECTS_DIR = $$PWD/build/obj/
UI_DIR = $$PWD/build/obj/
RCC_DIR = $$PWD/build/obj/
DESTDIR = $$PWD/build/

#QMAKE_CXXFLAGS_DEBUG += -D_GLIBCXX_DEBUG -fpermissive
QMAKE_CXXFLAGS_DEBUG += -fpermissive
QMAKE_CXXFLAGS_RELEASE += -fpermissive

## -------------------------------------------------
## specify external libraries and paths to be used
## -------------------------------------------------

GLEWPATH = /opt/local
#GLEWPATH = $$(HOME)/macports
QGLPATH = $$PWD/extlibs
GLEPATH = $$PWD/extlibs

# QGLViewer (as library)
INCLUDEPATH *= $$QGLPATH/include
LIBS *= $$QGLPATH/lib/libQGLViewer.a

# GLEW
INCLUDEPATH *= $$GLEWPATH/include
LIBS *= $$GLEWPATH/lib/libGLEW.a

# GLE
INCLUDEPATH *= $$GLEPATH/include
LIBS *= $$GLEPATH/lib/libgle.a
DEFINES += USE_GLE

LIBS += -framework Python
#CONFIG += c++11

## -------------------------------------------------

FORMS *= ./TrajectoryViewer.ui ./Statistics.ui

INCLUDEPATH *= ./include
HEADERS = include/Utils.h \
          include/RW_VASP.h \
          include/Atom.h \
          include/bootstrapping.h \
          include/ACData.h \
          include/qcustomplot.h \
          include/histogram.h \
          include/ACViewer.h \
          include/ACPlots.h \
          include/ACWindow.h \
          include/TrajectoryViewer.h \
          include/TrajectoryExplorerApp.h

SOURCES = src/Utils.cpp \
          src/RW_VASP.cpp \
          src/Atom.cpp \
          src/TrajectoryViewer.cpp \
          src/qcustomplot.cpp \
          src/bootstrapping.cpp \
          src/histogram.cpp \
          src/ACViewer.cpp \
          src/ACWindow.cpp \
          src/ACData.cpp \
          src/TrajectoryExplorerApp.cpp \
          src/main.cpp
