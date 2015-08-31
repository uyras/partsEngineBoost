#-------------------------------------------------
#
# Project created by QtCreator 2014-01-31T10:01:07
#
#-------------------------------------------------

TARGET = PartsEngineBoost
QT       += core
QT       -= gui

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

HEADERS += partarrayboost.h \
    PartArrayMPI.h \
    statemachinegmp.h

SOURCES += PartArrayMPI.cpp \
        statemachinegmp.cpp

OTHER_FILES += \
    README.md \
    .gitignore

LIBS+= -lboost_mpi -lboost_serialization

CONFIG(debug,debug|release) {
    SOURCES += main.cpp
    TEMPLATE = app
    CONFIG += console
}

CONFIG(release,debug|release){
    TEMPLATE = lib
    CONFIG += staticlib
    DESTDIR = $$PWD
}
