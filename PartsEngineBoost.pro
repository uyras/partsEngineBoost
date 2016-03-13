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

HEADERS += partarrayboost.h \
    PartArrayMPI.h \
    statemachinegmp.h \
    wanglandaumpi.h \
    gapmanager.h

SOURCES += PartArrayMPI.cpp \
        statemachinegmp.cpp \
    wanglandaumpi.cpp \
    gapmanager.cpp

OTHER_FILES += \
    README.md \
    .gitignore

LIBS+= -lboost_mpi -lboost_serialization

LIBS += -L$$PWD/../partsEngine/ -lPartsEngine
INCLUDEPATH += $$PWD/../partsEngine
DEPENDPATH += $$PWD/../partsEngine
PRE_TARGETDEPS += $$PWD/../partsEngine/libPartsEngine.a


TEMPLATE = lib
CONFIG += staticlib
DESTDIR = $$PWD
CONFIG(release,debug|release){
    DEFINES += QT_NO_DEBUG_OUTPUT
}
CONFIG += c++11
