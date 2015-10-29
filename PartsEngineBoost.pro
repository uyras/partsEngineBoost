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
    wanglandaumpi.h

SOURCES += PartArrayMPI.cpp \
        statemachinegmp.cpp \
    wanglandaumpi.cpp

OTHER_FILES += \
    README.md \
    .gitignore

LIBS+= -lboost_mpi -lboost_serialization

LIBS += -L$$PWD/../partsEngine/ -lPartsEngine
INCLUDEPATH += $$PWD/../partsEngine
DEPENDPATH += $$PWD/../partsEngine
PRE_TARGETDEPS += $$PWD/../partsEngine/libPartsEngine.a

CONFIG(debug,debug|release) {
    SOURCES += main.cpp
    TEMPLATE = app
    CONFIG += console
}

CONFIG(release,debug|release){
    TEMPLATE = lib
    CONFIG += staticlib
    DESTDIR = $$PWD
    DEFINES += QT_NO_DEBUG_OUTPUT
}
#DEFINES += QT_NO_DEBUG_OUTPUT
