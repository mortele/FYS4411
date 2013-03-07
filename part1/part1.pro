TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    functions.cpp \
    variationalmc.cpp\
    ../newfunctions.cpp \
    varmc.cpp \
    variationalloop.cpp \
    wavefunctions.cpp

HEADERS += \
    functions.h \
    variationalmc.h \
    varmc.h \
    variationalloop.h \
    wavefunctions.h

release {
    DEFINES += ARMA_NO_DEBUG
    QMAKE_LFLAGS -= -O1
    QMAKE_LFLAGS += -O3
    QMAKE_LFLAGS_RELEASE -= -O1
    QMAKE_LFLAGS_RELEASE += -O3
    QMAKE_CXXFLAGS -= -O2
    QMAKE_CXXFLAGS += -O3
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}
LIBS+= -larmadillo -lblas -llapack
