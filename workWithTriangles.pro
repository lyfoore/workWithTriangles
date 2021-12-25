TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    main.cpp \
    display/GUI.cpp \
    model.cpp \
    distance.cpp

HEADERS += \
    my_include/gl.h \
    my_include/glext.h \
    my_include/glu.h \
    my_include/glut.h \
    my_include/png.h \
    my_include/pngconf.h \
    display/GUI.h \
    model.h \
    distance.h


#QMAKE_CXXFLAGS += -O2
QMAKE_CXXFLAGS_RELEASE += -O3 -ffast-math  -msse -std=c++11


QMAKE_LFLAGS += -O3 -ffast-math  -msse -std=c++11
#QMAKE_LFLAGS += -O3  -msse -std=c++11
unix{
LIBS+=  -lGL -lGLU -lglut -lm
}
win32{
LIBS += -lopenGL32 -lGLU32 -lm
LIBS += -L$$PWD/my_lib -lglut32
}

