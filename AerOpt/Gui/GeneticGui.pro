#-------------------------------------------------
#
# Project created by QtCreator 2015-04-01T19:17:07
#
#-------------------------------------------------

QT  += \
    core \
    widgets \
    printsupport

CONFIG(debug, debug|release) {
    DESTDIR = build/debug-gcc
    OBJECTS_DIR = build/debug-gcc/obj
    MOC_DIR = build/debug-gcc/moc
    RCC_DIR = build/debug-gcc/rcc
    UI_DIR = build/debug-gcc/ui
    win32:QMAKE_CXXFLAGS += -static
    win32:QMAKE_CXXFLAGS += -static-libgcc
    win32:QMAKE_LFLAGS += -static
    win32:QMAKE_LFLAGS += -static-libgcc
} else {
    DESTDIR = ../
    OBJECTS_DIR = build/release-gcc/obj
    MOC_DIR = build/release-gcc/moc
    RCC_DIR = build/release-gcc/rcc
    UI_DIR = build/release-gcc/ui
    win32:QMAKE_CXXFLAGS += -static
    win32:QMAKE_CXXFLAGS += -static-libgcc
    win32:QMAKE_LFLAGS += -static
    win32:QMAKE_LFLAGS += -static-libgcc
}

#QMAKE_CC = /usr/local/Cellar/gcc@6/6.4.0/bin/gcc-6
#QMAKE_CXX = /usr/local/Cellar/gcc@6/6.4.0/bin/g++-6
#QMAKE_LINK = /usr/local/Cellar/gcc@6/6.4.0/bin/g++-6
QMAKE_MAC_SDK = macosx10.12

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3
QMAKE_CXXFLAGS += -std=c++11
QMAKE_CXXFLAGS += -Wno-unused-parameter
QMAKE_CXXFLAGS += -Wno-deprecated
QMAKE_CXXFLAGS += -Wno-unused-function
QMAKE_CXXFLAGS += -fmessage-length=0

TARGET = AerOptGui
TEMPLATE = app

SOURCES += \
    main.cpp \
    Gui/MainWindow.cpp \
    Gui/DebugOutput.cpp \
    Gui/TreeView.cpp \
    Gui/Menus/Menu.cpp \
    Gui/Canvas.cpp \
    Gui/GuiComponents/ProjectData.cpp \
    Gui/Dialogs/ConstraintsDialog.cpp \
    Gui/Dialogs/BoundaryDialog.cpp \
    Gui/Dialogs/OptimiserDialog.cpp \
    Gui/Dialogs/MeshDialog.cpp \
    Gui/Plotter/qcustomplot.cpp \
    Gui/Plotter/PlotterDialog.cpp \
    Gui/GuiComponents/Arrow.cpp \
    Gui/GuiComponents/BoundaryPoint.cpp

HEADERS  += \
    Gui/MainWindow.h \
    Gui/DebugOutput.h \
    Gui/TreeView.h \
    Gui/Menus/Menu.h \
    Gui/Canvas.h \
    Gui/GuiComponents/ProjectData.h \
    Gui/Dialogs/ConstraintsDialog.h \
    Gui/Enumerations.h \
    Gui/Dialogs/BoundaryDialog.h \
    Gui/Dialogs/OptimiserDialog.h \
    Gui/Dialogs/MeshDialog.h \
    Gui/Plotter/qcustomplot.h \
    Gui/Plotter/PlotterDialog.h \
    Gui/GuiComponents/Arrow.h \
    Gui/GuiComponents/BoundaryPoint.h

FORMS    += \
    Gui/MainWindow.ui \
    Gui/DebugOutput.ui \
    Gui/TreeView.ui \
    Gui/Dialogs/ConstraintsDialog.ui \
    Gui/Dialogs/ObjectiveDialog.ui \
    Gui/Dialogs/BoundaryDialog.ui \
    Gui/Dialogs/OptimiserDialog.ui \
    Gui/Dialogs/MeshDialog.ui \
    Gui/Plotter/PlotterDialog.ui

INCLUDEPATH += \
    Gui \
    Gui/Dialogs \
    Gui/Menus \
    Gui/GuiComponents \
    Gui/Plotter

RESOURCES += \
    Resourses.qrc
