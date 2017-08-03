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
    Gui/DebugOutput.cpp \
    Gui/Menus/Menu.cpp \
    Gui/Canvas.cpp \
    Gui/Dialogs/ConstraintsDialog.cpp \
    Gui/Dialogs/MeshDialog.cpp \
    Gui/Plotter/qcustomplot.cpp \
    Gui/Plotter/PlotterDialog.cpp \
    Gui/GuiComponents/Arrow.cpp \
    Gui/GuiComponents/BoundaryPoint.cpp \
    Gui/Dialogs/ConfigSimulationDialog.cpp \
    Gui/AppController.cpp \
    Gui/Dialogs/OptimisationManagerDialog.cpp \
    Gui/MainWindow.cpp \
    Gui/Models/OptimisationRun.cpp \
    Gui/Models/Profile.cpp \
    Gui/Models/Mesh.cpp \
    Core/ProfileModel.cpp \

HEADERS  += \
    Gui/DebugOutput.h \
    Gui/Menus/Menu.h \
    Gui/Canvas.h \
    Gui/Dialogs/ConstraintsDialog.h \
    Gui/Enumerations.h \
    Gui/Dialogs/MeshDialog.h \
    Gui/Plotter/qcustomplot.h \
    Gui/Plotter/PlotterDialog.h \
    Gui/GuiComponents/Arrow.h \
    Gui/GuiComponents/BoundaryPoint.h \
    Gui/Dialogs/ConfigSimulationDialog.h \
    Gui/AppController.h \
    Gui/Dialogs/OptimisationManagerDialog.h \
    Gui/MainWindow.h \
    Gui/Models/OptimisationRun.h \
    Gui/Models/Profile.h \
    Gui/Models/Mesh.h \
    Gui/CustomTypes.h \
    Core/ProfileModel.h \

FORMS    += \
    Gui/DebugOutput.ui \
    Gui/Dialogs/ConstraintsDialog.ui \
    Gui/Plotter/PlotterDialog.ui \
    Gui/Dialogs/ConfigSimulationDialog.ui \
    Gui/Dialogs/OptimisationManagerDialog.ui \
    Gui/MainWindow.ui \
    Gui/Dialogs/MeshDialog.ui

INCLUDEPATH += \
    Gui \
    Gui/Dialogs \
    Gui/Menus \
    Gui/GuiComponents \
    Gui/Models \
    Gui/Plotter

RESOURCES += \
    Resourses.qrc
