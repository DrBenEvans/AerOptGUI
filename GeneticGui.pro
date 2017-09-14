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
    #Gui/Dialogs/ConstraintsDialog.cpp \
    Gui/Dialogs/MeshDialog.cpp \
    Gui/Plotter/qcustomplot.cpp \
    Gui/GuiComponents/Arrow.cpp \
    Gui/Dialogs/ConfigSimulationDialog.cpp \
    #Gui/AppController.cpp \
    Gui/Dialogs/OptimisationManagerDialog.cpp \
    Gui/MainWindow.cpp \
    Core/Mesh.cpp \
    Core/Profile.cpp \
    #Gui/Dialogs/PlotterDialog.cpp \
    Core/Optimisation.cpp \
    #Core/Simulation.cpp \
    Gui/GuiComponents/ZoomPanView.cpp \
    Gui/GraphicsItems/BoundaryPointView.cpp \
    Gui/GraphicsItems/MeshView.cpp \
    Gui/GraphicsItems/ProfileView.cpp \
    Gui/GraphicsItems/ControlPointBoundingBox.cpp \
    Gui/GraphicsItems/ControlPointDragHandle.cpp \
    Core/BoundaryPoint.cpp \
    Gui/Models/ProfileModel.cpp \
    Gui/Models/OptimisationModel.cpp \
    Gui/GuiComponents/ControlPointView.cpp \
    Core/BoundaryPointModel.cpp \
    Core/MeshDialogModel.cpp \
    Gui/GraphicsItems/ViewScaler.cpp

HEADERS  += \
    Gui/DebugOutput.h \
    #Gui/Dialogs/ConstraintsDialog.h \
    Gui/Dialogs/MeshDialog.h \
    Gui/Plotter/qcustomplot.h \
    Gui/GuiComponents/Arrow.h \
    Gui/Dialogs/ConfigSimulationDialog.h \
    #Gui/AppController.h \
    Gui/Dialogs/OptimisationManagerDialog.h \
    Gui/MainWindow.h \
    Core/Mesh.h \
    Core/Profile.h \
    Core/Enumerations.h \
    Core/CustomTypes.h \
    #Gui/Dialogs/PlotterDialog.h \
    Core/Optimisation.h \
    #Core/Simulation.h \
    Gui/GuiComponents/ZoomPanView.h \
    Gui/GraphicsItems/MeshView.h \
    Gui/GraphicsItems/ProfileView.h \
    Gui/GraphicsItems/BoundaryPointView.h \
    Gui/GraphicsItems/ControlPointBoundingBox.h \
    Gui/GraphicsItems/ControlPointDragHandle.h \
    Core/BoundaryPoint.h \
    Gui/Models/ProfileModel.h \
    Gui/Models/OptimisationModel.h \
    Gui/GuiComponents/ControlPointView.h \
    Core/BoundaryPointModel.h \
    Core/MeshDialogModel.h \
    Gui/GraphicsItems/ViewScaler.h

FORMS    += \
    Gui/DebugOutput.ui \
#    Gui/Dialogs/ConstraintsDialog.ui \
#    Gui/Plotter/PlotterDialog.ui \
    Gui/Dialogs/ConfigSimulationDialog.ui \
    Gui/Dialogs/OptimisationManagerDialog.ui \
    Gui/MainWindow.ui \
    Gui/Dialogs/MeshDialog.ui \
    Gui/GuiComponents/ControlPointView.ui

INCLUDEPATH += \
    Gui \
    Gui/Dialogs \
    Gui/GuiComponents \
    Gui/GraphicsItems \
    Gui/Plotter \
    Gui/Models \
    Core

RESOURCES += \
    Resourses.qrc
