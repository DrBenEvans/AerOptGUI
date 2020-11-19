#-------------------------------------------------
#
# Project created by QtCreator 2015-04-01T19:17:07
#
#-------------------------------------------------

QT  += \
    core \
    widgets \
    charts \
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
    DESTDIR = release
    OBJECTS_DIR = build/release-gcc/obj
    MOC_DIR = build/release-gcc/moc
    RCC_DIR = build/release-gcc/rcc
    UI_DIR = build/release-gcc/ui
    win32:QMAKE_CXXFLAGS += -static
    win32:QMAKE_CXXFLAGS += -static-libgcc
    win32:QMAKE_LFLAGS += -static
    win32:QMAKE_LFLAGS += -static-libgcc
}

CONFIG+=static

QMAKE_MAC_SDK = macosx10.12

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3
QMAKE_CXXFLAGS += -std=c++14
QMAKE_CXXFLAGS += -Wno-unused-parameter
QMAKE_CXXFLAGS += -Wno-deprecated
QMAKE_CXXFLAGS += -Wno-unused-function
QMAKE_CXXFLAGS += -fmessage-length=0


TARGET = AerOptGui
TEMPLATE = app

SOURCES += \
	Gui/GuiComponents/ChartView.cpp \
    Gui/GuiComponents/ParallelCoordinatesWindow.cpp \
    Gui/GuiComponents/chart.cpp \
    main.cpp \
    Gui/DebugOutput.cpp \
    Gui/Dialogs/MeshDialog.cpp \
    Gui/Plotter/qcustomplot.cpp \
    Gui/Dialogs/ConfigSimulationDialog.cpp \
    Gui/MainWindow.cpp \
    Core/Mesh.cpp \
    Core/Profile.cpp \
    Core/Optimisation.cpp \
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
    Gui/GraphicsItems/ViewScaler.cpp \
    Core/FileManipulation.cpp \
    Core/ProcessManager.cpp \
    Gui/GuiComponents/PlotterWidget.cpp \
    Core/MeshDialogModel.cpp \
    Core/MeshModel.cpp \
    Gui/GraphicsItems/ColorMapper.cpp \
    Gui/Dialogs/windowsizing.cpp \
    Gui/Dialogs/ClusterLoginDialog.cpp \
    Core/clusterManager.cpp\
    Gui/Dialogs/ViewConfigureDialog.cpp
	
HEADERS  += \
    Gui/DebugOutput.h \
    Gui/Dialogs/MeshDialog.h \
    Gui/GuiComponents/ChartView.h \
    Gui/GuiComponents/ParallelCoordinatesWindow.h \
    Gui/GuiComponents/chart.h \
    Gui/Plotter/qcustomplot.h \
    Gui/Dialogs/ConfigSimulationDialog.h \
    Gui/MainWindow.h \
    Core/Mesh.h \
    Core/Profile.h \
    Core/Enumerations.h \
    Core/CustomTypes.h \
    Core/Optimisation.h \
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
    Gui/GraphicsItems/ViewScaler.h \
    Core/FileManipulation.h \
    Core/ProcessManager.h \
    Gui/GuiComponents/PlotterWidget.h \
    Core/MeshDialogModel.h \
    Core/MeshModel.h \
    Gui/GraphicsItems/ColorMapper.h \
    Gui/Dialogs/windowsizing.h \
    Gui/Dialogs/ClusterLoginDialog.h \
    Core/clusterManager.h \
    Gui/Dialogs/ViewConfigureDialog.h

FORMS    += \
    Gui/DebugOutput.ui \
    Gui/Dialogs/ConfigSimulationDialog.ui \
	Gui/GuiComponents/ParallelCoordinatesWindow.ui \
    Gui/MainWindow.ui \
    Gui/Dialogs/MeshDialog.ui \
    Gui/GuiComponents/ControlPointView.ui \
    Gui/Dialogs/ClusterLoginDialog.ui \
    Gui/Dialogs/ViewConfigureDialog.ui

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

LIBS += -L/opt/mxe/usr/x86_64-w64-mingw32.static/lib -lssh
INCLUDEPATH += /opt/mxe/usr/x86_64-w64-mingw32.static/include

LIBS += -L/opt/mxe/usr/x86_64-w64-mingw32.static/lib -llibcrypto  -llibssl
INCLUDEPATH += /opt/mxe/usr/x86_64-w64-mingw32.static/include/openssl

