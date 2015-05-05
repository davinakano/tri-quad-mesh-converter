DEFINES += DEBUG

DEFINES += DEBUGMODE

QMAKE_CXXFLAGS += -frounding-math

LIBS += -lumfpack \
        -lCGAL \
        -lgmp

HEADERS += \
    macros.h \
    types.h \
    dcel/edge.h \
    dcel/face.h \
    dcel/halfedge.h \
    dcel/operators.h \
    dcel/surface.h \
    dcel/vertex.h \
    qdcel/qedge.h \
    qdcel/qface.h \
    qdcel/qhalfedge.h \
    qdcel/qoperators.h \
    qdcel/qsurface.h \
    qdcel/qvertex.h \
    vtkReader.h \
    pack.h \
    geodesics/adg_curve_point.hpp \
    geodesics/adg_geodesics.hpp \
    geodesics/adg_mesh_interface.hpp \
    geodesics/geodesic_point.hpp \
    geodesics/geodesic_point_edge.hpp \
    geodesics/geodesic_point_face.hpp \
    geodesics/geodesic_point_vertex.hpp \
    geodesics/geometric.hpp \
    geodesics/mesh_interface.hpp \
    geodesics/priority_queue.hpp \
    geodesics/adg_from_dcel.hpp \
    meshref.h \
    myAdgGeodesics.h \
    meshpatch/patchvertex.hpp \
    meshpatch/patchedge.hpp \
    meshpatch/patch.hpp \
    meshpatch/mypair.hpp \
    numerics/ludcmp.hpp \
    geodesicMeshing/geodesicMeshing.h \
    floater/floaterpar_vertex_attributes.hpp \
    floater/floaterpar.hpp \
    bezier/tbezier.h \
    lssolver/tludcmp.h \
    tools.h \
    templateSimpleTest.h \
    tests/templateSimpleTest.h \
    vtkWriter.h \
    debug.h \
    cgalmethods.h

SOURCES += \
    main.cpp \
    vtkReader.cpp \
    myAdgGeodesics.cpp \
    numerics/ludcmp.cpp \
    geodesicMeshing/geodesicMeshing.cpp \
    floater/floaterpar.cpp \
    bezier/tbezier.cpp \
    lssolver/tludcmp.cpp \
    vtkWriter.cpp \
    cgalmethods.cpp
