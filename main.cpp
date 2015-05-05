#include "types.h"
#include "vtkReader.h"
#include "vtkWriter.h"
#include "meshref.h"
#include "geodesicMeshing/geodesicMeshing.h"
#include "floater/floaterpar.hpp"
#include "bezier/tbezier.h"
#include "tools.h"
#include "debug.h"
#include <sstream>

int main( int argc, char* argv[] ) {

    if( argc != 4 ) {
        return 3;
    }

    std::string filename = argv[1];

    int error_size = atoi(argv[2]);
    int function = atoi(argv[3]);

    vtkReader reader;
    reader.loadMesh( filename );

    TMesh* tm;
    tm = reader.getMeshTri();

    QMesh* qm;
    qm = reader.getMeshQuad();

    adgFDCel adg(tm);

    geodesicMeshing* geoMeshing = new geodesicMeshing(tm);

    QOperators qop(qm);

    QMesh::QEdgeIterator it = qm->edges_begin();

    QEdge* e = *it;

    tools T(error_size);

    T.initialDivision( tm, qm, geoMeshing );

    debug::iteration = 0;

    // Subdividing

    bool continua = T.templateMarkFaces(qm);

    while ( continua ){


        std::stringstream number;
        number << debug::iteration;

        vtkWriter::saveVTKQuad( filename + std::string("-quad-pre-") + number.str(), qm);

        switch(function) {
        case 1: T.correctList(qm);
            break;
        case 2: T.correctList2(qm);
            break;
        case 3: T.correctList3(qm);
            break;
        }



        vtkWriter::saveVTKQuad( filename + std::string("-quad-pos-") + number.str(), qm);

        T.subdivide( tm, qm, geoMeshing);

        vtkWriter::saveVTKTri(filename + std::string("-tri-")+ number.str(), tm, qm);

        debug::iteration++;

        continua = T.templateMarkFaces(qm);

    }

    std::stringstream number;
    number << debug::iteration;

    vtkWriter::saveVTKQuad( filename + std::string("-final-") + number.str(), qm);

    return 1;
}
