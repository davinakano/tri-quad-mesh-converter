#include <fstream>
#include <iostream>
#include <iomanip>

#include "macros.h"
#include "vtkReader.h"
#include "types.h"
#include "pack.h"

vtkReader::vtkReader()
{
    _MeshTri = 0;
    _MeshQuad = 0;
}

vtkReader::~vtkReader() {

    if( _MeshTri )
        delete _MeshTri;

    if( _MeshQuad )
        delete _MeshQuad;
}

TMesh* vtkReader::getMeshTri() {
    return _MeshTri;
}

QMesh* vtkReader::getMeshQuad() {
    return _MeshQuad;
}

bool vtkReader::loadMesh(std::string filename ) {

    ASSERT( !filename.empty() );
    ASSERT( _MeshTri == 0 );
    ASSERT( _MeshQuad == 0 );

    std::string nfilename = filename + std::string(".vtk");

    std::string auxstr;
    std::ifstream file( nfilename.c_str() );

    ASSERT( file.is_open() );

    file >> auxstr;
    while( auxstr.compare( "POINTS" ) != 0 )
        file >> auxstr;

    unsigned r_nv;
    file >> r_nv;

    ASSERT ( r_nv > 3 );

    // First step : Read Vertex -> r_coordsTri[]

    double* r_coordsTri = new double[ 3*r_nv ]; // Coordinates Vector

    file >> auxstr;

    for( unsigned i=0; i < r_nv; i++) {
        file >> r_coordsTri[3*i + 0];
        file >> r_coordsTri[3*i + 1];
        file >> r_coordsTri[3*i + 2];
    }

    file >> auxstr;
    while( auxstr.compare( "CELLS" ) != 0 )
        file >> auxstr;

    unsigned r_nc;
    int temp;
    file >> r_nc;

    ASSERT( r_nc >= 1 );

    file >> temp;

    //std::vector<Vertex*> vertexTri;
    //std::vector<Vertex*>::iterator itvTri;

    std::vector<int> vertexTri;
    std::vector<int>::iterator itvTri;

    std::vector<int> vertexQuad;
    //std::vector<int>::iterator itvQuad;

    // Second step : Read Cells -> Vector < Triangles > Vector < Quad >

    for( unsigned i=0; i<r_nc; i++) {

        file >> temp;
        if ( temp == 3 ){

            for ( int i = 1 ; i <= 3 ; i++ ){

                file >> temp;
                vertexTri.push_back(temp);
            }
        }

        else if ( temp == 4 ){

            for ( int j = 1 ; j <= 4 ; j++ ){

                file >> temp;
                vertexQuad.push_back(temp);
            }
        }

        else { ASSERT ( false ); };
    }

    // Third step : Transf. Vector < Triangles > -> rcells

    unsigned* r_cellsTri = new unsigned[ vertexTri.size() ];

    itvTri = vertexTri.begin();

    for( unsigned i=0 ; i < vertexTri.size() / 3 ; i++ ){

        r_cellsTri[3*i + 0] = *itvTri;
        itvTri++;
        r_cellsTri[3*i + 1] = *itvTri;
        itvTri++;
        r_cellsTri[3*i + 2] = *itvTri;
        itvTri++;
    }

    // Fourth step : Creation of the DCELL < Triangles >

    _MeshTri = new TMesh( r_nv, r_coordsTri, vertexTri.size() / 3, r_cellsTri );

    // Fifth step : Tranf. Vector < Quad > -> set pack2<int,int>

    std::set<pack2int> p2iSet;

    for( unsigned i=0 ; i < vertexQuad.size() ; i++ ){

        p2iSet.insert( pack2<int,int>(vertexQuad[i],p2iSet.size()) );
    }

    // Sixth step : Tranf. Inverted Set

    std::set<pack2int> p2iInvertedSet;

    for( std::set<pack2int>::iterator it = p2iSet.begin() ; it != p2iSet.end() ; it++ ){

        pack2int p = *it;

        p2iInvertedSet.insert( pack2<int,int>(p.second,p.first) );

    }

    // Seventh step : Tranf. Inverted Set -> r_coords

    double* r_coordsQuad = new double[ 3*p2iInvertedSet.size() ]; // Coordinates Vector of Quad's Vertex

    int i=0;

    for ( std::set<pack2int>::iterator it = p2iInvertedSet.begin() ; it != p2iInvertedSet.end(); it++ ){

        pack2int p = *it;

        ASSERT( i == p.first );

        ASSERT ( p.second < (int)r_nv );

        ASSERT ( p.second >= 0 );

        r_coordsQuad[3*i] = r_coordsTri[3*(p.second)];
        r_coordsQuad[3*i + 1] = r_coordsTri[3*(p.second) + 1];
        r_coordsQuad[3*i + 2] = r_coordsTri[3*(p.second) + 2];

        i++;
    }

    // Eighth step : Tranf. Vector < Quad > -> r_cells

    ASSERT ( vertexQuad.size() >= 4 );

    unsigned* r_cellsQuad = new unsigned[ vertexQuad.size() ];

    std::set<pack2int>::iterator it;

    for( unsigned i=0 ; i < vertexQuad.size() ; i++ ){

        it = p2iSet.find( pack2int(vertexQuad[i],0) );

        if ( it != p2iSet.end() ){

            pack2int p2i = *it;

            ASSERT ( p2i.first == vertexQuad[i] );

            r_cellsQuad[i] = p2i.second;

        }

        else { ASSERT ( false ); };
    }

    _MeshQuad = new QMesh( p2iSet.size(), r_coordsQuad, vertexQuad.size() / 4 , r_cellsQuad );

    // Making links between the Vertex and QVertex

    std::vector<Vertex*> TInvRef;
    std::vector<QVertex*> QInvRef;

    for( TMesh::VertexIterator it = _MeshTri->vertices_begin() ; it != _MeshTri->vertices_end() ; it++ ){

        Vertex* v = *it;
        TInvRef.push_back(v);

    }

    for( QMesh::QVertexIterator it = _MeshQuad->vertices_begin() ; it != _MeshQuad->vertices_end() ; it++ ){

        QVertex* v = *it;
        QInvRef.push_back(v);

    }

    for ( std::set<pack2int>::iterator it = p2iSet.begin() ; it != p2iSet.end() ; it++ ){

        pack2int p = *it;

        ASSERT ( p.first >= 0 );
        ASSERT ( p.first < (int)TInvRef.size() );
        ASSERT ( p.second >= 0 );
        ASSERT ( p.second < (int)QInvRef.size() );

        //TInvRef[p.first]->get_attributes().setQVertex(QInvRef[p.second]);
        //QInvRef[p.second]->get_attributes().setTVertex(TInvRef[p.first]);
        QInvRef[p.second]->get_attributes().setGP(new geodesics::GPV(TInvRef[p.first]));

    }

    i = 0;
    for( QMesh::QFaceIterator it = _MeshQuad->faces_begin() ; it != _MeshQuad->faces_end() ; it++ ){

        QFace* qf = *it;

        qf->get_attributes().id = i;
        i++;
    }


    file.close();

    return true;
}
