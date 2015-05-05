#include "vtkWriter.h"
#include <fstream>
#include <iostream>
#include <iomanip>

#include "macros.h"
#include "types.h"
#include "pack.h"

void vtkWriter::saveVTKTri( std::string filename, TMesh* tm, QMesh* qm ) {

    int i = 0;
    for( QMesh::QFaceIterator it = qm->faces_begin() ; it != qm->faces_end() ; it++ ){

        QFace* qf = *it;

        qf->get_attributes().id = i;
        i++;
    }


    filename += ".vtk";
    std::ofstream file( filename.c_str() );

    ASSERT( file.is_open() );

    unsigned nv = tm->get_number_of_vertices();

    file << "# vtk DataFile Version 1.0" << std::endl;
    file << "Mesh from of" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS " << nv << " float" << std::endl;

    i = 0; // sequencial index for the vertex

    for( TMesh::VertexIterator it = tm->vertices_begin(); it != tm->vertices_end(); ++it ) {
        TMesh::Vertex* v = *it;

        file << std::setprecision(15) << v->x();
        file << " " << std::setprecision(15) << v->y();
        file << " " << std::setprecision(15) << v->z();
        file << std::endl;
        v->get_attributes().setId(i); // setting an sequencial Id for the vertex, method setId() is located on types.h
        i++;
    }

    unsigned nf = tm->get_number_of_faces();
    file << std::endl << "CELLS " << nf << " " << 4*nf << std::endl;


    for( TMesh::FaceIterator it = tm->faces_begin(); it != tm->faces_end(); ++it ) {

        TMesh::Face* f = *it;

        file << "3";
        file << " " << f->get_halfedge()->get_origin()->get_attributes().getId();
        file << " " << f->get_halfedge()->get_next()->get_origin()->get_attributes().getId();
        file << " " << f->get_halfedge()->get_prev()->get_origin()->get_attributes().getId();
        file << std::endl;
    }

    file << std::endl << "CELL_TYPES " << nf << std::endl;
    for (unsigned i=0; i < nf; i++)
        file << "5 ";
    file << std::endl;
    file << std::endl;

    file << std::endl << "CELL_DATA " << nf << std::endl;
    file << "SCALARS region int 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for( TMesh::FaceIterator it = tm->faces_begin(); it != tm->faces_end(); ++it ) {

        TMesh::Face* f = *it;

        if( f->get_attributes().region.size() > 1 ) {
            file << "-1 ";
        } else if( f->get_attributes().region.size() == 1 ) {
            QFace* qf = (QFace*)(*f->get_attributes().region.begin());
            int r = qf->get_attributes().id;
            file << r << " ";
        } else {
            file << "-2 ";
        }
    }
    file << std::endl;
    file << std::endl;

    file << "SCALARS tri_type int 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for( TMesh::FaceIterator it = tm->faces_begin(); it != tm->faces_end(); ++it ) {

        TMesh::Face* f = *it;

        if( f->get_attributes().region.size() > 1 ) {
            file << "0 ";
        } else if( f->get_attributes().region.size() == 1 ) {
            file << "1 ";
        } else {
            file << "2 ";
        }
    }
    file << std::endl;
    file << std::endl;


    file.close();

}

void vtkWriter::saveVTKQuad( std::string filename, QMesh* _refinMesh ) {

    filename += ".vtk";
    std::ofstream file( filename.c_str() );

    ASSERT( file.is_open() );

    unsigned nv = _refinMesh->get_number_of_vertices();

    file << "# vtk DataFile Version 1.0" << std::endl;
    file << "Mesh from of" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS " << nv << " float" << std::endl;

    int i = 0; // sequencial index for the vertex

    for( QMesh::QVertexIterator it = _refinMesh->vertices_begin(); it != _refinMesh->vertices_end(); ++it ) {
        QMesh::QVertex* v = *it;

        file << std::setprecision(15) << v->x();
        file << " " << std::setprecision(15) << v->y();
        file << " " << std::setprecision(15) << v->z();
        file << std::endl;
        v->get_attributes().setId(i); // setting an sequencial Id for the vertex, method setId() is located on types.h
        i++;
    }

    unsigned nf = _refinMesh->get_number_of_faces();
    file << std::endl << "CELLS " << nf << " " << 5*nf << std::endl;

    for( QMesh::QFaceIterator it = _refinMesh->faces_begin(); it != _refinMesh->faces_end(); ++it ) {

        QMesh::QFace* f = *it;

        file << "4";
        file << " " << f->get_halfedge()->get_origin()->get_attributes().getId();
        file << " " << f->get_halfedge()->get_next()->get_origin()->get_attributes().getId();
        file << " " << f->get_halfedge()->get_next()->get_next()->get_origin()->get_attributes().getId();
        file << " " << f->get_halfedge()->get_prev()->get_origin()->get_attributes().getId();
        file << std::endl;
    }

    file << std::endl << "CELL_TYPES " << nf << std::endl;
    for (unsigned i=0; i < nf; i++)
        file << "9 ";
    file << std::endl;

    file << std::endl << "CELL_DATA " << nf << std::endl;
    file << "SCALARS subdivide int 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for( QMesh::QFaceIterator it = _refinMesh->faces_begin(); it != _refinMesh->faces_end(); ++it ) {

        QMesh::QFace* f = *it;


        file << f->get_attributes().getTemplateType() << " ";

    }

    file.close();

    file << std::endl;
    file << std::endl;
}


void vtkWriter::saveVTKTriPatch( std::string filename, unsigned nv, double* vset, unsigned nf, unsigned* fset, double* tvalue ) {

    filename += ".vtk";
    std::ofstream file( filename.c_str() );

    ASSERT( file.is_open() );

    file << "# vtk DataFile Version 1.0" << std::endl;
    file << "Mesh from of" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS " << nv << " float" << std::endl;


    for(unsigned i = 0 ; i< nv; i++){

        file << std::setprecision(15) << vset[3*i];
        file << " " << std::setprecision(15) << vset[3*i+1];
        file << " " << std::setprecision(15) << vset[3*i+2];
        file << std::endl;
    }


    file << std::endl << "CELLS " << nf << " " << 4*nf << std::endl;

    for( unsigned i = 0 ; i< nf; i++ ) {


        file << "3";
        file << " " << fset[3*i];
        file << " " << fset[3*i+1];
        file << " " << fset[3*i+2];
        file << std::endl;
    }

    file << std::endl << "CELL_TYPES " << nf << std::endl;
    for (unsigned i=0; i < nf; i++)
        file << "5 ";
    file << std::endl;
    file << std::endl;

    file << std::endl << "POINT_DATA " << nv << std::endl;
    file << "SCALARS tvalue double 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for( unsigned i = 0 ; i< nv; i++  ) {

        file << tvalue[i] << " ";

    }
    file << std::endl;
    file << std::endl;


    file.close();

}

void vtkWriter::saveVTKTriPatchWithoutTvalue( std::string filename, unsigned nv, double* vset, unsigned nf, unsigned* fset) {

    filename += ".vtk";
    std::ofstream file( filename.c_str() );

    ASSERT( file.is_open() );

    file << "# vtk DataFile Version 1.0" << std::endl;
    file << "Mesh from of" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS " << nv << " float" << std::endl;


    for(unsigned i = 0 ; i< nv; i++){

        file << std::setprecision(15) << vset[3*i];
        file << " " << std::setprecision(15) << vset[3*i+1];
        file << " " << std::setprecision(15) << vset[3*i+2];
        file << std::endl;
    }


    file << std::endl << "CELLS " << nf << " " << 4*nf << std::endl;

    for( unsigned i = 0 ; i< nf; i++ ) {


        file << "3";
        file << " " << fset[3*i];
        file << " " << fset[3*i+1];
        file << " " << fset[3*i+2];
        file << std::endl;
    }

    file << std::endl << "CELL_TYPES " << nf << std::endl;
    for (unsigned i=0; i < nf; i++)
        file << "5 ";
    file << std::endl;
    file << std::endl;


    file.close();

}


void vtkWriter::saveVTKGrid( std::string filename, unsigned dimensionX, unsigned dimensionY, unsigned dimensionZ, unsigned npoints, double* points) {

    filename += ".vtk";
    std::ofstream file( filename.c_str() );

    ASSERT( file.is_open() );

    file << "# vtk DataFile Version 1.0" << std::endl;
    file << "Mesh from of" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET STRUCTURED_GRID" << std::endl;
    file << "DIMENSIONS " << dimensionX << dimensionY << dimensionZ << std::endl;
    file << "POINTS " << npoints << " float" << std::endl;


    for(unsigned i = 0 ; i< npoints; i++){

        file << std::setprecision(15) << points[3*i];
        file << " " << std::setprecision(15) << points[3*i+1];
        file << " " << std::setprecision(15) << points[3*i+2];
        file << std::endl;
    }


    file.close();

}
