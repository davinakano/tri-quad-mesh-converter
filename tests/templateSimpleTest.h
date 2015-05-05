#ifndef TEMPLATESIMPLETEST_H
#define TEMPLATESIMPLETEST_H

#include "types.h"
#include "vtkReader.h"
#include "stdlib.h"
#include "meshref.h"

class templateSimpleTest
{
public :

    // ****************************************
    //              TESTE TEMPLATE 1
    // ****************************************

    /**
      * \fn void testeTemplate1(QMesh* qm)
      *
      * \brief Simple test of the first split template, running only on the six
      *        initial faces made by initialSubdivision method.
      *
      * \param qm Pointer of a quadrilateral mesh.
      */

    void testeTemplate1(QMesh* qm){

        QOperators qop(qm);

        std::list< QHalfedge* > _list_of_halfedges;

        QHalfedge* qh = (*qm->faces_begin())->get_halfedge();

        for( int i = 0; i <= 3; i++ ) {
            _list_of_halfedges.push_back(qh);
            qh = qh->get_next()->get_mate()->get_next();
        }

        for( std::list< QHalfedge* >::iterator it = _list_of_halfedges.begin(); it != _list_of_halfedges.end(); it++ ) {
            QHalfedge* qh = *it;
            //qop.split_template1(qh);
        }
    }

    // ****************************************
    //              TESTE TEMPLATE 2
    // ****************************************

    /**
      * \fn void testeTemplate2(QMesh* qm)
      *
      * \brief Simple test of the second split template, running only on the six
      *        initial faces made by initialSubdivision method.
      *
      * \param qm Pointer of a quadrilateral mesh.
      */

    void testeTemplate2(QMesh* qm){

        QOperators qop(qm);

        std::list< QFace* > _list_of_faces;

        for( QMesh::QFaceIterator it = qm->faces_begin(); it != qm->faces_end(); it++ ) {
            QFace* qf = *it;

            _list_of_faces.push_back(qf);
        }

        for( std::list< QFace* >::iterator it = _list_of_faces.begin(); it != _list_of_faces.end(); it++ ) {
            QFace* qhf = *it;
            //qop.split_template2(qhf->get_halfedge(), 0);
        }

    }

    // ****************************************
    //              TESTE TEMPLATE 3
    // ****************************************

    /**
      * \fn void testeTemplate3(QMesh* qm)
      *
      * \brief Simple test of the third split template, running only on the six
      *        initial faces made by initialSubdivision method.
      *
      * \param qm Pointer of a quadrilateral mesh.
      */

    void testeTemplate3(QMesh* qm){

        QOperators qop(qm);

        std::list< QHalfedge* > _list_of_halfedges;

        QHalfedge* qh = (*qm->faces_begin())->get_halfedge();

        for( int i = 0; i <= 2; i++ ) {
            _list_of_halfedges.push_back(qh);
            qh = qh->get_mate()->get_next();
        }

        for( std::list< QHalfedge* >::iterator it = _list_of_halfedges.begin(); it != _list_of_halfedges.end(); it++ ) {
            QHalfedge* qh = *it;
            qop.split_template3(qh);
            //break;
        }

    }
}

#endif // TEMPLATESIMPLETEST_H
