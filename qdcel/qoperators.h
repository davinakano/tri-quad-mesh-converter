/**
 * \file operators.h
 *
 * \brief Definition  of the class Operators, which performs some topological
 * operations in a Surface mesh.
 *
 * \author
 * Mario Augusto de Souza Lizier \n
 * Davi Yoshinori Cangussú Nakano \n
 * Universidade Federal de São Carlos, \n
 * Departamento de Computação, \n
 *
 * \version 1.0
 * \date Jan 2010
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

#ifndef QOPERATORS_H
#define	QOPERATORS_H

#include "qvertex.h"    // Vertex
#include "qhalfedge.h"  // Halfedge
#include "qedge.h"      // Edge
#include "qface.h"      // Face
#include "qsurface.h"   // Surface
#include "dcel/vertex.h"
#include <list>
#include <vector>

/**
 * \defgroup DCELNameSpace Namespace adcel.
 * @{
 */

/**
 * \namespace qdcel
 *
 * \brief   The   namespace   adcel   contains   the   definition   and
 * implementation of all classes of the data structure Doubly
 * Connected  Edge List  (DCEL), which  can be  used to  represent the
 * mesh.
 *
 */

namespace qdcel {

    /**
     * \class QOperators
     *
     * \brief This class performs some topological operations in a Surface mesh.
     */
    template <
            typename VAttrib = int,
            typename FAttrib = int,
            typename EAttrib = int,
            typename HAttrib = int
            >
            class QOperators {
    public:

        /**
         * \typedef Vertex
         *
         * \brief Defines  Vertex as  an alias for  dcel::Vertex< VAttrib,
         * FAttrib , EAttrib , HAttrib >.
         */
        typedef qdcel::QVertex< VAttrib, FAttrib, EAttrib, HAttrib > QVertex;

        /**
         * \typedef Face
         *
         * \brief  Defines  Face  as  an alias  for  dcel::Face<  VAttrib,
         * FAttrib , EAttrib , HAttrib >.
         */
        typedef qdcel::QFace< VAttrib, FAttrib, EAttrib, HAttrib > QFace;

        /**
         * \typedef Edge
         *
         * \brief  Defines  Edge  as  an alias  for  dcel::Edge<  VAttrib,
         * FAttrib , EAttrib , HAttrib >.
         */
        typedef qdcel::QEdge< VAttrib, FAttrib, EAttrib, HAttrib > QEdge;

        /**
         * \typedef Halfedge
         *
         * \brief  Defines  Halfedge   as  an  alias  for  dcel::Halfedge<
         * VAttrib, FAttrib , EAttrib , HAttrib >.
         */
        typedef qdcel::QHalfedge< VAttrib, FAttrib, EAttrib, HAttrib > QHalfedge;


        /**
         * \typedef Surface
         *
         * \brief Defines Surface as an alias for dcel::Surface < VAttrib,
         *  FAttrib , EAttrib , HAttrib >.
         */
        typedef qdcel::QSurface< VAttrib, FAttrib, EAttrib, HAttrib > QSurface;

        typedef dcel::Vertex< VAttrib, FAttrib, EAttrib, HAttrib > Vertex;

        QOperators(QSurface *_mesh);

        virtual ~QOperators();

        // bool split_triangle(Face *f, double x, double y, double z, Face *&f0, Face *&f1, Face *&f2);

        // void split_quad(Halfedge *h, double x, double y, double z, Face *&f0, Face *&f1, Face *&f2, Face *&f3);

        std::list<qdcel::QVertex< VAttrib, FAttrib, EAttrib, HAttrib >*> split_strip(QHalfedge *h, int nPoints, double *points);

        void split_template1(QHalfedge *he, std::vector<QVertex*> *createdVertex);

        void split_template2(QHalfedge *he, std::vector<QVertex*> *createdVertex);

        void split_template3(QHalfedge *he, std::vector<QVertex*> *createdVertex);



    private:

        void split_quad(QHalfedge *he, QVertex *nv1, QVertex *nv2, QEdge *e0, QEdge *e1);

        QVertex* add_vertex(double x, double y, double z, QHalfedge *h = 0) {

            QVertex* vertex = new QVertex(x, y, z, h);

            mesh->add_vertex(vertex);

            return vertex;
        }

        QFace* add_face(QHalfedge *h = 0) {

            QFace* face = new QFace(h);

            mesh->add_face(face);

            return face;
        }

        QEdge* add_edge(QHalfedge *h1 = 0, QHalfedge *h2 = 0) {
            
            QEdge* edge = new QEdge(h1, h2);

            mesh->add_edge(edge);

            return edge;
        }


        QSurface *mesh;
    };


    template <
            typename VAttrib,
            typename FAttrib,
            typename EAttrib,
            typename HAttrib
            >
    QOperators<VAttrib,FAttrib,EAttrib,HAttrib>::QOperators( QSurface *_mesh) : mesh(_mesh)
    {
    }

    template <
            typename VAttrib,
            typename FAttrib,
            typename EAttrib,
            typename HAttrib
            >
    QOperators<VAttrib,FAttrib,EAttrib,HAttrib>::~QOperators()
    {

    }

    template <
            typename VAttrib,
            typename FAttrib,
            typename EAttrib,
            typename HAttrib
            >
    void QOperators<VAttrib,FAttrib,EAttrib,HAttrib>::split_quad(QHalfedge *he, QVertex *nv1, QVertex *nv2, QEdge *e0, QEdge *e1) {

        // Half-Edges inside main quad
        QHalfedge *he1 = he->get_next();
        QHalfedge *he2 = he->get_next()->get_next();
        //QHalfedge *he3 = he->get_prev();

        // The four initial Vertex
        //QVertex *v0 = he->get_origin();
        QVertex *v1 = he1->get_origin();
        QVertex *v2 = he2->get_origin();
        //QVertex *v3 = he3->get_origin();

        QFace *f;

        // Creating the new Face - 'hNewFace' stands for Half Edges of New Face
        f = add_face();
        QHalfedge * hNewFace[4];
        hNewFace[0] = new QHalfedge(nv1, 0, f, 0, 0);
        hNewFace[1] = new QHalfedge(v1, 0, f, 0, hNewFace[0]);
        hNewFace[2] = new QHalfedge(v2, 0, f, 0, hNewFace[1]);
        hNewFace[3] = new QHalfedge(nv2, 0, f, hNewFace[0], hNewFace[2]);

        hNewFace[0]->set_next(hNewFace[1]);
        hNewFace[0]->set_prev(hNewFace[3]);

        hNewFace[1]->set_next(hNewFace[2]);

        hNewFace[2]->set_next(hNewFace[3]);

        f->set_halfedge(hNewFace[0]);

        // Creating the new Edge
        QEdge *newEdge;

        newEdge = add_edge(hNewFace[3], he1 );

        // Fixing the 'Half-Edge's to Edges' pointers
        hNewFace[0]->set_edge( e0 );
        hNewFace[1]->set_edge( he1->get_edge() );
        hNewFace[2]->set_edge( e1 );
        hNewFace[3]->set_edge( newEdge );

        // Fixing the bottom edge's first and/or second half edges
        if ( he1->get_edge()->get_first_halfedge() == he1 )
                he1->get_edge()->set_first_halfedge( hNewFace[1] );
        else
                he1->get_edge()->set_second_halfedge( hNewFace[1] );

        he1->set_edge( newEdge );

        if ( e0->get_first_halfedge() == 0 ) // e0 is a parameter of this function
            e0->set_first_halfedge( hNewFace[0] );
        else
            e0->set_second_halfedge( hNewFace[0] );

        if ( e1->get_first_halfedge() == 0 ) // e1 is a parameter of this function
            e1->set_first_halfedge( hNewFace[2] );
        else
            e1->set_second_halfedge( hNewFace[2] );

        nv1->set_halfedge( hNewFace[0] );
        v1->set_halfedge( hNewFace[1] );
        v2->set_halfedge( hNewFace[2] );
        nv2->set_halfedge( hNewFace[3] );

        he1->set_origin( nv1 );
        he2->set_origin( nv2 );

    }

    // Type 1 Template ( D - D - Straight)
    template <
            typename VAttrib,
            typename FAttrib,
            typename EAttrib,
            typename HAttrib
            >
    void QOperators<VAttrib,FAttrib,EAttrib,HAttrib>::split_template1(QHalfedge *he, std::vector<QVertex*> *createdVertex) {

        QVertex *v[4];
        QEdge   *newEdges[7];
        QEdge   *oldEdges[4];
        QHalfedge *oldQuad[4];
        QHalfedge *newQuad0[4];
        QHalfedge *newQuad1[4];
        QHalfedge *outsideH[6];
        QHalfedge *heSearch;
        QFace *newFace0;
        QFace *newFace1 = add_face();
        QFace *oldFace = he->get_face();
        int i = 0;

        newFace0 = oldFace;

        // Setting all old Vertex to v[]
        // Setting all old Half Edges inside of the Main Quad
        // Setting all old Edges inside of the Main QUad

        for (; i < 4 ; i++){
            v[i]        = he->get_origin();
            oldQuad[i]  = he;
            oldEdges[i] = he->get_edge();
            he          = he->get_next();
        }

        // Setting the old Half Edges to their respective Quads
        newQuad0[0] = oldQuad[0];
        newQuad0[0]->set_origin(v[0]);

        newQuad0[1] = oldQuad[1];
        newQuad0[1]->set_origin(v[1]);

        newQuad1[2] = oldQuad[2];
        newQuad1[2]->set_origin(v[2]);

        newQuad1[3] = oldQuad[3];
        newQuad1[3]->set_origin(v[3]);

        // Setting the old Edges to their respective Half Edges
        newEdges[0] = oldEdges[0];
        newEdges[1] = oldEdges[1];
        newEdges[3] = oldEdges[2];
        newEdges[4] = oldEdges[3];

        outsideH[0] = he->get_mate();
        outsideH[3] = he->get_next()->get_next()->get_mate();

        // Top Edge is subdivided?
        if (he->get_prev()->get_edge()->get_attributes().subdivide == true){

            heSearch = he->get_prev()->get_mate()->get_next();
            while (heSearch->get_next()->get_origin() != he->get_prev()->get_origin()) {
                heSearch = heSearch->get_mate();
                heSearch = heSearch->get_next();
            }

            newEdges[4]  = heSearch->get_edge();
            newEdges[5]  = he->get_prev()->get_edge();
            outsideH[4]  = heSearch;
            outsideH[5]  = he->get_prev()->get_mate();

        } else {

                    newEdges[4]  = he->get_prev()->get_edge();
                    newEdges[5]  = add_edge();
                    outsideH[4]  = he->get_prev()->get_mate();
                    outsideH[5]  = 0;
        }


        // Bottom Edge is subdivided?
        if (he->get_next()->get_edge()->get_attributes().subdivide == true){

            heSearch = he->get_next()->get_mate()->get_next();
            while (heSearch->get_next()->get_origin() != he->get_next()->get_origin()) {
                heSearch = heSearch->get_mate();
                heSearch = heSearch->get_next();
            }

            newEdges[1] = heSearch->get_edge();
            newEdges[2] = he->get_next()->get_edge();
            outsideH[1] = heSearch;
            outsideH[2] = he->get_next()->get_mate();

        } else {

                    newEdges[1] = he->get_next()->get_edge();
                    newEdges[2] = add_edge();
                    outsideH[1] = he->get_next()->get_mate();
                    outsideH[2] = 0;
               }

        newEdges[6]  = add_edge(); // Middle Edge

        // Left Quad
        // newQuad0[0] and newQuad0[1] were previously created
        newQuad0[2] = new QHalfedge(createdVertex->at(0), newEdges[6], newFace0, 0, newQuad0[1]);
        newQuad0[3] = new QHalfedge(createdVertex->at(1), newEdges[5], newFace0, newQuad0[0], newQuad0[2]);

        newQuad0[2]->set_next(newQuad0[3]);

        newQuad0[1]->set_prev(newQuad0[0]);
        newQuad0[1]->set_next(newQuad0[2]);

        newQuad0[0]->set_prev(newQuad0[3]);
        newQuad0[0]->set_next(newQuad0[1]);

        // Right Quad
        // newQuad1[2] and newQuad1[3] were previously created
        newQuad1[0] = new QHalfedge(createdVertex->at(1), newEdges[6], newFace1, 0, newQuad1[3]);
        newQuad1[1] = new QHalfedge(createdVertex->at(0), newEdges[2], newFace1, newQuad1[2], newQuad1[3]);

        newQuad1[0]->set_next(newQuad1[1]);
        newQuad1[0]->set_prev(newQuad1[3]);

        newQuad1[2]->set_next(newQuad1[3]);
        newQuad1[2]->set_prev(newQuad1[1]);

        newQuad1[3]->set_next(newQuad1[0]);
        newQuad1[3]->set_prev(newQuad1[2]);

        // ****** Edges to Half-Edges pointers *******************

        newEdges[0]->set_first_halfedge(newQuad0[0]);
        newEdges[0]->set_second_halfedge(outsideH[0]);

        newEdges[1]->set_first_halfedge(newQuad0[1]);
        newEdges[1]->set_second_halfedge(outsideH[1]);

        newEdges[2]->set_first_halfedge(newQuad1[1]);
        newEdges[2]->set_second_halfedge(outsideH[2]);

        newEdges[3]->set_first_halfedge(newQuad1[2]);
        newEdges[3]->set_second_halfedge(outsideH[3]);

        newEdges[4]->set_first_halfedge(newQuad1[3]);
        newEdges[4]->set_second_halfedge(outsideH[4]);

        newEdges[5]->set_first_halfedge(newQuad0[3]);
        newEdges[5]->set_second_halfedge(outsideH[5]);

        newEdges[6]->set_first_halfedge(newQuad0[2]);
        newEdges[6]->set_second_halfedge(newQuad1[0]);

        // ************************************************************

        // ********* Half-edges to Edges/Faces pointers ****************

        newQuad0[0]->set_face(newFace0);
        newQuad0[0]->set_edge(newEdges[0]);

        newQuad0[1]->set_face(newFace0);
        newQuad0[1]->set_edge(newEdges[1]);

        newQuad1[2]->set_face(newFace1);
        newQuad1[2]->set_edge(newEdges[3]);

        newQuad1[3]->set_face(newFace1);
        newQuad1[3]->set_edge(newEdges[4]);

        // ************************************************************

        // ********* Faces pointers ****************

        newFace0->set_halfedge(newQuad0[0]);
        newFace1->set_halfedge(newQuad1[0]);

        // ************************************************************

        // ********* Marking the edges ****************

        newEdges[1]->get_attributes().subdivide = true;
        newEdges[2]->get_attributes().subdivide = true;
        newEdges[4]->get_attributes().subdivide = true;
        newEdges[5]->get_attributes().subdivide = true;

        // ************************************************************

        outsideH[0]->set_edge(newEdges[0]);

        outsideH[1]->set_edge(newEdges[1]);

        if (outsideH[2])
            outsideH[2]->set_edge(newEdges[2]);

        outsideH[3]->set_edge(newEdges[3]);

        outsideH[4]->set_edge(newEdges[4]);

        if (outsideH[5])
            outsideH[5]->set_edge(newEdges[5]);

        createdVertex->at(0)->set_halfedge(newQuad0[2]);
        createdVertex->at(1)->set_halfedge(newQuad1[1]);

    }

    // Type 2 Template ( D - D - D - D )
    template <
            typename VAttrib,
            typename FAttrib,
            typename EAttrib,
            typename HAttrib
            >
    void QOperators<VAttrib,FAttrib,EAttrib,HAttrib>::split_template2(QHalfedge *he, std::vector<QVertex*> *createdVertex) {

        QVertex *v[4];
        QEdge   *newEdges[12];
        QEdge   *oldEdges[4];
        QHalfedge *oldQuad[4];
        QHalfedge *newQuad0[4];
        QHalfedge *newQuad1[4];
        QHalfedge *newQuad2[4];
        QHalfedge *newQuad3[4];
        QHalfedge *outsideH[8];
        QHalfedge *heSearch;
        QFace *newFace0;
        QFace *newFace1 = add_face();
        QFace *newFace2 = add_face();
        QFace *newFace3 = add_face();
        QFace *oldFace  = he->get_face();
        int i = 0;

        newFace0 = oldFace;

        // Setting all old Vertex to v[]
        // Setting all old Half Edges inside of the Main Quad
        // Setting all old Edges inside of the Main QUad

        for (; i < 4 ; i++){
            v[i]        = he->get_origin();
            oldQuad[i]  = he;
            oldEdges[i] = he->get_edge();
            he          = he->get_next();
        }

        // Left Edge is subdivided?
        if (he->get_edge()->get_attributes().subdivide == true){

            heSearch = he->get_mate()->get_next();
            while (heSearch->get_next()->get_origin() != he->get_origin()) {
                heSearch = heSearch->get_mate()->get_next();
                assert( heSearch != he );
            }

            newEdges[0]  = heSearch->get_edge();
            newEdges[1]  = he->get_edge();
            outsideH[0]  = heSearch;
            outsideH[1]  = he->get_mate();

        } else {

                    newEdges[0] = he->get_edge();
                    newEdges[1] = add_edge();
                    outsideH[0] = he->get_mate();
                    outsideH[1] = 0;
                }

        // Bottom Edge is subdivided?
        if (he->get_next()->get_edge()->get_attributes().subdivide == true){

            heSearch = he->get_next()->get_mate()->get_next();
            while (heSearch->get_next()->get_origin() != he->get_next()->get_origin()) {
                heSearch = heSearch->get_mate();
                heSearch = heSearch->get_next();
            }

            newEdges[2] = heSearch->get_edge();
            newEdges[3] = he->get_next()->get_edge();
            outsideH[2] = heSearch;
            outsideH[3] = he->get_next()->get_mate();

        } else {

                    newEdges[2] = he->get_next()->get_edge();
                    newEdges[3] = add_edge();
                    outsideH[2] = he->get_next()->get_mate();
                    outsideH[3] = 0;
               }

        // Right Edge is subdivided?
        if (he->get_next()->get_next()->get_edge()->get_attributes().subdivide == true){

            heSearch = he->get_next()->get_next()->get_mate()->get_next();
            while (heSearch->get_next()->get_origin() != he->get_next()->get_next()->get_origin()) {
                heSearch = heSearch->get_mate();
                heSearch = heSearch->get_next();
            }

            newEdges[4] = heSearch->get_edge();
            newEdges[5] = he->get_next()->get_next()->get_edge();
            outsideH[4] = heSearch;
            outsideH[5] = he->get_next()->get_next()->get_mate();

        }  else {

                    newEdges[4] = he->get_next()->get_next()->get_edge();
                    newEdges[5] = add_edge();
                    outsideH[4] = he->get_next()->get_next()->get_mate();
                    outsideH[5] = 0;
                }

        // Top Edge is subdivided?
        if (he->get_prev()->get_edge()->get_attributes().subdivide == true){

            heSearch = he->get_prev()->get_mate()->get_next();
            while (heSearch->get_next()->get_origin() != he->get_prev()->get_origin()) {
                heSearch = heSearch->get_mate();
                heSearch = heSearch->get_next();
            }

            newEdges[6]  = heSearch->get_edge();
            newEdges[7]  = he->get_prev()->get_edge();
            outsideH[6]  = heSearch;
            outsideH[7]  = he->get_prev()->get_mate();

        } else {

                    newEdges[6]  = he->get_prev()->get_edge();
                    newEdges[7]  = add_edge();
                    outsideH[6]  = he->get_prev()->get_mate();
                    outsideH[7]  = 0;
        }

        // Setting the old Half Edges to their respective Quads
        newQuad0[0] = oldQuad[0];
        newQuad0[0]->set_origin(v[0]);

        newQuad1[1] = oldQuad[1];
        newQuad1[1]->set_origin(v[1]);

        newQuad2[2] = oldQuad[2];
        newQuad2[2]->set_origin(v[2]);

        newQuad3[3] = oldQuad[3];
        newQuad3[3]->set_origin(v[3]);



        newEdges[8]  = add_edge();
        newEdges[9]  = add_edge();
        newEdges[10] = add_edge();
        newEdges[11] = add_edge();

        // Top Left Quad
        // newQuad0[0] was previously created
        newQuad0[1] = new QHalfedge(createdVertex->at(0), newEdges[8], newFace0, 0, 0);
        newQuad0[2] = new QHalfedge(createdVertex->at(4), newEdges[11], newFace0, 0, newQuad0[1]);
        newQuad0[3] = new QHalfedge(createdVertex->at(3), newEdges[7], newFace0, newQuad0[0], newQuad0[2]);

        newQuad0[2]->set_next(newQuad0[3]);
        newQuad0[1]->set_prev(newQuad0[0]);
        newQuad0[1]->set_next(newQuad0[2]);
        newQuad0[0]->set_prev(newQuad0[3]);
        newQuad0[0]->set_next(newQuad0[1]);

        // Bottom Left Quad
        // newQuad1[1] was previously created
        newQuad1[0] = new QHalfedge(createdVertex->at(0), newEdges[1], newFace1, newQuad1[1], 0);
        newQuad1[2] = new QHalfedge(createdVertex->at(1), newEdges[9], newFace1, 0, newQuad1[1]);
        newQuad1[3] = new QHalfedge(createdVertex->at(4), newEdges[8], newFace1, newQuad1[0], newQuad1[2]);

        newQuad1[0]->set_prev(newQuad1[3]);
        newQuad1[2]->set_next(newQuad1[3]);
        newQuad1[1]->set_prev(newQuad1[0]);
        newQuad1[1]->set_next(newQuad1[2]);

        // Bottom Right Quad
        // newQuad2[2] was previously created
        newQuad2[0] = new QHalfedge(createdVertex->at(4), newEdges[9], newFace2, 0, 0);
        newQuad2[1] = new QHalfedge(createdVertex->at(1), newEdges[3], newFace2, newQuad2[2], newQuad2[0]);
        newQuad2[3] = new QHalfedge(createdVertex->at(2), newEdges[10], newFace2, newQuad2[0], newQuad2[2]);

        newQuad2[0]->set_next(newQuad2[1]);
        newQuad2[0]->set_prev(newQuad2[3]);
        newQuad2[2]->set_next(newQuad2[3]);
        newQuad2[2]->set_prev(newQuad2[1]);

        // Top Right Quad
        // newQuad3[3] was previously created
        newQuad3[0] = new QHalfedge(createdVertex->at(3), newEdges[11], newFace3, 0, newQuad3[3]);
        newQuad3[1] = new QHalfedge(createdVertex->at(4), newEdges[10], newFace3, 0, newQuad3[0]);
        newQuad3[2] = new QHalfedge(createdVertex->at(2), newEdges[5], newFace3, newQuad3[3], newQuad3[1]);

        newQuad3[0]->set_next(newQuad3[1]);
        newQuad3[1]->set_next(newQuad3[2]);
        newQuad3[3]->set_next(newQuad3[0]);
        newQuad3[3]->set_prev(newQuad3[2]);

        // ****** Edges to Half-Edges pointers *******************

        newEdges[0]->set_first_halfedge(newQuad0[0]);
        newEdges[0]->set_second_halfedge(outsideH[0]);

        newEdges[1]->set_first_halfedge(newQuad1[0]);
        newEdges[1]->set_second_halfedge(outsideH[1]);

        newEdges[2]->set_first_halfedge(newQuad1[1]);
        newEdges[2]->set_second_halfedge(outsideH[2]);

        newEdges[3]->set_first_halfedge(newQuad2[1]);
        newEdges[3]->set_second_halfedge(outsideH[3]);

        newEdges[4]->set_first_halfedge(newQuad2[2]);
        newEdges[4]->set_second_halfedge(outsideH[4]);

        newEdges[5]->set_first_halfedge(newQuad3[2]);
        newEdges[5]->set_second_halfedge(outsideH[5]);

        newEdges[6]->set_first_halfedge(newQuad3[3]);
        newEdges[6]->set_second_halfedge(outsideH[6]);

        newEdges[7]->set_first_halfedge(newQuad0[3]);
        newEdges[7]->set_second_halfedge(outsideH[7]);

        newEdges[8]->set_first_halfedge(newQuad0[1]);
        newEdges[8]->set_second_halfedge(newQuad1[3]);

        newEdges[9]->set_first_halfedge(newQuad1[2]);
        newEdges[9]->set_second_halfedge(newQuad2[0]);

        newEdges[10]->set_first_halfedge(newQuad2[3]);
        newEdges[10]->set_second_halfedge(newQuad3[1]);

        newEdges[11]->set_first_halfedge(newQuad3[0]);
        newEdges[11]->set_second_halfedge(newQuad0[2]);

        // ************************************************************

        // ********* Half-edges to Edges/Faces pointers ****************

        newQuad0[0]->set_face(newFace0);
        newQuad0[0]->set_edge(newEdges[0]);

        newQuad1[1]->set_face(newFace1);
        newQuad1[1]->set_edge(newEdges[2]);

        newQuad2[2]->set_face(newFace2);
        newQuad2[2]->set_edge(newEdges[4]);

        newQuad3[3]->set_face(newFace3);
        newQuad3[3]->set_edge(newEdges[6]);

        // ************************************************************

        // ********* Faces pointers ****************

        newFace0->set_halfedge(newQuad0[0]);
        newFace1->set_halfedge(newQuad1[0]);
        newFace2->set_halfedge(newQuad2[0]);
        newFace3->set_halfedge(newQuad3[0]);

        // ************************************************************

        // ********* Marking the edges ****************

        for (i = 0; i < 8; i++)
          newEdges[i]->get_attributes().subdivide = true;

        // ************************************************************

        outsideH[0]->set_edge(newEdges[0]);

        if (outsideH[1])
            outsideH[1]->set_edge(newEdges[1]);

        outsideH[2]->set_edge(newEdges[2]);

        if (outsideH[3])
            outsideH[3]->set_edge(newEdges[3]);

        outsideH[4]->set_edge(newEdges[4]);

        if (outsideH[5])
            outsideH[5]->set_edge(newEdges[5]);

        outsideH[6]->set_edge(newEdges[6]);

        if (outsideH[7])
            outsideH[7]->set_edge(newEdges[7]);


        createdVertex->at(0)->set_halfedge( newQuad0[1] );
        createdVertex->at(1)->set_halfedge( newQuad1[2] );
        createdVertex->at(2)->set_halfedge( newQuad2[3] );
        createdVertex->at(3)->set_halfedge( newQuad3[0] );
        createdVertex->at(4)->set_halfedge( newQuad0[2] );
    }

    // Type 3 Template ( D - D twisted )
    template <
            typename VAttrib,
            typename FAttrib,
            typename EAttrib,
            typename HAttrib
            >
    void QOperators<VAttrib,FAttrib,EAttrib,HAttrib>::split_template3(QHalfedge *he, std::vector<QVertex*> *createdVertex) {

        QVertex *v[4];
        QEdge   *newEdges[5];
        QEdge   *oldEdges[4];
        QHalfedge *oldQuad[4];
        QHalfedge *newQuad1[4];
        QHalfedge *newQuad2[4];
        QHalfedge *outsideH[6];
        QHalfedge *heSearch;
        QFace *newFace1 = add_face();
        QFace *newFace2 = add_face();
        QFace *oldFace = he->get_face();
        int i = 0;

        // Setting all old Vertex to v[]

        for (; i < 4 ; i++){
            v[i] = he->get_origin();
            he   = he->get_next();
        }

        if (he->get_edge()->get_attributes().subdivide == true){

            heSearch = he->get_mate()->get_next();
            while (heSearch->get_next()->get_origin() != he->get_origin()) {
                heSearch = heSearch->get_mate();
                heSearch = heSearch->get_next();
            }

            oldEdges[0]  = heSearch->get_edge();
            newEdges[0]  = he->get_edge();
            outsideH[0]  = heSearch;
            outsideH[1]  = he->get_mate();

        } else {

            oldEdges[0] = he->get_edge();
            newEdges[0] = add_edge();
            outsideH[0] = he->get_mate();
            outsideH[1] = 0;
        }

        if (he->get_prev()->get_edge()->get_attributes().subdivide == true){

            // Searching for the e1[i] Edge
            heSearch = he->get_prev()->get_mate()->get_next();
            while (heSearch->get_next()->get_origin() != he->get_prev()->get_origin()) {
                heSearch = heSearch->get_mate();
                heSearch = heSearch->get_next();
            }

            oldEdges[3]  = he->get_prev()->get_edge();
            newEdges[4]  = heSearch->get_edge();
            outsideH[4]  = heSearch;
            outsideH[5]  = he->get_prev()->get_mate();

        } else {

            oldEdges[3]  = he->get_prev()->get_edge();
            newEdges[4]  = add_edge();
            outsideH[4]  = he->get_prev()->get_mate();
            outsideH[5]  = 0;
        }

        oldEdges[1] = he->get_next()->get_edge();
        oldEdges[2] = he->get_next()->get_next()->get_edge();

        newEdges[1] = add_edge();
        newEdges[2] = add_edge();
        newEdges[3] = add_edge();

        outsideH[2] = he->get_next()->get_mate();
        outsideH[3] = he->get_next()->get_next()->get_mate();

        oldQuad[0] = he;
        oldQuad[0]->set_edge(oldEdges[0]);
        oldQuad[0]->set_face(oldFace);
        oldQuad[0]->set_origin(v[0]);

        oldQuad[1] = he->get_next();
        oldQuad[1]->set_edge(oldEdges[1]);
        oldQuad[1]->set_face(oldFace);
        oldQuad[1]->set_origin(createdVertex->at(0));

        oldQuad[2] = he->get_next()->get_next();
        oldQuad[2]->set_edge(oldEdges[2]);
        oldQuad[2]->set_face(oldFace);
        oldQuad[2]->set_origin(createdVertex->at(2));

        oldQuad[3] = he->get_prev();
        oldQuad[3]->set_edge(oldEdges[3]);
        oldQuad[3]->set_face(oldFace);
        oldQuad[3]->set_origin(createdVertex->at(1));

        // newQuad1 (bottom quad)
        newQuad1[0] = new QHalfedge(createdVertex->at(0), newEdges[0], newFace1, 0, 0);
        newQuad1[1] = new QHalfedge(v[1], newEdges[1], newFace1, 0, newQuad1[0]);
        newQuad1[2] = new QHalfedge(v[2], newEdges[2], newFace1, 0, newQuad1[1]);
        newQuad1[3] = new QHalfedge(createdVertex->at(2), oldEdges[1], newFace1, newQuad1[0], newQuad1[2]);

        newQuad1[2]->set_next(newQuad1[3]);

        newQuad1[1]->set_next(newQuad1[2]);

        newQuad1[0]->set_next(newQuad1[1]);
        newQuad1[0]->set_prev(newQuad1[3]);

        // newQuad 2 (right quad)
        newQuad2[0] = new QHalfedge(createdVertex->at(1),oldEdges[2], newFace2, 0, 0);
        newQuad2[1] = new QHalfedge(createdVertex->at(2),newEdges[2], newFace2, 0, newQuad2[0]);
        newQuad2[2] = new QHalfedge(v[2], newEdges[3], newFace2, 0, newQuad2[1]);
        newQuad2[3] = new QHalfedge(v[3], newEdges[4], newFace2, newQuad2[0], newQuad2[2]);

        newQuad2[2]->set_next(newQuad2[3]);

        newQuad2[1]->set_next(newQuad2[2]);

        newQuad2[0]->set_next(newQuad2[1]);
        newQuad2[0]->set_prev(newQuad2[3]);

        // oldEdges
        oldEdges[0]->set_first_halfedge(oldQuad[0]);
        oldEdges[0]->set_second_halfedge(outsideH[0]);

        oldEdges[1]->set_first_halfedge(oldQuad[1]);
        oldEdges[1]->set_second_halfedge(newQuad1[3]);

        oldEdges[2]->set_first_halfedge(oldQuad[2]);
        oldEdges[2]->set_second_halfedge(newQuad2[0]);

        oldEdges[3]->set_first_halfedge(oldQuad[3]);
        oldEdges[3]->set_second_halfedge(outsideH[5]);

        // oldQuads
        oldQuad[0]->set_next(oldQuad[1]);
        oldQuad[0]->set_prev(oldQuad[3]);

        oldQuad[1]->set_next(oldQuad[2]);
        oldQuad[1]->set_prev(oldQuad[0]);

        oldQuad[2]->set_next(oldQuad[3]);
        oldQuad[2]->set_prev(oldQuad[1]);

        oldQuad[3]->set_next(oldQuad[0]);
        oldQuad[3]->set_prev(oldQuad[2]);

        oldQuad[0]->set_edge(oldEdges[0]);
        oldQuad[1]->set_edge(oldEdges[1]);
        oldQuad[2]->set_edge(oldEdges[2]);
        oldQuad[3]->set_edge(oldEdges[3]);

        oldFace->set_halfedge(oldQuad[0]);
        newFace1->set_halfedge(newQuad1[0]);
        newFace2->set_halfedge(newQuad2[0]);

        // PARAMETROS A CUIDAR
        oldEdges[0]->get_attributes().subdivide = true;
        oldEdges[3]->get_attributes().subdivide = true;
        newEdges[0]->get_attributes().subdivide = true;
        newEdges[4]->get_attributes().subdivide = true;
        // ***

        newEdges[0]->set_first_halfedge(newQuad1[0]);
        newEdges[0]->set_second_halfedge(outsideH[1]);

        newEdges[1]->set_first_halfedge(newQuad1[1]);
        newEdges[1]->set_second_halfedge(outsideH[2]);

        newEdges[2]->set_first_halfedge(newQuad1[2]);
        newEdges[2]->set_second_halfedge(newQuad2[1]);

        newEdges[3]->set_first_halfedge(newQuad2[2]);
        newEdges[3]->set_second_halfedge(outsideH[3]);

        newEdges[4]->set_first_halfedge(newQuad2[3]);
        newEdges[4]->set_second_halfedge(outsideH[4]);

        outsideH[0]->set_edge(oldEdges[0]);

        if (outsideH[1])
            outsideH[1]->set_edge(newEdges[0]);

        outsideH[2]->set_edge(newEdges[1]);
        outsideH[3]->set_edge(newEdges[3]);
        outsideH[4]->set_edge(newEdges[4]);

        if (outsideH[5])
            outsideH[5]->set_edge(oldEdges[3]);

        createdVertex->at(0)->set_halfedge( oldQuad[1] );
        createdVertex->at(1)->set_halfedge( newQuad2[0] );
        createdVertex->at(2)->set_halfedge( newQuad1[3] );

        newEdges[1]->get_attributes().geoline = oldEdges[1]->get_attributes().geoline;
        newEdges[3]->get_attributes().geoline = oldEdges[2]->get_attributes().geoline;

    }


    template <
            typename VAttrib,
            typename FAttrib,
            typename EAttrib,
            typename HAttrib
            >
    std::list<QVertex< VAttrib, FAttrib, EAttrib, HAttrib >*> QOperators<VAttrib,FAttrib,EAttrib,HAttrib>::split_strip(QHalfedge *h, int nPoints, double *points) {

        std::list<QEdge*> edges;
        std::list<QVertex*> v;

        for ( int i=0 ; i < nPoints ; i++ ){
            QVertex *vtx = add_vertex(points[3*i],points[3*i+1],points[3*i+2]);
            v.push_back( vtx );
        }

        for ( typename std::list<QVertex*>::iterator it = v.begin(); it != v.end(); ++it ){
            QEdge *e = add_edge();
            edges.push_back( e );
        }

        typename std::list<QEdge*>::iterator ite = edges.begin();
        typename std::list<QEdge*>::iterator iteAux = ite;
        iteAux++;

        typename std::list<QVertex*>::iterator itv = v.begin();
        typename std::list<QVertex*>::iterator itvAux = itv;
        itvAux++;

        for (; itv != v.end(); ++itvAux, ++itv, ++ite, ++iteAux ){

            if ( itvAux == v.end() )
               itvAux = v.begin();

            if ( iteAux == edges.end() )
                iteAux = edges.begin();

            split_quad( h, *itv, *itvAux, *ite, *iteAux );
            h = h->get_next()->get_next()->get_mate();

        }

        return v;

    }

}


#endif	/* QOPERATORS_H */

