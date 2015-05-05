/**
 * \file operators.h
 *
 * \brief Definition  of the class Operators, which performs some topological
 * operations in a Surface mesh.
 *
 * \author
 * Mario Augusto de Souza Lizier \n
 * Icaro Lins Leitao da Cunha \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
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

#ifndef OPERATORS_H
#define	OPERATORS_H

#include "vertex.h"    // Vertex
#include "halfedge.h"  // Halfedge
#include "edge.h"      // Edge
#include "face.h"      // Face
#include "surface.h"   // Surface

/**
 * \defgroup DCELNameSpace Namespace adcel.
 * @{
 */

/**
 * \namespace dcel
 *
 * \brief   The   namespace   adcel   contains   the   definition   and
 * implementation of all classes of the data structure Doubly
 * Connected  Edge List  (DCEL), which  can be  used to  represent the
 * mesh.
 *
 */

namespace dcel {

    /**
     * \class Operators
     *
     * \brief This class performs some topological operations in a Surface mesh.
     */
    template <
            typename VAttrib = int,
            typename FAttrib = int,
            typename EAttrib = int,
            typename HAttrib = int
            >
            class Operators {
    public:

        /**
         * \typedef Vertex
         *
         * \brief Defines  Vertex as  an alias for  dcel::Vertex< VAttrib,
         * FAttrib , EAttrib , HAttrib >.
         */
        typedef dcel::Vertex< VAttrib, FAttrib, EAttrib, HAttrib > Vertex;

        /**
         * \typedef Face
         *
         * \brief  Defines  Face  as  an alias  for  dcel::Face<  VAttrib,
         * FAttrib , EAttrib , HAttrib >.
         */
        typedef dcel::Face< VAttrib, FAttrib, EAttrib, HAttrib > Face;

        /**
         * \typedef Edge
         *
         * \brief  Defines  Edge  as  an alias  for  dcel::Edge<  VAttrib,
         * FAttrib , EAttrib , HAttrib >.
         */
        typedef dcel::Edge< VAttrib, FAttrib, EAttrib, HAttrib > Edge;

        /**
         * \typedef Halfedge
         *
         * \brief  Defines  Halfedge   as  an  alias  for  dcel::Halfedge<
         * VAttrib, FAttrib , EAttrib , HAttrib >.
         */
        typedef dcel::Halfedge< VAttrib, FAttrib, EAttrib, HAttrib > Halfedge;


        /**
         * \typedef Surface
         *
         * \brief Defines Surface as an alias for dcel::Surface < VAttrib,
         *  FAttrib , EAttrib , HAttrib >.
         */
        typedef dcel::Surface< VAttrib, FAttrib, EAttrib, HAttrib > Surface;



        Operators(Surface *_mesh);

        virtual ~Operators();

        bool split_triangle(Face *f, double x, double y, double z, Face *&f0, Face *&f1, Face *&f2);

        void split_quad(Halfedge *h, double x, double y, double z, Face *&f0, Face *&f1, Face *&f2, Face *&f3);

    private:

        void split_edge(Halfedge *h, double x, double y, double z);

        Vertex* add_vertex(double x, double y, double z, Halfedge *h = 0) {

            Vertex* vertex = new Vertex(x, y, z, h);

            mesh->add_vertex(vertex);

            return vertex;
        }

        Face* add_face(Halfedge *h = 0) {

            Face* face = new Face(h);

            mesh->add_face(face);

            return face;
        }

        Edge* add_edge(Halfedge *h1 = 0, Halfedge *h2 = 0) {
            
            Edge* edge = new Edge(h1, h2);

            mesh->add_edge(edge);

            return edge;
        }


        Surface *mesh;
    };


    template <
            typename VAttrib,
            typename FAttrib,
            typename EAttrib,
            typename HAttrib
            >
    Operators<VAttrib,FAttrib,EAttrib,HAttrib>::Operators( Surface *_mesh) : mesh(_mesh)
    {
    }

    template <
            typename VAttrib,
            typename FAttrib,
            typename EAttrib,
            typename HAttrib
            >
    Operators<VAttrib,FAttrib,EAttrib,HAttrib>::~Operators()
    {

    }

    template <
            typename VAttrib,
            typename FAttrib,
            typename EAttrib,
            typename HAttrib
            >
    bool Operators<VAttrib,FAttrib,EAttrib,HAttrib>::split_triangle(Face *f, double x, double y, double z, Face *&f0, Face *&f1, Face *&f2) {

        Halfedge *hv1 = f->get_halfedge()->get_next()->get_mate();
        Halfedge *hv2 = f->get_halfedge()->get_prev()->get_mate();
        Vertex *v0 = f->get_halfedge()->get_origin();
        Vertex *v1 = f->get_halfedge()->get_next()->get_origin();
        Vertex *v2 = f->get_halfedge()->get_prev()->get_origin();

        Vertex *v = add_vertex(x, y, z);

        f0 = f;

        f1 = add_face();
        Halfedge * h1[3];
        h1[0] = new Halfedge(v, 0, f1, 0, 0);
        h1[1] = new Halfedge(v1, 0, f1, 0, h1[0]);
        h1[2] = new Halfedge(v2, 0, f1, h1[0], h1[1]);
        h1[0]->set_next(h1[1]);
        h1[0]->set_prev(h1[2]);
        h1[1]->set_next(h1[2]);
        f1->set_halfedge(h1[0]);


        f2 = add_face();
        Halfedge * h2[3];
        h2[0] = new Halfedge(v, 0, f2, 0, 0);
        h2[1] = new Halfedge(v2, 0, f2, 0, h2[0]);
        h2[2] = new Halfedge(v0, 0, f2, h2[0], h2[1]);
        h2[0]->set_next(h2[1]);
        h2[0]->set_prev(h2[2]);
        h2[1]->set_next(h2[2]);
        f2->set_halfedge(h2[0]);

        Edge *e[3];

        e[0] = add_edge(h1[0], f->get_halfedge()->get_next() );
        e[1] = add_edge(h2[0], h1[2] );
        e[2] = add_edge(h2[2], f->get_halfedge()->get_prev() );

        if( hv1->get_edge()->get_first_halfedge() == hv1 )
            hv1->get_edge()->set_second_halfedge( h1[1] );
        else
            hv1->get_edge()->set_first_halfedge( h1[1] );

        if( hv2->get_edge()->get_first_halfedge() == hv2 )
            hv2->get_edge()->set_second_halfedge( h2[1] );
        else
            hv2->get_edge()->set_first_halfedge( h2[1] );

        h1[0]->set_edge( e[0] );
        h1[1]->set_edge( hv1->get_edge() );
        h1[2]->set_edge( e[1] );

        h2[0]->set_edge( e[1] );
        h2[1]->set_edge( hv2->get_edge() );
        h2[2]->set_edge( e[2] );

        f->get_halfedge()->get_next()->set_edge( e[0] );
        f->get_halfedge()->get_prev()->set_edge( e[2] );

        f->get_halfedge()->get_prev()->set_origin( v );

        v->set_halfedge(h1[0]);
        v0->set_halfedge(h2[2]);
        v1->set_halfedge(h1[1]);
        v2->set_halfedge(h2[1]);

        return true;
    }

    template <
            typename VAttrib,
            typename FAttrib,
            typename EAttrib,
            typename HAttrib
            >
    void Operators<VAttrib,FAttrib,EAttrib,HAttrib>::split_edge( Halfedge* h , double x , double y , double z ) {
        Vertex* vnew = add_vertex( x , y , z ) ;

        Face* f1 = h->get_face() ;
        Face* f2 = h->get_mate()->get_face() ;

        Halfedge* h_next = h->get_next() ;
        Halfedge* h_mate = h->get_mate() ;
        Halfedge* h_mate_next = h_mate->get_next() ;

        Halfedge* h1new = new Halfedge( vnew , 0 , f1 , h_next , h ) ;
        Halfedge* h2new = new Halfedge( vnew , 0 , f2 , h_mate_next , h_mate ) ;

        vnew->set_halfedge( h1new ) ;

        h->set_next( h1new ) ;
        h_next->set_prev( h1new ) ;

        h_mate->set_next( h2new ) ;
        h_mate_next->set_prev( h2new ) ;

        Edge* enew = add_edge( h1new , h_mate) ;

        h1new->set_edge( enew ) ;
        h_mate->set_edge( enew ) ;

        h->get_edge()->set_first_halfedge( h ) ;
        h->get_edge()->set_second_halfedge( h2new ) ;

        h2new->set_edge( h->get_edge() ) ;

        return ;
    }

    template <
            typename VAttrib,
            typename FAttrib,
            typename EAttrib,
            typename HAttrib
            >
    void Operators<VAttrib,FAttrib,EAttrib,HAttrib>::split_quad(Halfedge *h, double x, double y, double z, Face *&f0, Face *&f1, Face *&f2, Face *&f3) {
        split_edge( h , x , y , z ) ;

        Halfedge* h1new = new Halfedge( h->get_next()->get_origin() , 0 , 0 , h->get_prev() , h ) ;
        Halfedge* h2new = new Halfedge( h->get_prev()->get_origin() , 0 , 0 , h->get_next() , h->get_next()->get_next() ) ;

        Edge* e1new = add_edge( h1new , h2new ) ;

        h1new->set_edge( e1new ) ;
        h2new->set_edge( e1new ) ;

        Halfedge* h_mate = h->get_mate() ;

        Halfedge* h3new = new Halfedge( h_mate->get_next()->get_next()->get_origin() , 0 , 0 , h_mate , h_mate->get_next() ) ;
        Halfedge* h4new = new Halfedge( h_mate->get_origin() , 0 , 0 , h_mate->get_prev()->get_prev() , h_mate->get_prev() ) ;

        Edge* e2new = add_edge( h3new , h4new ) ;

        h3new->set_edge( e2new ) ;
        h4new->set_edge( e2new ) ;

        Face* f1new = add_face( h2new  ) ;
        Face* f2new = add_face( h4new ) ;

        h->get_face()->set_halfedge( h ) ;
        h1new->set_face( h->get_face() ) ;
        h->set_next( h1new ) ;
        h1new->get_next()->set_prev( h1new ) ;

        h_mate->get_face()->set_halfedge( h_mate ) ;
        h3new->set_face( h_mate->get_face() ) ;
        h_mate->set_prev( h3new ) ;
        h3new->get_prev()->set_next( h3new ) ;

        h2new->get_next()->set_prev( h2new ) ;
        h2new->get_prev()->set_next( h2new ) ;
        h2new->set_face( f1new ) ;
        h2new->get_prev()->set_face( f1new ) ;
        h2new->get_next()->set_face( f1new ) ;

        h4new->get_next()->set_prev( h4new ) ;
        h4new->get_prev()->set_next( h4new ) ;
        h4new->set_face( f2new ) ;
        h4new->get_prev()->set_face( f2new ) ;
        h4new->get_next()->set_face( f2new ) ;

        f0 = h2new->get_face();
        f1 = h1new->get_face();
        f2 = h3new->get_face();
        f3 = h4new->get_face();

        return ;
    }


}


#endif	/* OPERATORS_H */

