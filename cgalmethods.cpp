#include "cgalmethods.h"
#include "macros.h"


#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/parameterize.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Cartesian<double>             Kernel;
typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;
typedef CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> PolyhedronBuilder;
typedef CGAL::Triangulation_2<Kernel>         Triangulation;

template <class HDS>
class buildMesh : public CGAL::Modifier_base<HDS> {
public:
    buildMesh(std::set< TMesh::Vertex* >& _vertices, std::set< TMesh::Face* >& _faces)
        : vertices(_vertices), faces(_faces) { }
    void operator()( HDS& hds)  {
        CGAL::Polyhedron_incremental_builder_3<HDS> builder( hds, true);

        int nv = vertices.size();
        int nf = faces.size();

        builder.begin_surface(nv, nf );
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;

        int i = 0;
        for(std::set< TMesh::Vertex* >::iterator itv = vertices.begin(); itv!=vertices.end(); ++itv) {
            TMesh::Vertex* v = *itv;
            v->get_attributes().tempInt = i;
            builder.add_vertex( Point(v->x(), v->y(), v->z()) );
            i++;
        }

        for(std::set< TMesh::Face* >::iterator itf = faces.begin(); itf!=faces.end(); ++itf) {
            TMesh::Face* f = *itf;

            builder.begin_facet();
            builder.add_vertex_to_facet( f->get_halfedge()->get_origin()->get_attributes().tempInt );
            builder.add_vertex_to_facet( f->get_halfedge()->get_next()->get_origin()->get_attributes().tempInt );
            builder.add_vertex_to_facet( f->get_halfedge()->get_prev()->get_origin()->get_attributes().tempInt );
            builder.end_facet();
        }
        builder.end_surface();
    }

    std::set< TMesh::Vertex* >& vertices;
    std::set< TMesh::Face* >& faces;
};


void cgalMethods::parametrization(std::set< Vertex* >& vertices, std::set< Face* >& faces, std::set< pack3<Vertex*,double,double> >& vertices_uv) {

    Polyhedron mesh;
    buildMesh<HalfedgeDS> builder(vertices,faces);
    mesh.delegate(builder);

    //***************************************
    // Create Polyhedron adaptor
    // Note: no cutting => we support only
    // meshes that are topological disks
    //***************************************

    typedef CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron>
            Parameterization_polyhedron_adaptor;
    Parameterization_polyhedron_adaptor mesh_adaptor(mesh);

    //***************************************
    // Floater Mean Value Coordinates parameterization
    // (defaults are circular border and OpenNL solver)
    //***************************************

    typedef CGAL::Parameterizer_traits_3<Parameterization_polyhedron_adaptor>
            Parameterizer;  // Type that defines the error codes

    Parameterizer::Error_code err = CGAL::parameterize(mesh_adaptor);
    switch(err) {
    case Parameterizer::OK: // Success
        break;
    case Parameterizer::ERROR_EMPTY_MESH: // Input mesh not supported
    case Parameterizer::ERROR_NON_TRIANGULAR_MESH:
    case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:
    case Parameterizer::ERROR_BORDER_TOO_SHORT:
        std::cerr << "Input mesh not supported: " << Parameterizer::get_error_message(err) << std::endl;
        return;
        break;
    default: // Error
        std::cerr << "Error: " << Parameterizer::get_error_message(err) << std::endl;
        return;
        break;
    };

    //***************************************
    // Output
    //***************************************

    std::set< Vertex* >::iterator itv = vertices.begin();

    // Raw output: dump (u,v) pairs
    Polyhedron::Vertex_const_iterator pVertex;
    for (pVertex = mesh.vertices_begin();
         pVertex != mesh.vertices_end();
         pVertex++)
    {
        // (u,v) pair is stored in any halfedge
        double u = mesh_adaptor.info(pVertex->halfedge())->uv().x();
        double v = mesh_adaptor.info(pVertex->halfedge())->uv().y();
        vertices_uv.insert( pack3<Vertex*,double,double>( *itv, u, v));
        itv++;
    }


}

void cgalMethods::intersection( double s0u, double s0v, double s1u, double s1v, double r0u, double r0v, double r1u, double r1v, int &tipo, double& value) {

    CGAL::Point_2<K> s0(s0u, s0v);
    CGAL::Point_2<K> s1(s1u, s1v);
    CGAL::Point_2<K> r0(r0u, r0v);
    CGAL::Point_2<K> r1(r1u, r1v);

    CGAL::Segment_2<K> s( s0, s1 );
    CGAL::Segment_2<K> r( r0, r1 );

    CGAL::Object result = CGAL::intersection(s, r);
    if (const CGAL::Point_2<K> *ipoint = CGAL::object_cast<CGAL::Point_2<K> >(&result)) {

        // handle the point intersection case with *ipoint.
        if (*ipoint == s0 )
            tipo = 0;
        else if (*ipoint == s1 )
            tipo = 1;
        else {
            tipo = 2;

            double total = (s1.x()-s0.x())*(s1.x()-s0.x()) + (s1.y()-s0.y())*(s1.y()-s0.y());
            double part1 = (ipoint->x()-s0.x())*(ipoint->x()-s0.x()) + (ipoint->y()-s0.y())*(ipoint->y()-s0.y());
            value = sqrt( part1/total );

            //double vx = fabs( (value*s0.x() + (1.0 - value)*s1.x()) - ipoint->x());
            //double vy = fabs( (value*s0.y() + (1.0 - value)*s1.y()) - ipoint->y());
            //ASSERT( vx < 1e-10 );
            //ASSERT( vy < 1e-10 );

            if( value >= 1.0 ) {
                tipo = 1;
            } else if ( value <= 0.0 ) {
                tipo = 0;
            }

            // ASSERT( (value>0)&&(value<1.0) );
        }

    } else
        if (const CGAL::Segment_2<K> *iseg = CGAL::object_cast<CGAL::Segment_2<K> >(&result)) {

            // handle the segment intersection case with *iseg.
            tipo = 3;

        } else {

            // handle the no intersection case.
            tipo = 4;
        }

}
