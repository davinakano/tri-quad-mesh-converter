#ifndef GEODESICMESHING_H
#define GEODESICMESHING_H

#include <list>
#include <vector>
#include <set>

#include "types.h"
#include "myAdgGeodesics.h"

#include "meshref.h"
//#include "geodesics.h"


class geodesicMeshing {

private:
/*
    typedef geodesics::GP< TMesh > GP ;
    typedef geodesics::GPV< TMesh > GPV ;
    typedef geodesics::GPE< TMesh > GPE ;
    typedef geodesics::GPF< TMesh > GPF ;
    typedef geodesics::Geodesics< TMesh > Geodesics ;
*/

public:

    geodesicMeshing( TMesh* refinMesh );
    ~geodesicMeshing();

    bool getGeodesic( QFace* region, Vertex* vi, Vertex* ve, std::list< Vertex* >& vertices );

    bool remeshingGeodesics( std::list< geodesics::GP* >& pts, std::list< Vertex* >& vertices );

private:

    bool split2( geodesics::GPV* gpv_prev, geodesics::GPE* gpe_next, Halfedge*& he_left, Halfedge*& he_right);
    bool split2( geodesics::GPE* gpe_prev, geodesics::GPV* gpv_next, Halfedge*& he_left, Halfedge*& he_right);
    bool split3( geodesics::GPE* gpe_prev, geodesics::GPE* gpe_next, Halfedge*& he_left, Halfedge*& he_right);
    bool deleteGPs();

    std::set< geodesics::GP* > _todelete;

    geodesics* _geo;

    adgFDCel* _adgMesh;
    TMesh* _refinMesh;
};

#endif // OFGEODESICS_H
