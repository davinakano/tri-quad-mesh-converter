#ifndef TYPES_H
#define TYPES_H

#include <set>

#include "qdcel/qsurface.h"
#include "qdcel/qoperators.h"

#include "dcel/surface.h"
#include "dcel/operators.h"

class myVertexAttributes;
class myFaceAttributes;
class myEdgeAttributes;
class myHalfedgeAttributes;

class myQVertexAttributes;
class myQFaceAttributes;
class myQEdgeAttributes;
class myQHalfedgeAttributes;

typedef qdcel::QSurface<myQVertexAttributes,myQFaceAttributes,myQEdgeAttributes,myQHalfedgeAttributes> QMesh;
typedef qdcel::QOperators<myQVertexAttributes,myQFaceAttributes,myQEdgeAttributes,myQHalfedgeAttributes> QOperators;

typedef QMesh::QHalfedge QHalfedge;
typedef QMesh::QEdge QEdge;
typedef QMesh::QVertex QVertex;
typedef QMesh::QFace QFace;

typedef dcel::Surface<myVertexAttributes,myFaceAttributes,myEdgeAttributes,myHalfedgeAttributes> TMesh;
typedef dcel::Operators<myVertexAttributes,myFaceAttributes,myEdgeAttributes,myHalfedgeAttributes> TOperators;

typedef TMesh::Halfedge Halfedge;
typedef TMesh::Edge Edge;
typedef TMesh::Vertex Vertex;
typedef TMesh::Face Face;

class myVertexAttributes {
public:
    int id;
    bool status, cl;
    double distance;
    int tempInt;
    double tempDouble;

    int getId(){
        return id;
    }

    void setId(int newId){
        this->id = newId;
    }

    // Needed by adg_from_dcel.hpp

    bool alive(){
        return this->status;
    }

    void set_alive(bool a){
        this->status = a;
    }

    bool close(){
        return this->cl;
    }

    void set_close(bool a){
        this->cl = a;
    }

    double dist(){
        return this->distance;
    }

    void set_dist(double d){
        this->distance = d;
    }



};

class myFaceAttributes {
public:
    int id;
    bool status;

    std::set< QFace* > region;

    void activate(){
        this->status = true;
    }

    void deactivate(){
        this->status = false;
    }

    bool is_active(){
        return this->status;
    }

};

class myEdgeAttributes {
public:
    int id;
};

class myHalfedgeAttributes {
public:
    int id;
};

#include "myAdgGeodesics.h"

typedef myApproxGeodesics geodesics;

// Attributes of the QVertex

class myQVertexAttributes {
private:
    int id;
    geodesics::GP* gPoint;

public:

    myQVertexAttributes() {
        id = -1;
        gPoint = 0;
    }

    int getId(){
        return id;
    }

    void setId(int newId){
        this->id = newId;
    }

    // Pointer to the correspondent vertex of the Triangle Mesh

    void setGP(geodesics::GP* v){
        gPoint = v;
    }

    geodesics::GP* getGP(){
        return gPoint;
    }

};

class myQFaceAttributes {

public:

    enum OrientationCut {

         // Taking the qf->get_halfedge() as reference

         TEMPLATE_NNNN = 0000, // Template 0 - No subdivisions
         TEMPLATE_NNND = 0001, // Template 1 - Invalid
         TEMPLATE_NNDN = 0010, // Tempalte 1 - Invalid
         TEMPLATE_NNDD = 0011, // Template 3 - Top Right Subdivisions
         TEMPLATE_NDNN = 0100, // Template 1 - Invalid
         TEMPLATE_NDND = 0101, // Template 2 - Top Down Subdivisions
         TEMPLATE_NDDN = 0110, // Template 3 - Bottom Right Subdivisions
         TEMPLATE_NDDD = 0111, // Template 5 - Invalid
         TEMPLATE_DNNN = 1000, // Template 1 - Invalid
         TEMPLATE_DNND = 1001, // Template 3 - Left Top Subdivisions
         TEMPLATE_DNDN = 1010, // Template 2 - Left Right Subdivisions
         TEMPLATE_DNDD = 1011, // Template 5 - Invalid
         TEMPLATE_DDNN = 1100, // Template 3 - Left Bottom Subdivisions
         TEMPLATE_DDND = 1101, // Template 5 - Invalid
         TEMPLATE_DDDN = 1110, // Template 5 - Invalid
         TEMPLATE_DDDD = 1111  // Template 4 - Cross Subdivisions
    };

    int id;
    bool subdivide, border;
    OrientationCut orientation;

    std::set< Face* > faces;

    myQFaceAttributes(){
        subdivide = false;
        border = false;
        orientation = TEMPLATE_NNNN;
    }

    int getTemplateType(){

        switch(orientation){

            case TEMPLATE_NNNN: return 0;
            case TEMPLATE_NNND: return 1;
            case TEMPLATE_NNDN: return 1;
            case TEMPLATE_NNDD: return 3;
            case TEMPLATE_NDNN: return 1;
            case TEMPLATE_NDND: return 2;
            case TEMPLATE_NDDN: return 3;
            case TEMPLATE_NDDD: return 5;
            case TEMPLATE_DNNN: return 1;
            case TEMPLATE_DNND: return 3;
            case TEMPLATE_DNDN: return 2;
            case TEMPLATE_DNDD: return 5;
            case TEMPLATE_DDNN: return 3;
            case TEMPLATE_DDND: return 5;
            case TEMPLATE_DDDN: return 5;
            case TEMPLATE_DDDD: return 4;


        }
    }
};

class myQEdgeAttributes {
public:
    int id;
    bool subdivide;
    bool checkEdge;

    std::list< Vertex* > geoline;

    myQEdgeAttributes(){
        subdivide = false;
    }
};

class myQHalfedgeAttributes {
public:
    int id;
};

#endif
