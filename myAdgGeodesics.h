#ifndef MYADGGEODESICS_H
#define MYADGGEODESICS_H

#include "geodesics/adg_geodesics.hpp"
#include "geodesics/adg_from_dcel.hpp"

typedef adgfromdcel::adg_from_dcel adgFDCel;

class myApproxGeodesics : public dg::Approx_geodesics< TMesh > {

public:

    myApproxGeodesics( MI* m , double error , int niter ) : dg::Approx_geodesics< TMesh >(m, error, niter ) {

    }

    void setRegion( QFace* f ) {
        _region = f;
    }



private:



    QFace* _region;




    bool inRegion( Face* f, QFace* qf ) const;



    /**
     * \fn void get_restricted_adj_vertices( Vertex* v , std::vector< Vertex* >& a ) const
     *
     * \brief Finds all mesh vertices connected to a given mesh vertex
     * by an edge that is incident to a marked face.
     *
     * \param v A pointer to a geodesic vertex.
     * \param a A reference to a list of adjacent mesh vertices.
     */
    virtual void get_restricted_adj_vertices( Vertex* v , std::vector< Vertex* >& a ) const ;

};

#endif // MYADGGEODESICS_H
