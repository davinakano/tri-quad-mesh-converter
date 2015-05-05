#ifndef MESHREF_H
#define MESHREF_H

#include "types.h"
#include "macros.h"
#include "vtkReader.h"
#include "stdlib.h"
#include <sstream>

class MeshRef{

public:

    QMesh* mesh;
    geodesics* geo;
    adgFDCel* adg;

    MeshRef(QMesh* m, geodesics* g, adgFDCel* a){
        mesh = m;
        geo = g;
        adg = a;
    }

    ~MeshRef(){

    }

    geodesics::GP* findMiddlePoint( QFace* qf, geodesics::GP* v1, geodesics::GP* v2 ) {

        ASSERT( qf );
        ASSERT( v1 );
        ASSERT( v2 );

        std::list< geodesics::GP* > lv;

        geo->setRegion( qf );
        geo->compute(v1,v2,lv);

        ASSERT( lv.size() > 1 );

        geodesics::GP *g = NULL;

        if ( lv.size() == 2 ) {
            // pegar o meio
            // g = ?
        } else {

            std::list< geodesics::GP* >::iterator it = lv.begin();

            for ( unsigned i=0 ; i < lv.size()/2 ; i++)
                it++;

            g = *it;
        }

        return g;

    }


    void split(QHalfedge* Qh){

        QHalfedge* Qh2;
        Qh2 = Qh;

        std::list< geodesics::GP* > gpList;

        do {

            geodesics::GP* gpMiddle = findMiddlePoint( Qh2->get_face(), Qh2->get_origin()->get_attributes().getGP() ,
                                                       Qh2->get_next()->get_origin()->get_attributes().getGP() );

            ASSERT( gpMiddle );

            gpList.push_back(gpMiddle);
            Qh2 = Qh2->get_next()->get_next()->get_mate();

        } while (Qh != Qh2);

        double *points = new double[3*gpList.size()];
        int i = 0;

        std::list< geodesics::GP* >::iterator itgpList = gpList.begin();

        for ( ; itgpList != gpList.end() ; itgpList++ )
        {
            double x,y,z;
            geodesics::GP* gpoint = *itgpList;

            gpoint->get_coords(adg,x,y,z);

            points[3*i] = x;
            points[3*i+1] = y;
            points[3*i+2] = z;

            i++;
        }

        QOperators qop(mesh);

        std::list< QVertex* > qvList;

        qvList = qop.split_strip(Qh, gpList.size(), points);

        ASSERT( gpList.size() == qvList.size() );

        itgpList = gpList.begin();

        for ( std::list< QVertex* >::iterator itqvList = qvList.begin(); itqvList != qvList.end() ; itqvList++){

            QVertex* q = *itqvList;
            q->get_attributes().setGP(*itgpList);

            itgpList++;

        }

    }

};

#endif // MESHREF_H
