#ifndef CGALMETHODS_H
#define CGALMETHODS_H

#include "types.h"
#include <set>
#include "pack.h"

class cgalMethods
{
public:

    static void parametrization(std::set<Vertex *> &vertices, std::set<Face *> &faces, std::set<pack3<Vertex*,double,double> > &vertices_uv);
    static void intersection( double s0u, double s0v, double s1u, double s1v, double r0u, double r0v, double r1u, double r1v, int &tipo, double& value);
};

#endif // CGALMETHODS_H
