#ifndef VTKREADER_H
#define VTKREADER_H

#include <string>
#include <vector>
#include <list>
#include <set>

#include "types.h"

class vtkReader
{
public:
    vtkReader();
    ~vtkReader();

    TMesh* getMeshTri();
    QMesh* getMeshQuad();

    bool loadMesh( std::string filename );

private:

    TMesh* _MeshTri;
    QMesh* _MeshQuad;
};

#endif // VTKREADER_H
