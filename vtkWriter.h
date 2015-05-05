#ifndef VTKWRITER_H
#define VTKWRITER_H


#include <string>
#include <vector>
#include <list>
#include <set>

#include "types.h"


class vtkWriter
{
public:
    static void saveVTKTri( std::string filename, TMesh* tm, QMesh* qm );
    static void saveVTKQuad( std::string filename, QMesh* _refinMesh );
    static void saveVTKTriPatch( std::string filename, unsigned nv, double* vset, unsigned nf, unsigned* fset, double* tvalue );
    static void saveVTKGrid( std::string filename, unsigned dimensionX, unsigned dimensionY, unsigned dimensionZ, unsigned npoints, double* points);
    static void saveVTKTriPatchWithoutTvalue( std::string filename, unsigned nv, double* vset, unsigned nf, unsigned* fset);

};

#endif // VTKWRITER_H
