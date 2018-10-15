#ifndef _POINTDISTFIELD_H
#define _POINTDISTFIELD_H

// base code includes
#include<Wrapper.h>
#include<Geometry.h>
#include<Triangulation.h>


namespace ttk{

  class PointDistField : public Debug{

    public:
       PointDistField();
      ~PointDistField();
       int execute();

    protected:
    };
}


#endif // DISTANCEFIELD_H
