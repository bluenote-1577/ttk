/// \ingroup base
/// \class ttk::ScalarFieldFromPointCloud
/// \author Jim Shaw <jimshawster@gmail.com>
/// \date January 2019
///
/// \brief TTK processing package for converting point clouds into scalar fields. 
///
/// This package converts a set of points into a regular grid with a scalar field.
// The user can opt to obtain two types of scalar fields, a distance field or a kernel density estimation. 
///
/// \sa ttkScalarFieldFromPointCLoud.cpp %for a usage example.


#ifndef _SCALARFIELDFROMPOINTCLOUD_H
#define _SCALARFIELDFROMPOINTCLOUD_H

// base code includes
#include<Wrapper.h>
#include<Geometry.h>
#include<Triangulation.h>
#include<list>


namespace ttk{

    /**
     * Compute a scalar field an a regular grid from a point cloud. 
     */ 
  class ScalarFieldFromPointCloud: public Debug{

    public:
       ScalarFieldFromPointCloud();
      ~ScalarFieldFromPointCloud();

       int execute();

       //Must be called before execute(). 
       int preprocess();

       //If the point cloud is a 2-d (planar) point cloud. 
       inline bool isPlanar() const{
           return isPlanar_;
       }

       //The dimension n of the grid. So the grid has either n^3 points or n^2 points if planar. 
       inline int setNumberGridPoints(SimplexId numberGridPoints){
         numberGridPoints_ = numberGridPoints;
         return 0;
       }

       inline int setOutputScalarFieldPointer(void *data){
         outputScalarFieldPointer_=data;
         return 0;
       }

       inline int setNumberCloudPoints(SimplexId numberCloudPoints){
         numberCloudPoints_ = numberCloudPoints;

         return 0;
       }

       //The offset determines how closely the grid wraps around the point cloud.
       inline int setOffset(double offset){
         offset_ = offset;
         return 0;
       }

       inline int setGaussianKDE(bool gaussianKDE){
         gaussianKDE_ = gaussianKDE;
         return 0;
       }

       inline int setBandwidth(double bandwidth){
         bandwidth_ = bandwidth;
         return 0;
       }

       inline int setAutoBandwidth(bool autoBandwidth){
         autoBandwidth_ = autoBandwidth;
         return 0;
       }
       
       inline int setupTriangulation(Triangulation* triangulation){
         triangulation_=triangulation;
         return 0;
       }

       inline int setInputPointCloud(void *data){
         inputPointCloud_=data;
         return 0;
       }




    protected:
        Triangulation* triangulation_;
        void* inputPointCloud_;
        SimplexId numberCloudPoints_;
        SimplexId numberGridPoints_;
        double offset_;
        bool gaussianKDE_;
        double bandwidth_;
        bool autoBandwidth_;
        void* outputScalarFieldPointer_;
    private:
        double diagonalLength;
        double scaledBandwidth;
        bool donePreprocessing;
        SimplexId totalPointsGrid;
        double maxDistanceBound;
        double irrelevantExponent;
        double e;
        double pi;
        bool isPlanar_;
        std::vector<double> closestNeighbours
        (std::vector<double> *input) const;
    };
}

#endif // DISTANCEFIELD_H
