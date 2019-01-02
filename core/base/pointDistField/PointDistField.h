#ifndef _POINTDISTFIELD_H
#define _POINTDISTFIELD_H

// base code includes
#include<Wrapper.h>
#include<Geometry.h>
#include<Triangulation.h>
#include<list>


namespace ttk{

  class PointDistField : public Debug{

    public:
       PointDistField();
      ~PointDistField();

       int execute();

       int preprocess();

      inline bool isPlanar() const{
          return isPlanar_;
      }

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
        double maxDist;
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
