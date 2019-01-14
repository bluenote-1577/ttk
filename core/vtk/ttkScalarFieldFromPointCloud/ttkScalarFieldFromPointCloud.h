/// \ingroup vtk
/// \class ttkScalarFieldFromPointCloud
/// \author Jim Shaw <jimshawster@gmail.com>
/// \date January 2019
///
/// \brief TTK VTK-filter that wraps the scalarFieldFromFromPointCloud package.
///
/// VTK wrapping code for the @ScalarFieldFromPointCloud package.
/// 
/// \param Input Point Cloud (vtkDataSet)
/// \param Output Output scalar field and regular grid (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
///


#ifndef _TTK_SCALARFIELDFROMPOINTCLOUD_H
#define _TTK_SCALARFIELDFROMPOINTCLOUD_H

// VTK includes 
#include                  <vtkAppendPolyData.h>  
#include                  <vtkCharArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkIntArray.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkStructuredGrid.h>
#include                  <vtkImageData.h>
#include                  <vtkPointData.h>
#include                  <vtkSmartPointer.h>
#include                  <vtkType.h>
#include                  <vtkInformationVector.h>


// ttk code includes
#include                  <ttkWrapper.h>
#include                  <ScalarFieldFromPointCloud.h>
#include                  <Triangulation.h>
#include                  <cmath>
#include                  <vector>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkScalarFieldFromPointCloud
#else
class ttkScalarFieldFromPointCloud
#endif 
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:
      
    static ttkScalarFieldFromPointCloud* New();
   
    // macros
    vtkTypeMacro(ttkScalarFieldFromPointCloud, vtkDataSetAlgorithm);

    // default ttk setters
    vtkSetMacro(debugLevel_, int);
    vtkSetMacro(NumberGridPoints, int);
    vtkSetMacro(Offset, double);
    vtkSetMacro(Bandwidth, double);
    vtkSetMacro(GaussianKDE, bool);
    vtkSetMacro(Autobandwidth, bool);

    int FillOutputPortInformation(int port,
      vtkInformation *info){
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      return 1;
    }

    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else{
        threadNumber_ = ttk::OsCall::getNumberOfCores();
      }
      Modified();
    }
    
    void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }   
    
    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters
  
    int RequestData(vtkInformation *request,
            vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;


    
  protected:
    
    ttkScalarFieldFromPointCloud();
    
    ~ttkScalarFieldFromPointCloud();
    
        
    
  private:
    
    bool                  UseAllCores;
    bool                  GaussianKDE;
    bool                  Autobandwidth;
    int                   ThreadNumber;
    int                   NumberGridPoints;
    double                 Offset;
    double                 Bandwidth;
    double                 Mean;
    ttkTriangulation      triangulation;
    ttk::ScalarFieldFromPointCloud   scalarFieldFromPointCloud;

    // base code features
    int doIt(vtkDataSet *input, vtkUnstructuredGrid *output);
    int doIt_test(vtkDataSet *input, vtkUnstructuredGrid *output);

    std::vector<double> closestNeighbours(vtkDataSet* input, bool isPlanar);
    
    bool needsToAbort();
    
    int updateProgress(const float &progress);
   
};

#endif 
