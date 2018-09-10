#ifndef _TTK_PTLCOUDTEST_H
#define _TTK_PTCLOUDTEST_H

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
#include                  <vtkSphereSource.h>
#include                  <vtkType.h>
#include                  <vtkInformationVector.h>

// ttk code includes
#include                  <Wrapper.h>
#include                  <cmath>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkptcloudtest
#else
class ttkptcloudtest
#endif 
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:
      
    static ttkptcloudtest* New();
   
    // macros
    vtkTypeMacro(ttkptcloudtest, vtkDataSetAlgorithm);

    // default ttk setters
    vtkSetMacro(debugLevel_, int);

    int FillOutputPortInformation(int port,
      vtkInformation *info){
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
      info->Set(vtkDataObject::FIELD_NUMBER_OF_COMPONENTS(), 1);
      info->Set(vtkDataObject::FIELD_ATTRIBUTE_TYPE(), VTK_DOUBLE);
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
    
  protected:
    
    ttkptcloudtest();
    
    ~ttkptcloudtest();
    
    int RequestData(vtkInformation *request, 
      vtkInformationVector **inputVector, vtkInformationVector *outputVector);
    
    
  private:
    
    bool                  UseAllCores;
    int                   ThreadNumber;
    // base code features
    int doIt(vtkDataSet *input, vtkImageData *output);
    
    bool needsToAbort();
    
    int updateProgress(const float &progress);
   
};

#endif // _TTK_PTCLOUDTEST_H
