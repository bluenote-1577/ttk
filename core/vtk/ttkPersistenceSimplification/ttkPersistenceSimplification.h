/// \ingroup vtk
/// \class ttkPersistenceSimplification
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the persistenceSimplification processing package.
///
/// VTK wrapping code for the @PersistenceSimplification package.
/// 
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
///
/// \sa ttk::PersistenceSimplification
#pragma once

// VTK includes -- to adapt
#include                  <vtkCharArray.h>
#include                  <vtkDataArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkIntArray.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPointData.h>
#include                  <vtkSmartPointer.h>

// ttk code includes
#include                  <PersistenceSimplification.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a 
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK 
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkPersistenceSimplification
#else
class ttkPersistenceSimplification
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:
    
    static ttkPersistenceSimplification* New();
    vtkTypeMacro(ttkPersistenceSimplification, vtkDataSetAlgorithm)
    
    // default ttk setters
    vtkSetMacro(debugLevel_, int);
    
    void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }
    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters

    // set-getters macros to define from each variable you want to access from 
    // the outside (in particular from paraview) - to adapt.
    // Note that the XML file for the ParaView plug-in specification needs to be
    // edited accordingly.
    
    vtkSetMacro(DistinctOption, bool);
    vtkGetMacro(DistinctOption, bool);
    
    vtkSetMacro(AutoOption, bool);
    vtkGetMacro(AutoOption, bool);
    
    vtkSetMacro(UseMinOption, bool);
    vtkGetMacro(UseMinOption, bool);
    
    vtkSetMacro(UseMaxOption, bool);
    vtkGetMacro(UseMaxOption, bool);
    
    vtkSetMacro(CountMinOption, bool);
    vtkGetMacro(CountMinOption, bool);
    
    vtkSetMacro(CountMaxOption, bool);
    vtkGetMacro(CountMaxOption, bool);
    
    vtkSetMacro(CountAllOption, bool);
    vtkGetMacro(CountAllOption, bool);

    vtkSetMacro(CountMinArg, int);
    vtkGetMacro(CountMinArg, int);

    vtkSetMacro(CountMaxArg, int);
    vtkGetMacro(CountMaxArg, int);

    vtkSetMacro(CountAllArg, int);
    vtkGetMacro(CountAllArg, int);
   
    vtkSetMacro(PersistenceMaxArg, double);
    vtkGetMacro(PersistenceMaxArg, double);
   
    vtkSetMacro(PersistenceMinArg, double);
    vtkGetMacro(PersistenceMinArg, double);
   
    vtkSetMacro(PersistenceAllArg, double);
    vtkGetMacro(PersistenceAllArg, double);

    
    vtkSetMacro(ScalarField, std::string);
    vtkGetMacro(ScalarField, std::string);
    
  protected:
   
    ttkPersistenceSimplification(){

      CountMinArg = 2;
      CountMaxArg = 0;
      CountAllArg = 0;
      UseMaxOption = true;
      UseMinOption = false;
      PersistenceAllArg = 0.0;
      PersistenceMinArg = 0.0;
      PersistenceMaxArg = 0.0;
      DistinctOption = true;
      CountMinOption = true;
      CountMaxOption = true;
      CountAllOption = true;
      AutoOption = true;
      outputScalarField_ = NULL;
      
      UseAllCores = true;
      ThreadNumber = 1;
      debugLevel_ = 3;

    }
    
    ~ttkPersistenceSimplification(){};
    
    TTK_SETUP();
    
    
  private:
    // Thresholding options
    int                   CountMinArg, CountMaxArg, CountAllArg;
    double                PersistenceMaxArg, PersistenceMinArg, PersistenceAllArg;
    bool                  DistinctOption, AutoOption, CountMaxOption, CountMinOption, CountAllOption;
    bool                  UseMinOption, UseMaxOption;

    // default variables
    std::string           ScalarField;
    vtkDataArray          *outputScalarField_, *inputOffsets_;
    ttk::PersistenceSimplification            persistenceSimplification_;
};
