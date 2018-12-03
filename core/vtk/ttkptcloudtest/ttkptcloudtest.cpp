#include                  <ttkptcloudtest.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkptcloudtest)

ttkptcloudtest::ttkptcloudtest(){
    NumberGridPoints = 100;
    Offset = 0.1;
    GaussianKDE = false;
    Bandwidth = 0.01;
    Autobandwidth = false;
}

ttkptcloudtest::~ttkptcloudtest(){
}


// transmit abort signals -- to copy paste in other wrappers
bool ttkptcloudtest::needsToAbort(){
  return GetAbortExecute();
}


// transmit progress status -- to copy paste in other wrappers
int ttkptcloudtest::updateProgress(const float &progress){

  {
    stringstream msg;
    msg << "[ttkptcloudtest] " << progress*100 
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }
  
  UpdateProgress(progress);
  return 0;
}

int ttkptcloudtest::doIt(vtkDataSet *input, vtkUnstructuredGrid *output){

    int numpointsPointcloud = input->GetNumberOfPoints();
    std::vector<std::vector<double>> pointCloudCoordinates;

    //Detect the bounds for the grid. 
    for(int i = 0; i < numpointsPointcloud; i++){
        double x = input->GetPoint(i)[0];
        double y = input->GetPoint(i)[1];
        double z = input->GetPoint(i)[2];
        std::vector<double> xyz = {x,y,z};
        pointCloudCoordinates.push_back(xyz);
    }
    
    triangulation.setInputData(vtkImageData::New());
    triangulation.getTriangulation()->setWrapper(this);
    pointDistField.setGaussianKDE(GaussianKDE);
    pointDistField.setupTriangulation(triangulation.getTriangulation());
    pointDistField.setOffset(Offset);
    pointDistField.setAutoBandwidth(Autobandwidth);
    pointDistField.setBandwidth(Bandwidth);
    pointDistField.setNumberGridPoints(NumberGridPoints);
    pointDistField.setNumberCloudPoints(numpointsPointcloud);
    pointDistField.setInputPointCloud(static_cast<void*>(pointCloudCoordinates.data()));
    pointDistField.preprocess<double>();

    vtkDoubleArray* outputScalarField = vtkDoubleArray::New();
    outputScalarField->vtkDoubleArray::SetNumberOfComponents(1);
    outputScalarField->vtkDoubleArray::SetName("dist_field");

    if (pointDistField.isPlanar()){
        outputScalarField->vtkDoubleArray::SetNumberOfTuples(pow(NumberGridPoints,2));
    }
    else{
        outputScalarField->vtkDoubleArray::SetNumberOfTuples(pow(NumberGridPoints,3));
    }

    pointDistField.setOutputScalarFieldPointer(outputScalarField->GetVoidPointer(0));
    pointDistField.execute<double>();

    vtkSmartPointer<vtkUnstructuredGrid> triangle_unstruct_grid =
        triangulation.getVtkUnstructuredGrid();
    
    triangle_unstruct_grid->GetPointData()->AddArray(outputScalarField);
    output->ShallowCopy(triangle_unstruct_grid);

}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkptcloudtest::RequestData(vtkInformation *request, 
  vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;
  
  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
//  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), 1);
//  vtkImageData *output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::GetData(outputVector);
  
  doIt(input, output);
  
  {
    stringstream msg;
    msg << "[ttkptcloudtest] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  
  return 1;
}
