#include                  <ttkScalarFieldFromPointCloud.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkScalarFieldFromPointCloud)

ttkScalarFieldFromPointCloud::ttkScalarFieldFromPointCloud(){
    NumberGridPoints = 100;
    Offset = 0.1;
    GaussianKDE = false;
    Bandwidth = 0.01;
    Autobandwidth = false;
}

ttkScalarFieldFromPointCloud::~ttkScalarFieldFromPointCloud(){
}


// transmit abort signals -- to copy paste in other wrappers
bool ttkScalarFieldFromPointCloud::needsToAbort(){
  return GetAbortExecute();
}


// transmit progress status -- to copy paste in other wrappers
int ttkScalarFieldFromPointCloud::updateProgress(const float &progress){

  {
    stringstream msg;
    msg << "[ttkScalarFieldFromPointCloud] " << progress*100 
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }
  
  UpdateProgress(progress);
  return 0;
}

int ttkScalarFieldFromPointCloud::doIt(vtkDataSet *input, vtkUnstructuredGrid *output){

    int numpointsPointcloud = input->GetNumberOfPoints();
    std::vector<std::vector<double>> pointCloudCoordinates;

    for(int i = 0; i < numpointsPointcloud; i++){
        double x = input->GetPoint(i)[0];
        double y = input->GetPoint(i)[1];
        double z = input->GetPoint(i)[2];
        std::vector<double> xyz = {x,y,z};
        pointCloudCoordinates.push_back(xyz);
    }
    
    triangulation.setInputData(vtkImageData::New());
    triangulation.getTriangulation()->setWrapper(this);
    scalarFieldFromPointCloud.setGaussianKDE(GaussianKDE);
    scalarFieldFromPointCloud.setupTriangulation(triangulation.getTriangulation());
    scalarFieldFromPointCloud.setOffset(Offset);
    scalarFieldFromPointCloud.setAutoBandwidth(Autobandwidth);
    scalarFieldFromPointCloud.setBandwidth(Bandwidth);
    scalarFieldFromPointCloud.setNumberGridPoints(NumberGridPoints);
    scalarFieldFromPointCloud.setNumberCloudPoints(numpointsPointcloud);
    scalarFieldFromPointCloud.setInputPointCloud(static_cast<void*>(pointCloudCoordinates.data()));
    scalarFieldFromPointCloud.preprocess();

    vtkDoubleArray* outputScalarField = vtkDoubleArray::New();
    outputScalarField->vtkDoubleArray::SetNumberOfComponents(1);
    outputScalarField->vtkDoubleArray::SetName("PointCloudScalarField");

    if (scalarFieldFromPointCloud.isPlanar()){
        outputScalarField->vtkDoubleArray::SetNumberOfTuples(pow(NumberGridPoints,2));
    }
    else{
        outputScalarField->vtkDoubleArray::SetNumberOfTuples(pow(NumberGridPoints,3));
    }

    scalarFieldFromPointCloud.setOutputScalarFieldPointer(outputScalarField->GetVoidPointer(0));
    scalarFieldFromPointCloud.execute();

    //Get a triangulated unstructured grid. 
    vtkSmartPointer<vtkUnstructuredGrid> triangleUnstructGrid =
        triangulation.getVtkUnstructuredGrid();
    
    triangleUnstructGrid->GetPointData()->AddArray(outputScalarField);
    output->ShallowCopy(triangleUnstructGrid);

}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkScalarFieldFromPointCloud::RequestData(vtkInformation *request, 
  vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;
  
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::GetData(outputVector);
  
  doIt(input, output);
  
  {
    stringstream msg;
    msg << "[ttkScalarFieldFromPointCloud] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  
  return 1;
}
