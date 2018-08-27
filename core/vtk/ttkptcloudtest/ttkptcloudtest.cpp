#include                  <ttkptcloudtest.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkptcloudtest)

ttkptcloudtest::ttkptcloudtest(){
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

int ttkptcloudtest::doIt(vtkDataSet *input, vtkRectilinearGrid *output){

    double minx = 10000000;
    double miny = 10000000;
    double minz = 10000000;

    int numpoints = input->GetNumberOfPoints();
    std::vector<std::vector<double>> coords;
    for(int i = 0; i < numpoints; i++){
        double x = input->GetPoint(i)[0];
        double y = input->GetPoint(i)[1];
        double z = input->GetPoint(i)[2];
        std::vector<double> xyz = {x,y,z};
        coords.push_back(xyz);

        if (minx > x){
            minx = x;
        }

        if (miny > y){
            miny = y;
        }

        if (minz > z){
            minz = z;
        }
    }

    double diameter= 0;
    for(int i = 0; i < coords.size() - 1; i++){
        for(int k = i + 1; k < coords.size(); k++){
            double distx = pow((coords[i][0] - coords[k][0]),2);
            double disty = pow((coords[i][1] - coords[k][1]),2);
            double distz = pow((coords[i][2] - coords[k][2]),2);

            double dist = sqrt(distx+disty+distz);
            if(dist > diameter){
                diameter= dist;
            }
        }
    }

    diameter = diameter;

    std::cout << diameter << "diam" << '\n';

    vtkSmartPointer<vtkRectilinearGrid> grid = vtkSmartPointer<vtkRectilinearGrid>::New();

    int gridpoints = 100;
    double inc = diameter*3/gridpoints;

    grid->SetDimensions(gridpoints,gridpoints,gridpoints);

    vtkSmartPointer<vtkDoubleArray> xarray= vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> yarray= vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> zarray= vtkSmartPointer<vtkDoubleArray>::New();

    double x = minx - diameter;
    double y = miny - diameter;
    double z = minz - diameter;

    for(int i = 0; i < gridpoints; i++){
        x += inc;
        xarray->InsertNextValue(x); 
        y += inc;
        yarray->InsertNextValue(y); 
        z += inc;
        zarray->InsertNextValue(z); 
    }

    grid->SetXCoordinates(xarray);
    grid->SetYCoordinates(yarray);
    grid->SetZCoordinates(zarray);

    std::cout << "THERE ARE " << grid->GetNumberOfPoints() << '\n';
    std::cout << "THERE ARE " << grid->GetNumberOfCells() << '\n';

    vtkSmartPointer<vtkDoubleArray> dataarray= vtkSmartPointer<vtkDoubleArray>::New();
    dataarray->vtkDoubleArray::SetNumberOfComponents(1);
    dataarray->vtkDoubleArray::SetNumberOfTuples(pow(gridpoints,3));
    dataarray->vtkDoubleArray::SetName("dist_field");

    for(int i = 0; i < pow(gridpoints,3); i++){
        double* point = grid->GetPoint(i);
        double neighdist = 10000000;

        for(int k = 0; k < coords.size(); k++){ 
            double x,y,z;
            x = coords[k][0];
            y = coords[k][1];
            z = coords[k][2];

            double dist = sqrt(pow((point[0] - x),2) + pow((point[1] - y),2) + pow((point[2] - z),2));

            if(neighdist > dist){
                neighdist = dist;
            }
        }

        double* tup = new double[1];
        tup[0] = neighdist;
        dataarray->vtkDoubleArray::SetTuple(i,tup);
    }

    grid->GetPointData()->AddArray(dataarray);
    output->ShallowCopy(grid);
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkptcloudtest::RequestData(vtkInformation *request, 
  vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;
  
  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkRectilinearGrid *output = vtkRectilinearGrid::GetData(outputVector);
  
  doIt(input, output);
  
  {
    stringstream msg;
    msg << "[ttkptcloudtest] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  
  return 1;
}
