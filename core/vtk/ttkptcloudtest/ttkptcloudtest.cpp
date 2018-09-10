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

int ttkptcloudtest::doIt(vtkDataSet *input, vtkImageData *output){

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

    vtkSmartPointer<vtkImageData> grid = vtkSmartPointer<vtkImageData>::New();

    int gridpoints = 100;
    double inc = diameter*3/gridpoints;


//    vtkSmartPointer<vtkDoubleArray> xarray= vtkSmartPointer<vtkDoubleArray>::New();
//    vtkSmartPointer<vtkDoubleArray> yarray= vtkSmartPointer<vtkDoubleArray>::New();
//    vtkSmartPointer<vtkDoubleArray> zarray= vtkSmartPointer<vtkDoubleArray>::New();
//
    double orig_x = minx - diameter;
    double orig_y = miny - diameter;
    double orig_z = minz - diameter;
//
//    for(int i = 0; i < gridpoints; i++){
//        x += inc;
//        xarray->InsertNextValue(x); 
//        y += inc;
//        yarray->InsertNextValue(y); 
//        z += inc;
//        zarray->InsertNextValue(z); 
//    }
//
//    grid->SetXCoordinates(xarray);
//    grid->SetYCoordinates(yarray);
//    grid->SetZCoordinates(zarray);

    std::cout << "THERE ARE " << grid->GetNumberOfPoints() << '\n';
    std::cout << "THERE ARE " << grid->GetNumberOfCells() << '\n';

    int total_points = pow(gridpoints,3);
    vtkSmartPointer<vtkDoubleArray> dataarray= vtkSmartPointer<vtkDoubleArray>::New();
    dataarray->vtkDoubleArray::SetNumberOfComponents(1);
    dataarray->vtkDoubleArray::SetNumberOfTuples(total_points);
    dataarray->vtkDoubleArray::SetName("dist_field");
    grid->SetOrigin(orig_x,orig_y,orig_z);
    grid->SetSpacing(inc,inc,inc);
    grid->SetDimensions(gridpoints,gridpoints,gridpoints);
    grid->GetPointData()->SetNumberOfTuples(total_points);
    std::cerr << "origin, spacing\n";

    #pragma omp parallel for
    for(int i = 0; i < total_points; i++){
//        std::cerr << "point out";
        double* point;
        double xpoint,ypoint,zpoint;
        #pragma omp critical
        { 
            point = grid->GetPoint(i);
            xpoint = point[0];
            ypoint = point[1];
            zpoint = point[2];
        }
//        std::cerr << point[0] << point[1] << point[2] << '\n';
//        std::cerr << "point out";
        double neighdist = 10000000;

        for(int k = 0; k < coords.size(); k++){ 
            double x,y,z;
            x = coords[k][0];
            y = coords[k][1];
            z = coords[k][2];

            double dist = sqrt(pow((xpoint - x),2) + pow((ypoint - y),2) + pow((zpoint - z),2));

            if(neighdist > dist){
                neighdist = dist;
            }
        }
//        std::cerr << "calc done \n";

//        std::cerr << "scalar pointer\n";
        #pragma omp critical
        {
            double* tup = new double[1];
            tup[0] = neighdist;
            dataarray->SetTuple(i,tup);
        }
        std::cerr << i << '\n';
//        std::cerr << "neighdist\n";
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

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkImageData *output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  //vtkImageData *output = vtkImageData::GetData(outputVector);
  
  doIt(input, output);
  
  {
    stringstream msg;
    msg << "[ttkptcloudtest] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  
  return 1;
}
