#include                  <ttkptcloudtest.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkptcloudtest)

ttkptcloudtest::ttkptcloudtest(){
    NumberGridPoints = 100;
    Offset = 1;
    GaussianKDE = false;
    Bandwidth = 1;
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
    double maxx,maxy,maxz = -1;
    double e = 2.718281828459;
    double pi = 3.1415926535897;

    int numpoints = input->GetNumberOfPoints();
    std::vector<std::vector<double>> coords;
    for(int i = 0; i < numpoints; i++){
        double x = input->GetPoint(i)[0];
        double y = input->GetPoint(i)[1];
        double z = input->GetPoint(i)[2];
        std::vector<double> xyz = {x,y,z};
        coords.push_back(xyz);

        if(maxx < x){
            maxx = x;
        }

        if(maxy < y){
            maxy = y;
        }

        if(maxz < z){
            maxz = z;
        }

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

    double max_coord_dist;
    double xdist = maxx-minx;
    double ydist = maxy-miny;
    double zdist = maxz-minz;

    bool xplanar = false;
    bool yplanar = false;
    bool zplanar = false;

    if(xdist < 0.000001){
        xplanar = true;
    }

    if(ydist < 0.000001){
        yplanar = true;
    }

    if(zdist < 0.000001){
        zplanar = true;
    }

    vtkSmartPointer<vtkImageData> grid = vtkSmartPointer<vtkImageData>::New();

    int gridpoints = NumberGridPoints;
    double orig_x = minx - Offset;
    double orig_y = miny - Offset;
    double orig_z = minz - Offset;

    if(xplanar){
        orig_x = 0;
    }

    if(yplanar){
        orig_y = 0;
    }

    if(zplanar){
        orig_z = 0;
    }


    int total_points;
    total_points = pow(gridpoints,3);

    if(xplanar || yplanar || zplanar){
        total_points = pow(gridpoints,2);
    }

    std::cerr << total_points;


    vtkSmartPointer<vtkDoubleArray> dataarray= vtkSmartPointer<vtkDoubleArray>::New();
    dataarray->vtkDoubleArray::SetNumberOfComponents(1);
    dataarray->vtkDoubleArray::SetNumberOfTuples(total_points);
    dataarray->vtkDoubleArray::SetName("dist_field");

    grid->SetOrigin(orig_x,orig_y,orig_z);
    grid->SetSpacing((xdist+2*Offset)/NumberGridPoints,(ydist+2*Offset)/NumberGridPoints,(zdist+2*Offset)/NumberGridPoints);
    grid->SetDimensions(gridpoints,gridpoints,gridpoints);

    if(xplanar){
        grid->SetDimensions(1,gridpoints,gridpoints);
    }

    if(yplanar){
        grid->SetDimensions(gridpoints,1,gridpoints);
    }

    if(zplanar){
        grid->SetDimensions(gridpoints,gridpoints,1);
    }

    grid->GetPointData()->SetNumberOfTuples(total_points);
    int num_points_in_cloud = coords.size();
    double bandsquared = Bandwidth * Bandwidth;

    #pragma omp parallel for num_threads(threadNumber_)
    for(int i = 0; i < total_points; i++){

        double* point;
        double xpoint,ypoint,zpoint;

        #pragma omp critical
        { 
            point = grid->GetPoint(i);
            xpoint = point[0];
            ypoint = point[1];
            zpoint = point[2];
        }

        double neighdist = 10000000;
        double scalarValue = 0;

        for(int k = 0; k < coords.size(); k++){ 
            double x,y,z;
            x = coords[k][0];
            y = coords[k][1];
            z = coords[k][2];

            double distsquare = (pow((xpoint - x),2) + pow((ypoint - y),2) + pow((zpoint - z),2));

            if(GaussianKDE){
                if((distsquare)/(2 * bandsquared) < 10){
                    scalarValue += 1/(2.50662827 * Bandwidth) * pow(e,-distsquare/(2 * bandsquared));
                }
            }

            if(neighdist > distsquare){
                neighdist = distsquare;
            }
        }

        #pragma omp critical
        {
            double* tup = new double[1];
            if(GaussianKDE){
//                scalarValue /= num_points_in_cloud;
                tup[0] = scalarValue;
            }
            else{
                tup[0] = sqrt(neighdist);
            }
            dataarray->SetTuple(i,tup);
        }
    }

    grid->GetPointData()->AddArray(dataarray);
    output->ShallowCopy(grid);

}

int ttkptcloudtest::RequestInformation(
    vtkInformation *request,
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
  

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  //outInfo->Set(vtkDataObject::SPACING(), DataSpacing, 3);
  //outInfo->Set(vtkDataObject::ORIGIN(), DataOrigin, 3);
  //
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),0,NumberGridPoints-1,0,NumberGridPoints-1,0,NumberGridPoints-1);

  return 1;
}


// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkptcloudtest::RequestData(vtkInformation *request, 
  vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;
  
  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), 1);
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