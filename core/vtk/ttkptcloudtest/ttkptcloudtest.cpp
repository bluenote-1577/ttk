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

std::vector<double> ttkptcloudtest::closestNeighbours(vtkDataSet *input, bool isPlanar){
    std::vector<double> neighbourdistance;

//    Jan Orava(2016) - Nearest Neighbour Estimation - choice for optimal k.
    int num_neighbours;

    if (isPlanar){
        num_neighbours = floor(pow((0.587 * pow(input->GetNumberOfPoints(),0.8)),0.5));
    }
    else{
        num_neighbours = floor(pow((0.587 * pow(input->GetNumberOfPoints(),0.8)),1.0/3.0));
    }

    stringstream msg;
    msg << "[ttkptcloudtest] " << "Nearest neighbour k value = " << num_neighbours << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);

    if (num_neighbours == 0){
        num_neighbours = 1;
    }

    //For each point, populate a list sorted by closest neighbour. 
    for(int i = 0; i < input->GetNumberOfPoints(); i++){

        double x = input->GetPoint(i)[0];
        double y = input->GetPoint(i)[1];
        double z = input->GetPoint(i)[2];
        std::vector<double> nclosest_neighbours;

        for (int j = 0; j < input->GetNumberOfPoints(); j++){
            if(i != j){
                double x_check = input->GetPoint(j)[0];
                double y_check = input->GetPoint(j)[1];
                double z_check = input->GetPoint(j)[2];

                double dist = sqrt(pow(x-x_check,2)+pow(y-y_check,2)+pow(z-z_check,2));

                if(dist < 0.0000001){
                    //Don't want to divide by zero. 
                    continue;
                }

                auto loc = std::lower_bound(nclosest_neighbours.begin(), nclosest_neighbours.end(),dist);
                nclosest_neighbours.insert(loc,dist);
            }
        }

        /*A possible other way to calculate the bandwidth. 
         
        double cumsum = 0;
        for(int z = 0; z < num_neighbours; z++){
            cumsum += nclosest_neighbours[z];
        }
        cumsum = cumsum/num_neighbours;
        neighbourdistance.push_back(cumsum);

       */
        neighbourdistance.push_back(nclosest_neighbours[num_neighbours]);

    }

    /* The below section caps the bandwidth of points with points
     * with extremely close/far neighbours.
    Mean = 0;
    for(int i = 0; i < neighbourdistance.size(); i++){
        Mean += neighbourdistance[i];
    }
    Mean = Mean / neighbourdistance.size();
    std::cerr << "Mean distance " << Mean << '\n';

    for(int i = 0; i < neighbourdistance.size(); i++){
        if(neighbourdistance[i] <  Mean/3){

            std::cerr << neighbourdistance[i] << " 3/ Mean \n";
            neighbourdistance[i] = Mean/3;
        }

        else if (neighbourdistance[i] >  Mean*3){
            std::cerr << neighbourdistance[i] << " 3x Mean \n";
            neighbourdistance[i] = Mean*3;

        }
    }
    */

    return neighbourdistance;
}

int ttkptcloudtest::doIt(vtkDataSet *input, vtkUnstructuredGrid *output){

    double irrelevantExponent = 15;
    double maxDistanceBound = 100000000;
    double minx =maxDistanceBound;
    double miny =maxDistanceBound;
    double minz =maxDistanceBound;
    double maxx =-maxDistanceBound;
    double maxy =-maxDistanceBound;
    double maxz =-maxDistanceBound;
    double e = 2.718281828459;
    double pi = 3.1415926535897;

    int numpointsPointcloud = input->GetNumberOfPoints();
    std::vector<std::vector<double>> pointCloudCoordinates;

    //Detect the bounds for the grid. 
    for(int i = 0; i < numpointsPointcloud; i++){
        double x = input->GetPoint(i)[0];
        double y = input->GetPoint(i)[1];
        double z = input->GetPoint(i)[2];
        std::vector<double> xyz = {x,y,z};
        pointCloudCoordinates.push_back(xyz);

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

    double xDist = maxx-minx;
    double yDist = maxy-miny;
    double zDist = maxz-minz;

    double maxDist;

    if(yDist > zDist){
        if (yDist > xDist){
            maxDist = yDist;
        }
        else{
            maxDist = xDist;
        }
    }
    else{
        if(zDist > xDist){
            maxDist = zDist;
        }
        else{
            maxDist = xDist;
        }
    }

    bool xPlanar = false;
    bool yPlanar = false;
    bool zPlanar = false;

    if(xDist < 0.000001){
        xPlanar = true;
    }

    if(yDist < 0.000001){
        yPlanar = true;
    }

    if(zDist < 0.000001){
        zPlanar = true;
    }

    vtkSmartPointer<vtkImageData> grid = vtkSmartPointer<vtkImageData>::New();

    int gridPoints = NumberGridPoints;
    double origx = minx - Offset*xDist;
    double origy = miny - Offset*yDist;
    double origz = minz - Offset*zDist;
    double scaledBandwidth = Bandwidth * sqrt(xDist*xDist + yDist*yDist + zDist*zDist);

    if(xPlanar){
        origx = 0;
    }

    if(yPlanar){
        origy = 0;
    }

    if(zPlanar){
        origz = 0;
    }


    int totalPointsGrid = pow(gridPoints,3);

    if(xPlanar || yPlanar || zPlanar){
        totalPointsGrid = pow(gridPoints,2);
    }

//    std::cerr << totalPointsGrid;

    vtkSmartPointer<vtkDoubleArray> dataarray= vtkSmartPointer<vtkDoubleArray>::New();
    dataarray->vtkDoubleArray::SetNumberOfComponents(1);
    dataarray->vtkDoubleArray::SetNumberOfTuples(totalPointsGrid);
    dataarray->vtkDoubleArray::SetName("dist_field");

    grid->SetOrigin(origx,origy,origz);
    grid->SetSpacing((xDist*(1 + 2*Offset))/NumberGridPoints,(yDist*(1 + 2*Offset))/NumberGridPoints,(zDist*(1 + 2*Offset))/NumberGridPoints);
    
    double xSpace = (xDist*(1 + 2*Offset))/NumberGridPoints;
    double ySpace = (yDist*(1 + 2*Offset))/NumberGridPoints;
    double zSpace = (zDist*(1 + 2*Offset))/NumberGridPoints;

    grid->SetDimensions(gridPoints,gridPoints,gridPoints);
    if(xPlanar){
        grid->SetDimensions(1,gridPoints,gridPoints);
    }

    if(yPlanar){
        grid->SetDimensions(gridPoints,1,gridPoints);
    }

    if(zPlanar){
        grid->SetDimensions(gridPoints,gridPoints,1);
    }

    grid->GetPointData()->SetNumberOfTuples(totalPointsGrid);
    double bandwidthSquared = scaledBandwidth * scaledBandwidth;

    std::vector<double> kNeighbourDistances;
    if(Autobandwidth){
        kNeighbourDistances = closestNeighbours(input, xPlanar || yPlanar || zPlanar);
    }

    #pragma omp parallel for num_threads(threadNumber_)
    for(int i = 0; i < totalPointsGrid; i++){

        double* gridPoint;
        double xGridPoint,yGridPoint,zGridPoint;

        #pragma omp critical
        { 
            gridPoint= grid->GetPoint(i);
            xGridPoint = gridPoint[0];
            yGridPoint = gridPoint[1];
            zGridPoint = gridPoint[2];
        }

        double nearestNeighbourDistSquare = maxDistanceBound;
        double scalarValue = 0;

        for(int k = 0; k < numpointsPointcloud; k++){ 
            double x,y,z;
            x = pointCloudCoordinates[k][0];
            y = pointCloudCoordinates[k][1];
            z = pointCloudCoordinates[k][2];

            double distanceSquare = (pow((xGridPoint - x),2) + pow((yGridPoint - y),2) + pow((zGridPoint - z),2));

            if(GaussianKDE){
                double bandwidth = scaledBandwidth;
                if(Autobandwidth){
                    //By lowering the exponent, we get a more spread out distribution. 
                    double exponent = 1.0;
                    bandwidth = pow(kNeighbourDistances[k],exponent) * pow(maxDist,1-exponent);
                }

                if((distanceSquare)/(2 * bandwidth* bandwidth) < irrelevantExponent){

                    bandwidthSquared = bandwidth * bandwidth;

                    //2-D density. 
                    if(xPlanar || yPlanar || zPlanar){
                        scalarValue += 100.0/ bandwidthSquared / (2*pi) * pow(e,-distanceSquare/(2 * bandwidthSquared));
                    }
                    //3-D densties have a different multiplication factor.
                    else{
                        scalarValue += 100.0/ pow(bandwidth * 2 * pi,3.0/2.0) * pow(e,-distanceSquare/(2*bandwidthSquared));
                    }

                }
            }

            if(nearestNeighbourDistSquare > distanceSquare){
                nearestNeighbourDistSquare = distanceSquare;
            }
        }

        #pragma omp critical
        {
            double* tup = new double[1];
            if(GaussianKDE){
//                scalarValue /= numpointsPointcloud;
                tup[0] = scalarValue;
            }
            else{
                tup[0] = sqrt(nearestNeighbourDistSquare);
            }
            dataarray->SetTuple(i,tup);
        }
    }

    triangulation.setInputData(grid);
    vtkSmartPointer<vtkUnstructuredGrid> triangle_unstruct_grid =
        triangulation.getVtkUnstructuredGrid();
    triangle_unstruct_grid->GetPointData()->AddArray(dataarray);
    grid->GetPointData()->AddArray(dataarray);
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
