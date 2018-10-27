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

std::vector<double> ttkptcloudtest::closestNeighbours(vtkDataSet *input){
    std::vector<double> neighbourdistance;

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
                    //Dont' want to divide by zero. 
                    continue;
                }

                auto loc = std::lower_bound(nclosest_neighbours.begin(), nclosest_neighbours.end(),dist);
                nclosest_neighbours.insert(loc,dist);
            }
        }

        int size = nclosest_neighbours.size();
        int num_neighbours = input->GetNumberOfPoints() / 100;

        if (num_neighbours < 3){
            num_neighbours = 3;
        }

        if(size> num_neighbours){
            double cumsum = 0;
            for(int z = 0; z < num_neighbours; z++){
                cumsum += nclosest_neighbours[z];
            }
            neighbourdistance.push_back(cumsum/num_neighbours);

//            std::cerr << cumsum/num_neighbours<< '\n';
        }
        else{
            std::cerr << "Less than 10 points" << '\n';
        }
    }

    double mean = 0;
    for(int i = 0; i < neighbourdistance.size(); i++){
        mean += neighbourdistance[i];
    }
    mean = mean / neighbourdistance.size();
    std::cerr << "mean distance " << mean << '\n';

    for(int i = 0; i < neighbourdistance.size(); i++){
        if(neighbourdistance[i] <  mean/3){

            std::cerr << neighbourdistance[i] << " 3/ mean \n";
            neighbourdistance[i] = mean/3;
        }

        else if (neighbourdistance[i] >  mean*3){
            std::cerr << neighbourdistance[i] << " 3x mean \n";
            neighbourdistance[i] = mean*3;

        }
    }

    return neighbourdistance;
}

int ttkptcloudtest::doIt(vtkDataSet *input, vtkUnstructuredGrid *output){

    double minx = 10000000;
    double miny = 10000000;
    double minz = 10000000;
    double maxx = -10000000;
    double maxy = -10000000;
    double maxz = -10000000;
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
//    std::cerr << maxy << miny << '\n';
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
    double orig_x = minx - Offset*xdist;
    double orig_y = miny - Offset*ydist;
    double orig_z = minz - Offset*zdist;
    double scaled_Bandwidth = Bandwidth * sqrt(xdist*xdist + ydist*ydist + zdist*zdist);

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
    grid->SetSpacing((xdist*(1 + 2*Offset))/NumberGridPoints,(ydist*(1 + 2*Offset))/NumberGridPoints,(zdist*(1 + 2*Offset))/NumberGridPoints);
    
    double xspace = (xdist*(1 + 2*Offset))/NumberGridPoints;
    double yspace = (ydist*(1 + 2*Offset))/NumberGridPoints;
    double zspace = (zdist*(1 + 2*Offset))/NumberGridPoints;

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
    double bandsquared = scaled_Bandwidth * scaled_Bandwidth;

    std::vector<double> neighbourdistance;
    if(Autobandwidth){
        neighbourdistance = closestNeighbours(input);
    }

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
                if((distsquare)/(2 * scaled_Bandwidth * scaled_Bandwidth) < 10){
                    double bandwidth = scaled_Bandwidth;
                    if(Autobandwidth){
                        bandwidth = neighbourdistance[k];
                    }
                    bandsquared = bandwidth * bandwidth;
                    scalarValue += 1/(2.50662827 * bandwidth) * pow(e,-distsquare/(2 * bandsquared));
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

    triangulation.setInputData(grid);

    vtkSmartPointer<vtkUnstructuredGrid> triangle_unstruct_grid =
        triangulation.getVtkUnstructuredGrid();

    triangle_unstruct_grid->GetPointData()->AddArray(dataarray);

    grid->GetPointData()->AddArray(dataarray);
    //Convert to unstructured grid
//    vtkSmartPointer<vtkAppendFilter> appendfilter = 
//        vtkSmartPointer<vtkAppendFilter>::New();
//    appendfilter->AddInputData(grid);
//    appendfilter->Update();

//    vtkSmartPointer<vtkUnstructuredGrid> outgrid = appendfilter->GetOutput();
    

    output->ShallowCopy(triangle_unstruct_grid);

}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkptcloudtest::RequestData(vtkInformation *request, 
  vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;
  
  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), 1);
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
