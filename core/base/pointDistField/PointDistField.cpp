#include<PointDistField.h>

using namespace std;
using namespace ttk;

PointDistField::PointDistField(){
    irrelevantExponent = 15;
    maxDistanceBound = 100000000;
    e = 2.718281828459;
    pi = 3.1415926535897;
    donePreprocessing = false;

}
PointDistField::~PointDistField(){
}

std::vector<double> ttk::PointDistField::closestNeighbours(std::vector<double>* input) const {
    std::vector<double> neighbourdistance;

//    Jan Orava(2016) - Nearest Neighbour Estimation - choice for optimal k.
    SimplexId num_neighbours;

    if (isPlanar_){
        num_neighbours = floor(pow((0.587 * pow(numberCloudPoints_,0.8)),0.5));
    }
    else{
        num_neighbours = floor(pow((0.587 * pow(numberCloudPoints_,0.8)),1.0/3.0));
    }

    {
        std::stringstream msg;
        msg << "[PointDistField] " << "Nearest neighbour k value = " << num_neighbours << std::endl;
        dMsg(std::cout, msg.str(), infoMsg);
    }

    if (num_neighbours == 0){
        num_neighbours = 1;
    }


    Timer autoBandwidth_t;
    //For each point, populate a list sorted by closest neighbour. 
    #pragma omp parallel for num_threads(threadNumber_)
    for(SimplexId i = 0; i < numberCloudPoints_; i++){

        double x = input[i][0];
        double y = input[i][1];
        double z = input[i][2];
        std::list<double> nclosest_neighbours;
//        std::vector<double> nclosest_neighbours;

        for (SimplexId j = 0; j < numberCloudPoints_; j++){
            if(i != j){
                double x_check = input[j][0];
                double y_check = input[j][1];
                double z_check = input[j][2];

                double dist = (pow(x-x_check,2)+pow(y-y_check,2)+pow(z-z_check,2));

                if(dist < 0.0000001){
                    //Don't want to divide by zero. This came up while testing - Jim
                    continue;
                }

                auto loc = std::lower_bound(nclosest_neighbours.begin(), nclosest_neighbours.end(),dist);
                if(std::distance(nclosest_neighbours.begin(), loc) <= num_neighbours){
                    nclosest_neighbours.insert(loc,dist);
                    if(SimplexId(nclosest_neighbours.size()) > num_neighbours + 1){
                        nclosest_neighbours.pop_back();
                    }
                }
            }

        }

        #pragma omp critical
        {
        auto it = nclosest_neighbours.begin();
        std::advance(it, num_neighbours);
        neighbourdistance.push_back(sqrt(*it));
        }


    }

    {
      std::stringstream msg;
      msg << 
        "[PointDistField] Automatic bandwidth guess completed processing in "
        << autoBandwidth_t.getElapsedTime() << " s. " 
        << "(" << threadNumber_ << " thread(s)). " << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
    return neighbourdistance;
}

int ttk::PointDistField::preprocess(){
    double minx =maxDistanceBound;
    double miny =maxDistanceBound;
    double minz =maxDistanceBound;
    double maxx =-maxDistanceBound;
    double maxy =-maxDistanceBound;
    double maxz =-maxDistanceBound;

    std::vector<double>* inputPointCloud = static_cast<std::vector<double>*>(inputPointCloud_);

    for (SimplexId i = 0; i < numberCloudPoints_ ; i++){
        double x = inputPointCloud[i][0];
        double y = inputPointCloud[i][1];
        double z = inputPointCloud[i][2];

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

    SimplexId gridPoints = numberGridPoints_;
    double origx = minx - offset_*xDist;
    double origy = miny - offset_*yDist;
    double origz = minz - offset_*zDist;
    scaledBandwidth = bandwidth_ * sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
    diagonalLength = sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
    //TODO Debug message, remove in final iteration. 
    std::cout << "[PointDistField] diagonal length : " << diagonalLength << '\n';

    if(xPlanar){
        origx = 0;
    }

    if(yPlanar){
        origy = 0;
    }

    if(zPlanar){
        origz = 0;
    }

    totalPointsGrid = pow(gridPoints,3);

    if(xPlanar || yPlanar || zPlanar){
        totalPointsGrid = pow(gridPoints,2);
    }

    double xSpace = (xDist*(1 + 2*offset_))/numberGridPoints_;
    double ySpace = (yDist*(1 + 2*offset_))/numberGridPoints_;
    double zSpace = (zDist*(1 + 2*offset_))/numberGridPoints_;

    if(xPlanar){
        triangulation_->setInputGrid(origx,origy,origz,xSpace,ySpace,zSpace, 1,gridPoints,gridPoints);
    }
    else if(yPlanar){
        triangulation_->setInputGrid(origx,origy,origz,xSpace,ySpace,zSpace, gridPoints,1,gridPoints);
    }
    else if(zPlanar){
        triangulation_->setInputGrid(double(origx),double(origy),double(origz),double(xSpace),double(ySpace),double(zSpace), gridPoints,gridPoints,1);
    }
    else{
        triangulation_->setInputGrid(origx,origy,origz,xSpace,ySpace,zSpace, gridPoints,gridPoints,gridPoints);
    }

    isPlanar_ = zPlanar || xPlanar || yPlanar;
    donePreprocessing = true;
    return 0;
}


int ttk::PointDistField::execute(){

    if(!donePreprocessing){
        std::stringstream msg;
        msg << "[PointDistField] " << "Preprocesing not finished. Error! " << std::endl;
        dMsg(std::cout, msg.str(), fatalMsg);

        return -1;
    }

    std::vector<double>* inputPointCloud = static_cast<std::vector<double>*>(inputPointCloud_);
    double* scalars=static_cast<double*>(outputScalarFieldPointer_);
    std::vector<double> kNeighbourDistances;

    double mean_bandwidth = 0;
    double stdeviation = 0;

    if(autoBandwidth_){
        kNeighbourDistances = closestNeighbours(inputPointCloud);
        for (SimplexId i = 0; i < numberCloudPoints_ ; i++){
            mean_bandwidth += kNeighbourDistances[i];
            stdeviation += kNeighbourDistances[i] * kNeighbourDistances[i];
        }

        //TODO I tried using stdeviation for some calculations, should remove later.
        stdeviation /= numberCloudPoints_ - 1;
        stdeviation = sqrt(stdeviation);

        mean_bandwidth/= numberCloudPoints_;

        {
            std::stringstream msg;
            msg << "[PointDistField] Automatic Bandwidth set as : " << mean_bandwidth/diagonalLength << ".\n";
            dMsg(std::cout, msg.str(), fatalMsg);
        }

    }

    Timer scalarFieldTimer_t;

    #pragma omp parallel for num_threads(threadNumber_)
    for(SimplexId i = 0; i < totalPointsGrid; i++){

        float xGridPoint,yGridPoint,zGridPoint;
        triangulation_->getVertexPoint(i,xGridPoint,yGridPoint,zGridPoint);
        double nearestNeighbourDistSquare = maxDistanceBound;
        double scalarValue = 0;

        for(SimplexId k = 0; k < numberCloudPoints_; k++){ 

            double x,y,z;
            x = inputPointCloud[k][0];
            y = inputPointCloud[k][1];
            z = inputPointCloud[k][2];

            double distanceSquare = (pow((xGridPoint - x),2) + pow((yGridPoint - y),2) + pow((zGridPoint - z),2));

            if(gaussianKDE_){
                double bandwidth = scaledBandwidth;
                if(autoBandwidth_){
                    //By lowering the exponent, we get a more spread out distribution. 
                    //double exponent = 1;
//                    bandwidth = pow(kNeighbourDistances[k],exponent) * pow(maxDist,1-exponent);
//                    bandwidth = kNeighbourDistances[k];
                    bandwidth = mean_bandwidth;
                }

                if((distanceSquare)/(2 * bandwidth* bandwidth) < irrelevantExponent){

                    double bandwidthSquared = bandwidth * bandwidth;

                    //2-D density. 
                    if(isPlanar_){
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

        double scalarValToPush = scalarValue; 
        if(!gaussianKDE_){
            scalarValToPush = sqrt(nearestNeighbourDistSquare);
        }
        
        #pragma omp critical
        {
            scalars[i] = scalarValToPush;
        }
    }

    {
      std::stringstream msg;
      msg << 
        "[PointDistField] Scalar field computation completed processing in "
        << scalarFieldTimer_t.getElapsedTime() << " s. "
        << "(" << threadNumber_ << " thread(s)). " << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }

    return 0;
}



#
