#ifndef _POINTDISTFIELD_H
#define _POINTDISTFIELD_H

// base code includes
#include<Wrapper.h>
#include<Geometry.h>
#include<Triangulation.h>


namespace ttk{

  class PointDistField : public Debug{

    public:
       PointDistField();
      ~PointDistField();

      template <typename dataType>
       int execute();

      template <typename dataType>
       int preprocess();

      inline bool isPlanar() const{
          return isPlanar_;
      }

      inline int setNumberGridPoints(SimplexId numberGridPoints){
        numberGridPoints_ = numberGridPoints;
        return 0;
      }

      inline int setOutputScalarFieldPointer(void *data){
        outputScalarFieldPointer_=data;
        return 0;
      }

      inline int setNumberCloudPoints(SimplexId numberCloudPoints){
        numberCloudPoints_ = numberCloudPoints;
        return 0;
      }

      inline int setOffset(double offset){
        offset_ = offset;
        return 0;
      }

      inline int setGaussianKDE(bool gaussianKDE){
        gaussianKDE_ = gaussianKDE;
        return 0;
      }

      inline int setBandwidth(double bandwidth){
        bandwidth_ = bandwidth;
        return 0;
      }

      inline int setAutoBandwidth(bool autoBandwidth){
        autoBandwidth_ = autoBandwidth;
        return 0;
      }
      
      inline int setupTriangulation(Triangulation* triangulation){
        triangulation_=triangulation;
        return 0;
      }

      inline int setInputPointCloud(void *data){
        inputPointCloud_=data;
        return 0;
      }




    protected:
        Triangulation* triangulation_;
        void* inputPointCloud_;
        SimplexId numberCloudPoints_;
        SimplexId numberGridPoints_;
        double offset_;
        bool gaussianKDE_;
        double bandwidth_;
        bool autoBandwidth_;
        void* outputScalarFieldPointer_;
    private:
        double diagonalLength;
        double maxDist;
        double scaledBandwidth;
        bool donePreprocessing;
        SimplexId totalPointsGrid;
        double maxDistanceBound;
        double irrelevantExponent;
        double e;
        double pi;
        bool isPlanar_;
        template <typename dataType>
        std::vector<dataType> closestNeighbours
        (std::vector<dataType> *input) const;
    };
}

template <typename dataType>
int ttk::PointDistField::preprocess(){
    double minx =maxDistanceBound;
    double miny =maxDistanceBound;
    double minz =maxDistanceBound;
    double maxx =-maxDistanceBound;
    double maxy =-maxDistanceBound;
    double maxz =-maxDistanceBound;

    std::vector<dataType>* inputPointCloud = static_cast<std::vector<dataType>*>(inputPointCloud_);

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
    maxDist;

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
    std::cout << "[pointDistField] diagonal length : " << diagonalLength << '\n';

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
        triangulation_->setInputGrid(float(origx),float(origy),float(origz),float(xSpace),float(ySpace),float(zSpace), gridPoints,gridPoints,1);
    }
    else{
        triangulation_->setInputGrid(origx,origy,origz,xSpace,ySpace,zSpace, gridPoints,gridPoints,gridPoints);
    }

    isPlanar_ = zPlanar || xPlanar || yPlanar;
    donePreprocessing = true;
    return 0;
}


template <typename dataType>
int ttk::PointDistField::execute(){

    if(!donePreprocessing){
    std::cerr << "[PointDistField] Error : Preprocessing not completed." << std::endl;
        return -1;
    }

    std::vector<dataType>* inputPointCloud = static_cast<std::vector<dataType>*>(inputPointCloud_);
    dataType* scalars=static_cast<dataType*>(outputScalarFieldPointer_);
    std::vector<dataType> kNeighbourDistances;

    double mean_bandwidth = 0;
    double stdeviation = 0;
    if(autoBandwidth_){
        kNeighbourDistances = closestNeighbours(inputPointCloud);
        for (SimplexId i = 0; i < numberCloudPoints_ ; i++){
            mean_bandwidth += kNeighbourDistances[i];
            stdeviation += kNeighbourDistances[i] * kNeighbourDistances[i];
        }

        stdeviation /= numberCloudPoints_ - 1;
        stdeviation = sqrt(stdeviation);

        mean_bandwidth/= numberCloudPoints_;
        std::cout << "[PointDistField] Automatic Bandwidth set as : " << mean_bandwidth << ".\n";
        std::cout << "[PointDistField] Automatic+STDEV Bandwidth set as : " << mean_bandwidth + stdeviation/2<< ".\n";
    }

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
                    double exponent = 1;
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

        dataType scalarValToPush = scalarValue; 
        if(!gaussianKDE_){
            scalarValToPush = sqrt(nearestNeighbourDistSquare);
        }
        
        #pragma omp critical
        {
//            std::cerr << scalarValToPush<< '\n';
            scalars[i] = scalarValToPush;
        }
    }

    return 0;
}

template <typename dataType>
std::vector<dataType> ttk::PointDistField::closestNeighbours(std::vector<dataType>* input) const {
    std::vector<dataType> neighbourdistance;

//    Jan Orava(2016) - Nearest Neighbour Estimation - choice for optimal k.
    SimplexId num_neighbours;

    if (isPlanar_){
        num_neighbours = floor(pow((0.587 * pow(numberCloudPoints_,0.8)),0.5));
    }
    else{
        num_neighbours = floor(pow((0.587 * pow(numberCloudPoints_,0.8)),1.0/3.0));
    }

    std::stringstream msg;
    msg << "[ttkptcloudtest] " << "Nearest neighbour k value = " << num_neighbours << std::endl;
    dMsg(std::cout, msg.str(), advancedInfoMsg);

    if (num_neighbours == 0){
        num_neighbours = 1;
    }

    //For each point, populate a list sorted by closest neighbour. 
    for(SimplexId i = 0; i < numberCloudPoints_; i++){

        double x = input[i][0];
        double y = input[i][1];
        double z = input[i][2];
        std::vector<double> nclosest_neighbours;

        for (SimplexId j = 0; j < numberCloudPoints_; j++){
            if(i != j){
                double x_check = input[j][0];
                double y_check = input[j][1];
                double z_check = input[j][2];

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

    double sum = 0;
//    for(SimplexId i = 0; i < neighbourdistance.size(); i++){
//        std::cout << neighbourdistance[i] << ",";
//        sum += neighbourdistance[i];
//    }

    sum = sum/neighbourdistance.size();
    std::cout << "Mean neighdist = " << sum << " in paraview use : " << sum / diagonalLength;



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

#endif // DISTANCEFIELD_H
