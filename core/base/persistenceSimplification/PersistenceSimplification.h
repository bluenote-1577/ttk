/// \ingroup base
/// \class ttk::PersistenceSimplification 
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %persistenceSimplification processing package.
///
/// %PersistenceSimplification is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkPersistenceSimplification.cpp %for a usage example.

#pragma once

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>
#include                  <PersistenceDiagram.h>
#include                  <TopologicalSimplification.h>
#include                  <FTMTreePP.h>


namespace ttk{
  
  class PersistenceSimplification : public Debug{

    public:
        
      PersistenceSimplification();
      
      ~PersistenceSimplification();

      /// Execute the package.
      /// \pre If this TTK package uses ttk::Triangulation for fast mesh 
      /// traversals, the function setupTriangulation() must be called on this 
      /// object prior to this function, in a clearly distinct pre-processing 
      /// steps. An error will be returned otherwise.
      /// \note In such a case, it is recommended to exclude 
      /// setupTriangulation() from any time performance measurement.
      /// \param argment Dummy integer argument.
      /// \return Returns 0 upon success, negative values otherwise.
      template <typename scalarType, typename idType>
        int execute();

      template <typename scalarType>
        scalarType automaticThreshold(std::vector<std::tuple<SimplexId, SimplexId, scalarType>> pairs);
    
      /// Pass a pointer to an input array representing a scalarfield.
      /// The expected format for the array is the following:
      /// <vertex0-component0> <vertex0-component1> ... <vertex0-componentN>
      /// <vertex1-component0> <vertex1-component1> ... <vertex1-componentN>
      /// <vertexM-component0> <vertexM-component1> ... <vertexM-componentN>.
      /// The array is expected to be correctly allocated. 
      /// \param data Pointer to the data array.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa setVertexNumber() and setDimensionNumber().
      inline int setInputDataPointer(void *data){
        inputData_ = data;
        return 0;
      }

      inline int setInputOffsetDataPointer(void *data){
        offsetData_ = data;
        return 0;
      }

      inline int setOutputOffsetDataPointer(void *data){
        outputOffsetData_ = data;
        return 0;
      }

      // TODO_RC

      inline int setAutoOption(bool arg){
        AutoOption = arg;
        return 0;
      }

      inline int setDistinctOption(bool arg){
        DistinctOption = arg;
        return 0;
      }

      inline int setUseMinOption(int arg){
        UseMinOption = arg;
        return 0;
      }

      inline int setUseMaxOption(bool arg){
        UseMaxOption = arg;
        return 0;
      }

      inline int setCountMinOption(bool arg){
        CountMinOption = arg;
        return 0;
      }

      inline int setCountMaxOption(bool arg){
        CountMaxOption = arg;
        return 0;
      }

      inline int setCountAllOption(bool arg){
        CountAllOption = arg;
        return 0;
      }

      inline int setCountMinArg(int arg){
        CountMinArg = arg;
        return 0;
      }

      inline int setCountMaxArg(int arg){
        CountMaxArg = arg;
        return 0;
      }

      inline int setCountAllArg(int arg){
        CountAllArg = arg;
        return 0;
      }

      inline int setPersistenceMinArg(double arg){
        PersistenceMinArg = arg;
        return 0;
      }

      inline int setPersistenceMaxArg(double arg){
        PersistenceMaxArg = arg;
        return 0;
      }

      inline int setPersistenceAllArg(double arg){
        PersistenceAllArg = arg;
        return 0;
      }

      /// Pass a pointer to an output array representing a scalar field.
      /// The expected format for the array is the following:
      /// <vertex0-component0> <vertex0-component1> ... <vertex0-componentN>
      /// <vertex1-component0> <vertex1-component1> ... <vertex1-componentN>
      /// <vertexM-component0> <vertexM-component1> ... <vertexM-componentN>.
      /// The array is expected to be correctly allocated.
      /// \param data Pointer to the data array.
      /// \return Returns 0 upon success, negative values otherwise. 
      /// \sa setVertexNumber() and setDimensionNumber().
      inline int setOutputDataPointer(void *data){
        outputData_ = data;
        return 0;
      }
     
      // General documentation info:
      //
      /// Setup a (valid) triangulation object for this TTK base object.
      ///
      /// \pre This function should be called prior to any usage of this TTK 
      /// object, in a clearly distinct pre-processing step that involves no 
      /// traversal or computation at all. An error will be returned otherwise.
      ///
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement. Therefore, it is recommended to 
      /// call this function ONLY in the pre-processing steps of your program. 
      /// Note however, that your triangulation object must be valid when 
      /// calling this function (i.e. you should have filled it at this point, 
      /// see the setInput*() functions of ttk::Triangulation). See ttkPersistenceSimplification 
      /// for further examples.
      ///
      /// \param triangulation Pointer to a valid triangulation.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa ttk::Triangulation
      //
      //
      // Developer info:
      // ttk::Triangulation is a generic triangulation representation that 
      // enables fast mesh traversal, either on explicit triangulations (i.e.
      // tet-meshes) or implicit triangulations (i.e. low-memory footprint 
      // implicit triangulations obtained from regular grids).
      //
      // Not all TTK packages need such mesh traversal features. If your 
      // TTK package needs any mesh traversal procedure, we recommend to use 
      // ttk::Triangulation as described here.
      //
      // Each call to a traversal procedure of ttk::Triangulation 
      // must satisfy some pre-condition (see ttk::Triangulation for more 
      // details). Such pre-condition functions are typically called from this
      // function. 
      inline int setupTriangulation(Triangulation *triangulation){
        triangulation_ = triangulation;
       
        if(triangulation_){
          
          // TODO-1
          // Pre-condition functions.
          // Call all the required pre-condition functions here!
          // for example:
          triangulation_->preprocessVertexNeighbors();
          // end of TODO-1
          
        }
        
        return 0;
      }
    
    protected:
    
      void                  *inputData_, *outputData_, *offsetData_, *outputOffsetData_;
      Triangulation         *triangulation_;

      // TODO_RC
      int CountAllArg, CountMinArg, CountMaxArg;
      double PersistenceMaxArg, PersistenceMinArg, PersistenceAllArg;
      bool DistinctOption, AutoOption, CountMinOption, CountMaxOption, CountAllOption;
      bool UseMinOption, UseMaxOption;
  };
}

// if the package is a pure template class, uncomment the following line
// #include                  <PersistenceSimplification.cpp>

// template functions
template <typename scalarType, typename idType> int ttk::PersistenceSimplification::execute(){

  Timer t;
  
  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_)
    return -1;
  if(!inputData_)
    return -2;
  if(!outputData_)
    return -3;
#endif

  scalarType *outputData = (scalarType *) outputData_;
  scalarType *inputData = (scalarType *) inputData_;
  // SimplexId *offsetData = (SimplexId *) offsetData_;
  // SimplexId *outputOffsetData = (SimplexId *) outputOffsetData_;
  
  SimplexId vertexNumber = triangulation_->getNumberOfVertices();

  // init the output -- to adapt
  for(SimplexId i = 0; i < vertexNumber; i++){
    outputData[i] = inputData[i];
  }
  
  // the following open-mp processing is only relevant for embarrassingly 
  // parallel algorithms (such as smoothing) -- to adapt
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) 
#endif
  for(SimplexId i = 0; i < vertexNumber; i++){
    // TODO-2
    // processing here!
    // end of TODO-2
  }

  const SimplexId numberOfVertices=triangulation_->getNumberOfVertices();
  // convert offsets into a valid format for contour tree
  std::vector<SimplexId> voffsets(numberOfVertices);
  // std::copy(offsets,offsets+numberOfVertices,voffsets.begin()); // Was used in persistence curve

  // get contour tree
  ftm::FTMTreePP contourTree;
  contourTree.setupTriangulation(triangulation_, false);
  contourTree.setVertexScalars(inputData);
  contourTree.setTreeType(ftm::TreeType::Join_Split);
  contourTree.setVertexSoSoffsets(voffsets.data());
  contourTree.setSegmentation(false);
  contourTree.setThreadNumber(threadNumber_);
  contourTree.build<scalarType,idType>();

  // get persistence pairs
  std::vector<std::tuple<SimplexId, SimplexId, scalarType>> JTPairs;
  std::vector<std::tuple<SimplexId, SimplexId, scalarType>> STPairs;
  std::vector<SimplexId> vertexIdentifiers;
  contourTree.computePersistencePairs<scalarType>(JTPairs, true);
  contourTree.computePersistencePairs<scalarType>(STPairs, false);

  // merge pairs
  std::vector<std::tuple<SimplexId, SimplexId, scalarType>> CTPairs(JTPairs.size() +
                                                                        STPairs.size());
  const SimplexId JTSize = JTPairs.size();
  for (SimplexId i = 0; i < JTSize; ++i) {
     const auto& x = JTPairs[i];
     CTPairs[i]    = std::make_tuple(std::get<0>(x), std::get<1>(x), std::get<2>(x));
  }
  const SimplexId STSize = STPairs.size();
  for (SimplexId i = 0; i < STSize; ++i) {
     const auto& x       = STPairs[i];
     CTPairs[JTSize + i] = std::make_tuple(std::get<0>(x), std::get<1>(x), std::get<2>(x));
  }

  {
     auto cmp = [](const std::tuple<SimplexId, SimplexId, scalarType>& a,
                   const std::tuple<SimplexId, SimplexId, scalarType>& b) {
        return std::get<2>(a) < std::get<2>(b);
     };

     std::sort(CTPairs.begin(), CTPairs.end(), cmp);
     // remove the last pair which is present two times (global extrema pair)
     CTPairs.erase(CTPairs.end() - 1);
  }

  // Sort the STPairs so that the persistence thresholding is easy.
  {
     auto cmp = [](const std::tuple<SimplexId, SimplexId, scalarType>& a,
                   const std::tuple<SimplexId, SimplexId, scalarType>& b) {
        return std::get<2>(a) < std::get<2>(b);
     };

     std::sort(STPairs.begin(), STPairs.end(), cmp);
  }

  // Sort the JTPairs so that the persistence thresholding is easy.
  {
     auto cmp = [](const std::tuple<SimplexId, SimplexId, scalarType>& a,
                   const std::tuple<SimplexId, SimplexId, scalarType>& b) {
        return std::get<2>(a) < std::get<2>(b);
     };

     std::sort(JTPairs.begin(), JTPairs.end(), cmp);
  }

  // Sort the CTPairs so that the persistence thresholding is easy.
  {
     auto cmp = [](const std::tuple<SimplexId, SimplexId, scalarType>& a,
                   const std::tuple<SimplexId, SimplexId, scalarType>& b) {
        return std::get<2>(a) < std::get<2>(b);
     };

     std::sort(CTPairs.begin(), CTPairs.end(), cmp);
  }

  // for (SimplexId i = 0; i < int(STPairs.size()); i++){
  //   SimplexId a = std::get<0>(STPairs[i]);
  //   SimplexId b = std::get<1>(STPairs[i]);
  //   scalarType persistence = std::get<2>(STPairs[i]);
  //   std::cerr << "[PersistenceSimplification] a=" << a << ", b=" << b << ", p=" << persistence << ";\n";
  // }
  // std::cerr << todoArg << "\n";

  // ----------------------------------------
  // Algorithm for persistence thresholding
  // int clusterNumber;
  // scalarType persistenceThresh;
  // {
  //   int points_to_check = 2;
  //   float absthresh = 0.025 * std::get<2>(STPairs[STPairs.size()-1]);
  //   float relthresh = 0.1;
  //   int testing = 0;
  //   int against = 1;
  //   int m = 0;

  //   while(against < int(STPairs.size())) {
  //     m = absthresh + std::get<2>(STPairs[testing]) * relthresh;
  //     if(std::get<2>(STPairs[against]) < (m * (against - testing) + std::get<2>(STPairs[testing]))){
  //       testing+= 1;
  //       against = testing;
  //     }
  //     if((against - testing) == points_to_check){
  //       break;
  //     }
  //     against+= 1;
  //   }

  //   clusterNumber = STPairs.size() - testing - 1;
  //   persistenceThresh = std::get<2>(STPairs[testing]) + m;

  //   std::cerr << "[PersistenceSimplification] Found " << clusterNumber << " clusters\n";
  //   std::cerr << "[PersistenceSimplification] Thresholding at " << persistenceThresh << "\n";
  // }
  // ----------------------------------------
  
  if(DistinctOption){
    scalarType persistence;

    // MIN PAIRS
    if(UseMinOption){
      if(!AutoOption)
        persistence = PersistenceMinArg;
      else{
        std::cerr << "[PersistenceSimplification] Thresholding Min Pairs:\n";
        persistence = automaticThreshold(JTPairs);
      }
      if(AutoOption || !CountMinOption){
        for (SimplexId i = 0; i < int(JTPairs.size()); i++){
          if(std::get<2>(JTPairs[i]) > persistence){
            SimplexId a = std::get<0>(JTPairs[i]);
            SimplexId b = std::get<1>(JTPairs[i]);
            vertexIdentifiers.push_back(b);
            vertexIdentifiers.push_back(a);
          }
        }
      }
      else{
        for (SimplexId i = int(JTPairs.size()) - CountMinArg; i < int(JTPairs.size()); i++){
          SimplexId a = std::get<0>(JTPairs[i]);
          SimplexId b = std::get<1>(JTPairs[i]);
          vertexIdentifiers.push_back(b);
          vertexIdentifiers.push_back(a);
        }
      }
    }

    // MAX PAIRS
    if(UseMaxOption){
      if(!AutoOption)
        persistence = PersistenceMaxArg;
      else{
        std::cerr << "[PersistenceSimplification] Thresholding Max Pairs:\n";
        persistence = automaticThreshold(STPairs);
      }
      if(AutoOption || !CountMaxOption){
        for (SimplexId i = 0; i < int(STPairs.size()); i++){
          if(std::get<2>(STPairs[i]) > persistence){
            SimplexId a = std::get<0>(STPairs[i]);
            SimplexId b = std::get<1>(STPairs[i]);
            vertexIdentifiers.push_back(b);
            vertexIdentifiers.push_back(a);
          }
        }
      }
      else{
        for (SimplexId i = int(STPairs.size()) - CountMaxArg; i < int(STPairs.size()); i++){
          SimplexId a = std::get<0>(STPairs[i]);
          SimplexId b = std::get<1>(STPairs[i]);
          vertexIdentifiers.push_back(b);
          vertexIdentifiers.push_back(a);
        }
      }
    }
  }

  else{
    scalarType persistence;

    // COMBINED PAIRS
    if(!AutoOption)
      persistence = PersistenceAllArg;
    else{
      persistence = automaticThreshold(CTPairs);
    }
    if(AutoOption || !CountAllOption){
      for (SimplexId i = 0; i < int(CTPairs.size()); i++){
        if(std::get<2>(CTPairs[i]) > persistence){
          SimplexId a = std::get<0>(CTPairs[i]);
          SimplexId b = std::get<1>(CTPairs[i]);
          vertexIdentifiers.push_back(b);
          vertexIdentifiers.push_back(a);
        }
      }
    }
    else{
      for (SimplexId i = int(CTPairs.size()) - CountAllArg; i < int(CTPairs.size()); i++){
        SimplexId a = std::get<0>(CTPairs[i]);
        SimplexId b = std::get<1>(CTPairs[i]);
        vertexIdentifiers.push_back(b);
        vertexIdentifiers.push_back(a);
      }
    }
  }

  vertexIdentifiers.push_back(0);
  std::cerr << "[PersistenceSimplification] Using " << vertexIdentifiers.size() << " constraints\n";

  TopologicalSimplification topologicalSimplification;
  topologicalSimplification.setupTriangulation(triangulation_);
  // topologicalSimplification.setWrapper(this); // Don't need this
  topologicalSimplification.setInputScalarFieldPointer(inputData);
  topologicalSimplification.setVertexIdentifierScalarFieldPointer(&vertexIdentifiers[0]);
  topologicalSimplification.setInputOffsetScalarFieldPointer(offsetData_);
  topologicalSimplification.setOutputScalarFieldPointer(outputData);
  topologicalSimplification.setOutputOffsetScalarFieldPointer(outputOffsetData_);
  topologicalSimplification.setVertexNumber(vertexNumber);
  topologicalSimplification.setConstraintNumber(vertexIdentifiers.size());
  topologicalSimplification.setConsiderIdentifierAsBlackList(false);
  topologicalSimplification.setAddPerturbation(false);

  topologicalSimplification.execute<scalarType,idType>();

  {
    std::stringstream msg;
    msg << "[PersistenceSimplification] Data-set (" << vertexNumber
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }
  
  return 0;
}

template <typename scalarType> 
scalarType ttk::PersistenceSimplification::automaticThreshold(
  std::vector<std::tuple<SimplexId, SimplexId, scalarType>> pairs){
  int clusterNumber;
  scalarType persistenceThresh;

  int points_to_check = 2;
  float absthresh = 0.025 * std::get<2>(pairs[pairs.size()-1]);
  float relthresh = 0.1;
  int testing = 0;
  int against = 1;
  int m = 0;

  while(against < int(pairs.size())) {
    m = absthresh + std::get<2>(pairs[testing]) * relthresh;
    if(std::get<2>(pairs[against]) < (m * (against - testing) + std::get<2>(pairs[testing]))){
      testing+= 1;
      against = testing;
    }
    if((against - testing) == points_to_check){
      break;
    }
    against+= 1;
  }

  clusterNumber = pairs.size() - testing - 1;
  persistenceThresh = std::get<2>(pairs[testing]) + m;

  std::cerr << "[PersistenceSimplification] Found " << clusterNumber << " pairs\n";
  std::cerr << "[PersistenceSimplification] Thresholding at " << persistenceThresh << "\n";
  return persistenceThresh;
}