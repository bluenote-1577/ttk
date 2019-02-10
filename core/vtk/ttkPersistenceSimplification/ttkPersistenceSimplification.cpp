#include                  <ttkPersistenceSimplification.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPersistenceSimplification)

int ttkPersistenceSimplification::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  Memory m;
  
  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];
  
  triangulation_ = ttkTriangulation::getTriangulation(input);
 
  #ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_){
    stringstream msg;
    msg << "[ttkPersistenceSimplification] Error: "
      << "input triangulation pointer is NULL." << endl;
    dMsg(cerr, msg.str(), fatalMsg);
    return -1;
  }
#endif
  
  triangulation_->setWrapper(this);
  persistenceSimplification_.setWrapper(this);
  persistenceSimplification_.setupTriangulation(triangulation_);
  Modified();
  
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangulation_->isEmpty()){
    stringstream msg;
    msg << "[ttkPersistenceSimplification] Error: "
      << "ttkTriangulation allocation problem." << endl;
    dMsg(cerr, msg.str(), fatalMsg);
    return -1;
  }
#endif
  
  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);
  
  // in the following, the target scalar field of the input is replaced in the 
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you 
  // should proceed in the same way.
  vtkDataArray *inputScalarField = NULL;
  
#ifndef TTK_ENABLE_KAMIKAZE
  if(!input->GetPointData()){
    stringstream msg;
    msg << "[ttkPersistenceSimplification] Error: "
      << "input has no point data." << endl;
    dMsg(cerr, msg.str(), fatalMsg);
    return -1;
  }
#endif
  
  if(ScalarField.length()){
    inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
  }
  else{
    inputScalarField = input->GetPointData()->GetArray(0);
  }
  
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalarField){
    stringstream msg;
    msg << "[ttkPersistenceSimplification] Error: "
      << "input scalar field pointer is null." << endl;
    dMsg(cerr, msg.str(), fatalMsg);
    return -2;
  }
#endif
  
  // Set offsets
  const SimplexId numberOfVertices=input->GetNumberOfPoints();
  inputOffsets_ = NULL;
  
  if(ForceInputOffsetScalarField and InputOffsetScalarFieldName.length()){
    inputOffsets_=input->GetPointData()->GetArray(
      InputOffsetScalarFieldName.data());
  }
  
  if(tmpOffsets_){
    tmpOffsets_->Delete();
    tmpOffsets_=nullptr;
  }
  
  if(!tmpOffsets_){
    
    tmpOffsets_=ttkSimplexIdTypeArray::New();
    tmpOffsets_->SetNumberOfComponents(1);
    tmpOffsets_->SetNumberOfTuples(numberOfVertices);
    tmpOffsets_->SetName(ttk::OffsetScalarFieldName);
    for(SimplexId i=0; i<numberOfVertices; ++i)
      tmpOffsets_->SetTuple1(i,i);
  }
  if(!inputOffsets_)
    inputOffsets_ = tmpOffsets_;

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputOffsets_){
    stringstream msg;
    msg << "[ttkPersistenceSimplification] Error: "
      << "wrong input offset scalar field." << endl;
    dMsg(cerr, msg.str(), fatalMsg);
    return -1;
  }
#endif
  
#ifndef TTK_ENABLE_KAMIKAZE
  if(inputOffsets_->GetDataType()!=VTK_INT 
    and inputOffsets_->GetDataType()!=VTK_ID_TYPE){
    stringstream msg;
    msg << "[ttkPersistenceSimplification] Error: "
      << "input offset field type not supported." << endl;
    dMsg(cerr, msg.str(), fatalMsg);
    return -1;
  }
#endif
  
  // allocate the memory for the output scalar field
  if(!outputScalarField_){
    switch(inputScalarField->GetDataType()){
      
      case VTK_CHAR:
        outputScalarField_ = vtkCharArray::New();
        break;
        
      case VTK_DOUBLE:
        outputScalarField_ = vtkDoubleArray::New();
        break;

      case VTK_FLOAT:
        outputScalarField_ = vtkFloatArray::New();
        break;
       
      case VTK_INT:
        outputScalarField_ = vtkIntArray::New();
        break;

      case VTK_ID_TYPE:
        outputScalarField_ = vtkIdTypeArray::New();
        break;
        
      case VTK_SHORT:
        outputScalarField_=vtkShortArray::New();
        break;

      case VTK_UNSIGNED_CHAR:
        outputScalarField_=vtkUnsignedCharArray::New();
        break;
        
      case VTK_UNSIGNED_SHORT:
        outputScalarField_=vtkUnsignedShortArray::New();
        break;
        
      default:
        {
          stringstream msg;
          msg << "[ttkPersistenceSimplification] "
            << "Unsupported data type :(" << endl;
          dMsg(cerr, msg.str(), fatalMsg);
        }
    }
  }
  if(outputScalarField_){
    outputScalarField_->SetNumberOfTuples(input->GetNumberOfPoints());
    outputScalarField_->SetName(inputScalarField->GetName());
  }
  #ifndef TTK_ENABLE_KAMIKAZE
  else{
    stringstream msg;
    msg << "[ttkPersistenceSimplification] Error: "
      << "vtkDataArray allocation problem." << endl;
    dMsg(cerr, msg.str(), fatalMsg);
    return -9;
  }
  #endif
  
  if(OutputOffsetScalarFieldName.length()<=0)
    OutputOffsetScalarFieldName=ttk::OffsetScalarFieldName;
  

  // get offsets
  vtkSmartPointer<ttkSimplexIdTypeArray> 
    outputOffsets=vtkSmartPointer<ttkSimplexIdTypeArray>::New();
    
  if(outputOffsets){
    outputOffsets->SetNumberOfComponents(1);
    outputOffsets->SetNumberOfTuples(numberOfVertices);
    outputOffsets->SetName(OutputOffsetScalarFieldName.data());
  }
  
  // calling the executing package
  persistenceSimplification_.setInputDataPointer(
    inputScalarField->GetVoidPointer(0));
  persistenceSimplification_.setOutputDataPointer(
    outputScalarField_->GetVoidPointer(0));
  persistenceSimplification_.setInputOffsetDataPointer(
    inputOffsets_->GetVoidPointer(0));
  persistenceSimplification_.setOutputOffsetDataPointer(
    outputOffsets->GetVoidPointer(0));

  // Pass Paraview Thresholding options to the base code
  persistenceSimplification_.setDistinctOption(DistinctOption);
  persistenceSimplification_.setAutoOption(AutoOption);
  persistenceSimplification_.setCountMinOption(CountMinOption);
  persistenceSimplification_.setCountMaxOption(CountMaxOption);
  persistenceSimplification_.setCountAllOption(CountAllOption);
  persistenceSimplification_.setCountMinArg(CountMinArg);
  persistenceSimplification_.setCountMaxArg(CountMaxArg);
  persistenceSimplification_.setCountAllArg(CountAllArg);
  persistenceSimplification_.setPersistenceMinArg(PersistenceMinArg);
  persistenceSimplification_.setPersistenceMaxArg(PersistenceMaxArg);
  persistenceSimplification_.setPersistenceAllArg(PersistenceAllArg);
  persistenceSimplification_.setUseMinOption(UseMinOption);
  persistenceSimplification_.setUseMaxOption(UseMaxOption);
  

  int ret = 0;
  switch(inputScalarField->GetDataType()){
    ttkTemplateMacro({
      if(inputOffsets_->GetDataType()==VTK_INT)
        ret=persistenceSimplification_.execute<VTK_TT TTK_COMMA int>();
      if(inputOffsets_->GetDataType()==VTK_ID_TYPE)
        ret=persistenceSimplification_.execute<VTK_TT TTK_COMMA vtkIdType>();
    });
  }
  #ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in baseCode
  if(ret){
    stringstream msg;
    msg << "[ttkPersistenceSimplification] PersistenceSimplification.execute() "
      << "error code : " << ret << endl;
    dMsg(cerr, msg.str(), fatalMsg);
    return -12;
  }
#endif
  
  // on the output, replace the field array by a pointer to its processed
  // version
  if(ScalarField.length()){
    output->GetPointData()->RemoveArray(ScalarField.data());
  }
  else{
    output->GetPointData()->RemoveArray(0);
  }
  output->GetPointData()->AddArray(outputScalarField_);
  output->GetPointData()->AddArray(outputOffsets);
  
  {
    stringstream msg;
    msg << "[ttkPersistenceSimplification] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}
