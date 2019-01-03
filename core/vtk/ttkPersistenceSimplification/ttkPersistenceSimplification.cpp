#include                  <ttkPersistenceSimplification.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPersistenceSimplification)

int ttkPersistenceSimplification::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  Memory m;
  
  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];
  
  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);
 
  if(!triangulation)
    return -1;
  
  triangulation->setWrapper(this);
  persistenceSimplification_.setupTriangulation(triangulation);
  persistenceSimplification_.setWrapper(this);
 
  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);
  
  // in the following, the target scalar field of the input is replaced in the 
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you 
  // should proceed in the same way.
  vtkDataArray *inputScalarField = NULL;
  
  if(ScalarField.length()){
    inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
  }
  else{
    inputScalarField = input->GetPointData()->GetArray(0);
  }
  
  if(!inputScalarField)
    return -2;

  // Set offsets
  const SimplexId numberOfVertices=input->GetNumberOfPoints();

  inputOffsets_=ttkSimplexIdTypeArray::New();
  inputOffsets_->SetNumberOfComponents(1);
  inputOffsets_->SetNumberOfTuples(numberOfVertices);
  inputOffsets_->SetName(ttk::OffsetScalarFieldName);
  for(SimplexId i=0; i<numberOfVertices; ++i)
    inputOffsets_->SetTuple1(i,i);
  
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
        
      stringstream msg;
      msg << "[ttkPersistenceSimplification] Unsupported data type :(" << endl;
      dMsg(cerr, msg.str(), fatalMsg);
    }
  }
  outputScalarField_->SetNumberOfTuples(input->GetNumberOfPoints());
  outputScalarField_->SetName(inputScalarField->GetName());
  

  // get offsets
  inputOffsets_=ttkSimplexIdTypeArray::New();
  inputOffsets_->SetNumberOfComponents(1);
  inputOffsets_->SetNumberOfTuples(numberOfVertices);
  inputOffsets_->SetName(ttk::OffsetScalarFieldName);
  for(SimplexId i=0; i<numberOfVertices; ++i)
    inputOffsets_->SetTuple1(i,i);
  vtkSmartPointer<ttkSimplexIdTypeArray> outputOffsets=vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  if(outputOffsets){
    outputOffsets->SetNumberOfComponents(1);
    outputOffsets->SetNumberOfTuples(numberOfVertices);
    outputOffsets->SetName("OutputOffsetScalarField");
  }
  
  // calling the executing package
  persistenceSimplification_.setInputDataPointer(inputScalarField->GetVoidPointer(0));
  persistenceSimplification_.setOutputDataPointer(outputScalarField_->GetVoidPointer(0));
  persistenceSimplification_.setInputOffsetDataPointer(inputOffsets_->GetVoidPointer(0));
  persistenceSimplification_.setOutputOffsetDataPointer(outputOffsets->GetVoidPointer(0));

  // TODO_RC

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
  

  switch(inputScalarField->GetDataType()){
    ttkTemplateMacro(persistenceSimplification_.execute<VTK_TT TTK_COMMA int>());
  }
  
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
