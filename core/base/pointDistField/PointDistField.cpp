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

