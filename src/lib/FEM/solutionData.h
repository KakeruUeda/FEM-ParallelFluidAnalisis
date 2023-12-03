#include <iostream>
#include <limits.h>
#include <float.h>
#include <vector>
#include <assert.h>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <set>
#include <fstream>
#include <string.h>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>

#include <numeric>
#include <iterator>
#include <functional>

#include "eigen.h"
using namespace std;


class  SolutionData
{
  public:
    int  vecsize;

    VectorXd  solnDot, solnDotPrev, solnDotCur, solnExtrap;
    VectorXd  soln, solnPrev, solnPrev2, solnPrev3, solnPrev4, solnCur, solnInit, solnApplied;
    VectorXd  td;

    vector<int>  nodeMapPrev, nodeMap;

    SolutionData();
    ~SolutionData(){}


    void  initialise(int size1);

    void  setZero();
};