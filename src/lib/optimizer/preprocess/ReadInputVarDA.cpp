#include "VariationalDataAssimilation.h"

void VarDA::readInputVarDA()
{
  string str, base_label, label, method;

  base_label = "/VarDA";

  label = base_label + "/loopMAX";
  if(!tp.getInspectedValue(label, optMaxItr)){
    cout << label<< " is not set" << endl;
    exit(0);
  }

  label = base_label + "/output_itr";
  if(!tp.getInspectedValue(label, output_itr)){
    cout << label<< " is not set" << endl;
    exit(0);
  }

  label = base_label + "/alpha_cf";
  if(!tp.getInspectedValue(label, alpha_cf)){
    cout << label<< " is not set" << endl;
    exit(0);
  }

  label = base_label + "/beta1_cf";
  if(!tp.getInspectedValue(label, beta1_cf)){
    cout << label<< " is not set" << endl;
    exit(0);
  }

  label = base_label + "/beta2_cf";
  if(!tp.getInspectedValue(label, beta2_cf)){
    cout << label<< " is not set" << endl;
    exit(0);
  }

  label = base_label + "/controlBoundary";
  if(!tp.getInspectedValue(label, controlBoundary)){
    cout << label<< " is not set" << endl;
    exit(0);
  }
  if(    (controlBoundary != "top")   && (controlBoundary != "bottom")
      && (controlBoundary != "left")  && (controlBoundary != "right") 
      && (controlBoundary != "front") && (controlBoundary != "back"))
   {
    PetscPrintf(MPI_COMM_WORLD, "\n controlBoundary spelled wrong \n");
   }

  label = base_label + "/estimatedVelocity";
  if ( !tp.getInspectedValue(label, method)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  if(method == "CENTER"){
    estVel = EstimatedVelocity::CENTER;
  }else if(method == "AVERAGE"){
    estVel = EstimatedVelocity::AVERAGE;
  }else{
    cout << label << " is not set" << endl;
    exit(0);
  }

  return;
}
