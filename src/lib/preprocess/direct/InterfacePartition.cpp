#include "FEM.h"

void FEM::interfacePartition()
{
  vector<vector<double>> x_current(numOfNodeInElm,vector<double>(3,0e0));
  vector<double> sdf_current(numOfNodeInElm,0e0);
  
  vector<int> sub_cross_num(numOfElmGlobalFluid,0);

  for(int ic=0; ic<numOfElmGlobalFluid; ic++)
  {  
    if(elmFluid[ic]->getSubdomainId() == myId){
      if(phiFluid[ic]>0.999){
        continue;
      }
      
      for(int i=0;i<numOfNodeInElm;i++){
        sdf_current[i] = sdfFluid[elmFluid[ic]->nodeNumsPrevFluid[i]];
        for(int j=0;j<3;j++){
          x_current[i][j] = xFluid[elmFluid[ic]->nodeNumsPrevFluid[i]][j];
        }
      }  
      vector<vector<double>> sub_x_tmp;

      line_serch(sub_cross_num, sub_x_tmp, sdf_current, x_current, ic);
      double size = sub_x_tmp.size();
      
      if(size == 3) tetraPartition1(sub_x_tmp, sdf_current, x_current, ic);
      
      if(size == 4){
        int check = 0;
        for(int i=0; i<numOfNodeInElm; i++){
          if(sdf_current[i] < 0) check++;
        }
        if(check == 2){
          //tetraPartition2(sub_x_tmp, sdf_current, x_current, ic);
        }else if(check == 4 || check == 6){
          //tetraPartition3(sub_x_tmp, sdf_current, x_current, ic);  
        }else{
          cout << "error when size = 4" << endl;
          cout << "check = " << check << endl;
          exit(1);
        }
      }

    }
  }

  if(myId == 4){
    
    ofstream check1("check_cross_point.dat");
    for(int ic=0; ic<numOfElmGlobalFluid; ic++){
      if(phiFluid[ic]>0.999) continue;
      if(elmFluid[ic]->getSubdomainId() == 4){
        check1 << "ic = " << ic << " numOfCrossPoint = " << sub_cross_num[ic] << endl;
      }
    }
    check1.close();
  }

 if(myId == 4){
    
    ofstream check2("check_x.dat");
    for(int ic=0; ic<numOfElmGlobalFluid; ic++){
      if(phiFluid[ic]>0.999) continue;
      if(elmFluid[ic]->getSubdomainId() == 4){
        if(sub_cross_num[ic] == 3){
          check2 << "*** ic = " << ic << " ****" << endl;
          for(int i=0; i<11; i++){
            for(int j=0; j<4; j++){
              for(int k=0; k<3; k++){
                check2 << elmFluid[ic]->sub_x[i][j][k] << " " ;
              }
              check2 << endl;
            }
            check2 << endl; 
          }
        }
      }
    }
    check2.close();
  }

}


void FEM::tetraPartition1(vector<vector<double>> &sub_x_tmp, vector<double> &sdf_current, vector<vector<double>> &x_current, const int ic)
{
  vector<double> base_node(3);
  vector<double> isolate_node(3);
  
  int minus = 0; 
  int plus = 0;
  
  for(int i=0; i<numOfNodeInElm; i++){
    if(sdf_current[i] < 0){
      minus++;
    }
    if(sdf_current[i] >= 0){
      plus++;
    }
  }

  if(minus == 1){
    for(int i=0; i<numOfNodeInElm; i++){
      if(sdf_current[i] < 0){
        for(int k=0; k<3; k++){
          isolate_node[k] = x_current[i][k];
        }
      }
    }
  }else if(plus == 1){
    for(int i=0; i<numOfNodeInElm; i++){
      if(sdf_current[i] >= 0){
        for(int k=0; k<3; k++){
          isolate_node[k] = x_current[i][k];
        }
      }
    }
  }else{
    cout << " couldn't find isolated node" << endl;
    exit(1);
  }


  double max_length = 1e-12;
  double length;
  for(int i=0; i<numOfNodeInElm; i++){
    length = 0;
    for(int k=0; k<3; k++){
      length += pow(x_current[i][k] - isolate_node[k],2e0);
    }
    length = sqrt(length);
    
    if(length > max_length){
      for(int k=0; k<3; k++){
        base_node[k] = x_current[i][k];
      }
    }
  }
  
  for(int i=0; i<3; i++){
   
    int n1, n2, n3;
    if(i==0){n1=0; n2=1;} 
    if(i==1){n1=1; n2=2;} 
    if(i==2){n1=2; n2=0;} 

    int coordinate;
    for(int k=0; k<3; k++){
      double diff = fabs(sub_x_tmp[n1][k]- sub_x_tmp[n2][k]);
      if(diff < 1e-8){
        coordinate = k;
      }
    }

    int tmp=0;
    vector<vector<double>> surface_node(3,vector<double>(3,0e0));
    
    for(int j=0; j<numOfNodeInElm; j++){
      
      if(minus == 1){
        if(sdf_current[j] < 0) continue;
      }else if(plus == 1){
        if(sdf_current[j] >= 0) continue;
      }
      
      double diff = fabs(x_current[j][coordinate] - sub_x_tmp[n1][coordinate]);
      if(diff < 1e-8){
        for(int k=0; k<3; k++){
          surface_node[tmp][k] = x_current[j][k];
        }
        tmp++;
      }
    }
    vector<double> base_surface_node(3);
    vector<double> adjacent_node1(3);
    vector<double> adjacent_node2(3);

    double max_surface_length = 1e-12;
    double surface_length;

    double min_surface_length_ad1 = 1e12;
    double min_surface_length_ad2 = 1e12;
    double surface_length_ad1;
    double surface_length_ad2;

    for(int j=0; j<3; j++){
      surface_length = 0;
      for(int k=0; k<3; k++){
        surface_length += pow(surface_node[j][k] - isolate_node[k],2e0);
      }
      surface_length = sqrt(surface_length);
      
      if(surface_length > max_surface_length){
        for(int k=0; k<3; k++){
          base_surface_node[k] = x_current[i][k];
        }
      }

     surface_length_ad1 = 0;
     surface_length_ad2 = 0;
     for(int k=0; k<3; k++){
       surface_length_ad1 += pow(surface_node[j][k] - sub_x_tmp[n1][k],2e0);
       surface_length_ad2 += pow(surface_node[j][k] - sub_x_tmp[n2][k],2e0);
     }
     surface_length_ad1 = sqrt(surface_length_ad1);
     surface_length_ad2 = sqrt(surface_length_ad2);
     
     if(surface_length_ad1 < min_surface_length_ad1){
       for(int k=0; k<3; k++){
         adjacent_node1[k] = surface_node[j][k];
       }
     }

     if(surface_length_ad2 < min_surface_length_ad2){
       for(int k=0; k<3; k++){
         adjacent_node2[k] = surface_node[j][k];
       }
     }
   }

    for(int k=0; k<3; k++){

     elmFluid[ic]->setNumOfSubdomain(11);
     elmFluid[ic]->sub_x.resize(11, vector<vector<double>>(4,vector<double>(3,0e0)));
    
     elmFluid[ic]->sub_x[i*3][0][k] = sub_x_tmp[n1][k];
     elmFluid[ic]->sub_x[i*3][1][k] = base_surface_node[k];
     elmFluid[ic]->sub_x[i*3][2][k] = adjacent_node1[k];
     elmFluid[ic]->sub_x[i*3][3][k] = base_node[k];

     elmFluid[ic]->sub_x[i*3+1][0][k] = sub_x_tmp[n2][k];
     elmFluid[ic]->sub_x[i*3+1][1][k] = adjacent_node1[k];
     elmFluid[ic]->sub_x[i*3+1][2][k] = base_surface_node[k];
     elmFluid[ic]->sub_x[i*3+1][3][k] = base_node[k];

     elmFluid[ic]->sub_x[i*3+2][0][k] = sub_x_tmp[n1][k];
     elmFluid[ic]->sub_x[i*3+2][1][k] = sub_x_tmp[n2][k];
     elmFluid[ic]->sub_x[i*3+2][2][k] = base_surface_node[k];
     elmFluid[ic]->sub_x[i*3+2][3][k] = base_node[k];
    }

  }

   for(int k=0; k<3; k++){
     elmFluid[ic]->sub_x[9][0][k] = sub_x_tmp[0][k];
     elmFluid[ic]->sub_x[9][1][k] = sub_x_tmp[1][k];
     elmFluid[ic]->sub_x[9][2][k] = sub_x_tmp[2][k];
     elmFluid[ic]->sub_x[9][3][k] = base_node[k];

     elmFluid[ic]->sub_x[10][0][k] = sub_x_tmp[0][k];
     elmFluid[ic]->sub_x[10][1][k] = sub_x_tmp[1][k];
     elmFluid[ic]->sub_x[10][2][k] = sub_x_tmp[2][k];
     elmFluid[ic]->sub_x[10][3][k] = isolate_node[k];
   }
}

/*
void FEM::tetraPartition2(vector<vector<double>> &sub_x_tmp, vector<double> &sdf_current, vector<vector<double>> &x_current, const int ic)
{

  vector<double> adjacent_node1(3);
  vector<double> adjacent_node2(3);
  int a, b;
  double max_length = 1e-12;
  double length;
  
  for(int i=1; i<numOfNodeInElm; i++){  
    if(sdf_current[i] < 0) continue;
    length = 0;
    for(int k=0; k<3; k++){
      length += pow(sub_x_tmp[i][k] - sub_x_tmp[0][k],2e0);
    }
    length = sqrt(length);
    
    if(length > max_length){
      for(int k=0; k<3; k++){
        adjacent_node1[k] = sub_x_tmp[i][k];
        a = i;
      }
    }
  }

  double min_length = 1e12;

  for(int i=1; i<numOfNodeInElm; i++){  
    if(sdf_current[i] < 0) continue;
    length = 0;
    for(int k=0; k<3; k++){
      length += pow(sub_x_tmp[i][k] - sub_x_tmp[0][k],2e0);
    }
    length = sqrt(length);
    
    if(length < min_length){
      for(int k=0; k<3; k++){
        adjacent_node2[k] = sub_x_tmp[i][k];
        b = i;
      }
    }
  }

  for(int i=1; i<numOfNodeinElm; i++){
    if((i != a) && (i != b)){
      for(int k=0; k<3; k++){
        adjacent_node3[k] = sub_x_tmp[i][k];
      }
    }
  }

  vector<double> solid_node1(3);
  vector<double> solid_node2(3);









  
  vector<double> base_node(3);
  vector<double> isolate_node(3);
  
  int minus = 0; 
  int plus = 0;
  
  for(int i=0; i<numOfNodeInElm; i++){
    if(sdf_current[i] < 0){
      minus++;
    }
    if(sdf_current[i] >= 0){
      plus++;
    }
  }

  if(minus == 1){
    for(int i=0; i<numOfNodeInElm; i++){
      if(sdf_current[i] < 0){
        for(int k=0; k<3; k++){
          isolate_node[k] = x_current[i][k];
        }
      }
    }
  }else if(plus == 1){
    for(int i=0; i<numOfNodeInElm; i++){
      if(sdf_current[i] >= 0){
        for(int k=0; k<3; k++){
          isolate_node[k] = x_current[i][k];
        }
      }
    }
  }else{
    cout << " couldn't find isolated node" << endl;
    exit(1);
  }


  double max_length = 1e-12;
  double length;
  for(int i=0; i<numOfNodeInElm; i++){
    length = 0;
    for(int k=0; k<3; k++){
      length += pow(x_current[i][k] - isolate_node[k],2e0);
    }
    length = sqrt(length);
    
    if(length > max_length){
      for(int k=0; k<3; k++){
        base_node[k] = x_current[i][k];
      }
    }
  }
  
  for(int i=0; i<3; i++){
   
    int n1, n2, n3;
    if(i==0){n1=0; n2=1;} 
    if(i==1){n1=1; n2=2;} 
    if(i==2){n1=2; n2=0;} 

    int coordinate;
    for(int k=0; k<3; k++){
      double diff = fabs(sub_x_tmp[n1][k]- sub_x_tmp[n2][k]);
      if(diff < 1e-8){
        coordinate = k;
      }
    }

    int tmp=0;
    vector<vector<double>> surface_node(3,vector<double>(3,0e0));
    
    for(int j=0; j<numOfNodeInElm; j++){
      
      if(minus == 1){
        if(sdf_current[j] < 0) continue;
      }else if(plus == 1){
        if(sdf_current[j] >= 0) continue;
      }
      
      double diff = fabs(x_current[j][coordinate] - sub_x_tmp[n1][coordinate]);
      if(diff < 1e-8){
        for(int k=0; k<3; k++){
          surface_node[tmp][k] = x_current[j][k];
        }
        tmp++;
      }
    }
    vector<double> base_surface_node(3);
    vector<double> adjacent_node1(3);
    vector<double> adjacent_node2(3);

    double max_surface_length = 1e-12;
    double surface_length;

    double min_surface_length_ad1 = 1e12;
    double min_surface_length_ad2 = 1e12;
    double surface_length_ad1;
    double surface_length_ad2;

    for(int j=0; j<3; j++){
      surface_length = 0;
      for(int k=0; k<3; k++){
        surface_length += pow(surface_node[j][k] - isolate_node[k],2e0);
      }
      surface_length = sqrt(surface_length);
      
      if(surface_length > max_surface_length){
        for(int k=0; k<3; k++){
          base_surface_node[k] = x_current[i][k];
        }
      }

     surface_length_ad1 = 0;
     surface_length_ad2 = 0;
     for(int k=0; k<3; k++){
       surface_length_ad1 += pow(surface_node[j][k] - sub_x_tmp[n1][k],2e0);
       surface_length_ad2 += pow(surface_node[j][k] - sub_x_tmp[n2][k],2e0);
     }
     surface_length_ad1 = sqrt(surface_length_ad1);
     surface_length_ad2 = sqrt(surface_length_ad2);
     
     if(surface_length_ad1 < min_surface_length_ad1){
       for(int k=0; k<3; k++){
         adjacent_node1[k] = surface_node[j][k];
       }
     }

     if(surface_length_ad2 < min_surface_length_ad2){
       for(int k=0; k<3; k++){
         adjacent_node2[k] = surface_node[j][k];
       }
     }
   }

    for(int k=0; k<3; k++){

     elmFluid[ic]->setNumOfSubdomain(11);
     elmFluid[ic]->sub_x.resize(11, vector<vector<double>>(4,vector<double>(3,0e0)));
    
     elmFluid[ic]->sub_x[i*3][0][k] = sub_x_tmp[n1][k];
     elmFluid[ic]->sub_x[i*3][1][k] = base_surface_node[k];
     elmFluid[ic]->sub_x[i*3][2][k] = adjacent_node1[k];
     elmFluid[ic]->sub_x[i*3][3][k] = base_node[k];

     elmFluid[ic]->sub_x[i*3+1][0][k] = sub_x_tmp[n2][k];
     elmFluid[ic]->sub_x[i*3+1][1][k] = adjacent_node1[k];
     elmFluid[ic]->sub_x[i*3+1][2][k] = base_surface_node[k];
     elmFluid[ic]->sub_x[i*3+1][3][k] = base_node[k];

     elmFluid[ic]->sub_x[i*3+2][0][k] = sub_x_tmp[n1][k];
     elmFluid[ic]->sub_x[i*3+2][1][k] = sub_x_tmp[n2][k];
     elmFluid[ic]->sub_x[i*3+2][2][k] = base_surface_node[k];
     elmFluid[ic]->sub_x[i*3+2][3][k] = base_node[k];
    }

  }

   for(int k=0; k<3; k++){
     elmFluid[ic]->sub_x[9][0][k] = sub_x_tmp[0][k];
     elmFluid[ic]->sub_x[9][1][k] = sub_x_tmp[1][k];
     elmFluid[ic]->sub_x[9][2][k] = sub_x_tmp[2][k];
     elmFluid[ic]->sub_x[9][3][k] = base_node[k];

     elmFluid[ic]->sub_x[10][0][k] = sub_x_tmp[0][k];
     elmFluid[ic]->sub_x[10][1][k] = sub_x_tmp[1][k];
     elmFluid[ic]->sub_x[10][2][k] = sub_x_tmp[2][k];
     elmFluid[ic]->sub_x[10][3][k] = isolate_node[k];
   }
}
*/