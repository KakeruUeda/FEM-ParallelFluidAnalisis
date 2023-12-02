#include "FEM.h"
using namespace std;

void FEM::initialize()
{

  int  size1, size2, row, col;
  int  tempDOF, domTemp, ind;
  int  r, c, r1, c1, count=0, count1=0, count2=0, ii, jj, ee, dd, ind1, ind2, nsize;
  int  *tt, val1, val2, nnz, nnz_max_row, n1, n2, kk, e1, e2, ll, aa, bb, nn, dof;
  int  side, start1, start2, nr1, nr2, count_diag, count_offdiag, tempInt;
  
  readInput(); 
  setDomain();
  setBoundary();
  //setImage();
  node_map_get_old.resize(numOfNodeGlobal, 0);
  node_map_get_new.resize(numOfNodeGlobal, 0);
  node_proc_id.resize(numOfNodeGlobal, 0);

    // set sizes of some data arrays
  vector<bool>  vecBoolTempFalse(ndof, false);
	NodeTypeOld.resize(numOfNodeGlobal, vecBoolTempFalse);
  NodeTypeNew.resize(numOfNodeGlobal, vecBoolTempFalse);

  vector<int>  vecIntTempM1(ndof, -1);
  NodeDofArrayOld.resize(numOfNodeGlobal, vecIntTempM1);
  NodeDofArrayOld_withoutBd.resize(numOfNodeGlobal, vecIntTempM1);
  NodeDofArrayNew.resize(numOfNodeGlobal, vecIntTempM1);
  
  for(ii=0; ii<numOfBdNode; ii++){
    NodeTypeOld[DirichletBCs[ii][0]][DirichletBCs[ii][1]] = true;
  }

  ntotdofs_global = 0;
  for(int ii=0;ii<numOfNodeGlobal;++ii){
    for(int jj=0;jj<ndof;++jj){
      ////erase/////
      //if(!NodeTypeOld[ii][jj]){
        NodeDofArrayOld_withoutBd[ii][jj] = ntotdofs_global;
        NodeDofArrayOld[ii][jj] = ntotdofs_global++;
      //}

      if(NodeTypeOld[ii][jj]){
        NodeDofArrayOld[ii][jj] = -1;
      }
      ////erase////
    }
  }

  if(this_mpi_proc == 0)
  {
    cout << " Mesh statistics .....\n" << endl;
    cout << " numOfElmGlobal   = " << '\t' << numOfElmGlobal << endl;
    cout << " numOfNodeGlobal  = " << '\t' << numOfNodeGlobal  << endl;
    cout << " numOfNodeInElm   = " << '\t' << numOfNodeInElm << endl;
    cout << " ndof             = " << '\t' << ndof << endl;
    cout << " ntotdofs_local   = " << '\t' << ntotdofs_local  << endl;
    cout << " ntotdofs_global  = " << '\t' << ntotdofs_global << endl;
    cout << " n_mpi_procs      = " << '\t' << n_mpi_procs << endl;
  }

  PetscPrintf(MPI_COMM_WORLD, "\n\n    Before partitionMesh ... \n\n");
  divideMesh();
  PetscPrintf(MPI_COMM_WORLD, "\n    After partitionMesh ... \n\n"); 
  
  errpetsc = MPI_Barrier(MPI_COMM_WORLD);
  
  PetscPrintf(MPI_COMM_WORLD, "\n\n    Before prepareDataForParallel ... \n\n");
  prepareForParallel();
  PetscPrintf(MPI_COMM_WORLD, "\n    After prepareDataForParallel ... \n\n");
    
  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  SolnData.node_map_get_old = node_map_get_old;
  SolnData.node_map_get_new = node_map_get_new;

  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  for(ee=0; ee<numOfElmGlobal; ++ee)
  {
    if(elm[ee]->getSubdomainId() == this_mpi_proc)
    {
      nsize = ndof*numOfNodeInElm;

      elm[ee]->prepareElemData(nsize);

      numOfNodeInElm = elm[ee]->nodeNums.size();

      elm[ee]->forAssyVec.resize(nsize);
      elm[ee]->forAssyVec_withoutBd.resize(nsize);

      for(ii=0; ii<numOfNodeInElm; ++ii)
      {
        ind = ndof*ii;

        kk = elm[ee]->nodeNums[ii];

        for(jj=0;jj<ndof;++jj)
        {
          elm[ee]->forAssyVec[ind+jj] = NodeDofArrayNew[kk][jj];
          elm[ee]->forAssyVec_withoutBd[ind+jj] = NodeDofArrayOld_withoutBd[kk][jj];
        }
      }
    }
  }   
  /*
  if(this_mpi_proc == 1){
    ofstream Assy("forAssyVec1.dat"); 
    for(ee=0; ee<numOfElmGlobal; ee++){
      int size = elm[ee]->forAssyVec.size();
      for(ii=0; ii<size; ii++){
        if(ii != 0 && ii%4 == 0) Assy << " ";
        Assy << elm[ee]->forAssyVec[ii] << " ";
      } 
      Assy << endl;
    }
    Assy.close();
  }
  */


  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  assyForSoln.resize(ntotdofs_global);
  ind = 0;
  for(ii=0;ii<numOfNodeGlobal;++ii)
  {
    for(jj=0;jj<ndof;++jj)
    {
      ///erase////
      //if(NodeDofArrayNew[ii][jj] != -1)
      //{
        assyForSoln[ind++] = ii*ndof+jj;
      //}
    }
  }

  PetscPrintf(MPI_COMM_WORLD, "\n\n Element DOF values initialised \n\n");
  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  cout << " Total DOF   = " << '\t' << ntotdofs_local << '\t' << ntotdofs_global << endl;

  //vector<vector<int> > DDconn;
  vector<set<int> > forAssyMat;
  set<int>::iterator it;

  forAssyMat.resize(ntotdofs_global);

  //for(ii=row_start; ii<=row_end; ii++)
  //forAssyMat[ii].reserve(500);
  for(ee=0; ee<numOfElmGlobal; ee++)
  {
    if(elm[ee]->getSubdomainId() == this_mpi_proc)
    {
      tt = &(elm[ee]->forAssyVec_withoutBd[0]);
      nsize = elm[ee]->forAssyVec_withoutBd.size();

      for(ii=0;ii<nsize;ii++)
      {
          r = tt[ii];

          if(r != -1)
          {
            if(r >= row_start && r <= row_end)
            {
            for(jj=0;jj<nsize;jj++)
              {
                if(tt[jj] != -1)
                {
                  //printf("ii.... %5d \t %5d \t %5d \t %5d \n",ii, jj, r, tt[jj]);
                  forAssyMat[r].insert(tt[jj]);
                }
              }
            }
          }
      }
    }
  }

  errpetsc = MPI_Barrier(MPI_COMM_WORLD);
  PetscPrintf(MPI_COMM_WORLD, "\n\n Preparing matrix pattern DONE \n\n");


  PetscInt  *diag_nnz, *offdiag_nnz;

  errpetsc  = PetscMalloc1(ntotdofs_local,  &diag_nnz);
  errpetsc  = PetscMalloc1(ntotdofs_local,  &offdiag_nnz);


  kk = 0;
  nnz_max_row = 0;
  for(ii=row_start; ii<=row_end; ii++)
  {
    size1 = forAssyMat[ii].size();

    nnz_max_row = max(nnz_max_row, size1);

    count_diag=0, count_offdiag=0;
    for(it=forAssyMat[ii].begin(); it!=forAssyMat[ii].end(); ++it)
    {
      tempInt = *it;

      if(tempInt >= row_start && tempInt <= row_end)
        count_diag++;
      else
        count_offdiag++;
    }

    //if(this_mpi_proc == 1){
      //cout << " count_diag ..." << ii << '\t' << count_diag << '\t' << count_offdiag << endl;
    //}

    diag_nnz[kk]    = count_diag;
    offdiag_nnz[kk] = count_offdiag;
    kk++;
  }


  errpetsc = MPI_Barrier(MPI_COMM_WORLD);
  //exit(1);
  PetscPrintf(MPI_COMM_WORLD, "\n\n Initialising petsc solver \n\n");

  solverPetsc->initialise(ntotdofs_local, ntotdofs_global, diag_nnz, offdiag_nnz);
  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

    //Create parallel matrix, specifying only its global dimensions.
    //When using MatCreate(), the matrix format can be specified at
    //runtime. Also, the parallel partitioning of the matrix is
    //determined by Petsc at runtime.
    //Performance tuning note: For problems of substantial size,
    //preallocation of matrix memory is crucial for attaining good
    //performance. See the matrix chapter of the users manual for details.

  PetscPrintf(MPI_COMM_WORLD, " Initialise the Matrix pattern \n", errpetsc);

  ind = numOfNodeInElm*ndof;
  ind = ind*ind;
  PetscScalar  Klocal[ind];
  for(ii=0; ii<ind; ii++)  Klocal[ii] = 0.0;
  
  MatrixXd Klocal_tmp(numOfNodeInElm*ndof,numOfNodeInElm*ndof);

  vector<int>  vecIntTemp;
  for(ee=0; ee<numOfElmGlobal; ee++)
  {  
    Klocal_tmp.setZero();
    if(elm[ee]->getSubdomainId() == this_mpi_proc)
    {
      size1 = elm[ee]->forAssyVec.size();
      vecIntTemp = elm[ee]->forAssyVec_withoutBd;
      for(ii=0; ii<size1; ii++){
        if(elm[ee]->forAssyVec[ii] == -1){
          Klocal[ii*numOfNodeInElm*ndof+ii] = 1.0;
          Klocal_tmp(ii,ii) = 1.0;
          //printf("%5d \t %5d \t %5d \t %5d \n",ee, ii, vecIntTemp[ii], this_mpi_proc);
          //cout << "e = " << ee << " ii = " << ii << endl;
          //cout << "i = " << ii*numOfNodeInElm*ndof+ii << " Klocal = " << Klocal[ii*numOfNodeInElm*ndof+ii] << endl;
        }
      }
      //MatrixXdRM Klocal2_tmp = Klocal_tmp;
      errpetsc = MatSetValues(solverPetsc->mtx, size1, &vecIntTemp[0], size1, &vecIntTemp[0], &Klocal_tmp(0,0), INSERT_VALUES);
    }
  }
  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  //MatAssemblyBegin(solverPetsc->mtx,MAT_FINAL_ASSEMBLY);
  //MatAssemblyEnd(solverPetsc->mtx,MAT_FINAL_ASSEMBLY);
  //MatView(solverPetsc->mtx,PETSC_VIEWER_STDOUT_WORLD);
  //exit(1);
  /*
  if(this_mpi_proc == 1){
    ofstream Kb("Klocal_before.dat"); 
    for(int ic=0; ic<ind; ic++){
      Kb << Klocal[ic]  << " ";
      if((ic+1)%32 == 0){
        Kb << endl;
      } 
    }
    Kb.close();
    //exit(1);
  }
  */
  /*
  if(this_mpi_proc == 0){
    //MatView(mtx, PETSC_VIEWER_STDOUT_WORLD);
    //MatView(A,viewer);
    ofstream Assy("forAssyVec0_withoutBd.dat"); 
    for(int ic=0; ic<numOfElmGlobal; ic++){
      int size = elm[ic]->forAssyVec_withoutBd.size();
      for(ii=0; ii<size; ii++){
        if(ii != 0 && ii%4 == 0) Assy << " ";
        Assy << elm[ic]->forAssyVec_withoutBd[ii] << " ";
      } 
      Assy << endl;
    }
    Assy.close();
  }
    
  if(this_mpi_proc == 1){
    ofstream Assy1("forAssyVec1_withoutBd.dat"); 
    for(int ic=0; ic<numOfElmGlobal; ic++){
      int size = elm[ic]->forAssyVec_withoutBd.size();
      for(ii=0; ii<size; ii++){
        if(ii != 0 && ii%4 == 0) Assy1 << " ";
        Assy1 << elm[ic]->forAssyVec_withoutBd[ii] << " ";
      } 
      Assy1 << endl;
    }
    Assy1.close();
  }
  */
  
  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  solverPetsc->currentStatus = PATTERN_OK;

  PetscFree(diag_nnz); 
  PetscFree(offdiag_nnz);  

  //TODO: use vector<T>().swap[x) to empty the vector and also deallocate the memory
  for(ii=0; ii<NodeDofArrayNew.size(); ii++)
  {
    NodeDofArrayOld[ii].clear();
    NodeDofArrayNew[ii].clear();
  }
  NodeDofArrayOld.clear();
  NodeDofArrayNew.clear();

}

void FEM::readInput()
{ 
  string str,base_label,label;

  int tmp[3];
  double dummy[3];

  base_label = "/Domain";
  label = base_label + "/nx";
  if ( !tp.getInspectedVector(label,tmp,3)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  nx=tmp[0];
  ny=tmp[1];
  nz=tmp[2];

  label = base_label + "/Lx";
  if ( !tp.getInspectedVector(label,dummy,3)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  Lx=dummy[0];
  Ly=dummy[1];
  Lz=dummy[2];

  /////param////
  base_label = "/Parameter";
  label = base_label + "/mu";
  if ( !tp.getInspectedValue(label, mu)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  
}

void FEM::setDomain()
{
  dx=Lx/nx;
  dy=Ly/ny;
  dz=Lz/nz; 
  numOfNodeGlobal=(nx+1)*(ny+1)*(nz+1);
  numOfElmGlobal=nx*ny*nz;
  numOfNodeInElm=8;
  x.resize(numOfNodeGlobal,vector<double>(3, 0));
  element.resize(numOfElmGlobal);

  for(int ic=0;ic<numOfElmGlobal;ic++){
    element[ic].resize(8);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  int tmp2=0;
  for(int k=0; k<nz+1; k++){
    for(int i=0; i<ny+1; i++){
      for(int j=0; j<nx+1; j++){
        x[tmp2][0]=j*dx;
        x[tmp2][1]=i*dy;
        x[tmp2][2]=k*dz; 
        tmp2++;
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);


  tmp2=0;
  for(int k=0;k<nz;k++){
    for(int i=0;i<ny;i++){
      for(int j=0;j<nx;j++){
          element[tmp2][0]=j   +i*(nx+1) +k*(nx+1)*(ny+1);
          element[tmp2][1]=j+1 +i*(nx+1) +k*(nx+1)*(ny+1);
          element[tmp2][2]=j+1 +(i+1)*(nx+1) +k*(nx+1)*(ny+1);
          element[tmp2][3]=j   +(i+1)*(nx+1) +k*(nx+1)*(ny+1);
          element[tmp2][4]=j   +i*(nx+1) +(k+1)*(nx+1)*(ny+1);
          element[tmp2][5]=j+1 +i*(nx+1) +(k+1)*(nx+1)*(ny+1);
          element[tmp2][6]=j+1 +(i+1)*(nx+1) +(k+1)*(nx+1)*(ny+1);
          element[tmp2][7]=j   +(i+1)*(nx+1) +(k+1)*(nx+1)*(ny+1);
          tmp2++;
      }
    }
  }
  prepare();
}

void FEM::prepare()
{
  elm = new ElementBaseFEM* [numOfElmGlobal];

  for(int ee=0;ee<numOfElmGlobal;++ee){
    elm[ee] = new ElementBaseFEM();
    elm[ee]->SolnData = &(SolnData);
  }

  if(solverPetsc != NULL) delete solverPetsc;
  solverPetsc = NULL;

  solverPetsc = (PetscSolver*) new PetscSolver;

  SolnData.initialise(numOfNodeGlobal*ndof);
}



void FEM::setBoundary(){

  string str,base_label,label;

  int tmp[3];
  double dummy[3];
  
  bd_u.resize(numOfNodeGlobal,vector<double>(3,0));
  bd_iu.resize(numOfNodeGlobal,vector<int>(3,0));

  bd_p.resize(numOfNodeGlobal,0);
  bd_ip.resize(numOfNodeGlobal,0);

  for(int i=0;i<numOfNodeGlobal;i++){
    for(int j=0;j<3;j++){
      bd_iu[i][j] = 1;
      bd_u[i][j] = 0;
    }
  }

  for(int i=0;i<numOfNodeGlobal;i++){
    bd_ip[i] = 1;
    bd_p[i] = 0;
  }

  string bdType;
  base_label = "/Boundary";

  label = base_label + "/left/type";
  if ( !tp.getInspectedValue(label,bdType)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  
  if(bdType=="v"){
    double value[3];
    label = base_label + "/left/value";
    if ( !tp.getInspectedVector(label,value,3)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp=0;
    for(int i=0; i<nz+1; i++){
      for(int j=0; j<nx+1; j++){
        for(int k=0; k<3; k++){
          bd_iu[(nx+1)*j+(nx+1)*(ny+1)*i][k] = 0;
          bd_u[(nx+1)*j+(nx+1)*(ny+1)*i][k] = value[k];
          //cout << "letvalue = " << value[k] << endl; 
        }
      }
    }
  }else if(bdType=="p"){
    double value;
    label = base_label + "/left/value";
    if ( !tp.getInspectedValue(label,value)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = 0;
    for(int i=0;i<nz+1;i++){
      for(int j=0; j<nx+1; j++){
        bd_ip[(nx+1)*j+(nx+1)*(ny+1)*i] = 0;
        bd_p[(nx+1)*j+(nx+1)*(ny+1)*i] = value;
      }
    }
  }else if(bdType=="free"){
  }else{
    cout << bdType << " is undefined. Exit..." << endl;
  }

  label = base_label + "/right/type";
  if ( !tp.getInspectedValue(label,bdType)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  
  if(bdType=="v"){
    double value[3];
    label = base_label + "/right/value";
    if ( !tp.getInspectedVector(label,value,3)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = nx;
    for(int i=0;i<nz+1;i++){
      for(int j=0; j<nx+1; j++){
        for(int k=0;k<3;k++){
          bd_iu[(nx+1)*j+(nx+1)*(ny+1)*i+tmp][k] = 0;
          bd_u[(nx+1)*j+(nx+1)*(ny+1)*i+tmp][k] = value[k];
        }
      }
    }
  }else if(bdType=="v_poiseuille"){
    int tmp = nx;
    for(int i=0;i<ny+1;i++){
      for(int j=0;j<2;j++){
        bd_iu[(nx+1)*i+tmp][j] = 0;
        //valueはinputImageInfoで設定
      }
    }
  }else if(bdType=="p"){
    double value;
    label = base_label + "/right/value";
    if ( !tp.getInspectedValue(label,value)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = 0;
    for(int i=0; i<nz+1;i++){
      for(int j=0; j<nx+1; j++){
        bd_ip[(nx+1)*j+(nx+1)*(ny+1)*i+tmp] = 0;
        bd_p[(nx+1)*j+(nx+1)*(ny+1)*i+tmp] = value;
      }
    }
  }else if(bdType=="free"){
  }else{
    cout << bdType << " is undefined. Exit..." << endl;
  }

  label = base_label + "/bottom/type";
  if ( !tp.getInspectedValue(label,bdType)){
    cout << label << " is not set" << endl;
    exit(0);
  }

  if(bdType=="v"){
    double value[3];
    label = base_label + "/bottom/value";
    if ( !tp.getInspectedVector(label,value,3)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = 0;
    for(int i=0;i<nz+1;i++){
      for(int j=0;j<nx+1;j++){
        for(int k=0;k<3;k++){
          bd_iu[i*(nx+1)*(ny+1)+tmp+j][k] = 0;
          bd_u[i*(nx+1)*(ny+1)+tmp+j][k] = value[k];
        }
      }
    }
  }else if(bdType=="p"){
    double value;
    label = base_label + "/bottom/value";
    if ( !tp.getInspectedValue(label,value)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = 0;
    for(int i=0;i<nz+1;i++){
      for(int j=0; j<nx+1; j++){
        bd_ip[i*(nx+1)*(ny+1)+tmp+j] = 0;
        bd_p[i*(nx+1)*(ny+1)+tmp+j] = value;
      }
    }
  }else if(bdType=="free"){
  }else{
    cout << bdType << " is undefined. Exit..." << endl;
  }

  label = base_label + "/top/type";
  if ( !tp.getInspectedValue(label,bdType)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  //ofstream top["top.dat");

  if(bdType=="v"){
    double value[3];
    label = base_label + "/top/value";
    if ( !tp.getInspectedVector(label,value,3)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = (nx+1)*ny;
    for(int i=0;i<nz+1;i++){
      for(int j=0; j<nx+1; j++){
        for(int k=0;k<3;k++){
          bd_iu[i*(nx+1)*(ny+1)+tmp+j][k] = 0;
          bd_u[i*(nx+1)*(ny+1)+tmp+j][k]  = value[k];
        }
      }
    }
  }else if(bdType=="p"){
    double value;
    label = base_label + "/top/value";
    if ( !tp.getInspectedValue(label,value)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = (nx+1)*ny;
    for(int i=0;i<nz+1;i++){
      for(int j=0; j<nx+1; j++){
        bd_ip[i*(nx+1)*(ny+1)+tmp+j] = 0;
        bd_p[i*(nx+1)*(ny+1)+tmp+j] = value;
      }
    }
  }else if(bdType=="free"){
  }else{
    cout << bdType << " is undefined. Exit..." << endl;
  }

  
  label = base_label + "/front/type";
  if ( !tp.getInspectedValue(label,bdType)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  
  if(bdType=="v"){
    double value[3];
    label = base_label + "/front/value";
    if ( !tp.getInspectedVector(label,value,3)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp=0;
    for(int i=0; i<ny+1; i++){
      for(int j=0; j<nx+1; j++){
        for(int k=0; k<3; k++){
          bd_iu[j+(nx+1)*i][k] = 0;
          bd_u[j+(nx+1)*i][k] = value[k];
        }
      }
    }
  }else if(bdType=="p"){
    double value;
    label = base_label + "/front/value";
    if ( !tp.getInspectedValue(label,value)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = 0;
    for(int i=0;i<ny+1;i++){
      for(int j=0; j<nx+1; j++){
        bd_ip[j+(nx+1)*i] = 0;
        bd_p[j+(nx+1)*i] = value;
      }
    }
  }else if(bdType=="free"){
  }else{
    cout << bdType << " is undefined. Exit..." << endl;
  }


  
  label = base_label + "/back/type";
  if ( !tp.getInspectedValue(label,bdType)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  
  if(bdType=="v"){
    double value[3];
    label = base_label + "/back/value";
    if ( !tp.getInspectedVector(label,value,3)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp=0;
    for(int i=0; i<ny+1; i++){
      for(int j=0; j<nx+1; j++){
        for(int k=0; k<3; k++){
          bd_iu[j+(nx+1)*i+(nx+1)*(ny+1)*nz][k] = 0;
          bd_u[j+(nx+1)*i+(nx+1)*(ny+1)*nz][k] = value[k];
        }
      }
    }
  }else if(bdType=="p"){
    double value;
    label = base_label + "/back/value";
    if ( !tp.getInspectedValue(label,value)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = 0;
    for(int i=0;i<ny+1;i++){
      for(int j=0; j<nx+1; j++){
        bd_ip[j+(nx+1)*i+(nx+1)*(ny+1)*nz] = 0;
        bd_p[j+(nx+1)*i+(nx+1)*(ny+1)*nz] = value;
      }
    }
  }else if(bdType=="free"){
  }else{
    cout << bdType << " is undefined. Exit..." << endl;
  }

  int count = 0;
  vector<double> vecDbTmp(3,0);
  

  for(int i=0;i<numOfNodeGlobal;i++){
    for(int k=0; k<3; k++){
      if(bd_iu[i][k] == 0){
        vecDbTmp[0] = i;
        vecDbTmp[1] = k;
        vecDbTmp[2] = bd_u[i][k];
        DirichletBCs.push_back(vecDbTmp);
        count++;
      }
    }
  }
  
  ///pressure///
  vecDbTmp[0] = 0;
  vecDbTmp[1] = 3;
  vecDbTmp[2] = 1e0;
  DirichletBCs.push_back(vecDbTmp);
  count++;


  numOfBdNode = count;

  if(this_mpi_proc == 0){
    ofstream test("test.dat");
   
    test << "ndim  " << 3 << endl; 
    test << "npElm  " << numOfNodeInElm << endl; 
    test << "nNode  " << numOfNodeGlobal << endl; 
    test << "nElem  " << numOfElmGlobal << endl; 
    test << "nDBC  " << numOfBdNode << endl; 
    test << "nFBC  " << 0 << endl; 
    test << "nOutNodes  " << 0 << endl; 
    
    test << "nodes" <<  endl; 
    for(int i=0; i<numOfNodeGlobal; i++){
      test << i+1 << " " <<  x[i][0] << " " << x[i][1] << " " << x[i][2] << endl;
    }

    test << "elements" << endl;
    for(int i=0; i<numOfElmGlobal; i++){
      test << i+1 << " " << 1 << " " << 1 << " " << 1 << " " << element[i][0]+1 << " " <<  element[i][1]+1 << " " <<  element[i][2]+1 << " " <<  element[i][3]+1 << " " <<   element[i][4]+1 << " " <<  element[i][5]+1 << " " <<  element[i][6]+1 << " " <<  element[i][7]+1 << endl;
    }
    test << "prescribed boundary conditions" << endl;
    for(int i=0; i<numOfBdNode; i++){
      test <<  DirichletBCs[i][0]+1 << " " <<  DirichletBCs[i][1]+1 << " " <<  DirichletBCs[i][2] << endl;
    }
    test.close();
    //exit(1);
  }



}

int FEM::divideMesh()
{

  vector<int>  elm_proc_id(numOfElmGlobal, 0);

  if(this_mpi_proc == 0)
  {
    int  ee, ii, jj, kk, n2;
    int  nparts = n_mpi_procs, subdomain=0;

    /////////////////////////////////////////////////////////////////////////////
    //
    // Partition the mesh. Here, METIS is used.
    // 
    /////////////////////////////////////////////////////////////////////////////

    PetscInt  *eptr, *eind;

    errpetsc  = PetscMalloc1(numOfElmGlobal+1,  &eptr);CHKERRQ(errpetsc);
    errpetsc  = PetscMalloc1(numOfElmGlobal*numOfNodeInElm,  &eind);CHKERRQ(errpetsc);

    vector<int>  vecTemp2;

    eptr[0] = 0;
    kk = 0;
    for(ee=0; ee<numOfElmGlobal; ee++)
    {
      eptr[ee+1] = (ee+1)*numOfNodeInElm;
      //vecTemp2 = elems[ee]->nodeNums ;
      vecTemp2 = element[ee];
      //printVector(vecTemp2);
      for(ii=0; ii<numOfNodeInElm; ii++)
        eind[kk+ii] = vecTemp2[ii] ;
      kk += numOfNodeInElm;
    }

    int  ncommon_nodes;
    //if(ndim == 2)
    //  ncommon_nodes = 2;    // 3-noded tria or 4-noded quad
    //else
    //{
      //if(numOfNodeInElm == 4)       // 4-noded tetra element
        //ncommon_nodes = 3;
      //else
        ncommon_nodes = 4;  // 8-noded hexa element
    //}

    idx_t objval;
    idx_t options[METIS_NOPTIONS];

    METIS_SetDefaultOptions(options);

    // Specifies the partitioning method.
    //options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;          // Multilevel recursive bisectioning.
    options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;        // Multilevel k-way partitioning.

    //options[METIS_OPTION_NSEPS] = 10;

    //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;     // Edge-cut minimization
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;     // Total communication volume minimization

    options[METIS_OPTION_NUMBERING] = 0;  // C-style numbering is assumed that starts from 0.

    cout << " Executing METIS subroutine " << endl;

    // METIS partition routine
    //int ret = METIS_PartMeshNodal(&nElem_global, &nNode_global, eptr, eind, NULL, NULL, &nparts, NULL, options, &objval, elem_proc_id, node_proc_id);
    int ret = METIS_PartMeshDual(&numOfElmGlobal, &numOfNodeGlobal, eptr, eind, NULL, NULL, &ncommon_nodes, &nparts, NULL, options, &objval, &elm_proc_id[0], &node_proc_id[0]);

    if(ret == METIS_OK)
      cout << " METIS partition routine successful "  << endl;
    else
      cout << " METIS partition routine FAILED "  << endl;

    errpetsc = PetscFree(eptr); CHKERRQ(errpetsc);
    errpetsc = PetscFree(eind); CHKERRQ(errpetsc);
    
    //ofstream out_nodeId("check_nodeProcId.dat");
    //ofstream out_elmId("check_elmProcId.dat");
    //for(ee=0; ee<numOfNodeGlobal; ee++)
    //{
    //  out_nodeId << ee <<  " node_proc_id[ee] = " << node_proc_id[ee] << endl;
    //}
    //for(ee=0; ee<numOfElmGlobal; ee++)
    //{
    //  out_elmId << ee << " elm_proc_id[ee] = " << elm_proc_id[ee] << endl;
    //}
    //out_nodeId.close();
    //out_elmId.close();

    /// check ///
    string vtiFile;
    vtiFile = "check_MeshPartition.vti";
    export_vti(vtiFile, node_proc_id, elm_proc_id);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  errpetsc = MPI_Bcast(&elm_proc_id[0], numOfElmGlobal, MPI_INT, 0, MPI_COMM_WORLD); CHKERRQ(errpetsc);
  MPI_Barrier(MPI_COMM_WORLD);

  errpetsc = MPI_Bcast(&node_proc_id[0], numOfNodeGlobal, MPI_INT, 0, MPI_COMM_WORLD); CHKERRQ(errpetsc);
  MPI_Barrier(MPI_COMM_WORLD);
  
  for(int ee=0; ee<numOfElmGlobal; ee++)
  {
    elm[ee]->setSubdomainId(elm_proc_id[ee]);
  }

  if(this_mpi_proc == 0){
    ofstream getsubdomain("getSubDomainId.dat");
    for(int ee=0; ee<numOfElmGlobal; ee++){
      getsubdomain << elm[ee]->getSubdomainId() << endl;
    }
    getsubdomain.close();
  }
  //cout << "wwwww" << endl;
  //exit(1);

  MPI_Barrier(MPI_COMM_WORLD);
  numOfElmLocal = count(elm_proc_id.begin(), elm_proc_id.end(), this_mpi_proc);
  cout << " nunOfElmLocal =  " << numOfElmLocal << '\t' << this_mpi_proc << '\t' << n_mpi_procs << endl;
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  return 0;
}


int FEM::prepareForParallel()
{

  int numOfSubDomain = n_mpi_procs;
  int subdomain=0;
  int ee, ii, jj, kk, n1, n2, ind;

  nNode_owned = count(node_proc_id.begin(), node_proc_id.end(), this_mpi_proc);
  cout << " nNode_owned =  " << nNode_owned << '\t' << this_mpi_proc << '\t' << n_mpi_procs << endl;

  MPI_Barrier(MPI_COMM_WORLD);
  PetscPrintf(MPI_COMM_WORLD, "\n"); 
  MPI_Barrier(MPI_COMM_WORLD);

  vector<int>  nodelist_owned(nNode_owned);

  kk=0;
  for(ii=0; ii<numOfNodeGlobal; ii++)
    {
      if( node_proc_id[ii] == this_mpi_proc )
      {
        nodelist_owned[kk++] = ii;
      }
    }
  cout << " Locally owned nodes " << '\t' << this_mpi_proc << endl;
  //printVector(nodelist_owned);
  MPI_Barrier(MPI_COMM_WORLD);

  vector<int>  nNode_owned_vector(n_mpi_procs), nNode_owned_sum(n_mpi_procs);

  MPI_Allgather(&nNode_owned, 1, MPI_INT, &nNode_owned_vector[0], 1, MPI_INT, MPI_COMM_WORLD);

  nNode_owned_sum = nNode_owned_vector;
  for(ii=1; ii<n_mpi_procs; ii++)
  {
    nNode_owned_sum[ii] += nNode_owned_sum[ii-1];
  }

  node_start = 0;
  if(this_mpi_proc > 0) node_start = nNode_owned_sum[this_mpi_proc-1];
  node_end   = nNode_owned_sum[this_mpi_proc]-1;

  cout << " node_start =  " << node_start << '\t' << node_end << '\t' << this_mpi_proc << endl;
  
  MPI_Barrier(MPI_COMM_WORLD);

  vector<int>  displs(n_mpi_procs);

  displs[0] = 0;
  for(ii=0; ii<n_mpi_procs-1; ii++) displs[ii+1] = displs[ii] + nNode_owned_vector[ii];

  errpetsc = MPI_Allgatherv(&nodelist_owned[0], nNode_owned, MPI_INT, &node_map_get_old[0], &nNode_owned_vector[0], &displs[0], MPI_INT, MPI_COMM_WORLD);

  for(ii=0; ii<numOfNodeGlobal; ii++)
  {
    //n1 = node_map_get_old[ii];
    //node_map_get_new[n1] = ii;
    //for(jj=0; jj<ndof; jj++)
    //{
    //NodeTypeNew[ii][jj] = NodeTypeOld[n1][jj];
    //}
  }

  if(this_mpi_proc == 0){
    ofstream node_map_old("node_map_get_old.dat");
    for(ii=0; ii<numOfNodeGlobal; ii++)
    { 
      node_map_old << node_map_get_old[ii] << endl;
    }
    node_map_old.close();
  }

  for(ee=0; ee<numOfElmGlobal; ee++)
  {
    ////erase////
    //for(ii=0; ii<numOfNodeInElm; ii++) element[ee][ii] = node_map_get_new[element[ee][ii]];

      elm[ee]->nodeNums = element[ee];

      // element array is no longer needed.
      element[ee].clear();
  }

   ///////add////////
  for(ii=0; ii<numOfNodeGlobal; ii++)
  {
    node_map_get_new[ii] = node_map_get_old[ii];
    for(jj=0; jj<ndof; jj++)
    {
    NodeTypeNew[ii][jj] = NodeTypeOld[ii][jj];
    }
  }
  ////////add////////

  if(this_mpi_proc == 1){
    ofstream out2("this_mpi_proc=1.dat");
    for(ee=0; ee<numOfElmGlobal; ee++){
      for(ii=0; ii<numOfNodeInElm; ii++) out2 << elm[ee]->nodeNums[ii] << " " ;
      out2 << endl;
    }
    out2.close();
  }
  
  element.clear();

  // update Dirichlet BC information with new node numbers
  for(ii=0; ii<numOfBdNode; ii++)
  { 
    ////erase///
    //n1 = node_map_get_new[DirichletBCs[ii][0]];
    ////erase///
    //DirichletBCs[ii][0] = n1;
    ///erase////
    //NodeTypeNew[n1][DirichletBCs[ii][1]] = true;
  }
  MPI_Barrier(MPI_COMM_WORLD);

   //compute NodeDofArrayNew

  ind = 0;
  for(ii=0; ii<numOfNodeGlobal; ii++)
  {
    for(jj=0; jj<ndof; jj++)
    {
      ///erase///
      //if(NodeTypeNew[ii][jj] == false)
      //{
        ///NodeDofArrayNew[ii][jj] = ind++;
        NodeDofArrayNew[ii][jj] = NodeDofArrayOld[ii][jj];
        ind++;
      //}
    }
  }


  if(this_mpi_proc == 0){  
    ofstream NodeNew("NodeDofArrayNew.dat");
    for(ii=0; ii<numOfNodeGlobal; ii++)
    {
      for(jj=0; jj<ndof; jj++)
      {
        NodeNew << NodeDofArrayNew[ii][jj] << " ";
      }
      NodeNew << endl;
    }
    NodeNew.close();
  }
  
  if(ind != ntotdofs_global)
  {
    cerr << "Something wrong with NodeDofArrayNew " << endl;
    cout << "ind = " << ind << '\t' << ntotdofs_global << '\t' << this_mpi_proc << endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  PetscPrintf(MPI_COMM_WORLD, "\n"); 
  MPI_Barrier(MPI_COMM_WORLD);

  // compute first and last row indices of the rows owned by the local processor
  //row_start  =  1e9;
  //row_end    = -1e9;
  ntotdofs_local = 0;
  
  row_start = node_start*ndof;
  row_end = node_end*ndof+ndof-1;

  for(ii=node_start; ii<=node_end; ii++)
  {
    for(jj=0; jj<ndof; jj++)
    {
      ////erase////
      //if(NodeTypeNew[ii][jj] == false)
      //{
        //ind = NodeDofArrayNew[ii][jj];
        //row_start  = min(row_start, ind);
        //row_end    = max(row_end,   ind);
        ntotdofs_local++;
      //}
    }
  }

  cout << " ntotdofs_local = " << ntotdofs_local << '\t' << ntotdofs_global << '\t' << this_mpi_proc << endl;

  cout << " row_start  = " << row_start  << '\t' << row_end  << '\t' << this_mpi_proc << endl;
  MPI_Barrier(MPI_COMM_WORLD);
  PetscPrintf(MPI_COMM_WORLD, "\n"); 
  MPI_Barrier(MPI_COMM_WORLD);
  //check if the sum of local problem sizes is equal to that of global problem size
  ind=0;
  errpetsc = MPI_Allreduce(&ntotdofs_local, &ind, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ind != ntotdofs_global)
  {
    //cerr << " Sum of local problem sizes is not equal to global size" << endl;
    cout << " ind = " << ind << '\t' << ntotdofs_global << '\t' << this_mpi_proc << endl;
  }
  return 0;
  
}




FEM::FEM()
{
  MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &this_mpi_proc);

  ntotdofs_local = ntotdofs_global = 0;
  elm = nullptr;
  solverPetsc = nullptr;
}


FEM::~FEM()
{
  if(elm != nullptr)
  {
    for(unsigned int ii=0;ii<numOfElmGlobal;++ii) delete elm[ii];

    delete [] elm;
    elm = nullptr;
  }
}

