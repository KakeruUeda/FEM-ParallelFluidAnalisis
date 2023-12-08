#include "FEM.h"

void FEM::prepareMatrix(){

  int  size1;
  int  jpn;
  int  r, ii, jj, ee, kk, nsize;
  int  *tt, nnz_max_row;
  int  count_diag, count_offdiag, tempInt;
  
  nodeMapPrev.resize(numOfNodeGlobal, 0);
  nodeMap.resize(numOfNodeGlobal, 0);
  nodeId.resize(numOfNodeGlobal, 0);

  vector<bool>  vecBoolTempFalse(numOfDofsNode, false);
  nodeTypePrev.resize(numOfNodeGlobal, vecBoolTempFalse);
  nodeType.resize(numOfNodeGlobal, vecBoolTempFalse);

  vector<int>  vecIntTempM1(numOfDofsNode, -1);
  nodeDofArrayBCsPrev.resize(numOfNodeGlobal, vecIntTempM1);
  nodeDofArrayPrev.resize(numOfNodeGlobal, vecIntTempM1);
  nodeDofArrayBCs.resize(numOfNodeGlobal, vecIntTempM1);
  nodeDofArray.resize(numOfNodeGlobal, vecIntTempM1);
  
  for(ii=0; ii<numOfBdNode; ii++){
    nodeTypePrev[DirichletBCs_tmp[ii][0]][DirichletBCs_tmp[ii][1]] = true;
  }

  if(myId == 0){
    ofstream nodeTypeP("nodeTypePrev0.dat"); 
    for(int ic=0; ic<numOfNodeGlobal; ic++){
      for(int j=0; j<4; j++){
        nodeTypeP << "ic = " << ic << " nodeTypePrev[ic][j] = " << nodeTypePrev[ic][j] << " ";
      } 
      nodeTypeP << endl;
    }
    nodeTypeP.close();
  }

  
  numOfDofsGlobal = 0;
  for(int ii=0;ii<numOfNodeGlobal;ii++){
    for(int jj=0;jj<numOfDofsNode;jj++){
        nodeDofArrayPrev[ii][jj] = numOfDofsGlobal;
        nodeDofArrayBCsPrev[ii][jj] = numOfDofsGlobal++;
      if(nodeTypePrev[ii][jj]){
        nodeDofArrayBCsPrev[ii][jj] = -1;
      }
    }
  }

  if(myId == 0)
  {
    cout << "\n Mesh statistics\n" << endl;
    cout << " numOfElmGlobal   = " << '\t' << numOfElmGlobal << endl;
    cout << " numOfNodeGlobal  = " << '\t' << numOfNodeGlobal  << endl;
    cout << " numOfNodeInElm   = " << '\t' << numOfNodeInElm << endl;
    cout << " numOfDofsNode    = " << '\t' << numOfDofsNode << endl;
    cout << " numOfDofsLocal   = " << '\t' << numOfDofsLocal  << endl;
    cout << " numOfDofsGlobal  = " << '\t' << numOfDofsGlobal << endl;
    cout << " numOfId          = " << '\t' << numOfId << endl;
  }

  if(numOfId == 1)
  {
    elem_start = 0;
    elem_end   = numOfElmGlobal-1;
   
    numOfElmLocal = numOfElmGlobal;
    numOfNodeLocal = numOfNodeGlobal;
    numOfDofsLocal  = numOfDofsGlobal;
   
    row_start = 0;
    row_end   = numOfDofsGlobal-1;
   
    for(ii=0; ii<numOfNodeGlobal; ii++)
    {
      nodeMapPrev[ii] = ii;
      nodeMap[ii] = ii;
    }

    for(ee=0; ee<numOfElmGlobal; ee++){
      elm[ee]->nodeNumsPrev = element[ee];
      elm[ee]->nodeNums = element[ee];
      element[ee].clear();
    }
    element.clear();
   
    for(ii=0;ii<numOfNodeGlobal;++ii)
    {
      nodeType[ii] = nodeTypePrev[ii];
      nodeDofArray[ii] = nodeDofArrayPrev[ii];
      nodeDofArrayBCs[ii] = nodeDofArrayBCsPrev[ii];
    }
  }else{
    PetscPrintf(MPI_COMM_WORLD, "\n Mesh division start");
    divideMesh();
    PetscPrintf(MPI_COMM_WORLD, " Mesh division end \n\n");
    MPI_Barrier(MPI_COMM_WORLD);
  
    PetscPrintf(MPI_COMM_WORLD, " Parallel preparation start \n\n");
    prepareForParallel();
    PetscPrintf(MPI_COMM_WORLD, "\n Parallel preparation end \n\n");
    
    MPI_Barrier(MPI_COMM_WORLD);
  }

  SolnData.nodeMapPrev = nodeMapPrev;
  SolnData.nodeMap = nodeMap;

  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  for(ee=0; ee<numOfElmGlobal; ee++){
    if(elm[ee]->getSubdomainId() == myId){

      nsize = numOfDofsNode*numOfNodeInElm;
      elm[ee]->prepareElemData(nsize);
      numOfNodeInElm = elm[ee]->nodeNums.size();

      elm[ee]->nodeForAssyBCs.resize(nsize);
      elm[ee]->nodeForAssy.resize(nsize);
      for(ii=0; ii<numOfNodeInElm; ii++){
        jpn = numOfDofsNode*ii;
        kk = elm[ee]->nodeNums[ii];
        for(jj=0; jj<numOfDofsNode; jj++){
          elm[ee]->nodeForAssyBCs[jpn+jj] = nodeDofArrayBCs[kk][jj];
          elm[ee]->nodeForAssy[jpn+jj] =  nodeDofArray[kk][jj];
        }
      }
    }
  }   
  
  if(myId == 0){
    ofstream Assy("nodeForAssyBCs0.dat"); 
    for(int ic=0; ic<numOfElmGlobal; ic++){
      int size = elm[ic]->nodeForAssyBCs.size();
      for(int ii=0; ii<size; ii++){
        if(ii != 0 && ii%4 == 0) Assy << " ";
        Assy << elm[ic]->nodeForAssyBCs[ii] << " ";
      } 
      Assy << endl;
    }
    Assy.close();
  }
  if(myId == 1){
    ofstream Assy1("nodeForAssyBCs1.dat"); 
    for(int ic=0; ic<numOfElmGlobal; ic++){
      int size = elm[ic]->nodeForAssyBCs.size();
      for(int ii=0; ii<size; ii++){
        if(ii != 0 && ii%4 == 0) Assy1 << " ";
        Assy1 << elm[ic]->nodeForAssyBCs[ii] << " ";
      } 
      Assy1 << endl;
    }
    Assy1.close();
  }
  

  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  assyForSoln.resize(numOfDofsGlobal);
  jpn = 0;
  for(ii=0;ii<numOfNodeGlobal;++ii)
  {
    for(jj=0;jj<numOfDofsNode ;++jj)
    {
      assyForSoln[jpn++] = ii*numOfDofsNode+jj;
    }
  }

  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  printf(" numOfDofsLocal = %5d \t numOfDofsGlobal = %5d \t myId = %5d \n",numOfDofsLocal, numOfDofsGlobal, myId);
  
  errpetsc = MPI_Barrier(MPI_COMM_WORLD);
  vector<set<int> > forAssyMat;
  set<int>::iterator it;

  forAssyMat.resize(numOfDofsGlobal);

  for(ee=0; ee<numOfElmGlobal; ee++)
  {
    if(elm[ee]->getSubdomainId() == myId)
    {
      
      tt = &(elm[ee]->nodeForAssy[0]);
      nsize = elm[ee]->nodeForAssy.size();

      for(ii=0;ii<nsize;ii++)
      {
        r = tt[ii];
        //printf(" ee = %5d \t ii = %5d \t tt[ii] = %5d \t elm[ee]->getSubdomainId() %5d \t myId = %5d \n",ee, ii, r, elm[ee]->getSubdomainId(), myId);
        if(r >= row_start && r <= row_end)
        {
          //printf("IN... ee = %5d \t r = %5d \t elm[ee]->getSubdomainId() %5d \t myId = %5d \n",ee, r, elm[ee]->getSubdomainId(), myId);
          for(jj=0;jj<nsize;jj++)
          {
            forAssyMat[r].insert(tt[jj]);
          }
        }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  int c=0;

  for(ii=row_start; ii<=row_end; ii++)
  {
    nsize = elm[ii]->nodeForAssy.size();
    for(int r=0; r<numOfElmGlobal; r++)
    {
      for(it=forAssyMat[ii].begin(); it!=forAssyMat[ii].end(); ++it)
      {
        tempInt = *it;
        //printf(" r = %5d \t count = %5d \t tmpInt = %5d \t myId = %5d\n",r, c, tempInt, myId);
        c++;
      }
    }
  }

  PetscInt  *diag_nnz, *offdiag_nnz;

  errpetsc  = PetscMalloc1(numOfDofsLocal,  &diag_nnz);
  errpetsc  = PetscMalloc1(numOfDofsLocal,  &offdiag_nnz);


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
    diag_nnz[kk]    = count_diag;
    offdiag_nnz[kk] = count_offdiag;
    kk++;
  }

  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  PetscPrintf(MPI_COMM_WORLD, "\n Petsc solver initialize start\n");
  solverPetsc->initialise(numOfDofsLocal, numOfDofsGlobal, diag_nnz, offdiag_nnz, nnz_max_row);
  PetscPrintf(MPI_COMM_WORLD, " Petsc solver initialize done\n\n");
  
  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  jpn = numOfNodeInElm*numOfDofsNode;
  jpn = jpn*jpn;
  PetscScalar  Klocal[jpn];
  
  MatrixXd Klocal_tmp(numOfNodeInElm*numOfDofsNode,numOfNodeInElm*numOfDofsNode);

  vector<int>  vecIntTemp;
  for(ee=0; ee<numOfElmGlobal; ee++)
  {  
    for(ii=0; ii<jpn; ii++)  Klocal[ii] = 0.0;
    
    if(elm[ee]->getSubdomainId() == myId)
    {
      size1 = elm[ee]->nodeForAssyBCs.size();
      vecIntTemp = elm[ee]->nodeForAssy;
      for(ii=0; ii<size1; ii++){
        if(elm[ee]->nodeForAssyBCs[ii] == -1){
          Klocal[ii*numOfNodeInElm*numOfDofsNode+ii] = 1;
        }
      }
      //MatrixXdRM Klocal2_tmp = Klocal_tmp;
      errpetsc = MatSetValues(solverPetsc->mtx, size1, &vecIntTemp[0], size1, &vecIntTemp[0], Klocal, INSERT_VALUES);
    }
  }
  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  //errpetsc = MatAssemblyBegin(solverPetsc->mtx,MAT_FINAL_ASSEMBLY); 
  //errpetsc = MatAssemblyEnd(solverPetsc->mtx,MAT_FINAL_ASSEMBLY); 

  //errpetsc = MatView(solverPetsc->mtx,PETSC_VIEWER_STDOUT_WORLD); 
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
  for(ii=0; ii<nodeDofArrayBCs.size(); ii++)
  {
    nodeDofArrayBCsPrev[ii].clear();
    nodeDofArrayBCs[ii].clear();
    nodeDofArrayPrev[ii].clear();
    nodeDofArray[ii].clear();
  }
  nodeDofArrayBCsPrev.clear();
  nodeDofArrayBCs.clear();
  nodeDofArrayPrev.clear();
  nodeDofArray.clear();
}

int FEM::divideMesh()
{

  vector<int>  elmId(numOfElmGlobal, 0);

  if(myId == 0)
  {
    int  ee, ii, jj, kk, n2;
    int  nparts = numOfId, subdomain=0;

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

    ncommon_nodes = 4;  // 8-noded hexa element
  
    idx_t objval;
    idx_t options[METIS_NOPTIONS];

    METIS_SetDefaultOptions(options);

    // Specifies the partitioning method.
    //options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;          // Multilevel recursive bisectioning.
    options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;        

    //options[METIS_OPTION_NSEPS] = 10;

    //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;     // Edge-cut minimization
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;     // Total communication volume minimization

    options[METIS_OPTION_NUMBERING] = 0;  // C-style numbering is assumed that starts from 0.

    //cout << " Executing METIS subroutine " << endl;
    // METIS partition routine
    int ret = METIS_PartMeshDual(&numOfElmGlobal, &numOfNodeGlobal, eptr, eind, NULL, NULL, &ncommon_nodes, &nparts, NULL, options, &objval, &elmId[0], &nodeId[0]);

    if(ret == METIS_OK)
      cout << "\n\n METIS partition routine success "  << endl << endl;
    else
      cout << " METIS partition routine FAILED "  << endl << endl;

    errpetsc = PetscFree(eptr); CHKERRQ(errpetsc);
    errpetsc = PetscFree(eind); CHKERRQ(errpetsc);
    
    /// check ///
    string vtiFile;
    vtiFile = "check_MeshPartition.vti";
    export_vti(vtiFile, nodeId, elmId);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  errpetsc = MPI_Bcast(&elmId[0], numOfElmGlobal, MPI_INT, 0, MPI_COMM_WORLD); CHKERRQ(errpetsc);
  MPI_Barrier(MPI_COMM_WORLD);

  errpetsc = MPI_Bcast(&nodeId[0], numOfNodeGlobal, MPI_INT, 0, MPI_COMM_WORLD); CHKERRQ(errpetsc);
  MPI_Barrier(MPI_COMM_WORLD);
  
  for(int ee=0; ee<numOfElmGlobal; ee++)
  {
    elm[ee]->setSubdomainId(elmId[ee]);
  }

  /*
  if(myId == 0){
    ofstream getsubdomain("getSubDomainId.dat");
    for(int ee=0; ee<numOfElmGlobal; ee++){
      getsubdomain << elm[ee]->getSubdomainId() << endl;
    }
    getsubdomain.close();
  }
  */

  MPI_Barrier(MPI_COMM_WORLD);
  numOfElmLocal = count(elmId.begin(), elmId.end(), myId);
  //cout << " nunOfElmLocal =  " << numOfElmLocal << '\t' << myId << '\t' << numOfId << endl;
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  return 0;
}


int FEM::prepareForParallel()
{

  int numOfSubDomain = numOfId;
  int subdomain=0;
  int ee, ii, jj, kk, n1, n2, jpn;

  nNode_owned = count(nodeId.begin(), nodeId.end(), myId);
  printf(" nNode_owned =  %5d \t numOfId = %5d \t myId = %5d \n",node_start, numOfId, myId);
  
  MPI_Barrier(MPI_COMM_WORLD);
  PetscPrintf(MPI_COMM_WORLD, "\n"); 
  MPI_Barrier(MPI_COMM_WORLD);

  vector<int>  nodelist_owned(nNode_owned);

  kk=0;
  for(ii=0; ii<numOfNodeGlobal; ii++)
    {
      if( nodeId[ii] == myId )
      {
        nodelist_owned[kk++] = ii;
      }
    }

  MPI_Barrier(MPI_COMM_WORLD);

  vector<int>  nNode_owned_vector(numOfId), nNode_owned_sum(numOfId);

  MPI_Allgather(&nNode_owned, 1, MPI_INT, &nNode_owned_vector[0], 1, MPI_INT, MPI_COMM_WORLD);

  nNode_owned_sum = nNode_owned_vector;
  for(ii=1; ii<numOfId; ii++)
  {
    nNode_owned_sum[ii] += nNode_owned_sum[ii-1];
  }

  node_start = 0;
  if(myId > 0) node_start = nNode_owned_sum[myId-1];
  node_end   = nNode_owned_sum[myId]-1;

  //cout << " node_start =  " << node_start << '\t' << node_end << '\t' << myId << endl;
  printf(" node_start = %5d \t node_end = %5d \t myId = %5d \n",node_start, node_end, myId);
  
  MPI_Barrier(MPI_COMM_WORLD);

  vector<int>  displs(numOfId);

  displs[0] = 0;
  for(ii=0; ii<numOfId-1; ii++) displs[ii+1] = displs[ii] + nNode_owned_vector[ii];

  errpetsc = MPI_Allgatherv(&nodelist_owned[0], nNode_owned, MPI_INT, &nodeMapPrev[0], &nNode_owned_vector[0], &displs[0], MPI_INT, MPI_COMM_WORLD);

  for(ii=0; ii<numOfNodeGlobal; ii++){
    n1 = nodeMapPrev[ii];
    nodeMap[n1] = ii;
  }
  /*
  if(myId == 0){
    ofstream node_map_old("nodeMapPrev.dat");
    for(ii=0; ii<numOfNodeGlobal; ii++)
    { 
      node_map_old << nodeMapPrev[ii] << endl;
    }
    node_map_old.close();
  }
  */

  for(ee=0; ee<numOfElmGlobal; ee++){

    elm[ee]->nodeNumsPrev = element[ee];
    
    for(ii=0; ii<numOfNodeInElm; ii++) element[ee][ii] = nodeMap[element[ee][ii]];
    
    elm[ee]->nodeNums = element[ee];
    element[ee].clear();
  
  }
  element.clear();

  // update Dirichlet BC information with new node numbers
  for(ii=0; ii<numOfBdNode; ii++)
  { 
    n1 = nodeMap[DirichletBCs_tmp[ii][0]];
    DirichletBCs_tmp[ii][0] = n1;
    nodeType[n1][DirichletBCs_tmp[ii][1]] = true;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  jpn = 0;
  for(ii=0; ii<numOfNodeGlobal; ii++)
  {
    for(jj=0; jj<numOfDofsNode; jj++)
    {
      nodeDofArray[ii][jj] = jpn;
      nodeDofArrayBCs[ii][jj] = jpn++;
      if(nodeType[ii][jj]){
        nodeDofArrayBCs[ii][jj] = -1;
      }
    }
  }

 /*
  if(myId == 0){  
    ofstream NodeNew("nodeDofArrayBCs.dat");
    for(ii=0; ii<numOfNodeGlobal; ii++)
    {
      for(jj=0; jj<numOfDofsNode; jj++)
      {
        NodeNew << nodeDofArrayBCs[ii][jj] << " ";
      }
      NodeNew << endl;
    }
    NodeNew.close();
  }
  */
  
  if(jpn != numOfDofsGlobal)
  {
    cerr << "Something wrong with NodeDofArrayNew " << endl;
    cout << "jpn = " << jpn << '\t' << numOfDofsGlobal << '\t' << myId << endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  PetscPrintf(MPI_COMM_WORLD, "\n"); 
  MPI_Barrier(MPI_COMM_WORLD);

  // compute first and last row indices of the rows owned by the local processor
  //row_start  =  1e9;
  //row_end    = -1e9;
  numOfDofsLocal = 0;
  
  row_start = node_start*numOfDofsNode;
  row_end = node_end*numOfDofsNode+numOfDofsNode-1;

  for(ii=node_start; ii<=node_end; ii++)
  {
    for(jj=0; jj<numOfDofsNode; jj++)
    {
      ////erase////
      //if(NodeTypeNew[ii][jj] == false)
      //{
        //ind = NodeDofArrayNew[ii][jj];
        //row_start  = min(row_start, ind);
        //row_end    = max(row_end,   ind);
        numOfDofsLocal++;
      //}
    }
  }
  printf(" numOfDofsLocal = %5d/%5d \t row_start  = %5d \t row_end  = %5d \t myId  = %5d \n", numOfDofsLocal, numOfDofsGlobal, row_start, row_end, myId);
  
  MPI_Barrier(MPI_COMM_WORLD);
  //check if the sum of local problem sizes is equal to that of global problem size
  jpn=0;
  errpetsc = MPI_Allreduce(&numOfDofsLocal, &jpn, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(jpn != numOfDofsGlobal)
  {
    //cerr << " Sum of local problem sizes is not equal to global size" << endl;
    cout << " jpn = " << jpn << '\t' << numOfDofsGlobal << '\t' << myId << endl;
  }
  return 0;
  
}