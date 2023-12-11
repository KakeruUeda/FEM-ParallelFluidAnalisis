#include "FEM.h"

void FEM::prepareMatrix(){

  int  size1;
  int  jpn;
  int  r, ii, jj, ee, kk, nsize;
  int  *tt, nnz_max_row;
  int  count_diag, count_offdiag, tempInt;

  nodeMapPrevFluid.resize(numOfNodeGlobalFluid, 0);
  nodeMapFluid.resize(numOfNodeGlobalFluid, 0);
  nodeIdFluid.resize(numOfNodeGlobalFluid, 0);

  vector<bool>  vecBoolTempFalse(numOfDofsNode, false);
  nodeTypePrevFluid.resize(numOfNodeGlobalFluid, vecBoolTempFalse);
  nodeTypeFluid.resize(numOfNodeGlobalFluid, vecBoolTempFalse);

  vector<int>  vecIntTempM1(numOfDofsNode, -1);
  nodeDofArrayBCsPrevFluid.resize(numOfNodeGlobalFluid, vecIntTempM1);
  nodeDofArrayPrevFluid.resize(numOfNodeGlobalFluid, vecIntTempM1);
  nodeDofArrayBCsFluid.resize(numOfNodeGlobalFluid, vecIntTempM1);
  nodeDofArrayFluid.resize(numOfNodeGlobalFluid, vecIntTempM1);
  
  for(ii=0; ii<numOfBdNodeFluid; ii++){
    nodeTypePrevFluid[DirichletBCsFluid_tmp[ii][0]][DirichletBCsFluid_tmp[ii][1]] = true;
  }

  numOfDofsGlobalFluid = 0;
  for(int ii=0;ii<numOfNodeGlobalFluid;++ii){
    for(int jj=0;jj<numOfDofsNode;++jj){
        nodeDofArrayPrevFluid[ii][jj] = numOfDofsGlobalFluid;
        nodeDofArrayBCsPrevFluid[ii][jj] = numOfDofsGlobalFluid++;
      if(nodeTypePrevFluid[ii][jj]){
        nodeDofArrayBCsPrevFluid[ii][jj] = -1;
      }
    }
  }


  if(myId == 0)
  {
    cout << "\n Mesh statistics .....\n" << endl;
    cout << " numOfElmGlobalFluid   = " << '\t' << numOfElmGlobalFluid << endl;
    cout << " numOfNodeGlobalFluid  = " << '\t' << numOfNodeGlobalFluid << endl;
    cout << " numOfNodeInElm   = " << '\t' << numOfNodeInElm << endl;
    cout << " numOfDofsNode             = " << '\t' << numOfDofsNode << endl;
    cout << " numOfDofsLocalFluid   = " << '\t' << numOfDofsLocalFluid  << endl;
    cout << " numOfDofsGlobalFluid  = " << '\t' << numOfDofsGlobalFluid << endl;
    cout << " numOfId      = " << '\t' << numOfId << endl;
  }

  if(numOfId == 1)
  {
    elem_start = 0;
    elem_end   = numOfElmGlobal-1;
   
    numOfElmLocalFluid = numOfElmGlobalFluid;
    numOfNodeLocalFluid = numOfNodeGlobalFluid;
    numOfDofsLocalFluid  = numOfDofsGlobalFluid;
   
    row_start = 0;
    row_end   = numOfDofsGlobalFluid-1;
   
    for(ii=0; ii<numOfNodeGlobalFluid; ii++)
    {
      nodeMapPrevFluid[ii] = ii;
      nodeMapFluid[ii] = ii;
    }

    for(ee=0; ee<numOfElmGlobalFluid; ee++){
      elmFluid[ee]->nodeNumsPrevFluid = elementFluid[ee];
      elmFluid[ee]->nodeNumsFluid = elementFluid[ee];
      elementFluid[ee].clear();
    }
    elementFluid.clear();
   
    for(ii=0;ii<numOfNodeGlobalFluid;++ii)
    {
      nodeTypeFluid[ii] = nodeTypePrevFluid[ii];
      nodeDofArrayFluid[ii] = nodeDofArrayPrevFluid[ii];
      nodeDofArrayBCsFluid[ii] = nodeDofArrayBCsPrevFluid[ii];
    }
  }else{
    PetscPrintf(MPI_COMM_WORLD, "\n Mesh division start \n");
    divideMesh();
    PetscPrintf(MPI_COMM_WORLD, "\n Mesh division end \n");
    MPI_Barrier(MPI_COMM_WORLD);

    PetscPrintf(MPI_COMM_WORLD, "\n Parallel preparation start \n\n");
    prepareForParallel();
    PetscPrintf(MPI_COMM_WORLD, "\n Parallel preparation end \n");
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
 /// FLUID ONLY ///
  SolnDataFluid.nodeMapPrevFluid = nodeMapPrevFluid;
  SolnDataFluid.nodeMapFluid = nodeMapFluid;

  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  for(ee=0; ee<numOfElmGlobalFluid; ee++){
    if(elmFluid[ee]->getSubdomainId() == myId){

      nsize = numOfDofsNode*numOfNodeInElm;
      elmFluid[ee]->prepareElemData(nsize);
      numOfNodeInElm = elmFluid[ee]->nodeNumsFluid.size();

      elmFluid[ee]->nodeForAssyBCsFluid.resize(nsize);
      elmFluid[ee]->nodeForAssyFluid.resize(nsize);
      for(ii=0; ii<numOfNodeInElm; ++ii){
        jpn = numOfDofsNode *ii;
        kk = elmFluid[ee]->nodeNumsFluid[ii];
        for(jj=0;jj<numOfDofsNode ;++jj){
          elmFluid[ee]->nodeForAssyBCsFluid[jpn+jj] = nodeDofArrayBCsFluid[kk][jj];
          elmFluid[ee]->nodeForAssyFluid[jpn+jj] =  nodeDofArrayFluid[kk][jj];
        }
      }
    }
  }   
  /*
  if(myId == 0){
    ofstream Assy("nodeForAssyFluid0.dat"); 
    for(int ic=0; ic<numOfElmGlobalFluid; ic++){
      int size = elmFluid[ic]->nodeForAssyBCsFluid.size();
      for(int ii=0; ii<size; ii++){
        if(ii != 0 && ii%4 == 0) Assy << " ";
        Assy << elmFluid[ic]->nodeForAssyFluid[ii] << " ";
      } 
      Assy << endl;
    }
    Assy.close();
  }
  if(myId == 1){
    ofstream Assy1("nodeForAssyFluid1.dat"); 
    for(int ic=0; ic<numOfElmGlobalFluid; ic++){
      int size = elmFluid[ic]->nodeForAssyBCsFluid.size();
      for(int ii=0; ii<size; ii++){
        if(ii != 0 && ii%4 == 0) Assy1 << " ";
        Assy1 << elmFluid[ic]->nodeForAssyFluid[ii] << " ";
      } 
      Assy1 << endl;
    }
    Assy1.close();
  }
  */

  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  assyForSolnFluid.resize(numOfDofsGlobalFluid);
  jpn = 0;
  for(ii=0; ii<numOfNodeGlobalFluid; ii++)
  {
    for(jj=0; jj<numOfDofsNode; jj++)
    {
      assyForSolnFluid[jpn++] = ii*numOfDofsNode+jj;
    }
  }
  PetscPrintf(MPI_COMM_WORLD, "\n");
  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  printf(" numOfDofsLocalFluid = %5d \t numOfDofsGlobalFluid = %5d \t myId = %5d \n",numOfDofsLocalFluid, numOfDofsGlobalFluid, myId);
  errpetsc = MPI_Barrier(MPI_COMM_WORLD);
  
  vector<set<int>> forAssyMatFluid;
  set<int>::iterator it;

  forAssyMatFluid.resize(numOfDofsGlobalFluid);

  for(ee=0; ee<numOfElmGlobalFluid; ee++)
  {
    if(elmFluid[ee]->getSubdomainId() == myId)
    {
      tt = &(elmFluid[ee]->nodeForAssyFluid[0]);
      nsize = elmFluid[ee]->nodeForAssyFluid.size();

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
                forAssyMatFluid[r].insert(tt[jj]);
              }
            }
          }
        }
      }
    }
  }
  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  PetscInt  *diag_nnz, *offdiag_nnz;

  errpetsc  = PetscMalloc1(numOfDofsLocalFluid,  &diag_nnz);
  errpetsc  = PetscMalloc1(numOfDofsLocalFluid,  &offdiag_nnz);


  kk = 0;
  nnz_max_row = 0;
  for(ii=row_start; ii<=row_end; ii++)
  {
    size1 = forAssyMatFluid[ii].size();

    nnz_max_row = max(nnz_max_row, size1);

    count_diag=0, count_offdiag=0;
    for(it=forAssyMatFluid[ii].begin(); it!=forAssyMatFluid[ii].end(); ++it)
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
  solverPetscFluid->initialise(numOfDofsLocalFluid, numOfDofsGlobalFluid, diag_nnz, offdiag_nnz, nnz_max_row);
  PetscPrintf(MPI_COMM_WORLD, " Petsc solver initialize done\n\n");

  errpetsc = MPI_Barrier(MPI_COMM_WORLD);

  jpn = numOfNodeInElm*numOfDofsNode;
  jpn = jpn*jpn;
  PetscScalar  Klocal[jpn];
  
  MatrixXd Klocal_tmp(numOfNodeInElm*numOfDofsNode,numOfNodeInElm*numOfDofsNode);

  vector<int>  vecIntTemp;
  for(ee=0; ee<numOfElmGlobalFluid; ee++)
  {  
    for(ii=0; ii<jpn; ii++)  Klocal[ii] = 0.0;
    if(elmFluid[ee]->getSubdomainId() == myId)
    {
      size1 = elmFluid[ee]->nodeForAssyBCsFluid.size();
      vecIntTemp = elmFluid[ee]->nodeForAssyFluid;
      for(ii=0; ii<size1; ii++){
        if(elmFluid[ee]->nodeForAssyBCsFluid[ii] == -1){
          Klocal[ii*numOfNodeInElm*numOfDofsNode+ii] = 1;
        }
      }
      MatrixXdRM Klocal2_tmp = Klocal_tmp;
      errpetsc = MatSetValues(solverPetscFluid->mtx, size1, &vecIntTemp[0], size1, &vecIntTemp[0], Klocal, INSERT_VALUES);
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

  solverPetscFluid->currentStatus = PATTERN_OK;

  PetscFree(diag_nnz); 
  PetscFree(offdiag_nnz);  

  //TODO: use vector<T>().swap[x) to empty the vector and also deallocate the memory
  for(ii=0; ii<nodeDofArrayBCsFluid.size(); ii++)
  {
    nodeDofArrayBCsPrevFluid[ii].clear();
    nodeDofArrayBCsFluid[ii].clear();
    nodeDofArrayPrevFluid[ii].clear();
    nodeDofArrayFluid[ii].clear();
  }
  nodeDofArrayBCsPrevFluid.clear();
  nodeDofArrayBCsFluid.clear();
  nodeDofArrayPrevFluid.clear();
  nodeDofArrayFluid.clear();

}

int FEM::divideMesh()
{

  vector<int>  elmId(numOfElmGlobal, 0);
  vector<int>  elmIdFluid(numOfElmGlobalFluid, 0);

  if(myId == 0)
  {
    int  ee, ii, jj, kk, n2;
    int  nparts = numOfId, subdomain=0;

    PetscInt  *eptr, *eind;

    errpetsc  = PetscMalloc1(numOfElmGlobalFluid+1,  &eptr);CHKERRQ(errpetsc);
    errpetsc  = PetscMalloc1(numOfElmGlobalFluid*numOfNodeInElm,  &eind);CHKERRQ(errpetsc);

    vector<int>  vecTemp2;

    eptr[0] = 0;
    kk = 0;
    for(ee=0; ee<numOfElmGlobalFluid; ee++)
    {
      eptr[ee+1] = (ee+1)*numOfNodeInElm;
      //vecTemp2 = elems[ee]->nodeNums ;
      vecTemp2 = elementFluid[ee];
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


    // METIS partition routine
    int ret = METIS_PartMeshDual(&numOfElmGlobalFluid, &numOfNodeGlobalFluid, eptr, eind, NULL, NULL, &ncommon_nodes, &nparts, NULL, options, &objval, &elmIdFluid[0], &nodeIdFluid[0]);

    if(ret == METIS_OK)
      cout << "\n\n METIS partition routine success "  << endl << endl;
    else
      cout << " METIS partition routine FAILED "  << endl << endl;

    errpetsc = PetscFree(eptr); CHKERRQ(errpetsc);
    errpetsc = PetscFree(eind); CHKERRQ(errpetsc);
    
    /// check ///
    int tmp=0;
    for(int ii=0; ii<numOfElmGlobal; ii++){
      if(phi[ii]<1e-12){
        elmId[ii] = -1;
      }else{
        elmId[ii] = elmIdFluid[tmp++];
      }
    }
    nodeId.resize(numOfNodeGlobal,-1);
    
    int n1; 
    for(int ii=0; ii<numOfNodeGlobalFluid; ii++){
      n1 = nodeIdFluid[ii];
      nodeId[sortNode[ii]] = n1;
    }

    string vtiFile;
    vtiFile = outputDir + "/meshPartition.vti";
    export_vti(vtiFile, nodeId, elmId);
  }
  /**** prev
  MPI_Barrier(MPI_COMM_WORLD);
  exit(1);
  errpetsc = MPI_Bcast(&elmId[0], numOfElmGlobal, MPI_INT, 0, MPI_COMM_WORLD); CHKERRQ(errpetsc);
  MPI_Barrier(MPI_COMM_WORLD);

  errpetsc = MPI_Bcast(&nodeId[0], numOfNodeGlobal, MPI_INT, 0, MPI_COMM_WORLD); CHKERRQ(errpetsc);
  MPI_Barrier(MPI_COMM_WORLD);
  
  for(int ee=0; ee<numOfElmGlobal; ee++)
  {
    elm[ee]->setSubdomainId(elmId[ee]);
  }

  
  if(myId == 0){
    ofstream getsubdomain("getSubDomainId.dat");
    for(int ee=0; ee<numOfElmGlobal; ee++){
      getsubdomain << elm[ee]->getSubdomainId() << endl;
    }
    getsubdomain.close();
  }

  MPI_Barrier(MPI_COMM_WORLD);
  numOfElmLocal = count(elmId.begin(), elmId.end(), myId);
  //cout << " nunOfElmLocal =  " << numOfElmLocal << '\t' << myId << '\t' << numOfId << endl;
  
  MPI_Barrier(MPI_COMM_WORLD);
  */
  MPI_Barrier(MPI_COMM_WORLD);
  errpetsc = MPI_Bcast(&elmIdFluid[0], numOfElmGlobalFluid, MPI_INT, 0, MPI_COMM_WORLD); CHKERRQ(errpetsc);
  MPI_Barrier(MPI_COMM_WORLD);

  errpetsc = MPI_Bcast(&nodeIdFluid[0], numOfNodeGlobalFluid, MPI_INT, 0, MPI_COMM_WORLD); CHKERRQ(errpetsc);
  MPI_Barrier(MPI_COMM_WORLD);
  
  for(int ee=0; ee<numOfElmGlobalFluid; ee++)
  {
    elmFluid[ee]->setSubdomainId(elmIdFluid[ee]);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  numOfElmLocal = count(elmIdFluid.begin(), elmIdFluid.end(), myId);
  cout << " nunOfElmLocalFluid =  " << numOfElmLocalFluid << '\t' << myId << '\t' << numOfId << endl;
  
  MPI_Barrier(MPI_COMM_WORLD);


  
  return 0;
}


int FEM::prepareForParallel()
{

 int numOfSubDomain = numOfId;
 int subdomain=0;
 int ee, ii, jj, kk, n1, n2, jpn;

 nNode_owned = count(nodeIdFluid.begin(), nodeIdFluid.end(), myId);
 printf(" nNode_owned =  %5d \t numOfId = %5d \t myId = %5d \n",node_start, numOfId, myId);

 MPI_Barrier(MPI_COMM_WORLD);
 PetscPrintf(MPI_COMM_WORLD, "\n"); 
 MPI_Barrier(MPI_COMM_WORLD);

 vector<int>  nodelist_owned(nNode_owned);

 kk=0;
 for(ii=0; ii<numOfNodeGlobalFluid; ii++)
   {
     if( nodeIdFluid[ii] == myId )
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

 printf(" node_start = %5d \t node_end = %5d \t myId = %5d \n",node_start, node_end, myId);
 
 MPI_Barrier(MPI_COMM_WORLD);

 vector<int>  displs(numOfId);

 displs[0] = 0;
 for(ii=0; ii<numOfId-1; ii++) displs[ii+1] = displs[ii] + nNode_owned_vector[ii];

 errpetsc = MPI_Allgatherv(&nodelist_owned[0], nNode_owned, MPI_INT, &nodeMapPrevFluid[0], &nNode_owned_vector[0], &displs[0], MPI_INT, MPI_COMM_WORLD);

 for(ii=0; ii<numOfNodeGlobalFluid; ii++){
   n1 = nodeMapPrevFluid[ii];
   nodeMapFluid[n1] = ii;
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

 for(ee=0; ee<numOfElmGlobalFluid; ee++){

   elmFluid[ee]->nodeNumsPrevFluid = elementFluid[ee];
   
   for(ii=0; ii<numOfNodeInElm; ii++) elementFluid[ee][ii] = nodeMapFluid[elementFluid[ee][ii]];
   
   elmFluid[ee]->nodeNumsFluid = elementFluid[ee];
   elementFluid[ee].clear();
 
 }
 elementFluid.clear();

 // update Dirichlet BC information with new node numbers
 for(ii=0; ii<numOfBdNodeFluid; ii++)
 { 
   n1 = nodeMapFluid[DirichletBCsFluid_tmp[ii][0]];
   DirichletBCsFluid_tmp[ii][0] = n1;
   nodeTypeFluid[n1][DirichletBCsFluid_tmp[ii][1]] = true;
 }
 MPI_Barrier(MPI_COMM_WORLD);

 jpn = 0;
 for(ii=0; ii<numOfNodeGlobalFluid; ii++)
 {
   for(jj=0; jj<numOfDofsNode; jj++)
   {
     nodeDofArrayFluid[ii][jj] = jpn;
     nodeDofArrayBCsFluid[ii][jj] = jpn++;
     if(nodeTypeFluid[ii][jj]){
       nodeDofArrayBCsFluid[ii][jj] = -1;
     }
   }
 }

 /*
 if(myId == 0){  
   ofstream NodeNew("nodeDofArrayBCsFluid.dat");
   for(ii=0; ii<numOfNodeGlobalFluid; ii++)
   {
     for(jj=0; jj<numOfDofsNode; jj++)
     {
       NodeNew << nodeDofArrayBCsFluid[ii][jj] << " ";
     }
     NodeNew << endl;
   }
   NodeNew.close();
 }
 */
 
 if(jpn != numOfDofsGlobalFluid)
 {
   cerr << "Something wrong with NodeDofArrayNew " << endl;
   cout << "jpn = " << jpn << '\t' << numOfDofsGlobalFluid << '\t' << myId << endl;
 }

 MPI_Barrier(MPI_COMM_WORLD);
 PetscPrintf(MPI_COMM_WORLD, "\n");
 

 // compute first and last row indices of the rows owned by the local processor
 //row_start  =  1e9;
 //row_end    = -1e9;
 numOfDofsLocalFluid = 0;
 
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
       numOfDofsLocalFluid++;
     //}
   }
 }

 printf(" numOfDofsLocalFluid = %5d/%5d \t row_start  = %5d \t row_end  = %5d \t myId  = %5d \n", numOfDofsLocalFluid, numOfDofsGlobalFluid, row_start, row_end, myId);

 MPI_Barrier(MPI_COMM_WORLD);

 //check if the sum of local problem sizes is equal to that of global problem size
 jpn=0;
 errpetsc = MPI_Allreduce(&numOfDofsLocalFluid, &jpn, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

 if(jpn != numOfDofsGlobalFluid)
 {
   cerr << " Sum of local problem sizes is not equal to global size" << endl;
   cout << " jpn = " << jpn << '\t' << numOfDofsGlobalFluid << '\t' << myId << endl;
 }
 return 0;
  
}