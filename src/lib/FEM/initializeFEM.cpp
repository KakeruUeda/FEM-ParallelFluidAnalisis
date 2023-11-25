#include "FEM.h"
using namespace std;

void FEM::initialize()
{
  readInput();
  setDomain();
  //setBoundary();
  //setImage();

  partitionMesh();
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
          //elementfile << element[tmp2].node[0] << " " <<  element[tmp2].node[1] << " " <<  element[tmp2].node[2] << " " <<  element[tmp2].node[3]
          // <<" " <<   element[tmp2].node[4] << " " <<  element[tmp2].node[5] << " " <<  element[tmp2].node[6] << " " <<  element[tmp2].node[7] << endl;
          tmp2++;
      }
    }
  }
  prepare();
}

void FEM::prepare()
{
  elm = new ElementBaseFEM* [numOfElmGlobal];
}

int FEM::partitionMesh()
{
  
  PetscPrintf(MPI_COMM_WORLD, "\n     StabFEM::partitionMesh()  .... STARTED ...\n");

  node_proc_id.resize(numOfNodeGlobal, 0);

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

    ofstream out_nodeId("check_nodeProcId.dat");
    ofstream out_elmId("check_elmProcId.dat");
    for(ee=0; ee<numOfNodeGlobal; ee++)
    {
      out_nodeId << ee <<  " node_proc_id[ee] = " << node_proc_id[ee] << endl;
    }
    for(ee=0; ee<numOfElmGlobal; ee++)
    {
      out_elmId << ee << " elm_proc_id[ee] = " << elm_proc_id[ee] << endl;
    }
    out_nodeId.close();
    out_elmId.close();

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
  MPI_Barrier(MPI_COMM_WORLD);
  cout << "4" << endl;
  numOfElmLocal = count(elm_proc_id.begin(), elm_proc_id.end(), this_mpi_proc);
  cout << " nunOfElmLocal =  " << numOfElmLocal << '\t' << this_mpi_proc << '\t' << n_mpi_procs << endl;

  PetscPrintf(MPI_COMM_WORLD, "\n     StabFEM::partitionMesh()  .... FINISHED ...\n");
  
  return 0;
}


FEM::FEM()
{
  MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &this_mpi_proc);

  for (int i = 0; i < n_mpi_procs; i++) {
    if (i == this_mpi_proc) {
      printf("I'm rank %d\n", this_mpi_proc);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  elm = nullptr;
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
