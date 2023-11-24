#include "FEM.h"
using namespace std;

void FEM::initialize(string& domainFile)
{
  readInput(domainFile);
  setDomain();
  //setBoundary();
  //setImage();
  partitionMesh();
}

void FEM::readInput(string& domainFile)
{
  PetscPrintf(MPI_COMM_WORLD, "\nread MeshFile\n");
  
  ifstream  infile(domainFile);
  if(infile.fail())
  {
    cout << " Could not open 'domain.dat' file " << endl;
    exit(1);
  }

  string value;

  infile >> value >> nx;
  infile >> value >> ny;
  infile >> value >> nz;

  infile >> value >> Lx;
  infile >> value >> Ly;
  infile >> value >> Lz;
  
  infile.close();
  
  PetscPrintf(MPI_COMM_WORLD, "Control parameters are successfully read\n\n");
}

void FEM::setDomain()
{
  dx=Lx/nx;
  dy=Ly/ny;
  dz=Lz/nz; 
  numOfNode=(nx+1)*(ny+1)*(nz+1);
  numOfElm=nx*ny*nz;
  x.resize(numOfNode,vector<double>(2, 0));
  element.resize(numOfElm);
  
  if(this_mpi_proc == 0){
    cout << "nx = " << nx << endl;
    cout << "ny = " << ny << endl;
    cout << "nz = " << nz << endl;

    cout << "Lx = " << Lx << endl;
    cout << "Ly = " << Ly << endl;
    cout << "Lz = " << Lz << endl;
  }

}

void FEM::partitionMesh()
{
  /*
  PetscPrintf(MPI_COMM_WORLD, "\n     StabFEM::partitionMesh()  .... STARTED ...\n");

  vector<int>  elem_proc_id(nElem_global, 0);

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

    errpetsc  = PetscMalloc1(nElem_global+1,  &eptr);CHKERRQ(errpetsc);
    errpetsc  = PetscMalloc1(nElem_global*npElem,  &eind);CHKERRQ(errpetsc);

    vector<int>  vecTemp2;

    eptr[0] = 0;
    kk = 0;
    for(ee=0; ee<nElem_global; ee++)
    {
      eptr[ee+1] = (ee+1)*npElem;

      //vecTemp2 = elems[ee]->nodeNums ;
      vecTemp2 = elemConn[ee];
      //printVector(vecTemp2);

      for(ii=0; ii<npElem; ii++)
        eind[kk+ii] = vecTemp2[ii] ;

      kk += npElem;
    }

    int  ncommon_nodes;
    if(ndim == 2)
      ncommon_nodes = 2;    // 3-noded tria or 4-noded quad
    else
    {
      if(npElem == 4)       // 4-noded tetra element
        ncommon_nodes = 3;
      else
        ncommon_nodes = 4;  // 8-noded hexa element
    }

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
    int ret = METIS_PartMeshDual(&nElem_global, &nNode_global, eptr, eind, NULL, NULL, &ncommon_nodes, &nparts, NULL, options, &objval, &elem_proc_id[0], &node_proc_id[0]);

    if(ret == METIS_OK)
      cout << " METIS partition routine successful "  << endl;
    else
      cout << " METIS partition routine FAILED "  << endl;

    errpetsc = PetscFree(eptr); CHKERRQ(errpetsc);
    errpetsc = PetscFree(eind); CHKERRQ(errpetsc);

    if( 1 < 0)
    {
      for(ee=0; ee<nNode_global; ee++)
        cout << ee << '\t' << node_proc_id[ee] << endl;
      cout << endl;  cout << endl;  cout << endl;

      for(ee=0; ee<nElem_global; ee++)
        cout << ee << '\t' << elem_proc_id[ee] << endl;
      cout << endl;  cout << endl;  cout << endl;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  errpetsc = MPI_Bcast(&elem_proc_id[0], nElem_global, MPI_INT, 0, MPI_COMM_WORLD); CHKERRQ(errpetsc);
  MPI_Barrier(MPI_COMM_WORLD);

  errpetsc = MPI_Bcast(&node_proc_id[0], nNode_global, MPI_INT, 0, MPI_COMM_WORLD); CHKERRQ(errpetsc);
  MPI_Barrier(MPI_COMM_WORLD);

  for(int ee=0; ee<nElem_global; ee++)
  {
    elems[ee]->setSubdomainId(elem_proc_id[ee]);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  nElem_local = std::count(elem_proc_id.begin(), elem_proc_id.end(), this_mpi_proc);
  cout << " nElem_local =  " << nElem_local << '\t' << this_mpi_proc << '\t' << n_mpi_procs << endl;

  PetscPrintf(MPI_COMM_WORLD, "\n     StabFEM::partitionMesh()  .... FINISHED ...\n");
  */
  return;
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
}


FEM::~FEM()
{

}
