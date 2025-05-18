#include "FEM.h"

#include "PostProcess.h"

int main(int argc, char * argv[])
{
	MPI_Init(NULL, NULL);
	
	string petscfile = argv[2];
	PetscInitialize(NULL, NULL, petscfile.c_str(), NULL);

	FEM navier;
	PostProcess post;
	navier.solver = SOLVER::STEADY_NAVIERSTOKES;

	string input_file = argv[1];
	int ierror;
	if ((ierror = navier.tp.read(input_file)) != TP_NO_ERROR) {
		printf("\t Error at reading '%s' file\n", input_file.c_str());
		return 1;
	}
	string input_file_post = argv[1];
	int ierror_post;
	if ((ierror_post = post.tp_post.read(input_file_post)) != TP_NO_ERROR) {
		printf("\t Error at reading '%s' file\n", input_file_post.c_str());
		return 1;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	navier.initialize();

	MPI_Barrier(MPI_COMM_WORLD);
	post.readDAParam();

	int tmp = 0;

	/*
	MPI_Barrier(MPI_COMM_WORLD);
	for (auto it = post.numOfObsVoxels.begin(); it != post.numOfObsVoxels.end(); it++) {
		post.validateDADomain(navier, post.numOfObsVoxels[tmp][0],
		post.numOfObsVoxels[tmp][1], post.numOfObsVoxels[tmp][2]);
		tmp++;
	}
	*/

	MPI_Barrier(MPI_COMM_WORLD);
	navier.visualizeDomain();

	navier.nu = navier.mu / navier.rho;
	navier.Re = 1e0 / navier.nu;

	MPI_Barrier(MPI_COMM_WORLD);
	navier.SteadyNavierStokes();

	MPI_Barrier(MPI_COMM_WORLD);
	navier.visualizeResults();

	tmp = 0;

	MPI_Barrier(MPI_COMM_WORLD);
	for (auto it = post.numOfObsVoxels.begin(); it != post.numOfObsVoxels.end(); it++) {
		post.prepareForDataAssimilation(navier, post.numOfObsVoxels[tmp][0],
			post.numOfObsVoxels[tmp][1], post.numOfObsVoxels[tmp][2]);
		tmp++;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	post.extractDomain(navier);

	MPI_Barrier(MPI_COMM_WORLD);
	navier.solverDeallocate();

	MPI_Barrier(MPI_COMM_WORLD);
	PetscFinalize();
	MPI_Finalize();

	return 0;
}