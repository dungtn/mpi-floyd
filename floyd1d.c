/*
 ======================================================================
 Name        : floyd1d.c
 Author      : daisy
 Version     : 0.1
 Description : MPI Floyd-Warshall algorithm in C, Ansi-style
 Copyright   : Copyright 2012 Dung Thai
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
 ======================================================================
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define ROW 0
#define COL 1
#define INF 99

#define FILE_NOT_FOUND 404

inline int min(int a, int b) {
  return a < b ? a : b;
}

int floyd_all_pairs_sp_1d(int n, int nlocal, int *a) {
	int i, j, k, myrank;

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	int *krow = (int *)malloc(n * sizeof(int));
	int *kcol = (int *)malloc(nlocal * sizeof(int));

	for(k = 0; k < n; k++) {
		for(i = 0; i < n; i++) {
			krow[i] = a[(k % nlocal) * n + i];
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&krow[0], n, MPI_INT, k/nlocal, MPI_COMM_WORLD);

		for(i = 0; i < nlocal; i++) {
			kcol[i] = a[i * n + k];
		}
		for(i = 0; i < nlocal; i++) {
			for(j = 0; j < n; j++) {
				a[i * n + j] = min(a[i * n + j], kcol[i] + krow[j]);
//				if(myrank == 0) {
//					printf("k=%d kcol[%d]=%d krow[%d]=%d a[%d]=%d\n", k, i, kcol[i], j, krow[j], i * n + j, a[i * n + j]);
//				}
			}
		}
	}
	free(krow);
	free(kcol);
	return 0;
}

int main(int argc, char *argv[]) {

	int i, n, nlocal;
	int numprocs, myrank;
	MPI_File f; char* filename = "input/8";
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &f) != MPI_SUCCESS) {
		fprintf(stderr, "Cannot open file %s\n", filename);
		MPI_Abort(MPI_COMM_WORLD, FILE_NOT_FOUND);
		MPI_Finalize();
		return 1;
	}
	MPI_File_seek(f, 0, MPI_SEEK_SET);
	MPI_File_read(f, &n, 1, MPI_INT, &status);
	nlocal = n/numprocs; if(myrank == numprocs - 1) nlocal = nlocal + n % numprocs;

	int *a = (int *)malloc(nlocal * n * sizeof(int));
	MPI_File_seek(f, (myrank * nlocal * n + 1) * sizeof(int), MPI_SEEK_SET);
	MPI_File_read(f, &a[0], nlocal * n, MPI_INT, &status);
	MPI_File_close(&f);

// 	int j;
//	if(myrank == 3) {
//		for(i = 0; i < nlocal; i++) {
//			for(j = 0; j < n; j++) {
//				printf("%d ", a[i * n +j]);
//			}
//			printf("\n");
//		}
//	}

	double start = MPI_Wtime();
	floyd_all_pairs_sp_1d(n, nlocal, a);
	double stop = MPI_Wtime();
	printf("[%d] Completed in %1.3f seconds\n", myrank, stop-start);

//	if(myrank == 3) {
//		for(i = 0; i < nlocal; i++) {
//			for(j = 0; j < n; j++) {
//				printf("%d ", a[i * n +j]);
//			}
//			printf("\n");
//		}
//	}
	if(MPI_File_open(MPI_COMM_WORLD, "output/8", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &f) != MPI_SUCCESS) {
			printf("Cannot open file %s\n", "out");
			MPI_Abort(MPI_COMM_WORLD, FILE_NOT_FOUND);
			MPI_Finalize();
			return 1;
	}
	if(myrank == 0) {
		MPI_File_seek(f, 0, MPI_SEEK_SET);
		MPI_File_write(f, &n, 1, MPI_INT, &status);
	}
	for(i = 0; i < nlocal; i++) {
		MPI_File_seek(f, (myrank * nlocal * n + 1) * sizeof(int), MPI_SEEK_SET);
		MPI_File_write(f, &a[0], nlocal * n, MPI_INT, &status);
	}

	MPI_File_close(&f);
	free(a);

	MPI_Finalize();
	return 0;
}
