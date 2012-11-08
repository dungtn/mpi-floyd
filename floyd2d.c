/*
 ======================================================================================
 Name        : floyd2d.c
 Author      : daisy
 Version     : 0.1
 Description : MPI 2D all-pairs shortest path Floyd-Warshall algorithm in C, Ansi-style
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
 ======================================================================================
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

int floyd_all_pairs_sp_2d(int n, int nlocal, int *a, MPI_Comm comm_2d, MPI_Comm comm_row, MPI_Comm comm_col) {
	int i, j, k;
	int my2drank, mycoords[2];
	int row_root, col_root, coords[2];

	MPI_Comm_rank(comm_2d, &my2drank);
	MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);

	int *krow = (int *)malloc(nlocal * sizeof(int));
	int *kcol = (int *)malloc(nlocal * sizeof(int));

	for(k = 0; k < n; k++) {
		coords[ROW] = k / nlocal;
		coords[COL] = k / nlocal;
		if(k >= mycoords[ROW] * nlocal && k <= (mycoords[ROW] + 1) * nlocal) {
			for(i = 0; i < nlocal; i++) {
				krow[i] = a[(k % nlocal)*nlocal + i];
			}
		}
		MPI_Cart_rank(comm_col, coords, &col_root);
//		printf("[%d] column root %d\n", my2drank, col_root);
		MPI_Barrier(comm_col);
		MPI_Bcast(&krow[0], nlocal, MPI_INT, col_root, comm_col);
		if(k >= mycoords[COL] * nlocal && k <= (mycoords[COL] + 1) * nlocal) {
			for(i = 0; i < nlocal; i++) {
				kcol[i] = a[i*nlocal + (k % nlocal)];
			}
		}
		MPI_Cart_rank(comm_row, coords, &row_root);
//		printf("[%d] row root %d\n", my2drank, row_root);
		MPI_Barrier(comm_row);
		MPI_Bcast(&kcol[0], nlocal, MPI_INT, row_root, comm_row);

		for(i = 0; i < nlocal; i++) {
			for(j = 0; j < nlocal; j++) {
				a[i * nlocal + j] = min(a[i * nlocal + j], kcol[i] + krow[j]);
				if(my2drank == 3) {
					printf("k=%d kcol[%d]=%d krow[%d]=%d a[%d]=%d\n", k, i, kcol[i], j, krow[j], i * nlocal + j, a[i * nlocal + j]);
				}
			}
		}
	}
	free(krow);
	free(kcol);
	return 0;
}

int main(int argc, char *argv[]) {

	int i, n, nlocal;
	int numprocs, dims[2], periods[2], keep_dims[2];
	int myrank, my2drank, mycoords[2];
	MPI_File f; char* filename = "input/16";
	MPI_Comm comm_2d, comm_row, comm_col;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	dims[ROW] = dims[COL] = sqrt(numprocs);

	periods[ROW] = periods[COL] = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm_2d);

	MPI_Comm_rank(comm_2d, &my2drank);
	MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);

	keep_dims[ROW] = 0;
	keep_dims[COL] = 1;
	MPI_Cart_sub(comm_2d, keep_dims, &comm_row);

	keep_dims[ROW] = 1;
	keep_dims[COL] = 0;
	MPI_Cart_sub(comm_2d, keep_dims, &comm_col);

	if(MPI_File_open(comm_2d, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &f) != MPI_SUCCESS) {
		fprintf(stderr, "Cannot open file %s\n", filename);
		MPI_Abort(comm_2d, FILE_NOT_FOUND);
		MPI_Finalize();
		return 1;
	}
	MPI_File_seek(f, 0, MPI_SEEK_SET);
	MPI_File_read(f, &n, 1, MPI_INT, &status); nlocal = n/dims[ROW];

	int *a = (int *)malloc(nlocal * nlocal * sizeof(int));
	for(i = 0; i < nlocal; i++) {
		MPI_File_seek(f, ((mycoords[0] * nlocal  + i) * n + mycoords[1] * nlocal + 1) * sizeof(int), MPI_SEEK_SET);
		MPI_File_read(f, &a[i * nlocal], nlocal, MPI_INT, &status);
	}
	MPI_File_close(&f);

 	int j;
	if(my2drank == 3) {
		for(i = 0; i < nlocal; i++) {
			for(j = 0; j < nlocal; j++) {
				printf("%d ", a[i * nlocal +j]);
			}
			printf("\n");
		}
	}

	double start = MPI_Wtime();
	floyd_all_pairs_sp_2d(n, nlocal, a, comm_2d, comm_row, comm_col);
	double stop = MPI_Wtime();
	printf("[%d] Completed in %1.3f seconds\n", my2drank, stop-start);

	MPI_Comm_free(&comm_col);
	MPI_Comm_free(&comm_row);
	if(my2drank == 3) {
		for(i = 0; i < nlocal; i++) {
			for(j = 0; j < nlocal; j++) {
				printf("%d ", a[i * nlocal +j]);
			}
			printf("\n");
		}
	}
	if(MPI_File_open(comm_2d, "output/16", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &f) != MPI_SUCCESS) {
			printf("Cannot open file %s\n", "out");
			MPI_Abort(comm_2d, FILE_NOT_FOUND);
			MPI_Finalize();
			return 1;
	}
	if(my2drank == 0) {
		MPI_File_seek(f, 0, MPI_SEEK_SET);
		MPI_File_write(f, &n, 1, MPI_INT, &status);
	}
	for(i = 0; i < nlocal; i++) {
		MPI_File_seek(f, ((mycoords[0] * nlocal  + i) * n + mycoords[1] * nlocal + 1) * sizeof(int), MPI_SEEK_SET);
		MPI_File_write(f, &a[i * nlocal], nlocal, MPI_INT, &status);
	}

	MPI_File_close(&f);
	free(a);

	MPI_Comm_free(&comm_2d);
	MPI_Finalize();
	return 0;
}
