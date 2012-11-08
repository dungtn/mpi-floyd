/*
 ===============================================================================
 Name        : floyd.c
 Author      : dung.thai.ngoc at gmail dot com
 Version     : 0.1
 Description : All-pairs shortest path Floyd-Warshall algorithm in C, Ansi-style
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
 ===============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define INF 99

inline int min(int a, int b) {
  return a < b ? a : b;
}
inline int print(int d) {
	if(d == INF)
		printf("- ");
	else
		printf("%d ", d);
	return 0;
}
inline int fprint(FILE *f, int d) {
	if(d == INF)
		fprintf(f, "- ");
	else
		fprintf(f, "%d ", d);
	return 0;
}

int floyd_all_pairs_sp(int n, int **A) {
	int i, j, k;
	for(k = 0; k < n; k++) {
		for(i = 0; i < n; i++) {
			for(j = 0; j< n; j++) {
				if(i != j)  {
					if(A[i][j] == 0) A[i][j] = INF;
					A[i][j] = min(A[i][j], A[i][k] + A[k][j]);
				}
			}
		}
	}
	return 0;
}

double get_time() {
	struct timeval tv;
	struct timezone tz;

	gettimeofday(&tv, &tz);
	return (double)tv.tv_sec + ((double)tv.tv_usec/1000000.0);
}

int main(int argc, char * argv[]){

	if (argc != 2) {
		printf("Usage: %s filename", argv[0]);
		return 1;
	}

	FILE *f;
	if((f = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "Cannot open file %s\n", argv[1]);
		return 1;
	}

	int n; fscanf(f, "%d", &n);
	int i, j;
	int **A = (int **)malloc(n * sizeof(int*));
	for(i = 0; i < n; i++) {
		A[i] = (int *)malloc(n * sizeof(int));
		for(j = 0; j < n; j++) {
			fscanf(f, "%d ", &A[i][j]);
		}
	}
	fclose(f);

	double start = get_time();
	floyd_all_pairs_sp(n, A);
	double stop = get_time();
	printf("Completed in %1.3f seconds\n", stop-start);

	f = fopen("sp", "w");
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			fprint(f, A[i][j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);

	f = fopen("benchmark", "a");
	fprintf(f, "%d %d %1.3f", n, 1, stop-start);
	fclose(f);
	return EXIT_SUCCESS;
}
