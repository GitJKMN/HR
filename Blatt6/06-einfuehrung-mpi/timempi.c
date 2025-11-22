#define _DEFAULT_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include "mpi.h"

int main(int argc, char *argv[]) {


  int rank, size;
  MPI_Init( &argc, &argv ); 
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  
  int master_rank = size - 1;
  int is_master = rank == master_rank;

  if (is_master) {
    int ms[master_rank];
    int min_ms;
    int max_ms;
    for (int i = 0; i < master_rank; i++) {
      char recv[80];
      MPI_Recv(&recv, 80, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&ms[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("%s\n", recv);
    }
    min_ms = ms[0];
    max_ms = ms[0];
    for (int i = 0; i < master_rank; i++) {
      if (ms[i] < min_ms) {
        min_ms = ms[i];
      } else if (ms[i] > max_ms) {
        max_ms = ms[i];
      }
    }
    printf("[%d] Kleinster MS-Anteil: %d\n",rank, min_ms);
    printf("[%d] Größte Differenz: %d\n", rank, max_ms - min_ms);
  } else {
    struct timeval tv;
    time_t current_time;
    int micro_sec;
    char time_string[30];
    char output[80];
    char hostname[30];

    gettimeofday(&tv, NULL);
    gethostname(hostname, 30);

    current_time = tv.tv_sec;
    micro_sec = tv.tv_usec;

    strftime(time_string, 30, "%Y-%m-%d %T", localtime(&current_time));
    snprintf(output, 80, "[%d] %s // %s.%d", rank, hostname, time_string, micro_sec);
    MPI_Send(&output, 80, MPI_CHAR, master_rank, 0, MPI_COMM_WORLD);
    MPI_Send(&micro_sec, 1, MPI_INT, master_rank, 0, MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  printf("[%d] beendet jetzt!\n", rank);
  return 0;
}