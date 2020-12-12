// Author: Wes Kendall
// Copyright 2011 www.mpitutorial.com
// This code is provided freely with the tutorials on mpitutorial.com. Feel
// free to modify it for your own use. Any distribution of the code must
// either provide a link to www.mpitutorial.com or keep this header intact.
//
// MPI_Send, MPI_Recv example. Communicates the number -1 from process 0
// to process 1.
//
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <time.h>       /* time */
#include <vector>

using namespace std;

int main(int argc, char** argv)
{
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // We are assuming at least 2 processes for this task
  if (world_size < 4)
  {
    cerr << "World size must be 4 for " << argv[0] << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  unsigned int how, howrcv;

  // setup random number of vectors with random number of values
  unsigned int vec_rank = 15;

  srand (time(NULL));
  vector< vector<double> * > particles;

  for (unsigned int ps = 0; ps < 100; ++ps)
  {
    // int rnd_for_nodes = rand() % 4 + 1; //
    vector<double> *n = new vector<double>(vec_rank, 0);
    for (unsigned int v=0; v < vec_rank; ++v)
    {
      double rnd_for_nodes = (double)(rand() % vec_rank + 1); //
      (*n)[v] = rnd_for_nodes;
    }
    particles.push_back(n);
  }
  cout << "number of particles is " << particles.size() << endl;

  for (unsigned int nod = 0; nod < world_size; ++nod)
    if (nod != world_rank)
    {
      int rnd_for_send = rand() % 10 + 1;
      vector< vector<double> * > p_to_send;

      for (unsigned int i = 0; i < rnd_for_send; ++i)
        p_to_send.push_back(particles[i]);

      int how = p_to_send.size();

      MPI_Send (
        /* data         = */ &how,
        /* count        = */ 1,
        /* datatype     = */ MPI_UNSIGNED,
        /* destination  = */ nod,
        /* tag          = */ 0,
        /* communicator = */ MPI_COMM_WORLD);

      for (unsigned prtl = 0; prtl < how; ++prtl)
      {
        double dbuff[vec_rank];
        for (unsigned int val = 0; val < vec_rank; ++val)
          dbuff[val] = (*p_to_send[prtl])[val];

        // send then receive
        MPI_Send (
          /* data         = */ dbuff,
          /* count        = */ vec_rank,
          /* datatype     = */ MPI_DOUBLE,
          /* destination  = */ nod,
          /* tag          = */ 0,
          /* communicator = */ MPI_COMM_WORLD);
        cout << "Sent particle " << prtl << " to node " << nod << endl;
      }
    }

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  for (unsigned int nod = 0; nod < world_size; ++nod)
    if (nod != world_rank)
    {
      MPI_Recv(
        /* data         = */ &howrcv,
        /* count        = */ 1,
        /* datatype     = */ MPI_UNSIGNED,
        /* source       = */ nod,
        /* tag          = */ 0,
        /* communicator = */ MPI_COMM_WORLD,
        /* status       = */ MPI_STATUS_IGNORE);

      // for (unsigned prtl = 0; prtl < how; ++prtl)
      // {
      //   double dbuff[vec_rank];
      //   // for (unsigned int val = 0; val < vec_rank; ++val)
      //   //   dbuff[val] = (*p_to_send[prtl])[val];

      for (unsigned int i = 0; i < howrcv; ++i)
      {
        double dbuff[vec_rank];

        MPI_Recv(
          /* data         = */ dbuff,
          /* count        = */ vec_rank,
          /* datatype     = */ MPI_DOUBLE,
          /* source       = */ nod,
          /* tag          = */ 0,
          /* communicator = */ MPI_COMM_WORLD,
          /* status       = */ MPI_STATUS_IGNORE);

        cout << "Received particle " << i << " from node " << nod << endl;

        vector<double> *n = new vector<double>(vec_rank, 0);
        for (unsigned int v=0; v < vec_rank; ++v)
          (*n)[v] = dbuff[v];

        particles.push_back(n);
      }
    }
  MPI_Barrier(MPI_COMM_WORLD);
  cout << "Particles size is " << particles.size() << " for node " << world_rank << endl;

  // if (world_rank == 0)
  // {
  //   // howmuch
  //   how = pair_vec.size();
  //   double dbuff[how];
  //   for (unsigned int i=0; i < how; ++i)
  //     dbuff[i] = pair_vec[i];

  //   MPI_Send(
  //     /* data         = */ &how,
  //     /* count        = */ 1,
  //     /* datatype     = */ MPI_UNSIGNED,
  //     /* destination  = */ 1,
  //     /* tag          = */ 0,
  //     /* communicator = */ MPI_COMM_WORLD);

  //   cout << "sent number " << how << endl;

  //   // send then receive
  //   MPI_Send(
  //     /* data         = */ dbuff,
  //     /* count        = */ how,
  //     /* datatype     = */ MPI_DOUBLE,
  //     /* destination  = */ 1,
  //     /* tag          = */ 0,
  //     /* communicator = */ MPI_COMM_WORLD);

  //   // MPI_Barrier(MPI_COMM_WORLD);
  //   // cout << "-------------------- BOOM --------------------" << endl;

  //   // MPI_Recv(
  //   //   /* data         = */ &how,
  //   //   /* count        = */ 1,
  //   //   /* datatype     = */ MPI_DOUBLE,
  //   //   /* source       = */ 1,
  //   //   /* tag          = */ 0,
  //   //   /* communicator = */ MPI_COMM_WORLD,
  //   //   /* status       = */ MPI_STATUS_IGNORE);

  //   // cout << "received number  " << how << endl;

  //   // vector<double> recv_vec(how, 0);
  //   // for (unsigned int i = 0; i < how; ++i)
  //   // {
  //   //   MPI_Recv(
  //   //     /* data         = */ &recv_vec[i],
  //   //     /* count        = */ 1,
  //   //     /* datatype     = */ MPI_DOUBLE,
  //   //     /* source       = */ 1,
  //   //     /* tag          = */ 0,
  //   //     /* communicator = */ MPI_COMM_WORLD,
  //   //     /* status       = */ MPI_STATUS_IGNORE);

  //   //   cout << "received " << recv_vec[i] << endl;
  //   // }
  // }

  // else
  // {
  //   MPI_Recv (
  //     /* data         = */ &how,
  //     /* count        = */ 1,
  //     /* datatype     = */ MPI_DOUBLE,
  //     /* source       = */ 0,
  //     /* tag          = */ 0,
  //     /* communicator = */ MPI_COMM_WORLD,
  //     /* status       = */ MPI_STATUS_IGNORE);

  //   cout << "received number  " << how << endl;

  //   double dbuff[how];
  //   vector<double> recv_vec(how, 0);
  //   MPI_Recv (
  //     /* data         = */ dbuff,
  //     /* count        = */ how,
  //     /* datatype     = */ MPI_DOUBLE,
  //     /* source       = */ 0,
  //     /* tag          = */ 0,
  //     /* communicator = */ MPI_COMM_WORLD,
  //     /* status       = */ MPI_STATUS_IGNORE);


  //   for (unsigned int i=0; i < how; ++i)
  //     recv_vec[i] = dbuff[i];

  //   for (unsigned int i=0; i < how; ++i)
  //     cout << "received " << recv_vec[i] << endl;

  //   MPI_Barrier(MPI_COMM_WORLD);
  //   cout << "-------------------- MOOB --------------------" << endl;
  //   MPI_Barrier(MPI_COMM_WORLD);

  //   // // howmuch
  //   // how = unpair_vec.size();
  //   // MPI_Send (
  //   //   /* data         = */ &how,
  //   //   /* count        = */ 1,
  //   //   /* datatype     = */ MPI_UNSIGNED,
  //   //   /* destination  = */ 0,
  //   //   /* tag          = */ 0,
  //   //   /* communicator = */ MPI_COMM_WORLD);

  //   // cout << "sent number " << how << endl;

  //   // // send then receive
  //   // for (unsigned int i = 0; i < how; ++i)
  //   // {

  //   //   MPI_Send(
  //   //     /* data         = */ &unpair_vec[i],
  //   //     /* count        = */ 1,
  //   //     /* datatype     = */ MPI_DOUBLE,
  //   //     /* destination  = */ 0,
  //   //     /* tag          = */ 0,
  //   //     /* communicator = */ MPI_COMM_WORLD);
  //   //   cout << "sent  " << unpair_vec[i] << endl;
  //   // }
  // }
  MPI_Finalize();
}



//// EXAMPLE WITH STRUCTS
// #include <stdio.h>
// #include <stdlib.h>
// #include <mpi.h>
// #include <stddef.h>

// typedef struct car_s {
//         int shifts;
//         int topSpeed;
// } car;

// int main(int argc, char **argv) {

//     const int tag = 13;
//     int size, rank;

//     MPI_Init(&argc, &argv);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     if (size < 2) {
//         fprintf(stderr,"Requires at least two processes.\n");
//         exit(-1);
//     }

//     /* create a type for struct car */
//     const int nitems=2;
//     int          blocklengths[2] = {1,1};
//     MPI_Datatype types[2] = {MPI_INT, MPI_INT};
//     MPI_Datatype mpi_car_type;
//     MPI_Aint     offsets[2];

//     offsets[0] = offsetof(car, shifts);
//     offsets[1] = offsetof(car, topSpeed);

//     MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_car_type);
//     MPI_Type_commit(&mpi_car_type);

//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     if (rank == 0) {
//         car send;
//         send.shifts = 4;
//         send.topSpeed = 100;

//         const int dest = 1;
//         MPI_Send(&send,   1, mpi_car_type, dest, tag, MPI_COMM_WORLD);

//         printf("Rank %d: sent structure car\n", rank);
//     }
//     if (rank == 1) {
//         MPI_Status status;
//         const int src=0;

//         car recv;

//         MPI_Recv(&recv,   1, mpi_car_type, src, tag, MPI_COMM_WORLD, &status);
//         printf("Rank %d: Received: shifts = %d topSpeed = %d\n", rank,
//                  recv.shifts, recv.topSpeed);
//     }

//     MPI_Type_free(&mpi_car_type);
//     MPI_Finalize();

//     return 0;
// }


//// ANOTHER EXAMPLE
// #include "mpi.h"
// #include <stdio.h>

// typedef struct 
// {
//     char a;
//     int b;
//     short c;
// } my_struct;


// int main (int argc, char *argv[])
// {
//     int  numtasks, taskid;

//     MPI_Init(&argc, &argv);
//     MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
//     MPI_Comm_size(MPI_COMM_WORLD, &numtasks);


//     if (taskid == 0) 
//     {
//         my_struct m;
//         m.a = '!';
//         m.b = 1234;
//         m.c = 5678;

//         MPI_Send(&m, sizeof(my_struct), MPI_CHAR, 1, 0, MPI_COMM_WORLD);
//     }
//     else 
//     {
//         my_struct m;
//         MPI_Recv(&m, sizeof(my_struct), MPI_CHAR, 0, 0, MPI_COMM_WORLD, 
//                  MPI_STATUS_IGNORE);
//         printf("%c %d %d\n", m.a, m.b, m.c); 
//     }

//     MPI_Finalize();
// }



//// SEND VECTOR VIA MPI

// #include <iostream>
// #include <mpi.h>

// struct Point {
//   double x, y, z;
// };

// int main(int argc, char **argv) {
//   MPI_Init(&argc, &argv);

//   int rank, size;
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   MPI_Comm_size(MPI_COMM_WORLD, &size);

//   MPI_Datatype dt_point;
//   MPI_Type_contiguous(3, MPI_DOUBLE, &dt_point);
//   MPI_Type_commit(&dt_point);
  
//   constexpr int n_points = 10;
//   Point data[n_points];

//   // Process 0 sends the data
//   if (rank == 0) {
//     for (int i=0; i < n_points; ++i) {
//       data[i].x = (double)i;
//       data[i].y = (double)-i;
//       data[i].z = (double) i * i;
//     }

//     MPI_Send(data, n_points, dt_point, 1, 0, MPI_COMM_WORLD);
//   }
//   else { // Process 1 receives
//     MPI_Recv(data, n_points, dt_point, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

//     // Printing
//     for (int i=0; i < n_points; ++i) {
//       std::cout << "Point #" << i << " : (" << data[i].x << "; " << data[i].y << "; " << data[i].z << ")"
// 		<< std::endl;
//     }
//   }

//   MPI_Finalize();
// }



//// OPERATE WITH CUSTOM DATA TYPE
// #include <iostream>
// #include <cstdlib>
// #include <cmath>
// #include <mpi.h>

// constexpr int DOUBLE_MAX = 10;
// struct CustomData {
//   int n_values;
//   double values[DOUBLE_MAX];
// };

// int main(int argc, char **argv) {
  
//   MPI_Init(&argc, &argv);

//   int rank, size;

//   MPI_Comm_size(MPI_COMM_WORLD, &size);
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//   constexpr int n_structure_per_process = 5; // M = 5

//   // Random generator init
//   srand(rank * 10);
  
//   // Creating the dataset
//   CustomData data[n_structure_per_process];

//   // Generating the data
//   for (int i=0; i < n_structure_per_process; ++i) {
//     // Terrible way of generating random numbers, don't reproduce this at home
//     data[i].n_values = rand() % DOUBLE_MAX + 1;
//     for (int j=0; j < DOUBLE_MAX; ++j)
//       data[i].values[j] = (j < data[i].n_values ? (double)rand() / (double)RAND_MAX : 0.0);
//   }

//   // Copying the data to two different arrays
//   int int_send_buf[n_structure_per_process];
//   double dbl_send_buf[n_structure_per_process * DOUBLE_MAX];

//   for (int i=0; i < n_structure_per_process; ++i) {
//     int_send_buf[i] = data[i].n_values;
//     for (int j=0; j < data[i].n_values; ++j)
//       dbl_send_buf[i*DOUBLE_MAX + j] = data[i].values[j];
//   }

//   // Gathering everything on process 0
//   int *n_values = nullptr; 
//   double *dbl_values = nullptr; 

//   if (rank == 0) {
//     n_values = new int[n_structure_per_process * size];
//     dbl_values = new double[n_structure_per_process * size * DOUBLE_MAX];
//   }
  
//   MPI_Gather(int_send_buf, n_structure_per_process, MPI_INT,
// 	     n_values, n_structure_per_process, MPI_INT, 0, MPI_COMM_WORLD);
//   MPI_Gather(dbl_send_buf, n_structure_per_process * DOUBLE_MAX, MPI_DOUBLE,
// 	     dbl_values, n_structure_per_process * DOUBLE_MAX, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//   if (rank == 0) {
//     // Recopying the data and printing
//     CustomData gathered_data[n_structure_per_process * size];
//     memset(gathered_data, 0, n_structure_per_process * size * sizeof(CustomData));
//     for (int i=0; i < size; ++i) {
//       for (int j=0; j < n_structure_per_process; ++j) {
// 	int data_id = i * n_structure_per_process + j; // Linear index

// 	std::cout << "Data structure " << data_id << " : [";
	
// 	gathered_data[data_id].n_values = n_values[data_id];
// 	for (int k=0; k < n_values[data_id]; ++k) {
// 	  gathered_data[data_id].values[k] = dbl_values[i*n_structure_per_process*DOUBLE_MAX + j*DOUBLE_MAX + k];
// 	  std::cout << gathered_data[data_id].values[k] << (k == n_values[data_id]-1 ? "]" : "; ");
// 	}
// 	std::cout << std::endl;
//       }
//     }

//     // And freeing the memory
//     delete [] n_values;
//     delete [] dbl_values;
//   }
  
  
//   MPI_Finalize();
  
//   return 0;
// }

//// OPERATE WITH CUSTOM DATA TYPE 2
// #include <stdio.h>
// #include <stdlib.h>
// #include <mpi.h>
// #include <stddef.h>

// typedef struct car_s {
//         int shifts;
//         int topSpeed;
// } car;

// int main(int argc, char **argv) {

//     const int tag = 13;
//     int size, rank;

//     MPI_Init(&argc, &argv);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     if (size < 2) {
//         fprintf(stderr,"Requires at least two processes.\n");
//         exit(-1);
//     }

//     /* create a type for struct car */
//     const int nitems=2;
//     int          blocklengths[2] = {1,1};
//     MPI_Datatype types[2] = {MPI_INT, MPI_INT};
//     MPI_Datatype mpi_car_type;
//     MPI_Aint     offsets[2];

//     offsets[0] = offsetof(car, shifts);
//     offsets[1] = offsetof(car, topSpeed);

//     MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_car_type);
//     MPI_Type_commit(&mpi_car_type);

//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     if (rank == 0) {
//         car send;
//         send.shifts = 4;
//         send.topSpeed = 100;

//         const int dest = 1;
//         MPI_Send(&send,   1, mpi_car_type, dest, tag, MPI_COMM_WORLD);

//         printf("Rank %d: sent structure car\n", rank);
//     }
//     if (rank == 1) {
//         MPI_Status status;
//         const int src=0;

//         car recv;

//         MPI_Recv(&recv,   1, mpi_car_type, src, tag, MPI_COMM_WORLD, &status);
//         printf("Rank %d: Received: shifts = %d topSpeed = %d\n", rank,
//                  recv.shifts, recv.topSpeed);
//     }

//     MPI_Type_free(&mpi_car_type);
//     MPI_Finalize();

//     return 0;
// }


//// OPERATE WITH CUSTOM DATA TYPE 3
// #include <stdio.h>
// #include "mpi.h"

// int main( argc, argv )
// int argc;
// char **argv;
// {
//     int          rank;
//     struct { int a; double b;} value;
//     MPI_Datatype mystruct;
//     int          blocklens[2];
//     MPI_Aint     indices[2];
//     MPI_Datatype old_types[2];

//     MPI_Init( &argc, &argv );

//     MPI_Comm_rank( MPI_COMM_WORLD, &rank );

//     /* One value of each type */
//     blocklens[0] = 1;
//     blocklens[1] = 1;
//     /* The base types */
//     old_types[0] = MPI_INT;
//     old_types[1] = MPI_DOUBLE;
//     /* The locations of each element */
//     MPI_Address( &value.a, &indices[0] );
//     MPI_Address( &value.b, &indices[1] );
//     /* Make relative */
//     indices[1] = indices[1] - indices[0];
//     indices[0] = 0;
//     MPI_Type_struct( 2, blocklens, indices, old_types, &mystruct );
//     MPI_Type_commit( &mystruct );

//     do {
//         if (rank == 0)
//             scanf( "%d %lf", &value.a, &value.b );

//         MPI_Bcast( &value, 1, mystruct, 0, MPI_COMM_WORLD );

//         printf( "Process %d got %d and %lf
// ", rank, value.a, value.b );
//     } while (value.a >= 0);

//     /* Clean up the type */
//     MPI_Type_free( &mystruct );
//     MPI_Finalize( );
//     return 0;
// }

//// OPERATE WITH CUSTOM DATA TYPE 4
// #include "mpi.h"
// #include <stdio.h>

// struct Partstruct
// {
//     char c;
//     double d[6];
//     char b[7];
// };

// int main(int argc, char *argv[])
// {
//     struct Partstruct particle[1000];
//     int i, j, myrank;
//     MPI_Status status;
//     MPI_Datatype Particletype;
//     MPI_Datatype type[3] = { MPI_CHAR, MPI_DOUBLE, MPI_CHAR };
//     int blocklen[3] = { 1, 6, 7 };
//     MPI_Aint disp[3];
 
//     MPI_Init(&argc, &argv);
 
//     disp[0] = &particle[0].c - &particle[0];
//     disp[1] = &particle[0].d - &particle[0];
//     disp[2] = &particle[0].b - &particle[0];
//     MPI_Type_create_struct(3, blocklen, disp, type, &Particletype);
//     MPI_Type_commit(&Particletype);
 
//     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 
//     if (myrank == 0)
//     {
//         MPI_Send(particle, 1000, Particletype, 1, 123, MPI_COMM_WORLD);
//     }
//     else if (myrank == 1)
//     {
//         MPI_Recv(particle, 1000, Particletype, 0, 123, MPI_COMM_WORLD, &status);
//     }
//     MPI_Finalize();
//     return 0;
// }
