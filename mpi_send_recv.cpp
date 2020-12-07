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
