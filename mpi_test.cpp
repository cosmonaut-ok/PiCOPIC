#include <vector>
#include <mpi.h>

using namespace std;
using namespace MPI;

struct Child
{
  int foo;
  std::vector<float> bar;
  std::vector<int> baz;

  Child() : dtype(MPI_DATATYPE_NULL) {}
  // ~Child()
  // {
  //   if (dtype != MPI_DATATYPE_NULL)
  //     MPI_Datatype_free(&dtype);
  // }

  const MPI_Datatype mpi_dtype();
  void invalidate_dtype();

private:
  MPI_Datatype dtype;
  void make_dtype();
};

const MPI_Datatype Child::mpi_dtype()
{
  if (dtype == MPI_DATATYPE_NULL)
    make_dtype();
  return dtype;
}

void Child::invalidate_dtype()
{
  // if (dtype != MPI_DATATYPE_NULL)
  //   MPI_Datatype_free(&dtype);
}

void Child::make_dtype()
{
  const int nblock = 3;
  int block_count[nblock] = {1, bar.size(), baz.size()};
  MPI_Datatype block_type[nblock] = {MPI_INT, MPI_FLOAT, MPI_INT};
  MPI_Aint offset[nblock];
  MPI_Get_address(&foo, &offset[0]);
  MPI_Get_address(&bar[0], &offset[1]);
  MPI_Get_address(&baz[0], &offset[2]);

  MPI_Type_struct(nblock, block_count, offset, block_type, &dtype);
  MPI_Type_commit(&dtype);
}

int main ()
{
  int rank, size;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  Child kid;
  kid.foo = 5;
  kid.bar.resize(5);
  kid.baz.resize(10);

  Child kid2;

  if (rank == 0)
  {
    MPI_Send(MPI_BOTTOM, 1, kid.mpi_dtype(), 1, 0, MPI_COMM_WORLD);
  }

  if (rank == 1)
  {
    kid2.foo = 10;
    kid.bar.resize(5);
    kid.baz.resize(10);
    std::cout << kid2.foo << std::endl;
    MPI_Recv(MPI_BOTTOM, 1, kid2.mpi_dtype(), 0, 0, MPI_COMM_WORLD, NULL);
    std::cout << kid2.foo << std::endl;
  }

  MPI_Finalize();
}
