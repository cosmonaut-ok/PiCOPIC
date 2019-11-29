#include "lib.hpp"
#include "outEngine.hpp"
#include "H5Cpp.h"

using namespace std;

class OutEngineHDF5 : public OutEngine
{
public:
  OutEngineHDF5 () {};
  OutEngineHDF5 (string a_path, string a_subpath, unsigned int a_shape, int *a_size,
                  bool a_append, bool a_compress, unsigned int a_compress_level);

  void write_rec(string a_name, Grid<double> data);
  void write_vec(string a_name, Grid<double> data);
  void write_dot(string a_name, Grid<double> data);
  void write_1d_vector(string a_name, vector<double> data);
};
