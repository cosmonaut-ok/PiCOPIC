#include "msg.hpp"

using namespace std;

namespace msg
{
  int printed_step = 0;

  void print_probe_data(int probe_type, char* component, int step_number, int dump_interval, int* shape, double current_time)
  {
    int print_header_step = 20;
    if (printed_step % print_header_step == 0)
      {

	cout << endl
	     << left << setw(8) << "Step"
	     << left << setw(13) << "Saved Frame"
	     << left << setw(20) << "Dumping Probe Name"
	     << left << setw(21) << "Shape"
	     << left << setw(18) << "Model Time (sec)"
	     << left << setw(21) << "Simulation Duration"
	     << endl;
      }

    char type_comp[100];
    char probe_shape[100];

    switch (probe_type)
      {
      case 0:
	sprintf(type_comp, "frame/%s", component);
	sprintf(probe_shape, "[%i,%i,%i,%i]", shape[0], shape[1], shape[2], shape[3]);
	break;
      case 1:
	sprintf(type_comp, "column/%s", component);
	sprintf(probe_shape, "[%i]", shape[1]);
	break;
      case 2:
	sprintf(type_comp, "row/%s", component);
	sprintf(probe_shape, "[%i]", shape[0]);
	break;
      case 3:
	sprintf(type_comp, "dot/%s", component);
	sprintf(probe_shape, "[%i,%i]", shape[0], shape[1]);
	break;
      case 4:
	sprintf(type_comp, "mpframe/%s", component);
	sprintf(probe_shape, "[%i,%i,%i,%i]", shape[0], shape[1], shape[2], shape[3]);
	break;
      }

    cout << left << setw(8) << step_number * dump_interval
	 << left << setw(13) << step_number
	 << left << setw(20) << type_comp
	 << left << setw(21) << probe_shape
	 << left << setw(18) << current_time
	 << left << setw(18) << lib::get_simulation_duration()
	 << endl;
    ++printed_step;
  }
}
