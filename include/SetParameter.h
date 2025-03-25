#include <math.h>
#include <iostream>
#include <string>
#include <array>
#include <valarray>
#include <map>
#include <sstream>

#include "Input.h"


using namespace std;


class SetParameter {
    public:

        string output_folder;
        double user_thickness_of_inner_wafer;  // WAS t

        map<string, valarray<double>> profile_parameters;

        SetParameter(int argc, char* argv[]);
        ~SetParameter();

    private:

        Input input;
        int profile_layers_number = 5;
        int argv_num_of_input;

        array<const char*, 2> extra_names{"input_location", "thickness"};
        array<const char*, 3> parameter_names{"Ceiling", "Floor", "Width"};

        void get_from_ui();
        void get_from_file(string file);
        void get_from_args(int argc, char* argv[]);
        void check_output_folder(string folder);
};
