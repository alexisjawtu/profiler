#include <math.h>
#include <iostream>
#include <string>
#include <array>
#include <valarray>
#include <map>
#include <sstream>


using namespace std;


class SetParameter
{
    public:

    string output_folder;
    double user_thickness_of_inner_wafer;  // WAS t

    map<string, valarray<double>> profile_parameters;

    SetParameter(int argc, char *argv[]);
    ~SetParameter();
    
    private:
    
    int argv_num_of_input;

    void get_from_ui();
    void get_from_args(int argc, char *argv[]);
    void check_output_folder(string folder);

};
