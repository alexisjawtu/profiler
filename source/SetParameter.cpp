#include <sys/stat.h>

#include "SetParameter.h"


using namespace std;


struct stat statDirectory;


SetParameter::SetParameter(int argc, char* argv[]) {

    output_folder     = string(".");
    argv_num_of_input = 17;

    if (argc == 1)
        get_from_ui();
    else
        get_from_args(argc, argv);

    user_thickness_of_inner_wafer *= .0001;
    profile_parameters["Ceiling"] *= .0001;
    profile_parameters["Floor"]   *= .0001;
    profile_parameters["Width"]   *= .0001;

}


SetParameter::~SetParameter() {}


void SetParameter::get_from_ui() {

    string str_input;

    cout << "Input the Thickness of the wafer (unit: microns) ";
    cin >> user_thickness_of_inner_wafer;
    cout << endl << endl;

    cout << "Now, please input the parameters for the profile mesh.\n"
            "Please refer to pictures and examples in the documentation." << endl << endl;

    array<const char*, 3> parameter_names ({"Ceiling", "Floor", "Width"});

    profile_parameters["Ceiling"] = valarray<double> (5);
    profile_parameters["Floor"]   = valarray<double> (5);
    profile_parameters["Width"]   = valarray<double> (5);

    for (int i = 0; i < 5; i++) {
        for (auto &name: parameter_names) {
            printf("Set the profile %s %d (microns): ", name, i + 1);
            getline(cin, str_input);
            stringstream(str_input) >> profile_parameters[name][i];
        }
    }

    cout << endl;

}


void SetParameter::get_from_args(int argc, char* argv[]) {

    cout << "argc " << argc << endl;

    // Error in number of args
    if (argc != argv_num_of_input) {
        cerr << "ERROR: Wrong number of arguments. There should be " << argv_num_of_input
             << " arguments for the profile." << endl;
        exit(1);
    }
    
    user_thickness_of_inner_wafer = atof(argv[1]);

    profile_parameters["Ceiling"] = valarray<double> {
        atof(argv[2]), 
        atof(argv[5]),
        atof(argv[8]),
        atof(argv[11]),
        atof(argv[14])
    };

    profile_parameters["Floor"] = valarray<double> {
        atof(argv[3]),
        atof(argv[6]),
        atof(argv[9]),
        atof(argv[12]),
        atof(argv[15])
    };

    profile_parameters["Width"] = valarray<double> {
        atof(argv[4]),
        atof(argv[7]),
        atof(argv[10]),
        atof(argv[13]),
        atof(argv[16])
    };

}


void SetParameter::check_output_folder(string folder) {
    
    if (false) {

        // TODO check existence of the folder, to avoid overwriting
        // if (stat(prefix.c_str(), &statDirectory) == -1) {
        //     cerr << "Input directory does not exists!\n" ;
        //     exit(1);
        // };

    } else {

        output_folder = folder;
        mkdir(output_folder.data(), 0777);

    }

}
