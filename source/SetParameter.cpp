#include "SetParameter.h"


using namespace std;


struct stat statDirectory;


SetParameter::SetParameter(int argc, char* argv[]) {

    profile_parameters["Ceiling"] = valarray<double> (5);
    profile_parameters["Floor"]   = valarray<double> (5);
    profile_parameters["Width"]   = valarray<double> (5);

    output_folder = string(".");

    switch (argc)
    {
        case 1:
            get_from_ui();
            break;

        case 2:
            get_from_file(argv[1]);
            break;

        default:
            get_from_args(argc, argv);
            break;
    }

    user_thickness_of_inner_wafer *= .0001;
    profile_parameters["Ceiling"] *= .0001;
    profile_parameters["Floor"]   *= .0001;
    profile_parameters["Width"]   *= .0001;

}


SetParameter::~SetParameter() {}


void SetParameter::get_from_ui() {

    string str_input;

    cout << "Input the Thickness of the wafer (unit: microns) ";
    getline(cin, str_input);
    stringstream(str_input) >> user_thickness_of_inner_wafer;
    cout << "\n\v";

    cout << "Now, please input the parameters for the profile mesh.\n"
            "Please refer to pictures and examples in the documentation.\n\n";

    for (int i = 0; i < 5; i++) {
        for (auto& name: parameter_names) {
            printf("Set the profile %s %d (microns): ", name, i + 1);
            getline(cin, str_input);
            stringstream(str_input) >> profile_parameters[name][i];
        }
    }

    cout << endl;

}


void SetParameter::get_from_file(char* file) {
    input = Input(file);
    user_thickness_of_inner_wafer = input.thickness;
    cylinder_folder = input.cylinder_folder;

//    CONTINUE HERE: test the run

    for (int i = 0; i < input.levels; ++i)
    {
        //                                       TODO look for another professional 
        //                                            code example with large input.
        //
        //                                       TODO unify the maps<> 
        //                                            SetParameter::profile_params with
        //                                            Input::prof_params
        profile_parameters["Ceiling"][i] = input.prof_params["profile_ceil_" + num_to_str(i)];
        profile_parameters["Floor"][i] = input.prof_params["profile_floor_" + num_to_str(i)];
        profile_parameters["Width"][i] = input.prof_params["profile_width_" + num_to_str(i)];
    }
}


void SetParameter::get_from_args(int argc, char* argv[]) {
    // TODO: this method will dissapear. Do not use this anymore. 
    
    argv_num_of_input = parameter_names.size()
                        * profile_layers_number
                        + extra_names.size();

    if (argc - 1 != argv_num_of_input) {
        cerr << "ERROR: Wrong number of arguments. There should be " << argv_num_of_input
             << " arguments for the profile." << endl;
        exit(1);
    }
 
    cylinder_folder = argv[argc - 1];

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
