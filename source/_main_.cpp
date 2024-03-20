#include <iostream>
#include <typeinfo>
#include <fstream>
#include <string>
#include <valarray>
#include <array>
#include <vector>
#include <map>
#include <algorithm>


using namespace std;


int main(int argc, char* argv[])
{
    
    bool bit;
    ifstream diagonals_data ("diags.dat");

    if (!diagonals_data)
    {
        cout << "Something is wrong with the diagonals input file." << endl;
        exit(1);
    }

    diagonals_data.seekg(0, diagonals_data.end);  // set position at the end
    int diagonals = diagonals_data.tellg();       // tell which position is the end
    diagonals_data.seekg(0, diagonals_data.beg);  // set position back to beginning

    cout << "number " << diagonals - 1<< endl;
    
    //diagonals_data.getline(_raw_data_, diagonals + 1);

    cout << boolalpha << (diagonals_data.get() == '1') << endl;
    cout << boolalpha << (diagonals_data.get() == '1') << endl;
    cout << boolalpha << (diagonals_data.get() == '1') << endl;
    cout << boolalpha << (diagonals_data.get() == '1') << endl;
    
    cout << boolalpha << (diagonals_data.get() == '1') << endl;
    cout << boolalpha << (diagonals_data.get() == '1') << endl;
    cout << boolalpha << (diagonals_data.get() == '1') << endl;
    cout << boolalpha << (diagonals_data.get() == '1') << endl;

    cout << boolalpha << (diagonals_data.get() == '1') << endl;
    cout << boolalpha << (diagonals_data.get() == '1') << endl;
    cout << boolalpha << (diagonals_data.get() == '1') << endl;
    cout << boolalpha << (diagonals_data.get() == '1') << endl;

    cout << boolalpha << (diagonals_data.get() == '1') << endl;

    return 0;

}

/*
vector<size_t> triangles (triangles_0);
    auto print = [&](int id)
    {
        std::cout << "@" << id << ": ";
        for (int i : triangles_0)
            std::cout << i << ",\n";
    };

    sort(triangles.begin(), triangles.end());

    auto last = unique(triangles.begin(), triangles.end());

    triangles.erase(last, triangles.end());



    print(1);

// std::valarray<double> test = std::valarray<double> (5);
    // std::valarray<double> test2 = std::valarray<double> {0.31, 2, 3, -3.4};
    
    // std::array<double, 3> test3 = {3, 5, -3.2};
    // std::vector<double> test4 = {3, 5, -3.2, 9};

    // test2 *= 2;

    // for (int j = 0; j < 4; j++) {
    //     std::cout << test2[j] << std::endl;
    // }

    map<string, valarray<double>> profile_parameters;
    map<string, vector<string>> ui_messages;
    array<string, 3> parameter_names ({"ceilings", "floors", "widths"});

    profile_parameters["ceilings"] = valarray<double> (5);
    profile_parameters["floors"]   = valarray<double> (5);
    profile_parameters["widths"]   = valarray<double> (5);

    // TODO we could only iterate {"Ceiling", "Floor", "Width"} here!
    ui_messages["ceilings"] = vector<string> ({
        "Set the profile Ceiling 1 (microns): ",
        "Set the profile Ceiling 2 (microns): ",
        "Set the profile Ceiling 3 (microns): ",
        "Set the profile Ceiling 4 (microns): ",
        "Set the profile Ceiling 5 (microns): "
    });

    ui_messages["floors"] = vector<string> ({
        "Set the profile Floor 1 (microns): ",
        "Set the profile Floor 2 (microns): ",
        "Set the profile Floor 3 (microns): ",
        "Set the profile Floor 4 (microns): ",
        "Set the profile Floor 5 (microns): "
    });

    ui_messages["widths"] = vector<string> ({
        "Set the profile Width 1 (microns): ",
        "Set the profile Width 2 (microns): ",
        "Set the profile Width 3 (microns): ",
        "Set the profile Width 4 (microns): ",
        "Set the profile Width 5 (microns): "
    });

    for (int i = 0; i < 5; i++) {
        for (auto &name: parameter_names) {
            cout << ui_messages[name][i];
            cin >> profile_parameters[name][i];
        }
    }

    for (int i = 0; i < 5; i++) {
        for (auto &name: parameter_names) {
            cout << profile_parameters[name][i] << endl;
        }
    }



{
    char char_array[265];

    int index;

    std::ifstream data ("parameters.csv");

    data.getline(char_array, 256);
    data.getline(char_array, 256);

    
        take the first  --> folder
        take the second --> type
        then: parameters_tail

    std::string wafer_params (char_array);
    // parameters from the "scale" to the last 
    std::map<std::string, double> parameters_tail;  

    index = wafer_params.find_first_of(",");
    std::string folder  = wafer_params.substr(0, index);
    std::string numeric_parameters = wafer_params.substr(index + 1);

    index = numeric_parameters.find_first_of(",");
    int type = std::stoi(numeric_parameters.substr(0, index));
    numeric_parameters = numeric_parameters.substr(index + 1);

    index = numeric_parameters.find_first_of(",");
    parameters_tail.emplace("scale", std::stod(numeric_parameters.substr(0, index)));
    numeric_parameters = numeric_parameters.substr(index + 1);
    
    index = numeric_parameters.find_first_of(",");
    parameters_tail.emplace("diameter", std::stod( numeric_parameters.substr(0, index)  )  );
    numeric_parameters = numeric_parameters.substr(index + 1);

    index = numeric_parameters.find_first_of(",");
    parameters_tail.emplace("thickness", std::stod(numeric_parameters.substr(0, index)));
    numeric_parameters = numeric_parameters.substr(index + 1);

    index = numeric_parameters.find_first_of(",");
    parameters_tail.emplace("r_c", std::stod(numeric_parameters.substr(0, index)));
    numeric_parameters = numeric_parameters.substr(index + 1);

    index = numeric_parameters.find_first_of(",");
    parameters_tail.emplace("theta", std::stod(numeric_parameters.substr(0, index)));
    numeric_parameters = numeric_parameters.substr(index + 1);

    index = numeric_parameters.find_first_of(",");
    parameters_tail.emplace("h", std::stod(numeric_parameters.substr(0, index)));
    numeric_parameters = numeric_parameters.substr(index + 1);

    index = numeric_parameters.find_first_of(",");
    parameters_tail.emplace("dEA", std::stod(numeric_parameters.substr(0, index)));
    numeric_parameters = numeric_parameters.substr(index + 1);

    index = numeric_parameters.find_first_of(",");
    parameters_tail.emplace("dEB", std::stod(numeric_parameters.substr(0, index)));
    numeric_parameters = numeric_parameters.substr(index + 1);

    index = numeric_parameters.find_first_of(",");
    parameters_tail.emplace("dEC", std::stod(numeric_parameters.substr(0, index)));
    numeric_parameters = numeric_parameters.substr(index + 1);

    index = numeric_parameters.find_first_of(",");
    parameters_tail.emplace("theta_p", std::stod(numeric_parameters.substr(0, index)));
    numeric_parameters = numeric_parameters.substr(index + 1);

    std::cout << "folder " << folder << std::endl;
    std::cout << "type " << type << std::endl;
    std::cout << "parameters_tail" << std::endl;
    for (auto& i: parameters_tail) {
        std::cout << i.first << ", " << i.second << std::endl;
    }

    std::cout << index << ", " << numeric_parameters << std::endl;

    return 0;
}
*/
