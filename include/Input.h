#ifndef _INPUT_H
#define _INPUT_H

#include <string>
#include <map>
#include <iostream>
#include <fstream>


class Input {
    public:
        std::map<std::string, double> prof_params;
        std::string cylinder_folder;
        double thickness;
        double levels;

        void print()
        {
            int i {0};
            std::cout << i << ". " << cylinder_folder << '\n';

            for (auto& p: this -> prof_params)
            {
                i++;
                std::cout << i << ". " << p.first << " = " << p.second << '\n';
            }
        }

        Input() {};

        Input(std::string inputfile)
        {
            std::ifstream input(inputfile);
            std::string current_name;
            double current_value;

            /**
             * The first field is the name of the variable,
             * so we skip it.
            **/
            input >> current_name;  
            input >> cylinder_folder;  

            input >> current_name;
            input >> thickness;

            input >> current_name;
            input >> levels;

            while (!input.eof())
            {
                input >> current_name;
                input >> current_value;
                (this -> prof_params)[current_name] = current_value;
            }
            input.close();
        }

        ~Input() {}
};

#endif
