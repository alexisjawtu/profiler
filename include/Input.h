#ifndef _INPUT_H
#define _INPUT_H

#include <string>
#include <map>
#include <iostream>
#include <fstream>


class Input: public std::map<std::string, double> {
    public:
        std::string cylinder_folder;

        void print()
        {
            int i {0};
            std::cout << i << ". " << cylinder_folder << '\n';

            for (auto& p: *this)
            {
                i++;
                std::cout << i << ". " << p.first << " = " << p.second << '\n';
            }
        }

        Input()
        {
            std::map<std::string, double>();
        }

        Input(std::string inputfile)
        {
            std::map<std::string, double>();

            std::ifstream input(inputfile);
            std::string current_name;
            double current_value;

            /**
             * The first field is the name of the variable containing 
             * the cylinder output, so we skip it.
            **/

            input >> cylinder_folder;  
            input >> cylinder_folder;

            while (!input.eof())
            {
                input >> current_name;
                input >> current_value;
                (*this)[current_name] = current_value;
            }
        }

        ~Input() {}
};

#endif
