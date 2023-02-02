#pragma once 


#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <sys/stat.h>
#include "parameters.hpp"

//input
std::string pop_space(std::string rawString);
bool import_kv_string(std::string variableNameStr, std::string variableValueStr, Param& param);
bool import_file(Param& param, std::string filepath);
std::vector<std::vector<double>> import_model_mesh(std::string filepath);