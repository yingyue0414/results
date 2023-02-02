#include "io.hpp"

using namespace std;

/*
 * pop all spaces & tabs from the given string (\s and \t)
 * return popped string
 * the input string is unchanged
 */
std::string pop_space(std::string rawString) {
	std::string poppedString = rawString;
	while (poppedString.find(" ") != std::string::npos) {
		poppedString.erase(poppedString.find(" "), 1);
	}
    while (poppedString.find("\t") != std::string::npos) {
        poppedString.erase(poppedString.find("\t"), 1);
    }
	return poppedString;
}

/*
 * import parameters from key-value strings
 * store value in given Param object
 * pop space before match
 */
bool import_kv_string(std::string variableNameStr, std::string variableValueStr, Param& param) {
    if (variableNameStr.compare("numIterations") == 0) {
		param.numIterations = std::stoi(variableValueStr);
		std::cout << "NUMITERATIONS overridden: " << variableValueStr
				<< std::endl;
		return true;
	}
    if (variableNameStr.compare("meshpointOutput") == 0) {
		param.meshpointOutput = (variableValueStr.compare("true") == 0);
		std::cout << "MESHPOINTOUTPUT overridden: " << variableValueStr
				<< std::endl;
		return true;
	}
    if (variableNameStr.compare("xyzOutput") == 0) {
		param.xyzOutput = (variableValueStr.compare("true") == 0);
		std::cout << "XYZOUTPUT overridden: " << variableValueStr
				<< std::endl;
		return true;
	}
    if (variableNameStr.compare("sideX") == 0) {
		param.sideX = std::stod(variableValueStr);
		std::cout << "rSphere overridden: " << variableValueStr << std::endl;
		return true;
	}
    if (variableNameStr.compare("sideY") == 0) {
		param.sideY = std::stod(variableValueStr);
		std::cout << "rSphere overridden: " << variableValueStr << std::endl;
		return true;
	}
	if (variableNameStr.compare("lMeshSide") == 0) {
		param.l = std::stod(variableValueStr);
		std::cout << "LMESHSIDE overridden: " << variableValueStr << std::endl;
		return true;
	}
	if (variableNameStr.compare("c0Insertion") == 0) {
		param.C0 = std::stod(variableValueStr);
		std::cout << "C0INSERTION overridden: " << variableValueStr
				<< std::endl;
		return true;
	}
	if (variableNameStr.compare("c0Membrane") == 0) {
		param.c0 = std::stod(variableValueStr);
		std::cout << "C0MEMBRANE overridden: " << variableValueStr << std::endl;
		return true;
	}
	if (variableNameStr.compare("kcMembraneBending") == 0) {
		param.kc = std::stod(variableValueStr);
		std::cout << "KCMEMBRANEBENDING overridden: " << variableValueStr
				<< std::endl;
		return true;
	}
	if (variableNameStr.compare("usMembraneStretching") == 0) {
		param.us = std::stod(variableValueStr);
		std::cout << "USMEMBRANESTRETCHING overridden: " << variableValueStr
				<< std::endl;
		return true;
	}
	if (variableNameStr.compare("uvVolumeConstraint") == 0) {
		param.uv = std::stod(variableValueStr);
		std::cout << "UVVOLUMECONSTRAINT overridden: " << variableValueStr
				<< std::endl;
		return true;
	}
    if (variableNameStr.compare("timeStep") == 0) {
		param.timeStep = std::stod(variableValueStr);
		std::cout << "TIMESTEP overridden: " << variableValueStr
				<< std::endl;
		return true;
	}
    if (variableNameStr.compare("diffConst") == 0) {
		param.diffConst = std::stod(variableValueStr);
		std::cout << "DIFFCONST overridden: " << variableValueStr
				<< std::endl;
		return true;
	}
    if (variableNameStr.compare("KbT") == 0) {
		param.KbT = std::stod(variableValueStr);
		std::cout << "KBT overridden: " << variableValueStr
				<< std::endl;
		return true;
	}
	if (variableNameStr.compare("splinePointFileName") == 0) {
		param.splinePointFileName = variableValueStr;
		param.isEnergySplineIncluded = true;
		std::cout << "SPLINEPOINTFILENAME overriden: " << variableValueStr
				<< std::endl;
		return true;
	}
	std::cout << "VARIABLE NOT SUPPORTED: " << variableNameStr << std::endl;
	return false;
}

/*
 * load parameter file
 * input file path
 */
bool import_file(Param& param, std::string filepath) {
    //load in parameter file with the given filename
    std::ifstream ifile(filepath, std::ios::in);
    std::vector<std::string> parameters; //convert params file data to vector of rows

    //check to see that the file was opened correctly:
	if (!ifile.is_open()) {
		std::cout << "There was a problem opening the parameter file!\n" << endl;
		//exit(1); //exit or do additional error checking
	}

	std::string str = "";
	//keep reading line by line from the text file so long as data exists: 
    //(delete rows starting with # )
	while (getline(ifile, str)) {
		if (str[0]!='#') {
			parameters.push_back(str);
		}
	}

	//store kv pair
	std::string variableNameStr = "";
	std::string variableValueStr = "";

    //iterate over rows
    for (int i = 0; i < parameters.size(); ++i) {

		//use "=" and "#" as marks to find name-value pairs
		if (parameters[i].find("=") != std::string::npos) {

			variableNameStr = parameters[i].substr(0, parameters[i].find("=")); // raw variable name string! need to pop space!
			variableValueStr = parameters[i].substr(
					parameters[i].find("=") + 1); // need to delete comment
			//std::cout << "RAW VALUE: " << variableValueStr << std::endl; //for testing only
			//std::cout << variableValueStr.find("#") << std::endl;//for testing only
			if (variableValueStr.find("#") != std::string::npos) {
				variableValueStr = variableValueStr.substr(0,
						variableValueStr.find("#")); // raw value string! need to pop space!
			}
			//std::cout << "UNCOMMENTED VALUE: " << variableValueStr << std::endl; //for testing only

			// pop space
			//std::cout << variableNameStr << "::" << variableValueStr << std::endl;//for testing only
			variableNameStr = pop_space(variableNameStr);
			variableValueStr = pop_space(variableValueStr);
			//std::cout << variableNameStr << "::" << variableValueStr << std::endl;//for testing only

			// import kv string and load values to local variables
			import_kv_string(variableNameStr, variableValueStr, param);
		}
	}

	//end of import
	std::cout
			<< "============================END OF INPUT============================"
			<< std::endl;
    return true;
}

/*
 * load in model mesh file for adhesion of the triangular mesh
 * to the model mesh via adding an extra energy term -- E adhesion_geometry
 * and save the mesh infomation (vector<vector<double>> (n,3)) in model_mesh
 * 
 * input file path in the format of .csv (assuming no endline comma):
 * -- "x, y, z"
 * 
 */
vector<vector<double>> import_model_mesh(std::string filepath){
	//load in parameter file with the given filename
    std::ifstream ifile(filepath, std::ios::in);
    std::vector<std::string> meshdata; //convert params file data to vector of rows

	//check to see that the file was opened correctly:
	if (!ifile.is_open()) {
		std::cerr << "There was a problem opening the parameter file!\n";
		exit(1); //exit or do additional error checking
	}

	std::string str = "";
	//keep reading line by line from the text file so long as data exists: 
	//pop all spaces and tabs in the process
    //(delete rows starting with # )
	while (getline(ifile, str)) {
		if (str[0]!='#') {
			str = pop_space(str);
			meshdata.push_back(str);
		}
	}

	//iterate over rows of meshdata
	//initialize double vector
	vector<vector<double>> model_mesh (meshdata.size(), vector<double>(3));
	//strings of mesh data values along x, y, z-axis
	std::string x_str = "0";
	std::string y_str = "0";
	std::string z_str = "0";
	//index of comma
	int comma_index = 0;
	for (int i = 0; i < meshdata.size(); ++i){

		//search for comma (all spaces and tabs were popped) to find values
		//assuming no spaces and tabs in the data right now
		//stoi error if not enough values! (fewer than 3)
		//ignore the fourth value and all values afterwards
		
		//x
		comma_index = meshdata[i].find(",");
		x_str = meshdata[i].substr(0, comma_index);
		model_mesh[i][0] = std::stod(x_str);
		meshdata[i] = meshdata[i].substr(comma_index + 1);
		//y
		comma_index = meshdata[i].find(",");
		y_str = meshdata[i].substr(0, comma_index);
		model_mesh[i][1] = std::stod(y_str);
		//z
		z_str = meshdata[i].substr(comma_index + 1);
		model_mesh[i][2] = std::stod(y_str);
		
	}

	return model_mesh;
}