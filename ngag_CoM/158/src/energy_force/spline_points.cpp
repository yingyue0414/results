#include "energy_force.hpp"

using namespace std;

/*
 * This method takes in the vector of spline point and calculate the
 * average coordinates. Based on the difference between spline points
 * and mesh vertices, a difference vector is calculated and compared to
 * the target bond length. (Supposed only in Z direction). Afterwards,
 * all the mesh points are moved in the direction of the target difference
 * vector.
 * 
 * Input: (n, 3) vector for spline points;
 *        (m) vector for mesh vertices;
 *        double target bond length
 * Output: true if successfully moved mesh vertices;
 */
bool moveVerticesBasedOnSpline(vector<vector<double>>& splinePoints,
                               vector<Vertex>& vertices,
                               double lbond){
    //FOR TESTING ONLY => fix the move distance
    bool fix_dir_down = true;
    if (fix_dir_down){
        vector<double> moveVec {0.0, 0.0, 50.0};
        for (int i = 0; i < vertices.size(); i++){
            for (int j = 0; j < 3; j++) {
                vertices[i].Coord[j] += moveVec[j];
            }
        }
        return true;
    }
    //calculate average spline points and average vertices
    vector<double> avg_splinePoints(3, 0.0);
    getAverageVector(splinePoints, avg_splinePoints);
    vector<double> avg_vertices(3, 0.0);
    for (Vertex vertex : vertices){
        avg_vertices += vertex.Coord;
    }
    const_multiplication(avg_vertices, 1.0/(vertices.size()) , avg_vertices);
    //get target difference vector
    vector<double> target_diffVec {0.0, 0.0, lbond};
    vector<double> moveVec = avg_splinePoints - avg_vertices + target_diffVec;
    cout<<"MoveVec="<<moveVec[0]<<","<<moveVec[1]<<","<<moveVec[2]<<endl;
    //move all mesh vertices by moveVec
    for (int i = 0; i < vertices.size(); i++){
        for (int j = 0; j < 3; j++) {
            vertices[i].Coord[j] += moveVec[j];
        }
    }
    return true;
}

/*
 * Used in moveVerticesBasedOnSpline
 * Serves to get average vector of a 2D vector
 * Note that avgVec is assumed to be initialized
 * with zeros and has the same length as the first dimension of inVec2D.
 * Return true if succeeded
 */
bool getAverageVector(vector<vector<double>>& inVec2D, vector<double>& avgVec){
    for (vector<double> sp : inVec2D){
        avgVec += sp;
    }
    const_multiplication(avgVec, 1.0/(inVec2D.size()), avgVec);
    return true;
}

/*
 * Get a vector of indexes of vertices that are closest to
 * the splinePoints vector provided.
 * 
 * Input: (n, 3) vector for spline points;
 *        (m) vector for mesh vertices;
 * Output: (n) vector for vertices indexes that closest to
 *        each point in spline points respectively
 */
vector<int> getClosestVertexIndex(vector<vector<double>>& splinePoints,
                                  vector<Vertex>& vertices) {
    //initialize
    vector<int> splinePoints_correspondingVertexIndex(splinePoints.size(),0);
    //loop over spline points and find closest vertices for each
    for (int i = 0; i < splinePoints.size(); i++){
        //initialize with the first vertex
        double minDistance2 = getSquaredDistance(splinePoints[i], vertices[0]);
        double newDistance2 = minDistance2;
        int minVertexIndex = 0;
        //loop over vertices to search for min distance
        for (int j = 1; j < vertices.size(); j++){
            newDistance2 = getSquaredDistance(splinePoints[i], vertices[j]);
            //overide min and index if new d < min
            if (newDistance2 < minDistance2) {
                minDistance2 = newDistance2;
                minVertexIndex = j;
            }
        }
        //set corresponding vertex index to minVertexIndex
        splinePoints_correspondingVertexIndex[i] = minVertexIndex;
        //std::cout<<i<<","<<minVertexIndex<<endl;
    }
    return splinePoints_correspondingVertexIndex;
}

/*
 * Used in getCloesestVertexIndex
 * Calculate the squared distance between two points denoted by (3) vector
 * and Vertex respectively.
 * Return the squared distance
 */
double getSquaredDistance(vector<double>& splinePoint, Vertex& vertex){
    //get difference vector
    vector<double> diffVec = splinePoint - vertex.Coord;
    return (diffVec[0] * diffVec[0] + diffVec[1] * diffVec[1] + diffVec[2] * diffVec[2]);
}

/*
 * calculate spline energy and force and add to Energy and Force
 * instances of vertex
 * return total spline energy if energy and force successfully overridden
 * 
 * Input: (n,3) spline points
 *        (m) vertices
 *        (n) corresponding index of vertices for spline points
 *        bool doLocalSearch: if true, search for neighbors of corresponding vertices;
 *              override the old vertex with new minimum (neighbors + the original vertex)
 *        double splinePointsZcoordScaling: scales the z-coordinate of spline points
 * Ouput: total spline energy if succeeded
 */
double calculateSplineEnergyForce(vector<vector<double>>& splinePoints,
                                vector<Vertex>& vertices,
                                vector<int> splinePoints_correspondingVertexIndex,
                                double lbond, double springConst, bool doLocalSearch,
                                double splinePointsZcoordScaling){
    int index = -1; // current vertex index
    double totalEnergy = 0.0;
    bool verbose = false; //verbose mode
    if (verbose){
        std::cout<<"=========spline energy============"<<endl;
    }

    vector<vector<double>> scaledSplinePoints = splinePoints;

    //iterate over spline points
    for (int i = 0; i < scaledSplinePoints.size(); i++){
        index = splinePoints_correspondingVertexIndex[i];
        double distance = sqrt(getSquaredDistance(scaledSplinePoints[i], vertices[index]));
        //calculate energy = 0.5 * k * (r - l)^2
        totalEnergy += 0.5 * springConst * (distance - lbond) * (distance - lbond);
        //vertices[index].energy.EnergySpline = 0.5 * springConst * (distance - lbond) * (distance - lbond);
        //vertices[index].energy.EnergyTotal += vertices[index].energy.EnergySpline;
        //calculate force = - k * (r - l) * normalize(vertexvec- splinepointvec)
        vector<double> r_unit(3);
        getUnitVector(vertices[index].Coord - scaledSplinePoints[i], r_unit);
        vertices[index].force.ForceSpline = - springConst * (distance - lbond) * r_unit;
        if (verbose) {
            std::cout<<index<<","<<vertices[index].force.ForceSpline[0]<<","
                                 <<vertices[index].force.ForceSpline[1]<<","
                                 <<vertices[index].force.ForceSpline[2]<<endl;
        }
        vertices[index].force.ForceTotal += vertices[index].force.ForceSpline;
    }

    if (verbose){
        std::cout<<"==================================="<<endl;
    }
    //search for closest pt if do local search == true

    return totalEnergy;
}