#pragma once
// structure definition

#include <math.h>
#include <cmath>
#include <vector>
#include <string>

////////////////////////////

struct Force{
    std::vector<double> ForceCurvature {0.0, 0.0, 0.0};
    std::vector<double> ForceArea {0.0, 0.0, 0.0};
    std::vector<double> ForceVolume {0.0, 0.0, 0.0};
    std::vector<double> ForceThickness {0.0, 0.0, 0.0};
    std::vector<double> ForceTilt {0.0, 0.0, 0.0};
    std::vector<double> ForceRegularization {0.0, 0.0, 0.0};
    std::vector<double> ForceSpline {0.0, 0.0, 0.0};
    std::vector<double> ForceTotal {0.0, 0.0, 0.0};
};

struct Energy{
    double      energyCurvature = 0.0;
    double      energyArea = 0.0;
    double      energyVolume = 0.0;
    double      energyThickness = 0.0;
    double      energyTilt = 0.0;
    double      energyRegularization = 0.0;
    double      energySpline = 0.0; //see spline point in input.params
    double      energyTotal = 0.0;
};

struct Vertex{
    int            Index;
    bool           IsOutLayer;
    bool           IsBoundary = false;
    bool           IsGhost = false;
    std::vector<double> Coord {0.0, 0.0, 0.0};
    std::vector<double> CoordPrevious {0.0, 0.0, 0.0};
    std::vector<double> ReferenceCoord {0.0, 0.0, 0.0};
    std::vector<int>    AdjacentFace; // faces that have this vertex. there should be 5 or 6 faces.
    std::vector<int>    AdjacentVertex; // vertice that are nearby the vertex_i. There should be 6 or less.
    Force          forcePrevious;
    Force          force;
};

struct Face{
    int         Index;
    bool        IsOutLayer;
    bool        IsBoundary = false;
    bool        IsGhost = false;
    bool        IsInsertionPatch = false;
    std::vector<int> AdjacentVertex {0, 0, 0}; //{index_vertex0, index_vertex1, index_vertex2}
    std::vector<int> OneRingVertex; // there should be 12 or 11 vertices.
    std::vector<int> AdjacentFace; // faces that are adjacent to this face. There should be 12 or 11 faces.
    double      SpontCurvature = 0.0;
    double      MeanCurvature = 0.0;
    std::vector<double> normVector {0.0, 0.0, 0.0}; // norm vector of this face element
    double      ElementArea = 0.0; // local area of this face.
    double      ElementVolume = 0.0; // local volume of this face.
    Energy      energyPrevious;
    Energy      energy;
    std::vector<std::vector<double>> forceAtGaussQuad;
};

struct Deformation{
    int DeformShape = 0;
    int DeformArea  = 0;
    int Undeform    = 0;
};


// sructure of shape funcitons.
struct Shapefunctions{
    std::vector<std::vector<double>> sf;
};

struct SubMatrix {
    std::vector<std::vector<double>> irregM; 
    std::vector<std::vector<double>> irregM1; 
    std::vector<std::vector<double>> irregM2; 
    std::vector<std::vector<double>> irregM3; 
    std::vector<std::vector<double>> irregM4;
};

struct Particle{
    int Index;
    int FaceIndex; // face index that this particle is located on
    std::vector<double> vwu {1.0/3.0, 1.0/3.0, 1.0/3.0};
    std::vector<double> Coord {0.0, 0.0, 0.0};
    double D; // diffusion constant; nm2/us
};
struct CrossFace{
    int FaceIndex = -1;
    std::vector<std::vector<int>> PairNodes {{-1,-1},
                                   {-1,-1}}; // there should be only two pairs of nodes, because only two sides of the face can pass through the diffusing-plane
    std::vector<std::vector<double>> vwu{ {0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0} }; // there should be only two joint points
    double CrossDistance;
};

struct Param {
    double kc = 0.0;
    double us = 0.0;
    double uv = 0.0;
    double k;
    double K;
    double S0; //target area
    double S; // area, total area 
    double V0; //target volume
    double V; // volume, total volume 
    double C0; // spontaneous curvature of insertions
    double c0; // spontaneous curvature of membrane
    double meanL;
    double gama_shape;
    double gama_area;
    double sigma = 0.0;
    int    GaussQuadratureN;
    int    subDivideTimes;
    bool   isInsertionAreaConstraint = false;
    bool   isAdditiveScheme = false;
    bool   isGlobalConstraint;
    double s0;
    bool   isBoundaryFixed;
    bool   isBoundaryPeriodic;
    bool   isBoundaryFree;
    double sideX;
    double sideY;
    double Radius;
    double l;     // for flat subdivision 
    bool   usingNCG = true;
    bool   isNCGstucked = false;
    bool   usingRpi = true;
    
    bool   xyzOutput = true;
    bool   meshpointOutput = true;
    bool   isEnergySplineIncluded = true;
    std::vector<std::vector<double>> splinePoints;
    std::vector<int> splinePoints_correspondingVertexIndex;
    double springConst = 13.5; //for spline points
    double splinePointsZcoordScaling = 0.04; //for z scaling
    double lbond = 12.0;
    std::string splinePointFileName = "";
    int    numIterations = 1E5;
    double timeStep = 0.1; //us
    double diffConst = 0.01; //um^2/us
    double KbT = 4.17; //1KbT = 4.17 pN.nm
    std::vector<std::vector<int>> InsertionPatch;

    Energy energy;
    Energy energyPrevious;
    Deformation DeformationList;
}; 
