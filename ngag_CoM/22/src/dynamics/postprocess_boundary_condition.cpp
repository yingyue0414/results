#include "dynamics.hpp"

using namespace std;

//calculate real point relative to the given ghost point (index(real) - index(given)) in periodic boundary condition
//returns 0 if real point (not on 4th ring) given 
//an arbitrary real point is chosen if given 4th
//input: i = index of given pt; n = number of segments of columns; m = number of segments of rows
int get_relatve_pt_periodic(int i, int n, int m) {
    int indexRelative = 0;
    int colnum = i%(n+1);
    int rownum = i/(n+1);
    //get real point based on the position
    if (colnum < 4) {
        indexRelative += n-6;
    } else if (colnum > n-3) {
        indexRelative -= n-6;
    }
    if (rownum < 4) {
        indexRelative += (m-6)*(n+1);
    } else if (rownum > m-3) {
        indexRelative -= (m-6)*(n+1);
    }
    return indexRelative;
}

//post move step processing on triangular mesh - syncing ghost vertices and 4th ring real
//vertices with the corresponding mesh points
void  postprocess_ghost_periodic(Param& param, gsl_matrix *verticesOnMesh, vector<Vertex>& vertex) {
    double sidex = param.sideX; double sidey = param.sideY; double l = param.l;
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double dy = sqrt(3.0)/2.0 * dx; int m = round(sidey/dy);
    if  (pow(-1.0,m) < 0){
        m += 1;
    }
    
    double original = 0.0;
    //3.2 process real point
    for (int i = 0; i < vertex.size(); i++) {
        int indexRelative = 0;
        //calculate 4th point pair
        if (!(vertex[i].IsGhost)) {
            indexRelative = get_relatve_pt_periodic(i, n, m);
            if (indexRelative != 0) { //avoid redundant calculation
                //cout << "index relative:" << i << ", " << (i + indexRelative) << endl;

                double displacement = 0.0;
                for (int j = 0; j < 3; j++) {
                    //coord before displacement is stored in vertex[][]
                    //coord after displacement is stored in verticesOnMesh (gsl mat)
                    displacement = gsl_matrix_get(verticesOnMesh, i+indexRelative, j) - vertex[i+indexRelative].Coord[j];
                    original = vertex[i].Coord[j];
                    gsl_matrix_set(verticesOnMesh, i, j, original + displacement);
                }
            }
        }
    }

    //3.3 process ghost point
    for (int i = 0; i < vertex.size(); i++) {
        //calcualate real point for periodic boundary condition
        int indexRelative = 0;
        if (vertex[i].IsGhost) {
        //in order to sync up 4th ring pts
            indexRelative = get_relatve_pt_periodic(i, n, m);
            if (indexRelative != 0) { //avoid redundant calculation
                //cout << "index relative:" << i << ", " << (i + indexRelative) << endl;

                double displacement = 0.0;
                for (int j = 0; j < 3; j++) {
                    //coord before displacement is stored in vertex[][]
                    //coord after displacement is stored in verticesOnMesh (gsl mat)
                    displacement = gsl_matrix_get(verticesOnMesh, i+indexRelative, j) - vertex[i+indexRelative].Coord[j];
                    original = vertex[i].Coord[j];
                    gsl_matrix_set(verticesOnMesh, i, j, original + displacement);
                }
            }
            
        }
    }

    //sync vertex with verticeOnMesh
    for  (int i = 0; i < vertex.size(); i++) {
        //sync vertex with verticesOnMesh
        for (int j = 0; j < 3; j++) {
            vertex[i].Coord[j] = gsl_matrix_get(verticesOnMesh, i, j);
        }
    }
}

void postprocess_ghost_free (Param& param, gsl_matrix *verticesOnMesh,
                    vector<Vertex>& vertex, vector<Face>& face) {
    //output level for testing - >10 to enable; >100 to enable checkpoints
    int OUTPUTLEVEL = 0;
    
    //4.2.1 record all boundary vertex index; sync vertex with verticesOnMesh
    vector<bool> hasSet(vertex.size(), false);
    vector<int> boundPts;//record if the point position has been approximated
    for  (int i = 0; i < vertex.size(); i++) {
        //sync vertex with verticesOnMesh
        for (int j = 0; j < 3; j++) {
            vertex[i].Coord[j] = gsl_matrix_get(verticesOnMesh, i, j);
        }
        if (vertex[i].IsBoundary || vertex[i].IsGhost) {
            boundPts.push_back(i);
        } else {
            hasSet[i] = true;
        }
    }

    //adjacent vertex pair saves the ravelled list of adjacent vertices that form
    //a face with the given vertex; -1 denotes unused / undefined

    int itrwhile = 0; //test only, keep track of while iteration
    while (!(std::all_of(
            std::begin(hasSet), 
            std::end(hasSet), 
            [](bool i)
                    { 
                    return i;
                    } // if all of hasSet return true
        ))) {

        //for test iteration - PRESS TO CONTINUE
        //print_mat_contents(verticesOnMesh);
        if (OUTPUTLEVEL > 10) {
            cout << "Size remaining: " << boundPts.size() << endl; 

            cout << "Boundpts: ";
            for (auto j: boundPts)
                cout << j << " , ";
            cout << endl;
            cout << "hasSet: ";
            for (auto j: hasSet)
                cout << j << " , ";
            cout << endl;
        }
        if (OUTPUTLEVEL > 99) {
            cout << "Press Enter to Continue";
            cin.ignore(std::numeric_limits<streamsize>::max(),'\n'); 
        }

        //4.2.2 iterate through current boundary points and search for
        //neighbor triangle of each boundary point
        for (int bdIndex = 0; bdIndex < boundPts.size(); bdIndex++) {
            if (!hasSet[boundPts[bdIndex]]) {
                vector<int> * adjFace = &(vertex[boundPts[bdIndex]].AdjacentFace);
                //initialize parallelogram index
                int adjPt1 = -1;
                int adjPt2 = -1;
                int farPt = -1;
                for (int adjVPIndex = 0; adjVPIndex < adjFace->size(); adjVPIndex++) {
                    vector<int> AdjVertices = face[adjFace->at(adjVPIndex)].AdjacentVertex;
                    //erase-remove idiom to remove self from neighbor triangle vertices
                    AdjVertices.erase(std::remove(AdjVertices.begin(), AdjVertices.end(), 
                            boundPts[bdIndex]), AdjVertices.end());
                    if(hasSet[AdjVertices[0]] && hasSet[AdjVertices[1]] ) {

                        //4.2.3 get the point that makes a triangle with the given point
                        vector<int> * adjVP1 = &(vertex[AdjVertices[0]].AdjacentVertex);
                        std::sort(adjVP1->begin(), adjVP1->end());
                        vector<int> * adjVP2 = &(vertex[AdjVertices[1]].AdjacentVertex);
                        std::sort(adjVP2->begin(), adjVP2->end());
                        vector<int> adjVPCommon;
                        std::set_intersection(adjVP1->begin(), adjVP1->end(), adjVP2->begin(), adjVP2->end(), std::inserter(adjVPCommon, adjVPCommon.begin()));
                        for (int VPCommonInd = 0; VPCommonInd < adjVPCommon.size(); VPCommonInd++) {
                            if (adjVPCommon[VPCommonInd]!=boundPts[bdIndex] &&
                                        hasSet[adjVPCommon[VPCommonInd]]) {
                                adjPt1 = AdjVertices[0];
                                adjPt2 = AdjVertices[1];
                                farPt = adjVPCommon[VPCommonInd];
                            }
                        }

                        //output adjacent point vectors and their common elements
                        if (OUTPUTLEVEL > 10) {
                            cout << "adjVP1: ";
                            for (auto j: *adjVP1)
                                cout << j << " , ";
                            cout << endl;

                            cout << "adjVP2: ";
                            for (auto j: *adjVP2)
                                cout << j << " , ";
                            cout << endl;

                            cout << "adjVPCommon: ";
                            for (auto j: adjVPCommon)
                                cout << j << " , ";
                            cout << endl;
                        }

                        
                    }

                    //test only
                    if (OUTPUTLEVEL > 10) {
                        cout << "adjFace->at(adjVPIndex): " << adjFace->at(adjVPIndex) << endl;
                        cout << "Pt: " << boundPts[bdIndex] << endl;
                        cout << "AdjVertices(size): " << AdjVertices.size() << endl;
                        cout << "AdjVertices(ele): " << AdjVertices[0] <<" , " <<AdjVertices[1] << endl;
                        cout << "adjPt1: " << adjPt1 << endl;
                        cout << "adjPt2: " << adjPt2 << endl;
                        cout << "farPt: " << farPt << endl;
                    }
                    if (OUTPUTLEVEL > 99) {
                        cout << "Press Enter to Continue";
                        cin.ignore(std::numeric_limits<streamsize>::max(),'\n');   
                    }
                }
                //adjPt1, adjPt2, farPt set at this point; SKIP RESET if not
                //4.2.4 Parallelogram approximation for the boundary and ghost points
                //v(target) = v(a)+v(b)-v(far)
                if (farPt != -1) {
                    vertex[boundPts[bdIndex]].Coord = vertex[adjPt1].Coord + vertex[adjPt2].Coord - vertex[farPt].Coord;
                    hasSet[boundPts[bdIndex]] = true;
                    
                } 
            }
        }
        //output iteration number
        if (OUTPUTLEVEL > 10) {
            std::cout << itrwhile << endl;//test only
            itrwhile++;
        }   
    }
    //sync verticeOnMesh with vertex
    for  (int i = 0; i < vertex.size(); i++) {
        for (int j = 0; j < 3; j++) {
            gsl_matrix_set(verticesOnMesh, i, j, vertex[i].Coord[j]);
        }
    }
}