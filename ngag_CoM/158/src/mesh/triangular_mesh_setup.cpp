#include "mesh.hpp"
using namespace std;

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
vector<Vertex> setVertex_Loop_scheme(double sidex, double sidey, double l){ // vertex position
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double a = dx;
    double dy = sqrt(3.0)/2.0 * a; 
    int m = round(sidey/dy);      // y axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    double lx = n * dx;
    double ly = m * dy;

    int nodenum = (n+1)*(m+1);
    vector<Vertex> vertex(nodenum);

    #pragma omp parallel for
    for (int j = 0; j < m+1; j++){
        bool isEvenJ = false;
        if ( pow(-1.0,j) > 0.0 ){
            isEvenJ = true;
        }
        for (int i = 0; i < n+1; i++){
            int index = (n+1)*j + i;
            double x = i*dx;
            if ( isEvenJ == true ){
                x = x + a/2.0;
            }
            double y = j*dy;
            vertex[index].Index = index;
            vertex[index].Coord[0] = x - lx/2.0; 
            vertex[index].Coord[1] = y - ly/2.0; 
            vertex[index].Coord[2] = 0.0;
        }
    }

    return vertex;
}
vector<Face> setFace_Loop_scheme(double sidex, double sidey, double l){ // face and its surrounding vertex
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double a = dx;
    double dy = sqrt(3.0)/2.0 * a; 
    int m = round(sidey/dy);      // y axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    
    int facenum = m*n*2;
    vector<Face> face(facenum);
    #pragma omp parallel for
    for (int j = 0; j < m; j++){
        bool isEvenJ = false;
        if ( pow(-1.0,j) > 0.0 ){
            isEvenJ = true;
        }
        for (int i = 0; i < n; i++){
            int index = 2*n*j + i*2;
            int node1, node2, node3, node4;
            if ( isEvenJ == false ){
                node1 = (n+1)*j + i;
                node2 = (n+1)*(j+1) + i;
                node3 = (n+1)*j + (i+1);
                node4 = (n+1)*(j+1) + (i+1);
            }else{
                node1 = (n+1)*(j+1) + i;
                node2 = (n+1)*(j+1) + (i+1);
                node3 = (n+1)*j + i;
                node4 = (n+1)*j + (i+1);
            }
            face[index].Index = index;
            face[index].AdjacentVertex[0] = node1;
            face[index].AdjacentVertex[1] = node2;
            face[index].AdjacentVertex[2] = node3;
            index = index + 1;
            face[index].Index = index;
            face[index].AdjacentVertex[0] = node2;
            face[index].AdjacentVertex[1] = node4;
            face[index].AdjacentVertex[2] = node3;
        }
    }
    return face;
}

void determine_Boundary_vertex_face(double sidex, double sidey, double l, vector<Vertex>& vertex, vector<Face>& face){
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double a = dx;
    double dy = sqrt(3)/2 * a; 
    int m = round(sidey/dy);      // y axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }

    #pragma omp parallel for
    for (int j = 0; j < m+1; j++){
        for (int i = 0; i < n+1; i++){
            int index = (n+1)*j + i;
            if ( j == 0 || j == m || i == 0 || i == n ){
                vertex[index].IsBoundary = true; // if element is 1, then this vertex is on boundary
            }
        }
    }

    #pragma omp parallel for 
    for (int i = 0; i < face.size(); i++){
        int node1 = face[i].AdjacentVertex[0];
        int node2 = face[i].AdjacentVertex[1];
        int node3 = face[i].AdjacentVertex[2];
        if ( vertex[node1].IsBoundary == true || vertex[node2].IsBoundary == true || vertex[node3].IsBoundary == true ){
            face[i].IsBoundary = true; // this face is on boundary
        }
    }
}

void determine_Ghost_vertex_face(double sidex, double sidey, double l, bool isBoundaryFixed, bool isBoundaryPeriodic, bool isBoundaryFree, vector<Vertex>& vertex, vector<Face>& face){
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double a = dx;
    double dy = sqrt(3.0)/2.0 * a; 
    int m = round(sidey/dy);      // y axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    
    ///////////////////////////////////////// ghost vertex
    int vertexnum = (n+1)*(m+1);
    vector<int> TopBottom; 
    vector<int> LeftRight; 
    if (isBoundaryPeriodic == true && isBoundaryFree == false && isBoundaryFixed == false){
        TopBottom.insert(TopBottom.end(),{0, 1, 2, m-2, m-1, m});
        LeftRight.insert(LeftRight.end(),{0, 1, 2, n-2, n-1, n});
    }else if (isBoundaryFree == true && isBoundaryPeriodic == false && isBoundaryFixed == false){
        TopBottom.insert(TopBottom.end(),{0, m});
        LeftRight.insert(LeftRight.end(),{0, n});
    }
    // top and bottom ghost vertex
    for (int k = 0; k < TopBottom.size(); k++){
        int j = TopBottom[k];
        #pragma omp parallel for
        for (int i = 0; i < n+1; i++){
            int index = (n+1)*j + i;
            vertex[index].IsGhost = true; 
        }
    }
    // left and right ghost vertex
    for (int k = 0; k < LeftRight.size(); k++){
        int i = LeftRight[k];
        #pragma omp parallel for
        for (int j = 0; j < m+1; j++){
            int index = (n+1)*j + i;
            vertex[index].IsGhost = true; 
        }
    }
    //////////////////////////////////////// ghost face
    int facenum = m*n*2;
    TopBottom.clear();; 
    LeftRight.clear();; 
    if (isBoundaryPeriodic == true && isBoundaryFree == false && isBoundaryFixed == false){
        TopBottom.insert(TopBottom.end(),{0, 1, 2, m-3, m-2, m-1});
        LeftRight.insert(LeftRight.end(),{0, 1, 2, n-3, n-2, n-1});
    }else if (isBoundaryFree == true && isBoundaryPeriodic == false && isBoundaryFixed == false){
        TopBottom.insert(TopBottom.end(),{0, m-1});
        LeftRight.insert(LeftRight.end(),{0, n-1});
    }
    for (int k = 0; k < TopBottom.size(); k++){
        int j = TopBottom[k];
        #pragma omp parallel for
        for (int i = 0; i < n; i++){
            int index = 2*n*j + i*2;
            face[index].IsGhost = true; 
            face[index+1].IsGhost = true;
        }
    }
    for (int k = 0; k < LeftRight.size(); k++){
        int i = LeftRight[k];
        #pragma omp parallel for
        for (int j = 0; j < m; j++){
            int index = 2*n*j + i*2;
            face[index].IsGhost = true; 
            face[index+1].IsGhost = true;
        }
    }
}
// find the faces around the vertex, probaly 5 or 6 faces that has vertex_i
void determine_AdjacentFace_for_vertex(vector<Vertex>& vertex, vector<Face>& face){
    #pragma omp parallel for 
    for (int i = 0; i < vertex.size(); i++){
        vector<int> AdjacentFacetmp;
        for (int j = 0; j < face.size(); j++){
            for (int k = 0; k < 3; k++){ // AdjacentVertex.size() = 3
                if ( i == face[j].AdjacentVertex[k] ){
                    AdjacentFacetmp.push_back(j);
                }
            }
        }  
        vertex[i].AdjacentFace = AdjacentFacetmp;
    }
}

// find out the nearby vertices around the vertex_i. there could be 6 or less.
void determine_AdjacentVertex_for_vertex(vector<Vertex>& vertex, vector<Face>& face){  
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++){
        vector<int> AdjacentVertextmp;
        for (int j = 0; j < vertex[i].AdjacentFace.size(); j++){
            int FaceIndex = vertex[i].AdjacentFace[j]; 
            for (int k = 0; k < face[FaceIndex].AdjacentVertex.size(); k++){
                int VertexIndex = face[FaceIndex].AdjacentVertex[k];
                if ( VertexIndex != i ){
                    bool IsListed = false;
                    for (int m = 0; m < AdjacentVertextmp.size(); m++){
                        if ( VertexIndex == AdjacentVertextmp[m] ){
                            IsListed = true;
                        }
                    }
                    if ( IsListed == false ){
                        AdjacentVertextmp.push_back(VertexIndex);
                    }
                }
            }
        }
        vertex[i].AdjacentVertex = AdjacentVertextmp;
    }
}

// find the common vertex that node1 and node2 share, but not node3.
int Find_NodeIndex(int node1, int node2, int node3, vector<Vertex>& vertex){
    int node = -1;
    for (int i = 0; i < vertex[node1].AdjacentVertex.size(); i++){
        int nodetmp1 = vertex[node1].AdjacentVertex[i];
        for (int j = 0; j < vertex[node2].AdjacentVertex.size(); j++){
            int nodetmp2 = vertex[node2].AdjacentVertex[j];
            if ( nodetmp1 == nodetmp2 && nodetmp1 != node3 ){
                node = nodetmp1;
            }
        }
    }
    if ( node == -1 ){
        cout<<"Wrong! No efficent NodeIndex is found in Find_NodeIndex. Node1 = "<<node1<<", Node2 = "<<node2<<", Node3 = "<<node3<<endl;
        exit(0);
    }
    return node;
}

// To find out the one-ring vertices aound face_i. It should be 12 for the flat surface because we set it up only with regular patch.
// The boundary faces do not have complete one-ring, neither it will be called in the code, so no need to store their one-ring-vertex
void determine_OneRingVertex_for_face(vector<Vertex>& vertex, vector<Face>& face){ 
    // two types of patch: 1. regular patch with 12 one-ring vertices, each vertex has 6 closest nodes
    //                     2. irregular patch with 11 one-ring vertices, one vertex has 5 closest-nodes
    #pragma omp parallel for
    for (int i = 0; i < face.size(); i++){
        if ( face[i].IsBoundary == true ){ 
            continue;
        }
        face[i].OneRingVertex.clear();

        //int d1, d2, d3, d5, d6, d9, d10, d11, d12;
        int node0 = face[i].AdjacentVertex[0];  
        int node1 = face[i].AdjacentVertex[1];
        int node2 = face[i].AdjacentVertex[2];
        // regular patch, all three nodes have 6 neighbor faces or vertices.
        if ( vertex[node0].AdjacentFace.size() == 6 && vertex[node1].AdjacentFace.size() && vertex[node2].AdjacentFace.size() == 6 ){
            // note the order of the vertex
            int d4 = node0; int d7 = node1; int d8 = node2;
            // make sure d4, d7, d8 are in anti-clock-wise order
            vector<double> node4 = vertex[d4].Coord; vector<double> node7 = vertex[d7].Coord; vector<double> node8 = vertex[d8].Coord;
            vector<double> center = 1.0/3.0 * (node4 + node7 + node8);
            if ( dot(center, cross(node7-node4,node8-node4)) < 0 ){ // switch d7 and d8 position
                d7 = face[i].AdjacentVertex[2];
                d8 = face[i].AdjacentVertex[1];
                face[i].AdjacentVertex[1] = d7;
                face[i].AdjacentVertex[2] = d8;
            }
            int d3  = Find_NodeIndex(d4, d7, d8, vertex);
            int d11 = Find_NodeIndex(d7, d8, d4, vertex);
            int d5  = Find_NodeIndex(d4, d8, d7, vertex);
            int d1  = Find_NodeIndex(d3, d4, d7, vertex);
            int d2  = Find_NodeIndex(d4, d5, d8, vertex);
            int d6  = Find_NodeIndex(d3, d7, d4, vertex);
            int d9  = Find_NodeIndex(d8, d5, d4, vertex);
            int d10 = Find_NodeIndex(d7, d11, d8, vertex);
            int d12 = Find_NodeIndex(d8, d11, d7, vertex);
            vector<int> v { d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12 };
            face[i].OneRingVertex = v;
        // irregular patch, one node has 5 neighbors and the other two nodes have 6 neighbors.
        }else if ( vertex[node0].AdjacentFace.size() == 5 || vertex[node1].AdjacentFace.size() == 5 || vertex[node2].AdjacentFace.size() == 5  ){
            int d4, d7, d8;
            // make sure d4 is the one has 5 neighbors, and d4-d7-d8 are in anti-clock-wise order
            if ( vertex[node0].AdjacentFace.size() == 5 ){
                d4 = node0; d7 = node1; d8 = node2;
            }else if ( vertex[node1].AdjacentFace.size() == 5 ) {
                d4 = node1; d7 = node2; d8 = node0;
            }else if ( vertex[node2].AdjacentFace.size() == 5 ) {
                d4 = node2; d7 = node0; d8 = node1;
            }
            vector<double> node4 = vertex[d4].Coord; vector<double> node7 = vertex[d7].Coord; vector<double> node8 = vertex[d8].Coord;
            vector<double> center = 1.0/3.0 * (node4 + node7 + node8);
            if ( dot(center, cross(node7-node4,node8-node4)) < 0 ){ // switch d7 and d8 position
                int nodetmp = d7;
                d7 = d8;
                d8 = nodetmp;
                face[i].AdjacentVertex[0] = d4;
                face[i].AdjacentVertex[1] = d7;
                face[i].AdjacentVertex[2] = d8;
            }
            int d3  = Find_NodeIndex(d4, d7, d8, vertex);
            int d11 = Find_NodeIndex(d7, d8, d4, vertex);
            int d5  = Find_NodeIndex(d4, d8, d7, vertex);
            int d1  = Find_NodeIndex(d3, d4, d7, vertex);
            int d2  = Find_NodeIndex(d4, d5, d8, vertex);
            int d6  = Find_NodeIndex(d3, d7, d4, vertex);
            int d9  = Find_NodeIndex(d8, d5, d4, vertex);
            int d10 = Find_NodeIndex(d7, d11, d8, vertex);
            int d12 = Find_NodeIndex(d8, d11, d7, vertex);
            vector<int> v { d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12 };
            face[i].OneRingVertex = v;
        }
    }
}

void read_struture_vertex(vector<Vertex>& vertex, char* filename) {
    ifstream fin(filename);
    string line;
    int i = 0;
    while ( getline(fin,line) ) {
        istringstream sin(line);
        vector<string> positions;
        string info;
        while (getline(sin, info, ',')) {
            positions.push_back(info);
        }
        string xstr = positions[0];
        string ystr = positions[1];
        string zstr = positions[2];
        double x, y, z;
        stringstream sx, sy, sz;
        sx << xstr;
        sy << ystr;
        sz << zstr;
        sx >> x;
        sy >> y;
        sz >> z;
        vector<double> coord {x, y, z};
        vertex[i].Coord = coord;
        i++;
    }
    if ( i != vertex.size() ) {
        cout<< "Wrong! vertices number is "<<vertex.size()<<" but read number i = "<< i << endl;
    }
}

vector<vector<double>> setVMU(int GaussQuadratureN){ // To setup the Gaussian points for integral calculation. 'n' here is the Gausssian order.
    int n = GaussQuadratureN;
    vector<vector<double>> vmu;
    if (n==1){
        vector<vector<double> > vmutmp { {1.0/3.0, 1.0/3.0, 1.0/3.0} };
        vmu = vmutmp;
    }else if(n==2){
        vector<vector<double> > vmutmp { {1.0/6.0, 1.0/6.0, 4.0/6.0},
                                         {1.0/6.0, 4.0/6.0, 1.0/6.0},
                                         {4.0/6.0, 1.0/6.0, 1.0/6.0} };
        vmu = vmutmp;
    }else if(n==3){ // third order Guasssian may not converge. Weird!
        vector<vector<double> > vmutmp { {1.0/3.0, 1.0/3.0, 1.0/3.0},
                                         {1.0/5.0, 1.0/5.0, 3.0/5.0},
                                         {1.0/5.0, 3.0/5.0, 1.0/5.0},
                                         {3.0/5.0, 1.0/5.0, 1.0/5.0} };
        /*
        vector<vector<double> > vmutmp { {1.0/3.0, 1.0/3.0, 1.0/3.0},
                                         {2.0/15.0, 11.0/15.0, 2.0/15.0},
                                         {2.0/15.0, 2.0/15.0, 11.0/15.0},
                                         {11.0/15.0, 2.0/15.0, 2.0/15.0} };
        */
        vmu = vmutmp;
    }else if(n==4){
        vector<vector<double> > vmutmp { {0.44594849091597, 0.44594849091597, 0.10810301816807},
                                         {0.44594849091597, 0.10810301816807, 0.44594849091597},
                                         {0.10810301816807, 0.44594849091597, 0.44594849091597},
                                         {0.09157621350977, 0.09157621350977, 0.81684757298046},
                                         {0.09157621350977, 0.81684757298046, 0.09157621350977},
                                         {0.81684757298046, 0.09157621350977, 0.09157621350977} };
        vmu = vmutmp;
    }else if(n==5){
        vector<vector<double> > vmutmp { {0.33333333333333, 0.33333333333333, 0.33333333333333},
                                         {0.47014206410511, 0.47014206410511, 0.05971587178977},
                                         {0.47014206410511, 0.05971587178977, 0.47014206410511},
                                         {0.05971587178977, 0.47014206410511, 0.47014206410511},
                                         {0.10128650732346, 0.10128650732346, 0.79742698535309},
                                         {0.10128650732346, 0.79742698535309, 0.10128650732346},
                                         {0.79742698535309, 0.10128650732346, 0.10128650732346} };
        vmu = vmutmp;
    }else if(n==6){
        vector<vector<double> > vmutmp { {0.24928674517091, 0.24928674517091, 0.50142650965818},
                                         {0.24928674517091, 0.50142650965818, 0.24928674517091},
                                         {0.50142650965818, 0.24928674517091, 0.24928674517091},
                                         {0.06308901449150, 0.06308901449150, 0.87382197101700},
                                         {0.06308901449150, 0.87382197101700, 0.06308901449150},
                                         {0.87382197101700, 0.06308901449150, 0.06308901449150},
                                         {0.31035245103378, 0.63650249912140, 0.05314504984482},
                                         {0.63650249912140, 0.05314504984482, 0.31035245103378},
                                         {0.05314504984482, 0.31035245103378, 0.63650249912140},
                                         {0.63650249912140, 0.31035245103378, 0.05314504984482},
                                         {0.31035245103378, 0.05314504984482, 0.63650249912140},
                                         {0.05314504984482, 0.63650249912140, 0.31035245103378} };
        vmu = vmutmp;
    }

    return vmu;
}

vector<double> setVMUcoefficient(int GaussQuadratureN){ // coefficeints for Gaussian points. 'n' here is the Gausssian order.
    int n = GaussQuadratureN;
    vector<double> vmucoeff;
    if (n==1){
        vector<double> coefftmp { 1.0 };
        vmucoeff = coefftmp;
    }else if(n==2){
        vector<double> coefftmp { 1.0/3.0, 1.0/3.0, 1.0/3.0 };
        vmucoeff = coefftmp;
    }else if(n==3){
        vector<double> coefftmp { -0.56250000000000, 0.52083333333333, 0.52083333333333, 0.52083333333333 };
        vmucoeff = coefftmp;
    }else if(n==4){
        vector<double> coefftmp { 0.22338158967801, 0.22338158967801, 0.22338158967801, 0.10995174365532, 0.10995174365532, 0.10995174365532 };
        vmucoeff = coefftmp;
    }else if(n==5){
        vector<double> coefftmp { 0.22500000000000, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483 };
        vmucoeff = coefftmp;
    }else if(n==6){
        vector<double> coefftmp { 0.11678627572638, 0.11678627572638, 0.11678627572638, 0.05084490637021, 0.05084490637021, 0.05084490637021, 0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837 };
        vmucoeff = coefftmp;
    }

    return vmucoeff;
}

// shape functions on barycentric coordinates.
// 12 shape functions and their differential equations; 
/*
    shape_functions(0,:), shape functions;  
    shape_functions(1,:), differential to v; 
    shape_functions(2,:), differential to w; 
    shape_functions(3,:), double differential to v; 
    shape_functions(4,:), double differential to w;
    shape_functions(5,:), differential to v and w; 
    shape_functions(6,:), differential to w and v;
    // NOTE: The below function used the column as row, so I transpose it at the last line.  
*/
vector<vector<double>> determine_ShapeFunctions(vector<double>& vwu){
    double v = vwu[0]; double w = vwu[1]; double u = vwu[2];
    vector<vector<double>> sf(12,vector<double>(7)); // note, I transpose this matrix when output
    sf[0][0] = 1.0/12.0*(pow(u,4.0) + 2.0*pow(u,3.0)*v);
    sf[0][1] = 1.0/12.0*(-2.0*pow(u,3.0) - 6.0*pow(u,2.0)*v); 
    sf[0][2] = 1.0/12.0*(-4.0*pow(u,3.0) - 6.0*pow(u,2.0)*v); 
    sf[0][3] = u*v; 
    sf[0][4] = pow(u,2.0) + u*v; 
    sf[0][5] = 1.0/2.0*(pow(u,2.0) + 2.0*u*v); 
    sf[0][6] = 1.0/2.0*(pow(u,2.0) + 2.0*u*v); 
    sf[1][0] = 1.0/12.0*(pow(u,4.0) + 2.0*pow(u,3.0)*w); 
    sf[1][1] = 1.0/12.0*(-4.0*pow(u,3.0) - 6.0*pow(u,2.0)*w); 
    sf[1][2] = 1.0/12.0*(-2.0*pow(u,3.0) - 6.0*pow(u,2.0)*w);
    sf[1][3] = pow(u,2.0) + u*w; 
    sf[1][4] = u*w;
    sf[1][5] = 1.0/2.0*(pow(u,2.0) + 2.0*u*w); 
    sf[1][6] = 1.0/2.0*(pow(u,2.0) + 2.0*u*w); 
    sf[2][0] = 1.0/12.0*(pow(u,4.0) + 2.0*pow(u,3.0)*w + 6.0*pow(u,3.0)*v + 6.0*pow(u,2.0)*v*w + 12.0*pow(u,2.0)*pow(v,2.0) + 6.0*u*pow(v,2.0)*w + 6.0*u*pow(v,3.0) + 2.0*pow(v,3.0)*w + pow(v,4.0));
    sf[2][1] = 1.0/12.0*(2.0*pow(u,3.0) + 6.0*pow(u,2.0)*v - 6.0*u*pow(v,2.0) - 2.0*pow(v,3.0));
    sf[2][2] = 1.0/12.0*(-2.0*pow(u,3.0) - 6.0*pow(u,2.0)*w - 12.0*pow(u,2)*v - 12.0*u*v*w - 18.0*u*pow(v,2.0) - 6.0*pow(v,2.0)*w - 4.0*pow(v,3.0));
    sf[2][3] = -2.0*u*v;
    sf[2][4] = u*w + v*w + u*v + pow(v,2.0);
    sf[2][5] = 1.0/2.0*(-pow(u,2.0) - 2.0*u*v + pow(v,2.0));
    sf[2][6] = 1.0/2.0*(-pow(u,2.0) - 2.0*u*v + pow(v,2.0));
    sf[3][0] = 1.0/12.0*(6.0*pow(u,4.0) + 24.0*pow(u,3.0)*w + 24.0*pow(u,2.0)*pow(w,2.0) + 8.0*u*pow(w,3.0) + pow(w,4.0) + 24.0*pow(u,3.0)*v + 60.0*pow(u,2.0)*v*w + 36.0*u*v*pow(w,2.0) + 6.0*v*pow(w,3.0) + 24.0*pow(u,2.0)*pow(v,2.0) + 36.0*u*pow(v,2.0)*w + 12.0*pow(v,2.0)*pow(w,2.0) + 8.0*u*pow(v,3.0) + 6.0*pow(v,3.0)*w + pow(v,4.0));
    sf[3][1] = 1.0/12.0*(-12.0*pow(u,2.0)*w - 12.0*u*pow(w,2.0) - 2.0*pow(w,3.0) - 24.0*pow(u,2.0)*v - 48.0*u*v*w - 12.0*v*pow(w,2.0) -24.0*u*pow(v,2.0) - 18.0*pow(v,2.0)*w - 4.0*pow(v,3.0));
    sf[3][2] = 1.0/12.0*(-24.0*pow(u,2.0)*w - 24.0*u*pow(w,2.0) - 4.0*pow(w,3.0) - 12.0*pow(u,2.0)*v - 48.0*u*v*w - 18.0*v*pow(w,2.0) - 12.0*u*pow(v,2.0) - 12.0*pow(v,2.0)*w - 2.0*pow(v,3.0));
    sf[3][3] = -2.0*u*w - 2.0*pow(u,2.0) + v*w + pow(v,2.0);
    sf[3][4] = -2.0*pow(u,2.0) + pow(w,2.0) - 2.0*u*v + v*w;
    sf[3][5] = 1.0/2.0*(-2.0*pow(u,2.0) + pow(w,2.0) + 4.0*v*w + pow(v,2.0));
    sf[3][6] = 1.0/2.0*(-2.0*pow(u,2.0) + pow(w,2.0) + 4.0*v*w + pow(v,2.0));
    sf[4][0] = 1.0/12.0*(pow(u,4.0) + 6.0*pow(u,3.0)*w + 12.0*pow(u,2.0)*pow(w,2.0) + 6.0*u*pow(w,3.0) + pow(w,4.0) + 2.0*pow(u,3.0)*v + 6.0*pow(u,2.0)*v*w + 6.0*u*v*pow(w,2.0) + 2.0*v*pow(w,3.0));
    sf[4][1] = 1.0/12.0*(-2.0*pow(u,3.0) - 12.0*pow(u,2.0)*w - 18.0*u*pow(w,2.0) - 4.0*pow(w,3.0) - 6.0*pow(u,2.0)*v - 12.0*u*v*w - 6.0*v*pow(w,2.0));
    sf[4][2] = 1.0/12.0*(2.0*pow(u,3.0) + 6.0*pow(u,2.0)*w - 6.0*u*pow(w,2.0) - 2.0*pow(w,3.0));
    sf[4][3] = u*w + pow(w,2.0) + u*v + v*w;
    sf[4][4] = -2.0*u*w;
    sf[4][5] = 1.0/2.0*(-pow(u,2.0) - 2.0*u*w + pow(w,2.0));
    sf[4][6] = 1.0/2.0*(-pow(u,2.0) - 2.0*u*w + pow(w,2.0));
    sf[5][0] = 1.0/12.0*(2.0*u*pow(v,3.0) + pow(v,4.0)); 
    sf[5][1] = 1.0/12.0*(6.0*u*pow(v,2.0) + 2.0*pow(v,3.0)); 
    sf[5][2] = -1.0/6.0*pow(v,3.0);
    sf[5][3] = u*v; 
    sf[5][4] = 0.0;
    sf[5][5] = -1.0/2.0*pow(v,2.0); 
    sf[5][6] = -1.0/2.0*pow(v,2.0);
    sf[6][0] = 1.0/12.0*(pow(u,4.0) + 6.0*pow(u,3.0)*w + 12.0*pow(u,2.0)*pow(w,2.0) + 6.0*u*pow(w,3.0)+ pow(w,4.0) + 8.0*pow(u,3.0)*v + 36.0*pow(u,2.0)*v*w + 36.0*u*v*pow(w,2.0) + 8.0*v*pow(w,3.0) + 24.0*pow(u,2.0)*pow(v,2.0) + 60.0*u*pow(v,2.0)*w + 24.0*pow(v,2.0)*pow(w,2.0) + 24.0*u*pow(v,3.0) + 24.0*pow(v,3.0)*w + 6.0*pow(v,4.0));
    sf[6][1] = 1.0/12.0*(4.0*pow(u,3.0) + 18.0*pow(u,2.0)*w + 12.0*u*pow(w,2.0) + 2.0*pow(w,3.0) + 24.0*pow(u,2.0)*v + 48.0*u*v*w + 12.0*v*pow(w,2.0) + 24.0*u*pow(v,2.0) + 12.0*pow(v,2.0)*w);
    sf[6][2] = 1.0/12.0*(2.0*pow(u,3.0) + 6.0*pow(u,2.0)*w - 6.0*u*pow(w,2.0) - 2.0*pow(w,3.0) + 12.0*pow(u,2.0)*v - 12.0*v*pow(w,2.0) + 12.0*u*pow(v,2.0) - 12.0*pow(v,2.0)*w);
    sf[6][3] = pow(u,2.0) + u*w - 2.0*v*w - 2.0*pow(v,2.0);
    sf[6][4] = -2.0*u*w - 2.0*u*v - 2.0*v*w - 2.0*pow(v,2.0);
    sf[6][5] = 1.0/2.0*(pow(u,2.0) - 2.0*u*w - pow(w,2.0) - 4.0*v*w - 2.0*pow(v,2.0));
    sf[6][6] = 1.0/2.0*(pow(u,2.0) - 2.0*u*w - pow(w,2.0) - 4.0*v*w - 2.0*pow(v,2.0));
    sf[7][0] = 1.0/12.0*(pow(u,4.0) + 8.0*pow(u,3.0)*w + 24.0*pow(u,2.0)*pow(w,2.0) + 24.0*u*pow(w,3.0) + 6.0*pow(w,4.0) + 6.0*pow(u,3.0)*v + 36.0*pow(u,2.0)*v*w + 60.0*u*v*pow(w,2.0) + 24.0*v*pow(w,3.0) + 12.0*pow(u,2.0)*pow(v,2.0) + 36.0*u*pow(v,2.0)*w + 24.0*pow(v,2.0)*pow(w,2.0) + 6.0*u*pow(v,3.0) + 8.0*pow(v,3.0)*w + pow(v,4.0));
    sf[7][1] = 1.0/12.0*(2.0*pow(u,3.0) + 12.0*pow(u,2.0)*w + 12.0*u*pow(w,2.0) + 6.0*pow(u,2.0)*v - 12.0*v*pow(w,2.0) - 6.0*u*pow(v,2.0) - 12.0*pow(v,2.0)*w - 2.0*pow(v,3.0));
    sf[7][2] = 1.0/12.0*(4.0*pow(u,3.0) + 24.0*pow(u,2.0)*w + 24.0*u*pow(w,2.0) + 18.0*pow(u,2.0)*v + 48.0*u*v*w + 12.0*v*pow(w,2.0) + 12.0*u*pow(v,2.0) + 12.0*pow(v,2.0)*w + 2.0*pow(v,3.0));
    sf[7][3] = -2.0*u*w - 2.0*pow(w,2.0) - 2.0*u*v - 2.0*v*w;
    sf[7][4] = pow(u,2.0) - 2.0*pow(w,2.0) + u*v - 2.0*v*w;
    sf[7][5] = 1.0/2.0*(pow(u,2.0) - 2.0*pow(w,2.0) - 2.0*u*v - 4.0*v*w - pow(v,2.0));
    sf[7][6] = 1.0/2.0*(pow(u,2.0) - 2.0*pow(w,2.0) - 2.0*u*v - 4.0*v*w - pow(v,2.0));
    sf[8][0] = 1.0/12.0*(2.0*u*pow(w,3.0) + pow(w,4.0)); 
    sf[8][1] = -1.0/6.0*pow(w,3.0); 
    sf[8][2] = 1.0/12.0*(6.0*u*pow(w,2.0) + 2.0*pow(w,3.0));
    sf[8][3] = 0.0; 
    sf[8][4] = u*w;
    sf[8][5] = -1.0/2.0*pow(w,2.0); 
    sf[8][6] = -1.0/2.0*pow(w,2.0);
    sf[9][0] = 1.0/12.0*(2.0*pow(v,3.0)*w + pow(v,4.0));
    sf[9][1] = 1.0/12.0*(6.0*pow(v,2.0)*w + 4.0*pow(v,3.0)); 
    sf[9][2] = 1.0/6.0*pow(v,3.0);
    sf[9][3] = v*w + pow(v,2.0); 
    sf[9][4] = 0.0;
    sf[9][5] = 1.0/2.0*pow(v,2.0); 
    sf[9][6] = 1.0/2.0*pow(v,2.0);
    sf[10][0] = 1.0/12.0*(2.0*u*pow(w,3.0) + pow(w,4.0) + 6.0*u*v*pow(w,2.0) + 6.0*v*pow(w,3.0) + 6.0*u*pow(v,2.0)*w + 12.0*pow(v,2.0)*pow(w,2.0) + 2.0*u*pow(v,3.0) + 6.0*pow(v,3.0)*w + pow(v,4.0));
    sf[10][1] = 1.0/12.0*(4.0*pow(w,3.0) + 18.0*v*pow(w,2.0) + 6.0*u*pow(w,2.0) + 12.0*pow(v,2.0)*w + 12.0*u*v*w + 2.0*pow(v,3.0) + 6.0*u*pow(v,2.0));
    sf[10][2] = 1.0/12.0*(2.0*pow(w,3.0) + 6.0*u*pow(w,2.0) + 12.0*v*pow(w,2.0) + 12.0*u*v*w + 18.0*pow(v,2.0)*w + 6.0*u*pow(v,2.0) + 4.0*pow(v,3.0));
    sf[10][3] = pow(w,2.0) + v*w + u*w + u*v;
    sf[10][4] = u*w + v*w + u*v + pow(v,2.0);
    sf[10][5] = 1.0/2.0*(pow(w,2.0) + 4.0*v*w + 2.0*u*w + pow(v,2.0) + 2.0*u*v);
    sf[10][6] = 1.0/2.0*(pow(w,2.0) + 4.0*v*w + 2.0*u*w + pow(v,2.0) + 2.0*u*v);
    sf[11][0] = 1.0/12.0*(pow(w,4.0) + 2.0*v*pow(w,3.0)); 
    sf[11][1] = 1.0/6.0*pow(w,3.0); 
    sf[11][2] = 1.0/12.0*(4.0*pow(w,3.0) + 6.0*v*pow(w,2.0));
    sf[11][3] = 0.0; 
    sf[11][4] = pow(w,2.0) + v*w;
    sf[11][5] = 1.0/2.0*pow(w,2.0); 
    sf[11][6] = 1.0/2.0*pow(w,2.0);

    return transpose(sf);
}

// for irregular patch, more subdivision is needed. 
// for different sub-element, different new nodes are selected, select-matrix (SM)
// vertex here is the original 11 vertice, so vertex is 11*3 matrix
// M(17,11); mat M1(12,17); mat M2(12,17); mat M3(12,17); mat M4(11,17);
void determine_SubdivisionMatrix(vector<vector<double>>& M, vector<vector<double>>& SM1, vector<vector<double>>& SM2, vector<vector<double>>& SM3, vector<vector<double>>& SM4){
    int N = 6; double w = 3.0/8.0/N; // w=1/N*(5/8-(3/8+1/4*cos(2*pi/N))^2);
    int N1 = 5; double w1 = 3.0/8.0/N1; // w1=1/N1*(5/8-(3/8+1/4*cos(2*pi/N1))^2);
    double a = 3.0/8.0; double b = 1.0/8.0;
    vector<vector<double> > Mtmp { {a, b, a, b, 0, 0, 0, 0, 0, 0, 0},
                                   {b, a, a, 0, 0, b, 0, 0, 0, 0, 0},
                                   {w1, w1, 1.0-N1*w1, w1, 0, w1, w1, 0, 0, 0, 0},
                                   {b, 0, a, a, 0, 0, b, 0, 0, 0, 0},
                                   {0, a, b, 0, b, a, 0, 0, 0, 0, 0},
                                   {0, b, a, 0, 0, a, b, 0, 0, 0, 0},
                                   {0, 0, a, b, 0, b, a, 0, 0, 0, 0},
                                   {0, 0, b, a, 0, 0, a, b, 0, 0, 0},
                                   {0, b, 0, 0, a, a, 0, 0, b, 0, 0},
                                   {0, w, w, 0, w, 1.0-N*w,  w, 0, w, w, 0},
                                   {0, 0, b, 0, 0, a, a, 0, 0, b, 0},
                                   {0, 0, w, w, 0, w, 1.0-N*w, w, 0, w, w},
                                   {0, 0, 0, b, 0, 0, a, a, 0, 0, b},
                                   {0, 0, 0, 0, b, a, 0, 0, a, b, 0},
                                   {0, 0, 0, 0, 0, a, b, 0, b, a, 0},
                                   {0, 0, 0, 0, 0, b, a, 0, 0, a, b},
                                   {0, 0, 0, 0, 0, 0, a, b, 0, b, a} };
    M = Mtmp;
    vector<vector<double>> SM1tmp(12,vector<double>(17,0.0)); SM1 = SM1tmp;
    vector<int> element1 { 2, 3, 5, 6, 7, 9, 10, 11, 12, 14, 15, 16 };
    for (int i = 0; i < 12; i++){
        SM1[i][element1[i]] = 1.0;
    }
    vector<vector<double>> SM2tmp(12,vector<double>(17,0.0)); SM2 = SM2tmp;
    vector<int> element2 { 4, 1, 9, 5, 2, 14, 10, 6, 3, 15, 11, 7 };
    for (int i = 0; i < 12; i++){
        SM2[i][element2[i]] = 1.0;
    }
    vector<vector<double>> SM3tmp(12,vector<double>(17,0.0)); SM3 = SM3tmp;
    vector<int> element3 { 1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15 };
    for (int i = 0; i < 12; i++){ 
        SM3[i][element3[i]] = 1.0;
    }
    vector<vector<double>> SM4tmp(11,vector<double>(17,0.0)); SM4 = SM4tmp;
    vector<int> element4 { 0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11 };
    for (int i = 0; i < 11; i++){ 
        SM4[i][element4[i]] = 1.0;
    }
}

/////////////////////////////////
//sphere
void setsphere_Loop_scheme(vector<Vertex>&  vertex, vector<Face>&  face, double& meanl, double r, double l){
    double a = r*2.0*sin(M_PI/5.0)/(sqrt(4.0*pow(sin(M_PI/5.0),2.0)-1.0)+0.5*cos(M_PI/5.0)); // a is the side length of icosahedron
    int n = round(log(a/l)/log(2.0));                                // division times to make the side-length as l.
    seticosahedron(vertex, face, r);                               // set an icosahedron for division into a sphere   
    for (int j = 0; j < n; j++){
        vector<Vertex> oldvertex = vertex;
        vector<Face> oldface = face;
        determine_AdjacentFace_for_vertex(oldvertex, oldface);
        determine_AdjacentVertex_for_vertex(oldvertex, oldface);

        int facenumber = oldface.size();
        vector<Face> newface(facenumber*3);
        for (int i = 0; i < facenumber; i++){
            int vertexnumber = vertex.size();
            ///////////////////////////////////////////////////////
            //  new vertices
            int nodea = oldface[i].AdjacentVertex[0]; int nodeb = oldface[i].AdjacentVertex[1];
            vector<int> dot1 = finddot( nodea, nodeb, oldvertex, oldface );
            vector<double> newvertex1 = 3.0/8.0*(oldvertex[nodea].Coord + oldvertex[nodeb].Coord) + 1.0/8.0*(oldvertex[dot1[0]].Coord + oldvertex[dot1[1]].Coord);
            nodea = oldface[i].AdjacentVertex[1]; nodeb = oldface[i].AdjacentVertex[2];
            vector<int> dot2 = finddot( nodea, nodeb, oldvertex, oldface );
            vector<double> newvertex2 = 3.0/8.0*(oldvertex[nodea].Coord + oldvertex[nodeb].Coord) + 1.0/8.0*(oldvertex[dot2[0]].Coord + oldvertex[dot2[1]].Coord);        
            nodea = oldface[i].AdjacentVertex[2]; nodeb = oldface[i].AdjacentVertex[0];
            vector<int> dot3 = finddot( nodea, nodeb, oldvertex, oldface );
            vector<double> newvertex3 = 3.0/8.0*(oldvertex[nodea].Coord + oldvertex[nodeb].Coord) + 1.0/8.0*(oldvertex[dot3[0]].Coord + oldvertex[dot3[1]].Coord);
            // check whether the new vertexs are actually those existed ones
            int n1 = 1; int n2 = 1; int n3 = 1;
            int new1, new2, new3;
            for (int k = 0; k < vertexnumber; k++) { 
                if ( norm(vertex[k].Coord-newvertex1) < 1e-5 ){
                    n1 = 0;
                    new1 = k;
                }
                if ( norm(vertex[k].Coord-newvertex2) < 1e-5 ){
                    n2 = 0;
                    new2 = k;
                }
                if ( norm(vertex[k].Coord-newvertex3) < 1e-5 ){
                    n3 = 0;
                    new3 = k;
                }
            }
            if (n1 == 1){                     // if the new vertex is totally new,not same as existed one
                new1 = n1 + vertexnumber - 1;       // the new vertex's number 
                Vertex vertextmp; vertextmp.Coord = newvertex1;
                vertex.push_back(vertextmp); // add the new vertex 
            }
            if (n2 == 1){
                new2 = n1 + n2 + vertexnumber - 1;
                Vertex vertextmp; vertextmp.Coord = newvertex2;
                vertex.push_back(vertextmp); // add the new vertex 
            }
            if (n3 == 1){
                new3 = n1 + n2 + n3 + vertexnumber - 1;
                Vertex vertextmp; vertextmp.Coord = newvertex3;
                vertex.push_back(vertextmp); // add the new vertex 
            }
            /////////////////////////////////////////////////
            // new face
            vector<int> temp1 {oldface[i].AdjacentVertex[0], new1, new3}; 
            newface[3*i].AdjacentVertex = temp1 ;
            vector<int> temp2 {oldface[i].AdjacentVertex[1], new2, new1};
            newface[3*i+1].AdjacentVertex = temp2;
            vector<int> temp3 {oldface[i].AdjacentVertex[2], new3, new2};
            newface[3*i+2].AdjacentVertex = temp3;
            vector<int> temp { new1, new2, new3 };
            face[i].AdjacentVertex = temp;
        }
        // update the positions of oldvertices
        for (int i = 0; i < oldvertex.size(); i++){
            int N = oldvertex[i].AdjacentVertex.size();
            // w=1/N*(5/8-(3/8+1/4*cos(2*pi/N))^2);
            double w = 3.0/8.0/N; // w=1/(N+3/8/w);
            oldvertex[i].Coord = (1.0 - N*w) * oldvertex[i].Coord;
            for (int k = 0; k < N; k++){
                oldvertex[i].Coord = oldvertex[i].Coord + w * oldvertex[oldvertex[i].AdjacentVertex[k]].Coord;
            }
        }
        // update face and vertex.
        oldface = face;
        face.resize(oldface.size() + newface.size());
        //face.rows 0 to oldface.size()-1 are not changed.
        for ( int i = 0; i < newface.size(); i++ ){
            int facenumber = i + oldface.size();
            face[facenumber] = newface[i];
        }
        for ( int i = 0; i < oldvertex.size(); i++ ){
            vertex[i].Coord = oldvertex[i].Coord;
        }
    }
    //////////////////////////////////////////////////////////
    // calculate the smallest and largest side-length; out put smallest and largest 
    vector<double> sidelength(3*face.size()); 
    double sum = 0.0;
    #pragma omp parallel for
    for (int i = 0; i < face.size(); i++){
        int node0 = face[i].AdjacentVertex[0];
        int node1 = face[i].AdjacentVertex[1];
        int node2 = face[i].AdjacentVertex[2];
        sidelength[3*i] = norm(vertex[node0].Coord - vertex[node1].Coord);
        sidelength[3*i+1] = norm(vertex[node0].Coord - vertex[node2].Coord);
        sidelength[3*i+2] = norm(vertex[node2].Coord - vertex[node1].Coord);
        sum += sidelength[3*i] + sidelength[3*i+1] + sidelength[3*i+2];
    }
    meanl = sum / sidelength.size();
}

void seticosahedron(vector<Vertex>& vertex, vector<Face>& face, double r){
    // build the vertex
    vector<double> t(5);
    for (int i = 0; i < 5; i++){ t[i] = i * (2.0*M_PI/5.0); }
    vector<vector<double>> vertex1(5,vector<double>(3));  // vertex of up pentagon
    for (int i = 0; i < 5; i++){ 
        vector<double> v { cos(t[i]), sin(t[i]), 0.0 };
        vertex1[i] = v;
    }
    for (int i = 0; i < 5; i++){ t[i] += M_PI/5.0; }
    double a = 2.0*sin(M_PI/5.0);   // side length
    vector<vector<double>> vertex2(5,vector<double>(3));    // vertex of down pentagon
    for (int i = 0; i < 5; i++){ 
        vector<double> v { cos(t[i]), sin(t[i]), -a*sqrt(3.0)/2.0 };
        vertex2[i] = v;
    }
    double h = sqrt(a*a-1.0);          // distance between upest point and pentagon
    vector<double> vertex3 {0.0, 0.0, h}; // upest vertex
    vector<double> vertex4 {0.0, 0.0, -a*sqrt(3.0)/2.0-h }; // downest vertex
    //mat vertex(12,3);
    vertex[0].Coord = vertex3;
    vertex[1].Coord = vertex1[0];
    vertex[2].Coord = vertex1[1];
    vertex[3].Coord = vertex1[2];
    vertex[4].Coord = vertex1[3];
    vertex[5].Coord = vertex1[4];
    vertex[6].Coord = vertex2[0];
    vertex[7].Coord = vertex2[1];
    vertex[8].Coord = vertex2[2];
    vertex[9].Coord = vertex2[3];
    vertex[10].Coord = vertex2[4];
    vertex[11].Coord = vertex4;
    // move the center to the origin
    for ( int i = 0; i < 12; i++) {
        vertex[i].Coord[2] += a*sqrt(3.0)/4.0; 
    }
    // move the vertex to the surface of the sphere 
    for (int i = 0; i < 12; i++) {
        vertex[i].Coord = vertex[i].Coord /norm(vertex[i].Coord) * r;  
    }
    // build the face
    vector<vector<int>> facetmp { {0, 1, 2},
                                  {0, 2, 3},
                                  {0, 3, 4},
                                  {0, 4, 5},
                                  {0, 5, 1},
                                  {1, 6, 2},
                                  {2, 7, 3},
                                  {3, 8, 4},
                                  {4, 9, 5},
                                  {5, 10, 1},
                                  {2, 6, 7},
                                  {3, 7, 8},
                                  {4, 8, 9},
                                  {5, 9, 10},
                                  {1, 10, 6},
                                  {11, 7, 6},
                                  {11, 8, 7},
                                  {11, 9, 8},
                                  {11, 10, 9},
                                  {11, 6, 10} };
    for ( int i = 0; i < 20; i++){
        face[i].AdjacentVertex = facetmp[i];
    }
} 

vector<int> finddot(int a, int b, vector<Vertex>& vertex, vector<Face>& face){
    vector<int> dot;
    for (int i = 0; i < vertex[a].AdjacentVertex.size(); i++){
        for (int j = 0; j < vertex[b].AdjacentVertex.size(); j++){
            if ( vertex[a].AdjacentVertex[i] == vertex[b].AdjacentVertex[j] ){
                dot.push_back(vertex[a].AdjacentVertex[i]);
            }
        }
    }
    if ( dot.size() != 2 ){
        cout<<"Wrong during finddot, shouldn't find number of dots "<<dot.size()<<endl;
        exit(0);
    }
    return dot;
}

//output nearby triangles that forms a parallelogram with the given point
//@TODO
/*
vector<vector<int>> findNearByTriangle(int pt, vector<Vertex>& vertex, vector<Face>& face) {
    //get adjacent vertex
    vector<int> * neighbor_vertex = &(vertex[pt].AdjacentVertex);
    //iterate through adj vertex
    for (int i = 0; i < neighbor_vertex->size; i++) {
        
        vector<int> * adjVertex = &(vertex[i].AdjacentVertex);
        if (std::find(adjVertex->begin, adjVertex->end, pt) != adjVertex->end) { //has element
            neighbor_vertex.push_back(i);
        }
    }
}
*/