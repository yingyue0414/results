#include "single_particle_diffusion.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// detailed definition

std::vector<double> TwelveShapeFunctions(std::vector<double>& vwu){
    double v = vwu[0]; double w = vwu[1]; double u = vwu[2];
    std::vector<double> sf(12); 
    sf[0] = 1.0/12.0*(pow(u,4.0) + 2.0*pow(u,3.0)*v);
    sf[1] = 1.0/12.0*(pow(u,4.0) + 2.0*pow(u,3.0)*w); 
    sf[2] = 1.0/12.0*(pow(u,4.0) + 2.0*pow(u,3.0)*w + 6.0*pow(u,3.0)*v + 6.0*pow(u,2.0)*v*w + 12.0*pow(u,2.0)*pow(v,2.0) + 6.0*u*pow(v,2.0)*w + 6.0*u*pow(v,3.0) + 2.0*pow(v,3.0)*w + pow(v,4.0));
    sf[3] = 1.0/12.0*(6.0*pow(u,4.0) + 24.0*pow(u,3.0)*w + 24.0*pow(u,2.0)*pow(w,2.0) + 8.0*u*pow(w,3.0) + pow(w,4.0) + 24.0*pow(u,3.0)*v + 60.0*pow(u,2.0)*v*w + 36.0*u*v*pow(w,2.0) + 6.0*v*pow(w,3.0) + 24.0*pow(u,2.0)*pow(v,2.0) + 36.0*u*pow(v,2.0)*w + 12.0*pow(v,2.0)*pow(w,2.0) + 8.0*u*pow(v,3.0) + 6.0*pow(v,3.0)*w + pow(v,4.0));
    sf[4] = 1.0/12.0*(pow(u,4.0) + 6.0*pow(u,3.0)*w + 12.0*pow(u,2.0)*pow(w,2.0) + 6.0*u*pow(w,3.0) + pow(w,4.0) + 2.0*pow(u,3.0)*v + 6.0*pow(u,2.0)*v*w + 6.0*u*v*pow(w,2.0) + 2.0*v*pow(w,3.0));
    sf[5] = 1.0/12.0*(2.0*u*pow(v,3.0) + pow(v,4.0)); 
    sf[6] = 1.0/12.0*(pow(u,4.0) + 6.0*pow(u,3.0)*w + 12.0*pow(u,2.0)*pow(w,2.0) + 6.0*u*pow(w,3.0)+ pow(w,4.0) + 8.0*pow(u,3.0)*v + 36.0*pow(u,2.0)*v*w + 36.0*u*v*pow(w,2.0) + 8.0*v*pow(w,3.0) + 24.0*pow(u,2.0)*pow(v,2.0) + 60.0*u*pow(v,2.0)*w + 24.0*pow(v,2.0)*pow(w,2.0) + 24.0*u*pow(v,3.0) + 24.0*pow(v,3.0)*w + 6.0*pow(v,4.0));
    sf[7] = 1.0/12.0*(pow(u,4.0) + 8.0*pow(u,3.0)*w + 24.0*pow(u,2.0)*pow(w,2.0) + 24.0*u*pow(w,3.0) + 6.0*pow(w,4.0) + 6.0*pow(u,3.0)*v + 36.0*pow(u,2.0)*v*w + 60.0*u*v*pow(w,2.0) + 24.0*v*pow(w,3.0) + 12.0*pow(u,2.0)*pow(v,2.0) + 36.0*u*pow(v,2.0)*w + 24.0*pow(v,2.0)*pow(w,2.0) + 6.0*u*pow(v,3.0) + 8.0*pow(v,3.0)*w + pow(v,4.0));
    sf[8] = 1.0/12.0*(2.0*u*pow(w,3.0) + pow(w,4.0)); 
    sf[9] = 1.0/12.0*(2.0*pow(v,3.0)*w + pow(v,4.0));
    sf[10] = 1.0/12.0*(2.0*u*pow(w,3.0) + pow(w,4.0) + 6.0*u*v*pow(w,2.0) + 6.0*v*pow(w,3.0) + 6.0*u*pow(v,2.0)*w + 12.0*pow(v,2.0)*pow(w,2.0) + 2.0*u*pow(v,3.0) + 6.0*pow(v,3.0)*w + pow(v,4.0));
    sf[11] = 1.0/12.0*(pow(w,4.0) + 2.0*v*pow(w,3.0)); 
    // 12 shape functions 
    // shape_functions(1,:), shape functions;  
    return sf;
}

void determine_CrossFace_for_particle(const std::vector<double>& Direction, double dL, const std::vector<double>& Coord, const std::vector<double>& PlaneNorm, std::vector<Vertex>& vertex, std::vector<Face>& face, std::vector<CrossFace>& CrossFaces, bool& hitBoundary){
    bool StopDiffusing = false;
    int iteration = 0;
    while ( StopDiffusing == false ){
        iteration ++;
        
        CrossFace CrossFacetmp;

        int newFaceIndex = -1;
        int node0 = CrossFaces[iteration-1].PairNodes[1][0];
        int node1 = CrossFaces[iteration-1].PairNodes[1][1];
        for (int i = 0; i < vertex[node0].AdjacentFace.size(); i++){
            for (int j = 0; j < vertex[node1].AdjacentFace.size(); j++){
                if ( vertex[node0].AdjacentFace[i] == vertex[node1].AdjacentFace[j] && vertex[node0].AdjacentFace[i] != CrossFaces[iteration-1].FaceIndex ){
                    newFaceIndex = vertex[node0].AdjacentFace[i];
                }
            }
        }
        if ( newFaceIndex == -1 ){
            cout<<"Wrong! No efficient cross-face is found!"<<endl;
            exit(0); 
        }

        if (face[newFaceIndex].IsBoundary == true){
            cout<<"Warning! Diffusing on boundary faces!"<<endl;
            hitBoundary = true;
            StopDiffusing = true;
            break;
        }

        CrossFacetmp.FaceIndex = newFaceIndex;

            node0 = face[newFaceIndex].AdjacentVertex[0];
            node1 = face[newFaceIndex].AdjacentVertex[1];
        int node2 = face[newFaceIndex].AdjacentVertex[2];

        int PairNumber = 0;
        if ( dot(vertex[node0].Coord - Coord, PlaneNorm) * dot(vertex[node1].Coord - Coord, PlaneNorm) < 0.0 ){
            // node0 - node1 passthrough the diffusing-plane 
            std::vector<int> pairnodes {node0, node1};
            CrossFacetmp.PairNodes[PairNumber] = pairnodes;
            // found the joint point of node0-node1 line and the diffusing-plane
            // ' lamda*node0 + (1-lamda)*node1 - Coord ' should be perpendicular to plane-norm, so get lamda value
            double lamda = dot(Coord - vertex[node1].Coord, PlaneNorm) / dot(vertex[node0].Coord - vertex[node1].Coord, PlaneNorm);
            std::vector<double> vwu {1.0-lamda, 0.0, lamda}; //std::vector<double> vwu {lamda, 1.0-lamda, 0.0}; // the joint point of FaceIndex
            CrossFacetmp.vwu[PairNumber] = vwu;
            PairNumber ++; 
        }
        if ( dot(vertex[node1].Coord - Coord, PlaneNorm) * dot(vertex[node2].Coord - Coord, PlaneNorm) < 0.0 ){
            // node1 - node2 passthrough the diffusing-plane 
            std::vector<int> pairnodes {node1, node2};
            CrossFacetmp.PairNodes[PairNumber] = pairnodes;
            // found the joint point of node1-node2 line and the diffusing-plane
            // ' lamda*node1 + (1-lamda)*node2 - Coord ' should be perpendicular to plane-norm, so get lamda value
            double lamda = dot(Coord - vertex[node2].Coord, PlaneNorm) / dot(vertex[node1].Coord - vertex[node2].Coord, PlaneNorm);
            std::vector<double> vwu {lamda, 1.0-lamda, 0.0}; //std::vector<double> vwu {0.0, lamda, 1.0-lamda}; // the joint point of FaceIndex
            CrossFacetmp.vwu[PairNumber] = vwu;
            PairNumber ++; 
        }
        if ( dot(vertex[node2].Coord - Coord, PlaneNorm) * dot(vertex[node0].Coord - Coord, PlaneNorm) < 0.0 ){
            // node2 - node0 passthrough the diffusing-plane 
            std::vector<int> pairnodes {node2, node0};
            CrossFacetmp.PairNodes[PairNumber] = pairnodes;
            // found the joint point of node2-node0 line and the diffusing-plane
            // ' lamda*node2 + (1-lamda)*node0 - Coord ' should be perpendicular to plane-norm, so get lamda value
            double lamda = dot(Coord - vertex[node0].Coord, PlaneNorm) / dot(vertex[node2].Coord - vertex[node0].Coord, PlaneNorm);
            std::vector<double> vwu {0.0, lamda, 1.0-lamda}; //std::vector<double> vwu {1.0-lamda, 0.0, lamda}; // the joint point of FaceIndex
            CrossFacetmp.vwu[PairNumber] = vwu;
            PairNumber ++;
        }
        if ( PairNumber != 2 ){
            cout<<"Determine crossface. Wrong! Too many joints on one CrossFace! pairnumber = "<<PairNumber<<endl;
            exit(0);
        }
        // calcualte the arc-distance across this face
        std::vector<double> JointPoint1(3);
        std::vector<double> JointPoint2(3);
        int n = face[newFaceIndex].OneRingVertex.size();
        std::vector<std::vector<double>> dots(n,std::vector<double>(3,0.0)); // one ring vertices
        if ( n == 12 ){ // regular patch, use one-ring to calcualte JointPoints.
            for (int j = 0; j < 12; j++) {
                int nodenum = face[newFaceIndex].OneRingVertex[j];
                dots[j] = vertex[nodenum].Coord;
            }
            JointPoint1 = TwelveShapeFunctions(CrossFacetmp.vwu[0]) * dots;
            JointPoint2 = TwelveShapeFunctions(CrossFacetmp.vwu[1]) * dots;
         }else if ( n == 11 ){ // irregular patch, use three-adjacent-nodes to calcualte JointPoints.
            JointPoint1 = CrossFacetmp.vwu[0][2]*vertex[node0].Coord + CrossFacetmp.vwu[0][0]*vertex[node1].Coord + CrossFacetmp.vwu[0][1]*vertex[node2].Coord;
            JointPoint2 = CrossFacetmp.vwu[1][2]*vertex[node0].Coord + CrossFacetmp.vwu[1][0]*vertex[node1].Coord + CrossFacetmp.vwu[1][1]*vertex[node2].Coord;
        }
        double distance = norm(JointPoint1 - JointPoint2);
        double ArcDistance;
        if ( face[newFaceIndex].MeanCurvature < 1.0e-8 ){ // flate face, no curvature
            ArcDistance = distance;
        }else{
            double R = 1.0 / face[newFaceIndex].MeanCurvature;
            ArcDistance = R * 2.0*asin(distance/2.0/R);
        }
        CrossFacetmp.CrossDistance = ArcDistance;

        // make sure the the first pairnodes are closer to the particle!
        /*
        if ( norm(JointPoint1-Coord) > norm(JointPoint2-Coord) ){ // swap these two nodes
            swap(CrossFacetmp.PairNodes[0], CrossFacetmp.PairNodes[1]);
            swap(CrossFacetmp.vwu[0], CrossFacetmp.vwu[1]);
        }
        */
        if ( ( CrossFacetmp.PairNodes[1][0] == CrossFaces[iteration-1].PairNodes[1][0] &&  CrossFacetmp.PairNodes[1][1] == CrossFaces[iteration-1].PairNodes[1][1] ) 
            || ( CrossFacetmp.PairNodes[1][0] == CrossFaces[iteration-1].PairNodes[1][1] &&  CrossFacetmp.PairNodes[1][1] == CrossFaces[iteration-1].PairNodes[1][0] ) ){
            swap(CrossFacetmp.PairNodes[0], CrossFacetmp.PairNodes[1]);
            swap(CrossFacetmp.vwu[0], CrossFacetmp.vwu[1]);
        }

        CrossFaces.push_back(CrossFacetmp);

        // check whether diffusion far enough
        double sumDistance = 0.0;
        for (int j = 0; j < CrossFaces.size(); j++){
            sumDistance += CrossFaces[j].CrossDistance;
        }
        if ( sumDistance > dL ){
            StopDiffusing = true;
        }
    }
}

void single_particles_diffusing_on_surface(std::vector<Vertex>& vertex, std::vector<Face>& face, Param& param){
    int CopyNumber = 1;
    //std::vector<double> TitrationCenter {param.sideX, param.sideY, 0.0};
    int FaceIndex = 0;//625;//determine_faceindex(vertex, face, TitrationCenter);
    std::vector<double> vwuTitration {1.0/3.0, 1.0/3.0, 1.0/3.0};
    int TimeSteps = 1E5; // us
    double dt =  1; // us
    double D = 1.0;  // nm2/us
    ////////////////////////////////////////////////////////
    // initialize all the particles
    std::vector<Particle> particles(CopyNumber);
    {
        std::vector<std::vector<double>> dots(12,std::vector<double>(3,0.0)); // one ring vertices
        for (int j = 0; j < 12; j++) {
            int nodenum = face[FaceIndex].OneRingVertex[j];
            dots[j] = vertex[nodenum].Coord;
        }
        std::vector<double> Coord = TwelveShapeFunctions(vwuTitration) * dots;
        for (int i = 0; i < CopyNumber; i++){
            particles[i].Index = i;
            particles[i].D = D;
            particles[i].Coord = Coord;
            particles[i].FaceIndex = FaceIndex;
        }
    }
    
    ////////////////////////////////////////////////
    for (int iTime = 0; iTime < TimeSteps; iTime++){
        cout<<"Time step = "<<iTime<<endl;
        #pragma omp parallel for
        for(int i = 0; i < CopyNumber; i++){
            // define the diffusion displacement
            double D = particles[i].D;
            int FaceIndex = particles[i].FaceIndex;
            int node0 = face[FaceIndex].AdjacentVertex[0];
            int node1 = face[FaceIndex].AdjacentVertex[1];
            int node2 = face[FaceIndex].AdjacentVertex[2];
            // position on triangle, note it is not the same as particle.Coord
            std::vector<double> Coord = particles[i].vwu[2]*vertex[node0].Coord + particles[i].vwu[0]*vertex[node1].Coord + particles[i].vwu[1]*vertex[node2].Coord;
            
            double dx = sqrt(2.0*D*dt) * gr();
            double dy = sqrt(2.0*D*dt) * gr();
            double dL = sqrt(dx*dx + dy*dy);
            // define the diffusion direction
            
            std::vector<double> Norm = face[FaceIndex].normVector;
            std::vector<double> tmpv1 = cross(Coord,Norm); tmpv1 = tmpv1/norm(tmpv1);
            std::vector<double> tmpv2 = cross(Norm,tmpv1); tmpv2 = tmpv2/norm(tmpv2); // Norm, tmpv1, tmpv2 are local coordinates
            // rotate tmpv1 on tmpv1-tmpv2 plane, then take it as the diffusing direction.
            double theta = 2.0*M_PI*((double)rand()/RAND_MAX);
            std::vector<double> Direction = cos(theta)*tmpv1 + sin(theta)*tmpv2; Direction = Direction / norm(Direction);
            std::vector<double> PlaneNorm = cross(Norm, Direction); // norm vector of the plane
            // PlaneNorm, Direction, and Coord define the plane that the particle will diffuse on. This plane-surface joint curve-line is the trajectory.           
            ////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////
            // find out all the faces that the particle needs to diffuse through
            std::vector<CrossFace> CrossFaces; // store the faces that the particle needs to diffuse through

            // for each crossface, there should be two pairs of nodes, because each crossface will cross the diffusing-plane by two sides.
            // However, for the very first crossface, just need to store one pair of nodes. 
            CrossFace CrossFacetmp;
            CrossFacetmp.FaceIndex = FaceIndex;
            int PairNumber = 0;
            if ( dot(vertex[node0].Coord - Coord, PlaneNorm) * dot(vertex[node1].Coord - Coord, PlaneNorm) < 0.0 ){
                // node0 - node1 passthrough the diffusing-plane 
                std::vector<int> pairnodes {node0, node1};
                CrossFacetmp.PairNodes[PairNumber] = pairnodes;
                // found the joint point of node0-node1 line and the diffusing-plane
                // ' lamda*node0 + (1-lamda)*node1 - Coord ' should be perpendicular to plane-norm, so get lamda value
                double lamda = dot(Coord - vertex[node1].Coord, PlaneNorm) / dot(vertex[node0].Coord - vertex[node1].Coord, PlaneNorm);
                std::vector<double> vwu {1.0-lamda, 0.0, lamda}; // the joint point of FaceIndex. NOTE the vwu order!!!!!!!
                CrossFacetmp.vwu[PairNumber] = vwu;
                PairNumber ++;
            }
            if ( dot(vertex[node1].Coord - Coord, PlaneNorm) * dot(vertex[node2].Coord - Coord, PlaneNorm) < 0.0 ){              
                // node1 - node2 passthrough the diffusing-plane 
                std::vector<int> pairnodes {node1, node2};
                CrossFacetmp.PairNodes[PairNumber] = pairnodes;
                // found the joint point of node1-node2 line and the diffusing-plane
                // ' lamda*node1 + (1-lamda)*node2 - Coord ' should be perpendicular to plane-norm, so get lamda value
                double lamda = dot(Coord - vertex[node2].Coord, PlaneNorm) / dot(vertex[node1].Coord - vertex[node2].Coord, PlaneNorm);
                std::vector<double> vwu {lamda, 1.0-lamda, 0.0}; // the joint point of FaceIndex
                CrossFacetmp.vwu[PairNumber] = vwu;
                PairNumber ++;
            }
            if ( dot(vertex[node2].Coord - Coord, PlaneNorm) * dot(vertex[node0].Coord - Coord, PlaneNorm) < 0.0 ){        
                // node2 - node0 passthrough the diffusing-plane 
                std::vector<int> pairnodes {node2, node0};
                CrossFacetmp.PairNodes[PairNumber] = pairnodes;
                // found the joint point of node2-node0 line and the diffusing-plane
                // ' lamda*node2 + (1-lamda)*node0 - Coord ' should be perpendicular to plane-norm, so get lamda value
                double lamda = dot(Coord - vertex[node0].Coord, PlaneNorm) / dot(vertex[node2].Coord - vertex[node0].Coord, PlaneNorm);
                std::vector<double> vwu {0.0, lamda, 1.0-lamda}; // the joint point of FaceIndex
                CrossFacetmp.vwu[PairNumber] = vwu;
                PairNumber ++;
            }
            if ( PairNumber != 2 ){
                cout<<"Current face. Wrong! Too many joints on one CrossFace! pairnumber = "<<PairNumber<<endl;
                //cout<<"nodes = "<<node0<<", "<<node1<<", "<<node2<<endl;
                //cout<<"vwu = "<<particles[i].vwu[0]<<", "<<particles[i].vwu[1]<<", "<<particles[i].vwu[2]<<endl;
                //cout<<"face is regualr? "<<face[FaceIndex].OneRingVertex.size()<<endl;
                //cout<<dot(vertex[node0].Coord - Coord, PlaneNorm)<<","<<dot(vertex[node1].Coord - Coord, PlaneNorm)<<","<<dot(vertex[node2].Coord - Coord, PlaneNorm)<<endl;
                exit(0);
            }
            // determine which joint is toward the diffusing direction! And put it on the second position in PairNodes and vwu.
            std::vector<double> JointPoint1(3);
            std::vector<double> JointPoint2(3);
            int n = face[FaceIndex].OneRingVertex.size();
            std::vector<std::vector<double>> dots(12,std::vector<double>(3,0.0)); // one ring vertices
            if ( n == 12 ){ // regular patch, use one-ring to calcualte JointPoints.
                for (int j = 0; j < 12; j++) {
                    int nodenum = face[FaceIndex].OneRingVertex[j];
                    dots[j] = vertex[nodenum].Coord;
                }
                JointPoint1 = TwelveShapeFunctions(CrossFacetmp.vwu[0]) * dots;
                JointPoint2 = TwelveShapeFunctions(CrossFacetmp.vwu[1]) * dots;
            }else if ( n == 11 ){ // irregular patch, use three-adjacent-nodes to calcualte JointPoints.
                JointPoint1 = CrossFacetmp.vwu[0][2]*vertex[node0].Coord + CrossFacetmp.vwu[0][0]*vertex[node1].Coord + CrossFacetmp.vwu[0][1]*vertex[node2].Coord;
                JointPoint2 = CrossFacetmp.vwu[1][2]*vertex[node0].Coord + CrossFacetmp.vwu[1][0]*vertex[node1].Coord + CrossFacetmp.vwu[1][1]*vertex[node2].Coord;
            }
            if ( dot(JointPoint1-Coord,Direction) > 0.0 ){
                swap(CrossFacetmp.PairNodes[0], CrossFacetmp.PairNodes[1]);
                swap(CrossFacetmp.vwu[0], CrossFacetmp.vwu[1]);
                CrossFacetmp.PairNodes[0].clear();
                CrossFacetmp.vwu[0].clear();
            }else if ( dot(JointPoint2-Coord,Direction) > 0.0 ){
                CrossFacetmp.PairNodes[0].clear();
                CrossFacetmp.vwu[0].clear();
            }
            std::vector<double> JointPoint(3);
            if ( n == 12 ){
                JointPoint = TwelveShapeFunctions(CrossFacetmp.vwu[1]) * dots;
            }else if ( n == 11 ){
                JointPoint = CrossFacetmp.vwu[1][2]*vertex[node0].Coord + CrossFacetmp.vwu[1][0]*vertex[node1].Coord + CrossFacetmp.vwu[1][1]*vertex[node2].Coord;
            }

            double distance = norm(JointPoint - particles[i].Coord);
            double ArcDistance;
            if ( face[FaceIndex].MeanCurvature < 1e-8 ){ // flate face, no curvature
                ArcDistance = distance;
            }else{
                double R = 1.0 / face[FaceIndex].MeanCurvature;
                ArcDistance = R * 2.0*asin(distance/2.0/R);
            }
            CrossFacetmp.CrossDistance = ArcDistance;

            CrossFaces.push_back(CrossFacetmp);
            // find out all next cross-faces
            if ( CrossFacetmp.CrossDistance < dL ){
                bool hitBoundary = false;
                determine_CrossFace_for_particle(Direction, dL, Coord, PlaneNorm, vertex, face, CrossFaces, hitBoundary);
                if ( hitBoundary == true ){
                    //cout<<"Hit the boundary, Force to not diffuse in this time step!"<<endl;
                    CrossFaces.clear();
                    CrossFaces.push_back(CrossFacetmp);
                    dL = 0.0;
                }
            }
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            // determine the coord in the last cross-face.
            CrossFace CrossFaceTarget = CrossFaces.back();
            particles[i].FaceIndex = CrossFaceTarget.FaceIndex;
            node0 = face[CrossFaceTarget.FaceIndex].AdjacentVertex[0];
            node1 = face[CrossFaceTarget.FaceIndex].AdjacentVertex[1];
            node2 = face[CrossFaceTarget.FaceIndex].AdjacentVertex[2];
            std::vector<double> vwuTarget(3);
            std::vector<double> coordTarget(3);
            n = face[CrossFaceTarget.FaceIndex].OneRingVertex.size();
            if ( n == 12 ){ // one ring vertices
                for (int j = 0; j < 12; j++) {
                    int nodenum = face[CrossFaceTarget.FaceIndex].OneRingVertex[j];
                    dots[j] = vertex[nodenum].Coord;
                }
            } 
            if ( CrossFaces.size() == 1 ){ // particle doesn't diffuse out of the current face
                cout<<"crossface.size() = "<<CrossFaces.size()<<endl;
                cout<<"dL = "<<dL<<", lastdist = "<<CrossFaceTarget.CrossDistance<<endl;
                if ( CrossFaceTarget.CrossDistance > dL ){
                    double lamda;
                    if ( face[CrossFaceTarget.FaceIndex].MeanCurvature > 1e-5 ){
                        double rTarget =  1.0/ face[CrossFaceTarget.FaceIndex].MeanCurvature;
                        double thetaTarget = CrossFaceTarget.CrossDistance / rTarget;
                        double ratio = dL/ CrossFaceTarget.CrossDistance;
                        double angle = thetaTarget * ratio;
                        std::vector<double> JointPoint(3);
                        if ( n == 12 ) {
                            JointPoint = TwelveShapeFunctions(CrossFaceTarget.vwu[1]) * dots;
                        }else if ( n == 11 ){
                            JointPoint = CrossFaceTarget.vwu[1][2]*vertex[node0].Coord + CrossFacetmp.vwu[1][0]*vertex[node1].Coord + CrossFacetmp.vwu[1][1]*vertex[node2].Coord;
                        }
                        //lamda = rTarget * sin(angle) / norm(JointPoint-Coord); // angle from currect Coord;
                        lamda = sin(angle)/( sin(angle) + sin(thetaTarget-angle) );
                    }else{
                        lamda = dL/ CrossFaceTarget.CrossDistance;
                    }
                    // barycentric coordinates in the last cross-face. Note the order of nodes.
                    vwuTarget = lamda *  CrossFaceTarget.vwu[1] + (1.0-lamda)*particles[i].vwu;
                    if ( n == 12 ) {
                        coordTarget = TwelveShapeFunctions(vwuTarget) * dots;
                    }else if ( n == 11 ){
                        coordTarget = vwuTarget[2]*vertex[node0].Coord + vwuTarget[0]*vertex[node1].Coord + vwuTarget[1]*vertex[node2].Coord;
                    }
                }else{
                    cout<<"Wrong! CrossFaces haven't been calculated correctly, because the cross-distance is not correct!"<<endl;
                    exit(0);
                }
            }else if ( CrossFaces.size() > 1 ){
                double sumDistance = 0.0;
                for (int j = 0; j < CrossFaces.size(); j++){
                    sumDistance += CrossFaces[j].CrossDistance;
                }
                cout<<"crossface.size() = "<<CrossFaces.size()<<endl;
                cout<<"dL = "<<dL<<", sumdist = "<<sumDistance<<", lastdist = "<<CrossFaceTarget.CrossDistance<<endl;
                if ( sumDistance > dL && (sumDistance - CrossFaceTarget.CrossDistance < dL) ){
                    double lamda;
                    if ( face[CrossFaceTarget.FaceIndex].MeanCurvature > 1e-5 ){
                        double rTarget =  1.0/ face[CrossFaceTarget.FaceIndex].MeanCurvature;
                        double thetaTarget = CrossFaceTarget.CrossDistance / rTarget;
                        double ratio = ( dL-(sumDistance-CrossFaceTarget.CrossDistance) )/ CrossFaceTarget.CrossDistance;
                        double angle = thetaTarget * ratio;
                        std::vector<double> JointPoint1(3);
                        std::vector<double> JointPoint2(3);
                        if ( n == 12 ){
                            JointPoint1 = TwelveShapeFunctions(CrossFaceTarget.vwu[0]) * dots;
                            JointPoint2 = TwelveShapeFunctions(CrossFaceTarget.vwu[1]) * dots;
                        }else if ( n == 11 ){
                            JointPoint1 = CrossFaceTarget.vwu[0][2]*vertex[node0].Coord + CrossFaceTarget.vwu[0][0]*vertex[node1].Coord + CrossFaceTarget.vwu[0][1]*vertex[node2].Coord;
                            JointPoint2 = CrossFaceTarget.vwu[1][2]*vertex[node0].Coord + CrossFaceTarget.vwu[1][0]*vertex[node1].Coord + CrossFaceTarget.vwu[1][1]*vertex[node2].Coord;
                        }
                        //lamda = rTarget * sin(angle) / norm(JointPoint1-JointPoint2); // angle from joint point of 0;
                        lamda = sin(angle)/( sin(angle) + sin(thetaTarget-angle) );
                    }else{
                        lamda = ( dL-(sumDistance-CrossFaceTarget.CrossDistance) )/ CrossFaceTarget.CrossDistance;
                    }
                    // barycentric coordinates in the last cross-face. Note the order of nodes.
                    vwuTarget = lamda *  CrossFaceTarget.vwu[1] + (1.0-lamda)*CrossFaceTarget.vwu[0];
                    if ( n == 12 ){
                        coordTarget = TwelveShapeFunctions(vwuTarget) * dots;
                    }else if ( n == 11 ){
                        coordTarget = vwuTarget[2]*vertex[node0].Coord + vwuTarget[0]*vertex[node1].Coord + vwuTarget[1]*vertex[node2].Coord;
                    }
                }else{
                    cout<<"Wrong! CrossFaces haven't been calculated correctly, because the cross-distance is not correct!"<<endl;
                    exit(0);
                }
            }

            particles[i].vwu = vwuTarget;
            particles[i].Coord = coordTarget;
        }
        /*
        if ( iTime % 100 == 0 ){    
            int kk = iTime/100;
            char filename[20] = "positions%d.csv";
            sprintf(filename,"positions%d.csv",kk);
            ofstream outfile(filename);
            for (int j = 0; j < particles.size(); j++) {
                outfile << particles[j].Coord[0] <<','<<particles[j].Coord[1]<<','<<particles[j].Coord[2] << '\n';
            }
            outfile.close();
        }
        */
    }
}