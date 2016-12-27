#include "utlm.h"

using namespace std;

/*class mesh_face{

    double node_centre[2];
    double vertice[3];
    double area;
    double epsilon_relative;
    double mu_relative;
    double Z0;
};

class mesh_edge{

    double edge_vertice[2];
    int flip_index;
    double midpoint[2];
    double edgeLength;
    double linkLength;
    double Vlinki;
    double Vlinkr;
    double Vstub;
    double Ylink;
    double Ystub;


};*/
void creat_half_edge(const node_vec &mnode,
                     const faces &mface,
                     vector<edge> &edge_vec,
                     vector<int> &mesh_boundary){
    const int nEpf(3);
    const int no_faces(mface.eleVec.size());
    const int no_vet(mnode.nodex.size());
    const int no_edges(no_faces*nEpf);

    edge_vec.resize(no_edges);
    
    for(int inF=0;inF<no_faces;++inF){
        
        for(int inEpf=0;inEpf<nEpf;++inEpf){
            
            int v2((inEpf+1)%nEpf);
            int startVet(mface.eleVec[inF].ele_vet[inEpf]-1);
            int endVet(mface.eleVec[inF].ele_vet[v2]-1);
            
            int edgeIndex(inF*nEpf+inEpf);

//------------------Find faceId--and vertice index for each edge----
            edge_vec[edgeIndex].faceId=inF;
            edge_vec[edgeIndex].edgeVet[0]=startVet;
            edge_vec[edgeIndex].edgeVet[1]=endVet;
            
            int flip_edge_start(endVet);
            int flip_edge_end(startVet);

//------------------Find flip edge index-----------------------------
            for(int iFe=0;iFe<edgeIndex;++iFe){

                if(edge_vec[iFe].edgeVet[0]==flip_edge_start&&edge_vec[iFe].edgeVet[1]==flip_edge_end){
                    edge_vec[edgeIndex].edgeFlip=iFe;
                    edge_vec[iFe].edgeFlip=edgeIndex;        
                }
            }

//-------------------Find circumcentre for each face---------
            int v3((inEpf+2)%nEpf);
            int Vertice3(mface.eleVec[inF].ele_vet[v3]-1);

            double ax(mnode.nodex[startVet].node_vet[0]);
            double ay(mnode.nodex[startVet].node_vet[1]);

            double bx(mnode.nodex[endVet].node_vet[0]);
            double by(mnode.nodex[endVet].node_vet[1]);

            double cx(mnode.nodex[Vertice3].node_vet[0]);
            double cy(mnode.nodex[Vertice3].node_vet[1]);
            
            //cout<<"\nV1: "<<startVet<<"  V2: "<<endVet<<"  V3: "<<Vertice3;
            double tempD(2.0*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by)));

            double ccm_x(((ax*ax+ay*ay)*(by-cy)+(bx*bx+by*by)*(cy-ay)+(cx*cx+cy*cy)*(ay-by))/tempD);
            double ccm_y(((ax*ax+ay*ay)*(cx-bx)+(bx*bx+by*by)*(ax-cx)+(cx*cx+cy*cy)*(bx-ax))/tempD);

            edge_vec[edgeIndex].ccm[0]=ccm_x;
            edge_vec[edgeIndex].ccm[1]=ccm_y;

            //cout<<"\nFACE: "<<inF+1<<"     CCM: "<<edge_vec[edgeIndex].ccm[0]<<"    "<<edge_vec[edgeIndex].ccm[1];
            
//------------------Find Edge Length --------------
            double Ledge(sqrt((ax-bx)*(ax-bx)+(ay-by)*(ay-by)));
            edge_vec[edgeIndex].edgeLength=Ledge;
            
            double mpx((ax+bx)/2);
            double mpy((ay+by)/2);
            edge_vec[edgeIndex].midpoint[0]=mpx;
            edge_vec[edgeIndex].midpoint[1]=mpy;
        }    
    }
        map<int,int> _mesh_boundary;
    for(int iHe=0;iHe<no_edges;++iHe){
//--------------Create boundary mesh index vector--------
        
        if(edge_vec[iHe].edgeFlip==-1){
            double *vertice=new double[2];
            vertice=edge_vec[iHe].get_vertice();

            _mesh_boundary[vertice[0]]=iHe;

            delete[] vertice;
        }

//--------------Fill Link Length-------------------------
        int iFe(edge_vec[iHe].edgeFlip);
        double ccm_iHe_x(edge_vec[iHe].ccm[0]);
        double ccm_iHe_y(edge_vec[iHe].ccm[1]);

        double ccm_flip_x(edge_vec[iFe].ccm[0]);
        double ccm_flip_y(edge_vec[iFe].ccm[1]);

        double midpoint_x(edge_vec[iHe].midpoint[0]);
        double midpoint_y(edge_vec[iHe].midpoint[1]);

        double Llink(0);

        if(iFe==-1){
            double dx2((ccm_iHe_x-midpoint_x)*(ccm_iHe_x-midpoint_x));
            double dy2((ccm_iHe_y-midpoint_y)*(ccm_iHe_y-midpoint_y));

            Llink=sqrt(dx2+dy2);
        }
        else{
            double dx2((ccm_iHe_x-ccm_flip_x)*(ccm_iHe_x-ccm_flip_x));
            double dy2((ccm_iHe_y-ccm_flip_y)*(ccm_iHe_y-ccm_flip_y));

            Llink=sqrt(dx2+dy2)*0.5;
        }

        edge_vec[iHe].linkLength=Llink;

    }

    cout<<endl<<_mesh_boundary.size();

    for(map<int,int>::iterator it=_mesh_boundary.begin();it!=_mesh_boundary.end();++it){
            //cout<<"\nvertice index: "<<it->first<<"  edge index: "<<it->second;

        mesh_boundary.push_back(it->second);
        cout<<"\nBoundary index: "<<mesh_boundary.at(it->first);
    }

            
}

void min_edge_link_length(const vector<edge> &edge_vec,
                        double &_minEdge,
                        double &_minLink){
    _minEdge=1e9;
    _minLink=1e9;

    int shortestId(0);

    for(int it=0;it<edge_vec.size();++it){
        if(edge_vec[it].linkLength<_minLink){
            _minLink=edge_vec[it].linkLength;
            shortestId=it;
        }

        if(edge_vec[it].edgeLength<_minEdge){
            _minEdge=edge_vec[it].edgeLength;
        }
    }
}

void gaussian_wave_excite(const double &width, 
                        const double &delay, 
                        const double &t, 
                        double &Efield){
    
     const double t0(delay*width);
     double gamma(4*(t-t0)/(width));
     Efield=4/(width*constants::get_c0()*sqrt(constants::get_pi()))*exp(-gamma*gamma);

}


void calAdmittance( const double &dt,
                    faces &mface,
                    vector<edge> &edge_vector){
    
//---Calculate Addimittance--------------------
//---Y_LinkL, Y_stub_C, Y_Stub_L---------------------
//---Y_link_L=edge_length*del_t/(link_length*mu0*2)
//---Y_stub_C=epsilon*edge_length*link_length/del_t;
//---Y_stub_L=Y_link_L;
//---Y_stub=Y_stub_C - Y_stub_L
    const int no_edge(edge_vector.size());
    const int no_faces(mface.get_no_face());
    const int no_edge_pF(3);
    double mu0(constants::get_u0());
    double ep0(constants::get_e0());

    for(int i=0;i<no_faces;++i){
        
        double Ylink_total(0);
        double face_Z0(0);

        const double epsilonr(mface.get_epsilonr_with_id(i));

        for (int j=0;j<no_edge_pF;++j){
            const int iHe(i*3.0+j);
            double edge_length(edge_vector[iHe].edgeLength);
            double link_length(edge_vector[iHe].linkLength);

            double Y_link(edge_length*dt/(link_length*mu0*2));
            double Y_stub_C(edge_length*ep0*epsilonr*link_length/dt);
            double Y_stub(Y_stub_C - Y_link);

            edge_vector[iHe].Ylink=Y_link;
            edge_vector[iHe].Ystub=Y_stub;

            Ylink_total+=Y_link;
        }

        face_Z0=1/Ylink_total;

        mface.set_impedance(i,face_Z0);
    }
}

void set_inner_circle_different_face_number(
    faces &my_face,
    const vector<edge> &my_edges,
    const double &inner_circle_radius,
    const double inner_circle_centre[],
    const int &face_number){

    const int no_edge(my_edges.size());
    const int no_faces(my_face.eleVec.size());

    for(int iface=0;iface<no_faces;++iface){

        const int edge_id(iface*3+0);
        double* circm=new double[2];
        circm=my_edges[edge_id].get_ccm();
        double dist(sqrt((circm[0]-inner_circle_centre[0])*(circm[0]-inner_circle_centre[0])+(circm[1]-inner_circle_centre[1])*(circm[1]-inner_circle_centre[1])));

        if(dist<=inner_circle_radius){
            my_face.eleVec[iface].set_face_number(face_number);
        }
    }
}

void set_material_property(
    faces &my_face,
    const vector<double> &_epr,
    const vector<double> &_mur){
    
    const int no_face_material(my_face.find_no_material());
    const int no_faces(my_face.eleVec.size());
    for(int i=0;i<no_faces;++i){
        const int fnum_id(my_face.eleVec[i].get_fnum());
        my_face.eleVec[i].set_material(_epr[fnum_id-1],_mur[fnum_id-1]);
    }
}

void create_mesh_body_vector(
    const vector<edge> &my_edges,
    list<int> &mesh_body){

    const int no_edges(my_edges.size());
    for (int i=0;i<no_edges;++i){
        mesh_body.push_back(i);
    }

    list<int>::iterator ite=mesh_body.begin();
    while(ite!=mesh_body.end()){
        const int edge_id(*ite);
        const int flip_id(my_edges[edge_id].edgeFlip);
        mesh_body.remove(flip_id);

        ite++;

        if(flip_id==-1){mesh_body.remove(edge_id);}
    }
}

void create_reflection_coeff(
    const vector<edge> &my_edges,
    const vector<int> &my_bound_edges,
    const faces &my_faces,
    vector<double> &my_refl_coeff,
    vector<double> &Y_boundary,
    const vector<double> &my_condition){
//--boundary condition parameters:
//--matched boundary: 0;
//--short circuit (PEC): -1;
//--open circuit: 1;
    
    const int bound_edge_no(my_bound_edges.size());
    if(my_condition.size()!=bound_edge_no){
        cout<<"\nError: Boundary condtion size doesnot match boundary edge vector size"<<endl;

        return;
    }

    const int edge_no(my_edges.size());
    my_refl_coeff.reserve(edge_no);
    memset(&my_refl_coeff[0],0,sizeof(int)*edge_no);

    Y_boundary.reserve(edge_no);
    memset(&Y_boundary[0],0,sizeof(double)*edge_no);

    for(int i=0;i<bound_edge_no;++i){
        
        const double condition(my_condition[i]);
        const double edge_id(my_bound_edges[i]);

        if(condition==1.0)
        {
            my_refl_coeff[edge_id]=1.0;
            Y_boundary[edge_id]=0.0;
        }    
        else if(condition==-1.0)
        {
            my_refl_coeff[edge_id]=-1.0;
            Y_boundary[edge_id]=0.0;
        }
        else if(condition==0.0)
        {
            my_refl_coeff[edge_id]=0.0;
            const int face_id(my_edges[edge_id].get_face_id());

            const double epr(my_faces.get_epsilonr_with_id(face_id));

            const double Y0(sqrt(constants::get_e0()/constants::get_u0()));

            Y_boundary[edge_id]=Y0*(sqrt(epr*my_edges[edge_id].get_edge_length()));
        }
        else
        {
            my_refl_coeff[edge_id]=condition;
        }

    }
}

void scatter(const int &time_step, 
             vector<edge> &my_edges,                   
             faces &my_faces,
             map<int,int> &mesh_boundary)
{
    const int no_edge(my_edges.size());
    const int no_face(my_faces.get_no_face());
    const int no_edge_per_face(3);

    vector<double> nodeCurrent;
    nodeCurrent.reserve(no_face);
	memset(&nodeCurrent[0],0,sizeof(double)*no_face);

    vector<double> nodeVoltage;
    nodeVoltage.reserve(no_face);
	memset(&nodeVoltage[0],0,sizeof(double)*no_face);

    for(int kface=0;kface<no_face;++kface){

        for(int iepf=0;iepf<3;++iepf){
            const int edge_id(kface*3.0+iepf);
            nodeCurrent[kface]+=my_edges[edge_id].calc_Vlinki_times_Ylink();
        }
        
        nodeCurrent[kface]*=2.0;

        if(my_faces.get_epsilonr_with_id(kface)>1e8){
            nodeVoltage[kface]=0.0;
        }
        else{
            nodeVoltage[kface]=nodeCurrent[kface]*my_faces.get_Z0_with_id(kface);
        }

        my_edges[kface*3].set_Vlinkr(nodeVoltage[kface]);
        my_edges[kface*3+1].set_Vlinkr(nodeVoltage[kface]);
        my_edges[kface*3+2].set_Vlinkr(nodeVoltage[kface]);
    }

    //-----open circuit voltage and closed circuit current--

    for(map<int,int>::iterator it=mesh_boundary.begin();it!=mesh_boundary.end();++it){
        const int edge_id(it->second);
        const int face_id(my_edges[edge_id].get_face_id());
        const double face_epr(my_faces.get_epsilonr_with_id(face_id));
        
        double closeCurrent(0.0);
        double openVoltage(0.0);

        if(face_epr>1e8){
        //----PEC boundary----
            closeCurrent=0.0;
            openVoltage=my_edges[edge_id].get_Vlinkr();
        }
        else{
            const double linkVolt(my_edges[edge_id].get_Vlinkr());
            const double stubVolt(my_edges[edge_id].get_Vstub());
            const double linkY(my_edges[edge_id].get_Ylink());
            const double stubY(my_edges[edge_id].get_Ystub());

            closeCurrent=2.0*linkVolt*linkY+2.0*stubVolt*stubY;
            openVoltage=closeCurrent/(linkY+stubY);
        }
    }
}

/*
void connect(const int &time_step, vector<mesh_edge> &mesh_edges, vector<mesh_face> &mesh_body, const vector<int> &boundary, const vector<int> &mesh_internal)
{
    //how to decide an edge is on the boundary---------------------------
    const int no_edge(mesh_edges.size());
    const int no_face(mesh_body.size());
    const int no_bound_face(boundary.size());
    const int no_internal_nodes(mesh_internal.size());

    vector<double> Efield;
    Efield.reserve(no_edge);
    vector<double> Hfield;
    Hfield.reserve(no_edge);

    for(int it_bound_face=0;it_bound_face<no_bound_face;++it_bound_face){

        const int bound_edge_index(boundary[it_bound_face]);

        //how to define infinity Ystub or PEC?-------------------------
        if(mesh_edges[bound_edge_index].Ystub>1e8){
            Hfield[bound_edge_index]=0;
        }

    }

    for(int it_internal_nodes=0;it_internal_nodes<no_internal_nodes;++it_internal_nodes){

        const int internal_edge_index(mesh_internal[it_internal_nodes]);

        const int flip_edge_index(mesh_edges[internal_edge_index].flip_index);

        if(mesh_edges[internal_edge_index].Ystub>1e8||mesh_edges[flip_edge_index].Ystub>1e8){

            mesh_edges[internal_edge_index].Vlinki=0;
            mesh_edges[internal_edge_index].Vstub=0;
            mesh_edges[flip_edge_index].Vlinki=0;
            mesh_edges[flip_edge_index].Vstub=0;

        }
        else{
            double Ilinkr(mesh_edges[internal_edge_index].Vlinkr*mesh_edges[internal_edge_index].Ylink);
            double Ilinkr_flip(mesh_edges[flip_edge_index].Vlinkr*mesh_edges[flip_edge_index].Ylink);
            double Istub(mesh_edges[internal_edge_index].Vstub*mesh_edges[internal_edge_index].Ystub);
            double Istub_flip(mesh_edges[flip_edge_index].Vstub*mesh_edges[flip_edge_index].Ystub);

            double Iconnect(2*(Ilinkr+Ilinkr_flip+Istub+Istub_flip));
            double Ylink(mesh_edges[internal_edge_index].Ylink+mesh_edges[flip_edge_index].Ylink);
            double Ystub(mesh_edges[internal_edge_index].Ystub+mesh_edges[flip_edge_index].Ystub);
            double Ytotal(Ylink+Ystub);

            mesh_edges[internal_edge_index].Vlinki=Iconnect/Ytotal-mesh_edges[internal_edge_index].Vlinkr;
            mesh_edges[internal_edge_index].Vstub=Iconnect/Ytotal-mesh_edges[internal_edge_index].Vstub;
            mesh_edges[flip_edge_index].Vlinki=Iconnect/Ytotal-mesh_edges[flip_edge_index].Vlinkr;
            mesh_edges[flip_edge_index].Vstub=Iconnect/Ytotal-mesh_edges[flip_edge_index].Vstub;
        }
    }
}


void edge_excite(vector<mesh_edge> &mesh_edges, 
                const vector<int> &source_edge,
                const double &Vsource)
{

    int no_source_edge(source_edge.size());

    for(int sei=0;sei<no_source_edge;++sei){

        int source_edge_index(source_edge[sei]);

        int source_edge_flip(mesh_edges[source_edge_index].flip);

		double Ytotal(mesh_edges[source_edge_index].Ylink+mesh_edges[source_edge_index].Ystub);

		if(source_edge_flip!=0){

			Ytotal+=mesh_edges[source_edge_flip].Ylink+mesh_edges[source_edge_flip].Ystub;
	    }

		mesh_edges[source_edge_index].Vlinki+=Vsource/Ytotal;
		mesh_edges[source_edge_index].Vstub+=Vsource/Ytotal;
    }
}*/
