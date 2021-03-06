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

void gause_wave::evaluate(const vector<double> &time_array, vector<double> &Vsource_array){
    
    const int array_size(time_array.size());
    Vsource_array.reserve(array_size);

    const double t0(delay*width);
    const double pi(constants::get_pi());

    for(int i=0;i<array_size;++i){
        double t(time_array[i]);
        double gamma(4.0*(t-t0)/width);
        Vsource_array[i]=amplitude*(4.0/(width*wave_speed*sqrt(pi)))*exp(-gamma*gamma);
    }
}

void read_from_gmsh(const char filename[],
                    vector<node> &node_input,
                    vector<element> &face_input){

    cout<<"\nRead Data from gmsh file."<<endl;
    ifstream read_data(filename);

    if(!read_data){cerr<<"\nFile not open!";}

    string input_line;

    while(getline(read_data,input_line)){
       //cout<<"\nFile Open"<<endl; 
       //cout<<"\nInput: "<<input_line;
        
        if(input_line=="$MeshFormat"){
            cout<<"\nMesh Format read"<<endl;

            getline(read_data,input_line);
            double mesh_version, file_type, data_size;
            stringstream mesh_format(input_line);

            mesh_format>>mesh_version>>file_type>>data_size;
            cout<<"\nMesh Format Output \nversion: "<<mesh_version;
            cout<<"   file_type: "<<file_type;
            cout<<"   data_size: "<<data_size; 
            cout<<"\nMesh Format read finish.";
        }

        if(input_line=="$PhysicalNames"){
            cout<<"\nRead Face Property Number.";
            getline(read_data,input_line);
            int no_face_number;
            stringstream fnum_read(input_line);

            fnum_read>>no_face_number;
            
            for(int iFnum=0;iFnum<no_face_number;++iFnum){
                getline(read_data,input_line);
                int physical_dimension,face_number;
                string face_name;
                stringstream face_prop(input_line);
                face_prop>>physical_dimension>>face_number>>face_name;

                cout<<"\nFace Number Read. "<<physical_dimension<<"  "<<face_number<<"  "<<face_name;
            }
            
            getline(read_data,input_line);
            if(input_line!="$EndPhysicalNames"){
                cout<<"\nFace Number read not finish."<<input_line<<endl;
            }
            
            cout<<"\nFaceNumber Read Finish."<<endl;
        }

        if(input_line=="$Nodes"){
            
            getline(read_data,input_line);
            int node_size;
            istringstream read_node(input_line);
            read_node>>node_size;
            //cout<<"\nNode Size: "<<node_size;
            //node_input.reserve(node_size);

            for(int iNode=0;iNode<node_size;++iNode){
                getline(read_data,input_line);
                double node_x,node_y;
                int node_index;

                stringstream node_coord(input_line);
                node_coord>>node_index>>node_x>>node_y;

                //cout<<"\nNode Read -- ID: "<<node_index<<"  "<<node_x<<"  "<<node_y<<endl;

                node my_node(node_index,node_x,node_y);

                node_input.push_back(my_node);

            }

            getline(read_data,input_line);
            if(input_line!="$EndNodes"){
                cout<<"\nNode reading not Finish!"<<endl; 
                cout<<input_line<<endl;
            }
            
            cout<<"\nNode Read End."<<endl;
        }

        if(input_line=="$Elements"){
            
            int element_size;
            getline(read_data,input_line);
            stringstream read_element(input_line);
            read_element>>element_size;
            cout<<"\nElement Size: " <<element_size;

            //face_input.reserve(element_size);
            
            for(int iEle=0;iEle<element_size;++iEle){
                getline(read_data,input_line);
                //cout<<"\nInput line: "<<input_line;
                int ele_index,ele_type,no_tags,physical_entity,elementary_entity;
                int ele_vertex0,ele_vertex1,ele_vertex2;

                istringstream ele_input(input_line);
                ele_input>>ele_index>>ele_type>>no_tags>>physical_entity>>elementary_entity>>ele_vertex0>>ele_vertex1>>ele_vertex2;
                
                //cout<<"\nElement Read: "<<ele_index<<"  "<<ele_type<<"  "<<no_tags<<"  "<<physical_entity<<"  ";
                //cout<<elementary_entity<<"  "<<ele_vertex0<<"  "<<ele_vertex1<<"  "<<ele_vertex2<<endl;
                
                element my_element(ele_index,ele_vertex0,ele_vertex1,ele_vertex2);
                
                my_element.set_face_number(physical_entity);

                face_input.push_back(my_element);

                //cout<<"\nElements: "<<face_input[iEle];
            }

            getline(read_data,input_line);
            if(input_line!="$EndElements"){
                cout<<"\nElement reading not Finish!"<<endl;
                cout<<input_line<<endl;
            }
            
            cout<<"\nElement Read End."<<endl;
        }
    }
}

void create_half_edge(const node_vec &mnode,
                     const faces &mface,
                     vector<edge> &edge_vec,
                     vector<int> &mesh_boundary){
    cout<<"\nCreate edge vector from input Nodes and Faces"<<endl;
    
    if(mface.no_elements==0){
        cout<<"\nFace vector read fail."<<endl;
        return;
    }

    if(mnode.no_nodes==0){
        cout<<"\nNode vector read fail."<<endl;
        return;
    }

    const int nEpf(3);
    const int no_faces(mface.no_elements);
    const int no_vet(mnode.no_nodes);
    const int no_edges(no_faces*nEpf);

    edge_vec.resize(no_edges);
    
    for(int inF=0;inF<no_faces;++inF){
        
        for(int inEpf=0;inEpf<nEpf;++inEpf){
            
            int v2((inEpf+1)%nEpf);
            //vertex index starts from 0
            int startVet(mface.eleVec[inF].ele_vet[inEpf]-1);
            int endVet(mface.eleVec[inF].ele_vet[v2]-1);
            
            int edgeIndex(inF*nEpf+inEpf);
            //cout<<"\nEdge ID: "<<edgeIndex;
    //------------------Find faceId--and vertice index for each edge----
            edge_vec[edgeIndex].faceId=inF;   //starts from zero
            edge_vec[edgeIndex].edgeVet[0]=startVet;
            edge_vec[edgeIndex].edgeVet[1]=endVet;
            
            //cout<<"\nEdge vertex: "<<edge_vec[edgeIndex].edgeVet[0]<<"  "<<edge_vec[edgeIndex].edgeVet[0];

            int flip_edge_start(endVet);
            int flip_edge_end(startVet);

    //------------------Find flip edge index-----------------------------
            for(int iFe=0;iFe<edgeIndex;++iFe){

                if(edge_vec[iFe].edgeVet[0]==flip_edge_start&&edge_vec[iFe].edgeVet[1]==flip_edge_end){

                    //cout<<"\nEdge Index: "<<edgeIndex;//<<"  Flip ID: "<<iFe;
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

            //cout<<"\nBound Edge index: "<<iHe;
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

    cout<<endl<<"Mesh boundary size: "<<_mesh_boundary.size();

    for(map<int,int>::iterator it=_mesh_boundary.begin();it!=_mesh_boundary.end();++it){
            //cout<<"\nvertice index: "<<it->first<<"  edge index: "<<it->second;

        mesh_boundary.push_back(it->second);
        //cout<<"\nBoundary index: "<<mesh_boundary.at(it->first);
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
     double gamma(4.0*(t-t0)/(width));
     Efield=4.0/(width*constants::get_c0()*sqrt(constants::get_pi()))*exp(-gamma*gamma);
          
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
            const int iHe(i*3+j);
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

void create_PEC_bound_vector(
    const vector<edge> &my_edges,
    const faces &my_faces,
    vector<int> &pec_bound){

    const int no_edges(my_edges.size());
    const int no_faces(my_faces.get_no_face());

    for(int i=0;i<no_edges;++i){
        const int cur_face(my_edges[i].get_face_id());
        const int cur_fnum(my_faces.get_fnum_with_id(cur_face));
        const int flip_id(my_edges[i].get_flip_edge());
        const int flip_face(my_edges[flip_id].get_face_id());
        const int flip_fnum(my_faces.get_fnum_with_id(flip_face));

        const double face_epr(my_faces.get_epsilonr_with_id(cur_face));
        if(cur_fnum!=flip_fnum&&face_epr>1e8){
            pec_bound.push_back(i);
        }
    }
}

void create_mesh_body_vector(
    //remove edge ids that reverse the vertex id and belongs to the boundary 
    //(not material boundary) 
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

            Y_boundary[edge_id]=Y0*(sqrt(epr)*my_edges[edge_id].get_edge_length());
        }
        else
        {
            my_refl_coeff[edge_id]=condition;
        }

    }
}

void scatter(const int &time_step, 
             vector<edge> &my_edges,                   
             const faces &my_faces,
             const vector<int> &mesh_boundary){

    const int no_edge(my_edges.size());
    const int no_face(my_faces.get_no_face());
    const int no_edge_per_face(3);

    vector<double> nodeCurrent;
    nodeCurrent.reserve(no_face);
	memset(&nodeCurrent[0],0,sizeof(double)*no_face);

    vector<double> nodeVoltage;
    nodeVoltage.reserve(no_face);
	memset(&nodeVoltage[0],0,sizeof(double)*no_face);

    //-------------------------------------------------------------------
    //Node Voltage at the centre of the triangle-----
    //------V_node = 2.0*I_node*Z0_face (from my_faces_-----------------
    //------I_node = Sum(Vlink_i*Ylink) for each port of the node
    //-------------------------------------------------------------------
    for(int kface=0;kface<no_face;++kface){
        if(nodeCurrent[kface]!=0){
            cout<<"\nError - Node Current node initialised as 0.";
            return;
        }

        for(int iepf=0;iepf<3;++iepf){
            const int edge_id(kface*3.0+iepf);
            nodeCurrent[kface]+=my_edges[edge_id].calc_Vlinki_times_Ylink();
        }
        
        //nodeCurrent[kface]*=2.0;

        if(my_faces.get_epsilonr_with_id(kface)>1e8){
            nodeVoltage[kface]=0.0;
        }
        else{
            nodeVoltage[kface]=2.0*nodeCurrent[kface]*my_faces.get_Z0_with_id(kface);
        }

        //----reflected voltage for each port(edge)
        //----Vlinkr=V_node - Vlinki per edge
        my_edges[kface*3.0].set_Vlinkr(nodeVoltage[kface]-my_edges[kface*3].get_Vlinki());
        my_edges[kface*3+1].set_Vlinkr(nodeVoltage[kface]-my_edges[kface*3+1].get_Vlinki());
        my_edges[kface*3+2].set_Vlinkr(nodeVoltage[kface]-my_edges[kface*3.0+2].get_Vlinki());
    }

    //-----open circuit voltage and closed circuit current--
    double boundary_array_size(mesh_boundary.size()); //size of mesh boundary edge index vector

    for(int ibe=0;ibe<boundary_array_size;++ibe){   //ibe - iteration id for mesh boundary edge vector
        const int edge_id(mesh_boundary[ibe]);
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

void connect(const int &time_step, 
            vector<edge> &my_edges, 
            const faces &my_faces, 
            const vector<int> &my_bound_edges, 
            const vector<int> &my_pec_bounds,
            const vector<int> &my_inner_edges,
            const vector<double> &my_refl_coeff,
            const vector<double> &my_y_boundary){

    //how to decide an edge is on the boundary---------------------------
    const int no_edge(my_edges.size());
    const int no_face(my_faces.get_no_face());
    const int no_bound_edge(my_bound_edges.size());
    
    if(no_edge==0) {
        cout<<"\nEdge vector fails to send.";
        return;
    }
    if(no_face==0){
        cout<<"\nFace vector fails to send.";
        return;
    }
    if(no_bound_edge==0){
        cout<<"\nBoundary edge id vector fails to send.";
        return;
    }

    vector<double> Efield;
    Efield.reserve(no_edge);
    vector<double> Hfield;
    Hfield.reserve(no_edge);

    //connection at boundary edges ----------------------------------------------
    for(int it_bound_edge=0;it_bound_edge<no_bound_edge;++it_bound_edge){

        const int edge_id(my_bound_edges[it_bound_edge]);
        const int face_id(my_edges[edge_id].get_face_id());
        const double face_epr(my_faces.get_epsilonr_with_id(face_id));
        //how to define infinity Ystub or PEC?-------------------------
        const double reflection_coeff(my_refl_coeff[edge_id]);
        //--matched boundary - 
        const double vlinkr(my_edges[edge_id].get_Vlinkr());
        const double vlinki(my_edges[edge_id].get_Vlinki());
        const double vstub(my_edges[edge_id].get_Vstub());

        const double ylink(my_edges[edge_id].get_Ylink());
        const double ystub(my_edges[edge_id].get_Ystub());

        const double port_current(2.0*(vlinkr*ylink+vstub*ystub));

        if(reflection_coeff==0){
            //-----reflection coefficient zero => matched boundary admittance--
            //Ytotal=Ylink+Ystub+Yboundary---------------------------
            const double total_admit(ylink+ystub+my_y_boundary[edge_id]);

            const double port_voltage(port_current/total_admit);

            my_edges[edge_id].set_Vlinki(port_voltage-vlinkr);
            my_edges[edge_id].set_Vstub(port_voltage-vstub);
        }
        else{
            //speficied boundary condition or PEC(short circuit)
            //when PEC, refl_coeff=-1, simply reflect
            const double total_admit(ylink+ystub);

            const double port_voltage(port_current/total_admit);

            my_edges[edge_id].set_Vlinki(vlinkr*reflection_coeff);
            my_edges[edge_id].set_Vstub(vstub*reflection_coeff);
        }

    }

    //------PEC boundary connection process--------------
    if(my_pec_bounds.size()!=0){
        //pec boundary -> simply reflect
        cout<<"\nPEC Boundary connect";
        const int pec_bound_size(my_pec_bounds.size());
        for(int i=0;i<pec_bound_size;++i){

            int edge_id(my_pec_bounds[i]);
            const double v_linkr(my_edges[edge_id].get_Vlinkr());
            const double v_stub(my_edges[edge_id].get_Vstub());
            const double Ylink(my_edges[edge_id].get_Ylink());
            const double Ystub(my_edges[edge_id].get_Ystub());

            const double Iconnect(2.0*(v_linkr*Ylink+v_stub*Ystub));
            const double Ytotal(Ylink+Ystub);

            my_edges[edge_id].set_Vlinki(-1.0*v_linkr);
            my_edges[edge_id].set_Vstub(-1.0*v_stub);
        }
    }

    const int no_inner_edge(my_inner_edges.size());

    for(int it_inner_edge=0;it_inner_edge<no_inner_edge;++it_inner_edge){

        const int edge_id(my_inner_edges[it_inner_edge]);

        const int flip_edge_id(my_edges[edge_id].get_flip_edge());

        const int cur_face_id(my_edges[edge_id].get_face_id());
        const int flip_face_id(my_edges[flip_edge_id].get_face_id());

        const double cur_face_epr(my_faces.get_epsilonr_with_id(cur_face_id));
        const double flip_face_epr(my_faces.get_epsilonr_with_id(flip_face_id));

        const double edge_Vlinkr(my_edges[edge_id].get_Vlinkr());
        const double flip_Vlinkr(my_edges[flip_edge_id].get_Vlinkr());
        const double edge_Vstub(my_edges[edge_id].get_Vstub());
        const double flip_Vstub(my_edges[flip_edge_id].get_Vstub());

        if(cur_face_epr>1e8||flip_face_epr>1e8){
            //---PEC nodes----------------------------
            my_edges[edge_id].set_Vlinki(0.0);
            my_edges[edge_id].set_Vstub(0.0);
            my_edges[flip_edge_id].set_Vlinkr(0.0);
            my_edges[flip_edge_id].set_Vstub(0.0);

        }
        else{

            const double Ilinkr(my_edges[edge_id].get_Vlinkr()*my_edges[edge_id].get_Ylink());
            const double Ilinkr_flip(my_edges[flip_edge_id].get_Vlinkr()*my_edges[flip_edge_id].get_Ylink());
            const double Istub(my_edges[edge_id].get_Vstub()*my_edges[edge_id].get_Ystub());
            const double Istub_flip(my_edges[flip_edge_id].get_Vstub()*my_edges[flip_edge_id].get_Ystub());

            const double Iconnect(2.0*(Ilinkr+Ilinkr_flip+Istub+Istub_flip));
            const double Ylink(my_edges[edge_id].get_Ylink()+my_edges[flip_edge_id].get_Ylink());
            const double Ystub(my_edges[edge_id].get_Ystub()+my_edges[flip_edge_id].get_Ystub());
            const double Ytotal(Ylink+Ystub);

            const double connect_voltage(Iconnect/Ytotal);

            my_edges[edge_id].set_Vlinki(connect_voltage-edge_Vlinkr);
            my_edges[edge_id].set_Vstub(connect_voltage-edge_Vstub);
            my_edges[flip_edge_id].set_Vlinki(connect_voltage-flip_Vlinkr);
            my_edges[flip_edge_id].set_Vstub(connect_voltage-flip_Vstub);
        }
    }
}

void edge_excite(vector<edge> &mesh_edges, 
                const vector<int> &mesh_excite_edge_id,
                const double &Vsource){

    int no_source_edge(mesh_excite_edge_id.size());

    for(int sei=0;sei<no_source_edge;++sei){    //source edge id vector

        int source_edge_index(mesh_excite_edge_id[sei]);

        int source_edge_flip(mesh_edges[source_edge_index].get_flip_edge());

        double Ytotal(mesh_edges[source_edge_index].get_Ylink()+mesh_edges[source_edge_index].get_Ystub());

		if(source_edge_flip!=0){

            Ytotal+=mesh_edges[source_edge_flip].get_Ylink()+mesh_edges[source_edge_flip].get_Ystub();
	    }

        const double Vlinkr(mesh_edges[source_edge_index].get_Vlinkr());
        const double link_voltage_incident(mesh_edges[source_edge_index].get_Vlinki());
        const double stub_voltage(mesh_edges[source_edge_index].get_Vstub());
        mesh_edges[source_edge_index].set_Vlinki(link_voltage_incident+Vsource/Ytotal);
        mesh_edges[source_edge_index].set_Vstub(stub_voltage+Vsource/Ytotal);
    }
}

void run(vector<edge> &edge_vector,
         const faces &face_vector,
         const vector<int> &mesh_boundary_edge_id,
         const vector<int> &mesh_pec_boundary,
         const vector<int> &mesh_body,
         const vector<int> &mesh_excite_edge_id,
         const vector<double> &mesh_reflection_coefficient,
         const vector<double> &mesh_Y_boundary,
         const vector<double> &excite_field_td,
         const int &number_time_step){

    for(int k=0;k<number_time_step;++k){
        
        scatter(k,edge_vector,face_vector,mesh_boundary_edge_id);
        
        connect(k,edge_vector,face_vector,mesh_boundary_edge_id,mesh_pec_boundary,mesh_body,mesh_reflection_coefficient,mesh_Y_boundary);

        edge_excite(edge_vector,mesh_excite_edge_id,excite_field_td[k]);
        
    }
}
