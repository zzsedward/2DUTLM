#include "utlm.h"

int main(int argc, char* argv[])
{
    cout<<"Hello World!"<<endl;

    int timestep(2e3);
    
    double node_no_vertice;
    int node_dim;
    int node_attri;
    int node_bound_marker;
    double times(0.);
    double minLink(0.),minEdge(0.);
    double Efield(0.);
    double delay(1.);
    
    //ifstream fin("cyl_res20.node");
    //fin>>node_dim>>node_attri>>node_bound_marker;

    //node node1(fin);
    //cout<<node1;

    node_vec nodeVecMine("cyl_res20.node");
    //cout<<nodeVecMine;

    faces eleVecMine("cyl_res20.ele");
    //cout<<eleVecMine;

    vector<edge> edgeVector;   
    //cout<<edgeVector[edgeIndex].edgeVet[0];
    vector<int> boundaryVector;

    creat_half_edge(nodeVecMine,eleVecMine,edgeVector,boundaryVector);
    
    for(vector<int>::iterator it=boundaryVector.begin();
    
    it!=boundaryVector.end();++it){
            cout<<"\nbound edge index: "<<*it;
    }

    min_edge_link_length(edgeVector,minEdge,minLink);

    double del_t(minLink*sqrt(2*constants::get_e0()*constants::get_u0())/1.4);

    calAdmittance(del_t,eleVecMine,edgeVector);
    //cout<<"\nDel t: "<<del_t;

    const double r_in(0.5);
    const double origin[2]={0.0,0.0};
    const int inner_fnum(2);

    set_inner_circle_different_face_number(eleVecMine,edgeVector,r_in,origin,inner_fnum);
    //cout<<eleVecMine;

    vector<double>epsilon_r;
    epsilon_r.push_back(1.0);
    epsilon_r.push_back(1e10);

    vector<double>mu_r;
    mu_r.push_back(1.0);
    mu_r.push_back(1.0);

    set_material_property(eleVecMine,epsilon_r,mu_r);
    //cout<<eleVecMine;

    vector<int> pec_boundary_edges;
    create_PEC_bound_vector(edgeVector,eleVecMine,pec_boundary_edges);

    list<int> mesh_body_list;
    create_mesh_body_vector(edgeVector,mesh_body_list);

    /*for(list<int>::iterator it=mesh_body_list.begin();it!=mesh_body_list.end();++it){
        cout<<"\n"<<*it;
    }*/

    vector<double> reflect_coeff;
    vector<double> Y_boundary;

    vector<double> bound_condition;
    bound_condition.reserve(boundaryVector.size());
    memset(&bound_condition[0],-1.0,sizeof(double)*boundaryVector.size());

    create_reflection_coeff(edgeVector,boundaryVector,eleVecMine,reflect_coeff,Y_boundary,bound_condition);

    double width(3.*minEdge*10/constants::get_c0());
    //cout<<"\nWidth: "<<width;

    /*for(int it=0;it<timestep;++it){
        double excite(0);
        times=it*del_t;
        gaussian_wave_excite(width,delay,times,excite);
        cout<<"\ntime step: "<<it<<"    Vsource: "<<excite;
    }
    cout<<endl;*/

}


