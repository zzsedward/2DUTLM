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
    cout<<eleVecMine;

    /*vector<edge> edgeVector;
    //cout<<edgeVector[edgeIndex].edgeVet[0];
    creat_half_edge(nodeVecMine,eleVecMine,edgeVector);

    min_edge_link_length(edgeVector,minEdge,minLink);

    double del_t(minLink*sqrt(2*constants::get_e0()*constants::get_u0())/1.4);

    calAdmittance(del_t,eleVecMine,edgeVector);
    //cout<<"\nDel t: "<<del_t;

    double width(3.*minEdge*10/constants::get_c0());
    //cout<<"\nWidth: "<<width;

    for(int it=0;it<timestep;++it){
        double excite(0);
        times=it*del_t;
        gaussian_wave_excite(width,delay,times,excite);
        cout<<"\ntime step: "<<it<<"    Vsource: "<<excite;
    }
    cout<<endl;*/

}