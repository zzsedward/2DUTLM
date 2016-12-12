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
                     vector<edge> &edge_vec){
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
                map<int,int> mesh_boundary;
            for(int iHe=0;iHe<no_edges;++iHe){
//--------------Create boundary mesh index vector--------
                
                if(edge_vec[iHe].edgeFlip==-1){

                    mesh_boundary[edge_vec[iHe].edgeVet[0]]=iHe;
                }

                // for(int i=0;i<mesh_boundary.size();++i){
                //     cout<<"\nedge index: "<<mesh_boundary[i];
                // }
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
                cout<<mesh_boundary.size();
            for(map<int,int>::iterator it=mesh_boundary.begin();it!=mesh_boundary.end();++it){
                    cout<<"\nvertice index: "<<it->first<<"  edge index: "<<it->second;
                }
//--------------Fill Addimittance---------------------------
     /*     for(int iHe=0;iHe<no_edges;++iHe){
                cout<<"\nEdge: "<<iHe;
                cout<<"   Edge Length: "<<edge_vec[iHe].edgeLength;
                cout<<"   Link Length: "<<edge_vec[iHe].linkLength<<endl;
            }*/

//---------
            
}

int main(int argc, char* argv[])
{
    cout<<"Hello World!"<<endl;
    int timestep(1e4);
    double del_t(1e-10);
    double node_no_vertice;
    int node_dim;
    int node_attri;
    int node_bound_marker;
    
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
    creat_half_edge(nodeVecMine,eleVecMine,edgeVector);
    /*int eleNo,eleDim,eleAtt;
    eleIn>>eleNo>>eleDim>>eleAtt;

    element ele1(eleIn);
    cout<<ele1;
    for(int it=0;it<timestep;++it){

        scatter(it,&);

        connect(it);
    }*/


}



/*void scatter(const int &time_step, vector<mesh_edge> &mesh_edges, vector<mesh_face> &mesh_body)
{
    int no_edge=mesh_edges.size();
    int no_face=mesh_body.size();
    int no_edge_per_face(3);

    vector<double> nodeCurrent;
    nodeCurrent.reserve(no_face);
	memset(&nodeCurrent[0],0,sizeof(double)*no_face);

    vector<double> nodeVoltage;
    nodeVoltage.reserve(no_face);
	memset(&nodeVoltage[0],0,sizeof(double)*no_face);

    for(int kface=0;kface<no_face;++kface){

        nodeCurrent[kface]=2*(mesh_edges[kface*3+0].Vlinki*mesh_edges[kface*3+0].Ylink+mesh_edges[kface*3+1].Vlinki*mesh_edges[kface*3+1].Ylink+mesh_edges[kface*3+2].Vlinki*mesh_edges[kface*3+2].Ylink);
        nodeVoltage[kface]=nodeCurrent[kface]*mesh_body[kface].Z0;

        mesh_edges[kface*3+0].Vlinkr=nodeVoltage[kface]-mesh_edges[kface*3+0].Vlinki;

        mesh_edges[kface*3+1].Vlinkr=nodeVoltage[kface]-mesh_edges[kface*3+1].Vlinki;

        mesh_edges[kface*3+2].Vlinkr=nodeVoltage[kface]-mesh_edges[kface*3+2].Vlinki;

    }

    //-----open circuit voltage and closed circuit current------------------

}

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
