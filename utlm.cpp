#include <string.h>
#include <complex>
#include <math.h>
#include <iterator>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <time.h>
#include <functional>
#include <vector>

using namespace std;
class mesh_face{

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


};


int main(int argc, char* argv[])
{
    cout<<"Hello World!"<<endl;
    int timestep(1e4);
    double del_t(1e-10);

    for(int it=0;it<timestep;++it){

        scatter(it,&);

        connect(it);
    }


}

void scatter(const int &time_step, vector<mesh_edge> &mesh_edges, vector<mesh_face> &mesh_body)
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
}
