#include "utlm.h"

using namespace std;


int main(){

	cout<<"Hello World!"<<endl;
	
	vector<node> read_nodes;
	vector<element> read_faces;
	read_from_gmsh("ex4.msh",read_nodes,read_faces);

	node_vec mNodes(read_nodes);
	//cout<<"\nNodes: "<<mNodes;

	faces mFaces(read_faces);
	//cout<<"\nface size： "<<mFaces.no_elements;
	//cout<<"\nFaces: "<<mFaces<<endl<<endl;

	double times(0.),min_link(0.),min_edge(0.);
	double E_field(0.);
	double delay(1.);

	vector<edge> edgeVector;

	vector<int> boundaryVector;

	create_half_edge(mNodes,mFaces,edgeVector,boundaryVector);

	for(vector<int>::iterator it=boundaryVector.begin();it!=boundaryVector.end();++it){
		cout<<"\nBoundary Edge Index: "<<*it;
	}
	cout<<endl;

	min_edge_link_length(edgeVector,min_edge,min_link);

	double dt(min_link*sqrt(2*constants::get_e0*constants::get_u0())/1.4);

	cout<<"\nTime Step: "<<dt<<endl;

	calAdmittance(dt,mFaces,edgeVector);


}
