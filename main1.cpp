#include "utlm.h"

using namespace std;


int main(){

	cout<<"Hello World!"<<endl;
	node xx;
	node yy(1,1.1,1);
	
	vector<node> read_nodes;
	vector<element> read_faces;
	read_from_gmsh("ex4.msh",read_nodes,read_faces);

	node_vec mNodes(read_nodes);
	cout<<"\nNodes: "<<mNodes;

	faces mFaces(read_faces);
	//cout<<"\nface size： "<<mFaces.no_elements;
	cout<<"\nFaces: "<<mFaces<<endl<<endl;
	
	return 0;
}
