#include "utlm.h"

using namespace std;


int main(){

	cout<<"Hello World!"<<endl;
	node xx;
	node yy(1,1.1,1);

	cout<<xx<<endl;
	cout<<yy<<endl;
	
	vector<node> read_nodes;
	vector<element> read_faces;
	read_from_gmsh("ex4.msh",read_nodes,read_faces);
<<<<<<< HEAD
=======
	
	for(vector<node>::iterator iNode=read_nodes.begin();iNode!=read_nodes.end();++iNode){
		cout<<"\n"<<*iNode;}
	cout<<"\nEnd"<<endl;

	node_vec mNodes(read_nodes);
	cout<<"\nNodes: "<<endl<<mNodes;
>>>>>>> master
	return 0;

}
