#include "utlm.h"

using namespace std;


int main(){

	cout<<"Hello World!"<<endl;
	
	vector<node> read_nodes;
	vector<element> read_faces;
	read_from_gmsh("mesh1.msh",read_nodes,read_faces);

	node_vec mNodes(read_nodes);
	cout<<"\nNodes: "<<mNodes;

	faces mFaces(read_faces);
	cout<<"\nface size： "<<mFaces.eleVec.size();
	cout<<"\nFaces: "<<mFaces<<endl<<endl;

	double times(0.),min_link(0.),min_edge(0.);
	double E_field(0.);
	int total_time_steps(10);

	vector<edge> edgeVector;

	vector<int> boundaryVector;

	create_half_edge(mNodes,mFaces,edgeVector,boundaryVector);
	EdgeVector edge_vectors(edgeVector);
	//cout<<edge_vectors;

	for(vector<int>::iterator it=boundaryVector.begin();it!=boundaryVector.end();++it){
		cout<<"\nBoundary Edge Index: "<<*it;
	}
	cout<<endl;

	min_edge_link_length(edgeVector,min_edge,min_link);
	cout<<"\nminimum link length: "<<min_link;

	double dt(min_link*sqrt(2.0*constants::get_e0()*constants::get_u0())/1.4);

	cout<<"\nTime Step: "<<dt<<endl;

	calAdmittance(dt,mFaces,edgeVector);

	//set material properties - cylindrical resonator - everywhere air
	vector<double> mu_r,epsilon_r;
	mu_r.push_back(1.0);
	epsilon_r.push_back(1.0);

	set_material_property(mFaces,epsilon_r,mu_r);

	//create mesh body index list  -  remove the filp edge index 
	list<int> mesh_body_list;
	create_mesh_body_vector(edgeVector,mesh_body_list);
	for(list<int>::iterator iii=mesh_body_list.begin();iii!=mesh_body_list.end();++iii){
		cout<<"\nMesh Body id: "<<*iii;
	}

	//Set reflection coefficients for boundary edges--------------------
	//---in this case short circuit r=-1-----------------------
	vector<double> bound_condition;
	bound_condition.resize(boundaryVector.size());
	memset(&bound_condition[0],-1.0,sizeof(double)*boundaryVector.size());
						
	vector<double> reflection_coeff,Y_boundary;

	create_reflection_coeff(edgeVector,boundaryVector,mFaces,reflection_coeff,Y_boundary,bound_condition);	

	//Set excitation wave ---------------------------------------
	double width(25*min_edge/constants::get_c0());
	double delay(1.);

	gause_wave gause_excite(width,delay);
	cout<<"\nGause speed: "<<gause_excite.wave_speed;
	cout<<"\nWave direction: "<<gause_excite.direction[0]<<"  "<<gause_excite.direction[1];

	vector<double> Efield_vector_time;
	vector<double> time_array;
	for(int itime=0;itime<total_time_steps;++itime){ 	//itime - iteration index for integer total time steps
		time_array.push_back((double) itime*dt);
	}

	gause_excite.evaluate(time_array,Efield_vector_time);

	cout<<"\nExcitation field magnitude: "<<endl;
	for(int istep=0;istep<total_time_steps;++istep){
		cout<<Efield_vector_time[istep]<<endl;
	}
}
