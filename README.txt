Read Mesh -> set material property; set boundary; calculate admittance ->
-> set excitation -> Run TLM;

Read & load mesh:
	filetype: node, ele;
	node file format:
		vertices, boundary marker;
	element file format:
		vertex index for each face;

Convert to halfedge mesh:
	parameters:
		vertices (mesh point vertices coordinate);
		material boundaries;
		mesh_boundary (boundary edge indices generated);
		mesh_body (inner body edge indices generated);
		face(vertice indices for each face, face number, area);
		halfedge(face indices, vertices indices, flip edge 
				indices, circumcenter coordinate, midpoint coordinate,
				edgelength,linklength, Vlinki, Vlinkr, Vstub);
		
		no_faces;
		nVpf(number of vertice per face);
		triangulation;
		create_half_edge_function;
		setLength_function;
		find_boundary_edge_function;	

	create_half_edge function:
		for(int inF=0;inF<nF;++inF){
			for(int i_nVpf=0;i_nVpf<nVpf;++i_nVpf)
			{
				v2=(i_nVpf+1)%nVpf; //needs validate
				start_vertex_index=face[inF].vertice[i_nVpf];
				end_vertex_index=face[inF].vertice[v2];

				half_edge_index=(inF-1)*nVpf+i_nVpf;
				halfedge_vec[half_edge_index].vertice[0]=start_vertex_index;
				halfedge_vec[half_edge_index].vertice[1]=end_vertex_index;

				for(int iFe=0;iFe<half_edge_index;++iFe){
					flip_edge_vertice[0]=end_vertex_index;
					flip_edge_vertice[1]=start_vertex_index;
					
					if(halfedge_vec[iFe].vertice[0]=flip_edge_vertice[1]
						&&halfedge_vec[iFe].vertice[1]=flip_edge_vertice[0]){
							halfedge_vec[half_edge_index].flip=iFe;
							halfedge_vec[iFe].flip=half_edge_index;
						}
				}
			}
		}





vector<mesh> mesh_obj;

int main()
{
	vector<mesh> mesh1;
	mesh1.reserve(number_of_cells);

	for (int i_ts=0;i_ts<timesteps;++i_ts){
		
		scatter;
		
		connect;

		edge_Excite;
	
	}
}

class mesh{

	double node_centre[2];
	double edge[3];
	double neig[3];
	double edge_length[3];
	double link_length[3];
	double Vlinki[3];
	double Vstubr[3];
	double Impedance;
	double epsilonr;
	double mur;
	double area;
	
	mesh(int &_numbe_edge, vector<int> _vertice_index,)
	int number_edges;
	vector<int> vertice_index; 

	
}

 
	
void scatter()
{
	int no_edge(sizeof(vector(mesh_edge));
	int no_edge_per_face(3);
	
	for(int i_edge=0;i_edge<no_edge_per_face;++i_edge)
	{
		Io=2*Vlinki[n]*Ylink[n];
	}
	
	Vo=Io*Z;
	 
	
}

void connect()
{
	
	int index_edge, index_neigh;
	
	if(Y_stub>1e10)
	{	
		E_i[index_edge].Vlinki=0;

		H_i[index_edge].Vlinki=0;

		H_x[index_edge].Vlinki=0;

		H_y[index_edge].Vlinki=0;

		Vlinki[index_edge]=0;

		Vstub[index_edge]=0;

		Vlinki[index_flip]=0;

		Vstub[index_flip]=0;
}

//Excite source edge with Vsource at specific time step; add source effect
//Calculate total admittance at the edge and add up linkline and stub voltages
void edge_excite(vector<mesh_edge> &mesh_edges, vector<int> &source_edge, double Vsource)
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