Read Mesh -> set material property; set boundary; calculate admittance ->
-> set excitation -> Run TLM;

Read & load mesh:
	filetype: node, ele;
	node file format:
		first line: number_vertices; dimension; attribute; boundary_marker; 
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

	//create_half_edge function:
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

	//find circumcenter
	vertice a(x1,y1), b(x2,y2), c(x3,y3);
	d=2*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));

	cc.x=((x1*x1+y1*y1)(y2-y3)+(x2*x2+y2*y2)(y3-y1)+(x3*x3+y3*y3)(y1-y2))/d;
	cc.y=((x1*x1+y1*y1)(x3-x2)+(x2*x2+y2*y2)(x1-x3)+(x3*x3+y3*y3)(x2-x1))/d;
	
	//find_boundary_edge_function:
	//find boundary edge when flip edge=-1;
		
		vector<int> boundary_edge_vec;
		for(int iHe=0;iHe<no_half_edges;iHe++){
			if(halfedge_vec[iHe].flip==-1){
				boundary_edge_vec.push_back(iHe);
			}
			
	
	//find edge length and link length
	for(int iHe=0;iHe<no_half_edges;iHe++){
		vertex_index1=halfedge_vec[iHe].vertice[0];
		vertex_index2=halfedge_vec[iHe].vertice[1];
		
		vertex11=mesh_vertice[vertex_index1].vertex[0];
		vertex12=mesh_vertice[vertex_index1].vertex[1];
		vertex21=mesh_vertice[vertex_index2].vertex[0];
		vertex22=mesh_vertice[vertex_index2].vertex[1];
		
		edge_length=sqrt((vertex11-vertex21)^2+(vertex12-vertex22)^2);
		midpoint[0]=(vertex11+vertex21)/2.0;
		midpoint[1]=(vertex12+vertex22)/2.0;

		flip_edge_index=half_edge_vec[iHe].flip;
		CCM_iHe;
		CCM_flip;
		if(flip_edge_index==-1){
			linklength=distance(CCM_iHe,midpoint);}
		else
			{linklength=0.5*distance(CCM_iHe,CCM_flip);

	}    	

	//triangle area
	for (int inF=0;inF<no_faces;++inF){
		edgelength0=half_edge_vec[inF*3+0].edge_length;
		edgelength1=half_edge_vec[inF*3+1].edge_length;
		edgelength2=half_edge_vec[inF*3+2].edge_length;
		
		half_perimeter=(edgelength0+edgelength1+edgelength2)*0.5;
		arear=sqrt(half_perimeter*(half_perimeter-edgelength0)*(half_perimeter-
				edgelength1)*(half_perimeter-edgelength2));
	}

	//Calculate Adimittance----------------------------
	for (int iHe=0;iHe<no_half_edge;++iHe){
		edgeLength;
		linkLength;   
		deltaT;   //timestep
		YLlink=0.5*edgeLength*deltaT/(linkLength*mu);
		YCstub=edgeLength*linkLength*epsilon/deltaT;
		Ystub=YCstub-YLlink;

		//connect flag
		
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