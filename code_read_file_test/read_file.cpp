#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

int main(int argc, char* argv[]){

	cout<<"\nHello World!"<<endl;
	ifstream read_data("ex2.msh");
	
	if(!read_data){cerr<<"\nFile no open!";}
	
	cout<<"\nFile open successfully!";
	
	string input_line;

	while(getline(read_data,input_line)){
        //cout<<"\nInput line: "<<input_line;

        if(input_line=="$MeshFormat\r"){
            cout<<"\nmesh format read"<<endl;

            getline(read_data,input_line);
            double mesh_version, file_type, data_size;
            stringstream mesh_format(input_line);

            mesh_format>>mesh_version>>file_type>>data_size;
            cout<<"\nMesh Format Output \nversion: "<<mesh_version;
            cout<<"   file_type: "<<file_type;
            cout<<"   data_size: "<<data_size; 
            cout<<"\nMesh Format read finish.";
        }
        
        if(input_line=="$Nodes\r"){
            
            getline(read_data,input_line);
            double node_size;
            istringstream read_node(input_line);
            read_node>>node_size;
            cout<<"\nNode Size: "<<node_size;
            
            for(int iNode=0;iNode<node_size;++iNode){
                getline(read_data,input_line);
                double node_x,node_y,node_id;
                stringstream node_coord(input_line);
                node_coord>>node_id>>node_x>>node_y;

                cout<<"\nNode Read -- ID: "<<node_id<<"  "<<node_x<<"  "<<node_y<<endl;
            }

            getline(read_data,input_line);
            if(input_line!="$EndNodes\r"){
                cout<<"\nNode reading not Finish!"<<endl; 
                cout<<input_line<<endl;
            }
            
            cout<<"\nNode Read End."<<endl;
        }

        if(input_line=="$Elements\r"){
            
            double element_size;
            getline(read_data,input_line);
            stringstream read_element(input_line);
            read_element>>element_size;
            cout<<"\nElement Size: " <<element_size;
            
            for(int iEle=0;iEle<element_size;++iEle){
                getline(read_data,input_line);
                cout<<"\nInput line: "<<input_line;
                double ele_id,ele_type,no_tags,physical_entity,elementary_entity,ele_vertex0,ele_vertex1,ele_vertex2;
                istringstream ele_input(input_line);
                ele_input>>ele_id>>ele_type>>no_tags>>physical_entity>>elementary_entity>>ele_vertex0>>ele_vertex1>>ele_vertex2;
                
                cout<<"\nElement Read: "<<ele_id<<"  "<<ele_type<<"  "<<no_tags<<"  "<<physical_entity<<"  ";
                cout<<elementary_entity<<"  "<<ele_vertex0<<"  "<<ele_vertex1<<"  "<<ele_vertex2<<endl;

            }

            getline(read_data,input_line);
            if(input_line!="$EndElements\r"){
                cout<<"\nElement reading not Finish!"<<endl;
                cout<<input_line<<endl;
            }
            
            cout<<"\nElement Read End."<<endl;
        }

		//if(input_line[0]=='$') continue;
		//double data1,data2,data3;
		//stringstream is_data(input_line);
		
		//is_data>>data1>>data2>>data3;
		//cout<<"\nData input: \n"<<data1<<endl<<data2<<endl<<data3<<endl;
	}
    
   
	cout<<"\nPROGRAM ENDS"<<endl;
}
