#ifndef UTLM_H
#define UTLM_H
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
#include <list>
#include <map>

using namespace std;

struct constants{
    static const double get_pi(){return(M_PI);}
    static const double get_e0(){return(8.8541878176e-12);}
    static const double get_u0(){return(4.0e-7*get_pi());}
    static const double get_c0(){return(1./sqrt(get_e0()*get_u0()));}
};
struct node{

    int node_id;
    double node_vet[2];
    int node_bound_marker;

    node(const int _id=0, const double _x=0.0, const double _y=0.0, const int _bound_marker=0)
        :node_id(_id), node_bound_marker(_bound_marker){
            node_vet[0]=_x;
            node_vet[1]=_y;
        }

    node(const int _id, const double* const _vertex, const int _bound_marker)
        :node_id(_id), node_bound_marker(_bound_marker){
            memcpy(node_vet,_vertex,2.0*sizeof(double));
        }

    node(const node &_node)
        :node_id(_node.node_id), node_bound_marker(_node.node_bound_marker){
            memcpy(node_vet,_node.node_vet,2*sizeof(double));
        }

    node(const char filename[]){
        ifstream fin(filename);
        fin>>node_id>>node_vet[0]>>node_vet[1]>>node_bound_marker;
    }

    node(ifstream &fin){
        fin>>node_id>>node_vet[0]>>node_vet[1]>>node_bound_marker;
        cout<<"\nID: "<<node_id;
    }
    const node& operator=(const node& _node){
        if(this==&_node) return (*this);
    }

    friend istream& operator>>(istream &in, node &_node){
        in>>_node.node_id>>_node.node_vet[0]>>_node.node_vet[1]>>_node.node_bound_marker;
        return(in);
    }

    friend ostream& operator<<(ostream &out, const node &_node){
        out<<endl<<_node.node_id<<"  ";
        out<<_node.node_vet[0]<<"  "<<_node.node_vet[1]<<"  "<<_node.node_bound_marker<<endl;

        return(out);
    }

};

struct element{
    
    int ele_id;
    int ele_vet[3];
    int ele_attri;
    int fnum;
    double epr;
    double mur;
    double Z0;

    element(const int _id=0, const int _v1=0, 
            const int _v2=0, const int _v3=0, 
            const int _attr=0, int _fnum=1,
            double _epr=1.0, double _mur=1.0,
            double _Z0=0.0)
        :ele_id(_id), ele_attri(_attr),
         fnum(_fnum), epr(_epr), mur(_mur),Z0(_Z0){
            ele_vet[0]=_v1;
            ele_vet[1]=_v2;
            ele_vet[2]=_v3;
    }

    element(const int _id=0, const int* const _vertex=NULL, 
            const int _attr=0, int _fnum=1,
            double _epr=1.0, double _mur=1.0,
            double _Z0=0.0)
        :ele_id(_id), ele_attri(_attr), 
         fnum(_fnum), epr(_epr), mur(_mur),Z0(_Z0){
           
            memcpy(ele_vet,_vertex,3.0*sizeof(int));
    }

    element(const element &_element)
        :ele_id(_element.ele_id), ele_attri(_element.ele_attri), 
         fnum(_element.fnum),epr(_element.epr),mur(_element.mur),
         Z0(_element.Z0){
            
            memcpy(ele_vet,_element.ele_vet,3*sizeof(int));
    }

    //--Set face number -----------------------------------
    void set_face_number(const int face_number){
        fnum=face_number;
    }

    int get_fnum() const {
        return(fnum);
    }

    //--Set material properties----------------------------
    void set_material(
        const double &_epr,
        const double &_mur){
            epr=_epr;
            mur=_mur;
    }

    //--Get relative epsilon -----------------------
    double get_epsilonr() const {
        return(epr);
    }

    //--Set Z0 ------------------------------------
    void set_Z0(const double &_Z0){
        Z0=_Z0;
    }

    double get_Z0() const{
        return(Z0);
    }


    //--Input constructors---------------------------------
    element(const char filename[]){
        ifstream fin(filename);
        fin>>ele_id>>ele_vet[0]>>ele_vet[1]>>ele_vet[2]>>ele_attri;
    }

    element(ifstream &fin){
        fin>>ele_id>>ele_vet[0]>>ele_vet[1]>>ele_vet[2]>>ele_attri;
       
        cout<<"\nID: "<<ele_id;
    }

    const element& operator=(const element& _element){
        if(this==&_element) return (*this);
    }

    friend istream& operator>>(istream &in, element &_element){
        in>>_element.ele_id>>_element.ele_vet[0]>>_element.ele_vet[1];
        in>>_element.ele_vet[2]>>_element.ele_attri;

        return(in);
    }

    friend ostream& operator<<(ostream &out, const element &_element){
        out<<endl<<_element.ele_id<<"  ";
        out<<_element.ele_vet[0]<<"  "<<_element.ele_vet[1]<<"  "<<_element.ele_vet[2];
        out<<"  "<<_element.ele_attri;
        out<<"  "<<_element.fnum;
        out<<"  "<<_element.epr<<"  "<<_element.mur<<endl;

        return(out);
    }
};

//---------Struct for node vectors ----------------
struct node_vec{
    
    vector<node> nodex;
    int no_nodes;

    node_vec(const vector<node>& _node_vec){

        no_nodes=_node_vec.size();
	    nodex.reserve(no_nodes);

        memcpy(&nodex[0],&_node_vec[0],no_nodes*sizeof(node));
    }

    //--Input constructors----------------------------
    node_vec(const char filename[]){

        ifstream fin(filename);
        int nodeSize,nodeDim,nodeAttri,nodeBound,nodeId;
        double x,y;
        fin>>nodeSize>>nodeDim>>nodeAttri>>nodeBound;

        for(int iNode=0;iNode<nodeSize;++iNode){
            fin>>nodeId>>x>>y>>nodeBound;
            node nodeTemp(nodeId,x,y,nodeBound);
            nodex.push_back(nodeTemp);
        }
    }

    //--Output operator---------------------------------
    friend ostream& operator<<(ostream& out, const node_vec &_node_vec){
        int nodeSize(_node_vec.no_nodes);

        for(int iNode=0;iNode<nodeSize;++iNode){
            out<<_node_vec.nodex[iNode];
        }
        return(out);
    }
        

};

//--------- Struct for element vectors as faces-------
struct faces{
    
    vector<element> eleVec;
    int no_elements;
    int no_material;

    //--Constructor------------------------------------
    faces(const vector<element>& _eleVec, int _no_material=0):
        no_material(_no_material){

        no_elements=_eleVec.size();
        eleVec.reserve(no_elements);

        memcpy(&eleVec[0],&_eleVec[0],no_elements*sizeof(element));
    }

    //--Input constructor from ele file----------------------------
    faces(const char filename[]){

        ifstream fin(filename);
        int eleSize,eleDim,eleAttri,eleId;
        int v1,v2,v3;
        fin>>eleSize>>eleDim>>eleAttri;

            for(int iEle=0;iEle<eleSize;++iEle){
            fin>>eleId>>v1>>v2>>v3>>eleAttri;
            element eleTemp(eleId,v1,v2,v3,eleAttri);
            eleVec.push_back(eleTemp);
        }
    }

    //--Output Operator------------------------------------
    friend ostream& operator<<(ostream& out, const faces &_eleVec){

        int eleSize(_eleVec.no_elements);
        //cout<<"\nFaces size: "<<eleSize;
        for(int iEle=0;iEle<eleSize;++iEle){
            out<<_eleVec.eleVec[iEle];
        }
        return(out);
    }

    //--Get vector size------------------------------
    int get_no_face() const{
        return (eleVec.size());
    }

    //--Set Z0 ---------------------
    void set_impedance(const int &face_id,
                       const int &_Z0){
        const int no_face(eleVec.size());
        if(face_id>=0 && face_id<no_face){
            eleVec[face_id].set_Z0(_Z0);
        }
        else{cout<<"\nInvalid face id!"<<endl;}
    }

    double get_Z0_with_id(const int &face_id) const{
        return(eleVec[face_id].get_Z0());
    }

    //--Set face numbers-----------------------------
    void faces_face_number(const vector<int> face_id,
                            const int face_number){
         
         const int id_size(face_id.size());
         for(int i=0;i<id_size;++i){
             int faceId(face_id[i]);
             eleVec[faceId].set_face_number(face_number);
         }
    }

    //--Get face number with id------------------
    int get_fnum_with_id(const int &face_id) const{
        return(eleVec[face_id].get_fnum());
    }

    //--get no of material-----------------------------
    int find_no_material() const{
        const int no_faces(eleVec.size());
        int max_fnum(1);
        for(int i=0;i<no_faces;++i){
            
            const int fnum_face(eleVec[i].get_fnum());
            if(max_fnum<fnum_face){
                max_fnum=fnum_face;
            }
        }

        return(max_fnum);
    }

    //--get relative epsilon-----------------------
    double get_epsilonr_with_id(const int& face_id) const{
        return(eleVec[face_id].get_epsilonr());
    }

};


struct edge{

    int edgeVet[2];
    int faceId;
    int edgeFlip;
    double ccm[2];
    double midpoint[2];
    double edgeLength;
    double linkLength;
    double Vlinki;
    double Vlinkr;
    double Vstub;
    double Ylink;
    double Ystub;
    
    edge():faceId(0),edgeFlip(-1),
                edgeLength(0),linkLength(0),
                Vlinki(0),Vlinkr(0),Vstub(0),
                Ylink(0),Ystub(0){

            memset(&edgeVet[0],0,2*sizeof(int));
            memset(&ccm[0],0,2*sizeof(double));
            memset(&midpoint[0],0,2*sizeof(double));
    }

    edge(int *_edgeVet, 
            int &_faceId, 
            int &_edgeFlip,
            double *_ccm,
            double *_midpoint,
            double _edgeLength,
            double _linkLength,
            double _Vlinki,
            double _Vlinkr,
            double _Vstub,
            double _Ylink,
            double _Ystub)
            :faceId(_faceId),edgeFlip(_edgeFlip),
                edgeLength(_edgeLength),linkLength(_linkLength),
                Vlinki(_Vlinki),Vlinkr(_Vlinkr),Vstub(_Vstub),
                Ylink(_Ylink),Ystub(_Ystub){

        memcpy(&edgeVet[0],&_edgeVet[0],2*sizeof(int));
        memcpy(&ccm[0],&_ccm[0],2*sizeof(double));
        memcpy(&midpoint[0],&_midpoint[0],2*sizeof(double));
            
    }

    ~edge(){}

    friend ostream& operator<<(ostream &out, const edge &_edge){
        out<<endl<<_edge.edgeVet[0]<<"  "<<_edge.edgeVet[1];
        
        out<<" face ID: "<<_edge.faceId;

        return(out);
    }

    int get_face_id() const{
        return(faceId);
    }

    int get_flip_edge() const{
        return(edgeFlip);
    }

    double get_edge_length() const{
        return(edgeLength);
    }

    double* get_vertice() const{
        double *_vertice=new double[2];
        _vertice[0]=edgeVet[0];
        _vertice[1]=edgeVet[1];

        return(_vertice);
    }

    double* get_ccm() const{
        static double _ccm[2];
        _ccm[0]=ccm[0];
        _ccm[1]=ccm[1];
        return(_ccm);
    }

    double calc_Vlinki_times_Ylink() const{
        return(Vlinki*Ylink);
    }

    double get_Vlinki() const{
        return(Vlinki);
    }

    void set_Vlinki(const double &voltage){
        Vlinki=voltage;
    }

    void set_Vlinkr(const double &voltage){
        Vlinkr=voltage;
    }

    void set_Vstub(const double &voltage){
        Vstub=voltage;
    }

    double get_Vlinkr() const{
        return(Vlinkr);
    }

    double get_Vstub() const{
        return(Vstub);
    }

    double get_Ylink() const{
        return(Ylink);
    }

    double get_Ystub() const{
        return(Ystub);
    }

};

struct EdgeVector{

    vector<edge> edge_vector;
    
    //EdgeVector(){cout<<"Empty edge vector"<<endl<<edge_vector.size()<<endl;}

    EdgeVector(const vector<edge> &_edge_vector){

        cout<<"\nedge vector create."<<endl;
        const int no_edges(_edge_vector.size());
        //cout<<"\nnumber of edges: "<<no_edges;
        edge_vector.resize(no_edges);
        
        memcpy(&edge_vector[0],&_edge_vector[0],no_edges*sizeof(edge));
    }

    friend ostream& operator<<(ostream& out, const EdgeVector &_edge_vector){

        int no_edges(_edge_vector.edge_vector.size());
        out<<"\nEdge vector size: "<<no_edges<<endl;

        for(int iedge=0;iedge<no_edges;++iedge){
            out<<_edge_vector.edge_vector[iedge];
        }

        out<<endl;

        return(out);
    }


};

struct gause_wave{
    double width;
    double delay;
    double wave_speed;
    double amplitude;
    double direction[2];

    gause_wave(double _width=0.0,double _delay=1.0,double _wave_speed=constants::get_c0(),double _amplitude=1.0)
            :width(_width),delay(_delay),wave_speed(_wave_speed),amplitude(_amplitude){
                direction[0]=1;
                direction[1]=0;
            }
    
    void evaluate(const vector<double> &time_array, vector<double> &Vsource_array);
};


void read_from_gmsh(const char filename[],
		    vector<node> &node_input,
		    vector<element> &face_input);

void create_half_edge(const node_vec &mnode,
                     const faces &mface,
                     vector<edge> &edge_vec,
                     vector<int> &mesh_boundary);

void min_edge_link_length(const vector<edge> &edge_vec,
                        double &_minEdge,
                        double &_minLink);

void gaussian_wave_excite(const double &width, 
                        const double &delay, 
                        const double &t, 
                        double &Efield);

void calAdmittance(const double &dt,
                    faces &mface,
                    vector<edge> &edge_vector);

void set_inner_circle_different_face_number(
     faces &my_face,
     const vector<edge> &my_edges,
     const double &inner_circle_radius,
     const double inner_circle_centre[],
     const int &face_number);

void set_material_property(
    faces &my_face,
    const vector<double> &_epr,
    const vector<double> &_mur);

void create_PEC_bound_vector(
    const vector<edge> &my_edges,
    const faces &my_faces,
    vector<int> &pec_bound);

void create_mesh_body_vector(
    const vector<edge> &my_edges,
    list<int> &mesh_body);

void create_mat_bound_vector(
    const vector<edge> &my_edge);

void create_reflection_coeff(
    const vector<edge> &my_edges,
    const vector<int> &my_bound_edges,
    const faces &my_faces,
    vector<double> &my_refl_coeff,
    vector<double> &Y_boundary,
    const vector<double> &my_condition);

void scatter(const int &time_step,
             vector<edge> &edge_vector,
             const faces &my_faces,
             const vector<int> &mesh_boundary);

void connect(const int &time_step,
             vector<edge> &my_edges,
             const faces &my_faces,
             const vector<int> &my_bound_edges,
             const vector<int> &my_pec_bounds,
             const vector<int> &my_inner_edges,
             const vector<double> &my_refl_coeff,
             const vector<double> &my_y_boundary);

void edge_excite(vector<edge> &mesh_edges, 
                const vector<int> &source_edge,
                const double &Vsource);

void run(vector<edge> &edge_vector,
         const faces &face_vector,
         const vector<int> &mesh_boundary_edge_id,
         const vector<int> &mesh_pec_boundary,
         const vector<int> &mesh_body,
         const vector<int> &mesh_excite_edge_id,
         const vector<double> &mesh_reflection_coefficient,
         const vector<double> &mesh_Y_boundary,
         const vector<double> &excite_field_td,
         const int &number_time_step); 
        


#endif
