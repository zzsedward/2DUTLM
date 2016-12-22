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

    node(const int _id, const double _x, const double _y, const int _bound_marker)
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

    element(const int _id=0, const int _v1=0, 
            const int _v2=0, const int _v3=0, 
            const int _attr=0, int _fnum=1,
            double _epr=1.0, double _mur=1.0)
        :ele_id(_id), ele_attri(_attr),
         fnum(_fnum), epr(_epr), mur(_mur){
            ele_vet[0]=_v1;
            ele_vet[1]=_v2;
            ele_vet[2]=_v3;
    }

    element(const int _id=0, const int* const _vertex=NULL, 
            const int _attr=0, int _fnum=1,
            double _epr=1.0, double _mur=1.0)
        :ele_id(_id), ele_attri(_attr), 
         fnum(_fnum), epr(_epr), mur(_mur){
           
            memcpy(ele_vet,_vertex,3.0*sizeof(int));
    }

    element(const element &_element)
        :ele_id(_element.ele_id), ele_attri(_element.ele_attri), 
         fnum(_element.fnum),epr(_element.epr),mur(_element.mur){
            
            memcpy(ele_vet,_element.ele_vet,3*sizeof(int));
    }

//--Set face number -----------------------------------
    void element_face_number(const int face_number){
        fnum=face_number;
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

    node_vec(const vector<node>& _node_vec){

        int no_nodes(_node_vec.size());

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
    friend ostream& operator<<(ostream& out, node_vec &_node_vec){
        int nodeSize(_node_vec.nodex.size());

        for(int iNode=0;iNode<nodeSize;++iNode){
            out<<_node_vec.nodex[iNode];
        }
        return(out);
    }
        

};

//--------- Struct for element vectors as faces-------
struct faces{
    
    vector<element> eleVec;

//--Constructor------------------------------------
    faces(const vector<element>& _eleVec){

        int noEle(_eleVec.size());

        memcpy(&eleVec[0],&_eleVec[0],noEle*sizeof(element));
    }


//--Input constructor----------------------------
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
    friend ostream& operator<<(ostream& out, faces &_eleVec){
        int eleSize(_eleVec.eleVec.size());

        for(int iEle=0;iEle<eleSize;++iEle){
            out<<_eleVec.eleVec[iEle];
        }
        return(out);
    }

//--Set face numbers-----------------------------
    void faces_face_number(const vector<int> face_id,
                            const int face_number){
         
         const int id_size(face_id.size());
         for(int i=0;i<id_size;++i){
             int faceId(face_id[i]);
             eleVec[faceId].element_face_number(face_number);
         }
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

};

void creat_half_edge(const node_vec &mnode,
                     const faces &mface,
                     vector<edge> &edge_vec);

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

void scatter(const int &time_step,
             vector<edge> &edge_vector);

void set_inner_circle_different_face_number(
     const faces &my_face,
     const double &inner_circle_radius,
     const double &inner_circle_centre[2],
     const int &face_number);


#endif