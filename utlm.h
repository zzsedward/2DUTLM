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

using namespace std;

struct node{
    private:
        int id;
        double vertex[2];
        int bound_marker;
    public:
        node(const int _id, const double _x, const double _y, const int _bound_marker)
            :id(_id), bound_marker(_bound_marker){
                vertex[0]=_x;
                vertex[1]=_y;
            }
        node(const int _id, const double* const _vertex, const int _bound_marker)
            :id(_id), bound_marker(_bound_marker){
                memcpy(vertex,_vertex,2.0*sizeof(double));
            }

        node(const node &_node)
            :id(_node.id), bound_marker(_node.bound_marker){
                memcpy(vertex,_node.vertex,2*sizeof(double));
            }

        node(const char filename[]){
            ifstream fin(filename);
            fin>>id>>vertex[0]>>vertex[1]>>bound_marker;
        }

        node(ifstream &fin){
            fin>>id>>vertex[0]>>vertex[1]>>bound_marker;
            cout<<"\nID: "<<id;
        }
        const node& operator=(const node& _node){
            if(this==&_node) return (*this);
        }

        friend istream& operator>>(istream &in, node &_node){
            in>>_node.id>>_node.vertex[0]>>_node.vertex[1]>>_node.bound_marker;
            return(in);
        }

        friend ostream& operator<<(ostream &out, const node &_node){
            out<<endl<<_node.id<<"  ";
            out<<_node.vertex[0]<<"  "<<_node.vertex[1]<<"  "<<_node.bound_marker<<endl;

            return(out);
        }

};

struct element{
    private:
        int id;
        int vertex[3];
        int attribute;
    public:
        element(const int _id, const int _v1, const int _v2, const int _v3, const int _attr)
            :id(_id), attribute(_attr){
                vertex[0]=_v1;
                vertex[1]=_v2;
                vertex[2]=_v3;
            }
        element(const int _id, const int* const _vertex, const int _attr)
            :id(_id), attribute(_attr){
                memcpy(vertex,_vertex,3.0*sizeof(int));
            }

        element(const element &_element)
            :id(_element.id), attribute(_element.attribute){
                memcpy(vertex,_element.vertex,3*sizeof(int));
            }

        element(const char filename[]){
            ifstream fin(filename);
            fin>>id>>vertex[0]>>vertex[1]>>vertex[2]>>attribute;
        }

        element(ifstream &fin){
            fin>>id>>vertex[0]>>vertex[1]>>vertex[2]>>attribute;
            cout<<"\nID: "<<id;
        }
        const element& operator=(const element& _element){
            if(this==&_element) return (*this);
        }

        friend istream& operator>>(istream &in, element &_element){
            in>>_element.id>>_element.vertex[0]>>_element.vertex[1];
            in>>_element.vertex[2]>>_element.attribute;
            return(in);
        }

        friend ostream& operator<<(ostream &out, const element &_element){
            out<<endl<<_element.id<<"  ";
            out<<_element.vertex[0]<<"  "<<_element.vertex[1]<<"  "<<_element.vertex[2];
            out<<"  "<<_element.attribute<<endl;

            return(out);
        }
};

struct node_vec{
    private:
        vector<node> nodex;
    public:
        node_vec(const vector<node>& _node_vec){

            int no_nodes(_node_vec.size());

            memcpy(&nodex[0],&_node_vec[0],no_nodes*sizeof(node));
        }

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

        friend ostream& operator<<(ostream& out, node_vec &_node_vec){
            int nodeSize(_node_vec.nodex.size());

            for(int iNode=0;iNode<nodeSize;++iNode){
                out<<_node_vec.nodex[iNode];
            }
            return(out);
        }
        

};

struct ele_vec{
    private:
        vector<element> eleVec;
    public:
        ele_vec(const vector<element>& _eleVec){

            int noEle(_eleVec.size());

            memcpy(&eleVec[0],&_eleVec[0],noEle*sizeof(element));
        }

        ele_vec(const char filename[]){

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

        friend ostream& operator<<(ostream& out, ele_vec &_eleVec){
            int eleSize(_eleVec.eleVec.size());

            for(int iEle=0;iEle<eleSize;++iEle){
                out<<_eleVec.eleVec[iEle];
            }
            return(out);
        }
        
};

#endif