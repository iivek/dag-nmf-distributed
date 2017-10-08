/* 
 * File:   SigTRow.cpp
 * Author: ivek
 * 
 * Created on January 19, 2012, 11:27 AM
 */

#include "SigTRow.h"
#include "TmpCol.h"
#include "VCol.h"
#include "TRow.h"
#include "SupportXCol.h"
#include "VertexVisitor.h"

//#include <cblas.h>
#include <cblas.h>

#include <iostream>
#include <string>
#include <vector>

using namespace std;

SigTRow::SigTRow()
//        : values(VertexProxy::K->get_val())
{
//    std::cout << "in SigTRow constructor" << std::endl;
    this->values = new element_type[VertexProxy::K->get_val()];
    memset(this->values, 0, sizeof (element_type) * VertexProxy::K->get_val());
}

SigTRow::SigTRow(const SigTRow& orig)
//        : values(VertexProxy::K->get_val())
{
//    std::cout << "in SigTRow copy constructor" << std::endl;      
//    std::copy(orig.values.begin(), orig.values.end(), this->values.begin());
    this->values = new element_type[VertexProxy::K->get_val()];
    memcpy(this->values, orig.values, sizeof (element_type) * VertexProxy::K->get_val());
}

SigTRow::~SigTRow() {
    delete(this->values);
}

UpdateFunctionDelegator* SigTRow::clone() {
    UpdateFunctionDelegator* clone = new SigTRow(*this);
    return (UpdateFunctionDelegator*) clone;

}

void SigTRow::accept(VertexVisitor &v, gl::iscope& scope, gl::icallback& schedule)  {
    v.visit(this, scope, schedule);
}    

void SigTRow::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
    std::cout << "SigTRow update invoked, " << scope.color() <<std::endl;

    /* Get a reference to the vertex data and visible edges */
    //    vertex_data& us = scope.vertex_data();    

    gl::edge_list in_edges = scope.in_edge_ids();

    /* init accumulator to allzeros */
    memset(&this->values[0], 0, sizeof(element_type)*VertexProxy::K->get_val());    
    // Iterate through neighboring vertices
    /*
     * The order in which vertices are fetched is important here - as in
     * ModelNMF's initgraph()
     *           
     * (in_edges.size()-1)/3  VCol objects followed by
     * a TRow object followed by     
     * (in_edges.size()-1)/3 TmpCol objects followed by
     * (in_edges.size()-1)/3 SupportXCol objects followed by
     */
    unsigned int helper = (in_edges.size()-1)/3;
    unsigned int vColOffset = 0;
    unsigned int tRowOffset = vColOffset + helper;    
    unsigned int tmpColOffset = tRowOffset + 1 ;
    unsigned int supportXColOffset = tmpColOffset + helper;

    const element_type* buffer;
    /* Fetching triplets of SupportXCol, LvCols and TmpCols
     * TODO: use SIMD for elementwise multiplication
     */    
    unsigned int numPasses = (in_edges.size() - 1)/3;
    
    const VertexProxy& tneighbor = scope.neighbor_vertex_data(scope.source(in_edges[ tRowOffset ]));
    const TRow& tdelegator = (const TRow&) *tneighbor.getDelegator();
        
    for (size_t i = 0; i < numPasses; ++i) {
        // SupportXCol:
        const VertexProxy& neighbor = scope.neighbor_vertex_data(scope.source(in_edges[ i+supportXColOffset ]));
        const SupportXCol& delegator = (const SupportXCol&) *neighbor.getDelegator();
        unsigned int location = delegator.binarySearch(tdelegator.row);
        // VCol:       
        const VertexProxy& neighbor1 = scope.neighbor_vertex_data(scope.source(in_edges[ i+vColOffset ]));
        const VCol& delegator1 = (const VCol&) *neighbor1.getDelegator();
        buffer = &delegator1.expELogValues[0];
        
        // TmpCol:
        const VertexProxy& neighbor3 = scope.neighbor_vertex_data(scope.source(in_edges[ i+tmpColOffset ]));
        const TmpCol& delegator3 = (const TmpCol&) *neighbor3.getDelegator();
        const element_type& multiplier = delegator3.values[location];                               
        
        cblas_daxpy(VertexProxy::K->get_val(), multiplier, buffer, 1, &this->values[0], 1);
        //std::cout<<this->values[0]<<" "<<this->values[1]<<" "<<this->values[2] << " "<<std::endl;
        //std::cout<<buffer[0]<<" "<<buffer[1]<<" "<<buffer[2] << " "<<std::endl;
        
    }

    /* Fetching TRow
     * TODO: use SIMD for elementwise multiplication
     */
    buffer = &tdelegator.expELogValues[0];
    for (unsigned int col = 0; col < VertexProxy::K->get_val(); ++col) {
        this->values[col] *= buffer[col];
    } 
    //   scheduler.add_task(gl::update_task(scope.vertex(), update_function), 1.0);
}

