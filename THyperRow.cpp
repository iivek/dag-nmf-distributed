/* 
 * File:   THyperRow.cpp
 * Author: ivek
 * 
 * Created on January 21, 2012, 5:27 PM
 */

#include "THyperRow.h"
#include "VertexVisitor.h"
#include "TRow.h"

#include <iostream>
#include <string>
#include <vector>

THyperRow::THyperRow()
//    :aValues(VertexProxy::K->get_val())
//    ,bValues(VertexProxy::K->get_val())
{
//    std::cout << "in THyperRow constructor" << std::endl;
    this->aValues = new element_type[VertexProxy::K->get_val()];
    this->bValues = new element_type[VertexProxy::K->get_val()];
}

THyperRow::THyperRow(const THyperRow& orig)
//    :aValues(VertexProxy::K->get_val())
//    ,bValues(VertexProxy::K->get_val())
{
    
    //    std::cout << "in THyperRow copy constructor" << std::endl;
    // copying values member
//    std::copy(orig.aValues.begin(), orig.aValues.end(), this->aValues.begin());
//    std::copy(orig.bValues.begin(), orig.bValues.end(), this->bValues.begin());
    this->aValues = new element_type[VertexProxy::K->get_val()];
    this->bValues = new element_type[VertexProxy::K->get_val()];
    memcpy(this->aValues, orig.aValues, sizeof(element_type)*VertexProxy::K->get_val());
    memcpy(this->bValues, orig.bValues, sizeof(element_type)*VertexProxy::K->get_val());
};

THyperRow::~THyperRow() {
//    std::cout << "THyperRow destructor" << std::endl;
    delete(this->aValues);
    delete(this->bValues);
}

UpdateFunctionDelegator* THyperRow::clone() {
    UpdateFunctionDelegator* clone = new THyperRow(*this);
    return (UpdateFunctionDelegator*) clone;
}

void THyperRow::accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) {
    //        std::cout << "THyperRow accepted" << std::endl;
    v.visit(this, scope, schedule);
}

void THyperRow::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
//    std::cout << "THyperRow update invoked, " << scope.color() <<std::endl;
    
    /* TODO: use SIMD for elementwise operations
     */
    
/*    gl::edge_list in_edges = scope.in_edge_ids();
    
    const VertexProxy& neighborWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ 0 ]));
    const TRow* const neighbor = (const TRow* const) neighborWrapper.getDelegator();
    
    for(size_t col=0; col < VertexProxy::K->get_val(); ++col) {
        bValues[col] = neighbor->eValues[col];
    }
*/  
}

