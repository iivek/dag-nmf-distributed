/* 
 * File:   ScaleHyperCol.cpp
 * Author: ivek
 * 
 * Created on February 16, 2012, 10:35 AM
 */

#include "ScaleHyperCol.h"
#include "VCol.h"
//#include "MessageCollector.h"
#include "VertexVisitor.h"

#include <iostream>
#include <string>
#include <vector>

ScaleHyperCol::ScaleHyperCol() {
//    std::cout << "in ShapeHyperCol constructor" << std::endl;
}

ScaleHyperCol::ScaleHyperCol(const ScaleHyperCol& orig) : GenericScale(orig) {
//    std::cout << "in ScaleHyperCol copy constructor" << std::endl;    
}

ScaleHyperCol::~ScaleHyperCol() {
//    std::cout << "ScaleHyperCol destructor" << std::endl;
}

UpdateFunctionDelegator* ScaleHyperCol::clone() {
    UpdateFunctionDelegator* clone = new ScaleHyperCol(*this);
    return (UpdateFunctionDelegator*) clone;
}

void ScaleHyperCol::accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) {
    //        std::cout << "ScaleHyperCol accepted" << std::endl;
    v.visit(this, scope, schedule);
}

void ScaleHyperCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
//    std::cout << "ScaleHyperCol update invoked, " << scope.color() <<std::endl;
    /* TODO: use SIMD for elementwise operations
     */
/*        
    gl::edge_list in_edges = scope.in_edge_ids();
    
    
    const VertexProxy& neighborWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ 0 ]));
    const VCol* const neighbor = (const VCol* const) neighborWrapper.getDelegator();
    
    for(size_t col=0; col < VertexProxy::K->get_val(); ++col) {
        scaleRelatedValues[col] = neighbor->scaleRelatedValues[col];
    } 
 */ 
}