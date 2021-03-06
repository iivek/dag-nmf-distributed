/* 
 * File:   TRow.cpp
 * Author: ivek
 * 
 * Created on January 17, 2012, 10:29 AM
 */

#include "TRow.h"
#include "THyperRow.h"
#include "SigTRow.h"
#include "SupportXCol.h"
#include "VCol.h"

#include "VertexVisitor.h"

#include <cblas.h>
#include "math.h"
#include <boost/math/special_functions/digamma.hpp>

#include <iostream>
#include <string>
#include <vector>

TRow::TRow(unsigned int row)
{
    this->row = row;
    this->eValues = new element_type[VertexProxy::K->get_val()];
    this->expELogValues = new element_type[VertexProxy::K->get_val()];
}

TRow::TRow(const TRow& orig)
{
    this->eValues = new element_type[VertexProxy::K->get_val()];
    this->expELogValues = new element_type[VertexProxy::K->get_val()];
    this->row = orig.row;
    memcpy( this->eValues, orig.eValues, sizeof(element_type)*VertexProxy::K->get_val() );
    memcpy( this->expELogValues, orig.expELogValues, sizeof(element_type)*VertexProxy::K->get_val() );    
};

TRow::~TRow() {
    delete(this->eValues);
    delete(this->expELogValues);
}

UpdateFunctionDelegator* TRow::clone() {
    UpdateFunctionDelegator* clone = new TRow(*this);
    return (UpdateFunctionDelegator*) clone;
}

void TRow::accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) {
    v.visit(this, scope, schedule);
}

void TRow::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
    gl::edge_list in_edges = scope.in_edge_ids();

    /* The order in which vertices are fetched is important here - as listed in 
     * ModelNMF's initgraph()
     * 
     * (in_edges.size()-2)    VCol objects followed by
     * 1                      THyperRow object     
     * 1                      SigTRow object followed by     
     */
    unsigned int helper = in_edges.size()-2;
    unsigned int offsetVCol = 0;
    unsigned int offsetTHyperRow = offsetVCol+helper;    
    unsigned int offsetSigTRow = offsetTHyperRow + 1;            
    
    element_type alpha[VertexProxy::K->get_val()];
    element_type beta[VertexProxy::K->get_val()];
    
    /* TODO: use SIMD for elementwise operations
     */
    const VertexProxy& neighborTHyperRow = scope.neighbor_vertex_data(
        scope.source(in_edges[offsetTHyperRow]));
    const THyperRow& delegatorTHyperRow = (const THyperRow&) *neighborTHyperRow.getDelegator();
    const VertexProxy& neighborSigTRow = scope.neighbor_vertex_data(
        scope.source(in_edges[offsetSigTRow]));                        
    const SigTRow& delegatorSigTRow = (const SigTRow&) *neighborSigTRow.getDelegator();
    
    /* TODO: use SIMD addition
     */
    for(size_t i = 0; i < VertexProxy::K->get_val(); ++i) {
        alpha[i] = delegatorTHyperRow.aValues[i] + delegatorSigTRow.values[i];
    }
    
    /* TODO: SIMD addition
     */
    // initialize accumulator to zeros
    memset(beta, 0, sizeof(element_type) * VertexProxy::K->get_val());
    for (size_t i = 0; i < helper; ++i) {        
        // VCol:
        const VertexProxy& neighborVCol = scope.neighbor_vertex_data(scope.source(in_edges[offsetVCol+i]));
        const VCol& delegatorVCol = (const VCol&) *neighborVCol.getDelegator();
        
        for(size_t accross = 0; accross < VertexProxy::K->get_val(); ++accross) {
            beta[accross] += delegatorVCol.scaleRelatedValues[accross];
        }        
//        cblas_daxpy(VertexProxy::K->get_val(), 1.0, delegatorVCol.eLogValues, 1, beta, 1);
    }
    
    /* TODO: use SIMD addition and division
     */           
    for(size_t i = 0; i < VertexProxy::K->get_val(); ++i) {
        beta[i] = 1.0/(delegatorTHyperRow.aValues[i]*delegatorTHyperRow.bValues[i] + beta[i]);        
        this->eValues[i] = alpha[i]*beta[i];
        this->expELogValues[i] = exp(boost::math::digamma(alpha[i]))*beta[i];
    }
      
}