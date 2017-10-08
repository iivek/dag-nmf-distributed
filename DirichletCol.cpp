/* 
 * File:   DirichletCol.cpp
 * Author: ivek
 * 
 * Created on January 31, 2012, 5:50 PM
 */

#include "DirichletCol.h"
#include "DiscreteCol.h"
#include "DirichletHyperCol.h"
#include "VertexVisitor.h"

#include "math.h"
#include <boost/math/special_functions/digamma.hpp>

DirichletCol::DirichletCol(unsigned int order)
{
    this->order = order;
    this->values = new element_type[VertexProxy::K->get_val() * this->order];
}

DirichletCol::DirichletCol(const DirichletCol& orig)
{
    this->order = orig.order;
    this->values = new element_type[VertexProxy::K->get_val() * this->order];
    memcpy(this->values, orig.values, sizeof (element_type) * VertexProxy::K->get_val() * this->order);    
};

DirichletCol::~DirichletCol() {
    delete(values);
}

UpdateFunctionDelegator* DirichletCol::clone() {
    UpdateFunctionDelegator* clone = new DirichletCol(*this);
    return (UpdateFunctionDelegator*) clone;
}

void DirichletCol::accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) {
    v.visit(this, scope, schedule);
}

void DirichletCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
    /* Neighbors:
     * 1. Discrete
     * 2. hyperparameter
     */
    /* Accumulators for messages */
    element_type naturalParameters[getOrder()][VertexProxy::K->get_val()];
    //    memset(naturalParameters, 0, sizeof(element_type) *  mixtureSize * VertexProxy::K->get_val() );
    gl::edge_list in_edges = scope.in_edge_ids();    
    const VertexProxy& discreteWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ 0 ]));
    const DiscreteCol* const discrete = (const DiscreteCol* const) discreteWrapper.getDelegator();

    const VertexProxy& dirichletHyperWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ 1 ]));
    const DirichletHyperCol* const hyper = (const DirichletHyperCol* const) dirichletHyperWrapper.getDelegator();
    
    const double* discreteIter = &discrete->values[0];  // indesing-style notation remains from the times where arrays have been used instead of vectors,  
    const double* hyperIter = &hyper->values[0];
    double* thisIter = &this->values[0];
    
    double sumU[VertexProxy::K->get_val()];     // sum accross Dirichlet components
    memset(&sumU[0], 0, sizeof(double)*VertexProxy::K->get_val() );
    for(unsigned int col=0; col<VertexProxy::K->get_val(); ++col)   {
        for(unsigned int i=0; i<getOrder(); ++i)  {
            naturalParameters[i][col] = *discreteIter++ + *hyperIter++;
            sumU[col] += naturalParameters[i][col]; 
        }
    }
    for(unsigned int i=0; i<getOrder(); ++i)  {
        for(unsigned int row=0; row<VertexProxy::K->get_val(); ++row)   {
            *thisIter++ = boost::math::digamma(naturalParameters[i][row]) -
                    boost::math::digamma(sumU[row]);
        }
    }
}