/* 
 * File:   AuxCol.cpp
 * Author: ivek
 * 
 * Created on January 26, 2012, 4:59 PM
 */

#include "AuxCol.h"
#include "VCol.h"
#include "DiscreteCol.h"
#include "VertexVisitor.h"

#include "math.h"
#include "ShapeHyperCol.h"
#include "DiscreteCol.h"
#include <boost/math/special_functions/digamma.hpp>

using namespace graphlab;

AuxCol::AuxCol(unsigned int numChildren)
{
    this->numChildren = numChildren;
    this->indexInMixture = new unsigned int[this->getNumChildren()];
    this->childIsMixture = new bool[this->getNumChildren()];
    this->eLogValues = new element_type[VertexProxy::K->get_val()];
        
}

AuxCol::AuxCol(const AuxCol& orig)
        : GenericScale(orig)
{        
    this->numChildren = orig.numChildren;
    this->indexInMixture = new unsigned int[this->getNumChildren()];
    this->childIsMixture = new bool[this->getNumChildren()];

    memcpy(this->indexInMixture, orig.indexInMixture, sizeof(unsigned int) * getNumChildren());
    memcpy(this->childIsMixture, orig.childIsMixture, sizeof(bool) * getNumChildren());
    this->eLogValues = new element_type[VertexProxy::K->get_val()];
    memcpy(this->eLogValues, orig.eLogValues, sizeof (element_type) * VertexProxy::K->get_val());        
};

AuxCol::~AuxCol() {
    delete(this->scaleRelatedValues);
    delete(this->eLogValues);
    delete(this->indexInMixture);
    delete(this->childIsMixture);
}

UpdateFunctionDelegator* AuxCol::clone() {
    UpdateFunctionDelegator* clone = new AuxCol(*this);
    return (UpdateFunctionDelegator*) clone;
}


void AuxCol::accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) {
    v.visit(this, scope, schedule);
}

void AuxCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
    // TODO: SIMDs everywhere            
    /* Accumulators for messages */
    element_type naturalParameters[2][VertexProxy::K->get_val()];
    memset(naturalParameters, 0, sizeof (element_type)*2 * VertexProxy::K->get_val());

    gl::edge_list in_edges = scope.in_edge_ids();
    /* Traverse the neighbors and collect messages
     */
    const unsigned int hyperOffset = 1;   
        
    /* message from parent */
    const VertexProxy& shapeWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ 0 ]));
    const ShapeHyperCol* const shape = (const ShapeHyperCol* const) shapeWrapper.getDelegator();        // our shape hyperparameter
    const VertexProxy& parentWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ in_edges.size()-getNumChildren()-1 ]));
    const VCol* const parent = (const VCol* const) parentWrapper.getDelegator();                        // our VCol parent
    for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row)   {
        naturalParameters[0][row] -= parent->scaleRelatedValues[row] * shape->aValues[row];        
        naturalParameters[1][row] += shape->aValues[row]-1;
    }
    
    /* messages from children */
    unsigned int discreteIter = hyperOffset + getNumChildren();
    unsigned int childIter = in_edges.size() - getNumChildren();
    for(unsigned int iter=0; iter<getNumChildren(); ++iter)   {
        const VertexProxy& childShapeWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ hyperOffset+iter ]));
        const ShapeHyperCol* const childShape = (const ShapeHyperCol* const) childShapeWrapper.getDelegator();
        const VertexProxy& childWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ childIter++ ]));
        const VCol* const child = (const VCol* const) childWrapper.getDelegator();
        if( childIsMixture[iter] ) {
            const VertexProxy& discreteWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ discreteIter++ ]));
            const DiscreteCol* const discrete = (const DiscreteCol* const) discreteWrapper.getDelegator();
            for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row)   {
                naturalParameters[0][row] -= child->scaleRelatedValues[row] * childShape->aValues[row] *
                        discrete->values[ indexInMixture[iter]*VertexProxy::K->get_val() + row];
                naturalParameters[1][row] += childShape->aValues[row] *
                        discrete->values[ indexInMixture[iter]*VertexProxy::K->get_val() + row];
            }
        } else  {
             for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row)   {
                naturalParameters[0][row] -= child->scaleRelatedValues[row] * childShape->aValues[row];
                naturalParameters[1][row] += childShape->aValues[row];
            }
        }
    }  
    for (unsigned int i = 0; i < VertexProxy::K->get_val(); ++i) {
        element_type alpha = naturalParameters[1][i] + 1;
        element_type beta = -1.0 / naturalParameters[0][i];
        this->scaleRelatedValues[i] = alpha*beta;
        this->eLogValues[i] = boost::math::digamma(alpha) + log(beta);
    }
    
}