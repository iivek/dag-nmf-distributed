/* 
 * File:   DiscreteCol.cpp
 * Author: ivek
 * 
 * Created on January 31, 2012, 5:50 PM
 */

#include "DiscreteCol.h"
#include "VertexVisitor.h"

#include "math.h"
#include "DirichletCol.h"
#include "ShapeHyperCol.h"
#include "GenericScale.h"
#include "VCol.h"
#include "AuxCol.h"
#include "math.h"
#include <boost/math/special_functions/gamma.hpp>

DiscreteCol::DiscreteCol(unsigned int mixtureSize)
{
    this->mixtureSize = mixtureSize;
    this->values = new element_type[VertexProxy::K->get_val() * this->mixtureSize];
}

DiscreteCol::DiscreteCol(const DiscreteCol& orig)
{
    this->mixtureSize = orig.mixtureSize;
    this->values = new element_type[VertexProxy::K->get_val() * this->mixtureSize];
    memcpy(values, orig.values, sizeof (element_type) * VertexProxy::K->get_val() * mixtureSize);
}

DiscreteCol::~DiscreteCol() {
    delete(values);
}

UpdateFunctionDelegator* DiscreteCol::clone() {
    UpdateFunctionDelegator* clone = new DiscreteCol(*this);
    return (UpdateFunctionDelegator*) clone;
}

void DiscreteCol::accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) {
    v.visit(this, scope, schedule);
}

void DiscreteCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
    /* TODO: SIMDs
     */
    /* Accumulators for messages */
    element_type naturalParameters[getMixtureSize()][VertexProxy::K->get_val()];
    gl::edge_list in_edges = scope.in_edge_ids();

    const unsigned int hyperOffset = 1;

    /* Neighbors:
     * 1. ShapeHyperCol of the VCol
     * 2. DirichletCol
     * 3. AuxCols (parents of the VCol)
     * 4. the VCol
     */   

    const VertexProxy& shapeWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ 0 ]));
    const ShapeHyperCol * const shape = (const ShapeHyperCol* const) shapeWrapper.getDelegator();

    const VertexProxy& vColWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ in_edges.size() - 1 ]));
    const VCol* const v = (const VCol * const) vColWrapper.getDelegator();

    // additive part which is same accross all mixture elements - calculated
    // only once, then copied to all mixture elements
    unsigned int mix = 0;
    for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row) {
        naturalParameters[mix][row] = (shape->aValues[row] - 1) * log(v->expELogValues[row])
                - boost::math::lgamma(shape->aValues[row]);
    }
    
    for (unsigned int mix = 1; mix < getMixtureSize(); ++mix) {
        memcpy(&naturalParameters[mix][0], &naturalParameters[0][0], sizeof(element_type) * VertexProxy::K->get_val());
    }
    {
        const VertexProxy& dirichletWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ 1 ]));
        const DirichletCol* const dirichlet = (const DirichletCol * const) dirichletWrapper.getDelegator();
                
        // add dirichlet->values to natural parameters        
        const double* p = &dirichlet->values[0];
        for (mix = 0; mix < getMixtureSize(); ++mix) {
            for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row) {
                naturalParameters[mix][row] += *p++;
            }
        }
    }

    size_t current = 2;
    for (mix = 0; mix < getMixtureSize(); ++mix) {
        const VertexProxy& parentWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ current++ ]));
        const AuxCol* const parent = (const AuxCol* const) parentWrapper.getDelegator();
        for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row) {
            double temp = shape->aValues[row] * parent->scaleRelatedValues[row];
            naturalParameters[mix][row] += -temp * v->scaleRelatedValues[row] +
                    shape->aValues[row]*log(temp);
        }
    } 
    
    /* Turning the natural parameters into probabilities - normalization
     * Numericals issue with normalization stem from exp(of_a_large_negative_number)
     * in its naive form.
     */    
//    double normalizingFactors[VertexProxy::K->get_val()];
//    memset(normalizingFactors, 0, sizeof(element_type) * VertexProxy::K->get_val() );
//    for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row) {
//        for (mix = 0; mix<getMixtureSize(); ++mix)     {
//            naturalParameters[mix][row] = exp( naturalParameters[mix][row] );            
//            normalizingFactors[row] += naturalParameters[mix][row];
//        }                
//    }
//    double* p = &this->values[0];
//    for (mix = 0; mix < getMixtureSize(); ++mix)     {        
//        for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row) {
//            *p++ = naturalParameters[mix][row]/normalizingFactors[row];
//        }
//    }
    /* Naive form ends */
    
    /* Fixing this by using a log-sum-exp trick. See the Matlab prototype for
     * more comments and explanations.
     */
    // log-sum-exp trick    
    double logNormalizingFactors[VertexProxy::K->get_val()];
    memset(logNormalizingFactors, 0, sizeof(element_type) * VertexProxy::K->get_val() );
    /* First, finding maxima accross rows, for each col. These maxima will be
     * used to better utilize the numerical range before applying the exp function */
    element_type maxima[VertexProxy::K->get_val()];
    memset(maxima, -std::numeric_limits<float>::infinity(), sizeof(element_type) * VertexProxy::K->get_val());
    for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row) {
        for (mix = 0; mix<getMixtureSize(); ++mix)     {
            if( naturalParameters[mix][row] >= maxima[row] ) {
                maxima[row] = naturalParameters[mix][row];
            }
        }
    }
    // The actual calculation
    for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row) {
        for (mix = 0; mix<getMixtureSize(); ++mix)     {
            logNormalizingFactors[row] += exp (
                naturalParameters[mix][row]-maxima[row] );
        }
        logNormalizingFactors[row] = maxima[row] + log(logNormalizingFactors[row]);
    }
    double* p = &this->values[0];
    for (mix = 0; mix < getMixtureSize(); ++mix)     {        
        for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row) {
            *p++ = exp(naturalParameters[mix][row]-logNormalizingFactors[row]);
        }
    }
    // /log-sum-exp trick
}