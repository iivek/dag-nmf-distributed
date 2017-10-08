/* 
 * File:   VCol.cpp
 * Author: ivek
 * 
 * Created on January 17, 2012, 10:29 AM
 */

#include "VCol.h"

#include "VertexVisitor.h"

#include "SigVCol.h"
#include "TRow.h"
#include "DiscreteCol.h"

#include "math.h"
#include "ShapeHyperCol.h"
#include "AuxCol.h"
#include <boost/math/special_functions/digamma.hpp>

VCol::VCol(bool hasDiscrete, bool hasGammaChild)
{
    this->hasDiscrete = hasDiscrete;
    this->hasGammaChild = hasGammaChild;
    this->expELogValues = new element_type[VertexProxy::K->get_val()];
}

VCol::VCol(const VCol& orig) : GenericScale(orig)
{
    this->hasDiscrete = orig.hasDiscrete;
    this->hasGammaChild = orig.hasGammaChild;    
    this->expELogValues = new element_type[VertexProxy::K->get_val()];
    memcpy(this->expELogValues, orig.expELogValues, sizeof (element_type) * VertexProxy::K->get_val());
};

VCol::~VCol() {
    delete(this->expELogValues);
}

UpdateFunctionDelegator* VCol::clone() {
    UpdateFunctionDelegator* clone = new VCol(*this);
    return (UpdateFunctionDelegator*) clone;
}

void VCol::accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) {
    v.visit(this, scope, schedule);
}

/**
 * Topologically ordered vertices come into play at this point; such an ordering
 * is required to determine whether a neighbor is child or parent 
 * 
 * @param scope
 * @param scheduler
 */
void VCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
    /* Accumulators for messages */
    element_type naturalParameters[2][VertexProxy::K->get_val()];
    memset(naturalParameters, 0, sizeof (element_type)*2 * VertexProxy::K->get_val());

    gl::edge_list in_edges = scope.in_edge_ids();
    /* Traverse the neighbors and collect messages
     */
    gl::vertex_id ourID = scope.vertex();

    const ShapeHyperCol* coparent;
    const ShapeHyperCol* shapeDelegator;
    const DiscreteCol* discreteDelegator;

    unsigned int i; // neighbor iterator
    if (hasChild()) {
        // coparent is the first neighbor, followed by our shape
        const VertexProxy& coparentWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ 0 ]));
        coparent = (const ShapeHyperCol* const) coparentWrapper.getDelegator();

        const VertexProxy& shapeParent = scope.neighbor_vertex_data(scope.source(in_edges[ 1 ]));
        shapeDelegator = (const ShapeHyperCol * const) shapeParent.getDelegator();
        i = 2;
    } else {
        // no coparent - shape is first neighbor
        const VertexProxy& shapeParent = scope.neighbor_vertex_data(scope.source(in_edges[ 0 ]));
        shapeDelegator = (const ShapeHyperCol* const) shapeParent.getDelegator();
        i = 1;
    }

    if (this->isMixture()) {
        // a discrete parent exists
        const VertexProxy& discreteParent = scope.neighbor_vertex_data(scope.source(in_edges[ i++ ]));
        discreteDelegator = (const DiscreteCol * const) discreteParent.getDelegator();
        unsigned int discreteIter = 0;

        /* TODO: SIMD */
        for (i; i < in_edges.size() - 1; ++i) {
            /* Check if neighbor is a parent or a child */
            if (scope.source(in_edges[i]) < ourID) {
                /* Parent */
                const VertexProxy& currParent = scope.neighbor_vertex_data(scope.source(in_edges[ i ]));
                const GenericScale * const currDelegator = (const GenericScale * const) currParent.getDelegator();
                // TODO: SIMD
                for (int row = 0; row < VertexProxy::K->get_val(); ++row) {
                    // two message components                    
                    naturalParameters[0][row] += -shapeDelegator->aValues[row] * currDelegator->scaleRelatedValues[row] *
                            discreteDelegator->values[discreteIter];
                    naturalParameters[1][row] += (shapeDelegator->aValues[row] - 1) * discreteDelegator->values[discreteIter];
                   discreteIter++;
                }
            } else {
                // reached the children block - for all i'>i we'll have children
                break;
            }
        }
        
        if (hasChild()) {
            // what follows is an AuxCol
            const VertexProxy& aWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ i++ ]));
            const AuxCol * const auxDelegator = (const AuxCol * const) aWrapper.getDelegator();
            for (int row = 0; row < VertexProxy::K->get_val(); ++row) {
                naturalParameters[0][row] += -coparent->aValues[row] * auxDelegator->scaleRelatedValues[row];
                naturalParameters[1][row] += coparent->aValues[row];
            }
        }
        
        for (; i < in_edges.size() - 1; ++i) {
            /* Child is TRow*/
            const VertexProxy& tWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ i ]));
            const TRow * const tDelegator = (const TRow * const) tWrapper.getDelegator();
            for (int col = 0; col < VertexProxy::K->get_val(); ++col) {
                naturalParameters[0][col] -= tDelegator->eValues[col];
            }
        }
        // and the last neighbor is a SigVCol child
        const VertexProxy& sigVWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ in_edges.size() - 1 ]));
        const SigVCol * const sigVDelegator = (const SigVCol * const) sigVWrapper.getDelegator();
        for (int row = 0; row < VertexProxy::K->get_val(); ++row) {
            naturalParameters[1][row] += sigVDelegator->values[row];
        }
    } else {
        // isn't a mixture - we'll have only one gamma parent
        const VertexProxy& currParent = scope.neighbor_vertex_data(scope.source(in_edges[ i++ ]));
        const GenericScale * const currDelegator = (const GenericScale * const) currParent.getDelegator();
        // TODO: SIMD
        for (int row = 0; row < VertexProxy::K->get_val(); ++row) {
            // two message components                    
            naturalParameters[0][row] += -shapeDelegator->aValues[row] * currDelegator->scaleRelatedValues[row];
            naturalParameters[1][row] += (shapeDelegator->aValues[row] - 1);
        }
        if (hasChild()) {
            // what follows is an AuxCol
            const VertexProxy& aWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ i++ ]));
            const AuxCol * const auxDelegator = (const AuxCol * const) aWrapper.getDelegator();
            for (int row = 0; row < VertexProxy::K->get_val(); ++row) {
                naturalParameters[0][row] += -coparent->aValues[row] * auxDelegator->scaleRelatedValues[row];
                naturalParameters[1][row] += coparent->aValues[row];
            }
        }

        for (; i < in_edges.size()-1; ++i) {
            /* Child is TRow*/
            const VertexProxy& tWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ i ]));
            const TRow* const tDelegator = (const TRow* const) tWrapper.getDelegator();
            for (int col = 0; col < VertexProxy::K->get_val(); ++col) {
                naturalParameters[0][col] -= tDelegator->eValues[col];
            }
        }
        // and the last neighbor is a SigVCol child
        const VertexProxy& sigVWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ in_edges.size() - 1 ]));
        const SigVCol * const sigVDelegator = (const SigVCol * const) sigVWrapper.getDelegator();
        for (int row = 0; row < VertexProxy::K->get_val(); ++row) {
            naturalParameters[1][row] += sigVDelegator->values[row];
        }
        

    }
    /* TODO: SIMD
     */
    //    element_type alpha[VertexProxy::K->get_val()];
    //    element_type beta[VertexProxy::K->get_val()];
    for (i = 0; i < VertexProxy::K->get_val(); ++i) {

        element_type alpha = naturalParameters[1][i] + 1;
        element_type beta = -1.0 / naturalParameters[0][i];
        this->scaleRelatedValues[i] = alpha*beta;
        this->expELogValues[i] = exp(boost::math::digamma(alpha)) * beta;
    }
 
}