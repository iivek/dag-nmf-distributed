/* 
 * File:   VCol.cpp
 * Author: ivek
 * 
 * Created on January 17, 2012, 10:29 AM
 */

#include "VCol.h"

//#include "MessagePasser.h"
#include "VertexVisitor.h"

#include "SigVCol.h"
#include "TRow.h"
#include "DiscreteCol.h"

#include "math.h"
#include "ShapeHyperCol.h"
#include "AuxCol.h"
#include <boost/math/special_functions/digamma.hpp>

VCol::VCol(bool hasDiscrete, bool hasGammaChild)
//        :expELogValues(VertexProxy::K->get_val())
{
//    std::cout << "in VCol constructor" << std::endl;
//    this->hasDiscrete = hasDiscrete;
//    this->hasGammaChild = hasGammaChild;
    
    this->hasDiscrete = hasDiscrete;
    this->hasGammaChild = hasGammaChild;
    this->expELogValues = new element_type[VertexProxy::K->get_val()];
}

VCol::VCol(const VCol& orig) : GenericScale(orig)
//        ,expELogValues(VertexProxy::K->get_val())
{
    //    std::cout << "in VCol copy constructor" << std::endl;
//    this->hasDiscrete = orig.hasDiscrete;
//    this->hasGammaChild = orig.hasGammaChild;    
//    std::copy(orig.expELogValues.begin(), orig.expELogValues.end(), this->expELogValues.begin());
    
    this->hasDiscrete = orig.hasDiscrete;
    this->hasGammaChild = orig.hasGammaChild;    
    this->expELogValues = new element_type[VertexProxy::K->get_val()];
    memcpy(this->expELogValues, orig.expELogValues, sizeof (element_type) * VertexProxy::K->get_val());
};

VCol::~VCol() {
    //    std::cout << "in VCol destructor" << std::endl;
    delete(this->expELogValues);
}

UpdateFunctionDelegator* VCol::clone() {
    UpdateFunctionDelegator* clone = new VCol(*this);
    return (UpdateFunctionDelegator*) clone;
}

void VCol::accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) {
    //        std::cout << "VCol accepted" << std::endl;
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
    std::cout << "      VCol update invoked " << scope.color() << std::endl;

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
          std::cout << "shapeDelegator1=" << shapeDelegator->aValues[0] <<  "; "<< std::endl;
        i = 2;
    } else {
        // no coparent - shape is first neighbor
        const VertexProxy& shapeParent = scope.neighbor_vertex_data(scope.source(in_edges[ 0 ]));
        shapeDelegator = (const ShapeHyperCol* const) shapeParent.getDelegator();
        std::cout << "shapeDelegator2=" << shapeDelegator->aValues[0] <<  "; " << std::endl;
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
                    std::cout << "shape=" << shapeDelegator->aValues[row] << " scale= " << currDelegator->scaleRelatedValues[row] << "; discrete_iter=" << discreteDelegator->values[discreteIter] << "; ";
                    std::cout << -shapeDelegator->aValues[row] * currDelegator->scaleRelatedValues[row] *
                            discreteDelegator->values[discreteIter] << "; ";
                    discreteIter++;
                }
                std::cout << "out."<< std::endl;
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
                std::cout << "shape= " << auxDelegator->scaleRelatedValues[row] << "; coparent_shape =" << coparent->aValues[row] << "; ";
            }
            std::cout << std::endl;
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
              std::cout << "parent shape=" << shapeDelegator->aValues[row] <<  "; ";
                  
        }
        std::cout << std::endl;
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


    //    for (unsigned int i = 0; i < in_edges.size(); ++i) {
    //        /* Check if neighbor is a parent or a child */
    //        if (scope.source(in_edges[i]) < ourID) {
    //            std::cout<<"currentID of parent: "<<scope.source(in_edges[i])<<std::endl;
    //            /* Parent */
    //            const VertexProxy& parent = scope.neighbor_vertex_data(
    //                    scope.source(in_edges[ i ]));
    //            const MessagePasser& messageDelegator =
    //                    (const MessagePasser&) ((const VHyperCol&) *parent.getDelegator());
    //            messageDelegator.acceptAsParent((MessageCollector&) *this, *this, discrete, naturalParameters);
    //        } else {
    //            /* Child */
    //            const VertexProxy& child = scope.neighbor_vertex_data(
    //                    scope.source(in_edges[ i ]));
    //            const MessagePasser& messageDelegator =
    //                    (const MessagePasser&) ((const VHyperCol&) *child.getDelegator());
    //            messageDelegator.acceptAsChild((MessageCollector&) *this, *this, discrete, naturalParameters);            
    //        }
    //    }
    /* Collecting messages done */

    /* TODO: SIMD
     */
    //    element_type alpha[VertexProxy::K->get_val()];
    //    element_type beta[VertexProxy::K->get_val()];
    for (i = 0; i < VertexProxy::K->get_val(); ++i) {

        element_type alpha = naturalParameters[1][i] + 1;
        element_type beta = -1.0 / naturalParameters[0][i];
        this->scaleRelatedValues[i] = alpha*beta;
        this->expELogValues[i] = exp(boost::math::digamma(alpha)) * beta;
//        std::cout << this->scaleRelatedValues[i] << " ";
    }
//    std::cout << std::endl;
 
}

/*
void VCol::contributionToParent(const AuxCol * const parent,
        const DiscreteCol * const discrete, element_type* messageAccumulator) const {
//    std::cout << "      VCol::contributionToParent" << std::endl;
}

void VCol::contributionToChild(const SigVCol * const child,
        const DiscreteCol * const discrete, element_type* messageAccumulator) const {
//    std::cout << "      VCol::contributionToChild" << std::endl;
}
 */
/* MessagePasser stuff
 */
/*
void VCol::acceptAsParent(const MessageCollector& child_, const SigVCol& child,
        const DiscreteCol * const discrete, element_type* messageAccumulator) const {
    std::cout << "VCol::acceptAsParent, SigVCol" << std::endl;
}

void VCol::acceptAsParent(const MessageCollector& child_, const TRow& child,
        const DiscreteCol * const discrete, element_type* messageAccumulator) const {
    std::cout << "VCol::acceptAsParent, TRow" << std::endl;
}

void VCol::acceptAsParent(const MessageCollector& child_, const AuxCol& child,
        const DiscreteCol * const discrete, element_type* messageAccumulator) const {
    std::cout << "VCol::acceptAsParent, AuxCol" << std::endl;
}

void VCol::acceptAsChild(const MessageCollector& parent_, const VHyperCol& parent,
        const DiscreteCol * const discrete, element_type* messageAccumulator) const {
    std::cout << "VCol::acceptAsChild, VHyperCol" << std::endl;
}

void VCol::acceptAsChild(const MessageCollector& parent_, const ShapeHyperCol& parent,
        const DiscreteCol * const discrete, element_type* messageAccumulator) const {
    std::cout << "VCol::acceptAsChild, ShapeHyperCol" << std::endl;
}

void VCol::acceptAsChild(const MessageCollector& parent_, const AuxCol& parent,
        const DiscreteCol * const discrete, element_type* messageAccumulator) const {
    std::cout << "VCol::acceptAsChild, AuxCol" << std::endl;
}

void VCol::acceptAsChild(const MessageCollector& parent_, const DiscreteCol& parent,
        const DiscreteCol * const discrete, element_type* messageAccumulator) const {
    std::cout << "VCol::acceptAsChild, DiscreteCol" << std::endl;
}
 */
/* MessagePasser ends
 */