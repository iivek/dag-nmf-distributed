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
//        : indexInMixture(numChildren)
//        , childIsMixture(numChildren)
//        , eLogValues(VertexProxy::K->get_val())
{
    //    std::cout << "in AuxCol constructor" << std::endl;
    this->numChildren = numChildren;
    this->indexInMixture = new unsigned int[this->getNumChildren()];
    this->childIsMixture = new bool[this->getNumChildren()];
    this->eLogValues = new element_type[VertexProxy::K->get_val()];
        
}

AuxCol::AuxCol(const AuxCol& orig)
        : GenericScale(orig)
//        , indexInMixture(orig.getNumChildren())
//        , childIsMixture(orig.getNumChildren())
//        , eLogValues(VertexProxy::K->get_val())
{        
    //    std::cout << "in AuxCol copy constructor" << std::endl;
    
    this->numChildren = orig.numChildren;
    this->indexInMixture = new unsigned int[this->getNumChildren()];
    this->childIsMixture = new bool[this->getNumChildren()];
    
/*    std::copy(orig.indexInMixture.begin(), orig.indexInMixture.end(), this->indexInMixture.begin());
    std::copy(orig.childIsMixture.begin(), orig.childIsMixture.end(), this->childIsMixture.begin());
    std::copy(orig.eLogValues.begin(), orig.eLogValues.end(), this->eLogValues.begin());
*/        
    memcpy(this->indexInMixture, orig.indexInMixture, sizeof(unsigned int) * getNumChildren());
    memcpy(this->childIsMixture, orig.childIsMixture, sizeof(bool) * getNumChildren());
    this->eLogValues = new element_type[VertexProxy::K->get_val()];
    memcpy(this->eLogValues, orig.eLogValues, sizeof (element_type) * VertexProxy::K->get_val());        
};

AuxCol::~AuxCol() {
    //    std::cout << "in AuxCol destructor" << std::endl;
    delete(this->scaleRelatedValues);
    delete(this->eLogValues);
    delete(this->indexInMixture);
    delete(this->childIsMixture);
}

UpdateFunctionDelegator* AuxCol::clone() {
    UpdateFunctionDelegator* clone = new AuxCol(*this);
    return (UpdateFunctionDelegator*) clone;
}

/*void AuxCol::save(oarchive& oarc) const     {
    oarc << eLogValues << indexInMixture << childIsMixture;
}
    
void AuxCol::load(iarchive& iarc)   {
    iarc >> eLogValues >> indexInMixture >> childIsMixture;
}
 */
/*
void AuxColWrapper::save(oarchive& oarc) const     {
    oarc << ((AuxCol*)this->essence)->eLogValues << ((AuxCol*)this->essence)->indexInMixture << ((AuxCol*)this->essence)->childIsMixture;
}
    
void AuxColWrapper::load(iarchive& iarc)   {
    iarc >> ((AuxCol*)this->essence)->eLogValues;// >> ((AuxCol*)this->essence)indexInMixture >> ((AuxCol*)this->essence)childIsMixture;
}
*/
void AuxCol::accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) {
//    std::cout << "AuxCol accepted" << std::endl;
    v.visit(this, scope, schedule);
}

void AuxCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
    std::cout << "AuxCol update invoked, " << scope.color() << std::endl;
    
    // TODO: SIMDs everywhere        
    
    /* Accumulators for messages */
    element_type naturalParameters[2][VertexProxy::K->get_val()];
    memset(naturalParameters, 0, sizeof (element_type)*2 * VertexProxy::K->get_val());

    gl::edge_list in_edges = scope.in_edge_ids();
    /* Traverse the neighbors and collect messages
     */
//    gl::vertex_id ourID = scope.vertex();
    
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
    
//        std::cout << "as from shape parent" << std::endl;
//            for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row) {
//                    std::cout << shape->aValues[row] << " ";
//                std::cout << std::endl;
//            }
//        
//            std::cout << " scales from vcol parent" << std::endl;
//            for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row) {
//                for (size_t col = 0; col<2; ++col) {
//                    std::cout << parent->scaleRelatedValues[row] << " ";
//                }
//                std::cout << std::endl;
//            }
    
    /* messages from children */
    unsigned int discreteIter = hyperOffset + getNumChildren();
    unsigned int childIter = in_edges.size() - getNumChildren();
    for(unsigned int iter=0; iter<getNumChildren(); ++iter)   {
        const VertexProxy& childShapeWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ hyperOffset+iter ]));
        const ShapeHyperCol* const childShape = (const ShapeHyperCol* const) childShapeWrapper.getDelegator();
        const VertexProxy& childWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ childIter++ ]));
        const VCol* const child = (const VCol* const) childWrapper.getDelegator();
        if( childIsMixture[iter] ) {
//            std::cout << "child is mixture" << std::endl;
            const VertexProxy& discreteWrapper = scope.neighbor_vertex_data(scope.source(in_edges[ discreteIter++ ]));
            const DiscreteCol* const discrete = (const DiscreteCol* const) discreteWrapper.getDelegator();
            
//            double* p = &discrete->values[ indexInMixture[iter]*VertexProxy::K->get_val() ];
//            // points to the forst row in iter-th mixture. ++p points to the next in row
//            for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row)   {
//                naturalParameters[0][row] -= child->scaleRelatedValues[row] * childShape->aValues[row] * (*p);
//                naturalParameters[1][row] += childShape->aValues[row] * (*p);
//                ++p;
//            }
            
            for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row)   {
                naturalParameters[0][row] -= child->scaleRelatedValues[row] * childShape->aValues[row] *
                        discrete->values[ indexInMixture[iter]*VertexProxy::K->get_val() + row];
                naturalParameters[1][row] += childShape->aValues[row] *
                        discrete->values[ indexInMixture[iter]*VertexProxy::K->get_val() + row];
            }
           

        } else  {
//            std::cout << "child is single" << std::endl;
            // child is not a mixture, but a single gamma
             for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row)   {
                naturalParameters[0][row] -= child->scaleRelatedValues[row] * childShape->aValues[row];
                naturalParameters[1][row] += childShape->aValues[row];
            }
        }
        
//                    std::cout << "scales" << std::endl;
//
//            for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row) {
//                for (size_t col = 0; col<2; ++col) {
//                    std::cout << child->scaleRelatedValues[row] << " ";
//                }
//                std::cout << std::endl;
//            }
//            
//            std::cout << "as" << std::endl;
//            for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row) {
//                for (size_t col = 0; col<2; ++col) {
//                    std::cout << childShape->aValues[row] << " ";
//                }
//                std::cout << std::endl;
//            }
    }
   
    for (unsigned int i = 0; i < VertexProxy::K->get_val(); ++i) {
        element_type alpha = naturalParameters[1][i] + 1;
        element_type beta = -1.0 / naturalParameters[0][i];
        this->scaleRelatedValues[i] = alpha*beta;
        this->eLogValues[i] = boost::math::digamma(alpha) + log(beta);
    }
    
}
/*
void AuxCol::contributionToParent(const VCol * const parent,
        const DiscreteCol * const discrete, element_type* messageAccumulator) const {
//    std::cout << "        AuxCol::contributionToParent" << std::endl;

}

void AuxCol::contributionToChild(const VCol * const child,
        const DiscreteCol * const discrete, element_type* messageAccumulator) const {
//    std::cout << "        AuxCol::contributionToChild" << std::endl;
}

    
void AuxCol::acceptAsParent( const MessageCollector& child_, const VCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {
    std::cout << "AuxCol::acceptAsParent" << std::endl;
    //    v.visitParent(this, messageAccumulator);
}

void AuxCol::acceptAsChild( const MessageCollector& parent_, const VCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {
    std::cout << "AuxCol::acceptAsChild" << std::endl;    
}
 */ 