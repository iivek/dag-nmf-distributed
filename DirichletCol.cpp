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
//        : values(order*VertexProxy::K->get_val())
{
//    std::cout << "in DirichletCol constructor" << std::endl;
    
    this->order = order;
    this->values = new element_type[VertexProxy::K->get_val() * this->order];
}

DirichletCol::DirichletCol(const DirichletCol& orig)
//        : values(orig.getOrder()*VertexProxy::K->get_val())
{
//    std::cout << "in DirichletCol copy constructor" << std::endl;
    
//    std::copy(orig.values.begin(), orig.values.end(), this->values.begin());
    
    this->order = orig.order;
    this->values = new element_type[VertexProxy::K->get_val() * this->order];
    memcpy(this->values, orig.values, sizeof (element_type) * VertexProxy::K->get_val() * this->order);    

};

DirichletCol::~DirichletCol() {
//    std::cout << "in DirichletCol destructor" << std::endl;
    delete(values);
}

UpdateFunctionDelegator* DirichletCol::clone() {
    UpdateFunctionDelegator* clone = new DirichletCol(*this);
    return (UpdateFunctionDelegator*) clone;
}

void DirichletCol::accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) {
//    std::cout << "DirichletCol accepted" << std::endl;
    v.visit(this, scope, schedule);
}

void DirichletCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
    std::cout << "      DirichletCol update invoked, " << scope.color() <<std::endl;    
    
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
    
//    std::cout<<"Naturel parameters begin.";
    double sumU[VertexProxy::K->get_val()];     // sum accross Dirichlet components
    memset(&sumU[0], 0, sizeof(double)*VertexProxy::K->get_val() );
    // everything is organised accross rows (row by row)        
    for(unsigned int col=0; col<VertexProxy::K->get_val(); ++col)   {
        for(unsigned int i=0; i<getOrder(); ++i)  {
//            std::cout << *discreteIter << " " << *hyperIter << " ";
            naturalParameters[i][col] = *discreteIter++ + *hyperIter++;
            sumU[col] += naturalParameters[i][col]; 
//            std::cout << naturalParameters[i][col] << " ";
        }
//        std::cout<<sumU[col] << " " ;            
//        std::cout << std::endl;
    }
    
//    std::cout<<"Naturel parameters end.";
//    double sumU[getOrder()];     // sum accross Dirichlet components
//    memset(&sumU[0], 0, sizeof(double)*VertexProxy::K->get_val() );           
//    for(unsigned int i=0; i<getOrder(); ++i)  {
//        for(unsigned int col=0; col<VertexProxy::K->get_val(); ++col)   {
//            sumU[col] += naturalParameters[i][col];            
//        }
//    }
//    std::cout<<"sumu bergins";
//        for(unsigned int col=0; col<VertexProxy::K->get_val(); ++col)   {
//            std::cout<<sumU[col] << " " ;            
//        }
    std::cout<<std::endl;
    for(unsigned int i=0; i<getOrder(); ++i)  {
        for(unsigned int row=0; row<VertexProxy::K->get_val(); ++row)   {
            *thisIter++ = boost::math::digamma(naturalParameters[i][row]) -
                    boost::math::digamma(sumU[row]);
//            *thisIter = boost::math::digamma(naturalParameters[i][row]) -
//                    boost::math::digamma(sumU[row]);
//            std::cout << *thisIter++ << " ";
        }
//        std::cout<<std::endl;
    }
//    std::cout<<"New values end";
}

/*
void DirichletCol::contributionToParent( const DirichletHyperCol* const parent,
    const DiscreteCol* const discrete, element_type* messageAccumulator ) const {
//    std::cout<<"        DirichletCol::contributionToParent"<<std::endl;
    
}
void DirichletCol::contributionToChild( const DiscreteCol* const child,
    const DiscreteCol* const discrete, element_type* messageAccumulator ) const {
//    std::cout<<"        DirichletCol::contributionToChild"<<std::endl;
    
    
}

void DirichletCol::acceptAsParent( const MessageCollector& child_, const VCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {
    std::cout<<"DirichletCol::acceptAsParent"<<std::endl;
//    v.visitParent(this, messageAccumulator, messageComponents);
}


void DirichletCol::acceptAsChild( const MessageCollector& parent_, const DirichletHyperCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {
    std::cout<<"DirichletCol::acceptAsChild"<<std::endl;
    // not implemented
}
*/