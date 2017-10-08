/* 
 * File:   DirichletHyperCol.cpp
 * Author: ivek
 * 
 * Created on January 31, 2012, 6:25 PM
 */

#include "DirichletHyperCol.h"
#include "VertexVisitor.h"

#include "math.h"
#include <boost/math/special_functions/digamma.hpp>

DirichletHyperCol::DirichletHyperCol(unsigned int order)
//        : values(order*VertexProxy::K->get_val())
{
//    std::cout << "in DirichletHyperCol constructor" << std::endl;
    this->order = order;
    this->values = new element_type[VertexProxy::K->get_val() * this->order];
}

DirichletHyperCol::DirichletHyperCol(const DirichletHyperCol& orig)
//        : values(orig.getOrder()*VertexProxy::K->get_val())
{
//    std::cout << "in DirichletHyperCol copy constructor" << std::endl;
    
//    std::copy(orig.values.begin(), orig.values.end(), this->values.begin());
    
    this->order = orig.order;
    this->values = new element_type[VertexProxy::K->get_val() * this->order];
    memcpy(this->values, orig.values, sizeof (element_type) * VertexProxy::K->get_val() * this->order);
};

DirichletHyperCol::~DirichletHyperCol() {
//    std::cout << "in DirichletHyperCol destructor" << std::endl;
    delete(values);
}

UpdateFunctionDelegator* DirichletHyperCol::clone() {
    UpdateFunctionDelegator* clone = new DirichletHyperCol(*this);
    return (UpdateFunctionDelegator*) clone;
}

void DirichletHyperCol::accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) {
//    std::cout << "DirichletHyperCol accepted" << std::endl;
    v.visit(this, scope, schedule);
}

void DirichletHyperCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
    std::cout << "DirichletHyperCol update invoked, " << scope.color() <<std::endl;   
    
}

/*
void DirichletHyperCol::contributionToParent( const VCol* const parent,
    const DiscreteCol* const discrete, element_type* messageAccumulator ) const {
    std::cout<<"        DirichletHyperCol::contributionToParent"<<std::endl;
}

void DirichletHyperCol::contributionToChild( const DirichletCol* const child,
    const DiscreteCol* const discrete, element_type* messageAccumulator ) const {
//    std::cout<<"        DirichletHyperCol::contributionToChild"<<std::endl;          
}

void DirichletHyperCol::acceptAsParent( const MessageCollector& child_, const DirichletCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {
    std::cout<<"DirichletHyperCol::acceptAsParent"<<std::endl;
//    v.visitParent(this, messageAccumulator, messageComponents);
}
*/

