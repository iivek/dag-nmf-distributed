/* 
 * File:   ShapeHyperCol.cpp
 * Author: ivek
 * 
 * Created on January 24, 2012, 9:23 AM
 */

#include "ShapeHyperCol.h"
//#include "MessageCollector.h"
#include "VertexVisitor.h"

#include <iostream>
#include <string>
#include <vector>

ShapeHyperCol::ShapeHyperCol()
//        : aValues(VertexProxy::K->get_val())
{
//    std::cout << "in ShapeHyperCol constructor" << std::endl;
    this->aValues = new element_type[VertexProxy::K->get_val()];
}

ShapeHyperCol::ShapeHyperCol(const ShapeHyperCol& orig)
//        : aValues(VertexProxy::K->get_val())
{
    //    std::cout << "in ShapeHyperCol copy constructor" << std::endl;    <   
//    std::copy(orig.aValues.begin(), orig.aValues.end(), this->aValues.begin());
    
    this->aValues = new element_type[VertexProxy::K->get_val()];
    memcpy(this->aValues, orig.aValues, sizeof(element_type)*VertexProxy::K->get_val());
};

ShapeHyperCol::~ShapeHyperCol() {
//    std::cout << "ShapeHyperCol destructor" << std::endl;
    delete(this->aValues);
}

UpdateFunctionDelegator* ShapeHyperCol::clone() {
    UpdateFunctionDelegator* clone = new ShapeHyperCol(*this);
    return (UpdateFunctionDelegator*) clone;
}

void ShapeHyperCol::accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) {
    //        std::cout << "ShapeHyperCol accepted" << std::endl;
    v.visit(this, scope, schedule);
}

void ShapeHyperCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
//    std::cout << "ShapeHyperCol update invoked, " << scope.color() <<std::endl;
}

/*
void ShapeHyperCol::acceptAsParent( const MessageCollector& child_, const VCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {
    std::cout<<"ShapeHyperCol::acceptAsParent, VCol"<<std::endl;
//    return v.visitParent(this, messageAccumulator);
}

void ShapeHyperCol::acceptAsParent( const MessageCollector& child_, const AuxCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {
    std::cout<<"ShapeHyperCol::acceptAsParent, AuxCol"<<std::endl;
    // not implemented
}

void ShapeHyperCol::contributionToParent( const VCol* const parent,
    const DiscreteCol* const discrete, element_type* messageAccumulator ) const {
//    std::cout << "     ShapeHyperCol::contributionToParent"<<std::endl;    
}
void ShapeHyperCol::contributionToChild( const VCol* const child,
    const DiscreteCol* const discrete, element_type* messageAccumulator ) const {
//    std::cout << "      ShapeHyperCol::contributionToChild, child is VCol"<<std::endl;    
}
void ShapeHyperCol::contributionToChild( const AuxCol* const child,
    const DiscreteCol* const discrete, element_type* messageAccumulator ) const {
//    std::cout << "      ShapeHyperCol::contributionToChild, child is AuxCol"<<std::endl;    
}
*/
