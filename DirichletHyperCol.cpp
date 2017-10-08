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
{
    this->order = order;
    this->values = new element_type[VertexProxy::K->get_val() * this->order];
}

DirichletHyperCol::DirichletHyperCol(const DirichletHyperCol& orig)
{
    this->order = orig.order;
    this->values = new element_type[VertexProxy::K->get_val() * this->order];
    memcpy(this->values, orig.values, sizeof (element_type) * VertexProxy::K->get_val() * this->order);
};

DirichletHyperCol::~DirichletHyperCol() {
    delete(values);
}

UpdateFunctionDelegator* DirichletHyperCol::clone() {
    UpdateFunctionDelegator* clone = new DirichletHyperCol(*this);
    return (UpdateFunctionDelegator*) clone;
}

void DirichletHyperCol::accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) {
    v.visit(this, scope, schedule);
}

void DirichletHyperCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
}
