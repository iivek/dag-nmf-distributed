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
{
    this->aValues = new element_type[VertexProxy::K->get_val()];
}

ShapeHyperCol::ShapeHyperCol(const ShapeHyperCol& orig)
{
    this->aValues = new element_type[VertexProxy::K->get_val()];
    memcpy(this->aValues, orig.aValues, sizeof(element_type)*VertexProxy::K->get_val());
};

ShapeHyperCol::~ShapeHyperCol() {
    delete(this->aValues);
}

UpdateFunctionDelegator* ShapeHyperCol::clone() {
    UpdateFunctionDelegator* clone = new ShapeHyperCol(*this);
    return (UpdateFunctionDelegator*) clone;
}

void ShapeHyperCol::accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) {
    v.visit(this, scope, schedule);
}

void ShapeHyperCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
}