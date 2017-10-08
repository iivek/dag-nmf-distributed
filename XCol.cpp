/* 
 * File:   XCol.cpp
 * Author: ivek
 * 
 * Created on January 17, 2012, 10:28 AM
 */

#include "XCol.h"
#include "VertexVisitor.h"

XCol::XCol(unsigned int numElements)
{
    this->numElements = numElements;
    this->values = new element_type[this->numElements];
}

XCol::XCol(const XCol& orig)
{
    // deep copying the members
    this->numElements = orig.numElements;
    this->values = new element_type[this->numElements];
    memcpy(this->values, orig.values, sizeof (element_type) * this->numElements);
}

XCol::~XCol() {
    delete(values);
}

UpdateFunctionDelegator* XCol::clone() {
    UpdateFunctionDelegator* clone = new XCol(*this);
    return (UpdateFunctionDelegator*) clone;
}

void XCol::accept(VertexVisitor &v, gl::iscope& scope, gl::icallback& schedule) {
    v.visit(this, scope, schedule);
}

void XCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
}
