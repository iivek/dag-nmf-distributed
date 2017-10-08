/* 
 * File:   XCol.cpp
 * Author: ivek
 * 
 * Created on January 17, 2012, 10:28 AM
 */

#include "XCol.h"
#include "VertexVisitor.h"

XCol::XCol(unsigned int numElements)
//        :values(numElements)
{
//    std::cout << "XCol constructor" << std::endl;
    this->numElements = numElements;
    this->values = new element_type[this->numElements];
}

XCol::XCol(const XCol& orig)
//        :values(orig.getNumElements())
{
//    std::cout << "XCol copy constructor" << std::endl;
    // deep copying the members
    this->numElements = orig.numElements;
    this->values = new element_type[this->numElements];
    memcpy(this->values, orig.values, sizeof (element_type) * this->numElements);
    
//    std::copy(orig.values.begin(), orig.values.end(), this->values.begin());
}

XCol::~XCol() {
//    std::cout << "XCol destructor" << std::endl;
    delete(values);
}

UpdateFunctionDelegator* XCol::clone() {
    UpdateFunctionDelegator* clone = new XCol(*this);
    return (UpdateFunctionDelegator*) clone;
}

void XCol::accept(VertexVisitor &v, gl::iscope& scope, gl::icallback& schedule) {
    //        std::cout << "SparseCol accepted" << std::endl;
    v.visit(this, scope, schedule);
}

void XCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
//    std::cout << "XCol update invoked, " << scope.color() <<std::endl;
}
