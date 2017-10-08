/* 
 * File:   SupportXCol.cpp
 * Author: ivek
 * 
 * Created on January 20, 2012, 10:52 AM
 */

#include "SupportXCol.h"
#include "VertexVisitor.h"

#include <iostream>
#include <string>
#include <vector>


SupportXCol::SupportXCol(unsigned int numElements)
{
    this->numElements = numElements;
    this->support = new unsigned int[this->getNumElements()];
}

SupportXCol::SupportXCol(const SupportXCol& orig)
{
    // deep copying the members
    this->numElements = orig.getNumElements();
    this->support = new unsigned int[this->getNumElements()];
    memcpy(this->support, orig.support, sizeof (unsigned int) * this->numElements);
}

SupportXCol::~SupportXCol() {
    if( this->numElements > 0 )       {
        delete(support);
    }
}

UpdateFunctionDelegator* SupportXCol::clone() {
    UpdateFunctionDelegator* clone = new SupportXCol(*this);
    return (UpdateFunctionDelegator*) clone;
}

void SupportXCol::accept(VertexVisitor &v, gl::iscope& scope, gl::icallback& schedule) {
    v.visit(this, scope, schedule);
}

void SupportXCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
}

/* Note: if there are all zeros in this column, this method will crash - so make
 * sure this isn't the case before invoking it.
 * 
 * TODO: the -1 signal as well as datatypes: long & unsigned int
 */
unsigned int SupportXCol::binarySearch(const unsigned int& rowIndex) const {

    long mid;
    long first = 0;
    long last = getNumElements() - 1;

    while (first <= last) {
        mid = (first + last) / 2; //  split space in half
        if (rowIndex > support[mid]) {
            first = mid + 1; // proceed in the top half
        } else if (rowIndex < support[mid]) {           
            last = mid - 1; // proceed search in bottom half
        } else {
            return mid;
        }
    }
    return -1;
} 

