/* 
 * File:   THyperRow.cpp
 * Author: ivek
 * 
 * Created on January 21, 2012, 5:27 PM
 */

#include "THyperRow.h"
#include "VertexVisitor.h"
#include "TRow.h"

#include <iostream>
#include <string>
#include <vector>

THyperRow::THyperRow()
{
    this->aValues = new element_type[VertexProxy::K->get_val()];
    this->bValues = new element_type[VertexProxy::K->get_val()];
}

THyperRow::THyperRow(const THyperRow& orig)
{
    this->aValues = new element_type[VertexProxy::K->get_val()];
    this->bValues = new element_type[VertexProxy::K->get_val()];
    memcpy(this->aValues, orig.aValues, sizeof(element_type)*VertexProxy::K->get_val());
    memcpy(this->bValues, orig.bValues, sizeof(element_type)*VertexProxy::K->get_val());
};

THyperRow::~THyperRow() {
    delete(this->aValues);
    delete(this->bValues);
}

UpdateFunctionDelegator* THyperRow::clone() {
    UpdateFunctionDelegator* clone = new THyperRow(*this);
    return (UpdateFunctionDelegator*) clone;
}

void THyperRow::accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) {
    v.visit(this, scope, schedule);
}

void THyperRow::updateFunction(gl::iscope& scope, gl::icallback& scheduler) { 
}

