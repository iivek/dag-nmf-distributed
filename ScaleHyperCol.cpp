/* 
 * File:   ScaleHyperCol.cpp
 * Author: ivek
 * 
 * Created on February 16, 2012, 10:35 AM
 */

#include "ScaleHyperCol.h"
#include "VCol.h"
//#include "MessageCollector.h"
#include "VertexVisitor.h"

#include <iostream>
#include <string>
#include <vector>

ScaleHyperCol::ScaleHyperCol() {
}

ScaleHyperCol::ScaleHyperCol(const ScaleHyperCol& orig) : GenericScale(orig) {
}

ScaleHyperCol::~ScaleHyperCol() {
}

UpdateFunctionDelegator* ScaleHyperCol::clone() {
    UpdateFunctionDelegator* clone = new ScaleHyperCol(*this);
    return (UpdateFunctionDelegator*) clone;
}

void ScaleHyperCol::accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) {
    v.visit(this, scope, schedule);
}

void ScaleHyperCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
}