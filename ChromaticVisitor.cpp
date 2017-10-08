/* 
 * File:   ChromaticVisitor.cpp
 * Author: ivek
 * 
 * Created on January 9, 2012, 5:12 PM
 */

#include "ChromaticVisitor.h"
#include "TRow.h"
#include "VCol.h"
#include "TmpCol.h"
#include "XCol.h"
#include "SupportXCol.h"
#include "SigTRow.h"
#include "SigVCol.h"
#include "THyperRow.h"
#include "ShapeHyperCol.h"
#include "ScaleHyperCol.h"
#include "AuxCol.h"
#include "DirichletCol.h"
#include "DiscreteCol.h"
#include "DirichletHyperCol.h"

void ChromaticVisitor::visit(TRow *e, gl::iscope& scope, gl::icallback& schedule) {
    e->updateFunction(scope, schedule);
}

void ChromaticVisitor::visit(VCol *e, gl::iscope& scope, gl::icallback& schedule) {
    e->updateFunction(scope, schedule);
}

void ChromaticVisitor::visit(TmpCol *e, gl::iscope& scope, gl::icallback& schedule) {
    e->updateFunction(scope, schedule);
}

void ChromaticVisitor::visit(XCol *e, gl::iscope& scope, gl::icallback& schedule) {
    e->updateFunction(scope, schedule);
}

void ChromaticVisitor::visit(SupportXCol *e, gl::iscope& scope, gl::icallback& schedule) {
    e->updateFunction(scope, schedule);
}

void ChromaticVisitor::visit(SigTRow *e, gl::iscope& scope, gl::icallback& schedule) {
    e->updateFunction(scope, schedule);
}

void ChromaticVisitor::visit(SigVCol *e, gl::iscope& scope, gl::icallback& schedule) {
    e->updateFunction(scope, schedule);
}

void ChromaticVisitor::visit(THyperRow *e, gl::iscope& scope, gl::icallback& schedule) {
    e->updateFunction(scope, schedule);
}

void ChromaticVisitor::visit(ShapeHyperCol *e, gl::iscope& scope, gl::icallback& schedule) {
    e->updateFunction(scope, schedule);
}

void ChromaticVisitor::visit(ScaleHyperCol *e, gl::iscope& scope, gl::icallback& schedule) {
    e->updateFunction(scope, schedule);
}

void ChromaticVisitor::visit(AuxCol *e, gl::iscope& scope, gl::icallback& schedule) {
    e->updateFunction(scope, schedule);        
}

void ChromaticVisitor::visit(DirichletCol *e, gl::iscope& scope, gl::icallback& schedule) {
    e->updateFunction(scope, schedule);        
}

void ChromaticVisitor::visit(DiscreteCol *e, gl::iscope& scope, gl::icallback& schedule) {
    e->updateFunction(scope, schedule);        
}

void ChromaticVisitor::visit(DirichletHyperCol *e, gl::iscope& scope, gl::icallback& schedule) {
    e->updateFunction(scope, schedule);        
}