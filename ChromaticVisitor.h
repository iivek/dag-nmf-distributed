/* 
 * File:   ChromaticVisitor.h
 * Author: ivek
 *
 * Created on January 9, 2012, 5:12 PM
 */

#ifndef CHROMATICVISITOR_H
#define	CHROMATICVISITOR_H

#include "VertexVisitor.h"

class TRow;
class VCol;
class TmpCol;
class XCol;
class SupportXCol;
class SigTRow;
class SigVCol;
class THyperRow;
class VHyperCol;
class ShapeHyperCol;
class AuxCol;
class DirichletCol;
class DiscreteCol;
class DirichletHyperCol;

class ChromaticVisitor : public VertexVisitor {

    void visit(TRow *e, gl::iscope& scope, gl::icallback& schedule);
    void visit(VCol *e, gl::iscope& scope, gl::icallback& schedule);
    void visit(TmpCol *e, gl::iscope& scope, gl::icallback& schedule);
    void visit(XCol *e, gl::iscope& scope, gl::icallback& schedule);
    void visit(SupportXCol *e, gl::iscope& scope, gl::icallback& schedule);
    void visit(SigTRow *e, gl::iscope& scope, gl::icallback& schedule);
    void visit(SigVCol *e, gl::iscope& scope, gl::icallback& schedule);
    void visit(THyperRow *e, gl::iscope& scope, gl::icallback& schedule);    
    void visit(ShapeHyperCol *e, gl::iscope& scope, gl::icallback& schedule);
    void visit(ScaleHyperCol *e, gl::iscope& scope, gl::icallback& schedule);
    void visit(AuxCol *e, gl::iscope& scope, gl::icallback& schedule);
    void visit(DirichletCol *e, gl::iscope& scope, gl::icallback& schedule);
    void visit(DiscreteCol *e, gl::iscope& scope, gl::icallback& schedule);
    void visit(DirichletHyperCol *e, gl::iscope& scope, gl::icallback& schedule);
};

#endif	/* CHROMATICVISITOR_H */

