/* Visitor pattern constituent, update function
 * 
 * File:   VertexVisitor.h
 * Author: ivek
 *
 * Created on January 9, 2012, 4:25 PM
 */

#ifndef VERTEXVISITOR_H
#define	VERTEXVISITOR_H

#include "UpdateFunctionDelegator.h"

class TRow;
class VCol;
class TmpCol;
class XCol;
class SupportXCol;
class SigTRow;
class SigVCol;
class THyperRow;
class ShapeHyperCol;
class ScaleHyperCol;
class AuxCol;
class DirichletCol;
class DiscreteCol;
class DirichletHyperCol;

class VertexVisitor {
public:
    virtual void visit(TRow *e, gl::iscope& scope, gl::icallback& schedule) = 0;
    virtual void visit(VCol *e, gl::iscope& scope, gl::icallback& schedule) = 0;
    virtual void visit(TmpCol *e, gl::iscope& scope, gl::icallback& schedule) = 0;
    virtual void visit(XCol *e, gl::iscope& scope, gl::icallback& schedule) = 0;
    virtual void visit(SupportXCol *e, gl::iscope& scope, gl::icallback& schedule) = 0;
    virtual void visit(SigTRow *e, gl::iscope& scope, gl::icallback& schedule) = 0;
    virtual void visit(SigVCol *e, gl::iscope& scope, gl::icallback& schedule) = 0;
    virtual void visit(THyperRow *e, gl::iscope& scope, gl::icallback& schedule) = 0;    
    virtual void visit(ShapeHyperCol *e, gl::iscope& scope, gl::icallback& schedule) = 0;
    virtual void visit(ScaleHyperCol *e, gl::iscope& scope, gl::icallback& schedule) = 0;
    virtual void visit(AuxCol *e, gl::iscope& scope, gl::icallback& schedule) = 0;
    virtual void visit(DirichletCol *e, gl::iscope& scope, gl::icallback& schedule) = 0;
    virtual void visit(DiscreteCol *e, gl::iscope& scope, gl::icallback& schedule) = 0;
    virtual void visit(DirichletHyperCol *e, gl::iscope& scope, gl::icallback& schedule) = 0;    
};


#endif	/* VERTEXVISITOR_H */

