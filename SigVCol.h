/* 
 * File:   SigVCol.h
 * Author: ivek
 *
 * Created on January 19, 2012, 3:48 PM
 */

#ifndef SIGVCOL_H
#define	SIGVCOL_H

#include "UpdateFunctionDelegator.h"
#include "VertexProxy.h"
#include "VertexVisitor.h"
//#include "MessagePasser.h"
#include <vector>

class ParentVisitor;
class MessageCollector;

class SigVCol : public UpdateFunctionDelegator
{
public:

    SigVCol();
    SigVCol(const SigVCol& orig);
    ~SigVCol();
    UpdateFunctionDelegator* clone();

    /* elements of row this class describes, ordered by increasing indices */
    element_type* values;
    //std::vector<element_type> values;

    /* UpdateFunctionDelegator*/
    void accept(VertexVisitor &v, gl::iscope& scope, gl::icallback& schedule);
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);
    /* UpdateFunctionDelegator end */

private:

};

class SigVColWrapper : public VertexProxy {
public:

    SigVColWrapper() {
        this->essence = new SigVCol();
    }

    SigVColWrapper(const SigVColWrapper& orig) {
        this->essence = orig.essence->clone();
    };

    ~SigVColWrapper() {
        delete(this->essence);
    }
};

#endif	/* SIGVCOL_H */

