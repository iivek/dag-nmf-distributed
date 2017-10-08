/* 
 * File:   SigTRow.h
 * Author: ivek
 *
 * Created on January 19, 2012, 11:27 AM
 */

#ifndef SIGTROW_H
#define	SIGTROW_H

#include "UpdateFunctionDelegator.h"
#include <vector>
/* Wrapper related */
#include "VertexProxy.h"

class SigTRow : public UpdateFunctionDelegator  {
public:
    
    SigTRow();
    SigTRow(const SigTRow& orig);
    ~SigTRow();    
    UpdateFunctionDelegator* clone();
            
    /* elements of row this class describes, ordered by increasing indices */
    element_type* values;
    
    void accept(VertexVisitor &v, gl::iscope& scope, gl::icallback& schedule);
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);
    
private:
    
};

class SigTRowWrapper : public VertexProxy {
public:

    SigTRowWrapper(){        
        this->essence = new SigTRow();
    }
    SigTRowWrapper(const SigTRowWrapper& orig)      {        
        this->essence = orig.essence->clone();
    };
    ~SigTRowWrapper() {
        delete(this->essence);
    }
};

#endif	/* SIGTROW_H */

