/* 
 * File:   ATRow.h
 * Author: ivek
 *
 * Created on January 21, 2012, 5:27 PM
 */

#ifndef THYPERROW_H
#define	THYPERROW_H

#include "UpdateFunctionDelegator.h"
#include "VertexProxy.h"
#include <vector>

class THyperRow : public UpdateFunctionDelegator    {
public:            
    THyperRow();
    THyperRow(const THyperRow& orig);    
    ~THyperRow();
    UpdateFunctionDelegator* clone();

    void accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule);
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);
    
    element_type* aValues;       /* Order: ascending indices */
    element_type* bValues;       /* Order: ascending indices */
//    std::vector<element_type> aValues;       /* Order: ascending indices */
//    std::vector<element_type> bValues;       /* Order: ascending indices */

};

class THyperRowWrapper : public VertexProxy {
public:    
    THyperRowWrapper(){        
        this->essence = new THyperRow();
    }
    THyperRowWrapper(const THyperRowWrapper& orig)      {        
        this->essence = orig.essence->clone();
    };
    ~THyperRowWrapper() {
        delete(this->essence);
    }
};

#endif	/* THYPERROW_H */