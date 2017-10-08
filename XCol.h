/* 
 * File:   XCol.h
 * Author: ivek
 *
 * Created on January 17, 2012, 10:28 AM
 */

#ifndef XCOL_H
#define	XCOL_H

#include "UpdateFunctionDelegator.h"

/* Wrapper related */
#include "VertexProxy.h"

#include <iostream>
#include <string>
#include <vector>

class XCol : public UpdateFunctionDelegator     {
public:
    
    /* One col of sparse matrix. Support is not stored here, only values.
     */    
    element_type* values;
    unsigned int numElements;
    const unsigned int getNumElements() const {    return numElements;  }
//    std::vector<element_type> values;
//    const unsigned int getNumElements() const {    return values.size();  }
    //std::vector<element_type> row_index;        // support
    
    XCol(unsigned int numElements);
    XCol(const XCol& orig);
    ~XCol();
    
    UpdateFunctionDelegator* clone();
    void accept(VertexVisitor &v, gl::iscope& scope, gl::icallback& schedule);   
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);
    
    int binary_search(std::vector<unsigned int> sorted, int first,
        int last, const unsigned int& lookingFor);
    
private:
    XCol()      {}
};

class XColWrapper : public VertexProxy {
public:

    XColWrapper(unsigned int numElements) {
//        std::cout<<"XColWrapper constructor"<<std::endl;
        this->essence = new XCol(numElements);
    }
    XColWrapper(const XColWrapper& orig)      {
//        std::cout<<"XColWrapper copy constructor"<<std::endl;
        //essence = new XCol();
        this->essence = orig.essence->clone();
    }
    ~XColWrapper() {
        delete(this->essence);
    }

private:
    XColWrapper()       {}
};

#endif	/* XCOL_H */

