/* 
 * File:   SupportXCol.h
 * Author: ivek
 *
 * Created on January 20, 2012, 10:52 AM
 */

#ifndef SUPPORTXCOL_H
#define	SUPPORTXCOL_H

#include "UpdateFunctionDelegator.h"
#include <vector>
/* Wrapper related */
#include "VertexProxy.h"

class SupportXCol : public UpdateFunctionDelegator      {
public:
    /* One col of sparse matrix, Compressed Col Storage encoded
     */
    unsigned int* support;
    unsigned int numElements;           // number of nonzero elements in this col
    const unsigned int getNumElements() const      {   return numElements;     }
    
//    std::vector<unsigned int> support;
//    const unsigned int getNumElements() const      {   return this->support.size();     }
    
    SupportXCol(unsigned int numElements);
    SupportXCol(const SupportXCol& orig);
    ~SupportXCol();
    
    UpdateFunctionDelegator* clone();
    void accept(VertexVisitor &v, gl::iscope& scope, gl::icallback& schedule);   
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);
    
    unsigned int binarySearch(const unsigned int& rowIndex) const;
    
private:
    SupportXCol()       {}
};

class SupportXColWrapper : public VertexProxy {
public:

    SupportXColWrapper(unsigned int numElements) {
        essence = new SupportXCol(numElements);
    }
    SupportXColWrapper(const SupportXColWrapper& orig)      {
//        essence = new SupportXCol();
        essence = orig.essence->clone();
    }
    ~SupportXColWrapper() {
        delete(essence);
    }
    
private:
    SupportXColWrapper()        {}
};

#endif	/* SUPPORTXCOL_H */

