/* 
 * File:   TmpCol.h
 * Author: ivek
 *
 * Created on January 17, 2012, 10:29 AM
 */

#ifndef TMPCOL_H
#define	TMPCOL_H

#include "UpdateFunctionDelegator.h"
#include <vector>
/* Wrapper related */
#include "VertexProxy.h"


class TmpCol : public UpdateFunctionDelegator   {
public:
    /* One col of sparse matrix, Compressed Col Storage encoded
     */    
    element_type* values;
    unsigned int numElements;           // number of nonzero elements in this col
    const unsigned int getNumElements() const      {       return numElements;    }           // number of nonzero elements in this col
//    std::vector<element_type> values;
//    const unsigned int getNumElements() const      {       return values.size();    }           // number of nonzero elements in this col
    
    /* Support is stored in corresponding SupportXCol object. numElements is
     * actually redundant, (it's the same as in corresponding SupportXCol object)
     * but the design requires it (needed by the copy constructor)
     */
    
    //std::vector<element_type> values;
    //std::vector<element_type> row_index;        // support
    
    TmpCol(unsigned int numElements);
    TmpCol(const TmpCol& orig);
    ~TmpCol();
    
    UpdateFunctionDelegator* clone();
    void accept(VertexVisitor &v, gl::iscope& scope, gl::icallback& schedule);   
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);
    
private:
    TmpCol()    {}
};

class TmpColWrapper : public VertexProxy {
public:

    TmpColWrapper(unsigned int numElements) {
        essence = new TmpCol(numElements);
    }
    TmpColWrapper(const TmpColWrapper& orig)      {
//        essence = new TmpCol();
        essence = orig.essence->clone();
    }
    ~TmpColWrapper() {
        delete(essence);
    }
    
private:
    TmpColWrapper()     {}
};

#endif	/* TMPCOL_H */

