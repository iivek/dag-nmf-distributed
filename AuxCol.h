/* Vertex that describes the auxillary gamma variable
 * 
 * File:   AuxCol.h
 * Author: ivek
 *
 * Created on January 26, 2012, 4:59 PM
 */

#ifndef AUXCOL_H
#define	AUXCOL_H

#include "UpdateFunctionDelegator.h"
#include "GenericScale.h"
/* Wrapper related*/
#include "VertexProxy.h"

#include <iostream>
#include <string>
#include <vector>

class AuxCol : public GenericScale
{
    
public:            
    AuxCol(unsigned int numChildren);
    AuxCol(const AuxCol& orig);    
    ~AuxCol();
    UpdateFunctionDelegator* clone();
    
    /* Expectation of gamma moment is in the GenericScale superclass */
    /* Expectation of log gamma moment*/
    element_type* eLogValues;           /* Order: ascending indices */    
    unsigned int numChildren;           /* we'll have at least one child */    
    unsigned int* indexInMixture;
    bool* childIsMixture;
    const unsigned int getNumChildren() const {    return numChildren;   }
    
    /* Arrays to vectors: easier serialization but more memory used. In general. */    
    /* Graphlab's serialization */
//    void save(graphlab::oarchive& oarc) const;    
//    void load(graphlab::iarchive& iarc);
    
    /* UpdateFUnctionDelegator */
    void accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule);
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);   
    /* UpdateFUnctionDelegator ends */
    
private:
    AuxCol()    {}

};

class AuxColWrapper : public VertexProxy {
public:    
    AuxColWrapper(unsigned int numChildren){        
        this->essence = new AuxCol(numChildren);
    }
    AuxColWrapper(const AuxColWrapper& orig)      {        
        this->essence = orig.essence->clone();
    };
    ~AuxColWrapper() {
        delete(this->essence);
    }
    
private:
    AuxColWrapper()     {}
};

#endif	/* AUXCOL_H */

