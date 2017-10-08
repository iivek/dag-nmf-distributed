/* 
 * File:   DirichletCol.h
 * Author: ivek
 *
 * Created on January 31, 2012, 5:50 PM
 */

#ifndef DIRICHLETCOL_H
#define	DIRICHLETCOL_H

#include "UpdateFunctionDelegator.h"
/* Wrapper related*/
#include "VertexProxy.h"

#include <iostream>
#include <string>
#include <vector>

#include <graphlab/serialization/serialization_includes.hpp>



class DirichletCol : public UpdateFunctionDelegator
{
    
public:            
    DirichletCol(unsigned int order);
    DirichletCol(const DirichletCol& orig);    
    ~DirichletCol();
    UpdateFunctionDelegator* clone();
    
    /* Moments: */
    unsigned int order;
    element_type* values;
    const unsigned int getOrder() const {  return order;   }

/*    void save(oarchive& oarc) const     {
        oarc << eLogValues << indexInMixture;
    }
    
    void load(iarchive& iarc)   {
        iarc >> eLogValues >> indexInMixture;
    }
*/
    
    /* UpdateFUnctionDelegator */
    void accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule);
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);   
    /* UpdateFUnctionDelegator ends */
   
private:
    DirichletCol()      {}

};

class DirichletColWrapper : public VertexProxy {
public:    

    DirichletColWrapper(unsigned int order){        
        this->essence = new DirichletCol(order);
    }
    DirichletColWrapper(const DirichletColWrapper& orig)      {        
        this->essence = orig.essence->clone();
    };
    ~DirichletColWrapper() {
        delete(this->essence);
    }
private:
    DirichletColWrapper()       {}
};

#endif	/* DIRICHLETCOL_H */

