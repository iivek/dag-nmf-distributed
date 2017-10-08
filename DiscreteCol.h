/* 
 * File:   DiscreteCol.h
 * Author: ivek
 *
 * Created on January 31, 2012, 5:50 PM
 */

#ifndef DISCRETECOL_H
#define	DISCRETECOL_H

#include "UpdateFunctionDelegator.h"
/* Wrapper related*/
#include "VertexProxy.h"

#include <iostream>
#include <string>
#include <vector>


class DiscreteCol : public UpdateFunctionDelegator
{
    
public:            
    DiscreteCol(unsigned int mixtureSize);
    DiscreteCol(const DiscreteCol& orig);    
    ~DiscreteCol();
    UpdateFunctionDelegator* clone();
    
    element_type* values;     // organised per element of mixtures, accross rows
    unsigned int mixtureSize;
    const unsigned int getMixtureSize() const {  return mixtureSize;   }
    
//    std::vector<element_type> values;       /* Order: ascending indices */    
//    const unsigned int getMixtureSize() const {  return values.size()/VertexProxy::K->get_val();   }
    
    /* UpdateFUnctionDelegator */
    void accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule);
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);   
    /* UpdateFUnctionDelegator ends */
        
private:
    DiscreteCol()       {}

};

class DiscreteColWrapper : public VertexProxy {
public:    

    DiscreteColWrapper(unsigned int mixtureSize){        
        this->essence = new DiscreteCol(mixtureSize);
    }
    DiscreteColWrapper(const DiscreteColWrapper& orig)      {        
        this->essence = orig.essence->clone();
    };
    ~DiscreteColWrapper() {
        delete(this->essence);
    }
private:
    DiscreteColWrapper()        {}
};

#endif	/* DISCRETECOL_H */

