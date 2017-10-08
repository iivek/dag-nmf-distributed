/* 
 * File:   DirichletHyperCol.h
 * Author: ivek
 *
 * Created on January 31, 2012, 6:25 PM
 */

#ifndef DIRICHLETHYPERCOL_H
#define	DIRICHLETHYPERCOL_H

#include "UpdateFunctionDelegator.h"
/* Wrapper related*/
#include "VertexProxy.h"

#include <iostream>
#include <string>
#include <vector>


class DirichletHyperCol : public UpdateFunctionDelegator
{
    
public:            
    DirichletHyperCol(unsigned int order);
    DirichletHyperCol(const DirichletHyperCol& orig);    
    ~DirichletHyperCol();
    UpdateFunctionDelegator* clone();
    
    /* Moments: */
    double* values;
    unsigned int order;
    const unsigned int getOrder() const    {       return order;     }
    
//    std::vector<double> values;
//    const unsigned int getOrder() const    {       return values.size()/VertexProxy::K->get();     }
    
    /* UpdateFUnctionDelegator */
    void accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule);
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);   
    /* UpdateFUnctionDelegator ends */
    
private:
    DirichletHyperCol();
};

class DirichletHyperColWrapper : public VertexProxy {
public:    

    DirichletHyperColWrapper(unsigned int order){        
        this->essence = new DirichletHyperCol(order);
    }
    DirichletHyperColWrapper(const DirichletHyperColWrapper& orig)      {        
        this->essence = orig.essence->clone();
    };
    ~DirichletHyperColWrapper() {
        delete(this->essence);
    }
    
private:
    DirichletHyperColWrapper()  {}
};

#endif	/* DIRICHLETHYPERCOL_H */

