/* 
 * UpdateFunctionDelegator, as a visitor pattern backbone, has some abstract
 * methods. However, graphlab instantiates its vertices, which means
 * that vertex class cannot be abstract, so UpdateFunctionDelegator cannot be a
 * vertex class. As a remedy, VertexProxy's role is that of a non-abstract
 * mediator between graphlab and UpdateFunctionDelegator.
 * 
 * File:   VertexProxy.h
 * Author: ivek
 *
 * Created on January 9, 2012, 4:05 PM
 */

#ifndef VERTEXPROXY_H
#define	VERTEXPROXY_H

#include "UpdateFunctionDelegator.h"
#include "graphlab.hpp"

#include <graphlab/serialization/serialization_includes.hpp>

class VertexProxy {
public:
    
    /* Pointers to shared variables        */
    static gl::glshared_const<unsigned int>* M;
    static gl::glshared_const<unsigned int>* K;
    static gl::glshared_const<unsigned int>* N;
    
    VertexProxy(){
    };
    
    VertexProxy(const VertexProxy& orig)   {
        this->essence = orig.essence->clone();
    }
    
    ~VertexProxy(){
    };

   UpdateFunctionDelegator* getDelegator()    {
        return essence;
    }
        
    const UpdateFunctionDelegator* const getDelegator() const  {
        return essence;
    }       
    
    static void initSharedData(
        gl::glshared_const<unsigned int>* m, 
        gl::glshared_const<unsigned int>* n,
        gl::glshared_const<unsigned int>* k)        {
        
        M=m;
        N=n;
        K=k;
    }
    
    /* Serialization - only keep the wrapper, */
    void save(graphlab::oarchive& oarc) const   {};
    void load(graphlab::iarchive& iarc) {};

protected:
    UpdateFunctionDelegator* essence;
};

#endif	/* VERTEXPROXY_H */

