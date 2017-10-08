/* 
 * File:   DirichletHyperCol.h
 * Author: ivek
 *
 * Created on January 31, 2012, 6:25 PM
 */

#ifndef DIRICHLETHYPERCOL_H
#define	DIRICHLETHYPERCOL_H

#include "UpdateFunctionDelegator.h"
//#include "MessagePasser.h"
//#include "MessageCollector.h"
/* Wrapper related*/
#include "VertexProxy.h"

#include <iostream>
#include <string>
#include <vector>


class DirichletHyperCol : public UpdateFunctionDelegator
//        ,public MessagePasser,
//        ,public MessageCollector
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
        
    /* MessagePasser - message passing stuff, we are a parent of a Dirichlet col */
/*    void acceptAsParent( const MessageCollector& child_, const VCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsParent( const MessageCollector& child_, const SigVCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsParent( const MessageCollector& child_, const TRow& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsParent( const MessageCollector& child_, const AuxCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsParent( const MessageCollector& child_, const DiscreteCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsParent( const MessageCollector& child_, const DirichletCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;    //check

    void acceptAsChild( const MessageCollector& parent_, const VCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsChild( const MessageCollector& parent_, const VHyperCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsChild( const MessageCollector& parent_, const ShapeHyperCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsChild( const MessageCollector& parent_, const AuxCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsChild( const MessageCollector& parent_, const DiscreteCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsChild( const MessageCollector& parent_, const DirichletCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsChild( const MessageCollector& parent_, const DirichletHyperCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
*/    /* MessagePasser ends */
    
    /* Message passing related
     */
/*    void contributionToParent( const VCol* const parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;
    void contributionToChild( const DirichletCol* const child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;
*/    /* Message passing related ends*/
    
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

