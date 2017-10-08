/* 
 * File:   DirichletCol.h
 * Author: ivek
 *
 * Created on January 31, 2012, 5:50 PM
 */

#ifndef DIRICHLETCOL_H
#define	DIRICHLETCOL_H

#include "UpdateFunctionDelegator.h"
//#include "MessagePasser.h"
//#include "MessageCollector.h"
/* Wrapper related*/
#include "VertexProxy.h"

#include <iostream>
#include <string>
#include <vector>

#include <graphlab/serialization/serialization_includes.hpp>



class DirichletCol : public UpdateFunctionDelegator
//        ,public MessagePasser,
//        ,public MessageCollector
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
    
//    std::vector<element_type> values;       /* Order: ascending indices */    
//    const unsigned int getOrder() const {  return values.size()/VertexProxy::K->get();   }
    

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
        
    /* MessagePasser - message passing stuff, we are parent of VCol and child
     * of DirichletHyperCol
     */
/*    void acceptAsParent( const MessageCollector& child_, const VCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;    //check
    void acceptAsParent( const MessageCollector& child_, const SigVCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsParent( const MessageCollector& child_, const TRow& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsParent( const MessageCollector& child_, const AuxCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsParent( const MessageCollector& child_, const DiscreteCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsParent( const MessageCollector& child_, const DirichletCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}

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
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;    //check
*/    /* MessagePasser ends */
    
    /* Message passing related
     */
/*    void contributionToParent( const DirichletHyperCol* const parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;
    void contributionToChild( const DiscreteCol* const child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;
*/    /* Message passing related ends*/
    
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

