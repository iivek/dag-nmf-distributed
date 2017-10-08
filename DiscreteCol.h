/* 
 * File:   DiscreteCol.h
 * Author: ivek
 *
 * Created on January 31, 2012, 5:50 PM
 */

#ifndef DISCRETECOL_H
#define	DISCRETECOL_H

#include "UpdateFunctionDelegator.h"
//#include "MessagePasser.h"
//#include "MessageCollector.h"
/* Wrapper related*/
#include "VertexProxy.h"

#include <iostream>
#include <string>
#include <vector>


class DiscreteCol : public UpdateFunctionDelegator
//        ,public MessagePasser
//        ,public MessageCollector
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
        
    /* MessagePasser - message passing stuff, we are parent of VCol and child
     * of DirichletCol
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
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;    //check
    void acceptAsChild( const MessageCollector& parent_, const DirichletHyperCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
*/    /* MessagePasser ends */
    
    /* Message passing related
     */
/*    void contributionToParent( const DirichletCol* const parent,
        const DiscreteCol* const discrete,  element_type* messageAccumulator ) const;
*/    /* Message passing related ends*/
    
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

