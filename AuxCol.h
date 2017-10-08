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
//#include "MessagePasser.h"
//#include "MessageCollector.h"
/* Wrapper related*/
#include "VertexProxy.h"

#include <iostream>
#include <string>
#include <vector>

//#include <graphlab/serialization/serialization_includes.hpp>


//class AuxCol : public UpdateFunctionDelegator, public GenericScale
class AuxCol : public GenericScale
//        ,public MessagePasser,
//        ,public MessageCollector
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
//    std::vector<element_type> eLogValues;       /* Order: ascending indices */    
//    std::vector<unsigned int> indexInMixture;
//    std::vector<bool> childIsMixture;
//    const unsigned int numChildren() const {    return indexInMixture.size();   }
    
    /* Graphlab's serialization */
//    void save(graphlab::oarchive& oarc) const;    
//    void load(graphlab::iarchive& iarc);
    
    /* Alternative desing: break DiscreteCol into graphlab vertices, as many as
     * there are components in the mixture. Then AuxCol (us) would have vertices
     * only to one moment in the discrete distribution (which corresponds to it)
     * and above arrays would be obsolete in present form. Pros: we have access
     * not to entire discrete distro but only to the moment which we need - 
     * more precise memory mapping
     * In that case such broken down DiscreteCols should have vertices between
     * them to have access to all the moments when updating the complete
     * Discrete distribution.
     * A TODO task perhaps?
     */        
    
    /* UpdateFUnctionDelegator */
    void accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule);
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);   
    /* UpdateFUnctionDelegator ends */
        
    /* MessagePasser - message passing stuff, we are a parent/child of VCol */
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
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;    //check
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
    void contributionToChild( const VCol* const child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;
*/    /* Message passing related ends*/
    
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

