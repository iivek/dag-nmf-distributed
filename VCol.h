/* 
 * File:   VCol.h
 * Author: ivek
 *
 * Created on January 17, 2012, 10:29 AM
 */

#ifndef VCOL_H
#define	VCOL_H

#include "UpdateFunctionDelegator.h"
#include "GenericScale.h"
//#include "MessageCollector.h"
//#include "MessagePasser.h"

#include <iostream>
#include <string>
#include <vector>

/* Wrapper related*/
#include "VertexProxy.h"

class SigVCol;
class TRow;
class VHyperRow;

//class VCol :public UpdateFunctionDelegator, public GenericScale
class VCol : public GenericScale
//        ,public MessageCollector,
//        ,public MessagePasser
{
public:
    VCol(bool hasDiscrete, bool hasGammaChild);
    VCol(const VCol& orig);
    ~VCol();
    UpdateFunctionDelegator* clone();

    /* Moments: expectation of gamma and expectation of log gamma */
    /* Expectation of gamma moment is in the GenericScale superclass */
    element_type* expELogValues; /* Order: ascending indices */
//    std::vector<element_type> expELogValues;        /* Order: ascending indices */
    bool hasDiscrete;
    bool hasGammaChild;
    bool isMixture()  { return hasDiscrete; }
    bool hasChild()  { return hasGammaChild; }

    /* UpdateFUnctionDelegator */
    void accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule);
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);
    /* UpdateFunctionDelegator ends */

    /* MessagePasser - message passing stuff, we are parent of a AuxCol, TRow,
     * SigVCol and child of DiscreteCol and AuxCol*/
/*    void acceptAsParent( const MessageCollector& child_, const VCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsParent( const MessageCollector& child_, const SigVCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;    //check
    void acceptAsParent( const MessageCollector& child_, const TRow& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;    //check
    void acceptAsParent( const MessageCollector& child_, const AuxCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;    //check
    void acceptAsParent( const MessageCollector& child_, const DiscreteCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsParent( const MessageCollector& child_, const DirichletCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}   
    
    void acceptAsChild( const MessageCollector& parent_, const VCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsChild( const MessageCollector& parent_, const VHyperCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;    //check
    void acceptAsChild( const MessageCollector& parent_, const ShapeHyperCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;    //check
    void acceptAsChild( const MessageCollector& parent_, const AuxCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;    //check
    void acceptAsChild( const MessageCollector& parent_, const DiscreteCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;    //check
    void acceptAsChild( const MessageCollector& parent_, const DirichletCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsChild( const MessageCollector& parent_, const DirichletHyperCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
*/    /* MessagePasser ends */

    /* Message passing related
     */
/*    void contributionToParent(const AuxCol * const parent,
            const DiscreteCol * const discrete, element_type* messageAccumulator) const;
    void contributionToChild(const SigVCol * const child,
            const DiscreteCol * const discrete, element_type* messageAccumulator) const;
*/    /* Message passing related ends*/        

private:
    VCol();

};

class VColWrapper : public VertexProxy {
public:

    VColWrapper(bool hasDiscrete, bool hasGammaChild) {
        this->essence = new VCol(hasDiscrete, hasGammaChild);
    }

    VColWrapper(const VColWrapper& orig) {
        this->essence = orig.essence->clone();
    };

    ~VColWrapper() {
        delete(this->essence);
    }
    
private:
    VColWrapper()       {}
};

#endif	/* VCOL_H */

