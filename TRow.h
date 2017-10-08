/* 
 * File:   TRow.h
 * Author: ivek
 *
 * Created on January 17, 2012, 10:29 AM
 */

#ifndef TROW_H
#define	TROW_H

#include "UpdateFunctionDelegator.h"
//#include "MessagePasser.h"
#include <vector>

/* Wrapper related*/
#include "VertexProxy.h"


class MessageCollector;

class TRow : public UpdateFunctionDelegator
//        ,public MessagePasser
{
public:            
    TRow(unsigned int row);
    TRow(const TRow& orig);    
    ~TRow();
    UpdateFunctionDelegator* clone();
    /* our row index*/
    unsigned int row;
    /* Moments: expectation of gamma and expectation of log gamma */
    element_type* eValues;      /* Order: ascending indices */
    element_type* expELogValues;
//    std::vector<element_type> eValues;      /* Order: ascending indices */
//    std::vector<element_type> expELogValues;

    /* UpdateFunctionDelegator interface */
    void accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule);
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);
    /* UpdateFunctionDelegator ends */
    
    /* Message passing stuff - MessagePasser members
     * - graphlab's vertices defined in the implementation don't map one-to-one
     * to stochastic variables (see e.g. SigVCol, SigTRow) as a result,
     * messages from Poissonians will be collected from both TRows and SigVCols
     */   
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
*/    /* Message passing related ends*/

private:
    TRow()      {}
};

class TRowWrapper : public VertexProxy {
public:    

    TRowWrapper(unsigned int row)       {        
        this->essence = new TRow(row);
    }
    TRowWrapper(const TRowWrapper& orig)      {        
        this->essence = orig.essence->clone();
    };
    ~TRowWrapper() {
        delete(this->essence);
    }
    
private:
    TRowWrapper()       {}
};

#endif	/* TROW_H */

