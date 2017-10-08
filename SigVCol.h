/* 
 * File:   SigVCol.h
 * Author: ivek
 *
 * Created on January 19, 2012, 3:48 PM
 */

#ifndef SIGVCOL_H
#define	SIGVCOL_H

#include "UpdateFunctionDelegator.h"
#include "VertexProxy.h"
#include "VertexVisitor.h"
//#include "MessagePasser.h"
#include <vector>

class ParentVisitor;
class MessageCollector;

class SigVCol : public UpdateFunctionDelegator
        //,public MessagePasser
{
public:

    SigVCol();
    SigVCol(const SigVCol& orig);
    ~SigVCol();
    UpdateFunctionDelegator* clone();

    /* elements of row this class describes, ordered by increasing indices */
    element_type* values;
    //std::vector<element_type> values;

    /* UpdateFunctionDelegator*/
    void accept(VertexVisitor &v, gl::iscope& scope, gl::icallback& schedule);
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);
    /* UpdateFunctionDelegator end */

    /* MessagePasser - message passing stuff; we are a child of a VCol
     * 
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
/*    void contributionToParent(const VCol * const parent,
            const DiscreteCol * const discrete, element_type* messageAccumulator) const;
*/    /* Message passing related ends*/

private:

};

class SigVColWrapper : public VertexProxy {
public:

    SigVColWrapper() {
        this->essence = new SigVCol();
    }

    SigVColWrapper(const SigVColWrapper& orig) {
        this->essence = orig.essence->clone();
    };

    ~SigVColWrapper() {
        delete(this->essence);
    }
};

#endif	/* SIGVCOL_H */

