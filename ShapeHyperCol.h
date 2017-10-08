/* 
 * File:   ShapeHyperCol.h
 * Author: ivek
 *
 * Created on January 24, 2012, 9:23 AM
 */

#ifndef SHAPEHYPERCOL_H
#define	SHAPEHYPERCOL_H

#include "UpdateFunctionDelegator.h"
//#include "MessagePasser.h"
#include <vector>

/* Wrapper related*/
#include "VertexProxy.h"


class MessageCollector;

class ShapeHyperCol : public UpdateFunctionDelegator
//        ,public MessagePasser
{
public:
    ShapeHyperCol();
    ShapeHyperCol(const ShapeHyperCol& orig);
    ~ShapeHyperCol();
    UpdateFunctionDelegator* clone();

    element_type* aValues; /* Order: ascending indices */
    //std::vector<element_type> aValues; /* Order: ascending indices */

    /* UpdateFunctionDelegator */
    void accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule);
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);
    /* UpdateFunctionDelegator ends */

    /* MessagePasser - message passing stuff, we are a parent of a VCol/AuxCol */
/*    void acceptAsParent( const MessageCollector& child_, const VCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;    //check
    void acceptAsParent( const MessageCollector& child_, const SigVCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsParent( const MessageCollector& child_, const TRow& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
    void acceptAsParent( const MessageCollector& child_, const AuxCol& child,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const;    //check
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
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {}
*/    /* MessagePasser ends */

    /* Message passing related
     */
/*    void contributionToParent(const VCol * const parent,
            const DiscreteCol * const discrete, element_type* messageAccumulator) const;
    void contributionToChild(const VCol * const child,
            const DiscreteCol * const discrete, element_type* messageAccumulator) const;
    void contributionToChild(const AuxCol * const child,
            const DiscreteCol * const discrete, element_type* messageAccumulator) const;
*/

    /* Message passing related ends*/

};

class ShapeHyperColWrapper : public VertexProxy {
public:

    ShapeHyperColWrapper() {
        this->essence = new ShapeHyperCol();
    }

    ShapeHyperColWrapper(const ShapeHyperColWrapper& orig) {
        this->essence = orig.essence->clone();
    };

    ~ShapeHyperColWrapper() {
        delete(this->essence);
    }
};

#endif	/* SHAPEHYPERCOL_H */

