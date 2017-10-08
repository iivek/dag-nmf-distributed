/* 
 * File:   ShapeHyperCol.h
 * Author: ivek
 *
 * Created on January 24, 2012, 9:23 AM
 */

#ifndef SHAPEHYPERCOL_H
#define	SHAPEHYPERCOL_H

#include "UpdateFunctionDelegator.h"
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

