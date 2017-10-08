/* 
 * File:   ScaleHyperCol.h
 * Author: ivek
 *
 * Created on February 16, 2012, 10:35 AM
 */

#ifndef SCALEHYPERCOL_H
#define	SCALEHYPERCOL_H

#include "UpdateFunctionDelegator.h"
#include "GenericScale.h"

/* Wrapper related*/
#include "VertexProxy.h"


class ScaleHyperCol : public GenericScale
{
public:
    ScaleHyperCol();
    ScaleHyperCol(const ScaleHyperCol& orig);
    ~ScaleHyperCol();
    UpdateFunctionDelegator* clone();

    /* UpdateFunctionDelegator */
    void accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule);
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);
    /* UpdateFunctionDelegator ends */   

};

class ScaleHyperColWrapper : public VertexProxy {
public:

    ScaleHyperColWrapper() {
        this->essence = new ScaleHyperCol();
    }

    ScaleHyperColWrapper(const ScaleHyperColWrapper& orig) {
        this->essence = orig.essence->clone();
    };

    ~ScaleHyperColWrapper() {
        delete(this->essence);
    }
};

#endif	/* SCALEHYPERCOL_H */

