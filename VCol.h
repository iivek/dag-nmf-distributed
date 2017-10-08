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
#include <iostream>
#include <string>
#include <vector>

/* Wrapper related*/
#include "VertexProxy.h"

class SigVCol;
class TRow;
class VHyperRow;

class VCol : public GenericScale
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

