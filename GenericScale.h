/* 
 * File:   GenericScale.h
 * Author: ivek
 *
 * Created on February 16, 2012, 10:48 AM
 */

#ifndef GENERICSCALE_H
#define	GENERICSCALE_H

#include "UpdateFunctionDelegator.h"
#include "VertexProxy.h"
#include <vector>

/**
 * An interface through which vertices pass the scale hyperparameter to gamma
 * vertices;
 * 
 * @return 
 */
class GenericScale :public UpdateFunctionDelegator     {
public:
    element_type* scaleRelatedValues;
    
GenericScale()
{    
    this->scaleRelatedValues = new element_type[VertexProxy::K->get_val()];
}

GenericScale(const GenericScale& orig)
{   
//     std::copy(orig.scaleRelatedValues.begin(), orig.scaleRelatedValues.end(), this->scaleRelatedValues.begin());
     
    this->scaleRelatedValues = new element_type[VertexProxy::K->get_val()];
    memcpy(this->scaleRelatedValues, orig.scaleRelatedValues, sizeof(element_type)*VertexProxy::K->get_val());
};

~GenericScale() {
    delete(this->scaleRelatedValues);
}

};

#endif	/* GENERICSCALE_H */