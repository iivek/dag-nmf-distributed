/* 
 * File:   TmpCol.cpp
 * Author: ivek
 * 
 * Created on January 17, 2012, 10:29 AM
 */

#include "TmpCol.h"
#include "TRow.h"
#include "VCol.h"
#include "VertexVisitor.h"

#include "cblas.h"
#include "XCol.h"

#include <iostream>
#include <string>
#include <vector>

TmpCol::TmpCol(unsigned int numElements)
//        :values(numElements)
{
    //    std::cout << "TmpCol constructor" << std::endl;
    this->numElements = numElements;
    this->values = new element_type[this->numElements];
}

TmpCol::TmpCol(const TmpCol& orig)
//         :values(orig.getNumElements())
{
    //    std::cout << "TmpCol copy constructor" << std::endl;   
//    std::copy(orig.values.begin(), orig.values.end(), this->values.begin());   
    this->numElements = orig.numElements;
    this->values = new element_type[this->numElements];
    memcpy(this->values, orig.values, sizeof (element_type) * this->numElements);
}

TmpCol::~TmpCol() {
    //    std::cout << "in TmpCol destructor" << std::endl;
    if (this->numElements > 0) {
        delete(values);
    }
}

UpdateFunctionDelegator* TmpCol::clone() {
    UpdateFunctionDelegator* clone = new TmpCol(*this);
    return (UpdateFunctionDelegator*) clone;
}

void TmpCol::accept(VertexVisitor &v, gl::iscope& scope, gl::icallback& schedule) {
    //        std::cout << "SparseCol accepted" << std::endl;
    v.visit(this, scope, schedule);
}

void TmpCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
//    std::cout << "TmpCol update invoked, " << scope.color() << std::endl;
    /* Get a reference to the vertex data and visible edges */
    //    vertex_data& us = scope.vertex_data();    
    gl::edge_list in_edges = scope.in_edge_ids();

    /* Accessing neigboring XCol object - at numElements+1
     */
    const element_type* rawXcol;
    {
        size_t sourcev = scope.source(in_edges[ getNumElements() + 1 ]);
        const VertexProxy& neighbor = scope.neighbor_vertex_data(sourcev);
        const XCol& delegator = (const XCol&) *neighbor.getDelegator();
        rawXcol = &delegator.values[0];
    }
    
    /* Collecting neighboring TRow objects - from 1 to this->numElements
     */
    element_type left[this->getNumElements() * VertexProxy::K->get_val()];
    element_type* current = left;
    for (size_t i = 1; i < this->getNumElements() + 1; ++i) {
        size_t sourcev = scope.source(in_edges[i]);
        const VertexProxy& neighbor = scope.neighbor_vertex_data(sourcev);        
        const TRow& delegator = (const TRow&) *neighbor.getDelegator();
        memcpy(current, &delegator.expELogValues[0], sizeof(element_type) * VertexProxy::K->get_val());
        current += VertexProxy::K->get_val();        
    }

    /* Accessing neigboring VCol object - at 0
     */
    const element_type* right;
    {
        size_t sourcev = scope.source(in_edges[0]);
        const VertexProxy& neighbor = scope.neighbor_vertex_data(sourcev);
        const VCol& delegator = (const VCol&) *neighbor.getDelegator();
        right = &delegator.expELogValues[0];
//        std::cout<<"here VCOl    "<<delegator.scaleRelatedValues[0]<<std::endl;
    }

    /* Lt*Lv */
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
            this->getNumElements(), VertexProxy::K->get_val(),
            1.0, left, VertexProxy::K->get_val(),
            right, 1.0, 0.0, &this->values[0], 1.0);
    

    /* X./(Lt*Lv) */
    /* TODO: some kind of SIMD instruction instead of looping */
    //element_type* resultIter = this->values;
    for (int i = 0; i<this->getNumElements(); ++i) {
        this->values[i] = *rawXcol++ / this->values[i];
    }
//    std::cout<<std::endl;

}

