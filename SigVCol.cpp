/* 
 * File:   SigVCol.cpp
 * Author: ivek
 * 
 * Created on January 19, 2012, 3:48 PM
 */

#include "SigVCol.h"
#include "TRow.h"
#include "TmpCol.h"
#include "VCol.h"
#include "VertexVisitor.h"

//#include "MessageCollector.h"

//#include <cblas.h>
#include <cblas.h>

#include <iostream>
#include <string>
#include <vector>


using namespace std;

SigVCol::SigVCol()
//        : values( VertexProxy::K->get_val() )
{
    //    std::cout << "in SigVCol constructor" << std::endl;
    this->values = new element_type[VertexProxy::K->get_val()];
    memset(this->values, 0, sizeof (element_type) * VertexProxy::K->get_val());
}

SigVCol::SigVCol(const SigVCol& orig)
//        : values( VertexProxy::K->get_val() )
{
    //    std::cout << "in SigVCol copy constructor" << std::endl;
    
//    std::copy(orig.values.begin(), orig.values.end(), this->values.begin());
    this->values = new element_type[VertexProxy::K->get_val()];
    memcpy(this->values, orig.values, sizeof (element_type) * VertexProxy::K->get_val());
}

SigVCol::~SigVCol() {
    delete(this->values);
}

UpdateFunctionDelegator* SigVCol::clone() {
    UpdateFunctionDelegator* clone = new SigVCol(*this);
    return (UpdateFunctionDelegator*) clone;

}

void SigVCol::accept(VertexVisitor &v, gl::iscope& scope, gl::icallback& schedule) {
    v.visit(this, scope, schedule);
}

void SigVCol::updateFunction(gl::iscope& scope, gl::icallback& scheduler) {
    std::cout << "SigVCol update invoked, " << scope.color() << std::endl;

    /* Get a reference to the vertex data and visible edges */
    //    vertex_data& us = scope.vertex_data();    

    gl::edge_list in_edges = scope.in_edge_ids();

    /*
     * The order in which vertices are fetched is important here - as in
     * ModelNMF's initgraph()
     * 
     * one VCol object followed by
     * (in_edges.size()-2) TRow objects followed by     
     * one TmpCol object
     */
    unsigned int numPasses = in_edges.size() - 2; // #neigboring LtRows == sparsity of tmpCol neighbor
    unsigned int vColOffset = 0;
    unsigned int tRowOffset = vColOffset + 1;
    unsigned int tmpColOffset = tRowOffset + numPasses;

    // buffer to be filled with neigboring LtRows, rowwise, (corresponds to CblasRowMajor)
    element_type buffer[VertexProxy::K->get() * numPasses];
    element_type* location = buffer;
    /* Collecting TRows
     */
    for (size_t i = tRowOffset; i < tRowOffset + numPasses; ++i) {
        const VertexProxy& neighbor = scope.neighbor_vertex_data(scope.source(in_edges[ i ]));
        const TRow& delegator = (const TRow&) *neighbor.getDelegator();
        memcpy(location, &delegator.expELogValues[0], sizeof (element_type) * VertexProxy::K->get());
        location += VertexProxy::K->get();
    }
    /* Fetching TmpCol neighbor...
     */
    {
        const VertexProxy& neighbor = scope.neighbor_vertex_data(scope.source(in_edges[ tmpColOffset ]));
        const TmpCol& delegator = (const TmpCol&) *neighbor.getDelegator();
        /* ... and performing matrix multiplication
         */
        cblas_dgemv(CblasRowMajor, CblasTrans, numPasses, VertexProxy::K->get(), 1.0,
                buffer, VertexProxy::K->get(), &delegator.values[0], 1, 0.0, &this->values[0], 1);
    }
    /* Fetching VCol...
     * ( TODO: use SIMD for elementwise multiplication )
     */
    {
        const VertexProxy& neighbor = scope.neighbor_vertex_data(scope.source(in_edges[ vColOffset ]));
        const VCol& delegator = (const VCol&) *neighbor.getDelegator();
        /* ... and performing elementwise product with what we got previously
         */
        for (unsigned int row = 0; row < VertexProxy::K->get_val(); ++row) {
            this->values[row] *= delegator.expELogValues[row];
        }
    }

    //   scheduler.add_task(gl::update_task(scope.vertex(), update_function), 1.0);
}

/*
void SigVCol::contributionToParent(const VCol* const parent,
        const DiscreteCol * const discrete, element_type* messageAccumulator) const {
//    std::cout<<"SigVCol::contributionToParent"<<std::endl;
    element_type* second = &messageAccumulator[VertexProxy::K->get_val()];
    for (int row = 0; row < VertexProxy::K->get_val(); ++row) {
        second[row] += values[row];
    }     
}

void SigVCol::acceptAsChild( const MessageCollector& parent_, const VCol& parent,
        const DiscreteCol* const discrete, element_type* messageAccumulator ) const     {
    std::cout<<"SigVCol::acceptAsChild"<<std::endl;
    parent_.visitChild(this, &parent, discrete, messageAccumulator);
    
}
*/