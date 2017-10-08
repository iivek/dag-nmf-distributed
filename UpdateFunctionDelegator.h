/* Base class -inherit from it and wrap it into vertexProxy. Then pass it to 
 * graphlab's add_vertex
 * 
 * 
 * 
 * 
 *  * File:   UpdateFunctionDelegator.h
 * Author: ivek
 *
 * Created on January 9, 2012, 4:46 PM
 */

#ifndef UPDATEFUNCTIONDELEGATOR_H
#define	UPDATEFUNCTIONDELEGATOR_H

#include "graphlab.hpp"

class VertexProxy;
typedef char edge_data;
typedef graphlab::graph<VertexProxy, edge_data> graph_type;
typedef graphlab::types<graph_type> gl;
typedef double element_type;    /* Datatype of elements in our matrices*/

class VertexVisitor;

class UpdateFunctionDelegator {
public:        

    UpdateFunctionDelegator() {}        
//    UpdateFunctionDelegator(const UpdateFunctionDelegator& orig);
    ~UpdateFunctionDelegator() {        }

    virtual void accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule) = 0;
    
    /* If actions in the updateFunction depent on the classtype of vertices in
     * scope, a visitor pattern (similar to the one in
     * UpdateFunctionDelegator/VertexVisitor) can be added to support this type
     * of branching.
     * Alternatively, graphlab's guarantee of returning vertices in the same
     * order they have been added can be used to determine class of the vertex
     * in the scope.
     */ 
    virtual void updateFunction(gl::iscope& scope, gl::icallback& scheduler) = 0;
    virtual UpdateFunctionDelegator* clone() = 0;
private:
};

#endif	/* UPDATEFUNCTIONDELEGATOR_H */

