/* 
 * 
 * File:   ModelNMF.h
 * Author: ivek
 *
 * Created on January 5, 2012, 11:09 AM
 */

#ifndef MODELNMF_H
#define	MODELNMF_H


// TODO: try to move as much of the headers to .cpp
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
//#include "cblas.h"
#include <cblas.h>
#include <string>
#include <vector>

#include "graphlab.hpp"

#include "ChromaticVisitor.h"

#include "TRow.h"
#include "VCol.h"
#include "TmpCol.h"
#include "XCol.h"
#include "SupportXCol.h"
#include "SigTRow.h"
#include "SigVCol.h"
#include "THyperRow.h"
#include "ShapeHyperCol.h"
#include "ScaleHyperCol.h"
#include "AuxCol.h"
#include "DiscreteCol.h"
#include "DirichletCol.h"
#include "DirichletHyperCol.h"

class ModelNMF {
    

public:
    ModelNMF();
    ~ModelNMF();
    //ModelNMF(const ModelNMF& orig);       // no copy constructor
    
    /*  Here we create the shared data values
    */
    void initSharedData(unsigned int m, unsigned int n, unsigned int p);
    
    /* Throughout the code, we rely heavily on the order of vertices.
     */
    void initGraphMkII(graph_type& g, const char* const filename);
    void matOutput(graph_type& g, const char* const filename);
    
    /* Simply delegates the call to the vertex member method. Signature doesn't
     * match unless it's declared static (that is, i wasn't able to make it work
     * that way)
     */ 
    static void delegateUpdateFunction(gl::iscope& scope, gl::icallback& scheduler) ;
    
    /* helpers
     */  
    long int binary_search(unsigned int* sorted, int first, int last, const unsigned int& lookingFor);
            
    // TODO: make the following ones obsolete - only the graphlab's shared
    // variables should stay
    unsigned int m;
    unsigned int k;
    unsigned int n;
    
    static ChromaticVisitor chromaticVisitor;
    
    /* Graphlab's shared variables begin */
    static gl::glshared_const<unsigned int> M;
    static gl::glshared_const<unsigned int> K;
    static gl::glshared_const<unsigned int> N;
    /* Graphlab's shared variables end */
    
    gl::iscope* scope;
    gl::core glcore;
    
    /* Vertex group offsets, for easier vertex access from outside graphlab */    
    unsigned int xColOffset;   
    unsigned int tmpColOffset;
    unsigned int sigTRowOffset;
    unsigned int sigVColOffset;
    unsigned int tRowOffset;
    unsigned int tHyperRowOffset;
    unsigned int auxHyperCol;   // Aux hyperparameter ID
    /* And lists of IDs of (VCols together with AuxCols )and (VHyperCols) */
    
    std::vector<unsigned int> vColIDs;
    std::vector<unsigned int> auxColIDs;
    std::vector<unsigned int> vHyperColIDs;
    std::vector<unsigned int> shapeHyperColIDs;
    std::vector<unsigned int> discreteColIDs;
    std::vector<unsigned int> dirichletColIDs;
    std::vector<unsigned int> dirichletHyperColIDs;
    
    std::vector<unsigned int> forwardMapping;
    std::vector<unsigned int> numParents;
    unsigned int noparents;
    
private:
       
};

#endif	/* MODELNMF_H */

