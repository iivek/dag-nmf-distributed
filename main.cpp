/* 
 * File:   main.cpp
 * Author: ivek
 *
 * Created on January 5, 2012, 5:22 PM
 */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>

#include "graphlab.hpp"
//#include <cblas.h>
#include <cblas.h>
#include <boost/timer.hpp>

#include "math.h"     // debugging- for display only

#include "VertexProxy.h"
#include "VertexVisitor.h"
#include "ChromaticVisitor.h"
#include "ModelNMF.h"

using namespace std;

/*      Shared variables:
 * Dimensions
 * X ... MxN
 * T ... MxK
 * V ... KxN
 */
gl::glshared_const<unsigned int> ModelNMF::M;
gl::glshared_const<unsigned int> ModelNMF::K;
gl::glshared_const<unsigned int> ModelNMF::N;
/*      Static pointers to those variables      */
gl::glshared_const<unsigned int>* VertexProxy::M;
gl::glshared_const<unsigned int>* VertexProxy::K;
gl::glshared_const<unsigned int>* VertexProxy::N;


ChromaticVisitor ModelNMF::chromaticVisitor;

/**
 * TODO: Serialization would be useful to have, at least serialization of
 * the graph structure and VertexWrappers.
 * 
 * @param argc
 * @param argv
 * @return 
 */
int main(int argc, char** argv) {
   
    boost::timer t;        
    gl::core glcore;    // Create the core which contains the graph and engine    
    
    // Parse command line options
    graphlab::command_line_options opts;
    opts.scheduler_type = "chromatic";
    opts.scheduler_opts.add_option("max_iterations", 5);
    opts.scheduler_opts.add_option("color_graph", 0);
    opts.ncpus = 1;
    opts.scope_type = "edge";
    if (!opts.parse(argc, argv)) return EXIT_FAILURE;    
    // Initialize engine with command line options
    glcore.set_engine_options(opts);
    
    ModelNMF* nmf = new ModelNMF();
    nmf->initGraphMkII(glcore.graph(), "./resources/matlab/graphlab_test.mat");
    
//    glcore.graph().save("serialized.bin");
    
    /*
     *           Modifying coloring manually
     * What I aim at here is to have SigVCols, SigTRows of the same color, which
     * will guarantee that they are calculated from the V and T natural
     * parameters of the same generation.
     * 
     * I'll leave V and corresponding structure-related variables as returned
     * by compute_coloring() (it would be nice to have a compute_coloring method
     * for only a part of a graph) - Simple, not too many additional colors, and
     * it will work ok. But some colors may end up having no nodes assigned to
     * them.
     * 
     * If you modify the coloring, make sure that none of the neighboring
     * verteces do not have the same color -otherwise graphlab will automatically
     * color the graph and mess up the ordering of tmp and sigs and all...
     */
    size_t lastColor = 0;
    for (size_t v = nmf->xColOffset; v < nmf->xColOffset + 2*nmf->n; v=v+2) {
        glcore.graph().set_color(v, lastColor);
    }
    ++lastColor;
    for (size_t v = nmf->tmpColOffset; v < nmf->tmpColOffset + nmf->n; ++v) {
        glcore.graph().set_color(v, lastColor);
    }
    ++lastColor;
    for (size_t v = nmf->sigTRowOffset; v < nmf->sigTRowOffset + nmf->m; ++v) {
        glcore.graph().set_color(v, lastColor);
    }
//    ++lastColor;
    // sig_t and sig_v can have the same color
    for (size_t v = nmf->sigVColOffset; v < nmf->sigVColOffset + nmf->n; ++v) {
        glcore.graph().set_color(v, lastColor);
    }
    ++lastColor;
    for (size_t v = nmf->tRowOffset; v < nmf->tRowOffset + nmf->m; ++v) {
            ++lastColor;
        glcore.graph().set_color(v, lastColor);
    }
    ++lastColor;
    // to make the ordering as in MATLAB        
    for (size_t v = 0; v < nmf->vColIDs.size(); ++v) {        
            ++lastColor;
        glcore.graph().set_color(nmf->vColIDs[v], lastColor);
    }
    ++lastColor;
    for (size_t v = 0; v < nmf->auxColIDs.size(); ++v) {        
            ++lastColor;
        glcore.graph().set_color(nmf->auxColIDs[v], lastColor);
    }
    ++lastColor;
    for (size_t v = 0; v < nmf->discreteColIDs.size(); ++v) {     
            ++lastColor;
        glcore.graph().set_color(nmf->discreteColIDs[v], lastColor);
    }
    ++lastColor;
    for (size_t v = 0; v < nmf->dirichletColIDs.size(); ++v) {     
            ++lastColor;
        glcore.graph().set_color(nmf->dirichletColIDs[v], lastColor);
    }
    ++lastColor;  
    //hyperparameters. they can all have the same color
    for (size_t v = nmf->tHyperRowOffset; v < nmf->tHyperRowOffset+nmf->M.get_val(); ++v) {    
            ++lastColor;
        glcore.graph().set_color(v, lastColor);
    }
    ++lastColor;
    for (size_t v = 0; v < nmf->vHyperColIDs.size(); ++v) {        
            ++lastColor;
        glcore.graph().set_color(nmf->vHyperColIDs[v], lastColor);
    }
    ++lastColor;
    glcore.graph().set_color(nmf->auxHyperCol, lastColor);
    ++lastColor;            
    for (size_t v = 0; v < nmf->dirichletHyperColIDs.size(); ++v) {        
            ++lastColor;
        glcore.graph().set_color(nmf->dirichletHyperColIDs[v], lastColor);
    }
    ++lastColor;

    glcore.add_task_to_all(nmf->delegateUpdateFunction, 1.0); // used with chromatic scheduler, it simply sets the default update function            
    
    cout << endl << "Stuff before start() took " << t.elapsed()<<endl;
    
    double runtime = glcore.start();
    
    // output to a . mat
    nmf->matOutput(glcore.graph(), "./resources/matlab/graphlab_test_out.mat");

    return 0;
}