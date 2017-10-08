/* 
 * File:   ModelNMF.cpp
 * Author: ivek
 * 
 * Created on January 5, 2012, 11:09 AM
 */

#include "ModelNMF.h"

#include "math.h"
#include <mat.h>
#include <map>

using namespace std;


/* Something should be said on the input parameters, even though the example
 * should be self-explanatory.
 * Gamma pdfs are paraeterized with shapes and inverse *means*. 
 *
 */
ModelNMF::ModelNMF() {    
}

ModelNMF::~ModelNMF() {
}

void ModelNMF::initSharedData(unsigned int m, unsigned int n, unsigned int k) {

    M.set(m);
    N.set(n);
    K.set(k);
}

void ModelNMF::delegateUpdateFunction(
        gl::iscope& scope, gl::icallback& scheduler) {

    //const VertexProxy& current = scope.vertex_data();
    //((UpdateFunctionDelegator*) current.getDelegator())->accept(chromaticVisitor, scope, scheduler);
    VertexProxy& current = scope.vertex_data();
    ((UpdateFunctionDelegator*) current.getDelegator())->accept(chromaticVisitor, scope, scheduler);
}// update function ends

long int ModelNMF::binary_search(unsigned int* sorted, int first,
        int last, const unsigned int& lookingFor) {
    unsigned int mid;       

    while (first <= last) {
        mid = (first + last) / 2; //  split space in half
        if (lookingFor > sorted[mid]) {
            first = mid + 1; // proceed in the top half
        } else if (lookingFor < sorted[mid]) {
            last = mid - 1; // proceed search in bottom half
        } else {
            return mid;
        }
    }
    return -1; // not found
}

/*
 * TODO: handle empty rows/cols of X -> will crash if those exist, i guess
 * Usage of datatypes is inconsistent- sometimes it is double, sometimes
 * element_type. 
 * Input matrix needs to be reordered to match topologically sorted cols of right matrix
 */
void ModelNMF::initGraphMkII(graph_type& g, const char* const filename)     {       
    
    MATFile* pmat = matOpen(filename, "r");
    assert(pmat != NULL);
    /* Reminder: mxarrays are column major
     */    

    /* Note: all indexing we process here should be 0-based instead of Matlab
     * (Octave) styled 1-based
     */
    /* Processing the topologically sorted, already populated with auxillary
     * gammas, adjacency matrix over columns of the right, V matrix.
     */
    mxArray* w_mx = matGetVariable(pmat, "W");
    double* pM = (double*) mxGetData(w_mx);
    this->m = (unsigned int)* pM;
    mxArray* k_mx = matGetVariable(pmat, "K");
    double* pN = (double*) mxGetData(k_mx);
    this->n = (unsigned int)* pN;
    mxArray* i_mx = matGetVariable(pmat, "I");
    double* pI = (double*) mxGetData(i_mx);
    this->k = (unsigned int)* pI;
    std::cout << "m: " << this->m << " n: " << this->n <<" k: " << this->k <<std::endl;
    mxDestroyArray(w_mx);
    mxDestroyArray(k_mx);
    mxDestroyArray(i_mx);
    
    this->initSharedData(m, n, k);
    VertexProxy::initSharedData(&(this->M), &(this->N), &(this->K));        
    
    std::cout << "Checkpoint -3" <<  std::endl; //Debug
    
    mxArray* rowSupport_mx = matGetVariable(pmat, "rows_adjacency_colwise");
    unsigned int adjacencySparsity = mxGetM(rowSupport_mx); // number of nonzero entries    
    unsigned int* rowSupport = (unsigned int*) mxGetData(rowSupport_mx);    
    mxArray* colSupport_mx = matGetVariable(pmat, "cols_adjacency_colwise");
    
//    std::cout << "Checkpoint -2" <<  std::endl; //Debug

//    std::cout << "sparsity = " << adjacencySparsity << std::endl; //DEBUG

    assert(mxGetM(colSupport_mx) == adjacencySparsity);
    unsigned int* colSupport = (unsigned int*) mxGetData(colSupport_mx);
    /* Mapping of indices to get the original ordering of vertices (before
     * topological sorting)
     */
    /*
        for (unsigned int i = 0; i < adjacencySparsity; ++i) {
            std::cout << rowSupport[i] << " " << colSupport[i] << std::endl; //DEBUG        
        }
        std::cout <<
                std::endl; //DEBUG
     */
            
//    std::cout << "Checkpoint -1" <<  std::endl; //Debug 

    /* Number of parents and children (listed in topological order) that are
     * either main or auxillary gamma distributed variables. Note that if
     * numparents is 0, this doesn't mean that our graphlab vertex won't have
     * any parents; it simply means that there are no gamma variables as parents
     * (see code below). The same the same holds for numChildren.
     */
    mxArray* numChildren_mx = matGetVariable(pmat, "numchildren");
    unsigned int* numChildren = (unsigned int*) mxGetData(numChildren_mx);
    unsigned int numVertices = mxGetM(numChildren_mx);
    mxArray* numParents_mx = matGetVariable(pmat, "numparents");
    unsigned int* pNumParents = (unsigned int*) mxGetData(numParents_mx);
    numParents.resize(numVertices);
    memcpy(&numParents[0], pNumParents, sizeof(unsigned int) * numVertices);
    
    /* processing backwardMapping to forwardMapping.
     * 
     * What this all is actually about:
     * we started with adjacency matrix over V, columnwise, having K vCols. Then
     * the auxillary gamma variables got into play: their indices range from K
     * to N. Mapping from that ordering to the topological one is stored in
     * backward_mapping.
     * We'll parse backward_mapping to its inverse mapping - forward_mapping, but
     * only for the first K elements (auxillary cols excluded)
     */
/*    mxArray* backwardMapping_mx = matGetVariable(pmat, "backward_mapping");
    assert(mxGetM(backwardMapping_mx) == 1);
    unsigned int* backwardMapping = (unsigned int*) mxGetData(backwardMapping_mx);
*/
    mxArray* forwardMapping_mx = matGetVariable(pmat, "order");
    assert(mxGetM(forwardMapping_mx) == 1);
    //unsigned int* forwardMapping = (unsigned int*) mxGetData(forwardMapping_mx);
    unsigned int* pForwardMapping = ( unsigned int* )mxGetData(forwardMapping_mx);
    forwardMapping.resize(numVertices);
    memcpy(&forwardMapping[0], pForwardMapping, sizeof(unsigned int)*numVertices);       
    mxDestroyArray(forwardMapping_mx);    
    
    /* Gamma hyperparameters
     */
    mxArray* aVe_mx = matGetVariable(pmat, "a_ve");
    double* aVe = (double*) mxGetData(aVe_mx);
    mxArray* bVe_mx = matGetVariable(pmat, "b_ve");
    double* bVe = (double*) mxGetData(bVe_mx);                  
    
    /* Incorporating our topology into graphlab
     * Convention: first neighbor is always a hyperparameter parent
     */
    /* Aux col hyperparameter - all aux gammas will share a single one */
    mxArray* aVeAux_mx = matGetVariable(pmat, "a_ve_auxillary");
    double* aVeAux = (double*) mxGetData(aVeAux_mx);    
    if (numVertices > this->n) {
        ShapeHyperColWrapper vertex;
        ShapeHyperCol* shapeDelegator = (ShapeHyperCol*) vertex.getDelegator();
        /* initializing values */
        memcpy(&shapeDelegator->aValues[0], aVeAux, sizeof (element_type) * k);
        auxHyperCol = g.add_vertex(vertex);
        shapeHyperColIDs.push_back(auxHyperCol);
        std::cout << "AuxHyperCol ID: " << auxHyperCol << "  " << shapeDelegator->aValues[0] <<std::endl; //Debug
    }    
    mxDestroyArray(aVeAux_mx);

    /* First the gamma hyperparameters to have them at the beginning of the
     * neighbor list.
     */
    unsigned int shapeID;
    noparents = 0;
    for (unsigned int curr = 0; curr < numVertices; ++curr) {
     //           std::cout << "forwardMapping[curr] < N " << forwardMapping[curr] << " " << N << std::endl; //DEBUG
        if (forwardMapping[curr] < this->n) {
            /* Processing parents */
            if (numParents[curr] > 0) {
                /* Hyperparameter node
                 */
                ShapeHyperColWrapper shapeVertex;
                ShapeHyperCol* shapeDelegator = (ShapeHyperCol*) shapeVertex.getDelegator();
                memcpy(&shapeDelegator->aValues[0], aVe, sizeof (element_type) * this->k);
                aVe += this->k;     // Move the pointer
                shapeID = g.add_vertex(shapeVertex);
                shapeHyperColIDs.push_back(shapeID);
                std::cout << "shapeHyperCol ID: " << shapeID << std::endl; //Debug
                
                std::cout << "as that went into memcpy" << std::endl;
                for (unsigned int row = 0; row <  this->k; ++row) {
                    std::cout << shapeDelegator->aValues[0] << " ";
                    std::cout << std::endl;
                }                        
                        
                //if (numParents[curr] > 1) { }
                // else it has only auxillary parents
            } else {
                ++noparents;
                /* VCol vertex has no gammas above it -> adding parent
                 * hyperparameter vertices
                 */
                ShapeHyperColWrapper shapeVertex;
                ShapeHyperCol* shapeDelegator = (ShapeHyperCol*) shapeVertex.getDelegator();
                memcpy(&shapeDelegator->aValues[0], aVe, sizeof(element_type) * this->k);
                aVe += this->k;
                unsigned int hyperIndex = g.add_vertex(shapeVertex);
                vHyperColIDs.push_back(hyperIndex);
                std::cout << "vHyperCol ID: " << hyperIndex << std::endl; //Debug
                
                std::cout << "as that went into memcpy" << std::endl;
                for (unsigned int row = 0; row <  this->k; ++row) {
                    std::cout << shapeDelegator->aValues[0] << " ";
                    std::cout << std::endl;
                }

                ScaleHyperColWrapper scaleVertex;
                ScaleHyperCol* scaleDelegator = (ScaleHyperCol*) scaleVertex.getDelegator();
                memcpy(&scaleDelegator->scaleRelatedValues[0], bVe, sizeof(element_type) * this->k);
                bVe += this->k;

                hyperIndex = g.add_vertex(scaleVertex);
                vHyperColIDs.push_back(hyperIndex);
                std::cout << "vHyperCol ID: " << hyperIndex << std::endl; //Debug                                

            }

        }// else {}
    }
    shapeHyperColIDs.resize(shapeHyperColIDs.size());   
    
    /* Then the dirichlet priors, dirichlets and discretes, To have the
     * discretes at the second place (if they exist).
     */
    
    /* Mixture hyperparameters
     */      
    mxArray* expectationsDiscrete_mx = matGetVariable(pmat, "expectations_discrete");
    mxArray* expectationsDirichlet_mx = matGetVariable(pmat, "expectations_dirichlet");
    mxArray* uDirichlet_mx = matGetVariable(pmat, "u_dirichlet");
    unsigned int uCounter = 0;
    
    unsigned int discreteID;
    unsigned int dirichletID;
    unsigned int dirichletHyperID;
    for (unsigned int curr = 0; curr < numVertices; ++curr) {
        /*        std::cout << "forwardMapping[curr] < N " << forwardMapping[curr] << " " << N << std::endl; //DEBUG
         */
                std::cout << "numParents[curr] " << curr << " " <<  numParents[curr] << std::endl; //Debug
        if (forwardMapping[curr] < this->n) {
            /* Processing parents */
            if (numParents[curr] > 0) {
                if (numParents[curr] > 1) {             
                    /* add remaining ancestors, with edges:
                     * DiscreteCol, DirichletCol and DirichletHyperCol
                     */
                    DiscreteColWrapper discreteVertex(numParents[curr]);
                    DiscreteCol* discreteDelegator = (DiscreteCol*) discreteVertex.getDelegator();
                    
                    mxArray* expectationsDiscreteCell_mx = mxGetCell(expectationsDiscrete_mx, curr);
                    double* expectationsDiscreteContents = (double*) mxGetPr(expectationsDiscreteCell_mx);
                    memcpy(&discreteDelegator->values[0], expectationsDiscreteContents,
                            sizeof (element_type) * this->k * numParents[curr]);
                    discreteID = g.add_vertex(discreteVertex);
                    discreteColIDs.push_back(discreteID);                    
//                    expectationsDiscreteContents = NULL;
                    mxDestroyArray(expectationsDiscreteCell_mx);
//                    expectationsDiscreteCell_mx = NULL;
                    mxArray* dirichletHyperCell_mx = mxGetCell(uDirichlet_mx, uCounter++);
                    double* dirichletHyperContents = (double*) mxGetPr(dirichletHyperCell_mx);
                    DirichletHyperColWrapper dirichletHyperVertex( mxGetN(dirichletHyperCell_mx) );
                    DirichletHyperCol* dirichletHyperDelegator = (DirichletHyperCol*) dirichletHyperVertex.getDelegator();
                    memcpy(&dirichletHyperDelegator->values[0], dirichletHyperContents,
                            sizeof (element_type) * this->k * numParents[curr]);
                    dirichletHyperID = g.add_vertex(dirichletHyperVertex);
                    dirichletHyperColIDs.push_back(dirichletHyperID);                    
//                    dirichletHyperContents = NULL;
                    mxDestroyArray(dirichletHyperCell_mx);
//                    dirichletHyperCell_mx = NULL;
                    
                    mxArray* expectationsDirichletCell_mx = mxGetCell(expectationsDirichlet_mx, curr);
                    double* expectationsDirichletContents = (double*) mxGetPr(expectationsDirichletCell_mx);
                    DirichletColWrapper dirichletVertex(mxGetN(expectationsDirichletCell_mx));
                    DirichletCol* dirichletDelegator = (DirichletCol*) dirichletVertex.getDelegator();
                    memcpy(&dirichletDelegator->values[0], expectationsDirichletContents,
                            sizeof(element_type)*this -> k*numParents[curr]);
                    dirichletDelegator->order = numParents[curr];
                    dirichletID = g.add_vertex(dirichletVertex);
                    dirichletColIDs.push_back(dirichletID);
//                    expectationsDirichletContents = NULL;
                    mxDestroyArray(expectationsDirichletCell_mx);
//                    expectationsDirichletCell_mx = NULL;

                }// else it has only auxillary parents

            }// else {}

        }// else {}
    }
    std::cout<<"capacity before: "<< discreteColIDs.capacity()<<std::endl;
    /* 
    discreteColIDs.resize(discreteColIDs.size());
    dirichletHyperColIDs.resize(diri chletHyperColIDs.size());
    dirichletColIDs.resize(dirichletColIDs.size());
    std::cout<<"capacity after: "<< discreteColIDs.capacity()<<std::endl;
*/    
//    mxDestroyArray(expectationsDiscrete_mx);
//    mxDestroyArray(expectationsDirichlet_mx);
//    mxDestroyArray(uDirichlet_mx);    
  
    edge_data edata;
    unsigned int* vertexIndices; /* Storage for Graphlab-related indices of vertices - in the order vertices will have been added to graph (which is topological) */
    vertexIndices = new unsigned int[numVertices];
    unsigned int shapeIterator = 1; // skip the aux hyper
    unsigned int vIterator = 0;
    unsigned int discreteIterator = 0;
    unsigned int dirichletIterator = 0;
    unsigned int dirichletHyperIterator = 0;
    
    map<unsigned int, unsigned int> helper;        //  <ID_of_AuxCol, next_mixture> helper;

    for (unsigned int curr = 0; curr < numVertices; ++curr) {
        if (forwardMapping[curr] < this->n) {
            /* Processing parents */
            unsigned int ourID;
            if (numParents[curr] > 0) {
                
                bool hasDiscrete ;
                bool hasGammaChild;
                if (numParents[curr] > 1) {
                    hasDiscrete = true;
                } else {
                    hasDiscrete = false;
                }

                /* Current describes what should be embodied into VCol */
                VColWrapper vertex(hasDiscrete, false); // no gamma child
                VCol* vDelegator = (VCol*) vertex.getDelegator();
                ourID = g.add_vertex(vertex);
                vColIDs.push_back(ourID);
                std::cout << "VCol ID: " << ourID << std::endl; //Debug
                vertexIndices[curr] = ourID;

                /* edges to existing parents */
                shapeID = shapeHyperColIDs[shapeIterator++];
                g.add_edge(shapeID, ourID, edata);
                std::cout << shapeID << " -> " << ourID << std::endl;
                g.add_edge(ourID, shapeID, edata);
                std::cout << ourID << " -> " << shapeID << std::endl;

                bool isMixture = false;
                
                if (numParents[curr] > 1) {
//                    std::cout << "Checkpoint 102"; //Debug
                    
                    isMixture = true;
                    
                    discreteID = discreteColIDs[discreteIterator++];
                    dirichletID = dirichletColIDs[dirichletIterator++];
                    dirichletHyperID = dirichletHyperColIDs[dirichletHyperIterator++];

                    g.add_edge(discreteID, ourID, edata);
//                    std::cout << discreteID << " -> " << ourID << std::endl;
                    g.add_edge(ourID, discreteID, edata);
                    std::cout << ourID << " -> " << discreteID << std::endl;
                    
                    /* Shape hyperparameter of currentVCol to Discrete parent */
                    g.add_edge(shapeID, discreteID, edata);
//                    std::cout << "Shape coparent of Discrete: "<<shapeID << " -> " << discreteID << std::endl;

                    g.add_edge(dirichletHyperID, dirichletID, edata);
                    std::cout << dirichletHyperID << " -> " << dirichletID << std::endl;
                    g.add_edge(dirichletID, dirichletHyperID, edata);
                    std::cout << dirichletID << " -> " << dirichletHyperID << std::endl;
                           
                    g.add_edge(dirichletID, discreteID, edata);
                    std::cout << dirichletID << " -> " << discreteID << std::endl;
                    g.add_edge(discreteID, dirichletID, edata);
                    std::cout << discreteID << " -> " << dirichletID << std::endl;
                }

                /* vertex has auxillary parents. They have already been
                 * added to graph (we add vertices in topological order) so only
                 * edges to and from parents are to be added.
                 * 
                 * Important: at this point cols_mask is required to be sorted
                 * in increasing order.
                 */
                /* TODO: fix location, make it unsigned */
                long int location = binary_search(colSupport, 0, adjacencySparsity - 1, curr);
                // colSupport[location] == current holds;                
                /* The following loop sets location to first entry in colSupport
                 * such that colSupport[location] == current
                 * (linear time)
                 */
                while (colSupport[location] == curr) {
                    --location;
                    if (location < 0) {
                        break;
                    }
                }
                ++location;

//                std::cout << "curr = " << curr << " numparents = " << numParents[curr] << " location = " << location << " parents: " << std::endl; //DEBUG
                /* Now we're ready to fetch all curr's parents:  */
                for (unsigned int i = 0; i < numParents[curr]; ++i) {
                    // ID of current's parents in graphlab: vertexIndices[ rowSupport[location++] ]                                        

                    unsigned int parentID = vertexIndices[ rowSupport[location++] ];
                    AuxCol* parent = (AuxCol*) g.vertex_data(parentID).getDelegator();
//                    ++parent->numChildren;
                    parent->childIsMixture[ helper[parentID] ] = isMixture;
                    parent->indexInMixture[ helper[parentID] ] = i;
                    ++helper[parentID];

                    g.add_edge(parentID, ourID, edata);
                    std::cout << parentID << " -> " << ourID << std::endl;
                    g.add_edge(ourID, parentID, edata);
                    std::cout << ourID << " -> " << parentID << std::endl;
                    // parent's coparents; current's hyperparameters
                    g.add_edge(shapeID, parentID, edata);
//                    g.add_edge(parentID, shapeID, edata);
                    std::cout << "coparent, shape: " << shapeID << " -> " << parentID << std::endl;
                    // ovdje si stao. provjeri id koji se pojavljuje 
                    // parent's coparents; current's discrete distro
                    if (numParents[curr] > 1) {
                        g.add_edge(discreteID, parentID, edata);
//                        std::cout << "coparent, discrete: " << discreteID << " -> " << parentID << std::endl;
                        // and vice versa
                        g.add_edge(parentID, discreteID, edata);
//                        std::cout << "coparent, discrete: " << parentID << " -> " << discreteID << std::endl;
                    }
                }
            } else {
//                std::cout << "Checkpoint 103"; //Debug

//                                std::cout << "Adding VColWrapper, no parents" << std::endl; //DEBUG               

                /* Current describes what should be embodied into VCol */
                VColWrapper vertex(false, false);       // no discrete, no gamma child
                VCol* delegator = (VCol*) vertex.getDelegator();
                ourID = g.add_vertex(vertex);
//                std::cout << "befpre poushback"; //Debug
                vColIDs.push_back(ourID);
//                std::cout << "VCol ID: " << ourID << std::endl; //Debug
                vertexIndices[curr] = ourID;

                /* Adding edges */
                //shape
                unsigned int hyperIndex = vHyperColIDs[vIterator++];
                g.add_edge(hyperIndex, ourID, edata);
//                std::cout << hyperIndex << " -> " << ourID << std::endl;
                g.add_edge(ourID, hyperIndex, edata);
//                std::cout << ourID << " -> " << hyperIndex << std::endl;
                //scale
                hyperIndex = vHyperColIDs[vIterator++];
                g.add_edge(hyperIndex, ourID, edata);
//                std::cout << hyperIndex << " -> " << ourID << std::endl;
                g.add_edge(ourID, hyperIndex, edata);
//                std::cout << ourID << " -> " << hyperIndex << std::endl;

                //                std::cout << "adding edges - OK" << std::endl; //DEBUG
            }

        } else {
//            std::cout << "Adding AuxColWrapper" << std::endl; //DEBUG
            //            std::cout << "curr = " << curr << std::endl; //DEBUG
            /* current describes what should be an AuxCol */

            /* Adding vertex */
            AuxColWrapper vertex(numChildren[curr]);
            AuxCol* delegator = (AuxCol*) vertex.getDelegator();
            unsigned int ourID = g.add_vertex(vertex);
            auxColIDs.push_back(ourID);
//            std::cout << "AuxCol ID: " << ourID << std::endl; //Debug
            vertexIndices[curr] = ourID;
            helper[ourID] = 0;

            // we have only one parent - let's find it
            /* TODO: fix location, make it unsigned */
            
            long int location = binary_search(colSupport, 0, adjacencySparsity - 1, curr);
            std::cout << "location " << location << " sparsity " << adjacencySparsity << " curr "<< curr <<std::endl; //DEBUG

            unsigned int parentID = vertexIndices[ rowSupport[location] ];
            std::cout << "parentID " << parentID << std::endl; //DEBUG

            VColWrapper& vWrapper = (VColWrapper&) g.vertex_data(parentID);
            VCol* parentDelegator = (VCol*) vWrapper.getDelegator();
            parentDelegator->hasGammaChild = true;
            
//            std::cout << "aaaaa" <<std::endl; //DEBUG
     
            /* Adding edges */
            g.add_edge(ourID, parentID, edata);
//            std::cout << ourID << " -> " << parentID << std::endl;
            g.add_edge(parentID, ourID, edata);
//            std::cout << parentID << " -> " << ourID << std::endl;
            g.add_edge(auxHyperCol, ourID, edata);
//            std::cout << auxHyperCol << " -> " << ourID << std::endl;
            g.add_edge(ourID, auxHyperCol, edata);
//            std::cout << ourID << " -> " << auxHyperCol << std::endl;
            // and an edge from auxHyperCol to parentID (auxHyperCol is
            // coparent to parentID)
            g.add_edge(auxHyperCol, parentID, edata);
//            std::cout << "coparent: " << parentID << " -> " << auxHyperCol << std::endl;
        }
        
    }
    delete(vertexIndices);
//    vColIDs.resize(vColIDs.size());
//    auxColIDs.resize(vColIDs.size());

//        std::cout << "adjacency OK" << std::endl; //Debug

    /* Don't need those anymore */
    mxDestroyArray(rowSupport_mx);
    mxDestroyArray(colSupport_mx);
    mxDestroyArray(numChildren_mx);

    /* As we haven't initialized gamma statistics (memory requirements), we'll
     * do a second pass
     * Edit: mxGetVariable returns a pointer so initializing all the variables
     * at once shouldn't stress the memory out. TODO: put all initialization
     * code in a single loop
     */

    /* gamma statistics: taken in the order in which parentless VCols
     * appear when going from 0 to numVertices
     */
    mxArray* markov1_mx = matGetVariable(pmat, "markov_statistic_1");
    double* markov1 = (double*) mxGetData(markov1_mx);
    // TODO: check dimensions
    mxArray* markov2_mx = matGetVariable(pmat, "markov_statistic_2");
    double* markov2 = (double*) mxGetData(markov2_mx);

    size_t vColCnt = 0;
    size_t auxColCnt = 0;
    for (unsigned int curr = 0; curr < numVertices; ++curr) {

        if (forwardMapping[curr] < this->n) {

//            std::cout << "initializing a VCol, " << vColIDs[vColCnt] << endl;
            const VColWrapper& vertex = (VColWrapper&) g.vertex_data(vColIDs[vColCnt]);
            ++vColCnt;
            VCol* delegator = (VCol*) vertex.getDelegator();
            /* initializing values */
            memcpy(&delegator->scaleRelatedValues[0], markov1, sizeof (element_type) * this->k);
            markov1 += this->k;           
            memcpy(&delegator->expELogValues[0], markov2, sizeof (element_type) * this->k);
            markov2 += this->k;

        } else {
            // initializing AuxCols

//            std::cout << "initializing an AuxCol, " << auxColIDs[auxColCnt] << endl;

            // we have an AuxCol            
            AuxColWrapper& vertex = (AuxColWrapper&) g.vertex_data(auxColIDs[auxColCnt++]);
            AuxCol* delegator = (AuxCol*) vertex.getDelegator();
            /* initializing values */
            memcpy(&delegator->scaleRelatedValues[0], markov1, sizeof (element_type) * this->k);
            markov1 += this->k;
            memcpy(&delegator->eLogValues[0], markov2, sizeof (element_type) * this->k);
            markov2 += this->k;
        }
    }
    
            /* Don't need those anymore */
    mxDestroyArray(markov1_mx);
    mxDestroyArray(markov2_mx);
    mxDestroyArray(aVe_mx);
    mxDestroyArray(bVe_mx);
    mxDestroyArray(numParents_mx);

    mxArray* aTm_mx = matGetVariable(pmat, "a_tm_transposed");
    double* aTm = (double*) mxGetData(aTm_mx);
    mxArray* bTm_mx = matGetVariable(pmat, "b_tm_transposed");
    double* bTm = (double*) mxGetData(bTm_mx);

    tHyperRowOffset = vColIDs.size() + auxColIDs.size() +
            vHyperColIDs.size() + shapeHyperColIDs.size() +
            discreteColIDs.size() + dirichletColIDs.size() + dirichletHyperColIDs.size();
    for (unsigned int row = 0; row < this->m; ++row) {
        THyperRowWrapper vertex;
        THyperRow* delegator = (THyperRow*) vertex.getDelegator();
        memcpy(&delegator->aValues[0], aTm, sizeof (element_type) * this->k);
        memcpy(&delegator->bValues[0], bTm, sizeof (element_type) * this->k);
        aTm += this->k;
        bTm += this->k;
        unsigned int bla =
                g.add_vertex(vertex);
//        std::cout << "THyperRow ID: " << bla << std::endl; //Debug

        // debug
/*        for (unsigned int i = 0; i < this->k; ++i) {
            std::cout << delegator->aValues[i] << " ";
        }
        std::cout << endl;

        for (unsigned int i = 0; i < this->k; ++i) {
            std::cout << delegator->bValues[i] << " ";
        }
        std::cout << endl;
        // debug ends
*/    }
    std::cout << "Destroying arrays  aTm_mx, bTm_mx";
    mxDestroyArray(aTm_mx);
    mxDestroyArray(bTm_mx);

    // M komada

    std::cout << "Checkpoint 1";
          
    /* T, rowwise
     */
    mxArray* eT_mx = matGetVariable(pmat, "E_t_transposed");
    element_type* eT = (double*) mxGetData(eT_mx);
    mxArray* lT_mx = matGetVariable(pmat, "L_t_transposed");
    element_type* lT = (double*) mxGetData(lT_mx);

    tRowOffset = tHyperRowOffset + this->m;
    for (unsigned int row = 0; row < this->m; ++row) {
//        std::cout << "row: " << row; //Debug
        TRowWrapper vertex(row);
        TRow* delegator = (TRow*) vertex.getDelegator();
        delegator->row = row;
        
//        std::cout << " copying "<<*eT<<" "<<*lT<<std::endl;; //Debug

        memcpy(&delegator->eValues[0], eT, sizeof(element_type) * this->k);
        memcpy(&delegator->expELogValues[0], lT, sizeof(element_type) * this->k);
        eT += this->k;
        lT += this->k;

        unsigned int bla =
                g.add_vertex(vertex); // == tRowOffset+row
        
//        std::cout << " added "<<std::endl; //Debug
        
        const TRowWrapper& vertex2 = (TRowWrapper&) g.vertex_data(bla);
        const TRow& delegator2 = (const TRow&) *vertex2.getDelegator();
//        std::cout << "row: "<< row<< "; TRow ID: " << bla << "; " << delegator2.eValues[0]<<std::endl; //Debug
    }
    
    std::cout << "Destroying arrays  "<<std::endl;
    mxDestroyArray(eT_mx);
    mxDestroyArray(lT_mx);  
    // M komada
        
    std::cout << "Further loading matlab stuff  "<<std::endl;
        
    mxArray* colwiseNum_mx = matGetVariable(pmat, "colwise_num"); // number of nonzeros per column in the mask    
    unsigned int* colwiseNum = (unsigned int*) mxGetData(colwiseNum_mx);
    assert(mxGetN(colwiseNum_mx) == n);        
    
//    std::cout << "here:  "<<std::endl;

    /* tmp sparse, columnwise
     */
    tmpColOffset = tRowOffset + this->m;
    for (unsigned int col = 0; col < this->n; ++col) {        
        
        if (colwiseNum[col] > 0)        {
            TmpColWrapper vertex(colwiseNum[col]);
            TmpCol* delegator = (TmpCol*) vertex.getDelegator();           
            // values not initialized because of the order vertices get updated
            unsigned int bla =  g.add_vertex(vertex);
//            std::cout << "Col: " << col << "; TmpCol ID: " << bla << std::endl; //Debug
            //            //        }         // if(location != -1) ends
        } else {
            TmpColWrapper vertex(colwiseNum[col]);
            TmpCol* delegator = (TmpCol*) vertex.getDelegator();
//            delegator->values = new element_type[1];
        }
    }
    // N komada

    /* SigT, rowwise
     */
    std::cout << "SigTRows " << std::endl; //Debug
    sigTRowOffset = tmpColOffset + this->n;
    for (unsigned int row = 0; row < this->m; ++row) {
        SigTRowWrapper vertex;
        SigTRow* delegator = (SigTRow*) vertex.getDelegator();
        unsigned int bla = g.add_vertex(vertex);
//        std::cout << "Row: "<< row << "; SigTRow ID: " << bla << std::endl; //Debug
    }
    // M komada

    /* SigV, colwise
     */
    std::cout << "SigVCols " << std::endl; //Debug
    sigVColOffset = sigTRowOffset + this->m;
    for (unsigned int col = 0; col < this->n; ++col) {
        SigVColWrapper vertex;
        SigVCol* delegator = (SigVCol*) vertex.getDelegator();
        unsigned int bla =
                g.add_vertex(vertex);
//        std::cout << "Col: "<< col << "; SigVCol ID: " << bla << std::endl; //Debug
    }
    // N komada


    /* processing X - original data matrix, columnwise
     */
    std::cout << "x_raw: " << std::endl; //Debug
    mxArray* colwiseX_mx = matGetVariable(pmat, "x_raw"); // mask of x
    double* colwiseX = (double*) mxGetData(colwiseX_mx);
    mxArray* colwiseMaskCols_mx = matGetVariable(pmat, "colwise_mask_cols"); // mask of x
    unsigned int* colwiseMaskCols = (unsigned int*) mxGetData(colwiseMaskCols_mx);
    size_t maskSparsity = mxGetM(colwiseMaskCols_mx); // number of nonzero entries
    mxArray* colwiseMaskRows_mx = matGetVariable(pmat, "colwise_mask_rows"); // mask of x
    unsigned int* colwiseMaskRows = (unsigned int*) mxGetData(colwiseMaskRows_mx);

//    std::cout << "here 2 " << *colwiseMaskCols << " " << maskSparsity << std::endl; //Debug
     
    unsigned int* colwiseMaskRowsBreakPoints[this->n];
    xColOffset = sigVColOffset + this->n; // starting ID for both XCols and SupportXCols - every other is a SupportXCols;
//    std::cout << "here 2 " << std::endl; //Debug
    for (unsigned int col = 0; col < this->n; ++col) {
        if (colwiseNum[col] > 0) {
//            std::cout << "before bin search" << std::endl; //Debug
            long int location = binary_search(colwiseMaskCols, 0, maskSparsity - 1, col);
//            std::cout << "location = " << location << std::endl; //Debug
//std::cout << "after bin search" << std::endl; //Debug
            /* The following loop sets location to first entry in colSupport
             * such that colSupport[location] == current
             * (additional intervention linear in nonzeros per column)
             */
            //        if(location != -1)      {       // column is not empty
            bool broken = false;
            while (colwiseMaskCols[location] == col) {
                --location;
                if (location == 0) {
                    broken = true;
                    break;
                }
            }
            if (!broken) {
                ++location;
            }
            colwiseMaskRowsBreakPoints[col] = &colwiseMaskRows[location];
                        
            XColWrapper vertex(colwiseNum[col]);
            XCol* delegator = (XCol*) vertex.getDelegator();            
            // copy the local values and number of nonzeros to our node
            memcpy(&delegator->values[0], colwiseX, sizeof(element_type) * delegator->getNumElements());
            colwiseX += delegator->getNumElements(); // prepare pointer for the next iteration
            unsigned int bla = g.add_vertex(vertex);
            std::cout << "XCol ID: " << bla << std::endl; //Debug
            
            SupportXColWrapper vertexW(colwiseNum[col]);
            SupportXCol* delegatorW = (SupportXCol*) vertexW.getDelegator();            
            // copy the supports for this column
            memcpy(&delegatorW->support[0], &colwiseMaskRows[location], sizeof(unsigned int)*delegatorW->getNumElements());
            
            bla = g.add_vertex(vertexW);
            std::cout << "SupportXCol ID: " << bla << std::endl; //Debug

            //        }         // if(location != -1) ends
        } else {
            
            XColWrapper vertex(0);
            XCol* delegator = (XCol*) vertex.getDelegator();
//            delegator->values = new element_type[1]; // otherwise destructor wouldn't work
//            delegator->values[0] = -1;

            SupportXColWrapper vertexW(0);
            SupportXCol* delegatorW = (SupportXCol*) vertexW.getDelegator();    
//            delegatorW->support = new unsigned int[1];

            unsigned int bla =  g.add_vertex(vertex);
            std::cout << "XCol ID: " << bla << std::endl; //Debug
            bla = g.add_vertex(vertexW);
            std::cout << "SupportXCol ID: " << bla << std::endl; //Debug
        }
    }
    // N komada

    mxDestroyArray(colwiseX_mx);

    
    
    std::cout << "going for edges          " << std::endl; //Debug
    /********************************** Edges
     */

    /*    std::cout << "colwiseMaskRowsBreakPoints          " << std::endl; //Debug
        for (unsigned int col = 0; col < N; ++col) {
            std::cout << *colwiseMaskRowsBreakPoints[col] << std::endl;
        }
        std::cout << std::endl;
     */
    for (unsigned int col = 0; col < this->n; ++col) {

        //        std::cout << "col " << col << " " << colwiseNum[col]<<std::endl; //Debug              
        edge_data edata;

//dbg        std::cout << "    Us: Vcol" << std::endl;
        unsigned int us = vColIDs[col];

        if (colwiseNum[col] > 0) {
            /* T to us - we need rows of T where x has columwise support */
//dbg            std::cout << "T to us" << std::endl;
            for (unsigned int supportIter = 0; supportIter < colwiseNum[col]; ++supportIter) {

                /*                std::cout << "  bmapping: " << backwardMapping[col] << " ";
                                std::cout << *(colwiseMaskRowsBreakPoints[col] + supportIter) << " ";
                 */ g.add_edge(tRowOffset + *(colwiseMaskRowsBreakPoints[col] + supportIter), us, edata);
//dbg                std::cout << tRowOffset + *(colwiseMaskRowsBreakPoints[col] + supportIter) << " -> " << us << std::endl;

            }
            //            std::cout << std::endl;
        }

        /* sigVCol to us */
//dbg        std::cout << "sigVCol to us" << std::endl;
        g.add_edge(sigVColOffset + col, us, edata);
//dbg        std::cout << sigVColOffset + col << " -> " << us << std::endl;

        /* 
         * Us: tmpCol
         */
//dbg        std::cout << "    Us: tmpCol" << std::endl;
        us = tmpColOffset + col;
        /* XCol to us */
//dbg        std::cout << "XCol to us" << std::endl;
        g.add_edge(xColOffset + 2 * col, us, edata);
//dbg        std::cout << xColOffset + 2 * col << " -> " << us << std::endl;

        if (colwiseNum[col] > 0) {
            /* T to us */
            for (unsigned int supportIter = 0; supportIter < colwiseNum[col]; ++supportIter) {
//dbg                std::cout << "T to us" << std::endl;
                g.add_edge(tRowOffset + *(colwiseMaskRowsBreakPoints[col] + supportIter), us, edata);
//dbg                std::cout << tRowOffset + *(colwiseMaskRowsBreakPoints[col] + supportIter) << " -> " << us << std::endl;
            }
            /* V to us */
//dbg            std::cout << "V to us" << std::endl;
            g.add_edge(vColIDs[col], us, edata);
//dbg            std::cout << vColIDs[col] << " -> " << us << std::endl;
        }

        /* 
         * Us: sigV 
         * We'll need a VCol, a TmpCol and some TRows
         */
//dbg        std::cout << "    Us: sigV" << std::endl;
        us = sigVColOffset + col;
        /* V to us */
//dbg        std::cout << "V to us" << std::endl;
        g.add_edge(vColIDs[col], us, edata);
//dbg        std::cout << vColIDs[col] << " -> " << us << std::endl;
        /* T to us    - we need rows of T where x has columwise support */
//dbg        std::cout << "T to us" << std::endl;
        if (colwiseNum[col] > 0) {
            for (unsigned int supportIter = 0; supportIter < colwiseNum[col]; ++supportIter) {
                g.add_edge(tRowOffset + *(colwiseMaskRowsBreakPoints[col] + supportIter), us, edata);
//dbg                std::cout << tRowOffset + *(colwiseMaskRowsBreakPoints[col] + supportIter) << " -> " << us << std::endl;
            }
        }
        /* col of tmp to us */
//dbg        std::cout << "tmp to us" << std::endl;
        g.add_edge(tmpColOffset + col, us, edata);
//dbg        std::cout << tmpColOffset + col << " -> " << us << std::endl;
    }
    
        /* We'll also need the rowwise support. Unsorted is OK */
    vector< vector<unsigned int> > rowwiseSupport(this->m*this->n);
    unsigned int* pCol = colwiseMaskCols;
    unsigned int* pRow = colwiseMaskRows;
//    std::cout << "rowwise support:  "<<std::endl;    
    for (unsigned int col = 0; col < this->n; ++col) {
        for (unsigned int rowIter = 0; rowIter < colwiseNum[col]; ++rowIter) {
            rowwiseSupport[ *(pRow++) ].push_back(*(pCol++));
        }
    }
    mxDestroyArray(colwiseNum_mx);  
    for (unsigned int row = 0; row < this->m; ++row) {
        rowwiseSupport[row].resize(rowwiseSupport[row].size(), 0);
    }
    //DBG
/*    for (unsigned int row = 0; row < this->m; ++row) {
        std::sort(rowwiseSupport[row].begin(), rowwiseSupport[row].end()); // default operator is <
        for (int i = 0; i < rowwiseSupport[row].size(); ++i) {
            std::cout << rowwiseSupport[row][i] << " ";
        }
        std::cout << std::endl;
    }
*/

    for (unsigned int row = 0; row < this->m; ++row) {
        edge_data edata;
        /*
         * Us: T
         */
//        std::cout << "    us: T" << std::endl;
        unsigned int us = tRowOffset + row;
        /* V -> us:         
         * Determine which Vs we need using rowwiseSupport and
         * correspondingly add edges: */
//        std::cout << "V to us" << std::endl;
        for (unsigned int supportIter = 0; supportIter < rowwiseSupport[row].size(); ++supportIter) {
            // Note that we don't need support explicitly.
            //g.add_edge(vColOffset + backwardMapping[ rowwiseSupport[row][supportIter] ], us, edata);
            g.add_edge(vColIDs[ rowwiseSupport[row][supportIter] ], us, edata);
//            std::cout << vColIDs[ rowwiseSupport[row][supportIter] ] << " -> " << us << std::endl;
        }
        /* at sigTRow row us */
//        std::cout << "sigTRow to us" << std::endl;
        g.add_edge(sigTRowOffset + row, us, edata);
//        std::cout << sigTRowOffset + row << " -> " << us << std::endl;
        /* hyperparameter row to TRow */
//        std::cout << "hyperparameter row to TRow" << std::endl;
        g.add_edge(tHyperRowOffset + row, us, edata);
//        std::cout << tHyperRowOffset + row << " -> " << us << std::endl;

        /* 
         * Us: sigT
         * We'll need a LtRow, TmpCols and LvCols
         * Our ID is sigTRowOffset+row
         */
//        std::cout << "    us: sigT" << std::endl;
        us = sigTRowOffset + row;
        /* T -> us */
//        std::cout << "T -> us" << std::endl;
        g.add_edge(tRowOffset + row, us, edata);
//        std::cout << tRowOffset + row << " -> " << us << std::endl;
        /* V and tmp and support of x -> us
         * determine which Lvs we need using rowwiseSupport and correspondingly
         * add edges */
//        std::cout << "V and tmp and support of x -> us" << std::endl;
        for (unsigned int supportIter = 0; supportIter < rowwiseSupport[row].size(); ++supportIter) {
            //g.add_edge(vColOffset + forwardMapping[ rowwiseSupport[row][supportIter] ], us, edata);
//            std::cout << "asdfasdf" << std::endl;
            g.add_edge(vColIDs[ rowwiseSupport[row][supportIter] ], us, edata);
//            std::cout << "vcol" << vColIDs[ rowwiseSupport[row][supportIter] ] << " -> " << us << std::endl;
            g.add_edge(tmpColOffset + rowwiseSupport[row][supportIter], us, edata);
            std::cout << "tmpcol" << tmpColOffset + rowwiseSupport[row][supportIter] << " -> " << us << std::endl;
            g.add_edge(xColOffset + 1 + 2 *rowwiseSupport[row][supportIter], us, edata);
            std::cout << "xcol" << xColOffset + 1 + 2 * rowwiseSupport[row][supportIter] << " -> " << us << std::endl;
        }
    }
    g.finalize();

    mxDestroyArray(colwiseMaskCols_mx);
    mxDestroyArray(colwiseMaskRows_mx);    
    
    matClose(pmat);
  
    // backward_mapping still lives
}

/* Outputs the parameters and moments and everything to a .mat file.
 * Problems with serialization (pointer tracking). Until it is functional, 
 * the state and outputs of the algorithm will be stored this way.
 */
void ModelNMF::matOutput(graph_type& g, const char* const filename)     {
        
    //MATFile* pmat = matOpen(filename, "w");
    std::cout << filename << std::endl;
         
    MATFile* pmat = matOpen(filename, "w7.3");
    assert(pmat != NULL);
    int status;
    /* Reminder: mxarrays are column major
     */   
    mxArray* mx;
    mxArray* mx2;
    mxArray* mx3;
    mxArray* mx4;
    double* pDouble;
    double* pDouble2;
    double* pDouble3;
    double* pDouble4;
    
     std::cout << "Checkpoint" << std::endl;
    
    /*Discrete*/
    mx = mxCreateCellMatrix(this->numParents.size(), 1);
    unsigned int v = 0;
    for (size_t i = 0; i < this->numParents.size(); ++i)    {
        if(numParents[i]>1)     {
            const VertexProxy& wrapper = g.vertex_data(discreteColIDs[v]);
            const DiscreteCol* const d = (const DiscreteCol* const) wrapper.getDelegator();
            mxArray* contents_mx = mxCreateDoubleMatrix(VertexProxy::K->get(), d->getMixtureSize(), mxREAL);
            pDouble = (double*) mxGetData(contents_mx);
            memcpy(pDouble, &d->values[0], sizeof(element_type)*VertexProxy::K->get()*d->getMixtureSize());
            mxSetCell(mx, i, contents_mx);
            ++v;
        }
    }        
    status = matPutVariable(pmat, "expectations_discrete", mx);
//    std::cout << "expectations_discrete" << std::endl;
    
    /*Dirichlet*/
    mx = mxCreateCellMatrix(this->numParents.size(), 1);
    v=0;
    for (size_t i = 0; i < this->numParents.size(); ++i)    {
        if(numParents[i]>1)     {
            const VertexProxy& wrapper = g.vertex_data(dirichletColIDs[v]);
            const DirichletCol* const d = (const DirichletCol* const) wrapper.getDelegator();
            mxArray* contents_mx = mxCreateDoubleMatrix(VertexProxy::K->get(), d->getOrder(), mxREAL);
            pDouble = (double*) mxGetData(contents_mx);
            memcpy(pDouble, &d->values[0], sizeof(element_type)*VertexProxy::K->get()*d->getOrder());
            mxSetCell(mx, i, contents_mx);
            ++v;
        }
    }
    status = matPutVariable(pmat, "expectations_dirichlet", mx);
    
    /*Dirichlet hyper*/
    mx = mxCreateCellMatrix(this->numParents.size(), 1);
    for (v = 0; v < this->dirichletHyperColIDs.size(); ++v)     {
        const VertexProxy& wrapper = g.vertex_data(dirichletHyperColIDs[v]);
        const DirichletHyperCol* const d = (const DirichletHyperCol* const) wrapper.getDelegator();
        mxArray* contents_mx = mxCreateDoubleMatrix(VertexProxy::K->get(), d->getOrder(), mxREAL);
        pDouble = (double*) mxGetData(contents_mx);
        memcpy(pDouble, &d->values[0], sizeof(element_type)*VertexProxy::K->get()*d->getOrder());
        mxSetCell(mx, v, contents_mx);
    }
    status = matPutVariable(pmat, "u_dirichlet", mx);
    
//    std::cout << "Checkpoint" << std::endl;
    
    /*Markov expectations - right matrix & auxillary gammas*/
    /*shapes, scales*/           
    unsigned int shapeIndex = 0;
    if( forwardMapping.size() > N.get() )   {
        mx = mxCreateDoubleMatrix(VertexProxy::K->get(), 1, mxREAL);
        pDouble3 = (double*) mxGetData(mx);
//        std::cout << "before first shape" << std::endl;
        const VertexProxy& wrapper = g.vertex_data(shapeHyperColIDs[shapeIndex]);
        const ShapeHyperCol* const d = (const ShapeHyperCol* const) wrapper.getDelegator();
        memcpy(pDouble3, &d->aValues[0], sizeof(element_type)*VertexProxy::K->get());
        ++shapeIndex;
    }
    status = matPutVariable(pmat, "a_ve_auxillary", mx);    
    
    unsigned int vIndex = 0;
    unsigned int auxIndex = 0;
    unsigned int vHyperIndex = 0;
    mx = mxCreateDoubleMatrix(VertexProxy::K->get(), forwardMapping.size(), mxREAL);
    pDouble = (double*) mxGetData(mx);
    mx2 = mxCreateDoubleMatrix(VertexProxy::K->get(), forwardMapping.size(), mxREAL);    
    pDouble2 = (double*) mxGetData(mx2);
    mx3 = mxCreateDoubleMatrix(VertexProxy::K->get(), vHyperColIDs.size()/2+ shapeHyperColIDs.size()-1, mxREAL);
    pDouble3 = (double*) mxGetData(mx3);
    mx4 = mxCreateDoubleMatrix(VertexProxy::K->get(), noparents, mxREAL);
    pDouble4 = (double*) mxGetData(mx4);
    std::cout << "Before looping: loops = " << n << std::endl;
    for (unsigned int curr = 0; curr < forwardMapping.size(); ++curr) {        
//        std::cout << "          curr: " << curr <<std::endl;
        if (forwardMapping[curr] < n)   {
            // VCol related
            if (numParents[curr] > 0)   {
                // shape
//                std::cout << "here 1 " <<std::endl;
                const VertexProxy& wrapper = g.vertex_data(shapeHyperColIDs[shapeIndex]);
                const ShapeHyperCol* const d = (const ShapeHyperCol* const) wrapper.getDelegator();
                memcpy(pDouble3, &d->aValues[0], sizeof (element_type) * VertexProxy::K->get());
                pDouble3 += VertexProxy::K->get();
                ++shapeIndex;
            } else {
                // shape
//                std::cout << "here 1 " <<std::endl;
                const VertexProxy& wrapper = g.vertex_data(vHyperColIDs[vHyperIndex]);
                const ShapeHyperCol* const d = (const ShapeHyperCol* const) wrapper.getDelegator();
                memcpy(pDouble3, &d->aValues[0], sizeof (element_type) * VertexProxy::K->get());
                pDouble3 += VertexProxy::K->get();
                ++vHyperIndex;
                //scale hyperparameter
//                std::cout << "here 2 " <<std::endl;
                const VertexProxy& wrapper2 = g.vertex_data(vHyperColIDs[vHyperIndex]);
                const ScaleHyperCol* const d2 = (const ScaleHyperCol* const) wrapper2.getDelegator();                
                memcpy(pDouble4, &d2->scaleRelatedValues[0], sizeof(element_type)*VertexProxy::K->get());
                pDouble4 += VertexProxy::K->get();
                ++vHyperIndex;
            }
//            std::cout << "here 3 " <<std::endl;
            const VertexProxy& wrapper = g.vertex_data(vColIDs[vIndex]);
            const VCol* const d = (const VCol* const) wrapper.getDelegator();
            memcpy(pDouble, &d->scaleRelatedValues[0], sizeof(element_type)*VertexProxy::K->get());
            memcpy(pDouble2, &d->expELogValues[0], sizeof(element_type)*VertexProxy::K->get());
            pDouble     += VertexProxy::K->get();
            pDouble2    += VertexProxy::K->get();
            ++vIndex;
        } else {
            // current describes what should be an AuxCol related
//            std::cout << "here 4 " <<auxIndex<<" " << auxColIDs[auxIndex]<< std::endl;
            const VertexProxy& wrapper = g.vertex_data(auxColIDs[auxIndex]);
            const AuxCol* const d = (const AuxCol* const) wrapper.getDelegator();
            memcpy(pDouble, &d->scaleRelatedValues[0], sizeof(element_type)*VertexProxy::K->get());
            memcpy(pDouble2, &d->eLogValues[0], sizeof(element_type)*VertexProxy::K->get());
            pDouble     += VertexProxy::K->get();
            pDouble2    += VertexProxy::K->get();
            ++auxIndex;                                        
        }
    }    
    status = matPutVariable(pmat, "markov_statistic_1", mx);
    status = matPutVariable(pmat, "markov_statistic_2", mx2);
    status = matPutVariable(pmat, "a_ve", mx3);
    status = matPutVariable(pmat, "b_ve", mx4);
    
    /*Left matrix*/
    mx = mxCreateDoubleMatrix(VertexProxy::K->get(), VertexProxy::M->get(), mxREAL);
    pDouble = (double*) mxGetData(mx);
    mx2 = mxCreateDoubleMatrix(VertexProxy::K->get(), VertexProxy::M->get(), mxREAL);
    pDouble2 = (double*) mxGetData(mx2);
    for (size_t v = tRowOffset; v < tRowOffset + VertexProxy::M->get(); ++v)    {
        const VertexProxy& wrapper = g.vertex_data(v);
        const TRow* const t = (const TRow* const) wrapper.getDelegator();
        memcpy(pDouble, &t->eValues[0], sizeof(element_type)*VertexProxy::K->get());
        memcpy(pDouble2, &t->expELogValues[0], sizeof(element_type)*VertexProxy::K->get());
        pDouble  += VertexProxy::K->get();
        pDouble2 += VertexProxy::K->get();
    }
    status = matPutVariable(pmat, "E_t_transposed", mx);
    status = matPutVariable(pmat, "L_t_transposed", mx2);
    
    /* Left matrix - hyperparameters */
    mx = mxCreateDoubleMatrix(VertexProxy::K->get(), VertexProxy::M->get(), mxREAL);
    pDouble = (double*) mxGetData(mx);
    mx2 = mxCreateDoubleMatrix(VertexProxy::K->get(), VertexProxy::M->get(), mxREAL);
    pDouble2 = (double*) mxGetData(mx2);
    for (size_t v = this->tHyperRowOffset; v < tHyperRowOffset + VertexProxy::M->get(); ++v)    {
        const VertexProxy& wrapper = g.vertex_data(v);
        const THyperRow* const d = (const THyperRow* const) wrapper.getDelegator();
        memcpy(pDouble, &d->aValues[0], sizeof(element_type)*VertexProxy::K->get());
        memcpy(pDouble2, &d->bValues[0], sizeof(element_type)*VertexProxy::K->get());
        pDouble  += VertexProxy::K->get();
        pDouble2 += VertexProxy::K->get();
    }
    
    status = matPutVariable(pmat, "a_tm_transposed", mx);
    status = matPutVariable(pmat, "b_tm_transposed", mx2);
    
    matClose(pmat);
    
}