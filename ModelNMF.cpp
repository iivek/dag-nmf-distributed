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
    mxDestroyArray(w_mx);
    mxDestroyArray(k_mx);
    mxDestroyArray(i_mx);
    
    this->initSharedData(m, n, k);
    VertexProxy::initSharedData(&(this->M), &(this->N), &(this->K));        
    
    
    mxArray* rowSupport_mx = matGetVariable(pmat, "rows_adjacency_colwise");
    unsigned int adjacencySparsity = mxGetM(rowSupport_mx); // number of nonzero entries    
    unsigned int* rowSupport = (unsigned int*) mxGetData(rowSupport_mx);    
    mxArray* colSupport_mx = matGetVariable(pmat, "cols_adjacency_colwise");
    
    assert(mxGetM(colSupport_mx) == adjacencySparsity);
    unsigned int* colSupport = (unsigned int*) mxGetData(colSupport_mx);
    /* Mapping of indices to get the original ordering of vertices (before
     * topological sorting)
     */
    
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
    }    
    mxDestroyArray(aVeAux_mx);

    /* First the gamma hyperparameters to have them at the beginning of the
     * neighbor list.
     */
    unsigned int shapeID;
    noparents = 0;
    for (unsigned int curr = 0; curr < numVertices; ++curr) {
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

                ScaleHyperColWrapper scaleVertex;
                ScaleHyperCol* scaleDelegator = (ScaleHyperCol*) scaleVertex.getDelegator();
                memcpy(&scaleDelegator->scaleRelatedValues[0], bVe, sizeof(element_type) * this->k);
                bVe += this->k;

                hyperIndex = g.add_vertex(scaleVertex);
                vHyperColIDs.push_back(hyperIndex);
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
                    mxDestroyArray(expectationsDiscreteCell_mx);
                    mxArray* dirichletHyperCell_mx = mxGetCell(uDirichlet_mx, uCounter++);
                    double* dirichletHyperContents = (double*) mxGetPr(dirichletHyperCell_mx);
                    DirichletHyperColWrapper dirichletHyperVertex( mxGetN(dirichletHyperCell_mx) );
                    DirichletHyperCol* dirichletHyperDelegator = (DirichletHyperCol*) dirichletHyperVertex.getDelegator();
                    memcpy(&dirichletHyperDelegator->values[0], dirichletHyperContents,
                            sizeof (element_type) * this->k * numParents[curr]);
                    dirichletHyperID = g.add_vertex(dirichletHyperVertex);
                    dirichletHyperColIDs.push_back(dirichletHyperID);                    
                    mxDestroyArray(dirichletHyperCell_mx);
                    
                    mxArray* expectationsDirichletCell_mx = mxGetCell(expectationsDirichlet_mx, curr);
                    double* expectationsDirichletContents = (double*) mxGetPr(expectationsDirichletCell_mx);
                    DirichletColWrapper dirichletVertex(mxGetN(expectationsDirichletCell_mx));
                    DirichletCol* dirichletDelegator = (DirichletCol*) dirichletVertex.getDelegator();
                    memcpy(&dirichletDelegator->values[0], expectationsDirichletContents,
                            sizeof(element_type)*this -> k*numParents[curr]);
                    dirichletDelegator->order = numParents[curr];
                    dirichletID = g.add_vertex(dirichletVertex);
                    dirichletColIDs.push_back(dirichletID);
                    mxDestroyArray(expectationsDirichletCell_mx);

                }// else it has only auxillary parents

            }// else {}

        }// else {}
    }
  
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
                vertexIndices[curr] = ourID;

                /* edges to existing parents */
                shapeID = shapeHyperColIDs[shapeIterator++];
                g.add_edge(shapeID, ourID, edata);
                g.add_edge(ourID, shapeID, edata);

                bool isMixture = false;
                
                if (numParents[curr] > 1) {
                    isMixture = true;
                    
                    discreteID = discreteColIDs[discreteIterator++];
                    dirichletID = dirichletColIDs[dirichletIterator++];
                    dirichletHyperID = dirichletHyperColIDs[dirichletHyperIterator++];

                    g.add_edge(discreteID, ourID, edata);
                    g.add_edge(ourID, discreteID, edata);
                    
                    /* Shape hyperparameter of currentVCol to Discrete parent */
                    g.add_edge(shapeID, discreteID, edata);

                    g.add_edge(dirichletHyperID, dirichletID, edata);
                    g.add_edge(dirichletID, dirichletHyperID, edata);
                           
                    g.add_edge(dirichletID, discreteID, edata);
                    g.add_edge(discreteID, dirichletID, edata);
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

                /* Now we're ready to fetch all curr's parents:  */
                for (unsigned int i = 0; i < numParents[curr]; ++i) {
                    // ID of current's parents in graphlab: vertexIndices[ rowSupport[location++] ]                                        

                    unsigned int parentID = vertexIndices[ rowSupport[location++] ];
                    AuxCol* parent = (AuxCol*) g.vertex_data(parentID).getDelegator();
                    parent->childIsMixture[ helper[parentID] ] = isMixture;
                    parent->indexInMixture[ helper[parentID] ] = i;
                    ++helper[parentID];

                    g.add_edge(parentID, ourID, edata);
                    g.add_edge(ourID, parentID, edata);
                    // parent's coparents; current's hyperparameters
                    g.add_edge(shapeID, parentID, edata);
                    // parent's coparents; current's discrete distro
                    if (numParents[curr] > 1) {
                        g.add_edge(discreteID, parentID, edata);
                        // ...and vice versa
                        g.add_edge(parentID, discreteID, edata);
                    }
                }
            } else {
                /* Current describes what should be embodied into VCol */
                VColWrapper vertex(false, false);       // no discrete, no gamma child
                VCol* delegator = (VCol*) vertex.getDelegator();
                ourID = g.add_vertex(vertex);
                vColIDs.push_back(ourID);
                vertexIndices[curr] = ourID;

                /* Adding edges */
                //shape
                unsigned int hyperIndex = vHyperColIDs[vIterator++];
                g.add_edge(hyperIndex, ourID, edata);
                g.add_edge(ourID, hyperIndex, edata);
                //scale
                hyperIndex = vHyperColIDs[vIterator++];
                g.add_edge(hyperIndex, ourID, edata);
                g.add_edge(ourID, hyperIndex, edata);
            }

        } else {
            /* current describes what should be an AuxCol */

            /* Adding vertex */
            AuxColWrapper vertex(numChildren[curr]);
            AuxCol* delegator = (AuxCol*) vertex.getDelegator();
            unsigned int ourID = g.add_vertex(vertex);
            auxColIDs.push_back(ourID);
            vertexIndices[curr] = ourID;
            helper[ourID] = 0;

            // we have only one parent - let's find it
            /* TODO: fix location, make it unsigned */
            
            long int location = binary_search(colSupport, 0, adjacencySparsity - 1, curr);
            unsigned int parentID = vertexIndices[ rowSupport[location] ];
            VColWrapper& vWrapper = (VColWrapper&) g.vertex_data(parentID);
            VCol* parentDelegator = (VCol*) vWrapper.getDelegator();
            parentDelegator->hasGammaChild = true;
            
            /* Adding edges */
            g.add_edge(ourID, parentID, edata);
            g.add_edge(parentID, ourID, edata);
            g.add_edge(auxHyperCol, ourID, edata);
            g.add_edge(ourID, auxHyperCol, edata);
            // and an edge from auxHyperCol to parentID (auxHyperCol is
            // coparent to parentID)
            g.add_edge(auxHyperCol, parentID, edata);
        }
        
    }
    delete(vertexIndices);

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
   }
    mxDestroyArray(aTm_mx);
    mxDestroyArray(bTm_mx);
          
    /* T, rowwise
     */
    mxArray* eT_mx = matGetVariable(pmat, "E_t_transposed");
    element_type* eT = (double*) mxGetData(eT_mx);
    mxArray* lT_mx = matGetVariable(pmat, "L_t_transposed");
    element_type* lT = (double*) mxGetData(lT_mx);

    tRowOffset = tHyperRowOffset + this->m;
    for (unsigned int row = 0; row < this->m; ++row) {
        TRowWrapper vertex(row);
        TRow* delegator = (TRow*) vertex.getDelegator();
        delegator->row = row;

        memcpy(&delegator->eValues[0], eT, sizeof(element_type) * this->k);
        memcpy(&delegator->expELogValues[0], lT, sizeof(element_type) * this->k);
        eT += this->k;
        lT += this->k;

        unsigned int bla = g.add_vertex(vertex); // == tRowOffset+row
        
        const TRowWrapper& vertex2 = (TRowWrapper&) g.vertex_data(bla);
        const TRow& delegator2 = (const TRow&) *vertex2.getDelegator();
    }
    
    mxDestroyArray(eT_mx);
    mxDestroyArray(lT_mx);  
    // M komada
        
    mxArray* colwiseNum_mx = matGetVariable(pmat, "colwise_num"); // number of nonzeros per column in the mask    
    unsigned int* colwiseNum = (unsigned int*) mxGetData(colwiseNum_mx);
    assert(mxGetN(colwiseNum_mx) == n);        
    
    /* tmp sparse, columnwise
     */
    tmpColOffset = tRowOffset + this->m;
    for (unsigned int col = 0; col < this->n; ++col) {        
        
        if (colwiseNum[col] > 0)        {
            TmpColWrapper vertex(colwiseNum[col]);
            TmpCol* delegator = (TmpCol*) vertex.getDelegator();           
            // values not initialized because of the order vertices get updated
            unsigned int bla =  g.add_vertex(vertex);
        } else {
            TmpColWrapper vertex(colwiseNum[col]);
            TmpCol* delegator = (TmpCol*) vertex.getDelegator();
        }
    }
    // N is their count

    /* SigT, rowwise
     */
    sigTRowOffset = tmpColOffset + this->n;
    for (unsigned int row = 0; row < this->m; ++row) {
        SigTRowWrapper vertex;
        SigTRow* delegator = (SigTRow*) vertex.getDelegator();
        unsigned int bla = g.add_vertex(vertex);
    }
    // M is their count

    /* SigV, colwise
     */
    sigVColOffset = sigTRowOffset + this->m;
    for (unsigned int col = 0; col < this->n; ++col) {
        SigVColWrapper vertex;
        SigVCol* delegator = (SigVCol*) vertex.getDelegator();
        unsigned int bla =
                g.add_vertex(vertex);
    }
    // N is the count


    /* processing X - original data matrix, columnwise
     */
    mxArray* colwiseX_mx = matGetVariable(pmat, "x_raw"); // mask of x
    double* colwiseX = (double*) mxGetData(colwiseX_mx);
    mxArray* colwiseMaskCols_mx = matGetVariable(pmat, "colwise_mask_cols"); // mask of x
    unsigned int* colwiseMaskCols = (unsigned int*) mxGetData(colwiseMaskCols_mx);
    size_t maskSparsity = mxGetM(colwiseMaskCols_mx); // number of nonzero entries
    mxArray* colwiseMaskRows_mx = matGetVariable(pmat, "colwise_mask_rows"); // mask of x
    unsigned int* colwiseMaskRows = (unsigned int*) mxGetData(colwiseMaskRows_mx);

    unsigned int* colwiseMaskRowsBreakPoints[this->n];
    xColOffset = sigVColOffset + this->n; // starting ID for both XCols and SupportXCols - every other is a SupportXCols;
    for (unsigned int col = 0; col < this->n; ++col) {
        if (colwiseNum[col] > 0) {
            long int location = binary_search(colwiseMaskCols, 0, maskSparsity - 1, col);
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
            
            SupportXColWrapper vertexW(colwiseNum[col]);
            SupportXCol* delegatorW = (SupportXCol*) vertexW.getDelegator();            
            // copy the supports for this column
            memcpy(&delegatorW->support[0], &colwiseMaskRows[location], sizeof(unsigned int)*delegatorW->getNumElements());
            
            bla = g.add_vertex(vertexW);
            //        }         // if(location != -1) ends
        } else {
            
            XColWrapper vertex(0);
            XCol* delegator = (XCol*) vertex.getDelegator();

            SupportXColWrapper vertexW(0);
            SupportXCol* delegatorW = (SupportXCol*) vertexW.getDelegator();    

            unsigned int bla =  g.add_vertex(vertex);
            bla = g.add_vertex(vertexW);
        }
    }
    // N is the count

    mxDestroyArray(colwiseX_mx);

    
    
    /********************************** Edges */


    for (unsigned int col = 0; col < this->n; ++col) {
        edge_data edata;
        unsigned int us = vColIDs[col];

        if (colwiseNum[col] > 0) {
            /* T to us - we need rows of T where x has columwise support */
            for (unsigned int supportIter = 0; supportIter < colwiseNum[col]; ++supportIter) {
                g.add_edge(tRowOffset + *(colwiseMaskRowsBreakPoints[col] + supportIter), us, edata);
            }
        }

        /* sigVCol to us */
        g.add_edge(sigVColOffset + col, us, edata);
        /* 
         * Us: tmpCol
         */
        us = tmpColOffset + col;
        /* XCol to us */
        g.add_edge(xColOffset + 2 * col, us, edata);
        if (colwiseNum[col] > 0) {
            /* T to us */
            for (unsigned int supportIter = 0; supportIter < colwiseNum[col]; ++supportIter) {
                g.add_edge(tRowOffset + *(colwiseMaskRowsBreakPoints[col] + supportIter), us, edata);
            }
            /* V to us */
            g.add_edge(vColIDs[col], us, edata);
        }

        /* 
         * Us: sigV 
         * We'll need a VCol, a TmpCol and some TRows
         */
        us = sigVColOffset + col;
        /* V to us */
        g.add_edge(vColIDs[col], us, edata);
        /* T to us    - we need rows of T where x has columwise support */
        if (colwiseNum[col] > 0) {
            for (unsigned int supportIter = 0; supportIter < colwiseNum[col]; ++supportIter) {
                g.add_edge(tRowOffset + *(colwiseMaskRowsBreakPoints[col] + supportIter), us, edata);
            }
        }
        /* col of tmp to us */
        g.add_edge(tmpColOffset + col, us, edata);
    }
    
    /* We'll also need the rowwise support. Unsorted is OK */
    vector< vector<unsigned int> > rowwiseSupport(this->m*this->n);
    unsigned int* pCol = colwiseMaskCols;
    unsigned int* pRow = colwiseMaskRows;
    for (unsigned int col = 0; col < this->n; ++col) {
        for (unsigned int rowIter = 0; rowIter < colwiseNum[col]; ++rowIter) {
            rowwiseSupport[ *(pRow++) ].push_back(*(pCol++));
        }
    }
    mxDestroyArray(colwiseNum_mx);  
    for (unsigned int row = 0; row < this->m; ++row) {
        rowwiseSupport[row].resize(rowwiseSupport[row].size(), 0);
    }

    for (unsigned int row = 0; row < this->m; ++row) {
        edge_data edata;
        /*
         * Us: T
         */
        unsigned int us = tRowOffset + row;
        /* V -> us:         
         * Determine which Vs we need using rowwiseSupport and
         * correspondingly add edges: */
        for (unsigned int supportIter = 0; supportIter < rowwiseSupport[row].size(); ++supportIter) {
            // Note that we don't need support explicitly.
            g.add_edge(vColIDs[ rowwiseSupport[row][supportIter] ], us, edata);
        }
        /* at sigTRow row us */
        g.add_edge(sigTRowOffset + row, us, edata);
        /* hyperparameter row to TRow */
        g.add_edge(tHyperRowOffset + row, us, edata);

        /* 
         * Us: sigT
         * We'll need a LtRow, TmpCols and LvCols
         * Our ID is sigTRowOffset+row
         */
        us = sigTRowOffset + row;
        /* T -> us */
        g.add_edge(tRowOffset + row, us, edata);
        /* V and tmp and support of x -> us
         * determine which Lvs we need using rowwiseSupport and correspondingly
         * add edges */
        for (unsigned int supportIter = 0; supportIter < rowwiseSupport[row].size(); ++supportIter) {
            g.add_edge(vColIDs[ rowwiseSupport[row][supportIter] ], us, edata);
            g.add_edge(tmpColOffset + rowwiseSupport[row][supportIter], us, edata);
            g.add_edge(xColOffset + 1 + 2 *rowwiseSupport[row][supportIter], us, edata);
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
        
    /*Markov expectations - right matrix & auxillary gammas*/
    /*shapes, scales*/           
    unsigned int shapeIndex = 0;
    if( forwardMapping.size() > N.get() )   {
        mx = mxCreateDoubleMatrix(VertexProxy::K->get(), 1, mxREAL);
        pDouble3 = (double*) mxGetData(mx);
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
    for (unsigned int curr = 0; curr < forwardMapping.size(); ++curr) {        
        if (forwardMapping[curr] < n)   {
            // VCol related
            if (numParents[curr] > 0)   {
                // shape
                const VertexProxy& wrapper = g.vertex_data(shapeHyperColIDs[shapeIndex]);
                const ShapeHyperCol* const d = (const ShapeHyperCol* const) wrapper.getDelegator();
                memcpy(pDouble3, &d->aValues[0], sizeof (element_type) * VertexProxy::K->get());
                pDouble3 += VertexProxy::K->get();
                ++shapeIndex;
            } else {
                // shape
                const VertexProxy& wrapper = g.vertex_data(vHyperColIDs[vHyperIndex]);
                const ShapeHyperCol* const d = (const ShapeHyperCol* const) wrapper.getDelegator();
                memcpy(pDouble3, &d->aValues[0], sizeof (element_type) * VertexProxy::K->get());
                pDouble3 += VertexProxy::K->get();
                ++vHyperIndex;
                //scale hyperparameter
                const VertexProxy& wrapper2 = g.vertex_data(vHyperColIDs[vHyperIndex]);
                const ScaleHyperCol* const d2 = (const ScaleHyperCol* const) wrapper2.getDelegator();                
                memcpy(pDouble4, &d2->scaleRelatedValues[0], sizeof(element_type)*VertexProxy::K->get());
                pDouble4 += VertexProxy::K->get();
                ++vHyperIndex;
            }
            const VertexProxy& wrapper = g.vertex_data(vColIDs[vIndex]);
            const VCol* const d = (const VCol* const) wrapper.getDelegator();
            memcpy(pDouble, &d->scaleRelatedValues[0], sizeof(element_type)*VertexProxy::K->get());
            memcpy(pDouble2, &d->expELogValues[0], sizeof(element_type)*VertexProxy::K->get());
            pDouble     += VertexProxy::K->get();
            pDouble2    += VertexProxy::K->get();
            ++vIndex;
        } else {
            // current describes what should be an AuxCol related
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