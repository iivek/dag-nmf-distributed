/* 
 * File:   TRow.h
 * Author: ivek
 *
 * Created on January 17, 2012, 10:29 AM
 */

#ifndef TROW_H
#define	TROW_H

#include "UpdateFunctionDelegator.h"
//#include "MessagePasser.h"
#include <vector>

/* Wrapper related*/
#include "VertexProxy.h"


class MessageCollector;

class TRow : public UpdateFunctionDelegator
//        ,public MessagePasser
{
public:            
    TRow(unsigned int row);
    TRow(const TRow& orig);    
    ~TRow();
    UpdateFunctionDelegator* clone();
    /* our row index*/
    unsigned int row;
    /* Moments: expectation of gamma and expectation of log gamma */
    element_type* eValues;      /* Order: ascending indices */
    element_type* expELogValues;
//    std::vector<element_type> eValues;      /* Order: ascending indices */
//    std::vector<element_type> expELogValues;

    /* UpdateFunctionDelegator interface */
    void accept(VertexVisitor& v, gl::iscope& scope, gl::icallback& schedule);
    void updateFunction(gl::iscope& scope, gl::icallback& scheduler);
    /* UpdateFunctionDelegator ends */
    
private:
    TRow()      {}
};

class TRowWrapper : public VertexProxy {
public:    

    TRowWrapper(unsigned int row)       {        
        this->essence = new TRow(row);
    }
    TRowWrapper(const TRowWrapper& orig)      {        
        this->essence = orig.essence->clone();
    };
    ~TRowWrapper() {
        delete(this->essence);
    }
    
private:
    TRowWrapper()       {}
};

#endif	/* TROW_H */

