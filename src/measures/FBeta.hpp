/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FScore.hpp
 * Author: user
 *
 * Created on July 21, 2019, 10:16 AM
 */

#ifndef FBETA_HPP
#define FBETA_HPP
#include "Measure.hpp"

class FBeta: public Measure {
public:
    FBeta(Graph* G1, Graph* G2);    
    FBeta(Graph* G1, Graph* G2, string name);    
    virtual ~FBeta();
    double eval(const Alignment& A);
    virtual void setBeta(double beta);
    virtual double getBeta();
private:
    double beta = 1;
};

#endif /* FBETA_HPP */

