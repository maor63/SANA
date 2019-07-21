/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   WeightedAccuracy.hpp
 * Author: user
 *
 * Created on July 9, 2019, 5:05 PM
 */

#ifndef WEIGHTEDACCURACY_HPP
#define WEIGHTEDACCURACY_HPP
#include "Measure.hpp"

class WeightedAccuracy: public Measure {
public:
    WeightedAccuracy(Graph* G1, Graph* G2);
    virtual ~WeightedAccuracy();
    double eval(const Alignment& A);
    
    void setAlpha(double anAlpha);
    double getAlpha();
    
private:
    double alpha = 1;
};

#endif /* WEIGHTEDACCURACY_HPP */

