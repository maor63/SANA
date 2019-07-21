/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MatthewsCorrelationCoefficient.hpp
 * Author: user
 *
 * Created on July 18, 2019, 4:51 PM
 */

#ifndef MATTHEWSCORRELATIONCOEFFICIENT_HPP
#define MATTHEWSCORRELATIONCOEFFICIENT_HPP
#include "Measure.hpp"

class MatthewsCorrelationCoefficient: public Measure {
public:
    MatthewsCorrelationCoefficient(Graph* G1, Graph* G2);
    virtual ~MatthewsCorrelationCoefficient();
    double eval(const Alignment& A);
private:

};

#endif /* MATTHEWSCORRELATIONCOEFFICIENT_HPP */

