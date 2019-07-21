/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Markedness.hpp
 * Author: user
 *
 * Created on July 18, 2019, 4:52 PM
 */

#ifndef MARKEDNESS_HPP
#define MARKEDNESS_HPP
#include "Measure.hpp"

class Markedness: public Measure {
public:
    Markedness(Graph* G1, Graph* G2);
    virtual ~Markedness();
    double eval(const Alignment& A);
private:

};

#endif /* MARKEDNESS_HPP */

