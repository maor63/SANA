/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   IliyaMeasure.hpp
 * Author: user
 *
 * Created on September 19, 2019, 9:17 AM
 */

#ifndef ILIAMEASURE_HPP
#define ILIAMEASURE_HPP
#include "Measure.hpp"

class IliaMeasure: public Measure {
public:
    IliaMeasure(Graph* G1, Graph* G2);
    virtual ~IliaMeasure();
    double eval(const Alignment& A);
private:

};

#endif /* ILIAMEASURE_HPP */

