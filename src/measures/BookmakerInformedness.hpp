/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BookmakerInformedness.hpp
 * Author: user
 *
 * Created on July 18, 2019, 4:52 PM
 */

#ifndef BOOKMAKERINFORMEDNESS_HPP
#define BOOKMAKERINFORMEDNESS_HPP
#include "Measure.hpp"

class BookmakerInformedness: public Measure{
public:
    BookmakerInformedness(Graph* G1, Graph* G2);
    virtual ~BookmakerInformedness();
    double eval(const Alignment& A);
private:

};

#endif /* BOOKMAKERINFORMEDNESS_HPP */

