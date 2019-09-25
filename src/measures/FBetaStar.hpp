/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FBetaStar.hpp
 * Author: user
 *
 * Created on July 21, 2019, 11:07 AM
 */

#ifndef FBETASTAR_HPP
#define FBETASTAR_HPP
#include "FBeta.hpp"

class FBetaStar: public FBeta {
public:
    FBetaStar(Graph* G1, Graph* G2);
    virtual ~FBetaStar();
    void setBetaA(Alignment trueA);
private:

};

#endif /* FBETASTAR_HPP */

