/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FBetaStar.cpp
 * Author: user
 * 
 * Created on July 21, 2019, 11:07 AM
 */

#include "FBetaStar.hpp"
#include <math.h> 

FBetaStar::FBetaStar(Graph* G1, Graph* G2): FBeta(G1, G2, "fbetastar") {
}

FBetaStar::~FBetaStar() {
}

void FBetaStar::setBetaA(Alignment trueA){
    double E1 = G1->getNumEdges();    
    double Ea_hat = G2->numNodeInducedSubgraphEdges(trueA.getMapping());    
    FBeta::setBeta(sqrt(E1 / Ea_hat));
    
}
