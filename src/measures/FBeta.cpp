/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FScore.cpp
 * Author: user
 * 
 * Created on July 21, 2019, 10:16 AM
 */

#include "FBeta.hpp"

FBeta::FBeta(Graph* G1, Graph* G2): Measure(G1, G2, "fbeta") {
    
}

FBeta::FBeta(Graph* G1, Graph* G2, string name): Measure(G1, G2, name) {
    
}

FBeta::~FBeta() {
}

double FBeta::getBeta(){
    return beta;
}

void FBeta::setBeta(double beta){
    this->beta = beta;
}

double FBeta::eval(const Alignment& A) {
    double Ea = A.numAlignedEdges(*G1, *G2);
    double E1 = G1->getNumEdges();    
    double Ea_hat = G2->numNodeInducedSubgraphEdges(A.getMapping());
//    double omega = E1 * (E1-1) / 2.0;    
    return ((1 + beta * beta) * Ea) / (E1 + beta * beta * Ea_hat);
}

