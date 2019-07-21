/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MatthewsCorrelationCoefficient.cpp
 * Author: user
 * 
 * Created on July 18, 2019, 4:51 PM
 */

#include "MatthewsCorrelationCoefficient.hpp"

#include <math.h> 

MatthewsCorrelationCoefficient::MatthewsCorrelationCoefficient(Graph* G1, Graph* G2): Measure(G1, G2, "mcc") {
}

MatthewsCorrelationCoefficient::~MatthewsCorrelationCoefficient() {
}

double MatthewsCorrelationCoefficient::eval(const Alignment& A) {
    double Ea = A.numAlignedEdges(*G1, *G2);
    double E1 = G1->getNumEdges();
    
    double Ea_hat = G2->numNodeInducedSubgraphEdges(A.getMapping());
    double omega = E1 * (E1-1) / 2.0;
    double numerator = omega * Ea - E1 * Ea_hat;
    double denomerator = sqrt(E1 * Ea_hat * (omega - E1) * (omega - Ea_hat));
    return numerator / denomerator;
}

