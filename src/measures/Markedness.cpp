/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Markedness.cpp
 * Author: user
 * 
 * Created on July 18, 2019, 4:52 PM
 */

#include "Markedness.hpp"

Markedness::Markedness(Graph* G1, Graph* G2): Measure(G1, G2, "mk") {
}

Markedness::~Markedness() {
}

double Markedness::eval(const Alignment& A) {
    double Ea = A.numAlignedEdges(*G1, *G2);
    double E1 = G1->getNumEdges();
    
    double Ea_hat = G2->numNodeInducedSubgraphEdges(A.getMapping());
    double omega = E1 * (E1-1) / 2.0;
    return Ea / E1 - (Ea_hat - Ea) / (omega - E1);
}

