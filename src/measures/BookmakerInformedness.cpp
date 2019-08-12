/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BookmakerInformedness.cpp
 * Author: user
 * 
 * Created on July 18, 2019, 4:52 PM
 */

#include "BookmakerInformedness.hpp"

BookmakerInformedness::BookmakerInformedness(Graph* G1, Graph* G2): Measure(G1, G2, "bm"){
}

BookmakerInformedness::~BookmakerInformedness() {
}

double BookmakerInformedness::eval(const Alignment& A) {
    double Ea = A.numAlignedEdges(*G1, *G2);
    double E1 = G1->getNumEdges();    
    double Ea_hat = G2->numNodeInducedSubgraphEdges(A.getMapping());
    double V1 = G1->getNumNodes();
    double omega = V1 * (V1-1) / 2.0;
    return Ea / Ea_hat - (E1 - Ea) / (omega - Ea_hat);
}

