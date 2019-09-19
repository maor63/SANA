/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   IliyaMeasure.cpp
 * Author: user
 * 
 * Created on September 19, 2019, 9:17 AM
 */

#include "IliaMeasure.hpp"

IliaMeasure::IliaMeasure(Graph* G1, Graph* G2): Measure(G1, G2, "ilia") {
}

IliaMeasure::~IliaMeasure() {
}

double IliaMeasure::eval(const Alignment& A) {
    double Ea = A.numAlignedEdges(*G1, *G2);
    double E1 = G1->getNumEdges();
    double V1 = G1->getNumNodes();
    double omega = V1 * (V1-1) / 2.0;    
    double Ea_hat = G2->numNodeInducedSubgraphEdges(A.getMapping());
    double TP = Ea;
    double TN = omega - (E1 + Ea_hat - Ea);
    
    double positive_part = TP / E1;
    double negative_part = TN / (omega - E1);
    return (positive_part + negative_part) / 2.0;
}