/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   WeightedAccuracy.cpp
 * Author: user
 * 
 * Created on July 9, 2019, 5:05 PM
 */

#include "WeightedAccuracy.hpp"
#include <stdio.h>
using namespace std;

WeightedAccuracy::WeightedAccuracy(Graph* G1, Graph* G2): Measure(G1, G2, "wacc") {
}

WeightedAccuracy::~WeightedAccuracy() {
}

void WeightedAccuracy::setAlpha(double anAlpha){
    alpha = anAlpha; 
}

double WeightedAccuracy::getAlpha(){
    return alpha;
}

double WeightedAccuracy::eval(const Alignment& A) {
    double Ea = A.numAlignedEdges(*G1, *G2);
    double E1 = G1->getNumEdges();
    
    double Ea_hat = G2->numNodeInducedSubgraphEdges(A.getMapping());
    double omega = E1 * (E1-1) / 2.0;
    double numerator = 2 * ((alpha + 1) * Ea + omega - (E1 + Ea_hat));
    double denomerator = 2 * omega + (alpha - 1) * (E1 + Ea_hat);
//    cout <<"E1: "<<E1 << " ,Ea: " << Ea << " ,Ea_hat: "<< Ea_hat<<" ,Omega: " << omega << endl;
    return numerator / denomerator;
}

