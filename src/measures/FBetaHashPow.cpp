/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FBetaHashPow.cpp
 * Author: user
 * 
 * Created on September 23, 2019, 5:18 PM
 */

#include "FBetaHashPow.hpp"

FBetaHashPow::FBetaHashPow(Graph* G1, Graph* G2): FBeta(G1, G2, "fbetahashpow") {
    double E1 = G1->getNumEdges();
    double V1 = G1->getNumNodes();
    double omega = V1 * (V1-1) / 2.0;    
    FBeta::setBeta(E1 / (omega - E1));
}

FBetaHashPow::~FBetaHashPow() {
}

