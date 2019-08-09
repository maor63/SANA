/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FBetaHash.cpp
 * Author: user
 * 
 * Created on August 8, 2019, 6:31 PM
 */

#include "FBetaHash.hpp"
#include <math.h>

FBetaHash::FBetaHash(Graph* G1, Graph* G2): FBeta(G1, G2, "fbetahash") {
    double E1 = G1->getNumEdges();
    double omega = E1 * (E1-1) / 2.0;    
    FBeta::setBeta(sqrt(E1 / (omega - E1)));
}

FBetaHash::~FBetaHash() {
}

