#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdexcept>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <queue>
#include <set>
#include <cmath>
#include <fstream>
#include <time.h>
#include <cassert>
#include "Graph.hpp"
#include "utils/Timer.hpp"
#include "computeGraphlets.hpp"
#include "utils/utils.hpp"
using namespace std;

Graph Graph::loadGraph(string name) {
    Graph g;
    g.loadGwFile("networks/"+name+"/"+name+".gw");
    g.name = name;
    return g;
}

//transform format
void Graph::edgeList2gw(string fin, string fout) {
    vector<string> nodes = removeDuplicates(fileToStrings(fin));
    map<string,uint> nodeName2IndexMap;
    uint numNodes = nodes.size();
    for (uint i = 0; i < numNodes; i++) {
        nodeName2IndexMap[nodes[i]] = i;
    }
    vector<vector<string> > edges = fileToStringsByLines(fin);
    vector<vector<ushort>> edgeList(edges.size(), vector<ushort> (2));

    for (uint i = 0; i < edges.size(); i++) {
        if (edges[i].size() != 2) {
            throw runtime_error("File not in edge-list format: "+fin);
        }
        string node1 = edges[i][0];
        string node2 = edges[i][1];
        uint index1 = nodeName2IndexMap[node1];
        uint index2 = nodeName2IndexMap[node2];
        edgeList[i][0] = index1;
        edgeList[i][1] = index2;
    }
    saveInGWFormat(fout, nodes, edgeList);
}

string Graph::getName() const {
    return name;
}

Graph::Graph() :
    edgeList(vector<vector<ushort> > (0)),
    adjMatrix(vector<vector<bool> > (0)),
    adjLists(vector<vector<ushort> > (0)),
    decimal_representation(0),
    lockedList(vector<bool> (0)),
    lockedTo(vector<string>(0)),
    lockedCount(0)
    {}

Graph::Graph(const Graph& G) {
    edgeList = vector<vector<ushort> > (G.edgeList);
    adjMatrix = vector<vector<bool> > (G.adjMatrix);
    adjLists = vector<vector<ushort> > (G.adjLists);
    decimal_representation = G.decimal_representation;
    connectedComponents = vector<vector<ushort> > (G.connectedComponents);
    lockedList = vector<bool> (G.lockedList);
    lockedTo = vector<string> (G.lockedTo);
    lockedCount = G.lockedCount;
}

Graph::Graph(uint n, const vector<vector<ushort> > edges) {
    adjLists = vector<vector<ushort> > (n, vector<ushort> (0));
    adjMatrix = vector<vector<bool> > (n, vector<bool> (n, false));
    edgeList = edges;

    lockedList = vector<bool> (n, false);
    lockedTo = vector<string> (n, "");
    lockedCount = 0;

    //only add edges preserved by alignment
    for (const auto& edge: edges) {
        ushort node1 = edge[0], node2 = edge[1];
        adjLists[node1].push_back(node2);
        adjLists[node2].push_back(node1);
        adjMatrix[node1][node2] = adjMatrix[node2][node1] = true;
    }
    initConnectedComponents();
    construct_decimal_representation();
}

uint Graph::getNumNodes() const {
    return adjLists.size();
}

uint Graph::getNumEdges() const {
    return edgeList.size();
}

uint Graph::getNumConnectedComponents() const {
    return connectedComponents.size();
}

void Graph::getAdjMatrix(vector<vector<bool> >& adjMatrixCopy) const {
    uint n = getNumNodes();
    adjMatrixCopy = vector<vector<bool> > (n, vector<bool> (n));
    for (uint i = 0; i < n; i++) {
        for (uint j = 0; j < n; j++) adjMatrixCopy[i][j] = adjMatrix[i][j];
    }
}

void Graph::getAdjLists(vector<vector<ushort> >& adjListsCopy) const {
    uint n = getNumNodes();
    adjListsCopy = vector<vector<ushort> > (n, vector<ushort> (0));
    for (uint i = 0; i < n; i++) {
        adjListsCopy[i] = adjLists[i];
    }
}

void Graph::getEdgeList(vector<vector<ushort> >& edgeListCopy) const {
    edgeListCopy = edgeList;
}

const vector<vector<ushort> >& Graph::getConnectedComponents() const {
    return connectedComponents;
}

void Graph::loadGwFile(const string& fileName) {
    //this function could be improved to deal with blank lines and comments
    stringstream errorMsg;

    ifstream infile(fileName.c_str());
    string line;
    //ignore header
    for (int i = 0; i < 4; i++) getline(infile, line);
    //read number of nodes
    int n;
    getline(infile, line);
    istringstream iss(line);
    if (!(iss >> n) or n <= 0) {
        errorMsg << "Failed to read node number: " << line;
        throw runtime_error(errorMsg.str().c_str());
    }
    //read (and ditch) nodes
    string node;
    for (int i = 0; i < n; i++) {
        getline(infile, line);
        istringstream iss(line);
        if (!(iss >> node)) {
            errorMsg << "Failed to read node " << i << " of " << n << ": " << line << " (" << node << ")";
            throw runtime_error(errorMsg.str().c_str());
        }
    }
    //read number of edges
    int m;
    getline(infile, line);
    istringstream iss2(line);
    if (!(iss2 >> m)) {
        errorMsg << "Failed to read edge number: " << line;
        throw runtime_error(errorMsg.str().c_str());
    }

    adjLists = vector<vector<ushort> > (n, vector<ushort>(0));
    adjMatrix = vector<vector<bool> > (n, vector<bool>(n, false));
    edgeList = vector<vector<ushort> > (m, vector<ushort>(2));
    lockedList = vector<bool> (n, false);
    lockedTo = vector<string> (n, "");

    //read edges
    for (int i = 0; i < m; i++) {
        getline(infile, line);
        istringstream iss(line);
        ushort node1, node2;
        if (!(iss >> node1 >> node2)) {
            errorMsg << "Failed to read edge: " << line;
            throw runtime_error(errorMsg.str().c_str());
        }
        node1--; node2--; //-1 because of remapping

        edgeList[i][0] = node1;
        edgeList[i][1] = node2;

        adjMatrix[node1][node2] = adjMatrix[node2][node1] = true;
        adjLists[node1].push_back(node2);
        adjLists[node2].push_back(node1);
    }
    infile.close();
    initConnectedComponents();
}

bool comp_vectors(const vector<ushort>* a,const vector<ushort>* b) {
   return a->size() > b->size();
}

void Graph::initConnectedComponents() {
    //this function takes most of the initialization time
    ushort n = getNumNodes();
    vector<vector<ushort>* > aux(0);
    unordered_set<ushort> nodes;
    for (ushort i = 0; i < n; i++) nodes.insert(i);
    while (nodes.size() > 0) {
        ushort u = *nodes.begin();
        nodes.erase(nodes.begin());
        unordered_set<ushort> group;
        group.insert(u);
        queue<ushort> Q;
        Q.push(u);
        while (not Q.empty()) {
            ushort v = Q.front();
            Q.pop();
            unordered_set<ushort> vNeighbors(adjLists[v].begin(), adjLists[v].end());
            for (const auto& node: group) vNeighbors.erase(node);
            for (const auto& node: vNeighbors) nodes.erase(node);
            group.insert(vNeighbors.begin(), vNeighbors.end());
            for (const auto& node: vNeighbors) Q.push(node);
        }
        aux.push_back(new vector<ushort> (group.begin(), group.end()));
    }
    sort(aux.begin(), aux.end(), comp_vectors);
    connectedComponents = vector<vector<ushort> > (aux.size(), vector<ushort> (0));
    for (uint i = 0; i < connectedComponents.size(); i++) {
        connectedComponents[i] = vector<ushort> (aux[i]->size());
        for (uint j = 0; j < aux[i]->size(); j++) {
            connectedComponents[i][j] = (*aux[i])[j];
        }
    }
}

uint Graph::numNodeInducedSubgraphEdges(const vector<ushort>& subgraphNodes) const {
    unordered_set<ushort> nodeSet(subgraphNodes.begin(), subgraphNodes.end());
    uint count = 0;
    for (uint i = 0; i < subgraphNodes.size(); i++) {
        ushort node1 = subgraphNodes[i];
        for (uint j = 0; j < adjLists[node1].size(); j++) {
            ushort node2 = adjLists[node1][j];
            count += nodeSet.count(node2);
        }
    }
    return count/2;
}

Graph Graph::nodeInducedSubgraph(const vector<ushort>& nodes) const {
    uint n = nodes.size();
    vector<ushort> rev = reverseMapping(nodes, getNumNodes());
    unordered_set<ushort> nodeSet(nodes.begin(), nodes.end());
    Graph G;
    G.adjLists = vector<vector<ushort> > (n, vector<ushort> (0));
    G.adjMatrix = vector<vector<bool> > (n, vector<bool> (n, false));
    //only add edges between induced nodes
    for (const auto& edge: edgeList) {
        ushort node1 = edge[0], node2 = edge[1];
        if (nodeSet.count(node1) and nodeSet.count(node2)) {
            ushort newNode1 = rev[node1];
            ushort newNode2 = rev[node2];
            G.adjLists[newNode1].push_back(newNode2);
            G.adjLists[newNode2].push_back(newNode1);
            vector<ushort> newEdge(2);
            newEdge[0] = newNode1;
            newEdge[1] = newNode2;
            G.edgeList.push_back(newEdge);
            G.adjMatrix[newNode1][newNode2] = G.adjMatrix[newNode2][newNode1] = true;
        }
    }
    G.initConnectedComponents();
    return G;
}

void Graph::printStats(int numConnectedComponentsToPrint, ostream& stream) const {
    stream << "n    = " << getNumNodes() << endl;
    stream << "m    = " << getNumEdges() << endl;
    stream << "#connectedComponents = " << getNumConnectedComponents() << endl;
    stream << "Largest connectedComponents (nodes, edges) = ";
    for (int i = 0; i < min(numConnectedComponentsToPrint, getNumConnectedComponents()); i++) {
        const vector<ushort>& nodes = getConnectedComponents()[i];
        Graph H = nodeInducedSubgraph(nodes);
        stream << "(" << H.getNumNodes() << ", " << H.getNumEdges() << ") ";
    }
    stream << endl;
}

void Graph::getDistanceMatrix(vector<vector<short> >& dist) const {
    string distMatrixFile = "networks/"+name+"/autogenerated/"+name+"_distMatrix.txt";
    if (fileExists(distMatrixFile)) {
        readMatrixFromFile(dist, getNumNodes(), getNumNodes(), distMatrixFile);
    }
    else {
        Timer T;
        T.start();
        cerr << "Computing "+name+" distance matrix...";
        computeDistanceMatrix(dist);
        cerr << "done (" << T.elapsedString() << ")" << endl;
        writeMatrixToFile(dist, distMatrixFile);
    }
}

//a -1 indicates that the distance is infinite
//floyd-warshall algorithm
void Graph::computeDistanceMatrix(vector<vector<short> >& dist) const {
    ushort n = getNumNodes();
    dist = vector<vector<short> > (n, vector<short> (n, -1));
    for (ushort v = 0; v < n; v++) {
        dist[v][v] = 0;
        for (ushort i = 0; i < adjLists[v].size(); i++) {
            ushort u = adjLists[v][i];
            dist[u][v] = 1;
        }
    }
    for (ushort k = 0; k < n; k++) {
        for (ushort i = 0; i < n; i++) {
            for (ushort j = 0; j < n; j++) {
                if (dist[i][k] != -1 and dist[k][j] != -1) {
                    if (dist[i][j] == -1 or dist[i][j] > dist[i][k] + dist[k][j]) {
                        dist[i][j] = dist[i][k] + dist[k][j];
                    }
                }
            }
        }
    }
}

//The first line contains two integers n and e - the number of nodes and edges.
//The following e lines describe undirected edges with space-separated ids of their endpoints.
//Node ids should be between 0 and n-1
void Graph::writeGraphEdgeListFormat(const string& fileName) {
    ushort n = getNumNodes();
    ushort e = getNumEdges();
    ofstream outfile;
    outfile.open(fileName.c_str());
    outfile << n << " " << e << endl;
    for (uint i = 0; i < e; i++) {
        outfile << edgeList[i][0] << " " << edgeList[i][1] << endl;
    }
    outfile.close();
}

void Graph::writeGraphEdgeListFormatNETAL(const string& fileName) {
    ushort e = getNumEdges();
    ofstream outfile;
    outfile.open(fileName.c_str());
    for (uint i = 0; i < e; i++) {
        outfile << edgeList[i][0] << "\t" << edgeList[i][1] << endl;
    }
    outfile.close();
}

vector<ushort> Graph::numNodesAround(ushort node, ushort maxDist) const {
    uint n = getNumNodes();
    vector<ushort> distances(n, n);
    distances[node] = 0;
    queue<ushort> Q;
    Q.push(node);
    while (not Q.empty()) {
        ushort u = Q.front();
        Q.pop();
        ushort dist = distances[u];
        if (dist == maxDist) break;
        for (uint i = 0; i < adjLists[u].size(); i++) {
            ushort v = adjLists[u][i];
            if (distances[v] < n) continue;
            distances[v] = dist+1;
            Q.push(v);
        }
    }
    vector<ushort> result(maxDist, 0);
    for (uint i = 0; i < n; i++) {
        if (distances[i] < n and distances[i] > 0) result[distances[i]-1]++;
    }
    return result;
}

vector<ushort> Graph::numEdgesAround(ushort node, ushort maxDist) const {
    uint n = getNumNodes();
    vector<ushort> distances(n, n);
    vector<bool> visited(n, false);
    distances[node] = 0;
    queue<ushort> Q;
    Q.push(node);
    vector<ushort> result(maxDist, 0);
    while (not Q.empty()) {
        ushort u = Q.front();
        Q.pop();
        ushort dist = distances[u];
        if (dist == maxDist) break;
        for (uint i = 0; i < adjLists[u].size(); i++) {
            ushort v = adjLists[u][i];
            if (not visited[v]) result[dist]++;
            if (distances[v] < n) continue;
            distances[v] = dist+1;
            Q.push(v);
        }
        visited[u] = true;
    }
    return result;
}

ushort Graph::randomNode() {
    return randInt(0, getNumNodes()-1);
}

//note: does not update CCs
void Graph::addEdge(ushort node1, ushort node2) 
{
    adjMatrix[node1][node2] = adjMatrix[node2][node1] = true;
    vector<ushort> edge(2);
    edge[0] = node1;
    edge[1] = node2;
    edgeList.push_back(edge);
    adjLists[node1].push_back(node2);
    adjLists[node2].push_back(node1);
}

//note: does not update CCs
void Graph::removeEdge(ushort node1, ushort node2) {
    adjMatrix[node1][node2] = adjMatrix[node2][node1] = false;
    uint m = getNumEdges();
    //update edge list
    for (uint i = 0; i < m; i++) {
        if ((edgeList[i][0] == node1 and edgeList[i][1] == node2) or
            (edgeList[i][0] == node2 and edgeList[i][1] == node1)) {
            vector<ushort> lastEdge = edgeList[m-1];
            edgeList[i] = lastEdge;
            edgeList.pop_back();
            break;
        }
    }
    //update adjacency lists
    for (uint i = 0; i < adjLists[node1].size(); i++) {
        if (adjLists[node1][i] == node2) {
            ushort lastNeighbor = adjLists[node1][adjLists[node1].size()-1];
            adjLists[node1][i] = lastNeighbor;
            adjLists[node1].pop_back();
            break;
        }
    }
    for (uint i = 0; i < adjLists[node2].size(); i++) {
        if (adjLists[node2][i] == node1) {
            ushort lastNeighbor = adjLists[node2][adjLists[node2].size()-1];
            adjLists[node2][i] = lastNeighbor;
            adjLists[node2].pop_back();
            break;
        }
    }
}

//note: does not update CCs
void Graph::addRandomEdge() {
    ushort node1 = 0, node2 = 0;
    while (node1 == node2 or adjMatrix[node1][node2]) {
        node1 = randomNode();
        node2 = randomNode();
    }
    addEdge(node1, node2);
}

//note: does not update CCs
void Graph::removeRandomEdge() {
    ushort node1 = 0, node2 = 0;
    while (node1 == node2 or not adjMatrix[node1][node2]) {
        node1 = randomNode();
        node2 = randomNode();
    }
    removeEdge(node1, node2);
}

void Graph::addRandomEdges(double addedEdgesProportion) {
    uint n = (double) getNumEdges() * addedEdgesProportion;
    for (uint i = 0; i <= n; i++) {
        addRandomEdge();
    }
    initConnectedComponents();
}

void Graph::removeRandomEdges(double removedEdgesProportion) {
    uint n = (double) getNumEdges() * removedEdgesProportion;
    for (uint i = 0; i <= n; i++) {
        removeRandomEdge();
    }
    initConnectedComponents();
}

void Graph::rewireRandomEdges(double rewiredEdgesProportion) {
    uint n = (double) getNumEdges() * rewiredEdgesProportion;
    for (uint i = 0; i <= n; i++) {
        addRandomEdge();
        removeRandomEdge();
    }
    initConnectedComponents();
}

string Graph::autogenFilesFolder() {
    return "networks/"+name+"/autogenerated/";
}

vector<vector<uint> > Graph::loadGraphletDegreeVectors() {
    string gdvsFileName = autogenFilesFolder() + name + "_gdv.bin";
    uint n = getNumNodes();
    if (fileExists(gdvsFileName)) {
        vector<vector<uint> > gdvs(n, vector<uint> (73));
        readMatrixFromBinaryFile(gdvs, gdvsFileName);
        return gdvs;
    }
    cerr << "Computing " << gdvsFileName << " ... ";
    Timer T;
    T.start();
    vector<vector<uint> > gdvs = computeGraphletDegreeVectors();
    cerr << "done (" << T.elapsedString() << ")" << endl;
    writeMatrixToBinaryFile(gdvs, gdvsFileName);
    string readeableVersionFile = autogenFilesFolder() + name + "_gdv.txt";
    writeMatrixToFile(gdvs, readeableVersionFile);
    return gdvs;
}

vector<vector<uint> > Graph::computeGraphletDegreeVectors() {
    uint n = getNumNodes();
    uint m = getNumEdges();

    string fileName = "tmp/compute_dgvs.txt";
    ofstream fout(fileName.c_str());
    fout << n << " " << m << endl;
    for (uint i = 0; i < m; i++) {
        fout << edgeList[i][0] << " " << edgeList[i][1] << endl;
    }
    fout.close();

    vector<vector<uint> > gdvs = computeGraphlets(5, fileName);
    return gdvs;
}

map<string,ushort> Graph::getNodeNameToIndexMap() const {
    string networkFile = "networks/"+name+"/"+name+".gw";

    ifstream infile(networkFile.c_str());
    string line;
    //ignore header
    for (int i = 0; i < 4; i++) getline(infile, line);
    //read number of nodes
    int n;
    getline(infile, line);
    istringstream iss(line);
    if (!(iss >> n) or n <= 0) {
        throw runtime_error("Failed to read node number: " + line);
    }

    //read nodes
    map<string, ushort> res;
    string node;
    for (ushort i = 0; i < (ushort) n; i++) {
        getline(infile, line);
        istringstream iss(line);
        if (!(iss >> node)) {
            throw runtime_error("Failed to read node "+intToString(i)+" of "+intToString(n)+": "+line+" ("+node+")");
        }
        node = node.substr(2,node.size()-4); //strip |{ and }|
        res[node] = i;
    }
    infile.close();
    return res;
}

map<ushort,string> Graph::getIndexToNodeNameMap() const {
    map<string,ushort> reverse = getNodeNameToIndexMap();
    map<ushort,string> res;
    for (const auto &nameIndexPair : reverse ) {
        res[nameIndexPair.second] = nameIndexPair.first;
    }
    return res;
}

vector<string> Graph::getNodeNames() const {
    vector<string> result(getNumNodes());
    map<ushort,string> map = getIndexToNodeNameMap();
    for (ushort i = 0; i < getNumNodes(); i++) {
        result[i] = map[i];
    }
    return result;
}

vector<uint> Graph::degreeDistribution() const {
    uint n = getNumNodes();
    uint maxDegree = adjLists[0].size();
    for (uint i = 1; i < n; i++) {
        if (adjLists[i].size() > maxDegree) maxDegree = adjLists[i].size();
    }
    vector<uint> res(maxDegree+1);
    for (uint i = 0; i < n; i++) {
        res[adjLists[i].size()]++;
    }
    return res;
}

double Graph::getAverageDistance() const {
    vector<vector<short> > dists;
    getDistanceMatrix(dists);
    uint n = getNumNodes();
    double distSum = 0;
    double distCount = 0;
    for (uint i = 0; i < n; i++) {
        for (uint j = i+1; j < n; j++) {
            if (dists[i][j] > 0) {
                distSum += dists[i][j];
                distCount++;
            }
        }
    }
    return distSum/distCount;
}

struct Point {
    double x;
    double y;
    double z;

    Point(double x, double y, double z): x(x), y(y), z(z) {}

    double distOrigin() {
        return sqrt(x*x + y*y + z*z);
    }

    void normalize() {
        double l = distOrigin();
        x/=l; y/=l; z/=l;
    }

    static Point randomNormalizedPoint() {
        Point p(randDouble(), randDouble(), randDouble());
        p.normalize();
        return p;
    }

    double dist(const Point& p) {
        double a = x-p.x;
        double b = y-p.y;
        double c = z-p.z;
        return sqrt(a*a+b*b+c*c);
    }
};

struct Sphere {
    Point center;
    double radius;
    Sphere(Point center, double radius): center(center), radius(radius) {}

    Point randomContainedPoint() {
        Point p = Point::randomNormalizedPoint();
        double dist = radius*pow(randDouble(),1./3.);
        p.x *= dist; p.y *= dist; p.z *= dist;
        p.x += center.x; p.y += center.y; p.z += center.z;
        return p;
    }
};

//model GEO-GD expansion from paper:
//GEOMETRIC EVOLUTIONARY DYNAMICS OF PROTEIN INTERACTION NETWORKS
void Graph::GeoGeneDuplicationModel(uint numNodes, uint numEdges, string outputFile) {
    if (numNodes < 5) cerr << "The minimum number of nodes is 5";
    const double epsilon = 1;
    //initial network
    Point center = Point(0,0,0);
    Sphere sphere(center, epsilon/2);
    vector<Point> points(5, sphere.randomContainedPoint());
    //add nodes
    for (uint i = 5; i < numNodes; i++) {
        uint parentIndex = randMod(i);
        Point p = points[parentIndex];
        Point dir = Point::randomNormalizedPoint();
        double dist = randDouble()*2*epsilon;
        p.x += dist*dir.x;
        p.y += dist*dir.y;
        p.z += dist*dir.z;
        points.push_back(p);
    }
    //add edges
    long long unsigned int numDists = (numNodes*(numNodes-1))/2;
    vector<double> dists(numDists, -1);
    long long unsigned int count = 0;
    for (uint i = 0; i < numNodes; i++) {
        for (uint j = 0; j < i; j++) {
            dists[count] = points[i].dist(points[j]);
            count++;
        }
    }
    sort(dists.begin(), dists.end());
    double cutOffDist = dists[numEdges];

    ofstream outfile;
    outfile.open(outputFile.c_str());
    //header
    outfile << "LEDA.GRAPH" << endl;
    outfile << "string" << endl;
    outfile << "short" << endl;
    outfile << "-2" << endl;
    outfile << numNodes << endl;
    for (uint i = 0; i < numNodes; i++) {
        outfile << "|{node" << i+1 << "}|" << endl;
    }
    outfile << numEdges << endl;
    for (uint i = 0; i < numNodes; i++) {
        for (uint j = 0; j < i; j++) {
            if (points[i].dist(points[j]) < cutOffDist) {
                outfile << j+1 << " " << i+1 << " |{}|" << endl;
            }
        }
    }
    outfile << endl;
}

void writeGWHeader(ofstream& outfile) {
    outfile << "LEDA.GRAPH" << endl;
    outfile << "string" << endl;
    outfile << "short" << endl;
    outfile << "-2" << endl;
}

void writeGWNodes(ofstream& outfile, const vector<string>& nodeNames) {
    uint numNodes = nodeNames.size();
    outfile << numNodes << endl;
    for (uint i = 0; i < numNodes; i++) {
        outfile << "|{" << nodeNames[i] << "}|" << endl;
    }
}

void writeGWEdges(ofstream& outfile, const vector<vector<ushort>>& edgeList) {
    uint numEdges = edgeList.size();
    outfile << numEdges << endl;
    for (uint i = 0; i < numEdges; i++) {
        outfile << edgeList[i][0]+1 << " " << edgeList[i][1]+1 << " 0 |{}|" << endl;
    }
}

void Graph::saveInGWFormat(string outputFile, const vector<string>& nodeNames,
    const vector<vector<ushort>>& edgeList) {
    ofstream outfile;
    outfile.open(outputFile.c_str());

    writeGWHeader(outfile);
    writeGWNodes(outfile, nodeNames);
    writeGWEdges(outfile, edgeList);
}

void Graph::saveInGWFormatShuffled(string outputFile, const vector<string>& nodeNames,
    const vector<vector<ushort>>& edgeList) {

    uint n = nodeNames.size();
    vector<ushort> origPos2NewPos(n);
    for (uint i = 0; i < n; i++) origPos2NewPos[i] = i;
    randomShuffle(origPos2NewPos);

    map<ushort, ushort> newPos2OrigPos;
    for (uint i = 0; i < n; i++) {
        newPos2OrigPos[origPos2NewPos[i]] = i;
    }

    vector<string> newNodeNames(n);
    for (uint i = 0; i < n; i++) {
        newNodeNames[i] = nodeNames[newPos2OrigPos[i]];
    }

    uint m = edgeList.size();
    vector<vector<ushort>> newEdgeList(m, vector<ushort> (2));
    for (uint i = 0; i < m; i++) {
        newEdgeList[i][0] = origPos2NewPos[edgeList[i][0]];
        newEdgeList[i][1] = origPos2NewPos[edgeList[i][1]];
    }
    saveInGWFormat(outputFile, newNodeNames, newEdgeList);
}

//assumes the gw file for this graph does not already exist
//so we create new node names
void Graph::saveInGWFormat(string outputFile) {
    uint numNodes = getNumNodes();
    vector<string> nodeNames(numNodes);
    for (uint i = 0; i < numNodes; i++) {
        nodeNames[i] = "node"+intToString(i+1);
    }
    saveInGWFormat(outputFile, nodeNames, edgeList);
}

void Graph::saveInGWFormatWithNames(string outputFile) {
    saveInGWFormat(outputFile, getNodeNames(), edgeList);
}

void Graph::saveInShuffledOrder(string outputFile) {
    uint numNodes = getNumNodes();
    map<ushort,string> index2Name = getIndexToNodeNameMap();
    vector<string> nodeNames(numNodes);
    for (uint i = 0; i < numNodes; i++) {
        nodeNames[i] = index2Name[i];
    }
    saveInGWFormatShuffled(outputFile, nodeNames, edgeList);
}

void Graph::saveGraphletsAsSigs(string outputFile) {
    std::vector<std::vector<uint>> graphlets = computeGraphletDegreeVectors();
    vector<string> nodeNames = getNodeNames();
    ofstream out;
    out.open(outputFile.c_str());
     for (unsigned int i = 0; i < nodeNames.size(); i++) {
         out << nodeNames[i] << "\t";
         for (int j = 0; j < 73; j++) {
             out << graphlets[i][j] << "\t";
         }
         out << std::endl;
     }

     out.close();
}

Graph Graph::randomNodeInducedSubgraph(uint numNodes) {
    uint n = getNumNodes();
    if (numNodes > n) cerr << "the subgraph cannot have more nodes" << endl;
    vector<ushort> v(getNumNodes());
    for (uint i = 0; i < getNumNodes(); i++) v[i] = i;
    randomShuffle(v);
    v = vector<ushort> (v.begin(), v.begin()+numNodes);
    return nodeInducedSubgraph(v);
}

bool Graph::isWellDefined() {
    bool isWellDefined = true;
    //check that adj lists don't have repeated entries
    //check that all indices in adj lists are within bounds
    //check that every edge in the adj lists appears in the adj matrix
    //check that no node is neighbor to itself in the adj lists
    uint numNodes = adjMatrix.size();
    for (uint i = 0; i < adjLists.size(); i++) {
        for (uint j = 0; j < adjLists[i].size(); j++) {
            uint neighbor = adjLists[i][j];
            if (neighbor < 0 or neighbor >= numNodes) {
                cerr << "Adjacency list ill defined: node " << i << " adjacent to ";
                cerr << "node out of range " << neighbor << endl;
                isWellDefined = false;
            }
            if (neighbor == i) {
                cerr << "Adjacency list ill defined: node " << i << " adjacent to itself" << endl;
                isWellDefined = false;

            }
            if (not adjMatrix[i][neighbor] or not adjMatrix[i][neighbor]) {
                cerr << "node " << i << " adjacent to node " << adjMatrix[i][neighbor];
                cerr << " in the adj lists but not on the adj matrix" << endl;
                isWellDefined = false;
            }
            for (uint k = j+1; k < adjLists[i].size(); k++) {
                if (neighbor == adjLists[i][k]) {
                    cerr << "Adjacency list ill defined: node " << i << " adjacent to ";
                    cerr << neighbor << " twice" << endl;
                    isWellDefined = false;
                }
            }
        }
    }
    //check that no node is adj to itself in the adj matrix
    for (uint i = 0; i < numNodes; i++) {
        if (adjMatrix[i][i]) {
            cerr << "Adjacency matrix ill defined: node " << i << " adjacent to itself" << endl;
            isWellDefined = false;
        }
    }
    //check that adj matrix is symmetric
    //check that every edge in the adj matrix appears in the adj lists (twice)
    for (uint i = 0; i < numNodes; i++) {
        for (uint j = 0; j < i; j++) {
            if (adjMatrix[i][j] != adjMatrix[j][i]) {
                cerr << "Adjacency matrix ill defined: node (" << i << "," << j;
                cerr << ") is not symmetric" << endl;
                isWellDefined = false;
            }
            if (adjMatrix[i][j]) {
                bool found = false;
                for (uint k = 0; not found and k < adjLists[i].size(); k++) {
                    if (adjLists[i][k] == j) {
                        found = true;
                    }
                }
                if (not found) {
                    cerr << "Node " << i << " adjacent to node " << j << " in adjacency";
                    cerr << " matrix but not adjacency lists" << endl;
                    isWellDefined = false;
                }
                found = false;
                for (uint k = 0; not found and k < adjLists[j].size(); k++) {
                    if (adjLists[j][k] == i) {
                        found = true;
                    }
                }
                if (not found) {
                    cerr << "Node " << i << " adjacent to node " << j << " in adjacency";
                    cerr << " matrix but not adjacency lists" << endl;
                    isWellDefined = false;
                }
            }
        }
    }
    //check that edge list does not have repeated entries
    //check no node is adjacent to itself in the edge list
    //check all indices in range
    //check that every edge in edge list appears in the adj matrix
    for (uint i = 0; i < edgeList.size(); i++) {
        if (edgeList[i][0] == edgeList[i][1]) {
            cerr << "Edge list ill defined: node " << edgeList[i][0];
            cerr << " adjacent to itself" << endl;
            isWellDefined = false;
        }
        if (edgeList[i][0] < 0 or edgeList[i][0] >= numNodes or
            edgeList[i][1] < 0 or edgeList[i][1] >= numNodes) {
            cerr << "Edge list ill defined: node out of range" << edgeList[i][0] << endl;
            isWellDefined = false;
        }
        if (not adjMatrix[edgeList[i][0]][edgeList[i][1]]) {
            cerr << "Nodes " << edgeList[i][0] << " and " << edgeList[i][1];
            cerr << " adjacent in the edge list but not in the adjacency matrix" << endl;
            isWellDefined = false;
        }
        for (uint j = 0; j < i; j++) {
            if (edgeList[i][0] == edgeList[j][0] and edgeList[i][1] == edgeList[j][1]) {
                cerr << "Edge list ill defined: edge (" << edgeList[i][0] << ",";
                cerr << edgeList[i][1] << ") is repeated" << endl;
                isWellDefined = false;
            }
            if (edgeList[i][0] == edgeList[j][1] and edgeList[i][1] == edgeList[j][0]) {
                cerr << "edge list ill defined: edge (" << edgeList[i][0] << ",";
                cerr << edgeList[i][1] << ") is repeated" << endl;
                isWellDefined = false;
            }
        }
    }
    //check that every edge in adj matrix appears in the edge list
    for (uint i = 0; i < numNodes; i++) {
        for (uint j = 0; j < i; j++) {
            if (adjMatrix[i][j]) {
                bool found = false;
                for (uint k = 0; not found and k < edgeList.size(); k++) {
                    if ((edgeList[k][0] == i and edgeList[k][1] == j) or
                        (edgeList[k][1] == i and edgeList[k][0] == j))
                        found = true;
                }
                if (not found) {
                    cerr << "Node " << i << " adjacent to node " << j << " in adjacency";
                    cerr << " matrix but not in edge list" << endl;
                    isWellDefined = false;
                }
            }
        }
    }
    return isWellDefined;
}

//returns whether this graph and other have exactly the
//same set of node names
bool Graph::sameNodeNames(const Graph& other) const {
    vector<string> thisNames = getNodeNames();
    vector<string> otherNames = other.getNodeNames();
    set<string> s1(thisNames.begin(), thisNames.end());
    set<string> s2(otherNames.begin(), otherNames.end());
    return s1 == s2;
}


void Graph::setLockedList(vector<string>& nodes, vector<string> & pairs){
    map<string,ushort> nodeMap = getNodeNameToIndexMap();
    const int size = nodeMap.size();
    vector<bool> locked (size, false);
    vector<string> lockPairs (size, "");
    for(uint i = 0; i < nodes.size(); i++){
        ushort index = nodeMap[nodes[i]];
        locked[index] = true;
        lockPairs[index] = pairs[i];
    }
    lockedList = locked;
    lockedTo = lockPairs;
    lockedCount = nodes.size();
}

vector<bool>& Graph::getLockedList(){
    return lockedList;
}
bool Graph::isLocked(uint index)
{
    return lockedList[index];
}

string Graph::getLockedTo(uint index){
    return lockedTo[index];
}

int Graph::getLockedCount(){
    return lockedCount;
}


//My functions and constructor 
Graph::Graph(int n) :
    edgeList(vector<vector<ushort> > (0)),
    adjMatrix(vector<vector<bool>> (n, vector<bool>(n))),
    adjLists(vector<vector<ushort>> (n, vector<ushort>(0))),
    lockedList(vector<bool> (0)),
    lockedTo(vector<string>(0)),         
    lockedCount(0)
    {}

void Graph::setAdjMatrix(vector<bool>& v)
{
    int k = 0;
    for (unsigned int i = 0; i < this->getNumNodes(); i++)
    {
        for (unsigned int j = 0; j < this->getNumNodes(); j++)
        {
            if (i < j)
            {
                adjMatrix[i][j] = v[k];
                if (v[k])
                {
                    this->addEdge(i,j);
                }
                k++;
            }
        }    
    }        
}

void Graph::print_adjMatrix(bool upper)
{
    for (unsigned int i = 0; i < this->getNumNodes(); i++)
    {
        for (unsigned int j = 0; j < this->getNumNodes(); j++)
        {
            if (upper)
            {
                if (i < j)
                    std::cout << adjMatrix[i][j] << " ";
                else
                    std::cout << "  ";
            }
            else
                std::cout << adjMatrix[i][j] << " ";
        }
        std::cout << "\n";
    }
}

uint Graph::getDegree(uint node) const
{
    return adjLists[node].size();
}

void Graph::set_decimal_representation(int n)
{
    decimal_representation = n;
}

int Graph::get_decimal_representation() const
{
    return decimal_representation;
}

void Graph::construct_decimal_representation()
{
    int result = 0;
    int k = 0;

    for (int i = adjMatrix.size() - 1; i >= 0; i--)
    {   
        for (int j = adjMatrix.size() - 1; j >= 0; j--)   
        {
            if (i < j)
            {
                if (adjMatrix[i][j])
                {
                    int n = pow(2, k);
                    result += n;
                }
                k++;
            }
        }
    }
    
    decimal_representation = result;
}

void Graph::print_adjLists()
{
    for (unsigned int i = 0; i < adjLists.size(); i++)
    {
        std::cout << "Start node: " << i << " --> ";
        for (unsigned int j = 0; j < adjLists[i].size(); j++)
        {
            std::cout << adjLists[i][j] << " ";
        }
        std::cout << "\n";
    }

}