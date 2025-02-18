#ifndef SANA_HPP
#define SANA_HPP
#include "Method.hpp"
#include <map>
#include <tuple>
#include <random>
#include "../measures/localMeasures/LocalMeasure.hpp"
#include "../measures/Measure.hpp"
#include "../measures/MeasureCombination.hpp"
#include "../utils/randomSeed.hpp"
#include "../utils/ParetoFront.hpp"
#include "../measures/ExternalWeightedEdgeConservation.hpp"

#ifdef WEIGHTED
#define PARAMS int aligEdges, int g1Edges, int inducedEdges, int g2Edges, double TCSum, int localScoreSum, int n1, double wecSum, double ewecSum, int ncSum, unsigned int trueA_back, double g1WeightedEdges, double g2WeightedEdges, int squaredAligEdges
#else
#define PARAMS int aligEdges, int g1Edges, int inducedEdges, int g2Edges, double TCSum, int localScoreSum, int n1, double wecSum, double ewecSum, int ncSum, unsigned int trueA_back
#endif

class SANA: public Method {

public:
    SANA(Graph* G1, Graph* G2,
      double TInitial, double TDecay, double t, bool usingIterations, bool addHillClimbing, MeasureCombination* MC, string& objectiveScore
#ifdef WEIGHTED
	, string& startAligName
#endif
    );
    ~SANA(){}

    Alignment run();
    unordered_set<vector<uint>*>* paretoRun(const string &fileName);
    void describeParameters(ostream& stream);
    string fileNameSuffix(const Alignment& A);

    void enableRestartScheme(double minutesNewAlignments, uint iterationsPerStep,
        uint numCandidates, double minutesPerCandidate, double minutesFinalist);

    //set temperature schedule automatically
    void searchTemperaturesByLinearRegression();
    void searchTemperaturesByStatisticalTest();
    void setTDecay(double t);
    void setTDecayAutomatically();
    //to compute TDecay automatically
    //returns a value of lambda such that with this TInitial, temperature reaches
    //0 after a certain number of minutes
    double solveTDecay();

    void setAcceptableTInitial();
    void setAcceptableTFinal();
    void setAcceptableTFinalFromManualTInitial();
    double findAcceptableTInitial(double temperature);
    double findAcceptableTFinal(double temperature);
    double findAcceptableTFinalFromManualTInitial(double temperature);

    //set temperature decay dynamically
    void setDynamicTDecay();

    double simpleSearchTInitial();

    double elapsedEstimate = 0;
    int order = 0;

    //returns the number of iterations until it stagnates when not using temperture
    void hillClimbingIterations(long long int iterTarget);
    Alignment hillClimbingAlignment(Alignment startAlignment, long long int idleCountTarget);
    Alignment hillClimbingAlignment(long long int idleCountTarget);

    //returns an approximation of the the logarithm in base e of the size of the search space
    double searchSpaceSizeLog();
    string startAligName = "";
    void prune(string& startAligName);
#ifdef CORES
    Matrix<ulong> getCoreFreq();
    vector<ulong> getCoreCount();
#endif
    //to compute TDecay automatically
    //returns a value of lambda such that with this TInitial, temperature reaches
    //0 after a certain number of minutes
    double searchTDecay(double TInitial, double minutes);
    double searchTDecay(double TInitial, uint iterations);
    double getTInitial(void), getTFinal(void), getTDecay(void);

private:
    int maxTriangles = 0;

    //Temperature Boundaries. Use these after the tinitial has been determined
    double lowerTBound = 0;
    double upperTBound = 0;

    //data structures for the networks
    uint n1;
    uint n2;
    double g1Edges; //stored as double because it appears in division
    
    double omega;
    double waccAlpha;
    double beta;
    double betaStar;
    double betaHash;
    double betaHashPow;
#ifdef WEIGHTED
    double g1WeightedEdges;
    double g2WeightedEdges;
#endif
    double g2Edges; //stored as double because it appears in division
    Matrix<MATRIX_UNIT> G1Matrix;
    Matrix<MATRIX_UNIT> G2Matrix;
    vector<vector<uint> > G1AdjLists;
    vector<vector<uint> > G2AdjLists;

    void initTau(void);
    vector<uint> unLockedNodesG1;
    bool nodesHaveType = false;
    //random number generation
    mt19937 gen;
    uniform_int_distribution<> G1RandomNode;
    uniform_int_distribution<> G1RandomUnlockedNodeDist;
    uniform_int_distribution<> G2RandomUnassignedNode;
    uniform_real_distribution<> randomReal;
    uniform_int_distribution<> G1RandomUnlockedGeneDist;
    uniform_int_distribution<> G1RandomUnlockedmiRNADist;
    uniform_int_distribution<> G2RandomUnassignedGeneDist;
    uniform_int_distribution<> G2RandomUnassignedmiRNADist;


    uint G1RandomUnlockedNode();
    uint G1RandomUnlockedNode(uint source1); // used in nodes-have-type because
    uint G1RandomUnlockedNode_Fast();
    uint G2RandomUnlockedNode(uint target1);
    uint G2RandomUnlockedNode_Fast();

    //temperature schedule
    double TInitial;
    double TFinal;
    double TDecay;
    double minutes = 0;
    bool usingIterations;
    uint maxIterations = 0;
    uint iterationsPerformed = 0;
    const double TInitialScaling = 1;
    const double TDecayScaling = 1;
    //to compute TDecay dynamically
    //vector holds "ideal" temperature values at certain execution times
    bool dynamic_tdecay;
    vector<double> tau;
    double SANAtime;

    double Temperature;
    double temperatureFunction(long long int iter, double TInitial, double TDecay);
    double acceptingProbability(double energyInc, double Temperature);
    double trueAcceptingProbability();
    //to compute TInitial automatically
    //returns a value of TInitial such that the temperature is random
    double scoreForTInitial(double TInitial);
    bool isRandomTInitial(double TInitial, double highThresholdScore, double lowThresholdScore);
    double scoreRandom();

    bool initializedIterPerSecond;
    double iterPerSecond;
    double getIterPerSecond();
    void initIterPerSecond();

    vector<double> energyIncSample(double temp = 0.0);
    double expectedNumAccEInc(double temp, const vector<double>& energyIncSample);

    //data structures for the solution space search
    double changeProbability[2];
    vector<bool> *assignedNodesG2;
    vector<uint> *unassignedNodesG2;
    vector<uint> *unassignedmiRNAsG2;
    vector<uint> *unassignedgenesG2;
    vector<uint>* A;
    //initializes all the necessary data structures for a new run
    void initDataStructures(const Alignment& startA);

    bool addHillClimbing; //for post-run hill climbing

    //objective function
    MeasureCombination* MC;
    double eval(const Alignment& A);
    bool scoreComparison(double newAligEdges, double newInducedEdges, double newTCSum, double newLocalScoreSum, double newWecSum, double newNcSum, double& newCurrentScore, double newEwecSum, double newSquaredAligEdges);
    double ecWeight;
    double s3Weight;
    double waccWeight;
    double fbetaWeight;
    double fbetastarWeight;
    double fBetaHashWeight;
    double fBetaHashPowWeight;
    double mccWeight;
    double bmWeight;
    double mkWeight;
    double iliaWeight;
    double icsWeight;
    double wecWeight;
    double secWeight;
    double ncWeight;
    double localWeight;
    double mecWeight;
    double sesWeight;
    double ewecWeight;
    double TCWeight;    
    

    enum class Score{sum, product, inverse, max, min, maxFactor, pareto};
    Score score;

    //For pareto mode
    unsigned int paretoInitial;
    unsigned int paretoCapacity;
    unsigned int paretoIterations;

    //restart scheme
    bool restart;
    //parameters
    double minutesNewAlignments;
    uint iterationsPerStep;
    uint numCandidates;
    double minutesPerCandidate;
    double minutesFinalist;
    //data structures
    uint newAlignmentsCount;
    vector<Alignment> candidates;
    vector<double> candidatesScores;
    //functions
    Alignment runRestartPhases();
    uint getLowestIndex() const;
    uint getHighestIndex() const;


    //to evaluate EC incrementally
    bool needAligEdges;
    int aligEdges;
    int aligEdgesIncChangeOp(uint source, uint oldTarget, uint newTarget);
    int aligEdgesIncSwapOp(uint source1, uint source2, uint target1, uint target2);

    // to evaluate SES incrementally
    bool needSquaredAligEdges;
    int squaredAligEdges;
    int squaredAligEdgesIncChangeOp(uint source, uint oldTarget, uint newTarget);
    int squaredAligEdgesIncSwapOp(uint source1, uint source2, uint target1, uint target2);

    //to evaluate EC incrementally
    bool needSec;
    double secSum;

    //to evaluate S3 incrementally
    bool needInducedEdges;
    int inducedEdges;
    int inducedEdgesIncChangeOp(uint source, uint oldTarget, uint newTarget);

    bool needTC;
    double TCSum;
    double TCIncChangeOp(uint source, uint oldTarget, uint newTarget);
    double TCIncSwapOp(uint source1, uint source2, uint target1, uint target2);

    //to evaluate nc incrementally
    bool needNC;
    int ncSum;
    vector<uint> trueA;
    int ncIncChangeOp(uint source, uint oldTarget, uint newTarget);
    int ncIncSwapOp(uint source1, uint source2, uint target1, uint target2);

    //to evaluate wec incrementally
    bool needWec;
    double wecSum;
    vector<vector<float> > wecSims;
    double WECIncChangeOp(uint source, uint oldTarget, uint newTarget);
    double WECIncSwapOp(uint source1, uint source2, uint target1, uint target2);

    //to evaluate ewec incrementally
    bool needEwec;
    ExternalWeightedEdgeConservation* ewec;
    double ewecSum;
    double EWECIncChangeOp(uint source, uint oldTarget, uint newTarget);
    double EWECIncSwapOp(uint source1, uint source2, uint target1, uint target2);
    double EWECSimCombo(uint source, uint target);

    //to evaluate local measures incrementally
    bool needLocal;
    double localScoreSum;
    map<string, double>* localScoreSumMap = new map<string, double>;
    vector<vector<float> > sims;
#ifdef CORES
    Matrix<ulong> coreFreq;
    vector<ulong> coreCount; // number of times this node in g1 was sampled.
    Matrix<double> weightedCoreFreq; // weighted by pBad below
    vector<double> totalCoreWeight; // sum of all pBads, for each node in G1.
#endif
    map<string, vector<vector<float> > > localSimMatrixMap;
    double localScoreSumIncChangeOp(vector<vector<float> > const & sim, uint const & source, uint const & oldTarget, uint const & newTarget);
    double localScoreSumIncSwapOp(vector<vector<float> > const & sim, uint const & source1, uint const & source2, uint const & target1, uint const & target2);



    //other execution options
    bool constantTemp; //tempertare does not decrease as a function of iteration
    bool enableTrackProgress; //shows output periodically
    void trackProgress(long long int i, bool end = false);
    double avgEnergyInc;


    //algorithm
    Alignment simpleRun(const Alignment& A, double maxExecutionSeconds,
        long long int& iter);
    Alignment simpleRun(const Alignment& A, long long int maxExecutionIterations,
        long long int& iter);
    Alignment simpleRun(const Alignment& startA, double maxExecutionSeconds, long long int maxExecutionIterations,
        long long int& iter);
    unordered_set<vector<uint>*>* simpleParetoRun(const Alignment& A, double maxExecutionSeconds,
        long long int& iter, const string &fileName);
    unordered_set<vector<uint>*>* simpleParetoRun(const Alignment& A, long long int maxExecutionIterations,
        long long int& iter, const string &fileName);

    double currentScore;
    double previousScore;
    double energyInc;
    vector<double> sampledProbability;
    void SANAIteration();
    void performChange(int type);
    void performSwap(int type);


    //others
    Timer timer;
    void setInterruptSignal(); //allows the execution to be paused with Control+C

    // Used to support locking
    Alignment getStartingAlignment();
    bool implementsLocking(){ return true; }

    double pForTInitial(double TInitial);
    double getPforTInitial(const Alignment& startA, double maxExecutionSeconds,
        long long int& iter);
    string getFolder();
    string haveFolder();
    string mkdir(const std::string& file);
    tuple<int, double, int, double, double, double> regress(double start, double end, int amount);


    //Mostly for pareto front, to hold multiple alignments and scores
    unordered_map<string, int> mapScoresToIndexes();
    void prepareMeasureDataByAlignment();
    void insertCurrentAndPrepareNewMeasureDataByAlignment(vector<double> &addScores);
    vector<double> translateScoresToVector();
    void insertCurrentAlignmentAndData();
    void removeAlignmentData(vector<uint>* toRemove);
    void initializeParetoFront();
    vector<double> getMeasureScores(double newAligEdges, double newInducedEdges, double newTCSum,
                                     double newLocalScoreSum, double newWecSum, double newNcSum,
                                     double newEwecSum, double newSquaredAligEdges);
    bool dominates(vector<double> &left, vector<double> &right);
    void printParetoFront(const string &fileName);
    void deallocateParetoData();
    int numOfMeasures;
    vector<string> measureNames;
    int currentMeasure;
    vector<double> currentScores;
    ParetoFront paretoFront;
    //vector<bool>* newAN = new vector<bool>(0);
    //vector<uint>* newUAN = new vector<uint>(0);
    //vector<uint>* newUmiRNA = new vector<uint>(0);
    //vector<uint>* newUG = new vector<uint>(0);
    unordered_map<string, int> scoreNamesToIndexes;
    unordered_set<vector<uint>*>* storedAlignments = new unordered_set<vector<uint>*>;
    unordered_map<vector<uint>*, vector<bool>*> storedAssignedNodesG2;
    unordered_map<vector<uint>*, vector<uint>*> storedUnassignedNodesG2;
    unordered_map<vector<uint>*, vector<uint>*> storedUnassignedmiRNAsG2;
    unordered_map<vector<uint>*, vector<uint>*> storedUnassignedgenesG2;
    unordered_map<vector<uint>*, int> storedAligEdges;
    unordered_map<vector<uint>*, int> storedSquaredAligEdges;
    unordered_map<vector<uint>*, int> storedInducedEdges;
    unordered_map<vector<uint>*, double> storedLocalScoreSum;
    unordered_map<vector<uint>*, double> storedWecSum;
    unordered_map<vector<uint>*, double> storedEwecSum;
    unordered_map<vector<uint>*, int> storedNcSum;
    unordered_map<vector<uint>*, double> storedTCSum;
    unordered_map<vector<uint>*, double> storedCurrentScore;
    unordered_map<vector<uint>*, map<string, double>*> storedLocalScoreSumMap;
    typedef double (*calc)(PARAMS);
    unordered_map<string, calc> measureCalculation;
    unordered_set<string> localScoreNames = { "edgec", "edged", "esim", "go", "graphlet",
                                              "graphletcosine", "graphletlgraal", "importance",
                                              "nodec", "noded", "sequence" };
};

#endif
