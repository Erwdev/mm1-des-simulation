// ============================================================================
// DES_SIM.CPP - M/M/1 Queue Simulation
// Discrete Event Simulation with Multi-Replication & Confidence Interval
// ============================================================================

#include <iostream>
#include <queue>
#include <vector>
#include <random>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <limits>
const int MAX_ITER = 1000000;
// ============================================================================
// SECTION 1: ENUMS & CONSTANTS
// ============================================================================

enum EventType
{
    ARRIVAL,
    DEPARTURE
};

enum TerminationMode
{
    BY_SERVED, // Terminate after maxServed customers
    BY_TIME    // Terminate after horizonT time units
};

// ============================================================================
// SECTION 2: STRUCTS
// ============================================================================

struct Event
{
    EventType type;
    double time;
    int customerID;
    // bool: return type, operator overloading for <  comparison
    // passing by reference const Event
    bool operator<(const Event &other) const
    {
        return time > other.time;
    }
}; // struct with data structure priority queue min heap

struct State
{
    double clock;                    // Current simulation time
    int numInSystem;                 // Total customers in system (queue + service)
    bool serverBusy;                 // Server status
    double nextArrivalTime;          // Scheduled next arrival time
    std::queue<double> arrivalTimes; // Queue of customer arrival times
    // secra otomatis empty di awal queue
    State()
    {
        clock = 0.0;
        numInSystem = 0;
        serverBusy = false;
        nextArrivalTime = 0.0;
    }

    void reset()
    {
        clock = 0.0;
        numInSystem = 0;
        serverBusy = false;
        nextArrivalTime = 0.0;

        while (!arrivalTimes.empty())
        {
            arrivalTimes.pop();
        } // looping to clear arrival time queue
    }
};

// ---------------------------------------------------------------------------
// STATS STRUCT
// ----------------------------------------------------------------------------
struct Stats
{
    double totalDelay;    // Sum of time in system (W)
    double areaQ;         // Time-weighted queue length integral
    double areaB;         // Time-weighted server busy integral
    int numServed;        // Count of departed customers
    int totalServed;      // Total served customers including warmup
    int numArrived;       // Count of arrived customers
    int numRejected;      // count of rejected customer (if queue capacity is exhausted the server is full )
    double warmupEndTime; // Time when warmup period ends
    
    Stats() : totalDelay(0.0), areaQ(0.0), areaB(0.0), numServed(0), numArrived(0), warmupEndTime(0.0) {}
    
    void reset()
    {
        totalDelay = 0.0;
        areaQ = 0.0;
        areaB = 0.0;
        numServed = 0;
        numArrived = 0;
        numRejected = 0;
        warmupEndTime = 0.0;
    }

};

// ----------------------------------------------------------------------------
// PARAMS STRUCT
// ----------------------------------------------------------------------------
struct Params {
    double lambda;              // Arrival rate
    double mu;                  // Service rate
    int maxServed;              // Max customers (BY_SERVED mode)
    double horizonT;            // Time horizon (BY_TIME mode)
    int warmup;                 // Warmup period (customers)
    int reps;                    // Number of replications
    int seed;                   // Random seed
    int queueCap;               // Queue capacity (-1 = unlimited)
    TerminationMode termMode;   // Termination mode
    std::string outdir;         // Output directory

    bool helpFlag = false;      // Help flag
    
    Params() : lambda(0.9), mu(1.0), maxServed(20000), horizonT(20000.0), warmup(1000), reps(10), seed(12345), queueCap(-1), termMode(BY_SERVED), outdir("./"), helpFlag(false) {}

    bool isValid() {
        return lambda < mu;
    }
};

// ----------------------------------------------------------------------------
// REPLICATION RESULT STRUCT
// Owner: ANGGOTA 2
// ----------------------------------------------------------------------------
struct RepResult
{
    int repID;
    double avgQ;        // Average queue length (L)
    double utilization; // Server utilization (rho)
    double avgDelay;    // Average time in system (W)
    double avgWait;     // Average time in queue (Wq)
    int numServed;
    double simTime;

    RepResult() : repID(0), avgQ(0.0), utilization(0.0), avgDelay(0.0), avgWait(0.0), numServed(0), simTime(0.0) {}
};

// STRUCT SECTION ENDS HERE

// ============================================================================
// SECTION 3: RANDOM NUMBER GENERATOR CLASS
// Owner: ANGGOTA 1
// ============================================================================
class RNG
{
private:
    std::mt19937_64 generator;

    // using mt19937_64 as instructed
    // instancing generator with mt19937_64 type

public:
    RNG(int seed)
    {
        generator.seed(seed);
    }

    // Formula: -ln(U) / rate, where U ~ Uniform(0,1)
    double exponential(double rate)
    {
        // create uniform distribution 0, 1
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double u = dist(generator);
        if (rate <= 0)
        {
            throw std::invalid_argument("Rate of incoming must be >= 0");
        }
        if (u == 0.0)
            u = std::numeric_limits<double>::min(); // avoiding log(0)
        return -log(u) / rate;
    }
};

// ============================================================================
// SECTION 4: DES CLASS
// ============================================================================
class DES
{
private:
    State state;
    Stats stats;
    Params params;
    RNG rng;
    std::priority_queue<Event> FEL; // Future Event List
    int nextCustomerID;

public:
    // ------------------------------------------------------------------------
    // CONSTRUCTOR
    // Owner: ALL (call from main)
    // ------------------------------------------------------------------------
    DES(Params p) : params(p), rng(p.seed)
    {
        nextCustomerID = 1;
    }

    // ------------------------------------------------------------------------
    // INIT - Initialize simulation
    // Owner: ANGGOTA 1
    // TODO: Reset state, stats, schedule first arrival
    // ------------------------------------------------------------------------
    void init()
    {
        // Reset state

        state.reset();
        stats.reset();
        // clear future event list(berisi priority queue kedatangan )
        while (!FEL.empty())
        {
            FEL.pop();
        }
        std::cout << "[INIT] state, stats, FEL cleared" << std::endl;
        // first arrival time generated from exponential distribution
        double firstArrivalTime = rng.exponential(params.lambda);
        scheduleEvent(ARRIVAL, firstArrivalTime);
        nextCustomerID = 1;

        std::cout << "[INIT] Simulation initialized" << std::endl;
        std::cout << "[INIT] First arrival scheduled at t=" << firstArrivalTime << std::endl;
    }

    // ------------------------------------------------------------------------
    // UPDATE TIME INTEGRALS
    // Owner: ANGGOTA 3
    // ✅ COMPLETED - Update area under Q(t) and B(t) curves
    // ------------------------------------------------------------------------
    void updateTimeIntegrals(double previousTime)
    {
        double delta = state.clock - previousTime;

        // Check if past warmup period
        if (isWarmupComplete())
        {
            // Calculate number in queue (excluding server)
            int numInQueue = state.numInSystem - (state.serverBusy ? 1 : 0);

            // Update area under Q(t) curve (queue length over time)
            stats.areaQ += numInQueue * delta;

            // Update area under B(t) curve (server busy time)
            stats.areaB += (state.serverBusy ? 1.0 : 0.0) * delta;
        }
    }

    // ------------------------------------------------------------------------
    // HANDLE ARRIVAL
    // Owner: ANGGOTA 1
    // TODO: Process arrival event
    // ------------------------------------------------------------------------
    void handleArrival(Event e)
    {

        double prev_clock = state.clock;
        state.clock = e.time;
        updateTimeIntegrals(prev_clock);
        stats.numArrived++;

        double next_arrival = state.clock + rng.exponential(params.lambda);

        scheduleEvent(ARRIVAL, next_arrival);
        if (!state.serverBusy)
        {
            state.serverBusy = true;
            state.numInSystem++;
            state.arrivalTimes.push(state.clock); // record arrival time

            double departure_time = state.clock + rng.exponential(params.mu);
            scheduleEvent(DEPARTURE, departure_time);

            std::cout << "[ARRIVAL] Customer " << e.customerID << " at t=" << e.time << " -> SERVICE, depart at t=" << departure_time << std::endl;
        }
        else
        {
            if (params.queueCap == -1 || state.numInSystem < params.queueCap)
            {
                state.arrivalTimes.push(state.clock);
                state.numInSystem++;
                // adding arrival time into the state clock queue
                std::cout << "[ARRIVAL] Customer " << e.customerID << " at t=" << e.time << " -> QUEUE, number in system=" << state.numInSystem << std::endl;
            }
            else
            {
                // queue is full, rejecting customer
                stats.numRejected++;
                std::cout << "[ARRIVAL] Customer " << e.customerID << " at t=" << e.time << " -> REJECTED full both queue and server" << std::endl;
            }
        }
    }

    // ------------------------------------------------------------------------
    // HANDLE DEPARTURE
    // Owner: ANGGOTA 2
    // ------------------------------------------------------------------------
    void handleDeparture(Event e)
    {
        double previousTime = state.clock;
        state.clock = e.time;
        updateTimeIntegrals(previousTime);

        stats.totalServed++;

        if (isWarmupComplete()) {
            double arrivalTime = state.arrivalTimes.front();
            double delay = state.clock - arrivalTime;
            stats.totalDelay += delay;
            stats.numServed++;
        }

        state.arrivalTimes.pop();
        state.numInSystem--;

        // Check next in queue if exists
        if (!state.arrivalTimes.empty()) {
            double departureTime = state.clock + rng.exponential(params.mu);
            scheduleEvent(DEPARTURE, departureTime);
        } else {
            state.serverBusy = false;
        }
        
        std::cout << "[DEPARTURE] Customer " << e.customerID << " at t=" << e.time << std::endl;
    }

    // ------------------------------------------------------------------------
    // CHECK TERMINATION
    // Owner: ANGGOTA 2
    // ------------------------------------------------------------------------
    bool shouldTerminate()
    {
        if (params.termMode == BY_SERVED) {
            std::cout << stats.numServed << "|" << params.maxServed << std::endl;
            if (stats.numServed >= params.maxServed) {
                std::cout << "[TERMINATION] Reached max served: " << stats.numServed << std::endl;
                return true;
            }
        } else if (params.termMode == BY_TIME) {
            if (state.clock >= params.horizonT) {
                std::cout << "[TERMINATION] Reached time horizon: " << state.clock << std::endl;
                return true;
            }
        }
        return false;
    }

    // ------------------------------------------------------------------------
    // CHECK WARMUP
    // Owner: ANGGOTA 2
    // ------------------------------------------------------------------------
    bool isWarmupComplete()
    {
        if (params.termMode == BY_SERVED) {
            if (stats.totalServed >= params.warmup) {
                stats.warmupEndTime = state.clock;
                // std::cout << "[WARMUP COMPLETE] at t=" << state.clock << std::endl;
                return true;
            }
        } else if (params.termMode == BY_TIME) {
            if (state.clock >= params.warmup) {
                stats.warmupEndTime = params.warmup;
                // std::cout << "[WARMUP COMPLETE] at t=" << params.warmup << std::endl;
                return true;
            }
        }
        return false;
    }

    // ------------------------------------------------------------------------
    // RUN - Main simulation loop
    // Owner: ANGGOTA 1 + ANGGOTA 2 (integration)
    // TODO: Event processing loop
    // ------------------------------------------------------------------------
    RepResult run()
    {
        init();

        std::cout << "\n===STARTING SIMULATION===" << std::endl;

        int iter = 0;
        while (!shouldTerminate())
        {
            if (FEL.empty())
            {
                std::cerr << "[ERROR] FEL is empty but simulation not terminated " << std::endl; // error handling
                break;
            }

            Event currentEvent = FEL.top(); // query the top item on the list
            FEL.pop();

            bool wasWarmupComplete = isWarmupComplete();

            if (currentEvent.type == ARRIVAL)
            {
                handleArrival(currentEvent);
            }
            else if (currentEvent.type == DEPARTURE)
            {
                handleDeparture(currentEvent);
            }

            if (!wasWarmupComplete && isWarmupComplete())
            {
                stats.warmupEndTime = state.clock;
                std::cout << "[WARMUP COMPLETE] Warmup period ended at t=" << stats.warmupEndTime << std::endl;
            }

            if (stats.numArrived % 1000 == 0)
            {
                std::cout << "[PROGRESS]ITER:" << iter << " Arrived: " << stats.numArrived
                          << ", Served: " << stats.numServed
                          << ", Rejected: " << stats.numRejected
                          << ", Clock: " << state.clock << std::endl;
            }
            iter++;
            if (iter >= MAX_ITER)
            {
                std::cout << "[WARNING] Max iteration reached!" << std::endl;
                break;
            }
        }

        std::cout << "\n===SIMULATION ENDED===" << std::endl;
        std::cout << "Final Clock: " << state.clock << std::endl;
        std::cout << "Total Served: " << stats.numServed << std::endl;
        std::cout << "Total Rejected: " << stats.numRejected << std::endl;

        return computeResults();     
    }

    // ------------------------------------------------------------------------
    // COMPUTE RESULTS
    // Owner: ANGGOTA 3
    // ------------------------------------------------------------------------
    RepResult computeResults()
    {
        RepResult result;

        // Calculate effective time period (excluding warmup)
        double period = state.clock - stats.warmupEndTime;
        if (period <= 0)
            period = state.clock;

        // Calculate performance metrics
        result.avgQ = (period > 0) ? (stats.areaQ / period) : 0.0;
        result.utilization = (period > 0) ? (stats.areaB / period) : 0.0;
        result.avgDelay = (stats.numServed > 0) ? (stats.totalDelay / stats.numServed) : 0.0;
        result.avgWait = result.avgDelay - (1.0 / params.mu);
        result.numServed = stats.numServed;
        result.simTime = state.clock;
        result.repID = 0; // Will be set by caller

        return result;
    }

    // ------------------------------------------------------------------------
    // SCHEDULE EVENT (Helper)
    // Owner: ANGGOTA 1
    // ------------------------------------------------------------------------
    void scheduleEvent(EventType type, double time)
    {
        // TODO ANGGOTA 1: Create event and push to FEL
        Event newEvent;
        newEvent.type = type;
        newEvent.time = time;
        newEvent.customerID = nextCustomerID++;

        FEL.push(newEvent);

        std::cout << "[SCHEDULE] Event type= " << (type == ARRIVAL ? "ARRIVAL" : "DEPARTURE") << ", time= " << time << ", customerID= " << newEvent.customerID << std::endl;
    }
};

// ============================================================================
// SECTION 5: MULTI-REPLICATION CONTROLLER
// Owner: ANGGOTA 2
// ============================================================================

std::vector<RepResult> runReplications(Params params) {
    std::vector<RepResult> results;

    for (int i = 0; i < params.reps; i++) {
        Params repParams = params;
        repParams.seed = params.seed + i;

        DES sim(repParams);
        RepResult result = sim.run();
        result.repID = i + 1;
        results.push_back(result);

        std::cout << "[REPLICATION " << (i + 1) << "/" << params.reps << "] Completed" << std::endl;
    }
    

    return results;
}

// ============================================================================
// SECTION 6: STATISTICS & CONFIDENCE INTERVAL
// Owner: ANGGOTA 3
// ============================================================================

// ----------------------------------------------------------------------------
// SUMMARY STATISTICS STRUCT
// ----------------------------------------------------------------------------
struct Summary
{
    std::string metric;
    double mean;
    double stdDev;
    double ci_lower;
    double ci_upper;
    double ci_width;
};

// ----------------------------------------------------------------------------
// CALCULATE MEAN
// ----------------------------------------------------------------------------
double calculateMean(std::vector<double> &data)
{
    if (data.empty())
        return 0.0;

    double sum = 0.0;
    for (double value : data)
    {
        sum += value;
    }
    return sum / data.size();
}

// ----------------------------------------------------------------------------
// CALCULATE STANDARD DEVIATION
// ----------------------------------------------------------------------------
double calculateStdDev(std::vector<double> &data, double mean)
{
    if (data.size() <= 1)
        return 0.0;

    double sumSquaredDev = 0.0;
    for (double value : data)
    {
        double dev = value - mean;
        sumSquaredDev += dev * dev;
    }

    // Sample variance: divide by (n-1)
    double variance = sumSquaredDev / (data.size() - 1);
    return std::sqrt(variance);
}

// ----------------------------------------------------------------------------
// GET T-VALUE for 95% Confidence Interval
// ----------------------------------------------------------------------------
double getTValue(int n)
{
    // t-distribution values for 95% CI (two-tailed, alpha=0.05)
    // Degrees of freedom = n - 1
    if (n <= 1)
        return 12.706;
    if (n == 2)
        return 4.303;
    if (n == 3)
        return 3.182;
    if (n == 4)
        return 2.776;
    if (n == 5)
        return 2.571;
    if (n <= 7)
        return 2.447;
    if (n <= 9)
        return 2.306;
    if (n <= 11)
        return 2.228;
    if (n <= 14)
        return 2.160;
    if (n <= 17)
        return 2.120;
    if (n <= 21)
        return 2.086;
    if (n <= 26)
        return 2.060;
    if (n <= 31)
        return 2.042;
    if (n <= 41)
        return 2.021;
    if (n <= 61)
        return 2.000;
    return 1.96; // For large n (approximates normal distribution)
}

// ----------------------------------------------------------------------------
// CALCULATE CONFIDENCE INTERVAL
// ----------------------------------------------------------------------------
Summary calculateCI(std::vector<double> &data, std::string metricName)
{
    Summary summary;
    summary.metric = metricName;

    int n = data.size();
    if (n == 0)
    {
        summary.mean = 0.0;
        summary.stdDev = 0.0;
        summary.ci_lower = 0.0;
        summary.ci_upper = 0.0;
        summary.ci_width = 0.0;
        return summary;
    }

    // Calculate mean and standard deviation
    summary.mean = calculateMean(data);
    summary.stdDev = calculateStdDev(data, summary.mean);

    // Get t-value for 95% confidence interval
    double tValue = getTValue(n);

    // Calculate margin of error
    double marginOfError = tValue * (summary.stdDev / std::sqrt(n));

    // Calculate confidence interval
    summary.ci_lower = summary.mean - marginOfError;
    summary.ci_upper = summary.mean + marginOfError;
    summary.ci_width = 2.0 * marginOfError;

    return summary;
}

// ----------------------------------------------------------------------------
// COMPUTE ALL SUMMARIES
// ----------------------------------------------------------------------------
std::vector<Summary> computeSummaries(std::vector<RepResult> &results)
{
    std::vector<Summary> summaries;

    if (results.empty())
        return summaries;

    // Extract avgQ values
    std::vector<double> avgQ_data;
    for (auto &r : results)
    {
        avgQ_data.push_back(r.avgQ);
    }
    summaries.push_back(calculateCI(avgQ_data, "AvgQ"));

    // Extract utilization values
    std::vector<double> util_data;
    for (auto &r : results)
    {
        util_data.push_back(r.utilization);
    }
    summaries.push_back(calculateCI(util_data, "Utilization"));

    // Extract avgDelay values
    std::vector<double> delay_data;
    for (auto &r : results)
    {
        delay_data.push_back(r.avgDelay);
    }
    summaries.push_back(calculateCI(delay_data, "AvgDelay"));

    // Extract avgWait values
    std::vector<double> wait_data;
    for (auto &r : results)
    {
        wait_data.push_back(r.avgWait);
    }
    summaries.push_back(calculateCI(wait_data, "AvgWait"));

    return summaries;
}

// ============================================================================
// SECTION 7: CSV OUTPUT
// Owner: ANGGOTA 2
// ============================================================================

// ----------------------------------------------------------------------------
// WRITE PER-REPLICATION CSV
// TODO ANGGOTA 2: Export detailed results
// ----------------------------------------------------------------------------
void writePerRepCSV(std::vector<RepResult> &results, std::string filename)
{
    std::ofstream outFile(filename);
    if (!outFile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    // Write CSV header
    outFile << "RepID,AvgQ,Utilization,AvgDelay,AvgWait,NumServed,SimTime\n";

    // Write each RepResult
    for (const auto &res : results)
    {
        outFile << res.repID << ","
                << std::fixed << std::setprecision(4) << res.avgQ << ","
                << std::fixed << std::setprecision(4) << res.utilization << ","
                << std::fixed << std::setprecision(4) << res.avgDelay << ","
                << std::fixed << std::setprecision(4) << res.avgWait << ","
                << res.numServed << ","
                << std::fixed << std::setprecision(4) << res.simTime
                << "\n";
    }

    outFile.close();
    std::cout << "CSV file written to: " << filename << "\n";
}

// ----------------------------------------------------------------------------
// WRITE SUMMARY CSV
// TODO ANGGOTA 2: Export confidence intervals
// ----------------------------------------------------------------------------
void writeSummaryCSV(std::vector<Summary> &summaries, std::string filename)
{
    std::ofstream outFile(filename);
    if (!outFile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    // Write CSV header
    outFile << "Metric,Mean,StdDev,CI_Lower,CI_Upper,CI_Width\n";

    // Write each summary
    for (const auto &s : summaries)
    {
        outFile << s.metric << ","
                << std::fixed << std::setprecision(4) << s.mean << ","
                << std::fixed << std::setprecision(4) << s.stdDev << ","
                << std::fixed << std::setprecision(4) << s.ci_lower << ","
                << std::fixed << std::setprecision(4) << s.ci_upper << ","
                << std::fixed << std::setprecision(4) << s.ci_width
                << "\n";
    }

    outFile.close();
    std::cout << "Summary CSV file written to: " << filename << "\n";
}

// ============================================================================
// SECTION 8: CLI PARSER
// Owner: ANGGOTA 2
// ============================================================================

// ----------------------------------------------------------------------------
// PARSE COMMAND LINE ARGUMENTS
// TODO ANGGOTA 2: Extract parameters from argv
// ----------------------------------------------------------------------------
Params parseArguments(int argc, char *argv[])
{    
    Params params;

    for (int i = 1; i < argc; ++i){
        std::string arg = argv[i];

        // return early if --help argument is found
        if (arg == "--help"){
           params.helpFlag = true;
            return params;
        }

        // parse other parameters
        if (arg == "--lambda" && i + 1 < argc){
            params.lambda = std::stod(argv[i + 1]);
            ++i;
        } else if (arg == "--mu" && i + 1 < argc){
            params.mu = std::stod(argv[i + 1]);
            ++i;
        } else if (arg == "--maxServed" && i + 1 < argc){
            params.maxServed = std::stoi(argv[i + 1]);
            ++i;
        } else if (arg == "--horizonT" && i + 1 < argc){
            params.horizonT = std::stod(argv[i + 1]);
            ++i;
        } else if (arg == "--warmup" && i + 1 < argc){
            params.warmup = std::stoi(argv[i + 1]);
            ++i;
        } else if (arg == "--reps" && i + 1 < argc){
            params.reps = std::stoi(argv[i + 1]);
            ++i;
        } else if (arg == "--seed" && i + 1 < argc){
            params.seed = std::stoi(argv[i + 1]);
            ++i;
        } else if (arg == "--queueCap" && i + 1 < argc){
            params.queueCap = std::stoi(argv[i + 1]);
            ++i;
        } else if (arg == "--term" && i + 1 < argc){
            if (argv[i + 1] == std::string("served")){
                params.termMode = BY_SERVED;
            } else if (argv[i + 1] == std::string("time")){
                params.termMode = BY_TIME;
            }
            ++i;
        } else if (arg == "--outdir" && i + 1 < argc){
            params.outdir = argv[i + 1];
            ++i;
        }      
    }
    
    // Validate  
    if (!params.isValid()){
        std::cerr << "Error: lambda must be less than mu for stability." << std::endl;
        exit(1);
    }  

    return params;
}

// ----------------------------------------------------------------------------
// PRINT HELP MESSAGE
// TODO ANGGOTA 2: Display usage information
// ----------------------------------------------------------------------------
void printHelp()
{
    // TODO ANGGOTA 2: Print all available parameters and examples
    std::cout << "Usage: ./des_sim [OPTIONS]\n";
    std::cout << "\nOptions:\n";
    std::cout << "  --lambda <value>    Arrival rate (default: 0.9)\n";
    std::cout << "  --mu <value>        Service rate (default: 1.0)\n";
    std::cout << "  --maxServed <value> Max customers to serve (default: 20000)\n";
    std::cout << "  --horizonT <value>  Time horizon (default: 20000.0)\n";
    std::cout << "  --warmup <value>    Warmup period (default: 1000)\n";
    std::cout << "  --reps <value>      Number of replications (default: 10)\n";
    std::cout << "  --seed <value>      Random seed (default: 12345)\n";
    std::cout << "  --queueCap <value>  Queue capacity (-1 = unlimited, default: -1)\n";
    std::cout << "  --outdir <path>     Output directory (default: ./output/)\n";
    std::cout << "  --term <mode>      Termination mode: 'served' or 'time' (default: served)\n";
    std::cout << "\nExample:\n";
    std::cout << "   ./des_sim --lambda 0.9 --mu 1.0 --maxServed 20000 --warmup 1000 --reps 10 --seed 12345 --queueCap -1 --term served --outdir ./\n";
}

// ============================================================================
// SECTION 9: VALIDATION & ANALYSIS
// Owner: ANGGOTA 3
// ============================================================================

// ----------------------------------------------------------------------------
// VALIDATE WITH THEORETICAL RESULTS
// ----------------------------------------------------------------------------
void validateResults(Params p, std::vector<Summary> &summaries)
{

    std::cout << "\n==================================================" << std::endl;
    std::cout << "  THEORETICAL VS SIMULATION COMPARISON" << std::endl;
    std::cout << "==================================================" << std::endl;

    // Calculate theoretical M/M/1 results
    double rho = p.lambda / p.mu;
    double L_theory = rho / (1.0 - rho);        // Avg number in system
    double W_theory = 1.0 / (p.mu - p.lambda);  // Avg time in system
    double Wq_theory = rho / (p.mu - p.lambda); // Avg time in queue
    double Lq_theory = p.lambda * Wq_theory;    // Avg number in queue

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nTheoretical M/M/1 Results:" << std::endl;
    std::cout << "  Rho (Utilization) : " << rho << std::endl;
    std::cout << "  L (Avg in system) : " << L_theory << std::endl;
    std::cout << "  W (Avg time)      : " << W_theory << std::endl;
    std::cout << "  Lq (Avg in queue) : " << Lq_theory << std::endl;
    std::cout << "  Wq (Avg wait)     : " << Wq_theory << std::endl;

    std::cout << "\nSimulation Results (with 95% CI):" << std::endl;

    // Find metrics in summaries
    double sim_util = 0.0, sim_L = 0.0, sim_W = 0.0, sim_Wq = 0.0;
    for (auto &s : summaries)
    {
        if (s.metric == "Utilization")
        {
            sim_util = s.mean;
            std::cout << "  Rho (Utilization) : " << s.mean
                      << " [" << s.ci_lower << ", " << s.ci_upper << "]" << std::endl;
        }
        else if (s.metric == "AvgQ")
        {
            sim_L = s.mean;
            std::cout << "  L (Avg in system) : " << s.mean
                      << " [" << s.ci_lower << ", " << s.ci_upper << "]" << std::endl;
        }
        else if (s.metric == "AvgDelay")
        {
            sim_W = s.mean;
            std::cout << "  W (Avg time)      : " << s.mean
                      << " [" << s.ci_lower << ", " << s.ci_upper << "]" << std::endl;
        }
        else if (s.metric == "AvgWait")
        {
            sim_Wq = s.mean;
            std::cout << "  Wq (Avg wait)     : " << s.mean
                      << " [" << s.ci_lower << ", " << s.ci_upper << "]" << std::endl;
        }
    }

    // Calculate deviations
    std::cout << "\nDeviation from Theory:" << std::endl;
    std::cout << std::fixed << std::setprecision(2);

    double dev_util = std::abs(sim_util - rho) / rho * 100.0;
    double dev_L = std::abs(sim_L - L_theory) / L_theory * 100.0;
    double dev_W = std::abs(sim_W - W_theory) / W_theory * 100.0;
    double dev_Wq = std::abs(sim_Wq - Wq_theory) / Wq_theory * 100.0;

    std::cout << "  Rho deviation     : " << dev_util << "%" << std::endl;
    std::cout << "  L deviation       : " << dev_L << "%" << std::endl;
    std::cout << "  W deviation       : " << dev_W << "%" << std::endl;
    std::cout << "  Wq deviation      : " << dev_Wq << "%" << std::endl;

    // Check if simulation is within acceptable range (e.g., < 5%)
    std::cout << "\nValidation Status:" << std::endl;
    double threshold = 5.0; // 5% deviation threshold
    bool valid = (dev_util < threshold) && (dev_L < threshold) &&
                 (dev_W < threshold) && (dev_Wq < threshold);

    if (valid)
    {
        std::cout << "  ✓ PASSED - All metrics within " << threshold << "% of theory" << std::endl;
    }
    else
    {
        std::cout << "  ✗ WARNING - Some metrics exceed " << threshold << "% deviation" << std::endl;
        std::cout << "    (May need more replications or longer warmup)" << std::endl;
    }
}

// ----------------------------------------------------------------------------
// VERIFY LITTLE'S LAW
// ----------------------------------------------------------------------------
void verifyLittlesLaw(Params p, std::vector<Summary> &summaries)
{
    std::cout << "\n==================================================" << std::endl;
    std::cout << "  LITTLE'S LAW VERIFICATION" << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "\nLittle's Law: L = λ × W" << std::endl;

    // Find L and W from summaries
    double L_sim = 0.0, W_sim = 0.0;
    for (auto &s : summaries)
    {
        if (s.metric == "AvgQ")
        {
            L_sim = s.mean;
        }
        else if (s.metric == "AvgDelay")
        {
            W_sim = s.mean;
        }
    }

    // Calculate L using Little's Law
    double L_littles = p.lambda * W_sim;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nFrom Simulation:" << std::endl;
    std::cout << "  L (measured)      : " << L_sim << std::endl;
    std::cout << "  W (measured)      : " << W_sim << std::endl;
    std::cout << "  λ (arrival rate)  : " << p.lambda << std::endl;

    std::cout << "\nLittle's Law Calculation:" << std::endl;
    std::cout << "  L = λ × W         : " << L_littles << std::endl;

    // Calculate deviation
    double deviation = std::abs(L_sim - L_littles) / L_sim * 100.0;
    std::cout << "\nDeviation         : " << std::fixed << std::setprecision(2)
              << deviation << "%" << std::endl;

    // Verification status
    std::cout << "\nVerification Status:" << std::endl;
    if (deviation < 1.0)
    {
        std::cout << "  ✓ EXCELLENT - Little's Law holds (< 1% deviation)" << std::endl;
    }
    else if (deviation < 5.0)
    {
        std::cout << "  ✓ GOOD - Little's Law verified (< 5% deviation)" << std::endl;
    }
    else
    {
        std::cout << "  ✗ WARNING - Significant deviation from Little's Law" << std::endl;
    }

    std::cout << "\nNote: Small deviations are expected due to:" << std::endl;
    std::cout << "  - Finite simulation time" << std::endl;
    std::cout << "  - Warmup period effects" << std::endl;
    std::cout << "  - Random sampling variability" << std::endl;
}

// ============================================================================
// SECTION 10: MAIN FUNCTION
// Owner: ALL (integration point)
// ============================================================================

int main(int argc, char *argv[])
{
    std::cout << "==================================================" << std::endl;
    std::cout << "  M/M/1 Queue Discrete Event Simulation" << std::endl;
    std::cout << "==================================================" << std::endl;

    Params params = parseArguments(argc, argv);
    if (params.helpFlag)
    {
        printHelp();
        return 0;
    }
    
    std::cout << "\nConfiguration:" << std::endl;
    std::cout << "  Lambda: " << params.lambda << std::endl;
    std::cout << "  Mu: " << params.mu << std::endl;
    std::cout << "  Rho: " << (params.lambda / params.mu) << std::endl;
    std::cout << "  Replications: " << params.reps << std::endl;
    
    // TODO ANGGOTA 2: Run replications
    std::vector<RepResult> results = runReplications(params);
    
    // TODO ANGGOTA 3: Compute summaries
    std::vector<Summary> summaries = computeSummaries(results);

    // ✅ ANGGOTA 3: Print results
    std::cout << "\n=== SUMMARY STATISTICS ===" << std::endl;
    for (auto &s : summaries)
    {
        std::cout << s.metric << ": " << s.mean
                  << " [" << s.ci_lower << ", " << s.ci_upper << "]" << std::endl;
    }

    // TODO ANGGOTA 2: Write CSV files
    writePerRepCSV(results, params.outdir + "results_per_rep.csv");
    writeSummaryCSV(summaries, params.outdir + "summary.csv");

    // ✅ ANGGOTA 3: Validation
    validateResults(params, summaries);
    verifyLittlesLaw(params, summaries);

    std::cout << "\n=== SIMULATION COMPLETE ===" << std::endl;

    return 0;
}

// ============================================================================
// END OF FILE
// ============================================================================
