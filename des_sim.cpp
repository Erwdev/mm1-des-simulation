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
const int MAX_ITER = 1000000;
// ============================================================================
// SECTION 1: ENUMS & CONSTANTS
// ============================================================================

enum EventType {
    ARRIVAL,
    DEPARTURE
};

enum TerminationMode {
    BY_SERVED,  // Terminate after maxServed customers
    BY_TIME     // Terminate after horizonT time units
};

// ============================================================================
// SECTION 2: STRUCTS
// ============================================================================


struct Event {
    EventType type;
    double time;
    int customerID;
    
    // bool: return type, operator overloading for <  comparison 
    // passing by reference const Event 
    bool operator<(const Event& other) const {
        return time > other.time;
    }

};

struct State {
    double clock;                  // Current simulation time
    int numInSystem;               // Total customers in system (queue + service)
    bool serverBusy;               // Server status
    double nextArrivalTime;        // Scheduled next arrival time
    std::queue<double> arrivalTimes;    // Queue of customer arrival times
    // secra otomatis empty di awal queue
    State() {
        clock = 0.0;
        numInSystem = 0;
        serverBusy = false;
        nextArrivalTime = 0.0;
    }

    void reset(double clock = 0.0, int numInSystem = 0, bool serverBusy = false, double nextArrivalTime = 0.0) {
        this->clock = clock;
        this->numInSystem = numInSystem;
        this->serverBusy = serverBusy;
        this->nextArrivalTime = nextArrivalTime;

        while (!arrivalTimes.empty()) {
            arrivalTimes.pop();
        }//looping to clear arrival time queue
    }

};

// ----------------------------------------------------------------------------
// STATS STRUCT
// Owner: ANGGOTA 2
// TODO: Add tracking variables for statistics
// ----------------------------------------------------------------------------
struct Stats {
    double totalDelay;       // Sum of time in system (W)
    double areaQ;            // Time-weighted queue length integral
    double areaB;            // Time-weighted server busy integral
    int numServed;           // Count of departed customers
    int numArrived;     // Count of arrived customers       
    int numRejected;         //count of rejected customer (if queue capacity is exhausted the server is full )
    double warmupEndTime;    // Time when warmup period ends
    
    // TODO ANGGOTA 2: Add warmup-related tracking
    // TODO ANGGOTA 2: Initialize in constructor
    
    // anggota 2 left to implement the initializer

    void reset(double totalDelay = 0.0, double areaQ = 0.0, double areaB = 0.0, int numServed = 0, int numArrived = 0, double warmupEndTime = 0.0) {
        this->totalDelay = totalDelay;
        this->areaQ = areaQ;
        this->areaB = areaB;
        this->numServed = numServed;
        this->numArrived = numArrived;
        this->warmupEndTime = warmupEndTime;
    }
};

// ----------------------------------------------------------------------------
// PARAMS STRUCT
// Owner: ANGGOTA 2
// TODO: Add all simulation parameters
// ----------------------------------------------------------------------------
struct Params {
    double lambda;              // Arrival rate
    double mu;                  // Service rate
    int maxServed;              // Max customers (BY_SERVED mode)
    double horizonT;            // Time horizon (BY_TIME mode)
    int warmup;                 // Warmup period (customers)
    int seed;                   // Random seed
    int queueCap;               // Queue capacity (-1 = unlimited)
    TerminationMode termMode;   // Termination mode
    std::string outdir;              // Output directory
    
    // TODO ANGGOTA 2: Add constructor with default values
    // TODO ANGGOTA 2: Add validation method (lambda < mu, etc)
};

// ----------------------------------------------------------------------------
// REPLICATION RESULT STRUCT
// Owner: ANGGOTA 2
// TODO: Store results from single replication
// ----------------------------------------------------------------------------
struct RepResult {
    int repID;
    double avgQ;          // Average queue length (L)
    double utilization;   // Server utilization (rho)
    double avgDelay;      // Average time in system (W)
    double avgWait;       // Average time in queue (Wq)
    int numServed;
    double simTime;
    
    // TODO ANGGOTA 2: Add constructor
};

// STRUCT SECTION ENDS HERE

// ============================================================================
// SECTION 3: RANDOM NUMBER GENERATOR CLASS
// Owner: ANGGOTA 1
// ============================================================================
class RNG {
private:
    std::mt19937_64 generator;
    // using mt19937_64 as instructed 
    // instancing generator with mt19937_64 type
public:
    RNG(int seed) {
        generator.seed(seed);
    }

    // Formula: -ln(U) / rate, where U ~ Uniform(0,1)
    double exponential(double rate) {
        // create uniform distribution 0, 1
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double u = dist(generator);
        if ( rate <= 0){
            throw std::invalid_argument("Rate of incoming must be >= 0");
        }
        if ( u == 0.0 ) u = std::numeric_limits<double>::min(); //avoiding log(0)
        return -log(u) / rate;
    }
};

// ============================================================================
// SECTION 4: DES CLASS
// ============================================================================
class DES {
private:
    State state;
    Stats stats;
    Params params;
    RNG rng;
    std::priority_queue<Event> FEL;  // Future Event List
    int nextCustomerID;
    
public:
    // ------------------------------------------------------------------------
    // CONSTRUCTOR
    // Owner: ALL (call from main)
    // ------------------------------------------------------------------------
    DES(Params p) : params(p), rng(p.seed) {
        nextCustomerID = 1;
    }
    
    // ------------------------------------------------------------------------
    // INIT - Initialize simulation
    // Owner: ANGGOTA 1
    // TODO: Reset state, stats, schedule first arrival
    // ------------------------------------------------------------------------
    void init(State& state, Stats& stats) {
        // Reset state
        
        state.reset();
        stats.reset();
        // clear future event list(berisi priority queue kedatangan )
        while (!FEL.empty()) {
            FEL.pop();
        }
        std::cout << "[INIT] state, stats, FEL cleared" << std::endl;
        // first arrival time generated from exponential distribution
        double firstArrivalTime = rng.exponential(params.lambda);
        scheduleEvent(ARRIVAL, firstArrivalTime);
        nextCustomerID = 1;
        std::cout << "[INIT] Simulation initialized" << std::endl;
        std::cout << "[INIT] First arrival scheduled at t=" << firstArrivalTime <<std::endl;
        
    }
    
    // ------------------------------------------------------------------------
    // UPDATE TIME INTEGRALS
    // Owner: ANGGOTA 3
    // TODO: Update area under Q(t) and B(t) curves
    // ------------------------------------------------------------------------
    void updateTimeIntegrals(double previousTime) {
        // TODO ANGGOTA 3:
        // 1. Calculate time delta = state.clock - previousTime
        // 2. If past warmup: update stats.areaQ += numInQueue * delta
        // 3. If past warmup: update stats.areaB += (serverBusy ? 1 : 0) * delta
        // Note: numInQueue = numInSystem - (serverBusy ? 1 : 0)
    }
    
    // ------------------------------------------------------------------------
    // HANDLE ARRIVAL
    // Owner: ANGGOTA 1
    // TODO: Process arrival event
    // ------------------------------------------------------------------------
    void handleArrival(Event e) {
        double prev_clock = state.clock;
        state.clock = e.time;
        updateTimeIntegrals(prev_clock);
        stats.numArrived++;

        double next_arrival = state.clock + rng.exponential(params.lambda);

        scheduleEvent(ARRIVAL, next_arrival);
        if (!state.serverBusy ){
            state.serverBusy = true;
            state.numInSystem++;
            state.arrivalTimes.push(state.clock);//record arrival time
            
            double departure_time = state.clock + rng.exponential(params.mu);
            scheduleEvent(DEPARTURE, departure_time);

            std::cout << "[ARRIVAL] Customer " << e.customerID << " at t=" << e.time << " -> SERVICE, depart at t=" << departure_time << std::endl;
        }
        else {
            if (params.queueCap == -1 || state.numInSystem < params.queueCap) {
                state.arrivalTimes.push(state.clock);
                state.numInSystem++;
                // adding arrival time into the state clock queue
                std::cout << "[ARRIVAL] Customer " << e.customerID << " at t=" << e.time << " -> QUEUE, number in system=" << state.numInSystem << std::endl;
            }else{
                // queue is full, rejecting customer 
                stats.numRejected++;
                std::cout << "[ARRIVAL] Customer " << e.customerID << " at t=" << e.time << " -> REJECTED full both queue and server" << std::endl;
            }
        }
    }
    
    // ------------------------------------------------------------------------
    // HANDLE DEPARTURE
    // Owner: ANGGOTA 2
    // TODO: Process departure event
    // ------------------------------------------------------------------------
    void handleDeparture(Event e) {
        // TODO ANGGOTA 2:
        // 1. Save previous clock time
        // 2. Update clock to event time
        // 3. Call updateTimeIntegrals(previousTime)
        // 4. Calculate delay for this customer:
        //    - Get arrival time from front of queue (if queue not empty)
        //    - delay = clock - arrivalTime
        // 5. If past warmup: update stats (totalDelay, numServed)
        // 6. numInSystem--
        // 7. Check queue:
        //    a. If queue NOT empty:
        //       - Pop from queue
        //       - Schedule DEPARTURE (time = clock + exponential(mu))
        //    b. If queue empty:
        //       - Set serverBusy = false
        
        std::cout << "[DEPARTURE] Customer " << e.customerID << " at t=" << e.time << std::endl;
    }
    
    // ------------------------------------------------------------------------
    // CHECK TERMINATION
    // Owner: ANGGOTA 2
    // TODO: Check if simulation should terminate
    // ------------------------------------------------------------------------
    bool shouldTerminate() {
        // TODO ANGGOTA 2:
        // 1. If BY_SERVED mode: return numServed >= maxServed
        // 2. If BY_TIME mode: return clock >= horizonT
        return false; // PLACEHOLDER
    }
    
    // ------------------------------------------------------------------------
    // CHECK WARMUP
    // Owner: ANGGOTA 2
    // TODO: Check if warmup period has ended
    // ------------------------------------------------------------------------
    bool isWarmupComplete() {
        // TODO ANGGOTA 2:
        // 1. If BY_SERVED: return numServed >= warmup
        // 2. If BY_TIME: return clock >= warmup
        // 3. Mark warmup end time if just completed
        return true; // PLACEHOLDER - assume always warmed up for now
    }
    
    // ------------------------------------------------------------------------
    // RUN - Main simulation loop
    // Owner: ANGGOTA 1 + ANGGOTA 2 (integration)
    // TODO: Event processing loop
    // ------------------------------------------------------------------------
    RepResult run() {
        // TODO ANGGOTA 1 & 2:
        // 1. Call init()
        // 2. While NOT shouldTerminate():
        //    a. Get next event from FEL (top + pop)
        //    b. Check if warmup complete
        //    c. Process event:
        //       - if ARRIVAL: handleArrival(event)
        //       - if DEPARTURE: handleDeparture(event)
        // 3. Calculate final statistics
        // 4. Return RepResult
        
        init(state, stats);
        std::cout << "\n===STARTING SIMULATION===" << std::endl;
        bool MAX_ITER_REACHED = false;
        int iter = 0;
        while (!shouldTerminate()){
            if (FEL.empty()) {
                std::cerr << "[ERROR] FEL is empty but simulation not terminated " << std::endl;//error handling 
                break;
            }

            Event currentEvent = FEL.top(); //query the top item on the list 
            FEL.pop();
            if (!isWarmupComplete() && isWarmupComplete()) {
                stats.warmupEndTime = state.clock;
                std::cout << "[WARMUP COMPLETE] Warmup period ended at t=" << stats.warmupEndTime << std::endl;
            }
            if(currentEvent.type == ARRIVAL){
                handleArrival(currentEvent);
            }
            else if (currentEvent.type == DEPARTURE) {
                handleDeparture(currentEvent);
            }

            if(stats.numArrived % 1000 == 0) {
                std::cout << "[PROGRESS]ITER:"<< iter << " Arrived: " << stats.numArrived 
                          << ", Served: " << stats.numServed 
                          << ", Rejected: " << stats.numRejected 
                          << ", Clock: " << state.clock << std::endl;
            }
            iter++;
            if(iter >= MAX_ITER ){
                std::cout <<'[WARNING] Max iteration reached!' << std::endl;
                MAX_ITER_REACHED = true;
                break;
            }
        }
        
        // MAIN LOOP HERE
        // will display final stats after a while of simulation 
        std::cout << "\n=== SIMULATION ENDED ===" << std::endl;
        std::cout << "Final Clock: " << state.clock << std::endl;
        std::cout << "Total Arrived: " << stats.numArrived << std::endl;
        std::cout << "Total Served: " << stats.numServed << std::endl;
        std::cout << "Total Rejected: " << stats.numRejected << std::endl;

        return computeResults();
    }
    
    // ------------------------------------------------------------------------
    // COMPUTE RESULTS
    // Owner: ANGGOTA 3
    // TODO: Calculate performance metrics
    // ------------------------------------------------------------------------
    RepResult computeResults() {
        // TODO ANGGOTA 3:
        // 1. Calculate time period (exclude warmup)
        //    period = clock - warmupEndTime
        // 2. Calculate metrics:
        //    avgQ = areaQ / period
        //    utilization = areaB / period
        //    avgDelay = totalDelay / numServed
        //    avgWait = avgDelay - (1.0 / mu)
        // 3. Create and return RepResult

    }
    
    // ------------------------------------------------------------------------
    // SCHEDULE EVENT (Helper)
    // Owner: ANGGOTA 1
    // ------------------------------------------------------------------------
    void scheduleEvent(EventType type, double time) {
        // TODO ANGGOTA 1: Create event and push to FEL
  
    }
};

// ============================================================================
// SECTION 5: MULTI-REPLICATION CONTROLLER
// Owner: ANGGOTA 2
// ============================================================================

// ----------------------------------------------------------------------------
// RUN MULTIPLE REPLICATIONS
// TODO ANGGOTA 2: Loop over replications
// ----------------------------------------------------------------------------
std::vector<RepResult> runReplications(Params params, int numReps) {
    // TODO ANGGOTA 2:
    // 1. Create vector to store results
    // 2. Loop for i = 1 to numReps:
    //    a. Create new Params with seed = baseSeed + i
    //    b. Create DES object
    //    c. Run simulation and get result
    //    d. Store result
    //    e. Print progress
    // 3. Return vector of results
    

    // return results;
}

// ============================================================================
// SECTION 6: STATISTICS & CONFIDENCE INTERVAL
// Owner: ANGGOTA 3
// ============================================================================

// ----------------------------------------------------------------------------
// SUMMARY STATISTICS STRUCT
// ----------------------------------------------------------------------------
struct Summary {
    std::string metric;
    double mean;
    double stdDev;
    double ci_lower;
    double ci_upper;
    double ci_width;
};

// ----------------------------------------------------------------------------
// CALCULATE MEAN
// TODO ANGGOTA 3: Calculate average from vector
// ----------------------------------------------------------------------------
double calculateMean(std::vector<double>& data) {
    // TODO ANGGOTA 3: Sum all values, divide by count

}

// ----------------------------------------------------------------------------
// CALCULATE STANDARD DEVIATION
// TODO ANGGOTA 3: Calculate sample std dev
// ----------------------------------------------------------------------------
double calculateStdDev(std::vector<double>& data, double mean) {
    // TODO ANGGOTA 3:
    // 1. Calculate sum of squared deviations: Σ(xi - mean)²
    // 2. Divide by (n-1) for sample variance
    // 3. Return sqrt(variance)
    

}

// ----------------------------------------------------------------------------
// CALCULATE CONFIDENCE INTERVAL
// TODO ANGGOTA 3: Calculate 95% CI using t-distribution
// ----------------------------------------------------------------------------
Summary calculateCI(std::vector<double>& data, std::string metricName) {
    // TODO ANGGOTA 3:
    // 1. Calculate mean
    // 2. Calculate std dev
    // 3. Get t-value for 95% CI (hardcode for common n, or use lookup table)
    //    Example t-values: n=10→2.262, n=15→2.145, n=20→2.093, n=30→2.045
    // 4. Calculate margin of error: t * (stdDev / sqrt(n))
    // 5. CI = [mean - margin, mean + margin]
    

}

// ----------------------------------------------------------------------------
// COMPUTE ALL SUMMARIES
// TODO ANGGOTA 3: Generate summary for all metrics
// ----------------------------------------------------------------------------
std::vector<Summary> computeSummaries(std::vector<RepResult>& results) {
    // TODO ANGGOTA 3:
    // 1. Extract each metric into separate vector
    //    - avgQ, utilization, avgDelay, avgWait
    // 2. Calculate CI for each metric
    // 3. Return vector of summaries
    

    
    // Extract avgQ

    
    
    // TODO: Repeat for other metrics (utilization, avgDelay, avgWait)
    
    // return summaries;
}

// ============================================================================
// SECTION 7: CSV OUTPUT
// Owner: ANGGOTA 2
// ============================================================================

// ----------------------------------------------------------------------------
// WRITE PER-REPLICATION CSV
// TODO ANGGOTA 2: Export detailed results
// ----------------------------------------------------------------------------
void writePerRepCSV(std::vector<RepResult>& results, std::string filename) {
    // TODO ANGGOTA 2:
    // 1. Open file for writing
    // 2. Write header: RepID,AvgQ,Utilization,AvgDelay,AvgWait,NumServed,SimTime
    // 3. Write each result as CSV row
    // 4. Close file
    
    
    // Header
   
    
    // Data rows


}

// ----------------------------------------------------------------------------
// WRITE SUMMARY CSV
// TODO ANGGOTA 2: Export confidence intervals
// ----------------------------------------------------------------------------
void writeSummaryCSV(std::vector<Summary>& summaries, std::string filename) {
    // TODO ANGGOTA 2:
    // 1. Open file for writing
    // 2. Write header: Metric,Mean,StdDev,CI_Lower,CI_Upper,CI_Width
    // 3. Write each summary as CSV row
    // 4. Close file
    
    
    
    // // Header
    // file << "Metric,Mean,StdDev,CI_Lower,CI_Upper,CI_Width\n";
    
    // Data rows
    
    std::cout << "Written: " << filename << std::endl;
}

// ============================================================================
// SECTION 8: CLI PARSER
// Owner: ANGGOTA 2
// ============================================================================

// ----------------------------------------------------------------------------
// PARSE COMMAND LINE ARGUMENTS
// TODO ANGGOTA 2: Extract parameters from argv
// ----------------------------------------------------------------------------
Params parseArguments(int argc, char* argv[]) {
    // TODO ANGGOTA 2:
    // 1. Create Params with default values
    // 2. Loop through argv:
    //    - Check for --lambda, --mu, --term, --reps, etc
    //    - Parse value after flag
    //    - Update Params
    // 3. Validate parameters (lambda < mu for stability)
    // 4. Return Params
    
    Params p;
    // Default values
    
    // Parse arguments
    
    
    // Validate
  
    return p;
}

// ----------------------------------------------------------------------------
// PRINT HELP MESSAGE
// TODO ANGGOTA 2: Display usage information
// ----------------------------------------------------------------------------
void printHelp() {
    // TODO ANGGOTA 2: Print all available parameters and examples
    std::cout << "Usage: ./des_sim [OPTIONS]\n";
    std::cout << "\nOptions:\n";
    std::cout << "  --lambda <value>    Arrival rate (default: 0.9)\n";
    std::cout << "  --mu <value>        Service rate (default: 1.0)\n";
    // TODO: Add all other options
    std::cout << "\nExample:\n";
    std::cout << "  ./des_sim --lambda 0.9 --mu 1.0 --maxServed 20000 --reps 10\n";
}

// ============================================================================
// SECTION 9: VALIDATION & ANALYSIS
// Owner: ANGGOTA 3
// ============================================================================

// ----------------------------------------------------------------------------
// VALIDATE WITH THEORETICAL RESULTS
// TODO ANGGOTA 3: Compare simulation vs M/M/1 theory
// ----------------------------------------------------------------------------
void validateResults(Params p, std::vector<Summary>& summaries) {
    // TODO ANGGOTA 3:
    // 1. Calculate theoretical M/M/1 results:
    //    rho = lambda / mu
    //    L = rho / (1 - rho)
    //    W = 1 / (mu - lambda)
    //    Wq = rho / (mu - lambda)
    //    Lq = lambda * Wq
    // 2. Compare with simulation results
    // 3. Calculate % deviation
    // 4. Print comparison table
    
    // std::cout << "\n=== THEORETICAL VS SIMULATION ===" << std::endl;
    // std::cout << std::fixed << std::setprecision(4);
    // std::cout << "Rho (utilization): " << rho << std::endl;
    // std::cout << "L (theory): " << L_theory << std::endl;
    // // TODO: Print simulation results and deviation
}

// ----------------------------------------------------------------------------
// VERIFY LITTLE'S LAW
// TODO ANGGOTA 3: Check L = lambda * W
// ----------------------------------------------------------------------------
void verifyLittlesLaw(Params p, std::vector<Summary>& summaries) {
    // TODO ANGGOTA 3:
    // 1. Get L (avgQ) and W (avgDelay) from summaries
    // 2. Calculate L_calculated = lambda * W
    // 3. Compare with L from simulation
    // 4. Print verification result
    
    std::cout << "\n=== LITTLE'S LAW VERIFICATION ===" << std::endl;
    // TODO: Implement verification
}

// ============================================================================
// SECTION 10: MAIN FUNCTION
// Owner: ALL (integration point)
// ============================================================================

int main(int argc, char* argv[]) {
    std::cout << "==================================================" << std::endl;
    std::cout << "  M/M/1 Queue Discrete Event Simulation" << std::endl;
    std::cout << "==================================================" << std::endl;
    
    // TODO ANGGOTA 2: Parse command line arguments
    Params params = parseArguments(argc, argv);
    int numReps = 10;  // TODO: Get from CLI
    
    // Print configuration
    std::cout << "\nConfiguration:" << std::endl;
    std::cout << "  Lambda: " << params.lambda << std::endl;
    std::cout << "  Mu: " << params.mu << std::endl;
    std::cout << "  Rho: " << (params.lambda / params.mu) << std::endl;
    std::cout << "  Replications: " << numReps << std::endl;
    
    // TODO ANGGOTA 2: Run replications
    std::vector<RepResult> results = runReplications(params, numReps);
    
    // TODO ANGGOTA 3: Compute summaries
    std::vector<Summary> summaries = computeSummaries(results);
    
    // TODO ANGGOTA 3: Print results
    std::cout << "\n=== SUMMARY STATISTICS ===" << std::endl;
    for (auto& s : summaries) {
        std::cout << s.metric << ": " << s.mean 
             << " [" << s.ci_lower << ", " << s.ci_upper << "]" << std::endl;
    }
    
    // TODO ANGGOTA 2: Write CSV files
    writePerRepCSV(results, params.outdir + "results_per_rep.csv");
    writeSummaryCSV(summaries, params.outdir + "summary.csv");
    
    // TODO ANGGOTA 3: Validation
    validateResults(params, summaries);
    verifyLittlesLaw(params, summaries);
    
    std::cout << "\n=== SIMULATION COMPLETE ===" << std::endl;
    
    return 0;
}

// ============================================================================
// END OF FILE
// ============================================================================

/*
COMPILATION:
g++ -std=c++17 -O2 -Wall -o des_sim des_sim.cpp

EXECUTION EXAMPLES:
./des_sim --lambda 0.9 --mu 1.0 --maxServed 20000 --warmup 1000 --reps 10 --term served
./des_sim --lambda 0.5 --mu 1.0 --horizonT 50000 --warmup 1000 --reps 20 --term time
*/