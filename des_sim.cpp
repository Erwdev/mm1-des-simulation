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

// ----------------------------------------------------------------------------
// EVENT STRUCT
// Owner: ANGGOTA 1
// TODO: Implement operator< for priority_queue (min-heap)
// ----------------------------------------------------------------------------
struct Event {
    EventType type;
    double time;
    int customerID;

    // TODO ANGGOTA 1: Implement operator< (reversed for min-heap)
    // Format: bool operator<(const Event& other) const { ... }
    // Hint: return time > other.time; (reversed!)
};

// ----------------------------------------------------------------------------
// STATE STRUCT
// Owner: ANGGOTA 1 (arrival-related) + ANGGOTA 2 (service-related)
// TODO: Complete state variables
// ----------------------------------------------------------------------------
struct State {
    double clock;                  // Current simulation time
    int numInSystem;               // Total customers in system (queue + service)
    bool serverBusy;               // Server status
    double nextArrivalTime;        // Scheduled next arrival time
    std::queue<double> arrivalTimes;    // Queue of customer arrival times

    // TODO ANGGOTA 1: Initialize in constructor
    // TODO ANGGOTA 2: Add any additional state if needed
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
    int numArrived;          // Count of arrived customers
    double warmupEndTime;    // Time when warmup period ends

    // TODO ANGGOTA 2: Add warmup-related tracking
    // TODO ANGGOTA 2: Initialize in constructor
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

// ============================================================================
// SECTION 3: RANDOM NUMBER GENERATOR CLASS
// Owner: ANGGOTA 1
// ============================================================================
class RNG {
private:
    std::mt19937_64 generator;

public:
    // TODO ANGGOTA 1: Constructor with seed
    RNG(int seed) {
        // TODO: Initialize generator with seed
    }

    // TODO ANGGOTA 1: Generate exponential random variable
    // Formula: -ln(U) / rate, where U ~ Uniform(0,1)
    double exponential(double rate) {
        // TODO: Implement inverse transform method
        // uniform_real_distribution<double> dist(0.0, 1.0);
        // double u = dist(generator);
        // return -log(u) / rate;
        return 0.0; // PLACEHOLDER
    }

    // TODO ANGGOTA 1: (Optional) Test method to verify distribution
    void test(double rate, int samples = 1000) {
        // TODO: Generate samples and print mean
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
    void init() {
        // TODO ANGGOTA 1: 
        // 1. Reset state (clock=0, numInSystem=0, serverBusy=false, etc)
        // 2. Reset stats (all zeros)
        // 3. Clear FEL
        // 4. Schedule first arrival at time = exponential(lambda)

        std::cout << "[INIT] Simulation initialized" << std::endl;
    }

    // ------------------------------------------------------------------------
    // UPDATE TIME INTEGRALS
    // Owner: ANGGOTA 3
    // ✅ COMPLETED - Update area under Q(t) and B(t) curves
    // ------------------------------------------------------------------------
    void updateTimeIntegrals(double previousTime) {
        double delta = state.clock - previousTime;

        // Check if past warmup period
        if (isWarmupComplete()) {
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
    void handleArrival(Event e) {
        // TODO ANGGOTA 1:
        // 1. Save previous clock time
        // 2. Update clock to event time
        // 3. Call updateTimeIntegrals(previousTime)
        // 4. Increment numArrived
        // 5. Schedule NEXT arrival (time = clock + exponential(lambda))
        // 6. Check server status:
        //    a. If server IDLE:
        //       - Set serverBusy = true
        //       - Schedule DEPARTURE (time = clock + exponential(mu))
        //    b. If server BUSY:
        //       - Check queue capacity (if queueCap != -1)
        //       - If not full: add arrival time to queue
        //       - If full: reject (track rejections)
        // 7. numInSystem++

        std::cout << "[ARRIVAL] Customer " << e.customerID << " at t=" << e.time << std::endl;
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

        init();

        // MAIN LOOP HERE

        return computeResults();
    }

    // ------------------------------------------------------------------------
    // COMPUTE RESULTS
    // Owner: ANGGOTA 3
    // ------------------------------------------------------------------------
    RepResult computeResults() {
        RepResult result;

        // Calculate effective time period (excluding warmup)
        double period = state.clock - stats.warmupEndTime;
        if (period <= 0) period = state.clock;

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

    std::vector<RepResult> results;
    return results;
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
// ----------------------------------------------------------------------------
double calculateMean(std::vector<double>& data) {
    if (data.empty()) return 0.0;

    double sum = 0.0;
    for (double value : data) {
        sum += value;
    }
    return sum / data.size();
}

// ----------------------------------------------------------------------------
// CALCULATE STANDARD DEVIATION
// ----------------------------------------------------------------------------
double calculateStdDev(std::vector<double>& data, double mean) {
    if (data.size() <= 1) return 0.0;

    double sumSquaredDev = 0.0;
    for (double value : data) {
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
double getTValue(int n) {
    // t-distribution values for 95% CI (two-tailed, alpha=0.05)
    // Degrees of freedom = n - 1
    if (n <= 1) return 12.706;
    if (n == 2) return 4.303;
    if (n == 3) return 3.182;
    if (n == 4) return 2.776;
    if (n == 5) return 2.571;
    if (n <= 7) return 2.447;
    if (n <= 9) return 2.306;
    if (n <= 11) return 2.228;
    if (n <= 14) return 2.160;
    if (n <= 17) return 2.120;
    if (n <= 21) return 2.086;
    if (n <= 26) return 2.060;
    if (n <= 31) return 2.042;
    if (n <= 41) return 2.021;
    if (n <= 61) return 2.000;
    return 1.96; // For large n (approximates normal distribution)
}

// ----------------------------------------------------------------------------
// CALCULATE CONFIDENCE INTERVAL
// ----------------------------------------------------------------------------
Summary calculateCI(std::vector<double>& data, std::string metricName) {
    Summary summary;
    summary.metric = metricName;

    int n = data.size();
    if (n == 0) {
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
std::vector<Summary> computeSummaries(std::vector<RepResult>& results) {
    std::vector<Summary> summaries;

    if (results.empty()) return summaries;

    // Extract avgQ values
    std::vector<double> avgQ_data;
    for (auto& r : results) {
        avgQ_data.push_back(r.avgQ);
    }
    summaries.push_back(calculateCI(avgQ_data, "AvgQ"));

    // Extract utilization values
    std::vector<double> util_data;
    for (auto& r : results) {
        util_data.push_back(r.utilization);
    }
    summaries.push_back(calculateCI(util_data, "Utilization"));

    // Extract avgDelay values
    std::vector<double> delay_data;
    for (auto& r : results) {
        delay_data.push_back(r.avgDelay);
    }
    summaries.push_back(calculateCI(delay_data, "AvgDelay"));

    // Extract avgWait values
    std::vector<double> wait_data;
    for (auto& r : results) {
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



    // Header
    // file << "Metric,Mean,StdDev,CI_Lower,CI_Upper,CI_Width\n";

    // Data rows

    // std::cout << "Written: " << filename << std::endl;
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
// ----------------------------------------------------------------------------
void validateResults(Params p, std::vector<Summary>& summaries) {
    std::cout << "\n==================================================" << std::endl;
    std::cout << "  THEORETICAL VS SIMULATION COMPARISON" << std::endl;
    std::cout << "==================================================" << std::endl;

    // Calculate theoretical M/M/1 results
    double rho = p.lambda / p.mu;
    double L_theory = rho / (1.0 - rho);              // Avg number in system
    double W_theory = 1.0 / (p.mu - p.lambda);        // Avg time in system
    double Wq_theory = rho / (p.mu - p.lambda);       // Avg time in queue
    double Lq_theory = p.lambda * Wq_theory;          // Avg number in queue

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
    for (auto& s : summaries) {
        if (s.metric == "Utilization") {
            sim_util = s.mean;
            std::cout << "  Rho (Utilization) : " << s.mean
                << " [" << s.ci_lower << ", " << s.ci_upper << "]" << std::endl;
        }
        else if (s.metric == "AvgQ") {
            sim_L = s.mean;
            std::cout << "  L (Avg in system) : " << s.mean
                << " [" << s.ci_lower << ", " << s.ci_upper << "]" << std::endl;
        }
        else if (s.metric == "AvgDelay") {
            sim_W = s.mean;
            std::cout << "  W (Avg time)      : " << s.mean
                << " [" << s.ci_lower << ", " << s.ci_upper << "]" << std::endl;
        }
        else if (s.metric == "AvgWait") {
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

    if (valid) {
        std::cout << "  ✓ PASSED - All metrics within " << threshold << "% of theory" << std::endl;
    }
    else {
        std::cout << "  ✗ WARNING - Some metrics exceed " << threshold << "% deviation" << std::endl;
        std::cout << "    (May need more replications or longer warmup)" << std::endl;
    }
}

// ----------------------------------------------------------------------------
// VERIFY LITTLE'S LAW
// ----------------------------------------------------------------------------
void verifyLittlesLaw(Params p, std::vector<Summary>& summaries) {
    std::cout << "\n==================================================" << std::endl;
    std::cout << "  LITTLE'S LAW VERIFICATION" << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "\nLittle's Law: L = λ × W" << std::endl;

    // Find L and W from summaries
    double L_sim = 0.0, W_sim = 0.0;
    for (auto& s : summaries) {
        if (s.metric == "AvgQ") {
            L_sim = s.mean;
        }
        else if (s.metric == "AvgDelay") {
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
    if (deviation < 1.0) {
        std::cout << "  ✓ EXCELLENT - Little's Law holds (< 1% deviation)" << std::endl;
    }
    else if (deviation < 5.0) {
        std::cout << "  ✓ GOOD - Little's Law verified (< 5% deviation)" << std::endl;
    }
    else {
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

    // ✅ ANGGOTA 3: Compute summaries
    std::vector<Summary> summaries = computeSummaries(results);

    // ✅ ANGGOTA 3: Print results
    std::cout << "\n=== SUMMARY STATISTICS ===" << std::endl;
    for (auto& s : summaries) {
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
