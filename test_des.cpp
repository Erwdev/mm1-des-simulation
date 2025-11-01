// ============================================================================
// TEST_DES.CPP - Unit Tests for M/M/1 DES Simulation (ANGGOTA 1)
// Testing: RNG, Event, State, Arrival Handling
// ============================================================================

#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>
#include <queue>
#include <random>
#include <chrono>

// ============================================================================
// MINIMAL DEFINITIONS (Copy from des_sim.cpp)
// ============================================================================

enum EventType {
    ARRIVAL,
    DEPARTURE
};

struct Event {
    EventType type;
    double time;
    int customerID;
    
    // ANGGOTA 1: operator< for priority_queue
    bool operator<(const Event& other) const {
        return time > other.time;  // Min-heap
    }
};

struct State {
    double clock;
    int numInSystem;
    bool serverBusy;
    double nextArrivalTime;
    std::queue<double> arrivalTimes;
    
    // Constructor
    State() : clock(0.0), numInSystem(0), serverBusy(false), nextArrivalTime(0.0) {}
    
    void reset() {
        clock = 0.0;
        numInSystem = 0;
        serverBusy = false;
        nextArrivalTime = 0.0;
        while (!arrivalTimes.empty()) arrivalTimes.pop();
    }
};

class RNG {
private:
    std::mt19937_64 generator;
    
public:
    RNG(int seed) {
        generator.seed(seed);
    }
    
    double exponential(double rate) {
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double u = dist(generator);
        return -log(u) / rate;
    }
};

// ============================================================================
// TEST FRAMEWORK
// ============================================================================

int passedTests = 0;
int failedTests = 0;
int totalTests = 0;

#define TEST(name) void name()

#define RUN_TEST(test) do { \
    totalTests++; \
    std::cout << "[" << totalTests << "] " << #test << "... "; \
    try { \
        test(); \
        std::cout << "✓ PASSED" << std::endl; \
        passedTests++; \
    } catch (const std::exception& e) { \
        std::cout << "✗ FAILED" << std::endl; \
        std::cout << "    Error: " << e.what() << std::endl; \
        failedTests++; \
    } \
} while(0)

#define ASSERT_TRUE(condition) if (!(condition)) { \
    throw std::runtime_error("Assertion failed: " #condition); \
}

#define ASSERT_FALSE(condition) if (condition) { \
    throw std::runtime_error("Assertion failed: NOT(" #condition ")"); \
}

#define ASSERT_EQ(actual, expected) if ((actual) != (expected)) { \
    throw std::runtime_error("Expected " + std::to_string(expected) + \
                           " but got " + std::to_string(actual)); \
}

#define ASSERT_NEAR(actual, expected, tolerance) if (std::abs((actual) - (expected)) > (tolerance)) { \
    throw std::runtime_error("Expected " + std::to_string(expected) + \
                           " ± " + std::to_string(tolerance) + \
                           " but got " + std::to_string(actual)); \
}

#define ASSERT_GT(actual, value) if ((actual) <= (value)) { \
    throw std::runtime_error(std::to_string(actual) + " is not > " + std::to_string(value)); \
}

#define ASSERT_LT(actual, value) if ((actual) >= (value)) { \
    throw std::runtime_error(std::to_string(actual) + " is not < " + std::to_string(value)); \
}

// ============================================================================
// SECTION 1: EVENT STRUCT TESTS
// ============================================================================

TEST(test_event_min_heap_property) {
    std::priority_queue<Event> pq;
    
    Event e1 = {ARRIVAL, 3.0, 1};
    Event e2 = {ARRIVAL, 1.0, 2};
    Event e3 = {ARRIVAL, 2.0, 3};
    
    pq.push(e1);
    pq.push(e2);
    pq.push(e3);
    
    // Should pop in order: 1.0, 2.0, 3.0
    Event first = pq.top(); pq.pop();
    ASSERT_EQ(first.time, 1.0);
    
    Event second = pq.top(); pq.pop();
    ASSERT_EQ(second.time, 2.0);
    
    Event third = pq.top(); pq.pop();
    ASSERT_EQ(third.time, 3.0);
}

TEST(test_event_comparison_order) {
    Event early = {ARRIVAL, 1.0, 1};
    Event late = {ARRIVAL, 2.0, 2};
    
    // For min-heap: late < early should be TRUE
    // Because operator< is reversed
    ASSERT_TRUE(late < early);
}

TEST(test_event_same_time) {
    Event e1 = {ARRIVAL, 1.0, 1};
    Event e2 = {DEPARTURE, 1.0, 2};
    
    // Should not crash, order doesn't matter for same time
    std::priority_queue<Event> pq;
    pq.push(e1);
    pq.push(e2);
    
    Event top = pq.top();
    ASSERT_EQ(top.time, 1.0);
}

// ============================================================================
// SECTION 2: STATE STRUCT TESTS
// ============================================================================

TEST(test_state_default_constructor) {
    State s;
    
    ASSERT_EQ(s.clock, 0.0);
    ASSERT_EQ(s.numInSystem, 0);
    ASSERT_FALSE(s.serverBusy);
    ASSERT_EQ(s.nextArrivalTime, 0.0);
    ASSERT_TRUE(s.arrivalTimes.empty());
}

TEST(test_state_reset_after_modification) {
    State s;
    
    // Modify state
    s.clock = 100.0;
    s.numInSystem = 5;
    s.serverBusy = true;
    s.nextArrivalTime = 105.0;
    s.arrivalTimes.push(10.0);
    s.arrivalTimes.push(20.0);
    s.arrivalTimes.push(30.0);
    
    // Reset
    s.reset();
    
    // Check all back to default
    ASSERT_EQ(s.clock, 0.0);
    ASSERT_EQ(s.numInSystem, 0);
    ASSERT_FALSE(s.serverBusy);
    ASSERT_EQ(s.nextArrivalTime, 0.0);
    ASSERT_TRUE(s.arrivalTimes.empty());
}

TEST(test_state_arrival_queue_operations) {
    State s;
    
    // Push arrivals
    s.arrivalTimes.push(1.0);
    s.arrivalTimes.push(2.0);
    s.arrivalTimes.push(3.0);
    
    ASSERT_EQ(s.arrivalTimes.size(), 3);
    ASSERT_EQ(s.arrivalTimes.front(), 1.0);
    
    // Pop
    s.arrivalTimes.pop();
    ASSERT_EQ(s.arrivalTimes.size(), 2);
    ASSERT_EQ(s.arrivalTimes.front(), 2.0);
}

// ============================================================================
// SECTION 3: RNG CLASS TESTS
// ============================================================================

TEST(test_rng_constructor) {
    // Should not crash
    RNG rng(12345);
}

TEST(test_rng_deterministic_with_same_seed) {
    RNG rng1(12345);
    RNG rng2(12345);
    
    double val1 = rng1.exponential(1.0);
    double val2 = rng2.exponential(1.0);
    
    // Should produce EXACT same value
    ASSERT_EQ(val1, val2);
}

TEST(test_rng_different_seeds_different_values) {
    RNG rng1(12345);
    RNG rng2(54321);
    
    double val1 = rng1.exponential(1.0);
    double val2 = rng2.exponential(1.0);
    
    // Extremely unlikely to be equal
    ASSERT_TRUE(val1 != val2);
}

TEST(test_rng_exponential_always_positive) {
    RNG rng(12345);
    
    for (int i = 0; i < 100; i++) {
        double val = rng.exponential(1.0);
        ASSERT_GT(val, 0.0);
    }
}

TEST(test_rng_exponential_mean_rate_0_5) {
    RNG rng(12345);
    double rate = 0.5;
    int samples = 10000;
    
    double sum = 0.0;
    for (int i = 0; i < samples; i++) {
        sum += rng.exponential(rate);
    }
    
    double empiricalMean = sum / samples;
    double theoreticalMean = 1.0 / rate;  // = 2.0
    
    // Within 10% tolerance
    double tolerance = 0.1 * theoreticalMean;
    ASSERT_NEAR(empiricalMean, theoreticalMean, tolerance);
}

TEST(test_rng_exponential_mean_rate_1_0) {
    RNG rng(12345);
    double rate = 1.0;
    int samples = 10000;
    
    double sum = 0.0;
    for (int i = 0; i < samples; i++) {
        sum += rng.exponential(rate);
    }
    
    double empiricalMean = sum / samples;
    double theoreticalMean = 1.0;
    
    double tolerance = 0.1;
    ASSERT_NEAR(empiricalMean, theoreticalMean, tolerance);
}

TEST(test_rng_exponential_mean_rate_2_0) {
    RNG rng(12345);
    double rate = 2.0;
    int samples = 10000;
    
    double sum = 0.0;
    for (int i = 0; i < samples; i++) {
        sum += rng.exponential(rate);
    }
    
    double empiricalMean = sum / samples;
    double theoreticalMean = 0.5;
    
    double tolerance = 0.05;
    ASSERT_NEAR(empiricalMean, theoreticalMean, tolerance);
}

TEST(test_rng_exponential_higher_rate_smaller_mean) {
    RNG rng1(12345);
    RNG rng2(12345);
    
    double rate1 = 0.5;
    double rate2 = 2.0;
    
    double sum1 = 0.0;
    double sum2 = 0.0;
    
    int samples = 1000;
    for (int i = 0; i < samples; i++) {
        sum1 += rng1.exponential(rate1);
        sum2 += rng2.exponential(rate2);
    }
    
    double mean1 = sum1 / samples;
    double mean2 = sum2 / samples;
    
    // mean1 (rate=0.5) should be > mean2 (rate=2.0)
    ASSERT_GT(mean1, mean2);
}

TEST(test_rng_exponential_variance) {
    RNG rng(12345);
    double rate = 1.0;
    int samples = 10000;
    
    double sum = 0.0;
    double sumSq = 0.0;
    
    for (int i = 0; i < samples; i++) {
        double val = rng.exponential(rate);
        sum += val;
        sumSq += val * val;
    }
    
    double mean = sum / samples;
    double variance = (sumSq / samples) - (mean * mean);
    double stdDev = sqrt(variance);
    
    // For exponential: StdDev = 1/rate = 1.0
    double theoreticalStdDev = 1.0 / rate;
    double tolerance = 0.1;
    
    ASSERT_NEAR(stdDev, theoreticalStdDev, tolerance);
}

TEST(test_rng_exponential_sequence_independence) {
    RNG rng(12345);
    
    double val1 = rng.exponential(1.0);
    double val2 = rng.exponential(1.0);
    double val3 = rng.exponential(1.0);
    
    // Values should be different (independent)
    ASSERT_TRUE(val1 != val2);
    ASSERT_TRUE(val2 != val3);
    ASSERT_TRUE(val1 != val3);
}

// ============================================================================
// SECTION 4: PRIORITY QUEUE TESTS (FEL Simulation)
// ============================================================================

TEST(test_fel_basic_operations) {
    std::priority_queue<Event> FEL;
    
    ASSERT_TRUE(FEL.empty());
    
    Event e = {ARRIVAL, 1.0, 1};
    FEL.push(e);
    
    ASSERT_FALSE(FEL.empty());
    ASSERT_EQ(FEL.size(), 1);
}

TEST(test_fel_multiple_events_sorted) {
    std::priority_queue<Event> FEL;
    
    FEL.push({ARRIVAL, 5.0, 1});
    FEL.push({ARRIVAL, 2.0, 2});
    FEL.push({DEPARTURE, 3.0, 3});
    FEL.push({ARRIVAL, 1.0, 4});
    FEL.push({DEPARTURE, 4.0, 5});
    
    // Should pop in time order: 1.0, 2.0, 3.0, 4.0, 5.0
    std::vector<double> times;
    while (!FEL.empty()) {
        times.push_back(FEL.top().time);
        FEL.pop();
    }
    
    ASSERT_EQ(times.size(), 5);
    ASSERT_EQ(times[0], 1.0);
    ASSERT_EQ(times[1], 2.0);
    ASSERT_EQ(times[2], 3.0);
    ASSERT_EQ(times[3], 4.0);
    ASSERT_EQ(times[4], 5.0);
}

TEST(test_fel_clear_operation) {
    std::priority_queue<Event> FEL;
    
    FEL.push({ARRIVAL, 1.0, 1});
    FEL.push({ARRIVAL, 2.0, 2});
    FEL.push({ARRIVAL, 3.0, 3});
    
    ASSERT_EQ(FEL.size(), 3);
    
    // Clear FEL
    while (!FEL.empty()) {
        FEL.pop();
    }
    
    ASSERT_TRUE(FEL.empty());
}

// ============================================================================
// SECTION 5: EDGE CASES
// ============================================================================

TEST(test_rng_very_small_rate) {
    RNG rng(12345);
    double rate = 0.01;  // Very small
    
    double val = rng.exponential(rate);
    
    // Should produce large values (mean = 1/0.01 = 100)
    ASSERT_GT(val, 0.0);
}

TEST(test_rng_very_large_rate) {
    RNG rng(12345);
    double rate = 100.0;  // Very large
    
    double val = rng.exponential(rate);
    
    // Should produce small values (mean = 1/100 = 0.01)
    ASSERT_GT(val, 0.0);
    ASSERT_LT(val, 10.0);  // Sanity check
}

TEST(test_state_many_arrivals_in_queue) {
    State s;
    
    // Add 1000 arrivals
    for (int i = 1; i <= 1000; i++) {
        s.arrivalTimes.push(i * 0.1);
    }
    
    ASSERT_EQ(s.arrivalTimes.size(), 1000);
    
    // Pop 500
    for (int i = 0; i < 500; i++) {
        s.arrivalTimes.pop();
    }
    
    ASSERT_EQ(s.arrivalTimes.size(), 500);
}

// ============================================================================
// SECTION 6: PERFORMANCE TESTS
// ============================================================================

TEST(test_rng_performance_100k_samples) {
    RNG rng(12345);
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < 100000; i++) {
        rng.exponential(1.0);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // Should complete in reasonable time (< 1 second)
    std::cout << "\n    [100k samples in " << duration.count() << " ms]";
    ASSERT_LT(duration.count(), 1000);
}

TEST(test_fel_performance_large_queue) {
    std::priority_queue<Event> FEL;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Add 10k events
    RNG rng(12345);
    for (int i = 0; i < 10000; i++) {
        double time = rng.exponential(1.0);
        FEL.push({ARRIVAL, time, i});
    }
    
    // Pop all
    while (!FEL.empty()) {
        FEL.pop();
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "\n    [10k FEL operations in " << duration.count() << " ms]";
    ASSERT_LT(duration.count(), 1000);
}

// ============================================================================
// MAIN TEST RUNNER
// ============================================================================

int main() {
    std::cout << "\n";
    std::cout << "========================================\n";
    std::cout << "  M/M/1 DES - UNIT TESTS (ANGGOTA 1)\n";
    std::cout << "  Testing: RNG, Event, State, FEL\n";
    std::cout << "========================================\n\n";
    
    // Section 1: Event Tests
    std::cout << "=== SECTION 1: EVENT STRUCT ===\n";
    RUN_TEST(test_event_min_heap_property);
    RUN_TEST(test_event_comparison_order);
    RUN_TEST(test_event_same_time);
    
    // Section 2: State Tests
    std::cout << "\n=== SECTION 2: STATE STRUCT ===\n";
    RUN_TEST(test_state_default_constructor);
    RUN_TEST(test_state_reset_after_modification);
    RUN_TEST(test_state_arrival_queue_operations);
    
    // Section 3: RNG Tests
    std::cout << "\n=== SECTION 3: RNG CLASS ===\n";
    RUN_TEST(test_rng_constructor);
    RUN_TEST(test_rng_deterministic_with_same_seed);
    RUN_TEST(test_rng_different_seeds_different_values);
    RUN_TEST(test_rng_exponential_always_positive);
    RUN_TEST(test_rng_exponential_mean_rate_0_5);
    RUN_TEST(test_rng_exponential_mean_rate_1_0);
    RUN_TEST(test_rng_exponential_mean_rate_2_0);
    RUN_TEST(test_rng_exponential_higher_rate_smaller_mean);
    RUN_TEST(test_rng_exponential_variance);
    RUN_TEST(test_rng_exponential_sequence_independence);
    
    // Section 4: FEL Tests
    std::cout << "\n=== SECTION 4: FEL (Priority Queue) ===\n";
    RUN_TEST(test_fel_basic_operations);
    RUN_TEST(test_fel_multiple_events_sorted);
    RUN_TEST(test_fel_clear_operation);
    
    // Section 5: Edge Cases
    std::cout << "\n=== SECTION 5: EDGE CASES ===\n";
    RUN_TEST(test_rng_very_small_rate);
    RUN_TEST(test_rng_very_large_rate);
    RUN_TEST(test_state_many_arrivals_in_queue);
    
    // Section 6: Performance
    std::cout << "\n=== SECTION 6: PERFORMANCE ===\n";
    RUN_TEST(test_rng_performance_100k_samples);
    RUN_TEST(test_fel_performance_large_queue);
    
    // Summary
    std::cout << "\n========================================\n";
    std::cout << "  TEST SUMMARY\n";
    std::cout << "========================================\n";
    std::cout << "Total:  " << totalTests << "\n";
    std::cout << "Passed: " << passedTests << " (" 
              << (100.0 * passedTests / totalTests) << "%)\n";
    std::cout << "Failed: " << failedTests << "\n";
    
    if (failedTests == 0) {
        std::cout << "\n✓✓✓ ALL TESTS PASSED! ✓✓✓\n";
        std::cout << "ANGGOTA 1 implementation is CORRECT!\n\n";
        return 0;
    } else {
        std::cout << "\n✗✗✗ SOME TESTS FAILED ✗✗✗\n";
        std::cout << "Please fix the failed tests above.\n\n";
        return 1;
    }
}