// Nick Guevara and William Dappen
// CS415 Gill
// Project 1a
// 2-28-26

// Task 1a: Implement the recursive algorithm Fib(k) to compute the kth Fibonacci Number and compute the corresponding number of additions A(k) employed by the recursive algorithm. 
//          Produce a scatterplot of A(k) and indicate the algorithm’s asymptotic complexity. Note that A(k) is a function of k; your program must choose the different values of k to generate the scatterplot. 
// Task 1b: Compute the worst-case efficiency of Euclid's algorithm to compute GCD. Note that the worst-case inputs for Euclid’s algorithm happen to be the consecutive elements of the Fibonacci sequence
//          (you may store these from Task 1a). So, compute GCD(m, n) where m = Fib(k+1) and n = Fib(k) for several values of k >= 1. 
//          produce a scatterplot showing the number of modulo divisions D(n) taken to compute GCD as a function of n and indicate the algorithm’s asymptotic complexity. 
//          Note that D(n) is a function of n; your program must choose the different values of n (obtained by choosing different values of k) to generate the scatterplot.
// Task 2:  The problem of exponentiation exp(a, n) computes an , where a is a constant. Implement the following algorithms for solving the problem of exponentiation:
//          Decrease by one, decrease by constant factor, divide and conquer. 
//          Produce a scatterplot showing the number of multiplications M(n) made by each of the algorithms and indicate their asymptotic complexity. 
//          Note that M(n) is a function of n; your program must choose the different values of n to generate the scatterplot. You can choose any value of a.
// Task 3:  Implement and analyze the complexity of Selection Sort and Insertion Sort. 
//          Run your code on test data Download test data that is sorted, random, and reverse sorted, which corresponds to the best case, average case, and worst case, respectively.
//          Produce a scatterplot showing the number of comparisons C(n) made by each of the algorithms and indicate their asymptotic complexity.  Note that n here is the size of the list (that is being sorted). 

/* Other Requirements:
    User testing mode: This mode prompts the user for values of certain parameters and returns an output based on the goals of each task:

    Task 1: User is prompted for the value of k; program outputs the value of Fib(k) and GCD(m, n) where m = Fib(k+1) and n = Fib(k) .
    Task 2: User is prompted for the value of a and n.; program outputs the value an using each of the three algorithms.
    Task 3: User is prompted for the size of the list n (between 10 - 100 at an increment of 10), the program loads the appropriate file data{n}.txt from the folder "smallSet" in the test data Download test data 
            and outputs the sorted output based on the two algorithms on terminal.
    
    Scatter plot mode: This mode does NOT prompt the user for any values. It generates a scatter plot for each task, as described in the description of each of the tasks:

    Task 1: Fib(k) for different values of k generated using the recursive algorithm may be stored for this purpose. OR, for sake of speed, you may use an iterative algorithm to compute the Fibonacci sequence once. 
            Use consecutive elements as input to the GCD(m, n) function for computing and plotting D(n). (Think about what should be the upper bound of k and clearly indicate it in your report)? 
    Task 2: You may choose any value of the constant a. Compute and plot M(n) for the three different algorithms in the same plot (it will help you to compare them!).
    Task 3: For different sizes of the list n (range 10 - 10000), read data from the folder "testSet" in the test data Download test data that is sorted (data{n}_sorted.txt), random (data{n}.txt), and reverse sorted (data{n}_rSorted.txt). 
            Use that data to compute C(n) for each of the two sorting algorithms. Produce three scatter plots that compare the complexity of the two algorithms in i) Best-case, ii) Average-case, and iii) Worst-case.

*/

#include <iostream>
#include <utility>
#include <fstream>
#include <vector>
#include <numeric> // Required for std::iota
#include <cstddef> // Required for std::size_t
#include <iterator>
#include <sstream>

//computes Fib(k) with naive recursive algorithm
//input k is an integer between 0 and 91
//outputs the kth number in the fibonacci sequence
std::pair<long long, long long> Fib(int k){ 
    if (k == 0) {
        return {0, 0};
    }

    if (k == 1) {
        return {1, 0};
    }

    auto left = Fib(k-1);
    auto right = Fib(k-2);
    long long kthfib = left.first + right.first;
    long long additions = left.second + right.second + 1;
    return {kthfib,additions};
}

long long FibIt(int n){
    if (n == 0) return 0;
    if (n == 1) return 1;
    long long a=0, b=1, c;

    for (int i = 2; i <= n; ++i) {
        c = a + b;
        a = b;
        b = c;
    }
    return b;
}

//Euclid's Algorithm to compute GCD.
//for worst case analysis it needs to take in fibonacci numbers from the fib function
std::pair<long long, long long> GCD(long long m, long long n){
    if(n == 0){
        return {m, 0};
    }
    
    auto result = GCD(n, m % n);
    return {result.first, result.second + 1};
}

//Exponentiation Algorithm, implemented as Decrease-by-one
//input a is a constant base. n is the power.
//returns a pair, the first is the evaluation of the exponentiation, the second is the number of multiplications. 
std::pair<long long, long long> EXPI(long long a, long long n){
    if (n == 0){
        return {1, 0};
    }
    std::pair<long long, long long> recursiveCall = EXPI(a, n-1);
    return {recursiveCall.first * a, recursiveCall.second + 1}; //one multiplication added to count
}


//Exponentiation algorithm, decrease-by-constant-factor
//input a is a constant base. n is the power.
//returns a pair, the first is the evaluation of the exponentiation, the second is the number of multiplications. 
std::pair<long long, long long> EXPII(long long a, long long n){
    if (n == 0){
        return {1, 0};
    }
    std::pair<long long, long long> recursiveCall = EXPII(a, n/2);
    long long half = recursiveCall.first;

    if(n % 2){
        return {a * half * half, recursiveCall.second + 2}; //if n is odd, two multiplications added to count
    }

    return {half*half, recursiveCall.second +1}; //if n is even, one multiplication added to count

}

//Exponentiation algorithm, Divide-and-Conquer
//input a is a constant base. n is the power.
//returns a pair, the first is the evaluation of the exponentiation, the second is the number of multiplications. 
std::pair<long long, long long> EXPIII(long long a, long long n){
    if(n == 0){
        return {1, 0};
    }
    std::pair<long long, long long> recursiveCall1 = EXPIII(a, n/2);
    std::pair<long long, long long> recursiveCall2 = EXPIII(a, n/2);
    if(n % 2){
        return {a * recursiveCall1.first * recursiveCall2.first, recursiveCall1.second + recursiveCall2.second + 2};
    }
    return {recursiveCall1.first * recursiveCall2.first, recursiveCall1.second + recursiveCall2.second + 1};
}


long long InSort(std::vector<long long> &data){
    long long comparisons = 0; // to output
    for(size_t i = 1; i < data.size(); i++) {
        long long j = static_cast<long long>(i) - 1; //highest index of the subarray which is already sorted
        long long key = data.at(i);
        while (j >= 0 && data.at(j) > key) {
            comparisons++;
            data.at(j + 1) = data.at(j); //shift data forward if it is greater than the key
            j--;
        }
        if(j>=0){
            comparisons++; //failed comp
        }

        data.at(j + 1) = key;
    }
    return comparisons;
}

long long SelSort(std::vector<long long> &data){
    if (data.size() < 2) {
        return 0;
    }

    long long comparisons = 0;
    for(size_t i = 0; i<data.size()-1; i++){
        size_t minIdx = i;
        for(size_t j = i+1; j<data.size(); j++){
            if(data.at(j) < data.at(minIdx)){
                minIdx = j;
            }
            comparisons++;
        }
        std::swap(data.at(i), data.at(minIdx));
    }
    return comparisons;
}

void fibUser(){
    int k = -1;
    
    while (k < 0 || k > 91){
        std::cout << "Calculating the kth fibonacci number." << std::endl << "Please choose k, 0 <= k <= 91 (will be very slow for k > 40) : ";
        std::cin >> k;
        std::cout << std::endl;
    }

    std::pair<long long, long long> result1 = Fib(k); // fib(k)
    std::pair<long long, long long> result2 = Fib(k+1); // fib(k+1)
    std::pair<long long, long long> gcdResult = GCD(result2.first, result1.first); // GCD (m,n) 

    std::cout << "Results written to console: k, Fib(k), A(k)." << std::endl << k << "," << result1.first << "," << result1.second << std::endl << std::endl;
    std::cout << "Calculating GCD for worst case inputs: " << std::endl;
    std::cout << "m = Fib(k+1) = " << result2.first << std::endl;
    std::cout << "n = Fib(k) = " <<result1.first << std::endl;
    std::cout << "GCD(m,n) = " << gcdResult.first << std::endl;
    std::cout << "GCD basic modulo operations = " << gcdResult.second << std::endl << std::endl;
}

void expUser(){
    int a = 0;
    int n = 0;
    std::cout << "Calculating a^n using three methods." << std::endl << "Please enter base a: ";
    std::cin >> a;
    std::cout << "The result of a^n may overflow. For a=2, choose n less than 63. For higher bases, please choose n to be a small number. M(a,n) will be correct for any numbers. Please enter exp n: ";
    std::cin >> n;
    std::cout << std::endl;

    std::pair<long long, long long> expOne = EXPI(a,n);
    std::pair<long long, long long> expTwo = EXPII(a,n);
    std::pair<long long, long long> expThree= EXPIII(a,n);
    std::cout << "Exponentiation using Decrease-by-one. Result = " << expOne.first << std::endl << std::endl;
    std::cout << "Exponentiation using Decrease-by-constant-factor. Result = " << expTwo.first <<  std::endl << std::endl;
    std::cout << "Exponentiation using Divide-and-conquer. Result = " << expThree.first <<  std::endl << std::endl;
}


void sortUser(){
    int n = 0;
    while (n < 10 || n > 100 || (n % 10)){
        std::cout << "Sorting a list using insertion sort and selection sort." << std::endl << "Please choose n, the size of the list to be sorted, 10 <= n <= 100, in increments of 10: ";
        std::cin >> n;
        std::cout << std::endl;
    }

    //read data from input txt file to construct vector to be sorted
    const std::string filename = "smallSet/data" + std::to_string(n) + ".txt"; // create a path
    std::ifstream file(filename);
    if (!file) {
        std::cout << "Error: could not open file: " << filename << std::endl;
        return;
    }

    std::vector<long long> data;
    data.reserve(n);
    long long x = 0;
    while (file >> x) {
        data.push_back(x);
    }
    if (data.size() != static_cast<size_t>(n)){
        std::cout << "Warning: expected " << n << " values, read " << data.size() << std::endl;
    }

    file.close();

    std::vector<long long> insertionSortData = data;
    std::vector<long long> selectionSortData = data;

    InSort(insertionSortData);
    std::cout << "List after Insertion Sort: ";
    for (size_t i = 0; i < insertionSortData.size(); i++){
        std::cout << insertionSortData[i] << (i+1 == insertionSortData.size() ? '\n' : ' ');
    }
    std::cout << std::endl;

    SelSort(selectionSortData);
    std::cout << "List after Selection Sort: ";
    for (size_t i = 0; i < selectionSortData.size(); i++){
        std::cout << selectionSortData[i] << (i+1 == selectionSortData.size() ? '\n' : ' ');
    }
    std::cout << std::endl;
}

std::vector<long long> fileToVector(const std::string& filename){
    std::ifstream inputFile(filename + ".txt");

    if (!inputFile)
        throw std::runtime_error("Failed to open file: " + filename);

    std::vector<long long> integers;
    int value;

    while (inputFile >> value) {
        integers.push_back(value);
    }

    return integers;
}

using Matrix = std::vector<std::vector<long long>>; // matrix[col][row]
std::string implRowsToCSV(
        const Matrix& m,
        const std::string& implName,
        const std::vector<long long>& x
) {
    const size_t rows = x.size();
    const size_t cols = m.size();

    if (rows == 0) {
        std::cerr << "implRowsToCSV: x is empty.\n";
        std::exit(-100);
    }
    if (cols == 0) {
        std::cerr << "implRowsToCSV: matrix has no columns.\n";
        std::exit(-100);
    }

    for (size_t c = 0; c < cols; ++c) {
        if (m[c].size() != rows) {
            std::cerr << "implRowsToCSV: column " << c
                      << " size (" << m[c].size()
                      << ") != x size (" << rows << ").\n";
            std::exit(-100);
        }
    }

    std::ostringstream out;
    for (size_t i = 0; i < rows; ++i) {
        out << implName << "," << x[i];
        for (size_t c = 0; c < cols; ++c) {
            out << "," << m[c][i];
        }
        out << "\n";
    }
    return out.str();
}

void createCSV(
        const std::vector<std::string>& titles,
        const std::vector<Matrix>& data,
        const std::vector<std::string>& implementations,
        const std::vector<long long>& x,
        const std::string& filename
) {
    if (data.size() != implementations.size()) {
        std::cerr << "createCSV: data.size() must match implementations.size().\n";
        std::exit(-1000);
    }
    if (titles.empty()) {
        std::cerr << "createCSV: titles is empty.\n";
        std::exit(-1000);
    }
    if (x.empty()) {
        std::cerr << "createCSV: x is empty.\n";
        std::exit(-1000);
    }

    std::ofstream out(filename + ".csv");
    if (!out) {
        std::cerr << "createCSV: failed to open output file.\n";
        std::exit(-1000);
    }

    // Header: impl,x,metric1,metric2,...
    out << "impl,";
    for (size_t i = 0; i + 1 < titles.size(); ++i) out << titles[i] << ",";
    out << titles.back() << "\n";

    // Body
    for (size_t k = 0; k < data.size(); ++k) {
        out << implRowsToCSV(data[k], implementations[k], x);
    }
}




//Task 1: Fib(k) for different values of k generated using the recursive algorithm may be stored for this purpose. OR, for sake of speed, you may use an iterative algorithm to compute the Fibonacci sequence once.
//Use consecutive elements as input to the GCD(m, n) function for computing and plotting D(n). (Think about what should be the upper bound of k and clearly indicate it in your report)?
void fibScatter(){
    int k = 92;

    std::vector<long long> FibList;
    FibList.reserve(k);
    for(int i = 0; i<=k; i++){
        FibList.push_back(FibIt(i)); // produce a vector of Fibonacci numbers
    }
    std::vector<long long int> CompList;
    CompList.reserve(k);
    for(int j = 1; j < FibList.size(); j++){
        std::pair<long long,long long> result = GCD(FibList[j], FibList[j-1]);
        //std::cout << "GCD(" << FibList.at(j) << ", " << FibList.at(j-1) << ")" << ": " << result.first << ", " << result.second << std::endl;
        CompList.push_back(result.second);
    }
    std::vector<long long int> x(k);
    std::iota(std::begin(x), std::end(x), 0);

    Matrix m;
    m.push_back(CompList);

    std::vector<std::string> titles = { "GCD(fib(n) fib(n-1))", "Number of Comparisons" };

    std::vector<Matrix> data = { m };

    std::vector<std::string> implementations = { "fib_gcd" };

    std::string filename = "fib_scatter_data";

    createCSV(titles, data, implementations, x, filename);
    std::cout << "Wrote " << filename << ".csv\n";
}

//Task 2: You may choose any value of the constant a. Compute and plot M(n) for the three different algorithms in the same plot (it will help you to compare them!).
void expScatter() {
    int n=2000;
    // x column: 1..n
    std::vector<long long> x(static_cast<size_t>(n));
    std::iota(x.begin(), x.end(), 1);

    // Each implementation gets its own Matrix with ONE metric column
    Matrix expI(1), expII(1), expIII(1);
    expI[0].reserve(n);
    expII[0].reserve(n);
    expIII[0].reserve(n);

    for (int i = 1; i <= n; ++i) {
        auto result1 = EXPI(2, i);
        auto result2 = EXPII(2, i);
        auto result3 = EXPIII(2, i);

        expI[0].push_back(result1.second);
        expII[0].push_back(result2.second);
        expIII[0].push_back(result3.second);
    }

    std::vector<std::string> titles = {"Exponent n which a is raised to", "Number of Multiplications"};

    std::vector<Matrix> data = { expI, expII, expIII };

    std::vector<std::string> implementations = {
            "Decrease-by-one",
            "Decrease-by-constant-factor",
            "Divide-and-Conquer"
    };

    std::string filename = "exp_scatter_data";

    createCSV(titles, data, implementations, x, filename);

    std::cout << "Wrote " << filename << ".csv\n";
}


//Task 3: For different sizes of the list n (range 10 - 10000), read data from the folder "testSet" in the test data Download test data that is sorted (data{n}_sorted.txt), random (data{n}.txt), and reverse sorted (data{n}_rSorted.txt).
//Use that data to compute C(n) for each of the two sorting algorithms. Produce three scatter plots that compare the complexity of the two algorithms in i) Best-case, ii) Average-case, and iii) Worst-case.
void sortScatter() {
    std::vector<long long> x;
    x.reserve(100);
    for (int n = 100; n <= 10000; n += 100){
        x.push_back(n);
    }

    const size_t rows = x.size();

       Matrix best(2), avg(2), worst(2);
    for (auto& col : best)  col.reserve(rows);
    for (auto& col : avg)   col.reserve(rows);
    for (auto& col : worst) col.reserve(rows);

    for (int n = 100; n <= 10000; n += 100) {
        std::string fileNum = std::to_string(n);

        // If your files are in "testSet", include that prefix as shown:
        auto sorted  = fileToVector("testSet/data" + fileNum + "_sorted");
        auto random  = fileToVector("testSet/data" + fileNum);
        auto rsorted = fileToVector("testSet/data" + fileNum + "_rSorted");

        // Best case
        {
            auto a1 = sorted;
            auto a2 = sorted;
            long long c1 = InSort(a1);
            long long c2 = SelSort(a2);
            best[0].push_back(c1);
            best[1].push_back(c2);
        }

        // Average case
        {
            auto a1 = random;
            auto a2 = random;
            long long c1 = InSort(a1);
            long long c2 = SelSort(a2);
            avg[0].push_back(c1);
            avg[1].push_back(c2);
        }

        // Worst case
        {
            auto a1 = rsorted;
            auto a2 = rsorted;
            long long c1 = InSort(a1);
            long long c2 = SelSort(a2);
            worst[0].push_back(c1);
            worst[1].push_back(c2);
        }
    }

    std::vector<std::string> titles = {"Size of the input array",
                                       "Calculations for Insertion Sort",
                                       "Calculations for Selection Sort"};

    // three implementations (blocks)
    std::vector<std::string> implementations = {"best_case", "average_case", "worst_case"};
    std::vector<Matrix> data = {best, avg, worst};

    std::string filename = "sort_scatter_data";
    createCSV(titles, data, implementations, x, filename);

    std::cout << "Wrote " << filename << ".csv\n";
}


int main(){
    //Give the choice of two modes
    int choice = -1;
    while (choice > 1 || choice < 0){
        std::cout << "Project 1a, Asymptotic Analysis of melange of Algorithms" << std::endl << std::endl << "Enter 0 for user testing mode, 1 for scatter plot mode: ";
        std::cin >> choice;
        std::cout << std::endl;
    }

    
    if (!choice){ // User testing mode
        std::cout << "User Mode. Testing Fib, GCD, exponentiation, and sorts." << std::endl << std::endl; 
        fibUser();
        expUser();
        sortUser();
        std::cout << "Finished user testing mode! Please run the program again to test the other mode. Goodbye." << std::endl << std::endl;

       
    } else { // scatter plot mode
        std::cout << "Scatter plot mode. Testing Fib, GCD, exponentiation, and sorts." <<std::endl << std::endl;

                fibScatter();
                expScatter();
                sortScatter();

        std::cout << "Finished scatter plot mode! Please run the program again to test the other mode. Goodbye." << std::endl << std::endl;
    }

    return 0;
}

