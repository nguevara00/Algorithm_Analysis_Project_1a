// Nick Guevara and Will Dappen
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

void fibScatter(){}
void expScatter(){}
void sortScatter(){}


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
        std::cout << "Scatter plot mode. Testing Fib, GCD, exponetiation, and sorts." <<std::endl << std::endl;
        //fibScatter();
        //expScatter();
        //sortScatter();
        std::cout << "Finished scatter plot mode! Please run the program again to test the other mode. Goodbye." << std::endl << std::endl;
    }

    return 0;
}

