// Implement the recursive algorithm Fib(k) to computer the kth fib# and compute the number of additions A(k) employed by the algorithm.

#include <iostream>
#include <utility>
#include <fstream>

auto Fib = [](auto&& self, int k) -> std::pair<long long, long long>{
    if (k == 0) {
        return {0, 0};
    }

    if (k == 1) {
        return {1, 0};
    }

    auto left = self(self, k-1);
    auto right = self(self, k-2);
    long long kthfib = left.first + right.first;
    long long additions = left.second + right.second + 1;
    return {kthfib,additions};
};

template <typename T>
int Scatter(int k, T&& Func)
{
    std::ofstream file("Fib_results.csv");
    if (!file){
        std::cout << "Error opening file." << std::endl;
        return 1;
    }

    file << "k, F(k), A(k)\n";
    for (size_t i = 0; i <= k; i++){
        std::pair<long long, long long> result = Func(Func, i);
        file << i << "," << result.first << "," << result.second << "\n";
    }

    file.close();
    std::cout << "Results written to Fib_results.csv" << std::endl;
    return 0;
}

int main(){

    int choice = 0;
    int k;

    std::cout << "Fib(k). Enter 0 for user mode, 1 for scatterplot mode: ";
    std::cin >> choice;

    if (!choice){
        std::cout << "Calculating the fibonacci number k. Please choose k, 0 <= k <= 92 : ";
        std::cin >> k;

        if (k < 0 || k > 92){
            std::cout << "0 <= k <= 92 only. Terminating program." << std::endl;
            return 1;
        }

        std::pair<long long, long long> result = Fib(Fib, k);
        std::cout << "Results written to console: k, Fib(k), A(k)." << std::endl << k << "," << result.first << "," << result.second << "\n";
        return 0;
    } else {
        std::cout << "Calculating the fibonacci sequence up to k for scatterplot data. Please choose k, 0 <= k <= 92 : ";
        std::cin >> k;

        if (k < 0 || k > 92){
            std::cout << "0 <= k <= 92 only. Terminating program." << std::endl;
            return 1;
        }

        Scatter(k, Fib);

    }
}

