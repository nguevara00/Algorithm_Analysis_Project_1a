CS415 – Project 1A

Asymptotic Analysis of a Mélange of Algorithms
Nick Guevara & William Dappen
Professor Gill
February 28, 2026

Overview

This project analyzes the asymptotic complexity of several classical algorithms by:
 - Counting fundamental operations
 - Generating scatter plot data (CSV format)
 - Comparing empirical growth to theoretical recurrence solutions

The project consists of three main tasks:
 - Fibonacci & Euclid’s GCD 
 - Exponentiation Algorithms 
 - Sorting Algorithms (Insertion & Selection Sort)

The program supports two execution modes:
 - User Testing Mode 
 - Scatter Plot Mode

Compilation

Compile using a C++17-compatible compiler:
Run: g++ -std=c++17 project1a.cpp -o project1a

Then run:
./project1a

Required Folder Structure:

project/
│
├── project1a.cpp
├── smallSet/
│     data10.txt
│     ...
├── testSet/
│     data100.txt
│     data100_sorted.txt
│     data100_rSorted.txt
│     ...

Program Modes
1. User Testing Mode (Enter 0)
 - This mode prompts the user for specific inputs and prints results directly to the console. 

Task 1
 - Prompts for k (0 ≤ k ≤ 91)
Outputs:
 - Fib(k)
 - A(k) (number of additions in recursive Fib)
 - GCD(Fib(k+1), Fib(k))
 - Number of modulo operations

Task 2
 - Prompts for base a and exponent n
Outputs:
 - Result of a^n using:
   - Decrease-by-one 
   - Decrease-by-constant-factor 
   - Divide-and-conquer
Multiplication counts are internally tracked

Be warned, larger bases overflow earlier than x^63

However, multiplication counts remain correct

Task 3
 - Prompts for list size n (10–100, increments of 10)
 - Loads data from smallSet/
Outputs sorted arrays using:
 - Insertion Sort 
 - Selection Sort

Scatter Plot Mode (Enter 1)

This mode generates CSV files for plotting externally (Excel, Python, D3, etc.).

No user input is required.

**IMPORTANT**
    
    To view the produced CSV files, open the file:

        csv display.html

    in a Javascript compatible web browser interface, and import the csv files using
    the **Browse...** button. In the case of Task 3, different views can be observed 
    via the "Y metric" dropdown.

Generated CSV Files
    Task 1 – Fibonacci & GCD

    fib_scatter_data.csv

    gcd_scatter_data.csv

Task 2 – Exponentiation
    
    exp_scatter_data.csv

Task 3 – Sorting

    sort_scatter_data.csv

This file has three implementation groups:
 - best_case
 - average_case 
 - worst_case

 