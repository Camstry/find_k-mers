#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <cctype>

double calculateTm(const std::string& primer) {
    int g = 0, c = 0, a = 0, t = 0;
    
    for (char nucleotide : primer) {
        switch (std::tolower(nucleotide)) {
            case 'g': g++; break;
            case 'c': c++; break;
            case 'a': a++; break;
            case 't': t++; break;
        }
    }
    
    double Tm = 64.9 + ((41.0 * (g + c - 16.4)) / (g + a + c + t));
    return Tm;
}

void processGenome(const std::string& filename, int K, int P, double TL, double TU) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    std::string genome;
    std::string line;
    
    // Read the entire file into a single string
    while (std::getline(file, line)) {
        genome += line;
    }

    file.close();

    // Create a hash table to store k-mer counts
    std::unordered_map<std::string, int> kmerCounts;

    // First pass: count the occurrences of each k-mer
    for (size_t i = 0; i <= genome.length() - K; ++i) {
        std::string kmer = genome.substr(i, K);
        kmerCounts[kmer]++;
    }
    
    // Open output file
    std::string outputFilename = filename + ".kmer";
    std::ofstream outFile(outputFilename);
    if (!outFile.is_open()) {
        std::cerr << "Error opening output file: " << outputFilename << std::endl;
        return;
    }

    // Write search parameters to file    
    outFile << "K-mer length: " << K << " Primer length: " << P 
            << " Upper Tm: " << TU << " Lower Tm: " << TL << std::endl;
    
    // Valid k-mer counter
    int kmerFound = 0;

    // Second pass: process unique k-mers
    for (size_t i = 0; i <= genome.length() - K; ++i) {
        std::string kmer = genome.substr(i, K);

        // Skip non-unique k-mers
        if (kmerCounts[kmer] > 1) {
            continue;
        }
        
        // Check if the k-mer is shorter than 2*P
        if (kmer.length() < 2 * P) {
            continue;
        }

        size_t startPosition = i;
        size_t endPosition = i + K - 1;

        // Extract the iPrimer and fPrimer from the k-mer
        std::string iPrimer = kmer.substr(0, P);
        std::string fPrimer = kmer.substr(K - P);

        // Check if iPrimer or fPrimer starts or ends with 'a' or 't'
        char start_iPrimer = std::tolower(iPrimer.front());
        char end_iPrimer = std::tolower(iPrimer.back());
        char start_fPrimer = std::tolower(fPrimer.front());
        char end_fPrimer = std::tolower(fPrimer.back());

        if (start_iPrimer == 'a' || start_iPrimer == 't' ||
            end_iPrimer == 'a' || end_iPrimer == 't' ||
            start_fPrimer == 'a' || start_fPrimer == 't' ||
            end_fPrimer == 'a' || end_fPrimer == 't') {
            continue;
        }

        // Calculate Tm for iPrimer and fPrimer
        double Tm_iPrimer = calculateTm(iPrimer);
        double Tm_fPrimer = calculateTm(fPrimer);

        // Skip k-mers with a primer Tm outside the specified range
        if (Tm_iPrimer < TL || Tm_iPrimer > TU || Tm_fPrimer < TL || Tm_fPrimer > TU) {
            continue;
        }
        
        // Print results
        std::cout << Tm_iPrimer << " : " 
                  << startPosition << " - " << iPrimer << "..."
                  << fPrimer << " - " << endPosition << " : "
                  << Tm_fPrimer << std::endl;
        outFile << "Tm: " << Tm_iPrimer << " | " << startPosition << " - "
                << iPrimer << "..." << fPrimer << " - " << endPosition
                << " | Tm:" << Tm_fPrimer << std::endl;
        kmerFound++;
    }

    // Print number of k-mers found
    std::cout << "K-mers identified: " << kmerFound << std::endl;
    outFile << "K-mers identified: " << kmerFound << std::endl;
    
    outFile.close();
}

int main() {
    std::string filename;
    int K, P;
    double TL, TU;

    std::cout << "Enter the filename containing the genome: ";
    std::cin >> filename;
    
    std::cout << "Enter the integer K for the length of k-mers: ";
    std::cin >> K;

    std::cout << "Enter the integer P for the length of the primers to display: ";
    std::cin >> P;
    
    std::cout << "Enter the lower bound for the melting temperature TL: ";
    std::cin >> TL;

    std::cout << "Enter the upper bound for the melting temperature TU: ";
    std::cin >> TU;

    if (K <= 0 || P <= 0) {
        std::cerr << "Invalid K or P value. Please enter positive integers." << std::endl;
        return 1;
    }
    if (TL > TU) {
        std::cerr << "Invalid bounds for melting temperature. TL must be less than or equal to TU." << std::endl;
        return 1;
    }
    
    processGenome(filename, K, P, TL, TU);

    std::cout << "Processing complete. Press Enter to close the console." << std::endl;
    std::cin.ignore();
    std::cin.get();
    return 0;
}