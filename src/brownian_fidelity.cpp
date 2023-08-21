/*
 * File Name: brownian_fidelity
 * Description: calculation the fidelity of the ground state of brownian model
 * Created by Hao-Xin on 2023/07/25;
 * Usage: ./fidelity params.txt
 * For the form of params.txt please refer to fidelity_parameters.txt
 */



#include "brownian_type.h"
#include "params_case.h"
#include "myutil.h"
#include <filesystem> // C++17 header for directory manipulation
#include <dirent.h> // For POSIX directory manipulation

using namespace std;

// Function to convert a double to a string without trailing zeros
std::string double_to_string(double value) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(10) << value; // Use a high precision initially
  std::string str = oss.str();

  // Remove trailing zeros and the decimal point if there are no decimals
  size_t pos = str.find_last_not_of('0');
  if (pos != std::string::npos && str[pos] == '.')
    pos--; // Remove the decimal point if no decimals

  return str.substr(0, pos + 1);
}

int main(int argc, char *argv[]) {
  std::string path;
  size_t L;
  double Delta;
  double h;
  std::vector<std::pair<double, double>> gamma_pairs;

  // Check if the parameter file is provided as an argument
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <parameter_file_path>" << std::endl;
    return 1;
  }

  // Open the parameter file
  std::ifstream paramFile(argv[1]);
  if (!paramFile) {
    std::cerr << "Error opening parameter file: " << argv[1] << std::endl;
    return 1;
  }

  // Read the path from the first line of the parameter file
  if (!(paramFile >> path)) {
    std::cerr << "Error reading path from the parameter file." << std::endl;
    return 1;
  }

  // Read parameters from the parameter file
  std::string paramName;
  while (paramFile >> paramName) {
    if (paramName == "L") {
      paramFile >> L;
    } else if (paramName == "Delta") {
      paramFile >> Delta;
    } else if (paramName == "h") {
      paramFile >> h;
    } else if (paramName == "gamma_pairs") {
      double x, y;
      paramFile >> x >> y;
      gamma_pairs.emplace_back(x, y);
    } else {
      std::cerr << "Unknown parameter: " << paramName << std::endl;
    }
  }

  // Close the parameter file
  paramFile.close();

  SiteVecT sites(L, pb_out);
  std::vector<TenElemT> fidelity(gamma_pairs.size());
  for (size_t i = 0; i < gamma_pairs.size(); i++) {
    double gamma1 = gamma_pairs[i].first;
    double gamma2 = gamma_pairs[i].second;
    std::cout << "group " << i << ": (" << gamma1 << ", " << gamma2 << "), \t ";
    std::string mps_path1, mps_path2;
    std::string mps_path1_match, mps_path2_match;

    mps_path1_match = "mpsbrownianL" + std::to_string(L)
        + "Delta" + double_to_string(Delta) + "h" + double_to_string(h)
        + "gamma" + double_to_string(gamma1) + "D";

    mps_path2_match = "mpsbrownianL" + std::to_string(L)
        + "Delta" + double_to_string(Delta) + "h" + double_to_string(h)
        + "gamma" + double_to_string(gamma2) + "D";

    int max_D_mps1 = 0;
    int max_D_mps2 = 0;

    // Traverse the "path" directory to find the maximum D for mps1 and mps2
    // C++17 realization
    /*
    for (const auto &entry: std::filesystem::directory_iterator(path)) {
      if (entry.is_directory()) {
        std::string filename = entry.path().filename().string();
        if (filename.find(mps_path1_match) != std::string::npos) {
          // Extract the D from the filename
          int D = std::stoi(filename.substr(filename.find("D") + 1));
          max_D_mps1 = std::max(max_D_mps1, D);
        } else if (filename.find(mps_path2_match) != std::string::npos) {
          // Extract the D from the filename
          int D = std::stoi(filename.substr(filename.find("D") + 1));
          max_D_mps2 = std::max(max_D_mps2, D);
        }
      }
    }
     */
    //before C++17
    DIR *dir = opendir(path.c_str());
    if (!dir) {
      // Handle error opening the directory
    }

    dirent *entry;
    while ((entry = readdir(dir))) {
      if (entry->d_type == DT_DIR) {
        std::string filename = entry->d_name;
        if (filename != "." && filename != ".." && filename.find(mps_path1_match) != std::string::npos) {
          // Extract the D from the filename
          int D = std::stoi(filename.substr(filename.rfind("D") + 1));  //rfind: the last string "D"
          max_D_mps1 = std::max(max_D_mps1, D);
        } else if (filename != "." && filename != ".." && filename.find(mps_path2_match) != std::string::npos) {
          // Extract the D from the filename
          int D = std::stoi(filename.substr(filename.rfind("D") + 1));
          max_D_mps2 = std::max(max_D_mps2, D);
        }
      }
    }

    if (max_D_mps1 == 0) {
      std::cout << "warning: do not find the mps for gamma = " << gamma1 << std::endl;
    }

    if (max_D_mps2 == 0) {
      std::cout << "warning: do not find the mps for gamma = " << gamma2 << std::endl;
    }

    std::cout << "(D1, D2) = (" << max_D_mps1 << ", " << max_D_mps2 << "), \t ";

    closedir(dir);

    // Construct the final paths for mps1 and mps2 using the maximum D
    mps_path1 = path + "/" + mps_path1_match + std::to_string(max_D_mps1);
    mps_path2 = path + "/" + mps_path2_match + std::to_string(max_D_mps2);

    FiniteMPST mps1(sites), mps2(sites);
    mps1.Load(mps_path1);
    mps2.Load(mps_path2);
    TenElemT overlap = FiniteMPSInnerProd(mps1, mps2);
    fidelity[i] = std::norm(overlap);
    std::cout << "fidelity = " << fidelity[i] << std::endl;
  }


  // Open the file to store the results
  std::ofstream
      outFile("fidelityL" + to_string((L)) + "Delta" + double_to_string(Delta) + "h" + double_to_string(h) + ".csv");
  if (!outFile) {
    std::cerr << "Error opening output file: fidelity.csv" << std::endl;
    return 1;
  }

  // Write the header for the CSV file
  outFile << "L,Delta,h,gamma1,gamma2,Fidelity" << std::endl;

  // Store the common parameters (L, Delta, h) to be written only once
  outFile << L << "," << Delta << "," << h << std::endl;

  // Store the corresponding gamma pairs and fidelity values
  for (size_t i = 0; i < gamma_pairs.size(); i++) {
    double gamma1 = gamma_pairs[i].first;
    double gamma2 = gamma_pairs[i].second;

    // Calculate the squared fidelity
    double fidelity_squared = fidelity[i];

    // Write the data to the "fidelity" file
    outFile << gamma1 << "," << gamma2 << "," << fidelity_squared << std::endl;
  }

  // Close the "fidelity" file
  outFile.close();

  return 0;
}