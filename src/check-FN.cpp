#include <algorithm>
#include <fstream>
#include <sstream>

#include "utils.cpp"

int main(int argc, char *argv[]) {

  if(argc != 4 and argc != 5){
    cout << "Usage: " << argv[0] << " data-matrix tree-clades FN-count-file *annot-tree-file" << endl;
    return 0;
  }

  // Read in data matrix, scrapping the two-line header
  matrix m = getMatrixFromFile(argv[1], 2);

  // Read in tree clades
  matrix treeClades = getMatrixFromFile(argv[2], 0);
  cout << "input read!" << endl;

  // Was the tree even a solution? If so, get the number of FNs
  int FNCount = -1;
  treeLabelling treeLabel;
  if(treeClades[0][0] == -1){ 
    cout << "INVALID TREE: none provided" << endl; 
  }else{
    FNCount = checkValidSolution(m, treeClades, treeLabel, true);
  }
  cout << "Solution is valid" << endl;

  // Get the total number of 0s (possible FNs)
  int numZeros = 0;
  for(vector<int>& r : m){
    for(int& c : r){ numZeros += (1 - c); }
  }

  // Write results to the FN-counter file
  ofstream outFile; outFile.open(argv[3]); 

  outFile << FNCount << endl;
  outFile << ((double)(FNCount)) / numZeros << endl;
  cout << "Minimum FN count " << FNCount << endl;
  cout << "FN frequency is " << ((double)(FNCount)) / numZeros << endl;
  outFile.close();

  // Write to the tree-annotation file for graphviz (edge labels and node counts), if one was provided
  if(argc == 5){
    cout << "Writing to file:" << argv[4] << endl;

    ofstream treeAnnotOutFile; treeAnnotOutFile.open(argv[4]); 
    treeAnnotOutFile << "Digraph g {" << endl;
    for(const pair<int, int>& pair : treeLabel.toParent){
      int nodeID = pair.first; set<int>& taxa = treeLabel.nodeLabeling[nodeID];

      // Write out each node and the number of taxa in it
      treeAnnotOutFile << nodeID << " [ label = \" " << taxa.size() << " \" ]; " << endl;
    }

    for(const pair<int, int>& pair : treeLabel.toParent){
      int nodeID = pair.first; int parent = treeLabel.toParent[nodeID];

      // Write the basics of each edge
      treeAnnotOutFile << parent << " -> " << nodeID; 

      // Label the edge with the characters gained or lost
      treeAnnotOutFile << " [ label = \" ";
      set<int>& gainsLosses = treeLabel.edgeLabeling[nodeID];
      for(int gainLoss : gainsLosses){
        treeAnnotOutFile << gainLoss << " ";
      } 
      treeAnnotOutFile << "\" ]; " << endl;
    }
    treeAnnotOutFile << "}" << endl;
    treeAnnotOutFile.close();
  }


}

