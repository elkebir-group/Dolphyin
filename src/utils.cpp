#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <queue>
#include <set>
#include <algorithm>

#include "PQTree.cpp"

using namespace std;

/*
---GENERAL FUNCTIONS---
*/ 

vector<int> set2Vec(set<int>& s){
  vector<int> v;
  for(int i : s){
    v.push_back(i);
  }
  return v;
}

set<int> vec2Set(vector<int>& v){
  set<int> s;
  for(int i : v){
    s.insert(i);
  }
  return s;
}

set<int> getSetIntersection(const set<int> a, const set<int> b){
    vector<int> res;
    set_intersection(a.begin(), a.end(), b.begin(), b.end(), inserter(res, res.begin()));
    return set<int>(res.begin(), res.end());
}

set<int> getSetUnion(const set<int> a, const set<int> b){
    vector<int> res;
    set_union(a.begin(), a.end(), b.begin(), b.end(), inserter(res, res.begin()));
    return set<int>(res.begin(), res.end());
}

set<int> getSetDif(const set<int> a, const set<int> b){
    vector<int> res;
    set_difference(a.begin(), a.end(), b.begin(), b.end(), inserter(res, res.begin()));
    return set<int>(res.begin(), res.end());
}

/*
---EXP. ITERATOR FUNCTIONS---
*/

typedef vector<bool> expIterator;

// Make an iterator
expIterator getIterator(int n){
  vector<bool> it = vector<bool>(n + 1, false);
  return it;
}

// Increment iterator
void incrIterator(expIterator& it){
  int n = it.size();
  for(int i = 0; i < n; i++){
    if(!it[i]){
      it[i] = true;
      return;
    }
    it[i] = false;
  }
}

// Check if iterator is at end
bool iteratorDone(expIterator& it){
  int n = it.size();
  // cout << "Iterator done? " << it[n - 1] << endl;
  return it[n - 1];
}

// Get state of iterator
set<int> getIteratorState(expIterator& it, set<int> elements){
  int n = it.size() - 1;
  set<int> state;
  int i = 0;
  for(int elem : elements){
    if(it[i]){ state.insert(elem); } i++;
  }
  return state;
}

/*
---MATRIX FUNCTIONS---
*/ 

typedef vector<vector<int>> matrix;

// Read a matrix from input
matrix getMatrixFromFile(string inputFileName, int linesToBurn){

  // Find where to read
  ifstream fileData(inputFileName, ios::in);

  // Build B, dataMatrix
  matrix dataMatrix; string line, token;

  // Burn lines (additional info)
  for(int i = 0; i < linesToBurn; i++){ getline(fileData, line); }

  // Read in the important data
  while(getline(fileData, line)){

    vector<int> curRow;
    stringstream str(line);

    while(getline(str, token, ' ')){
      int nextVal = stoi(token);
      curRow.push_back(nextVal);
    }
    dataMatrix.push_back(curRow);
  }
  return dataMatrix;

}

// Undo the effect of a row mapping
matrix expandMatrix(matrix& m, map<int, int> rowMapping){
  matrix newMatrix;
  for(int r = 0; r < m.size(); r++){
    newMatrix.push_back(m[rowMapping[r]]);
  }
  return newMatrix;
}

// Write a matrix to output
void writeMatrixToFile(string outputFileName, matrix& m){

  ofstream outFile;
  outFile.open(outputFileName);

  for(vector<int>& r : m){
    for(int& c : r){
      outFile << c << " ";
    }
    outFile << endl;
  }
  outFile.close();
}

struct Submatrix {             
  set<int> rows;     
  set<int> columns; 
  matrix* matrix; 

  int getVal(int r, int c) const { 
    return (*(matrix))[r][c]; 
  }
};

pair<set<int>, set<int>> getRowCols(matrix m){
  int RSize = m.size();
  int CSize = m[0].size();
  set<int> R;
  for(int i = 0; i < RSize; i++){
    R.insert(i); 
  }
  set<int> C;
  for(int j = 0; j < CSize; j++){
    C.insert(j);
  }
  return {R, C};
}

set<int> extractClade(Submatrix& B, int cladeCol){
  set<int> clade;
  for(int r : B.rows){
    if(B.getVal(r, cladeCol)){
      clade.insert(r);
    }
  }
  return clade;
}

set<int> extractTaxa(Submatrix& B, int taxaRow){
  set<int> row;
  for(int c : B.columns){
    if(B.getVal(taxaRow, c)){
      row.insert(c);
    }
  }
  return row;
}

vector<set<int>> getClades(Submatrix& B, set<int> cladeCols){
  vector<set<int>> clades;
  for(int c : cladeCols){
    clades.push_back(extractClade(B, c));
  }
  return clades;
}

map<int, set<int>> getCladesWithMapping(Submatrix& B, set<int> cladeCols){
  map<int, set<int>> clades;
  for(int c : cladeCols){
    clades[c] = extractClade(B, c);
  }
  return clades;
}

map<int, set<int>> getTaxaWithMapping(Submatrix& B, set<int> taxaRows){
  map<int, set<int>> taxa;
  for(int r : taxaRows){
    taxa[r] = extractTaxa(B, r);
  }
  return taxa;
}

Submatrix extractSubmatrix(Submatrix& initial, set<int> newRows, set<int> newColumns){

  if(getSetDif(newRows, initial.rows).size() > 0 || getSetDif(newColumns, initial.columns).size() > 0){
    return {{-1}};
  }

  Submatrix newSubmatrix;
  newSubmatrix.rows = newRows;
  newSubmatrix.columns = newColumns;
  newSubmatrix.matrix = initial.matrix;
  return newSubmatrix;
}

// Get rid of copies of rows, return the old -> new mapping and additional clades
map<int, int> deleteDuplRows(Submatrix& B){

  map<int, int> oldToNewRows;
  map<int, set<int>> bTaxa = getTaxaWithMapping(B, B.rows);
  set<int> newTaxa;

  // For each row, 
  for(int r : B.rows){
    bool keepTaxa = true;

    // Remove rows that are copies and log the original
    for(int r1 : newTaxa){
      if(bTaxa[r] == bTaxa[r1]){
        keepTaxa = false;
        oldToNewRows[r] = r1;
        break;
      }
    }

    // Else, just keep them
    if(keepTaxa){
      oldToNewRows[r] = r;
      newTaxa.insert(r);
    }
  }

  // Replace rows
  B.rows = newTaxa;
  return oldToNewRows;
}

bool deleteBlankColumns(Submatrix& B){
  set<int> newColumns;
  map<int, set<int>> bClades = getCladesWithMapping(B, B.columns);
  int numTaxa = B.rows.size();

  for(int c : B.columns){
    bool keepColumn = true;

    // Remove columns with 0 1s
    if(bClades[c].size() == 0){
      keepColumn = false;
      continue;
    }
    if(keepColumn){ newColumns.insert(c); }
  }
  bool changed = (B.columns != newColumns);
  B.columns = newColumns;
  return changed;
}

bool deleteTrivialColumns(Submatrix& B){
  set<int> newColumns;
  map<int, set<int>> bClades = getCladesWithMapping(B, B.columns);
  int numTaxa = B.rows.size();

  for(int c : B.columns){
    bool keepColumn = true;

    // Remove columns with less than 2 1s or full columns
    if(bClades[c].size() <= 1 || bClades[c].size() >= numTaxa){
      keepColumn = false;
      continue;
    }
    if(keepColumn){ newColumns.insert(c); }
  }
  bool changed = (B.columns != newColumns);
  B.columns = newColumns;
  return changed;
}

// Delete columns that are duplicates 
bool deleteDuplColumns(Submatrix& B, set<int> cladeCols, map<int, int>& colMapping){

  set<int> newColumns;
  map<int, set<int>> bClades = getCladesWithMapping(B, B.columns);
  int numTaxa = B.rows.size();
  bool deleted = false;

  // Clade columns get priority - if it comes between a clade and non, use clade
  for(int c : cladeCols){
    bool keepColumn = true;

    // Remove columns with identical somewhere else
    for(int c1 : newColumns){
      if(bClades[c] == bClades[c1]){
        keepColumn = false;
        colMapping[c] = c1;
        break;
      }
    }
    if(keepColumn){ newColumns.insert(c); colMapping[c] = c; } else { deleted = true; }
  }

  for(int c : B.columns){
    bool keepColumn = true;

    // Remove columns with identical somewhere else
    for(int c1 : newColumns){
      if(bClades[c] == bClades[c1]){
        keepColumn = false;
        colMapping[c] = c1;
        break;
      }
    }
    if(keepColumn){ newColumns.insert(c); colMapping[c] = c; } else { deleted = true; }
  }

  bool changed = (B.columns != newColumns);
  B.columns = newColumns;
  return changed;
}

// Reduce columns and rows, retaining row mapping, and return new clades needed
pair<map<int, int>, map<int, int>> reduceMatrix(Submatrix& B, vector<set<int>>& cladesFromReduce, set<int> cladeCols){
  bool reductionComplete = false;

  // Set up initial row and column mapping
  set<int> initialRows = B.rows;
  set<int> initialCols = B.columns;
  map<int, int> finalRowMapping;
  map<int, int> finalColMapping;
  for(int i : B.rows){ finalRowMapping[i] = i; }
  for(int i : B.columns){ finalColMapping[i] = i; }

  while(!reductionComplete){

    // Reduce rows, getting the row mapping
    map<int, int> currRowMapping = deleteDuplRows(B);

    // Then, update the final row mapping from previous
    for(int r : initialRows){
      finalRowMapping[r] = currRowMapping[finalRowMapping[r]];
    }    

    // Figure out the resulting groupings and add as clades
    for(int grouping : B.rows){
      set<int> newClade;
      for(int r : initialRows){
        if(finalRowMapping[r] == grouping){ newClade.insert(r); }
      }
      cladesFromReduce.push_back(newClade);
    }

    // Delete out trivial columns (including clade columns)
    bool wasTrivialColRemoved = deleteTrivialColumns(B);
    cladeCols = getSetIntersection(B.columns, cladeCols);

    // Delete out duplicate columns, retaining all the clade columns
    map<int, int> colMapping;
    bool wasDuplColRemoved = deleteDuplColumns(B, cladeCols, colMapping);
    cladeCols = getSetIntersection(B.columns, cladeCols);

    // Update column mapping
    for(int c : initialCols){
      finalColMapping[c] = colMapping[finalColMapping[c]];
    }   

    bool wasAnyColRemoved = wasTrivialColRemoved || wasDuplColRemoved;
    reductionComplete = !wasAnyColRemoved;
    reductionComplete = reductionComplete || (B.columns.size() == 0);

  }
  return {finalRowMapping, finalColMapping};
}

// Invert the columns of cols2Invert for all ROWS in initial
void invertColumns(Submatrix& initial, set<int> cols2Invert){
  for(int r : initial.rows){
    for(int c : cols2Invert){
      (*(initial.matrix))[r][c] = 1 - initial.getVal(r, c);
    }
  }
  return;
}

// Determine the block decomposition of a matrix
vector<Submatrix> splitSubmatrix(Submatrix& m){

  set<int> RRemain = m.rows;
  set<int> CRemain = m.columns;
  vector<Submatrix> components;

  // For each row, see if it's part of a split
  while(RRemain.size() > 0){

    set<int> rowComp = {*(RRemain.begin())};
    set<int> colComp = {};

    while(true){
      bool added = false;

      // Check columns to add
      for(int r : rowComp){
        for(int c : CRemain){
          if(m.getVal(r, c) == 1){
            colComp.insert(c);
            added = true;
          }
        }
      }
      CRemain = getSetDif(CRemain, colComp);

      // Check rows to add
      for(int c : colComp){
        for(int r : RRemain){
          if(m.getVal(r, c) == 1){
            rowComp.insert(r);
            added = true;
          }
        }
      }
      RRemain = getSetDif(RRemain, rowComp);

      if(!added){
        Submatrix bi;
        bi.rows = rowComp;
        bi.columns = colComp;
        bi.matrix = m.matrix;
        components.push_back(bi);
        break;
      }
    }
  }
  return components;
}

bool check2CladesCompatible(set<int>& clade1, set<int>& clade2){
  bool disjoint = (getSetIntersection(clade1, clade2).size() == 0);
  bool subSets = (getSetDif(clade1, clade2).size() == 0) || (getSetDif(clade2, clade1).size() == 0);
  return disjoint || subSets;
}

bool checkCladesCompatible(vector<set<int>>& clades){
  for(int i = 0; i < clades.size(); i++){
    for(int j = 0; j < i; j++){
      bool compatible = check2CladesCompatible(clades[i], clades[j]);
      if(!compatible){
        return false;
      }
    }
  }
  return true;
}

bool checkCladeAgaisntSet(vector<set<int>>& clades, set<int> newClade){
  for(int i = 0; i < clades.size(); i++){
    bool compatible = check2CladesCompatible(clades[i], newClade);
    if(!compatible){
      return false;
    }
  }
  return true;
}

struct treeLabelling {  

  // Given a clade / node (rooting the clade), what's it's parent?
  map<int, int> toParent; 

  // Given a clade / node (rooting the clade), what characters (1-indexed) were gained (+) or lost (-) on the edge leading to it?
  map<int, set<int>> edgeLabeling;  

  // Given a clade / node (rooting the clade), what rows / taxa map to that node (were in that clade, but not in its children)?
  map<int, set<int>> nodeLabeling;

};

int checkValidSolution(matrix& data, matrix& treeMatrix, treeLabelling& treeLabel, bool verbose = false){

  // Get dimensions
  int numRows = data.size();
  int numColumnsData = data[0].size();
  int numColumnsTreeMatrix = treeMatrix[0].size();

  // Get clades as sets of taxa
  if(verbose){ cout << "Reading tree's clades..." << endl; }
  vector<set<int>> treeClades;
  for(int c = 0; c < numColumnsTreeMatrix; c++){
    set<int> colSet;
    for(int r = 0; r < numRows; r++){ if(treeMatrix[r][c] == 1){ colSet.insert(r); } }

    // Don't add the clade if it's a duplicate of another clade
    bool duplicateClade = false;
    for(set<int>& otherClade : treeClades){
      if(otherClade == colSet){ duplicateClade = true; break; }
    }
    if(!duplicateClade){ treeClades.push_back(colSet); }
  } numColumnsTreeMatrix = treeClades.size();

  // Sort the clades by decreasing size!
  sort(treeClades.begin(), treeClades.end(), [](const set<int>& first, const set<int>& second) { return first.size() > second.size(); });
  
  // Check all clades are valid
  if(!checkCladesCompatible(treeClades)){ if(verbose){ cout << "Tree clades were invalid." << endl; } return -1; }
  if(verbose){ cout << "Tree clades are valid!" << endl; }

  // Now, flesh out the parent part of the treeLabel: for each clade, determine the clade least before it of which it's a subset (its parent)
  for(int i = 0; i < numColumnsTreeMatrix; i++){
    treeLabel.toParent[i] = -1;
    for(int j = i - 1; j >= 0; j--){
      if(getSetUnion(treeClades[i], treeClades[j]) == treeClades[j]){
        treeLabel.toParent[i] = j;
        break;
      }
    }
  }

  // Get clade differences if one clade is a subset of the other; record a parallel array of what clades made that (1 indexed!!)
  vector<set<int>> cladeDifs;
  vector<pair<int, int>> cladePairsforDif;
  if(verbose){ cout << "Calculating clade differences..." << endl; }
  for(int c1 = 0; c1 < numColumnsTreeMatrix; c1++){ for(int c2 = c1 + 1; c2 < numColumnsTreeMatrix; c2++){

    // If clade c2 is a subset of clade c1, add c1 - c2
    set<int> setUnion = getSetUnion(treeClades[c1], treeClades[c2]);
    if(setUnion == treeClades[c1]){
      cladeDifs.push_back(getSetDif(treeClades[c1], treeClades[c2]));
      cladePairsforDif.push_back({c1, c2});
    }

    // If clade c1 is a subset of clade c2, add c2 - c1
    if(setUnion == treeClades[c2]){
      cladeDifs.push_back(getSetDif(treeClades[c2], treeClades[c1]));
      cladePairsforDif.push_back({c2, c1});
    }
  }}
  if(verbose){ cout << "Clade differences calculated." << endl; }

  // Set up total false negative count
  int minFN = 0;

  // Get each column in the initial data as a set and record the number of FNs needed to make this set the perfect difference of two clades
  for(int c = 0; c < numColumnsData; c++){
    if(verbose){ cout << "Examining data column: " << c << endl; }

    set<int> colSet; for(int r = 0; r < numRows; r++){ if(data[r][c] == 1){ colSet.insert(r); } }
    int minFNColumn = numRows;

    // Record what the best clades were for this character
    int bestC1 = -1; int bestC2 = -1;

    // Check every treeClade
    for(int i = 0; i < cladeDifs.size(); i++){
      set<int>& cladeDif = cladeDifs[i]; int c1 = cladePairsforDif[i].first; int c2 = cladePairsforDif[i].second;

      // If the character's a subset of this cladDif, record the difference between the two and the clades that made this
      if(getSetUnion(colSet, cladeDif) == cladeDif){ 
        if((int)(cladeDif.size() - colSet.size()) < minFNColumn){
          minFNColumn = (int)(cladeDif.size() - colSet.size());
          bestC1 = c1; bestC2 = c2; 
        }
      }  
    }

    // If the column set is a subset of any clade, the FN count is cardinality of the clade minus the the column set size, minus 1 (each taxa is also a trivial clade)
    for(int c1 = 0; c1 < numColumnsTreeMatrix; c1++){ if(getSetUnion(colSet, treeClades[c1]) == treeClades[c1]){ 

      if((int)(treeClades[c1].size() - colSet.size()) < minFNColumn){
        minFNColumn = (int)(treeClades[c1].size() - colSet.size());
        if(minFNColumn == 1){
          minFNColumn = 0;
        }
        bestC1 = c1;
        bestC2 = -1;
      }
    }}

    // Record the best clades to explain the character
    if(bestC1 != -1){ treeLabel.edgeLabeling[bestC1].insert(c + 1); }
    if(bestC2 != -1){ treeLabel.edgeLabeling[bestC2].insert(-(c + 1)); }

    // Record the false negatives
    if(verbose){ cout << "Minimum FNs in data column " << c << " is: " << minFNColumn << endl; }
    minFN = minFN + minFNColumn;
  }

  // Lastly, set up the membership of rows to each node in the treeLabel: What's in this clade that's not in it's children?
  for(int i = 0; i < numColumnsTreeMatrix; i++){
    treeLabel.nodeLabeling[i] = treeClades[i];

    // Subtract from the parent
    int parent = treeLabel.toParent[i];
    if(parent == -1){ continue; }
    treeLabel.nodeLabeling[parent] = getSetDif(treeLabel.nodeLabeling[parent], treeClades[i]);
  }
  return minFN;
}


/*
---LINEAR CHAIN 1-Dollo FUNCTIONS ---
*/ 

struct RowDAG {      
  vector<vector<int>> taxaInNode;
  vector<set<int>> childrenOfNode;
  vector<set<int>> ancestorsOfnode;
  Submatrix* submatrix;               
};

RowDAG makeRowDAG(Submatrix& B, int b){

  vector<int> branchingPoint = (*B.matrix)[b];

  // Get matrix rows as sets
  map<int, set<int>> rowSets;
  for(int r : B.rows){
    set<int> rowSet;
    for(int c : B.columns){
      if(B.getVal(r, c) == 1){
        rowSet.insert(c);
      }
    }
    rowSets[r] = rowSet;
  }

  // Custom sorter for rows
  auto compareRows = [&](pair<vector<int>, int> row1, pair<vector<int>, int> row2){
    int sum1 = 0; int sum2 = 0;
    for(int c : B.columns){ 
      sum1 += row1.first[c]; sum2 += row2.first[c];
    }
    return (sum1 < sum2);
  };

  // Set up row dag
  RowDAG rd;
  rd.taxaInNode = {};
  rd.childrenOfNode = {};
  rd.ancestorsOfnode = {};
  rd.submatrix = &B;

  // Order submatrix rows by number of entries
  vector<pair<vector<int>, int>> orderedRows;
  for(int r : B.rows){
    orderedRows.push_back({(*B.matrix)[r], r});
  }
  sort(orderedRows.begin(), orderedRows.end(), compareRows);

  // For each row in ascending order
  for(int i = 0; i < orderedRows.size(); i++){

    int rowIndex = orderedRows[i].second;
    set<int> rowI = rowSets[rowIndex];

    // Check if a same row is present
    bool foundSameRow = false;
    int numDAGNodes = rd.taxaInNode.size();

    for(int j = numDAGNodes - 1; j >= 0; j--){

      set<int> rowJ = rowSets[rd.taxaInNode[j][0]];
      if(getSetDif(rowI, rowJ).size() == 0 && getSetDif(rowJ, rowI).size() == 0){

        // i and j are the same - add and move to next rowIndex
        rd.taxaInNode[j].push_back(rowIndex);
        foundSameRow = true;
        break;
      }
    }
    if(foundSameRow){
      continue;
    }

    // Add a new node in DAG
    rd.taxaInNode.push_back({rowIndex});
    rd.childrenOfNode.push_back({});
    rd.ancestorsOfnode.push_back({});
    int newNodeI = rd.taxaInNode.size() - 1;

    // For each row in rowDAG between root and i, in descending order
    for(int j = numDAGNodes - 1; j >= 0; j--){

      set<int> rowJ = rowSets[rd.taxaInNode[j][0]];
      int numInJNotI = (getSetDif(rowJ, rowI)).size();

      // If j isn't a subset of I, continue searching j
      if(numInJNotI != 0){ continue; }

      // J is a subset of I: first, check if j's in i's ancestor set already - if so, continue searching j
      set<int> curIAncestors = rd.ancestorsOfnode[newNodeI];
      if(curIAncestors.find(j) != curIAncestors.end()){ continue; }
      
      // Add i as child of j
      rd.childrenOfNode[j].insert(newNodeI);

      // Add j and the ancestorsOfnode J to I's ancestors
      rd.ancestorsOfnode[newNodeI] = getSetUnion(rd.ancestorsOfnode[newNodeI], rd.ancestorsOfnode[j]);
      rd.ancestorsOfnode[newNodeI].insert(j);
    }

  }
  return rd;
}

// Get a list of all maximal node paths from a specific node
vector<vector<int>> getAllMaximalNodePathsFromNode(RowDAG& rd, int nodeID){

  // Base Case: Return own node
  set<int> childrenOfCurNode = rd.childrenOfNode[nodeID];
  if(childrenOfCurNode.size() == 0){
    return {{nodeID}};
  }

  //For each child, get all paths, add self at front, return
  vector<vector<int>> results;
  for(int child : childrenOfCurNode){
    vector<vector<int>> childPaths = getAllMaximalNodePathsFromNode(rd, child);
    for(vector<int> childPath : childPaths){
      vector<int> newPath = childPath;
      newPath.insert(newPath.begin(), nodeID);
      results.push_back(newPath);
    }
  }
  return results;
}

// Get a list of all maximal node paths from the RowDAG
vector<vector<int>> getAllMaximalNodePaths(RowDAG& rd){

  int numNodes = rd.taxaInNode.size();
  vector<vector<int>> results;
  for(int i = 0; i < numNodes; i++){

    // If a source node, add results
    if(rd.ancestorsOfnode[i].size() == 0){
      vector<vector<int>> curNodePaths = getAllMaximalNodePathsFromNode(rd, i);
      for(vector<int> path : curNodePaths){
        results.push_back(path);
      }
    }
  }
  return results;
}

vector<int> getPQOrdering(set<int> onElements, vector<set<int>> constraints){
  
  // Set up tracker
  PQTree tracker = PQTree(onElements.size());

  // Get Row <--> PQ Index Mappings
  map<int, int> rowsToPQIndex;
  map<int, int> PQIndexToRows;
  int pqIndex = 1;
  for(int element : onElements){
    rowsToPQIndex[element] = pqIndex;
    PQIndexToRows[pqIndex] = element;
    pqIndex++;
  }

  for(set<int> constraint : constraints){
    if(!tracker.addConstraint(constraint, rowsToPQIndex)){
      return {-1};
    };
  }
  return tracker.getOrdering(PQIndexToRows);

}

// Check if a submatrix has the consec 1s subject to a branching point and a taxapath (grouping)
// If b = -1, just check if the consec 1s exists at ALL
vector<int> hasLinear(Submatrix& B, int b, vector<vector<int>> taxaPath){
  vector<set<int>> constraints;
  set<int> curConstraint;

  // Add in the DAG-based constraints
  if(b != -1){
    for(vector<int> taxaSet : taxaPath){
      curConstraint = getSetUnion(curConstraint, vec2Set(taxaSet));
      constraints.push_back(curConstraint);
      constraints.push_back(getSetDif(B.rows, curConstraint));
    }
  }

  // Add in the non-1, b[] = 0 C1p constraints (or just add all if b = -1)
  for(int c : B.columns){
    if(b == -1 || B.getVal(b, c) == 0){
      set<int> constraint = extractClade(B, c);
      constraints.push_back(constraint);
    }
  }

  // Force b to be at one end
  if(b != -1){
    constraints.push_back(getSetDif(B.rows, {b}));
  }

  // Return (b ends linear chain)
  vector<int> linearChain = getPQOrdering(B.rows, constraints);
  if(b != -1 && linearChain[0] == b){
    reverse(linearChain.begin(), linearChain.end());
  }
  return linearChain;
}

vector<pair<int, vector<int>>> getMaximalLinears(Submatrix& B, set<int> cladeCols){

  // 1) For each possible branching point, get the 1-rows DAG 
  // 1a) Remove clade violations and 1b) check for con-sec 1s
  vector<pair<int, vector<int>>> maxLinearsNo1GroupCheck;
  for(int b : B.rows){

    // Divide columns by 0 and 1 in b
    set<int> cols0;
    set<int> cols1;
    for(int c : B.columns){
      if(B.getVal(b, c) == 0){
        cols0.insert(c);
      }else{
        cols1.insert(c);
      }
    }
    Submatrix rd1b = extractSubmatrix(B, B.rows, cols1);
    RowDAG rd = makeRowDAG(rd1b, b);

    // Get all maximal node paths in the DAG - some nodes will have multiple taxa
    vector<vector<int>> allNodePaths = getAllMaximalNodePaths(rd);

    // For every node path
    for(vector<int> nodePath : allNodePaths){
      vector<vector<int>> taxaPath;
      set<int> allTaxaInPath;

      // Get the set of all taxa in each node and add it to the taxa path
      for(int node : nodePath){
        taxaPath.push_back(rd.taxaInNode[node]);
      }

      // 1a) Remove taxa in taxa path that violate existing clades
      for(vector<int>& taxaInNode : taxaPath){
        vector<int> toRemove;

        for(int taxa : taxaInNode){

          // Is the taxa in a clade b isn't?
          bool compatible = true;
          for(int c : cladeCols){
            if(B.getVal(taxa, c) == 1 && B.getVal(b, c) == 0){ compatible = false; }
          }

          // If so, mark it to be removed remove it
          if(!compatible){ toRemove.push_back(taxa); }
        }

        // Remove everything marked
        for(int taxaToRemove : toRemove){
          taxaInNode.erase(remove(taxaInNode.begin(), taxaInNode.end(), taxaToRemove), taxaInNode.end());
        }
        allTaxaInPath = getSetUnion(allTaxaInPath, vec2Set(taxaInNode));
      }

      // 1b) Does this taxa path actually have C1P on other columns?
      Submatrix matrixOnPotentialLinear = extractSubmatrix(B, allTaxaInPath, B.columns);
      vector<int> maximalChain = hasLinear(matrixOnPotentialLinear, b, taxaPath);    

      // Add path
      if(maximalChain[0] != -1){
        maxLinearsNo1GroupCheck.push_back({b, maximalChain});
      }
    }
  }

  // Check every chain path prefix
  vector<pair<int, vector<int>>> maxLinears;
  for(pair<int, vector<int>>& maxLinear : maxLinearsNo1GroupCheck){

    int b = maxLinear.first;
    vector<int>& linearChain = maxLinear.second;

    // Else, it passes - add every prefix
    vector<int> curPrefix;
    
    for(int i : linearChain){

      // Set up prefix
      curPrefix.push_back(i);

      // 2) We've removed clade violations and checked consec 1s - 
      // Lastly, check if every chain path, for every 0 column, contain all taxa with a 1 if a 1 exists
      set<int> cols0;
      set<int> cols1;
      for(int c : B.columns){
        if(B.getVal(i, c) == 0){
          cols0.insert(c);
        }else{
          cols1.insert(c);
        }
      }
      bool onesGrouped = true;

      for(int c : cols0){
        set<int> taxaWith1 = extractClade(B, c);
        if(getSetIntersection(taxaWith1, vec2Set(curPrefix)).size() > 0 && getSetDif(taxaWith1, vec2Set(curPrefix)).size() > 0){
          onesGrouped = false;
        }
      }
      
      if(!onesGrouped){ continue; }
      maxLinears.push_back({i, curPrefix});
    }
  }

  // Remove duplicates
  vector<pair<int, vector<int>>> maxLinearsMinimal;
  for(int i = 0; i < maxLinears.size(); i++){
    pair<int, vector<int>>& maxLinearI = maxLinears[i];
    bool isDuplicate = false;
    for(int j = 0; j < i; j++){  
      pair<int, vector<int>>& maxLinearJ = maxLinears[j]; 
      if(maxLinearI.first == maxLinearJ.first && maxLinearI.second == maxLinearJ.second){
        isDuplicate = true; break;
      }
    }
    if(!isDuplicate){ maxLinearsMinimal.push_back(maxLinearI);}
  }

  // Sort chains by size before returning
  auto compareChains = [](pair<int, vector<int>> chain1, pair<int, vector<int>> chain2){
    return (chain1.second.size() > chain2.second.size());
  };
  sort(maxLinearsMinimal.begin(), maxLinearsMinimal.end(), compareChains);
  return maxLinearsMinimal;
}

typedef vector<pair<string, vector<int>>> tree;

void addClade(set<int> clade, tree& T){
  T.push_back({"Clade        \t", set2Vec(clade)});
}

void addTree(tree& nT, tree& T){
  for(pair<string, vector<int>> m : nT){
    T.push_back(m);
  }
}

// Given a tree and a old -> new row mapping, undo: 
// Replace each taxa with everything that maps to it
// numTaxaOrig is the number of things in the map
tree expandTree(tree& T, set<int> initialTaxa, map<int, int>& rowMapping){
  tree newTree;
  map<int, vector<int>> reverseMap;
  set<int> reverseMapKeys;
  for(int i : initialTaxa){
    reverseMap[rowMapping[i]] = {};
    reverseMapKeys.insert(rowMapping[i]);
  }
  for(int i : initialTaxa){
    reverseMap[rowMapping[i]].push_back(i);
  }

  // For each clade in the tree, substitute
  for(pair<string, vector<int>> m : T){
    set<int> newClade;
    for(int r : m.second){
      vector<int> rCorrespond = reverseMap[r];
      for(int r1 : rCorrespond){
        newClade.insert(r1);
      }
    }
    addClade(newClade, newTree);
  }
  return newTree;
}

tree invalidTree(){
  tree T;
  T.push_back({"Invalid      \t", {-1}});
  return T;
}

bool isTreeInvalid(tree& T){
  return (T[0].second)[0] == -1;
}

// Assumes tree has ONLY clades (output listed in order of Tree taxa in first clade)
pair<set<int>, matrix> treeCladesToMatrix(tree& T){
  if(isTreeInvalid(T)){
    cout << "Invalid Tree" << endl;
    return {{-1}, {{-1}}};
  }
  matrix m;
  set<int> taxaSet;
  set<int> rowIndices;
  map<int, int> rowMapping;
  for(int i = 0; i < (T[0].second).size(); i++){
    m.push_back({});
    taxaSet.insert(T[0].second[i]);
    rowIndices.insert(T[0].second[i]);
    rowMapping[T[0].second[i]] = i;
  }
  for(pair<string, vector<int>> curClade : T){
    set<int> cladeSet = vec2Set(curClade.second);
    set<int> notInClade = getSetDif(taxaSet, cladeSet);
    for(int i : cladeSet){
      m[rowMapping[i]].push_back(1);
    }
    for(int i : notInClade){
      m[rowMapping[i]].push_back(0);
    }
  }
  return {rowIndices, m};
}
