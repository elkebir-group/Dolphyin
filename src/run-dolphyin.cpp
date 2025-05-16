#include <algorithm>
#include <fstream>
#include <sstream>
#include <random>

#include "utils.cpp"

// Given input matrix B and cladeset, return a list of the clades
tree solve1Dollo(Submatrix& B, set<int>& cladeCols, int numBranch, double p, double k, int seed){

  bool verbose = false;

  // Set up solution
  tree T;
  addClade(B.rows, T);
  set<int> initialTaxa = B.rows;
  set<int> initialCols = B.columns;

  // Reduce the matrix: delete additional rows/columns and trivial columns
  vector<set<int>> cladesFromReduce;
  pair<map<int, int>, map<int, int>> mappings = reduceMatrix(B, cladesFromReduce, cladeCols);
  map<int, int> rowMapping = mappings.first; map<int, int> colMapping = mappings.second;

  cladeCols = getSetIntersection(cladeCols, B.columns);

  // If empty, add new clades to T from the row consolidation and return
  if(B.rows.size() == 0 || B.columns.size() == 0){
    for(set<int> cladeFromReduced : cladesFromReduce){ addClade(cladeFromReduced, T); }
    return T;
  }

  // Perform false negative correction if needed: Equalize non-clade column 1s between p% pairs of rows (selected with replacement) if the character difference is less than k%
  // Retain the switched positions to undo them as required (failure to find a tree, i.e.)
  int numRows = B.rows.size(); int numPairs = round((double)(numRows * numRows) * p); vector<pair<int, int>> fnLocations;
  if(numPairs > 0){

    // Set up a random generator with the current seed
    mt19937 randGen(seed);

    // Create a distribution over the possible 0 to (numRows - 1) indices for row selection, based on the total row mapping counts
    vector<double> indexWeights = {};
    vector<double> indexEqual = {};
    for(int i = 0; i < numRows; i++){
      int ithRow = *next(B.rows.begin(), i);
      double weight = 0;
      for(int r1Mapped : initialTaxa){
        if(rowMapping[r1Mapped] == ithRow){
          weight = weight + 1;
        }
      }
      indexWeights.push_back(1.0 / weight);
      indexEqual.push_back(1.0);
    }
    discrete_distribution<> distribution(indexWeights.begin(), indexWeights.end());
    discrete_distribution<> distrEqual(indexEqual.begin(), indexEqual.end());

    for(int i = 0; i < numPairs; i++){
      int randIndex1 = distribution(randGen);
      int randIndex2 = distrEqual(randGen);
      int r1 = *next(B.rows.begin(), randIndex1); int r2 = *next(B.rows.begin(), randIndex2);

      // Skip row pair if the difference between them is less than k%
      int numDifColumns = 0;
      for(int c : B.columns){ if(B.getVal(r1, c) != B.getVal(r2, c)){ 

        // How many columns mapped to c (how important is c?)
        for(int cOld : initialCols){
          if(colMapping[cOld] == c){ numDifColumns++; }
        }
      }}
      double dif = (double)(numDifColumns) / initialCols.size();
      if(dif > k){ continue; }

      // Switch out all false negatives where a 0 is seen only in 1 row, and track the number of 1s that this corresponds to prior to matrix reduction
      for(int c : getSetDif(B.columns, cladeCols)){
        if(B.getVal(r1, c) == 0 && B.getVal(r2, c) == 1){
          (*(B.matrix))[r1][c] = 1; fnLocations.push_back({r1, c});
          for(int r1Mapped : initialTaxa){ if(rowMapping[r1Mapped] == r1){ (*(B.matrix))[r1Mapped][c] = 1; }}
        }
        if(B.getVal(r2, c) == 0 && B.getVal(r1, c) == 1){
          (*(B.matrix))[r2][c] = 1; fnLocations.push_back({r2, c});
          for(int r2Mapped : initialTaxa){ if(rowMapping[r2Mapped] == r2){ (*(B.matrix))[r2Mapped][c] = 1; }}
        }
      }
    }
  }

  // Check the base case: if linear AND respects clades, simply find the ordering and return
  vector<int> fullLinearOrdering = hasLinear(B, -1, {{}});
  if(fullLinearOrdering[0] != -1){

    // Check if it respects clades
    bool respectsClades = true;
    for(int c : cladeCols){
      int lastEntry = B.getVal(fullLinearOrdering[fullLinearOrdering.size() - 1], c);
      if(lastEntry == 0){
        respectsClades = false;
      }
    }

    // Set up tree and add linear chain as clades
    if(respectsClades){

      // Add linear chain clades
      set<int> newLinClade = B.rows;
      addClade(newLinClade, T);

      int prevInOrder = -1;
      for(int i : fullLinearOrdering){
        newLinClade.erase(i);
        addClade(newLinClade, T);

       prevInOrder = i;
      }

      // Expand tree, add new clades to T from the row consolidation
      T = expandTree(T, initialTaxa, rowMapping);
      for(set<int> cladeFromReduce : cladesFromReduce){ addClade(cladeFromReduce, T); }
      return T;
    }
  }

  // If the number of branching is 0, we needed the base else - now, we're done
  if(numBranch == 0){
    for(pair<int, int> fnPair : fnLocations){ 
      int r = fnPair.first; int c = fnPair.second;
      for(int r1Mapped : initialTaxa){ if(rowMapping[r1Mapped] == r){ (*(B.matrix))[r1Mapped][c] = 0; }}
    } 
    return invalidTree();
  }

  // Get list of possible linear chains
  vector<pair<int, vector<int>>> maxLinears = getMaximalLinears(B, cladeCols);

  // Add the empty chain
  maxLinears.push_back({-1, {}});

  // For each linear chain
  for(pair<int, vector<int>>& maxLinearPair : maxLinears){

    int b = maxLinearPair.first;
    vector<int> maxLinear = maxLinearPair.second;

    // Get taxa in and outside linear chain 
    set<int> X0 = vec2Set(maxLinear);
    set<int> Xn = getSetDif(B.rows, X0);

    // Check clade compatibility of chain
    bool chainSuccessful = true;
    if(b != -1){
      for(int c : cladeCols){

        // If b -> 0, everything in the chain must be 0
        if(B.getVal(b, c) == 0){
          for(int r : X0){ if(B.getVal(r, c) == 1){ chainSuccessful = false; break; } }
        }

        // If b -> 1, everything in outside the chain must be 1
        if(B.getVal(b, c) == 1){
          for(int r : Xn){ if(B.getVal(r, c) == 0){ chainSuccessful = false; break; } }
        }
      }
    }
    if(!chainSuccessful){  continue; }

    // Add linear chain as clades
    tree T;
    set<int> newLinClade = B.rows;
    addClade(newLinClade, T);
    for(int i : maxLinear){
      newLinClade.erase(i);
      addClade(newLinClade, T);
    }

    // If there's no such taxa, we're done!!!
    if(Xn.size() == 0){ 

      // Expand tree, add new clades to T from the row consolidation
      T = expandTree(T, initialTaxa, rowMapping);
      for(set<int> cladeFromReduce : cladesFromReduce){ addClade(cladeFromReduce, T); }
      return T; 
    }

    // Get submatrix on these remaining taxa
    Submatrix Bt = extractSubmatrix(B, Xn, B.columns);

    // Get all columns gained in linear chain and lost below: X0+Xt- <- all c in Xt- where B[b,c] is 1
    set<int> X0pXtm;
    if(b != -1){
      for(int c : Bt.columns){
        if(B.getVal(b, c) == 1){ X0pXtm.insert(c); }
      }
    }

    // Invert all columns gained in X0 and lost after in matrix Bt
    invertColumns(Bt, X0pXtm);

    // What columns are known previous-losses?
    set<int> prevLosses = getSetUnion(X0pXtm, cladeCols);

    // Any other columns that could be gained right before the split?
    set<int> extraColumnsPotentLost = getSetDif(Bt.columns, prevLosses);
    set<int> extraColumnsLost;
    set<int> extraColumnsUnknown;

    // Get the subtree groupings on the known, already-gained columns X0pXtm
    Submatrix BtPrevLosses = extractSubmatrix(B, Bt.rows, prevLosses);
    vector<Submatrix> knownSplits = splitSubmatrix(BtPrevLosses);
    int numKnownSplits = knownSplits.size();

    // If there's not even more than 1 split on the known columns, quit
    if(numKnownSplits <= 1){
      invertColumns(Bt, X0pXtm); continue;
    }

    // For every other column: does flipping it make it fit nicely into one of these groups?
    for(int c : extraColumnsPotentLost){
      set<int> extraColClade = extractClade(Bt, c);
      set<int> extraInverse = getSetDif(Bt.rows, extraColClade);

      bool fitsInto = false;

      for(int i = 0; i < numKnownSplits; i++){

        // The column fits into one category!
        if(getSetDif(extraColClade, knownSplits[i].rows).size() == 0){
          fitsInto = true;
          break;
        }

        // Inverting the column fits into one category!
        if(getSetDif(extraInverse, knownSplits[i].rows).size() == 0){
          extraColumnsLost.insert(c);
          fitsInto = true;
          break;
        }
      }

      // If itself or it's inverse didn't fit, we need to consider it for all options
      if(!fitsInto){
        extraColumnsUnknown.insert(c);
      }
    }

    // Make the flips on invertable columns, add them into the list!
    invertColumns(Bt, extraColumnsLost);
    X0pXtm = getSetUnion(X0pXtm, extraColumnsLost);

    // Time to check all of the unknown columns, set up for split checking
    vector<Submatrix> splits;
    int t;

    // Try flipping every subset of them and check for splits
    int numPotSplitCols = extraColumnsUnknown.size();
    expIterator expIter = getIterator(numPotSplitCols);
    bool possibleSplit = false;
    while(!iteratorDone(expIter)){

      bool newGainsValid = true;

      // Get the current iterator state
      set<int> colTryFlip = getIteratorState(expIter, extraColumnsUnknown);

      // Flip columns and check
      invertColumns(Bt, colTryFlip);

      // Find the component / subtrees in extracted matrix Bt 
      splits = splitSubmatrix(Bt);
      t = splits.size();

      // There should actually be a branching point
      if(t > 1){ 

        // Add in these new columns
        set<int> possibleFullX0pXtm = getSetUnion(X0pXtm, colTryFlip);
        possibleSplit = true; 

        // Remove pointless columns
        deleteTrivialColumns(Bt);

        // Check clade compatibility of splits themselves
        vector<set<int>> clades = getClades(B, cladeCols);
        for(int i = 0; i < t; i++){
          set<int> splitClade = splits[i].rows;
          if(!checkCladeAgaisntSet(clades, splitClade)){
            newGainsValid = false; 
            break;         
          }
          addClade(splitClade, T);
        }

        if(!newGainsValid){ 
          invertColumns(Bt, colTryFlip);
          incrIterator(expIter);
          continue; 
        }

        // Check clade compatibility of previous with splits's gained-linear later-lost columns
        vector<set<int>> cladeColsT;
        for(int i = 0; i < t; i++){
          
          // Get all columns gained in linear and lost in specific subtree: X0+Xi- <- X0+Xt- intersection Ci 
          Submatrix Bi = splits[i];
          set<int> X0pXim = getSetIntersection(Bi.columns, possibleFullX0pXtm);

          // Save these clade columns and add to full list
          cladeColsT.push_back(X0pXim);
          vector<set<int>> cladesI = getClades(Bi, X0pXim);
          for(int j = 0; j < cladesI.size(); j++){
            clades.push_back(cladesI[j]);
          }
        }
      
        if(!checkCladesCompatible(clades)){ 
          invertColumns(Bt, colTryFlip);
          incrIterator(expIter);
          continue; 
        }

        // For each subtree found, recurse
        for(int i = 0; i < t; i++){

          // Get split submatrix, get all columns gained in linear and lost in specific subtree
          Submatrix Bi = splits[i];
          set<int> cladeColsI = cladeColsT[i];

          // Add existing clades to input variables
          set<int> fullCladeCols = getSetUnion(cladeColsI, cladeCols);
          set<int> fullCols = getSetUnion(splits[i].columns, cladeCols);
          Bi.columns = fullCols;

          // Recurse on Bt[Ri, Ci - X0+Xi-], S'
          tree ti = solve1Dollo(Bi, fullCladeCols, numBranch - 1, p, k, seed + 1);
          if(isTreeInvalid(ti)){ newGainsValid = false; break; }

          // Add Ti to T
          addTree(ti, T);
        }

        if(!newGainsValid){ 
          invertColumns(Bt, colTryFlip);
          incrIterator(expIter);
          continue; 
        }

        // Expand tree, add new clades to T from the row consolidation
        T = expandTree(T, initialTaxa, rowMapping);
        for(set<int> cladeFromReduced : cladesFromReduce){ addClade(cladeFromReduced, T); }
        return T;

      }

      // Try the next one
      invertColumns(Bt, colTryFlip);
      incrIterator(expIter);

    }

    invertColumns(Bt, X0pXtm);
   
  }

  // We failed to find a tree, undo false negatives and return
  for(pair<int, int> fnPair : fnLocations){
    int r = fnPair.first; int c = fnPair.second;
    for(int rMapped : initialTaxa){ if(rowMapping[rMapped] == r){ (*(B.matrix))[rMapped][c] = 0; }}
  }
  return invalidTree();
}

int main(int argc, char *argv[]) {

  // Read in input, output locations, and parameters
  if(argc != 6){
    cout << "Usage: " << argv[0] << " data-matrix output-treeClades p k seed" << endl;
    return 0;
  }

  // Read in input data matrix, scrapping the two-line header
  matrix m = getMatrixFromFile(argv[1], 2);
  pair<set<int>, set<int>> rc = getRowCols(m);
  Submatrix B; B.rows = rc.first; B.columns = rc.second; B.matrix = &m;

  // Get the p, k, seed parameters (and numBranch, which is always -1)
  double p = stod(argv[3]); double k = stod(argv[4]); int seed = stoi(argv[5]);
  int numBranch = -1;

  // Solve the 1-Dollo instance
  set<int> cladeCols = {};
  tree T = solve1Dollo(B, cladeCols, numBranch, p, k, seed);

  // Write out the resulting clades to the tree output file
  pair<set<int>, matrix> mSolPair = treeCladesToMatrix(T);
  matrix mSol = mSolPair.second;
  string outMatrix = argv[2];
  writeMatrixToFile(outMatrix, mSol);

}