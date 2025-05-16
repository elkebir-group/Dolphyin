#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <queue>

using namespace std;

class PQTree {   
   
  // Manage a single nodes, under which more nodes may be enclosed
  class CSCTNode{

    public:
    
      // Make a leaf node, carrying a value for grouping (1 - n)
      CSCTNode(int number){
        leafNumber = number;
        nodeType = 'C';
      };

      // Make a general internal node, of either a (C)ircle or (S)qaure type
      CSCTNode(char type){
        leafNumber = -1;
        nodeType = type;
      };

      // Add a child node to a parent (Added as pointer)
      void addChild(CSCTNode* child){ children.push_back(child); }

      // Retrieve general information on a node
      bool isCircle(){ return nodeType == 'C'; }
      bool isLeaf(){ return leafNumber >= 1; }
      char getLabel(){ return label; }
      int getNumChildren(){ return children.size(); }

      // Check if a ode only has one child and return it if so, else NULL
      CSCTNode* getOnlyChild(){
        if(getNumChildren() == 1){ return children[0]; }
        return NULL;
      }

      // Return a possible ordering
      vector<int> getOrdering(){
        if(isLeaf()){ return {leafNumber}; }

        vector<int> result;
        for(CSCTNode* child : children){ 
          vector<int> childResult = child->getOrdering();
          for(int leafVal : childResult){
            result.push_back(leafVal);
          }
        }
        return result;
      }

      // Return a string representation of a node by a [S] and (C) format
      string asString(){
        if(isLeaf()){ return to_string(leafNumber); }

        string ans;
        string startDelim = "["; string endDelim = "]"; 
        if(isCircle()){ startDelim = "("; endDelim = ")"; }
        ans += startDelim;

        for(CSCTNode* child : children){ ans += child->asString() + ", "; }
        ans.pop_back(); ans.pop_back();

        return ans + endDelim;
      }

      // Return a string representation of a node by a [S] and (C) format, 
      // with the node labelling by a constraint appended to the front of each node
      string asStringLabelled(){
        if(isLeaf()){ 
          string ans; ans += getLabel();
          return ans + to_string(leafNumber); 
        }

        string ans;
        string startDelim = "["; string endDelim = "]"; 
        if(isCircle()){ startDelim = "("; endDelim = ")"; }
        ans += getLabel(); ans += startDelim;

        for(CSCTNode* child : children){ ans += child->asStringLabelled() + ", "; }
        ans.pop_back(); ans.pop_back();

        return ans + endDelim;
      }

      // Reverse the order of every terminal (1 - n) in a node, useful for debugging
      void reverseFullOrder(){
        if(isLeaf()){ return; }
        for(CSCTNode* child : children){ child->reverseFullOrder(); }
        reverse(children.begin(),children.end());
        if(label == 'L'){ label = 'R'; } else if(label == 'R'){ label = 'L'; }
      }

      // Reverse the order of only said node's children, useful for debugging
      void reverseOwnOrder(){ reverse(children.begin(),children.end()); }

      // Set all constraint-specific labels to 0 in all nodes enclosed under this node
      void clearLabels(){
        label = 0;
        for(CSCTNode* child : children){ child->clearLabels(); }
      }

      // Eliminate unncessary nodes and structural irritations, does not affect constraints
      void cleanNode(){
        if(isLeaf()){ return; }

        int numChildren = getNumChildren();
        for(int i = 0; i < numChildren; i++){
          
          // For each child, erase it entirely from the array if it contains no nodes
          if(!(children[i]->isLeaf()) && children[i]->getNumChildren() == 0){
            children.erase(children.begin() + i);
            i--; numChildren--;
            continue;
          }

          // For each child, replace it with grandchild if it only has one child (redundant child)
          // Delete child while maintaining spot, HAKER's still here
          CSCTNode* onlyGrandChild = children[i]->getOnlyChild();
          while(onlyGrandChild){

            // TODO: Delete child, manage memory
            children[i] = onlyGrandChild;
            onlyGrandChild = children[i]->getOnlyChild();
          }
        }  

        // Recurse on children
        for(CSCTNode* child : children){
          child->cleanNode();
        }
      }

      // Label all nodes below tree with proper character label based on include status
      // (A)ll, (L)eft, (N)one, (M)iddle, (I)nvalid
      void labelWithConstraint(map<int, bool>& inConstraint){

        // If leaf, mark based on if contained in constraint (Base case)
        if(isLeaf()){
          if(inConstraint[leafNumber]){ label = 'A'; } else { label = 'N'; }
          return;
        }

        // Recurse, labelling children
        int numChildren = 0;
        for(CSCTNode* child : children){
          child->labelWithConstraint(inConstraint);
          numChildren++;
        }

        // Count stats of labelled children for current determination
        string childrenLabels;
        for(int i = 0; i < numChildren; i++){
          childrenLabels += children[i]->getLabel();
        };
        int numAlls = count(childrenLabels.begin(), childrenLabels.end(), 'A');
        int numLefts = count(childrenLabels.begin(), childrenLabels.end(), 'L');
        int numMiddle = count(childrenLabels.begin(), childrenLabels.end(), 'M');
        bool seenInvalid = count(childrenLabels.begin(), childrenLabels.end(), 'I') > 0;
        int numNone = count(childrenLabels.begin(), childrenLabels.end(), 'N');


        // Determine general node label for obvious cases
        if(seenInvalid || numLefts > 2 || (numMiddle > 0 && numChildren - numNone > 1)){
          label = 'I'; return; }
        if(numNone == numChildren){ 
          label = 'N'; return; }
        if(numAlls == numChildren){
          label = 'A'; return; }
        if(numMiddle > 0){
          label = 'M'; return; }

        // Determine more difficult, node type specific cases
        if(isCircle()){
          if(numLefts == 2){ label = 'M'; return; }  
          if(numLefts < 2){ label = 'L'; return; }
        } else {
          
          // Scan for sqaure L through the string forwards (A*(L+e)N*)
          char prevLetter = 'A';
          bool leftLabel = true;
          for(int i = 0; i < numChildren; i++){
            char currChar = childrenLabels[i];
            if((prevLetter == 'L' && currChar == 'L') || 
              (prevLetter == 'L' && currChar == 'A') ||
              (prevLetter == 'N' && currChar == 'A') || 
              (prevLetter == 'N' && currChar == 'L')
            ){ leftLabel = false; break; }
            prevLetter = currChar;
          }
          if(leftLabel) { label = 'L'; return; }

          // Scan for sqaure L through the string backwards (N*(L+e)A*)
          prevLetter = 'A';
          leftLabel = true;
          for(int i = numChildren; i >= 0; i--){
            char currChar = childrenLabels[i];
            if((prevLetter == 'L' && currChar == 'L') || 
              (prevLetter == 'L' && currChar == 'A') ||
              (prevLetter == 'N' && currChar == 'A') || 
              (prevLetter == 'N' && currChar == 'L')
            ){ leftLabel = false; break; }
            prevLetter = currChar;
          }
          if(leftLabel) { label = 'L'; return; }

          // Scan for sqaure M through the string (N*(L+e)A*(L+e)N*)
          int startLeftIndex = 0;
          int endLeftIndex = numChildren - 1;
          char currChar;
          while(startLeftIndex < numChildren){
            char currChar = childrenLabels[startLeftIndex];
            if(currChar != 'N'){
              if(currChar == 'L'){ startLeftIndex++; }
              break; 
            }
            startLeftIndex++;
          }
          while(endLeftIndex >= 0){
            char currChar = childrenLabels[endLeftIndex];
            if(currChar != 'N'){
              if(currChar == 'L'){ endLeftIndex--; }
              break; 
            }
            endLeftIndex--;
          }

          bool middleAs = true;
          for(int i = startLeftIndex; i < endLeftIndex + 1; i++){
            char currChar = childrenLabels[i];
            if(currChar != 'A'){ middleAs = false; break; }
          }

          // Make our assesment, defaulting to I otherwise
          if(middleAs){ label = 'M'; return; }      
          label = 'I'; return;     
        }
      }

      // Modify the structure to force the constraint on the tree
      void implementConstraint(){

        // Search for the responsible parent (2 nodes with non-N labels)
        int numChildren = getNumChildren();
        int recurseIndex = -1;
        for(int i = 0; i < numChildren; i++){
          if(children[i]->getLabel() != 'N'){
            if(recurseIndex >= 0){
              recurseIndex = -2;
              break;
            }
            recurseIndex = i;
          }
        }
        if(recurseIndex >= 0){
          children[recurseIndex]->implementConstraint();
          return;
        }

        // If responsible parent is Circle, build the 3-node C-S-C structure
        if(isCircle()){
          CSCTNode* replaceNode = new CSCTNode('C');
          CSCTNode* newNode = new CSCTNode('S');
          CSCTNode* newNodeSub = new CSCTNode('C');
          int numLefts = 0;

          // Add nones to first, A to last, and handle Ls by splitting to middle
          for(CSCTNode* child : children){
            if(child->getLabel() == 'N'){
              replaceNode->addChild(child);
            }
            if(child->getLabel() == 'A'){
              newNodeSub->addChild(child);
            }
            if(child->getLabel() == 'L'){
              if(numLefts == 1){ newNode->addChild(newNodeSub);}
              child->splitOnConstraint(newNode, numLefts < 1);
              numLefts++;
            }
          }   

          // If two lefts weren't seen, make sure to add the last node anyways
          if(numLefts < 2){ newNode->addChild(newNodeSub); }

          // Connect our nodes and move on
          replaceNode->addChild(newNode);
          children = {replaceNode};
        }   

        // If responsible parent is Square, built the single S node
        // the parent's children must be of order (N*(L+e)A*(L+e)N*)
        else {
          CSCTNode* replaceNode = new CSCTNode('S');
          bool expectingLLeftofAs = true;

          // Add nones and As in order and handle Ls by splitting
          for(CSCTNode* child : children){
            if(child->getLabel() == 'N'){
              replaceNode->addChild(child);
            }
            if(child->getLabel() == 'L'){
              child->splitOnConstraint(replaceNode, expectingLLeftofAs);
              expectingLLeftofAs = false;
            }
            if(child->getLabel() == 'A'){
              replaceNode->addChild(child);
              expectingLLeftofAs = false;
            }
          }
          children = {replaceNode};
        }     
      }

      // After we identify the responsible parent, if we see an L,
      // Only some of the nodes must be constrained
      // We recurse of the L, splitting nodes and adding As to our newNode constrainer
      // nonesToLeft, do we add Ns to the left side of newNode or not? Order matters, newNode is S
      void splitOnConstraint(CSCTNode* newNode, bool nonesToLeft){
      
        
        // If the split is on a circle, partition the As and Ns into separate Circle nodes and add
        if(isCircle()){
          CSCTNode* noneNode = new CSCTNode('C');
          CSCTNode* allNode = new CSCTNode('C');
          CSCTNode* childLabelledLeft = NULL;

          // Add nones and As to right list, and record an L to recurse on (guranteed only 0 - 1)
          for(CSCTNode* child : children){
            if(child->getLabel() == 'N'){ noneNode->addChild(child); }
            if(child->getLabel() == 'L'){ childLabelledLeft = child; }
            if(child->getLabel() == 'A'){ allNode->addChild(child); }
          }
          if(nonesToLeft){
            newNode->addChild(noneNode);
            if(childLabelledLeft){ childLabelledLeft->splitOnConstraint(newNode, nonesToLeft); }
            newNode->addChild(allNode);
          }else{
            newNode->addChild(allNode);
            if(childLabelledLeft){ childLabelledLeft->splitOnConstraint(newNode, nonesToLeft); }
            newNode->addChild(noneNode);
          }
        
        // If the split is on a square, childrens' order is either (N*(L+e)A*) or (A*(L+e)N*)
        // Add in the order needed to accurately maintain the requested additions to newNode
        }else{
          bool childNonesLeft = children[0]->getLabel() == 'N' || children[children.size() - 1]->getLabel() == 'A';
          bool iterateForward = (childNonesLeft == nonesToLeft);
          int numChildren = getNumChildren();

          // Add nones and As to right list, and recurse on an L instantly to maintain order (guranteed only 0 - 1)
          if(iterateForward){
            for(int i = 0; i < numChildren; i++){
              CSCTNode* child = children[i];
              if(child->getLabel() == 'L'){ child->splitOnConstraint(newNode, nonesToLeft); } 
              else{ newNode->addChild(child); }
            }
          }else{
            for(int i = numChildren - 1; i >= 0; i--){
              CSCTNode* child = children[i];
              if(child->getLabel() == 'L'){ child->splitOnConstraint(newNode, nonesToLeft); } 
              else{ newNode->addChild(child); }
            }
          }
        }
      }      

    private:

      // Node fields
      vector<CSCTNode*> children;
      int leafNumber;
      char nodeType;
      char label;

  };

  public:      

    // Instantiate a tree with n nodes (1 - n) and no order constraints (all under C)
    PQTree(int n){
      root = new CSCTNode('C');
      for(int i = 1; i <= n; i++){
        root->addChild(new CSCTNode(i));
      }
    };

    // Return a string representation of the tree by a [S] and (C) format
    string asString(){
      return root->asString();
    }

    // Return a string representation of the tree by a [S] and (C) format
    // with the node labelling by a constraint appended to the front of each tree node
    string asStringLabelled(){
      return root->asStringLabelled();
    }

    // Force a constraint onto the tree, returning false if impossible
    bool forceConstraint(vector<int> constraint){

      // If constraint has less than 2 elements, return
      // It's unnecessary, and will make it impossible to find a responsible parent 
      if(constraint.size() < 2){ return true; }

      // Label our nodes to ID what nodes to include
      bool possible = labelWithConstraint(constraint);
      if(!possible){ return false; }

      // Add the constraint, and clean the tree afterwards to resolve redundancies
      root->implementConstraint();
      cleanTree();
      return true;
    }

    // Clear all labels, useful after a constraint is attempted
    void clearLabels(){
      root->clearLabels();
    }

    // Remove redundancies from the tree structure 
    void cleanTree(){

      //Delete chains in the root, resetting if needed
      CSCTNode* onlyChild = root->getOnlyChild();
      while(onlyChild){

        // Delete child
        root = onlyChild;
        onlyChild = root->getOnlyChild();
      }

      // Recurse on the rest of the tree
      root->cleanNode();
    }

    // Label a tree's nodes based on the requested constraint, populating an int map, 
    // returning if constraint is possible
    bool labelWithConstraint(vector<int>& constraint){
      map<int, bool> inConstraint;
      for(int each : constraint){ inConstraint[each] = true; }
      root->labelWithConstraint(inConstraint);
      return (root->getLabel() != 'I');
    }


    // INTERFACE FOR DOLLO

    // Return a possible ordering
    vector<int> getOrdering(map<int, int>& PQIndexToRows){
      vector<int> origOrdering = root->getOrdering();
      vector<int> mappedOrdering;
      for(int i : origOrdering){
        mappedOrdering.push_back(PQIndexToRows[i]);
      }
      return mappedOrdering;
    }

    // Add a constraint onto the tree from a 0-indexed set, returning false if impossible
    bool addConstraint(set<int> constraint, map<int, int>& rowsToPQIndex){
      vector<int> vectorCons;
      for(int i : constraint){
        vectorCons.push_back(rowsToPQIndex[i]);
      }
      return forceConstraint(vectorCons);
    }

  private:

    // Tree fields
    CSCTNode* root;

};
