#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
#include<unordered_set>
#include<queue>
#include <algorithm>
#include <any>
using namespace std;


struct node
{
    string name;
    node *next;
    node() : next(NULL){};
    node(string name) : name(name), next(NULL){};
};

struct list{
    node *header;
    node *newNode;
    node *tail;
    int length;
    list() : header(NULL), newNode(NULL), tail(NULL){};
    
    void addNode(string nodeName){
        newNode = new node();
        newNode->name = nodeName;
        if(this->getListSize() != 0){
        tail->next = newNode;
        tail = tail->next;
        }else{
            header = newNode;
            tail = newNode;
        }
    };
    int getListSize(){
        node *curr;
        if(header){
            curr = header;
            int counter = 0;
            while(curr){
                counter++;
                curr = curr->next;
            }
            return counter;
        }else{
            return 0;
        }
    };
};

struct breadboard{
    /*
    "nodes" and "edges" make up our adjacency list for our graph structure. This structure has 
    a hash map,edges, that is used to keep track of edges/components and their connected nodes.
    The vector list contains nodes and their connected edges (but not information on the secondary
    nodes connected to these edges)

    every edge/component has the possiblity of being included in calculations hence there
    is a need for storing properties such as resistances, voltages... for every component
    this will be a seperate hashmap that accepts the edges names but is linked to the edges "value"
    */
     unordered_map<string, vector<int>> edges; // this contains information on the nodes and their related components (i.e. resistor R1 is connected to nodes 1 and 3)
     vector<list> nodes;
     unordered_map<string, int> edgesValue; // by default all edges will have a value of zero.
     int lineCounter = 0; // counts the number of lines "specially added" in the circuit  due to components in parallel
     void addEdge(string name, vector<int> edgeNodes){
        if(!edges.count(name)){
            edges[name] = edgeNodes;
        }
     };
     void setNewNode(list newNode){
         nodes.push_back(newNode);
     };
     /*
     Adding the first component to the bread board
     */
     void addComponent(string componentName){
         list newEdgeNode1 = list();
         list newEdgeNode2 = list();
        //  cout << " Current Size: "<< nodes.size() << endl;
         int newIndex = nodes.size();
         
         newEdgeNode1.addNode(componentName); //although it says "addNode" we are basically relating the component/edge to its node.
         newEdgeNode2.addNode(componentName); // if the component is "R4" then we have two nodes related to it.

         setNewNode(newEdgeNode1);
         setNewNode(newEdgeNode2);
         vector<int> edgeNodes = {newIndex, newIndex + 1}; // creating nodes for the edges/components that is being added
         addEdge(componentName,edgeNodes);    
         assignValue(componentName, 0);  // default value of zero 0   
     };
     /*
     Adding second(or any higher order) component to the breadboard but not closing a loop
     */
     void addComponent(string componentName, int mergeNode1){
         list newEdgeNode1 = list();
         int prevIndex = nodes.size() - 1;

         nodes[mergeNode1].addNode(componentName); // this is for when there are already components on the breadboard and we want to add a new component and only one new node.
         //although it says "addNode" we are basically relating the component/edge to its node.
         //if the component is "R4" then we have two nodes related to it.
         newEdgeNode1.addNode(componentName); 
         setNewNode(newEdgeNode1);
         addEdge(componentName, {mergeNode1, prevIndex + 1});
         assignValue(componentName, 0); // default value of zero 0
     };
     /*
     Adding higher order component but closing a loop.
     */
     void addComponent(string componentName, int mergeNode1, int mergeNode2){ // here we are "bridging a gap and there is no new node added."
     
        // nodes[mergeNode1].addNode(componentName);
        // nodes[mergeNode2].addNode(componentName);
        // addEdge(componentName, {mergeNode1, mergeNode2});
        /*
            Update: Need to make sure that if a component is added in parallel with another component, two lines (assuming 2 node component) is added such that components are in series with a 
            zero resistance line.
            Check if there exists an edge that is already at given node1 and node2.(i.e. they are already adjacent to each other)
        */
        vector<int> adj = getAdjacentNodes(mergeNode1);
        bool alreadyAdj = false;
        for(auto anyElem : adj){
            if(anyElem == mergeNode2){
                alreadyAdj = true;
            }
        }
        if(alreadyAdj){
            int lineIndex= nodes.size();                   // this is for assigning the edges of our component (since the new component will be in parllel but it is "displaced" by a line hence not the same nodes). Need to create new nodes for this component

            list newEdgeNode1 = list();                   // create end node for line 1
            setNewNode(newEdgeNode1);
            list newEdgeNode2 = list();                   // create end node for line 2
            setNewNode(newEdgeNode2);
            string line1 = "L" + to_string(lineCounter);  // name for line 1

            nodes[mergeNode1].addNode(line1);             // add line1 to start node
            addEdge(line1,{mergeNode1,lineIndex});        // create "component" line1 with start node and new end node

            lineCounter++;                                // we need another (unique) name for line 2
            string line2 = "L" + to_string(lineCounter);
             lineCounter++;                                // we need another (unique) name for line in the future
            nodes[mergeNode2].addNode(line2);             //add line2 to start node (other end of component we are parallel to) 
            addEdge(line2,{mergeNode2,lineIndex+1});      // actually create line2 "component" to start node and new end node.

            // Now we can add the component in parallel.
            int componentNode1 = lineIndex;
            int componentNode2 = lineIndex + 1;
            nodes[componentNode1].addNode(line1);  // Adding line1 node as adjacent nodes to component
            nodes[componentNode2].addNode(line2);  // Adding line2 node as adjacent nodes to component

            addComponent(componentName, lineIndex,lineIndex+1); // component itself is now added
            
        }
        else{
            nodes[mergeNode1].addNode(componentName);
            nodes[mergeNode2].addNode(componentName);
            addEdge(componentName, {mergeNode1, mergeNode2});
        }

        assignValue(componentName, 0); // default value of zero 0
        
     };

     void assignValue(string componentName, int value){
        if(edges.count(componentName)){
            edgesValue[componentName] = value;
        }
     }

     bool hasCycle(){

        /*
        Will use the BFS traversal method
        Make use of queue
        set a flag for all vertices/nodes:
        -1 node is unvisited
        0 node is in queue
        1 node is visited

        set all nodes to -1
        Can start from any node so I choose to start from node 0 (easy to start with i.e. array structure)
        insert first node in queue
        set flag to 0 for first node
        pop the queue.
        for the given node that was popped
        change node status to 1 and,
        add all the adjacent nodes of popped node to the queue. except the nodes that have a status of 1.
        change these adjacent nodes flags to 0, however if one of the nodes is already 0 then we have a cycle.
        
        repeat 
        */
        unordered_map<int,int> marker;
        queue<int> visited;
        for (int i = 0; i < nodes.size(); i++){
            marker[i] = -1;
        }
        visited.push(0);
        do{
            int node = visited.front();
            visited.pop();
            marker[node] = 1;
            vector<int> adj = getAdjacentNodes(node);
            for(auto node : adj){
                if(marker[node] == -1){
                    marker[node] = 0;
                    visited.push(node);
                }else if(marker[node] == 0){
                    return true;
                }
            }

        } while (!visited.empty());
        return false;
     }

     
     unordered_map<int,int> marker;
     vector<vector<int>> cycles;
     void getCycles(int lastNode, int parentNode, int par[], int &cyclenumber){

        if(marker[lastNode] == 1){
            return;
        }
        if (marker[lastNode] == 0){ // still visiting
            vector<int> v;
            cyclenumber++;

            int curr = parentNode;
            v.push_back(curr);

            while(curr != lastNode){
                curr = par[curr];
                v.push_back(curr);
            }
            cycles.push_back(v);
            return;
        }
        par[lastNode] = parentNode;
        marker[lastNode] = 0;
        for(int givenNode: getAdjacentNodes(lastNode)){
            if(givenNode == par[lastNode]){
                continue; //skip if we are relooking at a previously visited node.
            }else{
                getCycles(givenNode, lastNode, par, cyclenumber);
            }
        }
        marker[lastNode] = 1;
     }
     bool DFS_detect_cycle(int n){
        
        /*
        colour the vertices (make nodes -1), where 1, 0, 1 -> unvisited, visiting, visited)
        recall for DFS we can simply just make use of a recursive function and use the call stack
        instead of implementing one.
        - Starting with the first node

        for each node:
            if node is already visiting return true;
            otherwise mark as visiting.
            for each adjacent node perform DFS on it.
            perform dfs on it. 
            mark current node as visited.
        */
        marker[n] = 0;
        for(auto eachElem : getAdjacentNodes(n)){
            if(marker[eachElem] == 0)
                return true;
            if(marker[eachElem] == -1 && DFS_detect_cycle(eachElem))
                return true;
        }
         marker[n] = 1;
         return false;
     }
     vector<int> getCommonNodes( vector<int> cycle1, vector<int> cycle2){ //function to check whether two loops are adjacent (by returning the commonNodes)
         vector<int> commonNodes;
         vector<int> ctCycle;
         vector<int> cwCycle;
         int compare = cycle1.size() - cycle2.size();
          if (!isSubset(cycle1, cycle2))
          { // adjacent loops have some nodes in common but are not subsets of one another
             if(compare >= 0){
                ctCycle = cycle2;
                cwCycle = cycle1;
             }else{
                ctCycle = cycle1;
                cwCycle = cycle2; 
             }
                 for(auto eachElem: ctCycle){
                     if(checkElemExist(eachElem,cwCycle)){
                         commonNodes.push_back(eachElem);
                     }
                 }
             }
         
         return commonNodes; // we need at least two nodes to make up an edge// note that regardless of the common edges if cycle1 and cycle 2 have a subset-superset relation they will not be considered adjacent
     }
     vector<int> getAdjacentNodes(int n){ //recall nodes are indicies of the nodes array
        vector<int> adjN;
        node *curr = nodes[n].header;
        while(curr){
            for(auto elem : edges[curr->name]){
                if(elem != n){ //recall we do not want to count the node itself as an adjacent node.
                    bool insert = true;
                    /*
                    we need to avoid a double entry so check the array for if the node is already present. 
                    we could have used a hash map but the expected size of array elements
                    will probably never get to an extremely large number
                    */ 
                    for(auto anyElem: adjN){ 
                       if(anyElem == elem){
                            insert = false;
                        }
                    }
                    if(insert){
                        adjN.push_back(elem);  
                    }
                    
                } 
            }
            curr = curr->next;
        }
        return adjN;
     }  

      vector<string> getEdges(vector<int> edgeNodes){
        /*
            given a list of nodes representing an edge, this function returns the component names 
            of between the nodes of the given edgeNodes.
        */
            vector<string> allEdges;
            //     int ctptr = 0;
            //     int cwptr = 1;
            //  while(cwptr < edgeNodes.size()) {
              
               
            //     string foundedEdge = getPairEdge(edgeNodes[ctptr],edgeNodes[cwptr]);
            //     if (foundedEdge.size() > 0)
            //         allEdges.push_back(foundedEdge); // note if 2 or more components are attached in parallel then there is more than one edge for the pair
            //     cwptr++;
            //      if(cwptr >= edgeNodes.size() && ctptr < edgeNodes.size()-2){
            //      ctptr++;
            //      cwptr = ctptr + 1;
            //     }
            // }
            auto nodePairRep = getPairedNodesEdges(edgeNodes);
            for(auto pair : nodePairRep){
                auto foundedEdge = getPairEdge(pair[0], pair[1]);
                allEdges.push_back(foundedEdge);
            }
            return allEdges;
        }

    /*
     TODO: try changing the return type to just string. 

     Initially I thought multiple components could be connected to two nodes 
     However, I decided that was alot of work so lines (i.e. L1, L2...) were added as new 
     edges of zero resistance when components are added in parallel.
    */
    string getPairEdge(int a, int b){ //for a given node a and b get the edge in between them note we can only have one edge connected between 2 nodes connected.
        string foundedEdge;
        node *curr = nodes[a].header;
        while(curr){
            if(edges[curr->name][0] == b || edges[curr->name][1] == b){ // at a given edge if any of the nodes are linked to b then we know this is an edge for ab
                foundedEdge = curr->name;
                break;
            }
            curr = curr->next;
        }
        return foundedEdge;
    }
    
 
    void removeSuperSets(){
        int cwptr = 1; // Compare With pointer
        int ctptr = 0; // Compare To pointer
        if(cycles.size() > 1) // there needs to be at least 2 cycles
            while (cwptr < cycles.size())
            {
                getCycleRemainder(&cycles[ctptr], &cycles[cwptr]);
                cwptr++;
                if (cwptr >= cycles.size()){
                    ctptr++;
                    cwptr = ctptr + 1;
                }
                
        }
    }
    void reduceHiddenSuperLoops(vector<int>* givenCycle){

        //TODO: 
       /*
        Given a cycle, we already assume it is a hidden loop.
        We then iterate through the cycle to get all the pair of nodes (i.e. edges) that 
        exist in the loop:

        - Now we could perform dijkstra's algorithm to detect the second shortest path (the path apart from the direct edge) and then use
        the second shortest path as the closing end of the loop.

        - however I would like to perform dijkstra's algorithm only once so I need to determine
        what pair to do it on.
       */
  
      
    }

    

    void getCycleRemainder(vector<int> *cycle1, vector<int> *cycle2){ // This function is suppose to "reduce" any existing superset-subset relation into two independent sets. (This function does not take into account hidden superloops)
        // Note: because of dfs approach to cycle detection we have nice cycles that can be easily "disjointed" i.e if we have 3201 and 201 the remaining subcycle/subset is 321 just remove the anything that is 
        if(isSubset(*cycle1,*cycle2)){
            vector<int> replaceSet;
            vector<int> *superset;
            vector<int> *subset;
            if(cycle1->size() > cycle2->size()){ // for cycle 2 to be a subset all its elements have to be within cycle 1
                superset = cycle1;
                subset = cycle2;
            }
            else{
                superset = cycle2;
                subset = cycle1;
            }
           
            for(auto eachElem: *superset){
                bool insert = true;
                for (int i = 1; i < subset->size() - 1; i++){ // if the element of the super set is within any of the mid elements of the sub set then we should not include it. i.e. it is subset' or subset compliment
                    if (eachElem == (*subset)[i])
                        {
                            insert = false;
                        }
                }
                if(insert){
                        replaceSet.push_back(eachElem);
                }               
            }
            *superset = replaceSet; 
        }
    }

     bool isSubset(vector<int> cycle1, vector<int> cycle2){ // Function to check if two cycles have a subset-superset relation
        bool checkAgain = true;
        vector<int> superset;
        vector<int> subset;
       
        if(cycle1.size() == cycle2.size()){
            return false;
        }
        if(cycle1.size() >  cycle2.size()){ // for cycle 2 to be a subset all its elements have to be within cycle 1
            superset = cycle1;
            subset = cycle2;
        }
        else{
            superset = cycle2;
            subset = cycle1;
        }
        
        for(auto eachElem: subset){
            if(checkAgain){
                checkAgain = checkElemExist(eachElem, superset);
                if(eachElem == subset[subset.size() - 1]){
                    return checkAgain;
                }
            }
        }
        return false;
     }

     bool checkElemExist(int n, vector<int> arr){
        for (int i = 0; i < arr.size(); i++){
            if(n == arr[i]){
                return true;
            }
        }
        return false;
     }

     void printOut(){

         cout << "Building circuit... " << endl;
         for (int i = 0; i < nodes.size(); i++)
         {
             list elem = nodes[i];
             node *curr;
             curr = elem.header;
             cout << "components attached to node: " << i << " " << endl;
             while (curr)
             {
                 cout << curr->name << " ";
                 curr = curr->next;
                 if(!curr){
                     cout << "" << endl;
                 }
             }
         }
        cout << " " << endl;
        cout << " " << endl;
        cout << " " << endl;

        //  cout << "Testing adjacent node function" << endl;
         for (int i = 0; i < nodes.size(); i++){
             cout << "Here are the adjacent nodes for node " << i << " : " << endl;
             vector<int> adj= getAdjacentNodes(i);
             for(auto node : adj){
                 cout << node << " ";
             }
             cout << " " << endl;
         }
        cout << " " << endl;
        cout << " " << endl;
        cout << " " << endl;
         cout << "Attempt at detecting cycles..." << endl;
         cout << "BFS: Do we have cycles " << hasCycle() << " !" << endl;

        for (int i = 0; i < nodes.size(); i++){
            marker[i] = -1;
        }
         cout << "DFS: Do we have cycles " << DFS_detect_cycle(0) << " !" << endl;

         int par[nodes.size()];
         int cyclenumber = 0;
         for (int i = 0; i < nodes.size(); i++){
            marker[i] = -1;
        }
         getCycles(1, 0, par, cyclenumber);
         cout << "Now attempting to calculated the number of cycles in circuit : " << cyclenumber << " cycle(s)" << endl;
         for(auto cycle:cycles){
            printArray(cycle);         
            cout << " " << endl;
         }
         vector<int> arr1 = {4,3,2,0,1};
         vector<int> arr2 = {3,2, 0, 1 };
         cout << "Checking if one cycle is a subset of another:" << endl;
         printArray(arr1);
         cout << ", ";
         printArray(arr2);
         cout<< " , subset? : " << isSubset(arr1, arr2) <<  endl;

         cout << "Checking if other subset can be extracted and replace superset:" << endl;
         getCycleRemainder(&arr1, &arr2);
         printArray(arr1);
         cout << ", ";
         printArray(arr2);
         cout<< " , subset? :" << isSubset(arr1, arr2) << endl;

        //  cout << "We want to extract the edges but before that we need to extract the ";

         cout << "Time to extract the edges" << endl;
         displayEdges();

         cout << "Notice that we had a superset and subset, lets clean this up and try again" << endl;
         removeSuperSets();
         displayEdges();

            /*
            --------------------------------------------------------------------------------------------------------
            Testing New Adjacent Cycles function "setCycleAdjacents"
            --------------------------------------------------------------------------------------------------------
            */

            cout << "Here's a list of every cycle and their adjacent cycles: " << endl;
            unordered_map<string, vector<vector<int>>> adjacentLookUp;
            setCycleAdjacents(cycles, adjacentLookUp);
            int adjacencyCounter[cycles.size()] ={}; // integer array
            for (int i = 0; i < cycles.size();i++)
            {
                auto cycle = cycles[i];
                string cycleKey = getKeyForCycle(cycle);
                if(adjacentLookUp.count(cycleKey))
                    adjacencyCounter[i] = adjacentLookUp[cycleKey].size();
                else{
                    adjacencyCounter[i] = 0;
                }
                cout << "listing for : ";
                printArray(cycle);
                cout << endl;
                cout << endl;
                for(auto adjCycle : adjacentLookUp[getKeyForCycle(cycle)]){
                    
                    cout << " ";
                    printArray(adjCycle);
                    cout << " : ";
                    vector<int> commonEdgeNodes = getCommonNodes(cycle, adjCycle);
                    vector<string> commonEdgeNames = getEdges(commonEdgeNodes);
                    printArray(commonEdgeNames);
                    cout << " -> ";
                    printArray(commonEdgeNodes);

                    cout << endl;
                }
                cout << endl;
            }

            /*
            --------------------------------------
            Potential Solution to hidden loop problem.
            --------------------------------------
            */

            int maxCount = 0;
            int maxCountIndex = 0;
            for (int i = 0; i < cycles.size(); i++){
                if (adjacencyCounter[i] > cycles[i].size())
                {

                    if(maxCount < adjacencyCounter[i]){
                        maxCountIndex = i; // storing index of the max count
                        maxCount = adjacencyCounter[i];
                    }
                    printArray(cycles[i]); 
                    cout << " Is an invalid loop with " << adjacencyCounter[i] << " counts" << endl;
                }
            }

            cout << "MaxCount Index : " << adjacencyCounter[maxCountIndex] << endl;
            printArray(cycles[maxCountIndex]);
            cout << endl;

            /*
            TODO: work on reduceHiddenSuperLoops function later

              reduceHiddenSuperLoops(&cycles[maxCountIndex]);
            */

           /*
           ------------------------------------------------------------------------------
             Calculating loop equations
           ------------------------------------------------------------------------------
           */

           // Check values are assigned properly to components.
            unordered_map<string, int>::iterator it = edgesValue.begin();
            while(it != edgesValue.end())
            {
                string componentName = it->first;
                cout << "Component: " << componentName << " has a value of -> " << edgesValue[componentName] << endl;
                it++;
           }

           // Calculating equivalent resisitance of a single loop.
           // For each pair within a cycle get the edge and add it's value to a sum if(it is a resistor)
           // A cycle will have it's voltage so if we see a component that starts with v we save that value is a special variable.
           cout << endl;
           cout << "Let's get the properties of these loops: " << endl;
           cout << endl;
          vector<any> loopEquations = getLoopEquations(cycles);
          cout << endl;
          for(auto matrix : loopEquations){
                if(matrix.type() == typeid(vector<int>))
                    {
                        cout << "Voltage Source Column : " << endl;
                        printArrayVertical(any_cast<vector<int>>(matrix)); 
                    }
                   
                else
                    {
                        cout << "Resitance Matrix Column : " << endl;
                        for(auto row : any_cast<vector<vector<int>>>(matrix)){
                        printArray(row);
                        cout << endl;
                        }
                    }
                cout << endl;
          }
        }
        /*
        --------------------------------------------------------------------------
            printOut() HELPER FUNCTIONS
        --------------------------------------------------------------------------        
        */
        void displayEdges(){
            for(auto cycle : cycles){
            printArray(cycle);
            cout << " " << endl;
           
            vector<vector<int>> pairedNodes = getPairedNodes(cycle);

            for(auto pair : pairedNodes){
                 int node1 = pair[0];
                 int node2 = pair[1];

                auto pairEdges =  getPairEdge(node1, node2); //note if 2 or more components are attached in parallel then there is more than one edge for the pair... Nevermind I made it so that we will add lines i.e. L1 and L2 so now every pair would have only one edge between them
                cout << node1 << " and " << node2 << " ";
                cout << pairEdges << endl;
            }
         }
        }
        void printArrayVertical(vector<int> arr){
            for (auto elem : arr)
            {
                cout << elem << endl;
            }
        }
        string printArray(vector<int> arr){
            string array = "";
            for (auto elem : arr)
            {
                array += elem + " ";
                cout << elem << " ";
            }
            return array;
        }
        string printArray(vector<string> arr){
         string array = "";
         for (auto elem : arr)
         {
            array += elem + " ";
            cout << elem << " ";
         }
         return array;
        }

      
        /*
        --------------------------------------------------------------------------------
        LOOP Calculating Functions
        --------------------------------------------------------------------------------
        */
        
        vector<any> getLoopEquations(vector<vector<int>> givenCycles){
            /*
                given all the cycles, an empty cycleIndexLookUP, and filled adjacentCycleLookUp table that gives you the adjacent cycles for one cycle:

                first initiallize a 2d vector for the loop equations (every row should be the length of given cycles)
                assign the index of every cycle, hence their loop number, to cycleIndexLookup
                next loop through all cycles again and for every cycle (loop equation) add 

              
                cycleIndexLookUp : initially empty, will return the index/loop number of the cycle 
                givenCycles : this is a list of all the cycles in the circuit. (loop numbers are assigned from here)
            */
            unordered_map<string, int> cycleIndexLookUp;
            vector<vector<int>> resistanceMatrix;
            vector<int> voltageVector(givenCycles.size(),0);
            for (int i = 0; i < givenCycles.size(); i++){
                // setting index of every cycle
                vector<int> cycle = givenCycles[i];
                cycleIndexLookUp[getKeyForCycle(cycle)] = i;


                // setting up row equations for every cycle/loop/loopequation
                vector<int> equationRow(givenCycles.size(), 0); // initiallized with all zeros.
                resistanceMatrix.push_back(equationRow);

            }

            // Get list of adjacent cycles, adjacentCycleLookUp : returns the adjacent cycles to a cycle.
            unordered_map<string, vector<vector<int>>> adjacentCycleLookUp;
            setCycleAdjacents(givenCycles, adjacentCycleLookUp);

            // Now we have our indicies and we have allocated space for our loop equations let's generate the actual equations:
            for(int i = 0; i < givenCycles.size(); i++){
                vector<int> cycle = givenCycles[i];

                //Add the main loop properties to the loop equation
                printArray(cycle);
                    
                vector<int> sum = singleCycleSum(cycle);
                int resistorSum = sum[0];
                int voltageSum = sum[1];
                int mainRow = i;
                int mainCol = cycleIndexLookUp[getKeyForCycle(cycle)];
                // add the mainloop reistance
                resistanceMatrix[mainRow][mainCol] = resistorSum;
                // add mainloop voltage
                voltageVector[mainRow] = voltageSum;

                cout << " loop resistance: " << resistorSum  << " loop voltage: " << voltageSum<< endl;

                //Add the adjacent loop properties to the loop equation
                vector<vector<int>> adjacentCycles = adjacentCycleLookUp[getKeyForCycle(cycle)];
                for(auto adjCycles: adjacentCycles){
                    int row = i;
                    int col = cycleIndexLookUp[getKeyForCycle(adjCycles)];
                    resistanceMatrix[row][col] = getContributedResistance(cycle, adjCycles);
                }
            }
            vector<any> results;
            results.push_back(resistanceMatrix);
            results.push_back(voltageVector);
            return results;
        }

        int getContributedResistance(vector<int> mainCycle, vector<int> adjacentCycle) {

                    vector<int> commonEdgeNodes = getCommonNodes(mainCycle, adjacentCycle);
                    vector<string> commonEdgeNames = getEdges(commonEdgeNodes);
                    /*
                        I need to get the common resistor names first then I can get the polarity. There's no need to get the polarity for the loop voltage,
                        Only the voltage in the original cycle is counted for the loop equation.

                        We only need the compared polarity for one of the resistors along the common edge (the comparedPolarity, not the polarity, should be the same for the other components)
                        Ex: 
                            commonEdge: R1, R2 where R1 {0,1} , R2 {2, 1}
                            cycle1 : 0 1 2 4 5, cycle2 : 8 2 1 0 7 

                            note the polarity of R1 and R2 for cycle1 are 1 and -1 repectively
                            For cycle 2: -1 and 1 respectively.

                            However the compared polarity : polarity of cycle1 * polarity of cycle2 is -1 and -1 for both R1 and R2 (i.e. the assumed currents of each cycle are in the opposit directions)
                    */
                    cout << "  Contribution from : ";
                    printArray(adjacentCycle);
                    int comparedPolarity;
                    for(auto componentName : commonEdgeNames){
                        if(componentName[0] == 'R'){
                           comparedPolarity = getPolarityWithCycle(adjacentCycle,componentName) * getPolarityWithCycle(mainCycle,componentName); // -1 if different, 1 if the same, 0 if component does not exist in one of the cycles
                           break;
                        }
                    }                  
                    vector<int> Adjsum = sumEdges(commonEdgeNodes);
                    cout << " loop resistance: " << Adjsum[0] << ", with relative polarity: " << comparedPolarity << endl;
                    return comparedPolarity*Adjsum[0];
        }

        void setCycleIndex(vector<vector<int>> givenCycles, unordered_map<string,int> cycleIndexLookUp){
            /*
                When making loop equations I have to keep track of the certain loop that the equation is for.
                given a hashmap I can store the indices once and then assign the loop equations to their appropriate places. (not loop equation building NOT here)
            */
            for (int i = 0; i < givenCycles.size(); i++){
                vector<int> cycle = givenCycles[i];
                cycleIndexLookUp[getKeyForCycle(cycle)] = i;
            }
        }
        string cycleToString(vector<int> givenCycle){
         return getKeyForCycle(givenCycle);
        }

        string getKeyForCycle(vector<int> givenCycle){
            /*
            This fucntion is suppose to help me convert a cycle into a string so that I can easily use hash maps with strings
            */
         string key;
            for (int node : givenCycle)
            {
                key += to_string(node);
            }
            return key;
        }
        void setCycleAdjacents(vector<vector<int>> givenCycles,unordered_map <string, vector<vector<int>>> &adjacentCycleLookUp ){
            /*
            Function to update a hashmap of cycles where the key is a string version of a cycle and the value is an array of the other cycles
            */
            int ctptr = 0;                            // compare to pointer
            int cwptr = 1;                            // compare with pointer
            while(cwptr < cycles.size()){
                vector<int> commonEdgeNodes = getCommonNodes(givenCycles[ctptr], givenCycles[cwptr]);
               
                if(commonEdgeNodes.size() > 1){         // that means we have atleast an edge (we need 2 nodes to make an edge)
                    adjacentCycleLookUp[getKeyForCycle(givenCycles[ctptr])].push_back(givenCycles[cwptr]);
                    adjacentCycleLookUp[getKeyForCycle(givenCycles[cwptr])].push_back(givenCycles[ctptr]);
                }
                cwptr++;
                if(cwptr >= givenCycles.size() && ctptr < givenCycles.size()-2){
                    ctptr++;
                    cwptr = ctptr + 1;
                }
            }
        }
        vector<int> singleCycleSum(vector<int> givenCycle){
            /*
            Returns a sum of the resistace in the cycle and the voltage in the cycle.
            */
            vector<vector<int>> cyclePairs = getPairedNodes(givenCycle);
            int resistanceSum = 0;
            int voltageSum = 0;
           
            for(auto pair : cyclePairs){
                string componentName= getPairEdge(pair[0], pair[1]); // again only vector of size one.
                char indicator = componentName[0];
                bool resistor = indicator == 'R';
                bool source = indicator == 'V';
                int componentValue = edgesValue[componentName];
                if(resistor){
                    resistanceSum += componentValue;
                }
                else if (source){
                   /*
                   checking polarity:
                   pair: {0,3}
                   V_name: {0,3}, therefore component in the positive configuration (otherwise negative).
                   */
                    voltageSum += getPolarity(pair,componentName)*componentValue;
                    
                }
            }
           
         return {resistanceSum,voltageSum};
        }

        vector<int> sumEdges(vector<int> givenEdge){
            /*
            Returns a sum of the resistace and voltage along an Edge. "givenEdge" is a list of nodes that make up an edge.
            */
            vector<vector<int>> cyclePairs = getPairedNodesEdges(givenEdge);
            int resistanceSum = 0;
            int voltageSum = 0;
           
            for(auto pair : cyclePairs){
                string componentName= getPairEdge(pair[0], pair[1]); // get component that exist between pair of nodes
                char indicator = componentName[0];
                bool resistor = indicator == 'R';
                bool source = indicator == 'V';
                int componentValue = edgesValue[componentName];
                if(resistor){
                    resistanceSum += componentValue;
                }
                else if (source){
                   // TODO: have to account for polarities in voltages.
                   /*
                   checking polarity:
                   pair: {0,3}
                   V_name: {0,3}, therefore component in the positive configuration (otherwise negative).
                   */
                    voltageSum += getPolarity(pair,componentName)*componentValue;
                    
                }
            }
           
         return {resistanceSum,voltageSum};
        }

        int getPolarityWithCycle(vector<int> givenCycle, string componentName){
            /*
            same as getPolarity() but it accepts a given cycle instead of a pair of nodes
            We treat polarities as factors to multiply by, i.e. (1, -1, 0)
            */
         auto pairs = getPairedNodes(givenCycle);
         for(auto pair : pairs){
            if(getPairEdge(pair[0],pair[1]) == componentName){
                return getPolarity(pair, componentName);
            }
         }
         return 0; // if the component does not exist within the given cycle then return a factor of 0;
        }
        
        int getPolarity(vector<int> pair, string componentName){
            /*
            given an edge/compent name, and the "in-cycle-direction" of how the component/edge was traversed 
            (i.e. did it visit the left node first then the right or vice-versa)
            this function determines if the in-cycle-direction is aligned with the way the component was added to the edge unordered set 
            (when we add a component one node is on the left and the other is on the right)
            Ex: hw.addComponent("V2",6,3} , givencycle -> [0 1 3 6]  in this case V2 has a negative polarity because we connected to 6 first before connecting to 3, but the cycle
            visits 3 first before visiting 6
            */  
        //    edges[componentName];
            if(pair[0] == edges[componentName][0] && pair[1] == edges[componentName][1]){  //
                return 1;
            }
            else{
                return -1;
            }
        }

        vector<vector<int>> getPairedNodes(vector<int> givenCycle){
            /*
             getPairedNodes is a function that takes in a cycle and iterates through that cycle
             to get all the pairs of nodes that make up an edge in that cycle:

             givenCycle : {0, 1, 2, 3}
             cycleNodePairs: {{0,1},{1,2},{2,3},{3,0}}
            */
         vector<vector<int>> cycleNodePairs;
         for (int i = 0; i < givenCycle.size(); i++)
         {
             int node1 = givenCycle[i % (givenCycle.size())];
             int node2 = givenCycle[(i + 1) % (givenCycle.size())];

             vector<int> nodePair = {node1, node2};
             cycleNodePairs.push_back(nodePair);
            }
            return cycleNodePairs;
        }

        vector<vector<int>> getPairedNodesEdges(vector<int> givenEdge){
            /*
             getPairedNodesEdges is a function that takes in an Edge and iterates through that edge
             to get all the pairs of nodes that make up the edge:

             givenEdge : {0, 1, 2, 3}
             edgeNodePairs: {{0,1},{1,2},{2,3}}
            */
         vector<vector<int>> edgeNodePairs;
         for (int i = 1; i < givenEdge.size(); i++)
         {
             int node1 = givenEdge[i-1];
             int node2 = givenEdge[i];

             vector<int> nodePair = {node1, node2};
             edgeNodePairs.push_back(nodePair);
            }
            return edgeNodePairs;
        }




       
};

int main(){
    
    breadboard hw1 = breadboard();
    
/*
Circuit with parallel components.
*/
    // hw1.addComponent("R1");  
    // hw1.addComponent("R2",1); 
    // hw1.addComponent("R3",2); 
    // hw1.addComponent("R4",3,0); 
    // hw1.addComponent("R5",3); 
    // hw1.addComponent("R6",4); 
    // hw1.addComponent("R7",5); 
    // hw1.addComponent("R8",6); 
    // hw1.addComponent("R9",7,4);
    // hw1.addComponent("R10", 6);
    // hw1.addComponent("R11", 8);
    // hw1.addComponent("R12", 9, 7);
    // hw1.addComponent("R13", 8, 9);
    // hw1.addComponent("R14", 10, 9);

/*
    Generic circuit example
*/
    // hw1.addComponent("R1");  
    // hw1.addComponent("R2",1); 
    // hw1.addComponent("R3",2); 
    // hw1.addComponent("R4",3,0);
    // hw1.addComponent("R5", 0);
    // hw1.addComponent("R6", 4, 1);

/*
    Testing simple subset example.
*/
    // hw1.addComponent("R1");
    // hw1.addComponent("R2",1);
    // hw1.addComponent("R3", 2);
    // hw1.addComponent("R4", 3);
    // hw1.addComponent("R5", 4);
    // hw1.addComponent("R6", 5);
    // hw1.addComponent("R7", 6,0);
    // hw1.addComponent("R8", 1, 6);
/*
    Complex subset example
*/
    // hw1.addComponent("R1");
    // hw1.addComponent("R2", 1);
    // hw1.addComponent("R3", 0,2);
    // hw1.addComponent("R4", 1);
    // hw1.addComponent("R5", 2, 3);
    // hw1.addComponent("R6", 1);
    // hw1.addComponent("R7", 3, 4);

/*
    Real-life example : This was an assigment from ECE 2205A - failed
*/

    // hw1.addComponent("R1");    // 0,1
    // hw1.addComponent("R2", 1); // 1,2
    // hw1.addComponent("R4", 2); // 2,3
    // hw1.addComponent("R3", 2, 3); // 4,5 used for lines
    // hw1.addComponent("R5",2);     // 2,6
    
    // hw1.addComponent("R8", 3);   // 3,7
    // hw1.addComponent("R9", 7);   // 6,8
    // hw1.addComponent("R6", 0,2); // 9 10 used for lines
    // hw1.addComponent("R7", 6, 3);
    // hw1.addComponent("R12", 0, 6);
    // hw1.addComponent("R10", 8, 3);
    // hw1.addComponent("R11", 6, 8);
    // hw1.addComponent("R13", 0, 8);
    // hw1.addComponent("R14", 0, 8);

    /* Simple loop detection practice*/

    // hw1.addComponent("R1");
    // hw1.assignValue("R1", 5);
    // hw1.addComponent("R2",1);
    // hw1.assignValue("R2", 1);
    // hw1.addComponent("R3", 2);
    // hw1.assignValue("R3", 2);
    // hw1.addComponent("V1", 3); // for voltage sources we assume the first element in the edge vector is the negative side.
    // hw1.assignValue("V1", 1);
    // hw1.addComponent("Vab", 0, 4); // for voltage sources we assume the first element in the edge vector is the negative side.
    // hw1.assignValue("Vab", 5);


    /*  multicomponent adjacent edge - failed*/
    // hw1.addComponent("R1");
    // hw1.addComponent("R2", 1);
    // hw1.addComponent("R3", 2);
    // hw1.addComponent("R4", 3);
    // hw1.addComponent("R5", 4);
    // hw1.addComponent("R6", 5,0);
    // hw1.addComponent("R8", 4);
    // hw1.addComponent("R7", 6, 1);

    /* 3 adjacent loop case*/

    // hw1.addComponent("R1");
    // hw1.assignValue("R1", 7);
    // hw1.addComponent("V2", 1);
    // hw1.assignValue("V2", 15);
    // hw1.addComponent("R2", 2);
    // hw1.assignValue("R2", 5);
    // hw1.addComponent("V1", 3);
    // hw1.assignValue("V1", 7);
    // hw1.addComponent("R3", 4);
    // hw1.assignValue("R3", 4);
    // hw1.addComponent("V3", 5);
    // hw1.assignValue("V3", -10);
    // hw1.addComponent("R4", 6);
    // hw1.assignValue("R4", 1);
    // hw1.addComponent("R5", 7, 0);
    // hw1.assignValue("R5", 8);
    // hw1.addComponent("R6", 0);
    // hw1.assignValue("R6", 4);
    // hw1.addComponent("R7", 2, 8);
    // hw1.assignValue("R7", 1);
    // hw1.addComponent("R8", 4, 8);
    // hw1.assignValue("R8", 3);

    // Will be testing with this cirucit for a while (original circuit from notes in Notion)

    hw1.addComponent("Ve");
    hw1.assignValue("Ve", 24);
    hw1.addComponent("R1", 1);
    hw1.assignValue("R1", 200);
    hw1.addComponent("R2", 2, 0);
    hw1.assignValue("R2", 200);
    hw1.addComponent("R3", 2);
    hw1.assignValue("R3", 50);
    hw1.addComponent("R4", 3);
    hw1.assignValue("R4", 150);
    hw1.addComponent("L1", 0, 4);

    hw1.printOut();
}

