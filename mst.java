import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Scanner;


/**
 * The mst class can generate a random graph, implement Minimum spanning tree using simple and fibonacci schema.
 * @since Oct, 25, 2013 
 * @author Da Zhao
 *
 */
public class mst {
	
	private static String mode;
	private static int numVertex;
	private static double density;
	private static String filename;
	public final static int MAXINT = Integer.MAX_VALUE;
	
	/**
	 * The inner class Graph is used to generate a graph from user's input or from given number of vertex and density.
	 * Using this class to generate the graph, one can get a graph that not too much uniformity, say, 
	 * the degree of a node can be varies a lot. See details in createGraph().
	 * @author Da Zhao
	 *
	 */
	private class Graph {
		private int numEdges;
		// The adjacent list, the ith position means vertex i. 
		// In HashMap, the first integer j means i connecting to j, the second integer c means the cost.
		public ArrayList<HashMap<Integer, Integer>> list;
		
		/**
		 * Constructor of Graph, initialize the adjacent list, and get the number of edges.
		 */
		private Graph(){
			list = new ArrayList<HashMap<Integer, Integer>>();
			for(int i = 0; i < numVertex; i++){
				list.add(new HashMap<Integer, Integer>());
			}
			// The number of edges from input
			Double withDensity = numVertex * (numVertex - 1) * density / 200;
			// The number of edges should not be less than numVertex -1, so that it is a connected graph using Prim
			numEdges = Math.max(numVertex - 1, withDensity.intValue());
//			System.out.println("The number of edges is: " + numEdges);
		}
		
		/**
		 * Insert an edge with vertex and cost
		 * @param vertexA
		 * @param vertexB
		 * @param cost
		 */
		private void insert(int vertexA, int vertexB, int cost){
			list.get(vertexA).put(vertexB, cost);
			list.get(vertexB).put(vertexA, cost);
		}
		
		/**
		 * Get the cost of an edge
		 * @param vertexA
		 * @param vertexB
		 * @return The cost of the edge
		 */
		public int getWeight(int vertexA, int vertexB){
/*			if(vertexA == vertexB){
				throw new Exception("two vertex are identical!");
			}*/
			return list.get(vertexA).get(vertexB);
		}
		
		/**
		 * Check if edge is in the graph
		 * @param vertexA
		 * @param vertexB
		 * @return true/false for existence of the edge
		 */
		public boolean containsEdge(int vertexA, int vertexB){
			if(vertexA < vertexB){
				return list.get(vertexA).containsKey(vertexB);
			}
			else{
				return list.get(vertexB).containsKey(vertexA);
			}
		}
		
		/**
		 * Create the graph randomly
		 */
		private void createGraph(){
			Random random = new Random();
			// if density is 100%, we directly create the graph
			if(density == 100){
				for(int i = 0; i < numVertex; i++){
					for(int j = i + 1; j < numVertex; j++){
						int cost = random.nextInt(1000) + 1;
						insert(i,j,cost);
					}
				}
			}
			// else density < 100%
			else{
				// step 1: Initialize two arrays. One has all vertex from graph, one is empty.
				ArrayList<Integer> tempList = new ArrayList<Integer>();
				ArrayList<Integer> QList = new ArrayList<Integer>();
				for(int i = 0; i < numVertex; i++){
					QList.add(i);
				}
				
				// step 2: cut the node from QList to tempList, until all nodes are added to tempList.
				// Then a connected graph is generated.
				int currQNode = random.nextInt(numVertex);
				int currTNode;
				tempList.add(currQNode);
				QList.remove(currQNode);
				while(!QList.isEmpty()){
					currTNode = tempList.get(random.nextInt(tempList.size()));
					int QPos = random.nextInt(QList.size());
					currQNode = QList.get(QPos);
					int cost = random.nextInt(1000) + 1;
					insert(currTNode, currQNode, cost);
					tempList.add(currQNode);
					QList.remove(QPos);
				}
				Collections.sort(tempList);

				// step 3: for the remaining edges, randomly select any two nodes that not added previously.
				int remaining = numEdges - numVertex + 1;
				int threshold = numVertex - 1;
				for(int i = 0; i < remaining; i++){
					int posA = random.nextInt(tempList.size());
					int posB = random.nextInt(tempList.size());
					int vertexA = tempList.get(posA);
					int vertexB = tempList.get(posB);
					int counter = 0;
					// If more than one node in tempList, we randomly generate two different vertex
					if(tempList.size() != 1){
						// We keep posA so that the graph can not be so uniformly. Otherwise, more probability will give to vertex with less edges.
						while(posA == posB){
							posB = random.nextInt(tempList.size());
						}
						vertexB = tempList.get(posB);
						while(containsEdge(vertexA, vertexB)){
							do{
								posB = random.nextInt(tempList.size());
							}while(posA == posB);
							vertexB = tempList.get(posB);
							++counter;
							// If the random duplicated edges reaches threshold, it indicates the vertex is nearly full linked.
							// We use a loop to check which vertex it does not point to.
							if(counter == threshold){
								for(int j = 0; j < tempList.size(); j++){
									if(!containsEdge(vertexA, tempList.get(j))){
										if(posA != j){
											posB = j;
											vertexB = tempList.get(posB);
											break;
										}
									}
								}
								break;
							}
						}
						
						int cost = random.nextInt(1000) + 1;
						insert(vertexA, vertexB, cost);
						// Remove the vertex if it is full connected, will not random it later.
						int buffer = 0;
						if(list.get(vertexA).size() == numVertex - 1){
							tempList.remove(posA);
							buffer = 1;
						}
						if(list.get(vertexB).size() == numVertex - 1){
							// Using buffer when posA is deleted, since posB will be offset 
							if(posB > posA){
								tempList.remove(posB - buffer);
							}
							else{
								tempList.remove(posB);
							}
						}
					}
					// If tempList size is 1, we directly generate the edge. This is the rare case.
					else{
						for(int j = 0; j < numVertex; j++){
							if(!containsEdge(vertexA, j)){
								if(vertexA != j){
									vertexB = j;
									break;
								}
							}
						}
						int cost = random.nextInt(1000) + 1;
						insert(vertexA, vertexB, cost);
					}
				}
			}
		}
	}
	
	/**
	 * The node used for both simple and fibonacci schema
	 * @author Da Zhao
	 *
	 */
	private class Node {
		int id;
		int key;
		int parent;
	}
	
	/**
	 * The f-node for fibonacci schema
	 * @author Da Zhao
	 *
	 */
	private class FibonacciNode {
		int degree; // num of children
		boolean marked;
		FibonacciNode next;
		FibonacciNode prev;
		FibonacciNode parent;
		FibonacciNode child;
		Node node;	
		private FibonacciNode(){
			node = new Node();
		}
	}

	/**
	 * The inner class Simple used for simple schema
	 * @author Da Zhao
	 *
	 */
	private class Simple {
		Graph graph = null;
		ArrayList<Node> nodes = null;
		ArrayList<Node> result = null;
		HashMap<Integer, Node> qMap = null;
		/**
		 * Constructor of Simple, initialize all nodes and containers.
		 * @param aGraph
		 */
		private Simple(Graph aGraph) {
			graph = aGraph;
			qMap = new HashMap<Integer, Node>();
			nodes = new ArrayList<Node>();
			result = new ArrayList<Node>();
			// initialize the node in the graph
			for(int i = 0; i < numVertex; i++){
				Node node = new Node();
				node.id = i;
				node.key = MAXINT;
				node.parent = -1;
				nodes.add(node);
				qMap.put(i, node);
				result.add(node);
			}
		}
		
		/**
		 * The entrance of simple schema
		 */
		public void start(){
			Random random = new Random();
			int currNode = random.nextInt(numVertex);
			qMap.get(currNode).key = 0;
			while(qMap.size() != 0){
				// Extract the node(vertex) with minimum key
				Node minNode = extractMin();
				currNode = minNode.id;
				// For every neighbor of the node, decrease the key if it is less than that in the container.
				for (Map.Entry<Integer, Integer> entry : graph.list.get(currNode).entrySet()) {
				    int neighborNode = entry.getKey();
				    if(qMap.containsKey(neighborNode)){
					    int currWeight = graph.getWeight(neighborNode, currNode);
				    	if(currWeight < qMap.get(neighborNode).key){
					    	result.get(neighborNode).key = currWeight;
					    	result.get(neighborNode).parent = currNode;
				    	}
				    }
				}
			}
			int sum = 0;
			for(Node aNode : result){
				if(aNode.key != MAXINT){
					sum += aNode.key;
				}
			}
			System.out.println(sum);
			for(Node aNode : result){
				if(aNode.key != 0){
					System.out.println(aNode.parent + " " + aNode.id);
				}
			}
		}
		
		/**
		 * Check if id is in the array, not used since using HashMap. We need O(1) time in this, to better compare with fibonacci schema.
		 * @param id
		 * @param array
		 * @return The index of the node with the given nodeID
		 */
		private int arrayContain(int id, ArrayList<Node> array){
			for(Node aNode : array){
				if(aNode.id == id) {
					return array.indexOf(aNode);
				}
			}
			return -1;
		}
		
		/**
		 * To get the minimum key in the array, this will take O(V) time for V nodes
		 * @return The node with minumum key
		 */
		private Node extractMin(){
			int min = nodes.get(0).key;
			Node minNode = nodes.get(0);
			int pos = 0;
			for (int i = 1; i < nodes.size(); i++){
				if(min > nodes.get(i).key){
					min = nodes.get(i).key;
					minNode = nodes.get(i);
					pos = i;
				}
			}
			nodes.remove(pos);
			qMap.remove(minNode.id);
			return minNode;
		}
	}
	
	/**
	 * The inner class Fibonacci used for fibonacci schema
	 * @author Da Zhao
	 *
	 */
	private class Fibonacci {
		Graph graph = null;
		FibonacciNode minFNode = null;
		ArrayList<Node> result = null;
		HashMap<Integer, FibonacciNode> qMap = null;
		/**
		 * Constructor of Fibonacci, initialize all nodes and containers.
		 * @param aGraph
		 */
		private Fibonacci(Graph aGraph){
			qMap = new HashMap<Integer, FibonacciNode>();
			result = new ArrayList<Node>();
			graph = aGraph;
			for(int i = 0; i < numVertex; i++){
				FibonacciNode fNode = new FibonacciNode();
				fNode.node.id = i;
				fNode.node.key = MAXINT;
				fNode.node.parent = -1;
				fNode.parent = null;
				insert(fNode);

			}
		}
		
		/**
		 * Insert the node to container
		 * @param aFNode
		 */
		private void insert(FibonacciNode aFNode){
			if(aFNode == null) return;
			aFNode.marked = false;
			aFNode.degree = 0;
			aFNode.parent = null;
			aFNode.child = null;
			
			if(minFNode == null){
				minFNode = aFNode;
				minFNode.next = minFNode;
				minFNode.prev = minFNode;
			}
			else{
				aFNode.prev = minFNode.prev;
				aFNode.next = minFNode;
				aFNode.prev.next = aFNode;
				aFNode.next.prev = aFNode;
				if(aFNode.node.key < minFNode.node.key){
					minFNode = aFNode;
				}
			}
			qMap.put(aFNode.node.id, aFNode);
			result.add(aFNode.node);
		}
		
		/**
		 * To get the minimum key in the f-heap, this will take O(logV) amortized time for V nodes
		 * @return The f-node with minimum key
		 */
		private FibonacciNode extractMin(){
			if(minFNode == null) {System.out.println("minFNode is null..."); return null;}
			FibonacciNode resultFNode = minFNode;
			FibonacciNode startFNode = new FibonacciNode();
			FibonacciNode currFNode = new FibonacciNode();
			qMap.remove(minFNode.node.id); 
			
			// step1: remove the minNode, combine its children to top level
			// The min node has child(ren)
			if(minFNode.degree != 0){
				startFNode = minFNode.child;
				startFNode.parent = null;
				currFNode = startFNode.next;
				while(startFNode != currFNode){
					currFNode.parent = null;
					currFNode = currFNode.next;
				}
				currFNode = currFNode.prev;
				// add children to top level and combine with others
				if(minFNode.next != minFNode){
					startFNode.prev = minFNode.prev;
					minFNode.prev.next = startFNode;
					currFNode.next = minFNode.next;
					minFNode.next.prev = currFNode;
				}
			}
			// The min node has no child, but has sibling, extract minNode directly
			else if(minFNode.next != minFNode){
				startFNode = minFNode.prev;
				startFNode.next = minFNode.next;
				minFNode.next.prev = startFNode;
			}
			// There is only one node before extractMin, no child, no sibling
			else{ 
				startFNode = minFNode;
				minFNode = null;
				return startFNode; 
			}
			
			// step2: link roots with same degree
			startFNode = combineRoots(startFNode);
			
			// step3: find the minimum key in top level
			int min = startFNode.node.key;
			minFNode = startFNode;
			currFNode = startFNode.next;
			while(currFNode != startFNode){
				if(min > currFNode.node.key){
					min = currFNode.node.key;
					minFNode = currFNode;
				}
				currFNode = currFNode.next;
			}
			return resultFNode;
		}
		
		private FibonacciNode combineRoots(FibonacciNode startFNode){
			HashMap<Integer, FibonacciNode> degreeMap = new HashMap<Integer, FibonacciNode>();
			degreeMap.put(startFNode.degree, startFNode);
			FibonacciNode currFNode = startFNode.next;
			FibonacciNode toMatchNode = currFNode;
			while(currFNode != startFNode){
				// The degree is not shown before
				int toMatchDegree = toMatchNode.degree;
				if(!degreeMap.containsKey(toMatchDegree)){
					degreeMap.put(toMatchNode.degree, currFNode);
				}
				// There is a node with the same degree before. combine them!
				else{
					while(degreeMap.containsKey(toMatchDegree)){
						FibonacciNode matchedNode = degreeMap.get(currFNode.degree);
						degreeMap.remove(toMatchDegree);
						if(matchedNode.node.key < toMatchNode.node.key){
							// remove currFNode and add it to the child of matchedNode, keeping the currFNode and startFNode at top level
							if(toMatchNode == currFNode){
								currFNode = currFNode.prev;
							}
							if(toMatchNode == startFNode){
								startFNode = startFNode.next;
							}
							combineNode(matchedNode,toMatchNode);
							toMatchNode = matchedNode;
						}
						else{
							// remove matchedNode and add it to the child of currFNode, keeping the startFNode at top level
							if(matchedNode == startFNode){
								startFNode = startFNode.next;
							}
							combineNode(toMatchNode, matchedNode);
							toMatchDegree += 1;
						}
					}
					degreeMap.put(toMatchDegree, toMatchNode);
				}
				currFNode = currFNode.next;
				toMatchNode = currFNode;
			}
			return startFNode;
		}
		
		/**
		 * Remove nodeB and add it to the child of nodeA
		 * @param nodeA
		 * @param nodeB
		 */
		private void combineNode(FibonacciNode nodeA, FibonacciNode nodeB){
			// remove nodeB
			nodeB.prev.next = nodeB.next;
			nodeB.next.prev = nodeB.prev;

			// add nodeB to the child of nodeA
			if(nodeA.degree != 0){
				nodeB.prev = nodeA.child.prev;
				nodeB.next = nodeA.child;
				nodeB.parent = nodeA;
				nodeA.child.prev.next = nodeB;
				nodeA.child.prev = nodeB;
				
			}
			else{
				nodeA.child = nodeB;
				nodeB.parent = nodeA;
				nodeB.next = nodeB;
				nodeB.prev = nodeB;
			}
			nodeA.degree += 1;
		}
		
		/**
		 * Update the key of the node with new value
		 * @param aFNode
		 * @param newKey
		 */
		private void decreaseKey(FibonacciNode aFNode, int newKey){
			if (newKey < aFNode.node.key) {
				aFNode.node.key = newKey;
				FibonacciNode parent = aFNode.parent;
	            if (parent != null && aFNode.node.key < parent.node.key) {
	                cut(aFNode, parent);
	            }
	            if (aFNode.node.key < minFNode.node.key) {
	                minFNode = aFNode;
	            }
	        }
		}
		
		/**
		 * Cut the current node to top level, update the marked. Also, do the cascadeCut for parent if possible
		 * @param toCut
		 * @param parent
		 */
		private void cut(FibonacciNode toCut, FibonacciNode parent){
			if(toCut.next != toCut){
				toCut.prev.next = toCut.next;
				toCut.next.prev = toCut.prev;
			}
			if(parent.degree > 1 && parent.child == toCut){
				parent.child = toCut.next;
			}
			else if(parent.degree == 1){
				parent.child = null;
			}
			parent.degree -= 1;
			toCut.parent = null;
			
			toCut.prev = minFNode.prev;
			toCut.next = minFNode;
			toCut.prev.next = toCut;
			toCut.marked = false;
			minFNode.prev = toCut;
	        if (toCut.parent != null) {
	            if (!toCut.parent.marked) {
	            	toCut.parent.marked = true;
	            }
	            else {
	                cut(toCut.parent,parent);
	            }
	        }
		}
		
		/**
		 * Find the nodeID level by level. This is not used, and replaced by using HashMap, to better compare with simple schema.
		 * @param nodeID
		 * @param startFNode
		 * @return The f-node with given nodeID
		 */
		private FibonacciNode findFNode(int nodeID, FibonacciNode startFNode){
			if(startFNode == null) return null;
			if(startFNode.node.id == nodeID){
				return startFNode;
			}
			if(startFNode.child != null){
				FibonacciNode resultFNode = findFNode(nodeID, startFNode.child);
				if(resultFNode != null){
					return resultFNode;
				}
			}
			FibonacciNode currFNode = startFNode.next;
			while(currFNode != startFNode){
				if(currFNode.node.id == nodeID){
					return currFNode;
				}
				if(currFNode.child != null){
					FibonacciNode resultFNode = findFNode(nodeID, currFNode.child);
					if(resultFNode != null){
						return resultFNode;
					}
				}
				currFNode = currFNode.next;
			}
			return null;
		}
		

		/**
		 * The entrance of fibonacci schema
		 */
		public void start(){
			Random random = new Random();
			int currFNodeID = random.nextInt(numVertex);
			minFNode = qMap.get(currFNodeID);
			minFNode.node.key = 0;
			while(qMap.size() != 0){
				// Extract the node(vertex) with minimum key		
				FibonacciNode currFNode = extractMin();
				// For every neighbor of the node, decrease the key if it is less than that in the container.
				for (Map.Entry<Integer, Integer> entry : graph.list.get(currFNode.node.id).entrySet()) {
				    int neighborNode = entry.getKey();
				    if(qMap.containsKey(neighborNode)){
					    int currWeight = graph.getWeight(neighborNode, currFNode.node.id);
				    	if(currWeight < qMap.get(neighborNode).node.key){
				    		decreaseKey(qMap.get(neighborNode), currWeight);
					    	result.get(neighborNode).key = currWeight;
					    	result.get(neighborNode).parent = currFNode.node.id;
				    	}
				    }
				}
			}
			int sum = 0;
			for(Node aNode : result){
				if(aNode.key != MAXINT){
					sum += aNode.key;
				}
			}
			System.out.println(sum);
			for(Node aNode : result){
				if(aNode.key != 0){
					System.out.println(aNode.parent + " " + aNode.id);
				}
			}
		}
	}
	
	/**
	 * The starting point of processing the program
	 */
	private void process(){
		if("-r".equals(mode)){
			// generate the graph
			Graph graph = new Graph();
			graph.createGraph();
			
			// process Simple
			long start = System.currentTimeMillis();
			Simple simple = new Simple(graph);
			simple.start();
			System.out.println("Simple time cost: " + (System.currentTimeMillis() - start));
			
			// process Fibonacci
			start = System.currentTimeMillis();
			Fibonacci fibonacci = new Fibonacci(graph);
			fibonacci.start();
			System.out.println("Fibonacci time cost: " + (System.currentTimeMillis() - start));
		}
		
		else if("-s".equals(mode)){
			Graph graph = readFromFile(filename);
			// process Simple
			long start = System.currentTimeMillis();
			Simple simple = new Simple(graph);
			simple.start();
			System.out.println("Simple time cost: " + (System.currentTimeMillis() - start));
		}
		else if("-f".equals(mode)){
			Graph graph = readFromFile(filename);
			// process Fibonacci
			long start = System.currentTimeMillis();
			Fibonacci fibonacci = new Fibonacci(graph);
			fibonacci.start();
			System.out.println("Fibonacci time cost: " + (System.currentTimeMillis() - start));
		}
	}
	
	/**
	 * Using scanner to read from file with Integer values
	 * @param fileName
	 * @return the generated graph from input
	 */
	private Graph readFromFile(String fileName){
		Scanner sc = null;
		Graph graph = null;
		try {
			sc = new Scanner(new File(fileName));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		int numV = sc.nextInt();
		int numE = sc.nextInt();
		numVertex = numV;
		graph = new Graph();
		graph.numEdges = numE;
		while (sc.hasNextInt()) {
			int vertexA = sc.nextInt();
			int vertexB = sc.nextInt();
			int cost = sc.nextInt();
			graph.list.get(vertexA).put(vertexB, cost);
			graph.list.get(vertexB).put(vertexA, cost);
		}
		sc.close();
		return graph;
	}
	
	/**
	 * @param argv
	 * @throws Exception
	 */
	public static void main (String argv[]) throws Exception{
		
		// parse the input parameters
		parseInput(argv);
		
		// process the MST
		mst myMST = new mst();
		myMST.process();
	}
	
	/**
	 * Validate and parse the input
	 * @param argv
	 * @throws Exception
	 */
	private static void parseInput(String[] argv) throws Exception{
		if(argv.length != 2 && argv.length != 3){
			System.out.println("Please use the correct input format:" +
					" \"mst -r n d\" or \"mst -s file-name\" or \"mst -f file-name\"");
		}
		try{
			// random mode
			if("-r".equals(mode = argv[0])){
				numVertex = Integer.valueOf(argv[1]);
				density = Double.valueOf(argv[2]);
//				System.out.println("In random mode, the number of vertex is: " + numVertex
//						+ ", and the density is: " + density + "%");
			}
			// user input mode with simple schema
			else if("-s".equals(mode = argv[0])){
				filename = argv[1];
//				System.out.println("In user input mode with simple schema, the file name is: " + filename);
			}
			// user input mode with fibonacci heap schema
			else if("-f".equals(mode = argv[0])){
				filename = argv[1];
//				System.out.println("In user input mode with f-heap schema, the file name is: " + filename);
			}
			else{
				throw new Exception();
			}
		} catch(Exception e){
			System.out.println("Failure occurs when parsing input: " + e);
			throw e;
		}
	}
}