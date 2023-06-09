import java.util.ArrayList;
import java.util.ArrayDeque;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

import tester.*;
import javalib.impworld.*;
import java.awt.Color;
import javalib.worldimages.*;
import java.util.Random;

/* DOCUMENTATION:
 * Run this game using the run configurations used in every homework thus far, making sure
 * the jar files are imported and in the right directory.
 * 
 * Play this game by following the instructions on the screen.
 * Once you click run, a screen will pop up with the maze and some 
 * key-press instructions on the right hand side
 * Have fun playing!
 */

// Vertex class which represents a singular cell in our maze
class Vertex {
  // we are unable to make this final because this is a state of a Vertex,
  // which can be changed depending on whether we visited or not
  private boolean visited; // If this vertex has been visited
  // unable to make this final because a Vertex's state of having a player could
  // change
  // if a player decides to move onto it or not
  private boolean hasPlayer; // If this vertex currently has the player
  // unable to make this final because a vertex may or may not be part of the
  // final path
  private boolean partOfProductivePath; // If this vertex is part of final path

  // these fields can't be final because they are getting reassigned during the
  // initialization
  // process for computation of gradient levels
  private int gradientLevelStart;
  private int gradientLevelEnd;

  // y : represents column position, negative means more left, positive means more
  // "right"
  // x: row position (negative means more "up"-> positive means more "down")
  private final Posn posn; // Position in grid

  // the neighbor fields can't be final since when are changing them during
  // initialization
  private Vertex above;
  private Vertex below;
  private Vertex left;
  private Vertex right; // Neighbors

  // Convenience Constructor
  Vertex(Posn posn) {
    this.visited = false;
    this.hasPlayer = false;
    this.partOfProductivePath = false;
    this.gradientLevelStart = 10000;
    this.gradientLevelEnd = 10000;
    this.posn = posn;
    this.above = this;
    this.below = this;
    this.left = this;
    this.right = this;
  }

  // Convenience Constructor to set neighbors
  Vertex(Posn posn, Vertex above, Vertex below, Vertex left, Vertex right) {
    this.visited = false;
    this.hasPlayer = false;
    this.partOfProductivePath = false;
    this.gradientLevelStart = 10000;
    this.gradientLevelEnd = 10000;
    this.posn = posn;
    this.above = above;
    this.below = below;
    this.left = left;
    this.right = right;
  }
  
  // Big convenience constructor to facilitate testing. This way, we can create example vertices 
  // with predetermined attributes to then test behavior
  Vertex(boolean visited, boolean hasPlayer, boolean partOfProductivePath, 
      int gradientLevelStart, int gradientLevelEnd, Posn posn, Vertex above,
      Vertex below, Vertex left, Vertex right) {
    this.visited = visited;
    this.hasPlayer = hasPlayer;
    this.partOfProductivePath = partOfProductivePath;
    this.gradientLevelStart = gradientLevelStart;
    this.gradientLevelEnd = gradientLevelEnd;
    this.posn = posn;
    this.above = above;
    this.below = below;
    this.left = left;
    this.right = right;
  }

  // Moves player onto this vertex, as used by the move method in the maze class
  void acceptAPlayer() {
    this.visited = true;
    this.hasPlayer = true;
  }

  // Moves player off of this vertex, as used by the move method in the maze class
  void removeAPlayer() {
    this.hasPlayer = false;
  }

  // Accepts this vertex as a member of the productive path in our maze(it'll now
  // render dark blue)
  // We need this in the maze class to backtrack our solution.
  void addToProductivePath() {
    this.partOfProductivePath = true;
  }

  // Resets this vertex to its original standing
  // (where it isn't visited, doesn't have player, and isn't part of productive
  // path)
  void reset() {
    this.visited = false;
    this.hasPlayer = false;
    this.partOfProductivePath = false;
  }

  // Creates a connection between this vertex and the other vertex
  void createConnection(Vertex other) {
    other.createConnectionHelper(this, this.posn);
  }

  // Helper for above method: only mutates THIS vertex's neighbors
  void createConnectionHelper(Vertex v2, Posn otherPosition) {
    if (this.posn.y > otherPosition.y) {
      this.left = v2;
    }
    else if (this.posn.y < otherPosition.y) {
      this.right = v2;
    }
    else if (this.posn.x > otherPosition.x) {
      this.above = v2;
    }
    else if (this.posn.x < otherPosition.x) {
      this.below = v2;
    }
  }

  // Returns list of all valid neighbors (adjacent vertices that aren't blocked)
  ArrayList<Vertex> getValidNeighbors() {
    ArrayList<Vertex> validNeighbors = new ArrayList<>();

    if (this.above != this) {
      validNeighbors.add(this.above);
    }
    if (this.left != this) {
      validNeighbors.add(this.left);
    }

    if (this.right != this) {
      validNeighbors.add(this.right);
    }
    if (this.below != this) {
      validNeighbors.add(this.below);
    }

    return validNeighbors;
  }

  // Moves the player in the specified direction (off of this vertex if possible).
  Vertex movePlayer(String dir) {
    Vertex destination = this;
    if (dir.equals("down")) {
      destination = this.below;
    }
    else if (dir.equals("up")) {
      destination = this.above;
    }
    else if (dir.equals("right")) {
      destination = this.right;
    }
    else if (dir.equals("left")) {
      destination = this.left;
    }
    this.removeAPlayer();
    destination.acceptAPlayer();
    return destination;
  }

  // Updates distance from start/end if the new value is smaller
  void updateDistanceFrom(boolean fromStart, boolean initialVertex, Vertex neighbor) {
    int newDistance = neighbor.gradientLevelStart;
    if (!fromStart) {
      newDistance = neighbor.gradientLevelEnd;
    }
    // If we're at the initial vertex, we set distance to 0
    if (initialVertex) {
      newDistance = 0;
    }
    else {
      newDistance += 1;
    }
    // Update depending on if we want start or end
    if (fromStart) {
      this.gradientLevelStart = Math.min(this.gradientLevelStart, newDistance);
    }
    else {
      this.gradientLevelEnd = Math.min(this.gradientLevelEnd, newDistance);
    }
  }

  // Draws this vertex
  WorldImage draw(int rows, int cols, boolean hideVisited, boolean showProductivePath,
      int gradientStatus) {

    // Default color is gray.
    // These if statements override each other in order of precedence.
    // This design is clean and eliminates nested statements.
    Color cellColor = Color.GRAY;
    if (visited && !hideVisited) {
      cellColor = new Color(0, 0, 182, 155);
    }
    if (hasPlayer) {
      cellColor = Color.red;
    }
    if (showProductivePath && this.partOfProductivePath) {
      cellColor = Color.blue;
    }
    if (this.posn.x == 0 && this.posn.y == 0) {
      cellColor = Color.green;
    }
    if (this.posn.x == rows - 1 && this.posn.y == cols - 1) {
      cellColor = Color.MAGENTA;
    }

    // Create gradient of colors that scales based on graph size
    if (gradientStatus == -1) {
      int change = Math.max(1, 255 / Math.max(1, (rows * cols / 4))) * this.gradientLevelEnd;
      cellColor = new Color(Math.max(0, 255 - change), 0, Math.min(change, 255));
    }

    if (gradientStatus == 1) {
      int change = Math.max(1, 255 / Math.max(1, (rows * cols / 4))) * this.gradientLevelStart;
      cellColor = new Color(Math.max(0, 255 - change), 0, Math.min(change, 255));
    }

    // Scale cell based on size of maze
    int cellWidth = Math.min(800 / rows, 1200 / cols);

    // Draw borders if they're there
    WorldImage cellPicture = new RectangleImage(cellWidth, cellWidth, OutlineMode.SOLID, cellColor);
    if (this.left == this) {
      cellPicture = cellPicture.movePinhole(-1 * (cellWidth / 2 - 1), 0);
      cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
          new RectangleImage(2, cellWidth, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
          .movePinhole(cellWidth / 2 - 1, 0);
    }
    if (this.right == this) {
      cellPicture = cellPicture.movePinhole(cellWidth / 2 - 1, 0);
      cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
          new RectangleImage(2, cellWidth, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
          .movePinhole(-1 * (cellWidth / 2 - 1), 0);
    }
    if (this.above == this) {
      cellPicture = cellPicture.movePinhole(0, -1 * (cellWidth / 2 - 1));
      cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
          new RectangleImage(cellWidth, 2, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
          .movePinhole(0, cellWidth / 2 - 1);
    }
    if (this.below == this) {
      cellPicture = cellPicture.movePinhole(0, cellWidth / 2 - 1);
      cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
          new RectangleImage(cellWidth, 2, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
          .movePinhole(0, -1 * (cellWidth / 2 - 1));
    }
    return cellPicture;

  }

}

// Edge class represents an edge between two vertices with a weight
class Edge {
  // It doesn't make sense to make the fields of this class private.
  // The sole purpose of this class is to store two vertices and a weight.
  // Not to perform operations with the vertices.
  // That role is delegated to other classes.
  // The MazeBuilder class accesses the vertices of edges to initialize the maze.
  // Better to do them there instead of here.
  final Vertex v1;
  final Vertex v2;
  private final int weight;

  // Convenience constructor
  Edge(Vertex v1, Vertex v2, int weight) {
    this.v1 = v1;
    this.v2 = v2;
    this.weight = weight;
  }

  // Compare edge weight to another --> used in comparator
  int compareTo(Edge other) {
    return this.weight - other.weight;
  }
}

// Maze class
// Represents a maze with 2d arraylist of vertices and a few other valuable fields
class Maze {
  // row major order
  private final ArrayList<ArrayList<Vertex>> vertices;
  private final int rows; // Size of maze (rows x cols)
  private final int cols;
  private final Vertex start; // Starting vertex (top left)
  // currentSpot is unable to final because it represents which vertex a player is
  // on,
  // which could change
  private Vertex currentSpot; // Vertex the player (or search algo) is currently on
  private final Vertex end; // Ending vertex (bottom right)
  private final ArrayDeque<Vertex> allVisited; // Every vertex that has been visited
  private final ArrayDeque<Vertex> productivePath; // Vertices part of the final path

  // Convenience constructor
  Maze(ArrayList<ArrayList<Vertex>> vertices) {
    this.vertices = vertices;
    this.rows = vertices.size();
    this.cols = vertices.get(0).size();
    this.start = this.vertices.get(0).get(0);
    this.currentSpot = start;
    this.currentSpot.acceptAPlayer();
    this.end = this.vertices.get(rows - 1).get(cols - 1);
    this.allVisited = new ArrayDeque<Vertex>();
    this.allVisited.add(this.currentSpot);
    this.productivePath = new ArrayDeque<Vertex>();
    this.productivePath.add(this.currentSpot);
  }

  // Renders the maze
  WorldImage render(boolean hideVisited, boolean showProductivePath, int gradientStatus) {
    // Loops through and renders each vertex in each row
    WorldImage all = new CircleImage(0, OutlineMode.SOLID, Color.BLACK);
    for (int i = 0; i < this.vertices.size(); i += 1) {
      WorldImage row = new CircleImage(0, OutlineMode.SOLID, Color.BLACK);
      for (int j = 0; j < this.vertices.get(0).size(); j += 1) {
        Vertex curr = this.vertices.get(i).get(j);
        row = new BesideImage(row,
            curr.draw(rows, cols, hideVisited, showProductivePath, gradientStatus));
      }
      all = new AboveImage(all, row);
    }
    return all;
  }

  // Moves the player in the proper direction, updating the productive path list
  // in the process.
  // Returns 1 if the move was unproductive (already visited this cell)
  // and 0 if the move was productive.
  int movePlayer(String dir) {
    return this.handleNewSpot(this.currentSpot.movePlayer(dir));
  }

  // Handles movement to a new spot
  // Updates visited list and mutates productive path if this move is bad
  int handleNewSpot(Vertex newSpot) {
    if (newSpot != currentSpot) {
      this.currentSpot = newSpot;
      // If we've already visited this vertex, that means the last move was bad
      // Technically, it means every move between the last visit and this visit was
      // bad
      // But, if you go from X --> Y --> Z. You must retrace the same path
      // backwards...
      // ..to get back to X because there are no cycles. So removing last will
      // suffice.
      if (this.productivePath.contains(newSpot)) {
        this.productivePath.removeLast();
        return 1;
      }
      else {
        this.productivePath.add(newSpot);
      }
      this.allVisited.add(newSpot);
    }
    return 0;
  }

  // Returns if the maze is complete (we're at the end)
  boolean isCompleted() {
    return this.currentSpot == this.end;
  }

  // Shows another element of the final path (backtracking from end, one by one on
  // each tick).
  boolean showAnotherElementOfProductivePath() {
    if (this.productivePath.size() == 0) {
      return false;
    }
    else {
      Vertex last = this.productivePath.removeLast();
      last.addToProductivePath();
      return true;
    }

  }

  // Clears all progress made thus far (starts you back at the beginning).
  void clearProgress() {
    for (ArrayList<Vertex> row : this.vertices) {
      for (Vertex v : row) {
        v.reset();
      }
    }
    this.currentSpot = this.start;
    this.allVisited.clear();
    this.productivePath.clear();
  }

  // This method performs a search on our maze and mutates the visited and
  // productivePath deques.
  // It is only ever run with a cleared maze, so we assume visited and
  // productive path both initially contain just the start spot
  // If dfs = true then performs dfs, if dfs = false then performs bfs
  // Returns whether or not we found a path
  boolean findPath(boolean dfs) {
    ArrayDeque<Vertex> worklist = new ArrayDeque<>();
    // The keys of the hashmap represent each of the vertices in the graph
    // The value of the hashmap represent the a vertex that is CONNECTED to the
    // vertex
    // in the respective key location. This allows us to traverse backwrads to find
    // the
    // entire path to a given Vertex
    HashMap<Vertex, Vertex> connections = new HashMap<>();
    worklist.add(this.currentSpot);

    while (!worklist.isEmpty()) {
      Vertex curr = worklist.removeFirst();
      this.allVisited.add(curr);
      this.currentSpot = curr;
      if (this.currentSpot == end) {
        while (connections.containsKey(curr)) {
          this.productivePath.addFirst(curr);
          curr = connections.get(curr);
        }
        return true;
      }
      ArrayList<Vertex> validNeighbors = curr.getValidNeighbors();

      for (Vertex neighbor : validNeighbors) {
        if (!this.allVisited.contains(neighbor)) {
          if (dfs) {
            worklist.addFirst(neighbor);
          }
          else {
            worklist.addLast(neighbor);
          }
          connections.put(neighbor, curr);
        }
      }
    }
    return false;

  }

  // Visits the next vertex in the list of all visited vertices.
  // The vertex in question will now be colored on the next tick.
  boolean colorNextVisited() {
    if (this.allVisited.isEmpty()) {
      return false;
    }

    Vertex nextVisited = this.allVisited.removeFirst();
    // We simulate a player stepping on and leaving the vertex in question.
    // This is "visiting" that vertex.
    nextVisited.acceptAPlayer();
    nextVisited.removeAPlayer();

    return true;
  }

  // Initializes distance from start/end in each vertex
  void computeGradient(boolean fromStart) {
    Vertex curr = this.start;
    if (!fromStart) {
      curr = this.end;
    }

    curr.updateDistanceFrom(fromStart, true, curr);

    ArrayList<Vertex> visited = new ArrayList<>();
    ArrayDeque<Vertex> worklist = new ArrayDeque<>();
    worklist.add(curr);

    // Breadth first search:
    // the first time we encounter a vertex is via the closest path to the start/end
    while (!worklist.isEmpty()) {
      Vertex parent = worklist.removeFirst();
      // Loop through all neighbors and update distance
      for (Vertex neighbor : parent.getValidNeighbors()) {
        if (!visited.contains(neighbor)) {
          worklist.addLast(neighbor);
          neighbor.updateDistanceFrom(fromStart, false, parent);
          visited.add(neighbor);
        }
      }
    }
  }

}

// Class that builds a maze 
class MazeBuilder {
  private final int rows;
  private final int cols;
  private final ArrayList<ArrayList<Vertex>> vertices; // List of vertices
  private final ArrayList<Edge> allPossibleEdges;

  private final HashMap<Vertex, Vertex> connections; // Connections for Kruskal's algo
  // Every vertex in the graph is represented by a key
  // and the values are the representative Vertex for its group in kruskals
  // algorithm (union)
  // so we can find which representative the key Vertex points at
  private final ArrayList<Edge> finalizedEdges; // Edges present in the final graph

  // Convenience Constructor
  MazeBuilder(int rows, int cols) {
    this.rows = rows;
    this.cols = cols;
    this.vertices = new ArrayList<>();
    this.allPossibleEdges = new ArrayList<>();
    this.connections = new HashMap<>();
    this.finalizedEdges = new ArrayList<>();
  }

  // Generates a random maze.
  // If bias = -1, then horizontally favoring bias.
  // If bias = 0, no bias. If bias = 1, vertically favoring bias
  Maze generateRandomMaze(int bias) {
    this.initializeEdgeBetweenEveryCell(bias);
    this.performKruskal();
    for (Edge edge : this.finalizedEdges) {
      Vertex v1 = edge.v1;
      Vertex v2 = edge.v2;
      v1.createConnection(v2);
      v2.createConnection(v1);
    }
    return new Maze(this.vertices);

  }

  // generate random maze for randomness testing with a seed
  Maze generateRandomMazeSeed(int bias, int seed) {
    this.initializeEdgeBetweenEveryCellSeed(bias, seed);
    this.performKruskal();
    for (Edge edge : this.finalizedEdges) {
      Vertex v1 = edge.v1;
      Vertex v2 = edge.v2;
      v1.createConnection(v2);
      v2.createConnection(v1);
    }
    return new Maze(this.vertices);

  }

  // Generates a fully walled maze without any bias.
  Maze generateFullyWalledMaze() {
    this.initializeEdgeBetweenEveryCell(0);
    this.performKruskal();
    return new Maze(this.vertices);
  }

  // Knocks down one wall of the maze
  // (until there aren't any more walls to knock down while maintaining the MST)
  boolean knockDownWall() {
    if (this.finalizedEdges.isEmpty()) {
      return false;
    }
    else {
      Edge edge = this.finalizedEdges.get(this.finalizedEdges.size() - 1);
      this.finalizedEdges.remove(this.finalizedEdges.size() - 1);
      Vertex v1 = edge.v1;
      Vertex v2 = edge.v2;
      v1.createConnection(v2);
      v2.createConnection(v1);
      return true;
    }

  }

  // Performs kruskal's algo to generate list of final edges to form min spanning
  // tree
  ArrayList<Edge> performKruskal() {
    Collections.sort(this.allPossibleEdges, new LeastEdge());
    int index = 0;
    // As long as we have less than n - 1 edges (n = # of vertices).
    // And as long as we don't run out of all possible edges
    while (finalizedEdges.size() < rows * cols - 1 && index < allPossibleEdges.size()) {
      Edge edge = this.allPossibleEdges.get(index);

      Vertex ep1 = edge.v1;
      while (this.connections.get(ep1) != ep1) {
        ep1 = this.connections.get(ep1);
      }

      Vertex ep2 = edge.v2;
      while (this.connections.get(ep2) != ep2) {
        ep2 = this.connections.get(ep2);
      }

      if (!ep1.equals(ep2)) {
        this.connections.put(ep1, ep2); // Connect endpoints to create one unified bubble
        this.connections.put(edge.v1, ep2);
        // Pseudo Path compression to shorten searches
        // better than nothing

        this.finalizedEdges.add(edge);
      }
      index += 1;
    }
    return this.finalizedEdges;
  }

  // getter purely for testing initializeEdgeBetweenEveryCell ; we couldn't find a
  // better alternative since
  // the list of edges is private
  ArrayList<Edge> getAllPossible() {
    return this.allPossibleEdges;
  }

  // Initializes edge between every cell (with bias if specified)
  void initializeEdgeBetweenEveryCell(int bias) {
    for (int i = 0; i < this.rows; i += 1) {
      ArrayList<Vertex> row = new ArrayList<>();
      for (int j = 0; j < this.cols; j += 1) {
        Vertex vertex = new Vertex(new Posn(i, j));
        row.add(vertex);
        this.connections.put(vertex, vertex);
        // We create an edge between every cell and its neighbor using this method
        if (i != 0) {
          // Bias is applied to make the edge weight smaller/larger
          // based on what we're being biased towards
          this.allPossibleEdges.add(new Edge(row.get(j), this.vertices.get(i - 1).get(j),
              (int) (Math.random() * 100 + 1 - 10 * bias))); // bias added/subtracted on
        }
        if (j != 0) {
          this.allPossibleEdges.add(
              new Edge(row.get(j), row.get(j - 1), (int) (Math.random() * 100 + 1 + 10 * bias)));
        }
      }
      this.vertices.add(row);
    }

    // Now, we have exactly one edge connecting each cell to each of its neighbors.
    // # of edges here = ((rows - 1) * cols + (cols - 1) * rows))
  }

  // initializeEdgeBetweenEveryCell for testing randomness
  void initializeEdgeBetweenEveryCellSeed(int bias, int rand) {
    for (int i = 0; i < this.rows; i += 1) {
      ArrayList<Vertex> row = new ArrayList<>();
      for (int j = 0; j < this.cols; j += 1) {
        Vertex vertex = new Vertex(new Posn(i, j));
        row.add(vertex);
        this.connections.put(vertex, vertex);
        // We create an edge between every cell and its neighbor using this method
        if (i != 0) {
          // Bias is applied to make the edge weight smaller/larger
          // based on what we're being biased towards
          this.allPossibleEdges.add(new Edge(row.get(j), this.vertices.get(i - 1).get(j),
              (int) (rand * 100 + 1 - 10 * bias))); // bias added/subtracted on
        }
        if (j != 0) {
          this.allPossibleEdges
              .add(new Edge(row.get(j), row.get(j - 1), (int) (rand * 100 + 1 + 10 * bias)));
        }
      }
      this.vertices.add(row);
    }

    // Now, we have exactly one edge connecting each cell to each of its neighbors.
    // # of edges here = ((rows - 1) * cols + (cols - 1) * rows))
  }

}

// World class that manages the actual game output
class MazeWorld extends World {
  // maze and mazebuilder can't be final since we are sometimes generating
  // different mazes
  private Maze maze; // Our maze
  private MazeBuilder builder;
  // Maze builder object --> needed for knocking down 1 by 1 since we want to save
  // info
  private final int rows; // rows & cols in maze
  private final int cols;
  // the rest of these fields below also can't be final, since they all represent
  // states of
  // the big bang world, which is constantly changed ontick or via user input,
  private boolean hideVisited; // whether to hide visited (toggle)
  private boolean showProductivePath; // whether we are currently backtracking the final solution
  private boolean allDoneWithManual; // whether we have finished manually completing maze
  private int wrongManualMoves; // number of incorrect manual moves
  private boolean showSearchProgress; // whether we are currently showing dfs/bfs traversal
  private boolean knockingDownWalls; // whether we are currently knocking down walls to build maze
  private int gradientStatus; // -1 for end gradient, 1 for start gradient, 0 for no gradient

  // Convenience constructor
  MazeWorld(int rows, int cols) {
    if (rows > 100 || cols > 60) {
      throw new IllegalArgumentException(
          "Too big! Choose a " + "maze with at most 100 rows and 60 columns");
    }
    else if (rows < 1 || cols < 1) {
      throw new IllegalArgumentException("Must have at least 1 row & column");
    }
    this.maze = new MazeBuilder(rows, cols).generateRandomMaze(0);
    this.builder = new MazeBuilder(rows, cols);
    this.rows = rows;
    this.cols = cols;
    this.hideVisited = false;
    this.showProductivePath = false;
    this.allDoneWithManual = false;
    this.wrongManualMoves = 0;
    this.showSearchProgress = false;
    this.knockingDownWalls = false;
    this.gradientStatus = 0;
  }

  // Convenience constructor for testing randomness
  MazeWorld(Maze maze) {
    this.maze = maze;
    this.builder = new MazeBuilder(2, 2);
    this.rows = 2;
    this.cols = 2;
    this.hideVisited = false;
    this.showProductivePath = false;
    this.allDoneWithManual = false;
    this.wrongManualMoves = 0;
    this.showSearchProgress = false;
    this.knockingDownWalls = false;
    this.gradientStatus = 0;
  }

  // Make scene and render the maze along with instructions on the right side
  public WorldScene makeScene() {
    WorldScene test = new WorldScene(1500, 800);
    WorldImage moves = new TextImage("Wrong Manual Moves: " + this.wrongManualMoves, 25,
        Color.BLUE);
    WorldImage press = new TextImage("Press", 25, Color.BLACK);
    WorldImage dfsInstructions = new TextImage("'d' for DFS solution", 13, Color.BLACK);
    WorldImage bfsInstructions = new TextImage("'b' for BFS solution", 13, Color.BLACK);
    WorldImage toggleInstructions = new TextImage("'t' to toggle hiding visited", 13, Color.BLACK);
    WorldImage createInstructions = new TextImage("'c' to create new maze", 13, Color.BLACK);
    WorldImage restartInstructions = new TextImage("'r' to restart on same maze", 13, Color.BLACK);
    WorldImage hInstructions = new TextImage("'h' to build maze w/ " + "horiz. wall bias", 13,
        Color.BLACK);
    WorldImage vInstructions = new TextImage("'v' to build maze w/ " + "vert. wall bias", 13,
        Color.BLACK);
    WorldImage kInstructions = new TextImage("'k' to animate knocking down walls", 13, Color.BLACK);
    WorldImage sInstructions = new TextImage("'s' for gradient of distance from start", 13,
        Color.BLACK);
    WorldImage eInstructions = new TextImage("'e' for gradient of distance from end", 13,
        Color.BLACK);
    WorldImage doubleClick = new TextImage("click 's' or 'e' again to revert screen", 13,
        Color.BLACK);

    WorldImage allInstructions = new AboveImage(moves, press, dfsInstructions, bfsInstructions,
        toggleInstructions, createInstructions, restartInstructions, hInstructions, vInstructions,
        kInstructions, sInstructions, eInstructions, doubleClick);
    test.placeImageXY(
        this.maze.render(this.hideVisited, this.showProductivePath, this.gradientStatus), 600, 400);
    test.placeImageXY(allInstructions, 1350, 400);
    return test;

  }

  // On each tick, respond according to the current status of the game
  public void onTick() {
    this.makeScene();
    // if we are knocking down walls, we stop if our knock-down attempt is
    // unsuccessful
    if (this.knockingDownWalls) {
      boolean success = this.builder.knockDownWall();
      if (!success) {
        this.resetVariables();
      }
    }
    // During dfs/bfs mode, we are in one of two states:
    // We either have finished and are backtracking final solution
    // Or we are currently showing the path the algo took to the solution
    else if (this.showSearchProgress) {
      if (this.showProductivePath) {
        boolean success = this.maze.showAnotherElementOfProductivePath();
        if (!success) {
          this.resetVariables();
        }
      }
      else {
        boolean success = this.maze.colorNextVisited();
        if (!success) {
          this.showProductivePath = true;
        }
      }
    }
    // If we haven't already finished, we check maze for completion. If so, we
    // backtrack path.
    else if (!this.allDoneWithManual) {
      if (maze.isCompleted()) {
        this.showProductivePath = true;
      }

      if (this.showProductivePath) {
        boolean success = this.maze.showAnotherElementOfProductivePath();
        if (!success) {
          this.allDoneWithManual = true;
        }
      }
    }

  }

  // Respond to key events
  public void onKeyEvent(String key) {
    // Manual player movement
    if (key.equals("right") || key.equals("left") || key.equals("up") || key.equals("down")) {
      if (!this.allDoneWithManual && !this.knockingDownWalls && !this.showSearchProgress) {
        this.wrongManualMoves += this.maze.movePlayer(key);
      }
    }
    // Toggle between hiding and seeing visited cells
    else if (key.equals("t")) {
      this.hideVisited = !this.hideVisited;
    }
    // Create a new random maze without bias (with same dimensions)
    else if (key.equals("c")) {
      this.maze = new MazeBuilder(this.rows, this.cols).generateRandomMaze(0);
      this.resetVariables();
    }
    // Restart on this same maze
    else if (key.equals("r")) {
      if (!this.knockingDownWalls) {
        this.maze.clearProgress();
        this.resetVariables();
      }
    }
    // Perform DFS traversal
    else if (key.equals("d")) {
      if (!this.knockingDownWalls) {
        this.maze.clearProgress();
        this.resetVariables();
        this.maze.findPath(true);
        this.showSearchProgress = true;
      }
    }
    // Perform BFS traversal
    else if (key.equals("b")) {
      if (!this.knockingDownWalls) {
        this.maze.clearProgress();
        this.resetVariables();
        this.maze.findPath(false);
        this.showSearchProgress = true;
      }
    }
    // Create new maze with horizontal bias
    else if (key.equals("h")) {
      this.resetVariables();
      this.maze = new MazeBuilder(this.rows, this.cols).generateRandomMaze(-1);
    }
    // Create new maze with vertical bias
    else if (key.equals("v")) {
      this.resetVariables();
      this.maze = new MazeBuilder(this.rows, this.cols).generateRandomMaze(1);
    }
    // Create new maze by showing each wall being knocked down
    else if (key.equals("k")) {
      this.resetVariables();
      this.builder = new MazeBuilder(this.rows, this.cols);
      this.maze = this.builder.generateFullyWalledMaze();
      this.knockingDownWalls = true;
    }
    // Show gradient of distances from start
    else if (key.equals("s")) {
      if (!this.knockingDownWalls) {
        if (this.gradientStatus == 1) {
          this.gradientStatus = 0;
        }
        else {
          this.maze.computeGradient(true);
          this.gradientStatus = 1;
        }
      }
    }
    // Show gradient of distances from end
    else if (key.equals("e")) {
      if (!this.knockingDownWalls) {
        if (this.gradientStatus == -1) {
          this.gradientStatus = 0;
        }
        else {
          this.maze.computeGradient(false);
          this.gradientStatus = -1;
        }
      }

    }
  }

  // Reset each of the variables --> when we want to start from scratch
  void resetVariables() {
    this.showProductivePath = false;
    this.allDoneWithManual = false;
    this.wrongManualMoves = 0;
    this.showSearchProgress = false;
    this.knockingDownWalls = false;
    this.gradientStatus = 0;
  }
}

// Comparator that compares two edges by weight
class LeastEdge implements Comparator<Edge> {
  public int compare(Edge o1, Edge o2) {
    return o1.compareTo(o2);
  }
}

// Examples class with tests and big bang
class ExamplesMazes {
  // test for bigbang
  void testGame(Tester t) {
    MazeWorld g = new MazeWorld(20, 20);
    g.bigBang(1500, 800, 0.001);
  }

  // test for compare in the LeastEdge Comparator
  boolean testLeastEdgeFunctionObj(Tester t) {
    Comparator<Edge> least = new LeastEdge();
    Vertex v1 = new Vertex(new Posn(1, 1));
    Vertex v2 = new Vertex(new Posn(2, 2));
    Edge heavier = new Edge(v1, v2, 5);
    Vertex v3 = new Vertex(new Posn(1, 1));
    Vertex v4 = new Vertex(new Posn(2, 2));
    Edge lighter = new Edge(v3, v4, 10);

    return t.checkExpect(least.compare(heavier, lighter), -5)
        && t.checkExpect(least.compare(lighter, lighter), 0)
        && t.checkExpect(least.compare(lighter, heavier), 5);
  }

  // test for getVlaidNeighbors in vertex class
  boolean testGetValidNeighbors(Tester t) {
    Vertex v2 = new Vertex(new Posn(2, 2));
    Vertex v3 = new Vertex(new Posn(2, 2));
    Vertex v1 = new Vertex(new Posn(1, 1), v2, v3, v2, v3);

    ArrayList<Vertex> neighbors = new ArrayList<Vertex>();
    neighbors.add(v2);
    neighbors.add(v3);
    neighbors.add(v2);
    neighbors.add(v3);
    ArrayList<Vertex> noNeighbors = new ArrayList<Vertex>();
    return t.checkExpect(v1.getValidNeighbors(), neighbors)
        && t.checkExpect(v3.getValidNeighbors(), noNeighbors);

  }

  // test for onKeyEvent
  void testOnKeyEvent(Tester t) {
    ArrayList<Vertex> row1 = new ArrayList<Vertex>();
    ArrayList<Vertex> row2 = new ArrayList<Vertex>();
    Vertex v1 = new Vertex(new Posn(0, 0));
    Vertex v2 = new Vertex(new Posn(1, 1));
    Vertex v3 = new Vertex(new Posn(1, 1));
    Vertex v4 = new Vertex(new Posn(1, 1));
    row1.add(v1);
    row1.add(v2);
    row2.add(v3);
    row2.add(v4);
    ArrayList<ArrayList<Vertex>> twoD = new ArrayList<ArrayList<Vertex>>();
    twoD.add(row1);
    twoD.add(row2);
    Maze maze = new Maze(twoD);
    MazeWorld mazeWorld = new MazeWorld(maze);

    mazeWorld.onKeyEvent("d");
    WorldScene test = new WorldScene(1500, 800);
    WorldImage moves = new TextImage("Wrong Manual Moves: " + 0, 25, Color.BLUE);
    WorldImage press = new TextImage("Press", 25, Color.BLACK);
    WorldImage dfsInstructions = new TextImage("'d' for DFS solution", 13, Color.BLACK);
    WorldImage bfsInstructions = new TextImage("'b' for BFS solution", 13, Color.BLACK);
    WorldImage toggleInstructions = new TextImage("'t' to toggle hiding visited", 13, Color.BLACK);
    WorldImage createInstructions = new TextImage("'c' to create new maze", 13, Color.BLACK);
    WorldImage restartInstructions = new TextImage("'r' to restart on same maze", 13, Color.BLACK);
    WorldImage hInstructions = new TextImage("'h' to build maze w/ " + "horiz. wall bias", 13,
        Color.BLACK);
    WorldImage vInstructions = new TextImage("'v' to build maze w/ " + "vert. wall bias", 13,
        Color.BLACK);
    WorldImage kInstructions = new TextImage("'k' to animate knocking down walls", 13, Color.BLACK);
    WorldImage sInstructions = new TextImage("'s' for gradient of distance from start", 13,
        Color.BLACK);
    WorldImage eInstructions = new TextImage("'e' for gradient of distance from end", 13,
        Color.BLACK);
    WorldImage doubleClick = new TextImage("click 's' or 'e' again to revert screen", 13,
        Color.BLACK);

    WorldImage allInstructions = new AboveImage(moves, press, dfsInstructions, bfsInstructions,
        toggleInstructions, createInstructions, restartInstructions, hInstructions, vInstructions,
        kInstructions, sInstructions, eInstructions, doubleClick);
    test.placeImageXY(maze.render(false, false, 0), 600, 400);
    test.placeImageXY(allInstructions, 1350, 400);

    t.checkExpect(mazeWorld.makeScene(), test);

    mazeWorld.onKeyEvent("r");
    WorldScene resetWScene = new WorldScene(1500, 800);
    resetWScene.placeImageXY(maze.render(false, false, 0), 600, 400);
    resetWScene.placeImageXY(allInstructions, 1350, 400);
    t.checkExpect(mazeWorld.makeScene(), resetWScene);

    mazeWorld.onKeyEvent("s");
    WorldScene gradientScene = new WorldScene(1500, 800);
    // the third parameter, 1 means show gradient
    gradientScene.placeImageXY(maze.render(false, true, 1), 600, 400);
    gradientScene.placeImageXY(allInstructions, 1350, 400);
    t.checkExpect(mazeWorld.makeScene(), gradientScene);

    mazeWorld.onKeyEvent("e");
    WorldScene gradientScene2 = new WorldScene(1500, 800);
    // the third parameter, -1 means show gradient at end
    gradientScene2.placeImageXY(maze.render(false, true, -1), 600, 400);
    gradientScene2.placeImageXY(allInstructions, 1350, 400);
    t.checkExpect(mazeWorld.makeScene(), gradientScene2);

    mazeWorld.onKeyEvent("r");
    mazeWorld.onKeyEvent("b");
    WorldScene bfsScene = new WorldScene(1500, 800);
    bfsScene.placeImageXY(maze.render(false, true, 0), 600, 400);
    bfsScene.placeImageXY(allInstructions, 1350, 400);
    t.checkExpect(mazeWorld.makeScene(), bfsScene);

    mazeWorld.onKeyEvent("d");
    WorldScene dfsScene = new WorldScene(1500, 800);
    dfsScene.placeImageXY(maze.render(false, true, 0), 600, 400);
    dfsScene.placeImageXY(allInstructions, 1350, 400);
    t.checkExpect(mazeWorld.makeScene(), dfsScene);

    mazeWorld.onKeyEvent("right");
    WorldScene moveRight = new WorldScene(1500, 800);
    maze.movePlayer("right");
    moveRight.placeImageXY(maze.render(false, true, 0), 600, 400);
    moveRight.placeImageXY(allInstructions, 1350, 400);
    t.checkExpect(mazeWorld.makeScene(), moveRight);

    mazeWorld.onKeyEvent("down");
    WorldScene moveDown = new WorldScene(1500, 800);
    maze.movePlayer("down");
    moveDown.placeImageXY(maze.render(false, false, 0), 600, 400);
    moveDown.placeImageXY(allInstructions, 1350, 400);
    t.checkExpect(mazeWorld.makeScene(), moveDown);

    mazeWorld.onKeyEvent("left");
    WorldScene moveLeft = new WorldScene(1500, 800);
    maze.movePlayer("left");
    moveLeft.placeImageXY(maze.render(false, false, 0), 600, 400);
    moveLeft.placeImageXY(allInstructions, 1350, 400);
    t.checkExpect(mazeWorld.makeScene(), moveLeft);

    mazeWorld.onKeyEvent("left");
    WorldScene moveUp = new WorldScene(1500, 800);
    maze.movePlayer("up");
    moveUp.placeImageXY(maze.render(false, false, 0), 600, 400);
    moveUp.placeImageXY(allInstructions, 1350, 400);

    t.checkExpect(mazeWorld.makeScene(), moveUp);

    WorldScene hidePaths = new WorldScene(1500, 800);
    mazeWorld.onKeyEvent("t");
    hidePaths.placeImageXY(maze.render(true, false, 0), 600, 400);
    hidePaths.placeImageXY(allInstructions, 1350, 400);
    t.checkExpect(mazeWorld.makeScene(), hidePaths);

    WorldScene newSameMaze = new WorldScene(1500, 800);
    // creates a new random maze, the makeScene should be same as the world scene
    mazeWorld.onKeyEvent("c");
    newSameMaze.placeImageXY(maze.render(false, false, 0), 600, 400);
    newSameMaze.placeImageXY(allInstructions, 1350, 400);

  }

  boolean testMakeScene(Tester t) {
    ArrayList<Vertex> row1 = new ArrayList<Vertex>();
    ArrayList<Vertex> row2 = new ArrayList<Vertex>();
    Vertex v1 = new Vertex(new Posn(0, 0));
    Vertex v2 = new Vertex(new Posn(1, 1));
    Vertex v3 = new Vertex(new Posn(1, 1));
    Vertex v4 = new Vertex(new Posn(1, 1));
    row1.add(v1);
    row1.add(v2);
    row2.add(v3);
    row2.add(v4);

    ArrayList<ArrayList<Vertex>> twoD = new ArrayList<ArrayList<Vertex>>();
    twoD.add(row1);
    twoD.add(row2);
    Maze maze = new Maze(twoD);
    MazeWorld mazeWorld = new MazeWorld(maze);

    WorldScene test = new WorldScene(1500, 800);
    WorldImage moves = new TextImage("Wrong Manual Moves: " + 0, 25, Color.BLUE);
    WorldImage press = new TextImage("Press", 25, Color.BLACK);
    WorldImage dfsInstructions = new TextImage("'d' for DFS solution", 13, Color.BLACK);
    WorldImage bfsInstructions = new TextImage("'b' for BFS solution", 13, Color.BLACK);
    WorldImage toggleInstructions = new TextImage("'t' to toggle hiding visited", 13, Color.BLACK);
    WorldImage createInstructions = new TextImage("'c' to create new maze", 13, Color.BLACK);
    WorldImage restartInstructions = new TextImage("'r' to restart on same maze", 13, Color.BLACK);
    WorldImage hInstructions = new TextImage("'h' to build maze w/ " + "horiz. wall bias", 13,
        Color.BLACK);
    WorldImage vInstructions = new TextImage("'v' to build maze w/ " + "vert. wall bias", 13,
        Color.BLACK);
    WorldImage kInstructions = new TextImage("'k' to animate knocking down walls", 13, Color.BLACK);
    WorldImage sInstructions = new TextImage("'s' for gradient of distance from start", 13,
        Color.BLACK);
    WorldImage eInstructions = new TextImage("'e' for gradient of distance from end", 13,
        Color.BLACK);
    WorldImage doubleClick = new TextImage("click 's' or 'e' again to revert screen", 13,
        Color.BLACK);

    WorldImage allInstructions = new AboveImage(moves, press, dfsInstructions, bfsInstructions,
        toggleInstructions, createInstructions, restartInstructions, hInstructions, vInstructions,
        kInstructions, sInstructions, eInstructions, doubleClick);
    test.placeImageXY(maze.render(false, false, 0), 600, 400);
    test.placeImageXY(allInstructions, 1350, 400);
    return t.checkExpect(mazeWorld.makeScene(), test);
  }

  // test for compareTo in Edge class
  boolean testLeastEdge(Tester t) {
    Vertex v1 = new Vertex(new Posn(1, 1));
    Vertex v2 = new Vertex(new Posn(2, 2));
    Edge heavier = new Edge(v1, v2, 5);
    Vertex v3 = new Vertex(new Posn(1, 1));
    Vertex v4 = new Vertex(new Posn(2, 2));
    Edge lighter = new Edge(v3, v4, 10);

    return t.checkExpect(heavier.compareTo(lighter), -5)
        && t.checkExpect(lighter.compareTo(lighter), 0)
        && t.checkExpect(lighter.compareTo(heavier), 5);

  }

  // test for movePlayer in vertex class
  boolean testMovePlayer(Tester t) {
    Vertex v2 = new Vertex(new Posn(2, 2));
    Vertex v3 = new Vertex(new Posn(1, 1));
    Vertex v4 = new Vertex(new Posn(2, 2));
    Vertex v5 = new Vertex(new Posn(2, 2));
    Vertex v1 = new Vertex(new Posn(1, 1), v2, v3, v4, v5);
    return t.checkExpect(v1.movePlayer("up"), v2) && t.checkExpect(v1.movePlayer("right"), v5)
        && t.checkExpect(v1.movePlayer("left"), v4) && t.checkExpect(v1.movePlayer("down"), v3)
        && t.checkExpect(v5.movePlayer("right"), v5);

  }

  // test for acceptAPlayer AND draw: should draw a red cell since after accepting
  // a player, the player is red
  boolean acceptAPlayer(Tester t) {
    Vertex visitedVertex = new Vertex(new Posn(1, 1));
    int cellWidth = Math.min(800 / 1, 1200 / 1);
    WorldImage cellPicture = new RectangleImage(cellWidth, cellWidth, OutlineMode.SOLID, Color.RED);

    cellPicture = cellPicture.movePinhole(-1 * (cellWidth / 2 - 1), 0);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(2, cellWidth, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(cellWidth / 2 - 1, 0);

    cellPicture = cellPicture.movePinhole(cellWidth / 2 - 1, 0);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(2, cellWidth, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(-1 * (cellWidth / 2 - 1), 0);

    cellPicture = cellPicture.movePinhole(0, -1 * (cellWidth / 2 - 1));
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(cellWidth, 2, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(0, cellWidth / 2 - 1);

    cellPicture = cellPicture.movePinhole(0, cellWidth / 2 - 1);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(cellWidth, 2, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(0, -1 * (cellWidth / 2 - 1));

    visitedVertex.acceptAPlayer();

    return t.checkExpect(visitedVertex.draw(1, 1, false, false, 0), cellPicture);
  }

  // test for removeAPlayer AND draw, after removing a player, the color shold be
  // (0, 0, 182, 155)
  boolean testRemoveAPlayer(Tester t) {
    Vertex visitedVertex = new Vertex(new Posn(1, 1));
    int cellWidth = Math.min(800 / 1, 1200 / 1);
    WorldImage cellPicture = new RectangleImage(cellWidth, cellWidth, OutlineMode.SOLID,
        new Color(0, 0, 182, 155));

    cellPicture = cellPicture.movePinhole(-1 * (cellWidth / 2 - 1), 0);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(2, cellWidth, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(cellWidth / 2 - 1, 0);

    cellPicture = cellPicture.movePinhole(cellWidth / 2 - 1, 0);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(2, cellWidth, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(-1 * (cellWidth / 2 - 1), 0);

    cellPicture = cellPicture.movePinhole(0, -1 * (cellWidth / 2 - 1));
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(cellWidth, 2, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(0, cellWidth / 2 - 1);

    cellPicture = cellPicture.movePinhole(0, cellWidth / 2 - 1);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(cellWidth, 2, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(0, -1 * (cellWidth / 2 - 1));

    visitedVertex.acceptAPlayer();
    visitedVertex.removeAPlayer();
    return t.checkExpect(visitedVertex.draw(1, 1, false, false, 0), cellPicture);
  }

  // test for addToProductivePath AND draw, should render Color.BLUE
  boolean testAddToProductivePath(Tester t) {
    Vertex visitedVertex = new Vertex(new Posn(10, 10));
    int cellWidth = Math.min(800 / 1, 1200 / 1);
    WorldImage cellPicture = new RectangleImage(cellWidth, cellWidth, OutlineMode.SOLID,
        Color.BLUE);

    cellPicture = cellPicture.movePinhole(-1 * (cellWidth / 2 - 1), 0);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(2, cellWidth, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(cellWidth / 2 - 1, 0);

    cellPicture = cellPicture.movePinhole(cellWidth / 2 - 1, 0);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(2, cellWidth, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(-1 * (cellWidth / 2 - 1), 0);

    cellPicture = cellPicture.movePinhole(0, -1 * (cellWidth / 2 - 1));
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(cellWidth, 2, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(0, cellWidth / 2 - 1);

    cellPicture = cellPicture.movePinhole(0, cellWidth / 2 - 1);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(cellWidth, 2, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(0, -1 * (cellWidth / 2 - 1));

    visitedVertex.addToProductivePath();

    return t.checkExpect(visitedVertex.draw(1, 1, false, true, 0), cellPicture);
  }

  // test for colorNextVisited AND draw
  boolean testColorNextVisited(Tester t) {
    Vertex visitedVertex = new Vertex(new Posn(1, 1));
    ArrayList<Vertex> row = new ArrayList<Vertex>();
    ArrayList<ArrayList<Vertex>> twoD = new ArrayList<ArrayList<Vertex>>();
    row.add(visitedVertex);
    twoD.add(row);

    Maze maze = new Maze(twoD);
    int cellWidth = Math.min(800 / 5, 1200 / 5);
    WorldImage cellPicture = new RectangleImage(cellWidth, cellWidth, OutlineMode.SOLID,
        new Color(0, 0, 182, 155));

    cellPicture = cellPicture.movePinhole(-1 * (cellWidth / 2 - 1), 0);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(2, cellWidth, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(cellWidth / 2 - 1, 0);

    cellPicture = cellPicture.movePinhole(cellWidth / 2 - 1, 0);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(2, cellWidth, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(-1 * (cellWidth / 2 - 1), 0);

    cellPicture = cellPicture.movePinhole(0, -1 * (cellWidth / 2 - 1));
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(cellWidth, 2, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(0, cellWidth / 2 - 1);

    cellPicture = cellPicture.movePinhole(0, cellWidth / 2 - 1);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(cellWidth, 2, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(0, -1 * (cellWidth / 2 - 1));

    visitedVertex.removeAPlayer();

    Maze mazeNoProgress = new Maze(twoD);
    mazeNoProgress.clearProgress();
    return t.checkExpect(maze.colorNextVisited(), true)
        && t.checkExpect(visitedVertex.draw(5, 5, false, false, 0), cellPicture)
        && t.checkExpect(mazeNoProgress.colorNextVisited(), false);
  }

  // test for other draw parts not tested by above
  boolean testDrawStart(Tester t) {
    Vertex visitedVertex = new Vertex(new Posn(0, 0));
    int cellWidth = Math.min(800 / 5, 1200 / 5);
    WorldImage cellPicture = new RectangleImage(cellWidth, cellWidth, OutlineMode.SOLID,
        Color.GREEN);

    cellPicture = cellPicture.movePinhole(-1 * (cellWidth / 2 - 1), 0);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(2, cellWidth, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(cellWidth / 2 - 1, 0);

    cellPicture = cellPicture.movePinhole(cellWidth / 2 - 1, 0);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(2, cellWidth, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(-1 * (cellWidth / 2 - 1), 0);

    cellPicture = cellPicture.movePinhole(0, -1 * (cellWidth / 2 - 1));
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(cellWidth, 2, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(0, cellWidth / 2 - 1);

    cellPicture = cellPicture.movePinhole(0, cellWidth / 2 - 1);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(cellWidth, 2, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(0, -1 * (cellWidth / 2 - 1));

    return t.checkExpect(visitedVertex.draw(5, 5, false, false, 0), cellPicture);
  }

  // test for other draw parts not tested by above
  boolean testDrawEnd(Tester t) {
    Vertex visitedVertex = new Vertex(new Posn(4, 4));
    int cellWidth = Math.min(800 / 5, 1200 / 5);
    WorldImage cellPicture = new RectangleImage(cellWidth, cellWidth, OutlineMode.SOLID,
        Color.MAGENTA);

    cellPicture = cellPicture.movePinhole(-1 * (cellWidth / 2 - 1), 0);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(2, cellWidth, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(cellWidth / 2 - 1, 0);

    cellPicture = cellPicture.movePinhole(cellWidth / 2 - 1, 0);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(2, cellWidth, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(-1 * (cellWidth / 2 - 1), 0);

    cellPicture = cellPicture.movePinhole(0, -1 * (cellWidth / 2 - 1));
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(cellWidth, 2, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(0, cellWidth / 2 - 1);

    cellPicture = cellPicture.movePinhole(0, cellWidth / 2 - 1);
    cellPicture = new OverlayOffsetAlign(AlignModeX.PINHOLE, AlignModeY.PINHOLE,
        new RectangleImage(cellWidth, 2, OutlineMode.SOLID, Color.BLACK), 0, 0, cellPicture)
        .movePinhole(0, -1 * (cellWidth / 2 - 1));

    return t.checkExpect(visitedVertex.draw(5, 5, false, false, 0), cellPicture);
  }

  // test for isComplete in Maze Class
  boolean testIsComplete(Tester t) {
    ArrayList<Vertex> row1 = new ArrayList<Vertex>();
    ArrayList<Vertex> row2 = new ArrayList<Vertex>();
    Vertex v1 = new Vertex(new Posn(0, 0));
    Vertex v2 = new Vertex(new Posn(1, 1));
    Vertex v3 = new Vertex(new Posn(1, 1));
    Vertex v4 = new Vertex(new Posn(1, 1));
    row1.add(v1);
    row1.add(v2);
    row2.add(v3);
    row2.add(v4);

    ArrayList<ArrayList<Vertex>> twoD = new ArrayList<ArrayList<Vertex>>();
    twoD.add(row1);
    twoD.add(row2);
    Maze maze = new Maze(twoD);
    maze.handleNewSpot(v2);
    maze.handleNewSpot(v4);

    ArrayList<Vertex> row1Incomplete = new ArrayList<Vertex>();
    ArrayList<Vertex> row2Incomplete = new ArrayList<Vertex>();
    Vertex v1I = new Vertex(new Posn(0, 0));
    Vertex v2I = new Vertex(new Posn(1, 1));
    Vertex v3I = new Vertex(new Posn(1, 1));
    Vertex v4I = new Vertex(new Posn(1, 1));
    row1Incomplete.add(v1I);
    row1Incomplete.add(v2I);
    row2Incomplete.add(v3I);
    row2Incomplete.add(v4I);

    ArrayList<ArrayList<Vertex>> twoDIncomplete = new ArrayList<ArrayList<Vertex>>();
    twoDIncomplete.add(row1Incomplete);
    twoDIncomplete.add(row2Incomplete);
    Maze mazeIncomplete = new Maze(twoDIncomplete);
    return t.checkExpect(maze.isCompleted(), true)
        && t.checkExpect(mazeIncomplete.isCompleted(), false);
  }

  // test for onTick
  void testOnTick(Tester t) {
    ArrayList<Vertex> row1 = new ArrayList<Vertex>();
    ArrayList<Vertex> row2 = new ArrayList<Vertex>();
    Vertex v1 = new Vertex(new Posn(0, 0));
    Vertex v2 = new Vertex(new Posn(1, 1));
    Vertex v3 = new Vertex(new Posn(1, 1));
    Vertex v4 = new Vertex(new Posn(1, 1));
    row1.add(v1);
    row1.add(v2);
    row2.add(v3);
    row2.add(v4);

    ArrayList<ArrayList<Vertex>> twoD = new ArrayList<ArrayList<Vertex>>();
    twoD.add(row1);
    twoD.add(row2);
    Maze maze = new Maze(twoD);
    MazeWorld mazeWorld = new MazeWorld(maze);

    WorldScene test = new WorldScene(1500, 800);
    WorldImage moves = new TextImage("Wrong Manual Moves: " + 0, 25, Color.BLUE);
    WorldImage press = new TextImage("Press", 25, Color.BLACK);
    WorldImage dfsInstructions = new TextImage("'d' for DFS solution", 13, Color.BLACK);
    WorldImage bfsInstructions = new TextImage("'b' for BFS solution", 13, Color.BLACK);
    WorldImage toggleInstructions = new TextImage("'t' to toggle hiding visited", 13, Color.BLACK);
    WorldImage createInstructions = new TextImage("'c' to create new maze", 13, Color.BLACK);
    WorldImage restartInstructions = new TextImage("'r' to restart on same maze", 13, Color.BLACK);
    WorldImage hInstructions = new TextImage("'h' to build maze w/ " + "horiz. wall bias", 13,
        Color.BLACK);
    WorldImage vInstructions = new TextImage("'v' to build maze w/ " + "vert. wall bias", 13,
        Color.BLACK);
    WorldImage kInstructions = new TextImage("'k' to animate knocking down walls", 13, Color.BLACK);
    WorldImage sInstructions = new TextImage("'s' for gradient of distance from start", 13,
        Color.BLACK);
    WorldImage eInstructions = new TextImage("'e' for gradient of distance from end", 13,
        Color.BLACK);
    WorldImage doubleClick = new TextImage("click 's' or 'e' again to revert screen", 13,
        Color.BLACK);

    WorldImage allInstructions = new AboveImage(moves, press, dfsInstructions, bfsInstructions,
        toggleInstructions, createInstructions, restartInstructions, hInstructions, vInstructions,
        kInstructions, sInstructions, eInstructions, doubleClick);
    test.placeImageXY(maze.render(false, false, 0), 600, 400);
    test.placeImageXY(allInstructions, 1350, 400);

    t.checkExpect(mazeWorld.makeScene(), test);
    mazeWorld.onTick();
    mazeWorld.onTick();
    mazeWorld.onTick();
    t.checkExpect(mazeWorld.makeScene(), test);
  }

  // test for reset in Vertex class
  // after changing the state of the vOriginal, it is equal to the exact same
  // after resetting
  boolean testReset(Tester t) {
    Vertex v = new Vertex(new Posn(1, 1));
    Vertex vOriginal = new Vertex(new Posn(1, 1));

    v.acceptAPlayer();
    v.reset();
    return t.checkExpect(v, vOriginal);
  }

  // test for create connection
  // the created vertices initially have no connections, so after creating
  // connections, to
  // verify, we can use the movePlayer method to test the neighbors exist,
  // (returns the correct neighbor),
  // or doesn't exist (return the same vertex)
  boolean testCreateConnection(Tester t) {
    Vertex vCenter = new Vertex(new Posn(1, 1));
    Vertex vCenterNoConnection = new Vertex(new Posn(1, 1));
    Vertex vLeft = new Vertex(new Posn(1, 0));

    Vertex vRight = new Vertex(new Posn(1, 2));

    Vertex vDown = new Vertex(new Posn(2, 1));

    vCenter.createConnection(vLeft);
    vLeft.createConnection(vCenter);
    vCenter.createConnection(vRight);
    vRight.createConnection(vCenter);
    vDown.createConnection(vCenter);
    vCenter.createConnection(vDown);

    return t.checkExpect(vCenterNoConnection.movePlayer("left"), vCenterNoConnection)
        && t.checkExpect(vCenter.movePlayer("left"), vLeft)
        && t.checkExpect(vCenter.movePlayer("right"), vRight)
        && t.checkExpect(vCenter.movePlayer("down"), vDown);

  }

  // test for createConnectionHelper, same testing logic as above
  boolean testCreateConnectionHelp(Tester t) {
    Vertex vCenter = new Vertex(new Posn(1, 1));
    Vertex vCenterNoConnection = new Vertex(new Posn(1, 1));
    Vertex vLeft = new Vertex(new Posn(1, 0));
    Vertex vRight = new Vertex(new Posn(1, 2));
    Vertex vDown = new Vertex(new Posn(2, 1));

    vCenter.createConnectionHelper(vLeft, new Posn(2, 0));
    vCenter.createConnectionHelper(vRight, new Posn(2, 3));
    vCenter.createConnectionHelper(vDown, new Posn(2, 1));

    return t.checkExpect(vCenterNoConnection.movePlayer("left"), vCenterNoConnection)
        && t.checkExpect(vCenter.movePlayer("left"), vLeft)
        && t.checkExpect(vCenter.movePlayer("right"), vRight)
        && t.checkExpect(vCenter.movePlayer("up"), vCenter)
        && t.checkExpect(vCenter.movePlayer("down"), vDown);
  }

  // test for movePlayer in Maze Class
  boolean testMovePlayerMaze(Tester t) {
    ArrayList<Vertex> vertices = new ArrayList<Vertex>();
    Vertex v2 = new Vertex(new Posn(2, 2));
    Vertex v3 = new Vertex(new Posn(1, 1));

    vertices.add(v2);
    vertices.add(v3);

    ArrayList<ArrayList<Vertex>> twoD = new ArrayList<ArrayList<Vertex>>();
    twoD.add(vertices);
    ArrayList<Vertex> vertices2 = new ArrayList<Vertex>();
    Vertex v4 = new Vertex(new Posn(3, 3));
    Vertex v5 = new Vertex(new Posn(3, 3));
    vertices2.add(v4);
    vertices2.add(v4);
    twoD.add(vertices2);
    Maze maze = new Maze(twoD);

    v2.createConnection(v3);
    v3.createConnection(v2);
    v4.createConnection(v2);
    v2.createConnection(v4);
    // make it so we visited v3 already
    maze.handleNewSpot(v3);

    return t.checkExpect(maze.movePlayer("right"), 1) && t.checkExpect(maze.movePlayer("left"), 0)
        && t.checkExpect(maze.movePlayer("down"), 0) && t.checkExpect(maze.movePlayer("up"), 0);
  }

  // test for handleNewSpot
  boolean testHandleNewSpot(Tester t) {
    ArrayList<Vertex> vertices = new ArrayList<Vertex>();
    Vertex v2 = new Vertex(new Posn(2, 2));
    Vertex v3 = new Vertex(new Posn(1, 1));
    vertices.add(v2);
    vertices.add(v3);
    ArrayList<ArrayList<Vertex>> twoD = new ArrayList<ArrayList<Vertex>>();
    twoD.add(vertices);
    Maze maze = new Maze(twoD);

    return t.checkExpect(maze.handleNewSpot(v3), 0) 
        && t.checkExpect(maze.handleNewSpot(v2), 1);
  }

  // test for innitializeEdgeBetweenEveryCell using seeded method for randomness
  boolean testGenerateRandomMaze(Tester t) {
    Random rand = new Random(40);
    MazeBuilder maze = new MazeBuilder(2, 2);
    MazeBuilder maze2 = new MazeBuilder(2, 2);
    int seed = rand.nextInt();

    return t.checkExpect(maze.generateRandomMazeSeed(0, seed),
        maze2.generateRandomMazeSeed(0, seed))
        && t.checkExpect(maze2.generateRandomMazeSeed(1, seed), 
            maze.generateRandomMazeSeed(1, seed))
        && t.checkExpect(maze2.generateRandomMazeSeed(-1, seed), 
            maze.generateRandomMazeSeed(-1, seed));
  }

  // after pressing S and then reseting variables, the makeScene should
  // render the same thing
  void testResetVariable(Tester t) {
    ArrayList<Vertex> row1 = new ArrayList<Vertex>();
    ArrayList<Vertex> row2 = new ArrayList<Vertex>();
    Vertex v1 = new Vertex(new Posn(0, 0));
    Vertex v2 = new Vertex(new Posn(1, 1));
    Vertex v3 = new Vertex(new Posn(1, 1));
    Vertex v4 = new Vertex(new Posn(1, 1));
    row1.add(v1);
    row1.add(v2);
    row2.add(v3);
    row2.add(v4);

    ArrayList<ArrayList<Vertex>> twoD = new ArrayList<ArrayList<Vertex>>();
    twoD.add(row1);
    twoD.add(row2);
    Maze maze = new Maze(twoD);
    MazeWorld mazeWorld = new MazeWorld(maze);

    WorldScene test = new WorldScene(1500, 800);
    WorldImage moves = new TextImage("Wrong Manual Moves: " + 0, 25, Color.BLUE);
    WorldImage press = new TextImage("Press", 25, Color.BLACK);
    WorldImage dfsInstructions = new TextImage("'d' for DFS solution", 13, Color.BLACK);
    WorldImage bfsInstructions = new TextImage("'b' for BFS solution", 13, Color.BLACK);
    WorldImage toggleInstructions = new TextImage("'t' to toggle hiding visited", 13, Color.BLACK);
    WorldImage createInstructions = new TextImage("'c' to create new maze", 13, Color.BLACK);
    WorldImage restartInstructions = new TextImage("'r' to restart on same maze", 13, Color.BLACK);
    WorldImage hInstructions = new TextImage("'h' to build maze w/ " + "horiz. wall bias", 13,
        Color.BLACK);
    WorldImage vInstructions = new TextImage("'v' to build maze w/ " + "vert. wall bias", 13,
        Color.BLACK);
    WorldImage kInstructions = 
        new TextImage("'k' to animate knocking down walls", 13, Color.BLACK);
    WorldImage sInstructions = new TextImage("'s' for gradient of distance from start", 13,
        Color.BLACK);
    WorldImage eInstructions = new TextImage("'e' for gradient of distance from end", 13,
        Color.BLACK);
    WorldImage doubleClick = new TextImage("click 's' or 'e' again to revert screen", 13,
        Color.BLACK);

    WorldImage allInstructions = new AboveImage(moves, press, dfsInstructions, bfsInstructions,
        toggleInstructions, createInstructions, restartInstructions, hInstructions, vInstructions,
        kInstructions, sInstructions, eInstructions, doubleClick);
    test.placeImageXY(maze.render(false, false, 0), 600, 400);
    test.placeImageXY(allInstructions, 1350, 400);

    ArrayList<Vertex> row1G = new ArrayList<Vertex>();
    ArrayList<Vertex> row2G = new ArrayList<Vertex>();
    Vertex v1G = new Vertex(new Posn(0, 0));
    Vertex v2G = new Vertex(new Posn(1, 1));
    Vertex v3G = new Vertex(new Posn(1, 1));
    Vertex v4G = new Vertex(new Posn(1, 1));
    row1G.add(v1G);
    row1G.add(v2G);
    row2G.add(v3G);
    row2G.add(v4G);

    ArrayList<ArrayList<Vertex>> twoDG = new ArrayList<ArrayList<Vertex>>();
    twoDG.add(row1G);
    twoDG.add(row2G);
    Maze mazeG = new Maze(twoDG);

    WorldScene testG = new WorldScene(1500, 800);
    testG.placeImageXY(mazeG.render(false, false, 1), 600, 400);
    testG.placeImageXY(allInstructions, 1350, 400);

    mazeWorld.onKeyEvent("s");
    mazeWorld.resetVariables();
    t.checkExpect(mazeWorld.makeScene(), test);

  }

  // test for showAnotherElementOfProductivePath
  boolean testShowAnotherElement(Tester t) {
    ArrayList<Vertex> row1 = new ArrayList<Vertex>();

    Vertex v1 = new Vertex(new Posn(0, 0));
    row1.add(v1);
    ArrayList<ArrayList<Vertex>> twoD = new ArrayList<ArrayList<Vertex>>();
    twoD.add(row1);
    Maze maze = new Maze(twoD);

    // removes the last productive path, so the next one will be false since already
    // visited
    return t.checkExpect(maze.showAnotherElementOfProductivePath(), true)
        && t.checkExpect(maze.showAnotherElementOfProductivePath(), false);
  }

  // test for clearProgress
  boolean testClearProgress(Tester t) {
    ArrayList<Vertex> row1 = new ArrayList<Vertex>();
    ArrayList<Vertex> row2 = new ArrayList<Vertex>();
    Vertex v1 = new Vertex(new Posn(0, 0));
    Vertex v2 = new Vertex(new Posn(1, 1));
    Vertex v3 = new Vertex(new Posn(1, 1));
    Vertex v4 = new Vertex(new Posn(1, 1));
    row1.add(v1);
    row1.add(v2);
    row2.add(v3);
    row2.add(v4);
    ArrayList<ArrayList<Vertex>> twoD = new ArrayList<ArrayList<Vertex>>();
    twoD.add(row1);
    twoD.add(row2);
    Maze maze = new Maze(twoD);
    maze.handleNewSpot(v2);
    maze.handleNewSpot(v4);

    Maze mazeIncomplete = new Maze(twoD);
    mazeIncomplete.handleNewSpot(v2);
    mazeIncomplete.handleNewSpot(v4);

    // after clearing progress, the mazeIncomplete went from complete -> incomplete
    mazeIncomplete.clearProgress();

    return t.checkExpect(maze.isCompleted(), true)
        && t.checkExpect(mazeIncomplete.isCompleted(), false);

  }

  // test for compjteGradient
  // after computing gradient, the colors of the cells should change colors
  // and draw different vertices
  boolean testComputeGradient(Tester t) {
    ArrayList<Vertex> row1 = new ArrayList<Vertex>();
    ArrayList<Vertex> row2 = new ArrayList<Vertex>();
    Vertex v1 = new Vertex(new Posn(0, 0));
    Vertex v2 = new Vertex(new Posn(1, 1));
    Vertex v3 = new Vertex(new Posn(1, 1));
    Vertex v4 = new Vertex(new Posn(1, 1));
    row1.add(v1);
    row1.add(v2);
    row2.add(v3);
    row2.add(v4);
    ArrayList<ArrayList<Vertex>> twoD = new ArrayList<ArrayList<Vertex>>();
    twoD.add(row1);
    twoD.add(row2);
    Maze maze = new Maze(twoD);

    maze.computeGradient(true);

    WorldImage all = new CircleImage(0, OutlineMode.SOLID, Color.BLACK);

    WorldImage row1d = new CircleImage(0, OutlineMode.SOLID, Color.BLACK);

    Vertex curr = v1;
    row1d = new BesideImage(row1d, curr.draw(2, 2, false, false, 1));

    Vertex curr2 = v2;
    row1d = new BesideImage(row1d, curr2.draw(2, 2, false, false, 1));

    all = new AboveImage(all, row1d);

    Vertex curr3 = v3;
    Vertex curr4 = v4;

    WorldImage row2d = new CircleImage(0, OutlineMode.SOLID, Color.BLACK);
    row2d = new BesideImage(row2d, curr3.draw(2, 2, false, false, 1));
    row2d = new BesideImage(row2d, curr4.draw(2, 2, false, false, 1));

    all = new AboveImage(all, row2d);

    return t.checkExpect(maze.render(false, false, 1), all);

  }

  // test for render in Maze class
  boolean testRender(Tester t) {
    ArrayList<Vertex> row1 = new ArrayList<Vertex>();
    ArrayList<Vertex> row2 = new ArrayList<Vertex>();
    Vertex v1 = new Vertex(new Posn(0, 0));
    Vertex v2 = new Vertex(new Posn(1, 1));
    Vertex v3 = new Vertex(new Posn(1, 1));
    Vertex v4 = new Vertex(new Posn(1, 1));
    row1.add(v1);
    row1.add(v2);
    row2.add(v3);
    row2.add(v4);

    ArrayList<ArrayList<Vertex>> twoD = new ArrayList<ArrayList<Vertex>>();
    twoD.add(row1);
    twoD.add(row2);
    Maze maze = new Maze(twoD);

    WorldImage all = new CircleImage(0, OutlineMode.SOLID, Color.BLACK);

    WorldImage row1d = new CircleImage(0, OutlineMode.SOLID, Color.BLACK);

    Vertex curr = v1;
    row1d = new BesideImage(row1d, curr.draw(2, 2, false, false, 0));

    Vertex curr2 = v2;
    row1d = new BesideImage(row1d, curr2.draw(2, 2, false, false, 0));

    all = new AboveImage(all, row1d);

    Vertex curr3 = v3;
    Vertex curr4 = v4;

    WorldImage row2d = new CircleImage(0, OutlineMode.SOLID, Color.BLACK);
    row2d = new BesideImage(row2d, curr3.draw(2, 2, false, false, 0));
    row2d = new BesideImage(row2d, curr4.draw(2, 2, false, false, 0));

    all = new AboveImage(all, row2d);

    return t.checkExpect(maze.render(false, false, 0), all);
  }

  // test for fullyWalled
  // if everything is fully walled, then finding a path will return false for
  // bfs/dfs since no way to the end
  boolean testGenerateFullyWalledMaze(Tester t) {
    MazeBuilder mazee = new MazeBuilder(2, 2);
    Maze res = mazee.generateFullyWalledMaze();
    return t.checkExpect(res.findPath(true), false) && t.checkExpect(res.findPath(false), false);
  }

  // test knockDownWall
  boolean testKnockDown(Tester t) {
    MazeBuilder mazee = new MazeBuilder(2, 2);
    mazee.generateRandomMaze(0);
    // can only have 3 walls, since n-1, so the fourth call should return false
    // (no walls to be removed)
    return t.checkExpect(mazee.knockDownWall(), true) && t.checkExpect(mazee.knockDownWall(), true)
        && t.checkExpect(mazee.knockDownWall(), true)
        && t.checkExpect(mazee.knockDownWall(), false);

  }

  // test for findPath
  // if there is no connections, there exists no path, so return false
  // else, if a path is found -> return true
  void testFindPath(Tester t) {
    ArrayList<Vertex> row1 = new ArrayList<Vertex>();
    ArrayList<Vertex> row2 = new ArrayList<Vertex>();
    Vertex v1 = new Vertex(new Posn(0, 0));
    Vertex v2 = new Vertex(new Posn(1, 1));
    Vertex v3 = new Vertex(new Posn(1, 1));
    Vertex v4 = new Vertex(new Posn(1, 1));
    row1.add(v1);
    row1.add(v2);
    row2.add(v3);
    row2.add(v4);

    ArrayList<ArrayList<Vertex>> twoD = new ArrayList<ArrayList<Vertex>>();
    twoD.add(row1);
    twoD.add(row2);
    Maze maze = new Maze(twoD);

    MazeBuilder randomMaze = new MazeBuilder(2, 2);
    Maze res = randomMaze.generateRandomMaze(0);

    t.checkExpect(maze.findPath(true), false);
    t.checkExpect(res.findPath(false), true);
  }

  // test for innitializeEdgeBetweenEveryCell
  void testInitializeEdgeBetweenEveryCell(Tester t) {
    Random rand = new Random(40);
    MazeBuilder maze = new MazeBuilder(2, 2);
    MazeBuilder maze2 = new MazeBuilder(2, 2);
    int seed = rand.nextInt();
    maze.initializeEdgeBetweenEveryCellSeed(0, seed);
    maze2.initializeEdgeBetweenEveryCellSeed(0, seed);
    t.checkExpect(maze.getAllPossible(), maze2.getAllPossible());
    
    maze.initializeEdgeBetweenEveryCellSeed(-1, seed);
    maze2.initializeEdgeBetweenEveryCellSeed(-1, seed);
    t.checkExpect(maze.getAllPossible(), maze2.getAllPossible());
    
    maze.initializeEdgeBetweenEveryCellSeed(1, seed);
    maze2.initializeEdgeBetweenEveryCellSeed(1, seed);
    t.checkExpect(maze.getAllPossible(), maze2.getAllPossible());
  }

  // test for kruskal algorithm
  boolean testKruskal(Tester t) {
    Random rand = new Random(40);
    int seed = rand.nextInt();
    MazeBuilder maze = new MazeBuilder(2, 2);
    maze.initializeEdgeBetweenEveryCellSeed(0, seed);

    Vertex v3 = new Vertex(new Posn(0, 1));
    Vertex v4 = new Vertex(new Posn(0, 0));
    Vertex v8 = new Vertex(new Posn(1, 1));
    Vertex v6 = new Vertex(new Posn(1, 0));

    // vertex:3 ->

    Edge ed1 = new Edge(v3, v4, -1123336207);
    Edge ed2 = new Edge(v6, v4, -1123336207);
    Edge ed3 = new Edge(v8, v3, -1123336207);
    ArrayList<Edge> finalEdges1 = maze.performKruskal();

    ArrayList<Edge> customEdges = new ArrayList<Edge>();
    customEdges.add(ed1);
    customEdges.add(ed2);
    customEdges.add(ed3);

    return t.checkExpect(finalEdges1.size(), 3) && t.checkExpect(finalEdges1, customEdges);
  }

  // test for constructor exceptions with varying maze sizes
  boolean testExceptions(Tester t) {
    return t.checkConstructorException(
        new IllegalArgumentException("Must have at least 1 row & column"), "MazeWorld", 0, 0)
        && t.checkConstructorException(
            new IllegalArgumentException(
                "Too big! Choose a " + "maze with at most 100 rows and 60 columns"),
            "MazeWorld", 101, 61);
  }
  
  // Test for updateDistanceFrom method in the vertex class
  // Tests initializing distance for the home vertex, along with changing for neighbors
  // Tests these for start and end cases
  void testUpdateDistanceFrom(Tester t) {
    Vertex exBorder = new Vertex(new Posn(1, 2));
    Vertex zeroDist = new Vertex(false, false, false, 0, 0, 
        new Posn(1, 1), exBorder, exBorder, exBorder, exBorder);
    Vertex bigInitialized = new Vertex(false, false, false, 10000, 10000, 
        new Posn(1, 1), exBorder, exBorder, exBorder, exBorder);
    Vertex initial = new Vertex(false, false, false, 10000, 0, 
        new Posn(1, 1), exBorder, exBorder, exBorder, exBorder);
    Vertex updated = new Vertex(false, false, false, 1, 10000, 
        new Posn(1, 1), exBorder, exBorder, exBorder, exBorder);
    
    Vertex updatedEnd = new Vertex(false, false, false, 1, 1, 
        new Posn(1, 1), exBorder, exBorder, exBorder, exBorder);
    
    initial.updateDistanceFrom(true, true, initial);
    
    t.checkExpect(initial, zeroDist);
    
    bigInitialized.updateDistanceFrom(true, false, initial);
    
    t.checkExpect(bigInitialized, updated);
    
    bigInitialized.updateDistanceFrom(false, false, initial);
    
    t.checkExpect(bigInitialized, updatedEnd);
    
  }

}
