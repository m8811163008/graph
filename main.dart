import 'package:algo/queue/queue_double_stack.dart';
import 'package:algo/stack/stack.dart';

class Vertex<T> {
  final int index;
  final T data;

  const Vertex({
    required this.index,
    required this.data,
  });

  @override
  String toString() {
    return 'Vertex{index: $index, data: $data}';
  }
}

class Edge<T> {
  final Vertex<T> source;
  final Vertex<T> destination;
  final double? weight;

  Edge(this.source, this.destination, [this.weight]);
}

enum EdgeType { directed, undirected }

abstract class Graph<E> {
  /// Returns all of the vertices in a graph.
  Iterable<Vertex<E>> get vertices;

  /// Create a vertex and adds it to the graph.
  Vertex<E> createVertex(E data);

  /// Connect two vertices in the graph with either a directed
  /// or undirected edge.The weight is optional.
  void addEdge(Vertex<E> source, Vertex<E> destination,
      {EdgeType edgeType, double? weight});

  /// Returns a list of *outgoing* edges from a specific vertex.
  List<Edge<E>> edges(Vertex<E> source);

  /// Returns the weight of the edge between two vertices.
  double? weight(Vertex<E> source, Vertex<E> destination);
}

class AdjacencyList<E> implements Graph<E> {
  // Data list stores values of a vertex and the corresponding
  // list of edges.
  final Map<Vertex<E>, List<Edge<E>>> _connections = {};
  // Assigns a unique index to each new vertex.
  var _nextIndex = 0;

  @override
  Iterable<Vertex<E>> get vertices => _connections.keys;
  @override
  Vertex<E> createVertex(E data) {
    // create a new vertex with a unique index.
    final vertex = Vertex(index: _nextIndex, data: data);
    _nextIndex++;
    // Then add the vertex as a key in the _connections map.The
    // outgoing edges is empty.
    _connections[vertex] = [];
    return vertex;
  }

  @override
  void addEdge(
    Vertex<E> source,
    Vertex<E> destination, {
    EdgeType edgeType = EdgeType.undirected,
    double? weight,
  }) {
    // Add edge from source to destination.
    // Checks if it exists in the _connections map.If it does
    // creates a new directed edge from the source to the
    // destination. Then add it to the vertex's list of edges.
    _connections[source]?.add(
      Edge(source, destination, weight),
    );

    // If edgeType is undirected add Edge from destination to
    // source. Add an edge going the other direction.
    if (edgeType == EdgeType.undirected) {
      _connections[destination]?.add(
        Edge(destination, source, weight),
      );
    }
  }

  /// Gets the stored outgoing edges for the provided vertex.
  @override
  List<Edge<E>> edges(Vertex<E> source) {
    return _connections[source] ?? [];
  }

  @override
  double? weight(Vertex<E> source, Vertex<E> destination) {
    final match =
        _connections[source]?.where((edge) => edge.destination == destination);
    if (match == null) return null;
    return match.first.weight;
  }

  @override
  String toString() {
    final result = StringBuffer();
    // Loop throw every key value pair in _connections.
    _connections.forEach((vertex, edges) {
      // For every vertex, find all destinations and join
      // them into a single, comma-separated string.
      final destinations = edges.map((edge) => edge.destination).join(', ');
      // Put each vertex and destinations in a new line.
      result.writeln('$vertex --> $destinations');
    });
    return result.toString();
  }
}

void main() {
  final graph = AdjacencyList<String>();
  final a = graph.createVertex('A');
  final b = graph.createVertex('B');
  final c = graph.createVertex('C');
  final e = graph.createVertex('E');
  final f = graph.createVertex('F');
  graph.addEdge(a, e, edgeType: EdgeType.directed);
  graph.addEdge(e, f, edgeType: EdgeType.directed);
  graph.addEdge(a, b, edgeType: EdgeType.directed);
  graph.addEdge(b, c, edgeType: EdgeType.directed);
  graph.addEdge(c, a, edgeType: EdgeType.directed);
  print(graph);
  print(graph.isCyclic2(a));
}

/// Each vertex has its own row and column in the table.
/// The cell where rows and column are intersect hold the edge
/// weights. It any particular cell is empty, that is, if the
/// weight is null, then that means there is no edge between the
/// row vertex and the column vertex.
class AdjacencyMatrix<E> implements Graph<E> {
  final List<Vertex<E>> _vertices = [];
  // represents network in matrix from by giving each *Vertex*
  // a row and column in a table. Edge that don't exists between
  // two cities are shown with a weight of 0 in the cells where
  // the rows and columns intersect.
  // Every rows represents a source and each column represents a
  // destination.
  final List<List<double?>?> _weights = [];
  var _nextIndex = 0;
  @override
  Iterable<Vertex<E>> get vertices => _vertices;
  @override
  Vertex<E> createVertex(E data) {
    // Add a new vertex to the list.
    final vertex = Vertex(index: _nextIndex, data: data);
    _nextIndex++;
    _vertices.add(vertex);
    // Append a null value at the end of every row.This in effect
    // creates a new destination column in the matrix.
    for (var i = 0; i < _weights.length; i++) {
      _weights[i]?.add(null);
    }
    // 3
    final row = List.filled(_vertices.length, null, growable: true);
    _weights.add(row);
    return vertex;
  }

  @override
  void addEdge(Vertex<E> source, Vertex<E> destination,
      {EdgeType edgeType = EdgeType.undirected, double? weight}) {
    // Always add a directed edge.
    _weights[source.index]?[destination.index] = weight;
    // If edge is undirected then also add another edge
    // going from the destination to the source.
    if (edgeType == EdgeType.undirected) {
      _weights[destination.index]?[source.index] = weight;
    }
  }

  @override
  List<Edge<E>> edges(Vertex<E> source) {
    List<Edge<E>> edges = [];
    // 1
    for (var column = 0; column < _weights.length; column++) {
      final weight = _weights[source.index]?[column];
      // 2
      if (weight == null) continue;
      // 3
      final destination = _vertices[column];
      edges.add(Edge(source, destination, weight));
    }
    return edges;
  }

  @override
  double? weight(Vertex<E> source, Vertex<E> destination) {
    return _weights[source.index]?[destination.index];
  }

  @override
  String toString() {
    final output = StringBuffer();
    // 1
    for (final vertex in _vertices) {
      output.writeln('${vertex.index}: ${vertex.data}');
    }
    // 2
    for (int i = 0; i < _weights.length; i++) {
      for (int j = 0; j < _weights.length; j++) {
        final value = (_weights[i]?[j] ?? '.').toString();
        output.write(value.padRight(6));
      }
      output.writeln();
    }
    return output.toString();
  }
}

extension BreadthFirstSearch<E> on Graph<E> {
  /// Takes starting vertex and traverse breadth first.
  /// When traversing a graph using BreadFirstSearch, each
  /// vertex is enqueued once. This process has a time
  /// complexity of O(V). During this traversal, you also visit
  /// all the edges.
  /// The time it takes to visit all edges is O(E). Adding the
  /// two together means that the overall time complexity for
  /// BFS is O(V + E)
  /// The space complexity of BSF is O(V) since you have to
  /// store the vertices in three separate structures: queue,
  /// enqueued and visited.
  List<Vertex<E>> breadthFirstTraverse(Vertex<E> source) {
    // queue keeps track of the neighboring vertices to visit next.
    final queue = QueueStack<Vertex<E>>();
    // this is a list that stores the order in which the vertices
    // were explored.
    List<Vertex<E>> visited = [];
    // Keeps track of which vertices have been enqueued before.
    // We use [Set] so that lookup is cheap and only takes O(1)
    // and and a list require O(n).
    Set<Vertex<E>> enqueued = {};
    // a queue will help you keep track of which vertices to
    // visit next. The first in first out approach guarantees
    // that all og a vertex's neighbors are visited before you
    // traverse one level deeper.
    // initialize the breadth first algorithm with source.
    queue.enqueue(source);
    enqueued.add(source);
    // You continue to dequeue a vertex from the queue until the
    // queue is empty.
    while (true) {
      final vertex = queue.dequeue();
      if (vertex == null) break;
      // Every time you dequeue a vertex you add it to visited list.
      visited.add(vertex);
      // Find the edges and iterate over them.
      final neighborEdges = edges(vertex);
      for (var edge in neighborEdges) {
        // For each edge, you check to see if its destination
        // vertex has been enqueued before, and, if not, you add
        // it to the queue.
        if (!enqueued.contains(edge.destination)) {
          queue.enqueue(edge.destination);
          enqueued.add(edge.destination);
        }
      }
    }
    return visited;
  }

  List<Vertex<E>> breadthFirstTraverseRecursive(Vertex<E> source) {
    QueueStack<Vertex<E>> queue = QueueStack();
    List<Vertex<E>> visited = [];
    Set<Vertex<E>> enqueuedVertices = {};
    queue.enqueue(source);
    enqueuedVertices.add(source);
    visited.add(source);
    return _breadthFirstTraverseRecursive(
        queue: queue, visited: visited, enqueuedVertices: enqueuedVertices);
  }

  List<Vertex<E>> _breadthFirstTraverseRecursive(
      {required QueueStack<Vertex<E>> queue,
      required List<Vertex<E>> visited,
      required Set<Vertex<E>> enqueuedVertices}) {
    //Base case.
    final dequeuedVertex = queue.dequeue();
    if (dequeuedVertex == null) {
      return [...visited];
    }
    final outgoingEdges = edges(dequeuedVertex);
    for (var edge in outgoingEdges) {
      if (!enqueuedVertices.contains(edge.destination)) {
        queue.enqueue(edge.destination);
        enqueuedVertices.add(edge.destination);
        visited.add(edge.destination);
      }
    }
    return _breadthFirstTraverseRecursive(
        queue: queue, visited: visited, enqueuedVertices: enqueuedVertices);
  }

  int maximumQueueSize(Vertex<E> source) {
    var maxQueueSize = 0;
    final queue = QueueStack<Vertex<E>>();
    queue.enqueue(source);
    maxQueueSize = queue.length;

    Set<Vertex<E>> enqueued = {};
    enqueued.add(source);

    while (true) {
      if (maxQueueSize < queue.length) {
        maxQueueSize = queue.length;
      }
      final dequeueVertex = queue.dequeue();
      if (dequeueVertex == null) break;
      final neighborVertices = edges(dequeueVertex);
      for (var edge in neighborVertices) {
        if (!enqueued.contains(edge.destination)) {
          queue.enqueue(edge.destination);
          enqueued.add(edge.destination);
        }
      }
    }

    return maxQueueSize;
  }

  List<Vertex<E>> breadthFirstTraverseRecursiveV2(Vertex<E> source) {
    QueueStack<Vertex<E>> queue = QueueStack();
    List<Vertex<E>> visited = [];
    Set<Vertex<E>> enqueuedVertices = {};
    queue.enqueue(source);
    enqueuedVertices.add(source);
    visited.add(source);
    _breadthFirstTraverseRecursiveV2(
        queue: queue, visited: visited, enqueuedVertices: enqueuedVertices);
    return visited;
  }

  void _breadthFirstTraverseRecursiveV2(
      {required QueueStack<Vertex<E>> queue,
      required List<Vertex<E>> visited,
      required Set<Vertex<E>> enqueuedVertices}) {
    // Base case.
    final dequeuedVertex = queue.dequeue();
    //1
    if (dequeuedVertex == null) {
      return;
    }
    //2. Mark the vertex as visited.
    visited.add(dequeuedVertex);
    final outgoingEdges = edges(dequeuedVertex);
    //3. For every edge of the current vertex, check to see if
    // the adjacent vertices have
    // been visited before inserting them into the queue.
    for (var edge in outgoingEdges) {
      if (!enqueuedVertices.contains(edge.destination)) {
        queue.enqueue(edge.destination);
        enqueuedVertices.add(edge.destination);
      }
    }
    // Recursively call this function until the queue is empty.
    _breadthFirstTraverseRecursive(
        queue: queue, visited: visited, enqueuedVertices: enqueuedVertices);
  }

  bool get isConnected {
    // always check for edge cases.
    // If there are no vertices, treat the graph as connected.
    if (vertices.isEmpty) return true;
    //2
    final traverseResult = breadthFirstTraverseRecursive(vertices.first);
    // Go through every vertex in the graph and check if it has
    // been visited before.

    // return listEquals(traverseResult, vertices.toList());
    for (var vertex in vertices) {
      if (!traverseResult.contains(vertex)) {
        return false;
      }
    }
    return true;
  }
}

extension DepthFirstSearch<E> on Graph<E> {
  /// DFS will visit every single vertex at least once.This
  /// process has time complexity of O(v).
  /// When traversing a graph in DFS you have to check all
  /// neighboring vertices
  /// to find one available to visit next.The time complexity of
  /// this is O(E) because you have to visit every edge in graph
  /// in the worst case.
  /// Overall, the time complexity of depth-first search is O(V
  /// + E)
  /// Space complexity is O(V).
  List<Vertex<E>> depthFirstSearch(Vertex<E> source) {
    // Stack is used to store your path through the graph.
    final stack = Stack<Vertex<E>>();
    // Pushed ensure you don not visit a vertex twice and using
    // Set ensures fast O(1) lookup.
    final pushed = <Vertex<E>>{};
    final visited = <Vertex<E>>[];
    stack.push(source);
    pushed.add(source);
    visited.add(source);
    // Continue to check the stack until is is empty.
    // You have labeled this loop outerLoop so that you have a
    // to continue to the next vertex, even from within a nested
    // for loop.
    outerLoop:
    while (stack.isNotEmpty) {
      final vertex = stack.peek;
      // find all the neighboring edges for the current vertex.
      final outgoingEdges = edges(vertex);
      // You loop through every outgoing edges and check if the
      // neighboring vertex has been seen. If not, you push it
      // onto the stack and add it to the visited list.
      for (var edge in outgoingEdges) {
        if (!pushed.contains(edge.destination)) {
          stack.push(edge.destination);
          pushed.add(edge.destination);
          visited.add(edge.destination);
          // Now that you've found a neighbor to visit you
          // continue to otherLoop and peek at the newly pushed
          // neighbor.
          continue outerLoop;
        }
      }
      // If the current vertex didn't have any unvisited
      // neighbors, you know you've reached a dead end and can
      // pop it off the stack.
      stack.pop();
    }
    return visited;
  }
}

extension RecursiveDfs<E> on Graph<E> {
  void _dfs(Vertex<E> source, List<Vertex<E>> visited, Set<Vertex<E>> pushed) {
    // Mark the source vertex as visited.
    pushed.add(source);
    visited.add(source);
    // Visit every neighboring edges.
    final neighbors = edges(source);
    for (var edge in neighbors) {
      // As long as the adjacent vertex has not been visited
      // yet, continue to dive deeper down the branch
      // recursively.
      if (!pushed.contains(edge.destination)) {
        _dfs(edge.destination, visited, pushed);
      }
    }
  }

  List<Vertex<E>> dfs(Vertex<E> start) {
    // This method initializes visited and pushed and then
    // starts the recursion process.
    //
    List<Vertex<E>> visited = [];
    Set<Vertex<E>> pushed = {};
    _dfs(start, visited, pushed);
    return visited;
  }
}

extension CyclicGraph<E> on Graph<E> {
  /// Performs a depth-first graph traversal by recursively
  /// diving down one path until you find a cycle and
  /// back-tracking by popping off the stack to find another
  /// path.The time complexity is O(V+E).
  bool hasCycle(Vertex<E> source) {
    Set<Vertex<E>> pushed = {};
    return _hasCycle(source, pushed);
  }

  bool _hasCycle(Vertex<E> source, Set<Vertex<E>> pushed) {
    //1 Initialize the algorithm by adding the source vertex.
    pushed.add(source);
    // 2 visit every neighbor edge.
    final outgoingEdges = edges(source);
    // loop through edges
    for (var edge in outgoingEdges) {
      // 3 If the adjacent vertex has not been visited before,
      if (!pushed.contains(edge.destination)) {
        // recursively dive deeper down a branch to check for a
        // cycle.
        if (_hasCycle(edge.destination, pushed)) {
          return true;
        }
      } else {
        // 4 If the adjacent vertex has been visited before,
        // you've found a cycle.
        return true;
      }
    }

    // 5 remove the source vertex so you can continue to find
    // other paths with a potential cycle.
    pushed.remove(source);

    // If you've reached this far, then no cycle was found.
    return false;
  }
}

extension CycleGraph2<E> on Graph<E> {
  bool isCyclic2(Vertex<E> source) {
    // Keep track of pushed vertex.
    final Set<Vertex<E>> pushed = {};
    return _isCyclic2(source, pushed);
  }

  bool _isCyclic2(Vertex<E> source, Set<Vertex<E>> pushed) {
    // Add current source to pushed.
    pushed.add(source);
    // Find outgoing edges.
    final outgoingEdges = edges(source);
    for (final edge in outgoingEdges) {
      if (!pushed.contains(edge.destination)) {
        // If pushed doesn't contain edge destination we go deeper.
        if (_isCyclic2(edge.destination, pushed)) {
          return true;
        }
      } else {
        // If pushed contain edge destination we found a cycle
        return true;
      }
    }
    //remove source to add new adjacent vertex
    pushed.remove(source);

    return false;
  }
}
