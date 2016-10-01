public class ReduceGraph {

	static public void main( String arg[] ) throws java.io.IOException {

		java.util.Random r = new java.util.Random(Long.parseLong(arg[2]));
		final it.unimi.dsi.webgraph.ImmutableGraph graph = it.unimi.dsi.webgraph.ImmutableGraph.load(arg[0]);
		int n = Integer.parseInt(arg[1]);
		int nodes = graph.numNodes();
		long arcs = graph.numArcs();
		int startNode;
		do startNode = r.nextInt(nodes) + 1;
		while (graph.successorArray(startNode).length == 0);
		int i = 0;
		int e = 0;
		int[] g = new int[n * n];
		java.util.ArrayList<Integer> al = new java.util.ArrayList<Integer>(n);
		java.util.ArrayList<Integer> a = new java.util.ArrayList<Integer>(n * n);
		al.add(startNode);

		while (i < n) {
			int[] nbrs = graph.successorArray(al.get(i));
			for (int j = 0; j < nbrs.length; j++) {
				if (!al.contains(nbrs[j])) {
					if (al.size() < n) al.add(nbrs[j]);
					else continue;
				}
				int k = al.indexOf(nbrs[j]);
				if (g[i * n + k] == 0) {
					g[i * n + k] = g[k * n + i] = ++e;
					a.add(i);
					a.add(k);
				}
			}
			i++;
		}

		for (i = 0; i < e; i++)
			System.out.format("%d %d\n", a.get(2 * i), a.get(2 * i + 1));
	}
}
