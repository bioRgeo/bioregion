// Helpers for the Iterative Hierarchical Consensus Tree (IHCT), see R/IHCT.R.
//
// A tree is described the way stats::hclust does: a "merge" matrix with one
// row per node (n - 1 rows for n sites), each row giving the two things joined
// at that node (a negative number -i means site number i, a positive number k
// means the node created at row k), and a "height" vector giving the level at
// which the two things are joined. Rows are ordered so that a node is always
// created before it is used.
//
// These functions only walk that structure; they never build trees themselves.

#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// Number of sites under each node, and number of site pairs first joined at
// each node (= number of sites on the left x number of sites on the right).
// Used by the closed-form cophenetic correlation of UPGMA trees:
//   sum of squared errors = sum(d^2) - sum_k pairs_k * height_k^2
// [[Rcpp::export]]
List ihct_node_sizes(IntegerMatrix merge) {
  int n_nodes = merge.nrow();
  IntegerVector size(n_nodes);
  NumericVector pairs(n_nodes);
  for (int k = 0; k < n_nodes; k++) {
    int a = merge(k, 0), b = merge(k, 1);
    int size_a = a < 0 ? 1 : size[a - 1];
    int size_b = b < 0 ? 1 : size[b - 1];
    size[k] = size_a + size_b;
    pairs[k] = (double) size_a * (double) size_b;
  }
  return List::create(_["size"] = size, _["pairs"] = pairs);
}

// Cophenetic correlation of a tree with a dissimilarity matrix, computed by
// listing every pair of sites once. The cophenetic distance of two sites is
// the height of the node where they are first joined.
//
// `d` must be the dissimilarity matrix in the order of the tree's sites (site
// number i in `merge` is row/column i of `d`). Works for any linkage method,
// at a cost proportional to the number of pairs.
// [[Rcpp::export]]
double ihct_cophenetic_correlation(IntegerMatrix merge, NumericVector height,
                                   NumericMatrix d) {
  int n_nodes = merge.nrow();
  int n = n_nodes + 1;
  // Place the sites in a left-to-right order such that the sites of every
  // node form a contiguous block: node k covers positions lo[k]..hi[k].
  std::vector<int> lo(n_nodes), hi(n_nodes);
  std::vector<int> position(n);           // position of each site in that order
  {
    // iterative depth-first walk from the root (last row), left child first
    std::vector<int> stack;
    std::vector<int> state(n_nodes, 0);   // 0: not started, 1: left done, 2: finished
    int next_pos = 0;
    stack.push_back(n_nodes - 1);
    while (!stack.empty()) {
      int k = stack.back();
      if (state[k] == 0) {
        lo[k] = next_pos;
        int a = merge(k, 0);
        state[k] = 1;
        if (a < 0) position[-a - 1] = next_pos++; else stack.push_back(a - 1);
      } else if (state[k] == 1) {
        int b = merge(k, 1);
        state[k] = 2;
        if (b < 0) position[-b - 1] = next_pos++; else stack.push_back(b - 1);
      } else {
        hi[k] = next_pos - 1;
        stack.pop_back();
      }
    }
  }
  std::vector<int> site_at(n);            // site at each position
  for (int i = 0; i < n; i++) site_at[position[i]] = i;

  long double sum_d = 0, sum_d2 = 0, sum_c = 0, sum_c2 = 0, sum_dc = 0;
  double n_pairs = 0;
  for (int k = 0; k < n_nodes; k++) {
    double h = height[k];
    int a = merge(k, 0), b = merge(k, 1);
    int a_lo = a < 0 ? position[-a - 1] : lo[a - 1];
    int a_hi = a < 0 ? position[-a - 1] : hi[a - 1];
    int b_lo = b < 0 ? position[-b - 1] : lo[b - 1];
    int b_hi = b < 0 ? position[-b - 1] : hi[b - 1];
    for (int p = a_lo; p <= a_hi; p++) {
      int i = site_at[p];
      for (int q = b_lo; q <= b_hi; q++) {
        double dij = d(i, site_at[q]);
        sum_d += dij; sum_d2 += dij * dij;
        sum_c += h; sum_c2 += h * h; sum_dc += dij * h;
        n_pairs += 1;
      }
    }
  }
  long double cov = sum_dc - sum_d * sum_c / n_pairs;
  long double var_d = sum_d2 - sum_d * sum_d / n_pairs;
  long double var_c = sum_c2 - sum_c * sum_c / n_pairs;
  if (var_d <= 0 || var_c <= 0) return NA_REAL;
  return (double) (cov / std::sqrt(var_d * var_c));
}
