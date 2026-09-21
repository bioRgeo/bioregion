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

// Depth-first walk of a tree, left child first. It places the sites in a
// left-to-right order in which the sites of every node form a contiguous block:
// node k covers the positions block_lo[k] to block_hi[k], leaf j sits at
// position leaf_pos[j], and site_at gives the leaf at each position.
static void tree_blocks(const IntegerMatrix& merge, std::vector<int>& block_lo,
                        std::vector<int>& block_hi, std::vector<int>& leaf_pos,
                        std::vector<int>& site_at) {
  int n_nodes = merge.nrow();
  std::vector<int> stack;
  std::vector<int> state(n_nodes, 0);   // 0: not started, 1: left done, 2: finished
  int next_pos = 0;
  stack.push_back(n_nodes - 1);         // the root is the last row
  while (!stack.empty()) {
    int k = stack.back();
    if (state[k] == 0) {
      block_lo[k] = next_pos;
      state[k] = 1;
      int a = merge(k, 0);
      if (a < 0) { leaf_pos[-a - 1] = next_pos; site_at[next_pos++] = -a - 1; }
      else stack.push_back(a - 1);
    } else if (state[k] == 1) {
      state[k] = 2;
      int b = merge(k, 1);
      if (b < 0) { leaf_pos[-b - 1] = next_pos; site_at[next_pos++] = -b - 1; }
      else stack.push_back(b - 1);
    } else {
      block_hi[k] = next_pos - 1;
      stack.pop_back();
    }
  }
}

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
  // the sites of every node are then the block of positions lo[k]..hi[k]
  std::vector<int> lo(n_nodes), hi(n_nodes), position(n), site_at(n);
  tree_blocks(merge, lo, hi, position, site_at);

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

// Remove sites from a tree and recompute the heights that change.
//
// The tree is given as `merge` and `height` (as above), the number of site
// pairs joined at each node (`pairs`, i.e. sites on the left x sites on the
// right) and `leaf_site`, the site of the dissimilarity matrix `d` that every
// leaf of the tree stands for (a 1-based row of `d`). `keep` says, for every
// site of `d`, whether it stays.
//
// The leaves that go are dropped, the nodes left with a single child are
// contracted, and every node that lost sites gets the mean dissimilarity
// between the two groups it still joins, which is the height UPGMA would give
// it. Nodes that lost nothing keep their height untouched.
//
// That height is obtained by subtracting the dissimilarities of the lost pairs
// from the node's total, so a node only costs the pairs it lost: removing one
// site from a tree of m sites costs m, not m^2, because the sites it was paired
// with at its successive ancestors are all different.
//
// Returns the pruned tree in the same four pieces. Its root is the last row of
// the `merge` it returns.
// [[Rcpp::export]]
List ihct_prune_tree(IntegerMatrix merge, NumericVector height, NumericVector pairs,
                     IntegerVector leaf_site, LogicalVector keep, NumericMatrix d) {
  int n_nodes = merge.nrow();
  int m = n_nodes + 1;
  std::vector<int> lo(n_nodes), hi(n_nodes), leaf_pos(m), site_at(m);
  tree_blocks(merge, lo, hi, leaf_pos, site_at);

  // The sites in the left-to-right order of the tree: the row of `d` they sit
  // in, whether they stay, and how many stay (or go) before each position, so
  // that the sites a node keeps, or loses, can be listed without scanning its
  // whole block.
  std::vector<int> row_at(m);
  std::vector<int> n_kept_before(m + 1, 0), n_gone_before(m + 1, 0);
  std::vector<int> kept_at, gone_at;
  for (int p = 0; p < m; p++) {
    int site = leaf_site[site_at[p]] - 1;
    row_at[p] = site;
    bool stays = keep[site];
    n_kept_before[p + 1] = n_kept_before[p] + (stays ? 1 : 0);
    n_gone_before[p + 1] = n_gone_before[p] + (stays ? 0 : 1);
    if (stays) kept_at.push_back(p); else gone_at.push_back(p);
  }

  // new leaf numbers, in the order the leaves had in the tree
  std::vector<int> new_leaf(m, 0);
  int n_left = 0;
  for (int j = 0; j < m; j++) if (keep[leaf_site[j] - 1]) new_leaf[j] = ++n_left;
  IntegerVector out_site(n_left);
  for (int j = 0; j < m; j++) if (new_leaf[j]) out_site[new_leaf[j] - 1] = leaf_site[j];

  int n_out = n_left > 0 ? n_left - 1 : 0;
  IntegerMatrix out_merge(n_out, 2);
  NumericVector out_height(n_out), out_pairs(n_out);
  std::vector<int> ref(n_nodes, 0);   // what each node became: 0 nothing left,
  int created = 0;                    // -j a single site, k the new node k
  for (int k = 0; k < n_nodes; k++) {
    int a = merge(k, 0), b = merge(k, 1);
    int ref_a = a < 0 ? -new_leaf[-a - 1] : ref[a - 1];
    int ref_b = b < 0 ? -new_leaf[-b - 1] : ref[b - 1];
    // one side gone: the node is contracted and stands for the other side
    if (ref_a == 0 || ref_b == 0) { ref[k] = ref_a != 0 ? ref_a : ref_b; continue; }

    int a_lo = a < 0 ? leaf_pos[-a - 1] : lo[a - 1];
    int a_hi = a < 0 ? leaf_pos[-a - 1] : hi[a - 1];
    int b_lo = b < 0 ? leaf_pos[-b - 1] : lo[b - 1];
    int b_hi = b < 0 ? leaf_pos[-b - 1] : hi[b - 1];
    int kept_a = n_kept_before[a_hi + 1] - n_kept_before[a_lo];
    int kept_b = n_kept_before[b_hi + 1] - n_kept_before[b_lo];
    double h = height[k];
    if (kept_a != a_hi - a_lo + 1 || kept_b != b_hi - b_lo + 1) {
      // the dissimilarities of the pairs this node no longer joins: a lost site
      // on the left against the whole right, then a site kept on the left
      // against a lost site on the right
      long double lost = 0;
      for (int t = n_gone_before[a_lo]; t < n_gone_before[a_hi + 1]; t++) {
        int row = row_at[gone_at[t]];
        for (int q = b_lo; q <= b_hi; q++) lost += d(row, row_at[q]);
      }
      for (int t = n_kept_before[a_lo]; t < n_kept_before[a_hi + 1]; t++) {
        int row = row_at[kept_at[t]];
        for (int u = n_gone_before[b_lo]; u < n_gone_before[b_hi + 1]; u++)
          lost += d(row, row_at[gone_at[u]]);
      }
      long double total = (long double) pairs[k] * (long double) height[k] - lost;
      h = (double) (total / ((long double) kept_a * (long double) kept_b));
      if (h < 0) h = 0;               // rounding only, the sum cannot be negative
    }
    created++;
    out_merge(created - 1, 0) = ref_a;
    out_merge(created - 1, 1) = ref_b;
    out_height[created - 1] = h;
    out_pairs[created - 1] = (double) kept_a * (double) kept_b;
    ref[k] = created;
  }

  return List::create(_["merge"] = out_merge, _["height"] = out_height,
                      _["pairs"] = out_pairs, _["leaf_site"] = out_site);
}

// The two groups of sites the root of a tree joins, as rows of the
// dissimilarity matrix (see ihct_prune_tree for `leaf_site`).
// [[Rcpp::export]]
List ihct_top_division(IntegerMatrix merge, IntegerVector leaf_site) {
  int n_nodes = merge.nrow();
  int m = n_nodes + 1;
  std::vector<int> lo(n_nodes), hi(n_nodes), leaf_pos(m), site_at(m);
  tree_blocks(merge, lo, hi, leaf_pos, site_at);

  List out(2);
  for (int side = 0; side < 2; side++) {
    int child = merge(n_nodes - 1, side);
    int child_lo = child < 0 ? leaf_pos[-child - 1] : lo[child - 1];
    int child_hi = child < 0 ? leaf_pos[-child - 1] : hi[child - 1];
    IntegerVector sites(child_hi - child_lo + 1);
    for (int p = child_lo; p <= child_hi; p++) sites[p - child_lo] = leaf_site[site_at[p]];
    out[side] = sites;
  }
  return out;
}
