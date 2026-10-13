// Helpers for the Iterative Hierarchical Consensus Tree (IHCT), see R/IHCT.R,
// and for tree_eval() in R/utils.R.
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
#include <algorithm>
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

// The dissimilarities of a group of sites, in the order they are given, as the
// vector stats::dist uses: the lower triangle read column by column.
//
// `sites` holds the rows of `dist_mat` the group is made of, already shuffled
// by the caller, so the result is exactly what stats::as.dist(dist_mat[sites,
// sites]) returns. Going straight from the full matrix to that vector saves
// the two square copies the R route makes (the sub-matrix, then the index
// matrices as.dist() builds with row() and col()), which together cost about
// twice as much as the clustering they feed.
// [[Rcpp::export]]
NumericVector ihct_shuffled_dist(NumericMatrix dist_mat, IntegerVector sites) {
  int m = sites.size();
  int n = dist_mat.nrow();
  std::vector<int> row(m);
  for (int i = 0; i < m; i++) row[i] = sites[i] - 1;   // 1-based rows from R

  NumericVector out(((R_xlen_t) m * (m - 1)) / 2);
  const double* src = &dist_mat[0];
  R_xlen_t k = 0;
  for (int j = 0; j < m - 1; j++) {
    const double* column = src + (R_xlen_t) row[j] * n;
    for (int i = j + 1; i < m; i++) out[k++] = column[row[i]];
  }
  out.attr("Size") = m;
  out.attr("Diag") = false;
  out.attr("Upper") = false;
  out.attr("method") = "user";
  out.attr("class") = "dist";
  return out;
}

// Fit of a tree to a dissimilarity matrix, over every pair of sites:
// - cophcor: the cophenetic correlation, i.e. the Pearson correlation between
//   the dissimilarities and the cophenetic distances (Sokal & Rohlf 1962);
// - msd: the mean squared difference between them (Maire et al. 2015).
// The cophenetic distance of two sites is the height of the node where they
// are first joined, so the pairs are listed node by node and no cophenetic
// matrix is ever built. Works for any linkage method, at a cost proportional
// to the number of pairs.
//
// By default `d` is the dissimilarity matrix in the order of the tree's sites
// (site number i in `merge` is row/column i of `d`). Give `leaf_site` -- the
// row of `d` each site of the tree stands for, as ihct_prune_tree() takes it --
// to read the dissimilarities from a larger matrix instead, which saves
// building the sub-matrix of the group.
//
// The correlation is computed in two passes, as stats::cor() does: the means
// first, then the sums of products of the deviations from the means. The
// one-pass formula (sum of squares minus squared sum over n) loses most of its
// digits when the dissimilarities are close to each other compared to their
// mean, e.g. when many are near 1, and all the more so where long double is
// no wider than double (macOS arm64). The first pass is cheap: the mean of
// the cophenetic distances comes from the heights and node sizes alone, and
// the mean of the dissimilarities is read straight through the matrix.
// For the same reason the msd sums the squared differences themselves.
//
// cophcor is NA when the dissimilarities or the heights are all equal.
// [[Rcpp::export]]
NumericVector tree_eval_cpp(IntegerMatrix merge, NumericVector height,
                            NumericMatrix d,
                            Nullable<IntegerVector> leaf_site = R_NilValue) {
  int n_nodes = merge.nrow();
  int n = n_nodes + 1;
  // the sites of every node are then the block of positions lo[k]..hi[k]
  std::vector<int> lo(n_nodes), hi(n_nodes), position(n), site_at(n);
  tree_blocks(merge, lo, hi, position, site_at);
  // the two children of node k are the blocks a_lo[k]..a_hi[k] and
  // b_lo[k]..b_hi[k], and node k joins every site of one to every site of the
  // other
  std::vector<int> a_lo(n_nodes), a_hi(n_nodes), b_lo(n_nodes), b_hi(n_nodes);
  for (int k = 0; k < n_nodes; k++) {
    int a = merge(k, 0), b = merge(k, 1);
    a_lo[k] = a < 0 ? position[-a - 1] : lo[a - 1];
    a_hi[k] = a < 0 ? position[-a - 1] : hi[a - 1];
    b_lo[k] = b < 0 ? position[-b - 1] : lo[b - 1];
    b_hi[k] = b < 0 ? position[-b - 1] : hi[b - 1];
  }

  // where each site of the tree sits in `d`
  std::vector<int> row(n);
  if (leaf_site.isNull()) {
    for (int i = 0; i < n; i++) row[i] = i;
  } else {
    IntegerVector site(leaf_site);
    for (int i = 0; i < n; i++) row[i] = site[i] - 1;
  }

  // first pass: the means
  long double n_pairs = 0, sum_c = 0, sum_d = 0;
  bool same_c = true, same_d = true;
  for (int k = 0; k < n_nodes; k++) {
    double pairs = (double) (a_hi[k] - a_lo[k] + 1) * (b_hi[k] - b_lo[k] + 1);
    n_pairs += pairs;
    sum_c += pairs * height[k];
    same_c = same_c && height[k] == height[0];
  }
  // the rows sorted, so that each column of `d` is read downwards
  std::vector<int> sorted_row(row);
  std::sort(sorted_row.begin(), sorted_row.end());
  const double* src = &d[0];
  R_xlen_t n_row = d.nrow();
  double d_first = n > 1 ? src[(R_xlen_t) sorted_row[0] * n_row + sorted_row[1]] : 0;
  for (int j = 0; j < n - 1; j++) {
    const double* column = src + (R_xlen_t) sorted_row[j] * n_row;
    for (int i = j + 1; i < n; i++) {
      double dij = column[sorted_row[i]];
      sum_d += dij;
      same_d = same_d && dij == d_first;
    }
  }
  long double mean_c = sum_c / n_pairs, mean_d = sum_d / n_pairs;
  // When all the values are equal, the sum divided by n can still be a hair
  // off that value, and every deviation from the mean would then be the same
  // tiny non-zero number: cophcor would come out as some arbitrary value
  // instead of NA. stats::cor() avoids this by refining the mean; here the
  // mean is set to the value itself.
  if (same_c && n_nodes > 0) mean_c = height[0];
  if (same_d && n > 1) mean_d = d_first;

  // second pass: the deviations from the means
  long double sum_dev_dc = 0, sum_dev_d2 = 0, sum_dev_c2 = 0, sum_err2 = 0;
  for (int k = 0; k < n_nodes; k++) {
    double h = height[k];
    long double dev_c = h - mean_c;   // the same for every pair of the node
    long double node_dev_d = 0;
    for (int p = a_lo[k]; p <= a_hi[k]; p++) {
      int i = row[site_at[p]];
      for (int q = b_lo[k]; q <= b_hi[k]; q++) {
        double dij = d(i, row[site_at[q]]);
        long double dev_d = dij - mean_d;
        node_dev_d += dev_d;
        sum_dev_d2 += dev_d * dev_d;
        double err = dij - h;
        sum_err2 += err * err;
      }
    }
    double pairs = (double) (a_hi[k] - a_lo[k] + 1) * (b_hi[k] - b_lo[k] + 1);
    sum_dev_dc += dev_c * node_dev_d;
    sum_dev_c2 += pairs * dev_c * dev_c;
  }

  double cophcor = NA_REAL;
  if (sum_dev_d2 > 0 && sum_dev_c2 > 0) {
    cophcor = (double) (sum_dev_dc / std::sqrt(sum_dev_d2 * sum_dev_c2));
    cophcor = std::min(1.0, std::max(-1.0, cophcor));   // as stats::cor()
  }
  return NumericVector::create(_["cophcor"] = cophcor,
                               _["msd"] = (double) (sum_err2 / n_pairs));
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
