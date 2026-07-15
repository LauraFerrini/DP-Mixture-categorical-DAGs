#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <cstdint>
#include <limits>

using namespace Rcpp;

// mode: 0 = parents, 1 = family
// base: 0 => values are 0..L-1, 1 => values are 1..L

// [[Rcpp::export]]
IntegerMatrix match_counts_one_cluster_all_rows_fast(IntegerMatrix Y,
                                                     IntegerVector rows_k,  // 1-based
                                                     NumericMatrix dag,     // q x q, dag[i,j]!=0 => i -> j
                                                     IntegerVector levels,  // length q
                                                     int mode = 1,          // 0=parents, 1=family
                                                     int base = 0) {        // 0 or 1
  const int n  = Y.nrow();
  const int q  = Y.ncol();
  const int nk = rows_k.size();
  
  // --- Minimal sanity checks (keep these; they are cheap and prevent UB) ---
  if (dag.nrow() != q || dag.ncol() != q) stop("dag must be q x q");
  if (levels.size() != q) stop("levels length must be q");
  if (mode != 0 && mode != 1) stop("mode must be 0 (parents) or 1 (family)");
  if (base != 0 && base != 1) stop("base must be 0 or 1");
 
  // Convert cluster rows to 0-based (NO bounds checks for speed; assume valid)
  std::vector<int> rk0(nk);
  for (int t = 0; t < nk; ++t) rk0[t] = rows_k[t] - 1;
 
  // Precompute column sets per node (0-based)
  std::vector< std::vector<int> > node_cols(q);
  for (int j0 = 0; j0 < q; ++j0) {
    auto &cols = node_cols[j0];
    cols.reserve(q);
    if (mode == 1) cols.push_back(j0);            // family includes child
    for (int i = 0; i < q; ++i) if (dag(i, j0) != 0.0) cols.push_back(i);
  }
  
  IntegerMatrix out(n, q);

  // Main loop over nodes
  for (int j0 = 0; j0 < q; ++j0) {
    const auto &cols0 = node_cols[j0];
   
    // parents-mode with no parents: empty tuple => all nk match
    if (mode == 0 && cols0.empty()) {
      for (int i = 0; i < n; ++i) out(i, j0) = nk;
      continue;
    }
    // ---------- Fast packed uint64 key path (NO NA/range checks) ----------
    // NOTE: This assumes:
    //   - no NA in Y
    //   - values are in-range for given levels/base
    //   - the radix product doesn't overflow uint64
    // If any of those can happen, add checks back or use a safe fallback.
  
    std::unordered_map<std::uint64_t, int> counts;
    counts.reserve((std::size_t)(nk * 1.3));
  
    // Count keys among cluster rows
    for (int t = 0; t < nk; ++t) {
      const int row = rk0[t];
    
      std::uint64_t key  = 0;
      std::uint64_t mult = 1;
     
      // pack key
      for (int c : cols0) {
        const int L   = levels[c];
        const int val = Y(row, c) - base;   // digit in 0..L-1 (assumed)
        key  += (std::uint64_t)val * mult;
        mult *= (std::uint64_t)L;
      }
      
      ++counts[key];
   }
    
    // Fill all rows (0 if configuration not present in cluster)
    for (int i = 0; i < n; ++i) {
      std::uint64_t key  = 0;
      std::uint64_t mult = 1;
     
      for (int c : cols0) {
        const int L   = levels[c];
        const int val = Y(i, c) - base;
        key  += (std::uint64_t)val * mult;
        mult *= (std::uint64_t)L;
      }
      
      auto it = counts.find(key);
      out(i, j0) = (it == counts.end()) ? 0 : it->second;
   }
  }
  
  return out;
} 



// [[Rcpp::export]]
NumericVector fa_prod_div(double a,
                          NumericVector Ical,
                          List fa_list) {
  
  int q = fa_list.size();
  NumericVector out(q);
  
  for (int j = 0; j < q; ++j) {
    IntegerVector idx = fa_list[j];
   
    double prod = 1.0;
    int m = idx.size();
  
    for (int k = 0; k < m; ++k) {
      int ii = idx[k] - 1;  // R -> C++ index
      prod *= Ical[ii];
   }
    
    out[j] = a / prod;
  }
  return out;
}


