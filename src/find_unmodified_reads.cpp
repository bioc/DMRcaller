// annotating the list of unmodif and coverage using rcpp script (find_unmod_and_coverage)

#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;

// [[Rcpp::export]]
List findUnmodCoverageCpp(
  List cov_reads,      // each element: IntegerVector
  List cov_strands,    // each element: CharacterVector (strand info)
  List mod_reads,      // each element: IntegerVector (already index)
  CharacterVector ref_strands
) {
int n = cov_reads.size();
List unmod_list(n);
IntegerVector readsN(n);

for (int i = 0; i < n; ++i) {
  IntegerVector reads = cov_reads[i];
  CharacterVector strands = cov_strands[i];
  IntegerVector mods = mod_reads[i];
  std::string ref_str_val = as<std::string>(ref_strands[i]);
  std::unordered_set<int> mod_set(mods.begin(), mods.end());

  std::vector<int> unmod_indices;
  std::unordered_set<int> same_strand_reads;

  for (int j = 0; j < reads.size(); ++j) {
    int this_read = reads[j];
    std::string this_strand = as<std::string>(strands[j]);
    if (this_strand == ref_str_val) {
      same_strand_reads.insert(this_read);
      if (mod_set.find(this_read) == mod_set.end()) {
        unmod_indices.push_back(this_read);
      }
    }
  }
  unmod_list[i] = wrap(unmod_indices);
  readsN[i] = same_strand_reads.size();
}
return List::create(
  _["unmod_list"] = unmod_list,
  _["readsN"] = readsN
);
}