#include "scater.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "utils.h"

#include <stdexcept>
#include <algorithm>
#include <vector>

/* Functions to compute the cumulative sum. */

void check_topset(const Rcpp::IntegerVector& top) {
    for (size_t t=1; t<top.size(); ++t) {
        if (top[t] < top[t-1]) { 
            throw std::runtime_error("numbers of top genes must be sorted"); 
        }
    }
    return;
}

template<typename T, class V, class IT>
void compute_cumsum (typename V::iterator it, size_t ngenes, const Rcpp::IntegerVector& top, IT out) {
    const size_t ntop=top.size();
    if (ntop==0) {
        return;
    }

    std::partial_sort(it, it + std::min(ngenes, size_t(top[ntop-1])), it + ngenes, std::greater<T>());
    int x=0;
    T accumulated=0;

    for (auto target_index : top) {
        while (x<target_index && x<ngenes) { // '<' as top contains 1-based indices.
            accumulated+=*(it+x);
            ++x;
        }
        (*out)=accumulated;
        ++out;
    }

    return;
}

/* Class to compute and store the per-cell statistics. */

template <typename T, class V>
struct per_cell_statistics {
    per_cell_statistics() : limit(0), counter(0) {}
    
    // Constructor and fill command for no subsetting, i.e., statistics computed on all genes.
    per_cell_statistics(size_t ncells, T detection_limit, size_t ngenes, Rcpp::IntegerVector TOP) : 
            top(TOP), limit(detection_limit), counter(0), temporary(ngenes),
            totals(ncells), detected(ncells), percentages(top.size(), ncells) {
        check_topset(top);
        return;
    }

    void fill(typename V::iterator it) {
        const size_t ngenes=temporary.size();
        std::copy(it, it+ngenes, temporary.begin());
        compute_summaries(temporary.begin(), ngenes);
    }
 
    // Constructor and fill command for statistics to be computed on a subset of genes.
    per_cell_statistics(size_t ncells, T detection_limit, Rcpp::IntegerVector sub, Rcpp::IntegerVector TOP) :
            per_cell_statistics(ncells, detection_limit, sub.size(), TOP) {
        subset=sub;
        return;
    }
   
    void fill_subset(typename V::iterator it) {
        auto tmp=temporary.begin();
        for (auto sIt=subset.begin(); sIt!=subset.end(); ++sIt, ++tmp) {
            (*tmp)=*(it + *sIt); 
        }
        compute_summaries(temporary.begin(), subset.size());
        return;
    }
private:
    Rcpp::IntegerVector top;
    T limit;
    size_t counter;

    Rcpp::IntegerVector subset;
    V temporary;

    // This function *will* modify the memory pointed to by 'it', so copying should be done beforehand.
    void compute_summaries (typename V::iterator it, size_t ngenes) {
        auto& total=totals[counter];
        auto& nfeat=detected[counter];
        for (size_t i=0; i<ngenes; ++i) {
            const auto& curval=*(it+i);
            total+=curval;
            if (curval > limit) {
                ++nfeat;
            }
        }

        auto pct_col=percentages.column(counter);
        compute_cumsum<T, V>(it, ngenes, top, pct_col.begin());
        for (auto& p : pct_col) { 
            p/=total;
            p*=100; 
        }
        ++counter;
        return;
    }
public:
    V totals;
    Rcpp::IntegerVector detected;
    Rcpp::NumericMatrix percentages;
};

/* Class to compute per-gene statistics. */

template <typename T, class V>
struct per_gene_statistics {
    per_gene_statistics() : limit(0) {}

    per_gene_statistics(size_t ngenes, T detection_limit) : totals(ngenes), detected(ngenes), limit(detection_limit) {}

    void compute_summaries (typename V::iterator it) {
        auto dIt=detected.begin();
        for (auto tIt=totals.begin(); tIt!=totals.end(); ++tIt, ++it, ++dIt) {
            const auto& curval=(*it);
            (*tIt)+=curval;
            if (curval > limit) {
                ++(*dIt);
            }
        }
        return;
    }

    V totals;
    Rcpp::IntegerVector detected;
private:
    T limit;
};

/* Functions to fill output lists. */

template <typename T, class V>
Rcpp::List create_output_per_cell(const per_cell_statistics<T, V>& PCS) {
    return Rcpp::List::create(PCS.totals, PCS.detected, PCS.percentages);
}

template <typename T, class V>
Rcpp::List create_output_per_gene(const per_gene_statistics<T, V>& PGS) {
    return Rcpp::List::create(PGS.totals, PGS.detected);
}

/* Main loop for combined_qc. */

template <typename T, class V, class M>
SEXP combined_qc_internal(M mat, Rcpp::IntegerVector start, Rcpp::IntegerVector end, Rcpp::List featcon, Rcpp::List cellcon, Rcpp::IntegerVector topset, V detection_limit) {
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();

    if (detection_limit.size()!=1) {
        throw std::runtime_error("detection limit should be a scalar");
    }
    T limit=detection_limit[0];

    // Defining the subset of cells for which to perform the calculation.
    const size_t firstcell=check_integer_scalar(start, "first cell index");
    const size_t lastcell=check_integer_scalar(end, "last cell index");
    if (firstcell > lastcell || lastcell > ncells) {
        throw std::runtime_error("cell indices for parallel execution are out of range");
    }
    const size_t n_usedcells=lastcell - firstcell;

    // Setting up per-cell statistics (for each feature control set).
    const size_t nfcontrols=featcon.size();
    per_cell_statistics<T, V> all_PCS(n_usedcells, limit, ngenes, topset);
    std::vector<per_cell_statistics<T, V> > control_PCS(nfcontrols);

    for (size_t fx=0; fx<nfcontrols; ++fx) {
        Rcpp::IntegerVector current=process_subset_vector(featcon[fx], ngenes, false); // converts to zero-index.
        control_PCS[fx]=per_cell_statistics<T, V>(n_usedcells, limit, current, topset);
    }

    // Setting up per-feature statistics (for each cell control set).
    const size_t nccontrols=cellcon.size();
    std::vector<std::vector<size_t> > chosen_ccs(n_usedcells); 

    per_gene_statistics<T, V> all_PGS(ngenes, limit);
    std::vector<per_gene_statistics<T, V> > control_PGS(nccontrols);
    
    for (size_t cx=0; cx<nccontrols; ++cx) {
        Rcpp::IntegerVector current=process_subset_vector(cellcon[cx], ncells, false); // converts to zero-index.
        for (auto curcell : current) {
            size_t cur_index=curcell;
            if (firstcell <= cur_index && cur_index < lastcell) {
                chosen_ccs[cur_index - firstcell].push_back(cx); 
            }
        }
        control_PGS[cx]=per_gene_statistics<T, V>(ngenes, limit);
    }

    // Running through the requested stretch of cells.
    V holder(ngenes);
    for (size_t c=0; c<n_usedcells; ++c) {
        mat->get_col(c + firstcell, holder.begin());

        all_PCS.fill(holder.begin());
        for (size_t fx=0; fx<nfcontrols; ++fx) {
            control_PCS[fx].fill_subset(holder.begin());
        }

        all_PGS.compute_summaries(holder.begin());
        auto& chosen_cc=chosen_ccs[c];
        for (auto& cx : chosen_cc) {
            control_PGS[cx].compute_summaries(holder.begin());
        }
    }

    // Creating a list for all per-cell statistics, and again for all per-feature statistics.
    Rcpp::List output_per_cell(1+nfcontrols);
    output_per_cell[0]=create_output_per_cell(all_PCS);
    for (size_t fx=0; fx<nfcontrols; ++fx) { 
        output_per_cell[fx+1]=create_output_per_cell(control_PCS[fx]);
    }

    Rcpp::List output_per_gene(1+nccontrols);
    output_per_gene[0]=create_output_per_gene(all_PGS);
    for (size_t cx=0; cx<nccontrols; ++cx) { 
        output_per_gene[cx+1]=create_output_per_gene(control_PGS[cx]);
    }

    return Rcpp::List::create(output_per_cell, output_per_gene);
}

SEXP combined_qc(SEXP matrix, SEXP start, SEXP end, SEXP featcon, SEXP cellcon, SEXP top, SEXP limit) {
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(matrix);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(matrix);
        return combined_qc_internal<int, Rcpp::IntegerVector>(mat.get(), start, end, featcon, cellcon, top, limit);
    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(matrix);
        return combined_qc_internal<double, Rcpp::NumericVector>(mat.get(), start, end, featcon, cellcon, top, limit);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}

/* Main loop for top_cumprop. */

template <typename T, class V, class M>
SEXP top_cumprop_internal(M mat, Rcpp::IntegerVector topset) {
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();

    check_topset(topset);
    Rcpp::NumericMatrix percentages(topset.size(), ncells);
    V holder(ngenes);

    for (size_t c=0; c<ncells; ++c) {
        mat->get_col(c, holder.begin());
        double totals=std::accumulate(holder.begin(), holder.end(), T(0));

        auto cur_col=percentages.column(c);
        compute_cumsum<T, V>(holder.begin(), ngenes, topset, cur_col.begin());
        for (auto& p : cur_col) { p/=totals; }
    }

    return percentages;
}

SEXP top_cumprop(SEXP matrix, SEXP top) {
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(matrix);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(matrix);
        return top_cumprop_internal<int, Rcpp::IntegerVector>(mat.get(), top);
    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(matrix);
        return top_cumprop_internal<double, Rcpp::NumericVector>(mat.get(), top);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}
