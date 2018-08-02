#include "scater.h"

/* Class to compute and store the per-cell statistics. */

template <typename T, class V>
struct per_cell_statistics {
    per_cell_statistics() : counter(0), limit(0) {}
    
    per_cell_statistics(size_t ncells, T detection_limit, Rcpp::IntegerVector TOP, size_t ngenes) : 
            top(TOP), limit(detection_limit), counter(0), 
            totals(ncells), detected(ncells), percentages(top.size(), ncells),
            temporary(ngenes) {
        for (size_t t=1; t<top.size(); ++t) {
            if (top[t] < top[t-1]) { 
                throw std::runtime_error("numbers of top genes must be sorted"); 
            }
        }
        return;
    }

    per_cell_statistics(size_t ncells, T detection_limit, Rcpp::IntegerVector sub, Rcpp::IntegerVector TOP) :
            per_cell_statistics(ncells, detection_limit, TOP, sub.size()) {
        subset=sub;
        return;
    }

    per_cell_statistics(size_t ncells, T detection_limit, size_t ngenes, Rcpp::IntegerVector TOP) : 
            per_cell_statistics(ncells, detection_limit, TOP, ngenes) {}

    void fill(typename V::iterator it) {
        const size_t ngenes=temporary.size();
        compute_summaries(it, ngenes);
        std::copy(it, it+ngenes, temporary.begin());
        compute_percentages(temporary.begin(), ngenes);
        ++counter;
        return;
    }

    void fill_subset(typename V::iterator it) {
        auto tmp=temporary.begin();
        for (auto sIt=subset.begin(); sIt!=subset.end(); ++sIt, ++tmp) {
            (*tmp)=*(it + *sIt - 1); // convert to zero indexing.
        }
        compute_summaries(temporary.begin(), subset.size());
        compute_percentages(temporary.begin(), subset.size());
        ++counter;
        return;
    }
private:
    Rcpp::IntegerVector top;
    T limit;
    size_t counter;

    Rcpp::IntegerVector subset;
    V temporary;

    void compute_summaries (typename V::iterator it, size_t ngenes) {
        auto& total=totals[counter];
        auto& nfeat=detected[counter];
        for (size_t i=0; i<ngenes; ++i, ++it) {
            const auto& curval=(*it);
            total+=curval;
            if (curval > limit) {
                ++nfeat;
            }
        }
        return;
    }

    void compute_percentages (typename V::iterator it, size_t ngenes) {
        if (top.size()==0) {
            return;
        }

        std::partial_sort(it, 
                it + std::min(ngenes, size_t(top[top.size()-1])), 
                it + ngenes, 
                std::greater<T>());

        int x=0;
        T accumulated=0;
        double total=totals[counter];
        auto pct_col=percentages.column(counter);

        for (size_t t=0; t<top.size(); ++t) { 
            int target_index=top[t] - 1; // 0-based index.
            while (x<=target_index && x<ngenes) {
                accumulated+=*(it+x);
                ++x;
            }
            pct_col[t]=double(accumulated)/total * 100;
        }
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

/* Main loop. */

template <typename T, class V, class M>
SEXP combined_qc_internal(M mat, Rcpp::List featcon, Rcpp::List cellcon, Rcpp::IntegerVector topset, V detection_limit) {
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();

    if (detection_limit.size()!=1) {
        throw std::runtime_error("detection limit should be a scalar");
    }
    T limit=detection_limit[0];

    // Processing feature controls.
    const size_t nfcontrols=featcon.size();
    std::vector<Rcpp::IntegerVector> feat_controls(nfcontrols);

    per_cell_statistics<T, V> all_PCS(ncells, limit, ngenes, topset);
    std::vector<per_cell_statistics<T, V> > control_PCS(nfcontrols);

    for (size_t fx=0; fx<nfcontrols; ++fx) {
        feat_controls[fx]=featcon[fx];
        control_PCS[fx]=per_cell_statistics<T, V>(ncells, limit, feat_controls[fx], topset);
    }

    // Processing cell controls.
    const size_t nccontrols=cellcon.size();
    std::vector<std::vector<size_t> > chosen_ccs(ncells);

    per_gene_statistics<T, V> all_PGS(ngenes, limit);
    std::vector<per_gene_statistics<T, V> > control_PGS(nccontrols);
    
    for (size_t cx=0; cx<nccontrols; ++cx) {
        Rcpp::IntegerVector current=cellcon[cx];
        for (auto g : current) {
            chosen_ccs[g-1].push_back(cx); // converted to zero-based indexing.
        }
        control_PGS[cx]=per_gene_statistics<T, V>(ngenes, limit);
    }

    // Running through all cells.
    V holder(ngenes);
    for (size_t c=0; c<ncells; ++c) {
        auto cIt=mat->get_const_col(c, holder.begin());

        all_PCS.fill(cIt);
        for (size_t fx=0; fx<nfcontrols; ++fx) {
            control_PCS[fx].fill_subset(cIt);
        }

        all_PGS.compute_summaries(cIt);
        auto& chosen_cc=chosen_ccs[c];
        for (auto& cx : chosen_cc) {
            control_PGS[cx].compute_summaries(cIt);
        }
    }

    // Creating a list for all per-cell statistics.
    Rcpp::List output_per_cell(1+nfcontrols);
    output_per_cell[0]=create_output_per_cell(all_PCS);
    for (size_t fx=0; fx<nfcontrols; ++fx) { 
        output_per_cell[fx+1]=create_output_per_cell<T, V>(control_PCS[fx]);
    }

    Rcpp::List output_per_gene(1+nccontrols);
    output_per_gene[0]=create_output_per_gene(all_PGS);
    for (size_t cx=0; cx<nccontrols; ++cx) { 
        output_per_gene[cx+1]=create_output_per_gene<T, V>(control_PGS[cx]);
    }

    return Rcpp::List::create(output_per_cell, output_per_gene);
}

SEXP combined_qc(SEXP matrix, SEXP featcon, SEXP cellcon, SEXP top, SEXP limit) {
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(matrix);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(matrix);
        return combined_qc_internal<int, Rcpp::IntegerVector>(mat.get(), featcon, cellcon, top, limit);
    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(matrix);
        return combined_qc_internal<double, Rcpp::NumericVector>(mat.get(), featcon, cellcon, top, limit);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}

