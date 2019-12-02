#include "Rcpp.h"

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
    size_t x=0;
    T accumulated=0;

    for (size_t target_index : top) {
        while (x<target_index && x<ngenes) { // '<' as top contains 1-based indices.
            accumulated+=*(it+x);
            ++x;
        }
        (*out)=accumulated;
        ++out;
    }

    return;
}

/*******************************************************
 * Class to compute and store the per-cell statistics. *
 *******************************************************/

template <typename T, class V>
struct per_cell_statistics {
    per_cell_statistics() : limit(0), counter(0) {}

    per_cell_statistics(size_t ncells, T detection_limit) :
        limit(detection_limit), counter(0), totals(ncells), detected(ncells) {}
    
    per_cell_statistics(size_t ncells, T detection_limit, size_t ngenes, Rcpp::IntegerVector TOP) : 
        limit(detection_limit), counter(0), totals(ncells), detected(ncells),
        top(TOP), temporary(ngenes), percentages(top.size(), ncells) 
    {
        check_topset(top);
        return;
    }

    void fill(typename V::iterator it, size_t ngenes) {
        auto& total=totals[counter];
        auto& nfeat=detected[counter];
        auto copy=it;

        for (size_t i=0; i<ngenes; ++i, ++copy) {
            const auto& curval=*copy;
            total+=curval;
            if (curval > limit) {
                ++nfeat;
            }
        }

        if (top.size()) {            
            std::copy(it, it+ngenes, temporary.begin());
            auto pct_col=percentages.column(counter);
            compute_cumsum<T, V>(temporary.begin(), ngenes, top, pct_col.begin());
            for (auto& p : pct_col) { 
                p/=total;
                p*=100; 
            }
        }

        ++counter;
        return;
    }
 
    void fill(typename V::iterator it, Rcpp::IntegerVector subset) {
        auto& total=totals[counter];
        auto& nfeat=detected[counter];
        for (auto s : subset) {
            const auto& curval=*(it+s);
            total+=curval;
            if (curval > limit) {
                ++nfeat;
            }
        }
        ++counter;
        return;
    }

    Rcpp::List yield() const {
        return Rcpp::List::create(totals, detected, percentages);
    }
private:
    T limit;
    size_t counter;
public:
    V totals;
    Rcpp::IntegerVector detected;
private:    
    Rcpp::IntegerVector top;
    V temporary;
public:
    Rcpp::NumericMatrix percentages;
};

template <class M>
Rcpp::List cell_qc_internal(Rcpp::RObject input, Rcpp::List featcon, Rcpp::IntegerVector topset, 
    typename M::vector detection_limit) 
{
    auto mat=beachmat::create_matrix<M>(input);
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();
    typename M::vector holder(ngenes);

    if (detection_limit.size()!=1) {
        throw std::runtime_error("detection limit should be a scalar");
    }
    typename M::type limit=detection_limit[0];

    // Setting up per-cell statistics (for each feature control set).
    typedef per_cell_statistics<typename M::type, typename M::vector> cell_stats;
    cell_stats all_PCS(ncells, limit, ngenes, topset);

    const size_t nfcontrols=featcon.size();
    std::vector<Rcpp::IntegerVector> all_subsets(nfcontrols);
    std::vector<cell_stats> control_PCS(nfcontrols);
    for (size_t fx=0; fx<nfcontrols; ++fx) {
        all_subsets[fx]=process_subset_vector(featcon[fx], ngenes, false); // converts to zero-index.
        control_PCS[fx]=cell_stats(ncells, limit);
    }

    // Running through the requested stretch of cells.
    for (size_t c=0; c<ncells; ++c) {
        mat->get_col(c, holder.begin());
        all_PCS.fill(holder.begin(), ngenes);
        for (size_t fx=0; fx<nfcontrols; ++fx) {
            control_PCS[fx].fill(holder.begin(), all_subsets[fx]);
        }
    }

    // Creating a list for all per-cell statistics.
    Rcpp::List output_per_cell(nfcontrols);
    for (size_t fx=0; fx<nfcontrols; ++fx) { 
        output_per_cell[fx]=control_PCS[fx].yield();
    }

    return Rcpp::List::create(all_PCS.yield(), output_per_cell);
}

// [[Rcpp::export(rng=false)]]
Rcpp::RObject per_cell_qc(Rcpp::RObject matrix, Rcpp::List featcon, Rcpp::IntegerVector top, SEXP limit) {
    auto mattype=beachmat::find_sexp_type(matrix);
    if (mattype==INTSXP) {
        return cell_qc_internal<beachmat::integer_matrix>(matrix, featcon, top, limit);
    } else if (mattype==REALSXP) {
        return cell_qc_internal<beachmat::numeric_matrix>(matrix, featcon, top, limit);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
}

/*******************************************************
 * Class to compute and store the per-gene statistics. *
 *******************************************************/

template <typename T, class V>
struct per_gene_statistics {
    per_gene_statistics() : limit(0) {}

    per_gene_statistics(size_t ngenes, T detection_limit) : 
        limit(detection_limit), counter(0), totals(ngenes), detected(ngenes) {}

    void fill(typename V::iterator it, size_t ncells) {
        auto& total=totals[counter];
        auto& nfeat=detected[counter];
        auto copy=it;

        for (size_t i=0; i<ncells; ++i, ++copy) {
            const auto& curval=(*copy);
            total+=curval;
            if (curval > limit) {
                ++nfeat;
            }
        }

        total/=ncells;
        nfeat/=ncells;
        ++counter; 
        return;
    }

    void fill(typename V::iterator it, Rcpp::IntegerVector subset) {
        auto& total=totals[counter];
        auto& nfeat=detected[counter];

        for (auto s : subset){ 
            const auto& curval=*(it + s);
            total+=curval;
            if (curval > limit) {
                ++nfeat;
            }
        }

        total/=subset.size();
        nfeat/=subset.size();
        ++counter;
        return;
    }

    Rcpp::List yield() {
        return Rcpp::List::create(totals, detected);
    }
private:
    T limit;
    size_t counter;
public:
    Rcpp::NumericVector totals;
    Rcpp::NumericVector detected;
};

template <class M>
Rcpp::RObject feature_qc_internal(Rcpp::RObject input, Rcpp::List cellcon, typename M::vector detection_limit) 
{
    auto mat=beachmat::create_matrix<M>(input);
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();
    typename M::vector holder(ncells);

    if (detection_limit.size()!=1) {
        throw std::runtime_error("detection limit should be a scalar");
    }
    typename M::type limit=detection_limit[0];

    // Setting up per-feature statistics (for each cell control set).
    typedef per_gene_statistics<typename M::type, typename M::vector> gene_stats;
    gene_stats all_PGS(ngenes, limit);

    const size_t nccontrols=cellcon.size();
    std::vector<Rcpp::IntegerVector> all_subsets(nccontrols);
    std::vector<gene_stats> control_PGS(nccontrols);
    for (size_t cx=0; cx<nccontrols; ++cx) {
        all_subsets[cx]=process_subset_vector(cellcon[cx], ngenes, false); // converts to zero-index.
        control_PGS[cx]=gene_stats(ngenes, limit);
    }
    
    // Running through the requested stretch of genes.
    for (size_t g=0; g<ngenes; ++g) {
        mat->get_row(g, holder.begin());
        all_PGS.fill(holder.begin(), ncells);
        for (size_t cx=0; cx<nccontrols; ++cx) {
            control_PGS[cx].fill(holder.begin(), all_subsets[cx]);
        }
    }

    // Creating a list for all per-gene statistics.
    Rcpp::List output_per_gene(nccontrols);
    for (size_t cx=0; cx<nccontrols; ++cx) { 
        output_per_gene[cx]=control_PGS[cx].yield();
    }

    return Rcpp::List::create(all_PGS.yield(), output_per_gene);
}

// [[Rcpp::export(rng=false)]]
Rcpp::RObject per_feature_qc(Rcpp::RObject matrix, Rcpp::List cellcon, SEXP limit) {
    auto mattype=beachmat::find_sexp_type(matrix);
    if (mattype==INTSXP) {
        return feature_qc_internal<beachmat::integer_matrix>(matrix, cellcon, limit);
    } else if (mattype==REALSXP) {
        return feature_qc_internal<beachmat::numeric_matrix>(matrix, cellcon, limit);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
}

/*******************************************************
 * Separate calculation of the top proportion per cell *
 *******************************************************/

template <class M>
Rcpp::NumericMatrix top_cumprop_internal(Rcpp::RObject incoming, Rcpp::IntegerVector topset) {
    auto mat=beachmat::create_matrix<M>(incoming);
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();

    check_topset(topset);
    Rcpp::NumericMatrix percentages(topset.size(), ncells);
    typename M::vector holder(ngenes); 

    for (size_t c=0; c<ncells; ++c) {
        mat->get_col(c, holder.begin()); // need to copy as cumsum will change ordering.
        double totals=std::accumulate(holder.begin(), holder.end(), static_cast<typename M::type>(0));

        auto cur_col=percentages.column(c);
        compute_cumsum<typename M::type, typename M::vector>(holder.begin(), ngenes, topset, cur_col.begin());
        for (auto& p : cur_col) { p/=totals; }
    }

    return percentages;
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericMatrix top_cumprop(Rcpp::RObject matrix, Rcpp::IntegerVector top) {
    auto mattype=beachmat::find_sexp_type(matrix);
    if (mattype==INTSXP) {
        return top_cumprop_internal<beachmat::integer_matrix>(matrix, top);
    } else if (mattype==REALSXP) {
        return top_cumprop_internal<beachmat::numeric_matrix>(matrix, top);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
}
