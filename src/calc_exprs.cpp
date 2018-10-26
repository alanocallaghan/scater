#include "scater.h"

template <class V, class M>
class normalizer {
public:
    normalizer(M mat, Rcpp::List sf_list, Rcpp::IntegerVector sf_to_use, Rcpp::RObject genes_sub) : 
            ptr(mat), vec(mat->get_nrow()), 
            size_factors(sf_list.size()), current_sfs(sf_list.size()), set_id(sf_to_use), 
            subset(process_subset_vector(genes_sub, mat->get_nrow())), smallest(0), largest(0) {

        const size_t nsets=sf_list.size();
        std::vector<int> detected(nsets);
        for (auto i : set_id) {
            if (i < 0 || i >= nsets) {
                throw std::runtime_error("size factor set ID out of range");
            }
            detected[i]=1;
        }

        const size_t& ncells=ptr->get_ncol();
        for (size_t i=0; i<nsets; ++i) {
            // We only error check SFs that get used; this allows invalid but 
            // unused SFs to exist, e.g., lib sizes of zero for no-row inputs.
            if (!detected[i]) { 
                continue;
            }

            Rcpp::NumericVector current(sf_list[i]);
            if (current.size() != ncells) { 
                throw std::runtime_error("length of 'size_fac' does not equal number of columns");
            }
            size_factors[i]=current;

            for (auto x : current) {
                if (ISNAN(x) || x<=0) {
                    throw std::runtime_error("size factors should be positive real numbers");
                }
            }
        }

        if (set_id.size()!=ptr->get_nrow()) { 
            throw std::runtime_error("size factor index vector must be equal to number of genes");
        }

        if (subset.size()) { 
            smallest=*std::min_element(subset.begin(), subset.end());
            largest=*std::max_element(subset.begin(), subset.end()) + 1;
        }
        return;
    }

    void get_cell(size_t j, Rcpp::NumericVector::iterator output) {
        // Moving the size factors to their own vector to avoid cache misses.
        for (size_t i=0; i<size_factors.size(); ++i) {
            current_sfs[i]=size_factors[i][j];
        }

        auto inIt=ptr->get_const_col(j, vec.begin(), smallest, largest);
        for (auto s : subset) {
            (*output) = *(inIt + s - smallest) / current_sfs[set_id[s]];
            ++output;
        }

        return;
    }

    size_t get_output_nrow() const {
        return subset.size();
    }
private:        
    M ptr;
    V vec;

    std::vector<Rcpp::NumericVector> size_factors;
    std::vector<double> current_sfs;
    Rcpp::IntegerVector set_id;

    Rcpp::IntegerVector subset;
    size_t smallest, largest;
};

/* Function to compute normalized expression values. */

template <class V, class M>
Rcpp::RObject norm_exprs_internal(M mat, Rcpp::List size_fac_list, Rcpp::IntegerVector sf_to_use, Rcpp::RObject prior_count, Rcpp::RObject log, Rcpp::RObject subset) {
    normalizer<V, M> norm(mat, size_fac_list, sf_to_use, subset);
    const size_t slen=norm.get_output_nrow();
    const size_t ncells=mat->get_ncol();

    // Pulling out the scalars.
    const double prior=check_numeric_scalar(prior_count, "prior count");
    const bool dolog=check_logical_scalar(log, "log specification");

    // Deciding whether or not to preserve sparsity in the output.
    const bool preserve_sparse=(prior==1 || !dolog); 
    beachmat::output_param OPARAM(mat->get_matrix_type(), true, preserve_sparse);
    OPARAM.optimize_chunk_dims(slen, ncells);
    auto optr=beachmat::create_numeric_output(slen, ncells, OPARAM);

    // Computing normalized expression values for each cell (plus a prior, if log-transforming).
    Rcpp::NumericVector current_cell(slen);
    for (size_t c=0; c<ncells; ++c) {
        norm.get_cell(c, current_cell.begin());

        if (dolog) {
            for (auto& val : current_cell) {
                if (!val && preserve_sparse) { // Ensure that it doesn't get turned into some slightly non-zero value.
                    ;
                } else { 
                    val=std::log(val + prior)/M_LN2;
                }
            }
        }
          
        optr->set_col(c, current_cell.begin());
    }

    return optr->yield();
}

SEXP norm_exprs(SEXP counts, SEXP size_fac, SEXP sf_use, SEXP prior_count, SEXP log, SEXP subset) {
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(counts);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(counts);
        return norm_exprs_internal<Rcpp::IntegerVector>(mat.get(), size_fac, sf_use, prior_count, log, subset);
    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(counts);
        return norm_exprs_internal<Rcpp::NumericVector>(mat.get(), size_fac, sf_use, prior_count, log, subset);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}

/* Function to compute average normalized expression values. */

template <class V, class M>
Rcpp::RObject ave_exprs_internal(M mat, Rcpp::List size_fac_list, Rcpp::IntegerVector sf_to_use, Rcpp::RObject subset, Rcpp::IntegerVector firstcell, Rcpp::IntegerVector lastcell) {
    normalizer<V, M> norm(mat, size_fac_list, sf_to_use, subset);
    const size_t slen=norm.get_output_nrow();
    const size_t first=check_integer_scalar(firstcell, "first cell");
    const size_t last=check_integer_scalar(lastcell, "last cell");

    // Summing normalized expression values for each gene (averaging in R, due to parallelization with BiocParallelParam).
    Rcpp::NumericVector current_cell(slen), output(slen);
    for (size_t c=first; c<last; ++c) {
        norm.get_cell(c, current_cell.begin());

        auto oIt=output.begin();
        for (auto ccIt=current_cell.begin(); ccIt!=current_cell.end(); ++ccIt, ++oIt) {
            (*oIt)+=(*ccIt);
        }
    }

    return output;
}

SEXP ave_exprs(SEXP counts, SEXP size_fac, SEXP sf_use, SEXP subset, SEXP first, SEXP last) {
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(counts);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(counts);
        return ave_exprs_internal<Rcpp::IntegerVector>(mat.get(), size_fac, sf_use, subset, first, last);
    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(counts);
        return ave_exprs_internal<Rcpp::NumericVector>(mat.get(), size_fac, sf_use, subset, first, last);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}
