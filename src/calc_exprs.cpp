#include "scater.h"

template <class V, class M>
class normalizer {
public:
    normalizer(M mat, Rcpp::List sf_list, Rcpp::IntegerVector sf_to_use, Rcpp::RObject genes_sub) : 
            ptr(mat), vec(mat->get_nrow()), 
            size_factors(sf_list.size()), current_sfs(sf_list.size()), set_id(sf_to_use), 
            subset(process_subset_vector(genes_sub, mat->get_nrow())), output(subset.size()) {

        const size_t& ncells=ptr->get_ncol();
        for (size_t i=0; i<sf_list.size(); ++i) {
            Rcpp::NumericVector current(sf_list[i]);
            if (current.size() != ncells) { 
                throw std::runtime_error("length of 'size_fac' does not equal number of columns");
            }
            size_factors[i]=current;
        }

        if (set_id.size()!=ptr->get_nrow()) { 
            throw std::runtime_error("size factor index vector must be equal to number of genes");
        }
        return;
    }

    Rcpp::NumericVector& get_cell(size_t j) {
        // Moving the size factors to their own vector to avoid cache misses.
        for (size_t i=0; i<size_factors.size(); ++i) {
            current_sfs[i]=size_factors[i][j];
        }

        auto inIt=ptr->get_const_col(j, vec.begin());
        auto oIt=output.begin();
        for (auto s : subset) {
            (*oIt) = *(inIt + s) / current_sfs[set_id[s]];
            ++oIt;
        }

        return output;
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
    Rcpp::NumericVector output;
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
    for (size_t c=0; c<ncells; ++c) {
        auto current_cell=norm.get_cell(c);

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
Rcpp::RObject ave_exprs_internal(M mat, Rcpp::List size_fac_list, Rcpp::IntegerVector sf_to_use, Rcpp::RObject subset) {
    normalizer<V, M> norm(mat, size_fac_list, sf_to_use, subset);
    const size_t slen=norm.get_output_nrow();
    const size_t ncells=mat->get_ncol();

    // Averaging normalized expression values for each gene.
    Rcpp::NumericVector output(slen);
    for (size_t c=0; c<ncells; ++c) {
        auto current_cell=norm.get_cell(c);

        auto oIt=output.begin();
        for (auto ccIt=current_cell.begin(); ccIt!=current_cell.end(); ++ccIt, ++oIt) {
            (*oIt)+=(*ccIt);
        }
    }

    for (auto& o : output) {
        o/=ncells;
    }
    return output;
}

SEXP ave_exprs(SEXP counts, SEXP size_fac, SEXP sf_use, SEXP subset) {
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(counts);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(counts);
        return ave_exprs_internal<Rcpp::IntegerVector>(mat.get(), size_fac, sf_use, subset);
    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(counts);
        return ave_exprs_internal<Rcpp::NumericVector>(mat.get(), size_fac, sf_use, subset);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}
