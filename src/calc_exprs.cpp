#include "scater.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/utils/const_column.h"
#include "utils.h"

#include <vector>
#include <stdexcept>
#include <algorithm>

template <class M>
class normalizer {
public:
    normalizer(M* mat, Rcpp::List sf_list, Rcpp::IntegerVector sf_to_use, Rcpp::RObject genes_sub) : 
        ptr(mat), col_holder(mat, false),
        size_factors(sf_list.size()), current_sfs(sf_list.size()), set_id(sf_to_use), 
        subset(process_subset_vector(genes_sub, mat->get_nrow())), smallest(0), largest(0) 
    {

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

        /* NOTE: very difficult to exploit sparsity due to the need
         * to consider subsets (that may also be duplicated or unordered),
         * hence the false in the const_column constructor.
         */
        col_holder.fill(j, smallest, largest);
        auto it=col_holder.get_values();
        for (auto s : subset) {
            (*output) = *(it + s - smallest) / current_sfs[set_id[s]];
            ++output;
        }

        return;
    }

    size_t get_output_nrow() const {
        return subset.size();
    }
private: 
    M* ptr;
    typename M::vector vec;
    beachmat::const_column<M> col_holder;

    std::vector<Rcpp::NumericVector> size_factors;
    std::vector<double> current_sfs;
    Rcpp::IntegerVector set_id;

    Rcpp::IntegerVector subset;
    size_t smallest, largest;
};

/*****************************************************
 * Function to compute normalized expression values. *
 *****************************************************/

template <class M>
Rcpp::RObject norm_exprs_internal(Rcpp::RObject input, Rcpp::List size_fac_list, Rcpp::IntegerVector sf_to_use, Rcpp::RObject prior_count, Rcpp::RObject log, Rcpp::RObject subset) {
    auto mat=beachmat::create_matrix<M>(input);
    normalizer<M> norm(mat.get(), size_fac_list, sf_to_use, subset);
    const size_t slen=norm.get_output_nrow();
    const size_t ncells=mat->get_ncol();

    // Pulling out the scalars.
    const double prior=check_numeric_scalar(prior_count, "prior count");
    const bool dolog=check_logical_scalar(log, "log specification");

    // Deciding whether or not to preserve sparsity in the output.
    beachmat::output_param OPARAM(mat.get());
    const bool preserve_sparse=(prior==1 || !dolog); 
    if (mat->get_class()=="dgCMatrix" && mat->get_package()=="Matrix" && !preserve_sparse) {
        OPARAM=beachmat::output_param("matrix", "base");
    }
   auto optr=beachmat::create_numeric_output(slen, ncells, OPARAM);

    // Computing normalized expression values for each cell (plus a prior, if log-transforming).
    Rcpp::NumericVector current_cell(slen);
    for (size_t c=0; c<ncells; ++c) {
        norm.get_cell(c, current_cell.begin());

        if (dolog) {
            for (auto& val : current_cell) {
                if (val || !preserve_sparse) { // Ensure that it doesn't get turned into some slightly non-zero value.
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
        return norm_exprs_internal<beachmat::integer_matrix>(counts, size_fac, sf_use, prior_count, log, subset);
    } else if (mattype==REALSXP) {
        return norm_exprs_internal<beachmat::numeric_matrix>(counts, size_fac, sf_use, prior_count, log, subset);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}

/*************************************************************
 * Function to compute average normalized expression values. *
 *************************************************************/

template <class M>
Rcpp::RObject ave_exprs_internal(Rcpp::RObject input, Rcpp::List size_fac_list, Rcpp::IntegerVector sf_to_use, Rcpp::RObject subset) {
    auto mat=beachmat::create_matrix<M>(input);
    normalizer<M> norm(mat.get(), size_fac_list, sf_to_use, subset);
    const size_t slen=norm.get_output_nrow();

    Rcpp::NumericVector current_cell(slen), output(slen);
    for (size_t c=0; c<mat->get_ncol(); ++c) {
        norm.get_cell(c, current_cell.begin());
        auto oIt=output.begin();
        for (auto cc : current_cell) {
            (*oIt)+=cc;
            ++oIt;
        }
    }

    return output;
}

SEXP ave_exprs(SEXP counts, SEXP size_fac, SEXP sf_use, SEXP subset) {
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(counts);
    if (mattype==INTSXP) {
        return ave_exprs_internal<beachmat::integer_matrix>(counts, size_fac, sf_use, subset);
    } else if (mattype==REALSXP) {
        return ave_exprs_internal<beachmat::numeric_matrix>(counts, size_fac, sf_use, subset);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}
