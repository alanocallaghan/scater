#include "scater.h"

template <typename T, class V, class M>
Rcpp::RObject calc_exprs_internal (M mat, 
        Rcpp::List size_fac_list, Rcpp::IntegerVector sf_to_use, 
        Rcpp::RObject prior_count, Rcpp::RObject log, 
        Rcpp::RObject sum, Rcpp::RObject subset) {

    // Checking dimensions.
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();
    Rcpp::IntegerVector subout=process_subset_vector(subset, mat, true);
    const size_t slen=subout.size();

    // Checking size factors.
    const size_t nsfsets=size_fac_list.size();
    std::vector<Rcpp::NumericVector> size_factors(nsfsets);
    std::vector<Rcpp::NumericVector::iterator> size_factors_it(nsfsets);
    auto sfIt=size_factors.begin();
    auto sfiIt=size_factors_it.begin();
    for (auto sflIt=size_fac_list.begin(); sflIt!=size_fac_list.end(); ++sflIt, ++sfIt, ++sfiIt) {
        (*sfIt)=(*sflIt);
        if (sfIt->size() != ncells) { 
            throw std::runtime_error("length of 'size_fac' does not equal number of columns");
        }
        (*sfiIt)=sfIt->begin();
    }
    
    // Checking size factor choices and subsetting them.
    if (sf_to_use.size()!=ngenes) { 
        throw std::runtime_error("size factor index vector must be equal to number of genes");
    }
    std::vector<size_t> chosen_sf_dex(slen);
    auto csdIt=chosen_sf_dex.begin();
    for (auto s : subout) { 
        const size_t& current=((*csdIt)=sf_to_use[s]-1);
        if (current < 0 || current > nsfsets) { 
            throw std::runtime_error("size factor set index is out of range");
        }
        ++csdIt;
    }

    // Checking scalars.
    if (prior_count.sexp_type()!=REALSXP || LENGTH(prior_count)!=1) { 
        throw std::runtime_error("'prior_count' should be a numeric scalar");
    }
    const double prior=Rcpp::NumericVector(prior_count)[0];
    if (log.sexp_type()!=LGLSXP || LENGTH(log)!=1) {
        throw std::runtime_error("log specification should be a logical scalar"); 
    }
    const bool dolog=Rcpp::LogicalVector(log)[0];
    if (sum.sexp_type()!=LGLSXP || LENGTH(sum)!=1) {
        throw std::runtime_error("sum specification should be a sumical scalar"); 
    }
    const bool dosum=Rcpp::LogicalVector(sum)[0];

    // Setting up output object (complicated setup to avoid initializing the matrix if summing).
    V input(ngenes);
    Rcpp::NumericVector output(slen);
    const bool preserve_sparse=(prior==1 && dolog) || (prior==0 && !dolog); // Deciding whether or not to preserve sparsity.

    std::vector<std::unique_ptr<beachmat::numeric_output> > holder;
    beachmat::numeric_output* optr=NULL;
    if (!dosum) { 
        holder.push_back(
                beachmat::create_numeric_output(slen, ncells,
                    beachmat::output_param(mat->get_matrix_type(), true, preserve_sparse))
                );
        optr=holder.front().get();
    }

    /* Computing normalized expression values for each cell, plus a prior.
     * We may or may not log-transform, and we may or may not sum across genes.
     */
    for (size_t c=0; c<ncells; ++c) {
        auto inIt=mat->get_const_col(c, input.begin());
        auto oIt=output.begin();
        auto csdIt=chosen_sf_dex.begin();

        for (const auto& s : subout) {
            double tmp=*(inIt+s);

            if (!tmp && preserve_sparse) { // Ensure that it doesn't get turned into some slightly non-zero value.
                if (!dosum) {
                    (*oIt)=0;
                }
            } else {            
                tmp/=*size_factors_it[*csdIt];
                tmp+=prior;
                
                if (dosum) { 
                    (*oIt)+=tmp;
                } else if (dolog) { 
                    (*oIt)=std::log(tmp)/M_LN2;
                } else {
                    (*oIt)=tmp;
                }
            }

            ++oIt;
            ++csdIt;
        }
          
        // Incrementing all the size factor iterators.
        for (auto&& sIt : size_factors_it) {
            ++sIt;
        }

        // Adding the result.
        if (!dosum) {
            optr->set_col(c, output.begin());
        }
    }

    // Cleaning up expected output.
    if (dosum) {
        if (dolog) {
            for (auto&& o : output) {
                o=std::log(o)/M_LN2;
            }
        }
        return output;
    } else {
        return optr->yield();
    }
}

SEXP calc_exprs(SEXP counts, SEXP size_fac, SEXP sf_use, SEXP prior_count, SEXP log, SEXP sum, SEXP subset) {
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(counts);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(counts);
        return calc_exprs_internal<int, Rcpp::IntegerVector>(mat.get(), size_fac, sf_use, prior_count, log, sum, subset);
    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(counts);
        return calc_exprs_internal<double, Rcpp::NumericVector>(mat.get(), size_fac, sf_use, prior_count, log, sum, subset);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}
