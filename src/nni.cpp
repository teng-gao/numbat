// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// Modified from R-package `ape' by Emmanuel Paradis and Klaus Schliep
static int iii;

void bar_reorderRcpp(int node, int nTips, const arma::Col<int> & e1,
    const arma::Col<int> & e2, std::vector<int> & neworder, const arma::Col<int> & L,
    const arma::Col<int> & xi, const arma::Col<int> & xj)
{
    int i = node - nTips - 1, j, k;

    for (j = xj[i] -1; j >= 0; j--)
        neworder[iii--] = L[xi[i] + j ] + 1;

    for (j = 0; j < xj[i]; j++) {
        k = e2[L[xi[i] + j ]];
        if (k > nTips)
            bar_reorderRcpp(k, nTips, e1, e2, neworder, L, xi, xj);
    }
}

// // [[Rcpp::export]]
arma::Mat<int> reorder_rows(arma::Mat<int> x, arma::Col<int> y) {

    // Create an output matrix
    arma::Mat<int> out = x;

    // Loop through each row and copy the data. 
    for (int i = 0; i < y.n_elem; ++i) {
        out.row(i) = x.row(y[i]-1);
    }

    return out;
}

// [[Rcpp::export]]
arma::Mat<int> reorderRcpp(arma::Mat<int> E) {

    int n = E.n_rows;
    int nTips = n/2 + 1;
    int root = nTips + 1;

    arma::Col<int> e1 = E.col(0);
    arma::Col<int> e2 = E.col(1);
    int m = max(e1), k, j;
    int nnode = m - nTips;
    
    arma::Col<int> L(n);
    std::vector<int> neworder(n);
    arma::Col<int> pos(nnode);
    arma::Col<int> xi(nnode);
    arma::Col<int> xj(nnode);
    for (int i = 0; i < n; i++) {
        xj[e1[i] - nTips - 1]++;
    }
    for (int i = 1; i < nnode; i++) {
        xi[i] = xi[i-1] + xj[i - 1];
    }
    for (int i = 0; i < n; i++) {
        k = e1[i] - nTips - 1;
        j = pos[k]; /* the current 'column' position corresponding to k */
        L[xi[k] + j] = i;
        pos[k]++;
    }

    iii = n - 1;

    bar_reorderRcpp(root, nTips, e1, e2, neworder, L, xi, xj);

    E = reorder_rows(E, neworder);

    return E;
}

// Modified from R-package `phangorn' by Klaus Schliep
// [[Rcpp::export]]
std::vector<arma::Mat<int>> nnin_cpp(const arma::Mat<int> E, const int n) {

    arma::Mat<int> E1 = E;
    arma::Mat<int> E2 = E;
    arma::Col<int> parent = E.col(0);
    arma::Col<int> child = E.col(1);
    int k = min(parent) - 1;
    arma::uvec indvec = find(child > k);
    int ind = indvec[n-1];
    int p1 = parent[ind];
    int p2 = child[ind];
    arma::uvec ind1_vec = find(parent == p1);
    ind1_vec = ind1_vec.elem(find(ind1_vec != ind));
    int ind1 = ind1_vec[0];
    arma::uvec ind2 = find(parent == p2);
    
    int e1 = child[ind1];
    int e2 = child[ind2[0]];
    int e3 = child[ind2[1]];

    E1(ind1, 1) = e2;
    E1(ind2[0], 1) = e1;
    E2(ind1, 1) = e3;
    E2(ind2[1], 1) = e1;

    std::vector<arma::Mat<int>> res(2);

    res[0] = reorderRcpp(E1);
    res[1] = reorderRcpp(E2);

    return res;
}

// Serial version
// [[Rcpp::export]]
List nni_cpp(const List tree) {
    
    arma::Mat<int> E = tree["edge"];
    int Nnode = tree["Nnode"];

    int n = E.n_rows/2 - 1;

    List res(2*n);

    for (int i = 0; i < n; i++) {
        std::vector<arma::Mat<int>> trees = nnin_cpp(E, i+1);
        arma::Mat<int> E1 = trees[0];
        arma::Mat<int> E2 = trees[1];
        
        List tree1 = List::create(Named("edge") = E1, Named("Nnode") = Nnode);
        List tree2 = List::create(Named("edge") = E2, Named("Nnode") = Nnode);

        res[2*i] = tree1;
        res[2*i+1] = tree2;
    }

    return res;

}


// struct get_nni : public Worker {

//     const arma::Mat<int> E;

//     std::vector<arma::Mat<int>> res;

//     // initialize with source and destination
//     get_nni(const arma::Mat<int> E, std::vector<arma::Mat<int>> res): E(E), res(res) {}

//     // take the square root of the range of elements requested
//     void operator()(std::size_t begin, std::size_t end) {
//         for (std::size_t i = begin; i < end; i++) {
//             std::vector<arma::Mat<int>> trees = nnin_cpp(E, i+1);
//             arma::Mat<int> tree1 = trees[0];
//             arma::Mat<int> tree2 = trees[1];
//             res[2*i] = tree1;
//             res[2*i+1] = tree2;
//         }
//     }
// };

// // [[Rcpp::export]]
// std::vector<arma::Mat<int>> nni_cpp_parallel(const List tree) {
    
//     arma::Mat<int> E = tree["edge"];

//     int n = E.n_rows/2 - 1;

//     std::vector<arma::Mat<int>> res(2*n);

//     get_nni get_nni(E, res);

//     parallelFor(0, n, get_nni);

//     return res;

// }
