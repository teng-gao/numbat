// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Modified from R-package `ape' by Emmanuel Paradis and Klaus Schliep
static int iii;

void bar_reorderRcpp(int node, int nTips, const IntegerVector & e1,
    const IntegerVector & e2, IntegerVector neworder, const IntegerVector & L,
    const IntegerVector & xi, const IntegerVector & xj)
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

// [[Rcpp::export]]
IntegerMatrix reorder_rows(IntegerMatrix x, IntegerVector y) {

    // Create an output matrix
    IntegerMatrix out(x.nrow(), x.ncol());

    // Loop through each row and copy the data. 
    for (int i = 0; i < y.length(); ++i) {
        out(i,_) = x(y[i]-1,_);
    }

    return out;
}

// [[Rcpp::export]]
IntegerMatrix reorderRcpp(IntegerMatrix E) {

    int n = E.nrow();
    int nTips = n/2 + 1;
    int root = nTips + 1;

    IntegerVector e1 = E( _, 0);
    IntegerVector e2 = E( _, 1);
    int m = max(e1), k, j;
    int nnode = m - nTips;
    
    IntegerVector L(n);
    IntegerVector neworder(n);
    IntegerVector pos(nnode);
    IntegerVector xi(nnode);
    IntegerVector xj(nnode);
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

IntegerVector cwhich(LogicalVector x) {
    IntegerVector idx = seq(0, x.size()-1);
    return idx[x];
}

// Modified from R-package `phangorn' by Klaus Schliep
// [[Rcpp::export]]
List nnin_cpp(const List tree, const int n) {

    List tree1 = clone(tree);
    List tree2 = clone(tree);
    IntegerMatrix edge = tree["edge"];
    IntegerVector parent = edge(_, 0);
    IntegerVector child = edge(_, 1);
    int k = min(parent) - 1;
    int ind = cwhich(wrap(child > k))[n-1];
    int p1 = parent[ind];
    int p2 = child[ind];
    IntegerVector ind1_vec = cwhich(wrap(parent == p1));
    ind1_vec = ind1_vec[ind1_vec != ind];
    int ind1 = ind1_vec[0];
    IntegerVector ind2 = cwhich(wrap(parent == p2));
    int e1 = child[ind1];
    int e2 = child[ind2[0]];
    int e3 = child[ind2[1]];

    IntegerMatrix edge1 = tree1["edge"];
    edge1(ind1, 1) = e2;
    edge1(ind2[0], 1) = e1;

    IntegerMatrix edge2 = tree2["edge"];
    edge2(ind1, 1) = e3;
    edge2(ind2[1], 1) = e1;

    tree1["edge"] = reorderRcpp(edge1);
    tree2["edge"] = reorderRcpp(edge2);

    tree1["tip.label"] = R_NilValue;
    tree2["tip.label"] = R_NilValue;

    List res = List::create(tree1, tree2);
    
    return res;
}

// [[Rcpp::export]]
std::vector<List> nni_cpp(const List tree) {
    
    IntegerMatrix E = tree["edge"];

    int n = E.nrow()/2 - 1;

    std::vector<List> res(2*n);

    for (int i = 0; i < n; i++) {
        List trees = nnin_cpp(tree, i+1);
        res[2*i] = trees[0];
        res[2*i+1] = trees[1];
    }

    return res;

}