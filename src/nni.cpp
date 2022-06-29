#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector cwhich(LogicalVector x) {
    IntegerVector idx = seq(0, x.size()-1);
    return idx[x];
}

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

    tree1["tip.label"] = R_NilValue;
    tree2["tip.label"] = R_NilValue;

    List res = List::create(tree1, tree2);
    
    return res;
}