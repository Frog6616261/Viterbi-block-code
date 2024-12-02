#include <iostream>
#include <vector>

#define vect_b std::vector<bool>
#define vect_vect_b std::vector<std::vector<bool>>

class Viterbi{
private:

//params matrix
int k_;
int n_;

//Graph Viterbi
std::vector<std::vector<int>> graph;

//Generate and check matrixs
vect_vect_b G_;
vect_vect_b Gt_;
vect_vect_b H_;
vect_vect_b Ht_;

std::vector<int> numers_of_NLD_col;



void set_min_span_form(){
    //convert the matrix to a stepwise view and remember not linear dependent columns
    for (int row = 0, col = 0; row < k_ && col < n_; row++, col++){
        int cur_row = row;
        while (cur_row < k_ && col < n_){
            if (G_[cur_row][col]){
                numers_of_NLD_col.push_back(col);
                break;
                }
            else{
                cur_row++;
                if (cur_row == k_){
                    cur_row = row;
                    col++;
                }
            }
        }
        if (col == n_) break;
        if (cur_row != row) add_vect_to_vect_(G_[row], G_[cur_row], n_);
        for (int change_row = row + 1; change_row < k_; change_row++){
            if (G_[change_row][col]) add_vect_to_vect_(G_[change_row], G_[row], n_);
        }              
    }
    //set the unit matrix
    int bound = (k_-1);
    for (int row = k_-2; row >= 0; row--){
        for (int col_num = k_-1; col_num > row; col_num--){
            int col = numers_of_NLD_col[col_num];

            if (G_[row][col]) add_vect_to_vect_(G_[row],G_[col_num], n_);
        }            
    }  
}

void set_H_and_Ht(){

}

void set_graph(){

}

void add_vect_to_vect_(vect_b& changing_vect, vect_b& added_vect, int long_vect){
    for (int i = 0; i < long_vect; i++) changing_vect[i] = changing_vect[i] ^ added_vect[i];
}

public:

Viterbi(int k, int n, vect_vect_b& G){
    k_ = k;
    n_ = n;
    //vect_vect_b new_G{G};    
    G_ = G;

    set_min_span_form();
    set_H_and_Ht();
    set_graph();
}

vect_b decode(std::vector<double> params){
    return vect_b{0};
}

vect_b encode(vect_b word_in){
    return vect_b{0};
}

double simulate(std::vector<double> params){
    return 0.0;
}

void print_G_(){
    std::cout<<std::endl;
    for (int row = 0; row < k_; row++){
        for (int col = 0; col < n_; col++) std::cout<<G_[row][col]<<" ";
        std::cout<<std::endl;
    }
}

void print_LND_vect(){
    std::cout<<std::endl;
    for (int col = 0; col < k_; col++) std::cout<<numers_of_NLD_col[col]<<" ";
    std::cout<<std::endl;
}

};
