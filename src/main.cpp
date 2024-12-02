#include <iostream>
#include <vector>
#include <set>
#include <tuple>
#include <cmath>
#include <bitset>
#include <random>
#include <algorithm>
#include <numeric>
#include <functional>


#define vect_b std::vector<bool>
#define vect_vect_b std::vector<std::vector<bool>>

const int SIZE_BITS = 32;

vect_vect_b read_generate_matrix(int k, int n){
    vect_vect_b G(k);
    for (int row = 0; row < k; row++){
        vect_b cur_row(n);
        bool a;
        for (int col = 0; col < n; col++){
            std::cin >> a;
            cur_row[col] = a;
        } 

        G[row] = cur_row;                
    }

    return G;
}

template<typename T>
inline std::vector<T> read_vector(int size) {
    std::vector<T> result(size);
    for (int i = 0; i < size; i++) {
        T word;
        std::cin >> word;
        result[i] = word;
    }

    return result;
}

class Viterbi{
private:

//params matrix
int k_;
int n_;

//Graph Viterbi
std::vector<std::vector<std::pair<bool, int>>> graph_;

//Count nodes in layers
std::vector<int> nodes_;

//Count edges
std::vector<int> edges_;

//Generate and check matrixs
vect_vect_b G_;
vect_vect_b Gs_; // start generate matrix


void set_min_span_form_(){
    //set start stepwise 
    for (int row = 0, col = 0; row < k_ && col < n_; row++, col++){
        int cur_row = row;
        while (cur_row < k_ && col < n_){
            if (G_[cur_row][col]) break;
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
    
    //set end stepwise
    vect_b is_use_row(k_, 0);
    int using_row = 0;
    for (int col = n_-1; col >= 0; col--){
        int cur_col = col;
        int num_cur_using_row = -1;
        for (int row = k_ - 1; row >= 0; row--){
            if (G_[row][col] && !is_use_row[row]){
                if (num_cur_using_row < 0){
                    num_cur_using_row = row;
                    is_use_row[row] = true;
                    using_row++;
                    continue;
                }
                else add_vect_to_vect_(G_[row], G_[num_cur_using_row], n_);
            }
        if (using_row == k_) break;             
        }
    }
     
}

void build_graph_(){
    //find activeting span bounderes
    std::vector<int> new_active(n_, -1); 
    std::vector<int> end_active(n_, -1); 
    for (int row = 0; row < k_; row++) {
        int cur_start = 0;
        int cur_end = n_ - 1;
        while (!G_[row][cur_start]) cur_start++;
        while (!G_[row][cur_end]) cur_end--;
        new_active[cur_start] = row;
        end_active[cur_end] = row;
    }

    //build graph
    std::set<int> past_rows;
    std::set<int> active_rows;
    std::set<int> cur_rows;
    std::set<int> future_rows;
    int numb_edges = 0;
    nodes_.push_back(1);


    for(int lay_num = 0; lay_num < n_; lay_num++){
        //indentify astive, current, past rows
        if (new_active[lay_num] != -1){
            active_rows.insert(new_active[lay_num]);
            future_rows.insert(new_active[lay_num]);
        } 
        if (end_active[lay_num] != -1){
            cur_rows.erase(end_active[lay_num]);
            future_rows.erase(end_active[lay_num]); 
        }    
        
        //count of nodes, edges
        nodes_.push_back(std::pow(2, future_rows.size()));
        numb_edges = std::pow(2, active_rows.size());
        graph_.push_back(std::vector<std::pair<bool, int>>(numb_edges));

        //full node layer weight and end of edges
        for (int num_edge = 0; num_edge < numb_edges; num_edge++){
            bool cur_weight = false;
            int cur_end_edge = get_end_way_numer(num_edge, past_rows, cur_rows, future_rows);
            std::vector<bool> active_rows_val = get_active_vector(num_edge, active_rows);

            //solve current weight
            int num_row_for_arr = 0;
            for (int row: active_rows){
                cur_weight = cur_weight ^ (G_[row][lay_num] && active_rows_val[num_row_for_arr]);
                num_row_for_arr++;               
            }

            //add value of edge in graph
            graph_[(lay_num + 1)][num_edge] = {cur_weight, cur_end_edge};
        }

    //redefining rows for next step
    past_rows = future_rows;
    active_rows = future_rows;
    cur_rows = future_rows;
    }
}

void add_vect_to_vect_(vect_b& changing_vect, vect_b& added_vect, int long_vect){
    for (int i = 0; i < long_vect; i++) changing_vect[i] = changing_vect[i] ^ added_vect[i];
}

int get_end_way_numer(
    int num_edge, std::set<int>& past_rows, std::set<int>& cur_rows, std::set<int>& future_rows){
    if (past_rows.size() == 0) return num_edge;
    if (cur_rows.size() == 0){
        if (future_rows.size() == 0) return 0;
        else return num_edge % 2;
    }

    int cur_mul = 1 << (cur_rows.size() - 1);
    int cur_pos = past_rows.size() - 1;
    std::vector<int> multipliers;
    std::vector<int> bit_pos;
    for (int row: past_rows){
        if (cur_rows.count(row)){
            multipliers.push_back(cur_mul);
            bit_pos.push_back(cur_pos);
            cur_mul >>= 1;
        }         
        cur_pos--;
    }

    int cur_bits_numb;
    if (future_rows.size() > cur_rows.size()) cur_bits_numb = num_edge / 2;
    else cur_bits_numb = num_edge; 

    std::bitset<SIZE_BITS> bit_arr(cur_bits_numb);
    int result = 0;
    int size = bit_pos.size();
    for(int i = 0; i < size; i++){
        result += bit_arr[bit_pos[i]] * multipliers[i];
    } 
    
    if (future_rows.size() > cur_rows.size()) return (result * 2) + (num_edge % 2);
    else return result;
}

std::vector<bool> get_active_vector(int num_edge, std::set<int>& active_rows){
    int size = active_rows.size();
    std::bitset<SIZE_BITS> bit_numb(num_edge);
    std::vector<bool> result;
    
    for (int i = (size - 1); i >= 0; i--) result.push_back(bit_numb[i]);

    return result;
}

public:

Viterbi(int k, int n, vect_vect_b& G){
    k_ = k;
    n_ = n;
    //vect_vect_b new_G{G};    
    G_ = G;
    Gs_ = G_;
    //graph_.emplace_back();
    graph_.emplace_back();

    set_min_span_form_();
    build_graph_();
}

vect_b encode(std::vector<bool> word_in){
    vect_b word_out(n_);
    for (int col = 0; col < n_; col++){
        for (int row = 0; row < k_; row++){
            word_out[col] = word_out[col] ^ (Gs_[row][col] && word_in[row]);
        }
    }

    return word_out;
}

vect_b decode(std::vector<double> word_in){
    //We go through the graph in reverse order.
    //In places where two edges enter one node, we eliminate one path.
    //We fill the vector with tuples with the accumulated error, the current value of the edge, and the position of the passed tuple for a given path.
    //At the end, we collect the vector starting from the end of the vector of tuples.

    std::vector<std::tuple<double, bool, int>> nodes_arr{{0.0, 1, -1}};
    int sp = 0; // start position for nodes array
    
    for (int lay = n_ - 1; lay >= 0; lay--){
        int edges = graph_[lay + 1].size();
        for (int num_edge = 0; num_edge < edges; num_edge++){   
            int past_pos_node = sp + graph_[lay + 1][num_edge].second;         
            bool cur_val = graph_[lay + 1][num_edge].first;
            double cur_err = std::get<0>(nodes_arr[past_pos_node])
             + abs(word_in[lay] - (cur_val ? -1 : 1));            

            if ((nodes_[lay] < nodes_[lay + 1]) && (num_edge % 2)){
                int pos = sp + nodes_[lay + 1] + (num_edge / 2);
                double past_err = std::get<0>(nodes_arr[pos]);
                if (cur_err < past_err){
                    std::get<0>(nodes_arr[pos]) = cur_err;
                    std::get<1>(nodes_arr[pos]) = cur_val;
                    std::get<2>(nodes_arr[pos]) = past_pos_node;
                }
            }
            else nodes_arr.emplace_back(cur_err, cur_val, past_pos_node);
        }

        sp += nodes_[lay + 1];
    }

    int cur = nodes_arr.size() - 1;
    std::vector<bool> word_out(n_);
        for (int i = 0; i < n_; i++) {
            word_out[i] = std::get<1>(nodes_arr[cur]);
            cur = std::get<2>(nodes_arr[cur]);
        }
    return word_out;
}

void print_G_(){
    std::cout<<std::endl;
    for (int row = 0; row < k_; row++){
        for (int col = 0; col < n_; col++) std::cout<<G_[row][col]<<" ";
        std::cout<<std::endl;
    }
}

void print_nodes_(){
    std::cout<<std::endl;
    for (int i = 0; i < nodes_.size(); i++){
        std::cout<<nodes_[i]<<" ";
    }
    std::cout<<std::endl;
}

void print_graph_(){
    std::cout<<std::endl;
    for (int layer = 0; layer < graph_.size(); layer++){
        for (int edge = 0; edge < graph_[layer].size(); edge++){
            std::cout<<"["<<graph_[layer][edge].first<<","<<graph_[layer][edge].second<<"]"<<"  ";
        }
        std::cout<<std::endl;
    }    
}

};

int main(int argc, char* kwarg[]){

    //read data
    std::freopen("input.txt", "r", stdin);
    std::freopen("output.txt", "w", stdout);

    int n, k;
    std::cin>>n;
    std::cin>>k;
    vect_vect_b input_generate_matrix{read_generate_matrix(k,n)};

    //create Viterbi object
    Viterbi code{k, n, input_generate_matrix};
    code.print_nodes_();

    //parametrs for Simulate
    std::random_device rd{};
    std::mt19937 gen{rd()};
    auto gen_b = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());

    //do the command
    std::string command;
    while(std::cin>>command){
        if (command == "Encode"){
            std::vector<bool> word_in = read_vector<bool>(k);
            word_in = code.encode(word_in);
            for (bool i : word_in) std::cout << i << " ";
            std::cout << std::endl;
        }

        if (command == "Decode"){
            std::vector<double> word_in = read_vector<double>(n);
            std::vector<bool> word_out = code.decode(word_in);
            for (bool i : word_out) std::cout << i << " ";
            std::cout << std::endl;
        } 

        if (command == "Simulate"){
            double snrb;
            int num_of_simulations, max_error;
            std::cin >> snrb >> num_of_simulations >> max_error;

            // calculate sigma for normal distribution
            double snr = pow(10, snrb / 10);
            snr = snr * k / n;
            double sigma = sqrt(1.0 / 2.0 / snr);
            std::normal_distribution<> nd{0, sigma};

            int errs = 0;
            int iters = 0;
            for (int i = 0; i < num_of_simulations; i++) {
                std::vector<bool> b = std::vector<bool>(k);
                for (int j = 0; j < k; j++) b[j] = gen_b();
                b = code.encode(b);
                // add noise
                std::vector<double> d;
                transform(b.begin(), b.end(), back_inserter(d),
                          [&nd, &gen](bool bb) -> double {
                              return 1 - 2 * (bb ? 1 : 0) + nd(gen);
                          });
                std::vector<bool> b2 = code.decode(d);
                if (b != b2) errs++;
                iters = i + 1;
                if (errs >= max_error) break;
            }
            std::cout << (double) errs / iters << std::endl;
        } 
    }
}