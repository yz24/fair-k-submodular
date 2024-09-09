// Copyright 2020, Grigorios Loukides
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Takuya Akiba nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <time.h> // clock
#include <math.h>
#include <random>
#include <memory> // auto_ptr
#include "jlog.h"
#include "tools.h"

using namespace std;
using namespace jlog_internal;
//typedef long long LL;

int n_;
int k_;
int tick;


double entropy(vector<vector<vector<int> > > &vals, vector<pair<int, int> > &S, int n_sensors, int n_ticks, int n_pos, double delta) {
	map<vector<int>, int> C;
	for (int t = 0; t < n_ticks; t++) {
		vector<int> x;
		for (auto item:S) {
			int p = item.first, s = item.second;
			x.push_back(vals[s][t][p]);
		}
		C[x]++;
	}
	double H = 0;
	for (auto it = C.begin(); it != C.end(); it++) {
		double pr = 1.0 * it->second / n_ticks;
		H += -pr * log(pr);
	}

    int random = rand()%2000;
    double div = H * delta * 1000 / (random-1000);
    double H1 = H + div;
	return H1;
}

// [s][t][p] = val
// <pos,sensor>
double mutual_information(vector<vector<vector<int> > > &vals,
		vector<pair<int, int> > &S, int n_sensors, int n_ticks, int n_pos, double delta) {
	// <p,s>
	// U = S \cup O
	set<pair<int, int> > O, U;
	set<pair<int, int> > Sset(S.begin(), S.end());
	for (int p = 0; p < n_pos; p++) {
		for (int s = 0; s < n_sensors; s++) {
			pair<int, int> ps = make_pair(p, s);
			if (!Sset.count(ps)) {
				O.insert(ps);
			}
			U.insert(ps);
		}
	}

	return entropy(vals, S, n_sensors, n_ticks, n_pos, delta);
}


void lazy_greedy_( vector<vector<vector<int> > > &es, int budget, double delta) {
    cout<<"\n----------TS greedy with B="<<budget<<"----------"<<endl;
    bool TS = true;
    // < gain, < <node,item>, tick > >
    priority_queue<pair<double, pair<pair<int, int>, int> > > que;
    for (int v = 0; v < n_; v++) {
        for (int z = 0; z < k_; z++) {
            que.push(make_pair(1e12, make_pair(make_pair(v, z), -1)));
        }
    }
    if(!TS){
        vector<int> budgets(k_);
//    int B = 0;
        for (int j = 0; j < budgets.size(); j++) {
            budgets[j] = budget;
        }
    }


    double current_f_value = 0.0;
    double candidate_f_value = 0.0;
    int num_evaluations = 0;
    vector<bool> used(n_);
    vector<pair<int, int> > S;
    for (int j = 0; j < budget; j++) {
//        printf(INFO "|S| = %d" DEF, j+1);
        double gain;
        pair<int, int> next;
        double start = get_current_time_sec();
        for (int num = 0;; num++) {
            double now = get_current_time_sec();

            pair<double, pair<pair<int, int>, int> > pp = que.top();
            que.pop();
            pair<int, int> s = pp.second.first;
            int last = pp.second.second;
//            if (used[s.first] || (!TS && budgets[s.second] == 0)) {
            if (used[s.first]) {
                continue;
            }

            if (last == j) {
                next = s;
                gain = pp.first;
                JLOG_ADD("num-of-eval", num);
                num_evaluations += num;
                break;
            }
            vector<pair<int, int> > SS(S);
            SS.push_back(s);
//            candidate_f_value = simulate(es, SS, beta);
            candidate_f_value = mutual_information(es, SS, k_, tick, n_, delta);

            que.push(make_pair(candidate_f_value - current_f_value, make_pair(s, j)));

        }
        current_f_value = candidate_f_value;
        S.push_back(next);
        used[next.first] = true;
//        if(!TS) budgets[next.second]--;

        if((j+1)%k_==0){
            printf(INFO "|S| = %d" DEF, j+1);
            JLOG_ADD("seed", next.first);
            JLOG_ADD("item", next.second);
            JLOG_ADD("gain", gain);
            JLOG_ADD("objective", current_f_value);
            JLOG_ADD("num-of-evaluation as of now", num_evaluations);
            cout<<endl;

            cout<<"seed set: ";
            vector<int> stats(k_, 0);
            for(auto s: S){
                stats[s.second]++;
                cout<<"<"<<s.first<<", "<<s.second<<">, ";
            }
            cout<<endl;
            for(auto st: stats){
                cout<< st<<" ";
            }
            cout<<endl;
            double f_value = mutual_information(es, S, k_, tick, n_, 0);
            cout<<f_value<<" "<<num_evaluations<<endl;
        }
    }

    cout<<"seed set: ";
    vector<int> stats(k_, 0);
    for(auto s: S){
        stats[s.second]++;
        cout<<"<"<<s.first<<", "<<s.second<<">, ";
    }
    cout<<endl;
    for(auto st: stats){
        cout<< st<<" ";
    }
    cout<<endl;
    double f_value = mutual_information(es, S, k_, tick, n_, 0);
    cout<<f_value<<" "<<num_evaluations<<endl;

}


void lazy_greedy_IS_( vector<vector<vector<int> > > &es, vector<int> uppers, int budget, double delta) {
    cout<<"\n----------IS greedy with B="<<budget<<"----------"<<endl;
    bool TS = true;
    // < gain, < <node,item>, tick > >
    priority_queue<pair<double, pair<pair<int, int>, int> > > que;
    for (int v = 0; v < n_; v++) {
        for (int z = 0; z < k_; z++) {
            que.push(make_pair(1e12, make_pair(make_pair(v, z), -1)));
        }
    }
    vector<int> budgets = std::move(uppers);


    double current_f_value = 0.0;
    double candidate_f_value = 0.0;
    int num_evaluations = 0;
    vector<bool> used(n_);
    vector<pair<int, int> > S;
    for (int j = 0; j < budget; j++) {
//        printf(INFO "|S| = %d" DEF, j+1);
        double gain;
        pair<int, int> next;
        double start = get_current_time_sec();
        for (int num = 0;; num++) {
            double now = get_current_time_sec();

            pair<double, pair<pair<int, int>, int> > pp = que.top();
            que.pop();
            pair<int, int> s = pp.second.first;
            int last = pp.second.second;
//            if (used[s.first] || (!TS && budgets[s.second] == 0)) {
            if (used[s.first] || budgets[s.second]<=0) {
                continue;
            }

            if (last == j) {
                next = s;
                gain = pp.first;
                JLOG_ADD("num-of-eval", num);
                num_evaluations += num;
                break;
            }
            vector<pair<int, int> > SS(S);
            SS.push_back(s);
//            candidate_f_value = simulate(es, SS, beta);
            candidate_f_value = mutual_information(es, SS, k_, tick, n_, delta);

            que.push(make_pair(candidate_f_value - current_f_value, make_pair(s, j)));

        }
        current_f_value = candidate_f_value;
        S.push_back(next);
        used[next.first] = true;
        budgets[next.second]--;
//        if(!TS) budgets[next.second]--;

        if((j+1)%k_==0){
            printf(INFO "|S| = %d" DEF, j+1);
            JLOG_ADD("seed", next.first);
            JLOG_ADD("item", next.second);
            JLOG_ADD("gain", gain);
            JLOG_ADD("objective", current_f_value);
            JLOG_ADD("num-of-evaluation as of now", num_evaluations);
            cout<<endl;

            cout<<"seed set: ";
            vector<int> stats(k_, 0);
            for(auto s: S){
                stats[s.second]++;
                cout<<"<"<<s.first<<", "<<s.second<<">, ";
            }
            cout<<endl;
            for(auto st: stats){
                cout<< st<<" ";
            }
            cout<<endl;
            double f_value = mutual_information(es, S, k_, tick, n_, 0);
            cout<<f_value<<" "<<num_evaluations<<endl;
        }
    }

    cout<<"seed set: ";
    vector<int> stats(k_, 0);
    for(auto s: S){
        stats[s.second]++;
        cout<<"<"<<s.first<<", "<<s.second<<">, ";
    }
    cout<<endl;
    for(auto st: stats){
        cout<< st<<" ";
    }
    cout<<endl;
    double f_value = mutual_information(es, S, k_, tick, n_, 0);
    cout<<f_value<<" "<<num_evaluations<<endl;

}



bool extendable(vector<pair<int, int> > S, pair<int, int> s, vector<int> lowers, vector<int> uppers, int B){
    vector<int> supp(k_, 0);
    supp[s.second]++;
    for(auto p : S){
        supp[p.second]++;
    }
    // check upper bounds, only check the s's type
    if(supp[s.second] > uppers[s.second]){
        return false;
    }
    // check lower bounds.
    int sum = 0;
    for(int i = 0; i<k_; i++){
        sum += max(lowers[i], supp[i]);
    }
    if(sum > B){
        return false;
    }
    return true;
}


// [u][i] = <v, p[1..k]>
void fair_greedy_(vector<vector<vector<int >> > &es, vector<int> lowers, vector<int> uppers, int B, double delta) {
    cout<<"\n---------- fair greedy with B="<<B<<"----------"<<endl;
    // < gain, < <node,item>, tick > >
    priority_queue<pair<double, pair<pair<int, int>, int> > > que;
    for (int v = 0; v < n_; v++) {
        for (int z = 0; z < k_; z++) {
            que.push(make_pair(1e12, make_pair(make_pair(v, z), -1)));
        }
    }
//    vector<int> budgets(k);
//    int B = 0;
//    for (int j = 0; j < budgets.size(); j++) {
//        budgets[j] = budget;
//    }

    double current_f_value = 0.0;
    double candidate_f_value = 0.0;
    int num_evaluations = 0;

    vector<bool> used(n_);
    vector<pair<int, int> > S;
    for (int j = 0; j < B; j++) {
        pair<double, pair<pair<int, int>, int> > pp = que.top();
//        if(pp.first < 0){
//            break;
//        }
        printf(INFO "|S| = %d" DEF, j+1);
        double gain;
        pair<int, int> next;
        for (int num = 0;; num++) {

            pair<double, pair<pair<int, int>, int> > pp = que.top();
            que.pop();
            pair<int, int> s = pp.second.first;
            int last = pp.second.second;

            if (used[s.first] || !extendable(S, s, lowers, uppers, B)) {
//                if(budgets[s.second] == 0)
//                    cout<<"type of "<<s.second<<" used up."<<endl;
                continue;
            }

            if (last == j) {
                next = s;
                gain = pp.first;
                JLOG_ADD("num-of-eval", num);
                num_evaluations += num;
                break;
            }
            vector<pair<int, int> > SS(S);
            SS.push_back(s);
//			double sigma = simulate(n, es, k, SS, beta);
//			double psigma = simulate(n, es, k, S, beta);
//            candidate_f_value = simulate(es, SS, beta);
            candidate_f_value = mutual_information(es, SS, k_, tick, n_, delta);

            que.push(make_pair(candidate_f_value - current_f_value, make_pair(s, j)));

        }
        current_f_value = candidate_f_value;
        S.push_back(next);
        used[next.first] = true;
//        if(!TS) budgets[next.second]--;

        if((j+1)%k_==0){
            printf(INFO "|S| = %d" DEF, j+1);
            JLOG_ADD("seed", next.first);
            JLOG_ADD("item", next.second);
            JLOG_ADD("gain", gain);
            JLOG_ADD("objective", current_f_value);
            JLOG_ADD("num-of-evaluation as of now", num_evaluations);
            cout<<endl;
        }
    }
    cout<<"seed set: ";
    vector<int> stats(k_,0);
    for(auto s: S){
        stats[s.second]++;
        cout<<"<"<<s.first<<", "<<s.second<<">, ";
    }
    cout<<endl;
    for(auto st: stats){
        cout<< st<<" ";
    }
    cout<<endl;
    double f_value =  mutual_information(es, S, k_, tick, n_, 0);
    cout<<f_value<<" "<<num_evaluations<<endl;

}

// [u][i] = <v, p[1..k]>
void fair_threshold_(vector<vector<vector<int > > > &es, vector<int> lowers, vector<int> uppers, int B, double eps, double delta) {
    cout<<"\n---------- fair threshold with eps="<<eps<<", B="<<B<<"----------"<<endl;
    // < gain, < <node,item>, tick > >
    map<pair<int, int>, int> evaluations;

    priority_queue<pair<double, pair<int, int>>> que;
    for (int v = 0; v < n_; v++) {
        for (int z = 0; z < k_; z++) {
            que.push(make_pair(1e12, make_pair(v, z)));
            evaluations[{v, z}] = 0;
        }
    }

    double current_f_value = 0.0;
    double candidate_f_value = 0.0;
    int num_evaluations = 0;

    vector<bool> used(n_);
    vector<pair<int, int> > S;
    for (int j = 0; j < B; j++) {

        if(que.empty()){
            break;
        }

        pair<int, int> next;
        double gain;
        int num = 0;
        while (true) {
            pair<double, pair<int, int> > pp = que.top();
            que.pop();
            pair<int, int> s = pp.second;

            if (used[s.first] || !extendable(S, s, lowers, uppers, B)) {
                continue;
            }

            vector<pair<int, int> > SS(S);
            SS.push_back(s);
//            candidate_f_value = simulate(es, SS, beta);
            candidate_f_value = mutual_information(es, SS, k_, tick, n_, delta);
            num ++;
            evaluations[s] ++;
            gain = candidate_f_value - current_f_value;

            if(gain >= que.top().first || gain >= pp.first * (1-eps)){
                S.push_back(s);
                next = s;
                used[s.first] = true;
                num_evaluations += num;

                break;
            }

            if(evaluations[s] <= int(log(B/(2*eps)))/eps){
                que.push(make_pair(gain, s));
            }
//          que.push(make_pair(candidate_f_value - current_f_value, make_pair(s, j)));

        }
        current_f_value = candidate_f_value;
        if((j+1)%k_==0){
            printf(INFO "|S| = %d" DEF, j+1);
            JLOG_ADD("seed", next.first);
            JLOG_ADD("item", next.second);
            JLOG_ADD("gain", gain);
            JLOG_ADD("objective", current_f_value);
            JLOG_ADD("num-of-evaluation as of now", num_evaluations);
            cout<<endl;
        }

    }
    cout<<"seed set: ";
    vector<int> stats(k_, 0);
    for(auto s: S){
        stats[s.second]++;
        cout<<"<"<<s.first<<", "<<s.second<<">, ";
    }
    cout<<endl;
    for(auto st: stats){
        cout<< st<<" ";
    }
    cout<<endl;
    double f_value = mutual_information(es, S, k_, tick, n_, 0);
    cout<<f_value<<", "<<num_evaluations<<endl;

}


vector<pair<int, int>> generate_random(int B, vector<int> lowers, vector<int> uppers){
//    vector<pair<int, int>> results(k);
    set<pair<int, int>> uniquePairs;
    vector<int> vec(n_);
    vector<int> types(k_);
    std::iota (std::begin(vec), std::end(vec), 0);
    std::iota (std::begin(types), std::end(types), 0);
    vector<int> type_count(k_, 0);
    // Generate B unique pairs

    for(int i=0; i<k_; i++){
        while(type_count[i]<lowers[i]){
            int randomVec = rand() % vec.size();
            if(uniquePairs.count({vec[randomVec], i}) > 0){
                continue;
            }
            uniquePairs.insert({vec[randomVec], i});
//            cout<<uniquePairs.size()<<": "<<vec[randomVec]<<" "<< i<<endl;

            type_count[i]++;
        }
    }

    int remaining = B - std::accumulate(lowers.begin(), lowers.end(), 0);
    while(remaining > 0){
        int randomVec = rand() % vec.size();
        int randomType = rand() % types.size();
        if(uniquePairs.count({vec[randomVec], types[randomType]}) > 0){
            continue;
        }
        uniquePairs.insert({vec[randomVec], types[randomType]});
//        cout<<uniquePairs.size()<<": "<<vec[randomVec]<<" "<< types[randomType]<<endl;
        remaining --;
        type_count[randomType] ++;
    }

    vector<pair<int, int>> results(uniquePairs.begin(), uniquePairs.end());

    return results;
}

vector<pair<int, int>> generate_random(int B){
//    vector<pair<int, int>> results(k);
    set<pair<int, int>> uniquePairs;
    vector<int> vec(n_);
    vector<int> types(k_);
    std::iota (std::begin(vec), std::end(vec), 0);
    std::iota (std::begin(types), std::end(types), 0);

    int c = 0;
    // Generate B unique pairs
    while (uniquePairs.size() < B) {
        int randomVec = rand() % vec.size();
        int randomType = rand() % types.size();
        uniquePairs.insert({vec[randomVec], types[randomType]});
        cout<<c<<": "<<uniquePairs.size()<<", "<<vec[randomVec]<<" "<< types[randomType]<<endl;
        c++;
    }
    vector<pair<int, int>> results(uniquePairs.begin(), uniquePairs.end());

    return results;
}



