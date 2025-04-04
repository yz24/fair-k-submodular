## Fairness in Monotone $k$-submodular Maximization
Submodular optimization has become increasingly prominent in machine learning and fairness has drawn much attention. In this paper, we propose to study the fair $k$-submodular maximization problem and develop a $\frac{1}{3}$-approximation greedy algorithm with a running time of $\mathcal{O}(knB)$. To the best of our knowledge, our work is the first to incorporate fairness in the context of $k$-submodular maximization, and our theoretical guarantee matches the best-known $k$-submodular maximization results without fairness constraints. In addition, we have developed a faster threshold-based algorithm that achieves a $(\frac{1}{3} - \epsilon)$ approximation with $\mathcal{O}(\frac{kn}{\epsilon} \log \frac{B}{\epsilon})$ evaluations of the function $f$. Furthermore, for both algorithms, we provide approximation guarantees when the $k$-submodular function is not accessible but only can be approximately accessed. We have extensively validated our theoretical findings through empirical research and examined the practical implications of fairness. Specifically, we have addressed the question: "What is the price of fairness?" through case studies on influence maximization with $k$ topics and sensor placement with $k$ types. The experimental results show that the fairness constraints do not significantly undermine the quality of solutions.

If you use this code, please cite
```@inproceedings{zhu2024fairness,
  title={Fairness in monotone k-submodular maximization: Algorithms and applications},
  author={Zhu, Yanhui and Basu, Samik and Pavan, A},
  booktitle={2024 IEEE International Conference on Big Data (BigData)},
  pages={173--178},
  year={2024},
  organization={IEEE}
}
```
