This repository goes with the paper:

**Deep Learning Evolved: Overcoming Sub-Optimal Local Minima with $(\mu /  \rho + \lambda)$-Evolution Strategies**

Abstract:
Integrating Evolution Strategies (ES) and Backpropagation (BP) within a deep neural network framework presents a significant challenge, as ES has previously only been shown to perform comparably to BP for smaller problems. In this study, we extend the application of ES to high-dimensional problems, using it to initialize the weights of a Deep Neural Network (DNN). Our experiments demonstrate that this novel ES approach can effectively overcome local minima and converge towards near-optimal global solutions. Following ES initialization, we employ traditional BP gradient methods to further refine the weights based on the initial set provided by ES.
A key finding of our research is the potential for Evolution Strategies to reduce the computational time required by traditional gradient-based backpropagation learning methods. This efficiency is achieved by providing an initial set of weights close to the global optimum, enabling the network to converge more rapidly than with random or zero-weight initialization approaches. Our approach offers a promising direction for future research in the efficient training of deep neural networks and opens up new possibilities for tackling high-dimensional problems with these networks.

```
@inproceedings{rivas2023deep,
  title={Deep Learning Evolved: Overcoming Sub-Optimal Local Minima with $(\mu /  \rho + \lambda)$-Evolution Strategies},
  author={Rivas, Pablo},
  booktitle={The 25th International Conference on Artificial Intelligence (ICAI 2023)},
  year={2023}
}
```
