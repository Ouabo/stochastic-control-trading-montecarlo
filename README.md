This repository contains a complete R implementation of a stochastic optimal execution model using dynamic programming and Monte Carlo simulation.
The project illustrates how an algorithm can learn when and how much to trade under uncertainty using a fully simulated market environment.

The price process is entirely simulated, allowing full control over volatility, drift, and noise. This makes the framework self-contained, reproducible, and ideal for experimentation or teaching.

The repository includes:

generation of simulated price paths based on a stochastic model

a dynamic-programming solver to compute the optimal execution policy

Monte Carlo routines to approximate conditional expectations

forward simulations showing how the optimal strategy behaves in new scenarios

visualization tools, including static figures and animated GIFs of
price, inventory, actions, and cumulative reward

three example market scenarios (rising, trending, falling) illustrating how
the policy adapts to different volatility patterns

This project serves as an accessible, end-to-end introduction to
stochastic optimal control in finance, showing how trading algorithms
can make optimal decisions in a simulated but realistic environment.
