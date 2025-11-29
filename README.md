# Droop Parameter-Based Economic Dispatch for Islanded Microgrid Considering Demand Response

This repository provides the full implementation, data, and replication framework for the research work:

> **Droop Parameter-Based Economic Dispatch for Islanded Microgrid Considering Demand Response**  
> Authors: *Ganesh Kumar Pandey*, *Debapriya Das*  
> International Journal of Ambient Energy, 2024  
> DOI: 10.1080/01430750.2024.2305326

---

## 🔍 Overview

This project presents a **day-ahead optimisation** framework for an **islanded AC microgrid** consisting of Droop-controlled Distributed Generators (DGs), renewable units, and incentive-based Demand Response (DR).

The optimisation integrates:

| Component | Purpose |
|----------|---------|
| **P–f & Q–V droop control** | Decentralised sharing of active/reactive power |
| **Economic Load Dispatch (ELD)** | Equal incremental cost principle |
| **Modified Newton–Raphson (MNR) Load Flow** | Frequency-dependent power flow in islanded MG |
| **Particle Swarm Optimisation (PSO)** | Solving multi-objective nonlinear problem |
| **Incentive-based DR** | Peak load reduction with customer utility maximisation |

The result:  
• Reduced operating cost  
• Improved renewable utilisation  
• Better fairness in DG cost allocation  
• Loading relief during peak hours

---

## 🧠 Mathematical Formulation

### 🎯 Objective
Minimise:

\[
\sum_{t=1}^{24}
\left[
\sum_{i=1}^{N_c} (y_i^t - \gamma_i^t x_i^t)
+ \frac{\sum_{j=1}^{N_{DG}} (\lambda_j^t - \bar{\lambda}^t)^2}
{N_{DG}-1}
\right]
\]

Where:  
- \( \lambda_j^t \) = incremental cost of DG \(j\) at hour \(t\)  
- \( x_i^t, y_i^t \) = curtailed power & incentive for customer \(i\)  
- \( \gamma_i^t \) = interruptibility factor  

### 🔒 Constraints
- DG generation limits  
- Droop coefficient bounds  
- Power balance:  
  \[
  \sum P_{DG} = P_{load} - P_{DR} + P_{loss}
  \]
- DR feasibility: rationality, incentive compatibility  
- Utility daily budget \( \leq \$1000 \)

---

## 📂 Repository Structure

