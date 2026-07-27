# Economic Dispatch with Demand Response in Islanded Microgrids

This repository contains the official MATLAB implementation of my published research on droop parameter-based economic dispatch with incentive-based demand response in islanded microgrids.

## Highlights

- Droop parameter optimization
- Economic load dispatch
- Incentive-based demand response
- Particle Swarm Optimization (PSO)
- Modified Newton-Raphson load flow
- IEEE 33-bus islanded microgrid

## Repository Structure

```
run_economic_dispatch.m     Main program
droop_load_flow.m           Load flow solver
objective_function.m        Objective function
initialize_droop.m          Droop initialization
initialize_customers.m      Customer initialization
initialize_curtailment.m    Curtailment initialization
bus_data.m                  Bus data
line_data.m                 Line data
build_ybus.m                Y-bus formation
polar_to_rect.m             Utility function
```

## Requirements

- MATLAB

## Run

```matlab
run_economic_dispatch
```

## License

MIT License