# STL_STT_CBF: Control Barrier Functions for the Full Class of Signal Temporal Logic Tasks Using Spatiotemporal Tubes_

This repository provides the implementation for synthesizing spatiotemporal tubes (STTs) and designing controllers to satisfy Signal Temporal Logic (STL) specifications. 
This has been presented in the paper "Control Barrier Functions for the Full Class of Signal Temporal Logic Tasks Using Spatiotemporal Tubes"

It includes two case studies:

- **Differential Drive Robot**
- **Drone**

## Repository Structure

### Tube Generation (Python)
- `DifferentialDriveRobot.py`: Generates STTs for the differential drive robot example.
- `Drone.py`: Generates STTs for the drone example.

### Control and Visualization (MATLAB)
- `DifferentialDriveRobot_control.m`: Synthesizes the closed-form controller and simulates the differential drive robot trajectory using the generated tube.
- `Drone_control.m`: Synthesizes the closed-form controller and simulates the drone trajectory using the generated tube.

### Data and Results
- `Drone_1.csv`: Tube and trajectory data for the drone.
- The `Figures` folder contains the `.png` files for the examples.
- The `Compare` folder contains the comparative results for cascaded and fragmented STL specifications, using both our previous and current approaches.

## Requirements

- Python 3.7+ with `z3-solver` installed (`pip install z3-solver`)
- MATLAB R2024B or later (for running `.m` files)

## How to Run

1. **Generate STTs**:
   - Run `DifferentialDriveRobot.py` or `Drone.py` in Python to generate STTs.

2. **Run Controller and Simulation**:
   - Open MATLAB.
   - Run `DifferentialDriveRobot_control.m` for the differential drive robot example or `Drone_control.m` for the drone.
   - This will simulate the system within the tube and visualize the results.

## Citation

If you use this code in your research, please cite the following paper:

```bibtex
@article{STT_STL_arxiv,
  title={Control Barrier Functions for the Full Class of Signal Temporal Logic Tasks Using Spatiotemporal Tubes},
  author={Das, Ratnangshu and Choudhury, Subhodeep and Jagtap, Pushpak},
  journal={arXiv preprint arXiv:2505.05323},
  year={2025}
}