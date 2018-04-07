A fault simulator for TRAX faults and transition faults (TF) using CUDA for Nvidia GPUs.

Copyright (c) 2012-2018, Matthew L. Beckler

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

A good high-level introduction to what this actually is is found in the file [project_report.pdf](./project_report.pdf).

First, circuit designs must be converted into the right "easy" format for use by the GPU fault simulator. Use `python convert_to_easy_format.py c432`, for example. Then you can use `python ref/fault_sim.py c432` to run the very slow reference implementation. Then use `./build.sh` to build the three variants of the GPU code. Then TODO.

For more information about this GPU fault simulator:

 * M. Beckler and R. D. S. Blanton, "GPU-accelerated fault dictionary generation for the TRAX fault model," 2017 International Test Conference in Asia (ITC-Asia), Taipei, 2017, pp. 34-39. - doi: 10.1109/ITC-ASIA.2017.8097107 - https://doi.org/10.1109/ITC-ASIA.2017.8097107

 * M. Beckler and R. D. Blanton, "Fault simulation acceleration for TRAX dictionary construction using GPUs," 2017 IEEE International Test Conference (ITC), Fort Worth, TX, 2017, pp. 1-9. - doi: 10.1109/TEST.2017.8242078 - https://doi.org/10.1109/TEST.2017.8242078

For more information about the TRAX fault model:

 * M. Beckler and R. D. Blanton, "On-Chip Diagnosis of Generalized Delay Failures using Compact Fault Dictionaries," in IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems, vol. PP, no. 99, pp. 1-1. - doi: 10.1109/TCAD.2018.2803621 - https://doi.org/10.1109/TCAD.2018.2803621

 * Beckler, M. and Blanton, R.D., "On-Chip Diagnosis for Early-Life and Wear-Out Failures," IEEE International Test Conference, Nov. 2012. - https://doi.org/10.1109/TEST.2012.6401580

