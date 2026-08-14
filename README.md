# This is a series of projects detailing my progress of simulations of general relativistic magnetohydrodynamics (GRMHD).

## Project no.1: The Brio Wu shock tube

### Non technical explanation:
Imagine two sealed rooms filled with air. One room is very hot and under high pressure, while the other is cold and under low pressure, with a thin wall separating them. If that wall suddenly disappears, you might expect the air to calmly mix and balance itself out. Instead, the air reorganizes itself violently. A shock wave pushes into the low-pressure side, while another wave spreads back through the high-pressure side. This is the basic idea behind a classic physics experiment called a shock tube.

The Brio-Wu shock tube takes this scenario and, instead of air, deals with plasma interacting with a magnetic field. Plasma is a state of matter made up of just charged particles floating around, and its interaction with magnetic fields makes its behavior much more complicated. This kind of physics is important to study when dealing with extreme environments in plasma such as stars, accretion disks, and astrophysical jets.

The Brio-Wu test is essentially a standard stress test for computer simulations of plasma. Scientists already have a good understanding of what should happen in this experiment, so if a computer program can reproduce the expected shock waves and other behavior correctly, it gives us confidence that the simulation is working as intended.

My project goes a step beyond simply reproducing the standard test. I am intentionally pushing the simulation into increasingly extreme conditions to find out where the numerical methods begin to fail; more specifically, I can change the density, pressure, velocity, or magnetic field and observe what happens. At some point, the computer may struggle to calculate a physically meaningful solution, or it may fail to find an answer or produce a result that doesn't make physical sense. This brings me to my goal: I want to understand which conditions cause these failures, how extreme the conditions need to become before they happen, and why the numerical method struggles at those points.

I'm currently studying this using Special Relativistic Magnetohydrodynamics (SRMHD), which models magnetized plasma moving at extremely high speeds without including strong gravity. The longer-term goal is to use what I learn here as a foundation for General Relativistic Magnetohydrodynamics (GRMHD), where the simulation also has to account for the extreme gravity found near objects such as black holes.

In simple terms, I'm not just asking “Can the simulation solve this problem?” I'm asking “How far can I push the simulation before it breaks, and what can that teach us about making these simulations more reliable?”
