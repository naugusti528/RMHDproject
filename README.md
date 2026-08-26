# General Relativistic Magnetohydrodynamics (GRMHD)

## Project no.1: The Brio Wu shock tube

### No jargon explanation:
Imagine two sealed rooms filled with air. One room is very hot and under high pressure, while the other is cold and under low pressure, with a thin wall separating them. If that wall suddenly disappears, you might expect the air to calmly mix and balance itself out. Instead, the air reorganizes itself violently. A shock wave pushes into the low-pressure side, while another wave spreads back through the high-pressure side. This is the basic idea behind a classic physics experiment called a shock tube.

The Brio-Wu shock tube takes this scenario and, instead of air, deals with plasma interacting with a magnetic field. Plasma is a state of matter made up of just charged particles floating around, and its interaction with magnetic fields makes its behavior much more complicated. This kind of physics is important to study when dealing with extreme environments in plasma such as stars, accretion disks, and astrophysical jets.

I am intentionally pushing the Brio Wu shock tube simulation into increasingly extreme conditions to find out where the numerical methods begin to fail. At some point, the computer may struggle to calculate a physically meaningful solution, which brings me to my goal: I want to understand which conditions cause these failures, how extreme the conditions need to become before they happen, and why the numerical method struggles at those points.

I'm currently studying this using Special Relativistic Magnetohydrodynamics (SRMHD), which models magnetized plasma moving at extremely high speeds without including strong gravity. The longer-term goal is to use what I build here as a foundation for General Relativistic Magnetohydrodynamics (GRMHD), where the simulation also has to account for the extreme gravity found near objects such as black holes. In simple terms, I'm not just asking “Can the simulation solve this problem?” I'm asking “How far can I push the simulation before it breaks, and what can that teach us about making these simulations more reliable?”

### Academic explanation:
The Brio Wu shock tube takes the Sod shock tube problem (a standard 1D hydrodynamic simulation benchmark) and adds magnetism to it. The nonlinearity of fluid dynamics is already complex as is, but adding magnetism turns the problem into a more complex, nonlinear, magnetohydrodynamic (MHD) problem. The software must ensure that the divergence of the magnetic field must strictly remain zero. Additionally, the compound wave is determined specifically by the slow magnetosonic eigenvalue because it's an intermediate-state structure arising from the out-of-phase coupling of thermal and magnetic pressures. Altogether, it means it originates in the center of the tube, so the simulation must produce it accurately without excessive numerical smoothing or hallucinated/unphysical oscillations.

Because there is no strong gravity present, we deal with Minkowski spacetime, which essentially leaves us with special relativistic MHD (SRMHD). I choose to incorporate relativity because down the line I hope to simulate . Examples of these environments include active galactic nuclei (AGNs), specifically their accretion disks and collimated jets, and I mention these phenomena in particular because the gaps in our knowledge of them are what I wish to research in depth. Before I can develop general relativistic MHD (GRMHD) for AGNs, I must first develop functioning and deployable SRMHD code, which is why I'm starting with the Brio Wu shock tube simulation.

The Brio Wu test is vulnerable to breaking from certain values of variables. My current goal is to find out what values of physical properties are needed to make the solver produce unphysical results, because if I'm going to move on to GRMHD, I need to make sure my SRMHD software is robust in extreme environments, like AGNs. As I investigate the fault tolerance of the Brio Wu test, I will also gain insights into the physical meanings of the breaking points, and then how I can guard my software against these breaking points, thereby improving durability for more severe scenarios.

### Update - August 19th 2026
I finished the solver and under standard conditions it produces the shock, compound, and rarefaction waves properly. My next goal is to fault test it with each variable condition present. To have some structure present in how I'm doing this, I'm going to focus on 2 distinct categories of the quantities I'm fault testing: one consists of the physical parameters (i.e. magnetic field, velocity, pressure, plasma beta, etc) and the other consists of numerical parameters (resolution, CFL, Newton-Raphson tolerance, etc). 

### Update - August 26h 2026
I tried to break the solver with increasingly low magnitudes, and nothing happened. I did a sweep where the pressure went across 8 orders of magnitude, and the solver stayed robust and functional.
Important to note though: this was when the velocity was 0. Next, I'm going to keep the beta sweep while implementing a velocity sweep. Maybe once I approach relativistic velocities, the solver might crash.
