# EM Antenna Field Simulator

An interactive browser-based electromagnetic field simulator for exploring antenna radiation, interference, polarization, and phased-array behavior.

**Live demo:** [www.antennasim.com](https://www.antennasim.com)

The goal of this project is to provide a fast and intuitive way to experiment with antenna configurations and directly visualize how their electromagnetic fields evolve in space and time.

## Features

### Antenna sources

The simulator currently supports several antenna models.

#### Center-fed dipole

A finite straight-wire dipole with a prescribed sinusoidal current distribution.

You can control:

- Current amplitude
- Excitation phase
- Antenna length and orientation through placement
- Feed gap

The inspector also displays the antenna's electrical length, `L / λ`.

#### Hertzian dipole

An ideal infinitesimal current element represented by its current moment, `Iℓ`.

You can control:

- Current moment
- Phase
- Dipole orientation

This is particularly useful for studying elementary dipole radiation, interference between sources, and polarization.

#### Small circular loop

A circular loop with a prescribed uniform circulating current.

You can control:

- Current amplitude
- Phase
- Loop radius

The simulator displays quantities including circumference, `C / λ`, and magnetic moment `IA`.

The uniform-current model is primarily intended for electrically small loops. The interface warns when the loop becomes too large for the small-loop approximation to be appropriate.

#### Linear antenna array

Create a complete linear array as a single editable component.

Array elements can currently be:

- Hertzian dipoles
- Small circular loops

The array parameters include:

- Number of elements
- Element spacing in wavelengths
- Global amplitude
- Global phase
- Progressive phase step
- Array-axis direction
- Dipole orientation, for Hertzian-dipole arrays
- Loop radius, for small-loop arrays

Each element has the same amplitude and a phase given by:

`φₙ = φ₀ + nΔφ`

where `Δφ` is the phase difference between neighboring elements.

This makes it possible to explore:

- Broadside arrays
- End-fire arrays
- Beam steering
- Constructive and destructive interference
- Grating lobes
- The effect of element spacing

The physical spacing automatically follows the wavelength when the simulation frequency or wave speed is changed.

## Field visualization

The fields are visualized in the XY plane containing the antenna geometry.

Three display modes are available.

### Electric field

Displays the instantaneous electric field.

- The heatmap represents field strength
- Arrows indicate the local field direction
- Animation shows the time evolution of the field

### Magnetic field

Displays the magnetic-field component perpendicular to the simulation plane.

- Red and blue indicate opposite field directions
- Intensity represents field magnitude

### Energy flux

Displays an energy-flow visualization based on the local relationship:

`S ∝ E × B`

Arrows show the local direction of energy flow.

This mode should be regarded as an intuitive energy-flux visualization rather than a calibrated power-density measurement.

## Interactive workspace

The simulator uses a component-based interface designed for quickly constructing and modifying antenna configurations.

### Add

Open the component palette and choose an antenna type.

Depending on the component, placement is performed by clicking its position, endpoints, or axis direction.

### Select

Click an antenna in the workspace or select it from the scene list.

Its parameters can then be edited in the **Properties** panel.

### Pan and zoom

The simulation area is not limited to the initially visible region.

- Drag using the Pan tool to move around the scene
- Use the mouse wheel to zoom around the cursor
- Use the zoom controls at the bottom of the interface
- Reset zoom to return to the default view

Zoom changes the view, not the underlying physical dimensions of the simulation.

### Quality

Low, Medium, and High quality modes trade computational cost for spatial and source-discretization resolution.

Higher resolution can be useful when examining smaller geometries or detailed near-field structure.

## Global controls

The top toolbar provides controls for:

- Field visualization mode
- Frequency
- Wave propagation speed
- Simulation quality
- Pause / resume

The wavelength used by the simulation is:

`λ = c / f`

where:

- `λ` is wavelength
- `c` is the selected propagation speed
- `f` is frequency

Changing frequency or wave speed therefore changes the electrical dimensions of antennas and wavelength-relative array spacing.

## Keyboard shortcuts

| Key | Action |
| --- | --- |
| `V` | Select tool |
| `H` | Pan tool |
| `A` | Open Add menu |
| `Space` | Pause / resume |
| `Delete` / `Backspace` | Delete selected component |
| `Esc` | Cancel the current operation |
| `0` | Reset zoom |

Middle-button dragging can also be used for panning.

## How the simulation works

Antenna current distributions are represented as collections of small current elements.

The contributions from these elements are combined as complex phasors using a propagating free-space Green-function model. The resulting vector potential is sampled over the simulation grid, from which the electric and magnetic fields are calculated.

Because the sources are represented by prescribed current distributions, arbitrary combinations of antennas can interfere coherently according to their amplitudes and phases.

The implementation also caches source-field calculations where possible so that operations such as changing excitation phase or amplitude, panning, and revisiting previously calculated regions can be handled efficiently.

## Physical scope and limitations

This project is intended primarily as an **educational and visualization tool**.

There are several important limitations:

- Antenna currents are prescribed rather than solved self-consistently.
- Mutual coupling between antennas or array elements is currently not modeled.
- Feed impedance, matching, input power, and radiation resistance are not calculated.
- Passive conductors and PEC boundaries are not currently modeled.
- The small-loop model assumes a uniform circulating current.
- The center-fed dipole uses an assumed sinusoidal current distribution.
- A small source-distance smoothing term is used to avoid singular behavior, so fields extremely close to a source should be interpreted qualitatively.
- The simulation displays a two-dimensional slice of the electromagnetic field rather than a full 3D visualization.
- It is not a replacement for full-wave methods such as FDTD, FEM, MoM, or commercial electromagnetic solvers.

These simplifications make it possible to interact with antenna configurations in real time while retaining many of the important wave phenomena the simulator is intended to demonstrate.

## Things to explore

Some interesting experiments include:

- Change the length of a dipole relative to wavelength
- Place two dipoles with different phase offsets
- Create circular or elliptical polarization using perpendicular Hertzian dipoles
- Compare electric-dipole and small-loop fields
- Build broadside and end-fire arrays
- Vary array spacing and observe grating lobes
- Sweep the progressive phase of an array to steer its beam
- Compare near-field and far-field behavior
- Observe interference between multiple independent antennas

For example, for a linear array with spacing `d`, the progressive phase can be varied to steer the direction in which contributions from neighboring elements add constructively.

## Technology

The simulator runs entirely in the browser and is implemented in JavaScript using **p5.js** for rendering and interaction.

No server-side computation is required.

## Project status

The simulator is actively being developed.

Possible future extensions include:

- Additional antenna models
- Passive conducting structures and PEC boundaries
- Additional array geometries
- More advanced field-analysis tools
- Improved far-field and radiation-pattern visualization

Feedback, bug reports, and suggestions are welcome.
