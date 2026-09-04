# Deadpanorama

Deadpanorama is an ImageJ2/Fiji plugin for 3D segmentation and analysis of
*[deadpan](https://flybase.org/reports/FBgn0010109)*-labelled neuroblast nuclei in confocal image stacks of the
*Drosophila* larval central nervous system. It maps nuclei from a nuclear-marker
channel, reconstructs them in 3D, and measures signal in the image channels for
each reconstructed nucleus.

The Maven coordinates declared by the project are `uk.ac.cam.gurdon:Deadpanorama:0.10.2`.
`pom.xml` and the Java source live under `src/main/java/`, alongside legacy
compiled `.class` files (Java 8) kept for reference; see Installation below for the current, versioned JAR.

## Installation

A pre-built plugin, `target/Deadpanorama-0.10.2.jar` (Java 21; with a
companion `target/Deadpanorama-0.10.2-sources.jar`), is versioned in this
repository and can be installed directly: copy it into `Fiji.app/plugins/`,
or install it through **Plugins > Install…** in Fiji, then skip to step 4
below.

To build from source instead:

1. Install Maven and a JDK.
2. From this directory, run:

   ```bash
   mvn -B clean package
   ```

3. If the build succeeds, copy
   `target/Deadpanorama-0.10.2.jar` into `Fiji.app/plugins/`, or install it
   through **Plugins > Install…** in Fiji.
4. Restart Fiji and run **Plugins > Deadpanorama**.

The POM inherits from `pom-scijava` 13.0.0 and declares ImageJ/SciJava
dependencies.

The plugin is registered through the SciJava `@Plugin` annotation in
[`Deadpanorama.java`](./src/main/java/uk/ac/cam/gurdon/deadpanorama/Deadpanorama.java)
with the menu path `Plugins>Deadpanorama`; no `plugins.config` file is required.

## Input and limitations

The plugin expects an open multi-channel Z-stack in Fiji. One channel is used
for nuclear mapping and the selected analysis operations can measure any valid
image channel, including the mapping channel.

Spatial calibration is used for the Gaussian scale, object area, centroid
coordinates, and volume. The GUI labels these values in µm, µm², and µm³, but
the code uses the units in the ImageJ image calibration. If the image is not
calibrated, the results are pixel-based. If an image has multiple time points,
the implementation reads time point 1.

An optional image ROI restricts mapping to 2D components whose centroid lies
inside the ROI.

## Segmentation pipeline

`Segmenter.getNuclei()` performs the following steps:

1. **3D Difference of Gaussians** on the selected mapping channel. The two
   Gaussian scales are `Sigma` and `Sigma × Sigma Factor`; the X/Y/Z pixel
   calibration is used to convert the physical `Sigma` value to pixel units.
2. **Histogram clustering** on the DoG stack. An iterative K-means-like
   clustering and merging procedure starts from `Starting K` intensity levels
   and returns the levels used for mapping.
3. **Per-slice 2D segmentation** for the selected Z range. Each intensity level
   is thresholded into ROIs; watershed splitting can optionally be applied.
   Slices are processed in parallel.
4. **2D component filtering** by calibrated area, mean intensity normalized by
   the maximum stack intensity, and circularity. Components already captured at
   another intensity level are suppressed.
5. **Optional ROI filtering** by component centroid.
6. **3D reconstruction** by joining components whose Z positions are within
   `Max Slice Distance` and whose projected overlap fraction reaches
   `Min Join Overlap`.
7. **Flatness filtering** after reconstruction. A reconstructed nucleus is
   removed if any of its 2D components has a flatness below `Min Flatness`.
   Flatness is the reciprocal square root of density variation, where density
   variation is variance divided by integrated density.
8. **Overlay and inspection.** Reconstructed nuclei are drawn on the source
   image. **To Cell Counter** can create a Cell Counter image containing the
   reconstructed centroids for manual inspection or correction.
9. **Analysis operations** run on the final nucleus list and write to an ImageJ
   Results table.

## Parameters

Parameters are set in the plugin GUI. Advanced controls are hidden until
**Advanced** is selected. When Advanced is not selected, the following values
are forced to their defaults: Sigma Factor, Starting K, Watershed, Min
Flatness, Max Slice Distance, and Min Join Overlap. Other values remain
editable. The parameter log button writes the current values to the ImageJ log.

| Field | GUI label | Default | Meaning |
|---|---|---:|---|
| `deadpanC` | Deadpan Channel | 1 | Channel used for nucleus mapping |
| `startZ`, `endZ` | Slice Range | 1–100 | Z-slices used for mapping |
| `sigma` | Sigma | 1.0 | Gaussian standard deviation; the GUI treats this as µm |
| `K` | Sigma Factor | 4.0 | Multiplier for the second DoG Gaussian scale |
| `startK` | Starting K | 16 | Initial number of histogram intensity levels |
| `watershed` | Watershed | off | Split thresholded 2D objects before filtering |
| `minA`, `maxA` | Area Range | 10–100 | Minimum and maximum calibrated 2D component area; the GUI displays µm² |
| `threshold` | Threshold | 0.0 | Minimum mean intensity divided by the maximum stack intensity |
| `minCirc` | Min Circularity | 0.6 | Minimum 2D circularity; 1 is a circle |
| `minFlatness` | Min Flatness | 0.2 | Minimum post-reconstruction flatness |
| `joinZ` | Max Slice Distance | 1 | Maximum Z-slice separation for joining components |
| `joinOverlap` | Min Join Overlap | 0.2 | Minimum projected overlap fraction; 1 is complete overlap |

The defaults are loaded and saved through ImageJ `Prefs` keys beginning with
`DeadParams.`. The selected operation types and channel numbers are saved in
`Deadpanorama.cfg` in the user's home directory. The operation list's saved
codes do not include the `outer` checkbox state, so those choices are not
restored.

## Analysis operations

Operations are added or removed in the **Analysis** panel. Each operation adds
columns to the Results table, designated by the zero-based `Deadpan Nucleus` index.
For operations that measure signal, `outer` requests measurements from the
surrounding ring rather than from pixels inside the reconstructed nucleus.

| Operation | Channels | Behavior |
|---|---:|---|
| Statistics | 1 | Reports mean, population standard deviation, minimum, and maximum intensity per nucleus |
| Correlation | 2 | Reports Pearson correlation, Manders overlap coefficient, and Li intensity correlation quotient per nucleus; one `outer` control applies to both channels |
| Scatter Plot | 2 | Plots normalized per-nucleus mean intensities and applies the bundled OPTICS clustering routine to the points |
| Count Positive | 1 | Uses a MaxEntropy threshold on the per-nucleus mean-intensity histogram when the observed range is sufficiently wide; otherwise no nucleus is called positive |
| Signal Ratio | 2 | Reports both channel means and the channel-A/channel-B mean ratio; the current implementation uses channel A's `outer` state for both channels |
| 3D Render | 1 | Creates a composite signal/mask image and renders a rotating 3D projection using the bundled `Projectile` implementation |

## Output

For every reconstructed nucleus, the Results table contains:

| Column | Meaning |
|---|---|
| `Deadpan Nucleus` | Zero-based nucleus index |
| `X`, `Y`, `Z` | Calibrated centroid coordinates |
| `Volume (<unit>³)` | Calibrated volume using the image's current calibration unit |

The Results table is shown with the title `<image title> Deadpanorama` after
all selected operations finish. Operation-specific columns are appended to the
same table. The reconstructed nuclei remain available as an overlay on the
source image.

## Build details

The module is a Maven JAR project under the `pom-scijava` 13.0.0 parent. The
child POM declares these dependencies; versions not shown there are inherited
from the parent:

| Dependency | Role |
|---|---|
| `net.imagej:imagej` | ImageJ2 APIs and runtime components |
| `net.imagej:ij` | ImageJ 1.x compatibility APIs |
| `net.imagej:imagej-legacy` | ImageJ 1.x legacy integration |
| `org.scijava:scijava-common` | SciJava plugin framework |
| `sc.fiji:imagescience:3.0.0` | Image-science routines |
| `org.jfree:jfreechart` | Scatter-plot and histogram charts |


## Bundled third-party code

| Path | Origin and licensing information |
|---|---|
| [`src/main/java/vos/de/kurt/`](./src/main/java/vos/de/kurt/) | Cell Counter code by Kurt De Vos and the Board of Regents of the University of Wisconsin–Madison. |
| [`src/main/java/uk/ac/cam/gurdon/deadpanorama/Projectile.java`](./src/main/java/uk/ac/cam/gurdon/deadpanorama/Projectile.java) | Adapted from the ImageJ Projector implementation and attributed in the source to Pascal code contributed by Michael Castle of the University of Michigan Mental Health Research Institute. |

## License

Deadpanorama is released under the **GNU General Public License v3.0 or later**;
see [`LICENSE`](./LICENSE).

## Citation

For the associated study citation, see the parent repository's
[`README.md`](../README.md#citation).
