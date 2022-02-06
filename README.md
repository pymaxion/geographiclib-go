# geographiclib-go
### A Go port of [GeographicLib](https://geographiclib.sourceforge.io/)

# Description
This is a Go implementation of the geodesic algorithms from Charles F. F. Karney's [GeographicLib](https://geographiclib.sourceforge.io/). Though not an official implementation of GeographicLib, geographiclib-go has feature parity with the officially-maintained [Java](https://geographiclib.sourceforge.io/html/java/) and [Python](https://geographiclib.sourceforge.io/html/python/) implementations and additionally includes a utility to validate the implementation against the official 500k-line [GeodTest.dat](find link) repository of geodesic test data.

More information about the GeographicLib project can be found at https://geographiclib.sourceforge.io. This README lifts heavily from GeographicLib's excellent project documentation, especially the section [Geodesics on an Ellipsoid](#geodesics-on-an-ellipsoid), which copies the documentation verbatim to introduce the direct and inverse geodesic problems.
           
# Contents
1. [Setup](#setup)
    1. [Installation](#installation)
    1. [Testing](#testing)
1. [Geodesics on an Ellipsoid](#geodesics-on-an-ellipsoid)
    1. [Introduction](#introduction)
    1. [Additional properties](#additional-properties)
    1. [Multiple shortest geodesics](#multiple-shortest-geodesics)
    1. [Background](#background)
    1. [References](#references)

# Setup
## Installation
< TODO >

## Testing
< TODO >

# Geodesics on an Ellipsoid
## Introduction
Consider a ellipsoid of revolution with equatorial radius <img src="https://render.githubusercontent.com/render/math?math=a"/>, polar semi-axis <img src="https://render.githubusercontent.com/render/math?math=b"/>, and flattening <img src="https://render.githubusercontent.com/render/math?math=f = \\frac{a − b}{a}"/>. Points on the surface of the ellipsoid are characterized by their latitude <img src="https://render.githubusercontent.com/render/math?math=φ"/> and longitude <img src="https://render.githubusercontent.com/render/math?math=λ"/>. (Note that latitude here means the geographical latitude, the angle between the normal to the ellipsoid and the equatorial plane).

The shortest path between two points on the ellipsoid at (φ1, λ1) and (φ2, λ2) is called the geodesic. Its length is s12 and the geodesic from point 1 to point 2 has forward azimuths α1 and α2 at the two end points. In this figure, we have λ12 = λ2 − λ1.

![](img/Geodesic_problem_on_an_ellipsoid.svg)

A geodesic can be extended indefinitely by requiring that any sufficiently small segment is a shortest path; geodesics are also the straightest curves on the surface.

Traditionally two geodesic problems are considered:

* the direct problem — given φ1, λ1, α1, s12, determine φ2, λ2, and α2; this is solved by Geodesic.Direct.
* the inverse problem — given φ1, λ1, φ2, λ2, determine s12, α1, and α2; this is solved by Geodesic.Inverse.

## Additional properties
The routines also calculate several other quantities of interest

* S12 is the area between the geodesic from point 1 to point 2 and the equator; i.e., it is the area, measured counter-clockwise, of the quadrilateral with corners (φ1,λ1), (0,λ1), (0,λ2), and (φ2,λ2). It is given in meters2.
* m12, the reduced length of the geodesic is defined such that if the initial azimuth is perturbed by dα1 (radians) then the second point is displaced by m12 dα1 in the direction perpendicular to the geodesic. m12 is given in meters. On a curved surface the reduced length obeys a symmetry relation, m12 + m21 = 0. On a flat surface, we have m12 = s12. 
* M12 and M21 are geodesic scales. If two geodesics are parallel at point 1 and separated by a small distance dt, then they are separated by a distance M12 dt at point 2. M21 is defined similarly (with the geodesics being parallel to one another at point 2). M12 and M21 are dimensionless quantities. On a flat surface, we have M12 = M21 = 1.
* σ12 is the arc length on the auxiliary sphere. This is a construct for converting the problem to one in spherical trigonometry. The spherical arc length from one equator crossing to the next is always 180°.

If points 1, 2, and 3 lie on a single geodesic, then the following addition rules hold:

* s13 = s12 + s23
* σ13 = σ12 + σ23
* S13 = S12 + S23
* m13 = m12M23 + m23M21
* M13 = M12M23 − (1 − M12M21) m23/m12
* M31 = M32M21 − (1 − M23M32) m12/m23

## Multiple shortest geodesics
The shortest distance found by solving the inverse problem is (obviously) uniquely defined. However, in a few special cases there are multiple azimuths which yield the same shortest distance. Here is a catalog of those cases:

* φ1 = −φ2 (with neither point at a pole). If α1 = α2, the geodesic is unique. Otherwise there are two geodesics and the second one is obtained by setting [α1,α2] ← [α2,α1], [M12,M21] ← [M21,M12], S12 ← −S12. (This occurs when the longitude difference is near ±180° for oblate ellipsoids.)
* λ2 = λ1 ± 180° (with neither point at a pole). If α1 = 0° or ±180°, the geodesic is unique. Otherwise there are two geodesics and the second one is obtained by setting [α1,α2] ← [−α1,−α2], S12 ← −S12. (This occurs when φ2 is near −φ1 for prolate ellipsoids.)
* Points 1 and 2 at opposite poles. There are infinitely many geodesics which can be generated by setting [α1,α2] ← [α1,α2] + [δ,−δ], for arbitrary δ. (For spheres, this prescription applies when points 1 and 2 are antipodal.)
* s12 = 0 (coincident points). There are infinitely many geodesics which can be generated by setting [α1,α2] ← [α1,α2] + [δ,δ], for arbitrary δ.

## Background
The algorithms implemented by this package are given in Karney (2013) and are based on Bessel (1825) and Helmert (1880); the algorithm for areas is based on Danielsen (1989). These improve on the work of Vincenty (1975) in the following respects:

* The results are accurate to round-off for terrestrial ellipsoids (the error in the distance is less than 15 nanometers, compared to 0.1 mm for Vincenty).
* The solution of the inverse problem is always found. (Vincenty’s method fails to converge for nearly antipodal points.)
* The routines calculate differential and integral properties of a geodesic. This allows, for example, the area of a geodesic polygon to be computed.

## References
* F. W. Bessel, [The calculation of longitude and latitude from geodesic measurements (1825)](https://arxiv.org/abs/0908.1824), Astron. Nachr. 331(8), 852–861 (2010), translated by C. F. F. Karney and R. E. Deakin.
* F. R. Helmert, [Mathematical and Physical Theories of Higher Geodesy, Vol 1](https://doi.org/10.5281/zenodo.32050), (Teubner, Leipzig, 1880), Chaps. 5–7.
* T. Vincenty, [Direct and inverse solutions of geodesics on the ellipsoid with application of nested equations](http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf), Survey Review 23(176), 88–93 (1975).
* J. Danielsen, [The area under the geodesic](https://doi.org/10.1179/003962689791474267), Survey Review 30(232), 61–66 (1989).
* C. F. F. Karney, [Algorithms for geodesics](https://doi.org/10.1007/s00190-012-0578-z), J. Geodesy 87(1) 43–55 (2013); [addenda](https://geographiclib.sourceforge.io/geod-addenda.html).
* C. F. F. Karney, [Geodesics on an ellipsoid of revolution](https://arxiv.org/abs/1102.1215v1), Feb. 2011; [errata](https://geographiclib.sourceforge.io/geod-addenda.html#geod-errata).
* [A geodesic bibliography](https://geographiclib.sourceforge.io/geodesic-papers/biblio.html).
* The wikipedia page, [Geodesics on an ellipsoid](https://en.wikipedia.org/wiki/Geodesics_on_an_ellipsoid).
