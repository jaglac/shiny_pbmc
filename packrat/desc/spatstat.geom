Package: spatstat.geom
Version: 3.0-6
Date: 2023-01-30
Title: Geometrical Functionality of the 'spatstat' Family
Authors@R: c(person("Adrian", "Baddeley", 
                    role = c("aut", "cre", "cph"),
       	            email = "Adrian.Baddeley@curtin.edu.au",
		    comment = c(ORCID="0000-0001-9499-8382")),
	     person("Rolf", "Turner", 
                    role = c("aut", "cph"),
 	            email="r.turner@auckland.ac.nz",
		    comment=c(ORCID="0000-0001-5521-5218")),
	     person("Ege",   "Rubak", 
                    role = c("aut", "cph"),
		    email = "rubak@math.aau.dk",
		    comment=c(ORCID="0000-0002-6675-533X")),
	     person("Tilman",    "Davies",        role = "ctb"),
	     person("Ute",       "Hahn",          role = "ctb"),
	     person("Abdollah",  "Jalilian",      role = "ctb"),
 	     person("Greg",      "McSwiggan",     role = c("ctb", "cph")),
 	     person("Sebastian", "Meyer",         role = c("ctb", "cph")),
	     person("Jens",      "Oehlschlaegel", role = c("ctb", "cph")),
	     person("Suman",     "Rakshit",       role = "ctb"),
	     person("Dominic",   "Schuhmacher",   role = "ctb"),
	     person("Rasmus",    "Waagepetersen", role = "ctb"))
Maintainer: Adrian Baddeley <Adrian.Baddeley@curtin.edu.au>
Depends: R (>= 3.5.0), spatstat.data (>= 3.0), stats, graphics,
        grDevices, utils, methods
Imports: spatstat.utils (>= 3.0), deldir (>= 1.0-2), polyclip (>=
        1.10-0)
Suggests: spatstat.random, spatstat.explore, spatstat.model,
        spatstat.linnet, maptools (>= 0.9-9), spatial, fftwtools (>=
        0.9-8), spatstat (>= 2.0-0)
Additional_repositories: https://spatstat.r-universe.dev
Description: Defines spatial data types and supports geometrical operations
	     on them. Data types include point patterns, windows (domains),
	     pixel images, line segment patterns, tessellations and hyperframes.
	     Capabilities include creation and manipulation of data
	     (using command line or graphical interaction),
	     plotting, geometrical operations (rotation, shift, rescale,
	     affine transformation), convex hull, discretisation and
	     pixellation, Dirichlet tessellation, Delaunay triangulation,
	     pairwise distances, nearest-neighbour distances,
	     distance transform, morphological operations
	     (erosion, dilation, closing, opening), quadrat counting,
	     geometrical measurement, geometrical covariance,
	     colour maps, calculus on spatial domains,
	     Gaussian blur, level sets of images, transects of images,
	     intersections between objects, minimum distance matching.
	     (Excludes spatial data on a network, which are supported by
	     the package 'spatstat.linnet'.)
License: GPL (>= 2)
URL: http://spatstat.org/
NeedsCompilation: yes
ByteCompile: true
BugReports: https://github.com/spatstat/spatstat.geom/issues
Packaged: 2023-01-30 00:51:47 UTC; adrian
Author: Adrian Baddeley [aut, cre, cph]
    (<https://orcid.org/0000-0001-9499-8382>),
  Rolf Turner [aut, cph] (<https://orcid.org/0000-0001-5521-5218>),
  Ege Rubak [aut, cph] (<https://orcid.org/0000-0002-6675-533X>),
  Tilman Davies [ctb],
  Ute Hahn [ctb],
  Abdollah Jalilian [ctb],
  Greg McSwiggan [ctb, cph],
  Sebastian Meyer [ctb, cph],
  Jens Oehlschlaegel [ctb, cph],
  Suman Rakshit [ctb],
  Dominic Schuhmacher [ctb],
  Rasmus Waagepetersen [ctb]
Repository: CRAN
Date/Publication: 2023-01-30 10:20:02 UTC
Built: R 4.2.2; x86_64-w64-mingw32; 2023-03-07 03:45:35 UTC; windows
ExperimentalWindowsRuntime: ucrt
Archs: x64
