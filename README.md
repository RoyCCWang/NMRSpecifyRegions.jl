# NMRSpecifyRegions

Given a set of frequency locations (e.g., the frequency locations from a Discrete Fourier transform of the NMR FID data) and a list of metabolites, this package filters the set such that the locations near a resonance component remains. This package relies on data from the GISSMO database for the resonance component location prediction.


# Install
From a Julia REPL or script, add the custom registry before adding the package:

```
using Pkg

Pkg.Registry.add(url = "https://github.com/RoyCCWang/RWPublicJuliaRegistry")

Pkt.add("NMRSpecifyRegions)
```

# Citation
Our work is undergoing peer review.

# License
This project is licensed under the Mozilla Public License v2.0; see the LICENSE file for details. Individual source files may contain the following tag instead of the full license text:
```
SPDX-License-Identifier: MPL-2.0
```

Using SPDX enables machine processing of license information based on the SPDX License Identifiers and makes it easier for developers to see at a glance which license they are dealing with.
